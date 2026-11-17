"""P8: LLM safety rationales for 24 targets, audited."""
import json
import sys
from pathlib import Path

import anthropic
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from tssc import RES, report  # noqa: E402

PRICE_IN, PRICE_OUT = 4.0, 20.0  # USD per million tokens


def sample(oof):
    o = oof.sort_values(["M3 M1 + curated events", "target"], ascending=[False, True])
    top = o.head(12)
    rest = o.iloc[12:]
    rnd = rest.sample(12, random_state=7)
    return pd.concat([top.assign(arm="top risk"), rnd.assign(arm="random")])


def main():
    t = pd.read_csv(RES / "targets.csv")
    drugs = pd.read_csv(RES / "drugs.csv")
    links = pd.read_csv(RES / "drug_targets.csv")
    warns = pd.read_csv(RES / "drug_warnings.csv")
    organ = pd.read_csv(RES / "organ_scores.csv").set_index("target")
    s = sample(pd.read_csv(RES / "oof_predictions.csv"))
    packer = report.Packer(t, drugs, links, warns, organ)

    out = RES / "llm_reports.jsonl"
    done = {}
    if out.exists():
        for line in out.read_text(encoding="utf-8").splitlines():
            j = json.loads(line)
            done[j["target"]] = j
    report.load_key()
    client = anthropic.Anthropic(max_retries=4)
    for _, r in s.iterrows():
        if r["target"] in done:
            continue
        rows = packer.pack(r["target"])
        try:
            parsed, use = report.ask(client, r["symbol"], rows)
        except anthropic.APIStatusError as e:
            print(f"{r['symbol']}: API error {e.status_code}", file=sys.stderr)
            continue
        rec = {"target": r["target"], "symbol": r["symbol"], "arm": r["arm"], "y": int(r["y"]),
               "pack": rows, "report": parsed, "usage": use}
        done[r["target"]] = rec
        with open(out, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec) + "\n")
        print(f"{r['symbol']:10s} {use['stop_reason']:9s} in {use['input_tokens']} out {use['output_tokens']}")

    aud, summ = [], []
    for rec in done.values():
        if rec["report"] is None:
            summ.append({"symbol": rec["symbol"], "arm": rec["arm"], "y": rec["y"], "risk_level": None,
                         "stop_reason": rec["usage"]["stop_reason"]})
            continue
        a = report.audit(rec["report"], [tuple(x) for x in rec["pack"]])
        aud.append(a.assign(symbol=rec["symbol"], arm=rec["arm"]))
        summ.append({"symbol": rec["symbol"], "arm": rec["arm"], "y": rec["y"],
                     "risk_level": rec["report"]["risk_level"], "sentences": len(a),
                     "stop_reason": rec["usage"]["stop_reason"],
                     "input_tokens": rec["usage"]["input_tokens"],
                     "output_tokens": rec["usage"]["output_tokens"]})
    aud = pd.concat(aud, ignore_index=True)
    aud.to_csv(RES / "p8_sentence_audit.csv", index=False)
    summ = pd.DataFrame(summ)
    summ.to_csv(RES / "p8_reports.csv", index=False)
    ok = aud["has_citation"] & aud["ids_valid"] & (aud["classes_unsupported"] == 0)
    fact = aud[aud["kind"] == "fact"]
    cost = (summ["input_tokens"].sum() * PRICE_IN + summ["output_tokens"].sum() * PRICE_OUT) / 1e6
    res = pd.DataFrame([{
        "reports": int(summ["risk_level"].notna().sum()), "sentences": len(aud),
        "fact": int((aud["kind"] == "fact").sum()), "hypothesis": int((aud["kind"] == "hypothesis").sum()),
        "gap": int((aud["kind"] == "gap").sum()),
        "share_sentences_pass": ok.mean(),
        "numbers_in_facts": int(fact["numbers"].sum()),
        "share_numbers_supported": 1 - fact["numbers_unsupported"].sum() / max(fact["numbers"].sum(), 1),
        "class_mentions": int(fact["classes"].sum()),
        "class_mentions_unsupported": int(fact["classes_unsupported"].sum()),
        "cost_usd": round(cost, 2)}])
    res.to_csv(RES / "p8_audit.csv", index=False)
    print(res.T.to_string())
    print(pd.crosstab(summ["risk_level"], summ["y"]).to_string())
    print(np.round(summ.groupby("arm")["y"].mean(), 2).to_string())


if __name__ == "__main__":
    main()
