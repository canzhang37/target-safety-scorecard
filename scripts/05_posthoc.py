"""Post hoc checks, not in the plan."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from tssc import RES, data, report, stats  # noqa: E402


def faers_over_count(t):
    g = t["class_l1"].values
    y = t["y"].values
    base = t[["log_n"]].values
    m0 = stats.oof_logistic(base, y, g)
    mf = stats.oof_logistic(np.hstack([base, t[["faers_score"]].values]), y, g)
    pt, dr = stats.boot_auc(y, {"M0": m0, "M0 + FAERS": mf}, reps=1000, seed=3)
    d = dr["M0 + FAERS"] - dr["M0"]
    return pd.DataFrame([{"auroc_M0": pt["M0"], "auroc_M0_faers": pt["M0 + FAERS"],
                          "diff": pt["M0 + FAERS"] - pt["M0"],
                          "low": stats.ci(d)[0], "high": stats.ci(d)[1]}])


def shared_profiles():
    """FAERS profiles attached to more than one parent molecule."""
    mol = data.read("drug_molecule", ["id", "parentId", "childChemblIds", "maximumClinicalStage",
                                      "drugType", "name"])
    pmap = data.parent_map(mol)
    adr = data.faers_counts(data.read("openfda_significant_adverse_drug_reactions"), pmap)
    prof = (adr.groupby("chembl_id")[["event", "count"]]
               .apply(lambda g: tuple(sorted(zip(g["event"], g["count"])))))
    df = pd.DataFrame({"chembl": prof.index, "parent": [pmap.get(c, c) for c in prof.index],
                       "profile": prof.values})
    grp = df.groupby("profile")["parent"].agg(lambda s: sorted(set(s)))
    multi = grp[grp.map(len) > 1]
    names = mol.set_index("id")["name"]
    cohort = set(pd.read_csv(RES / "drugs.csv")["drug"])
    rows = [{"parents": "; ".join(p), "names": "; ".join(str(names.get(x, x)) for x in p),
             "terms": len(prof_), "cohort_drugs": len(set(p) & cohort)}
            for prof_, p in multi.items()]
    return pd.DataFrame(rows), df["parent"].nunique()


def numbers_outside_facts():
    tot = bad = 0
    for line in (RES / "llm_reports.jsonl").read_text(encoding="utf-8").splitlines():
        r = json.loads(line)
        ids = {i: t for i, t in r["pack"]}
        for s in r["report"]["sentences"]:
            if s["kind"] == "fact":
                continue
            ct = " ".join(ids.get(i, "") for i in s["evidence_ids"]).lower()
            for n in report.numbers(s["text"]):
                tot += 1
                bad += not report.supported(n, ct)
    return tot, bad


def main():
    t = pd.read_csv(RES / "targets.csv")
    f = faers_over_count(t)
    f.to_csv(RES / "posthoc_faers_over_count.csv", index=False)
    print(f.round(3).to_string(index=False))
    sp, n_par = shared_profiles()
    sp.to_csv(RES / "posthoc_faers_shared_profiles.csv", index=False)
    print(f"{len(sp)} profiles shared by {sp['parents'].str.count(';').add(1).sum()} of {n_par} "
          f"parent molecules; cohort drugs involved {sp['cohort_drugs'].sum()}")
    tot, bad = numbers_outside_facts()
    pd.DataFrame([{"numbers": tot, "unsupported": bad}]).to_csv(
        RES / "posthoc_numbers_non_fact.csv", index=False)
    print(f"numbers in hypothesis and gap sentences: {tot}, unsupported {bad}")


if __name__ == "__main__":
    main()
