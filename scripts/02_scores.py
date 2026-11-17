"""P1-P6: shared warnings, scores, weighting, circularity."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from tssc import PRE_DRUG, RES, data, stats  # noqa: E402

RISK = ["risk_" + f for f in PRE_DRUG]


def models(t, y, n_col):
    """Out-of-fold scores for M0-M3 plus raw scores."""
    g = t["class_l1"].values
    base = np.log(t[[n_col]].values)
    eq = data.equal_weight_score(t, RISK).values[:, None]
    flags = stats.impute_with_flags(t, RISK).values
    cur = t[["event_clinical", "event_experimental"]].values
    out = {
        "M0 drug count": stats.oof_logistic(base, y, g),
        "M1 count + equal-weight evidence": stats.oof_logistic(np.hstack([base, eq]), y, g),
        "M2 count + learned weights": stats.oof_logistic(np.hstack([base, flags]), y, g),
        "M3 M1 + curated events": stats.oof_logistic(np.hstack([base, eq, cur]), y, g),
    }
    return out, eq.ravel()


def main():
    t = pd.read_csv(RES / "targets.csv")
    links = pd.read_csv(RES / "drug_targets.csv")
    drugs = pd.read_csv(RES / "drugs.csv")
    print(f"targets {len(t)}, with warning {t['y'].sum()}, classes {t['class_l1'].nunique()}")

    # P1 within-target concordance, single-target drugs
    s = links.merge(drugs[["drug", "warned", "n_targets"]], on="drug")
    s = s[s["n_targets"] == 1]
    k = s.groupby("target")["drug"].transform("size")
    s = s[k >= 2].sort_values(["target", "drug"])
    share, null, p, pairs = stats.concordance_test(s["warned"].astype(int), s["target"])
    p1 = pd.DataFrame([{"drugs": len(s), "targets": s["target"].nunique(), "pairs": pairs,
                        "concordant_share": share, "null_share": null, "p_one_sided": p}])
    p1.to_csv(RES / "p1_concordance.csv", index=False)
    print(p1.round(4).to_string(index=False))

    y = t["y"].values
    oof, eq = models(t, y, "n_drugs")
    scores = dict(oof)
    scores["FAERS terms (raw)"] = t["faers_score"].values
    scores["drug count (raw)"] = t["log_n"].values
    scores["equal-weight evidence (raw)"] = eq
    scores["curated event, clinical sources"] = t["event_clinical"].values
    scores["curated event, experimental sources"] = t["event_experimental"].values
    for c in RISK:
        scores[c] = t[c].fillna(t[c].median()).values
    point, draws = stats.boot_auc(y, scores, reps=1000, seed=1)
    rows = [{"score": k, "auroc": point[k], "low": stats.ci(draws[k])[0], "high": stats.ci(draws[k])[1]}
            for k in scores]
    auc = pd.DataFrame(rows)
    auc.to_csv(RES / "auroc.csv", index=False)
    print(auc.round(3).to_string(index=False))

    def diff(a, b):
        d = draws[a] - draws[b]
        lo, hi = stats.ci(d)
        return {"comparison": f"{a} minus {b}", "diff": point[a] - point[b], "low": lo, "high": hi}

    dd = pd.DataFrame([
        diff("M1 count + equal-weight evidence", "M0 drug count"),
        diff("M2 count + learned weights", "M1 count + equal-weight evidence"),
        diff("M3 M1 + curated events", "M1 count + equal-weight evidence"),
        diff("curated event, clinical sources", "curated event, experimental sources"),
    ])
    dd.to_csv(RES / "auroc_differences.csv", index=False)
    print(dd.round(3).to_string(index=False))
    pd.DataFrame({"target": t["target"], "symbol": t["symbol"], "y": y, **oof,
                  "equal_weight": eq}).to_csv(RES / "oof_predictions.csv", index=False)

    # single-target sensitivity cohort
    ts = t[t["n_single"] > 0].reset_index(drop=True)
    ys = ts["y_single"].astype(int).values
    oofs, _ = models(ts, ys, "n_single")
    sc = {"drug count (raw)": np.log(ts["n_single"].values), **oofs}
    ps, ds = stats.boot_auc(ys, sc, reps=1000, seed=2)
    d = ds["M1 count + equal-weight evidence"] - ds["M0 drug count"]
    sens = pd.DataFrame([{"targets": len(ts), "with_warning": int(ys.sum()),
                          "auroc_count": ps["drug count (raw)"],
                          "auroc_M0": ps["M0 drug count"],
                          "auroc_M1": ps["M1 count + equal-weight evidence"],
                          "M1_minus_M0": ps["M1 count + equal-weight evidence"] - ps["M0 drug count"],
                          "low": stats.ci(d)[0], "high": stats.ci(d)[1]}])
    sens.to_csv(RES / "sensitivity_single_target.csv", index=False)
    print(sens.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
