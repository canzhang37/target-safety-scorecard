"""P7: are warnings in organs where the target is expressed?"""
import sys
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from tssc import ORGANS, RES, organ, stats  # noqa: E402


def main():
    links = pd.read_csv(RES / "drug_targets.csv")
    warns = pd.read_csv(RES / "drug_warnings.csv")
    t = pd.read_csv(RES / "targets.csv")
    w = warns[warns["toxicityClass"].isin(ORGANS)].merge(links[["drug", "target"]], on="drug")
    classes = w.groupby("target")["toxicityClass"].agg(lambda s: sorted(set(s)))
    ids = sorted(classes.index)
    print(f"targets with an organ-mapped warning: {len(ids)}")

    raw = RES / "raw"
    raw.mkdir(exist_ok=True)
    cache = raw / "gtex_rows.parquet"
    have = pd.read_parquet(cache) if cache.exists() else pd.DataFrame(columns=["target"])
    todo = [i for i in ids if i not in set(have["target"])]
    with requests.Session() as s:
        parts = [organ.fetch_gtex(i, s) for i in todo]
    gtex = pd.concat([have] + parts, ignore_index=True) if parts else have
    gtex.to_parquet(cache, index=False)
    gtex = gtex[gtex["target"].isin(ids)]
    print(f"GTEx rows {len(gtex)}, tissues {gtex['tissue'].nunique()}")

    z = pd.DataFrame({c: organ.organ_z(gtex, kw) for c, kw in ORGANS.items()}).reindex(ids)
    z.index.name = "target"
    z.round(4).to_csv(RES / "organ_scores.csv")

    long = []
    for c in ORGANS:
        for i in ids:
            if pd.notna(z.loc[i, c]):
                long.append({"target": i, "toxicityClass": c, "member": c in classes[i], "z": z.loc[i, c]})
    long = pd.DataFrame(long)
    obs, p, levels = stats.stratified_diff_test(long["z"], long["member"], long["toxicityClass"])
    per = (long.groupby(["toxicityClass", "member"])["z"].agg(["size", "mean"]).unstack()
               .set_axis(["n_other", "n_class", "z_other", "z_class"], axis=1)
               .reset_index())
    per["difference"] = per["z_class"] - per["z_other"]
    per.to_csv(RES / "p7_organ_by_class.csv", index=False)
    print(per.round(3).to_string(index=False))
    sym = t.set_index("target")["symbol"]
    out = pd.DataFrame([{"classes_tested": len(levels), "weighted_mean_difference": obs, "p_one_sided": p}])
    out.to_csv(RES / "p7_organ_match.csv", index=False)
    print(out.round(4).to_string(index=False))
    long.assign(symbol=long["target"].map(sym)).to_csv(RES / "p7_organ_long.csv", index=False)


if __name__ == "__main__":
    main()
