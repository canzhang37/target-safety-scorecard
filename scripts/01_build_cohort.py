"""Cohort of drug targets, outcome and evidence."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from tssc import RES, data  # noqa: E402


def main():
    RES.mkdir(exist_ok=True)
    drugs, links, warns, t = data.build()
    drugs.to_csv(RES / "drugs.csv", index=False)
    links.to_csv(RES / "drug_targets.csv", index=False)
    warns.to_csv(RES / "drug_warnings.csv", index=False)
    t.to_csv(RES / "targets.csv", index=False)
    s = {
        "approved_drugs": len(drugs),
        "drugs_with_warning": int(drugs["warned"].sum()),
        "single_target_drugs": int((drugs["n_targets"] == 1).sum()),
        "targets": len(t),
        "targets_with_warning": int(t["y"].sum()),
        "targets_single_cohort": int((t["n_single"] > 0).sum()),
        "warning_rows": len(warns),
    }
    for k, v in s.items():
        print(f"{k:24s} {v}")
    print(drugs["drugType"].value_counts().to_string())
    print(warns["warningType"].value_counts().to_string())


if __name__ == "__main__":
    main()
