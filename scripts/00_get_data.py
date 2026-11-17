"""Download the Open Targets Platform tables used here (26.09, about 37 MB)."""
import sys
from pathlib import Path

import requests

RELEASE = "26.09"
BASE = f"https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{RELEASE}/output/"
TABLES = ["drug_molecule", "drug_mechanism_of_action", "drug_warning", "target",
          "target_safety_event", "mouse_phenotype", "openfda_significant_adverse_drug_reactions",
          "reactome", "target_prioritisation"]
OUT = Path(__file__).resolve().parents[1] / "data" / "raw" / RELEASE


def files(table):
    import re
    r = requests.get(BASE + table + "/", timeout=60)
    r.raise_for_status()
    return sorted(set(re.findall(r'href="([^"/]+\.parquet)"', r.text)))


def main():
    for t in TABLES:
        d = OUT / t
        d.mkdir(parents=True, exist_ok=True)
        for f in files(t):
            p = d / f
            if p.exists():
                continue
            with requests.get(BASE + t + "/" + f, timeout=300, stream=True) as r:
                r.raise_for_status()
                tmp = p.with_suffix(".part")
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
                tmp.rename(p)
        print(t, sum(x.stat().st_size for x in d.glob("*.parquet")) / 1e6, "MB", file=sys.stderr)


if __name__ == "__main__":
    main()
