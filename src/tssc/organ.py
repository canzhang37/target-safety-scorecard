"""GTEx expression from the Open Targets API and organ scores."""
import time

import numpy as np
import pandas as pd
import requests

API = "https://api.platform.opentargets.org/api/v4/graphql"
QUERY = """query($id: String!) { target(ensemblId: $id) {
  baselineExpression(page: {index: 0, size: 3000}) { count rows {
    datasourceId unit median
    tissueBiosample { biosampleName } tissueBiosampleParent { biosampleName } } } } }"""


def fetch_gtex(target_id, session=None, tries=4):
    s = session or requests
    for k in range(tries):
        try:
            r = s.post(API, json={"query": QUERY, "variables": {"id": target_id}}, timeout=120)
            r.raise_for_status()
            be = r.json()["data"]["target"]["baselineExpression"]
            break
        except (requests.RequestException, KeyError, TypeError):
            if k == tries - 1:
                raise
            time.sleep(2 ** k)
    if be["count"] > len(be["rows"]):
        raise RuntimeError(f"{target_id}: page too small ({be['count']} rows)")
    rows = [{"target": target_id, "tissue": x["tissueBiosample"]["biosampleName"],
             "parent": (x["tissueBiosampleParent"] or {}).get("biosampleName"),
             "tpm": x["median"]} for x in be["rows"] if x["datasourceId"] == "gtex"]
    return pd.DataFrame(rows)


def matches(tissue, parent, keywords):
    t, p = (tissue or "").lower(), (parent or "").lower()
    return any(k == p or k in t for k in keywords)


def organ_z(gtex, keywords):
    """Highest z-score among matching tissues, per target."""
    g = gtex.copy()
    g["v"] = np.log1p(g["tpm"])
    st = g.groupby("target")["v"].agg(["mean", "std"])
    g = g.join(st, on="target")
    g["z"] = (g["v"] - g["mean"]) / g["std"].replace(0, np.nan)
    hit = [matches(a, b, keywords) for a, b in zip(g["tissue"], g["parent"])]
    return g[hit].groupby("target")["z"].max()
