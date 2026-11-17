"""Approved drugs, their targets and warnings, and target evidence."""
import numpy as np
import pandas as pd

from . import (CLINICAL_SOURCES, DROP_ACTIONS, EXPERIMENTAL_SOURCES, KEEP_DRUG_TYPES, PRE_DRUG,
               RAW, WARNING_TYPES)


def read(table, columns=None):
    return pd.read_parquet(RAW / table, columns=columns)


def parent_map(mol):
    """Salt or child ChEMBL id -> parent id."""
    par = mol["parentId"].where(mol["parentId"].notna(), mol["id"])
    m = dict(zip(mol["id"], par))
    for kids, p in zip(mol["childChemblIds"], mol["id"]):
        if kids is not None:
            for k in kids:
                m.setdefault(k, m.get(p, p))
    return m


def approved_drugs(mol, pmap):
    """Approved parent molecules of kept types."""
    a = mol[mol["maximumClinicalStage"] == "APPROVAL"].copy()
    a["drug"] = a["id"].map(pmap)
    # type from the parent row, else a child
    types = mol.assign(drug=mol["id"].map(pmap)).sort_values(
        "parentId", na_position="first").drop_duplicates("drug").set_index("drug")
    out = pd.DataFrame({"drug": sorted(set(a["drug"]))})
    out["drugType"] = out["drug"].map(types["drugType"])
    out["name"] = out["drug"].map(types["name"])
    return out[out["drugType"].isin(KEEP_DRUG_TYPES)].reset_index(drop=True)


def drug_targets(moa, pmap, drugs, protein_coding):
    """One row per (drug, target) from qualifying mechanisms."""
    m = moa[(moa["targetType"] == "single protein") & ~moa["actionType"].isin(DROP_ACTIONS)]
    m = m.explode("chemblIds").explode("targets").dropna(subset=["chemblIds", "targets"])
    m = m.assign(drug=m["chemblIds"].map(lambda c: pmap.get(c, c)), target=m["targets"])
    m = m[m["drug"].isin(set(drugs["drug"])) & m["target"].isin(protein_coding)]
    return (m.groupby(["drug", "target"])
             .agg(actionType=("actionType", "first"), mechanism=("mechanismOfAction", "first"))
             .reset_index())


def drug_warnings(warn, pmap):
    w = warn[warn["warningType"].isin(WARNING_TYPES)].explode("chemblIds").dropna(subset=["chemblIds"])
    w = w.assign(drug=w["chemblIds"].map(lambda c: pmap.get(c, c)))
    cols = ["drug", "warningType", "toxicityClass", "country", "year", "description", "efoTerm"]
    return w[cols].drop_duplicates().reset_index(drop=True)


def faers_counts(adr, pmap):
    a = adr[adr["llr"] > adr["critval"]]
    a = a.assign(drug=a["chembl_id"].map(lambda c: pmap.get(c, c)))
    return a


def level1_class(tc):
    if tc is None or len(tc) == 0:
        return "unclassified"
    l1 = [c["label"] for c in tc if c["level"] == "l1"]
    return l1[0] if l1 else "unclassified"


def build():
    """Cohort tables: drugs, drug-target links, warnings, target features."""
    mol = read("drug_molecule", ["id", "parentId", "childChemblIds", "maximumClinicalStage",
                                 "drugType", "name"])
    pmap = parent_map(mol)
    tgt = read("target", ["id", "approvedSymbol", "biotype", "targetClass", "pathways"])
    pc = set(tgt.loc[tgt["biotype"] == "protein_coding", "id"])
    drugs = approved_drugs(mol, pmap)
    links = drug_targets(read("drug_mechanism_of_action"), pmap, drugs, pc)
    drugs = drugs[drugs["drug"].isin(set(links["drug"]))].reset_index(drop=True)
    warns = drug_warnings(read("drug_warning"), pmap)
    warns = warns[warns["drug"].isin(set(drugs["drug"]))].reset_index(drop=True)
    drugs["warned"] = drugs["drug"].isin(set(warns["drug"]))
    ntarg = links.groupby("drug")["target"].nunique()
    drugs["n_targets"] = drugs["drug"].map(ntarg)

    faers = faers_counts(read("openfda_significant_adverse_drug_reactions"), pmap)
    drugs["faers_terms"] = drugs["drug"].map(faers.groupby("drug")["event"].nunique()).fillna(0)

    lk = links.merge(drugs[["drug", "warned", "n_targets", "faers_terms"]], on="drug")
    t = lk.groupby("target").agg(n_drugs=("drug", "nunique"), y=("warned", "max"),
                                 faers=("faers_terms", "max")).reset_index()
    single = lk[lk["n_targets"] == 1].groupby("target").agg(
        n_single=("drug", "nunique"), y_single=("warned", "max"))
    t = t.join(single, on="target")
    t["n_single"] = t["n_single"].fillna(0).astype(int)
    t["y"] = t["y"].astype(int)
    t["log_n"] = np.log(t["n_drugs"])
    t["faers_score"] = np.log1p(t["faers"])

    info = tgt.set_index("id")
    t["symbol"] = t["target"].map(info["approvedSymbol"])
    t["class_l1"] = t["target"].map(info["targetClass"]).map(level1_class)

    pri = read("target_prioritisation").set_index("targetId")
    for f in PRE_DRUG:
        v = t["target"].map(pri[f])
        if f == "isCancerDriverGene":
            v = v.fillna(0)
        t["risk_" + f] = -v

    ev = read("target_safety_event", ["targetId", "datasource"])
    clin = set(ev.loc[ev["datasource"].isin(CLINICAL_SOURCES), "targetId"])
    expm = set(ev.loc[ev["datasource"].isin(EXPERIMENTAL_SOURCES), "targetId"])
    t["event_clinical"] = t["target"].isin(clin).astype(int)
    t["event_experimental"] = t["target"].isin(expm).astype(int)
    return drugs, links, warns, t.sort_values("target").reset_index(drop=True)


def equal_weight_score(t, cols):
    """Mean percentile rank; missing gets the middle rank."""
    r = t[cols].rank(pct=True, method="average")
    return r.fillna(0.5).mean(axis=1)
