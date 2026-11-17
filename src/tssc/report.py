"""Evidence packs, LLM safety rationales and their audit."""
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from . import ROOT, data

MODEL = "claude-opus-5-5"
TOX_CLASSES = ["carcinogenicity", "cardiotoxicity", "dermatological toxicity", "drug misuse",
               "gastrointestinal toxicity", "hematological toxicity", "hepatotoxicity",
               "immune system toxicity", "infectious disease", "metabolic toxicity",
               "musculoskeletal toxicity", "nephrotoxicity", "neurotoxicity", "ocular toxicity",
               "occular toxicity", "ototoxicity", "psychiatric toxicity", "respiratory toxicity",
               "teratogenicity", "vascular toxicity"]

ORGAN_NAME = {"hepatotoxicity": "liver", "cardiotoxicity": "heart", "nephrotoxicity": "kidney",
              "neurotoxicity": "brain", "psychiatric toxicity": "brain",
              "dermatological toxicity": "skin", "gastrointestinal toxicity": "gut",
              "respiratory toxicity": "lung", "hematological toxicity": "blood",
              "musculoskeletal toxicity": "skeletal muscle", "vascular toxicity": "artery"}

SYSTEM = """You are a target safety scientist writing a short safety assessment of one drug target for an internal review.
You receive an evidence pack: rows with IDs. Use only these rows.
Rules:
- Write 6 to 10 sentences. Each sentence must list the IDs of the rows it rests on.
- Mark a sentence "fact" if it restates what the cited rows say, "hypothesis" if it proposes a mechanism (on-target or off-target) that the rows suggest but do not state, and "gap" if it names missing or weak evidence; a gap sentence cites the row that shows the absence.
- Use only numbers that appear in the cited rows, written as they appear there. Do not compute new numbers or counts.
- Name a toxicity class only if a cited row names it.
- Give an overall risk level (low, medium or high) for developing a new drug against this target."""

SCHEMA = {
    "type": "object",
    "properties": {
        "risk_level": {"type": "string", "enum": ["low", "medium", "high"]},
        "sentences": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "text": {"type": "string"},
                    "kind": {"type": "string", "enum": ["fact", "hypothesis", "gap"]},
                    "evidence_ids": {"type": "array", "items": {"type": "string"}},
                },
                "required": ["text", "kind", "evidence_ids"],
                "additionalProperties": False,
            },
        },
    },
    "required": ["risk_level", "sentences"],
    "additionalProperties": False,
}


def _fmt(x, d=3):
    return "missing" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{d}f}"


class Packer:
    """Builds the evidence pack of one target from the raw tables."""

    def __init__(self, targets, drugs, links, warns, organ=None):
        self.t = targets.set_index("target")
        self.drugs = drugs.set_index("drug")
        self.links, self.warns = links, warns
        self.organ = organ
        tg = data.read("target", ["id", "constraint", "pathways"]).set_index("id")
        self.tg = tg[tg.index.isin(self.t.index)]
        self.pri = data.read("target_prioritisation").set_index("targetId")
        mp = data.read("mouse_phenotype", ["targetFromSourceId", "modelPhenotypeLabel",
                                           "modelPhenotypeClasses"])
        self.mp = mp[mp["targetFromSourceId"].isin(self.t.index)]
        ev = data.read("target_safety_event")
        self.ev = ev[ev["targetId"].isin(self.t.index)]
        mol = data.read("drug_molecule", ["id", "parentId", "childChemblIds", "maximumClinicalStage",
                                          "drugType", "name"])
        pmap = data.parent_map(mol)
        self.adr = data.faers_counts(data.read("openfda_significant_adverse_drug_reactions"), pmap)

    def pack(self, tid):
        rows = []
        r = self.t.loc[tid]
        p = self.pri.loc[tid] if tid in self.pri.index else pd.Series(dtype=float)
        loeuf = None
        cons = self.tg.loc[tid, "constraint"] if tid in self.tg.index else None
        if cons is not None:
            for c in cons:
                if c["constraintType"] == "lof":
                    loeuf = c["oeUpper"]
        rows.append(("G1", f"Genetic constraint: gnomAD loss-of-function LOEUF {_fmt(loeuf, 2)}; "
                           f"Open Targets constraint factor {_fmt(p.get('geneticConstraint'))} "
                           "(-1 most constrained, 1 least)"))
        mp = self.mp[self.mp["targetFromSourceId"] == tid]
        if len(mp):
            cl = pd.Series([c["label"] for cs in mp["modelPhenotypeClasses"] for c in cs]).value_counts()
            top = "; ".join(f"{k} ({v})" for k, v in cl.head(6).items())
            rows.append(("G2", f"Mouse knockouts: {len(mp)} phenotypes; most frequent classes: {top}; "
                               f"Open Targets mouse KO score {_fmt(p.get('mouseKOScore'))} (-1 most severe, 0 none)"))
        else:
            rows.append(("G2", "Mouse knockouts: no phenotype data in Open Targets"))
        rows.append(("G3", f"Tissue expression: specificity factor {_fmt(p.get('tissueSpecificity'))}, "
                           f"distribution factor {_fmt(p.get('tissueDistribution'))} "
                           "(-1 least favourable, 1 most favourable)"))
        rows.append(("G4", "Cancer driver gene: " + ("yes" if r["risk_isCancerDriverGene"] > 0 else "no")))
        rows.append(("G5", f"Closest paralog factor {_fmt(p.get('paralogMaxIdentityPercentage'))} "
                           "(-1 identical paralog, 0 none close)"))
        ev = self.ev[self.ev["targetId"] == tid]
        if len(ev):
            items = {}
            for _, e in ev.iterrows():
                eff = ", ".join(sorted({x["direction"] for x in (e["effects"] if e["effects"] is not None else [])}))
                tis = ", ".join(sorted({b["tissueLabel"] for b in (e["biosamples"] if e["biosamples"] is not None else [])
                                        if b["tissueLabel"]}))
                key = (e["event"], e["datasource"], eff, tis)
                items[key] = items.get(key, 0) + 1
            for k, ((evn, src, eff, tis), n) in enumerate(list(items.items())[:10], 1):
                rows.append((f"S{k}", f"Curated safety event: {evn}; source {src}"
                                      + (f"; effect {eff}" if eff else "") + (f"; tissues {tis}" if tis else "")
                                      + (f"; listed {n} times" if n > 1 else "")))
        else:
            rows.append(("S0", "Curated safety events: none recorded in Open Targets"))
        pw = self.tg.loc[tid, "pathways"] if tid in self.tg.index else None
        if pw is not None and len(pw):
            top = sorted({x["topLevelTerm"] for x in pw if x["topLevelTerm"]})
            rows.append(("P1", "Reactome top-level pathways: " + "; ".join(top[:10])))
            rows.append(("P2", "Reactome pathways (first 6): " + "; ".join(x["pathway"] for x in pw[:6])))
        else:
            rows.append(("P1", "Reactome pathways: none annotated"))
        if self.organ is not None and tid in self.organ.index:
            z = self.organ.loc[tid].dropna()
            z = z.rename(ORGAN_NAME).groupby(level=0).first()
            rows.append(("X1", "GTEx expression of this target, z-score across tissues (0 = its average "
                               "tissue): " + "; ".join(f"{k} {v:.2f}" for k, v in z.items())))
        ds = self.links[self.links["target"] == tid]["drug"].tolist()
        ds = sorted(ds, key=lambda d: (not self.drugs.loc[d, "warned"], str(self.drugs.loc[d, "name"])))
        for k, d in enumerate(ds[:6], 1):
            w = self.warns[self.warns["drug"] == d]
            wt = "; ".join(sorted({f"{a} ({b if isinstance(b, str) else 'class not given'}, {c}, {int(y) if y == y else 'year not given'})"
                                   for a, b, c, y in zip(w["warningType"], w["toxicityClass"], w["country"], w["year"])}))
            # salts repeat the parent's terms
            a = (self.adr[self.adr["drug"] == d].sort_values("llr", ascending=False)
                 .drop_duplicates("event").head(5))
            ft = "; ".join(f"{e} (reports {c}, LLR {l:.1f})" for e, c, l in zip(a["event"], a["count"], a["llr"]))
            act = self.links[(self.links["target"] == tid) & (self.links["drug"] == d)]["actionType"].iloc[0]
            rows.append((f"D{k}", f"Approved drug {self.drugs.loc[d, 'name']} ({self.drugs.loc[d, 'drugType']}, "
                                  f"{act.lower()}): warnings: {wt or 'none'}; top FAERS signals: {ft or 'none'}"))
        if len(ds) > 6:
            rows.append(("D0", f"{len(ds)} approved drugs act on this target; 6 are listed"))
        return rows


def prompt(symbol, rows):
    pack = "\n".join(f"[{i}] {txt}" for i, txt in rows)
    return f"Target: {symbol}\n\nEvidence pack:\n{pack}"


def load_key(path=None):
    if "ANTHROPIC_API_KEY" in os.environ:
        return
    p = Path(path or os.environ.get("TSSC_ENV_FILE", ROOT / ".env"))
    for line in p.read_text(encoding="utf-8-sig").splitlines():
        k, _, v = line.partition("=")
        if k.strip() == "ANTHROPIC_API_KEY" and v.strip():
            os.environ["ANTHROPIC_API_KEY"] = v.strip().strip('"')


def ask(client, symbol, rows):
    """One structured call; returns parsed output and usage."""
    r = client.beta.messages.create(
        model=MODEL, max_tokens=16000, system=SYSTEM,
        betas=["server-side-fallback-2026-07-01"], fallbacks="default",
        output_config={"effort": "medium", "format": {"type": "json_schema", "schema": SCHEMA}},
        messages=[{"role": "user", "content": prompt(symbol, rows)}],
    )
    use = {"input_tokens": r.usage.input_tokens, "output_tokens": r.usage.output_tokens,
           "stop_reason": r.stop_reason, "model": r.model}
    if r.stop_reason != "end_turn":
        return None, use
    text = next(b.text for b in r.content if b.type == "text")
    return json.loads(text), use


# skips digits inside names such as CYP3A4 or IL-6
NUM = re.compile(r"(?<![\w\-.])-?\d+(?:,\d{3})*(?:\.\d+)?(?![\w\-])")


def numbers(s):
    s = s.replace("−", "-")
    return [m.replace(",", "") for m in NUM.findall(s)]


def supported(num, cited_text):
    """A number is supported if a cited number rounds to it."""
    d = len(num.split(".")[1]) if "." in num else 0
    x = float(num)
    return any(round(float(c), d) == x for c in numbers(cited_text))


def audit(report, rows):
    ids = {i: txt for i, txt in rows}
    out = []
    for k, s in enumerate(report["sentences"]):
        cited = [i for i in s["evidence_ids"] if i in ids]
        bad = [i for i in s["evidence_ids"] if i not in ids]
        ctext = " ".join(ids[i] for i in cited).lower()
        nums = numbers(s["text"]) if s["kind"] == "fact" else []
        unsup = [n for n in nums if not supported(n, ctext)]
        classes = [c for c in TOX_CLASSES if c in s["text"].lower()] if s["kind"] == "fact" else []
        cls_bad = [c for c in classes if c not in ctext]
        out.append({"sentence": k, "kind": s["kind"], "n_ids": len(s["evidence_ids"]),
                    "has_citation": len(s["evidence_ids"]) > 0, "ids_valid": not bad,
                    "numbers": len(nums), "numbers_unsupported": len(unsup),
                    "unsupported": "; ".join(unsup), "classes": len(classes),
                    "classes_unsupported": len(cls_bad), "text": s["text"]})
    return pd.DataFrame(out)
