import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tssc import data, organ, report, stats  # noqa: E402


def test_concordant_pairs_counts():
    labels = [1, 1, 0, 0, 0, 1]
    groups = ["a", "a", "a", "b", "b", "c"]
    same, total = stats.concordant_pairs(labels, groups)
    # a: 1 of 3 pairs agree; b: 1 of 1
    assert (same, total) == (2, 4)


def test_concordance_test_detects_shared_labels():
    groups = np.repeat(np.arange(40), 3)
    labels = np.repeat(np.arange(40) % 2, 3)
    share, null, p, _ = stats.concordance_test(labels, groups, perms=500)
    assert share == 1.0 and null < 0.6 and p < 0.01


def test_equal_weight_middle_rank_for_missing():
    t = pd.DataFrame({"a": [1.0, 2.0, np.nan], "b": [3.0, 1.0, 2.0]})
    s = data.equal_weight_score(t, ["a", "b"])
    assert np.isclose(s.iloc[2], (0.5 + 2 / 3) / 2)


def test_stratified_difference_finds_planted_shift():
    rng = np.random.default_rng(0)
    strata = np.repeat(["x", "y"], 60)
    member = np.tile(np.r_[np.ones(20, bool), np.zeros(40, bool)], 2)
    z = rng.normal(size=120) + member * 1.0
    obs, p, levels = stats.stratified_diff_test(z, member, strata, perms=500)
    assert obs > 0.5 and p < 0.01 and len(levels) == 2


def test_organ_matching_rule():
    assert organ.matches("venous blood", "blood", ["blood"])
    assert not organ.matches("tibial artery", "blood vessel", ["blood"])
    assert organ.matches("gastrocnemius medialis", "muscle tissue", ["muscle tissue"])
    assert organ.matches("coronary artery", "heart", ["artery", "aorta"])


def test_organ_z_takes_best_matching_tissue():
    g = pd.DataFrame({"target": "T", "tissue": ["liver", "lung", "skin", "brain"],
                      "parent": ["liver", "lung", "skin", "brain"], "tpm": [100.0, 1.0, 1.0, 1.0]})
    z = organ.organ_z(g, ["liver"])
    assert z["T"] > 1


def test_numbers_skip_gene_names():
    s = "CYP3A4 and IL-6 rose; LOEUF 0.34, z -0.52, 1,234 reports in 2004."
    assert report.numbers(s) == ["0.34", "-0.52", "1234", "2004"]


def test_supported_uses_rounding():
    assert report.supported("0.82", "factor -0.819 (-1 most constrained)") is False
    assert report.supported("-0.82", "factor -0.819 (-1 most constrained)")
    assert report.supported("4797", "reports 2907, llr 4797.3")
    assert not report.supported("12", "reports 13")


def test_audit_flags_bad_ids_numbers_and_classes():
    rows = [("G1", "LOEUF 0.34"), ("D1", "warnings: Withdrawn (hepatotoxicity, 2004)")]
    rep = {"risk_level": "high", "sentences": [
        {"text": "LOEUF is 0.34.", "kind": "fact", "evidence_ids": ["G1"]},
        {"text": "LOEUF is 0.5.", "kind": "fact", "evidence_ids": ["G1"]},
        {"text": "It was withdrawn for cardiotoxicity.", "kind": "fact", "evidence_ids": ["D1"]},
        {"text": "Maybe off-target.", "kind": "hypothesis", "evidence_ids": ["X9"]},
    ]}
    a = report.audit(rep, rows)
    assert a["numbers_unsupported"].tolist() == [0, 1, 0, 0]
    assert a["classes_unsupported"].tolist() == [0, 0, 1, 0]
    assert a["ids_valid"].tolist() == [True, True, True, False]


def test_parent_map_sends_salts_to_parent():
    mol = pd.DataFrame({"id": ["P", "S", "Q"], "parentId": [None, "P", None],
                        "childChemblIds": [["S"], None, None]})
    m = data.parent_map(mol)
    assert m["S"] == "P" and m["P"] == "P" and m["Q"] == "Q"


def test_oof_logistic_keeps_groups_out_of_training():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(100, 2))
    y = (x[:, 0] + rng.normal(size=100) > 0).astype(int)
    g = np.repeat(np.arange(10), 10)
    p = stats.oof_logistic(x, y, g)
    assert p.shape == (100,) and ((p > 0) & (p < 1)).all()
