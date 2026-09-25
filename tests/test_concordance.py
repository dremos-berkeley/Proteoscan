import json
import os
import re

import pandas as pd
import pytest

from conftest import SAMPLES, ref_dir, require_ref
from peaksparity import concordance as cc
from peaksparity.peaks_export import read_export, write_export

LEVELS_ALWAYS = ["psm", "scan", "peptide", "peptide_sequence", "peptide_IL_collapsed",
                 "protein_group", "protein_top_accession", "protein_accession"]


def _assert_perfect(m):
    for name, lv in m["levels"].items():
        assert lv["recall"] == 1.0, name
        assert lv["precision"] == 1.0, name
        assert lv["jaccard"] == 1.0, name
        assert lv["ref_only"] == 0 and lv["test_only"] == 0, name
        assert lv["shared"] == lv["n_ref"] == lv["n_test"], name


@pytest.mark.parametrize("k", SAMPLES)
def test_self_concordance_is_perfect(k, tmp_path):
    d = require_ref(k)
    res = cc.compare(d, d, features=True)
    m = res["metrics"]
    assert set(LEVELS_ALWAYS + ["feature"]) <= set(m["levels"])
    _assert_perfect(m)
    assert m["headline"]["psm_recall"] == 1.0 and m["headline"]["psm_precision"] == 1.0
    assert m["psm_same_scan_different_peptide"] == {"ref_only": 0, "test_only": 0}
    assert not m["warnings"]
    # PSM keys de-duplicate the per-feature repeats PEAKS writes
    assert m["levels"]["psm"]["n_ref_rows"] >= m["levels"]["psm"]["n_ref"]
    # PEAKS' own "#Peptides" in the summary is the I/L-collapsed count
    stats = m["ref_statistics"][f"Sample {k}"]
    assert m["levels"]["peptide_IL_collapsed"]["n_ref"] == stats["peptides"]
    assert m["levels"]["protein_group"]["n_ref"] == stats["protein_groups"]
    assert m["levels"]["protein_accession"]["n_ref"] == stats["proteins_all"]
    assert m["levels"]["protein_top_accession"]["n_ref"] == stats["proteins_top"]
    assert m["levels"]["scan"]["n_ref"] == stats["identified_scans"]
    cc.write_report(res, tmp_path)
    for f in ["metrics.json", "summary.md", "psm_ref_only.csv", "psm_test_only.csv",
              "peptide_ref_only.csv", "peptide_test_only.csv", "protein_group_ref_only.csv",
              "protein_group_test_only.csv", "feature_ref_only.csv", "feature_test_only.csv"]:
        assert (tmp_path / f).exists(), f
    assert "100.00%" in (tmp_path / "summary.md").read_text()


def test_cli_self(tmp_path, capsys):
    d = require_ref(1)
    rc = cc.main(["--ref", d, "--test", d, "--out", str(tmp_path / "out"), "--features"])
    assert rc == 0
    m = json.loads((tmp_path / "out" / "metrics.json").read_text())
    assert m["headline"]["psm_recall"] == 1.0
    assert "PSM recall 100.00%" in capsys.readouterr().out


@pytest.mark.parametrize("a,b", [(1, 2), (3, 4)])
def test_different_samples_low_overlap(a, b, tmp_path):
    da, db = require_ref(a), require_ref(b)
    res = cc.compare(da, db, features=True)
    m = res["metrics"]
    # different raw files -> no PSM can match
    assert m["levels"]["psm"]["shared"] == 0
    assert m["headline"]["psm_recall"] == 0.0
    assert m["levels"]["feature"]["shared"] == 0
    # ignoring the file name, scan numbers still almost never coincide
    res2 = cc.compare(da, db, features=True, ignore_source=True)
    m2 = res2["metrics"]
    assert m2["levels"]["psm"]["jaccard"] < 0.05
    assert m2["levels"]["feature"]["jaccard"] < 0.10
    assert m2["levels"]["peptide"]["jaccard"] < 0.5
    assert m2["psm_ref_only_breakdown"]["same_scan_same_peptide_other_charge"] >= 0
    assert sum(m2["psm_ref_only_breakdown"].values()) == m2["levels"]["psm"]["ref_only"]
    assert sum(m2["psm_test_only_breakdown"].values()) == m2["levels"]["psm"]["test_only"]
    cc.write_report(res2, tmp_path)
    assert (tmp_path / "summary.md").read_text().startswith("# Concordance report")


def _to_proforma(s):
    names = {"+15.99": "Oxidation", "+0.98": "Deamidated", "+57.02": "Carbamidomethyl"}
    return re.sub(r"\(([+-][0-9.]+)\)", lambda m: f"[{names[m.group(1)]}]", s)


def test_proforma_test_side_is_equivalent(tmp_path):
    d = require_ref(3)
    exp = read_export(d)
    exp.psms["Peptide"] = exp.psms["Peptide"].map(_to_proforma)
    exp.peptides["Peptide"] = exp.peptides["Peptide"].map(_to_proforma)
    exp.psms["Source File"] = exp.psms["Source File"].str.replace(".raw", ".mzML", regex=False)
    assert exp.psms["Peptide"].str.contains(r"\[").any()
    write_export(exp, tmp_path)
    m = cc.compare(d, str(tmp_path))["metrics"]
    _assert_perfect(m)


def test_perturbations_are_classified(tmp_path):
    """Edit a copy of Sample 1 in known ways and check every count."""
    d = require_ref(1)
    exp = read_export(d)
    p = exp.psms
    key = ["Source File", "Scan", "Peptide", "z"]
    single = p[~p.duplicated(key, keep=False)]
    # scans holding exactly one PSM row, so edits do not collide with chimeric partners
    one_per_scan = single[~single.duplicated("Scan", keep=False)]
    plain = one_per_scan[~one_per_scan["Peptide"].str.contains(r"\(")]
    has_i = plain[plain["Peptide"].str.contains("I")]
    has_nq = plain[plain["Peptide"].str.contains("[NQ]") & ~plain["Peptide"].str.contains("I")]
    others = plain.drop(has_i.index[:1]).drop(has_nq.index[:1])

    i_charge, i_il, i_mod, i_drop, i_other = (
        others.index[0], has_i.index[0], has_nq.index[0], others.index[1], others.index[2])
    t = p.copy()
    t.loc[i_charge, "z"] = t.loc[i_charge, "z"] + 1
    t.loc[i_il, "Peptide"] = t.loc[i_il, "Peptide"].replace("I", "L", 1)
    t.loc[i_mod, "Peptide"] = re.sub(r"([NQ])", r"\1(+0.98)", t.loc[i_mod, "Peptide"], count=1)
    t.loc[i_other, "Peptide"] = "GGGGGGGGGGK"
    t = t.drop(index=i_drop)
    exp.psms = t
    write_export(exp, tmp_path, tables=["psms", "peptides", "proteins"])

    m = cc.compare(d, str(tmp_path))["metrics"]
    psm = m["levels"]["psm"]
    n = psm["n_ref"]
    assert psm["ref_only"] == 5 and psm["test_only"] == 4
    assert psm["shared"] == n - 5
    assert psm["n_test"] == n - 1
    rb = m["psm_ref_only_breakdown"]
    assert rb["same_scan_same_peptide_other_charge"] == 1
    assert rb["same_scan_IL_variant"] == 1
    assert rb["same_scan_same_sequence_other_mods"] == 1
    assert rb["same_scan_other_peptide"] == 1
    # the dropped PSM's peptide is still identified at other scans or not at all
    assert rb["other_scan_same_peptide"] + rb["not_found"] + rb["other_scan_same_sequence"] == 1
    assert m["psm_same_scan_different_peptide"]["ref_only"] == 3
    tb = m["psm_test_only_breakdown"]
    assert tb["same_scan_same_peptide_other_charge"] == 1
    assert tb["same_scan_IL_variant"] == 1
    assert tb["same_scan_same_sequence_other_mods"] == 1
    assert tb["same_scan_other_peptide"] == 1
    assert m["levels"]["peptide"]["recall"] == 1.0  # peptides.csv unchanged
    assert m["levels"]["protein_group"]["recall"] == 1.0
    assert abs(m["headline"]["psm_recall"] - (n - 5) / n) < 1e-12
    assert abs(m["headline"]["psm_precision"] - (n - 5) / (n - 1)) < 1e-12


def test_protein_groups_match_by_any_accession(tmp_path):
    d = require_ref(1)
    exp = read_export(d, tables=["psms", "peptides", "proteins"])
    pr = exp.proteins.copy()
    # split group 1 into two test groups; drop one non-top member elsewhere
    g1 = pr["Protein Group"] == 1
    idx = pr.index[g1]
    pr.loc[idx[: len(idx) // 2], "Protein Group"] = 999
    # rename the only accession of one singleton group -> ref-only + test-only
    single = pr.groupby("Protein Group").filter(lambda g: len(g) == 1)
    j = single.index[0]
    old_acc = pr.loc[j, "Accession"]
    pr.loc[j, "Accession"] = "X00000|FAKE_HUMAN"
    exp.proteins = pr
    write_export(exp, tmp_path)
    m = cc.compare(d, str(tmp_path))["metrics"]
    pg = m["levels"]["protein_group"]
    n = pg["n_ref"]
    assert pg["n_test"] == n + 1
    assert pg["ref_only"] == 1 and pg["test_only"] == 1
    assert pg["shared"] == n - 1  # the split group is one component
    assert pg["recall"] == (n - 1) / n
    assert pg["precision"] == n / (n + 1)
    top = m["levels"]["protein_top_accession"]
    was_top = bool(exp.proteins.loc[j, "Top"])
    assert top["ref_only"] == int(was_top)
    assert old_acc not in set(pd.read_csv(tmp_path / "db.proteins.csv")["Accession"])


def test_psm_min_score_filter():
    d = require_ref(1)
    m_all = cc.compare(d, d)["metrics"]
    m_thr = cc.compare(d, d, min_10lgp="ref")["metrics"]
    assert m_thr["options"]["min_10lgp_applied"] == 28.36
    assert m_thr["levels"]["psm"]["n_ref"] < m_all["levels"]["psm"]["n_ref"]
    assert m_thr["headline"]["psm_recall"] == 1.0


def test_norm_helpers():
    assert cc.norm_source("C:\\data\\DrewR_Adam_91526_DDA_OT_1.raw") == "DrewR_Adam_91526_DDA_OT_1"
    assert cc.norm_source("/x/DrewR_Adam_91526_DDA_OT_1.mzML.gz") == "DrewR_Adam_91526_DDA_OT_1"
    assert cc.norm_accession("sp|P02768|ALBU_HUMAN") == "P02768|ALBU_HUMAN"
    assert cc.norm_accession("P02768|ALBU_HUMAN") == "P02768|ALBU_HUMAN"
    s = cc.set_metrics({1, 2, 3, 4}, {3, 4, 5})
    assert s == {"n_ref": 4, "n_test": 3, "shared": 2, "ref_only": 2, "test_only": 1,
                 "recall": 0.5, "precision": 2 / 3, "jaccard": 0.4}
    assert cc.set_metrics(set(), set())["recall"] is None
