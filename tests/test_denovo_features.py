"""Tests for proteoscan.denovo_features (CAA / tag features against de novo calls)."""
import numpy as np
import pandas as pd
import pytest

from proteoscan.denovo_features import DENOVO_FEATURE_COLUMNS, add_denovo_features, denovo_features


def _dn(rows):
    return pd.DataFrame([{"query_id": q, "sequence": s, "per_residue_scores": np.asarray(c, float), "picked": p}
                         for q, s, c, p in rows])


def test_exact_and_il():
    dn = _dn([(1, "PEPTLDEK", [0.9] * 8, True)])
    c = pd.DataFrame({"query_id": [1], "peptide_sage": ["PEPTIDEK"]})
    f = denovo_features(c, dn).iloc[0]
    assert f.caa == 8 and f.tag_len == 8 and f.caa_pct == 100 and f.denovo_exact == 1
    assert f.denovo_conf_mean == pytest.approx(0.9)


def test_low_confidence_collapses_and_breaks_tag():
    dn = _dn([(1, "GGPEPTIDEK", [.9, .9, .9, .2, .9, .9, .9, .9, .9, .9], True)])
    c = pd.DataFrame({"query_id": [1], "peptide_sage": ["GGPEPTIDEK"]})
    f = denovo_features(c, dn).iloc[0]
    assert f.caa == 9 and f.tag_len == 6


def test_mass_blocks_and_c_anchor():
    # N (114.043) vs GG (114.043): same block mass, different residues -> not common
    dn = _dn([(1, "NPEPTIDEK", [0.9] * 9, True),
              # W vs M(ox): different mass inside; N-anchored AAA + C-anchored PEPGGGK
              (2, "AAAWPEPGGGK", [0.9] * 11, True)])
    c = pd.DataFrame({"query_id": [1, 2], "peptide_sage": ["GGPEPTIDEK", "AAAM[+15.9949]PEPGGGK"]})
    f = denovo_features(c, dn)
    assert list(f.caa) == [8, 10]
    assert list(f.tag_len) == [8, 7]


def test_no_double_counting_and_beams():
    dn = _dn([(1, "PEPTLDEKK", [0.95] * 9, True), (1, "PEPTLDEK", [0.9] * 8, False)])
    c = pd.DataFrame({"query_id": [1, 2], "peptide_sage": ["PEPTIDEK", "PEPTIDEK"]})
    f = denovo_features(c, dn)
    assert f.caa.iloc[0] <= 8 and f.max_caa_any_beam.iloc[0] == 8
    assert f.denovo_exact.iloc[0] == 0
    assert f.iloc[1].isna().all()  # no de novo call for query 2


def test_nterm_and_mods_folded():
    dn = _dn([(1, "[Acetyl]-M[Oxidation]PEPTK", [0.9, 0.8, 0.9, 0.9, 0.9, 0.9, 0.2], True)])
    c = pd.DataFrame({"query_id": [1, 1], "peptide_sage": ["[+42.0106]-M[+15.9949]PEPTK", "M[+15.9949]PEPTK"]})
    f = denovo_features(c, dn)
    assert f.denovo_exact.tolist() == [1, 0]
    assert f.caa.iloc[0] == 5 and f.denovo_conf_mean.iloc[0] == pytest.approx(0.88)


def test_add_columns():
    dn = _dn([(1, "PEPTIDEK", [0.9] * 8, True)])
    c = pd.DataFrame({"query_id": [1], "peptide_sage": ["PEPTIDEK"], "x": [3]})
    out = add_denovo_features(c, dn)
    assert set(DENOVO_FEATURE_COLUMNS) <= set(out.columns) and out["x"].iloc[0] == 3
