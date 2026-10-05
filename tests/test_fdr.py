"""Tests for proteoscan.fdr (PEAKS-style target/decoy thresholds)."""
import math

import numpy as np
import pandas as pd
import pytest

from conftest import SAMPLES, require_ref
from proteoscan.fdr import (
    fdr_curve,
    peptide_fdr,
    protein_group_fdr,
    psm_fdr,
    q_values,
    sequence_key,
    threshold_at_fdr,
    threshold_summary,
)

S = [50, 40, 40, 30, 25, 20, 10]
D = [0, 0, 1, 0, 0, 1, 1]


def test_fdr_curve_is_decoys_over_targets():
    c = fdr_curve(S, D)
    assert list(c.threshold) == [50, 40, 30, 25, 20, 10]  # tied scores form one threshold
    assert list(c.n_target) == [1, 2, 3, 4, 4, 4]
    assert list(c.n_decoy) == [0, 1, 1, 1, 2, 3]
    assert c.fdr.tolist() == pytest.approx([0, 0.5, 1 / 3, 0.25, 0.5, 0.75])  # D/T, not 2D/(T+D)
    assert c.q_value.tolist() == pytest.approx([0, 0.25, 0.25, 0.25, 0.5, 0.75])


def test_threshold_lowest_passing_not_monotonised():
    # FDR at 40 is 0.5 > 0.34 but the lowest threshold with D/T <= 0.34 is 25
    assert threshold_at_fdr(S, D, 0.34) == 25
    assert threshold_at_fdr(S, D, 0.0) == 50
    assert threshold_at_fdr([10.0], [True], 0.01) == math.inf
    s = threshold_summary(S, D, 0.25)
    assert s == {"threshold": 25.0, "n_target": 4, "n_decoy": 1, "fdr": 0.25}


def test_q_values():
    q = q_values(S + [np.nan], D + [0])
    assert q[:7].tolist() == pytest.approx([0, 0.25, 0.25, 0.25, 0.25, 0.5, 0.75])
    assert np.isnan(q[7])


def test_psm_fdr_one_hit_per_query():
    psms = pd.DataFrame({
        "Scan": [1, 1, 2, 3, 3, 4],
        "Feature Id": [7, 8, 9, 10, 10, 11],  # scan 1 is chimeric: two queries
        "Peptide": ["AILK", "PEPK", "DECK", "LLLK", "ILLK", "XXK"],
        "-10LgP": [40, 35, 30, 25, 25, 20],
        "decoy": [False, False, True, False, False, True],
    })
    r = psm_fdr(psms, alpha=0.34)
    assert len(r["table"]) == 5  # (scan 3, feature 10) I/L variants -> one PSM
    assert r["threshold"] == 25 and r["n_target"] == 3 and r["n_decoy"] == 1
    arr = psm_fdr((np.array([3.0, 2.0, 1.0]), np.array([False, True, False])), alpha=0.5)
    assert arr["threshold"] == 1.0


def test_peptide_fdr_unique_sequences():
    assert sequence_key("LIM(+15.99)K") == "LLMK"
    psms = pd.DataFrame({
        "Peptide": ["PEPTIDEK", "PEPTLDEK", "PEPTIDEK", "M(+15.99)AAK", "MAAK", "DECOYK"],
        "-10LgP": [30, 45, 20, 25, 28, 26],
        "decoy": [False, False, False, False, False, True],
    })
    r = peptide_fdr(psms, alpha=1.0)
    t = r["table"].set_index("key")
    assert len(t) == 3  # PEPTLDEK (I/L) and PTM forms collapse
    assert t.loc["PEPTLDEK", "-10LgP"] == 45 and t.loc["MAAK", "-10LgP"] == 28
    assert peptide_fdr(psms, alpha=0.0)["threshold"] == 28


def test_protein_group_fdr_top_hit():
    prot = pd.DataFrame({
        "Protein Group": [1, 1, 2, 3, 3],
        "Top": [True, False, True, True, True],
        "-10LgP": [100.0, 120.0, 50.0, 40.0, 40.0],  # non-top rows never define the group
        "decoy": [False, True, True, False, False],
    })
    r = protein_group_fdr(prot, alpha=0.5)
    assert len(r["table"]) == 3
    assert r["n_decoy"] == 1 and r["n_target"] == 2 and r["threshold"] == 40.0
    assert protein_group_fdr(prot, alpha=0.0)["threshold"] == 100.0


@pytest.mark.parametrize("k", SAMPLES)
def test_reference_counting_units(k):
    """db.summary.html counts: #PSMs = distinct (Scan, Feature Id) pairs of db.psms.csv;
    #Sequences = distinct I/L-collapsed stripped sequences of db.peptides.csv."""
    d = require_ref(k)
    from proteoscan.peaks_export import read_export

    exp = read_export(d, tables=["psms", "peptides"], verify=False)
    stats = exp.params["statistics"]["Total"]
    assert len(exp.psms.drop_duplicates(["Scan", "Feature Id"])) == stats["psms"]
    assert exp.peptides.Peptide.map(sequence_key).nunique() == stats["sequences"]
    assert exp.peptides.Peptide.str.replace("I", "L").nunique() == stats["peptides"]


def test_pseudocount():
    c = fdr_curve(S, D, pseudocount=1)
    assert c.fdr.tolist() == pytest.approx([1.0, 1.0, 2 / 3, 0.5, 0.75, 1.0])  # (D+1)/T
    assert threshold_summary(S, D, 0.5, pseudocount=1) == {"threshold": 25.0, "n_target": 4, "n_decoy": 1, "fdr": 0.5}
    assert threshold_at_fdr(S, D, 0.34, pseudocount=1) == math.inf
    assert q_values(S, D, pseudocount=1).tolist() == pytest.approx([0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 1.0])
    assert threshold_at_fdr(S, D, 0.34) == threshold_at_fdr(S, D, 0.34, pseudocount=0) == 25


@pytest.mark.parametrize("seed", range(20))
def test_threshold_equals_qvalue_thresholding(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(5, 400))
    s = np.round(rng.normal(0, 1, n), 1)  # rounding creates ties
    d = rng.random(n) < rng.uniform(0.05, 0.6)
    s[~d] += rng.uniform(0, 2)
    for alpha in (0.0, 0.01, 0.05, 0.2):
        for pc in (0, 1):
            thr = threshold_at_fdr(s, d, alpha, pseudocount=pc)
            q = q_values(s, d, pseudocount=pc)
            passing = s[q <= alpha]
            assert thr == (passing.min() if len(passing) else math.inf)
