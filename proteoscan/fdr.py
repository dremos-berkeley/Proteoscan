"""PEAKS-style target/decoy FDR and threshold selection (stage 6).

Rules (PEAKS DB paper, Zhang et al. 2012 MCP; PEAKS Studio 12 manual sections 6.1 / 6.3.6;
docs/PEAKS_ALGORITHMS.md section 2.2; counting units checked against the four
DrewR_Adam_91526 exports, see docs/PROTEIN_INFERENCE.md):

* FDR at a score threshold t is ``D(>= t) / T(>= t)``: decoy hits over *target* hits
  scoring at least t.  It is not ``2D / (T + D)`` and not ``(D + 1) / T``.  With PEAKS's
  decoy fusion the decoy half of every fused entry receives the same protein-level bonuses
  as the target half, which is what makes the plain ratio unbiased.
* PSM level: one hit per query.  A query is an MS2 spectrum x precursor feature (chimeric
  spectra give several queries); I/L variants of the same hit are one PSM.  (The
  "#PSMs" of db.summary.html equals the number of distinct (Scan, Feature Id) pairs in
  db.psms.csv in all four references.)
* Peptide level: "Peptides are scored based on their top PSM"; the FDR counts *unique
  sequences*: residue sequences with I = L and modifications ignored (db.summary.html
  "#Sequences" = distinct I/L-collapsed stripped sequences of db.peptides.csv, and its
  filter line reports "FDR (Peptide Sequences)").
* Protein-group level: a group counts as decoy when its top protein is a decoy; the group
  score is the top protein's -10lgP.
* Threshold: PEAKS reports the score threshold chosen for a target FDR (e.g. peptide
  -10lgP >= 23.02 at 1 % peptide FDR).  We return the *lowest* observed score t such that
  D(>= t) / T(>= t) <= alpha.  This is the same accepted set as conventional q-value
  thresholding: with q(x) = min over thresholds t <= x of FDR(t), a hit has q <= alpha exactly
  when some threshold at or below its score passes, i.e. when its score >= the lowest passing
  threshold.  So there is no separate "monotonised" variant; ``q_values`` gives the per-hit
  q-values for convenience (tests/test_fdr.py checks the equivalence).
* Pseudocount: every function takes a keyword-only ``pseudocount`` (default 0) and then uses
  (D + pseudocount) / T.  The accepted-set criterion stays D/T (PEAKS parity; validated by
  entrapment, docs/STATUS.md); (D + 1)/T is the conservative estimate reported alongside it,
  which matters when few or no decoys pass (small defined mixtures, README).
* Ties: all hits with the same score are accepted or rejected together.
* Decoys are excluded from every reported count.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

__all__ = [
    "fdr_curve",
    "threshold_at_fdr",
    "threshold_summary",
    "q_values",
    "psm_fdr",
    "peptide_fdr",
    "protein_group_fdr",
    "sequence_key",
]


def _arrays(scores, is_decoy):
    s = np.asarray(scores, dtype=float)
    d = np.asarray(is_decoy, dtype=bool)
    if s.shape != d.shape:
        raise ValueError("scores and is_decoy must have the same length")
    ok = ~np.isnan(s)
    return s[ok], d[ok]


def fdr_curve(scores, is_decoy, *, pseudocount: int = 0):
    """FDR at every distinct score threshold.

    Returns a DataFrame sorted by threshold descending with columns ``threshold``,
    ``n_target`` = T(>= t), ``n_decoy`` = D(>= t), ``fdr`` = (D + pseudocount)/T (inf when
    T = 0), ``q_value`` = min over thresholds t' <= t of fdr(t').
    """
    s, d = _arrays(scores, is_decoy)
    if len(s) == 0:
        return pd.DataFrame(columns=["threshold", "n_target", "n_decoy", "fdr", "q_value"])
    order = np.argsort(-s, kind="mergesort")
    s, d = s[order], d[order]
    cum_d = np.cumsum(d)
    cum_t = np.cumsum(~d)
    # last index of each distinct score (ties accepted together)
    last = np.r_[np.nonzero(np.diff(s))[0], len(s) - 1]
    thr = s[last]
    nt = cum_t[last]
    nd = cum_d[last]
    with np.errstate(divide="ignore", invalid="ignore"):
        fdr = np.where(nt > 0, (nd + pseudocount) / np.maximum(nt, 1), np.inf)
    q = np.minimum.accumulate(fdr[::-1])[::-1]
    return pd.DataFrame({"threshold": thr, "n_target": nt, "n_decoy": nd, "fdr": fdr, "q_value": q})


def threshold_summary(scores, is_decoy, alpha, *, pseudocount: int = 0):
    """Lowest score threshold t with (D(>= t) + pseudocount) / T(>= t) <= alpha, with its counts.

    Returns a dict ``threshold`` (``inf`` when no threshold qualifies, so that nothing
    passes), ``n_target``, ``n_decoy`` and ``fdr`` at that threshold.  Hits pass when
    ``score >= threshold``.
    """
    c = fdr_curve(scores, is_decoy, pseudocount=pseudocount)
    ok = c[(c.fdr <= alpha) & (c.n_target > 0)] if len(c) else c
    if ok.empty:
        return {"threshold": math.inf, "n_target": 0, "n_decoy": 0, "fdr": 0.0}
    r = ok.iloc[-1]  # sorted descending -> last = lowest threshold
    return {"threshold": float(r.threshold), "n_target": int(r.n_target),
            "n_decoy": int(r.n_decoy), "fdr": float(r.fdr)}


def threshold_at_fdr(scores, is_decoy, alpha, *, pseudocount: int = 0):
    """Lowest observed score t with (D(>= t) + pseudocount) / T(>= t) <= alpha (``inf`` if none).

    Accept hits with ``score >= threshold_at_fdr(...)``; equivalent to accepting hits with
    ``q_values(...) <= alpha`` (module docstring).
    """
    return threshold_summary(scores, is_decoy, alpha, pseudocount=pseudocount)["threshold"]


def q_values(scores, is_decoy, *, pseudocount: int = 0):
    """q-value per hit: min over thresholds t <= score of (D(>= t) + pseudocount)/T(>= t)."""
    s = np.asarray(scores, dtype=float)
    c = fdr_curve(s, is_decoy, pseudocount=pseudocount)
    if c.empty:
        return np.full(len(s), np.nan)
    thr = c.threshold.to_numpy()[::-1]  # ascending
    q = c.q_value.to_numpy()[::-1]
    idx = np.searchsorted(thr, s, side="right") - 1
    out = np.where(idx >= 0, q[np.clip(idx, 0, None)], np.nan)
    out[np.isnan(s)] = np.nan
    return out


def _result(df, score_col, decoy_col, alpha, extra=None):
    res = threshold_summary(df[score_col], df[decoy_col], alpha)
    table = df.copy()
    table["q_value"] = q_values(table[score_col], table[decoy_col])
    table["pass"] = (table[score_col] >= res["threshold"]) & (~table[decoy_col].astype(bool))
    res["table"] = table
    if extra:
        res.update(extra)
    return res


def sequence_key(peptide):
    """PEAKS peptide-FDR unit: residue sequence without modifications, I -> L."""
    from .seqnorm import strip_mods

    return strip_mods(peptide).replace("I", "L")


def psm_fdr(psms, alpha=0.01, score_col="-10LgP", decoy_col="decoy", query_cols=None):
    """PSM-level FDR.

    ``psms``: one row per candidate PSM (target and decoy), or a (scores, is_decoy) pair of
    arrays.  Only the best hit per query is kept; ``query_cols`` defaults to
    ("Source File", "Scan", "Feature Id") restricted to the columns present (PEAKS's query =
    spectrum x precursor feature).  I/L variants of one query collapse to one PSM.
    Returns :func:`threshold_summary`'s dict plus ``table`` (the per-query best PSMs with
    ``q_value`` and ``pass``).
    """
    if isinstance(psms, tuple):
        df = pd.DataFrame({score_col: psms[0], decoy_col: psms[1]})
        return _result(df, score_col, decoy_col, alpha)
    df = psms
    if query_cols is None:
        query_cols = [c for c in ("Source File", "Scan", "Feature Id") if c in df.columns]
    if query_cols:
        df = df.sort_values(score_col, ascending=False, kind="mergesort") \
               .drop_duplicates(list(query_cols), keep="first")
    return _result(df.reset_index(drop=True), score_col, decoy_col, alpha)


def peptide_fdr(psms, alpha=0.01, score_col="-10LgP", decoy_col="decoy", peptide_col="Peptide",
                key=sequence_key):
    """Peptide-level FDR over unique sequences scored by their best PSM.

    ``psms``: PSM (or peptide) rows with ``peptide_col``, ``score_col`` and ``decoy_col``;
    or a (scores, is_decoy) pair of arrays already at peptide level.  ``key`` maps a
    peptide string to the counting unit (default: stripped, I -> L).  A sequence seen as
    both target and decoy counts as target.  Returns the threshold dict plus ``table``
    (one row per sequence: ``key``, best ``score_col``, ``decoy``, ``q_value``, ``pass``).
    Peptidoforms pass when their own best PSM scores >= threshold.
    """
    if isinstance(psms, tuple):
        df = pd.DataFrame({score_col: psms[0], decoy_col: psms[1]})
        return _result(df, score_col, decoy_col, alpha)
    df = psms[[peptide_col, score_col, decoy_col]].copy()
    df["key"] = df[peptide_col].map(key)
    g = df.groupby("key", sort=False)
    per = pd.DataFrame({
        score_col: g[score_col].max(),
        decoy_col: ~g[decoy_col].apply(lambda v: (~v.astype(bool)).any()),
    }).reset_index()
    return _result(per, score_col, decoy_col, alpha)


def protein_group_fdr(groups, alpha=0.01, score_col="-10LgP", decoy_col="decoy",
                      group_col="Protein Group", top_col="Top"):
    """Protein-group FDR: a group is decoy when its top protein is a decoy.

    ``groups``: protein rows (target and decoy) with ``group_col``, ``score_col``,
    ``decoy_col`` and optionally ``top_col``; the group's top hit is its best-scoring Top
    protein.  Or a (group scores, top-is-decoy) pair of arrays.  Returns the threshold dict
    plus ``table`` (one row per group).
    """
    if isinstance(groups, tuple):
        df = pd.DataFrame({score_col: groups[0], decoy_col: groups[1]})
        return _result(df, score_col, decoy_col, alpha)
    df = groups
    if top_col in df.columns:
        df = df[df[top_col].astype(bool)]
    df = df.sort_values(score_col, ascending=False, kind="mergesort") \
           .drop_duplicates(group_col, keep="first")
    return _result(df[[group_col, score_col, decoy_col]].reset_index(drop=True),
                   score_col, decoy_col, alpha)
