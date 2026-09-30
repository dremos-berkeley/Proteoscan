"""Two-pass search with a protein-support feature (PEAKS DB's "protein feature", decoy-fused).

PEAKS DB shortlists proteins, then scores peptides with a bonus for peptides of proteins that
have other good peptides; its decoy-fusion design gives each decoy peptide the bonus of the
target protein it was derived from, so the bonus cannot separate targets from decoys by itself.

Here: pass 1 = full-database search + rescoring + filter. Pass 2 = search of every query against
only the proteins accepted in pass 1 (plus their decoys). Candidates of both passes are merged and
each gets `prot_support`: the number of distinct pass-1 peptide sequences of its best-supported
protein, where a decoy protein inherits the count of its target (decoy fusion).
"""
from __future__ import annotations

import re
from collections import defaultdict

import numpy as np
import pandas as pd

DECOY_PREFIXES = ("#DECOY#", "DECOY_", "rev_")


def _acc(protein_id: str) -> tuple[str, bool]:
    p, decoy = protein_id, False
    for t in DECOY_PREFIXES:
        if p.startswith(t):
            p, decoy = p[len(t):], True
    parts = p.split("|")
    return (parts[1] if len(parts) >= 3 and parts[0] in ("sp", "tr") else parts[0]), decoy


def protein_support(accepted: pd.DataFrame) -> tuple[dict[str, int], dict[str, set]]:
    """accepted: pass-1 accepted target PSM candidates (columns stripped, proteins).
    Returns ({accession: n distinct sequences}, {accession: set(sequences)})."""
    seqs: dict[str, set] = defaultdict(set)
    for stripped, proteins in zip(accepted["stripped"], accepted["proteins"]):
        key = stripped.replace("I", "L")
        for p in str(proteins).split(";"):
            a, decoy = _acc(p)
            if not decoy:
                seqs[a].add(key)
    return {a: len(s) for a, s in seqs.items()}, seqs


def add_protein_features(cands: pd.DataFrame, support: dict[str, int], seqs: dict[str, set]) -> pd.DataFrame:
    out = cands.copy()
    vals = np.zeros(len(out))
    for i, (stripped, proteins) in enumerate(zip(out["stripped"], out["proteins"])):
        key = stripped.replace("I", "L")
        best = 0
        for p in str(proteins).split(";"):
            a, decoy = _acc(p)
            n = support.get(a, 0)
            if not decoy and key in seqs.get(a, ()):
                n -= 1  # do not let a peptide support itself
            best = max(best, n)
        vals[i] = best
    out["prot_support"] = vals
    out["log_prot_support"] = np.log1p(vals)
    return out


def merge_passes(pass1: pd.DataFrame, pass2: pd.DataFrame, keep: int = 5,
                 query_col: str = "query_id") -> pd.DataFrame:
    """Union of candidates per query; identical (query, peptide) rows are collapsed."""
    p1 = pass1.copy(); p1["from_pass2"] = 0
    p2 = pass2.copy(); p2["from_pass2"] = 1
    both = pd.concat([p1, p2], ignore_index=True, sort=False)
    both = both.sort_values(["from_pass2"]).drop_duplicates([query_col, "peptide_sage"], keep="first")
    both = both.sort_values([query_col, "hyperscore"], ascending=[True, False])
    both["rank"] = both.groupby(query_col).cumcount() + 1
    both = both[both["rank"] <= keep].copy()
    best = both.groupby(query_col)["hyperscore"].transform("max")
    nxt = both.groupby(query_col)["hyperscore"].shift(-1)
    both["delta_best"] = best - both["hyperscore"]
    both["delta_next"] = (both["hyperscore"] - nxt).fillna(both["hyperscore"])
    return both.reset_index(drop=True)
