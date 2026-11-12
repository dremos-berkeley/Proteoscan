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


def _acc(protein_id: str, twin_prefix: str = "") -> tuple[str, bool]:
    """(accession, is_decoy). With `twin_prefix` an entrapment twin keeps the prefix in front of its accession."""
    p, decoy = protein_id, False
    for t in DECOY_PREFIXES:
        if p.startswith(t):
            p, decoy = p[len(t):], True
    twin = bool(twin_prefix) and p.startswith(twin_prefix)
    if twin:
        p = p[len(twin_prefix):]
    parts = p.split("|")
    acc = parts[1] if len(parts) >= 3 and parts[0] in ("sp", "tr") else parts[0]
    return (twin_prefix + acc if twin else acc), decoy


def protein_support(accepted: pd.DataFrame, twin_prefix: str = "") -> tuple[dict[str, int], dict[str, set]]:
    """accepted: pass-1 accepted target PSM candidates (columns stripped, proteins).
    Returns ({accession: n distinct sequences}, {accession: set(sequences)}).

    Entrapment validation (`twin_prefix`): only real targets collect support, and the twin of each protein
    inherits its count, exactly as a decoy inherits the support of the target it was derived from, so the
    feature cannot separate targets from entrapment sequences by itself."""
    seqs: dict[str, set] = defaultdict(set)
    for stripped, proteins in zip(accepted["stripped"], accepted["proteins"]):
        key = stripped.replace("I", "L")
        for p in str(proteins).split(";"):
            a, decoy = _acc(p, twin_prefix)
            if not decoy and not (twin_prefix and a.startswith(twin_prefix)):
                seqs[a].add(key)
    support = {a: len(s) for a, s in seqs.items()}
    if twin_prefix:
        for a in list(support):
            support[twin_prefix + a] = support[a]
    return support, seqs


def add_protein_features(cands: pd.DataFrame, support: dict[str, int], seqs: dict[str, set],
                         twin_prefix: str = "") -> pd.DataFrame:
    out = cands.copy()
    vals = np.zeros(len(out))
    for i, (stripped, proteins) in enumerate(zip(out["stripped"], out["proteins"])):
        key = stripped.replace("I", "L")
        best = 0
        for p in str(proteins).split(";"):
            a, decoy = _acc(p, twin_prefix)
            n = support.get(a, 0)
            if not decoy and key in seqs.get(a, ()):
                n -= 1  # do not let a peptide support itself
            best = max(best, n)
        vals[i] = best
    out["prot_support"] = vals
    out["log_prot_support"] = np.log1p(vals)
    return out


def rerank(both: pd.DataFrame, keep: int = 5, query_col: str = "query_id") -> pd.DataFrame:
    """Recompute rank and the delta columns per query after candidates were merged."""
    both = both.sort_values([query_col, "hyperscore"], ascending=[True, False])
    both["rank"] = both.groupby(query_col).cumcount() + 1
    both = both[both["rank"] <= keep].copy()
    best = both.groupby(query_col)["hyperscore"].transform("max")
    nxt = both.groupby(query_col)["hyperscore"].shift(-1)
    both["delta_best"] = best - both["hyperscore"]
    both["delta_next"] = (both["hyperscore"] - nxt).fillna(both["hyperscore"])
    return both.reset_index(drop=True)


def merge_passes(pass1: pd.DataFrame, pass2: pd.DataFrame, keep: int = 5,
                 query_col: str = "query_id") -> pd.DataFrame:
    """Union of candidates per query; identical (query, peptide) rows are collapsed, pass-1 rows
    (which may already carry features) win."""
    p1 = pass1.copy(); p1["from_pass2"] = 0
    p2 = pass2.copy(); p2["from_pass2"] = 1
    both = pd.concat([p1, p2], ignore_index=True, sort=False)
    both = both.sort_values(["from_pass2"]).drop_duplicates([query_col, "peptide_sage"], keep="first")
    return rerank(both, keep, query_col)
