"""Protein inference imitating PEAKS Studio (stage 7): protein -10lgP, grouping, Top flag,
#Peptides / #Unique, coverage, #Spec, Area, Average Mass, Description and Gene.

Input is a table of identified peptides in the PEAKS ``db.peptides.csv`` layout (one row per
peptidoform, best PSM) and the protein FASTA.  Output is a pair of DataFrames with exactly
the columns of PEAKS ``db.proteins.csv`` and ``db.protein-peptides.csv``
(:func:`peaksparity.peaks_export.canonical_columns`).

Every rule below was established on the four DrewR_Adam_91526 reference exports; see
docs/PROTEIN_INFERENCE.md for the evidence.  Summary:

Peptide -> protein mapping
    A peptide supports a protein when its residue sequence occurs in the protein with I and L
    treated as identical.  PEAKS writes one peptide row per I/L variant that exists in the
    database (e.g. AILFLPR for CCL21 and ALLFIPR for HS90B, same spectrum and score), and a
    protein lists only the variant(s) that occur verbatim in it.  Domination and uniqueness
    work on the I/L-collapsed sequence ("peptide family").  Start/End are the 1-based first
    occurrence; coverage uses only that occurrence.

Protein -10lgP
    sum_k s_k / k over the protein's distinct residue sequences (not peptidoforms: all PTM
    forms and positions of one sequence count once, at their best score), s_1 >= s_2 >= ...,
    in double precision, reported as a Java float.  (This differs from the published rule,
    which counts different PTM sets separately; the exports show all PTM forms collapse.)

Grouping
    Protein X dominates Y when every peptide family of Y with score >= T (T = 15) is also a
    family of X.  See :func:`group_proteins` for the exact greedy procedure, the Top rule and
    tie handling.

Unique, counts, Area, #Spec
    A peptide family is unique when all the proteins it maps to lie in one group.
    #Peptides / #Unique count distinct (Start, End) of the protein's supporting / unique
    rows (so PTM forms of one sequence count once).  #Spec counts the PSMs of all supporting
    peptidoforms.  Area is the float32 sum of the areas of the distinct MS1 features of the
    unique supporting peptidoforms (empty when there is none).

Filter
    Proteins with -10lgP < 15 or #Unique < 1 are removed after grouping; peptides keep the
    uniqueness computed before the filter.

Average Mass
    Sum of average residue masses (classic table, see ``AVERAGE_RESIDUE``) plus the
    *monoisotopic* mass of water (18.01056), cast to float32.  X = mean of the 20 residues.

Description / Gene
    Description is the FASTA header after the first token; Gene the ``GN=`` value.
"""

from __future__ import annotations

import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd

from .peaks_export import canonical_columns
from .seqnorm import parse

__all__ = [
    "AVERAGE_RESIDUE",
    "ProteinDB",
    "read_protein_db",
    "average_mass",
    "digest_positions",
    "map_peptides",
    "Enzyme",
    "extended_peptide",
    "protein_score",
    "group_proteins",
    "infer_proteins",
    "read_mzid_psm_scores",
]

#: average residue masses PEAKS uses (fitted on 1336 exported proteins; the classic
#: Expasy-style table).  X is the mean of the 20 standard residues, which reproduces
#: PEAKS for a protein containing X; U, B, Z, O are our own choices (no evidence).
AVERAGE_RESIDUE = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
AVERAGE_RESIDUE["X"] = sum(AVERAGE_RESIDUE.values()) / 20.0
AVERAGE_RESIDUE.update({"U": 150.0388, "O": 237.3018,
                        "B": (114.1038 + 115.0886) / 2, "Z": (128.1307 + 129.1155) / 2,
                        "J": 113.1594})
WATER_MONO = 18.010565


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------

_GN_RE = re.compile(r"\bGN=(\S+)")


def _accession(token):
    """PEAKS accession from the first header token: ``sp|P14678|RSMB_HUMAN`` -> ``P14678|RSMB_HUMAN``."""
    parts = token.split("|")
    if len(parts) >= 3 and parts[0] in ("sp", "tr", "gi", "ref", "gb", "emb", "dbj", "pdb", "pir", "prf"):
        return "|".join(parts[1:])
    return token


class ProteinDB:
    """Proteins of a FASTA file: parallel lists ``accession``, ``description``, ``gene``,
    ``sequence`` plus ``index`` (accession -> position)."""

    def __init__(self, entries):
        self.accession, self.description, self.gene, self.sequence = [], [], [], []
        self.index = {}
        for header, seq in entries:
            head = header[1:] if header.startswith(">") else header
            tok, _, rest = head.strip().partition(" ")
            acc = _accession(tok)
            if acc in self.index:
                continue
            m = _GN_RE.search(rest)
            self.index[acc] = len(self.accession)
            self.accession.append(acc)
            self.description.append(rest.strip())
            self.gene.append(m.group(1) if m else "")
            self.sequence.append(seq.upper().replace("*", ""))
        self._collapsed = None
        self._automaton = None

    def __len__(self):
        return len(self.accession)

    @property
    def collapsed(self):
        if self._collapsed is None:
            self._collapsed = [s.replace("I", "L") for s in self.sequence]
        return self._collapsed

    def find(self, peptides):
        """All occurrences of stripped peptide sequences, with I == L.

        Returns {peptide: [(protein index, [0-based starts, ascending]), ...]} with proteins in
        database order."""
        peps = sorted(set(peptides))
        by_key = defaultdict(list)
        for p in peps:
            by_key[p.replace("I", "L")].append(p)
        hits = defaultdict(lambda: defaultdict(list))  # key -> protein -> starts
        try:
            import ahocorasick
        except ImportError:  # pragma: no cover - slow fallback
            ahocorasick = None
        if ahocorasick is not None and by_key:
            A = ahocorasick.Automaton()
            for k in by_key:
                A.add_word(k, k)
            A.make_automaton()
            for i, s in enumerate(self.collapsed):
                for end, k in A.iter(s):
                    hits[k][i].append(end - len(k) + 1)
        else:
            for k in by_key:
                for i, s in enumerate(self.collapsed):
                    st = s.find(k)
                    while st >= 0:
                        hits[k][i].append(st)
                        st = s.find(k, st + 1)
        out = {}
        for k, ps in by_key.items():
            lst = [(i, sorted(v)) for i, v in sorted(hits.get(k, {}).items())]
            for p in ps:
                out[p] = lst
        return out


class Enzyme:
    """Digestion rule used to decide whether an occurrence of a peptide in a protein counts.

    PEAKS assigns a peptide only to proteins in which it is a valid digestion product under
    the search settings (e.g. semi-specific trypsin, <= 2 missed cleavages): the exported
    Accession lists never contain proteins where the occurrence has two non-enzymatic
    termini, and they follow the "not before P" trypsin rule (AVFPSIVGR is not assigned to
    ACTBM where it is followed by P...).  A protein N-terminus, and position 2 after an
    initiator Met, count as enzymatic termini.
    """

    def __init__(self, cleave_after="KR", not_before="P", specificity="semi", max_missed=2):
        self.cleave_after = cleave_after
        self.not_before = not_before or ""
        self.specificity = specificity  # "specific", "semi", "unspecific"/"none"
        self.max_missed = max_missed

    @classmethod
    def from_params(cls, params):
        """From :func:`peaksparity.peaks_export.parse_summary_html` output."""
        params = params or {}
        enz = str(params.get("enzyme", "Trypsin")).lower()
        mode = str(params.get("digest_mode", "Semi-Specific")).lower()
        spec = "semi" if "semi" in mode else ("specific" if "specific" in mode else "unspecific")
        if enz in ("trypsin", "trypsin/p") or "trypsin" in enz:
            ca, nb = "KR", ("" if enz.endswith("/p") else "P")
        elif enz.startswith("lys-c") or enz.startswith("lysc"):
            ca, nb = "K", "P"
        else:
            spec = "unspecific"
            ca, nb = "", ""
        mm = params.get("max_missed_cleavages", 2)
        return cls(ca, nb, spec, 2 if mm is None else int(mm))

    def _site(self, seq, j):
        """True if the protein can be cleaved between seq[j-1] and seq[j]."""
        if j <= 0 or j >= len(seq):
            return True
        return seq[j - 1] in self.cleave_after and seq[j] not in self.not_before

    def valid(self, protein, start, length):
        if self.specificity in ("unspecific", "none") or not self.cleave_after:
            return True
        end = start + length
        nt = self._site(protein, start) or (start == 1 and protein[0] == "M")
        ct = self._site(protein, end)
        if self.specificity == "specific" and not (nt and ct):
            return False
        if self.specificity == "semi" and not (nt or ct):
            return False
        missed = sum(1 for j in range(start + 1, end) if self._site(protein, j))
        return missed <= self.max_missed


DEFAULT_ENZYME = Enzyme()


def _read_fasta(path):
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line, []
            elif line:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def read_protein_db(fasta, exclude_prefix="#DECOY#"):
    """Read a FASTA file (or pass through a :class:`ProteinDB`)."""
    if isinstance(fasta, ProteinDB):
        return fasta
    entries = [(h, s) for h, s in _read_fasta(fasta)
               if not (exclude_prefix and h[1:].startswith(exclude_prefix))]
    return ProteinDB(entries)


def average_mass(seq):
    """PEAKS "Average Mass": average residue masses + monoisotopic water, as float32."""
    return float(np.float32(sum(AVERAGE_RESIDUE.get(a, 0.0) for a in seq) + WATER_MONO))


def extended_peptide(peptide, protein_seq, start0, length):
    """PEAKS "Extended Peptide": ``K.PEPTIDE.R``; the flank and its dot are omitted at a
    protein terminus (``MAGSGAGVR.C``, ``R.TLYGFGG``)."""
    pre = protein_seq[start0 - 1] + "." if start0 > 0 else ""
    end = start0 + length
    post = "." + protein_seq[end] if end < len(protein_seq) else ""
    return f"{pre}{peptide}{post}"


def map_peptides(db, stripped, enzyme=DEFAULT_ENZYME):
    """Map stripped peptide sequences to proteins of ``db``.

    Returns {peptide: [(protein index, 0-based start, exact), ...]}.  Only occurrences that
    are valid digestion products under ``enzyme`` count.  ``exact`` is True when the protein
    has the peptide verbatim (then start is the first valid verbatim occurrence), False when
    it matches only with I == L (first valid I/L occurrence).
    """
    enzyme = enzyme or Enzyme(specificity="unspecific")
    out = {}
    for pep, lst in db.find(stripped).items():
        n = len(pep)
        res = []
        for i, starts in lst:
            seq = db.sequence[i]
            first_il = None
            first_exact = None
            for st in starts:
                if not enzyme.valid(seq, st, n):
                    continue
                if first_il is None:
                    first_il = st
                if seq.startswith(pep, st):
                    first_exact = st
                    break
            if first_exact is not None:
                res.append((i, first_exact, True))
            elif first_il is not None:
                res.append((i, first_il, False))
        out[pep] = res
    return out


def digest_positions(fasta, peptides, enzyme=DEFAULT_ENZYME):
    """Positions of peptides in a FASTA (or :class:`ProteinDB`).

    ``peptides``: iterable of (modified) peptide strings.  Returns a DataFrame with one row per
    (peptide, protein) in which the peptide is a valid digestion product: Peptide, Accession,
    Start, End (1-based, first occurrence), Extended Peptide, and IL_only (True when the
    protein contains the sequence only with I/L swapped).
    """
    db = read_protein_db(fasta)
    peptides = list(dict.fromkeys(peptides))
    stripped = {p: parse(p)[0] for p in peptides}
    hits = map_peptides(db, set(stripped.values()), enzyme)
    rows = []
    for p in peptides:
        s = stripped[p]
        for i, st, exact in hits[s]:
            seq = db.sequence[i]
            rows.append((p, db.accession[i], st + 1, st + len(s),
                         extended_peptide(p, seq, st, len(s)), not exact))
    return pd.DataFrame(rows, columns=["Peptide", "Accession", "Start", "End",
                                       "Extended Peptide", "IL_only"])


# --------------------------------------------------------------------------
# Scores
# --------------------------------------------------------------------------


def protein_score(scores):
    """PEAKS protein -10lgP from the best scores of its distinct peptide sequences:
    sum_k s_k / k with s_1 >= s_2 >= ... (double precision)."""
    v = sorted((float(s) for s in scores), reverse=True)
    return float(sum(s / (k + 1) for k, s in enumerate(v)))


def read_mzid_psm_scores(path):
    """Full-precision PSM scores from a PEAKS ``db.mzidentml_1_2_0.xml``.

    Returns a DataFrame (Peptide, RT, score) with one row per SpectrumIdentificationItem;
    ``Peptide`` is the PEAKS modified sequence and ``RT`` the result's retention time text.
    The CSV exports round -10lgP to 2 decimals, which is not enough to reproduce protein
    -10lgP to 2 decimals; the mzIdentML carries the float value.
    """
    import xml.etree.ElementTree as ET

    rows = []
    for _ev, el in ET.iterparse(path):
        if not el.tag.endswith("SpectrumIdentificationResult"):
            continue
        rt = None
        for c in el:
            if c.tag.endswith("cvParam") and c.get("name") == "retention time":
                rt = c.get("value")
        for sii in el:
            if not sii.tag.endswith("SpectrumIdentificationItem"):
                continue
            pep = sii.get("peptide_ref", "")
            if pep.startswith("PEPTIDE_"):
                pep = pep[len("PEPTIDE_"):]
            sc = None
            for c in sii:
                if c.tag.endswith("cvParam") and c.get("name") == "PEAKS:peptideScore":
                    sc = float(c.get("value"))
            rows.append((pep, rt, sc))
        el.clear()
    return pd.DataFrame(rows, columns=["Peptide", "RT", "score"])
