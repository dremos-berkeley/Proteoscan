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
    treated as identical *and* that occurrence is a valid digestion product under the search
    settings (semi-specific trypsin, K/R not before P, <= 2 missed cleavages; a protein
    N-terminus or position 2 after an initiator Met is an enzymatic terminus).  PEAKS writes
    one peptide row per I/L variant that exists in the database (e.g. AILFLPR for CCL21 and
    ALLFIPR for HS90B, same spectrum and score); a protein lists only the variant(s) that
    occur verbatim in it.  Domination and uniqueness work on the I/L-collapsed sequence
    ("peptide family").  Start/End are the 1-based first occurrence; coverage uses only that
    occurrence.

Protein -10lgP
    sum_k s_k / k over the protein's distinct residue sequences (not peptidoforms: all PTM
    forms and positions of one sequence count once, at their best score), s_1 >= s_2 >= ...,
    accumulated in float32.  (The published rule counts different PTM sets separately; the
    exports show that all PTM forms of a sequence collapse.)

Grouping and numbering
    Proteins are processed in descending order of their number of supporting peptide rows
    (peptidoforms).  An unassigned protein X opens group n+1; every unassigned protein Y whose
    families with score >= T (T = 15) are a subset of X's joins it (X dominates Y).  Top =
    same family set as X.  Group numbers are the creation order *before* the protein filter,
    so groups whose proteins are all filtered out leave gaps.  PEAKS's order among proteins
    with equal row counts is not reproducible (it looks hash-based); we break ties by
    -10lgP, then accession.

Unique, counts, Area, #Spec
    A peptide family is unique when all the proteins it maps to (before filtering) lie in one
    group.  #Peptides / #Unique count distinct (Start, End) of the protein's supporting /
    unique rows (so PTM forms of one sequence count once).  #Spec counts the PSM rows of all
    supporting peptidoforms.  Area is the float32 sum of the areas of the distinct MS1
    features of the unique supporting peptidoforms (empty when there is none), printed with
    3 significant digits rounding HALF_UP.  A feature counts for the protein tables only if
    one of its PSMs has -10lgP >= 20.

Filter
    Proteins with -10lgP < 15 or #Unique < 1 are removed after grouping.

Average Mass
    Sum of average residue masses (classic table, see ``AVERAGE_RESIDUE``) plus the
    monoisotopic mass of water (18.01056), cast to float32.  X = mean of the 20 residues.

Description / Gene
    Description is the FASTA header after the first token; Gene the ``GN=`` value.
"""

from __future__ import annotations

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
    "protein_peptide_mass",
    "MONO_RESIDUE_5",
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
# PEAKS's database stores selenocysteine as X (SELH_HUMAN: local U44, PEAKS X44), and
# ProteinDB converts U -> X on reading; U is mapped to the X mass for callers that do not.
AVERAGE_RESIDUE.update({"U": AVERAGE_RESIDUE["X"], "O": 237.3018,
                        "B": (114.1038 + 115.0886) / 2, "Z": (128.1307 + 129.1155) / 2,
                        "J": 113.1594})
WATER_MONO = 18.01056  # PEAKS uses the 5-decimal monoisotopic water here (1260/1260 exact)


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
    ``sequence`` plus ``index`` (accession -> position).  Selenocysteine (U) is stored as X,
    as in PEAKS's own database."""

    def __init__(self, entries, u_to_x=True):
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
            seq = seq.upper().replace("*", "")
            if u_to_x:
                seq = seq.replace("U", "X")
            self.sequence.append(seq)
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
    sum_k s_k / k with s_1 >= s_2 >= ..., evaluated in float32 exactly as
    ``t = (float)(t + (float)(s_k / k))`` (reproduces all 1336 reference values bit for bit;
    double precision differs in the last printed digit for ~18% of proteins)."""
    f32 = np.float32
    t = f32(0.0)
    for k, s in enumerate(sorted((float(x) for x in scores), reverse=True), start=1):
        t = f32(t + f32(f32(s) / f32(k)))
    return float(t)


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


# --------------------------------------------------------------------------
# Peptide mass as written in db.protein-peptides.csv
# --------------------------------------------------------------------------

#: monoisotopic residue masses with 5 decimals (the table that reproduces the "Mass" column of
#: db.protein-peptides.csv for all 3932 distinct reference peptides)
MONO_RESIDUE_5 = {
    "A": "71.03711", "C": "103.00919", "D": "115.02694", "E": "129.04259", "F": "147.06841",
    "G": "57.02146", "H": "137.05891", "I": "113.08406", "L": "113.08406", "K": "128.09496",
    "M": "131.04049", "N": "114.04293", "P": "97.05276", "Q": "128.05858", "R": "156.10111",
    "S": "87.03203", "T": "101.04768", "V": "99.06841", "W": "186.07931", "Y": "163.06333",
    "U": "150.95364", "O": "237.14773",
}
_MONO5 = None


def _exact_delta(delta):
    """Exact monoisotopic delta for a PEAKS 2-decimal delta (15.99 -> 15.994915)."""
    from .seqnorm import MOD_MASSES

    r = round(delta, 2)
    best = None
    for v in MOD_MASSES.values():
        if round(v, 2) == r:
            if best is None or abs(v - delta) < abs(best - delta):
                best = v
    return delta if best is None else best


def protein_peptide_mass(peptide):
    """"Mass" of db.protein-peptides.csv: exact decimal sum of 5-decimal residue masses,
    water 18.01056 and exact modification deltas, rounded HALF_DOWN to 4 decimals.

    (db.peptides.csv prints the float32 search mass instead, so the two files disagree in
    the 4th decimal for ~30% of peptides; this reproduces the protein-peptides value.)"""
    from decimal import ROUND_HALF_DOWN, Decimal

    global _MONO5
    if _MONO5 is None:
        _MONO5 = {k: Decimal(v) for k, v in MONO_RESIDUE_5.items()}
    stripped, mods = parse(peptide)
    total = Decimal("18.01056")
    for a in stripped:
        total += _MONO5.get(a, Decimal(0))
    for _pos, d in mods:
        total += Decimal(repr(_exact_delta(d)))
    return float(total.quantize(Decimal("0.0001"), rounding=ROUND_HALF_DOWN))


# --------------------------------------------------------------------------
# Grouping
# --------------------------------------------------------------------------


def group_proteins(prot_fams, dom_fams, order):
    """Greedy PEAKS grouping.

    ``prot_fams``: {protein: set of peptide families}; ``dom_fams``: {protein: frozenset of
    families with score >= T}; ``order``: proteins in processing order.

    Proteins are processed in ``order`` (PEAKS: most supporting peptide rows first).  An
    unassigned protein X opens a new group; every still unassigned protein Y whose
    families (score >= T) are a subset of X's joins X's group.  Y is Top when its set equals
    X's (mutual domination; e.g. ACTB/ACTG).  Returns (group_of, groups, top) with groups a
    list of member lists (leader first) numbered in creation order.
    """
    by_fam = defaultdict(list)
    for p, fs in prot_fams.items():
        for f in fs:
            by_fam[f].append(p)
    rank = {p: r for r, p in enumerate(order)}
    group_of, groups, top = {}, [], {}
    for x in order:
        if x in group_of:
            continue
        g = len(groups)
        group_of[x] = g
        members = [x]
        top[x] = True
        dx = dom_fams[x]
        cand = set()
        for f in prot_fams[x]:
            cand.update(by_fam[f])
        for y in sorted(cand, key=rank.__getitem__):
            if y in group_of:
                continue
            dy = dom_fams[y]
            if dy <= dx:
                group_of[y] = g
                members.append(y)
                top[y] = dy == dx
        groups.append(members)
    return group_of, groups, top


# --------------------------------------------------------------------------
# Main entry point
# --------------------------------------------------------------------------


def _col(df, *names):
    for n in names:
        if n in df.columns:
            return n
    return None


def _sample_from_columns(df):
    for c in df.columns:
        if c.startswith("Area ") and not c.endswith(" Mod"):
            return c[5:]
    return "Sample 1"


def _text(v):
    if v is None:
        return ""
    try:
        if v != v:
            return ""
    except (TypeError, ValueError):
        pass
    return str(v)


def _f32_sum(values):
    """Sum in float32 (PEAKS accumulates protein areas as Java floats)."""
    t = np.float32(0.0)
    for v in values:
        t = np.float32(t + np.float32(v))
    return float(t)


def _half_up_sci(x):
    """PEAKS prints the protein Area with 3 significant digits rounding HALF_UP (one exact
    tie in the references: ERF3A 1065000.0f -> "1.07E6"), whereas
    :func:`peaksparity.peaks_export.format_value` rounds HALF_EVEN.  For an exact tie
    return the next double up so that either formatter prints the PEAKS value."""
    if x != x or x == 0:
        return x
    from decimal import Decimal

    d = abs(Decimal(float(x)))
    m = d.scaleb(-d.adjusted()) * 100  # 3 significant digits -> integer part
    if m == m.to_integral_value() + Decimal("0.5"):
        return float(np.nextafter(float(x), np.inf if x > 0 else -np.inf))
    return x


def infer_proteins(peptides_df, fasta_path, score_threshold=15.0, min_unique=1, *,
                   psms_df=None, features_df=None, sample=None, enzyme=DEFAULT_ENZYME,
                   mapping="fasta", score_col=None, domination_threshold=None,
                   feature_min_psm_score=20.0, return_details=False, psms=None, sample_name=None):
    """PEAKS-style protein inference.

    Parameters
    ----------
    peptides_df : DataFrame in the ``db.peptides.csv`` layout (one row per peptidoform).
        Required: ``Peptide`` and ``-10LgP`` (or ``-10lgP``).  Used when present: ``Mass``,
        ``Length``, ``ppm``, ``m/z``, ``z``, ``RT``, ``Scan``, ``Source File``,
        ``Area <sample>``, ``#Feature``, ``PTM``, ``AScore``, ``Found By``, ``Accession``.
    fasta_path : FASTA path or :class:`ProteinDB` (target proteins; ``#DECOY#`` entries are
        skipped).
    score_threshold : protein -10lgP filter (PEAKS default 15); also the domination
        threshold T unless ``domination_threshold`` is given.
    min_unique : minimum #Unique per reported protein.
    psms_df : optional ``db.psms.csv``-style table (Peptide, -10LgP, Feature Id, Area) for
        #Spec (PSM rows per peptide) and for the feature-level Area.
    features_df : optional ``<sample>.db.features.csv`` (Feature Id, Area) with unrounded
        feature areas; without it the (rounded) PSM or peptide areas are summed.
    sample : sample name for the per-sample columns (default: taken from the peptides'
        ``Area <name>`` column, else "Sample 1").
    enzyme : :class:`Enzyme` deciding which occurrences count (default semi-specific
        trypsin, 2 missed cleavages); None = any occurrence.
    mapping : "fasta" (map by sequence) or "accession" (use the peptides' ``Accession``
        column, ``:``/``;`` separated, as the protein list; positions still from the FASTA).
    score_col : column with full-precision peptide scores (default: ``score`` if present,
        else the -10LgP column).  PEAKS's CSV rounds -10lgP to 2 decimals, which limits the
        protein -10lgP to about +-0.01.
    feature_min_psm_score : a feature counts for a peptide's protein-level Area / #Feature
        only if it carries a PSM with at least this -10lgP (PEAKS: 20, see
        docs/PROTEIN_INFERENCE.md).
    return_details : also return a dict with the unfiltered groups and per-protein data.
    psms, sample_name : aliases of ``psms_df`` and ``sample``.

    Returns
    -------
    (proteins_df, protein_peptides_df) with exactly the PEAKS columns, or
    (proteins_df, protein_peptides_df, details) with ``return_details``.
    """
    if psms_df is None:
        psms_df = psms
    if sample is None:
        sample = sample_name
    db = read_protein_db(fasta_path)
    T = score_threshold if domination_threshold is None else domination_threshold
    pep = peptides_df.reset_index(drop=True)
    if sample is None:
        sample = _sample_from_columns(pep)
    lgp_col = _col(pep, "-10LgP", "-10lgP")
    if score_col is None:
        score_col = "score" if "score" in pep.columns else lgp_col
    if score_col is None:
        raise ValueError("peptides_df needs a -10LgP (or score) column")
    peptides = pep["Peptide"].astype(str).tolist()
    scores = pep[score_col].astype(float).to_numpy()
    stripped = [parse(p)[0] for p in peptides]
    fams = [s.replace("I", "L") for s in stripped]
    n = len(pep)

    # ---- peptide -> protein occurrences ----------------------------------
    if mapping == "fasta":
        hits = map_peptides(db, set(stripped), enzyme)
    elif mapping == "accession":
        acc_col = _col(pep, "Accession")
        if acc_col is None:
            raise ValueError("mapping='accession' needs an Accession column")
        hits = {}
        for s, a in zip(stripped, pep[acc_col].astype(str)):
            if s in hits:
                continue
            lst = []
            for acc in re.split(r"[:;]", a):
                i = db.index.get(acc.strip())
                if i is None:
                    continue
                seq = db.sequence[i]
                st = -1
                pos = seq.find(s)
                while pos >= 0:
                    if enzyme is None or enzyme.valid(seq, pos, len(s)):
                        st = pos
                        break
                    pos = seq.find(s, pos + 1)
                if st < 0:
                    st = seq.find(s)
                exact = st >= 0
                if not exact:
                    st = seq.replace("I", "L").find(s.replace("I", "L"))
                if st >= 0:
                    lst.append((i, st, exact))
            hits[s] = sorted(lst)
    else:
        raise ValueError("mapping must be 'fasta' or 'accession'")

    fam_best = defaultdict(float)
    for f, sc in zip(fams, scores):
        fam_best[f] = max(fam_best[f], sc)

    # protein -> {row: start}; families; exact families
    prot_rows = defaultdict(dict)
    prot_fams = defaultdict(set)
    il_only = defaultdict(dict)  # protein -> {fam: start} for families matched only by I/L
    for r, s in enumerate(stripped):
        for i, st, exact in hits.get(s, ()):
            prot_fams[i].add(fams[r])
            if exact:
                prot_rows[i][r] = st
            else:
                il_only[i].setdefault(fams[r], []).append((r, st))
    for i, d in il_only.items():
        have = {fams[r] for r in prot_rows[i]}
        for f, lst in d.items():
            if f not in have:  # no verbatim variant among the input peptides
                for r, st in lst:
                    prot_rows[i][r] = st

    proteins = sorted(prot_fams)
    pscore, nrows = {}, {}
    for i in proteins:
        best = {}
        for r in prot_rows[i]:
            best[stripped[r]] = max(best.get(stripped[r], -np.inf), scores[r])
        pscore[i] = protein_score(best.values())
        nrows[i] = len(prot_rows[i])
    dom = {i: frozenset(f for f in prot_fams[i] if fam_best[f] >= T) for i in proteins}
    order = sorted(proteins, key=lambda i: (-nrows[i], -pscore[i], db.accession[i]))
    group_of, groups, top = group_proteins(prot_fams, dom, order)

    fam_groups = defaultdict(set)
    for i in proteins:
        for f in prot_fams[i]:
            fam_groups[f].add(group_of[i])
    unique_fam = {f: len(g) == 1 for f, g in fam_groups.items()}

    # ---- PSM counts and features ------------------------------------------
    n_spec = np.ones(n)
    pep_features = [None] * n
    feat_area = {}
    if psms_df is not None and len(psms_df):
        cnt = psms_df.groupby("Peptide").size()
        n_spec = pep["Peptide"].map(cnt).fillna(0).to_numpy(float)
        if "Feature Id" in psms_df.columns:
            psc = _col(psms_df, "-10LgP", "-10lgP")
            ok = psms_df["Feature Id"].fillna(0).astype(int) > 0
            if psc is not None and feature_min_psm_score is not None:
                ok &= psms_df[psc].astype(float) >= feature_min_psm_score
            sub = psms_df[ok]
            fmap = sub.groupby("Peptide")["Feature Id"].agg(lambda v: sorted(set(int(x) for x in v)))
            pep_features = [fmap.get(p, []) for p in peptides]
            if features_df is not None and len(features_df):
                feat_area = dict(zip(features_df["Feature Id"].astype(int),
                                     features_df["Area"].astype(float)))
            if "Area" in sub.columns:
                for fid, a in zip(sub["Feature Id"].astype(int), sub["Area"].astype(float)):
                    if fid not in feat_area and a == a:
                        feat_area[fid] = a
    area_col = _col(pep, f"Area {sample}", "Area")

    def pep_area(r):
        """(area, n_features) for protein-peptides."""
        fs = pep_features[r]
        if fs is None:
            a = pep[area_col].iloc[r] if area_col else np.nan
            nf = pep["#Feature"].iloc[r] if "#Feature" in pep.columns else (0 if a != a else 1)
            return a, nf
        if not fs:
            return np.nan, 0
        return _f32_sum(feat_area.get(f, 0.0) for f in fs), len(fs)

    # ---- per-protein summary -----------------------------------------------
    ptm_col = _col(pep, "PTM")
    rows_out, pp_out = [], []
    kept = {}
    for i in proteins:
        seq = db.sequence[i]
        rows = prot_rows[i]
        se = {(st + 1, st + len(stripped[r])) for r, st in rows.items()}
        uniq_rows = [r for r in rows if unique_fam[fams[r]]]
        use = {(rows[r] + 1, rows[r] + len(stripped[r])) for r in uniq_rows}
        n_unique = len(use)
        if pscore[i] < score_threshold or n_unique < min_unique or not rows:
            continue
        cov = np.zeros(len(seq), bool)
        for s0, e0 in se:
            cov[s0 - 1:e0] = True
        coverage = 100.0 * cov.mean() if len(seq) else 0.0
        spec = float(sum(n_spec[r] for r in rows))
        # Area: distinct features of the unique supporting peptides (float32 sum)
        if all(pep_features[r] is not None for r in uniq_rows):
            fs = sorted({f for r in uniq_rows for f in pep_features[r]})
            area = _f32_sum(feat_area.get(f, 0.0) for f in fs) if fs else np.nan
        else:
            vals = [pep_area(r)[0] for r in uniq_rows]
            vals = [v for v in vals if v == v]
            area = _f32_sum(vals) if vals else np.nan
        area = _half_up_sci(area)
        ptms = set()
        if ptm_col:
            for r in rows:
                for t in _text(pep[ptm_col].iloc[r]).split(";"):
                    if t.strip():
                        ptms.add(t.strip())
        kept[i] = True
        rows_out.append({
            "Protein Group": group_of[i] + 1,
            "Top": bool(top[i]),
            "Accession": db.accession[i],
            "Gene": db.gene[i],
            "-10LgP": pscore[i],
            "Coverage(%)": coverage,
            f"Coverage(%) {sample}": coverage,
            "#Peptides": len(se),
            f"Area {sample}": area,
            "#Unique": n_unique,
            f"#Spec {sample}": spec,
            "PTM": "; ".join(sorted(ptms)),
            "Average Mass": average_mass(seq),
            "Description": db.description[i],
        })
        for r, st in sorted(rows.items(), key=lambda kv: (kv[1], peptides[kv[0]])):
            a, nf = pep_area(r)
            L = len(stripped[r])

            def g(c, default=""):
                cc = _col(pep, c)
                return pep[cc].iloc[r] if cc else default

            pp_out.append({
                "Protein Group": group_of[i] + 1,
                "Accession": db.accession[i],
                "Gene": db.gene[i],
                "Top": bool(top[i]),
                "Peptide": peptides[r],
                "Extended Peptide": extended_peptide(peptides[r], seq, st, L),
                "Unique": bool(unique_fam[fams[r]]),
                "-10LgP": float(pep[lgp_col].iloc[r]) if lgp_col else float(scores[r]),
                "Mass": protein_peptide_mass(peptides[r]),
                "Length": L,
                "ppm": g("ppm", np.nan),
                "m/z": g("m/z", np.nan),
                "z": g("z", pd.NA),
                "RT": g("RT", np.nan),
                "Scan": g("Scan", pd.NA),
                "Source File": _text(g("Source File")),
                f"Area {sample}": a,
                "#Feature": nf,
                "Start": st + 1,
                "End": st + L,
                "PTM": _text(g("PTM")),
                "AScore": _text(g("AScore")),
                "Found By": _text(g("Found By", "DB Search")),
            })

    pcols = canonical_columns("proteins", (sample,))
    ppcols = canonical_columns("protein_peptides", (sample,))
    proteins_df = pd.DataFrame(rows_out, columns=pcols)
    if len(proteins_df):
        proteins_df = proteins_df.sort_values(
            ["Protein Group", "Top", "-10LgP", "Accession"],
            ascending=[True, False, False, True], kind="mergesort").reset_index(drop=True)
    prot_rank = {a: k for k, a in enumerate(proteins_df["Accession"])}
    pp_df = pd.DataFrame(pp_out, columns=ppcols)
    if len(pp_df):
        pp_df["_r"] = pp_df["Accession"].map(prot_rank)
        pp_df = pp_df.sort_values(["_r", "Start", "End", "Peptide"], kind="mergesort") \
                     .drop(columns="_r").reset_index(drop=True)
    for c in ("Protein Group", "#Peptides", "#Unique"):
        proteins_df[c] = proteins_df[c].astype("Int64")
    for c in ("Protein Group", "Length", "Start", "End", "#Feature", "z", "Scan"):
        pp_df[c] = pd.to_numeric(pp_df[c], errors="coerce").astype("Int64")
    proteins_df["Top"] = proteins_df["Top"].astype(bool)
    pp_df["Top"] = pp_df["Top"].astype(bool)
    pp_df["Unique"] = pp_df["Unique"].astype(bool)
    if not return_details:
        return proteins_df, pp_df
    details = {
        "db": db,
        "groups": [[db.accession[i] for i in mem] for mem in groups],
        "group_of": {db.accession[i]: g + 1 for i, g in group_of.items()},
        "top": {db.accession[i]: t for i, t in top.items()},
        "score": {db.accession[i]: s for i, s in pscore.items()},
        "n_rows": {db.accession[i]: v for i, v in nrows.items()},
        "unique_family": unique_fam,
        "order": [db.accession[i] for i in order],
    }
    return proteins_df, pp_df, details
