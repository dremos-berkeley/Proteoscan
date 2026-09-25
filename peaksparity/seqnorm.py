"""Canonical representation of modified peptide sequences.

Canonical form = PEAKS style: residues with mass-delta annotations rounded to
two decimals, always signed and with a leading zero::

    PEPTM(+15.99)IDE          residue modification
    (+42.01)PEPTIDE           N-terminal modification
    PEPTIDE-(-0.98)           C-terminal modification (no PEAKS example seen;
                              this is our own convention)

Accepted inputs:

* PEAKS:      ``N(+.98)``, ``M(+15.99)``, ``C(+57.02)``, ``(+42.01)PEPTIDE``
* ProForma:   ``PEPTM[Oxidation]IDE``, ``PEPTM[+15.9949]IDE``,
              ``PEPTM[UNIMOD:35]IDE``, ``[Acetyl]-PEPTIDE``, ``PEPTIDE-[Amidated]``
* OpenMS:     ``PEPTM(Oxidation)IDE``, ``.(Acetyl)PEPTIDE``, ``PEPTIDE.(Amidated)``
* Comet/TPP:  ``n[42.0106]PEPTIDE`` (unsigned numbers are taken as deltas)
* MaxQuant:   ``_PEPTM(Oxidation (M))IDE_``
* flanked:    ``K.PEPTIDE.R`` (PEAKS "Extended Peptide"; flanks are dropped)

Several modifications on one site are summed into one delta.  Residue letters
are kept verbatim: I and L are distinct and the site of each modification is
part of the identity.  Deltas are rounded HALF_UP from the shortest decimal
representation of the value (15.9949 -> +15.99, 0.984016 -> +0.98,
-17.026549 -> -17.03).
"""

from __future__ import annotations

import re
from decimal import ROUND_HALF_UP, Decimal
from functools import lru_cache

__all__ = [
    "MOD_MASSES",
    "UNIMOD_IDS",
    "canonical",
    "parse",
    "strip_mods",
    "mods_of",
    "equal_modseq",
    "format_delta",
    "mod_delta",
    "il_collapse",
]

#: monoisotopic mass deltas by (lower-case) name / alias
MOD_MASSES = {}
#: Unimod accession number -> delta
UNIMOD_IDS = {}


def _add(delta, unimod_id, *names):
    for n in names:
        MOD_MASSES[n.lower()] = delta
    if unimod_id is not None:
        UNIMOD_IDS[unimod_id] = delta


_add(57.021464, 4, "Carbamidomethyl", "Carbamidomethylation", "CAM", "Iodoacetamide derivative")
_add(15.994915, 35, "Oxidation", "Ox", "Hydroxylation")
_add(0.984016, 7, "Deamidated", "Deamidation", "Deamidation (NQ)", "Deam", "de")
_add(42.010565, 1, "Acetyl", "Acetylation", "ac", "Acetylation (N-term)", "Acetylation (Protein N-term)",
     "Acetylation (K)")
_add(71.037114, 24, "Propionamide", "Propionamidation", "Acrylamide adduct")
_add(-17.026549, 28, "Gln->pyro-Glu", "Pyro-glu from Q", "Pyro-glu (Q)", "Pyroglu from Q", "Pyro_glu_Q",
     "Glu(pyro)Q")
_add(-18.010565, 27, "Glu->pyro-Glu", "Pyro-glu from E", "Pyro-glu (E)", "Pyroglu from E", "Pyro_glu_E")
_add(43.005814, 5, "Carbamyl", "Carbamylation")
_add(79.966331, 21, "Phospho", "Phosphorylation", "ph")
_add(14.01565, 34, "Methyl", "Methylation")
_add(28.0313, 36, "Dimethyl", "Dimethylation")
_add(42.04695, 37, "Trimethyl", "Trimethylation")
_add(-0.984016, 2, "Amidated", "Amidation")
_add(-17.026549, 385, "Ammonia-loss")
_add(-18.010565, 23, "Dehydrated", "Dehydration", "Water loss")
_add(114.042927, 121, "GG", "GlyGly", "Ubiquitination", "GlyGly (K)")
_add(27.994915, 122, "Formyl", "Formylation")
_add(229.162932, 737, "TMT6plex", "TMT", "TMT10plex", "TMT11plex")
_add(304.207146, 2016, "TMTpro", "TMTpro16plex")
_add(-131.040485, 765, "Met-loss")
_add(-89.02998, 766, "Met-loss+Acetyl")
_add(39.994915, 26, "Pyro-carbamidomethyl")
_add(203.079373, 43, "HexNAc")
_add(31.989829, 425, "Dioxidation")
_add(21.981943, 30, "Cation:Na", "Sodium")

_NUM_RE = re.compile(r"^([+-]?)(\d+\.?\d*|\.\d+)$")
_Q = Decimal("0.01")


def format_delta(delta):
    """Round a mass delta to 2 decimals, PEAKS style: ``+15.99``, ``-17.03``, ``+0.98``."""
    d = delta if isinstance(delta, Decimal) else Decimal(repr(float(delta)))
    q = d.quantize(_Q, rounding=ROUND_HALF_UP)
    if q == 0:
        q = Decimal("0.00")
    s = format(q, "f")
    return s if s.startswith("-") else "+" + s


@lru_cache(maxsize=4096)
def mod_delta(text):
    """Mass delta (float) of a modification given by number, Unimod name/alias or UNIMOD:n."""
    t = text.strip()
    m = _NUM_RE.match(t)
    if m:
        v = float(m.group(2))
        return -v if m.group(1) == "-" else v
    low = t.lower()
    for pre in ("u:", "unimod:"):
        if low.startswith(pre):
            rest = t[len(pre):].strip()
            if rest.isdigit():
                if int(rest) in UNIMOD_IDS:
                    return UNIMOD_IDS[int(rest)]
                raise ValueError(f"unknown Unimod accession {text!r}")
            low = rest.lower()
            break
    if low in MOD_MASSES:
        return MOD_MASSES[low]
    base = re.sub(r"\s*[(@].*$", "", low).strip()  # "Oxidation (M)", "Oxidation@M"
    if base in MOD_MASSES:
        return MOD_MASSES[base]
    raise ValueError(f"unknown modification {text!r}")


def _read_group(s, i):
    """Read a balanced (...) or [...] group starting at s[i]; return (content, next index)."""
    open_ch = s[i]
    close_ch = ")" if open_ch == "(" else "]"
    depth = 0
    for j in range(i, len(s)):
        if s[j] == open_ch:
            depth += 1
        elif s[j] == close_ch:
            depth -= 1
            if depth == 0:
                return s[i + 1: j], j + 1
    raise ValueError(f"unbalanced {open_ch!r} in {s!r}")


_FLANK_RE = re.compile(r"^(?:[A-Z]|-)\.(.+)\.(?:[A-Z]|-)$")


@lru_cache(maxsize=200000)
def parse(seq):
    """Parse a modified peptide.

    Returns ``(stripped, mods)`` with ``mods`` a tuple of ``(position, delta)``
    sorted by position, where position 0 is the N-terminus, 1..n the residues
    and n+1 the C-terminus; deltas on the same position are summed.
    """
    s = seq.strip()
    if len(s) >= 2 and s[0] == "_" and s[-1] == "_":  # MaxQuant
        s = s[1:-1]
    m = _FLANK_RE.match(s)
    if m and not s.startswith(".(") and not s.endswith(")"):
        s = m.group(1)
    residues = []
    mods = {}

    def add(pos, content):
        mods[pos] = mods.get(pos, 0.0) + mod_delta(content)

    i, n = 0, len(s)
    # ---- N-terminal prefix
    if s.startswith("n[") or s.startswith("n("):  # Comet / TPP
        content, i = _read_group(s, 1)
        add(0, content)
    elif s.startswith(".(") or s.startswith(".["):  # OpenMS
        content, i = _read_group(s, 1)
        add(0, content)
    else:
        # ProForma "[X]-" (possibly repeated) or PEAKS leading "(X)"
        while i < n and s[i] in "[(":
            content, j = _read_group(s, i)
            add(0, content)
            i = j
            if i < n and s[i] == "-":
                i += 1
    # ---- residues
    while i < n:
        c = s[i]
        if "A" <= c <= "Z":
            residues.append(c)
            i += 1
            while i < n and s[i] in "[(":
                content, i = _read_group(s, i)
                add(len(residues), content)
            continue
        if c in "-." and i + 1 < n and s[i + 1] in "[(":  # C-term: "-[X]", "-(X)", ".(X)"
            content, i = _read_group(s, i + 1)
            add(len(residues) + 1, content)
            continue
        if c == "c" and i + 1 < n and s[i + 1] in "[(":  # Comet C-term
            content, i = _read_group(s, i + 1)
            add(len(residues) + 1, content)
            continue
        raise ValueError(f"cannot parse {seq!r} at position {i} ({c!r})")
    if not residues:
        raise ValueError(f"no residues in {seq!r}")
    return "".join(residues), tuple(sorted(mods.items()))


def _render(stripped, mods):
    n = len(stripped)
    by = {p: format_delta(d) for p, d in mods}
    out = []
    if 0 in by:
        out.append(f"({by[0]})")
    for k, aa in enumerate(stripped, start=1):
        out.append(aa)
        if k in by:
            out.append(f"({by[k]})")
    if n + 1 in by:
        out.append(f"-({by[n + 1]})")
    return "".join(out)


@lru_cache(maxsize=200000)
def _canonical(seq, fixed_items, nterm_to_first_residue):
    stripped, mods = parse(seq)
    md = dict(mods)
    if nterm_to_first_residue and 0 in md:
        md[1] = md.get(1, 0.0) + md.pop(0)
    for aa, delta in fixed_items:
        for k, r in enumerate(stripped, start=1):
            if r == aa and k not in md:
                md[k] = delta
    return _render(stripped, tuple(sorted(md.items())))


def canonical(seq, fixed=None, nterm_to_first_residue=False):
    """Canonical PEAKS-style modified sequence, e.g. ``"PEPTM[Oxidation]IDE" -> "PEPTM(+15.99)IDE"``.

    ``fixed``: optional {residue: delta or mod name} of fixed modifications to
    add to every unmodified matching residue (for tools that omit them).
    ``nterm_to_first_residue``: fold N-terminal mods into the first residue
    (off by default; PEAKS distinguishes them).
    """
    items = ()
    if fixed:
        items = tuple(sorted(
            (aa, float(v) if isinstance(v, (int, float)) else mod_delta(str(v))) for aa, v in fixed.items()
        ))
    return _canonical(seq, items, bool(nterm_to_first_residue))


def strip_mods(seq):
    """Plain residue sequence (no modifications, no flanks)."""
    return parse(seq)[0]


def mods_of(seq):
    """List of ``(position, site, delta_str)``; site is the residue letter,
    ``"N-term"`` (position 0) or ``"C-term"`` (position n+1)."""
    stripped, mods = parse(seq)
    n = len(stripped)
    out = []
    for p, d in mods:
        site = "N-term" if p == 0 else ("C-term" if p == n + 1 else stripped[p - 1])
        out.append((p, site, format_delta(d)))
    return out


def equal_modseq(a, b, **kw):
    """True if both strings denote exactly the same modified peptide
    (same residues, I != L, same modification sites and rounded deltas)."""
    return canonical(a, **kw) == canonical(b, **kw)


def il_collapse(seq):
    """Replace I by L (for I/L-ambiguity diagnostics only; not used for matching)."""
    return seq.replace("I", "L")
