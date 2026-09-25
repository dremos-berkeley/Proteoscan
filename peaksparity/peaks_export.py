"""Reader and byte-faithful writer for PEAKS Studio DB-search export directories.

A PEAKS "Analysis N - Sample k" export directory contains (among HTML/figures):

    db.psms.csv                  PSMs of accepted peptides
    db.peptides.csv              one row per accepted (modified) peptide
    db.proteins.csv              one row per protein (grouped by "Protein Group")
    db.protein-peptides.csv      protein x peptide mapping
    <Sample>.db.features.csv     all LC-MS features of the sample
    <Sample>.db.denovoOnly.csv   de novo only spectra (ALC above threshold)
    db.ptmProfile.csv            PTM summary
    db.modifiedPeptide.csv       modified peptide vs. unmodified counterpart
    db.summary.html              search parameters and filter thresholds

``read_export(dir)`` returns a :class:`PeaksExport` holding one pandas DataFrame
per CSV (``None`` when the file is absent) plus ``params`` parsed from
``db.summary.html``.  ``write_export(export, dir)`` writes the CSVs back with
the same column order, per-column quoting and number formatting as PEAKS, so
that read -> write reproduces the reference files byte for byte.

PEAKS writes CSVs with a hand-rolled writer rather than a CSV library, so the
quoting is a property of the *column*, not of the value (see
docs/PEAKS_EXPORT_FORMAT.md).  Three quoting modes occur:

    always    value is always quoted, including empty values ("")
    never     value is never quoted (we still quote if the value contains
              a comma, quote or newline, so output stays valid CSV)
    nonempty  quoted when non-empty, bare when empty (only AScore)

and the following number formats (``kind``):

    int       plain integer
    fixedN    fixed N decimals (Java ``String.format("%.Nf")``: HALF_UP on the
              shortest decimal representation)
    sci       Java ``DecimalFormat("0.##E0")`` e.g. 6.58E9, 1.1E5, 3E7
    jfloat    Java ``Float.toString`` (JDK <= 18 FloatingDecimal) e.g.
              129882.984, 3.0000496E7, 4531940.0, 2.95233664E8
    dec2      Java ``DecimalFormat("#.##")`` e.g. .98, 57.02
    bool      true / false
    str       text

Missing numbers are written as an empty field, except "MS2 Correlation"
which uses ``NaN``.  Missing text is represented as ``""`` (empty string) in
the DataFrames, never NaN.

No dataclasses are used; tables are plain DataFrames.
"""

from __future__ import annotations

import glob
import html
import os
import re
import struct
from decimal import ROUND_HALF_EVEN, ROUND_HALF_UP, Decimal

import numpy as np
import pandas as pd

__all__ = [
    "TABLES",
    "PeaksExport",
    "TableLayout",
    "read_export",
    "write_export",
    "read_table",
    "write_table",
    "parse_summary_html",
    "canonical_columns",
    "empty_table",
    "format_value",
    "java_float_to_string",
]

# --------------------------------------------------------------------------
# File names
# --------------------------------------------------------------------------

#: logical table name -> file name template ({sample} = sample name)
TABLES = {
    "psms": "db.psms.csv",
    "peptides": "db.peptides.csv",
    "proteins": "db.proteins.csv",
    "protein_peptides": "db.protein-peptides.csv",
    "features": "{sample}.db.features.csv",
    "denovo_only": "{sample}.db.denovoOnly.csv",
    "ptm_profile": "db.ptmProfile.csv",
    "modified_peptides": "db.modifiedPeptide.csv",
}
SUMMARY_FILE = "db.summary.html"
_PER_SAMPLE = {"features", "denovo_only"}

# --------------------------------------------------------------------------
# Column specifications (the PEAKS layout observed in PEAKS Studio 12/13
# DB-search exports).  "{S}" is a sample-name placeholder ("Sample 1").
# Entries: (name, kind, quote, extras) with extras keys
#   na      : text written for missing values (default "")
#   hq      : header cell quoted (default True)
#   default : column name used by canonical_columns() for pattern names
# --------------------------------------------------------------------------


def _c(name, kind, quote="never", **extras):
    return (name, kind, quote, extras)


_TAG = "tag(>={S}%)"  # tag(>=0.0%) in psms/peptides, tag(>=0%) in denovoOnly

_SPECS = {
    "psms": [
        _c("Peptide", "str"),
        _c("-10LgP", "fixed2"),
        _c("Mass", "fixed4"),
        _c("Tag Length", "int"),
        _c("CAA (%)", "fixed1"),
        _c("Length", "int"),
        _c("Delta RT", "fixed4"),
        _c("MS2 Correlation", "fixed2", na="NaN"),
        _c("ppm", "fixed1"),
        _c("m/z", "fixed5"),
        _c("raw m/z", "fixed5"),
        _c("z", "int"),
        _c("RT", "fixed4"),
        _c("raw RT", "fixed4"),
        _c("Area", "sci"),
        _c("Scan", "int"),
        _c("Source File", "str", "always"),
        _c("Feature Id", "int"),
        _c("Accession", "str", "always"),
        _c("PTM", "str", "always"),
        _c("AScore", "str", "nonempty"),
        _c("Positional Confidence", "str"),
        _c(_TAG, "str", default="tag(>=0.0%)"),
    ],
    "peptides": [
        _c("Peptide", "str"),
        _c("-10LgP", "fixed2"),
        _c("Mass", "fixed4"),
        _c("Tag Length", "int"),
        _c("CAA (%)", "fixed1"),
        _c("Length", "int"),
        _c("Delta RT", "fixed4"),
        _c("MS2 Correlation", "fixed2", na="NaN"),
        _c("ppm", "fixed1"),
        _c("m/z", "fixed5"),
        _c("z", "int"),
        _c("RT", "fixed4"),
        _c("Area {S}", "sci"),
        _c("Scan", "int"),
        _c("Source File", "str"),
        _c("#Feature", "int"),
        _c("#Feature {S}", "int"),
        _c("Accession", "str"),
        _c("PTM", "str", "always"),
        _c("AScore", "str", "nonempty"),
        _c("Positional Confidence", "str"),
        _c(_TAG, "str", default="tag(>=0.0%)"),
        _c("Found By", "str"),
    ],
    "proteins": [
        _c("Protein Group", "int"),
        _c("Top", "bool"),
        _c("Accession", "str", "always"),
        _c("Gene", "str", "always"),
        _c("-10LgP", "jfloat"),
        _c("Coverage(%)", "fixed2"),
        _c("Coverage(%) {S}", "fixed2"),
        _c("#Peptides", "int"),
        _c("Area {S}", "sci"),
        _c("#Unique", "int"),
        _c("#Spec {S}", "jfloat"),
        _c("PTM", "str", "always"),
        _c("Average Mass", "jfloat"),
        _c("Description", "str", "always"),
    ],
    "protein_peptides": [
        _c("Protein Group", "int"),
        _c("Accession", "str", "always"),
        _c("Gene", "str", "always"),
        _c("Top", "bool"),
        _c("Peptide", "str"),
        _c("Extended Peptide", "str"),
        _c("Unique", "bool"),
        _c("-10LgP", "fixed2"),
        _c("Mass", "fixed4"),
        _c("Length", "int"),
        _c("ppm", "fixed1"),
        _c("m/z", "fixed5"),
        _c("z", "int"),
        _c("RT", "fixed4"),
        _c("Scan", "int"),
        _c("Source File", "str"),
        _c("Area {S}", "sci"),
        _c("#Feature", "int"),
        _c("Start", "int"),
        _c("End", "int"),
        _c("PTM", "str", "always"),
        _c("AScore", "str", "nonempty"),
        _c("Found By", "str"),
    ],
    "features": [
        _c("Source File", "str", "always"),
        _c("Feature Id", "int"),
        _c("m/z", "fixed5"),
        _c("z", "int"),
        _c("RT", "fixed4"),
        _c("RT Begin", "fixed4"),
        _c("RT End", "fixed4"),
        _c("Area", "jfloat"),
        _c("Quality", "fixed2"),
        _c("DB Peptide", "str"),
        _c("-10LgP", "fixed2"),
        _c("Denovo Peptide", "str"),
        _c("ALC (%)", "fixed1"),
        _c("ppm", "fixed1"),
        _c("Accession", "str"),
    ],
    "denovo_only": [
        _c("Source File", "str", "always"),
        _c("Scan", "int"),
        _c("Peptide", "str"),
        _c("Tag length", "int"),
        _c("ALC (%)", "fixed1"),
        _c("Length", "int"),
        _c("m/z", "fixed5"),
        _c("z", "int"),
        _c("RT", "fixed4"),
        _c("Area", "jfloat"),
        _c("Mass", "fixed4"),
        _c("ppm", "fixed1"),
        _c("PTM", "str", "always"),
        _c("local confidence (%)", "str"),
        _c("mode", "str"),
        _c(_TAG, "str", hq=False, default="tag(>=0%)"),  # header NOT quoted
        _c("Accession", "str"),
        _c("Feature Id", "int"),
    ],
    "ptm_profile": [
        _c("Name", "str", "always"),
        _c("Mass (Da)", "dec2"),
        _c("Position", "str", "always"),
        _c("# peptides", "int"),
        _c("# proteins", "int"),
        _c("# confident modification sites", "int"),
        _c("# modification sites", "int"),
    ],
    "modified_peptides": [
        _c("Peptide Mod", "str"),
        _c("-10LgP Mod", "fixed2"),
        _c("Mass Mod", "fixed4"),
        _c("Length Mod", "int"),
        _c("PPM Mod", "fixed1"),
        _c("m/z Mod", "fixed5"),
        _c("z Mod", "int"),
        _c("RT Mod", "fixed4"),
        _c("Scan Mod", "int"),
        _c("Source File Mod", "str"),
        _c("Area {S} Mod", "sci"),
        _c("#Feature {S} Mod", "int"),
        _c("Found By Mod", "str"),
        _c("PTM", "str", "always"),
        _c("AScore", "str", "nonempty"),
        _c("Accession", "str"),
        _c("RT shift", "fixed4"),
        _c("Peptide", "str"),
        _c("-10LgP", "fixed2"),
        _c("Mass", "fixed4"),
        _c("Length", "int"),
        _c("PPM", "fixed1"),
        _c("m/z", "fixed5"),
        _c("z", "int"),
        _c("RT", "fixed4"),
        _c("Scan", "int"),
        _c("Source File", "str"),
        _c("Area {S}", "sci"),
        _c("#Feature {S}", "int"),
        _c("Found By", "str"),
    ],
}

_FLOAT_KINDS = {"fixed1", "fixed2", "fixed3", "fixed4", "fixed5", "fixed6", "sci", "jfloat", "dec2"}
_KINDS = _FLOAT_KINDS | {"int", "bool", "str"}
_QUOTES = ("always", "never", "nonempty")


def _pattern_regex(name):
    return re.compile(re.escape(name).replace(re.escape("{S}"), "(.+)"))


_SPEC_LOOKUP = {}
for _t, _entries in _SPECS.items():
    _lit = {}
    _pat = []
    for _e in _entries:
        if "{S}" in _e[0]:
            _pat.append((_pattern_regex(_e[0]), _e))
        else:
            _lit[_e[0]] = _e
    _pat.sort(key=lambda p: -len(p[1][0]))  # most specific ("Area {S} Mod") first
    _SPEC_LOOKUP[_t] = (_lit, _pat)


def _spec_for(table, column):
    """Return (kind, quote, na, header_quoted) for a column, or None if unknown."""
    if table not in _SPEC_LOOKUP:
        return None
    lit, pat = _SPEC_LOOKUP[table]
    e = lit.get(column)
    if e is None:
        for rx, entry in pat:
            if rx.fullmatch(column):
                e = entry
                break
    if e is None:
        return None
    name, kind, quote, ex = e
    return kind, quote, ex.get("na", ""), ex.get("hq", True)


def canonical_columns(table, samples=("Sample 1",)):
    """PEAKS column order for ``table``; per-sample columns are expanded for each sample."""
    if isinstance(samples, str):
        samples = (samples,)
    cols = []
    for name, _kind, _q, ex in _SPECS[table]:
        if "default" in ex:
            cols.append(ex["default"])
        elif "{S}" in name:
            cols.extend(name.replace("{S}", s) for s in samples)
        else:
            cols.append(name)
    return cols


def _dtype_for(kind):
    if kind in _FLOAT_KINDS:
        return "float64"
    if kind == "int":
        return "Int64"
    if kind == "bool":
        return "bool"
    return object


def empty_table(table, samples=("Sample 1",)):
    """An empty DataFrame with the PEAKS columns (in order) and matching dtypes."""
    cols = canonical_columns(table, samples)
    return pd.DataFrame({c: pd.Series([], dtype=_dtype_for(_spec_for(table, c)[0])) for c in cols})


# --------------------------------------------------------------------------
# Number formatting (Java emulation)
# --------------------------------------------------------------------------

# FloatingDecimal.insignificantDigitsNumber (JDK 8..18)
_INSIG = [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8,
          9, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15,
          15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19]


def java_float_to_string(x):
    """Emulate ``Float.toString((float) x)`` of JDK <= 18.

    Integer-valued floats go through FloatingDecimal's exact-integer path, which
    prints all integer digits (e.g. 2.95233664E8, not the shortest 2.9523366E8);
    other values use the shortest repr that round-trips to float32 (verified
    against ~130k values from PEAKS exports).
    """
    f = np.float32(x)
    if np.isnan(f):
        return "NaN"
    if np.isinf(f):
        return "Infinity" if f > 0 else "-Infinity"
    if f == 0:
        return "-0.0" if np.signbit(f) else "0.0"
    neg = bool(f < 0)
    a = abs(float(f))
    bits = struct.unpack(">I", struct.pack(">f", a))[0]
    exp = (bits >> 23) & 0xFF
    if exp != 0 and a.is_integer() and exp - 127 <= 62:
        bin_exp = exp - 127
        p2 = bin_exp - 24 - 1
        insig = _INSIG[p2] if (bin_exp > 24 and 1 < p2 < len(_INSIG)) else 0
        v = int(a)
        dec_shift = 0
        if insig:
            p = 10 ** insig
            r = v % p
            v //= p
            dec_shift = insig
            if r >= p // 2:
                v += 1
        sv = str(v)
        digits = sv.rstrip("0") or "0"
        dec_exp = len(sv) + dec_shift
    else:
        s = np.format_float_scientific(np.float32(a), unique=True, trim="-")
        mant, e = s.split("e")
        digits = mant.replace(".", "").rstrip("0") or "0"
        dec_exp = int(e) + 1
    if 0 < dec_exp < 8:
        ip = digits[:dec_exp].ljust(dec_exp, "0")
        fp = digits[dec_exp:] or "0"
        out = ip + "." + fp
    elif -3 < dec_exp <= 0:
        out = "0." + "0" * (-dec_exp) + digits
    else:
        out = digits[0] + "." + (digits[1:] or "0") + "E" + str(dec_exp - 1)
    return ("-" if neg else "") + out


def _fmt_fixed(x, n):
    """Java String.format("%.nf"): HALF_UP applied to the shortest decimal repr."""
    s = f"{x:.{n}f}"
    r = repr(float(x))
    if "e" in r or "inf" in r or "nan" in r:
        if r in ("nan", "inf", "-inf"):
            return {"nan": "NaN", "inf": "Infinity", "-inf": "-Infinity"}[r]
        return str(Decimal(r).quantize(Decimal(1).scaleb(-n), rounding=ROUND_HALF_UP))
    dot = r.find(".")
    if dot >= 0 and len(r) - dot - 1 == n + 1 and r[-1] == "5":
        return str(Decimal(r).quantize(Decimal(1).scaleb(-n), rounding=ROUND_HALF_UP))
    return s


_Q2 = Decimal("0.01")


def _fmt_sci(x):
    """Java DecimalFormat("0.##E0") (HALF_EVEN on the exact binary value)."""
    x = float(x)
    if x == 0:
        return "0E0"
    d = Decimal(x)
    sign = "-" if d < 0 else ""
    d = abs(d)
    e = d.adjusted()
    m = d.scaleb(-e).quantize(_Q2, rounding=ROUND_HALF_EVEN)
    if m >= 10:
        e += 1
        m = d.scaleb(-e).quantize(_Q2, rounding=ROUND_HALF_EVEN)
    ms = format(m, "f").rstrip("0").rstrip(".")
    return f"{sign}{ms}E{e}"


def _fmt_dec2(x):
    """Java DecimalFormat("#.##"): no leading zero, no trailing zeros."""
    q = Decimal(float(x)).quantize(_Q2, rounding=ROUND_HALF_EVEN)
    s = format(q, "f").rstrip("0").rstrip(".")
    neg = s.startswith("-")
    if neg:
        s = s[1:]
    if s.startswith("0."):
        s = s[1:]
    if s == "":
        s = "0"
    return ("-" if neg and s != "0" else "") + s


def format_value(kind, v):
    """Format one non-missing value of the given kind the way PEAKS does."""
    if kind == "str":
        return str(v)
    if kind == "int":
        return str(int(v))
    if kind == "bool":
        if isinstance(v, str):
            return v.lower()
        return "true" if bool(v) else "false"
    if kind.startswith("fixed"):
        return _fmt_fixed(float(v), int(kind[5:]))
    if kind == "sci":
        return _fmt_sci(v)
    if kind == "jfloat":
        return java_float_to_string(float(v))
    if kind == "dec2":
        return _fmt_dec2(v)
    raise ValueError(f"unknown kind {kind!r}")


def _isna(v):
    if v is None or v is pd.NA:
        return True
    try:
        return bool(v != v)  # NaN / NaT
    except (TypeError, ValueError):
        return False


# --------------------------------------------------------------------------
# Layout
# --------------------------------------------------------------------------


class TableLayout:
    """How a table is laid out on disk: column order, kinds, quoting, NA text."""

    __slots__ = ("table", "columns", "kinds", "quotes", "na_reps", "header_quoted",
                 "line_terminator", "notes")

    def __init__(self, table, columns, kinds, quotes, na_reps, header_quoted,
                 line_terminator="\n", notes=None):
        self.table = table
        self.columns = list(columns)
        self.kinds = list(kinds)
        self.quotes = list(quotes)
        self.na_reps = list(na_reps)
        self.header_quoted = list(header_quoted)
        self.line_terminator = line_terminator
        self.notes = list(notes or [])

    @classmethod
    def for_columns(cls, table, columns, df=None):
        """Layout from the static PEAKS spec (unknown columns inferred from dtype)."""
        kinds, quotes, nas, hqs = [], [], [], []
        for c in columns:
            sp = _spec_for(table, c) if table else None
            if sp is None:
                kind = "str"
                if df is not None and c in df.columns:
                    dt = df[c].dtype
                    if pd.api.types.is_bool_dtype(dt):
                        kind = "bool"
                    elif pd.api.types.is_integer_dtype(dt):
                        kind = "int"
                    elif pd.api.types.is_float_dtype(dt):
                        kind = "jfloat"
                sp = (kind, "never", "", True)
            kinds.append(sp[0])
            quotes.append(sp[1])
            nas.append(sp[2])
            hqs.append(sp[3])
        return cls(table, columns, kinds, quotes, nas, hqs)

    def entry(self, column):
        i = self.columns.index(column)
        return self.kinds[i], self.quotes[i], self.na_reps[i], self.header_quoted[i]

    def _key(self):
        return (self.table, self.columns, self.kinds, self.quotes, self.na_reps,
                self.header_quoted, self.line_terminator)

    def __eq__(self, other):
        return isinstance(other, TableLayout) and self._key() == other._key()

    def __repr__(self):
        return f"TableLayout({self.table!r}, {len(self.columns)} columns)"

    def to_dict(self):
        return {
            "table": self.table,
            "columns": [
                {"name": c, "kind": k, "quote": q, "na": n, "header_quoted": h}
                for c, k, q, n, h in zip(self.columns, self.kinds, self.quotes,
                                         self.na_reps, self.header_quoted)
            ],
            "line_terminator": self.line_terminator,
            "notes": self.notes,
        }


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------

_FIELD = re.compile(r'"((?:[^"]|"")*)"(?=,|$)|([^,]*)')
_NEEDS_QUOTE = re.compile(r'[",\n\r]')


def _split_line(line):
    vals, quoted = [], []
    pos, n = 0, len(line)
    while True:
        m = _FIELD.match(line, pos)
        g1 = m.group(1)
        if g1 is not None:
            vals.append(g1.replace('""', '"') if '""' in g1 else g1)
            quoted.append(True)
        else:
            vals.append(m.group(2))
            quoted.append(False)
        pos = m.end()
        if pos >= n:
            return vals, quoted
        if line[pos] != ",":  # pragma: no cover - malformed
            raise ValueError(f"malformed CSV line near char {pos}: {line[:120]!r}")
        pos += 1
        if pos == n:
            vals.append("")
            quoted.append(False)
            return vals, quoted


_INT_RE = re.compile(r"-?(0|[1-9]\d*)")
_FIX_RE = re.compile(r"-?\d+\.(\d+)")
_SCI_RE = re.compile(r"-?\d(\.\d{1,2})?E-?\d+")


def _infer_kind(values):
    nonempty = [v for v in values if v != "" and v != "NaN"]
    if not nonempty:
        return "str"
    if all(_INT_RE.fullmatch(v) for v in nonempty):
        return "int"
    if all(v in ("true", "false") for v in nonempty):
        return "bool"
    decs = set()
    for v in nonempty:
        m = _FIX_RE.fullmatch(v)
        if not m:
            break
        decs.add(len(m.group(1)))
    else:
        if len(decs) == 1 and 1 <= next(iter(decs)) <= 6:
            return f"fixed{next(iter(decs))}"
        return "jfloat"
    if all(_SCI_RE.fullmatch(v) for v in nonempty):
        return "sci"
    return "str"


def _infer_quote(values, quoted, preferred):
    """Pick a quoting mode consistent with the observed file, preferring the spec."""
    ok = {
        "always": all(quoted),
        "never": all((not q) or bool(_NEEDS_QUOTE.search(v)) for v, q in zip(values, quoted)),
        "nonempty": all(q == (v != "") or (q and _NEEDS_QUOTE.search(v)) for v, q in zip(values, quoted)),
    }
    if preferred in ok and ok[preferred]:
        return preferred, True
    for mode in _QUOTES:
        if ok[mode]:
            return mode, True
    nq = sum(quoted)
    return ("always" if nq * 2 >= len(quoted) else "never"), False


def _convert(values, kind, na_rep, verify):
    """Convert column strings to a Series; returns (series, kind, na_rep) with fallback to str."""
    if kind == "str":
        return pd.Series(values, dtype=object), "str", ""
    if kind == "int":
        if any(v != "" and not _INT_RE.fullmatch(v) for v in values):
            return _convert(values, "str", "", verify)
        if "" in values:
            s = pd.array([int(v) if v != "" else None for v in values], dtype="Int64")
            return pd.Series(s), "int", ""
        return pd.Series(np.array([int(v) for v in values], dtype=np.int64)), "int", ""
    if kind == "bool":
        if any(v not in ("true", "false") for v in values):
            return _convert(values, "str", "", verify)
        return pd.Series(np.array([v == "true" for v in values], dtype=bool)), "bool", ""
    # floats
    has_empty = "" in values
    has_nan = "NaN" in values
    if has_empty and has_nan:
        return _convert(values, "str", "", verify)
    na = "NaN" if has_nan else ("" if has_empty else na_rep)
    try:
        arr = np.array([np.nan if (v == "" or v == "NaN") else float(v) for v in values], dtype=np.float64)
    except ValueError:
        return _convert(values, "str", "", verify)
    if verify:
        for v, x in zip(values, arr):
            if v == "" or v == "NaN":
                continue
            if format_value(kind, x) != v:
                return _convert(values, "str", "", verify)
    return pd.Series(arr), kind, na


def read_table(path, table=None, verify=True):
    """Read one PEAKS CSV.  Returns ``(DataFrame, TableLayout)``.

    ``verify=True`` checks that every numeric cell re-formats to the exact text
    in the file; columns failing that are kept as text so writing stays exact.
    """
    with open(path, "r", encoding="utf-8", newline="") as fh:
        text = fh.read()
    lt = "\r\n" if "\r\n" in text[:10000] else "\n"
    lines = text.split(lt)
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        raise ValueError(f"{path}: empty file")
    header, hq = _split_line(lines[0])
    ncol = len(header)
    rows_v = []
    rows_q = []
    for i, ln in enumerate(lines[1:], start=2):
        v, q = _split_line(ln)
        if len(v) != ncol:
            raise ValueError(f"{path}:{i}: {len(v)} fields, header has {ncol}")
        rows_v.append(v)
        rows_q.append(q)
    cols_v = list(zip(*rows_v)) if rows_v else [()] * ncol
    cols_q = list(zip(*rows_q)) if rows_q else [()] * ncol
    if len(set(header)) != len(header):
        raise ValueError(f"{path}: duplicate column names")

    data = {}
    kinds, quotes, nas, notes = [], [], [], []
    for j, col in enumerate(header):
        vals = list(cols_v[j])
        sp = _spec_for(table, col) if table else None
        if sp is None:
            kind, pref_q, na = _infer_kind(vals), "never", ""
            notes.append(f"column {col!r} not in PEAKS spec; inferred kind {kind}")
        else:
            kind, pref_q, na, _ = sp
        quote, exact = _infer_quote(vals, cols_q[j], pref_q)
        if not exact:
            notes.append(f"column {col!r}: inconsistent quoting, write will not be byte-exact")
        series, kind2, na2 = _convert(vals, kind, na, verify)
        if kind2 != kind:
            notes.append(f"column {col!r}: values do not follow kind {kind}; kept as text")
        data[col] = series
        kinds.append(kind2)
        quotes.append(quote)
        nas.append(na2)
    df = pd.DataFrame(data, columns=header)
    if not rows_v:
        df = pd.DataFrame({c: pd.Series([], dtype=_dtype_for(k)) for c, k in zip(header, kinds)})
    layout = TableLayout(table, header, kinds, quotes, nas, hq, lt, notes)
    return df, layout


# --------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------


def _quote(s, mode):
    if mode == "always" or (mode == "nonempty" and s != "") or _NEEDS_QUOTE.search(s):
        return '"' + s.replace('"', '""') + '"'
    return s


def _format_column(series, kind, quote, na_rep):
    out = []
    append = out.append
    if kind == "str":
        for v in series.tolist():
            s = "" if _isna(v) else str(v)
            append(_quote(s, quote))
        return out
    for v in series.tolist():
        s = na_rep if _isna(v) else format_value(kind, v)
        append(_quote(s, quote))
    return out


def _reconcile_layout(df, table, layout):
    cols = list(df.columns)
    if layout is not None and layout.columns == cols:
        return layout
    spec_layout = TableLayout.for_columns(table, cols, df)
    if layout is None:
        return spec_layout
    kinds, quotes, nas, hqs = [], [], [], []
    for c in cols:
        src = layout if c in layout.columns else spec_layout
        k, q, n, h = src.entry(c)
        kinds.append(k)
        quotes.append(q)
        nas.append(n)
        hqs.append(h)
    return TableLayout(table, cols, kinds, quotes, nas, hqs, layout.line_terminator)


def write_table(df, path, table=None, layout=None):
    """Write a DataFrame in PEAKS CSV format.  Returns the layout used.

    Column order is the DataFrame's column order; use :func:`canonical_columns`
    or :func:`empty_table` to get PEAKS order for a pipeline-produced table.
    """
    layout = _reconcile_layout(df, table, layout)
    header = ",".join(
        _quote(c, "always" if h else "never") for c, h in zip(layout.columns, layout.header_quoted)
    )
    formatted = [
        _format_column(df[c], k, q, n)
        for c, k, q, n in zip(layout.columns, layout.kinds, layout.quotes, layout.na_reps)
    ]
    lt = layout.line_terminator
    parts = [header]
    parts.extend(",".join(row) for row in zip(*formatted))
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(lt.join(parts) + lt)
    return layout


# --------------------------------------------------------------------------
# db.summary.html
# --------------------------------------------------------------------------


def _clean(s):
    s = re.sub(r"<[^>]+>", " ", s)
    s = html.unescape(s).replace("\xa0", " ")
    return re.sub(r"\s+", " ", s).strip()


def _table_after(text, caption):
    m = re.search(re.escape(caption) + r".*?(<table.*?</table>)", text, re.S)
    if not m:
        return None
    rows = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", m.group(1), re.S):
        rows.append([_clean(c) for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", tr, re.S)])
    return rows


def _num(s):
    s = s.replace("≥", "").replace(">=", "").replace("%", "").replace(",", "").strip()
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


def _snake(s):
    return re.sub(r"[^0-9a-z]+", "_", s.lower()).strip("_")


def _yes(s):
    return s.strip().lower() in ("yes", "enabled", "true", "on")


_T6_SIMPLE = {
    "Search Engine Name": ("search_engine", str),
    "Precursor Mass Search Type": ("precursor_mass_type", str),
    "Enzyme": ("enzyme", str),
    "Max Missed Cleavages": ("max_missed_cleavages", int),
    "Digest Mode": ("digest_mode", str),
    "Max Variable PTM Per Peptide": ("max_variable_ptm_per_peptide", int),
    "Database": ("database", str),
    "Taxon": ("taxon", str),
    "Searched Entries": ("searched_entries", int),
    "Deep Learning Boost": ("deep_learning_boost", _yes),
    "FDR Estimation": ("fdr_estimation", _yes),
    "Find Unspecified PTMs with PEAKS PTM": ("peaks_ptm", _yes),
    "Find Mutations and Sequence Variants with SPIDER": ("spider", _yes),
}
_T2 = {
    "Peptide -10lgP": "peptide_min_10lgp",
    "Protein -10lgP": "protein_min_10lgp",
    "Proteins unique peptides": "protein_min_unique_peptides",
    "De novo score(%)": "denovo_min_alc_pct",
}
_T3 = {
    "FDR (Peptide-Spectrum Matches)": "fdr_psm_pct",
    "FDR (Peptide Sequences)": "fdr_peptide_pct",
    "FDR (Protein Group)": "fdr_protein_group_pct",
    "De Novo Only Spectra": "denovo_only_spectra",
}
_T1_KEYS = ["ms1_scans", "ms2_scans", "chimera_scans", "features", "psms", "identified_scans",
            "identified_features", "peptides", "sequences", "protein_groups", "proteins_all",
            "proteins_top"]
_MOD_RE = re.compile(r"^(.*?)\s*\(([+-]?\d*\.?\d+)\)$")


def parse_summary_html(path_or_text):
    """Parse search parameters (Table 6), thresholds (Table 2), FDRs (Table 3),
    data statistics (Table 1) and missed cleavages (Table 5) from db.summary.html."""
    if "<" not in path_or_text[:200] and os.path.exists(path_or_text):
        with open(path_or_text, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    else:
        text = path_or_text
    p = {"raw": {}}

    # ---- Table 6: search parameters
    m = re.search(r"Table 6\.</b>.*?<div class=\"align-left\">(.*?)</div>", text, re.S)
    raw6 = {}
    if m:
        current_list = None
        for item in re.split(r"<br\s*/?>", m.group(1)):
            it = html.unescape(item).replace("\xa0", " ")
            if not it.strip():
                continue
            is_sub = bool(re.match(r"\s*(&nbsp;?|\xa0)", item))
            it = _clean(it)
            if is_sub and current_list is not None:
                raw6[current_list].append(it)
                continue
            if ":" not in it:
                continue
            label, val = [x.strip() for x in it.split(":", 1)]
            if val == "" and label.endswith("Modifications"):
                raw6[label] = []
                current_list = label
            else:
                raw6[label] = val
                current_list = None
    p["raw"]["search_parameters"] = raw6
    for label, val in raw6.items():
        if label in _T6_SIMPLE:
            key, conv = _T6_SIMPLE[label]
            try:
                p[key] = conv(val)
            except (TypeError, ValueError):
                p[key] = val
        elif label in ("Parent Mass Error Tolerance", "Fragment Mass Error Tolerance"):
            key = "precursor_tolerance" if label.startswith("Parent") else "fragment_tolerance"
            mm = re.match(r"([\d.]+)\s*(\S+)", val)
            if mm:
                p[key] = float(mm.group(1))
                p[key + "_unit"] = mm.group(2)
            else:
                p[key] = val
        elif label == "Peptide Length Range":
            mm = re.match(r"(\d+)\s*-\s*(\d+)", val)
            if mm:
                p["peptide_length_min"], p["peptide_length_max"] = int(mm.group(1)), int(mm.group(2))
        elif label == "Confident Amino Acid Threshold":
            p["confident_aa_threshold_pct"] = _num(val)
        elif label.endswith("Modifications") and isinstance(val, list):
            mods = []
            for s in val:
                mm = _MOD_RE.match(s)
                mods.append({"name": mm.group(1), "delta": float(mm.group(2))} if mm else {"name": s, "delta": None})
            p["fixed_modifications" if label.startswith("Fixed") else "variable_modifications"] = mods
        else:
            p[_snake(label)] = val

    # ---- Tables 2 and 3
    for cap, mapping, rk in (("Table 2.</b>", _T2, "filtration"), ("Table 3.</b>", _T3, "filtered_statistics")):
        rows = _table_after(text, cap) or []
        raw = {r[0]: r[1] for r in rows if len(r) >= 2}
        p["raw"][rk] = raw
        for label, val in raw.items():
            p[mapping.get(label, _snake(label))] = _num(val)

    # ---- Table 1: statistics of data
    rows = _table_after(text, "Table 1.</b>") or []
    stats = {}
    for r in rows:
        if len(r) == len(_T1_KEYS) + 1 and all(isinstance(_num(x), int) for x in r[1:]):
            stats[r[0]] = dict(zip(_T1_KEYS, (_num(x) for x in r[1:])))
    p["statistics"] = stats

    # ---- Table 5: missed cleavages
    rows = _table_after(text, "Table 5.</b>") or []
    if rows:
        hdr = rows[0][1:]
        p["missed_cleavages"] = {r[0]: dict(zip(hdr, (_num(x) for x in r[1:]))) for r in rows[1:]}
    return p


# --------------------------------------------------------------------------
# Export directory
# --------------------------------------------------------------------------


class PeaksExport:
    """In-memory PEAKS export: one DataFrame per table (None if absent)."""

    TABLE_NAMES = tuple(TABLES)

    def __init__(self, psms=None, peptides=None, proteins=None, protein_peptides=None,
                 features=None, denovo_only=None, ptm_profile=None, modified_peptides=None,
                 params=None, sample_name=None, layouts=None, source_dir=None):
        self.psms = psms
        self.peptides = peptides
        self.proteins = proteins
        self.protein_peptides = protein_peptides
        self.features = features
        self.denovo_only = denovo_only
        self.ptm_profile = ptm_profile
        self.modified_peptides = modified_peptides
        self.params = params if params is not None else {}
        self.sample_name = sample_name
        self.layouts = layouts if layouts is not None else {}
        self.source_dir = source_dir

    def tables(self):
        """Dict of table name -> DataFrame for present tables."""
        return {n: getattr(self, n) for n in self.TABLE_NAMES if getattr(self, n) is not None}

    def __repr__(self):
        parts = [f"{n}={len(df)}" for n, df in self.tables().items()]
        return f"PeaksExport(sample={self.sample_name!r}, {', '.join(parts)})"


def _find_sample_file(directory, suffix, sample_name):
    if sample_name:
        p = os.path.join(directory, f"{sample_name}{suffix}")
        return (p, sample_name) if os.path.exists(p) else (None, None)
    hits = sorted(glob.glob(os.path.join(glob.escape(directory), "*" + suffix)))
    if not hits:
        return None, None
    if len(hits) > 1:
        raise ValueError(f"{directory}: several *{suffix} files; pass sample_name=")
    base = os.path.basename(hits[0])
    return hits[0], base[: -len(suffix)]


def read_export(directory, tables=None, verify=True, sample_name=None):
    """Read a PEAKS export directory into a :class:`PeaksExport`.

    ``tables``: optional iterable of table names to load (default all).
    """
    if not os.path.isdir(directory):
        raise FileNotFoundError(directory)
    want = set(TABLES) if tables is None else set(tables)
    unknown = want - set(TABLES)
    if unknown:
        raise ValueError(f"unknown tables {sorted(unknown)}")
    exp = PeaksExport(source_dir=os.path.abspath(directory))
    for name, tmpl in TABLES.items():
        if name in _PER_SAMPLE:
            path, sname = _find_sample_file(directory, tmpl.replace("{sample}", ""), sample_name)
            if sname and exp.sample_name is None:
                exp.sample_name = sname
        else:
            path = os.path.join(directory, tmpl)
            if not os.path.exists(path):
                path = None
        if path is None or name not in want:
            continue
        df, layout = read_table(path, name, verify=verify)
        setattr(exp, name, df)
        exp.layouts[name] = layout
    if exp.sample_name is None and exp.peptides is not None:
        for c in exp.peptides.columns:
            if c.startswith("Area "):
                exp.sample_name = c[5:]
                break
    summ = os.path.join(directory, SUMMARY_FILE)
    if os.path.exists(summ):
        exp.params = parse_summary_html(summ)
    return exp


def write_export(export, directory, tables=None):
    """Write the tables of ``export`` as PEAKS CSVs into ``directory``.

    Uses the layouts captured by :func:`read_export` when available, otherwise
    the static PEAKS spec.  Returns {table name: path}.
    """
    os.makedirs(directory, exist_ok=True)
    out = {}
    for name, df in export.tables().items():
        if tables is not None and name not in tables:
            continue
        fname = TABLES[name]
        if "{sample}" in fname:
            fname = fname.format(sample=export.sample_name or "Sample 1")
        path = os.path.join(directory, fname)
        write_table(df, path, name, export.layouts.get(name))
        out[name] = path
    return out
