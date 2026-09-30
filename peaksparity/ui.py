"""Self-contained HTML results viewer for a PEAKS-format DB-search export (stage 10, `ui`).

Usage::

    python -m peaksparity.ui --export <export dir> [--mzml run.mzML ...] [--fasta db.fasta]
        [--ref <PEAKS export dir>] --out report.html

Writes one HTML file with inline CSS, JavaScript and data (no network access, no external
libraries) that opens from ``file://``.  Tabs: Summary (counts, filters, search parameters,
score / precursor-error plots, missed cleavages, PTM profile, and a Concordance panel when
``--ref`` is given), Proteins (group table, sequence coverage map, supporting peptides),
Peptides and PSMs (sortable / filterable / paginated tables) and an annotated MS2 spectrum
viewer (b/y ions, 2+ fragments for precursor charge >= 3, optional H2O/NH3 losses, fragment
table with mass errors).

Peak lists are read from ``--mzml`` with pyopenms (random access by raw scan number: the
native id ``... scan=N``) and only for the best PSM of each peptide, keeping the top
``--top-peaks`` peaks per spectrum and stopping when the embedded peak data reach
``--max-spectra-mb``.  See docs/UI.md.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from . import __version__
from . import seqnorm
from .concordance import compare as concordance_compare
from .concordance import norm_accession, norm_source
from .peaks_export import read_export
from .proteins import Enzyme, read_protein_db

__all__ = ["build_report", "collect_report_data", "read_spectra", "render_html", "main"]

_ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ui_assets")
_TABLES = ("psms", "peptides", "proteins", "protein_peptides", "ptm_profile")
_SCAN_RE = re.compile(r"\bscan=(\d+)\b")
DEFAULT_FRAGMENT_TOL = (0.5, "Da")


def _log(msg):
    print(f"[ui] {msg}", file=sys.stderr, flush=True)


# --------------------------------------------------------------------------
# JSON helpers
# --------------------------------------------------------------------------


def _py(v):
    """Plain JSON-able Python value (NaN/NA -> None)."""
    if v is None or v is pd.NA:
        return None
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        f = float(v)
        return None if (math.isnan(f) or math.isinf(f)) else f
    if isinstance(v, dict):
        return {str(k): _py(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, set, frozenset)):
        return [_py(x) for x in v]
    if isinstance(v, np.ndarray):
        return [_py(x) for x in v.tolist()]
    return str(v)


def _column_values(s):
    if pd.api.types.is_bool_dtype(s.dtype):
        return [bool(x) for x in s.tolist()]
    if pd.api.types.is_integer_dtype(s.dtype):
        return [None if x is pd.NA or x is None else int(x) for x in s.tolist()]
    if pd.api.types.is_float_dtype(s.dtype):
        return [None if x != x else float(x) for x in s.tolist()]
    return [_py(x) if not isinstance(x, str) else x for x in s.tolist()]


def _guess_kind(s):
    if pd.api.types.is_bool_dtype(s.dtype):
        return "bool"
    if pd.api.types.is_integer_dtype(s.dtype):
        return "int"
    if pd.api.types.is_float_dtype(s.dtype):
        return "jfloat"
    return "str"


def table_payload(df, layout=None, extra=None):
    """``{"columns", "kinds", "rows"}`` for a DataFrame; ``extra`` = [(name, kind, values)]."""
    cols = list(df.columns)
    if layout is not None and list(layout.columns) == cols:
        kinds = list(layout.kinds)
    else:
        kinds = [_guess_kind(df[c]) for c in cols]
    data = [_column_values(df[c]) for c in cols]
    for name, kind, values in reversed(list(extra or ())):  # extra columns go first
        cols.insert(0, name)
        kinds.insert(0, kind)
        data.insert(0, [_py(v) for v in values])
    rows = [list(r) for r in zip(*data)] if data and len(df) else []
    return {"columns": cols, "kinds": kinds, "rows": rows}


def _dumps(obj):
    """Compact JSON safe to embed in a <script type="application/json"> element."""
    s = json.dumps(obj, separators=(",", ":"), allow_nan=False, ensure_ascii=False)
    # "<" only occurs inside JSON strings, where \u003c is an equivalent escape; this rules
    # out "</script>" and "<!--" inside the element.  U+2028/9 are escaped for old parsers.
    return s.replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")


# --------------------------------------------------------------------------
# Peptides and modifications
# --------------------------------------------------------------------------


def _exact_delta_map():
    """{"+15.99": 15.994915, ...} from seqnorm.MOD_MASSES (first name wins)."""
    out = {}
    for v in seqnorm.MOD_MASSES.values():
        out.setdefault(seqnorm.format_delta(v), float(v))
    return out


_EXACT = _exact_delta_map()


def peptide_info(peptides):
    """{peptide string: [stripped, [[position, exact delta], ...]]}; position 0 = N-term,
    1..n residues, n+1 = C-term.  Rounded PEAKS deltas are replaced by the exact
    monoisotopic delta of the known modification with that rounding."""
    out = {}
    for p in peptides:
        if not isinstance(p, str) or p in out or not p:
            continue
        try:
            stripped, mods = seqnorm.parse(p)
        except ValueError:
            continue
        ml = []
        for pos, d in mods:
            key = seqnorm.format_delta(d)
            ml.append([int(pos), _EXACT.get(key, float(d))])
        out[p] = [stripped, ml]
    return out


def _missed_cleavages(stripped_list, enzyme):
    """Histogram {"0", "1", "2", "3", "4+"} of internal cleavage sites per peptide."""
    c = Counter()
    for s in stripped_list:
        n = 0
        if enzyme.cleave_after:
            for j in range(1, len(s)):
                if s[j - 1] in enzyme.cleave_after and s[j] not in enzyme.not_before:
                    n += 1
        c["4+" if n >= 4 else str(n)] += 1
    return {k: int(c.get(k, 0)) for k in ("0", "1", "2", "3", "4+")}


def _split(v, seps):
    if not isinstance(v, str) or not v:
        return []
    return [x.strip() for x in re.split(seps, v) if x.strip()]


def _mod_names(params, frames):
    """{rounded delta string: PTM name} from the search parameters and PTM columns."""
    names = {}
    for key in ("fixed_modifications", "variable_modifications"):
        for m in params.get(key) or []:
            if m.get("delta") is not None:
                names.setdefault(seqnorm.format_delta(m["delta"]), m["name"])
    for df in frames:
        if df is None or "PTM" not in df.columns or "Peptide" not in df.columns:
            continue
        for pep, ptm in zip(df["Peptide"].tolist()[:20000], df["PTM"].tolist()[:20000]):
            nm = _split(ptm, r";")
            if not nm:
                continue
            try:
                _st, mods = seqnorm.parse(pep)
            except (ValueError, TypeError):
                continue
            if len(mods) == len(nm):
                for (_p, d), n in zip(mods, nm):
                    names.setdefault(seqnorm.format_delta(d), n)
    return names


def computed_ptm_profile(peptides, protein_peptides, pepinfo, names):
    """PTM profile (Name, Mass, Position, # peptides, # proteins, # modification sites)
    computed from the peptide tables when db.ptmProfile.csv is absent."""
    agg = {}

    def ent(delta_s, name):
        return agg.setdefault(name, {"mass": delta_s, "res": set(), "peps": set(),
                                     "prots": set(), "sites": set()})

    if protein_peptides is not None and {"Peptide", "Accession"} <= set(protein_peptides.columns):
        starts = protein_peptides["Start"].tolist() if "Start" in protein_peptides.columns else [None] * len(protein_peptides)
        for pep, acc, st in zip(protein_peptides["Peptide"], protein_peptides["Accession"], starts):
            info = pepinfo.get(pep)
            if not info:
                continue
            for pos, d in info[1]:
                ds = seqnorm.format_delta(d)
                e = ent(ds, names.get(ds, ds))
                res = "N-term" if pos == 0 else ("C-term" if pos > len(info[0]) else info[0][pos - 1])
                e["res"].add(res)
                e["peps"].add(pep)
                e["prots"].add(acc)
                if st is not None and st == st:
                    e["sites"].add((acc, int(st) + max(pos, 1) - 1))
    elif peptides is not None:
        for pep, accs in zip(peptides["Peptide"], peptides.get("Accession", [""] * len(peptides))):
            info = pepinfo.get(pep)
            if not info:
                continue
            for pos, d in info[1]:
                ds = seqnorm.format_delta(d)
                e = ent(ds, names.get(ds, ds))
                e["res"].add(info[0][pos - 1] if 1 <= pos <= len(info[0]) else "term")
                e["peps"].add(pep)
                for a in _split(accs, r"[:;]"):
                    e["prots"].add(a)
    rows = []
    for name, e in sorted(agg.items(), key=lambda kv: -len(kv[1]["peps"])):
        rows.append([name, e["mass"], "".join(sorted(e["res"])), len(e["peps"]), len(e["prots"]),
                     len(e["sites"]) if e["sites"] else None])
    return {"columns": ["Name", "Mass (Da)", "Position", "# peptides", "# proteins", "# modification sites"],
            "kinds": ["str", "str", "str", "int", "int", "int"], "rows": rows}


# --------------------------------------------------------------------------
# Spectra
# --------------------------------------------------------------------------


def _scan_of(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        m = re.search(r"(\d+)$", str(v))
        return int(m.group(1)) if m else None


def wanted_spectra(exp):
    """Ordered (source key, scan) of the best PSM of every peptide (best score first)."""
    df = exp.peptides
    if df is None or not {"Scan", "Source File"} <= set(df.columns):
        df = exp.psms
        if df is None:
            return []
        df = df.sort_values("-10LgP", ascending=False, kind="stable").drop_duplicates("Peptide")
    else:
        df = df.sort_values("-10LgP", ascending=False, kind="stable")
    out, seen = [], set()
    for src, scan in zip(df["Source File"].tolist(), df["Scan"].tolist()):
        sc = _scan_of(scan)
        if sc is None:
            continue
        k = (norm_source(src), sc)
        if k not in seen:
            seen.add(k)
            out.append(k)
    return out


class _MzmlReader:
    """Random access to MS2 spectra of an indexed mzML by raw scan number (pyopenms)."""

    def __init__(self, path):
        import pyopenms as oms
        self.path = path
        self.od = oms.OnDiscMSExperiment()
        if not self.od.openFile(path):
            raise IOError(f"pyopenms could not open {path} (needs an indexed mzML)")
        self.n = self.od.getNrSpectra()
        self.index = None  # scan -> spectrum index, built only if index != scan - 1

    def _build_index(self):
        _log(f"{os.path.basename(self.path)}: native ids are not in scan order; indexing all spectra")
        meta = self.od.getMetaData()
        self.index = {}
        for i in range(meta.getNrSpectra()):
            m = _SCAN_RE.search(meta.getSpectrum(i).getNativeID())
            if m:
                self.index[int(m.group(1))] = i

    def get(self, scan):
        if self.index is None:
            i = scan - 1
            if 0 <= i < self.n:
                s = self.od.getSpectrum(i)
                m = _SCAN_RE.search(s.getNativeID())
                if m and int(m.group(1)) == scan:
                    return s
            self._build_index()
        i = self.index.get(scan)
        return None if i is None else self.od.getSpectrum(i)


def _spectrum_entry(s, top_n):
    mz, inten = s.get_peaks()
    mz = np.asarray(mz, dtype=np.float64)
    inten = np.asarray(inten, dtype=np.float64)
    keep = inten > 0
    mz, inten = mz[keep], inten[keep]
    if len(mz) > top_n:
        idx = np.argpartition(-inten, top_n - 1)[:top_n]
        mz, inten = mz[idx], inten[idx]
    order = np.argsort(mz)
    mz, inten = mz[order], inten[order]
    bp = float(inten.max()) if len(inten) else 0.0
    rel = np.round(100.0 * inten / bp, 2) if bp > 0 else inten
    e = {
        "mz": [round(float(x), 4) for x in mz],
        "i": [float(x) for x in rel],
        "bp": float(f"{bp:.4g}"),
        "rt": round(s.getRT() / 60.0, 4),
        "ms": int(s.getMSLevel()),
        "id": s.getNativeID(),
    }
    prec = s.getPrecursors()
    if prec:
        e["pmz"] = round(prec[0].getMZ(), 5)
        e["pz"] = int(prec[0].getCharge())
    return e


def read_spectra(mzml_paths, wanted, top_n=150, max_bytes=40e6):
    """Peak lists for ``wanted`` [(source key, scan), ...] in priority order.

    Returns ``(spectra, info)`` with spectra = {"<source key>|<scan>": entry}.  Stops adding
    spectra once their JSON exceeds ``max_bytes``.
    """
    info = {"requested": len(wanted), "embedded": 0, "not_found": 0, "no_mzml": 0,
            "skipped_size_cap": 0, "bytes": 0, "top_n": top_n,
            "max_bytes": int(max_bytes), "files": {}}
    spectra = {}
    if not mzml_paths or not wanted:
        info["no_mzml"] = len(wanted) if not mzml_paths else 0
        return spectra, info
    by_key = {norm_source(os.path.basename(p)): p for p in mzml_paths}
    sources = {s for s, _ in wanted}
    if len(mzml_paths) == 1 and len(sources) == 1 and next(iter(sources)) not in by_key:
        by_key[next(iter(sources))] = mzml_paths[0]
        _log(f"single mzML: mapping source {next(iter(sources))!r} to {mzml_paths[0]}")
    readers = {}
    total = 0
    for src, scan in wanted:
        path = by_key.get(src)
        if path is None:
            info["no_mzml"] += 1
            continue
        if total >= max_bytes:
            info["skipped_size_cap"] += 1
            continue
        if path not in readers:
            _log(f"opening {path}")
            readers[path] = _MzmlReader(path)
            info["files"][src] = path
        s = readers[path].get(scan)
        if s is None or s.size() == 0:
            info["not_found"] += 1
            continue
        e = _spectrum_entry(s, top_n)
        size = len(json.dumps(e, separators=(",", ":"))) + len(src) + 12
        if total + size > max_bytes:
            info["skipped_size_cap"] += 1
            total = max_bytes  # stop here: keep priority order strict
            continue
        spectra[f"{src}|{scan}"] = e
        total += size
    info["embedded"] = len(spectra)
    info["bytes"] = int(min(total, max_bytes)) if spectra else 0
    return spectra, info


# --------------------------------------------------------------------------
# Report data
# --------------------------------------------------------------------------


def _fragment_tolerance(params, override=None, unit=None):
    if override is not None:
        return {"value": float(override), "unit": unit or "Da", "source": "command line"}
    v = params.get("fragment_tolerance")
    if isinstance(v, (int, float)):
        u = str(params.get("fragment_tolerance_unit") or "Da")
        return {"value": float(v), "unit": "ppm" if u.lower() == "ppm" else "Da",
                "source": "export search parameters"}
    return {"value": DEFAULT_FRAGMENT_TOL[0], "unit": DEFAULT_FRAGMENT_TOL[1],
            "source": "default (no fragment tolerance in the export)"}


def _counts(exp):
    c = {}
    ps = exp.psms
    if ps is not None:
        c["psm_rows"] = len(ps)
        if {"Scan", "Source File"} <= set(ps.columns):
            src = [norm_source(x) for x in ps["Source File"].tolist()]
            scans = ps["Scan"].tolist()
            c["identified_scans"] = len(set(zip(src, scans)))
            if "Feature Id" in ps.columns:
                c["psms"] = len(set(zip(src, scans, ps["Feature Id"].tolist())))
            keys = set()
            for s, sc, p, z in zip(src, scans, ps["Peptide"].tolist(), ps["z"].tolist()):
                try:
                    p = seqnorm.canonical(p)
                except ValueError:
                    pass
                keys.add((s, sc, p, z))
            c["psm_keys"] = len(keys)
    pep = exp.peptides if exp.peptides is not None else None
    peps = pep["Peptide"].tolist() if pep is not None else (
        sorted(set(ps["Peptide"].tolist())) if ps is not None else [])
    c["peptide_rows"] = len(peps)
    c["peptides"] = len({seqnorm.il_collapse(p) for p in peps})
    stripped = set()
    for p in peps:
        try:
            stripped.add(seqnorm.il_collapse(seqnorm.strip_mods(p)))
        except ValueError:
            stripped.add(p)
    c["sequences"] = len(stripped)
    pr = exp.proteins
    if pr is not None:
        c["proteins_all"] = len(pr)
        if "Protein Group" in pr.columns:
            c["protein_groups"] = int(pr["Protein Group"].nunique())
        if "Top" in pr.columns:
            c["proteins_top"] = int(pr["Top"].astype(str).str.lower().eq("true").sum())
    return c


def _concordance(ref, export_dir, exp):
    """Run peaksparity.concordance.compare and derive per-row "In ref" flags for the tables."""
    _log(f"concordance vs {ref}")
    res = concordance_compare(ref, export_dir)
    m, t = res["metrics"], res["tables"]
    flags = {}
    # PSM rows: in ref unless their key is in the test-only table
    tonly = t.get("psm_test_only")
    if exp.psms is not None and tonly is not None:
        only = set(zip(tonly["source"], tonly["scan"].astype(int), tonly["modseq"], tonly["z"].astype(int)))
        vals = []
        for src, scan, pep, z in zip(exp.psms["Source File"], exp.psms["Scan"], exp.psms["Peptide"], exp.psms["z"]):
            try:
                ms = seqnorm.canonical(pep)
            except ValueError:
                ms = "UNPARSED:" + str(pep)
            key = (norm_source(src), _scan_of(scan), ms, int(z) if z == z else 0)
            vals.append(key not in only)
        flags["psms"] = vals
    ponly = t.get("peptide_test_only")
    if exp.peptides is not None and ponly is not None:
        only = set(ponly["peptide"].tolist())
        flags["peptides"] = [p not in only for p in exp.peptides["Peptide"].tolist()]
    gonly = t.get("protein_group_test_only")
    if exp.proteins is not None and gonly is not None and "Protein Group" in exp.proteins.columns:
        only = set(gonly["protein_group"].tolist())
        flags["proteins"] = [g not in only for g in exp.proteins["Protein Group"].tolist()]
    out_tables = {}
    for name, df in t.items():
        df = df.copy()
        for c in df.columns:  # sets / lists -> text
            if df[c].dtype == object:
                df[c] = [";".join(sorted(map(str, v))) if isinstance(v, (set, frozenset, list, tuple)) else v
                         for v in df[c].tolist()]
        out_tables[name] = table_payload(df)
    return {"metrics": _py(m), "tables": out_tables}, flags


def collect_report_data(export_dir, mzml=(), fasta=None, ref=None, top_n=150, max_spectra_mb=40.0,
                        frag_tol=None, frag_unit=None, title=None):
    """Everything the page needs.  Returns ``(data, spectra)`` (two JSON-able dicts)."""
    export_dir = os.path.abspath(export_dir)
    _log(f"reading export {export_dir}")
    exp = read_export(export_dir, tables=_TABLES, verify=False)
    params = exp.params or {}
    notes = []
    if not params:
        notes.append("The export has no db.summary.html: search parameters, filter thresholds and "
                     "FDRs are not available; counts are computed from the CSV files.")

    # peptide parsing (all peptide strings the page may annotate)
    peps = []
    for df in (exp.psms, exp.peptides, exp.protein_peptides):
        if df is not None and "Peptide" in df.columns:
            peps.extend(df["Peptide"].tolist())

    conc, flags = None, {}
    if ref:
        conc, flags = _concordance(ref, export_dir, exp)
        for name in ("psm_ref_only", "psm_test_only", "peptide_ref_only", "peptide_test_only"):
            tp = conc["tables"].get(name)
            if tp and "peptide" in tp["columns"]:
                j = tp["columns"].index("peptide")
                peps.extend(r[j] for r in tp["rows"])
    pepinfo = peptide_info(peps)

    # tables
    tables = {}
    for name in ("psms", "peptides", "proteins", "protein_peptides"):
        df = getattr(exp, name)
        if df is None:
            continue
        extra = [("In ref", "bool", flags[name])] if name in flags else None
        tables[name] = table_payload(df, exp.layouts.get(name), extra)

    # missed cleavages
    enzyme = Enzyme.from_params(params)
    pep_list = exp.peptides["Peptide"].tolist() if exp.peptides is not None else (
        sorted(set(exp.psms["Peptide"])) if exp.psms is not None else [])
    mc = {"computed": _missed_cleavages([pepinfo[p][0] for p in pep_list if p in pepinfo], enzyme),
          "enzyme": {"cleave_after": enzyme.cleave_after, "not_before": enzyme.not_before,
                     "specificity": enzyme.specificity, "max_missed": enzyme.max_missed}}
    if params.get("missed_cleavages"):
        mc["summary"] = _py(params["missed_cleavages"])

    # PTM profile
    if exp.ptm_profile is not None:
        ptm = table_payload(exp.ptm_profile, exp.layouts.get("ptm_profile"))
        ptm["source"] = "db.ptmProfile.csv"
    else:
        names = _mod_names(params, (exp.psms, exp.peptides, exp.protein_peptides))
        ptm = computed_ptm_profile(exp.peptides, exp.protein_peptides, pepinfo, names)
        ptm["source"] = "computed from the peptide tables (db.ptmProfile.csv absent)"

    # protein sequences
    fasta_used = fasta
    if not fasta_used and os.path.exists(os.path.join(export_dir, "db.proteins.fasta")):
        fasta_used = os.path.join(export_dir, "db.proteins.fasta")
    sequences = {}
    if fasta_used and exp.proteins is not None:
        _log(f"reading FASTA {fasta_used}")
        db = read_protein_db(fasta_used)
        missing = 0
        for acc in dict.fromkeys(exp.proteins["Accession"].tolist()):
            i = db.index.get(acc)
            if i is None:
                i = db.index.get(norm_accession(acc))
            if i is None:
                missing += 1
            else:
                sequences[acc] = db.sequence[i]
        if missing:
            notes.append(f"{missing} protein accession(s) were not found in {fasta_used}; "
                         "their coverage maps are unavailable.")
    elif exp.proteins is not None:
        notes.append("No FASTA given (--fasta): sequence coverage maps are unavailable.")

    # spectra
    mzml = [os.path.abspath(p) for p in (mzml or [])]
    wanted = wanted_spectra(exp)
    spectra, sinfo = read_spectra(mzml, wanted, top_n=top_n, max_bytes=max_spectra_mb * 1e6)
    if not mzml:
        notes.append("No mzML given (--mzml): spectra are not embedded; the spectrum view shows "
                     "theoretical fragments only.")
    else:
        if sinfo["skipped_size_cap"]:
            notes.append(f"{sinfo['skipped_size_cap']} spectra were not embedded because the "
                         f"{max_spectra_mb:g} MB cap was reached (lowest-scoring peptides dropped first).")
        if sinfo["not_found"]:
            notes.append(f"{sinfo['not_found']} scans were not found in the mzML file(s).")
        if sinfo["no_mzml"]:
            notes.append(f"{sinfo['no_mzml']} spectra belong to source files without a matching --mzml.")

    raw_params = (params.get("raw") or {}).get("search_parameters") or {}
    p_clean = {k: v for k, v in params.items() if k != "raw"}
    stats = params.get("statistics") or {}
    data = {
        "meta": {
            "title": title or f"{os.path.basename(export_dir)}",
            "export_dir": export_dir,
            "sample": exp.sample_name,
            "generated": _dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "tool": f"peaks-parity {__version__}",
            "fasta": os.path.abspath(fasta_used) if fasta_used else None,
            "mzml": mzml,
            "ref": os.path.abspath(ref) if ref else None,
        },
        "params": _py(p_clean),
        "raw_params": _py(raw_params),
        "statistics": _py(stats.get("Total") or (next(iter(stats.values())) if stats else {})),
        "frag_tol": _fragment_tolerance(params, frag_tol, frag_unit),
        "counts": _counts(exp),
        "tables": tables,
        "pepinfo": pepinfo,
        "sequences": sequences,
        "missed_cleavages": mc,
        "ptm_profile": ptm,
        "spectra_info": sinfo,
        "concordance": conc,
        "notes": notes,
    }
    return data, spectra


def _asset(name):
    with open(os.path.join(_ASSETS, name), encoding="utf-8") as fh:
        return fh.read()


def render_html(data, spectra):
    """The complete HTML document (str)."""
    title = f"Results · {data['meta']['title']}"
    css = _asset("report.css")
    js = _asset("report.js")
    if "</script" in js.lower() or "</style" in css.lower():
        raise ValueError("asset contains a closing tag")
    parts = [
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
        "<meta name=\"generator\" content=\"peaks-parity ui\">\n",
        f"<title>{_html_escape(title)}</title>\n<style>\n", css, "\n</style>\n</head>\n<body>\n",
        "<div id=\"app\"><noscript>This report needs JavaScript.</noscript></div>\n",
        "<script type=\"application/json\" id=\"pp-data\">", _dumps(data), "</script>\n",
        "<script type=\"application/json\" id=\"pp-spectra\">", _dumps(spectra), "</script>\n",
        "<script>\n", js, "\n</script>\n</body>\n</html>\n",
    ]
    return "".join(parts)


def _html_escape(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;"))


def build_report(export_dir, out, mzml=(), fasta=None, ref=None, top_n=150, max_spectra_mb=40.0,
                 frag_tol=None, frag_unit=None, title=None):
    """Write the HTML report to ``out``; returns a small dict of sizes/counts."""
    data, spectra = collect_report_data(export_dir, mzml=mzml, fasta=fasta, ref=ref, top_n=top_n,
                                        max_spectra_mb=max_spectra_mb, frag_tol=frag_tol,
                                        frag_unit=frag_unit, title=title)
    doc = render_html(data, spectra)
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    size = os.path.getsize(out)
    info = {"out": os.path.abspath(out), "bytes": size, "spectra": len(spectra),
            "spectra_bytes": data["spectra_info"]["bytes"],
            "rows": {k: len(v["rows"]) for k, v in data["tables"].items()}}
    _log(f"wrote {out} ({size / 1e6:.1f} MB, {len(spectra)} spectra)")
    return info


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m peaksparity.ui",
                                 description="Write a self-contained HTML results viewer for a PEAKS-format export.")
    ap.add_argument("--export", required=True, help="export directory (PEAKS DB-search CSV layout)")
    ap.add_argument("--mzml", nargs="*", default=[], help="indexed mzML file(s) for the spectrum view")
    ap.add_argument("--fasta", help="protein FASTA for coverage maps (default: <export>/db.proteins.fasta if present)")
    ap.add_argument("--ref", help="reference PEAKS export directory: adds a Concordance panel")
    ap.add_argument("--out", required=True, help="output .html")
    ap.add_argument("--top-peaks", type=int, default=150, help="peaks kept per embedded spectrum")
    ap.add_argument("--max-spectra-mb", type=float, default=40.0, help="cap on embedded peak data (MB)")
    ap.add_argument("--frag-tol", type=float, help="fragment tolerance (default: export params, else 0.5 Da)")
    ap.add_argument("--frag-unit", choices=["Da", "ppm"], help="unit of --frag-tol (default Da)")
    ap.add_argument("--title", help="report title (default: export directory name)")
    a = ap.parse_args(argv)
    info = build_report(a.export, a.out, mzml=a.mzml, fasta=a.fasta, ref=a.ref, top_n=a.top_peaks,
                        max_spectra_mb=a.max_spectra_mb, frag_tol=a.frag_tol, frag_unit=a.frag_unit,
                        title=a.title)
    print(json.dumps(info))
    return 0


if __name__ == "__main__":
    sys.exit(main())
