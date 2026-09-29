"""PEAKS ``db.ptmProfile.csv`` and ``db.modifiedPeptide.csv`` from the psms/peptides/protein tables.

Every rule here was inferred from the four DrewR_Adam_91526 reference exports and is
documented, with the evidence, in docs/PTM_TABLES.md.  Short version:

* A PTM site is **confident** when the PSM's positional confidence at the residue is >= 2 and
  its AScore for the site is >= 15.  Fixed modifications (no AScore) are always confident.
* ``db.modifiedPeptide.csv`` has one row per peptide carrying a variable modification.  The
  "Mod" half comes from the *last* PSM (psms.csv row order) among the peptide's PSMs with the
  best full-precision score; ``db.peptides.csv`` uses the *first*.  ``m/z`` columns hold the
  MS2 scan's isolation-window target m/z, not a PSM m/z.  The unmodified half is the peptide
  with the variable modifications removed, taken from db.peptides.csv; its ``#Feature``
  counts distinct Feature Ids *including 0*.  Rows are sorted by best score descending, then
  stripped sequence descending, then modification pattern (see ``_mod_rank``).
* ``db.ptmProfile.csv`` has one row per PTM name: ``# peptides`` = distinct stripped sequences
  in db.peptides.csv carrying it; proteins and sites are counted over db.protein-peptides.csv
  rows with ``Top`` and ``Unique`` both true, the site being ``Start + position - 1``;
  a site is confident when *any* PSM of a supporting peptide is confident there.

Full precision is needed in three places that the 2-decimal / 4-decimal CSV values cannot
always resolve: the choice of the "Mod" representative PSM and the row order (exact score
ties are common because PEAKS scores are quantised), and ``RT shift``.  Pass ``scores`` /
``rts`` (aligned with ``psms`` rows) when available; otherwise the CSV values are used with
the fall-backs described in docs/PTM_TABLES.md.
"""
from __future__ import annotations

import os
import re
from collections import defaultdict
from decimal import ROUND_HALF_EVEN, Decimal

import numpy as np
import pandas as pd

from .peaks_export import canonical_columns
from .seqnorm import mods_of, strip_mods

__all__ = [
    "ASCORE_MIN",
    "POSITIONAL_CONFIDENCE_MIN",
    "build_modified_peptides",
    "build_ptm_profile",
    "read_isolation_targets",
    "peaks_reference_scores",
]

#: a site is confident when AScore >= ASCORE_MIN (evidence: false max 14.63, true min 15.62) ...
ASCORE_MIN = 15.0
#: ... and the PSM's positional confidence at the residue >= this (false max 1, true min 2)
POSITIONAL_CONFIDENCE_MIN = 2
#: AScore printed in db.modifiedPeptide.csv / db.peptides.csv is capped at 1000 (psms: 1001)
ASCORE_CAP = 1000.0

_ASCORE_RE = re.compile(r"^(?:[A-Z]|N-term|C-term)?(\d+):(.+):(-?\d+(?:\.\d+)?)$")
_DEAMIDATION = {("N", "+0.98"), ("Q", "+0.98")}


# --------------------------------------------------------------------------
# small parsers
# --------------------------------------------------------------------------


def _text(v):
    if v is None or v is pd.NA or (isinstance(v, float) and v != v):
        return ""
    return str(v)


def _parse_ascore(s):
    """'M3:Oxidation (M):1001.00;N15:...' -> {position: (name, score)}."""
    out = {}
    for tok in _text(s).split(";"):
        if not tok:
            continue
        m = _ASCORE_RE.match(tok)
        if m:
            out[int(m.group(1))] = (m.group(2), float(m.group(3)))
    return out


def _parse_pc(s):
    return [int(float(x)) for x in _text(s).split()]


def _sites(peptide, ptm):
    """[(position, residue, delta, name)] for the modified sites of a PEAKS peptide string, in
    position order; names come from the ';'-joined PTM column (one entry per site)."""
    ms = mods_of(peptide)
    names = [n for n in _text(ptm).split(";") if n]
    if len(names) != len(ms):
        names = (names + [f"Unknown ({d})" for _, _, d in ms])[: len(ms)]
    return [(p, aa, d, n) for (p, aa, d), n in zip(ms, names)]


def _variable_names(*tables):
    """PTM names that ever carry an AScore (fixed modifications never do)."""
    names = set()
    for t in tables:
        if t is None or "AScore" not in t:
            continue
        for s in t["AScore"].dropna().unique():
            for name, _ in _parse_ascore(s).values():
                names.add(name)
    return names


def _confident(pc, ascore, fixed):
    if fixed:
        return True
    return pc >= POSITIONAL_CONFIDENCE_MIN and ascore is not None and ascore >= ASCORE_MIN


def _site_confidence(psm_row_pc, psm_row_asc, pos, fixed):
    pc = psm_row_pc[pos - 1] if 0 < pos <= len(psm_row_pc) else 0
    asc = psm_row_asc.get(pos, (None, None))[1]
    return _confident(pc, asc, fixed)


def _half_even(x, n):
    """Round the exact binary value of ``x`` HALF_EVEN to n decimals (PEAKS prints the
    isolation m/z this way: float 848.765625 -> 848.76562)."""
    q = Decimal(float(x)).quantize(Decimal(1).scaleb(-n), rounding=ROUND_HALF_EVEN)
    return float(q)


def _dec_diff(a, b, n=4):
    """a - b for values that carry n decimals (avoids binary noise before fixedN formatting)."""
    q = Decimal(1).scaleb(-n)
    return float(Decimal(repr(float(a))).quantize(q) - Decimal(repr(float(b))).quantize(q))


# --------------------------------------------------------------------------
# isolation targets (the m/z columns of db.modifiedPeptide.csv)
# --------------------------------------------------------------------------

_SPEC_RE = re.compile(rb"<spectrum\b[^>]*\bid=\"[^\"]*?scan=(\d+)")
_ISO_RE = re.compile(rb"name=\"isolation window target m/z\"")
_VAL_RE = re.compile(rb"value=\"([^\"]+)\"")


def read_isolation_targets(mzml_path):
    """{raw scan number: isolation window target m/z} for every MS2 spectrum of an mzML file
    (ThermoRawFileParser writes the Thermo float32 target, e.g. ``702.0225219726562``)."""
    out = {}
    scan = None
    with open(mzml_path, "rb") as fh:
        for line in fh:
            if b"<spectrum " in line:
                m = _SPEC_RE.search(line)
                scan = int(m.group(1)) if m else None
            elif scan is not None and b"isolation window target" in line and _ISO_RE.search(line):
                v = _VAL_RE.search(line)
                if v:
                    out[scan] = float(v.group(1))
    return out


def _iso_lookup(isolation_mz):
    if isolation_mz is None:
        return None
    if callable(isolation_mz):
        return isolation_mz

    def f(source_file, scan):
        v = isolation_mz.get((source_file, scan))
        if v is None:
            v = isolation_mz.get(scan)
        if v is None:
            v = isolation_mz.get(int(scan))
        return v
    return f


# --------------------------------------------------------------------------
# representative PSMs
# --------------------------------------------------------------------------


def _psm_index(psms, scores):
    """{peptide: [row positions in psms order]} and the score array used for ranking."""
    by_pep = defaultdict(list)
    for i, p in enumerate(psms["Peptide"].tolist()):
        by_pep[p].append(i)
    if scores is None:
        sc = psms["-10LgP"].to_numpy(dtype=float)
    else:
        sc = np.asarray(pd.Series(scores).to_numpy(dtype=float))
        if len(sc) != len(psms):
            raise ValueError("scores must be aligned with psms rows")
        sc = np.where(np.isnan(sc), psms["-10LgP"].to_numpy(dtype=float), sc)
    return by_pep, sc


def _rep_rows(rows, sc, full, ceiling, pep_row=None):
    """(first, last) psms row among the best-scoring PSMs of one peptide.

    With full-precision scores the tie set is exact.  With 2-decimal scores only, a tie is
    trusted when it is at the sample's score ceiling (PEAKS caps -10lgP; every capped PSM has
    the same float); otherwise the peptides.csv representative (``pep_row``) is used for both."""
    best = max(sc[i] for i in rows)
    ties = [i for i in rows if sc[i] == best]
    if full or best >= ceiling:
        return ties[0], ties[-1]
    first = pep_row if pep_row is not None else ties[0]
    return first, first


def _match_peptide_row(psms, rows, prow):
    """The psms row that db.peptides.csv reports (same Scan, z, m/z, RT), or None."""
    for i in rows:
        r = psms.iloc[i]
        if (r["Scan"] == prow["Scan"] and r["z"] == prow["z"] and r["RT"] == prow["RT"]
                and r["m/z"] == prow["m/z"]):
            return i
    return None


# --------------------------------------------------------------------------
# db.modifiedPeptide.csv
# --------------------------------------------------------------------------


def _mod_rank(peptide, variable):
    """Per-residue rank used to order peptides with the same stripped sequence (descending):
    deamidated N/Q sort below the unmodified residue, other variable modifications above it.
    (7/7 reference cases; equivalent to comparing N(+.98) as D, Q(+.98) as E, M(+15.99) as M~.)"""
    stripped = strip_mods(peptide)
    rank = [0] * len(stripped)
    for pos, aa, delta in mods_of(peptide):
        if not 1 <= pos <= len(stripped):
            continue
        if (aa, delta) in _DEAMIDATION:
            rank[pos - 1] = -1
        else:
            rank[pos - 1] = 1
    return tuple(rank)


def build_modified_peptides(psms, peptides, proteins=None, protein_peptides=None,
                            sample_name="Sample 1", *, scores=None, rts=None, isolation_mz=None,
                            variable_ptms=None):
    """PEAKS ``db.modifiedPeptide.csv``.

    psms, peptides: PEAKS-format tables (``peaks_export`` column names).
    proteins, protein_peptides: accepted for a uniform signature; not used.
    scores: optional full-precision -10lgP per psms row (same order).  Needed to reproduce
        exact ties (representative PSM and row order); falls back to the 2-decimal ``-10LgP``.
    rts: optional full-precision RT (minutes) per psms row, used for ``RT shift``.
    isolation_mz: {(source file, scan) or scan: isolation-window target m/z} or a callable
        (source_file, scan) -> m/z.  See :func:`read_isolation_targets`.  Without it the
        representative PSM's ``raw m/z`` (or ``m/z`` when 0) is written, which is not PEAKS.
    variable_ptms: names of variable modifications; default: every name seen in an AScore.
    """
    S = sample_name
    cols = canonical_columns("modified_peptides", (S,))
    psms = psms.reset_index(drop=True)
    peptides = peptides.reset_index(drop=True)
    variable = set(variable_ptms) if variable_ptms is not None else _variable_names(psms, peptides)
    by_pep, sc = _psm_index(psms, scores)
    full = scores is not None
    ceiling = float(psms["-10LgP"].max()) if len(psms) else np.inf
    rt_arr = None if rts is None else np.asarray(pd.Series(rts).to_numpy(dtype=float))
    iso = _iso_lookup(isolation_mz)
    pep_by_seq = {p: i for i, p in enumerate(peptides["Peptide"].tolist())}
    area_col, nfeat_col = f"Area {S}", f"#Feature {S}"

    def rep_of(pepstr):
        prow = peptides.iloc[pep_by_seq[pepstr]]
        rows = by_pep.get(pepstr, [])
        if not rows:
            return prow, None, None, None
        pep_row = _match_peptide_row(psms, rows, prow)
        first, last = _rep_rows(rows, sc, full, ceiling, pep_row)
        return prow, first, last, max(sc[i] for i in rows)

    def mz_of(source_file, scan, fallback_row):
        v = iso(source_file, scan) if iso is not None else None
        if v is None and fallback_row is not None:
            r = psms.iloc[fallback_row]
            v = r["raw m/z"] if "raw m/z" in psms and r["raw m/z"] else r["m/z"]
        return np.nan if v is None else _half_even(np.float32(v), 5)

    def rt_of(row_i, fallback):
        if rt_arr is not None and row_i is not None and not np.isnan(rt_arr[row_i]):
            return rt_arr[row_i]
        return fallback

    records = []
    for pepstr, ptm in zip(peptides["Peptide"].tolist(), peptides["PTM"].tolist()):
        sites = _sites(pepstr, ptm)
        if not any(n in variable for *_, n in sites):
            continue
        prow, first, last, best = rep_of(pepstr)
        rrow = psms.iloc[last] if last is not None else prow
        asc = _parse_ascore(rrow["AScore"])
        pc = _parse_pc(rrow["Positional Confidence"])
        flags, ascores = [], []
        for pos, aa, delta, name in sites:
            fixed = name not in variable
            ok = _site_confidence(pc, asc, pos, fixed)
            flags.append(f"{aa}{pos}:{name}:{'true' if ok else 'false'}")
        stripped = strip_mods(pepstr)
        for pos in sorted(asc):
            name, s = asc[pos]
            aa = stripped[pos - 1] if 0 < pos <= len(stripped) else ""
            ascores.append(f"{aa}{pos}:{name}:{min(s, ASCORE_CAP):.2f}")
        rec = {
            "Peptide Mod": pepstr,
            "-10LgP Mod": prow["-10LgP"],
            "Mass Mod": prow["Mass"],
            "Length Mod": prow["Length"],
            "PPM Mod": rrow["ppm"],
            "m/z Mod": mz_of(rrow["Source File"], rrow["Scan"], last),
            "z Mod": rrow["z"],
            "RT Mod": rrow["RT"],
            "Scan Mod": rrow["Scan"],
            "Source File Mod": rrow["Source File"],
            f"Area {S} Mod": prow[area_col],
            f"#Feature {S} Mod": prow[nfeat_col],
            "Found By Mod": prow["Found By"],
            "PTM": ";".join(flags),
            "AScore": ";".join(ascores),
            "Accession": prow["Accession"],
            "_score": best if best is not None else prow["-10LgP"],
            "_stripped": stripped,
            "_rank": _mod_rank(pepstr, variable),
        }
        # unmodified counterpart: same residues, variable modifications removed, fixed kept
        um = _unmodified(pepstr, sites, variable)
        if um in pep_by_seq:
            urow_p, ufirst, _ulast, _ = rep_of(um)
            urows = by_pep.get(um, [])
            fids = {int(psms.iloc[i]["Feature Id"]) for i in urows} if "Feature Id" in psms else set()
            rec.update({
                "Peptide": um,
                "-10LgP": urow_p["-10LgP"],
                "Mass": urow_p["Mass"],
                "Length": urow_p["Length"],
                "PPM": urow_p["ppm"],
                "m/z": mz_of(urow_p["Source File"], urow_p["Scan"], ufirst),
                "z": urow_p["z"],
                "RT": urow_p["RT"],
                "Scan": urow_p["Scan"],
                "Source File": urow_p["Source File"],
                f"Area {S}": urow_p[area_col],
                f"#Feature {S}": len(fids) if urows else urow_p[nfeat_col],
                "Found By": urow_p["Found By"],
            })
            ra = rt_of(last, rrow["RT"])
            rb = rt_of(ufirst, urow_p["RT"])
            rec["RT shift"] = (ra - rb) if rt_arr is not None else _dec_diff(ra, rb)
        records.append(rec)

    records.sort(key=lambda r: (r["_score"], r["_stripped"], r["_rank"]), reverse=True)
    df = pd.DataFrame.from_records(records, columns=cols + ["_score", "_stripped", "_rank"])
    df = df[cols]
    for c in cols:
        if c in ("Length Mod", "z Mod", "Scan Mod", f"#Feature {S} Mod", "Length", "z", "Scan",
                 f"#Feature {S}"):
            df[c] = pd.array(df[c].tolist(), dtype="Int64") if len(df) else pd.Series([], dtype="Int64")
        elif c in ("Peptide Mod", "Source File Mod", "Found By Mod", "PTM", "AScore", "Accession",
                   "Peptide", "Source File", "Found By"):
            df[c] = df[c].map(_text).astype(object)
        else:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype(float)
    return df.reset_index(drop=True)


def _unmodified(pepstr, sites, variable):
    """Render the peptide with its variable modifications removed (fixed ones kept)."""
    stripped = strip_mods(pepstr)
    keep = {pos: delta for pos, _aa, delta, name in sites if name not in variable}
    out = []
    if 0 in keep:
        out.append(f"({keep[0]})")
    for k, aa in enumerate(stripped, start=1):
        out.append(aa)
        if k in keep:
            out.append(f"({keep[k]})")
    if len(stripped) + 1 in keep:
        out.append(f"-({keep[len(stripped) + 1]})")
    return "".join(out)


# --------------------------------------------------------------------------
# db.ptmProfile.csv
# --------------------------------------------------------------------------


def _first_start(fasta, protein_peptides):
    from .proteins import read_protein_db
    db = read_protein_db(fasta)
    starts = []
    for acc, pep in zip(protein_peptides["Accession"], protein_peptides["Peptide"]):
        i = db.index.get(acc)
        st = strip_mods(pep)
        starts.append(db.sequence[i].find(st) + 1 if i is not None else 0)
    return starts


def build_ptm_profile(psms, peptides, proteins=None, protein_peptides=None, fasta=None,
                      sample_name="Sample 1", *, variable_ptms=None):
    """PEAKS ``db.ptmProfile.csv`` (the 7-column PTM summary; = summary Table 4).

    Protein positions come from ``protein_peptides['Start']`` (PEAKS's first valid occurrence);
    ``fasta`` is only read when that column is missing, taking the first verbatim occurrence.
    ``proteins`` is accepted for a uniform signature and not used.
    """
    cols = canonical_columns("ptm_profile", (sample_name,))
    variable = set(variable_ptms) if variable_ptms is not None else _variable_names(psms, peptides)

    # per-PTM: stripped sequences, delta, residues
    seqs = defaultdict(set)
    delta_of, residues_of = {}, defaultdict(set)
    for pepstr, ptm in zip(peptides["Peptide"].tolist(), peptides["PTM"].tolist()):
        st = strip_mods(pepstr)
        for pos, aa, delta, name in _sites(pepstr, ptm):
            seqs[name].add(st)
            delta_of.setdefault(name, delta)
            residues_of[name].add(aa)

    # confidence per (peptide, position): any PSM confident there
    conf = defaultdict(bool)
    for pepstr, ptm, a, pcs in zip(psms["Peptide"].tolist(), psms["PTM"].tolist(),
                                   psms["AScore"].tolist(), psms["Positional Confidence"].tolist()):
        asc, pc = _parse_ascore(a), _parse_pc(pcs)
        for pos, aa, delta, name in _sites(pepstr, ptm):
            key = (pepstr, pos)
            if not conf[key]:
                conf[key] = _site_confidence(pc, asc, pos, name not in variable)

    # sites on Top & Unique protein-peptide rows
    prots, sites, conf_sites = defaultdict(set), defaultdict(set), defaultdict(set)
    if protein_peptides is not None and len(protein_peptides):
        pp = protein_peptides
        starts = pp["Start"].tolist() if "Start" in pp else _first_start(fasta, pp)
        ptm_by_pep = dict(zip(peptides["Peptide"], peptides["PTM"]))
        for acc, pepstr, top, uniq, start in zip(pp["Accession"], pp["Peptide"], pp["Top"],
                                                 pp["Unique"], starts):
            if not (_true(top) and _true(uniq)):
                continue
            ptm = ptm_by_pep.get(pepstr, pp.loc[pp["Peptide"] == pepstr, "PTM"].iloc[0])
            for pos, aa, delta, name in _sites(pepstr, ptm):
                site = (acc, int(start) + pos - 1)
                prots[name].add(acc)
                sites[name].add(site)
                fixed = name not in variable
                if fixed or conf.get((pepstr, pos), False):
                    conf_sites[name].add(site)

    rows = []
    for name in seqs:
        m = re.search(r"\(([A-Za-z-]+)\)\s*$", name)
        position = m.group(1) if m else "".join(sorted(residues_of[name]))
        rows.append({
            "Name": name,
            "Mass (Da)": float(delta_of[name]),
            "Position": position,
            "# peptides": len(seqs[name]),
            "# proteins": len(prots[name]),
            "# confident modification sites": len(conf_sites[name]),
            "# modification sites": len(sites[name]),
        })
    rows.sort(key=lambda r: (-r["# peptides"], 0 if r["Name"] in variable else 1, r["Name"]))
    df = pd.DataFrame.from_records(rows, columns=cols)
    for c in ("# peptides", "# proteins", "# confident modification sites", "# modification sites"):
        df[c] = df[c].astype("int64")
    df["Mass (Da)"] = df["Mass (Da)"].astype(float)
    return df


def _true(v):
    if isinstance(v, str):
        return v.strip().lower() == "true"
    return bool(v)


# --------------------------------------------------------------------------
# reference helpers (feeding PEAKS's own exports back in)
# --------------------------------------------------------------------------


def _mgf_entries(mgf_path):
    """{SCANS index: (raw scan, PEPMASS text, charge)} from a PEAKS MGF export (headers only)."""
    out = {}
    pm = ch = rs = None
    with open(mgf_path, "rb") as fh:
        for line in fh:
            c = line[:1]
            if c not in (b"P", b"C", b"R", b"S"):
                continue
            if line.startswith(b"PEPMASS="):
                pm = float(line[8:].split()[0])
            elif line.startswith(b"CHARGE="):
                ch = int(line[7:].strip().rstrip(b"+-"))
            elif line.startswith(b"RAWSCANS="):
                rs = int(line[9:].strip().split(b",")[0])
            elif line.startswith(b"SCANS="):
                out[int(line[6:].strip())] = (rs, pm, ch)
    return out


def peaks_reference_scores(export_dir, psms, cache_dir=None):
    """Full-precision -10lgP for each row of a PEAKS ``db.psms.csv`` (NaN where the mzIdentML
    has no item, i.e. PSMs below the peptide threshold).

    ``db.mzidentml_1_2_0.xml`` holds ``PEAKS:peptideScore`` per SpectrumIdentificationItem;
    ``spectrumID="index=N"`` is MGF entry ``SCANS=N+1``, whose RAWSCANS/PEPMASS/CHARGE identify
    the psms row (Scan, m/z, z) together with the peptide.
    """
    import glob
    import xml.etree.ElementTree as ET

    key_cache = None
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        key_cache = os.path.join(cache_dir, os.path.basename(os.path.normpath(export_dir)).replace(" ", "_")
                                 + ".mzid_scores.tsv")
    if key_cache and os.path.exists(key_cache):
        m = pd.read_csv(key_cache, sep="\t")
    else:
        mgf = sorted(glob.glob(os.path.join(glob.escape(export_dir), "*.mgf")))
        if not mgf:
            raise FileNotFoundError(f"no MGF in {export_dir}")
        entries = _mgf_entries(mgf[0])
        rows = []
        for _ev, el in ET.iterparse(os.path.join(export_dir, "db.mzidentml_1_2_0.xml")):
            if not el.tag.endswith("SpectrumIdentificationResult"):
                continue
            idx = int(el.get("spectrumID").split("=")[1])
            raw, pm, ch = entries.get(idx + 1, (None, None, None))
            for sii in el:
                if not sii.tag.endswith("SpectrumIdentificationItem"):
                    continue
                pep = sii.get("peptide_ref", "")
                pep = pep[len("PEPTIDE_"):] if pep.startswith("PEPTIDE_") else pep
                sc = None
                for c in sii:
                    if c.tag.endswith("cvParam") and c.get("name") == "PEAKS:peptideScore":
                        sc = float(c.get("value"))
                rows.append((pep, raw, pm, ch, sc))
            el.clear()
        m = pd.DataFrame(rows, columns=["Peptide", "Scan", "mz", "z", "score"])
        if key_cache:
            m.to_csv(key_cache, sep="\t", index=False)
    d = {(p, int(s), round(float(mz), 5), int(z)): sc
         for p, s, mz, z, sc in zip(m["Peptide"], m["Scan"], m["mz"], m["z"], m["score"])
         if s == s}
    return pd.Series([d.get((p, int(s), round(float(mz), 5), int(z)), np.nan)
                      for p, s, mz, z in zip(psms["Peptide"], psms["Scan"], psms["m/z"], psms["z"])],
                     index=psms.index, dtype=float)
