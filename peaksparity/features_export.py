"""Build PEAKS' ``<Sample>.db.features.csv`` from our MS1 feature table and PSMs.

Rules (checked against the four DrewR_Adam_91526_DDA_OT reference exports, see
docs/FEATURES_EXPORT.md):

* Only features that carry a DB PSM (any row of db.psms.csv with that Feature Id) or a de novo
  result (ALC >= 50) are listed, sorted by Feature Id; `Feature Id` 0 means "no feature" and
  never appears. A feature with a DB PSM gets no de novo peptide.
* m/z, z, RT, RT Begin, RT End and Area are the feature's; `Quality` is always 0.00.
* `-10LgP` is the maximum over the feature's PSMs (3386/3386 in sample 4). `DB Peptide` is the
  set of peptides of the PSM rows at that maximum, i.e. the I/L variants that
  `apply_protein_mapping` split into rows, `;`-joined (3371/3386; the rest list extra
  non-I/L peptides by a rule we could not pin down). PEAKS' order of the variants looks
  arbitrary (hash order); we sort them.
* `Accession` holds the proteins of the lexicographically first variant only (66/66 multi-variant
  features in samples 1-4), in their psms order, `:`-joined.
* `ppm` is the *uncalibrated* error of the feature m/z against the peptide mass,
  ``((m/z - proton) * z - M) / M * 1e6`` (psms' ppm is after recalibration; the two differ by
  ~8 ppm in sample 4). With the exact peptide mass this lands within 0.1 ppm of PEAKS' value for
  99.8% of DB rows but reproduces the printed digit for only ~74%: PEAKS' value sits ~0.03 ppm
  higher on average for a reason we have not found.
* De novo rows: best-ALC de novo result per feature (99.9% peptide, 100% ALC agreement against
  denovoOnly.csv), its own ppm; `-10LgP` 0.00, `DB Peptide` empty; `ALC (%)` is 0.0 on DB rows.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .peaks_export import canonical_columns

PROTON = 1.007276466


def _pick(df: pd.DataFrame, *names):
    for n in names:
        if n in df.columns:
            return df[n]
    return None


def _normalise_features(features_df: pd.DataFrame, default_source: str | None) -> pd.DataFrame:
    """Accept features.detect_features output (feature_id, mz, charge, rt, rt_begin, rt_end, area)
    or a PEAKS-style table; returns a frame with PEAKS column names."""
    f = features_df
    fid = _pick(f, "Feature Id", "feature_id")
    if fid is None:
        raise ValueError("features_df needs a feature_id / 'Feature Id' column")
    rt = _pick(f, "RT", "rt").astype(float)
    src = _pick(f, "Source File", "source_file")
    if src is None:
        src = pd.Series([default_source or ""] * len(f), index=f.index)
    else:
        src = src.astype(str).str.replace(r"\.mzML$", ".raw", regex=True)
    rb = _pick(f, "RT Begin", "rt_begin")
    re_ = _pick(f, "RT End", "rt_end")
    area = _pick(f, "Area", "area")
    out = pd.DataFrame({
        "Source File": src.to_numpy(),
        "Feature Id": fid.astype(int).to_numpy(),
        "m/z": _pick(f, "m/z", "mz").astype(float).to_numpy(),
        "z": _pick(f, "z", "charge").astype(int).to_numpy(),
        "RT": rt.to_numpy(),
        "RT Begin": (rb.astype(float) if rb is not None else rt).to_numpy(),
        "RT End": (re_.astype(float) if re_ is not None else rt).to_numpy(),
        # PEAKS stores areas as Java floats; keep the float32 value so sums match the file
        "Area": (area.astype(np.float32).astype(float).to_numpy() if area is not None
                 else np.full(len(f), np.nan)),
    })
    return out.drop_duplicates(["Source File", "Feature Id"])


def _ppm(mz, z, mass):
    mz, z, mass = (np.asarray(v, dtype=float) for v in (mz, z, mass))
    with np.errstate(invalid="ignore", divide="ignore"):
        return ((mz - PROTON) * z - mass) / mass * 1e6


def _peptide_mass(peptide: str, fallback: float) -> float:
    """Unrounded monoisotopic mass of a PEAKS-notation peptide (psms' Mass has 4 decimals)."""
    try:
        from pyteomics import mass as pm

        from .proteins import _exact_delta
        from .seqnorm import parse
        s, mods = parse(peptide)
        return float(pm.fast_mass(s) + sum(_exact_delta(d) for _, d in mods))
    except Exception:
        return fallback


def _db_assignments(psms_tbl: pd.DataFrame, default_source: str | None) -> pd.DataFrame:
    """One row per (Source File, Feature Id) with DB Peptide, -10LgP, Mass, Accession."""
    p = psms_tbl
    fid = p["Feature Id"].fillna(0).astype(int)
    p = p.loc[fid > 0].copy()
    if p.empty:
        return pd.DataFrame(columns=["Source File", "Feature Id", "DB Peptide", "-10LgP", "_mass", "_acc"])
    p["Feature Id"] = p["Feature Id"].astype(int)
    if "Source File" not in p.columns:
        p["Source File"] = default_source or ""
    p["_lgp"] = p["-10LgP"].astype(float)
    p["_row"] = np.arange(len(p))
    mx = p.groupby(["Source File", "Feature Id"])["_lgp"].transform("max")
    best = p[p["_lgp"] >= mx - 1e-9]
    rows = []
    for (src, fid_), g in best.groupby(["Source File", "Feature Id"], sort=False):
        peps = sorted(set(g["Peptide"].astype(str)))
        first = g[g["Peptide"].astype(str) == peps[0]].sort_values("_row")
        accs: list[str] = []
        for a in first["Accession"].fillna("").astype(str):
            for x in a.split(";"):
                if x and x not in accs:
                    accs.append(x)
        rows.append({"Source File": src, "Feature Id": fid_, "DB Peptide": ";".join(peps),
                     "-10LgP": float(g["_lgp"].max()),
                     "_mass": _peptide_mass(peps[0], float(first["Mass"].iloc[0])),
                     "_acc": ":".join(accs)})
    return pd.DataFrame(rows)


def _denovo_assignments(denovo_df: pd.DataFrame, default_source: str | None, min_alc: float) -> pd.DataFrame:
    """Best-ALC de novo result per feature. Accepts a PEAKS denovoOnly-style table
    (Feature Id, Peptide, ALC (%), Mass, Accession) or snake_case equivalents
    (feature_id, peptide, alc, mass, accession)."""
    d = denovo_df
    fid = _pick(d, "Feature Id", "feature_id")
    pep = _pick(d, "Peptide", "Denovo Peptide", "peptide")
    alc = _pick(d, "ALC (%)", "alc")
    if fid is None or pep is None or alc is None:
        raise ValueError("denovo_df needs Feature Id, Peptide and ALC (%) columns")
    mass = _pick(d, "Mass", "mass")
    ppm = _pick(d, "ppm")
    acc = _pick(d, "Accession", "accession")
    src = _pick(d, "Source File", "source_file")
    x = pd.DataFrame({
        "Source File": (src.astype(str).str.replace(r"\.mzML$", ".raw", regex=True).to_numpy()
                        if src is not None else default_source or ""),
        "Feature Id": pd.to_numeric(fid, errors="coerce").fillna(0).astype(int).to_numpy(),
        "Denovo Peptide": pep.astype(str).to_numpy(),
        "ALC (%)": alc.astype(float).to_numpy(),
        "_mass": mass.astype(float).to_numpy() if mass is not None else np.nan,
        "_ppm": ppm.astype(float).to_numpy() if ppm is not None else np.nan,
        "_acc": (acc.fillna("").astype(str).str.replace(";", ":").to_numpy() if acc is not None else ""),
    })
    x = x[(x["Feature Id"] > 0) & (x["ALC (%)"] >= min_alc)]
    return (x.sort_values("ALC (%)", ascending=False, kind="stable")
             .drop_duplicates(["Source File", "Feature Id"]))


def build_features_table(features_df: pd.DataFrame, psms_tbl: pd.DataFrame, denovo_df: pd.DataFrame | None = None,
                         sample_name: str = "Sample 1", min_alc: float = 50.0) -> pd.DataFrame:
    """PEAKS ``<Sample>.db.features.csv`` rows.

    features_df: our feature table (features.detect_features / ms1_features.parquet) or a
        PEAKS-style one; may carry a source_file / 'Source File' column (needed for multi-file
        samples, otherwise the PSMs' single source file is used).
    psms_tbl: db.psms.csv-style table (after apply_protein_mapping) with Feature Id, Peptide,
        -10LgP, Mass, Accession, Source File.
    denovo_df: optional de novo results with a Feature Id (see _denovo_assignments).
    sample_name is accepted for symmetry with the other builders (the file name carries it).
    """
    srcs = psms_tbl["Source File"].dropna().unique() if "Source File" in psms_tbl.columns else []
    default_source = str(srcs[0]) if len(srcs) == 1 else None
    feats = _normalise_features(features_df, default_source)
    db = _db_assignments(psms_tbl, default_source)
    cols = canonical_columns("features")
    parts = []
    if len(db):
        m = feats.merge(db, on=["Source File", "Feature Id"], how="inner")
        m["Quality"] = 0.0
        m["Denovo Peptide"] = ""
        m["ALC (%)"] = 0.0
        m["ppm"] = _ppm(m["m/z"], m["z"], m["_mass"])
        m["Accession"] = m["_acc"]
        parts.append(m)
    if denovo_df is not None and len(denovo_df):
        dn = _denovo_assignments(denovo_df, default_source, min_alc)
        if len(db):
            have = set(zip(db["Source File"], db["Feature Id"]))
            dn = dn[[k not in have for k in zip(dn["Source File"], dn["Feature Id"])]]
        m = feats.merge(dn, on=["Source File", "Feature Id"], how="inner")
        m["Quality"] = 0.0
        m["DB Peptide"] = ""
        m["-10LgP"] = 0.0
        # the de novo result's own ppm (against the feature precursor) when given, else recompute
        m["ppm"] = np.where(np.isfinite(m["_ppm"].astype(float)), m["_ppm"],
                            _ppm(m["m/z"], m["z"], m["_mass"]))
        m["Accession"] = m["_acc"]
        parts.append(m)
    if not parts:
        return pd.DataFrame(columns=cols)
    out = pd.concat(parts, ignore_index=True, sort=False)
    out = out.sort_values(["Source File", "Feature Id"], kind="stable").reset_index(drop=True)
    return out[cols]


def features_from_entries(entries: pd.DataFrame, source_file: str | None = None) -> pd.DataFrame:
    """Fallback feature table from precursors.parquet entries (no elution bounds: RT Begin/End =
    apex). Used when the run was given a precomputed feature table that is no longer at hand."""
    e = entries[entries["feature_id"].fillna(0).astype(int) > 0]
    f = e.drop_duplicates("feature_id")[["feature_id", "mz", "charge", "rt", "area"]].copy()
    if source_file is not None:
        f["source_file"] = source_file
    return f


def feature_concordance(ref_features: pd.DataFrame, test_features: pd.DataFrame, ppm_tol: float = 10.0,
                        rt_tol: float = 0.5, db_only: bool = True, ignore_source: bool = False) -> dict:
    """Match features one-to-one with concordance's rule (same source and z, m/z within ppm_tol,
    RT within rt_tol min) and report recall/precision plus area and peptide agreement on the
    matched pairs. db_only restricts both sides to rows with a DB Peptide."""
    from scipy.stats import spearmanr

    from .concordance import _feature_frame, _match_features

    def db(df):
        if not db_only or "DB Peptide" not in df.columns:
            return df.reset_index(drop=True)
        s = df["DB Peptide"].fillna("").astype(str)
        return df[s != ""].reset_index(drop=True)

    R, T = db(ref_features), db(test_features)
    rf, tf = _feature_frame(R, ignore_source), _feature_frame(T, ignore_source)
    pairs = _match_features(rf, tf, ppm_tol, rt_tol)
    ri = np.array([p[0] for p in pairs], dtype=int)
    ti = np.array([p[1] for p in pairs], dtype=int)
    res = {"n_ref": len(rf), "n_test": len(tf), "shared": len(pairs),
           "recall": len(pairs) / len(rf) if len(rf) else None,
           "precision": len(pairs) / len(tf) if len(tf) else None}
    if len(pairs):
        ra, ta = rf["area"].to_numpy()[ri], tf["area"].to_numpy()[ti]
        ok = (ra > 0) & (ta > 0)
        lr = np.log10(ta[ok] / ra[ok])
        res.update(
            n_area_pairs=int(ok.sum()),
            area_spearman_log=float(spearmanr(np.log10(ra[ok]), np.log10(ta[ok])).correlation),
            area_pearson_log=float(np.corrcoef(np.log10(ra[ok]), np.log10(ta[ok]))[0, 1]),
            median_log10_ratio=float(np.median(lr)),
            iqr_log10_ratio=float(np.subtract(*np.percentile(lr, [75, 25]))),
            frac_within_2x=float(np.mean(np.abs(lr) <= np.log10(2))),
            median_ppm_diff=float(np.median([p[2] for p in pairs])),
            median_rt_diff_min=float(np.median([p[3] for p in pairs])),
        )
        if "DB Peptide" in R.columns and "DB Peptide" in T.columns:
            def il(s):
                return {x.replace("I", "L") for x in str(s).split(";") if x}
            rp = R["DB Peptide"].fillna("").to_numpy()[ri]
            tp = T["DB Peptide"].fillna("").to_numpy()[ti]
            both = [(a, b) for a, b in zip(rp, tp) if a and b]
            res["peptide_agree_IL"] = (float(np.mean([bool(il(a) & il(b)) for a, b in both]))
                                       if both else None)
    return res
