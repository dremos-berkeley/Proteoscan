#!/usr/bin/env python
"""Evaluate proteoscan.monoiso.correct_precursor against PEAKS's featureless precursors.

Usage::

    python scripts/eval_monoiso.py prepare --label s4 --mzml <run.mzML> --mgf <PEAKS run.mgf> \
        --out <workdir>
    python scripts/eval_monoiso.py eval --label s4 --out <workdir> [--psms <db.psms.csv>] \
        [--param min_ratio=0.3 ...]
    python scripts/eval_monoiso.py grid --label s4 --out <workdir> --psms <db.psms.csv> [--threads 8]
    python scripts/eval_monoiso.py speed --label s4 --out <workdir>

prepare  reads the mzML (MS1 centroids, MS2 metadata: scan RT, instrument charge, isolation
         target from the filter string, master MS1 from precursor/@spectrumRef) and the PEAKS
         MGF headers, and caches <label>_ms1.npz and <label>_featureless.parquet in --out.
         Featureless entries are MGF entries whose RTINSECONDS is within 0.01 s of the scan RT.
eval     runs correct_precursor(master MS1, target, instrument charge) on every featureless
         entry and compares with PEAKS's PEPMASS (5 ppm) and the isotope shift class k
         (relative to the target peak; k = round((m/z - target) * z / 1.00335)).  With --psms,
         PSMs with Feature Id 0 give a "truth" k from the peptide's calculated mass.
grid     small parameter grid over min_ratio / margin / tol_ppm / heavy-arithmetic threshold.
speed    entries per second of correct_precursor (single thread).
"""
from __future__ import annotations

import argparse
import itertools
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from proteoscan import monoiso  # noqa: E402
from proteoscan.monoiso import C13, PROTON, correct_precursor  # noqa: E402

KS = [-4, -3, -2, -1, 0]
MBINS = [0, 1000, 1500, 2000, 2500, 3000, 4000, 6000, 20000]


# --------------------------------------------------------------------------------------
# prepare
# --------------------------------------------------------------------------------------

def read_mgf_headers(path) -> pd.DataFrame:
    """(entry, rawscan, pepmass, charge, rt_s) from a PEAKS MGF (header lines only, via grep)."""
    out = subprocess.run(["grep", "-E", "^(BEGIN IONS|PEPMASS=|CHARGE=|RAWSCANS=|RTINSECONDS=)", str(path)],
                         check=True, capture_output=True, text=True).stdout.splitlines()
    rows, cur = [], None
    for line in out:
        if line.startswith("BEGIN"):
            if cur is not None:
                rows.append(cur)
            cur = {}
        else:
            k, v = line.split("=", 1)
            cur[k] = v
    if cur is not None:
        rows.append(cur)
    d = pd.DataFrame(rows)
    return pd.DataFrame({
        "entry": np.arange(len(d)),
        "rawscan": d["RAWSCANS"].astype(int).to_numpy(),
        "pepmass": d["PEPMASS"].str.split().str[0].astype(float).to_numpy(),
        "charge": d["CHARGE"].str.rstrip("+-").astype(int).to_numpy(),
        "rt_s": d["RTINSECONDS"].astype(float).to_numpy(),
    })


def read_mzml(path):
    """MS1 centroids (flat arrays + per-scan offsets) and MS2 metadata."""
    import pyopenms as oms

    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(path), exp)
    scanre = re.compile(r"scan=(\d+)")
    filtre = re.compile(r"ms2 ([\d.]+)@")
    ms1_scan, mzs, ints, rows = [], [], [], []
    for s in exp:
        nid = s.getNativeID()
        nid = nid.decode() if isinstance(nid, bytes) else nid
        scan = int(scanre.search(nid).group(1))
        if s.getMSLevel() == 1:
            m, i = s.get_peaks()
            o = np.argsort(m, kind="stable")
            ms1_scan.append(scan)
            mzs.append(np.asarray(m, np.float64)[o])
            ints.append(np.asarray(i, np.float32)[o])
            continue
        p = s.getPrecursors()[0]
        fs = s.getMetaValue("filter string") if s.metaValueExists("filter string") else None
        if fs is None and s.getAcquisitionInfo().size():
            a = s.getAcquisitionInfo()[0]
            fs = a.getMetaValue("filter string") if a.metaValueExists("filter string") else None
        fs = fs.decode() if isinstance(fs, bytes) else fs
        ref = p.getMetaValue("spectrum_ref") if p.metaValueExists("spectrum_ref") else ""
        ref = ref.decode() if isinstance(ref, bytes) else ref
        mref = scanre.search(ref or "")
        mf = filtre.search(fs or "")
        rows.append((scan, s.getRT(), p.getMZ(), p.getCharge(), float(mf.group(1)) if mf else p.getMZ(),
                     int(mref.group(1)) if mref else -1, p.getIsolationWindowLowerOffset(),
                     p.getIsolationWindowUpperOffset()))
    ms2 = pd.DataFrame(rows, columns=["rawscan", "scan_rt", "prec_mz", "prec_z", "target_mz", "master",
                                      "iso_lo", "iso_hi"])
    offs = np.concatenate([[0], np.cumsum([len(m) for m in mzs])])
    ms1 = {"scan": np.asarray(ms1_scan), "offs": offs, "mz": np.concatenate(mzs), "int": np.concatenate(ints)}
    return ms1, ms2


def cmd_prepare(a):
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    g = read_mgf_headers(a.mgf)
    print(f"MGF: {len(g)} entries ({time.time() - t0:.0f} s)", flush=True)
    ms1, ms2 = read_mzml(a.mzml)
    print(f"mzML: {len(ms1['scan'])} MS1, {len(ms2)} MS2 ({time.time() - t0:.0f} s)", flush=True)
    np.savez(out / f"{a.label}_ms1.npz", **ms1)
    d = g.merge(ms2, on="rawscan", how="left")
    d["featureless"] = (d.rt_s - d.scan_rt).abs() <= 0.01
    d.to_parquet(out / f"{a.label}_entries.parquet", index=False)
    fl = d[d.featureless].reset_index(drop=True)
    fl.to_parquet(out / f"{a.label}_featureless.parquet", index=False)
    print(f"featureless: {len(fl)} of {len(d)}; master MS1 missing: {(fl.master < 0).sum()}")


# --------------------------------------------------------------------------------------
# evaluation
# --------------------------------------------------------------------------------------

_MS1 = None
_FL = None


def _load(out, label):
    global _MS1, _FL
    z = np.load(Path(out) / f"{label}_ms1.npz")
    ms1 = {k: z[k] for k in z.files}
    ms1["pos"] = {int(s): i for i, s in enumerate(ms1["scan"])}
    _MS1 = ms1
    _FL = pd.read_parquet(Path(out) / f"{label}_featureless.parquet")
    monoiso._envelope_table()   # build once before worker processes fork
    return ms1, _FL


def _spec(ms1, scan):
    i = ms1["pos"][int(scan)]
    lo, hi = ms1["offs"][i], ms1["offs"][i + 1]
    return ms1["mz"][lo:hi], ms1["int"][lo:hi]


_IMPL = correct_precursor


def use_module(path):
    """Evaluate correct_precursor from another file (e.g. an older monoiso.py) instead."""
    global _IMPL
    import importlib.util
    spec = importlib.util.spec_from_file_location("monoiso_alt", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["monoiso_alt"] = mod
    spec.loader.exec_module(mod)
    _IMPL = mod.correct_precursor


def run_correction(ms1, fl, params: dict, rows=None) -> pd.DataFrame:
    """Apply correct_precursor to featureless entries; returns our m/z, k, source."""
    idx = np.arange(len(fl)) if rows is None else rows
    master = fl.master.to_numpy()
    tgt = fl.target_mz.to_numpy()
    z = fl.prec_z.to_numpy()
    res_mz = np.empty(len(idx))
    res_k = np.empty(len(idx), int)
    res_src = np.empty(len(idx), object)
    for n, i in enumerate(idx):
        m, it = _spec(ms1, master[i])
        c = _IMPL(m, it, float(tgt[i]), int(z[i]), **params)
        res_mz[n], res_k[n], res_src[n] = c.mz, c.shift, c.source
    return pd.DataFrame({"row": idx, "our_mz": res_mz, "our_shift": res_k, "our_src": res_src})


def _worker(args):
    params, rows = args
    return run_correction(_MS1, _FL, params, rows)


def run_parallel(params, threads, n):
    chunks = np.array_split(np.arange(n), threads * 4)
    with ProcessPoolExecutor(threads) as ex:
        parts = list(ex.map(_worker, [(params, c) for c in chunks]))
    return pd.concat(parts, ignore_index=True)


def truth_table(psms_path, fl) -> pd.DataFrame:
    """PSMs with Feature Id 0 matched to featureless entries by (scan, m/z); truth mono m/z.

    Biased towards PEAKS: a PSM exists only where PEAKS's PEPMASS matched the peptide."""
    p = pd.read_csv(psms_path)
    p = p[p["Feature Id"] == 0]
    t = pd.DataFrame({"rawscan": p["Scan"].astype(int), "psm_mz": p["m/z"].astype(float),
                      "psm_z": p["z"].astype(int), "calc_mass": p["Mass"].astype(float)})
    t["key"] = (t.psm_mz * 1e4).round().astype(np.int64)
    f = fl[["rawscan", "pepmass"]].reset_index().rename(columns={"index": "row"})
    f["key"] = (f.pepmass * 1e4).round().astype(np.int64)
    t = t.merge(f[["row", "rawscan", "key"]], on=["rawscan", "key"], how="inner")
    t = t.drop_duplicates("row")
    t["true_mz"] = t.calc_mass / t.psm_z + PROTON
    return t[["row", "psm_z", "calc_mass", "true_mz"]]


def anchors(ms1, fl, tol_ppm=10.0) -> np.ndarray:
    """m/z of the master-MS1 centroid nearest the isolation target (nan if none within tol)."""
    out = np.full(len(fl), np.nan)
    for n, (sc, t) in enumerate(zip(fl.master.to_numpy(), fl.target_mz.to_numpy())):
        m, _ = _spec(ms1, sc)
        i = np.searchsorted(m, t)
        c = [j for j in (i - 1, i) if 0 <= j < len(m)]
        if c:
            j = min(c, key=lambda q: abs(m[q] - t))
            if abs(m[j] - t) <= t * tol_ppm * 1e-6:
                out[n] = m[j]
    return out


def truth_other_scans(psms_path, fl, ms1, ppm=3.0, rt_min=1.0, max_j=3, control_ppm=0.0) -> pd.DataFrame:
    """Less biased truth: a peptide identified in ANOTHER scan (any PSM, with or without feature)
    whose isotope j (0..max_j), at the entry's charge, matches the target peak's centroid (the
    anchor) within `ppm` of the run's median raw PSM mass error, eluting within `rt_min` minutes
    of the scan.  Then the true mono is j isotopes below the target peak, whatever PEAKS chose
    for this scan.  Entries with conflicting candidates are dropped.  Random matches are not
    excluded; control_ppm != 0 moves the window off the calibration to count them."""
    p = pd.read_csv(psms_path)
    # raw MS1 mass offset of identified precursors (PEAKS's own "ppm" column is recalibrated)
    obs = (p["m/z"] - PROTON) * p["z"]
    d = obs - p["Mass"]
    d = d - np.round(d / C13) * C13
    offset = float(np.median(d / p["Mass"] * 1e6)) + control_ppm
    mass = p["Mass"].to_numpy(float)
    z = p["z"].to_numpy(int)
    rt = p["raw RT"].to_numpy(float)
    rt2 = p["RT"].to_numpy(float)
    scan = p["Scan"].to_numpy(int)
    anc = anchors(ms1, fl)
    rows = []
    fz = fl.prec_z.to_numpy()
    frt = fl.scan_rt.to_numpy() / 60.0
    fscan = fl.rawscan.to_numpy()
    for zz in np.unique(fz):
        sel = np.where(z == zz)[0]
        if not len(sel):
            continue
        iso = (mass[sel, None] / zz + PROTON + np.arange(max_j + 1)[None, :] * C13 / zz).ravel()
        iso = iso * (1 + offset * 1e-6)
        jj = np.tile(np.arange(max_j + 1), len(sel))
        src = np.repeat(sel, max_j + 1)
        o = np.argsort(iso)
        iso, jj, src = iso[o], jj[o], src[o]
        ents = np.where((fz == zz) & np.isfinite(anc))[0]
        lo = np.searchsorted(iso, anc[ents] * (1 - ppm * 1e-6))
        hi = np.searchsorted(iso, anc[ents] * (1 + ppm * 1e-6))
        for e, a, b in zip(ents, lo, hi):
            if a == b:
                continue
            c = np.arange(a, b)
            s = src[c]
            keep = (scan[s] != fscan[e]) & ((np.abs(rt[s] - frt[e]) <= rt_min) | (np.abs(rt2[s] - frt[e]) <= rt_min))
            if not keep.any():
                continue
            js = np.unique(jj[c[keep]])
            if len(js) == 1:
                rows.append((e, int(js[0])))
    t = pd.DataFrame(rows, columns=["row", "j"])
    t["true_k"] = -t.j
    return t[["row", "true_k"]]


def annotate(fl, res) -> pd.DataFrame:
    d = fl.iloc[res.row.to_numpy()].reset_index(drop=True).copy()
    for c in ("row", "our_mz", "our_shift", "our_src"):
        d[c] = res[c].to_numpy()
    z = d.prec_z.to_numpy()
    d["peaks_k"] = np.round((d.pepmass - d.target_mz) * z / C13).astype(int)
    d["our_k"] = np.round((d.our_mz - d.target_mz) * z / C13).astype(int)
    d["mz_agree"] = (d.our_mz - d.pepmass).abs() <= d.pepmass * 5e-6
    d["k_agree"] = d.our_k == d.peaks_k
    d["M"] = (d.target_mz - PROTON) * z
    d["Mbin"] = pd.cut(d.M, MBINS)
    return d


def _with_true_k(d, truth):
    t = d.merge(truth, on="row")
    if "true_k" not in t:
        t["true_k"] = np.round((t.true_mz - t.target_mz) * t.prec_z.to_numpy() / C13).astype(int)
    return t


def summarize(d, truth=None, truth2=None) -> dict:
    s = {"n": len(d), "mz_agree": d.mz_agree.mean(), "k_agree": d.k_agree.mean()}
    for tag, tr in (("psm", truth), ("other", truth2)):
        if tr is not None and len(tr):
            t = _with_true_k(d, tr)
            s[f"{tag}_n"] = len(t)
            s[f"{tag}_peaks"] = (t.peaks_k == t.true_k).mean()
            s[f"{tag}_ours"] = (t.our_k == t.true_k).mean()
            if tag == "other":
                tt = t[t.true_k >= -1]
                s["other01_peaks"] = (tt.peaks_k == tt.true_k).mean()
                s["other01_ours"] = (tt.our_k == tt.true_k).mean()
    return s


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.3f}"
    return str(v)


def md(df: pd.DataFrame) -> str:
    cols = [str(c) for c in df.columns]
    lines = ["| " + str(df.index.name or "") + " | " + " | ".join(cols) + " |", "|" + "---|" * (len(cols) + 1)]
    for i, row in zip(df.index, df.itertuples(index=False)):
        lines.append("| " + str(i) + " | " + " | ".join(_fmt(v) for v in row) + " |")
    return "\n".join(lines)


def report_tables(d, truth=None, truth2=None) -> str:
    o = []
    kk = lambda x: x.clip(-4, 1)  # noqa: E731
    o.append(f"entries: {len(d)}; m/z agreement (5 ppm) {d.mz_agree.mean():.4f}; "
             f"k agreement {d.k_agree.mean():.4f}")
    o.append("\nConfusion, rows = ours k, cols = PEAKS k (relative to the target peak):\n")
    o.append(md(pd.crosstab(kk(d.our_k).rename("ours"), kk(d.peaks_k).rename("PEAKS"), margins=True)))
    g = d.groupby("Mbin", observed=False)
    t = pd.DataFrame({"n": g.size(), "k agree": g.k_agree.mean(), "m/z agree": g.mz_agree.mean()})
    for k in KS:
        t[f"PEAKS {k}"] = g.peaks_k.apply(lambda x, k=k: (x.clip(-4, 1) == k).mean())
    for k in KS:
        t[f"ours {k}"] = g.our_k.apply(lambda x, k=k: (x.clip(-4, 1) == k).mean())
    o.append("\nBy neutral mass (from the target m/z), fractions:\n")
    o.append(md(t))
    o.append("\nm/z agreement by PEAKS k (rows) and our source:\n")
    o.append(md(pd.crosstab(kk(d.peaks_k), d.our_src, values=d.mz_agree, aggfunc="mean").round(3)
                .join(kk(d.peaks_k).value_counts().rename("n"))))
    for title, tr in (("PSM truth: Feature Id 0 PSMs on the entry itself (biased towards PEAKS: a PSM "
                       "exists only where PEAKS's PEPMASS matched the peptide)", truth),
                      ("Other-scan truth: a peptide identified in another scan within 1 min whose isotope j <= 3 "
                       "matches the target-peak centroid (3 ppm around the median raw PSM error, same z)", truth2)):
        if tr is None or not len(tr):
            continue
        t = _with_true_k(d, tr)
        o.append(f"\n{title}: n = {len(t)}; PEAKS right {np.mean(t.peaks_k == t.true_k):.4f}; "
                 f"ours right {np.mean(t.our_k == t.true_k):.4f}\n")
        if "control" in tr.attrs:
            n_k = t.true_k.value_counts()
            o.append("Random matches expected (same window moved +-12 ppm off the calibration), by truth k: "
                     + ", ".join(f"k={k}: {int(v)} of {int(n_k.get(k, 0))}" for k, v in sorted(tr.attrs["control"].items())) + "\n")
            for kmax in (0, -1):
                tt = t[t.true_k >= kmax]
                o.append(f"  restricted to truth k >= {kmax}: n = {len(tt)}; PEAKS right "
                         f"{np.mean(tt.peaks_k == tt.true_k):.4f}; ours right {np.mean(tt.our_k == tt.true_k):.4f}\n")
        o.append(md(pd.crosstab(kk(t.our_k).rename("ours"), kk(t.true_k).rename("truth"), margins=True)))
        o.append("\nPEAKS k vs truth k:\n")
        o.append(md(pd.crosstab(kk(t.peaks_k).rename("PEAKS"), kk(t.true_k).rename("truth"), margins=True)))
    return "\n".join(o)


_INT_PARAMS = {"max_shift", "min_peaks"}


def _cast(k, v):
    if k in _INT_PARAMS:
        return int(v)
    return float(v)


def parse_params(items):
    p = {}
    for it in items or []:
        k, v = it.split("=", 1)
        p[k] = _cast(k, v)
    return p


def _truths(a, fl, ms1):
    if not a.psms:
        return None, None
    t2 = truth_other_scans(a.psms, fl, ms1)
    ctl = pd.concat([truth_other_scans(a.psms, fl, ms1, control_ppm=c) for c in (-12.0, 12.0)])
    t2.attrs["control"] = (ctl.true_k.value_counts() / 2).round(0).to_dict()
    return truth_table(a.psms, fl), t2


def cmd_eval(a):
    ms1, fl = _load(a.out, a.label)
    params = parse_params(a.param)
    t0 = time.time()
    res = run_parallel(params, a.threads, len(fl))
    print(f"[{a.label}] {len(fl)} entries in {time.time() - t0:.1f} s on {a.threads} workers; params {params}")
    d = annotate(fl, res)
    truth, truth2 = _truths(a, fl, ms1)
    txt = report_tables(d, truth, truth2)
    print(txt)
    if a.save:
        d.drop(columns=["Mbin"]).to_parquet(Path(a.out) / f"{a.label}_eval.parquet", index=False)
        (Path(a.out) / f"{a.label}_eval.md").write_text(txt + "\n")


GRID = {
    "min_ratio": [0.0, 0.05, 0.1, 0.2, 0.3],
    "margin": [0.0, 0.03, 0.06],
    "tol_da": [0.009, 0.01, 0.011],
    "arith_min_mass": [0.0, 2000.0, 1e9],
}


def cmd_grid(a):
    ms1, fl = _load(a.out, a.label)
    truth, truth2 = _truths(a, fl, ms1)
    grid = json_grid(a.grid) if a.grid else GRID
    keys = list(grid)
    out = []
    for vals in itertools.product(*[grid[k] for k in keys]):
        params = {k: _cast(k, v) for k, v in zip(keys, vals)}
        res = run_parallel(params, a.threads, len(fl))
        d = annotate(fl, res)
        s = summarize(d, truth, truth2)
        s.update(params)
        out.append(s)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()}, flush=True)
    r = pd.DataFrame(out).sort_values("mz_agree", ascending=False)
    r.to_csv(Path(a.out) / f"{a.label}_grid.csv", index=False)
    print(r.head(20).to_string())


def json_grid(s):
    import json
    return json.loads(s)


def cmd_speed(a):
    ms1, fl = _load(a.out, a.label)
    params = parse_params(a.param)
    n = min(a.sample, len(fl))
    t0 = time.perf_counter()
    run_correction(ms1, fl, params, np.arange(n))
    dt = time.perf_counter() - t0
    print(f"[{a.label}] {n} entries in {dt:.2f} s single-threaded: {n / dt:.0f} entries/s "
          f"(~{200000 / (n / dt):.0f} s per 200k scans)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["prepare", "eval", "grid", "speed"])
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mzml")
    ap.add_argument("--mgf")
    ap.add_argument("--psms")
    ap.add_argument("--param", nargs="*", help="correct_precursor keyword overrides, key=value")
    ap.add_argument("--grid", help="JSON dict of parameter lists (default: built-in grid)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--sample", type=int, default=30000, help="speed: number of entries timed")
    ap.add_argument("--save", action="store_true")
    ap.add_argument("--module", help="evaluate correct_precursor from this .py file instead of proteoscan.monoiso")
    a = ap.parse_args(argv)
    if a.module:
        use_module(a.module)
    {"prepare": cmd_prepare, "eval": cmd_eval, "grid": cmd_grid, "speed": cmd_speed}[a.cmd](a)


if __name__ == "__main__":
    main()
