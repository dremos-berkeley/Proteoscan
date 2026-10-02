#!/usr/bin/env python
"""Check a proteoscan installation (macOS or Linux) and run a small end-to-end smoke test.

    python scripts/check_install.py                 # checks + smoke test
    python scripts/check_install.py --skip-smoke    # tools and imports only
    python scripts/check_install.py --work DIR --keep

1. Platform: OS, CPU architecture, Rosetta translation, RAM.
2. Python imports (required and optional) with versions; OpenMP runtimes loaded in-process
   (two different libomp copies in one process abort with "OMP: Error #15" on macOS).
3. External tools: comet, ThermoRawFileParser (plus sage, percolator, mono, dotnet if present).
   Executed, not just located, so a wrong-architecture or unsigned binary is caught here.
4. Smoke test on synthetic data generated here (tiny FASTA + mzML written with pyopenms):
   the real `proteoscan.pipeline.run` path, i.e. MS1 features -> precursor list ->
   Comet -> read_comet -> MS2PIP + DeepLC features -> mokapot -> PEAKS-style filter ->
   protein inference -> CSV export, with per-stage timings, then checks the export.

The first run downloads the MS2PIP HCD2021 model (about 955 MB) into ~/.ms2pip; later runs
take well under a minute. Exit status 0 = everything required works.
"""
from __future__ import annotations

import argparse
import ctypes
import importlib
import inspect
import json
import multiprocessing as mp
import os
import platform
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:  # allow running from a checkout without `pip install -e .`
    sys.path.insert(0, str(REPO))

OK, WARN, FAIL = "ok", "WARN", "FAIL"
_results: list[tuple[str, str, str]] = []


def report(status: str, what: str, detail: str = "") -> None:
    _results.append((status, what, detail))
    print(f"  [{status:4s}] {what:34s} {detail}", flush=True)


# --------------------------------------------------------------------------------------
# 1. platform
# --------------------------------------------------------------------------------------

def _sysctl(name: str) -> str | None:
    try:
        return subprocess.run(["sysctl", "-n", name], capture_output=True, text=True, timeout=5).stdout.strip()
    except Exception:
        return None


def total_ram_gb() -> float | None:
    if sys.platform == "darwin":
        v = _sysctl("hw.memsize")
        return int(v) / 2**30 if v and v.isdigit() else None
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2**30
    except (ValueError, OSError, AttributeError):
        return None


def peak_rss_gb(who: int) -> float:
    import resource
    r = resource.getrusage(who).ru_maxrss
    return r / 2**30 if sys.platform == "darwin" else r / 2**20  # bytes on macOS, KiB on Linux


def check_platform() -> None:
    print("\n== Platform")
    mach = platform.machine()
    info = f"{platform.system()} {platform.release()} {mach}"
    if sys.platform == "darwin":
        info = f"macOS {platform.mac_ver()[0]} {mach}"
        if _sysctl("sysctl.proc_translated") == "1":
            report(WARN, "Rosetta 2", "this Python is x86_64 running under Rosetta on Apple Silicon "
                   "(expected only for an osx-64 env; see docs/INSTALL_MACOS.md)")
        brand = _sysctl("machdep.cpu.brand_string")
        if brand:
            info += f" ({brand})"
    report(OK, "system", info)
    report(OK, "python", f"{sys.version.split()[0]} {sys.executable}")
    ram = total_ram_gb()
    if ram is not None:
        st = OK if ram >= 32 else WARN
        report(st, "RAM", f"{ram:.0f} GB" + ("" if ram >= 32 else
               " (a 3 GB raw file currently peaks at ~34 GB in the MS2PIP stage and will swap; "
               "see docs/INSTALL_MACOS.md, 'Run time and memory')"))
    report(OK, "CPUs", str(os.cpu_count()))
    report(OK, "multiprocessing start method", mp.get_start_method(allow_none=True) or "(default)")


# --------------------------------------------------------------------------------------
# 2. imports
# --------------------------------------------------------------------------------------

REQUIRED = [
    ("numpy", "numpy"), ("pandas", "pandas"), ("scipy", "scipy"), ("scikit-learn", "sklearn"),
    ("lightgbm", "lightgbm"), ("xgboost", "xgboost"), ("pyarrow", "pyarrow"), ("numba", "numba"),
    ("pyopenms", "pyopenms"), ("pyteomics", "pyteomics"), ("psm-utils", "psm_utils"),
    ("ms2pip", "ms2pip"), ("deeplc", "deeplc"), ("mokapot", "mokapot"), ("torch", "torch"),
    ("pyyaml", "yaml"), ("click", "click"), ("pyahocorasick", "ahocorasick"),
    ("proteoscan", "proteoscan"),
]
OPTIONAL = [("pytest", "pytest"), ("casanovo", "casanovo"), ("peptdeep", "peptdeep"),
            ("matplotlib", "matplotlib"), ("polars", "polars")]


def _version(dist: str, mod) -> str:
    try:
        from importlib.metadata import version
        return version(dist)
    except Exception:
        return str(getattr(mod, "__version__", "?"))


def check_imports() -> bool:
    print("\n== Python packages")
    import warnings
    ok = True
    for dist, name in REQUIRED + OPTIONAL:
        required = (dist, name) in REQUIRED
        t = time.time()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mod = importlib.import_module(name)
            extra = ""
            if name == "torch":
                mps = getattr(getattr(mod, "backends", None), "mps", None)
                extra = (f" cuda={mod.cuda.is_available()} mps={bool(mps and mps.is_available())} "
                         f"threads={mod.get_num_threads()}")
            if name == "proteoscan":
                extra = f" ({Path(mod.__file__).parent})"
            report(OK, dist, f"{_version(dist, mod)}{extra}  [{time.time() - t:.1f}s]")
        except Exception as e:  # ImportError, OSError from a broken dylib, ...
            msg = f"{type(e).__name__}: {str(e).splitlines()[0][:150] if str(e) else ''}"
            report(FAIL if required else WARN, dist, msg if required else f"optional, not usable ({msg})")
            ok &= not required
    check_openmp()
    return ok


def loaded_libraries() -> list[str]:
    if sys.platform == "darwin":
        libc = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
        libc._dyld_image_count.restype = ctypes.c_uint32
        libc._dyld_get_image_name.restype = ctypes.c_char_p
        libc._dyld_get_image_name.argtypes = [ctypes.c_uint32]
        return [libc._dyld_get_image_name(i).decode() for i in range(libc._dyld_image_count())]
    maps = Path("/proc/self/maps")
    if maps.exists():
        return sorted({ln.split()[-1] for ln in maps.read_text().splitlines() if "/" in ln})
    return []


def check_openmp() -> None:
    libs = [p for p in loaded_libraries() if re.search(r"lib(g|i|)omp[0-9.]*\.(dylib|so)", os.path.basename(p))]
    real = sorted({os.path.realpath(p) for p in libs})
    names = ", ".join(real) if real else "none loaded"
    llvm = [p for p in real if re.search(r"libomp[0-9.]*\.(dylib|so)", os.path.basename(p))]
    if len(llvm) > 1:
        report(WARN, "OpenMP runtimes", f"{len(llvm)} different libomp copies: {names}. On macOS this can abort with "
               "'OMP: Error #15'; install torch/lightgbm/xgboost/scikit-learn from conda-forge "
               "(env/environment-macos.yml) or, as a last resort, export KMP_DUPLICATE_LIB_OK=TRUE")
    else:
        report(OK, "OpenMP runtimes", names)


# --------------------------------------------------------------------------------------
# 3. external tools
# --------------------------------------------------------------------------------------

def macho_arch(path: str) -> str | None:
    """CPU type of a Mach-O / ELF binary, or None for scripts."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(20)
    except OSError:
        return None
    if head[:4] == b"\xcf\xfa\xed\xfe":
        cpu = struct.unpack_from("<i", head, 4)[0]
        return {0x100000C: "arm64", 0x1000007: "x86_64"}.get(cpu, hex(cpu))
    if head[:4] in (b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca"):
        return "universal"
    if head[:4] == b"\x7fELF":
        em = struct.unpack_from("<H", head, 18)[0]
        return {0x3E: "x86_64", 0xB7: "aarch64"}.get(em, hex(em))
    return None


def run_tool(label: str, cmd: list[str], pattern: str, required: bool, hint: str = "") -> bool:
    exe = shutil.which(cmd[0])
    if not exe:
        report(FAIL if required else WARN, label, ("not on PATH. " if required else "optional, not installed. ") + hint)
        return not required
    arch = macho_arch(os.path.realpath(exe))
    where = f"{exe}" + (f" [{arch}]" if arch else "")
    try:
        p = subprocess.run([exe] + cmd[1:], capture_output=True, text=True, timeout=60)
        out = (p.stdout or "") + (p.stderr or "")
    except OSError as e:
        why = str(e)
        if getattr(e, "errno", None) == 86 or "Bad CPU type" in why:
            why += " -> binary is for another CPU; on Intel Macs Comet 2026.x must be built from source (docs/INSTALL_MACOS.md)"
        report(FAIL if required else WARN, label, f"{where}: cannot execute: {why}")
        return not required
    except subprocess.TimeoutExpired:
        report(FAIL if required else WARN, label, f"{where}: timed out")
        return not required
    if p.returncode in (-9, 137):
        report(FAIL if required else WARN, label, f"{where}: killed by the OS (SIGKILL). On Apple Silicon this means an "
               "unsigned binary: run `codesign --force --sign - <binary>` (scripts/install_trfp.sh does this)")
        return not required
    m = re.search(pattern, out)
    if not m:
        first = out.strip().splitlines()[:2]
        report(FAIL if required else WARN, label, f"{where}: unexpected output (rc={p.returncode}): {first}")
        return not required
    report(OK, label, f"{m.group(1).strip()}  {where}")
    return True


def check_tools() -> bool:
    print("\n== External tools")
    ok = run_tool("comet", ["comet"], r"Comet version\s+\"?([^\"\n]+)", True,
                  "conda install -c bioconda comet-ms=2026011 (Apple Silicon) or build from source (Intel)")
    trfp = shutil.which("ThermoRawFileParser") and "ThermoRawFileParser" or "thermorawfileparser"
    ok &= run_tool("ThermoRawFileParser", [trfp, "--version"], r"(?m)^\s*([0-9]+\.[0-9]+\.[0-9]+\S*)\s*$",
                   False, "needed only for .raw input: bash scripts/install_trfp.sh")
    run_tool("sage", ["sage", "--version"], r"(\d+\.\d+\.\d+)", False, "optional search engine")
    run_tool("percolator", ["percolator", "--help"], r"Percolator version ([^\s,]+)", False,
             "optional; mokapot is used instead")
    run_tool("dotnet", ["dotnet", "--list-runtimes"], r"(Microsoft\.NETCore\.App \S+)", False,
             "only for the framework-dependent ThermoRawFileParser route")
    run_tool("mono", ["mono", "--version"], r"version (\S+)", False, "only for ThermoRawFileParser 1.4.x")
    return ok


# --------------------------------------------------------------------------------------
# 4. synthetic data
# --------------------------------------------------------------------------------------

AA_FREQ = {"A": 7.0, "R": 5.6, "N": 3.6, "D": 4.7, "C": 2.3, "Q": 4.7, "E": 7.1, "G": 6.6, "H": 2.6,
           "I": 4.3, "L": 10.0, "K": 5.7, "M": 2.1, "F": 3.7, "P": 6.3, "S": 8.3, "T": 5.3, "W": 1.2,
           "Y": 2.7, "V": 6.0}
HYDRO = {"A": 1.1, "R": -0.4, "N": -0.9, "D": -0.8, "C": 0.5, "Q": -0.7, "E": -0.6, "G": -0.2, "H": -0.6,
         "I": 3.5, "L": 3.8, "K": -1.1, "M": 2.2, "F": 4.1, "P": 0.2, "S": -0.4, "T": 0.0, "W": 4.5,
         "Y": 2.0, "V": 2.6}  # crude retention coefficients; DeepLC is calibrated on them
PROTON = 1.00727646688
H2O = 18.0105646837
C13 = 1.0033548378
CAM = 57.021464


def _residue_masses():
    from pyteomics import mass
    m = {a: mass.std_aa_mass[a] for a in AA_FREQ}
    m["C"] += CAM  # fixed carbamidomethylation
    return m


def _tryptic(seq: str) -> list[str]:
    return [p for p in re.sub(r"(?<=[KR])(?!P)", "\n", seq).split("\n") if p]


def make_synthetic(work: Path, n_proteins: int = 80, n_planted: int = 450, n_foreign: int = 250, seed: int = 7):
    """Write synthetic.fasta and synthetic.mzML; return the planted peptide set.

    mokapot needs well over 100 target PSMs per cross-validation fold to reach 1% FDR
    (its q-values count one extra decoy), hence several hundred planted peptides."""
    import numpy as np
    import pyopenms as oms

    rng = np.random.default_rng(seed)
    aas = np.array(list(AA_FREQ))
    p = np.array(list(AA_FREQ.values()))
    p /= p.sum()

    def protein(n):
        return "M" + "".join(rng.choice(aas, size=n - 1, p=p))

    prots = [protein(int(rng.integers(300, 700))) for _ in range(n_proteins)]
    foreign_prots = [protein(500) for _ in range(n_proteins)]
    fasta = work / "synthetic.fasta"
    with open(fasta, "w") as fh:
        for i, s in enumerate(prots, 1):
            fh.write(f">sp|SYN{i:03d}|SYN{i:03d}_HUMAN Synthetic protein {i} OS=Homo sapiens OX=9606 GN=SYN{i} PE=1 SV=1\n")
            for k in range(0, len(s), 60):
                fh.write(s[k:k + 60] + "\n")

    def pick(proteins, n, exclude=frozenset()):
        peps = sorted({q for s in proteins for q in _tryptic(s) if 7 <= len(q) <= 20 and q not in exclude})
        idx = rng.choice(len(peps), size=min(n, len(peps)), replace=False)
        return [peps[i] for i in sorted(idx)]

    target_peps = set(q for s in prots for q in _tryptic(s))
    planted = pick(prots, n_planted)
    foreign = pick(foreign_prots, n_foreign, exclude=frozenset(target_peps))
    res = _residue_masses()

    peps = []
    for kind, seqs in (("planted", planted), ("foreign", foreign)):
        for s in seqs:
            M = sum(res[a] for a in s) + H2O
            z = 2 if M < 1800 else 3
            h = np.mean([HYDRO[a] for a in s])
            peps.append(dict(seq=s, kind=kind, M=M, z=z, mz=(M + z * PROTON) / z, h=h,
                             abund=float(10 ** rng.normal(6.5, 0.4))))
    hs = np.array([q["h"] for q in peps])
    rt_lo, rt_hi = 120.0, 2520.0  # seconds
    for q in peps:  # retention time: monotone in hydrophobicity, plus noise
        q["rt"] = float(np.clip(rt_lo + (q["h"] - hs.min()) / (np.ptp(hs) or 1) * (rt_hi - rt_lo)
                                + rng.normal(0, 20), rt_lo, rt_hi))

    def envelope(M, z, k=5):
        lam = M * 0.000594
        w = np.exp(-lam) * np.cumprod(np.r_[1.0, lam / np.arange(1, k)])
        return (M + z * PROTON) / z + np.arange(k) * C13 / z, w / w.max()

    sigma = 6.0
    ms1_rts = np.arange(0.0, rt_hi + 150.0, 2.0)
    exp = oms.MSExperiment()
    scan = 0
    by_ms1: dict[int, list[dict]] = {}
    for q in peps:  # MS2 right after the MS1 closest to the apex
        by_ms1.setdefault(int(np.argmin(np.abs(ms1_rts - q["rt"]))), []).append(q)
    centroid = getattr(getattr(oms.SpectrumSettings, "SpectrumType", None), "CENTROID", None)
    by_rt = sorted(peps, key=lambda q: q["rt"])
    rts_sorted = np.array([q["rt"] for q in by_rt])
    for q in peps:
        q["env"] = envelope(q["M"], q["z"])
    for i, t in enumerate(ms1_rts):
        scan += 1
        mzs, ints = [rng.uniform(350, 1500, 60)], [10 ** rng.normal(4.3, 0.3, 60)]
        lo, hi = np.searchsorted(rts_sorted, [t - 4 * sigma, t + 4 * sigma])
        for q in by_rt[lo:hi]:
            f = np.exp(-0.5 * ((t - q["rt"]) / sigma) ** 2)
            m, w = q["env"]
            mzs.append(m + rng.normal(0, m * 1e-6))
            ints.append(q["abund"] * f * w)
        mz = np.concatenate(mzs)
        it = np.concatenate(ints)
        o = np.argsort(mz)
        sp = oms.MSSpectrum()
        sp.setMSLevel(1)
        sp.setRT(float(t))
        sp.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
        if centroid is not None:
            sp.setType(centroid)
        sp.set_peaks((mz[o].astype(np.float64), it[o].astype(np.float32)))
        exp.addSpectrum(sp)
        for j, q in enumerate(by_ms1.get(i, [])):
            scan += 1
            q["scan"] = scan
            s = q["seq"]
            cum = np.cumsum([res[a] for a in s])
            b = cum[:-1] + PROTON
            y = (cum[-1] - cum[:-1]) + H2O + PROTON
            n = len(s)
            pos = np.arange(1, n)
            yi = 0.3 + np.exp(-0.5 * ((pos - n / 2) / (n / 4)) ** 2) * rng.uniform(0.5, 1.0, n - 1)
            bi = 0.4 * yi[::-1] * rng.uniform(0.3, 1.0, n - 1)
            if "P" in s[1:]:  # enhanced cleavage N-terminal to proline
                for k in [k for k in range(1, n) if s[k] == "P"]:
                    yi[k - 1] *= 3.0
            keep_b = rng.random(n - 1) < 0.8
            keep_y = rng.random(n - 1) < 0.9
            nm = rng.uniform(150, min(2000, q["M"]), 30)
            mz = np.concatenate([b[keep_b], y[keep_y], nm])
            it = np.concatenate([bi[keep_b], yi[::-1][keep_y], rng.uniform(0.01, 0.15, 30)]) * 1e5
            o = np.argsort(mz)
            sp = oms.MSSpectrum()
            sp.setMSLevel(2)
            sp.setRT(float(t + 0.1 + 0.15 * j))
            sp.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
            if centroid is not None:
                sp.setType(centroid)
            pc = oms.Precursor()
            pc.setMZ(float(q["mz"]))
            pc.setCharge(int(q["z"]))
            pc.setIntensity(float(q["abund"]))
            pc.setIsolationWindowLowerOffset(0.8)
            pc.setIsolationWindowUpperOffset(0.8)
            sp.setPrecursors([pc])
            sp.set_peaks((mz[o].astype(np.float64), it[o].astype(np.float32)))
            exp.addSpectrum(sp)
    exp.updateRanges()
    mzml = work / "synthetic.mzML"
    oms.MzMLFile().store(str(mzml), exp)
    return fasta, mzml, peps


# --------------------------------------------------------------------------------------
# 4b. smoke test through proteoscan.pipeline.run
# --------------------------------------------------------------------------------------

def _timed(stage: str, fn, timings: dict):
    def wrapper(*a, **k):
        t = time.time()
        try:
            return fn(*a, **k)
        finally:
            timings[stage] = timings.get(stage, 0.0) + time.time() - t
    return wrapper


def smoke_test(work: Path, threads: int, processes: int, two_pass: bool = False) -> bool:
    print(f"\n== Smoke test (work dir {work})")
    import numpy as np
    import pandas as pd

    ms2pip_dir = Path.home() / ".ms2pip"
    if not (ms2pip_dir / "model_20210416_HCD2021_Y.xgboost").exists():
        report(WARN, "MS2PIP model", f"HCD2021 not in {ms2pip_dir}; this run downloads ~955 MB first")

    t0 = time.time()
    try:
        fasta, mzml, peps = make_synthetic(work)
    except Exception as e:
        report(FAIL, "synthetic data", f"{type(e).__name__}: {e}")
        traceback.print_exc()
        return False
    n_ms2 = len(peps)
    report(OK, "synthetic data", f"{mzml.name} ({mzml.stat().st_size / 1e6:.1f} MB, {n_ms2} MS2), "
           f"{fasta.name}  [{time.time() - t0:.1f}s]")

    from proteoscan import pipeline, rescore_features
    from proteoscan.config import SearchConfig

    cfg_path = REPO / "configs" / "peaks_facility_human.yml"
    cfg = SearchConfig.from_yaml(cfg_path) if cfg_path.exists() else SearchConfig()
    cfg.database = str(fasta)
    cfg.threads = threads

    timings: dict[str, float] = {}
    saved = {}
    patches = {
        (pipeline, "stage_precursors"): "MS1 features + precursor list",
        (pipeline, "run_comet"): "Comet search",
        (pipeline, "read_comet"): "read_comet",
        (rescore_features, "ms2pip_correlations"): "MS2PIP features",
        (rescore_features, "deeplc_predictions"): "DeepLC features",
        (pipeline, "rescore_tda"): "mokapot rescoring",
        (pipeline, "stage_two_pass"): "two-pass search",
        (pipeline, "write_export"): "PEAKS CSV export",
    }
    patches = {k: v for k, v in patches.items() if hasattr(*k)}
    for (mod, name), stage in patches.items():
        saved[(mod, name)] = getattr(mod, name)
        setattr(mod, name, _timed(stage, getattr(mod, name), timings))
    out = work / "run"
    t = time.time()
    try:
        kw = {}
        if "two_pass" in inspect.signature(pipeline.run).parameters:
            # The second pass re-scores only candidates new in pass 2; on a data set this small
            # there are fewer than the 50 confident PSMs DeepLC needs to calibrate them, so the
            # pass is skipped unless --two-pass is given (full-size runs are unaffected).
            kw["two_pass"] = two_pass
            if not two_pass:
                report(OK, "two-pass stage", "skipped on the synthetic set (use --two-pass to include it)")
        stats = pipeline.run(cfg, [mzml], out, mode="tda", processes=processes, sample_name="Sample 1",
                             candidates_per_spectrum=3, **kw)
    except Exception as e:
        report(FAIL, "pipeline.run", f"{type(e).__name__}: {e}")
        traceback.print_exc()
        return False
    finally:
        for (mod, name), fn in saved.items():
            setattr(mod, name, fn)
    total = time.time() - t
    for stage in patches.values():
        if stage in timings:
            report(OK, f"stage: {stage}", f"{timings[stage]:.1f}s")
    report(OK, "pipeline.run total", f"{total:.1f}s, peak RSS {peak_rss_gb(0):.2f} GB (python), "
           f"{peak_rss_gb(-1):.2f} GB (largest child, e.g. comet)")

    ok = True
    feat = pd.read_parquet(out / "candidates_features.parquet")
    for col in ("ms2_corr", "pred_rt", "hyperscore"):
        v = feat[col].to_numpy(dtype=float) if col in feat else np.array([])
        good = len(v) and np.isfinite(v).all() and np.nanstd(v) > 0
        report(OK if good else FAIL, f"feature column {col}",
               f"n={len(v)} mean={np.nanmean(v):.3f} sd={np.nanstd(v):.3f}" if len(v) else "missing")
        ok &= bool(good)

    export = Path(stats["export_dir"])
    files = sorted(p.name for p in export.glob("*.csv"))
    need = {"db.psms.csv", "db.peptides.csv", "db.proteins.csv", "db.protein-peptides.csv"}
    miss = need - set(files)
    report(FAIL if miss else OK, "export files", f"missing {sorted(miss)}" if miss else ", ".join(files))
    ok &= not miss
    if not miss:
        psms = pd.read_csv(export / "db.psms.csv")
        found = set(psms["Peptide"].astype(str).str.replace(r"\(.*?\)|\[.*?\]", "", regex=True)
                    .str.replace(r"^[A-Z-]\.|\.[A-Z-]$", "", regex=True))
        planted = {q["seq"] for q in peps if q["kind"] == "planted"}
        foreign = {q["seq"] for q in peps if q["kind"] == "foreign"}
        rec = len(planted & found) / len(planted)
        fp = len(foreign & found)
        st = OK if rec >= 0.8 else FAIL
        report(st, "planted peptides identified", f"{len(planted & found)}/{len(planted)} ({rec:.0%}); "
               f"PSM rows {len(psms)}; foreign peptides reported {fp} (must be 0)")
        ok &= rec >= 0.8 and fp == 0
    report(OK if ok else FAIL, "smoke test", json.dumps({k: stats[k] for k in
           ("n_target_psms", "n_decoy_psms", "n_peptides", "n_protein_groups") if k in stats}))
    return ok


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--work", help="work directory for the smoke test (default: a new temporary directory)")
    ap.add_argument("--keep", action="store_true", help="keep the temporary work directory")
    ap.add_argument("--skip-smoke", action="store_true", help="only check tools and imports")
    ap.add_argument("--threads", type=int, default=min(4, os.cpu_count() or 1), help="Comet threads")
    ap.add_argument("--processes", type=int, default=min(4, os.cpu_count() or 1), help="MS2PIP processes")
    ap.add_argument("--two-pass", action="store_true",
                    help="also run the pipeline's second (protein-restricted) search pass")
    ap.add_argument("--start-method", choices=["spawn", "fork", "forkserver"], default="spawn",
                    help="multiprocessing start method (spawn = the macOS default; used on every platform "
                         "so that a Linux run exercises the macOS code path)")
    a = ap.parse_args(argv)
    mp.set_start_method(a.start_method, force=True)
    import warnings
    # DeepLC 4.5 ignores mass-only ProForma modifications (C[+57.0215] etc.) and warns per batch
    warnings.filterwarnings("ignore", message="Skipping modification without known composition")
    import logging
    logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s %(message)s")

    t0 = time.time()
    check_platform()
    ok_imports = check_imports()
    ok_tools = check_tools()
    ok_smoke = True
    if not a.skip_smoke:
        if not ok_imports or not shutil.which("comet"):
            report(FAIL, "smoke test", "skipped: fix the failures above first")
            ok_smoke = False
        else:
            work = Path(a.work) if a.work else Path(tempfile.mkdtemp(prefix="proteoscan_check_"))
            work.mkdir(parents=True, exist_ok=True)
            try:
                ok_smoke = smoke_test(work, a.threads, a.processes, a.two_pass)
            finally:
                if not a.work and not a.keep:
                    shutil.rmtree(work, ignore_errors=True)
                else:
                    print(f"  work directory kept: {work}")
    n_fail = sum(s == FAIL for s, _, _ in _results)
    n_warn = sum(s == WARN for s, _, _ in _results)
    verdict = "PASS" if ok_imports and ok_tools and ok_smoke and n_fail == 0 else "FAIL"
    print(f"\n== {verdict}: {n_fail} failure(s), {n_warn} warning(s), {time.time() - t0:.0f}s total")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
