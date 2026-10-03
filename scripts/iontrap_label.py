"""Label ProteomeTools ion-trap MS2 spectra by targeted matching against the known pool content.

    python scripts/iontrap_label.py --pool TUM_first_pool_30 --split-dir <run>/split --out labels.parquet

For each pool peptide (Supplementary Table 1 of Zolg et al. 2017, plus the spiked QC peptides; C as
carbamidomethyl, optionally each single M oxidised) and charge 1-4, candidate scans are the MS2 scans
of each ion-trap scan type (hcd_it, cid_it) whose precursor m/z matches within 15 ppm (isotope offset
0 or +1, charge equal to the scan's charge when the scan reports one). Each (candidate, scan) is scored
with a binomial peak-match score over b/y ions (1+, and 2+ for precursors >= 3+) against the top 10
peaks per 100 m/z, as in Andromeda/MS Amanda. Every target is paired with its reversed decoy (same
precursor mass, C-terminal residue kept) scored on the same scans; the score threshold per scan type
is set where the best-per-(peptide, charge) decoys reach 1% of the targets. Kept: the best scan per
(modified peptide, charge, scan type) above that threshold with at least 6 matched b/y ions, a scan
being assigned to one peptide only.

Output columns: pool, run, scan_type, scan, modseq, sequence, charge, iso, score, n_by, n_theo,
precursor_mz, rt, nce, activation, scan_mode, lo_mz, hi_mz, mzml.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # single-threaded numpy/scipy: many runs go in parallel

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from proteoscan.fragpred import MONO, H2O, PROTON, fragment_mz  # noqa: E402

TRUTH = Path("/drive_a/drew/public-data/PXD004732/truth")
C13 = 1.0033548
FILTER = re.compile(r"^(FTMS|ITMS) \+ \w+ \w+ ?(?P<flags>(?:[a-z] )*)Full ms2 [0-9.]+@(?P<act>[a-z]+)(?P<nce>[0-9.]+)")
MODE = {"r": "rapid", "t": "turbo", "n": "normal", "z": "zoom", "e": "enhanced"}
OX, CAM = 15.994915, 57.021464


def parse_filter(f: str):
    m = FILTER.search(f or "")
    if not m:
        return None, np.nan, ""
    flags = m.group("flags").split()
    mode = next((MODE[x] for x in flags if x in MODE), "normal")
    return m.group("act"), float(m.group("nce")), mode


def ms2_table(mzml: Path):
    import pyopenms as oms
    od = oms.OnDiscMSExperiment()
    od.openFile(str(mzml))
    meta = od.getMetaData()
    rows = []
    for i in range(meta.getNrSpectra()):
        s = meta[i]
        if s.getMSLevel() != 2:
            continue
        mm = re.search(r"scan=(\d+)", s.getNativeID())
        p = s.getPrecursors()[0]
        fs = str(s.getMetaValue("filter string")) if s.metaValueExists("filter string") else ""
        act, nce, mode = parse_filter(fs)
        w = s.getInstrumentSettings().getScanWindows()
        lo, hi = (w[0].begin, w[0].end) if len(w) else (0.0, 1e9)
        an = "ft" if fs.startswith("FTMS") else "it"
        rows.append((i, int(mm.group(1)), p.getMZ(), p.getCharge(), s.getRT() / 60.0, act, nce, mode, lo, hi, an))
    return od, pd.DataFrame(rows, columns=["idx", "scan", "precursor_mz", "scan_charge", "rt", "activation", "nce",
                                           "scan_mode", "lo_mz", "hi_mz", "analyzer"])


def top_peaks(mz, it, q=10, window=100.0):
    if len(mz) == 0:
        return mz, it
    b = np.floor(mz / window).astype(int)
    keep = np.zeros(len(mz), bool)
    order = np.lexsort((-it, b))
    bs = b[order]
    start = np.r_[0, np.nonzero(np.diff(bs))[0] + 1]
    rank = np.arange(len(order)) - np.repeat(start, np.diff(np.r_[start, len(order)]))
    keep[order[rank < q]] = True
    return mz[keep], it[keep]


def forms_for(seq: str, max_ox: int = 1):
    """(modseq, residue masses) for C carbamidomethyl, plus single M oxidations."""
    base = np.array([MONO[a] + (CAM if a == "C" else 0.0) for a in seq])
    def render(ox_pos):
        return "".join(a + ("[+57.0215]" if a == "C" else "") + ("[+15.9949]" if i in ox_pos else "") for i, a in enumerate(seq))
    out = [(render(()), base)]
    if max_ox:
        for i, a in enumerate(seq):
            if a == "M":
                r = base.copy(); r[i] += OX
                out.append((render((i,)), r))
    return out


def score_match(res: np.ndarray, charge: int, mz: np.ndarray, lo: float, hi: float, tol: float):
    types = ["b_z1", "y_z1"] + (["b_z2", "y_z2"] if charge >= 3 else [])
    theo = fragment_mz(res, types).ravel()
    theo = theo[(theo >= lo) & (theo <= hi)]
    n = len(theo)
    if n == 0 or len(mz) == 0:
        return 0.0, 0, n
    l = np.searchsorted(mz, theo - tol, "left"); h = np.searchsorted(mz, theo + tol, "right")
    hit = h > l
    k = int(hit.sum())
    by1 = int(hit[: 2 * (len(res) - 1)].sum()) if len(hit) else 0
    inwin = mz[(mz >= lo) & (mz <= hi)]
    p = min(0.5, len(inwin) * 2 * tol / max(hi - lo, 1.0))
    sc = -np.log10(max(binom.sf(k - 1, n, p), 1e-300)) if k > 0 else 0.0
    return float(sc), by1, n


def label_run(mzml: Path, pool: str, peptides: list[str], ppm: float = 15.0, tol: float = 0.5, fdr: float = 0.01,
              min_by: int = 6, scan_type: str = "", max_ox: int = 1):
    od, ms2 = ms2_table(mzml)
    ms2 = ms2.sort_values("precursor_mz").reset_index(drop=True)
    pmz = ms2["precursor_mz"].to_numpy()
    zc = ms2["scan_charge"].to_numpy(); sidx = ms2["idx"].to_numpy(); lo_a = ms2["lo_mz"].to_numpy(); hi_a = ms2["hi_mz"].to_numpy()
    cache: dict[int, tuple] = {}
    rows = []
    for decoy in (False, True):
        for seq in peptides:
            s = seq[:-1][::-1] + seq[-1] if decoy else seq
            if len(s) < 5 or any(a not in MONO for a in s):
                continue
            for modseq, res in forms_for(s, max_ox):
                mass = res.sum() + H2O
                for z in (1, 2, 3, 4):
                    for iso in (0, 1):
                        target = (mass + iso * C13 + z * PROTON) / z
                        d = target * ppm * 1e-6
                        a, b = np.searchsorted(pmz, target - d), np.searchsorted(pmz, target + d, "right")
                        for j in range(a, b):
                            if zc[j] and int(zc[j]) != z:
                                continue
                            if j not in cache:
                                mz, it = od.getSpectrum(int(sidx[j])).get_peaks()
                                mz = np.asarray(mz, float); it = np.asarray(it, float)
                                o = np.argsort(mz); cache[j] = top_peaks(mz[o], it[o])
                            sc, nby, nth = score_match(res, z, cache[j][0], lo_a[j], hi_a[j], tol)
                            rows.append((decoy, seq, modseq, s, z, iso, j, sc, nby, nth))
    if not rows:
        return pd.DataFrame(), {}
    m = pd.DataFrame(rows, columns=["decoy", "target_seq", "modseq", "sequence", "charge", "iso", "j", "score", "n_by", "n_theo"])
    m = m.sort_values("score", ascending=False)
    best = m.drop_duplicates(["decoy", "modseq", "charge"])
    best = best[best["n_by"] >= min_by]
    t = best.sort_values("score", ascending=False)
    nd = np.cumsum(t["decoy"].to_numpy()); nt = np.cumsum(~t["decoy"].to_numpy())
    q = nd / np.maximum(nt, 1)
    ok = np.nonzero(q <= fdr)[0]
    thr = float(t["score"].to_numpy()[ok[-1]]) if len(ok) else np.inf
    acc = best[(~best["decoy"]) & (best["score"] >= thr)].sort_values("score", ascending=False)
    acc = acc.drop_duplicates("j")  # a scan carries one peptide
    stats = dict(scan_type=scan_type, n_ms2=len(ms2), threshold=thr, n_target_best=int((~best["decoy"]).sum()),
                 n_decoy_best=int(best["decoy"].sum()), n_decoy_above=int((best["decoy"] & (best["score"] >= thr)).sum()),
                 n_labelled=len(acc), n_peptides=acc["target_seq"].nunique())
    out = acc.merge(ms2, left_on="j", right_index=True).drop(columns=["j", "decoy", "target_seq", "idx"])
    out["pool"] = pool; out["scan_type"] = scan_type; out["mzml"] = str(mzml); out["run"] = mzml.name.split(".")[0]
    return out, stats


def label_events(mzml: Path, pool: str, peptides: list[str], ppm: float = 15.0, tol_it: float = 0.5,
                 tol_ft: float = 0.02, fdr: float = 0.01, min_by: int = 6, max_ox: int = 1,
                 label_types=("hcd_it", "cid_it")):
    """Event mode (default): the instrument acquires every precursor as a block of consecutive MS2
    scans (CID-IT, HCD-IT, HCD-OT x2). A candidate is scored on the whole block (sum of the per-scan
    binomial scores; Orbitrap scans at 0.02 Da, ion-trap scans at 0.5 Da), the decoy-calibrated 1%
    threshold is applied to the best block per (modified peptide, charge), and the accepted block's
    ion-trap scans are labelled when they themselves match >= min_by b/y ions. The scan's own charge
    annotation is not required to agree (some runs annotate nearly every precursor as 2+)."""
    od, ms2 = ms2_table(mzml)
    ms2 = ms2.sort_values("scan").reset_index(drop=True)
    ms2["scan_type"] = ms2["activation"].astype(str) + "_" + ms2["analyzer"]
    sc_ = ms2["scan"].to_numpy(); pm_ = ms2["precursor_mz"].to_numpy()
    new_ev = np.r_[True, (np.diff(sc_) != 1) | (np.abs(np.diff(pm_)) > 0.01)]
    ms2["event"] = np.cumsum(new_ev) - 1
    ev = ms2.groupby("event").agg(precursor_mz=("precursor_mz", "first"), first=("scan", "first")).reset_index()
    ev = ev.sort_values("precursor_mz").reset_index(drop=True)
    members = ms2.groupby("event").indices
    epmz = ev["precursor_mz"].to_numpy(); eid = ev["event"].to_numpy()
    sidx = ms2["idx"].to_numpy(); lo_a = ms2["lo_mz"].to_numpy(); hi_a = ms2["hi_mz"].to_numpy()
    tol_a = np.where(ms2["analyzer"].to_numpy() == "ft", tol_ft, tol_it); types_a = ms2["scan_type"].to_numpy()
    cache: dict[int, tuple] = {}
    rows = []
    for decoy in (False, True):
        for seq in peptides:
            s = seq[:-1][::-1] + seq[-1] if decoy else seq
            if len(s) < 5 or any(a not in MONO for a in s):
                continue
            for modseq, res in forms_for(s, max_ox):
                mass = res.sum() + H2O
                for z in (1, 2, 3, 4):
                    for iso in (0, 1):
                        target = (mass + iso * C13 + z * PROTON) / z
                        d = target * ppm * 1e-6
                        a, b = np.searchsorted(epmz, target - d), np.searchsorted(epmz, target + d, "right")
                        for k in range(a, b):
                            total = 0.0; per = {}
                            for j in members[eid[k]]:
                                if j not in cache:
                                    mz, it = od.getSpectrum(int(sidx[j])).get_peaks()
                                    mz = np.asarray(mz, float); it = np.asarray(it, float)
                                    o = np.argsort(mz); cache[j] = top_peaks(mz[o], it[o])
                                sc, nby, nth = score_match(res, z, cache[j][0], lo_a[j], hi_a[j], tol_a[j])
                                total += sc
                                if types_a[j] in label_types:
                                    per[types_a[j]] = (j, sc, nby, nth)
                            rows.append((decoy, seq, modseq, s, z, iso, int(eid[k]), total, per))
    if not rows:
        return pd.DataFrame(), {}
    m = pd.DataFrame(rows, columns=["decoy", "target_seq", "modseq", "sequence", "charge", "iso", "event", "event_score", "per"])
    best = m.sort_values("event_score", ascending=False).drop_duplicates(["decoy", "modseq", "charge"])
    nd = np.cumsum(best["decoy"].to_numpy()); nt = np.cumsum(~best["decoy"].to_numpy())
    ok = np.nonzero(nd / np.maximum(nt, 1) <= fdr)[0]
    thr = float(best["event_score"].to_numpy()[ok[-1]]) if len(ok) else np.inf
    acc = best[(~best["decoy"]) & (best["event_score"] >= thr)].drop_duplicates("event")  # one peptide per block
    out = []
    for r in acc.itertuples(index=False):
        for t, (j, sc, nby, nth) in r.per.items():
            if nby >= min_by:
                out.append(dict(modseq=r.modseq, sequence=r.sequence, charge=r.charge, iso=r.iso, event=r.event,
                                event_score=r.event_score, score=sc, n_by=nby, n_theo=nth, j=j))
    out = pd.DataFrame(out)
    stats = dict(n_ms2=len(ms2), n_events=len(ev), threshold=thr, n_target_best=int((~best["decoy"]).sum()),
                 n_decoy_above=int((best["decoy"] & (best["event_score"] >= thr)).sum()), n_events_accepted=len(acc),
                 n_peptides=int(acc["target_seq"].nunique()))
    if not len(out):
        return out, stats
    out = out.merge(ms2.drop(columns=["event"]), left_on="j", right_index=True).drop(columns=["j", "idx"])
    for t, g in out.groupby("scan_type"):
        stats[f"n_labelled_{t}"] = len(g)
    out["pool"] = pool; out["mzml"] = str(mzml); out["run"] = mzml.name.split(".")[0]
    return out, stats


def pool_peptides(pool: str) -> list[str]:
    t = pd.read_parquet(TRUTH / "pool_peptides.parquet")
    s = list(dict.fromkeys(t.loc[t["pool"] == pool, "sequence"]))
    s += [x for x in pd.read_parquet(TRUTH / "qc_peptides.parquet")["sequence"] if x not in set(s)]
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--split-dir", help="per-type mode: directory of split mzMLs")
    ap.add_argument("--mzml", help="event mode: the unsplit mzML")
    ap.add_argument("--out", required=True)
    ap.add_argument("--types", nargs="+", default=["hcd_it", "cid_it"])
    a = ap.parse_args()
    bench = set(Path("/drive_a/drew/public-data/PXD004732/pools.txt").read_text().split())
    peps = pool_peptides(a.pool)
    if a.mzml:
        df, st = label_events(Path(a.mzml), a.pool, peps, label_types=tuple(a.types))
        st.update(pool=a.pool, benchmark_pool=a.pool in bench, mode="events")
        print(st, flush=True)
        df.to_parquet(a.out)
        pd.DataFrame([st]).to_csv(Path(a.out).with_suffix(".stats.csv"), index=False)
        return
    parts, stats = [], []
    for t in a.types:
        mz = sorted(Path(a.split_dir).glob(f"*.{t}.mzML"))
        if not mz:
            continue
        df, st = label_run(mz[0], a.pool, peps, scan_type=t)
        st.update(pool=a.pool, benchmark_pool=a.pool in bench)
        stats.append(st); parts.append(df)
        print(st, flush=True)
    out = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    out.to_parquet(a.out)
    pd.DataFrame(stats).to_csv(Path(a.out).with_suffix(".stats.csv"), index=False)


if __name__ == "__main__":
    main()
