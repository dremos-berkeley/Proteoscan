"""Trace-based MS1 feature detection (``algorithm: traces`` in peaksparity.features).

A small, dependency-light detector written for centroided high-resolution MS1 (Orbitrap),
built to imitate the behaviour visible in PEAKS feature exports (short features that end
where the signal disappears or at a valley; charge 1-9; monoisotopic m/z).  Steps:

1. link centroids of consecutive MS1 scans into mass traces (nearest m/z within
   ``trace_ppm``, gaps of up to ``trace_gap`` scans);
2. split each trace into elution peaks at valleys deeper than ``valley_ratio`` x the
   smaller neighbouring maximum (3-point smoothed profile), keeping peaks with at least
   ``min_points`` centroids;
3. for every elution peak as monoisotopic candidate and every charge, chain isotope
   *traces* at mono + k * 1.00286/z (``iso_ppm``): within the candidate's scan range the
   trace must have >= ``iso_min_points`` centroids and correlate with it (Pearson >=
   ``min_corr``); isotopes are compared on whole traces so that differently split
   isotope traces still group.  The chain's summed intensities are scored against an
   averagine (Poisson) envelope (cosine >= ``min_fit``);
4. accept hypotheses greedily by score (explained intensity x fit), each elution peak used
   by at most one feature.

Everything runs in numba on the flat MS1 arrays of :func:`peaksparity.features._ms1_arrays`.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np
from numba import njit

ISO_SPACING = 1.00286864  # average isotope spacing of peptides (averagine)
C13 = 1.0033548378
PROTON = 1.00727646688

DEFAULTS: dict[str, Any] = {
    "trace_ppm": 6.0,
    "trace_gap": 1,
    "valley_ratio": 0.5,
    "smooth": 1,
    "min_points": 2,
    "iso_ppm": 8.0,
    "min_corr": 0.5,
    "min_fit": 0.6,
    "max_isotopes": 6,
    "charge_low": 1,
    "charge_high": 8,
    "min_isotopes": 2,
    "min_intensity": 0.0,
    "score_power": 1.0,
    "iso_min_points": 2,
}


# --------------------------------------------------------------------------------------
# 1. mass traces
# --------------------------------------------------------------------------------------

@njit(cache=True)
def _link(mz, inten, offs, ppm, max_gap):
    n = len(mz)
    nscan = len(offs) - 1
    tid = np.full(n, -1, np.int64)
    cap = 1 << 16
    act_mz = np.empty(cap)
    act_w = np.empty(cap)
    act_id = np.empty(cap, np.int64)
    act_last = np.empty(cap, np.int64)
    n_act = 0
    next_id = 0
    for s in range(nscan):
        a, b = offs[s], offs[s + 1]
        npk = b - a
        # candidate pairs (distance, peak, active)
        cd = np.empty(2 * npk)
        cp = np.empty(2 * npk, np.int64)
        ca = np.empty(2 * npk, np.int64)
        nc = 0
        if n_act > 0:
            for p in range(a, b):
                j = np.searchsorted(act_mz[:n_act], mz[p])
                for q in (j - 1, j):
                    if 0 <= q < n_act:
                        d = abs(act_mz[q] - mz[p]) / mz[p] * 1e6
                        if d <= ppm:
                            cd[nc] = d
                            cp[nc] = p
                            ca[nc] = q
                            nc += 1
        order = np.argsort(cd[:nc])
        peak_used = np.zeros(npk, np.bool_)
        act_used = np.zeros(n_act, np.bool_)
        for o in order:
            p = cp[o]
            q = ca[o]
            if peak_used[p - a] or act_used[q]:
                continue
            peak_used[p - a] = True
            act_used[q] = True
            tid[p] = act_id[q]
            w = act_w[q] + inten[p]
            act_mz[q] = (act_mz[q] * act_w[q] + mz[p] * inten[p]) / w
            act_w[q] = w
            act_last[q] = s
        # rebuild active set: survivors + new traces
        m = 0
        new_mz = np.empty(n_act + npk)
        new_w = np.empty(n_act + npk)
        new_id = np.empty(n_act + npk, np.int64)
        new_last = np.empty(n_act + npk, np.int64)
        for q in range(n_act):
            if s - act_last[q] <= max_gap:
                new_mz[m] = act_mz[q]
                new_w[m] = act_w[q]
                new_id[m] = act_id[q]
                new_last[m] = act_last[q]
                m += 1
        for p in range(a, b):
            if not peak_used[p - a]:
                tid[p] = next_id
                new_mz[m] = mz[p]
                new_w[m] = inten[p]
                new_id[m] = next_id
                new_last[m] = s
                next_id += 1
                m += 1
        o2 = np.argsort(new_mz[:m])
        if m > cap:
            cap = m * 2
            act_mz = np.empty(cap)
            act_w = np.empty(cap)
            act_id = np.empty(cap, np.int64)
            act_last = np.empty(cap, np.int64)
        for i in range(m):
            act_mz[i] = new_mz[o2[i]]
            act_w[i] = new_w[o2[i]]
            act_id[i] = new_id[o2[i]]
            act_last[i] = new_last[o2[i]]
        n_act = m
    return tid, next_id


# --------------------------------------------------------------------------------------
# 2. elution peaks
# --------------------------------------------------------------------------------------

@njit(cache=True)
def _segments(order, tid, scan, mz, inten, n_traces, valley_ratio, min_points, smooth):
    """Split traces into elution peaks.  Returns segment table + dense profiles."""
    n = len(order)
    # trace boundaries in `order`
    starts = np.empty(n_traces + 1, np.int64)
    starts[:] = n
    for i in range(n - 1, -1, -1):
        starts[tid[order[i]]] = i
    starts[n_traces] = n
    for t in range(n_traces - 1, -1, -1):
        if starts[t] > starts[t + 1]:
            starts[t] = starts[t + 1]
    cap = n // 2 + 16
    s_first = np.empty(cap, np.int64)
    s_last = np.empty(cap, np.int64)
    s_apex = np.empty(cap, np.int64)
    s_mz = np.empty(cap)
    s_apexint = np.empty(cap)
    s_sum = np.empty(cap)
    s_npts = np.empty(cap, np.int64)
    s_poff = np.empty(cap + 1, np.int64)
    s_trace = np.empty(cap, np.int64)
    prof = np.empty(n * 2 + 16)
    # compact trace table (traces with >= min_points centroids)
    tcap = n_traces
    t_first = np.empty(tcap, np.int64)
    t_last = np.empty(tcap, np.int64)
    t_mz = np.empty(tcap)
    t_poff = np.empty(tcap + 1, np.int64)
    t_seg0 = np.empty(tcap, np.int64)
    t_seg1 = np.empty(tcap, np.int64)
    t_prof = np.empty(n * 2 + 16)
    nt = 0
    ntp = 0
    t_poff[0] = 0
    ns = 0
    npf = 0
    s_poff[0] = 0
    for t in range(n_traces):
        a, b = starts[t], starts[t + 1]
        if b - a < min_points:
            continue
        sc0 = scan[order[a]]
        sc1 = scan[order[b - 1]]
        L = sc1 - sc0 + 1
        raw = np.zeros(L)
        rmz = np.zeros(L)
        tw = 0.0
        twm = 0.0
        for i in range(a, b):
            k = order[i]
            raw[scan[k] - sc0] = inten[k]
            rmz[scan[k] - sc0] = mz[k]
            tw += inten[k]
            twm += inten[k] * mz[k]
        if ntp + L > len(t_prof):
            t_prof = np.concatenate((t_prof, np.empty(len(t_prof) + L)))
        for i in range(L):
            t_prof[ntp + i] = raw[i]
        ntp += L
        t_first[nt] = sc0
        t_last[nt] = sc1
        t_mz[nt] = twm / tw
        t_poff[nt + 1] = ntp
        t_seg0[nt] = ns
        sm = np.empty(L)
        for i in range(L):
            lo = max(0, i - smooth)
            hi = min(L - 1, i + smooth)
            acc = 0.0
            wt = 0.0
            for j in range(lo, hi + 1):
                w = smooth + 1.0 - abs(j - i)  # triangular kernel
                acc += raw[j] * w
                wt += w
            sm[i] = acc / wt
        # cut points at deep valleys between local maxima
        cuts = np.empty(L + 1, np.int64)
        nc = 0
        cuts[nc] = 0
        nc += 1
        prev_max = -1
        left_peak = 0.0
        for i in range(L):
            up = i == 0 or sm[i] > sm[i - 1]
            down = i == L - 1 or sm[i] >= sm[i + 1]
            if not (up and down):
                continue
            if prev_max >= 0:
                vi = prev_max
                for j in range(prev_max, i + 1):
                    if sm[j] < sm[vi]:
                        vi = j
                if sm[vi] < valley_ratio * min(left_peak, sm[i]):
                    cuts[nc] = vi
                    nc += 1
                    left_peak = sm[i]
                elif sm[i] > left_peak:
                    left_peak = sm[i]
            else:
                left_peak = sm[i]
            prev_max = i
        cuts[nc] = L
        nc += 1
        for c in range(nc - 1):
            lo = cuts[c]
            hi = cuts[c + 1]  # exclusive; the valley scan belongs to the right part
            # trim zeros at both ends
            while lo < hi and raw[lo] == 0.0:
                lo += 1
            while hi > lo and raw[hi - 1] == 0.0:
                hi -= 1
            npts = 0
            for i in range(lo, hi):
                if raw[i] > 0:
                    npts += 1
            if npts < min_points:
                continue
            if ns >= cap:
                cap2 = cap * 2
                s_first = np.concatenate((s_first, np.empty(cap2 - cap, np.int64)))
                s_last = np.concatenate((s_last, np.empty(cap2 - cap, np.int64)))
                s_apex = np.concatenate((s_apex, np.empty(cap2 - cap, np.int64)))
                s_mz = np.concatenate((s_mz, np.empty(cap2 - cap)))
                s_apexint = np.concatenate((s_apexint, np.empty(cap2 - cap)))
                s_sum = np.concatenate((s_sum, np.empty(cap2 - cap)))
                s_npts = np.concatenate((s_npts, np.empty(cap2 - cap, np.int64)))
                s_poff = np.concatenate((s_poff, np.empty(cap2 - cap, np.int64)))
                s_trace = np.concatenate((s_trace, np.empty(cap2 - cap, np.int64)))
                cap = cap2
            if npf + (hi - lo) > len(prof):
                prof = np.concatenate((prof, np.empty(len(prof))))
            best = lo
            wsum = 0.0
            wmz = 0.0
            tot = 0.0
            for i in range(lo, hi):
                prof[npf + i - lo] = raw[i]
                if sm[i] > sm[best]:
                    best = i
                if raw[i] > 0:
                    tot += raw[i]
                    wsum += raw[i]
                    wmz += raw[i] * rmz[i]
            npf += hi - lo
            s_first[ns] = sc0 + lo
            s_last[ns] = sc0 + hi - 1
            s_apex[ns] = sc0 + best
            s_mz[ns] = wmz / wsum
            s_apexint[ns] = raw[best] if raw[best] > 0 else sm[best]
            s_sum[ns] = tot
            s_npts[ns] = npts
            s_poff[ns + 1] = npf
            s_trace[ns] = nt
            ns += 1
        t_seg1[nt] = ns
        nt += 1
    return (s_first[:ns], s_last[:ns], s_apex[:ns], s_mz[:ns], s_apexint[:ns], s_sum[:ns],
            s_npts[:ns], s_poff[:ns + 1], prof[:npf], s_trace[:ns],
            t_first[:nt], t_last[:nt], t_mz[:nt], t_poff[:nt + 1], t_prof[:ntp],
            t_seg0[:nt], t_seg1[:nt])


# --------------------------------------------------------------------------------------
# 3./4. isotope grouping
# --------------------------------------------------------------------------------------

@njit(cache=True)
def _averagine(mass, kmax):
    lam = mass * 0.000594  # Poisson approximation of the averagine isotope envelope
    p = np.empty(kmax)
    v = math.exp(-lam)
    for k in range(kmax):
        p[k] = v
        v = v * lam / (k + 1)
    return p


@njit(cache=True)
def _iso_match(i, tr, s_first, s_last, s_poff, prof, t_first, t_last, t_poff, t_prof):
    """Correlation, summed intensity and #points of trace ``tr`` within segment i's range."""
    lo = max(s_first[i], t_first[tr])
    hi = min(s_last[i], t_last[tr])
    if hi < lo:
        return -1.0, 0.0, 0
    n = s_last[i] - s_first[i] + 1
    sx = 0.0
    sy = 0.0
    sxx = 0.0
    syy = 0.0
    sxy = 0.0
    npt = 0
    for sc in range(s_first[i], s_last[i] + 1):
        x = prof[s_poff[i] + sc - s_first[i]]
        y = 0.0
        if lo <= sc <= hi:
            y = t_prof[t_poff[tr] + sc - t_first[tr]]
        if y > 0:
            npt += 1
        sx += x
        sy += y
        sxx += x * x
        syy += y * y
        sxy += x * y
    if n < 2:
        return 1.0, sy, npt
    vx = sxx - sx * sx / n
    vy = syy - sy * sy / n
    if vx <= 0 or vy <= 0:
        return 0.0, sy, npt
    return (sxy - sx * sy / n) / math.sqrt(vx * vy), sy, npt


@njit(cache=True)
def _seg_at(tr, sc, t_seg0, t_seg1, s_first, s_last):
    """Segment of trace tr that contains scan sc (or -1)."""
    for g in range(t_seg0[tr], t_seg1[tr]):
        if s_first[g] <= sc <= s_last[g]:
            return g
    return -1


@njit(cache=True)
def _group(s_first, s_last, s_apex, s_mz, s_sum, s_npts, s_poff, prof, t_first, t_last, t_mz,
           t_poff, t_prof, t_seg0, t_seg1, zlo, zhi, kmax, iso_ppm, min_corr, min_fit,
           min_iso, score_power, iso_min_points):
    ns = len(s_mz)
    nt = len(t_mz)
    torder = np.argsort(t_mz)
    tmz = t_mz[torder]
    cap = ns * 2 + 16
    h_mono = np.empty(cap, np.int64)
    h_z = np.empty(cap, np.int64)
    h_n = np.empty(cap, np.int64)
    h_score = np.empty(cap)
    h_fit = np.empty(cap)
    h_iso = np.empty((cap, kmax), np.int64)
    nh = 0
    chain = np.empty(kmax, np.int64)
    obs = np.empty(kmax)
    for i in range(ns):
        x0 = 0.0
        for sc in range(s_first[i], s_last[i] + 1):
            x0 += prof[s_poff[i] + sc - s_first[i]]
        for z in range(zlo, zhi + 1):
            n = 1
            obs[0] = x0
            chain[0] = -1
            for k in range(1, kmax):
                target = s_mz[i] + k * ISO_SPACING / z
                tol = target * iso_ppm * 1e-6 + k * 0.0006 / z
                lo = np.searchsorted(tmz, target - tol)
                best = -1
                best_c = min_corr
                best_y = 0.0
                for jj in range(lo, nt):
                    if tmz[jj] > target + tol:
                        break
                    tr = torder[jj]
                    c, y, npt = _iso_match(i, tr, s_first, s_last, s_poff, prof, t_first,
                                           t_last, t_poff, t_prof)
                    if npt < iso_min_points:
                        continue
                    if c >= best_c:
                        best_c = c
                        best = tr
                        best_y = y
                if best < 0:
                    break
                chain[k] = best
                obs[k] = best_y
                n += 1
            if n < min_iso:
                continue
            mass = (s_mz[i] - PROTON) * z
            th = _averagine(mass, n)
            num = 0.0
            dx = 0.0
            dy = 0.0
            tot = 0.0
            for k in range(n):
                num += obs[k] * th[k]
                dx += obs[k] * obs[k]
                dy += th[k] * th[k]
                tot += obs[k]
            fit = num / math.sqrt(dx * dy)
            if fit < min_fit:
                continue
            if nh >= cap:
                cap2 = cap * 2
                h_mono = np.concatenate((h_mono, np.empty(cap2 - cap, np.int64)))
                h_z = np.concatenate((h_z, np.empty(cap2 - cap, np.int64)))
                h_n = np.concatenate((h_n, np.empty(cap2 - cap, np.int64)))
                h_score = np.concatenate((h_score, np.empty(cap2 - cap)))
                h_fit = np.concatenate((h_fit, np.empty(cap2 - cap)))
                h_iso = np.concatenate((h_iso, np.empty((cap2 - cap, kmax), np.int64)))
                cap = cap2
            h_mono[nh] = i
            h_z[nh] = z
            h_n[nh] = n
            h_fit[nh] = fit
            h_score[nh] = tot * fit ** score_power
            for k in range(n):
                h_iso[nh, k] = chain[k]
            nh += 1
    # greedy acceptance; an elution peak is either a monoisotopic peak or (at the
    # feature's apex) an isotope of one accepted feature
    used = np.zeros(ns, np.bool_)
    owner = np.full(ns, -1, np.int64)
    owner_k = np.full(ns, -1, np.int64)
    ho = np.argsort(-h_score[:nh])
    f_mono = np.empty(nh, np.int64)
    f_z = np.empty(nh, np.int64)
    f_n = np.empty(nh, np.int64)
    f_fit = np.empty(nh)
    f_iso1 = np.empty(nh, np.int64)
    nf = 0
    iso_seg = np.empty(kmax, np.int64)
    for o in ho:
        i = h_mono[o]
        if used[i]:
            continue
        n = 1
        for k in range(1, h_n[o]):
            g = _seg_at(h_iso[o, k], s_apex[i], t_seg0, t_seg1, s_first, s_last)
            if g >= 0 and used[g]:
                break
            iso_seg[k] = g
            n += 1
        if n < min_iso:
            continue
        used[i] = True
        owner[i] = nf
        owner_k[i] = 0
        for k in range(1, n):
            g = iso_seg[k]
            if g >= 0:
                used[g] = True
                owner[g] = nf
                owner_k[g] = k
        f_mono[nf] = i
        f_z[nf] = h_z[o]
        f_n[nf] = n
        f_fit[nf] = h_fit[o]
        f_iso1[nf] = iso_seg[1] if n > 1 else -1
        nf += 1
    return f_mono[:nf], f_z[:nf], f_n[:nf], f_fit[:nf], f_iso1[:nf], owner, owner_k


def run_traces(ms1: dict[str, np.ndarray], params: dict[str, Any] | None = None,
               verbose: bool = True, debug: dict | None = None) -> list[dict[str, Any]]:
    """Detect features on flat MS1 arrays; returns records like _featuremap_records.

    If ``debug`` is a dict it receives the elution-peak table and segment ownership.
    """
    import time

    p = dict(DEFAULTS)
    p.update(params or {})
    t0 = time.time()
    mz, inten, offs, rt = ms1["mz"], ms1["int"], ms1["offs"], ms1["rt"]
    if p["min_intensity"] > 0:
        keep = inten >= p["min_intensity"]
        scan_all = np.repeat(np.arange(len(offs) - 1), np.diff(offs))
        mz, inten, scan_all = mz[keep], inten[keep], scan_all[keep]
        offs = np.r_[0, np.cumsum(np.bincount(scan_all, minlength=len(offs) - 1))]
    tid, ntr = _link(mz, inten, offs, float(p["trace_ppm"]), int(p["trace_gap"]))
    scan = np.repeat(np.arange(len(offs) - 1), np.diff(offs))
    order = np.lexsort((scan, tid))
    t1 = time.time()
    seg = _segments(order, tid, scan, mz, inten, ntr, float(p["valley_ratio"]), int(p["min_points"]),
                    int(p["smooth"]))
    (s_first, s_last, s_apex, s_mz, s_apexint, s_sum, s_npts, s_poff, prof, s_trace,
     t_first, t_last, t_mz, t_poff, t_prof, t_seg0, t_seg1) = seg
    t2 = time.time()
    f_mono, f_z, f_n, f_fit, f_iso1, owner, owner_k = _group(
        s_first, s_last, s_apex, s_mz, s_sum, s_npts, s_poff, prof, t_first, t_last, t_mz,
        t_poff, t_prof, t_seg0, t_seg1, int(p["charge_low"]), int(p["charge_high"]),
        int(p["max_isotopes"]), float(p["iso_ppm"]), float(p["min_corr"]), float(p["min_fit"]),
        int(p["min_isotopes"]), float(p["score_power"]), int(p["iso_min_points"]))
    t3 = time.time()
    if verbose:
        print(f"[features] traces: {ntr} traces ({t1 - t0:.0f}s), {len(s_mz)} elution peaks "
              f"({t2 - t1:.0f}s), {len(f_mono)} features ({t3 - t2:.0f}s)", flush=True)
    if debug is not None:
        debug.update(s_first=s_first, s_last=s_last, s_apex=s_apex, s_mz=s_mz,
                     s_apexint=s_apexint, s_sum=s_sum, s_npts=s_npts, owner=owner,
                     owner_k=owner_k, f_mono=f_mono, f_z=f_z, f_n=f_n, f_fit=f_fit)
    i1 = np.where(f_iso1 >= 0, f_iso1, f_mono)
    first = np.minimum(s_first[f_mono], s_first[i1])
    last = np.maximum(s_last[f_mono], s_last[i1])
    recs = []
    for k in range(len(f_mono)):
        m = f_mono[k]
        recs.append({
            "mz": float(s_mz[m]),
            "charge": int(f_z[k]),
            "rt_model": float(rt[s_apex[m]] if s_apexint[m] >= s_apexint[i1[k]] else rt[s_apex[i1[k]]]),
            "hull_rt_begin": float(rt[first[k]]),
            "hull_rt_end": float(rt[last[k]]),
            "mono_rt_begin": float(rt[s_first[m]]),
            "mono_rt_end": float(rt[s_last[m]]),
            "intensity": float(s_sum[m]),
            "quality": float(f_fit[k]),
            "n_isotopes": int(f_n[k]),
        })
    return recs
