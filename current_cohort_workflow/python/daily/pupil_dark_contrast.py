#!/usr/bin/env python3
"""
pupil_dark_contrast.py - repair of the package's pupil where it jumps to another dark blob, found with a strong
dark-range contrast and the glint (added 2026-11-18; user, on the day-5 QC sheet of 9670f green2: "this drawing is
pupil" (6 frames drawn next to the glint), "why not increase contrast a lot!!!!!! to get easier for pupil", then
"try this for day1,2,3,4,5 all"). An option for listed sessions, not a rule; the package's output is not changed.

Why a repair and not a new outline (pupil_dark_calibrate.py, 2026-11-18, the user's page drawings): where the
package sits on the right blob its outline matches the drawings (diameter ratio 1.01, centre 1.8 px); the
dark-contrast outline is smaller (ratio 0.75 at the half level). So the package's own detection is kept wherever it is
on the pupil, and only frames where it sits on another blob are changed.

Pass 1, per video frame, inside the session's eye box (crop as the package: imcrop [x y w h] = rows y..y+h,
cols x..x+w, 1-based):
  1. glint = the largest saturated blob (grey > 230, 4-900 px, not on the box edge) = the corneal reflection;
     it and other small bright specks are painted over (inpaint) and thin bright whiskers removed (grey opening 5 px);
  2. strong dark-range stretch: grey levels between the box's 0.5th and 35th percentile spread over 0-255;
  3. pupil candidates = round dark blobs within 50 px of the glint (75 px of the box centre without a glint):
     disc-minus-ring darkness at radii 5-15 px, local maxima, best 6 per frame, each outlined at the half level
     between its dark core and its ring.
Pass 2 (from the saved candidates; --rebuild redoes only this):
  4. the dark-contrast pupil of a frame = the candidate nearest to glint + the session's pupil-glint offset (median
     offset of the strongest candidates), within 14 px; without a glint, nearest to the previous one; else none;
  5. the package's per-frame detection = its accepted blob (hybrid_detections.csv, accepted_src > 0: the U-Net blob
     of highest probability, >= 0.9);
  6. per frame: both on the same blob (centres <= 12 px apart) -> the package's detection ('package'); no
     dark-contrast pupil -> the package's detection ('package_unchecked'); the package on another blob -> the
     dark-contrast pupil, its diameter x the session's size factor (median package / dark-contrast diameter over the
     same-blob frames) ('replaced'), only if that size is within +-25 % of the package's own same-blob diameter
     around that time (rolling median +-15 s; else no detection, 'other_blob_size_gate_none'); no package detection
     -> none ('none', the package's own rule; with --fill the dark-contrast pupil x factor, 'filled');
  7. the package's step-3 tracking (python copy of pupil_track_session.m: search radius 30 px around the last pupil
     or the warm-up anchor, reject a jump > 18 px or a size change > 18 %, hold the last pupil, EMA centre 0.18 /
     diameter 0.12, re-lock after 15 frames = this cohort's reacquire_frames) and its smoothing (linear gap fill +
     Savitzky-Golay 3, 11; edges held). The copy is checked on the package's own detections against its results2.csv
     (port_check in source.json). Frames whose tracked value is the same as the copy's on the package's own
     detections (= not touched by a replacement) keep the package's own values exactly (results2 and CLEAN).
Outputs: <pupil tree>/dayN/<cage>/<mouse>/pupil_outputs_darkcontrast/
  pupil_dark_candidates.csv, pupil_dark_glint.csv, pass1.json   pass 1
  pupil_dark_frames.csv   per frame: package detection, dark-contrast pupil, distance, action, merged detection
  pupil_dark_CLEAN.csv    package CLEAN columns (the package's own timestamps), from the merged detections
  source.json             parameters, offset, factor, % per action, port check, usable (below)
  dark_vs_package.png / dark_vs_package_minutes.csv
  check_sheet.png         8 replaced frames + 4 frames spread over the session (+ check_frames_*.png with --frames)
usable = glint in >= 80 % of frames, dark-contrast pupil in >= 80 %, and the two on the same blob in >= 50 % of the
frames where both exist, over the session AND in every minute with >= 300 such frames (else the dark-contrast pupil
is not the package's pupil in most frames or drifts to another dark area in part of the session: look first).
usage: python pupil_dark_contrast.py --day 5 --cage 9670f --mouse green2 [--all | --rebuild] [--fill] [--frames ...]
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.io import loadmat
from scipy.signal import savgol_filter

T = Path(__file__).resolve().parent
P = T.parent
DARK_LO, DARK_HI, GLINT_T = 0.5, 35.0, 230
R_LIST = (5, 7, 9, 11, 13, 15)
SEARCH_PX, PRIOR_PX, TOPK = 50, 14, 6
MIN_SCORE, MAX_CORE, MIN_CONTRAST = 20.0, 140.0, 25.0
FRAC = 0.5          # outline level between the blob's core and its ring (pupil_dark_calibrate.py)
# the package's step-3 tracking (pupil_track_session.m) and this cohort's reacquire_frames (qc.json)
TRACK = dict(search=30.0, warmup=30, max_jump=18.0, max_dfrac=0.18, a_c=0.18, a_d=0.12, reacq=15)
SAME_BLOB_PX = 12.0
MIN_SAME_FRAMES = 100
GATE, GATE_WIN = 0.25, 1801       # replacement size within +-25 % of the package's same-blob diameter (+-15 s)
USABLE = dict(glint_pct=80.0, found_pct=80.0, same_pct=50.0)
GEN = "pupil_outputs_general"


def _disk(r):
    k = np.zeros((2 * r + 1, 2 * r + 1), np.float32)
    cv2.circle(k, (r, r), r, 1, -1)
    return k


KERN = {}
for _r in R_LIST:
    _r2 = int(round(1.6 * _r)) + 3
    _D, _A = _disk(_r), _disk(_r2)
    _A[_r2 - (_r + 2):_r2 + (_r + 2) + 1, _r2 - (_r + 2):_r2 + (_r + 2) + 1] -= _disk(_r + 2)
    KERN[_r] = (_D / _D.sum(), _A / _A.sum())


def pupil_root() -> Path:
    sys.path.insert(0, str(T))
    import concat_choice as CC
    return CC.active_pupil_root(P)


def video_of(day, cage, mouse) -> Path:
    vs = sorted((P / "video" / f"day{day}" / cage / mouse).glob("*.mp4"))
    if not vs:
        raise SystemExit(f"no video for day{day} {cage} {mouse}")
    return vs[0]


def crop_of(day, cage, mouse):
    """package convention: imcrop(frame, [x y w h]) -> rows y..y+h, cols x..x+w (1-based) = (h+1) x (w+1)"""
    x, y, w, h = [int(round(v)) for v in loadmat(P / "video" / f"day{day}" / cage / mouse / "roi.mat")["roi"].ravel()]
    return x - 1, y - 1, w + 1, h + 1


def read_frames(video, crop, frames=None):
    """(1-based frame number, grey crop uint8), all frames or only the listed ones (read in order)"""
    x0, y0, cw, ch = crop
    cap = cv2.VideoCapture(str(video))
    want = set(frames) if frames is not None else None
    last = max(frames) if frames is not None else None
    k = 0
    while True:
        if last is not None and k >= last:
            break
        if not cap.grab():
            break
        k += 1
        if want is not None and k not in want:
            continue
        ok, fr = cap.retrieve()
        if ok:
            yield k, cv2.cvtColor(fr[y0:y0 + ch, x0:x0 + cw], cv2.COLOR_BGR2GRAY)
    cap.release()


def prep(C):
    """(stretched float image, glint (x, y, r) 0-based or None)"""
    H, W = C.shape
    n, lab, st, cen = cv2.connectedComponentsWithStats((C > GLINT_T).astype(np.uint8), connectivity=8)
    best = None
    for i in range(1, n):
        x, y, w, h, a = st[i]
        if x <= 0 or y <= 0 or x + w >= W or y + h >= H or a < 4 or a > 900:
            continue
        if best is None or a > st[best][4]:
            best = i
    glint = None
    if best is not None:
        glint = (float(cen[best][0]), float(cen[best][1]), float(np.sqrt(st[best][4] / np.pi)))
    n2, lab2, st2, _ = cv2.connectedComponentsWithStats((C > 200).astype(np.uint8), connectivity=8)
    small = np.isin(lab2, [i for i in range(1, n2) if st2[i][4] <= 900])
    paint = cv2.dilate(small.astype(np.uint8), np.ones((5, 5), np.uint8))
    Ci = cv2.inpaint(C, paint, 3, cv2.INPAINT_TELEA) if paint.any() else C
    Co = cv2.morphologyEx(Ci, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)))
    Cs = cv2.GaussianBlur(Co, (0, 0), 1.0).astype(np.float32)
    lo, hi = np.percentile(Cs, [DARK_LO, DARK_HI])
    return np.clip((Cs - lo) / max(1.0, hi - lo) * 255, 0, 255), glint


def segment(s, cx, cy, r, frac=None):
    """outline of the dark blob at (cx, cy) (0-based) of scale r: connected pixels darker than core + frac x
    (ring - core), holes filled; a stricter level (frac - 0.15) if it leaks to the window edge"""
    frac = FRAC if frac is None else frac
    H, W = s.shape
    rr = int(min(2.2 * r + 4, 32))
    x0, x1, y0, y1 = max(0, cx - rr), min(W, cx + rr + 1), max(0, cy - rr), min(H, cy + rr + 1)
    win = s[y0:y1, x0:x1]
    yy, xx = np.mgrid[y0:y1, x0:x1]
    dd = np.hypot(xx - cx, yy - cy)
    near = dd <= max(2.0, 0.5 * r)
    core = float(np.percentile(win[near], 20))
    ring = win[(dd >= r + 2) & (dd <= 1.6 * r + 3)]
    if ring.size < 10:
        return None
    iris = float(np.median(ring))
    if iris - core < MIN_CONTRAST:
        return None
    leak = 0
    for f in (frac, frac - 0.15):
        t = core + f * (iris - core)
        lab, n = ndi.label((win < t) & (dd <= rr))
        seed = np.argmin(np.where(near, win, np.inf))
        k = lab.ravel()[seed]
        if k == 0:
            return None
        M = ndi.binary_fill_holes(lab == k)
        leak = int((M & (dd >= rr - 1)).any())
        if not leak:
            break
    if leak:
        M = M & (dd <= 1.6 * r + 3)
    ys, xs = np.nonzero(M)
    a = int(M.sum())
    return dict(cx=float(xs.mean()) + x0 + 1, cy=float(ys.mean()) + y0 + 1, diam=2 * np.sqrt(a / np.pi), area=a,
                core=core, ring=iris, leak=leak, mask=(M, x0, y0))


def peaks(s, ref, search=SEARCH_PX):
    """round dark blobs within `search` px of ref (0-based): [(cx, cy, r, score)], strongest first"""
    H, W = s.shape
    m = int(search) + 30
    bx0, bx1 = max(0, int(ref[0]) - m), min(W, int(ref[0]) + m + 1)
    by0, by1 = max(0, int(ref[1]) - m), min(H, int(ref[1]) + m + 1)
    sb = s[by0:by1, bx0:bx1]
    S = np.full(sb.shape, -1e9, np.float32)
    Rm = np.zeros(sb.shape, np.int16)
    Dm = np.zeros(sb.shape, np.float32)
    for r, (D, A) in KERN.items():
        d = cv2.filter2D(sb, -1, D, borderType=cv2.BORDER_REPLICATE)
        sc = cv2.filter2D(sb, -1, A, borderType=cv2.BORDER_REPLICATE) - d
        u = sc > S
        S[u], Rm[u], Dm[u] = sc[u], r, d[u]
    yy, xx = np.mgrid[by0:by1, bx0:bx1]
    ok = ((xx - ref[0]) ** 2 + (yy - ref[1]) ** 2 <= search ** 2) & (S > MIN_SCORE) & (Dm < MAX_CORE)
    peak = ok & (S >= cv2.dilate(S, np.ones((7, 7), np.uint8)))
    ys, xs = np.nonzero(peak)
    return [(int(xs[j]) + bx0, int(ys[j]) + by0, int(Rm[ys[j], xs[j]]), float(S[ys[j], xs[j]]))
            for j in np.argsort(-S[ys, xs])]


def candidates(s, ref, search=SEARCH_PX):
    out = []
    for cx, cy, r, sc in peaks(s, ref, search)[:TOPK]:
        seg = segment(s, cx, cy, r)
        if seg:
            out.append(dict(score=sc, r=r, **seg))
    return out


def frame_candidates(C):
    s, gl = prep(C)
    ref = (gl[0], gl[1]) if gl else (s.shape[1] / 2, s.shape[0] / 2)
    return s, gl, candidates(s, ref, SEARCH_PX if gl else 1.5 * SEARCH_PX)


def choose(cands, glint1, offset, prev=None):
    """the candidate nearest to glint + offset (glint1 = 1-based glint centre; or nearest to the previous pupil when
    there is no glint), within PRIOR_PX -> (candidate or None, distance)"""
    if glint1 is not None and offset is not None:
        ref = (glint1[0] + offset[0], glint1[1] + offset[1])
    elif prev is not None:
        ref = prev
    else:
        return None, np.nan
    best, bd = None, np.inf
    for c in cands:
        d = float(np.hypot(c["cx"] - ref[0], c["cy"] - ref[1]))
        if d < bd:
            best, bd = c, d
    return (best, bd) if best is not None and bd <= PRIOR_PX else (None, bd)


def track_like_package(d, x, y, s=TRACK):
    """python copy of the step-3 loop of pupil_track_session.m: per-frame detection (diameter, centre; NaN = no blob)
    -> the held, EMA'd trace (Diameter_pixels, CenterX_px, CenterY_px) and the rejected frames. One blob per frame
    (maxp detector); the search radius drops a blob whose centre is > 30 px from the reference (the package cuts the
    blob at the radius, so a blob straddling it can differ slightly)."""
    n = len(d)
    od, ox, oy = np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan)
    rej = np.zeros(n, bool)
    nan2 = np.array([np.nan, np.nan])
    lc, ld, anchor, warm, miss = nan2.copy(), np.nan, nan2.copy(), [], 0
    nwarm = max(5, int(round(0.5 * s["warmup"])))
    for k in range(n):
        has = bool(np.isfinite(d[k]) and np.isfinite(x[k]) and np.isfinite(y[k]))
        c0 = np.array([x[k], y[k]], float)
        if has:
            if len(warm) < s["warmup"]:
                warm.append(c0)
            if len(warm) >= nwarm:
                anchor = np.median(np.array(warm), axis=0)
        ref = lc if np.all(np.isfinite(lc)) else anchor
        if not np.all(np.isfinite(ref)) and has:
            ref = c0
        found = has and not (np.all(np.isfinite(ref)) and np.hypot(*(c0 - ref)) > s["search"])
        if found:
            c, dk = c0, float(d[k])
            acc = not (np.all(np.isfinite(lc)) and np.hypot(*(c - lc)) > s["max_jump"])
            if acc and np.isfinite(ld) and abs(dk - ld) > s["max_dfrac"] * max(1.0, ld):
                acc = False
            if not acc:
                rej[k] = True
                c, dk = (lc.copy(), ld) if np.isfinite(ld) else (nan2.copy(), np.nan)
        else:
            c, dk = (lc.copy(), ld) if np.isfinite(ld) else (nan2.copy(), np.nan)
        if np.all(np.isfinite(c)):
            if np.all(np.isfinite(lc)):
                c = (1 - s["a_c"]) * lc + s["a_c"] * c
            lc = c
        if np.isfinite(dk):
            if np.isfinite(ld):
                dk = (1 - s["a_d"]) * ld + s["a_d"] * dk
            ld = dk
        od[k] = ld
        ox[k], oy[k] = lc
        if s["reacq"] > 0:
            miss = 0 if (found and not rej[k]) else miss + 1
            if miss >= s["reacq"]:
                lc, ld, miss, anchor, warm = nan2.copy(), np.nan, 0, nan2.copy(), []
    return od, ox, oy, rej


def smooth_like_package(d):
    """fillmissing(raw, 'linear') + sgolayfilt(., 3, 11) of pupil_track_session.m (edges held, not extrapolated)"""
    v = pd.Series(d).interpolate(limit_area="inside").ffill().bfill().to_numpy()
    if not np.isfinite(v).any():
        v = np.zeros_like(v)
    return savgol_filter(v, 11, 3) if len(v) >= 11 else v


def session_dir_of(a):
    return pupil_root() / f"day{a.day}" / a.cage / a.mouse


def pass1(video, crop, out_dir, quiet=False):
    rows, glints, nfr = [], [], 0
    for k, C in read_frames(video, crop):
        nfr = k
        _, gl, cs = frame_candidates(C)
        glints.append((k,) + ((gl[0] + 1, gl[1] + 1, gl[2]) if gl else (np.nan, np.nan, np.nan)))
        for j, c in enumerate(cs):
            c.pop("mask")
            rows.append(dict(frame=k, k=j, glint_x=gl[0] + 1 if gl else np.nan, glint_y=gl[1] + 1 if gl else np.nan,
                             **c))
        if k % 5000 == 0 and not quiet:
            print(f"  frame {k}", flush=True)
    pd.DataFrame(rows).to_csv(out_dir / "pupil_dark_candidates.csv", index=False, float_format="%.6g")
    pd.DataFrame(glints, columns=["frame", "glint_x", "glint_y", "glint_r"]).to_csv(
        out_dir / "pupil_dark_glint.csv", index=False, float_format="%.6g")
    (out_dir / "pass1.json").write_text(json.dumps(dict(frames=nfr, video=str(video), crop_x0y0wh=list(crop),
                                                        made=datetime.now().isoformat(timespec="seconds")), indent=1),
                                        encoding="utf-8")
    return nfr


def pass2(a, sd, out_dir, fill=False, quiet=False):
    p1 = json.loads((out_dir / "pass1.json").read_text(encoding="utf-8"))
    nfr = int(p1["frames"])
    gen = sd / GEN
    R2 = pd.read_csv(gen / "pupil_diameter_FINAL_U-Net_results2.csv")
    PC = pd.read_csv(gen / "pupil_diameter_FINAL_U-Net_results_CLEAN.csv")
    HD = pd.read_csv(gen / "hybrid_detections.csv")
    if len(R2) != nfr or len(HD) != nfr or len(PC) != nfr or \
            not np.array_equal(HD.frame.to_numpy(), np.arange(1, nfr + 1)):
        raise SystemExit(f"frame count {nfr} differs from the package's rows (results2 {len(R2)}, detections "
                         f"{len(HD)}, CLEAN {len(PC)}) - stop (the frames would not line up)")
    acc = HD["accepted_src_1unet_2dark_0none"].to_numpy() > 0
    pdm, pxm, pym = (np.where(acc, HD[c].to_numpy(float), np.nan) for c in ("best_diam_px", "best_cx_px", "best_cy_px"))
    # the python copy of step 3 on the package's own detections must give the package's own trace
    od0, ox0, oy0, _ = track_like_package(pdm, pxm, pym)
    e_d = np.abs(od0 - R2.Diameter_pixels.to_numpy(float))
    e_c = np.hypot(ox0 - R2.CenterX_px.to_numpy(float), oy0 - R2.CenterY_px.to_numpy(float))
    both = np.isfinite(e_d)
    port = dict(frames_diam_within_0p01px_pct=round(100 * float(np.mean(e_d[both] <= 0.01)), 2) if both.any() else None,
                diam_median_abs_diff_px=round(float(np.median(e_d[both])), 4) if both.any() else None,
                diam_p99_abs_diff_px=round(float(np.percentile(e_d[both], 99)), 3) if both.any() else None,
                centre_median_diff_px=round(float(np.nanmedian(e_c)), 4),
                nan_mismatch_frames=int(np.sum(np.isfinite(od0) != np.isfinite(R2.Diameter_pixels.to_numpy(float)))))
    # the dark-contrast pupil per frame
    Cd = pd.read_csv(out_dir / "pupil_dark_candidates.csv")
    gp = out_dir / "pupil_dark_glint.csv"
    if gp.exists():
        G = pd.read_csv(gp)
        gl_x, gl_y = G.glint_x.to_numpy(float), G.glint_y.to_numpy(float)
        glint_pct = 100 * float(np.isfinite(gl_x).mean())
        glint_src = "every frame (pupil_dark_glint.csv)"
    else:
        gl_x = np.full(nfr, np.nan)
        gl_y = np.full(nfr, np.nan)
        g0 = Cd.drop_duplicates("frame")
        gl_x[g0.frame.to_numpy() - 1], gl_y[g0.frame.to_numpy() - 1] = g0.glint_x, g0.glint_y
        glint_pct = 100 * float(np.isfinite(gl_x[np.isin(np.arange(1, nfr + 1), g0.frame)]).mean())
        glint_src = "frames with a candidate only (pass 1 without pupil_dark_glint.csv)"
    top = Cd[(Cd.k == 0) & Cd.glint_x.notna()]
    top = top[top.score >= top.score.median()]
    offset = (float((top.cx - top.glint_x).median()), float((top.cy - top.glint_y).median())) if len(top) else None
    by = {k: g.to_dict("records") for k, g in Cd.groupby("frame")}
    nd, nx, ny, dprior = (np.full(nfr, np.nan) for _ in range(4))
    prev = None
    for k in range(1, nfr + 1):
        g1 = (gl_x[k - 1], gl_y[k - 1]) if np.isfinite(gl_x[k - 1]) else None
        c, dist = choose(by.get(k, []), g1, offset, prev)
        dprior[k - 1] = dist
        if c is not None:
            prev = (c["cx"], c["cy"])
            nd[k - 1], nx[k - 1], ny[k - 1] = c["diam"], c["cx"], c["cy"]
    found = np.isfinite(nd)
    dist = np.hypot(pxm - nx, pym - ny)
    same = acc & found & (dist <= SAME_BLOB_PX)
    other = acc & found & (dist > SAME_BLOB_PX)
    ratio = pdm[same] / nd[same]
    factor = float(np.median(ratio)) if same.sum() >= MIN_SAME_FRAMES else np.nan
    act = np.where(same, "package", np.where(acc & ~found, "package_unchecked", np.where(other, "replaced", "none")))
    md, mx, my = pdm.copy(), pxm.copy(), pym.copy()
    # size gate for a replacement: within +-GATE of the package's own same-blob diameter around that time (rolling
    # median over +-15 s); a dark-contrast blob cut by a whisker would otherwise put a dip into the trace
    ref_d = pd.Series(np.where(same, pdm, np.nan)).rolling(GATE_WIN, center=True, min_periods=30).median() \
        .interpolate(limit_area="inside").ffill().bfill().to_numpy()
    if np.isfinite(factor):
        cand = nd * factor
        ok_size = np.abs(cand / ref_d - 1) <= GATE
        rep = other & ok_size
        md[rep], mx[rep], my[rep] = cand[rep], nx[rep], ny[rep]
        bad = other & ~ok_size
        md[bad], mx[bad], my[bad] = np.nan, np.nan, np.nan
        act = np.where(bad, "other_blob_size_gate_none", act)
        if fill:
            fl = ~acc & found & ok_size
            md[fl], mx[fl], my[fl] = cand[fl], nx[fl], ny[fl]
            act = np.where(fl, "filled", act)
    else:
        act = np.where(other, "package_other_blob_kept", act)
    act = np.where(~acc & found & (act == "none"), "none_dark_found", act)
    od, ox, oy, rej = track_like_package(md, mx, my)
    # frames whose tracked value differs from the copy run on the package's own detections = changed by the repair;
    # every other frame keeps the package's own value exactly (so the copy's small differences never enter)
    chg = ~(np.isclose(od, od0, equal_nan=True, rtol=0, atol=1e-9) & np.isclose(ox, ox0, equal_nan=True, rtol=0,
                                                                                    atol=1e-9)
            & np.isclose(oy, oy0, equal_nan=True, rtol=0, atol=1e-9))
    od = np.where(chg, od, R2.Diameter_pixels.to_numpy(float))
    ox = np.where(chg, ox, R2.CenterX_px.to_numpy(float))
    oy = np.where(chg, oy, R2.CenterY_px.to_numpy(float))
    sm = smooth_like_package(od)
    near = np.convolve(chg.astype(float), np.ones(11), mode="same") > 0      # within 5 frames of a changed frame
    sm = np.where(near, sm, PC.SmoothedDiameter_px.to_numpy(float))
    F = pd.DataFrame(dict(frame=np.arange(1, nfr + 1), glint_x=gl_x, glint_y=gl_y, pkg_d=pdm, pkg_cx=pxm, pkg_cy=pym,
                          dark_d=nd, dark_cx=nx, dark_cy=ny, dark_dist_prior_px=dprior, dist_pkg_dark_px=dist,
                          action=act, merged_d=md, merged_cx=mx, merged_cy=my, tracked_rejected=rej.astype(int),
                          trace_changed=chg.astype(int)))
    F.to_csv(out_dir / "pupil_dark_frames.csv", index=False, float_format="%.6g")
    C = pd.DataFrame(dict(Timestamp_seconds=R2.Timestamp_seconds.to_numpy(), Diameter_pixels=od, CenterX_px=ox,
                          CenterY_px=oy, UsedScoreFallback=R2.UsedScoreFallback.to_numpy(), SmoothedDiameter_px=sm))
    C.to_csv(out_dir / "pupil_dark_CLEAN.csv", index=False, float_format="%.10g")
    pct = {f"{k}_pct": round(100 * float(np.mean(act == k)), 2) for k in
           ("package", "package_unchecked", "replaced", "other_blob_size_gate_none", "filled", "none",
            "none_dark_found", "package_other_blob_kept")}
    both_n = int((acc & found).sum())
    same_pct = round(100 * float(same.sum()) / both_n, 2) if both_n else 0.0
    found_pct = round(100 * float(found.mean()), 2)
    # the same check minute by minute: a drift to another dark area in part of the session (e.g. a second glint
    # taking over) must not pass because the other part agrees (day 2 9670f white, 2026-11-18)
    minute = (R2.Timestamp_seconds.to_numpy(float) // 60).astype(int)
    mins = pd.DataFrame(dict(m=minute, both=acc & found, same=same)).groupby("m").sum()
    mins = mins[mins.both >= 300]
    min_same = round(float((100 * mins.same / mins.both).min()), 2) if len(mins) else 0.0
    usable = bool(glint_pct >= USABLE["glint_pct"] and found_pct >= USABLE["found_pct"] and
                  same_pct >= USABLE["same_pct"] and min_same >= USABLE["same_pct"] and np.isfinite(factor))
    why = [w for w, bad in (
        (f"glint in only {glint_pct:.0f} % of frames", glint_pct < USABLE["glint_pct"]),
        (f"dark-contrast pupil in only {found_pct:.0f} % of frames", found_pct < USABLE["found_pct"]),
        (f"same blob as the package in only {same_pct:.0f} % of frames", same_pct < USABLE["same_pct"]),
        (f"same blob in only {min_same:.0f} % of frames in its worst minute", min_same < USABLE["same_pct"]),
        ("too few same-blob frames for the size factor", not np.isfinite(factor))) if bad]
    pkc = PC.SmoothedDiameter_px.to_numpy(float)
    new = C.SmoothedDiameter_px.to_numpy(float)
    info = dict(tool="pupil_dark_contrast.py", made=datetime.now().isoformat(timespec="seconds"),
                session=f"day{a.day}_{a.cage}_{a.mouse}", video=p1.get("video"), frames=nfr,
                usable=usable, not_usable_because=why, glint_pct=round(glint_pct, 2), glint_from=glint_src,
                dark_found_pct=found_pct, same_blob_pct_of_both=same_pct, same_blob_pct_worst_minute=min_same,
                offset_px=[round(offset[0], 2), round(offset[1], 2)] if offset else None,
                size_factor_package_over_dark=round(factor, 4) if np.isfinite(factor) else None,
                size_factor_iqr=[round(float(np.percentile(ratio, 25)), 3), round(float(np.percentile(ratio, 75)), 3)]
                if len(ratio) else None,
                actions=pct, package_accepted_pct=round(100 * float(acc.mean()), 2),
                tracked_rejected_pct=round(100 * float(rej.mean()), 2),
                median_diam_px=dict(package=round(float(np.nanmedian(pkc)), 2), repaired=round(float(np.nanmedian(new)), 2)),
                frames_changed_gt_1px_pct=round(100 * float(np.nanmean(np.abs(new - pkc) > 1.0)), 2),
                trace_changed_pct=round(100 * float(chg.mean()), 2),
                unchanged_frames_identical_to_package=bool(np.array_equal(
                    np.where(chg, 0, od), np.where(chg, 0, R2.Diameter_pixels.to_numpy(float)), equal_nan=True)),
                port_check=port, fill=bool(fill),
                settings=dict(dark_range_pct=[DARK_LO, DARK_HI], glint_t=GLINT_T, radii=list(R_LIST),
                              search_px=SEARCH_PX, prior_px=PRIOR_PX, frac=FRAC, same_blob_px=SAME_BLOB_PX,
                              tracking=TRACK, usable_rule=USABLE),
                note="option for listed sessions (user 2026-11-18: pupil next to the glint; strong contrast)")
    (out_dir / "source.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
    compare(a, out_dir, C, PC, F)
    if not quiet:
        print(f"{info['session']}: usable {usable} {why}; replaced {pct['replaced_pct']} % of frames, same blob "
              f"{same_pct} %, factor {info['size_factor_package_over_dark']}, median {info['median_diam_px']}; "
              f"port check {port}")
    return info


def compare(a, out_dir, C, pk, F):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = pk.Timestamp_seconds.to_numpy()
    act = F.action.to_numpy()
    per = pd.DataFrame(dict(minute=(t // 60).astype(int), rep=(act == "replaced") | (act == "filled"),
                            same=act == "package", unchecked=act == "package_unchecked",
                            pk=pk.SmoothedDiameter_px, nw=C.SmoothedDiameter_px))
    g = per.groupby("minute").agg(replaced_pct=("rep", lambda v: 100 * v.mean()),
                                  same_blob_pct=("same", lambda v: 100 * v.mean()),
                                  unchecked_pct=("unchecked", lambda v: 100 * v.mean()),
                                  package_median_px=("pk", "median"), repaired_median_px=("nw", "median"))
    g.to_csv(out_dir / "dark_vs_package_minutes.csv", float_format="%.4g")
    s10 = pd.DataFrame(dict(b=(t // 10).astype(int), pk=pk.SmoothedDiameter_px, nw=C.SmoothedDiameter_px)).groupby(
        "b").median()
    fig, ax = plt.subplots(3, 1, figsize=(13, 8.5), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1.3, 1.3]))
    ax[0].plot(t / 60, pk.SmoothedDiameter_px, color="0.75", lw=0.4, label="package per frame")
    ax[0].plot(t / 60, C.SmoothedDiameter_px, color="#8fd3ea", lw=0.4, label="repaired per frame")
    ax[0].plot((s10.index * 10 + 5) / 60, s10.pk, color="k", lw=1.4, label="package (10-s median)")
    ax[0].plot((s10.index * 10 + 5) / 60, s10.nw, color="#0088b8", lw=1.6, label="repaired (10-s median)")
    ax[0].set_ylabel("pupil diameter (px)")
    ax[0].legend(loc="upper right", fontsize=8, ncol=4)
    ax[0].set_title(f"{a.cage} {a.mouse} day {a.day}: package pupil vs repaired (package kept where it is on the pupil; "
                    f"replaced where it sat on another blob: {100 * np.mean(act == 'replaced'):.1f} % of frames)",
                    fontsize=10)
    ax[1].plot(t / 60, F.dist_pkg_dark_px, color="#7a4fb0", lw=0.4)
    ax[1].axhline(SAME_BLOB_PX, color="r", lw=0.8, ls="--")
    ax[1].set_ylabel("package - dark-contrast\ncentre distance (px)")
    ax[1].set_ylim(0, 80)
    ax[2].bar(g.index + 0.5, g.replaced_pct, width=0.9, color="#e08a2c", label="replaced (% frames)")
    ax[2].plot(g.index + 0.5, g.same_blob_pct, "o-", color="#0088b8", ms=3, label="package on the same blob (% frames)")
    ax[2].plot(g.index + 0.5, g.unchecked_pct, "s-", color="0.5", ms=3, label="no dark-contrast pupil (% frames)")
    ax[2].set_ylim(0, 105)
    ax[2].set_ylabel("% of frames")
    ax[2].set_xlabel("session time (min)")
    ax[2].legend(loc="center right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "dark_vs_package.png", dpi=110)
    plt.close(fig)


def sheet(a, video, crop, frames, out_png, pkg, offset, clean=None):
    tiles, prev = [], None
    for k, Cg in read_frames(video, crop, sorted(frames)):
        s, gl, cs = frame_candidates(Cg)
        if offset is not None:
            c, _ = choose(cs, (gl[0] + 1, gl[1] + 1) if gl else None, offset, prev)
        else:
            c = cs[0] if cs else None
        lo, hi = np.percentile(Cg, [2, 98])
        plain = np.clip((Cg.astype(np.float32) - lo) / max(1.0, hi - lo) * 255, 0, 255).astype(np.uint8)
        pair = []
        for im in (plain, s.astype(np.uint8)):
            v = cv2.cvtColor(cv2.resize(im, (im.shape[1] * 2, im.shape[0] * 2), interpolation=cv2.INTER_CUBIC),
                             cv2.COLOR_GRAY2BGR)
            for tab, col in ((pkg, (0, 230, 255)), (clean, (60, 220, 60))):
                if tab is not None and k <= len(tab):
                    r = tab.iloc[k - 1]
                    if np.isfinite(r.Diameter_pixels) and np.isfinite(r.CenterX_px):
                        cv2.circle(v, (int(round((r.CenterX_px - 1) * 2)), int(round((r.CenterY_px - 1) * 2))),
                                   int(round(r.Diameter_pixels)), col, 1, lineType=cv2.LINE_AA)
            if c is not None:
                M, x0, y0 = c["mask"]
                big = np.zeros(im.shape, np.uint8)
                big[y0:y0 + M.shape[0], x0:x0 + M.shape[1]] = M
                cnt, _ = cv2.findContours(cv2.resize(big * 255, (v.shape[1], v.shape[0]),
                                                     interpolation=cv2.INTER_NEAREST),
                                          cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
                cv2.drawContours(v, cnt, -1, (255, 220, 0), 1)
            pair.append(v)
        txt = f"f{k}"
        if pkg is not None and k <= len(pkg):
            txt += f"  package {pkg.iloc[k - 1].Diameter_pixels:.1f}"
        if clean is not None and k <= len(clean):
            txt += f"  repaired {clean.iloc[k - 1].Diameter_pixels:.1f} px"
        t = np.hstack([pair[0], np.full((pair[0].shape[0], 6, 3), 255, np.uint8), pair[1]])
        cv2.putText(t, txt, (6, 22), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 0, 255), 2)
        tiles.append(t)
        if c is not None:
            prev = (c["cx"], c["cy"])
    if not tiles:
        return
    cols = 2
    while len(tiles) % cols:
        tiles.append(np.full_like(tiles[0], 255))
    gapc = np.full((tiles[0].shape[0], 16, 3), 255, np.uint8)
    rows = [np.hstack([tiles[i], gapc, tiles[i + 1]]) for i in range(0, len(tiles), cols)]
    gapr = np.full((12, rows[0].shape[1], 3), 255, np.uint8)
    img = np.vstack([x for r in rows for x in (r, gapr)])
    bar = np.full((36, img.shape[1], 3), 255, np.uint8)
    cv2.putText(bar, f"{a.cage} {a.mouse} day {a.day}: left = as the QC sheet, right = strong dark contrast;  "
                     "yellow = package, green = repaired, cyan = dark-contrast blob", (6, 25),
                cv2.FONT_HERSHEY_SIMPLEX, 0.65, (0, 0, 0), 2)
    cv2.imwrite(str(out_png), np.vstack([bar, img]))


def sheet_frames(F, n_rep=8, n_spread=4):
    """frames for the check sheet: replaced frames spread over the session (one per replaced stretch) + spread frames"""
    rep = F.frame[F.action == "replaced"].to_numpy()
    ks = []
    if len(rep):
        starts = rep[np.r_[True, np.diff(rep) > 30]]
        ks += [int(v) for v in starts[np.linspace(0, len(starts) - 1, min(n_rep, len(starts))).astype(int)] + 3]
    ks += [int(v) for v in np.linspace(1, len(F), n_spread + 2)[1:-1]]
    return sorted(set(min(max(1, k), len(F)) for k in ks))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--day", type=int, required=True)
    ap.add_argument("--cage", required=True)
    ap.add_argument("--mouse", required=True)
    ap.add_argument("--frames", type=int, nargs="*", default=None)
    ap.add_argument("--all", action="store_true", help="pass 1 (every frame) + pass 2")
    ap.add_argument("--rebuild", action="store_true", help="pass 2 only, from the saved candidates")
    ap.add_argument("--fill", action="store_true", help="also fill frames without a package detection")
    ap.add_argument("--no-sheet", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    cv2.setNumThreads(2)
    video, crop = video_of(a.day, a.cage, a.mouse), crop_of(a.day, a.cage, a.mouse)
    sd = session_dir_of(a)
    out_dir = sd / "pupil_outputs_darkcontrast"
    out_dir.mkdir(parents=True, exist_ok=True)
    if a.all:
        pass1(video, crop, out_dir, a.quiet)
    info = pass2(a, sd, out_dir, a.fill, a.quiet) if (a.all or a.rebuild) else None
    sj = out_dir / "source.json"
    if info is None and sj.exists():
        info = json.loads(sj.read_text(encoding="utf-8"))
    offset = tuple(info["offset_px"]) if info and info.get("offset_px") else None
    pkg = pd.read_csv(sd / GEN / "pupil_diameter_FINAL_U-Net_results2.csv")
    clean = pd.read_csv(out_dir / "pupil_dark_CLEAN.csv") if (out_dir / "pupil_dark_CLEAN.csv").exists() else None
    if (a.all or a.rebuild) and not a.no_sheet:
        F = pd.read_csv(out_dir / "pupil_dark_frames.csv", usecols=["frame", "action"])
        sheet(a, video, crop, sheet_frames(F), out_dir / "check_sheet.png", pkg, offset, clean)
    if a.frames:
        png = out_dir / f"check_frames_{'_'.join(map(str, a.frames[:6]))}.png"
        sheet(a, video, crop, a.frames, png, pkg, offset, clean)
        print(png)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
