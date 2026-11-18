#!/usr/bin/env python3
"""
tst_mouse_lock.py - TST automatic scoring of the FIRST 3 MIN that keeps following the MOUSE when the camera view
changes inside a video (added 2026-11-18). The user: "try to fix it... always only mouse tracking" (9669m green2 day 5:
the camera moved in the first 20 s and the package's scorer followed another dark blob there, scoring those seconds as
non-moving); "TST should be 3min so just use first 3min"; "always mouse is more less upper middle".

The package scorer (tools/auto_pupil_tst_pipeline_20261114_iris_on/tst) is NOT changed. It finds the mouse with
per-video maps made from 300 frames spread over the video (dark threshold, always-dark mask, 'where the mouse usually
is'). When the camera is moved or zoomed during a recording, the frames before / during the move do not fit those
maps. This tool runs the package's own scoring code with three things around it:
  1. FIRST 3 MIN: the scorer is only given the first --max-s seconds (default 180) of the video: frames after that are
     not decoded, and the frame count / numFrames / duration it uses are those of the first 180 s. Shorter videos are
     used whole (the package's QC then warns about the length).
  2. CAMERA VIEW FIX, ONLY WHEN NEEDED: the package scorer first runs on the untouched first 3 min and the mouse QC
     (3.) checks every frame. Only when more than 1 % of the frames, or most frames of any second, are not the mouse
     is the view fix tried (a zoomed view with a plain white background has too few key points for a reliable
     alignment, so it is never applied without need - day-5 test 2026-11-18: 9669m white/black were tracked 100 %
     correctly untouched, while forced alignment moved white's result by 4 pp from noise):
     ORB key points of frames sampled at 2 fps are matched to one reference frame of the main view (RANSAC; rotation
     + scale + shift). A frame is 'moved' when it was aligned reliably (>= 30 inliers) and its corners land more than
     2 px (at the scorer's 320 px width) from the reference view; samples that could not be aligned count only next
     to moved ones. Every frame of a moved stretch (+- 0.5 s) is aligned one by one (frames that cannot be aligned get
     the transform interpolated from their neighbours) and warped into the main view (uncovered border = the light
     background grey). The fixed run is used only if its QC passes and is better; otherwise the untouched result is
     kept and the video is reported as NEEDS A LOOK.
  3. MOUSE QC on every frame, from the scorer's own per-frame features. 'Not mouse' when the tracked blob's area is
     outside 0.4-2.5 x the video's reference mouse area (a_ref), its centre jumps more than one mouse length (l1)
     from its 3-s running median, or its centre is outside the upper middle of the frame (x 10-80 % of the width,
     y in the top 60 %). Calibrated 2026-11-18 on 3 day-5 videos (mouse centre x 27-37 %, y 11-34 %): normal videos
     0.0-0.3 % of frames fail; the broken first 20 s of 9669m green2 fail 76 % on area/jump alone. Failing frames are
     listed, never hidden or dropped.
  A check sheet of the corrected stretch shows the blob the package maps pick (red, original frame) above the blob
  after the fix (green, aligned frame).
Outputs in --out:
  AutoScoringResults_<video>.mat / _summary.csv   the CHOSEN run (package format, 0 = Non-moving, 1 = Moving);
                         untouched_<max_s>s/ and (when tried) viewfix_<max_s>s/ hold both runs
  <video>_lock.json      max_s, chosen run and why, % non-moving and mouse QC of each run, the view fix details
  <video>_view.csv       per aligned frame: frame, t_s, moved, scale, rot_deg, dx_px, dy_px (640 px units), inliers
  <video>_mouse_qc.csv   per frame: frame, t_s, valid, area_ratio, jump_l1, cx_frac, cy_frac, not_mouse
  <video>_lock_check.png only when the camera moved
usage: python tst_mouse_lock.py --video <mp4> --out <folder> [--max-s 180] [--compare-plain]
  --compare-plain  also score the first 3 min with untouched frames (plain package maps) for the before/after %.
Runs at Windows idle priority, one video at a time, 2 threads (like one package worker). Feature caches go to
<out>/cache_*_<max_s>s and are never shared with a run of other settings.
"""
from __future__ import annotations

import argparse
import ctypes
import json
import math
import sys
import types
from pathlib import Path

import cv2
import numpy as np
import pandas as pd

T = Path(__file__).resolve().parent
TST = T / "auto_pupil_tst_pipeline_20261114_iris_on" / "tst"
sys.path.insert(0, str(TST))
import tst_features as TF  # noqa: E402
import tst_auto_score as TAS  # noqa: E402

ORB_W = 640                        # alignment works on 640 px wide grey frames
MOVE_PX = 4.0                      # corner shift at 640 px width = 2 px at the scorer's 320 px
FIX_IF_NOT_MOUSE_PCT = 1.0         # the view fix is tried only when the untouched run has more not-mouse frames than
                                   # this, or any second in which most frames are not the mouse
SAMPLE_S = 0.5                     # 2 fps sampling for the view scan
MIN_INL = 30                       # fewer RANSAC inliers = alignment not trusted
AREA_LO, AREA_HI = 0.4, 2.5        # QC: tracked blob area vs the video's reference mouse area
X_LO, X_HI, Y_HI = 0.10, 0.80, 0.60   # QC: 'upper middle' (fractions of width / height)
REAL_VC = cv2.VideoCapture
ORIG_NUM_FRAMES = TAS.vm.manual_num_frames


def idle():
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x40)
    except Exception:  # noqa: BLE001 - not Windows
        pass


def small(frame):
    h, w = frame.shape[:2]
    s = ORB_W / w
    return cv2.cvtColor(cv2.resize(frame, (ORB_W, int(round(h * s))), interpolation=cv2.INTER_AREA),
                        cv2.COLOR_BGR2GRAY), s


ORB = cv2.ORB_create(nfeatures=2000)
BF = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)


def keypoints(g):
    k, d = ORB.detectAndCompute(g, None)
    return (np.float32([p.pt for p in k]) if k else np.zeros((0, 2), np.float32)), d


def align(kd, ref_kd):
    """2x3 similarity mapping frame pixels -> reference pixels (640 px units) and the inlier count; None if unsure"""
    (p, d), (p0, d0) = kd, ref_kd
    if d is None or d0 is None or len(p) < MIN_INL:
        return None, 0
    m = BF.match(d, d0)
    if len(m) < MIN_INL:
        return None, len(m)
    src = p[[x.queryIdx for x in m]]
    dst = p0[[x.trainIdx for x in m]]
    A, inl = cv2.estimateAffinePartial2D(src, dst, method=cv2.RANSAC, ransacReprojThreshold=2.0, maxIters=3000,
                                         confidence=0.995)
    n = int(inl.sum()) if inl is not None else 0
    if A is None or n < MIN_INL:
        return None, n
    return A, n


def is_moved(A, shape):
    h, w = shape
    c = np.float32([[0, 0], [w, 0], [0, h], [w, h]])
    return float(np.abs(c @ A[:, :2].T + A[:, 2] - c).max()) > MOVE_PX


def params(A):
    s = float(np.sqrt(abs(np.linalg.det(A[:, :2]))))
    return s, float(np.degrees(np.arctan2(A[1, 0], A[0, 0]))), float(A[0, 2]), float(A[1, 2])


def scan_views(video: Path, max_frames: int):
    """2-fps samples of the first max_frames -> reference sample, per-sample transforms; then every frame of the
    moved stretches"""
    cap = REAL_VC(str(video))
    fps = cap.get(cv2.CAP_PROP_FPS) or 30.0
    step = max(1, int(round(fps * SAMPLE_S)))
    samples, kds, k, shape, scale = {}, {}, 0, None, 1.0
    while k < max_frames:
        if not cap.grab():
            break
        k += 1
        if (k - 1) % step == 0:
            ok, fr = cap.retrieve()
            if ok:
                g, scale = small(fr)
                shape = g.shape
                samples[k] = g
                kds[k] = keypoints(g)
    cap.release()
    n_frames = k
    ks = sorted(samples)
    best = None
    for frac in (0.5, 0.25, 0.75, 0.1, 0.9):              # reference = the sample that most others match
        kr = ks[int(frac * (len(ks) - 1))]
        res = {kk: align(kds[kk], kds[kr]) for kk in ks}
        n_same = sum(1 for A, _ in res.values() if A is not None and not is_moved(A, shape))
        if best is None or n_same > best[0]:
            best = (n_same, kr, res)
        if n_same >= 0.5 * len(ks):
            break
    n_same, kr, res = best
    # a sample is 'moved' only when it was aligned reliably AND its view differs; samples that could not be aligned
    # (too few key points) count only when they touch moved samples (e.g. the blurred / zoomed first seconds)
    status = {kk: ("unknown" if res[kk][0] is None else "moved" if is_moved(res[kk][0], shape) else "same")
              for kk in ks}
    moved_s = set(kk for kk in ks if status[kk] == "moved")
    grow = True
    while grow:
        grow = False
        for i, kk in enumerate(ks):
            if status[kk] == "unknown" and kk not in moved_s and (
                    (i > 0 and ks[i - 1] in moved_s) or (i + 1 < len(ks) and ks[i + 1] in moved_s)):
                moved_s.add(kk)
                grow = True
    moved_s = sorted(moved_s)
    stretches = []
    for kk in moved_s:
        a, b = max(1, kk - step), min(n_frames, kk + step)
        if stretches and a <= stretches[-1][1] + 1:
            stretches[-1][1] = max(stretches[-1][1], b)
        else:
            stretches.append([a, b])
    per = {}
    if stretches:
        cap = REAL_VC(str(video))
        k, last = 0, stretches[-1][1]
        while k < last:
            if not cap.grab():
                break
            k += 1
            if any(a <= k <= b for a, b in stretches):
                ok, fr = cap.retrieve()
                if ok:
                    per[k] = align(keypoints(small(fr)[0]), kds[kr])
        cap.release()
    return dict(fps=fps, n_frames=n_frames, step=step, ref_frame=kr, ref_gray=samples[kr], n_samples=len(ks),
                n_same=n_same, stretches=stretches, per=per, scale=scale, shape=shape)


def full_res_warps(v):
    """frame -> 2x3 warp in full-resolution pixels for the frames that differ from the reference view"""
    out, rows = {}, []
    for a, b in v["stretches"]:
        ks = list(range(a, b + 1))
        got = {k: v["per"].get(k, (None, 0)) for k in ks}
        good = [k for k in ks if got[k][0] is not None]
        for k in ks:
            A, n = got[k]
            interp = False
            if A is None and good:                      # cannot be aligned (blur while the camera moves): neighbours
                lo = max((g for g in good if g < k), default=None)
                hi = min((g for g in good if g > k), default=None)
                if lo is not None and hi is not None:
                    w = (k - lo) / (hi - lo)
                    A = (1 - w) * got[lo][0] + w * got[hi][0]
                else:
                    A = got[lo if lo is not None else hi][0]
                interp = True
            if A is None:
                continue
            mv = is_moved(A, v["shape"])
            sc, rot, dx, dy = params(A)
            rows.append(dict(frame=k, t_s=round((k - 1) / v["fps"], 3), moved=int(mv), scale=round(sc, 4),
                             rot_deg=round(rot, 3), dx_px=round(dx, 2), dy_px=round(dy, 2), inliers=n,
                             interpolated=int(interp)))
            if mv:
                Af = A.astype(np.float64).copy()
                Af[:, 2] /= v["scale"]                  # 640 px units -> full resolution
                out[k] = Af
    return out, pd.DataFrame(rows, columns=["frame", "t_s", "moved", "scale", "rot_deg", "dx_px", "dy_px", "inliers",
                                            "interpolated"])


class WarpCapture:
    """cv2.VideoCapture for the scorer: only the first MAX frames, the listed frames warped into the main view."""
    WARPS: dict = {}
    FILL = (140, 140, 140)
    MAX = 10 ** 9

    def __init__(self, path, *a):
        self._c = REAL_VC(path, *a)
        self._k = 0
        self._w = WarpCapture.WARPS.get(str(Path(str(path)).resolve()).lower(), {})

    def _warp(self, fr):
        A = self._w.get(self._k)
        if A is None or fr is None:
            return fr
        h, w = fr.shape[:2]
        return cv2.warpAffine(fr, A, (w, h), flags=cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT,
                              borderValue=WarpCapture.FILL)

    def grab(self):
        if self._k >= WarpCapture.MAX:
            return False
        ok = self._c.grab()
        if ok:
            self._k += 1
        return ok

    def retrieve(self, *a):
        ok, fr = self._c.retrieve(*a)
        return ok, (self._warp(fr) if ok else fr)

    def read(self, *a):
        if self._k >= WarpCapture.MAX:
            return False, None
        ok, fr = self._c.read(*a)
        if ok:
            self._k += 1
            fr = self._warp(fr)
        return ok, fr

    def get(self, prop):
        v = self._c.get(prop)
        return min(v, WarpCapture.MAX) if prop == cv2.CAP_PROP_FRAME_COUNT else v

    def set(self, prop, val):
        if prop in (cv2.CAP_PROP_POS_FRAMES, cv2.CAP_PROP_POS_MSEC, cv2.CAP_PROP_POS_AVI_RATIO):
            raise RuntimeError("tst_mouse_lock: seeking would break the frame count of the warps")
        return self._c.set(prop, val)

    def __getattr__(self, n):
        return getattr(self._c, n)


class _CV2(types.ModuleType):
    def __getattr__(self, n):
        return getattr(cv2, n)


STATS = {}


def hook(warps: dict, video: Path, fill, max_frames: int, max_s: float):
    WarpCapture.WARPS = {str(video.resolve()).lower(): warps}
    WarpCapture.FILL = fill
    WarpCapture.MAX = max_frames
    proxy = _CV2("cv2")
    proxy.VideoCapture = WarpCapture
    TF.cv2 = proxy
    orig = getattr(TF.video_stats, "__wrapped__", TF.video_stats)

    def vs(v):
        STATS[str(v)] = orig(v)
        return STATS[str(v)]
    vs.__wrapped__ = orig
    TF.video_stats = vs

    def num_frames(video_, n_decoded, fps_opencv):
        nf = ORIG_NUM_FRAMES(video_, n_decoded, fps_opencv)
        if n_decoded >= max_frames:                     # the video was cut at max_s: numFrames = frames of max_s
            fr = nf["frame_rate"]
            return dict(num_frames=int(math.floor(max_s * fr + 1e-6)), duration_s=float(max_s), frame_rate=fr,
                        source=f"first {max_s:g} s of the video (tst_mouse_lock; user 2026-11-18: TST = first 3 min)")
        return nf
    TAS.vm.manual_num_frames = num_frames


def unhook():
    TF.cv2 = cv2
    TF.video_stats = getattr(TF.video_stats, "__wrapped__", TF.video_stats)
    TAS.vm.manual_num_frames = ORIG_NUM_FRAMES
    WarpCapture.MAX = 10 ** 9


def score(video: Path, out: Path, cache: Path):
    out.mkdir(parents=True, exist_ok=True)
    return TAS._job((str(video), str(out), str(TAS.DEFAULT_MODEL), str(cache / f"{video.stem}.npz"), "manual"))


def mouse_qc(feat_npz: Path):
    z = np.load(feat_npz, allow_pickle=False)
    F, names = z["F"], [str(x) for x in z["names"]]
    col = lambda c: F[:, names.index(c)]  # noqa: E731
    fps, a_ref, l1 = float(z["fps"]), float(z["a_ref"]), float(z["l1"])
    W = TF.FEAT_WIDTH
    H = int(round(int(z["height"]) * W / int(z["width"])))
    valid = col("valid") > 0
    ratio = col("area") / max(a_ref, 1.0)
    rm = lambda x: pd.Series(np.where(valid, x, np.nan)).rolling(int(3 * fps), center=True,  # noqa: E731
                                                                 min_periods=10).median().to_numpy()
    jump = np.hypot(col("cx") - rm(col("cx")), col("cy") - rm(col("cy"))) / max(l1, 1.0)
    xf, yf = col("cx") / W, col("cy") / H
    with np.errstate(invalid="ignore"):
        bad_size = valid & ((ratio < AREA_LO) | (ratio > AREA_HI))
        bad_jump = valid & (jump > 1.0)
        bad_pos = valid & ((xf < X_LO) | (xf > X_HI) | (yf > Y_HI))
    bad = bad_size | bad_jump | bad_pos
    t = np.arange(len(F)) / fps
    Q = pd.DataFrame(dict(frame=np.arange(1, len(F) + 1), t_s=np.round(t, 3), valid=valid.astype(int),
                          area_ratio=np.round(ratio, 3), jump_l1=np.round(jump, 3), cx_frac=np.round(xf, 3),
                          cy_frac=np.round(yf, 3), not_mouse=bad.astype(int)))
    sec = pd.Series(bad).groupby(t.astype(int)).mean()
    return Q, dict(not_mouse_pct=round(100 * float(bad.mean()), 2), size_pct=round(100 * float(bad_size.mean()), 2),
                   jump_pct=round(100 * float(bad_jump.mean()), 2), position_pct=round(100 * float(bad_pos.mean()), 2),
                   seconds_mostly_not_mouse=[int(s) for s, v in sec.items() if v > 0.5], a_ref=round(a_ref, 1),
                   l1=round(l1, 2), frames=int(len(F)), seconds=round(len(F) / fps, 2))


def check_sheet(video: Path, warps: dict, stats_fixed: dict, stats_plain, out_png: Path, fps: float,
                plain_bad=None):
    """5 frames: where the package maps tracked something else (plain_bad, if known) or spread over the aligned
    frames, plus one frame after the last aligned one; original + package blob above, aligned + fixed blob below"""
    ks_all = sorted(warps)
    src = sorted(set(plain_bad) & set(ks_all)) if plain_bad is not None and len(plain_bad) else ks_all
    if not src:
        src = ks_all
    pick = sorted({src[int(i)] for i in np.linspace(0, len(src) - 1, min(4, len(src)))} |
                  {max(ks_all) + int(2 * fps)})
    got, cap, k = {}, REAL_VC(str(video)), 0
    while k < max(pick):
        if not cap.grab():
            break
        k += 1
        if k in pick:
            ok, fr = cap.retrieve()
            if ok:
                got[k] = fr
    cap.release()
    top, bot = [], []
    for k in pick:
        if k not in got:
            continue
        fr = got[k]
        A = warps.get(k)
        wf = fr if A is None else cv2.warpAffine(fr, A, (fr.shape[1], fr.shape[0]), borderMode=cv2.BORDER_CONSTANT,
                                                 borderValue=WarpCapture.FILL)
        for img, st, colr, dst in ((fr, stats_plain, (0, 0, 255), top), (wf, stats_fixed, (0, 200, 0), bot)):
            g = TF.gray_small(img)
            vis = cv2.resize(img, (int(g.shape[1] * 1.5), int(g.shape[0] * 1.5)))
            if st is not None:
                m, _ = TF.mouse_mask(g, st["thr"], st["static"], st["freq"])
                cs, _ = cv2.findContours(cv2.resize(m * 255, (vis.shape[1], vis.shape[0]),
                                                    interpolation=cv2.INTER_NEAREST), cv2.RETR_EXTERNAL,
                                         cv2.CHAIN_APPROX_NONE)
                cv2.drawContours(vis, cs, -1, colr, 2)
            cv2.putText(vis, f"{(k - 1) / fps:.1f} s", (6, 22), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 0, 255), 2)
            dst.append(vis)
    if top:
        sheet = np.vstack([np.hstack(top), np.hstack(bot)])
        bar = np.full((30, sheet.shape[1], 3), 255, np.uint8)
        cv2.putText(bar, "top: original frames, blob picked with the package maps (red) | bottom: frames aligned to "
                         "the main view, blob after the fix (green)", (6, 21), cv2.FONT_HERSHEY_SIMPLEX, 0.5,
                    (0, 0, 0), 1)
        sheet = np.vstack([sheet, bar])
        cv2.imwrite(str(out_png), sheet)


def qc_ok(q: dict) -> bool:
    return q["not_mouse_pct"] <= FIX_IF_NOT_MOUSE_PCT and not q["seconds_mostly_not_mouse"]


def run_scored(video, out_sub, cache, warps, fill, max_frames, max_s):
    hook(warps, video, fill, max_frames, max_s)
    try:
        row = score(video, out_sub, cache)
        stats = STATS.get(str(video))
    finally:
        unhook()
    Q, q = mouse_qc(cache / f"{video.stem}.npz")
    return row, stats, Q, q


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--video", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-s", type=float, default=180.0)
    ap.add_argument("--force-fix", action="store_true", help="try the view fix even when the untouched run passes")
    ap.add_argument("--no-fix", action="store_true", help="never try the view fix")
    ap.add_argument("--compare-plain", action="store_true", help="(kept for old calls; the untouched run is always made)")
    a = ap.parse_args(argv)
    idle()
    cv2.setNumThreads(1)
    video, out = Path(a.video), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    stem, tag = video.stem, f"{a.max_s:g}s"
    cap = REAL_VC(str(video))
    fps0 = cap.get(cv2.CAP_PROP_FPS) or 30.0
    ok, fr0 = cap.read()
    cap.release()
    max_frames = int(math.floor(a.max_s * fps0 + 1e-6))
    fill = (int(np.median(cv2.cvtColor(fr0, cv2.COLOR_BGR2GRAY))),) * 3 if ok else (140, 140, 140)
    # 1) the package scorer on the untouched first max_s seconds
    prow, stats_p, Qp, qp = run_scored(video, out / f"untouched_{tag}", out / f"cache_untouched_{tag}", {}, fill,
                                       max_frames, a.max_s)
    info = dict(video=str(video), max_s=a.max_s, fps=fps0,
                untouched=dict(pct_nonmoving=prow.get("pct_nonmoving"), qc_status=prow.get("qc_status"),
                               qc_messages=prow.get("qc_messages"), numFrames=prow.get("numFrames"),
                               numFrames_source=prow.get("numFrames_source"), mouse_qc=qp),
                mouse_qc_rule=f"not mouse = area outside {AREA_LO}-{AREA_HI} x a_ref, or centre > 1 l1 from its 3-s "
                              f"median, or centre outside x {X_LO:.0%}-{X_HI:.0%} / top {Y_HI:.0%} of the frame; the "
                              f"view fix is tried when > {FIX_IF_NOT_MOUSE_PCT} % of frames are not the mouse or any "
                              f"second is mostly not the mouse")
    chosen, Q, why = "untouched", Qp, "the tracked blob is the mouse in the untouched run"
    need = (not qc_ok(qp) or a.force_fix) and not a.no_fix
    if need:
        # 2) camera view fix, only for videos whose untouched run lost the mouse
        v = scan_views(video, max_frames)
        warps, view_df = full_res_warps(v)
        view_df.to_csv(out / f"{stem}_view.csv", index=False)
        info["view"] = dict(reference_frame=v["ref_frame"], samples=v["n_samples"], samples_in_main_view=v["n_same"],
                            moved_stretches_s=[[round((x - 1) / v["fps"], 2), round((y - 1) / v["fps"], 2)]
                                               for x, y in v["stretches"]],
                            warped_frames=len(warps), interpolated_frames=int(view_df["interpolated"].sum()),
                            move_threshold="corner shift > 2 px at the scorer's 320 px width")
        print(f"{video.name}: untouched run lost the mouse in {qp['not_mouse_pct']} % of frames -> view fix: "
              f"main view in {v['n_same']}/{v['n_samples']} samples, moved {info['view']['moved_stretches_s']} s, "
              f"{len(warps)} frames aligned", flush=True)
        if warps:
            frow, stats_f, Qf, qf = run_scored(video, out / f"viewfix_{tag}", out / f"cache_viewfix_{tag}", warps, fill,
                                               max_frames, a.max_s)
            info["viewfix"] = dict(pct_nonmoving=frow.get("pct_nonmoving"), qc_status=frow.get("qc_status"),
                                   qc_messages=frow.get("qc_messages"), numFrames=frow.get("numFrames"), mouse_qc=qf)
            check_sheet(video, warps, stats_f, stats_p, out / f"{stem}_lock_check.png", v["fps"],
                        Qp.frame[Qp.not_mouse == 1].to_numpy())
            if qf["not_mouse_pct"] < qp["not_mouse_pct"] and qc_ok(qf):
                chosen, Q, why = "viewfix", Qf, "the untouched run lost the mouse; after aligning the moved frames it " \
                                                 "is the mouse again"
            else:
                why = "NEEDS A LOOK: the view fix did not make the tracked blob the mouse - untouched result kept"
        else:
            why = "NEEDS A LOOK: the untouched run lost the mouse but no camera move was found - untouched result kept"
    src = out / f"{chosen}_{tag}"
    for f in (f"AutoScoringResults_{stem}.mat", f"AutoScoringResults_{stem}_summary.csv"):
        if (src / f).exists():
            (out / f).write_bytes((src / f).read_bytes())
    Q.to_csv(out / f"{stem}_mouse_qc.csv", index=False)
    final = info[chosen]
    info.update(chosen=chosen, why=why, pct_nonmoving=final["pct_nonmoving"], mouse_qc=final["mouse_qc"])
    (out / f"{stem}_lock.json").write_text(json.dumps(info, indent=1, default=float), encoding="utf-8")
    msg = (f"  {video.name}: first {a.max_s:g} s {final['pct_nonmoving']:.1f} % non-moving [{chosen}] - not-mouse "
           f"frames {final['mouse_qc']['not_mouse_pct']} %")
    if chosen != "untouched" or need:
        msg += f" (untouched: {info['untouched']['pct_nonmoving']:.1f} %, not-mouse {qp['not_mouse_pct']} %)"
    print(msg + f"; {why}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
