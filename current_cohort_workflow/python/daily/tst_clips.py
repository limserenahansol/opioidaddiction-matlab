#!/usr/bin/env python3
"""
tst_clips.py - a few short annotated clips of the automatic TST detection (added 2026-11-18; user: "i wanna see short
video your non moving and moving just show me a few short video (very a few)").

Each clip: --len seconds of the video around a chosen time, real speed (30 fps), cropped around where the mouse goes
(the same crop as tst_examples.py), with in every frame:
  top bar   MOVING (orange) / NON-MOVING (blue) = the scorer's call for THIS frame (behaviour time, 'pred_frame'),
            the smoothed P(moving) and the time in the video;
  green     the mouse found in this frame (the scorer's mouse_mask, following the previous frame as the scorer does);
  bottom    where the clip is inside the first 3 min (orange / blue = the calls).
Frames of a 'viewfix' run are shown aligned to the main view, as the scorer saw them. Encoded with the ffmpeg of
tools/.venv (imageio-ffmpeg; H.264, yuv420p) so that it plays everywhere. Idle priority.
usage: python tst_clips.py --lock-root analysis/tst_test/lock2 --out analysis/tst_test/clips_day5
         --clip 9669mwhite_TST_day5:16.1 --clip 9669mwhite_TST_day5:95.5 [--len 5]
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import subprocess
import sys
from pathlib import Path

import cv2
import numpy as np
import scipy.io as sio

T = Path(__file__).resolve().parent
sys.path.insert(0, str(T))
import tst_mouse_lock as TML  # noqa: E402
import tst_examples as TE  # noqa: E402

TF = TML.TF
OUT_W, OUT_H = 640, 360
FFMPEG = sorted(glob.glob(str(T / ".venv" / "Lib" / "site-packages" / "imageio_ffmpeg" / "binaries" / "ffmpeg*.exe")))


def crop_box(F, names, max_frames, W0, H0):
    col = lambda c: F[:, names.index(c)]  # noqa: E731
    k = W0 / TF.FEAT_WIDTH
    ok = col("valid")[:max_frames] > 0
    bl, bt = col("bb_left")[:max_frames][ok], col("bb_top")[:max_frames][ok]
    br, bb = bl + col("bb_w")[:max_frames][ok], bt + col("bb_h")[:max_frames][ok]
    x0, x1 = np.percentile(bl, 1) * k, np.percentile(br, 99) * k
    y0, y1 = np.percentile(bt, 1) * k, np.percentile(bb, 99) * k
    mg = 0.35 * max(x1 - x0, y1 - y0)
    x0, x1, y0, y1 = x0 - mg, x1 + mg, y0 - mg, y1 + mg
    w, h = x1 - x0, y1 - y0
    if w / h < 16 / 9:
        x0, x1 = (x0 + x1) / 2 - h * 8 / 9, (x0 + x1) / 2 + h * 8 / 9
    else:
        y0, y1 = (y0 + y1) / 2 - w * 9 / 32, (y0 + y1) / 2 + w * 9 / 32
    w, h = min(W0, x1 - x0), min(H0, y1 - y0)
    x0 = int(round(min(max(0, x0), W0 - w)))
    y0 = int(round(min(max(0, y0), H0 - h)))
    return x0, y0, int(round(w)), int(round(h))


def one_clip(folder: Path, t_mid: float, length: float, out_mp4: Path, cache: dict):
    stem = folder.name
    info = json.loads((folder / f"{stem}_lock.json").read_text(encoding="utf-8"))
    video, chosen, max_s, fps = Path(info["video"]), info["chosen"], float(info["max_s"]), float(info["fps"])
    M = sio.loadmat(folder / f"AutoScoringResults_{stem}.mat", squeeze_me=True)
    y = np.asarray(M["pred_frame"]).astype(int).ravel()
    p = np.asarray(M["prob_moving_smoothed"]).astype(float).ravel()
    z = np.load(folder / f"cache_{chosen}_{max_s:g}s" / f"{stem}.npz", allow_pickle=False)
    F, names = z["F"], [str(x) for x in z["names"]]
    W0, H0 = int(z["width"]), int(z["height"])
    max_frames = int(math.floor(max_s * fps + 1e-6))
    warps = TE.warps_from_view(folder, stem, W0) if chosen == "viewfix" else {}
    if stem not in cache:
        TML.hook(warps, video, (140, 140, 140), max_frames, max_s)
        try:
            cache[stem] = TF.video_stats(str(video))
        finally:
            TML.unhook()
    st = cache[stem]
    cx, cy, cw, ch = crop_box(F, names, max_frames, W0, H0)
    j0 = max(0, int(round((t_mid - length / 2) * fps)))
    j1 = min(len(y) - 1, j0 + int(round(length * fps)) - 1)
    proc = subprocess.Popen([FFMPEG[-1], "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "bgr24", "-s",
                             f"{OUT_W}x{OUT_H}", "-r", f"{fps:g}", "-i", "-", "-c:v", "libx264", "-pix_fmt", "yuv420p",
                             "-crf", "20", "-preset", "veryfast", str(out_mp4)], stdin=subprocess.PIPE,
                            creationflags=0x40 if sys.platform == "win32" else 0)
    cap = TML.REAL_VC(str(video))
    j, prev = -1, None
    k_open = np.ones((3, 3), np.uint8)
    n_mov = 0
    while j < j1:
        ok, fr = cap.read()
        if not ok:
            break
        j += 1
        if j < j0 - int(fps):                      # warm-up second: the mask follows the previous frame
            continue
        A = warps.get(j + 1)
        if A is not None:
            fr = cv2.warpAffine(fr, A, (W0, H0), borderMode=cv2.BORDER_CONSTANT, borderValue=(140, 140, 140))
        m, _ = TF.mouse_mask(TF.gray_small(fr), st["thr"], st["static"], st["freq"], prev, k_open)
        prev = m
        if j < j0:
            continue
        vis = fr.copy()
        cs, _ = cv2.findContours(cv2.resize(m * 255, (W0, H0), interpolation=cv2.INTER_NEAREST), cv2.RETR_EXTERNAL,
                                 cv2.CHAIN_APPROX_NONE)
        cv2.drawContours(vis, cs, -1, (0, 210, 0), max(2, int(round(cw / OUT_W * 2))))
        vis = cv2.resize(vis[cy:cy + ch, cx:cx + cw], (OUT_W, OUT_H), interpolation=cv2.INTER_AREA)
        mv = int(y[j]) == 1
        n_mov += mv
        cv2.rectangle(vis, (0, 0), (OUT_W, 34), (0, 140, 255) if mv else (200, 90, 20), -1)
        cv2.putText(vis, f"{'MOVING' if mv else 'NON-MOVING'}   P(moving) {p[j]:.2f}", (10, 24),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.75, (255, 255, 255), 2)
        cv2.putText(vis, f"t = {j / fps:5.1f} s", (OUT_W - 150, 24), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 255, 255), 2)
        cv2.putText(vis, stem.replace("_TST_day5", " TST day 5"), (10, OUT_H - 34), cv2.FONT_HERSHEY_SIMPLEX, 0.55,
                    (0, 0, 0), 2)
        bar = vis[OUT_H - 16:OUT_H]
        n = min(len(y), max_frames)
        xs = (np.arange(n) * OUT_W / n).astype(int)
        bar[:] = 255
        for xx, vv in zip(xs, y[:n]):
            bar[:, xx] = (0, 140, 255) if vv == 1 else (200, 90, 20)
        xx = int(j * OUT_W / n)
        cv2.rectangle(vis, (max(0, xx - 2), OUT_H - 18), (min(OUT_W - 1, xx + 2), OUT_H - 1), (0, 0, 0), -1)
        proc.stdin.write(vis.tobytes())
    cap.release()
    proc.stdin.close()
    proc.wait()
    return j1 - j0 + 1, n_mov


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lock-root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--clip", action="append", required=True, help="<video folder name>:<centre time s>")
    ap.add_argument("--len", type=float, default=5.0)
    a = ap.parse_args(argv)
    if not FFMPEG:
        raise SystemExit("ffmpeg of tools/.venv (imageio-ffmpeg) not found")
    TML.idle()
    cv2.setNumThreads(1)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cache = {}
    for c in a.clip:
        name, t = c.rsplit(":", 1)
        mp4 = out / f"{name}_{float(t):.0f}s.mp4"
        n, n_mov = one_clip(Path(a.lock_root) / name, float(t), a.len, mp4, cache)
        print(f"  {mp4.name}: {n} frames, called moving in {100 * n_mov / max(1, n):.0f} %")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
