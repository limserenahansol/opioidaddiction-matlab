#!/usr/bin/env python3
"""
tst_examples.py - example frames of the automatic TST detection, for the user to confirm (added 2026-11-18; user: "now
8 mice TST is done can you do auto and give me examples of detection for me to confirm").

For every video scored by tst_mouse_lock.py (<lock-root>/<video>/ with <video>_lock.json and the chosen
AutoScoringResults_<video>.mat), one sheet:
  * 6 example frames: the middle of the 3 longest MOVING and the 3 longest NON-MOVING stretches of the prediction
    (behaviour time, 'pred_frame'), first 3 min only, in time order;
  * on each frame: green = the mouse found in that frame, magenta = the mouse 0.5 s earlier (if the two outlines
    differ, the mouse moved), white box = the box the scorer tracked (its own per-frame features); label = the call
    (MOVING / NON-MOVING) and the smoothed P(moving);
  * a strip of the whole first 3 min: orange = moving, blue = non-moving, black ticks = the example frames;
  * the header: % non-moving, which run was used (untouched / camera-view fix) and the mouse QC.
Frames of a 'viewfix' run are shown aligned to the main view, as the scorer saw them. The outlines are drawn again
with the scorer's own mouse_mask (without its 'previous frame' term), so in rare frames they can differ slightly from
the tracked blob; the white box is the tracked blob itself.
Outputs: <out>/<video>_examples.png and <out>/TST_examples_overview.png (all videos, smaller).
usage: python tst_examples.py --lock-root analysis/tst_test/lock2 --out analysis/tst_test/examples_day5
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import scipy.io as sio

T = Path(__file__).resolve().parent
sys.path.insert(0, str(T))
import tst_mouse_lock as TML  # noqa: E402  (also puts the package tst folder on the path)

TF = TML.TF
TILE_W, TILE_H = 480, 270
ORANGE, BLUE, GREEN, MAGENTA, WHITE = (0, 140, 255), (200, 90, 20), (0, 210, 0), (230, 0, 230), (255, 255, 255)


def runs(y, val, min_len):
    """(start, end) 0-based frame runs where y == val, at least min_len frames"""
    out, s = [], None
    for j, v in enumerate(np.append(y, -1)):
        if v == val and s is None:
            s = j
        elif v != val and s is not None:
            if j - s >= min_len:
                out.append((s, j - 1))
            s = None
    return out


def warps_from_view(folder: Path, stem: str, width: int):
    f = folder / f"{stem}_view.csv"
    if not f.exists():
        return {}
    d = pd.read_csv(f)
    d = d[d.moved == 1]
    s = TML.ORB_W / width
    out = {}
    for r in d.itertuples():
        th = math.radians(r.rot_deg)
        a, b = r.scale * math.cos(th), r.scale * math.sin(th)
        out[int(r.frame)] = np.array([[a, -b, r.dx_px / s], [b, a, r.dy_px / s]], np.float64)
    return out


def one(folder: Path, out_dir: Path):
    stem = folder.name
    info = json.loads((folder / f"{stem}_lock.json").read_text(encoding="utf-8"))
    video = Path(info["video"])
    chosen, max_s = info["chosen"], float(info["max_s"])
    M = sio.loadmat(folder / f"AutoScoringResults_{stem}.mat", squeeze_me=True)
    y = np.asarray(M["pred_frame"]).astype(int).ravel()
    p = np.asarray(M["prob_moving_smoothed"]).astype(float).ravel()
    fps = float(info["fps"])
    z = np.load(folder / f"cache_{chosen}_{max_s:g}s" / f"{stem}.npz", allow_pickle=False)
    F, names = z["F"], [str(x) for x in z["names"]]
    col = lambda c: F[:, names.index(c)]  # noqa: E731
    W = int(z["width"])
    warps = warps_from_view(folder, stem, W) if chosen == "viewfix" else {}
    max_frames = int(math.floor(max_s * fps + 1e-6))
    TML.hook(warps, video, (140, 140, 140), max_frames, max_s)
    try:
        st = TF.video_stats(str(video))
    finally:
        TML.unhook()
    # examples: middle of the 3 longest moving and 3 longest non-moving runs (>= 0.5 s)
    pick = []
    for val in (1, 0):
        r = sorted(runs(y, val, int(0.5 * fps)), key=lambda x: x[0] - x[1])[:3]
        pick += [((a + b) // 2, val) for a, b in r]
    pick.sort()
    need = sorted({j for j, _ in pick} | {max(0, j - int(0.5 * fps)) for j, _ in pick})
    frames, cap, j = {}, TML.REAL_VC(str(video)), -1
    while j < max(need):
        ok, fr = cap.read()
        if not ok:
            break
        j += 1
        if j in need:
            A = warps.get(j + 1)
            if A is not None:
                fr = cv2.warpAffine(fr, A, (fr.shape[1], fr.shape[0]), borderMode=cv2.BORDER_CONSTANT,
                                    borderValue=(140, 140, 140))
            frames[j] = fr
    cap.release()
    # one crop per video around where the tracked mouse goes (1st-99th percentile of its boxes + margin, 16:9), so the
    # mouse is large enough to judge; the same crop for all 6 frames
    fr0 = next(iter(frames.values()))
    H0, W0 = fr0.shape[:2]
    k = W0 / TF.FEAT_WIDTH                                   # 320 px feature units -> full-resolution pixels
    ok_rows = col("valid")[:max_frames] > 0
    bl, bt = col("bb_left")[:max_frames][ok_rows], col("bb_top")[:max_frames][ok_rows]
    br, bb = bl + col("bb_w")[:max_frames][ok_rows], bt + col("bb_h")[:max_frames][ok_rows]
    cx0, cx1 = np.percentile(bl, 1) * k, np.percentile(br, 99) * k
    cy0, cy1 = np.percentile(bt, 1) * k, np.percentile(bb, 99) * k
    mg = 0.35 * max(cx1 - cx0, cy1 - cy0)
    cx0, cx1, cy0, cy1 = cx0 - mg, cx1 + mg, cy0 - mg, cy1 + mg
    w, h = cx1 - cx0, cy1 - cy0
    if w / h < 16 / 9:
        cx0, cx1 = (cx0 + cx1) / 2 - h * 8 / 9, (cx0 + cx1) / 2 + h * 8 / 9
    else:
        cy0, cy1 = (cy0 + cy1) / 2 - w * 9 / 32, (cy0 + cy1) / 2 + w * 9 / 32
    w, h = min(W0, cx1 - cx0), min(H0, cy1 - cy0)
    cx0 = int(round(min(max(0, cx0), W0 - w)))
    cy0 = int(round(min(max(0, cy0), H0 - h)))
    cw, ch = int(round(w)), int(round(h))
    lw = max(1, int(round(cw / TILE_W)))
    tiles = []
    for j, val in pick:
        if j not in frames:
            continue
        vis = frames[j].copy()
        for jj, colr, th in ((max(0, j - int(0.5 * fps)), MAGENTA, 1), (j, GREEN, 2)):
            if jj not in frames:
                continue
            m, _ = TF.mouse_mask(TF.gray_small(frames[jj]), st["thr"], st["static"], st["freq"])
            cs, _ = cv2.findContours(cv2.resize(m * 255, (W0, H0), interpolation=cv2.INTER_NEAREST),
                                     cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            cv2.drawContours(vis, cs, -1, colr, th * lw)
        if j < len(F):
            x0, y0 = col("bb_left")[j] * k, col("bb_top")[j] * k
            cv2.rectangle(vis, (int(x0), int(y0)), (int(x0 + col("bb_w")[j] * k), int(y0 + col("bb_h")[j] * k)), WHITE,
                          lw)
        vis = cv2.resize(vis[cy0:cy0 + ch, cx0:cx0 + cw], (TILE_W, TILE_H), interpolation=cv2.INTER_AREA)
        lab = "MOVING" if val == 1 else "NON-MOVING"
        cv2.rectangle(vis, (0, 0), (TILE_W, 30), ORANGE if val == 1 else BLUE, -1)
        cv2.putText(vis, f"{lab}  P(moving) {p[j]:.2f}   t = {j / fps:.1f} s", (8, 21), cv2.FONT_HERSHEY_SIMPLEX,
                    0.6, WHITE, 2)
        tiles.append(vis)
    while len(tiles) < 6:
        tiles.append(np.full((TILE_H, TILE_W, 3), 230, np.uint8))
    grid = np.vstack([np.hstack(tiles[:3]), np.hstack(tiles[3:6])])
    # strip of the first 3 min
    n = min(len(y), max_frames)
    strip = np.zeros((40, grid.shape[1], 3), np.uint8)
    xs = (np.arange(n) * grid.shape[1] / max(1, n)).astype(int)
    for xx, v in zip(xs, y[:n]):
        strip[:, xx] = ORANGE if v == 1 else BLUE
    for j, _ in pick:
        xx = int(j * grid.shape[1] / max(1, n))
        cv2.line(strip, (xx, 0), (xx, 39), (0, 0, 0), 3)
    scale = np.full((22, grid.shape[1], 3), 255, np.uint8)
    for s in range(0, int(max_s) + 1, 30):
        xx = int(s * fps * grid.shape[1] / max(1, n))
        cv2.putText(scale, f"{s}s", (min(xx, grid.shape[1] - 40), 16), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 0), 1)
    q = info["mouse_qc"]
    head = np.full((70, grid.shape[1], 3), 255, np.uint8)
    title = (f"{stem}: first {max_s:g} s  {info['pct_nonmoving']:.1f} % non-moving   (run: "
             f"{'camera-view fix' if chosen == 'viewfix' else 'untouched'}; mouse QC: {q['not_mouse_pct']} % frames "
             f"not the mouse)")
    cv2.putText(head, title, (8, 26), cv2.FONT_HERSHEY_SIMPLEX, 0.65, (0, 0, 0), 2)
    cv2.putText(head, "green = mouse now, magenta = mouse 0.5 s before (different outline = it moved), white box = "
                      "tracked box; strip: orange = moving, blue = non-moving, black = examples",
                (8, 56), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (60, 60, 60), 1)
    sheet = np.vstack([head, grid, strip, scale])
    out_dir.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(out_dir / f"{stem}_examples.png"), sheet)
    return sheet, dict(video=stem, pct_nonmoving=round(info["pct_nonmoving"], 1), run=chosen,
                       not_mouse_pct=q["not_mouse_pct"], examples_s=[round(j / fps, 1) for j, _ in pick],
                       examples_call=["moving" if v else "non-moving" for _, v in pick])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lock-root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", nargs="*", default=None, help="video names (folder names) to do; default all")
    a = ap.parse_args(argv)
    TML.idle()
    cv2.setNumThreads(1)
    out = Path(a.out)
    sheets, rows = [], []
    for folder in sorted(p for p in Path(a.lock_root).iterdir() if p.is_dir() and (p / f"{p.name}_lock.json").exists()
                         and (not a.only or p.name in a.only)):
        sheet, row = one(folder, out)
        sheets.append(sheet)
        rows.append(row)
        print(f"  {row['video']}: {row['pct_nonmoving']} % non-moving [{row['run']}], examples at {row['examples_s']} s")
    if sheets:
        w = max(s.shape[1] for s in sheets)
        small = [cv2.resize(s, (w // 2, int(s.shape[0] * (w // 2) / s.shape[1]))) for s in sheets]
        cols = 2
        rows_img = []
        for i in range(0, len(small), cols):
            part = small[i:i + cols]
            h = max(x.shape[0] for x in part)
            part = [np.vstack([x, np.full((h - x.shape[0], x.shape[1], 3), 255, np.uint8)]) for x in part]
            while len(part) < cols:
                part.append(np.full((h, w // 2, 3), 255, np.uint8))
            rows_img.append(np.hstack(part))
        cv2.imwrite(str(out / "TST_examples_overview.png"), np.vstack(rows_img))
        pd.DataFrame(rows).to_csv(out / "TST_examples_index.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
