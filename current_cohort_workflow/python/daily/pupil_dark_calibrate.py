#!/usr/bin/env python3
"""
pupil_dark_calibrate.py - how big the dark-contrast pupil outline (pupil_dark_contrast.py) is compared with the user's
own pupil drawings (added 2026-11-18). The outline level (FRAC, between the blob's dark core and the ring around it)
is chosen so that the new diameter matches the user's drawings; the package's raw detection is shown alongside.

Drawings: analysis/pupil_teach/marks.json (choice 'draw'; the label belongs to image_frame, a 1-based video frame;
ellipse cx, cy, w, h in crop pixels, 0-based pixel centres). by = None: drawn by the user on the drawing page (exact);
by = 'claude-read-from-screenshot': read from a screenshot (position reliable, size approximate) - reported apart.
User diameter = sqrt(w*h) (equivalent diameter of the ellipse, as the package's EquivDiameter).
Per drawing: the frame is decoded exactly (pupil_example_sheets.grab), the dark-contrast image made (prep), the round
dark blob nearest to the user's centre taken (within 12 px), and outlined at FRAC 0.35 ... 0.8.
Package: raw detection of that frame (hybrid_detections.csv, the U-Net blob before tracking), counted only where its
centre lies within max(4, 0.35 d_user) px of the user's centre (= the same blob).
Output: analysis/pupil_dark/calibration/ calibration_frames.csv, calibration_summary.csv, calibration_sheet.png
usage: python pupil_dark_calibrate.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import cv2
import numpy as np
import pandas as pd

T = Path(__file__).resolve().parent
P = T.parent
sys.path.insert(0, str(T))
import pupil_dark_contrast as D  # noqa: E402
import pupil_example_sheets as PES  # noqa: E402

FRACS = (0.35, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8)
OUT = P / "analysis" / "pupil_dark" / "calibration"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    marks = json.loads((P / "analysis" / "pupil_teach" / "marks.json").read_text(encoding="utf-8"))
    root = D.pupil_root()
    rows, tiles = [], []
    det_cache = {}
    for key, m in sorted(marks.items()):
        if m.get("choice") != "draw" or not m.get("ellipse"):
            continue
        day, cage, mouse = int(m["day"]), m["cage"], m["mouse"]
        k = int(m.get("image_frame") or m["frame"])
        try:
            video, crop = D.video_of(day, cage, mouse), D.crop_of(day, cage, mouse)
        except SystemExit as exc:
            print(f"  skip {key}: {exc}")
            continue
        if m.get("crop") and (int(m["crop"]["w"]), int(m["crop"]["h"])) != (crop[2], crop[3]):
            print(f"  skip {key}: eye box changed since the drawing ({m['crop']} vs {crop[2:]})")
            continue
        C = PES.grab(video, crop, k)
        if C is None:
            print(f"  skip {key}: frame {k} not decoded")
            continue
        e = m["ellipse"]
        ux, uy, ud = float(e["cx"]), float(e["cy"]), float(np.sqrt(e["w"] * e["h"]))   # 0-based
        s, gl = D.prep(C)
        pk = D.peaks(s, (ux, uy), search=12)
        best = min(pk, key=lambda p: np.hypot(p[0] - ux, p[1] - uy)) if pk else None
        r = dict(key=key, session=m["session"], day=day, cage=cage, mouse=mouse, frame=k,
                 by="page" if not m.get("by") else "screenshot", user_d=ud, user_cx=ux + 1, user_cy=uy + 1,
                 glint=int(gl is not None))
        segs = {}
        if best is not None:
            r.update(peak_dist=float(np.hypot(best[0] - ux, best[1] - uy)), peak_r=best[2], peak_score=best[3])
            for f in FRACS:
                sg = D.segment(s, best[0], best[1], best[2], frac=f)
                segs[f] = sg
                r[f"new_d_{f}"] = sg["diam"] if sg else np.nan
                r[f"new_dist_{f}"] = float(np.hypot(sg["cx"] - 1 - ux, sg["cy"] - 1 - uy)) if sg else np.nan
        sd = root / f"day{day}" / cage / mouse / "pupil_outputs_general" / "hybrid_detections.csv"
        if sd not in det_cache:
            det_cache[sd] = pd.read_csv(sd) if sd.exists() else None
        H = det_cache[sd]
        if H is not None:
            h = H[H.frame == k]
            if len(h):
                h = h.iloc[0]
                r.update(pkg_d=float(h.best_diam_px), pkg_dist=float(np.hypot(h.best_cx_px - 1 - ux,
                                                                             h.best_cy_px - 1 - uy)))
        rows.append(r)
        if len(tiles) < 24 and best is not None:
            lo, hi = np.percentile(C, [2, 98])
            v = np.clip((C.astype(np.float32) - lo) / max(1.0, hi - lo) * 255, 0, 255).astype(np.uint8)
            v = cv2.cvtColor(cv2.resize(v, (v.shape[1] * 2, v.shape[0] * 2), interpolation=cv2.INTER_CUBIC),
                             cv2.COLOR_GRAY2BGR)
            cv2.ellipse(v, (int(round(ux * 2)), int(round(uy * 2))), (int(round(e["w"])), int(round(e["h"]))), 0, 0,
                        360, (0, 140, 255), 2, lineType=cv2.LINE_AA)
            for f, col in ((0.5, (255, 220, 0)), (0.65, (60, 255, 60))):
                sg = segs.get(f)
                if sg:
                    M, x0, y0 = sg["mask"]
                    big = np.zeros(C.shape, np.uint8)
                    big[y0:y0 + M.shape[0], x0:x0 + M.shape[1]] = M
                    cnt, _ = cv2.findContours(cv2.resize(big * 255, (v.shape[1], v.shape[0]),
                                                         interpolation=cv2.INTER_NEAREST),
                                              cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
                    cv2.drawContours(v, cnt, -1, col, 1)
            cv2.putText(v, f"{m['session']} f{k} ({r['by']})", (5, 18), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 255), 2)
            tiles.append(cv2.resize(v, (360, int(360 * v.shape[0] / v.shape[1]))))
    F = pd.DataFrame(rows)
    F.to_csv(OUT / "calibration_frames.csv", index=False, float_format="%.4g")
    summ = []
    for by in ("page", "screenshot", "all"):
        G = F if by == "all" else F[F.by == by]
        G = G[G.peak_dist <= 12] if "peak_dist" in G else G
        for f in FRACS:
            ratio = (G[f"new_d_{f}"] / G.user_d).dropna()
            summ.append(dict(drawings=by, what=f"new FRAC {f}", n=len(ratio), median_ratio_to_user=ratio.median(),
                             q25=ratio.quantile(0.25), q75=ratio.quantile(0.75),
                             median_centre_err_px=G[f"new_dist_{f}"].median()))
        same = G[G.pkg_dist <= np.maximum(4, 0.35 * G.user_d)] if "pkg_dist" in G else G.iloc[0:0]
        ratio = (same.pkg_d / same.user_d).dropna()
        summ.append(dict(drawings=by, what="package raw (same blob only)", n=len(ratio),
                         median_ratio_to_user=ratio.median(), q25=ratio.quantile(0.25), q75=ratio.quantile(0.75),
                         median_centre_err_px=same.pkg_dist.median()))
    S = pd.DataFrame(summ)
    S.to_csv(OUT / "calibration_summary.csv", index=False, float_format="%.3g")
    print(S.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    if tiles:
        h = max(t.shape[0] for t in tiles)
        tiles = [np.vstack([t, np.full((h - t.shape[0], t.shape[1], 3), 255, np.uint8)]) for t in tiles]
        while len(tiles) % 6:
            tiles.append(np.full_like(tiles[0], 255))
        img = np.vstack([np.hstack(tiles[i:i + 6]) for i in range(0, len(tiles), 6)])
        bar = np.full((30, img.shape[1], 3), 255, np.uint8)
        cv2.putText(bar, "orange = your drawing,  cyan = new outline at FRAC 0.5,  green = FRAC 0.65", (6, 21),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 0, 0), 1)
        cv2.imwrite(str(OUT / "calibration_sheet.png"), np.vstack([bar, img]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
