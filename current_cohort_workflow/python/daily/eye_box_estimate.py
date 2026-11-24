#!/usr/bin/env python3
"""
eye_box_estimate.py - proposes the day's eye-box centres for eye_box_fixed.py (added 2026-11-23; the method Claude used
by hand on day 9: the darkest blob of eye size nearest the previous day's centre, median over 9 frames from 10 % to
90 % of the video). Report only: writes a check sheet with the mouse's FIXED box (video/eye_box_sizes.csv) drawn at the
proposal on 3 frames, and prints the eye_box_fixed.py command. Nothing is set until that command is run with --apply.
usage: python eye_box_estimate.py --day N [--prev M]   (M default N-1; the previous centres come from video/roi_log.csv)
Output: analysis/eye_box_day<N>/day<N>_eye_box_proposal.jpg + .csv
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import cv2
import numpy as np

T = Path(__file__).resolve().parent
P = T.parent
MICE = [(c, m) for c in ("9669m", "9670f") for m in ("black", "green1", "green2", "white")]


def prev_centres(day: int) -> dict:
    out = {}
    with open(P / "video" / "roi_log.csv", newline="", encoding="utf-8-sig") as fh:
        for r in csv.reader(fh):
            if len(r) >= 6 and r[1] == f"day{day}":
                m = re.search(r"centred on the eye opening \((\d+), (\d+)\)", r[5])
                if m:
                    out[(r[2], r[3])] = (int(m.group(1)), int(m.group(2)))
    return out


def sizes() -> dict:
    with open(P / "video" / "eye_box_sizes.csv", newline="", encoding="utf-8-sig") as fh:
        return {(r["cage"], r["mouse"]): (int(r["w"]), int(r["h"])) for r in csv.DictReader(fh)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--day", type=int, required=True)
    ap.add_argument("--prev", type=int)
    a = ap.parse_args(argv)
    prev = prev_centres(a.prev or a.day - 1)
    sz = sizes()
    out = P / "analysis" / f"eye_box_day{a.day}"
    out.mkdir(parents=True, exist_ok=True)
    rows, tiles, sets = [], [], []
    for c, m in MICE:
        vids = sorted((P / "video" / f"day{a.day}" / c / m).glob("*.mp4"))
        if not vids or (c, m) not in prev:
            print(f"  {c} {m}: no video or no previous centre - skipped")
            continue
        cx, cy = prev[(c, m)]
        cap = cv2.VideoCapture(str(vids[0]))
        n = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
        pts = []
        for fr in np.linspace(n * 0.1, n * 0.9, 9).astype(int):
            cap.set(cv2.CAP_PROP_POS_FRAMES, int(fr))
            ok, f = cap.read()
            if not ok:
                continue
            g = cv2.cvtColor(f, cv2.COLOR_BGR2GRAY)
            x0, y0 = max(0, cx - 170), max(0, cy - 130)
            win = cv2.GaussianBlur(g[y0:cy + 130, x0:cx + 170], (9, 9), 0)
            msk = (win <= np.percentile(win, 8)).astype(np.uint8)     # darkest 8 % = eye opening (and dark fur)
            k, lab, st, cen = cv2.connectedComponentsWithStats(msk, 8)
            cand = [(np.hypot(cen[i][0] + x0 - cx, cen[i][1] + y0 - cy), i) for i in range(1, k) if 300 < st[i, 4] < 20000]
            if cand:
                _, i = min(cand)
                pts.append((cen[i][0] + x0, cen[i][1] + y0))
        ex, ey = (int(round(np.median([p[0] for p in pts]))), int(round(np.median([p[1] for p in pts])))) if pts else (cx, cy)
        spread = float(np.median([np.hypot(p[0] - ex, p[1] - ey) for p in pts])) if pts else float("nan")
        w, h = sz.get((c, m), (200, 150))
        rows.append(dict(cage=c, mouse=m, prev_cx=cx, prev_cy=cy, cx=ex, cy=ey, shift_px=round(float(np.hypot(ex - cx, ey - cy)), 1),
                         spread_px=round(spread, 1), n_frames=len(pts), w=w, h=h))
        sets.append(f"--set {c} {m} {ex} {ey}")
        for fr in (int(n * 0.2), int(n * 0.5), int(n * 0.8)):
            cap.set(cv2.CAP_PROP_POS_FRAMES, fr)
            ok, f = cap.read()
            if not ok:
                continue
            x0, y0 = max(0, ex - 260), max(0, ey - 190)
            t = f[y0:y0 + 380, x0:x0 + 520].copy()
            cv2.rectangle(t, (ex - w // 2 - x0, ey - h // 2 - y0), (ex + w // 2 - x0, ey + h // 2 - y0), (0, 255, 255), 2)
            cv2.drawMarker(t, (cx - x0, cy - y0), (255, 0, 255), cv2.MARKER_CROSS, 16, 2)
            cv2.putText(t, f"{c} {m} f{fr}", (5, 22), cv2.FONT_HERSHEY_SIMPLEX, 0.65, (0, 0, 255), 2)
            tiles.append(cv2.resize(t, (390, 285)))
    if tiles:
        while len(tiles) % 6:
            tiles.append(np.full_like(tiles[0], 255))
        sheet = np.vstack([np.hstack(tiles[i:i + 6]) for i in range(0, len(tiles), 6)])
        cv2.imwrite(str(out / f"day{a.day}_eye_box_proposal.jpg"), sheet, [cv2.IMWRITE_JPEG_QUALITY, 85])
    with open(out / f"day{a.day}_eye_box_proposal.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["cage"])
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print(f"  {r['cage']} {r['mouse']}: ({r['cx']}, {r['cy']})  shift from day {a.prev or a.day - 1} {r['shift_px']} px, "
              f"spread {r['spread_px']} px over {r['n_frames']} frames")
    print(out / f"day{a.day}_eye_box_proposal.jpg")
    print("python tools/eye_box_fixed.py --day", a.day, " ".join(sets), "--apply")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
