#!/usr/bin/env python3
"""
eye_box_fixed.py - one eye-box size per mouse for all days (added 2026-11-20; the user answered "네, 모두 다시" to
"fix one eye-box size per mouse and re-track days 1-7").

Why: the user's package crops the eye box and resizes it to the U-Net's square 128 x 128 input, so the box SIZE changes
the pupil size it reads (same video, d5 9670f green2: a 240 -> 185 px wide box reads -3.0 px = -12 %; d6 9669m green2
-2.0 px), while moving a box of the same size onto the eye does not (-0.02 / +0.3 px) and cuts the errors. So every
day uses the mouse's size from video/eye_box_sizes.csv, and only the position follows the eye: centred on the eye
opening, which Claude reads on the day's frames before the pupil step (mice sit differently every day).

usage:
  python eye_box_fixed.py --day N --set CAGE MOUSE CX CY [--set ...] [--why TEXT] [--apply]
        writes video/dayN/<cage>/<mouse>/roi.mat = the mouse's fixed size centred at (CX, CY), the eye opening's
        centre in full-frame pixels; an existing roi.mat is kept as roi_before_<date>.mat; one line per box goes to
        video/roi_log.csv. Dry run unless --apply.
  python eye_box_fixed.py --check [--day N]
        lists every session (video folder with a video) whose roi.mat is missing or not its mouse's fixed size
        (CHECK lines, exit 1 if any). Mice without a fixed size (not in the csv) are listed as info only.
The daily routine (daily_summary.py) runs --check for the day before the pupil step (report only).
eye_box_sizes.csv columns: cage,mouse,w,h,since,reason
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

T = Path(__file__).resolve().parent
P = T.parent
FRAME_W, FRAME_H = 1920, 1080
VIDEO_EXT = (".mp4", ".avi")


def sizes(project: Path) -> dict:
    f = project / "video" / "eye_box_sizes.csv"
    if not f.exists():
        return {}
    with open(f, newline="", encoding="utf-8-sig") as h:
        return {(r["cage"].strip(), r["mouse"].strip()): (int(r["w"]), int(r["h"])) for r in csv.DictReader(h)}


def box_at(cx: float, cy: float, w: int, h: int) -> list:
    """MATLAB imcrop box [x y w h] (1-based) of size w x h centred at full-frame pixel (cx, cy), kept in the frame"""
    x = int(round(cx - w / 2)) + 1
    y = int(round(cy - h / 2)) + 1
    return [min(max(1, x), FRAME_W - w), min(max(1, y), FRAME_H - h), w, h]


def read_roi(f: Path):
    try:
        return [int(round(v)) for v in np.asarray(loadmat(f)["roi"], float).ravel()[:4]]
    except Exception:
        return None


def sessions(project: Path, day: int | None):
    days = [day] if day else sorted(int(p.name[3:]) for p in (project / "video").glob("day*") if p.name[3:].isdigit())
    for d in days:
        for vd in sorted((project / "video" / f"day{d}").glob("*/*")):
            if vd.is_dir() and any(p.suffix.lower() in VIDEO_EXT for p in vd.iterdir()):
                yield d, vd.parent.name, vd.name, vd


def check(project: Path, day: int | None) -> int:
    S = sizes(project)
    bad = 0
    print(f"eye boxes vs the fixed size per mouse ({project / 'video' / 'eye_box_sizes.csv'})")
    if not S:
        print("  CHECK: no eye_box_sizes.csv - no fixed sizes yet")
        return 1
    for d, c, m, vd in sessions(project, day):
        roi = read_roi(vd / "roi.mat") if (vd / "roi.mat").exists() else None
        want = S.get((c, m))
        if want is None:
            print(f"  info  day{d} {c} {m}: no fixed size for this mouse (box {roi})")
        elif roi is None:
            bad += 1
            print(f"  CHECK day{d} {c} {m}: no roi.mat - set it with --set before the pupil step")
        elif tuple(roi[2:4]) != want:
            bad += 1
            print(f"  CHECK day{d} {c} {m}: box {roi[2]}x{roi[3]} is not the mouse's fixed size {want[0]}x{want[1]}")
        else:
            print(f"  ok    day{d} {c} {m}: box {' '.join(map(str, roi))}")
    return 1 if bad else 0


def set_boxes(project: Path, day: int, items: list, why: str, apply: bool) -> int:
    S = sizes(project)
    stamp = datetime.now().strftime("%Y%m%d")
    rc = 0
    for c, m, cx, cy in items:
        if (c, m) not in S:
            print(f"  !! {c} {m}: no fixed size in eye_box_sizes.csv")
            rc = 1
            continue
        w, h = S[(c, m)]
        vd = project / "video" / f"day{day}" / c / m
        if not vd.is_dir():
            print(f"  !! day{day} {c} {m}: no video folder {vd}")
            rc = 1
            continue
        box = box_at(float(cx), float(cy), w, h)
        rf = vd / "roi.mat"
        before = read_roi(rf) if rf.exists() else None
        print(f"  day{day} {c} {m}: eye centre ({cx}, {cy}) -> box {box}" + (f" (before {before})" if before else ""))
        if not apply:
            continue
        if rf.exists():
            keep = vd / f"roi_before_{stamp}.mat"
            if not keep.exists():
                shutil.copy2(rf, keep)
        savemat(str(rf), {"roi": np.array([box], dtype=float)})
        with open(project / "video" / "roi_log.csv", "a", newline="", encoding="utf-8") as hh:
            csv.writer(hh).writerow([datetime.now().isoformat(timespec="seconds"), f"day{day}", c, m,
                                     " ".join(map(str, box)),
                                     f"fixed size {w}x{h} (eye_box_sizes.csv) centred on the eye opening ({cx}, {cy}) "
                                     f"by eye_box_fixed.py; {why}" + (f"; before: {' '.join(map(str, before))}"
                                                                      if before else "")])
    if not apply:
        print("  (dry run: add --apply to write)")
    return rc


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", default=str(P))
    ap.add_argument("--day", type=int)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--set", nargs=4, action="append", metavar=("CAGE", "MOUSE", "CX", "CY"))
    ap.add_argument("--why", default="eye opening read on the day's frames by Claude")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args(argv)
    project = Path(a.project)
    if a.set:
        if not a.day:
            ap.error("--set needs --day")
        return set_boxes(project, a.day, a.set, a.why, a.apply)
    return check(project, a.day)


if __name__ == "__main__":
    sys.exit(main())
