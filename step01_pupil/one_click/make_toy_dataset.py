"""make_toy_dataset.py (2026-11-12) - a tiny copy of real sessions in the one-click folder layout, for quick tests.

For each toy session: the first NT frames of BOTH camera videos (so the runner has to pick the pupil one by its frame
count), roi.mat (left out for one session to test ROIMode 'previous'), and digital.csv cut to exactly NT camera TTL rising
edges. One TTL file is put in a 'ttl' sub-folder of the video day instead of the session folder, to test the TTL search.
Also saves the matching rows of the existing (per-video U-Net) CLEAN pupil CSV as a reference.
Old step 3 skipped video frame 1, so reference row k = video frame k+1.

Usage: python -X utf8 make_toy_dataset.py [out_dir] [n_frames]     (reads video_root, writes only out_dir)
"""
import shutil
import sys
from pathlib import Path

import cv2
import pandas as pd

import paths_config as pc  # noqa: E402
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else pc.get("work_dir") / "toy"
NT = int(sys.argv[2]) if len(sys.argv) > 2 else 900
H = pc.get("video_root")
S = pd.read_csv(pc.get("work_dir") / "sessions.csv", dtype=str)
S["cage"] = S["cage"].str.zfill(4)
# (day, cage, color, keep roi.mat, TTL location)   - mice 0910_red and 6100_orange were NOT used to train the hold-out model
TOY = [("day12", "0910", "red", True, "session"), ("day13", "0910", "red", False, "video_ttl"),
       ("day12", "6100", "orange", True, "session")]


def cut_video(src, dst, n, fps):
    v = cv2.VideoCapture(str(src))
    w, h = int(v.get(3)), int(v.get(4))
    out = cv2.VideoWriter(str(dst), cv2.VideoWriter_fourcc(*"MJPG"), fps, (w, h))
    k = 0
    while k < n:
        ok, f = v.read()
        if not ok:
            break
        out.write(f)
        k += 1
    v.release()
    out.release()
    return k


for day, cage, color, keep_roi, ttl_where in TOY:
    r = S[(S.day == day) & (S.cage == cage) & (S.color == color)].iloc[0]
    vdir = OUT / "videos" / day / cage / color
    vdir.mkdir(parents=True, exist_ok=True)
    src_dir = H / day / cage / color
    for vid in sorted(src_dir.glob("*.avi")):
        is_pupil = str(vid) == r.pupil_video
        n = cut_video(vid, vdir / vid.name, NT if is_pupil else NT + 60, 40.0)   # other camera: 60 extra frames
        print(f"{day} {cage}_{color}: {vid.name} -> {n} frames {'(pupil)' if is_pupil else '(other camera)'}")
    if keep_roi:
        shutil.copy2(src_dir / "roi.mat", vdir / "roi.mat")
    D = pd.read_csv(r.digital)
    cam = D["Channel 6"].to_numpy() > 0.5
    rise = [i for i in range(1, len(cam)) if cam[i] and not cam[i - 1]] if not cam[0] else [0] + \
           [i for i in range(1, len(cam)) if cam[i] and not cam[i - 1]]
    cut = rise[NT] if len(rise) > NT else len(D)
    tdir = (OUT / "sessions" / day / cage / color) if ttl_where == "session" else (OUT / "videos" / day / "ttl" / cage / color)
    tdir.mkdir(parents=True, exist_ok=True)
    D.iloc[:cut].to_csv(tdir / "digital.csv", index=False, float_format="%.9f")
    ref = pd.read_csv(r.pupil_csv).iloc[: NT - 1]
    (OUT / "reference").mkdir(exist_ok=True)
    ref.to_csv(OUT / "reference" / f"{day}_{cage}_{color}_reference_CLEAN.csv", index=False)
    print(f"   digital.csv: {min(NT, len(rise))} camera edges -> {tdir}; roi.mat {'copied' if keep_roi else 'LEFT OUT'}")
print("toy dataset ->", OUT)
