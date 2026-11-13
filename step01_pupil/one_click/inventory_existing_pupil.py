"""inventory_existing_pupil.py (2026-11-12) - map every finished session to the files it was made from.

For each <session_root>\\dayN\\<cage>\\<color>\\concat_out_*\\summary.txt (the alignment that Step 02 used):
  digital.csv path, final pupil CSV path, nPupil, TTL edges   (from summary.txt)
  the pupil VIDEO = the .avi in <video_root>\\dayN\\<cage>\\<color>\\ whose frame count is closest to the TTL camera edges
  ROI = <video_root>\\dayN\\<cage>\\<color>\\roi.mat
And every U-Net training run <video_root>\\dayN\\<cage>\\training_data1\\run_* with its ROI, pair count, validation Dice, and the
mouse it belongs to (ROI identical to that colour's roi.mat - runs do not record their video).
Writes <work>\\sessions.csv and <work>\\runs.csv. Read-only on both roots.
Usage: python -X utf8 inventory_existing_pupil.py [session_root] [video_root] [work_dir]   (defaults: one_click_paths.json)
"""
import glob
import os
import re
import sys
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import scipy.io as sio

import paths_config as pc  # noqa: E402
K = pc.get("session_root", 1)
H = pc.get("video_root", 2)
W = pc.get("work_dir", 3)
W.mkdir(parents=True, exist_ok=True)


def read_roi(p):
    try:
        return [int(v) for v in np.ravel(sio.loadmat(str(p))["roi"])]
    except Exception:
        return None


def nframes(p):
    v = cv2.VideoCapture(str(p))
    n = int(v.get(cv2.CAP_PROP_FRAME_COUNT))
    v.release()
    return n


rows = []
for s in sorted(K.glob("day*/*/*/concat_out_*/summary.txt")):
    day, cage, color = s.parts[-5], s.parts[-4], s.parts[-3]
    txt = s.read_text(encoding="utf-8", errors="replace")
    g = lambda pat: (re.search(pat, txt).group(1).strip() if re.search(pat, txt) else "")
    r = dict(day=day, day_index=int(day[3:]), cage=cage, color=color, summary=str(s),
             digital=g(r"Digital file:\s*(.+)"), pupil_csv=g(r"Pupil file:\s*(.+)"),
             n_pupil=int(g(r"nPupil rows:\s*(\d+)") or -1), ttl_edges=int(g(r"TTL rising edges \(raw\)=(\d+)") or -1),
             strategy=g(r"Chosen strategy:\s*(.+)"))
    vids = sorted(glob.glob(str(H / day / cage / color / "*.avi")))
    fr = [nframes(v) for v in vids]
    r["videos"] = ";".join(os.path.basename(v) for v in vids)
    r["video_frames"] = ";".join(map(str, fr))
    if vids and r["ttl_edges"] > 0:
        k = int(np.argmin([abs(f - r["ttl_edges"]) for f in fr]))
        r["pupil_video"], r["pupil_video_frames"] = vids[k], fr[k]
        r["frame_minus_ttl"] = fr[k] - r["ttl_edges"]
    roi = read_roi(H / day / cage / color / "roi.mat")
    r["roi"] = "" if roi is None else " ".join(map(str, roi))
    r["pupil_csv_exists"] = os.path.exists(r["pupil_csv"])
    rows.append(r)
    print(day, cage, color, r.get("pupil_video_frames"), r["ttl_edges"], r["roi"], flush=True)
S = pd.DataFrame(rows)
S.to_csv(W / "sessions.csv", index=False)

runs = []
for rd in sorted(H.glob("day*/*/training_data1/run_*")):
    day, cage = rd.parts[-4], rd.parts[-3]
    if not re.fullmatch(r"day\d+", day):
        continue
    imgs = len(list((rd / "images").glob("frame_*.png"))) if (rd / "images").is_dir() else 0
    labs = len(list((rd / "labels").glob("frame_*_mask.png"))) if (rd / "labels").is_dir() else 0
    dice = np.nan
    vs = rd / "step2_val_summary.txt"
    if vs.exists():
        m = re.search(r"Dice on \d+ pairs = ([\d.]+)", vs.read_text(errors="replace"))
        dice = float(m.group(1)) if m else np.nan
    roi = None
    try:
        c = sio.loadmat(str(rd / "step1_config.mat"), variable_names=["roi"])
        roi = [int(v) for v in np.ravel(c["roi"])]
    except Exception:
        pass
    color = ""
    if roi is not None:
        for cdir in (H / day / cage).iterdir():
            if cdir.is_dir() and (cdir / "roi.mat").exists() and read_roi(cdir / "roi.mat") == roi:
                color = cdir.name
                break
    runs.append(dict(day=day, day_index=int(day[3:]), cage=cage, color=color, run=str(rd), n_images=imgs, n_labels=labs,
                     val_dice=dice, roi="" if roi is None else " ".join(map(str, roi)),
                     has_model=(rd / "trainedPupilNet.mat").exists()))
R = pd.DataFrame(runs)
R.to_csv(W / "runs.csv", index=False)
print(f"\n{len(S)} sessions, {len(R)} runs -> {W}")
print("sessions: pupil csv exists", int(S.pupil_csv_exists.sum()), "| video frame-TTL diff counts:",
      S["frame_minus_ttl"].value_counts().head(8).to_dict())
print("runs with colour matched:", int((R.color != "").sum()), "| total label pairs:", int(R.n_labels.sum()),
      "| runs Dice>=0.9 & >=100 pairs:", int(((R.val_dice >= 0.9) & (R.n_labels >= 100)).sum()))
