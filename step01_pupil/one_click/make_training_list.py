"""make_training_list.py (2026-11-12) - pick the image/mask pairs for ONE general pupil U-Net from the per-video runs.

Input: runs.csv from inventory_existing_pupil.py. A run is used when its own validation Dice >= 0.90, it has >= 100 pairs
and its mouse is known (ROI matches that colour's roi.mat). Up to PER_RUN pairs are taken per run, evenly spread over its
frames, so every mouse-day weighs about the same.
Split: HOLDOUT mice are left out of training completely (their runs -> test_pairs.csv) to test the model on animals it has
never seen, like a new cohort. Write --all to also make the list for the final model that uses all 14 mice.
Outputs (work dir): train_pairs.csv, test_pairs.csv (or all_pairs.csv)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import paths_config as pc  # noqa: E402
W = pc.get("work_dir")
HOLDOUT = {("0910", "red"), ("0911", "white"), ("6100", "orange")}   # one per cage except 6099; Active+Passive mix
PER_RUN, PER_TEST_RUN = 120, 30
R = pd.read_csv(W / "runs.csv", dtype={"cage": str, "color": str})
R["cage"] = R["cage"].str.zfill(4)
good = R[(R.val_dice >= 0.90) & (R.n_labels >= 100) & R.color.notna() & (R.color != "")].copy()


def pairs(run, k):
    lab = sorted((Path(run) / "labels").glob("frame_*_mask.png"))
    out = []
    for p in [lab[i] for i in np.unique(np.linspace(0, len(lab) - 1, min(k, len(lab))).round().astype(int))]:
        img = p.with_name(p.name.replace("_mask", "")).parent.parent / "images" / p.name.replace("_mask", "")
        if img.exists():
            out.append((str(img), str(p)))
    return out


def build(df, k):
    rows = []
    for _, r in df.iterrows():
        for img, lab in pairs(r.run, k):
            rows.append(dict(image=img, label=lab, day_index=r.day_index, cage=r.cage, color=r.color, run=r.run))
    return pd.DataFrame(rows)


if "--all" in sys.argv:
    A = build(good, PER_RUN)
    A.to_csv(W / "all_pairs.csv", index=False)
    print(f"all_pairs: {len(A)} pairs from {A.run.nunique()} runs, {A.groupby(['cage','color']).ngroups} mice")
else:
    ho = good.apply(lambda r: (r.cage, r.color) in HOLDOUT, axis=1)
    T = build(good[~ho], PER_RUN)
    V = build(good[ho], PER_TEST_RUN)
    T.to_csv(W / "train_pairs.csv", index=False)
    V.to_csv(W / "test_pairs.csv", index=False)
    print(f"train: {len(T)} pairs, {T.run.nunique()} runs, {T.groupby(['cage','color']).ngroups} mice")
    print(f"test (held-out mice {sorted(HOLDOUT)}): {len(V)} pairs, {V.run.nunique()} runs")
