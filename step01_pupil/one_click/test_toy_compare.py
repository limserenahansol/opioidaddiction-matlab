"""test_toy_compare.py (2026-11-12) - compare one-click toy outputs with the existing per-video U-Net results.

For each toy session: new  <toy>/out/dayN/<cage>/<color>/pupil_outputs_general/pupil_diameter_FINAL_U-Net_results_CLEAN.csv
vs reference <toy>/reference/dayN_<cage>_<color>_reference_CLEAN.csv (made by step1+step3 for that video).
Old step 3 skipped video frame 1, so new frame k pairs with reference row k-1. The first 60 frames are left out of the
summary (both trackers are still in their 30-frame warm-up / EMA start). Also checks the concat_out (rows, alignment, licks).
Usage: python -X utf8 test_toy_compare.py [toy_dir]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import paths_config as pc  # noqa: E402
T = Path(sys.argv[1]) if len(sys.argv) > 1 else pc.get("work_dir") / "toy"
rows = []
for ref_p in sorted((T / "reference").glob("*_reference_CLEAN.csv")):
    day, cage, color = ref_p.name.split("_")[:3]
    od = T / "out" / day / cage / color
    new_p = od / "pupil_outputs_general" / "pupil_diameter_FINAL_U-Net_results_CLEAN.csv"
    if not new_p.exists():
        rows.append(dict(session=f"{day} {cage}_{color}", status="no output"))
        continue
    new = pd.read_csv(new_p)["SmoothedDiameter_px"].to_numpy()[1:]      # frames 2..N
    ref = pd.read_csv(ref_p)["SmoothedDiameter_px"].to_numpy()
    n = min(len(new), len(ref))
    a, b = new[60:n], ref[60:n]
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    cc = sorted(od.glob("concat_out_*"))
    comb = pd.read_csv(cc[-1] / "combined_pupil_digital.csv") if cc else None
    summ = (cc[-1] / "summary.txt").read_text(encoding="utf-8", errors="replace") if cc else ""
    dig = None
    for cand in [T / "sessions" / day / cage / color / "digital.csv", T / "videos" / day / "ttl" / cage / color / "digital.csv"]:
        if cand.exists():
            dig = pd.read_csv(cand)
    lick_edges = int(((dig["Channel 7"] > 0.5) & ~(dig["Channel 7"].shift(1, fill_value=0) > 0.5)).sum()) if dig is not None else -1
    rows.append(dict(
        session=f"{day} {cage}_{color}", frames_compared=int(ok.sum()),
        pearson_r=round(float(np.corrcoef(a, b)[0, 1]), 3),
        median_new_px=round(float(np.median(a)), 2), median_ref_px=round(float(np.median(b)), 2),
        median_abs_diff_pct=round(float(np.median(np.abs(a - b) / b) * 100), 1),
        within_10pct=round(float(np.mean(np.abs(a - b) / b <= 0.10) * 100), 1),
        concat_rows=(len(comb) if comb is not None else 0),
        align=(summ.split("Alignment mode:")[1].split("\n")[0].strip() if "Alignment mode:" in summ else ""),
        lick_onsets=(int(((comb["Lick_TTL"].fillna(0) > 0.5) & ~(comb["Lick_TTL"].fillna(0).shift(1, fill_value=0) > 0.5)).sum())
                     if comb is not None else -1),
        lick_edges_in_ttl=lick_edges))
R = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(R.to_string(index=False))
R.to_csv(T / "toy_compare.csv", index=False)
print("->", T / "toy_compare.csv")
