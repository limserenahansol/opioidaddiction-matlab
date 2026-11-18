#!/usr/bin/env python3
"""
wholeeye_clean.py - the whole-eye trace in the package's CLEAN format, for mice with eye damage (user 2026-11-17:
"9670 black has eye damage ... use whole eye as alternate", "9669 black also has eye damages so ... both";
"네, 넣어 주세요" = into the tables).

For every session that eye_track_all.m tracked (analysis/wholeeye/day<N>/<cage>/<mouse>/wholeeye_frames.csv) and
whose mouse is approved in <pupil tree>/wholeeye_mice.csv (or the session in wholeeye_sessions.csv):
  <pupil tree>/day<N>/<cage>/<mouse>/pupil_outputs_wholeeye/wholeeye_CLEAN.csv with the package CLEAN columns
    Timestamp_seconds   (video frame - 1) / 60, as the package's CLEAN file (the mp4 stamps frames 1/60 s apart)
    Diameter_pixels     equivalent diameter of the eye opening = 2 * sqrt(area / pi), NaN where no eye was found
    CenterX_px/_Y_px    the eye's centre, 1-based in the package's eye-box crop of that day (as the package writes)
    UsedScoreFallback   0
    SmoothedDiameter_px Diameter_pixels with every gap interpolated (edges held) and a 5-frame running median - gap-
                        filled like the package's SmoothedDiameter_px, which its Step 02 table carries
  plus EyeArea_px2, EyeMajor_px, EyeMinor_px (for reference; not used by the tables)
  and wholeeye_source.json (SHA-1 of the input, model, crop, frames, % found). An up-to-date file is not rewritten.
One row per video frame. The pupil tree = concat_choice.active_pupil_root (pupil_active_root.txt).
usage: python wholeeye_clean.py --project <project> [--dry-run]
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import concat_choice as CC  # noqa: E402


def one(src: Path, sd: Path, dry: bool) -> str:
    out_dir = sd / CC.WHOLE_DIR
    out, js = out_dir / CC.WHOLE_CSV, out_dir / "wholeeye_source.json"
    h = CC.sha1(src)
    if out.exists() and js.exists():
        try:
            if json.loads(js.read_text(encoding="utf-8")).get("frames_csv_sha1") == h:
                return "up to date"
        except ValueError:
            pass
    qc = sd / "pupil_outputs_general" / "qc.json"
    if not qc.exists():
        raise RuntimeError("no package qc.json (eye box) in this session folder")
    roi = [float(v) for v in json.loads(qc.read_text(encoding="utf-8"))["roi"]]
    E = pd.read_csv(src)
    if dry:
        return f"would write {CC.WHOLE_DIR}/{CC.WHOLE_CSV} ({len(E)} frames)"
    E = E.sort_values("frame").reset_index(drop=True)
    if not np.array_equal(E["frame"].to_numpy(int), np.arange(1, len(E) + 1)):
        raise RuntimeError("frames are not 1..N in wholeeye_frames.csv")
    found = E["found"].to_numpy(float) > 0
    area = np.where(found, E["area_px"].to_numpy(float), np.nan)
    d = 2.0 * np.sqrt(area / np.pi)
    sm = pd.Series(d).interpolate(limit_area="inside").ffill().bfill()
    sm = sm.rolling(5, center=True, min_periods=1).median()
    # package convention: CenterX_px is 1-based in the eye-box crop, whose 0-based origin is (roi_x - 1, roi_y - 1)
    cx = np.where(found, E["cx"].to_numpy(float) - (roi[0] - 1) + 1, np.nan)
    cy = np.where(found, E["cy"].to_numpy(float) - (roi[1] - 1) + 1, np.nan)
    T = pd.DataFrame(dict(Timestamp_seconds=(E["frame"].to_numpy(float) - 1) / 60.0, Diameter_pixels=d,
                          CenterX_px=cx, CenterY_px=cy, UsedScoreFallback=0, SmoothedDiameter_px=sm.to_numpy(),
                          EyeArea_px2=area, EyeMajor_px=np.where(found, E["major_px"], np.nan),
                          EyeMinor_px=np.where(found, E["minor_px"], np.nan)))
    pkg = sd / "pupil_outputs_general" / "pupil_diameter_FINAL_U-Net_results_CLEAN.csv"
    n_pkg = len(pd.read_csv(pkg, usecols=[0])) if pkg.exists() else None
    out_dir.mkdir(exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    T.to_csv(tmp, index=False, float_format="%.6g")
    tmp.replace(out)
    q = {}
    qj = src.parent / "wholeeye_qc.json"
    if qj.exists():
        try:
            q = json.loads(qj.read_text(encoding="utf-8"))
        except ValueError:
            q = {}
    js.write_text(json.dumps(dict(
        tool="wholeeye_clean.py", made=datetime.now().isoformat(timespec="seconds"), frames_csv=str(src),
        frames_csv_sha1=h, clean_sha1=CC.sha1(out), model=q.get("model"), crop=q.get("crop"), package_roi=roi,
        frames=int(len(T)), package_clean_rows=n_pkg, found_pct=round(100 * float(found.mean()), 2),
        diameter="equivalent diameter of the eye opening, 2*sqrt(area/pi), px",
        note="whole eye instead of the pupil (eye damage); user 2026-11-17"), indent=1), encoding="utf-8")
    warn = "" if n_pkg in (None, len(T)) else f"  !! {len(T)} frames, the package's CLEAN file has {n_pkg} rows"
    return f"wrote {CC.WHOLE_CSV} ({len(T)} frames, eye found {100 * found.mean():.1f} %){warn}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    project = Path(a.project)
    root = CC.active_pupil_root(project)
    keys = CC.approved_sessions(root, CC.LISTS[CC.WHOLE])
    bad = 0
    for src in sorted((project / "analysis" / "wholeeye").glob("day*/*/*/wholeeye_frames.csv")):
        day_s, cage, mouse = src.parts[-4], src.parts[-3], src.parts[-2]
        if not day_s[3:].isdigit():
            continue
        k = CC.norm_key(day_s[3:], cage, mouse)
        sd = CC.session_dir(root, k[0], cage, mouse)
        if k not in keys or sd is None:
            print(f"  {day_s} {cage}_{mouse}: not approved for the whole eye in {root.name} (or no session folder) - skipped")
            continue
        try:
            msg = one(src, sd, a.dry_run)
        except Exception as exc:  # noqa: BLE001
            msg, bad = f"FAILED ({exc})", bad + 1
        print(f"  whole eye {day_s} {cage}_{mouse}: {msg}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
