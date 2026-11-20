#!/usr/bin/env python3
"""
pupil_box_swap.py - put a better pupil tracking of a session into the real pupil tree (added 2026-11-20; user rule:
"always you need to use the best way (less error making tracker) to get pupil diameter tracking so if you find better
you can do use better for real data and always use best data for real data for analysis").

The better tracking is the user's package itself (same code, model and settings), run on the same video with a better
eye box in the test tree analysis/pupil_boxtest (box_test runs: run_all_config_teach_boxtest_*.json). Claude checks it
first, three ways: fewer frames out, the same pupil where both runs are fine (centre in full-frame pixels), and the
check pictures. Then the session gets a row in <pupil tree>/box_swaps.csv and this tool, per row not yet applied:
  1. moves the session's pupil_outputs_general, pupil_outputs_darkcontrast and every concat_out_* folder to
     <session>/_superseded_<date>_box/ (nothing is deleted)
  2. copies the test run's pupil_outputs_general and concat_out_* into the session folder; in qc.json the test video
     path (a hard link) becomes the real video path (video/dayN/<cage>/<mouse>/<same name>)
  3. keeps the old eye box as video/dayN/<cage>/<mouse>/roi_before_<date>_box.mat, writes the test box to roi.mat and
     adds a line to video/roi_log.csv
  4. writes <session>/box_swap.json (what, when, why, the numbers) and marks the row applied
Derived copies (alignment fix, dark-contrast repair, ...) are made again by the daily routine on the new table where
their list still approves the session. Afterwards run the daily routine (daily_summary.py --no-arrange) so every
table, the pupil check pictures, step02_base and Step 02 follow.

box_swaps.csv columns: day,cage,mouse,box,test_root,out_before_pct,out_after_pct,centre_diff_px,decided,reason,applied
usage: python pupil_box_swap.py [--apply]        (dry run unless --apply)
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

T = Path(__file__).resolve().parent
P = T.parent
sys.path.insert(0, str(T))
import concat_choice as CC  # noqa: E402

LIST = "box_swaps.csv"
COLS = ["day", "cage", "mouse", "box", "test_root", "out_before_pct", "out_after_pct", "centre_diff_px", "decided",
        "reason", "applied"]


def read_rows(root: Path):
    f = root / LIST
    if not f.exists():
        return []
    with open(f, newline="", encoding="utf-8-sig") as h:
        return list(csv.DictReader(h))


def write_rows(root: Path, rows):
    with open(root / LIST, "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=COLS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in COLS})


def one(root: Path, r: dict, apply: bool, stamp: str):
    d, c, m = int(r["day"]), r["cage"], r["mouse"]
    box = [float(v) for v in str(r["box"]).split()]
    test = Path(r["test_root"]) if Path(r["test_root"]).is_absolute() else P / r["test_root"]
    sd = CC.session_dir(root, d, c, m)
    st = CC.session_dir(test, d, c, m)
    tag = f"day{d} {c} {m}"
    if sd is None or st is None:
        print(f"  !! {tag}: session folder missing (real {sd}, test {st})")
        return False
    tq = st / "pupil_outputs_general" / "qc.json"
    truns = [p for p in st.glob("concat_out_*") if p.is_dir() and (p / "combined_pupil_digital.csv").exists()]
    if not tq.exists() or not truns:
        print(f"  !! {tag}: the test run is not complete ({st})")
        return False
    q = json.loads(tq.read_text(encoding="utf-8"))
    tvid = Path(q.get("video", ""))
    vdir = P / "video" / f"day{d}" / c / m
    rvid = vdir / tvid.name
    if not rvid.exists():
        print(f"  !! {tag}: real video {rvid} not found")
        return False
    if [round(v) for v in q.get("roi", [])] != [round(v) for v in box]:
        print(f"  !! {tag}: the test run's box {q.get('roi')} is not the listed box {box}")
        return False
    sup = sd / f"_superseded_{stamp}_box"
    old = [p for p in sd.iterdir() if p.is_dir() and (p.name.startswith("concat_out_") or
                                                      p.name in ("pupil_outputs_general", "pupil_outputs_darkcontrast"))]
    print(f"  {tag}: move {len(old)} old folder(s) -> {sup.name}; copy {len(truns)} run(s) + pupil_outputs_general from "
          f"the test tree; roi.mat -> {[int(v) for v in box]}")
    if not apply:
        return True
    sup.mkdir(exist_ok=True)
    for p in old:
        shutil.move(str(p), str(sup / p.name))
    shutil.copytree(st / "pupil_outputs_general", sd / "pupil_outputs_general")
    for p in truns:
        shutil.copytree(p, sd / p.name)
    q["video"] = str(rvid)
    (sd / "pupil_outputs_general" / "qc.json").write_text(json.dumps(q, indent=2), encoding="utf-8")
    rf = vdir / "roi.mat"
    keep = vdir / f"roi_before_{stamp}_box.mat"
    before = " ".join(str(int(v)) for v in loadmat(rf)["roi"].ravel()) if rf.exists() else ""
    if rf.exists() and not keep.exists():
        shutil.copy2(rf, keep)
    savemat(str(rf), {"roi": np.array([box], dtype=float)})
    with open(P / "video" / "roi_log.csv", "a", newline="", encoding="utf-8") as h:
        csv.writer(h).writerow([datetime.now().isoformat(timespec="seconds"), f"day{d}", c, m,
                                " ".join(str(int(v)) for v in box),
                                f"better pupil tracking with this eye-centred box (box test {r.get('test_root')}: out "
                                f"{r.get('out_before_pct')} -> {r.get('out_after_pct')} %, same pupil); swapped into the "
                                f"real tree by pupil_box_swap.py ({r.get('decided')}); before: {before} "
                                f"(kept as {keep.name})"])
    note = dict(when=datetime.now().isoformat(timespec="seconds"), tool="tools/pupil_box_swap.py",
                rule="user 2026-11-20: always use the best (least-error) pupil tracking for the real data",
                box_before=before, box_after=" ".join(str(int(v)) for v in box), test_session=str(st),
                test_video=str(tvid), real_video=str(rvid), superseded=str(sup),
                out_before_pct=r.get("out_before_pct"), out_after_pct=r.get("out_after_pct"),
                centre_diff_px=r.get("centre_diff_px"), reason=r.get("reason"), decided=r.get("decided"))
    (sd / "box_swap.json").write_text(json.dumps(note, indent=2), encoding="utf-8")
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--project", default=str(P))
    a = ap.parse_args(argv)
    root = CC.active_pupil_root(Path(a.project))
    rows = read_rows(root)
    todo = [r for r in rows if not str(r.get("applied", "")).strip()]
    print(f"box swaps: {len(rows)} row(s) in {root / LIST}, {len(todo)} to apply" + ("" if a.apply else " (dry run)"))
    stamp = datetime.now().strftime("%Y%m%d")
    ok_all = True
    for r in todo:
        ok = one(root, r, a.apply, stamp)
        ok_all &= ok
        if ok and a.apply:
            r["applied"] = datetime.now().isoformat(timespec="seconds")
    if a.apply and todo:
        write_rows(root, rows)
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
