#!/usr/bin/env python3
"""
straub_day.py - Straub tail videos of one day, scored by hand like the hot plate (added 2026-11-21; user: "ill also do
straub tail analysis on monday using similar as manual scoring for hotplate what we used").

How the previous cohort did it (checked 2026-11-21 in the archive):
  - days 10 (morphine) and 16 (withdrawal); one video per mouse, about 5 min (8608-10098 frames at 30 fps), named like
    0910black2straub.mp4;
  - scored with the lab's manual scorer (step01_manual_assays/manual_scoring_video_multibehavior_batch2.m): key 'a'
    HELD while the tail is up (Straub tail = tail raised > 30 deg from the body axis), no key otherwise. The scorer names
    key 'a' 'Licking', so in Step 02: STRAUB_Pct_Licking = Straub tail %, STRAUB_Pct_Non_moving = Straub-FREE %
    (e.g. 0910 black: day 10 Licking 52.95 % / Non-moving 47.05 %; day 16 2.00 % / 98.00 %).
  - Step 02 reads assays/day<N>/<cage>/<mouse>/ManualScoringResults*straub*.mat (prepare_step02_base.py copies it).

For day N, like the HOT step of assays_day.py (that file is not changed):
 1. ARRANGE: videos in --src (default: assay_video_dirs of tools/local_paths.json, else Logitech Capture's folder; 2
    levels) whose name has 'straub', a mouse of the cohort ('<cage><sex><colour>', e.g. 9670fgreen1) and 'day<N>' are
    COPIED (originals stay) to <project>/straub_video/day<N>/<cage>/<mouse>/. Files written in the last 2 min wait.
 2. SCORING COPIES: the first --seconds (300 = the 5-min observation) of each video, no audio, constant frame rate, in
    <project>/straub_video/day<N>_first5min/. Score them at real time (key 'a' held while the tail is up):
      cd('<project>/tools/opioidaddiction-matlab/current_cohort_workflow/matlab')
      manual_scoring_video_multibehavior_speed('<project>/straub_video/day<N>_first5min', 1)
 3. FILING: output/ManualScoringResults_*.mat -> assays/day<N>/<cage>/<mouse>/ (the name keeps 'straub').
 4. AUTOMATIC (added 2026-11-21 after the user confirmed it on an annotated video; skip with --no-auto): straub_auto_v8.py
    since 2026-11-22 (v8: body = dark fur hull, so a reflection is not body; curved tails kept; a tail leaving the top of
    the body = Straub; a vertical body = rearing / grooming is NOT COUNTED unless the tail is clearly raised above it -
    user rules and 5 answered frames of 2026-11-22; v6b = straub_auto.py / straub_angle.py, kept)
    on every scoring copy -> analysis/straub_auto/day<N>/ (per-frame tail angle; Straub = > 25 deg from the body axis,
    or the tail raised over the back; 0.5-s majority). Shown next to the hand score, never instead of it.
 5. TABLE: analysis/days/day<N>/straub_day<N>.csv: seconds scored, straub_pct (key 'a'), straub_free_pct (no key),
    other keys % (should be 0; '!!' if not), auto_straub_pct, auto_measured_pct, tail-immersion latency of the day.
Dry run unless --apply. Windows idle priority.
usage: python straub_day.py --day 10 [--src <folder> ...] [--seconds 300] [--apply]
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import time
import shutil
from pathlib import Path

T = Path(__file__).resolve().parent
P = T.parent
sys.path.insert(0, str(T))
import assays_day as AD  # noqa: E402  (MOUSE, VIDEO_EXT, default_src, idle, known_mice, ffmpeg, video_fps, ...)

KEY_STRAUB = 1          # score code of key 'a' (behaviorNames{1} = 'Licking' in the lab scorer)


def parse(name: str, mice: dict):
    """(is_straub, cage folder, mouse folder, day) from a file name; None parts when not understood"""
    low = name.lower()
    d = re.search(r"day(\d+)(?!\d)", low)
    m = AD.MOUSE.search(name)
    day = int(d.group(1)) if d else None
    if m is None:
        return "straub" in low, None, None, day
    cage4, color = m.group("cage"), m.group("color").lower()
    cages = [c for c in mice if c[:4] == cage4]
    cage = cages[0] if len(cages) == 1 else None
    if cage is not None and color not in mice[cage]:
        color = None
    return "straub" in low, cage, color, day


NAME_FIXES = P / "straub_video" / "name_fixes.csv"   # added 2026-11-22: opt-in list of recordings saved under a wrong name


def name_fixes() -> dict:
    """{original file name (lower case): (cage folder, mouse folder, day)} for the rows of straub_video/name_fixes.csv
    with an 'approved' cell (a recording the user said was saved under the wrong mouse name; the original is never
    renamed, its copy gets the right name). Added 2026-11-22 (user: '9670f_green2Straub_prday9.mp4 ... this is green1')"""
    if not NAME_FIXES.exists():
        return {}
    out = {}
    with open(NAME_FIXES, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if (r.get("approved") or "").strip() and (r.get("file") or "").strip():
                out[r["file"].strip().lower()] = (r["cage"].strip(), r["mouse"].strip(), int(r["day"]))
    return out


def arrange(day: int, srcs, mice, apply: bool, log: list):
    fixes = name_fixes()
    for src in srcs:
        if not src.is_dir():
            continue
        cands = [p for p in src.glob("*") if p.is_file()] + [p for p in src.glob("*/*") if p.is_file()]
        for v in sorted(cands):
            if v.suffix.lower() not in AD.VIDEO_EXT:
                continue
            fx = fixes.get(v.name.lower())
            if fx is not None:                            # listed in name_fixes.csv: the user's mouse, copy renamed
                ok, (cage, mouse, d) = True, fx
            else:
                ok, cage, mouse, d = parse(v.stem, mice)
            if not ok or d != day:
                continue
            if cage is None or mouse is None:
                log.append(f"  !! {v.name}: mouse not understood (known: {sorted(mice)}) - not copied")
                continue
            if time.time() - v.stat().st_mtime < 120:
                log.append(f"  {v.name}: written less than 2 min ago (still recording?) - next run")
                continue
            dst = P / "straub_video" / f"day{day}" / cage / mouse / (
                f"{cage}{mouse}_straub_day{d}{v.suffix.lower()}" if fx is not None else v.name)
            if fx is not None:
                log.append(f"  {v.name}: name fixed by name_fixes.csv -> {cage} {mouse} day {d}")
            if AD.same_file(v, dst):
                continue
            if dst.exists():
                log.append(f"  !! {v.name}: a different file of that name is already in {dst.parent} - not copied")
                continue
            log.append(f"  {'copy' if apply else 'would copy'} {v.name} -> {dst.relative_to(P).as_posix()}")
            if apply:
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(v, dst)


def copies(day: int, seconds: float, apply: bool, log: list):
    vids = sorted(p for p in (P / "straub_video" / f"day{day}").glob("*/*/*")
                  if p.is_file() and p.suffix.lower() in AD.VIDEO_EXT)
    sdir = P / "straub_video" / f"day{day}_first5min"
    ff = AD.ffmpeg()
    for v in vids:
        dst = sdir / f"{v.stem}_first5min.mp4"
        if dst.exists():
            continue
        if ff is None:
            log.append("  !! ffmpeg (tools/.venv imageio-ffmpeg) not found - no scoring copies")
            break
        fps = AD.video_fps(v)
        nfr = int(round(seconds * fps))
        log.append(f"  {'make' if apply else 'would make'} {dst.relative_to(P).as_posix()} (first {seconds:g} s = "
                   f"{nfr} frames at {fps:g} fps{'; CHECK: 60 fps in the header' if fps > 45 else ''})")
        if apply:
            sdir.mkdir(parents=True, exist_ok=True)
            subprocess.run([str(ff), "-hide_banner", "-loglevel", "error", "-y", "-i", str(v), "-frames:v", str(nfr),
                            "-vf", f"setpts=N/{fps:g}/TB", "-r", f"{fps:g}", "-an", "-c:v", "libx264", "-crf", "16",
                            "-preset", "veryfast", "-pix_fmt", "yuv420p", str(dst)],
                           creationflags=0x40 if sys.platform == "win32" else 0)
    return vids, sdir


def file_scores(day: int, sdir: Path, mice, apply: bool, log: list):
    out = sdir / "output"
    for m in sorted(out.glob("ManualScoringResults_*.mat")) if out.exists() else []:
        ok, cage, mouse, d = parse(m.stem, mice)
        if not ok or cage is None or mouse is None or d != day:
            log.append(f"  !! {m.name}: not a straub file of a known mouse / day {day} - not copied")
            continue
        AD.copy_if_changed(m, P / "assays" / f"day{day}" / cage / mouse / m.name, apply, log,
                           "Straub tail hand score, first 5 min")


def auto_step(day: int, sdir: Path, apply: bool, log: list):
    """automatic Straub tail (straub_auto.py, confirmed by the user on 2026-11-21) on every scoring copy: per-frame tail
    angle, Straub = > 25 deg from the body axis (or raised over the back), 0.5-s majority -> analysis/straub_auto/day<N>/"""
    out = P / "analysis" / "straub_auto" / f"day{day}"
    for v in sorted(sdir.glob("*_first5min.mp4")) if sdir.exists() else []:
        js = out / f"{v.stem}_summary.json"
        if js.exists() and js.stat().st_mtime >= v.stat().st_mtime:
            continue
        if not apply:
            log.append(f"  would run the automatic Straub tail on {v.name}")
            continue
        import straub_auto_v8 as SAU      # v8 since 2026-11-22 (user: correct body line for a mouse seen from behind,
        S = SAU.run(v, out, None)         # vertical body = rearing / grooming not counted; 5 frames answered)
        log.append(f"  automatic Straub tail {v.name}: {S['straub_pct']} % (angle measured in {S['measured_pct']} % of "
                   f"frames) -> {js.relative_to(P).as_posix()}")


def table(day: int, mice, log: list):
    import numpy as np
    import scipy.io as sio
    tail = {}
    tf = P / "assays" / "tail_immersion_all.csv"
    if tf.exists():
        with open(tf, newline="", encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                if int(r["day"]) == day:
                    tail[(r["cage"], r["mouse"])] = r["latency_s"]
    autos = {}
    for js in sorted((P / "analysis" / "straub_auto" / f"day{day}").glob("*_summary.json")):
        _, ac, am, ad = parse(js.stem, mice)
        if ac and am and ad == day:
            autos[(ac, am)] = json.loads(js.read_text(encoding="utf-8"))
    rows = []
    for cage in sorted(mice):
        for mouse in sorted(mice[cage]):
            r = dict(day=day, cage=cage, mouse=mouse)
            au = autos.get((cage, mouse))
            if au:
                r["auto_straub_pct"] = au.get("straub_pct")
                r["auto_measured_pct"] = au.get("measured_pct")
                r["auto_vertical_not_counted_pct"] = au.get("excluded_vertical_pct")
            adir = P / "assays" / f"day{day}" / cage / mouse
            f = sorted(adir.glob("ManualScoringResults*straub*.mat")) if adir.exists() else []
            r["straub_file"] = f[-1].name if f else ""
            if f:
                S = sio.loadmat(f[-1], squeeze_me=True)
                sc = np.asarray(S["score"]).ravel()
                n = int(np.asarray(S.get("numFrames", len(sc))).ravel()[0])
                sc = sc[:n]
                fr = float(np.asarray(S.get("frameRate", 30.0)).ravel()[0])
                r["seconds_scored"] = round(len(sc) / fr, 1)
                r["straub_pct"] = round(100 * float(np.mean(sc == KEY_STRAUB)), 2)
                r["straub_free_pct"] = round(100 * float(np.mean(sc == 0)), 2)
                other = 100 * float(np.mean((sc != 0) & (sc != KEY_STRAUB)))
                r["other_keys_pct"] = round(other, 2)
                if other > 0:
                    log.append(f"  !! {cage} {mouse}: {other:.2f} % of frames have a key other than 'a' - check")
                if len(sc) / fr < 290:
                    log.append(f"  !! {cage} {mouse}: only {len(sc) / fr:.0f} s scored (< 5 min)")
            r["tail_immersion_s"] = tail.get((cage, mouse), "")
            rows.append(r)
    out = P / "analysis" / "days" / f"day{day}" / f"straub_day{day}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = []
    for r in rows:
        cols += [c for c in r if c not in cols]
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return out, rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--day", type=int, required=True)
    ap.add_argument("--src", nargs="*", help="folders with the day's videos (default: assay_video_dirs)")
    ap.add_argument("--seconds", type=float, default=300.0, help="length of the scoring copy (default 300 s = 5 min)")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--no-auto", action="store_true", help="skip the automatic Straub tail (straub_auto.py)")
    a = ap.parse_args(argv)
    AD.idle()
    mice = AD.known_mice(P)
    log = []
    print(f"=== STRAUB TAIL day {a.day} (first {a.seconds:g} s, scored by hand: key 'a' held while the tail is up) ===")
    arrange(a.day, [Path(s) for s in a.src] if a.src else AD.default_src(), mice, a.apply, log)
    vids, sdir = copies(a.day, a.seconds, a.apply, log)
    file_scores(a.day, sdir, mice, a.apply, log)
    if not a.no_auto:
        auto_step(a.day, sdir, a.apply, log)
    if vids:
        mdir = T / "opioidaddiction-matlab" / "current_cohort_workflow" / "matlab"
        n_done = len(list((sdir / "output").glob("ManualScoringResults_*.mat"))) if (sdir / "output").exists() else 0
        log.append(f"  {len(vids)} video(s), {n_done} scored. To score at real time (hold 'a' while the tail is up): "
                   f"in MATLAB cd('{mdir.as_posix()}'); manual_scoring_video_multibehavior_speed('{sdir.as_posix()}', 1)")
    else:
        log.append(f"  no Straub video for day {a.day} yet (name needs 'straub', the mouse and 'day{a.day}', e.g. "
                   f"9670fgreen1_straub_day{a.day}.mp4)")
    if a.apply and vids:                       # no table on days without Straub videos
        out, rows = table(a.day, mice, log)
        log.append(f"  table: {out.relative_to(P).as_posix()} ({sum(1 for r in rows if r['straub_file'])} scored)")
    print("\n".join(log) + ("" if a.apply else "\n(dry run - nothing written)"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
