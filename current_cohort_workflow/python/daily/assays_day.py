#!/usr/bin/env python3
"""
assays_day.py - the TST and hot-plate (HOT) videos of one experiment day into the cohort's files and Step 02 (added
2026-11-18, day 5). The user: "TST should be 3min so just use first 3min", "all HOT should be 30s long", "TST to step2
you can do whatever but for us to use this consistently with previous cohort and the other days".

For day N (TST + HOT on days 5, 10, 13, 16, 18; the tail-immersion latencies go in with tail_immersion_entry.py):
 1. ARRANGE: videos in --src (default: assay_video_dirs of tools/local_paths.json, else <home>/Videos/Logitech/
    LogiCapture; 2 folder levels) whose
    name has 'tst' or 'hot', a mouse of the cohort ('<cage><sex><colour>', e.g. 9670fgreen1) and 'day<N>' are COPIED
    (the originals stay) to <project>/tst_video/day<N>/<cage>/<mouse>/ and <project>/hp_video/day<N>/<cage>/<mouse>/.
    A file written in the last 2 min (maybe still recording) waits for the next run.
 2. TST: every TST video of the day is scored by tst_mouse_lock.py (the package's scorer, not changed; first 180 s;
    per-frame mouse check; camera-view fix only when the untouched run lost the mouse) into
    <project>/analysis/tst/day<N>/<video>/. The chosen result goes to Step 02 as
      assays/day<N>/<cage>/<mouse>/ManualScoringResults_<video>_AUTO.mat
    = the file name the ORIGINAL Step 02 reads ('ManualScoringResults*tst*.mat'), so the TST_* columns are made exactly
    as for the previous cohort; '_AUTO' marks the automatic score and shows in Step 02's TST_File column. A hand score
    of that mouse-day (ManualScoringResults*tst*.mat without '_AUTO' in that folder) wins: the _AUTO file is then not
    written. Example sheets for checking: analysis/tst/day<N>/examples/ (tst_examples.py).
 3. HOT: the FIRST 30 s of every HOT video (30 s x its frame rate, e.g. 900 frames at 30 fps; no audio) as a scoring
    copy in <project>/hp_video/day<N>_first30s/ for the lab's manual scorer:
      cd('<project>/tools/opioidaddiction-matlab/step01_manual_assays')
      manual_scoring_video_multibehavior_batch2('<project>/hp_video/day<N>_first30s')
    Its ManualScoringResults_*.mat (written to that folder's output/) are copied to assays/day<N>/<cage>/<mouse>/
    (the name keeps 'HOT', which Step 02 reads as 'ManualScoringResults*hot*.mat').
 4. TABLE: analysis/days/day<N>/assays_day<N>.csv, one row per mouse of the cohort: TST % non-moving (first 3 min), run
    used, % frames that were not the mouse, HOT scored or not with % per behaviour, tail-immersion latency
    (assays/tail_immersion_all.csv).
The daily routine's step02 base update (prepare_step02_base.py --assays) then copies the assay files next to the task
logs, and the original Step 02 fills TST_*, HOT_* and Immersion_Latency_s. A file already in step02_base is not
replaced there: a changed result is printed with '!!'.
Dry run unless --apply. Runs at Windows idle priority; tst_mouse_lock.py runs one video at a time with the Python of
tools/.venv_pkg313 (scikit-learn).
usage: python assays_day.py --day 5 [--src <folder> ...] [--apply] [--no-tst] [--no-examples]
"""
from __future__ import annotations

import argparse
import csv
import ctypes
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

T = Path(__file__).resolve().parent
P = T.parent
PKG_PY = T / ".venv_pkg313" / "Scripts" / "python.exe"
VIDEO_EXT = {".mp4", ".avi", ".mov", ".mkv", ".wmv", ".m4v"}
MOUSE = re.compile(r"(?P<cage>\d{4})(?P<sex>[mf])?[\s_-]*(?P<color>black|white|green\d|orange|red|blue|yellow|pink)",
                   re.I)
def default_src() -> list:
    """assay_video_dirs of tools/local_paths.json (this PC's folders, not in the repository), else Logitech Capture's
    default folder"""
    try:
        got = json.loads((T / "local_paths.json").read_text(encoding="utf-8")).get("assay_video_dirs", [])
    except (OSError, ValueError):
        got = []
    return [Path(p) for p in got] or [Path.home() / "Videos" / "Logitech" / "LogiCapture"]
HOT_S = 30.0


def idle():
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x40)
    except Exception:  # noqa: BLE001
        pass


def known_mice(project: Path) -> dict:
    out = {}
    for p in (project / "python").glob("day*/*/*"):
        if p.is_dir():
            out.setdefault(p.parent.name.lower(), set()).add(p.name.lower())
    return out


def parse(name: str, mice: dict):
    """(assay, cage folder, mouse folder, day) from a video name; None parts when not understood"""
    low = name.lower()
    assay = "tst" if "tst" in low else ("hot" if "hot" in low else None)
    d = re.search(r"day(\d+)(?!\d)", low)
    m = MOUSE.search(name)
    if m is None:
        return assay, None, None, int(d.group(1)) if d else None
    cage4, color = m.group("cage"), m.group("color").lower()
    cages = [c for c in mice if c[:4] == cage4]
    cage = cages[0] if len(cages) == 1 else None
    if cage is not None and color not in mice[cage]:
        color = None
    return assay, cage, color, int(d.group(1)) if d else None


def ffmpeg():
    c = sorted((T / ".venv" / "Lib" / "site-packages" / "imageio_ffmpeg" / "binaries").glob("ffmpeg*.exe"))
    return c[-1] if c else None


def video_fps(path: Path) -> float:
    import cv2
    cap = cv2.VideoCapture(str(path))
    f = cap.get(cv2.CAP_PROP_FPS) or 30.0
    cap.release()
    return round(f, 3)


def same_file(a: Path, b: Path) -> bool:
    return b.exists() and a.stat().st_size == b.stat().st_size


def arrange(day: int, srcs, mice, apply: bool, log: list):
    n = 0
    for src in srcs:
        if not src.is_dir():
            continue
        cands = [p for p in src.glob("*") if p.is_file()] + [p for p in src.glob("*/*") if p.is_file()]
        for v in sorted(cands):
            if v.suffix.lower() not in VIDEO_EXT:
                continue
            assay, cage, mouse, d = parse(v.stem, mice)
            if assay is None or d != day:
                continue
            if cage is None or mouse is None:
                log.append(f"  !! {v.name}: mouse not understood (known: {sorted(mice)}) - not copied")
                continue
            if time.time() - v.stat().st_mtime < 120:
                log.append(f"  {v.name}: written less than 2 min ago (still recording?) - next run")
                continue
            dst = P / ("tst_video" if assay == "tst" else "hp_video") / f"day{day}" / cage / mouse / v.name
            if same_file(v, dst):
                continue
            if dst.exists():
                log.append(f"  !! {v.name}: a different file of that name is already in {dst.parent} - not copied")
                continue
            log.append(f"  {'copy' if apply else 'would copy'} {v.name} -> {dst.relative_to(P).as_posix()}")
            if apply:
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(v, dst)
            n += 1
    return n


def day_videos(kind: str, day: int):
    root = P / ("tst_video" if kind == "tst" else "hp_video") / f"day{day}"
    return sorted(p for p in root.glob("*/*/*") if p.is_file() and p.suffix.lower() in VIDEO_EXT)


def copy_if_changed(src: Path, dst: Path, apply: bool, log: list, what: str):
    if dst.exists() and dst.read_bytes() == src.read_bytes():
        return False
    base = P / "step02_base" / dst.relative_to(P / "assays")
    if dst.exists() and base.exists():
        log.append(f"  !! {what}: {dst.name} changed, but step02_base still has the old copy ({base.parent}) - "
                   f"tell Claude")
    log.append(f"  {'write' if apply else 'would write'} {dst.relative_to(P).as_posix()} ({what})")
    if apply:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    return True


def tst_step(day: int, apply: bool, log: list, examples: bool):
    vids = day_videos("tst", day)
    out_root = P / "analysis" / "tst" / f"day{day}"
    for v in vids:
        cage, mouse, stem = v.parent.parent.name, v.parent.name, v.stem
        out = out_root / stem
        js = out / f"{stem}_lock.json"
        done = False
        if js.exists():
            try:
                done = Path(json.loads(js.read_text(encoding="utf-8"))["video"]).resolve() == v.resolve()
            except (ValueError, KeyError, OSError):
                done = False
        if not done:
            if not apply:
                log.append(f"  would score {v.name} (first 180 s) -> {out.relative_to(P).as_posix()}")
                continue
            r = subprocess.run([str(PKG_PY), "-B", str(T / "tst_mouse_lock.py"), "--video", str(v), "--out", str(out)],
                               capture_output=True, text=True)
            for line in (r.stdout or "").splitlines():
                if line.startswith("  ") or "view fix" in line or "FAILED" in line:
                    log.append("  " + line.strip())
            if r.returncode or not js.exists():
                log.append(f"  !! {v.name}: tst_mouse_lock.py failed ({(r.stderr or '').strip()[-300:]})")
                continue
        mat = out / f"AutoScoringResults_{stem}.mat"
        if not mat.exists():
            log.append(f"  !! {v.name}: no AutoScoringResults (QC error?) - see {js}")
            continue
        adir = P / "assays" / f"day{day}" / cage / mouse
        hand = [f for f in adir.glob("ManualScoringResults*tst*.mat") if "_auto" not in f.name.lower()] \
            if adir.exists() else []
        if hand:
            log.append(f"  {cage} {mouse}: hand TST score {hand[0].name} is used; the automatic one is not written")
            continue
        copy_if_changed(mat, adir / f"ManualScoringResults_{stem}_AUTO.mat", apply, log, "TST automatic, first 180 s")
    if examples and apply and vids:
        subprocess.run([str(PKG_PY), "-B", str(T / "tst_examples.py"), "--lock-root", str(out_root), "--out",
                        str(out_root / "examples")], capture_output=True, text=True)
        log.append(f"  example sheets: {(out_root / 'examples').relative_to(P).as_posix()}")


def hot_step(day: int, apply: bool, log: list):
    vids = day_videos("hot", day)
    sdir = P / "hp_video" / f"day{day}_first30s"
    ff = ffmpeg()
    for v in vids:
        dst = sdir / f"{v.stem}_first30s.mp4"
        if dst.exists():
            continue
        if ff is None:
            log.append("  !! ffmpeg (tools/.venv imageio-ffmpeg) not found - no HOT scoring copies")
            break
        fps = video_fps(v)
        nfr = int(round(HOT_S * fps))
        log.append(f"  {'make' if apply else 'would make'} {dst.relative_to(P).as_posix()} (first {HOT_S:g} s = "
                   f"{nfr} frames at {fps:g} fps)")
        if apply:
            sdir.mkdir(parents=True, exist_ok=True)
            subprocess.run([str(ff), "-hide_banner", "-loglevel", "error", "-y", "-i", str(v), "-frames:v", str(nfr),
                            "-vf", f"setpts=N/{fps:g}/TB", "-r", f"{fps:g}", "-an", "-c:v", "libx264", "-crf", "16",
                            "-preset", "veryfast", "-pix_fmt", "yuv420p", str(dst)],
                           creationflags=0x40 if sys.platform == "win32" else 0)
    mice = known_mice(P)
    for m in sorted((sdir / "output").glob("ManualScoringResults_*.mat")) if (sdir / "output").exists() else []:
        assay, cage, mouse, d = parse(m.stem, mice)
        if cage is None or mouse is None or d != day:
            log.append(f"  !! {m.name}: mouse/day not understood - not copied")
            continue
        copy_if_changed(m, P / "assays" / f"day{day}" / cage / mouse / m.name, apply, log, "HOT hand score, first 30 s")
    if vids:
        n_done = len(list((sdir / "output").glob("ManualScoringResults_*.mat"))) if (sdir / "output").exists() else 0
        log.append(f"  HOT: {len(vids)} video(s), {n_done} scored by hand. To score: in MATLAB "
                   f"cd('{(T / 'opioidaddiction-matlab' / 'step01_manual_assays').as_posix()}'); "
                   f"manual_scoring_video_multibehavior_batch2('{sdir.as_posix()}')")


def table(day: int):
    import numpy as np
    import scipy.io as sio
    mice = known_mice(P)
    tail = {}
    tf = P / "assays" / "tail_immersion_all.csv"
    if tf.exists():
        with open(tf, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if int(r["day"]) == day:
                    tail[(r["cage"], r["mouse"])] = r["latency_s"]
    rows = []
    for cage in sorted(mice):
        for mouse in sorted(mice[cage]):
            r = dict(day=day, cage=cage, mouse=mouse)
            adir = P / "assays" / f"day{day}" / cage / mouse
            tst = sorted(adir.glob("ManualScoringResults*tst*.mat")) if adir.exists() else []
            hand = [f for f in tst if "_auto" not in f.name.lower()]
            used = (hand or tst or [None])[0]
            r["tst_file"] = used.name if used else ""
            r["tst_source"] = "" if used is None else ("hand" if used in hand else "automatic, first 180 s")
            if used is not None:
                S = sio.loadmat(used, squeeze_me=True)
                sc = np.asarray(S["score"]).ravel()
                n = int(np.asarray(S.get("numFrames", len(sc))).ravel()[0])
                r["tst_pct_nonmoving"] = round(100 * float(np.mean(sc[:n] == 0)), 2)
                stem = used.stem.replace("ManualScoringResults_", "").replace("_AUTO", "")
                js = P / "analysis" / "tst" / f"day{day}" / stem / f"{stem}_lock.json"
                if js.exists() and used not in hand:
                    info = json.loads(js.read_text(encoding="utf-8"))
                    r["tst_run"] = info.get("chosen")
                    r["tst_not_mouse_pct"] = info.get("mouse_qc", {}).get("not_mouse_pct")
            hot = [f for f in adir.glob("ManualScoringResults*hot*.mat")] if adir.exists() else []
            r["hot_file"] = hot[0].name if hot else ""
            if hot:
                S = sio.loadmat(hot[0], squeeze_me=True)
                sc = np.asarray(S["score"]).ravel()
                names = [str(x) for x in np.atleast_1d(S["behaviorNames"])]
                for i, nm in enumerate(names):
                    r[f"hot_pct_{nm.replace(' ', '_')}"] = round(100 * float(np.mean(sc == i)), 2)
                r["hot_seconds"] = round(len(sc) / float(S.get("frameRate", 30.0)), 2)
            r["tail_immersion_s"] = tail.get((cage, mouse), "")
            rows.append(r)
    out = P / "analysis" / "days" / f"day{day}" / f"assays_day{day}.csv"
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
    ap.add_argument("--src", nargs="*", default=None)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--no-tst", action="store_true")
    ap.add_argument("--no-examples", action="store_true")
    a = ap.parse_args(argv)
    idle()
    mice = known_mice(P)
    srcs = [Path(s) for s in a.src] if a.src else default_src()
    log = []
    print(f"=== ASSAYS day {a.day} (TST first 180 s automatic, HOT first 30 s by hand) ===")
    arrange(a.day, srcs, mice, a.apply, log)
    if not a.no_tst:
        tst_step(a.day, a.apply, log, not a.no_examples)
    hot_step(a.day, a.apply, log)
    for line in log:
        print(line)
    if a.apply:
        out, rows = table(a.day)
        n_tst = sum(1 for r in rows if r.get("tst_file"))
        n_hot = sum(1 for r in rows if r.get("hot_file"))
        print(f"  table: {out.relative_to(P).as_posix()} ({n_tst} TST, {n_hot} HOT scored)")
    else:
        print("  dry run - nothing written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
