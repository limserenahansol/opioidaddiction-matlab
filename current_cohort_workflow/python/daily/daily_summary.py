#!/usr/bin/env python3
"""
daily_summary.py - end-of-day routine ("summary of today behavior result").

  1) ARRANGE  find the day's Python logs, Saleae captures (.sal + digital.csv) and camera videos wherever
              they were saved, pair them by START TIME (not by file name), and put them into
              <project>/python|salae|video/day<N>/<cage>/<mouse>/   (cage = e.g. 9670f, mouse = e.g. green1)
              Python logs are copied (the originals stay), Saleae and video files are moved.
              Every action is written to <project>/_arrangement_log.csv. Files still open (recording) are skipped.
  2) SYNC     sync_sessions.py   (Python + Saleae + video on one clock)
  3) ANALYSIS quick_analysis.py  (per-day summary figures + across-day master tables)
  4) PUPIL    the user's pupil/TST package for this day (TST too when tst_video\\dayN has videos); ~10 min/session
  5) STEP 02  copy this day's concat_out + task log (+ assay files) into step02_base (input of the original Step 02)
  6) (--longitudinal) original Step 02 -> ALL_mice_longitudinal.csv, then the cohort figures

Usage:
  python daily_summary.py                    today (all steps, pupil included)
  python daily_summary.py --date 2026-11-13  another day
  python daily_summary.py --dry-run          only show what would be moved
  python daily_summary.py --no-arrange       only sync + analysis (+ pupil)
  python daily_summary.py --no-pupil         skip the slow pupil step today (run it later; done sessions are skipped)
  python daily_summary.py --longitudinal     also rebuild ALL_mice_longitudinal.csv + cohort figures
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
import zipfile
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
PROJECT = TOOLS.parent
HOME = Path.home()
# where Python / Logic / IC Capture may have saved things (searched for files of the given day)
SEARCH = [
    (HOME / "Downloads" / "behavior_task_windows", 1),
    (HOME / "Downloads", 2),
    (HOME / "Documents", 3),
    (HOME / "Desktop", 3),
    (HOME / "Videos", 3),
    (HOME, 2),
    (PROJECT, 6),
]
SKIP_PARTS = {"AppData", ".claude", ".venv", "sync", "analysis", "_tests_and_aborted", "node_modules", "tools",
              "pupil", "$Recycle.Bin", "OneDrive"}
NAME_RE = re.compile(r"^(?P<cage>\d{4}[mfMF]?)_*(?P<mouse>[A-Za-z]+\d*)_*(?P<task>fr\d+|pr\d*[a-z0-9]*)?_*day(?P<day>\d+)",
                     re.I)
MIN_SESSION_S = 120          # shorter Python runs are tests / aborted
DEEP_MAX = 25                # max folder depth followed for folders changed on the searched day
ALIAS_FILE = PROJECT / "session_aliases.csv"   # optional: session,cage,mouse (exact Python log name -> folder)
CAM_HZ = 30.26


# ----------------------------------------------------------------------------- small helpers
def in_use(p: Path) -> bool:
    """True if another program (IC Capture, Logic, Python) still has the file open."""
    if sys.platform != "win32":
        return False
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.CreateFileW.restype = ctypes.c_void_p
    h = k32.CreateFileW(str(p), 0x80000000 | 0x40000000, 0, None, 3, 0x80, None)
    if h is None or h == ctypes.c_void_p(-1).value:
        return True
    k32.CloseHandle(ctypes.c_void_p(h))
    return False


def ffmpeg_exe():
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:
        return shutil.which("ffmpeg")


def scan(day: date):
    """All candidate files touched on that day (local time)."""
    t0 = datetime.combine(day, datetime.min.time()).timestamp()
    t1 = t0 + 86400 + 6 * 3600          # sessions that run past midnight
    found, seen = {"jsonl": [], "sal": [], "csv": [], "video": []}, set()
    for root, depth in SEARCH:
        if not root.is_dir():
            continue
        stack = [(root, 0)]
        while stack:
            d, lvl = stack.pop()
            try:
                entries = list(d.iterdir())
            except OSError:
                continue
            for p in entries:
                if p.name in SKIP_PARTS or p.name.startswith("."):
                    continue
                if p.is_dir():
                    if lvl + 1 < depth:
                        stack.append((p, lvl + 1))
                    elif lvl + 1 < DEEP_MAX:
                        # Logic 2 often saves each capture inside the previous export folder (a chain 8+ levels
                        # deep): follow folders that were created / changed on that day, however deep
                        try:
                            if t0 <= p.stat().st_mtime <= t1:
                                stack.append((p, lvl + 1))
                        except OSError:
                            pass
                    continue
                key = str(p).lower()
                if key in seen:
                    continue
                try:
                    mt = p.stat().st_mtime
                except OSError:
                    continue
                if not (t0 <= mt <= t1):
                    continue
                ext = p.suffix.lower()
                kind = {".jsonl": "jsonl", ".sal": "sal", ".csv": "csv", ".mp4": "video", ".avi": "video"}.get(ext)
                if kind is None and p.suffix == "" and zipfile.is_zipfile(p):
                    kind = "sal"                               # Logic capture saved without extension
                if kind:
                    seen.add(key)
                    found[kind].append(p)
    return found


def read_python(p: Path):
    try:
        with open(p, encoding="utf-8-sig", errors="replace") as f:
            lines = [ln for ln in f if ln.strip()]
        first = json.loads(lines[0])
        last = None
        for ln in reversed(lines):
            try:
                last = json.loads(ln); break
            except json.JSONDecodeError:
                continue
        if first.get("tag") != "meta" or "t" not in first:
            return None
        start = datetime.fromisoformat(first["t"])
        return dict(path=p, name=p.stem, start=start.timestamp(), start_str=first["t"][:19],
                    dur=float(last.get("elapsed", 0)), ended=last.get("tag") == "experiment_end",
                    rewards=sum('"reward_cmd"' in ln for ln in lines), paradigm=first.get("paradigm", ""))
    except Exception:
        return None


def read_sal(p: Path):
    try:
        with zipfile.ZipFile(p) as z:
            m = json.loads(z.read("meta.json"))
        st = m["data"]["captureStartTime"]
        return dict(path=p, start=(float(st["unixTimeMilliseconds"]) + float(st.get("fractionalMilliseconds", 0))) / 1000)
    except Exception:
        return None


def read_csv(p: Path):
    """Saleae digital export: first camera rising edge and last time (relative to capture start)."""
    try:
        with open(p, encoding="utf-8-sig", errors="replace") as f:
            hdr = [h.strip().lower() for h in f.readline().split(",")]
            if not hdr or not hdr[0].startswith("time") or "camera" not in hdr:
                return None
            ci = hdr.index("camera")
            ri = hdr.index("reward") if "reward" in hdr else None
            prev, first_cam, last_t, n_cam = "0", None, None, 0
            rprev, r_fall, rises = "0", -1.0, []
            for ln in f:
                parts = ln.split(",")
                if len(parts) <= ci:
                    continue
                t, v = float(parts[0]), parts[ci].strip()
                last_t = t
                if v == "1" and prev != "1":
                    n_cam += 1
                    if first_cam is None:
                        first_cam = t
                prev = v
                if ri is not None and len(parts) > ri:
                    rv = parts[ri].strip()
                    if rv == "1" and rprev != "1" and t - r_fall > 0.001:     # ignore valve contact bounce
                        rises.append(t)
                    if rv == "0" and rprev == "1":
                        r_fall = t
                    rprev = rv
        return dict(path=p, first_cam=first_cam, last=last_t, n_cam=n_cam, rewards=rises)
    except Exception:
        return None


def read_video(p: Path, ff):
    out = dict(path=p, created=None)
    if not ff:
        return out
    try:
        head = subprocess.run([ff, "-hide_banner", "-i", str(p)], capture_output=True, text=True,
                              errors="replace", timeout=60).stderr
        m = re.search(r"creation_time\s*:\s*(\S+)", head)
        if m:
            out["created"] = datetime.fromisoformat(m.group(1).replace("Z", "+00:00")).timestamp()
        out["readable"] = "moov atom not found" not in head and "Invalid data" not in head
    except Exception:
        pass
    return out


# ----------------------------------------------------------------------------- pairing
def known_mice():
    """{cage: set(mouse)} from the folders already filed in python/dayN/<cage>/<mouse>/."""
    out = {}
    for d in (PROJECT / "python").glob("day*/*/*"):
        if d.is_dir():
            out.setdefault(d.parent.name.lower(), set()).add(d.name.lower())
    return out


def resolve_mouse_names(real, log):
    """Make the mouse folder the same as on earlier days.
    1) session_aliases.csv (session,cage,mouse) wins.  2) If the name is not a known mouse of that cage, but exactly
    one known mouse of the cage starts with it and has no session yet that day (e.g. 'green' -> 'green1' when
    'green2' already ran), use that one and print a CHECK line.  Otherwise keep the name and warn."""
    alias = {}
    if ALIAS_FILE.exists():
        import csv
        with open(ALIAS_FILE, encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                if r.get("session"):
                    alias[r["session"].strip().lower()] = (r["cage"].strip().lower(), r["mouse"].strip().lower())
    known = known_mice()
    for s in real:
        key = s["py"]["name"].lower()
        if key in alias:
            s["cage"], s["mouse"] = alias[key]
            s["problems"].append(f"folder set by session_aliases.csv -> {s['cage']}/{s['mouse']}")
    for s in real:
        cage, mouse = s["cage"], s["mouse"]
        km = known.get(cage, set())
        if not km or mouse in km:
            continue
        today = {x["mouse"] for x in real if x["cage"] == cage and x is not s and x["day"] == s["day"]}
        cand = sorted(k for k in km if k.startswith(mouse) and k not in today)
        if len(cand) == 1:
            s["problems"].append(f"mouse name '{mouse}' is not used before in {cage}; filed as '{cand[0]}' "
                                 f"(the only {cage} mouse starting with '{mouse}' without a session that day). "
                                 f"If wrong, add a line to session_aliases.csv")
            s["mouse"] = cand[0]
        else:
            s["problems"].append(f"new mouse name '{mouse}' in {cage} (known: {', '.join(sorted(km))}) - "
                                 f"if it is one of these, add a line to session_aliases.csv")


def pair_day(found, ff, log):
    pys = [x for x in (read_python(p) for p in found["jsonl"]) if x]
    sals = [x for x in (read_sal(p) for p in found["sal"]) if x]
    csvs = [x for x in (read_csv(p) for p in found["csv"]) if x]
    vids = [read_video(p, ff) for p in found["video"] if p.stat().st_size > 1_000_000]
    sessions, used = [], set()
    # the same run can exist several times (Downloads copy, renamed copy, project copy): one per START time,
    # preferring the copy already in the project; runs already filed under _tests_and_aborted are ignored
    known_tests = set()
    for t in (PROJECT / "_tests_and_aborted" / "python").glob("*.jsonl"):
        x = read_python(t)
        if x:
            known_tests.add(round(x["start"]))
    uniq = {}
    for p in sorted(pys, key=lambda x: (PROJECT not in x["path"].parents, len(x["name"]))):
        k = round(p["start"])
        if k in known_tests and PROJECT not in p["path"].parents:
            continue
        uniq.setdefault(k, p)
    for p in sorted(uniq.values(), key=lambda x: x["start"]):
        m = NAME_RE.match(p["name"])
        s = dict(py=p, sal=None, csv=None, video=None, problems=[])
        if not p["ended"] and (datetime.now().timestamp() - p["path"].stat().st_mtime) < 300:
            s["status"] = "RUNNING (skipped)"
        elif p["name"].lower().startswith("test") or p["dur"] < MIN_SESSION_S or not m:
            s["status"] = "test/aborted" if m or p["name"].lower().startswith("test") or p["dur"] < MIN_SESSION_S \
                else "name not understood"
        else:
            s["status"] = "session"
            s["cage"], s["mouse"] = m.group("cage").lower(), m.group("mouse").lower()
            s["day"] = int(m.group("day"))
        sessions.append(s)
    real = [s for s in sessions if s["status"] == "session"]
    resolve_mouse_names(real, log)
    # Saleae capture: started 0-10 min before Python
    for s in real:
        c = [(s["py"]["start"] - x["start"], x) for x in sals if id(x) not in used and -5 <= s["py"]["start"] - x["start"] <= 600]
        if c:
            lag, x = min(c, key=lambda z: z[0]); s["sal"] = x; used.add(id(x)); s["sal_lag"] = lag
    # Saleae CSV: (1) same number of valve openings during the session as Python reward commands,
    # (2) camera pulses for the session length, (3) first camera pulse at (python start - capture start) if the
    # .sal is known, (4) prefer the CSV saved in the same folder as that .sal
    for s in real:
        best = None
        n_rew, dur = s["py"]["rewards"], s["py"]["dur"]
        for x in csvs:
            if id(x) in used or x["first_cam"] is None:
                continue
            in_sess = [r for r in x["rewards"] if x["first_cam"] - 0.1 <= r <= x["first_cam"] + dur + 2]
            if abs(len(in_sess) - n_rew) > max(2, 0.03 * n_rew):
                continue
            if abs(x["n_cam"] - dur * CAM_HZ) > max(60, 0.03 * dur * CAM_HZ):
                continue
            if s["sal"] is not None:
                err = abs(x["first_cam"] - s["sal_lag"])
                if err > 1.0:
                    continue
                score = (x["path"].parent != s["sal"]["path"].parent, abs(len(in_sess) - n_rew), err)
            else:
                py_end = s["py"]["start"] + dur
                gap = x["path"].stat().st_mtime - py_end
                if not 0 <= gap <= 3600:
                    continue
                score = (False, abs(len(in_sess) - n_rew), gap)
            if best is None or score < best[0]:
                best = (score, x)
        if best:
            s["csv"] = best[1]; used.add(id(best[1]))
            if s["sal"] is None:
                s["problems"].append("no .sal found - CSV matched by reward count, pulse count and save time")
    # video: IC Capture stamps the file with the first frame time (= Python start, whole seconds)
    for s in real:
        c = [(abs(v["created"] - s["py"]["start"]), v) for v in vids
             if id(v) not in used and v.get("created") is not None and -1.5 <= v["created"] - s["py"]["start"] <= 3]
        if not c:
            py_end = s["py"]["start"] + s["py"]["dur"]
            c = [(abs(v["path"].stat().st_mtime - py_end), v) for v in vids if id(v) not in used and
                 v.get("created") is None and v.get("readable") and 0 <= v["path"].stat().st_mtime - py_end <= 900]
            if c:
                s["problems"].append("video matched by save time (no creation time)")
        if c:
            _, v = min(c, key=lambda z: z[0]); s["video"] = v; used.add(id(v))
        else:
            unfinished = [v for v in vids if not v.get("readable", True) and id(v) not in used]
            if unfinished:
                s["problems"].append("an unfinished/unreadable video exists - stop the IC Capture recording")
    for s in real:
        for k in ("sal", "csv", "video"):
            if s[k] is None:
                s["problems"].append(f"no {k} found")
    return sessions


# ----------------------------------------------------------------------------- moving
def target_dir(kind, s):
    return PROJECT / kind / f"day{s['day']}" / s["cage"] / s["mouse"]


def place(action, src: Path, dst: Path, dry, log):
    if src.resolve() == dst.resolve():
        return "already in place"
    if dst.exists():
        if dst.stat().st_size == src.stat().st_size:
            return "already there (same size)"
        k = 2
        while dst.with_name(f"{dst.stem}_{k}{dst.suffix}").exists():
            k += 1
        dst = dst.with_name(f"{dst.stem}_{k}{dst.suffix}")
    if in_use(src):
        return "IN USE (still recording/open) - left in place"
    if dry:
        return f"would {action} -> {dst}"
    dst.parent.mkdir(parents=True, exist_ok=True)
    if action == "copy":
        shutil.copy2(src, dst)
    else:
        shutil.move(str(src), str(dst))
    log.append([datetime.now().isoformat(timespec="seconds"), action, str(src), str(dst), "ok (daily_summary)"])
    return f"{action} -> {dst.relative_to(PROJECT)}"


def arrange(day: date, dry: bool):
    ff = ffmpeg_exe()
    found = scan(day)
    log = []
    sessions = pair_day(found, ff, log)
    print(f"\n=== ARRANGE {day} : {len(sessions)} Python log(s) found ===")
    rows = []
    for s in sessions:
        p = s["py"]
        if s["status"] != "session":
            res = ""
            if s["status"] == "test/aborted":
                res = place("copy", p["path"], PROJECT / "_tests_and_aborted" / "python" / p["path"].name, dry, log)
            print(f"  {p['name']:30s} {p['start_str'][11:]}  {p['dur']/60:5.1f} min  {s['status']}  {res}")
            rows.append(dict(session=p["name"], start=p["start_str"], status=s["status"], result=res))
            continue
        res = [place("copy", p["path"], target_dir("python", s) / f"{p['name']}.jsonl", dry, log)]
        if s["sal"]:
            res.append(place("move", s["sal"]["path"], target_dir("salae", s) / f"{p['name']}.sal", dry, log))
        if s["csv"]:
            res.append(place("move", s["csv"]["path"], target_dir("salae", s) / "digital.csv", dry, log))
        if s["video"]:
            v = s["video"]["path"]
            res.append(place("move", v, target_dir("video", s) / f"{p['name']}{v.suffix.lower()}", dry, log))
        flag = "OK" if not s["problems"] else "CHECK: " + "; ".join(s["problems"])
        print(f"  {p['name']:30s} {p['start_str'][11:]}  {p['dur']/60:5.1f} min  day{s['day']} {s['cage']}/{s['mouse']}  {flag}")
        for r in res:
            if not r.startswith("already"):
                print(f"      {r}")
        rows.append(dict(session=p["name"], start=p["start_str"], status=flag, result=" | ".join(res)))
    if log and not dry:
        lp = PROJECT / "_arrangement_log.csv"
        new = not lp.exists()
        with open(lp, "a", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            if new:
                w.writerow(["when", "action", "from", "to", "result"])
            w.writerows(log)
    # clean up empty folders that the moves left behind (only outside the project's data folders)
    if not dry:
        for x in [c["path"].parent for c in [*(s["csv"] for s in sessions if s.get("csv"))]]:
            try:
                if x != HOME and PROJECT not in x.parents and not any(x.iterdir()):
                    x.rmdir()
            except OSError:
                pass
    return rows


PACKAGE = TOOLS / "auto_pupil_tst_pipeline_20261114_iris_on"     # the user's pupil + TST package (RUN_ALL_18DAYS)


WORKFLOW = TOOLS / "opioidaddiction-matlab" / "current_cohort_workflow"   # git repo copy (Step 02 wrapper etc.)


def matlab_exe():
    return shutil.which("matlab") or next((str(p) for p in sorted(Path(r"C:\Program Files\MATLAB").glob("R*/bin/matlab.exe"))), None)


def run_package_day(day_n: int, tst: bool) -> bool:
    """Pupil (and TST when that day has TST videos) for ONE day with the user's pupil/TST package. Runs every day
    by default since 2026-11-14 (finished sessions are skipped, so the day-18 RUN_ALL_18DAYS.bat run only has to do
    the QC page and the sessions marked redo). Writes a one-day copy of run_all_config.json (days = [N]; TST videos
    of tst_video\\dayN only, results in pupil_tst_final\\TST_auto\\dayN) and runs MATLAB -batch (ROI: previous day's)."""
    cfg = PACKAGE / "run_all_config.json"
    if not cfg.exists():
        print("\n(pupil/TST step skipped: package or run_all_config.json missing)")
        return False
    c = json.loads(cfg.read_text(encoding="utf-8"))
    c.update(days=[day_n], dry_run=False, run_pupil=True, run_tst=bool(tst), roi_mode="previous",
             _note=f"one-day run from daily_summary.py (day {day_n}); the main config is run_all_config.json")
    if tst:
        c.update(tst_video_root=str(PROJECT / "tst_video" / f"day{day_n}"),
                 tst_out_root=str(Path(c.get("pupil_out_root") or PROJECT / "pupil_tst_final") / "TST_auto" / f"day{day_n}"))
    day_cfg = PACKAGE / f"run_all_config_day{day_n}.json"
    day_cfg.write_text(json.dumps(c, indent=2), encoding="utf-8")
    matlab = matlab_exe()
    if not matlab:
        print("\n(pupil/TST step skipped: MATLAB not found)")
        return False
    print(f"\n=== PUPIL{' + TST' if tst else ''} (package, day {day_n}; about 10 min per session) ===")
    cmd = f"cd('{PACKAGE}'); RUN_ALL_18DAYS('{day_cfg}')"
    subprocess.run([matlab, "-batch", cmd])
    return True


def update_step02_base(day_n: int) -> bool:
    """Copy this day's newest concat_out + task log (+ assay files) into <project>\\step02_base (what the ORIGINAL
    Step 02 reads), so the longitudinal table can be built at any time with run_step02_on_base.m."""
    script = WORKFLOW / "python" / "prepare_step02_base.py"
    if not script.exists():
        print("\n(step02_base not updated: current_cohort_workflow not found in tools)")
        return False
    print(f"\n=== STEP 02 BASE (day {day_n}) ===")
    subprocess.run([sys.executable, str(script), "--pupil-out", str(PROJECT / "pupil_tst_final"), "--python-root",
                    str(PROJECT / "python"), "--assays", str(PROJECT / "assays"), "--base", str(PROJECT / "step02_base"),
                    "--days", str(day_n), "--apply"])
    return True


def run_longitudinal() -> bool:
    """ORIGINAL Step 02 on <project>\\step02_base (new longitudinal_outputs\\run_###), then the cohort figures."""
    matlab = matlab_exe()
    base = PROJECT / "step02_base"
    if not matlab or not base.exists():
        print("\n(longitudinal step skipped: MATLAB or step02_base missing)")
        return False
    print("\n=== STEP 02 (original, via run_step02_on_base) ===")
    subprocess.run([matlab, "-batch", f"cd('{WORKFLOW / 'matlab'}'); run_step02_on_base('{base}');"])
    runs = sorted((base / "longitudinal_outputs").glob("run_*"))
    if not runs:
        return False
    csvp = runs[-1] / "ALL_mice_longitudinal.csv"
    man = PROJECT / "cohort_manifest.csv"
    print(f"\n=== COHORT FIGURES ({runs[-1].name}) ===")
    cmd = [sys.executable, str(WORKFLOW / "python" / "cohort_plots.py"), "--longitudinal", str(csvp),
           "--out", str(PROJECT / "analysis" / "cohort_figs" / runs[-1].name)]
    if man.exists():
        cmd += ["--manifest", str(man)]
    subprocess.run(cmd)
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--date", help="YYYY-MM-DD (default: today)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-arrange", action="store_true")
    ap.add_argument("--no-pupil", action="store_true",
                    help="skip the pupil/TST package today (default: run it for this day, ~10 min per session)")
    ap.add_argument("--pupil", action="store_true", help="(old flag; pupil now runs by default)")
    ap.add_argument("--longitudinal", action="store_true",
                    help="also rebuild ALL_mice_longitudinal.csv with the original Step 02 and the cohort figures")
    args = ap.parse_args(argv)
    day = date.fromisoformat(args.date) if args.date else date.today()
    if not args.no_arrange:
        arrange(day, args.dry_run)
    if args.dry_run:
        return 0
    sys.path.insert(0, str(TOOLS))
    import sync_sessions, quick_analysis
    print("\n=== SYNC ===")
    rc = sync_sessions.main(["--root", str(PROJECT)])
    print("\n=== ANALYSIS ===")
    quick_analysis.main(["--root", str(PROJECT)])
    day_folders = sorted(int(p.name[3:]) for p in (PROJECT / "python").glob("day*") if p.name[3:].isdigit()) \
        if (PROJECT / "python").exists() else []
    if day_folders:
        # the day folder of this date = the newest day folder that has a log from this date
        dn = next((d for d in reversed(day_folders)
                   if any(datetime.fromtimestamp(p.stat().st_mtime).date() == day
                          for p in (PROJECT / "python" / f"day{d}").rglob("*.jsonl"))), day_folders[-1])
        if not args.no_pupil:
            tst_dir = PROJECT / "tst_video" / f"day{dn}"
            has_tst = tst_dir.exists() and any(p.suffix.lower() in (".mp4", ".avi", ".mov", ".mkv", ".wmv", ".m4v")
                                               for p in tst_dir.rglob("*"))
            if run_package_day(dn, has_tst):
                print("\n=== ANALYSIS (again, now with pupil) ===")
                quick_analysis.main(["--root", str(PROJECT)])
        else:
            print("\n(pupil/TST: skipped today with --no-pupil; it can be run later, finished sessions are skipped)")
        update_step02_base(dn)
        if args.longitudinal:
            run_longitudinal()
    days = sorted((PROJECT / "analysis" / "days").glob("day*")) if (PROJECT / "analysis" / "days").exists() else []
    print("\nDone. Per-day summaries:", ", ".join(d.name for d in days) or "-")
    print(f"  {PROJECT / 'analysis' / 'days'}\\dayN\\figs\\00_day_dashboard.png")
    print(f"  across days: {PROJECT / 'analysis' / 'master_day_level.csv'} and analysis\\figs_across_days\\")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
