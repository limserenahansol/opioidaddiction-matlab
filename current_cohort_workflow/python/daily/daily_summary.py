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
  4) PUPIL    first (2026-11-20) the eye-box size check: every day-N box must have its mouse's fixed size from
              video\\eye_box_sizes.csv (eye_box_fixed.py --check; report only; set the boxes with eye_box_fixed.py
              --set); then the user's pupil/TST package for this day (pupil only since 2026-11-18; TST: step 5b);
              ~10 min/session, then the options for rare sessions - the alignment fixes of alignment_fixes.csv
              (fix_concat_alignment.py) and the full-pupil traces of fullpupil_sessions.csv (fullpupil_sessions.py;
              its screen is not run daily) - then pupil grades + alignment check
              (pupil_session_quality.py), guess traces and the eye-box check sheet; then (2026-11-18) the pupil check
              pictures of all days (pupil_check_examples.py -> analysis\\pupil_check\\: one picture per session of the
              pupil the tables use, with the frames out; Claude writes the verdicts after looking)
  5) CHECKS   (report only, added 2026-11-15) session_check.py --all for this date (protocol, valve time, licks,
              camera frames vs pulses, Saleae settings vs tools\\protocol_expected.csv) and pr_daily_checks.py
              (per-mouse PR table: analysis\\days\\dayN\\pr_checks_dayN.txt/.csv and figs\\12_pr_checks.png);
              then (2026-11-18) pr_table_days.py: the PR of every mouse on every day in one table
              (analysis\\days\\pr_table_days1-N.csv/.png), and (2026-11-20) ti_table_days.py: the tail-immersion
              latency of every mouse on every day from day 5 (analysis\\days\\ti_table_days5-N.csv/.png + plot)
  5b) ASSAYS (added 2026-11-18; --no-assays skips it) assays_day.py: the day's TST / hot-plate videos (Logitech
              Capture folder, names <cage><mouse>_TST_dayN / _HOT_dayN) copied to tst_video / hp_video; TST scored
              automatically on the FIRST 3 MIN (tst_mouse_lock.py: the package scorer unchanged + per-frame mouse
              check + camera-view fix only when needed) -> assays\\dayN\\<cage>\\<mouse>\\ManualScoringResults_
              <video>_AUTO.mat; HOT: first-30-s copies for the hand scorer (hp_video\\dayN_first30s), its .mat files
              filed; table analysis\\days\\dayN\\assays_dayN.csv (TST, HOT, tail immersion)
  5c) STRAUB TAIL (added 2026-11-21; skipped with --no-assays) straub_day.py: <cage><mouse>_straub_dayN videos ->
              straub_video, first-5-min copies for the hand scorer (key 'a' held while the tail is up), its .mat files
              filed for Step 02, and the AUTOMATIC score of every copy (straub_auto.py) -> analysis\\straub_auto\\dayN;
              table analysis\\days\\dayN\\straub_dayN.csv (hand % next to automatic %)
  6) STEP 02  copy this day's concat_out + task log (+ assay files) into step02_base (input of the original Step 02)
  7) STEP 02 (every day since 2026-11-17; --no-longitudinal skips it) original Step 02 -> ALL_mice_longitudinal.csv,
              then the previous cohort's figures (cohort_plots.py; whole-eye mice left out of the pupil figures)
  8) QC SUMMARY (added 2026-11-17) qc_summary.py: one picture of all days so far (rewards Python = Saleae, PR, valve,
              session length, licks Saleae vs Python, frames vs pulses, pupil / whole eye); every CHECK item is
              explained (yellow) or red; open issues in analysis\\days\\qc_open_issues.csv
  9) ROUTINE PPT (added 2026-11-17) routine_ppt.py -> analysis\\ppt\\routine_dayN_<date>.pptx: overview, QC per
              session, data problems and fixes, routine problems fixed (analysis\\days\\routine_log.csv), where all
              plots are, basic plots (previous-cohort figures + today's figures), what to do better next
              (analysis\\days\\next_steps.csv)
 10) DATA INDEX (added 2026-11-18) data_index.py -> analysis\\data_index.csv + analysis\\days\\dayN\\data_index_dayN.csv:
              every data file saved so far (relative path, rows, size, which one is final), so the analysis after
              day 18 reads them and never re-runs the daily steps; report only (MISSING lines)

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
sys.path.insert(0, str(TOOLS))
import concat_choice as CC   # noqa: E402  (2026-11-17: which pupil tree the tables use; whole-eye approvals)


def pupil_root() -> Path:
    """the pupil/TST package tree the tables use (concat_choice.active_pupil_root: <project>/pupil_active_root.txt;
    since 2026-11-17 pupil_tst_teach = the package with the model retrained on the user's drawings), else
    pupil_tst_final (the original package run)"""
    return CC.active_pupil_root(PROJECT)
def local_paths(key: str) -> list:
    """this PC's folders from tools/local_paths.json (kept out of the repository, which stores no absolute paths):
    raw_dirs (where IC Capture / Logic save; searched), skip_roots (old project copies never searched),
    assay_video_dirs (TST / HOT recordings, assays_day.py). [] when the file or key is missing."""
    try:
        return [Path(p) for p in json.loads((TOOLS / "local_paths.json").read_text(encoding="utf-8")).get(key, [])]
    except (OSError, ValueError):
        return []


# where Python / Logic / IC Capture may have saved things (searched for files of the given day); raw_dirs of
# local_paths.json = the local raw folder (not OneDrive) IC Capture / Logic save to from day 5 (user 2026-11-18)
SEARCH = [
    (HOME / "Downloads" / "behavior_task_windows", 1),
    (HOME / "Downloads", 2),
    (HOME / "Documents", 3),
    (HOME / "Desktop", 3),
    (HOME / "Videos", 3),
    *[(p, 3) for p in local_paths("raw_dirs")],
    (HOME, 2),
    (PROJECT, 6),
]


def shell_folder(name: str):
    """Where Windows puts a user folder now (read only). OneDrive moved Documents and Desktop on 2026-11-16, so IC
    Capture / Logic may now save there instead of HOME\\Documents. None if unknown."""
    try:
        import os
        import winreg
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER,
                            r"Software\Microsoft\Windows\CurrentVersion\Explorer\User Shell Folders") as k:
            v, _ = winreg.QueryValueEx(k, name)
        return Path(os.path.expandvars(v))
    except (OSError, ImportError):
        return None


for _name in ("Personal", "Desktop", "My Video"):          # Documents, Desktop, Videos as Windows sees them now
    _p = shell_folder(_name)
    if _p is not None and _p.is_dir() and all(_p != r for r, _ in SEARCH):
        SEARCH.insert(len(SEARCH) - 2, (_p, 3))
SKIP_PARTS = {"AppData", ".claude", ".venv", "sync", "analysis", "_tests_and_aborted", "node_modules", "tools",
              "pupil", "$Recycle.Bin", "OneDrive"}
# 2026-11-18 the project was copied out of OneDrive to a local drive (user: "please use local drive not one drive");
# old copies (skip_roots of local_paths.json) stay as backups and must never be searched for a day's files (Documents
# still is, because IC Capture / Logic may still save there)
SKIP_ROOTS = {str(p).lower() for p in local_paths("skip_roots") if p.resolve() != PROJECT.resolve()}
NAME_RE = re.compile(r"^(?P<cage>\d{4}[mfMF]?)_*(?P<mouse>[A-Za-z]+\d*)_*"
                     r"(?P<task>fr\d+|pr\d*[a-z0-9]*|e?passive|active)?_*day(?P<day>\d+)",   # passive: days 6-10
                     re.I)
# 2026-11-19 (user, day 6): "you can make your rule more flexible if I made also a bit different name you can still
# find in file 'black' and 'Passive' no matter of position". The day-6 passive logs 9670fgreen1passive_prday6 and
# 9670fblackpassive_prday6 were 'name not understood' / mouse 'blackpassive' under NAME_RE alone. parse_name() keeps
# NAME_RE's answer whenever it gives a mouse already filed in that cage (so every name that worked before is read the
# same way) and otherwise finds the cage, the mouse colour and 'day<N>' anywhere in the name, in any order and case.
COLOURS = ("green1", "green2", "green3", "green4", "black", "white", "orange", "red", "blue", "yellow", "pink",
           "green")
CAGE_ANY = re.compile(r"(?<!\d)(\d{4})([mf])?(?!\d)", re.I)
DAY_ANY = re.compile(r"day[\s_-]*(\d{1,2})(?!\d)", re.I)


class NameMatch:
    """What NAME_RE.match gives (m.group('cage' | 'mouse' | 'task' | 'day'), m['...']), plus .flexible and .notes."""

    def __init__(self, cage, mouse, task, day, flexible=False, notes=()):
        self._g = {"cage": cage, "mouse": mouse, "task": task, "day": str(day)}
        self.flexible = flexible
        self.notes = list(notes)

    def group(self, key):
        return self._g[key]

    __getitem__ = group


def parse_name(name: str, known: dict | None = None):
    """cage / mouse / task / day of a Python log name (a NameMatch), or None when it cannot be read.
    1) NAME_RE (<cage><mouse>_<task>day<N>) when it gives a mouse already filed in that cage, the start of one (e.g.
       'green' -> green1 / green2, settled later by resolve_mouse_names), or when no mice are filed yet.
    2) Otherwise the cage (a filed cage folder such as 9670f, or its 4 digits, else any stand-alone 4 digits + m/f), the
       first filed mouse name / colour after it, and 'day<N>', anywhere in the name ('passive' anywhere too). Several
       mouse names -> the first one after the cage, with a CHECK note.
    3) Nothing found by 2) -> NAME_RE's answer if there was one (a new mouse name, as before), else None."""
    known = known or {}
    low = name.lower()
    m = NAME_RE.match(name)
    strict = None
    if m:
        cage, mouse = m.group("cage").lower(), m.group("mouse").lower()
        notes = []
        if known and cage not in known:                # e.g. '9669' without m/f -> the one filed cage 9669m
            kc = [k for k in known if k[:4] == cage[:4]]
            if len(kc) == 1:
                notes.append(f"cage '{cage}' read as the filed cage {kc[0]}")
                cage = kc[0]
        strict = NameMatch(cage, mouse, (m.group("task") or "").lower(), m.group("day"), notes=notes)
        km = known.get(cage, set())
        if not km or mouse in km or any(k.startswith(mouse) for k in km):
            return strict
    d = DAY_ANY.search(low)
    if not d:
        return strict
    cage, cpos = None, None
    for kc in sorted(known, key=len, reverse=True):
        i = low.find(kc)
        if i < 0 and re.fullmatch(r"\d{4}[mf]", kc):
            j = [mm.start() for mm in CAGE_ANY.finditer(low) if mm.group(1) == kc[:4]]
            i = j[0] if j else -1
        if i >= 0 and (cpos is None or i < cpos):
            cage, cpos = kc, i
    if cage is None:
        mm = CAGE_ANY.search(low)
        if not mm:
            return strict
        cage, cpos = (mm.group(1) + (mm.group(2) or "")).lower(), mm.start()
    km = known.get(cage, set())
    names = (km | {re.sub(r"\d+$", "", k) for k in km}) if km else set(COLOURS)
    hits = []
    for nm in sorted(names, key=len, reverse=True):          # green1 before green
        for mm in re.finditer(re.escape(nm), low):
            if not any(a <= mm.start() < a + len(b) for a, b in hits):
                hits.append((mm.start(), nm))
    if not hits:
        return strict
    after = sorted(h for h in hits if h[0] >= cpos)
    _, mouse = (after or sorted(hits))[0]
    day_n = int(d.group(1))
    notes = [f"name read by the flexible rule (cage, mouse and day<N> anywhere in the name): cage {cage}, mouse "
             f"{mouse}, day {day_n}"]
    others = sorted({b for _, b in hits} - {mouse})
    if others:
        notes.append(f"CHECK: the name also contains {', '.join(others)}; read as {mouse} (the first after the cage). "
                     f"If wrong, add a line to session_aliases.csv")
    return NameMatch(cage, mouse, "passive" if "passive" in low else "", day_n, flexible=True, notes=notes)


LICK_TAG = re.compile(r'"tag":\s*"(?:observed_)?lick"')   # 2026-11-19: Python lick records, for the CSV choice
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
                    if str(p).lower() in SKIP_ROOTS:          # the old project copy left in OneDrive (2026-11-18)
                        continue
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
                    rewards=sum('"reward_cmd"' in ln for ln in lines), paradigm=first.get("paradigm", ""),
                    licks=sum(bool(LICK_TAG.search(ln)) for ln in lines))   # lick / observed_lick records
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
            li = hdr.index("lick") if "lick" in hdr else None          # 2026-11-19: lick onsets, for the CSV choice
            prev, first_cam, last_t, n_cam = "0", None, None, 0
            rprev, r_fall, rises = "0", -1.0, []
            lprev, lick_on = "0", []
            for ln in f:
                parts = ln.split(",")
                if len(parts) <= ci:
                    continue
                t, v = float(parts[0]), parts[ci].strip()
                last_t = t
                if li is not None and len(parts) > li:
                    lv = parts[li].strip()
                    if lv == "1" and lprev != "1":
                        lick_on.append(t)
                    lprev = lv
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
        return dict(path=p, first_cam=first_cam, last=last_t, n_cam=n_cam, rewards=rises,
                    licks=lick_on if li is not None else None)
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
    known_now = known_mice()
    for p in sorted(uniq.values(), key=lambda x: x["start"]):
        m = parse_name(p["name"], known_now)           # 2026-11-19: flexible names (see parse_name)
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
            s["problems"].extend(m.notes)
            said, is_replay = "passive" in p["name"].lower(), "replay" in str(p.get("paradigm", "")).lower()
            if said != is_replay:                      # 2026-11-19: the name and the log disagree -> say so
                s["problems"].append(f"CHECK: the name {'says' if said else 'does not say'} passive, but the log's "
                                     f"paradigm is {p.get('paradigm') or 'not given'}")
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
    # .sal is known, (4) prefer the CSV saved in the same folder as that .sal.
    # 2026-11-19 (user: "네, 고쳐 주세요"): a passive replay has the SAME reward times as its active partner and every
    # export is called digital.csv, so (1)-(3) cannot tell the two apart (day 6: 9669m black and white got each
    # other's CSV). Before (4) the choice now prefers (a) a CSV whose path names the session (an export folder named
    # like the log, or the session's own project folder salae/dayN/<cage>/<mouse>/), (b) a CSV whose lick onsets in
    # the session agree with the Python log's lick / observed_lick records, (c) the CSV saved soonest after the
    # session ended; a CSV saved more than 60 s before the session ended comes last.
    for s in real:
        best = None
        n_rew, dur = s["py"]["rewards"], s["py"]["dur"]
        py_end = s["py"]["start"] + dur
        nm = s["py"]["name"].lower()
        try:
            own = (PROJECT / "salae" / f"day{s['day']}" / s["cage"] / s["mouse"]).resolve()
        except Exception:   # noqa: BLE001
            own = None
        py_lk = s["py"].get("licks")
        for x in csvs:
            if id(x) in used or x["first_cam"] is None:
                continue
            in_sess = [r for r in x["rewards"] if x["first_cam"] - 0.1 <= r <= x["first_cam"] + dur + 2]
            if abs(len(in_sess) - n_rew) > max(2, 0.03 * n_rew):
                continue
            if abs(x["n_cam"] - dur * CAM_HZ) > max(60, 0.03 * dur * CAM_HZ):
                continue
            named = nm in str(x["path"]).lower() or (own is not None and x["path"].parent.resolve() == own)
            lk = x.get("licks")
            if lk is not None and py_lk is not None:
                n_lk = sum(1 for t in lk if x["first_cam"] - 0.1 <= t <= x["first_cam"] + dur + 2)
                lick_bad = abs(n_lk - py_lk) > max(20, 0.5 * max(n_lk, py_lk))
            else:
                n_lk, lick_bad = None, False
            try:
                gap = x["path"].stat().st_mtime - py_end
            except OSError:
                gap = 0.0
            late = gap if gap >= -60 else 1e9
            if s["sal"] is not None:
                err = abs(x["first_cam"] - s["sal_lag"])
                if err > 1.0:
                    continue
                score = (not named, lick_bad, x["path"].parent != s["sal"]["path"].parent, abs(len(in_sess) - n_rew),
                         late, err)
            else:
                if not 0 <= gap <= 3600:
                    continue
                score = (not named, lick_bad, False, abs(len(in_sess) - n_rew), late, gap)
            if best is None or score < best[0]:
                best = (score, x, n_lk, lick_bad)
        if best:
            s["csv"] = best[1]; used.add(id(best[1]))
            if best[3]:
                s["problems"].append(f"CHECK: the chosen Saleae CSV has {best[2]} lick onsets in the session but the "
                                     f"Python log has {py_lk} lick records - it may belong to another session")
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


def remap_moved_paths(c: dict) -> list:
    """The project folder can move (OneDrive moved Documents, 2026-11-16). A config path inside a project folder that
    no longer exists is pointed at the same place under PROJECT. Only the one-day copy is changed; run_all_config.json
    itself is left as it is. Paths that exist are never touched."""
    changed = []
    for k, v in list(c.items()):
        if not isinstance(v, str) or not re.match(r"^[A-Za-z]:[\\/]", v):
            continue
        p = Path(v)
        if p.exists() or PROJECT.name not in p.parts:
            continue
        q = PROJECT.joinpath(*p.parts[p.parts.index(PROJECT.name) + 1:])
        if q.exists():
            c[k] = str(q)
            changed.append(f"{k}: {v} -> {q}")
    return changed


def run_package_day(day_n: int, tst: bool) -> bool:
    """Pupil (and TST when that day has TST videos) for ONE day with the user's pupil/TST package. Runs every day
    by default since 2026-11-14 (finished sessions are skipped, so the day-18 RUN_ALL_18DAYS.bat run only has to do
    the QC page and the sessions marked redo). Writes a one-day copy of run_all_config.json (days = [N]; TST videos
    of tst_video\\dayN only, results in pupil_tst_final\\TST_auto\\dayN) and runs MATLAB -batch (ROI: previous day's)."""
    # 2026-11-17: when the tables use the new package's tree (pupil_tst_teach), the one-day config is made from
    # run_all_config_teach_day1.json (same package code; the model retrained on the user's drawings; output in
    # pupil_tst_teach). run_all_config.json (the original package run) is not changed.
    teach = PACKAGE / "run_all_config_teach_day1.json"
    use_teach = pupil_root().name != "pupil_tst_final" and teach.exists()
    cfg = teach if use_teach else PACKAGE / "run_all_config.json"
    if not cfg.exists():
        print("\n(pupil/TST step skipped: package or run_all_config.json missing)")
        return False
    c = json.loads(cfg.read_text(encoding="utf-8"))
    c.update(days=[day_n], dry_run=False, run_pupil=True, run_tst=bool(tst), roi_mode="previous",
             _note=f"one-day run from daily_summary.py (day {day_n}); made from {cfg.name}")
    if use_teach and Path(str(c.get("pupil_out_root", ""))).resolve() != pupil_root().resolve():
        print(f"  !! {cfg.name} writes to {c.get('pupil_out_root')}, but the tables use {pupil_root()}")
    for m in remap_moved_paths(c):
        print(f"  (project folder moved) {m}")
    if tst:
        c.update(tst_video_root=str(PROJECT / "tst_video" / f"day{day_n}"),
                 tst_out_root=str(Path(c.get("pupil_out_root") or PROJECT / "pupil_tst_final") / "TST_auto" / f"day{day_n}"))
    day_cfg = PACKAGE / (f"run_all_config_teach_day{day_n}.json" if use_teach else f"run_all_config_day{day_n}.json")
    day_cfg.write_text(json.dumps(c, indent=2), encoding="utf-8")
    matlab = matlab_exe()
    if not matlab:
        print("\n(pupil/TST step skipped: MATLAB not found)")
        return False
    print(f"\n=== PUPIL{' + TST' if tst else ''} (package, day {day_n}; about 10 min per session) ===")
    cmd = f"cd('{PACKAGE}'); RUN_ALL_18DAYS('{day_cfg}')"
    subprocess.run([matlab, "-batch", cmd])
    return True


def run_wholeeye() -> bool:
    """the whole eye of the mice approved in <pupil tree>/wholeeye_mice.csv (2026-11-17: 9670f black, 9669m black - eye
    damage, the pupil cannot be seen): tracking crops (wholeeye_build_set.py --crops-only; a new day needs its eye box
    video\\dayN\\<cage>\\<mouse>\\roi.mat first) -> eye_track_all.m with the one-time eye model (sessions already done
    are skipped) -> wholeeye_clean.py (pupil_outputs_wholeeye\\wholeeye_CLEAN.csv) -> wholeeye_sessions.py (Step 02
    copy concat_out_*_wholeeye)."""
    root = pupil_root()
    if not CC.read_rows(root / CC.WHOLE_MICE, warn=None, need=("cage", "mouse")):
        print("  (no mouse in wholeeye_mice.csv)")
        return False
    # eye_unet_black3 (2026-11-18): 21 user outlines incl. the white/cloudy part of the 9670f black eye (user: the white
    # part IS eye); every whole-eye session uses this one model (black2 = 18 outlines, its outputs kept in
    # analysis\\wholeeye\\_superseded_20261118_black2)
    model = PROJECT / "analysis" / "wholeeye" / "model" / "eye_unet_black3.mat"
    matlab = matlab_exe()
    if not model.exists() or not matlab:
        print("  (whole-eye step skipped: the eye model or MATLAB is missing)")
        return False
    subprocess.run([sys.executable, str(TOOLS / "wholeeye_build_set.py"), "--crops-only"])
    ft = TOOLS / "opioidaddiction-matlab" / "current_cohort_workflow" / "pupil_finetune"
    crops, out = PROJECT / "analysis" / "wholeeye" / "crops.json", PROJECT / "analysis" / "wholeeye"
    subprocess.run([matlab, "-batch", f"maxNumCompThreads(8); cd('{ft.as_posix()}'); "
                                      f"eye_track_all('{crops.as_posix()}', '{model.as_posix()}', '{out.as_posix()}')"])
    subprocess.run([sys.executable, str(TOOLS / "wholeeye_clean.py"), "--project", str(PROJECT)])
    rc = subprocess.run([sys.executable, str(TOOLS / "wholeeye_sessions.py"), "--project", str(PROJECT)]).returncode
    if rc:
        print("  !! wholeeye_sessions.py reported a problem (see above)")
    return True


def update_step02_base(day_n: int) -> bool:
    """Copy the newest concat_out + task log (+ assay files) of every day into <project>\\step02_base (what the
    ORIGINAL Step 02 reads), so the longitudinal table can be built at any time with run_step02_on_base.m. All days,
    not only day_n: a session redone or fixed later replaces its old copy (moved to step02_base\\_superseded_<date>);
    copies that are already there are left alone."""
    script = WORKFLOW / "python" / "prepare_step02_base.py"
    if not script.exists():
        print("\n(step02_base not updated: current_cohort_workflow not found in tools)")
        return False
    print(f"\n=== STEP 02 BASE (all days; new today: day {day_n}) ===")
    subprocess.run([sys.executable, str(script), "--pupil-out", str(pupil_root()), "--python-root",
                    str(PROJECT / "python"), "--assays", str(PROJECT / "assays"), "--base", str(PROJECT / "step02_base"),
                    "--apply"])
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
    # whole-eye mice (2026-11-17): their Diameter_px is the whole eye, not the pupil -> not in the pupil figures
    try:
        import concat_choice as CC
        rows = CC.read_rows(pupil_root() / CC.WHOLE_MICE, need=("cage", "mouse"))
        keys = [f"{re.sub(r'[^0-9]', '', r['cage'])}_{r['mouse'].strip().lower()}" for r in rows
                if CC.is_approved(r.get("approved"))]
    except Exception as e:                       # never stop the figures for this
        print(f"  (whole-eye mice not read: {e})")
        keys = []
    if keys:
        cmd += ["--no-pupil-mice", *sorted(set(keys))]
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
                    help="(old flag; Step 02 + the cohort figures now run by default for the routine PPT)")
    ap.add_argument("--no-longitudinal", action="store_true",
                    help="skip Step 02 + the previous-cohort figures today (the PPT then shows the newest older run)")
    ap.add_argument("--no-ppt", action="store_true", help="skip the routine PowerPoint")
    ap.add_argument("--no-assays", action="store_true", help="skip the TST / HOT step (assays_day.py) today")
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
            # TST is no longer scored by the package's own TST step (2026-11-18): assays_day.py below runs the same
            # package scorer on the FIRST 3 MIN with the per-frame mouse check (user: "TST should be 3min so just use
            # first 3min", "always only mouse tracking"); the package run here does the pupil only
            # fixed eye-box size per mouse (added 2026-11-20, user: "네, 모두 다시" = one box size per mouse for all
            # days, because the box size changes the pupil size the package reads): report only, before the package
            print("\n=== EYE BOX SIZE CHECK (fixed size per mouse, video\\eye_box_sizes.csv; report only) ===")
            if subprocess.run([sys.executable, str(TOOLS / "eye_box_fixed.py"), "--project", str(PROJECT), "--check",
                               "--day", str(dn)]).returncode:
                print("  !! CHECK: set the day's boxes with eye_box_fixed.py --day N --set CAGE MOUSE CX CY --apply\n"
                      "     (eye opening centre read on the day's frames) before the pupil step")
            if run_package_day(dn, False):
                # approved alignment fixes (added 2026-11-15): sync-aligned copy of the package's concat_out for the
                # sessions in pupil_tst_final\alignment_fixes.csv (the user's approval); any other mismatch is only
                # listed by pupil_session_quality.py below (ask the user first)
                # ALL approved sessions, not only today's: a session redone later (RUN_ALL_18DAYS redo) gets its
                # fix again here; up-to-date sessions are skipped in a second
                print("\n=== ALIGNMENT FIXES (approved sessions only) ===")
                rc_fix = subprocess.run([sys.executable, str(TOOLS / "fix_concat_alignment.py"), "--project",
                                         str(PROJECT)]).returncode
                if rc_fix:
                    print("  !! fix_concat_alignment.py reported a problem (see above and alignment_fixes_log.csv)")
                # fully dilated pupils (added 2026-11-15): an OPTION for rare sessions - only the sessions the user
                # listed in fullpupil_sessions.csv get the full-pupil trace + Step 02 copy; nothing else changes.
                # (The screen for such sessions is not run daily: fullpupil_sessions.py --screen-only when wanted.)
                print("\n=== FULL PUPIL (listed sessions only) ===")
                rc_fp = subprocess.run([sys.executable, str(TOOLS / "fullpupil_sessions.py"), "--project",
                                        str(PROJECT)]).returncode
                if rc_fp:
                    print("  !! fullpupil_sessions.py reported a problem (see above)")
                # best-effort pupil of hard eyes (added 2026-11-16): an OPTION for the sessions the user approved in
                # besteffort_sessions.csv - a Step 02 copy with only the pupil columns from pupil_besteffort.py's trace
                print("\n=== BEST-EFFORT PUPIL (listed sessions only) ===")
                rc_be = subprocess.run([sys.executable, str(TOOLS / "besteffort_sessions.py"), "--project",
                                        str(PROJECT)]).returncode
                if rc_be:
                    print("  !! besteffort_sessions.py reported a problem (see above)")
                # sessions without video that the user approved for a behavior-only Step 02 table (2026-11-16)
                print("\n=== BEHAVIOR-ONLY TABLES (listed no-video sessions only) ===")
                rc_bo = subprocess.run([sys.executable, str(TOOLS / "behavioronly_sessions.py"), "--project",
                                        str(PROJECT)]).returncode
                if rc_bo:
                    print("  !! behavioronly_sessions.py reported a problem (see above)")
                # whole eye of mice with eye damage (2026-11-17; wholeeye_mice.csv in the pupil tree): the one-time eye
                # model tracks the new sessions, then their Step 02 copy gets the whole eye in the pupil columns
                print("\n=== WHOLE EYE (mice in wholeeye_mice.csv) ===")
                run_wholeeye()
                # dark-contrast repair of the package's pupil (2026-11-18; darkcontrast_sessions.csv in the pupil tree):
                # only the listed sessions get a Step 02 copy whose pupil columns are repaired where the package sat
                # on another dark blob (pupil_dark_contrast.py; the user: only 9670f green2 day 5)
                print("\n=== DARK-CONTRAST PUPIL REPAIR (listed sessions only) ===")
                rc_dc = subprocess.run([sys.executable, str(TOOLS / "darkcontrast_sessions.py"), "--project",
                                        str(PROJECT)]).returncode
                if rc_dc:
                    print("  !! darkcontrast_sessions.py reported a problem (see above)")
                # per-session pupil grade (added 2026-11-15): manual pupil_grades.csv, else automatic stability grade.
                # All days are recomputed (a redone day's alignment check stays current); grade lines: today only
                print("\n=== PUPIL QUALITY ===")
                subprocess.run([sys.executable, str(TOOLS / "pupil_session_quality.py"),
                                "--pupil-root", str(pupil_root()), "--report-days", str(dn)])
                subprocess.run([sys.executable, str(TOOLS / "pupil_guess_trace.py"),
                                "--pupil-root", str(pupil_root()), "--days", str(dn)])
                print("\n=== ANALYSIS (again, now with pupil) ===")
                quick_analysis.main(["--root", str(PROJECT)])
                # pupil check pictures (added 2026-11-18, user: "show me png files as example ... your confirmation this
                # session is okay if it is noisy how many frames are out"): one picture per session of the pupil the
                # tables use + the frames out; report only. Claude then looks at the new day's pictures and writes the
                # verdicts in analysis\\pupil_check\\pupil_check_verdicts.csv (the deck shows LOOK until then)
                print("\n=== PUPIL CHECK PICTURES (all days so far) ===")
                if subprocess.run([sys.executable, "-B", str(TOOLS / "pupil_check_examples.py"), "--days",
                                   *map(str, range(1, dn + 1))]).returncode:
                    print("  !! pupil_check_examples.py could not finish (see above)")
            # eye-box check (added 2026-11-15): the package re-uses the previous day's ROI; mice sit differently daily
            print("\n=== PUPIL ROI CHECK ===")
            subprocess.run([sys.executable, str(TOOLS / "pupil_roi_check.py"), "--root", str(PROJECT), "--day", str(dn),
                            "--pupil-root", str(pupil_root())])
            print("  If an eye is not inside its yellow box: set \"redraw_roi\" for that session in run_all_config.json\n"
                  "  and run RUN_ALL_18DAYS.bat (it opens the ROI window and redoes the session).")
        else:
            print("\n(pupil/TST: skipped today with --no-pupil; it can be run later, finished sessions are skipped)")
        # report-only checks (added 2026-11-15). Exit code 1 only means "there are CHECK lines", not a failure.
        # session_check: every session of this date vs tools\protocol_expected.csv (also run by hand after each
        # session with SESSION_CHECK.bat); pr_daily_checks: one row per mouse for the day (valve time, measured
        # lock-out, Python-only licks the task counted, session end, PR vs the previous cohort's days 3-5).
        print("\n=== SESSION CHECKS (report only) ===")
        rc_sc = subprocess.run([sys.executable, "-B", str(TOOLS / "session_check.py"), "--date", day.isoformat(),
                                "--all"]).returncode
        if rc_sc == 1:
            print("  CHECK lines above; reports in analysis\\session_checks\\")
        elif rc_sc not in (0, 2):
            print(f"  !! session_check.py ended with code {rc_sc} (see above)")
        print("\n=== PR CHECKS (report only) ===")
        rc_pr = subprocess.run([sys.executable, "-B", str(TOOLS / "pr_daily_checks.py"), "--project", str(PROJECT),
                                "--day", str(dn)]).returncode
        if rc_pr == 1:
            print(f"  CHECK lines: analysis\\days\\day{dn}\\pr_checks_day{dn}.txt")
        elif rc_pr:
            print("  !! pr_daily_checks.py could not write its outputs (file open in Excel?) - see above")
        # PR score of every mouse on every day so far, one table (added 2026-11-18, user: "all days' PR score table")
        print("\n=== PR TABLE (all days so far) ===")
        rc_pt = subprocess.run([sys.executable, "-B", str(TOOLS / "pr_table_days.py")]).returncode
        if rc_pt == 1:
            print("  !! rewards Python != Saleae in a session (see above)")
        elif rc_pt:
            print("  !! pr_table_days.py could not finish (see above)")
        # tail immersion of every mouse on every day from day 5, one table + plot (added 2026-11-20, user: "in summary
        # we should also show TI data aswell ofc in all days"); the values are the ones the user sent in chat
        print("\n=== TI TABLE (tail immersion, all days so far) ===")
        if subprocess.run([sys.executable, "-B", str(TOOLS / "ti_table_days.py")]).returncode:
            print("  !! ti_table_days.py could not finish (see above)")
        # TST + HOT of the day (added 2026-11-18; TST/HP days 5, 10, 13, 16, 18): videos copied into tst_video /
        # hp_video, TST scored automatically on the first 3 min -> assays\dayN\...\ManualScoringResults_<video>_AUTO.mat
        # (the name the original Step 02 reads), HOT first-30-s copies for the hand scorer and its .mat files filed;
        # nothing happens on days without assay videos. Before the step02 base update, which copies the assay files.
        if not args.no_assays:
            pkg_py = TOOLS / ".venv_pkg313" / "Scripts" / "python.exe"
            if pkg_py.exists():
                print()
                if subprocess.run([str(pkg_py), "-B", str(TOOLS / "assays_day.py"), "--day", str(dn),
                                   "--apply"]).returncode:
                    print("  !! assays_day.py could not finish (see above)")
            else:
                print("\n(assays step skipped: tools\\.venv_pkg313 missing)")
        # STRAUB TAIL (added 2026-11-21; user: "from now on you can use this for my current cohort from tmr"):
        # straub_day.py arranges <cage><mouse>_straub_dayN videos, makes first-5-min copies for the hand scorer, files
        # its .mat files and scores every copy AUTOMATICALLY (straub_auto.py: tail > 25 deg from the body axis or
        # raised over the back, 0.5-s majority; confirmed by the user on an annotated video 2026-11-21) ->
        # analysis\straub_auto\dayN + analysis\days\dayN\straub_dayN.csv. Nothing happens on days without Straub videos.
        if not args.no_assays:
            print()
            if subprocess.run([sys.executable, "-B", str(TOOLS / "straub_day.py"), "--day", str(dn), "--apply"]).returncode:
                print("  !! straub_day.py could not finish (see above)")
        update_step02_base(dn)
        # Step 02 + the previous cohort's figures every day since 2026-11-17 (the routine PPT shows them as the basic
        # plots); --no-longitudinal skips it
        if not args.no_longitudinal:
            run_longitudinal()
        # QC summary of ALL days so far (added 2026-11-17): the daily report is QC first (user: "daily routine need to
        # focus on QC ... later I will do analysis after 18 days"). One picture + csv; report only, changes nothing.
        print("\n=== QC SUMMARY (all days so far) ===")
        rc_qc = subprocess.run([sys.executable, "-B", str(TOOLS / "qc_summary.py")]).returncode
        if rc_qc:
            print("  !! qc_summary.py could not finish (see above)")
        # data index (added 2026-11-18, user: "we are making csv or output data per day so in order to use and recall
        # this for 18days later"): every saved data file with its relative path and which one is final
        print("\n=== DATA INDEX ===")
        if subprocess.run([sys.executable, "-B", str(TOOLS / "data_index.py"), "--quiet"]).returncode:
            print("  !! data_index.py could not finish (see above)")
        # independent EVALUATOR (added 2026-11-22, user: "I sometimes find you read wrongly file so please alsway every day
        # please make evaluator whether you read file correctly, and QC correctly and analysiscorrectly"): the raw files
        # are read again with its own small readers and compared with the routine's outputs (files, QC, analysis);
        # report only. Before the deck, which shows it right after the overview (FAIL items in 'Needs your decision')
        print("\n=== EVALUATOR (independent check: files read right? QC right? analysis right?) ===")
        rc_ev = subprocess.run([sys.executable, "-B", str(TOOLS / "daily_evaluator.py"), "--day", str(dn)]).returncode
        if rc_ev == 1:
            print(f"  !! FAIL items above (analysis\\days\\day{dn}\\evaluator_day{dn}.txt): look before using the day")
        elif rc_ev:
            print("  !! daily_evaluator.py could not finish (see above)")
        # routine PowerPoint (added 2026-11-17, user: "PPT: QC, problem fixed, folder address ... representative plots
        # as fundamental, what to do better next"): analysis\\ppt\\routine_dayN_<date>.pptx
        if not args.no_ppt:
            print("\n=== ROUTINE PPT ===")
            rc_ppt = subprocess.run([sys.executable, "-B", str(TOOLS / "routine_ppt.py")]).returncode
            if rc_ppt:
                print("  !! routine_ppt.py could not finish (see above)")
    days = sorted((PROJECT / "analysis" / "days").glob("day*")) if (PROJECT / "analysis" / "days").exists() else []
    print("\nDone. QC first:", PROJECT / "analysis" / "days" / "qc_summary_days1-N.png", "(red = needs a look)")
    print("  Per-day summaries:", ", ".join(d.name for d in days) or "-")
    print(f"  {PROJECT / 'analysis' / 'days'}\\dayN\\figs\\00_day_dashboard.png")
    print(f"  across days: {PROJECT / 'analysis' / 'master_day_level.csv'} and analysis\\figs_across_days\\")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
