#!/usr/bin/env python3
"""
session_check.py - read-only check to run right after each session: double-click SESSION_CHECK.bat once
IC Capture is stopped and the Saleae capture is saved (and exported to CSV), so recording mistakes are found
the same day, while the next mice can still be fixed.

The session's files are found the way the evening run finds them (daily_summary.scan + pair_day; nothing is
moved). The evening sync (sync_sessions.process_session) is run on them into analysis/session_checks/ (never
the real sync folder), and its numbers are compared with the row of tools/protocol_expected.csv for the
session's day:
  A files      Python log, .sal, Saleae CSV, video (a video whose recording started more than 3 s after
               Python is named: the evening run will not file it)
  B name       cage / mouse / day from the file name; the same cage/mouse as another session of the same date
               (CHECK); a '_N' at the end of the name (the task adds it when the output file already exists);
               day number already used on an earlier date, different from the other sessions of the same
               date, or not the next day (CHECK when no other session of that date uses the number). A full
               run whose name cannot be read is CHECK (the evening run would not file it); it is still checked
               with a provisional name.
  C protocol   paradigm, requirement sequence, lock-out (reward_cmd -> next trial_start), injector_run,
               session length, experiment_end. Several rows of one day (days 6-10: active PR and passive
               replay) = the row whose paradigm is the log's. A passive (replay) log also gets the log it
               replayed (replay_source): found, another mouse, same date, finished before this session, not a
               replay itself; the pair in cohort_manifest.csv (when filled); the replayed rewards vs the
               source's (count, and times from the first trial_start within 0.5 s)
  D valve      Saleae reward pulses (closer than 1 ms merged): open time, count vs Python reward_cmd, and
               whether the sync matched every reward_cmd to an opening
  E licks      Saleae vs Python licks; Python-only licks the task counted as valid (CHECK on PR days)
  F camera     clock: the camera trigger runs from 'S' (Python start) to 'T' (sent with experiment_end), so
               its pulse train must start at the Python start, end at experiment_end and last as long as
               Python says (PC clock changes are also looked up in the Windows event log). When they disagree,
               the train is the session window, and reward/lick matching is judged only on the side of the
               clock change where the sync's single clock line fits; camera pulses vs video frames; first
               pulse; pulses after the end (the session's own train or another run's); video closed before
               the end; other Python runs while the video was recorded. Without a Saleae CSV: frames vs
               Python duration x 30.26 Hz (estimate).
  G Saleae     sample rate, threshold, channel names, 1-us glitch filter (from the .sal meta.json)
  H sync       status and warnings of the evening sync
Each line is PASS / CHECK / INFO. INFO = context, or no expectation for this day (empty cell or no row).
Python logs that are not checked (test runs, aborted runs, empty 0-byte logs) are listed; an empty log or an
aborted run with a session name counts as something to check.
Test run (a log whose name starts with or contains 'test', e.g. the morning valve test): when it is the newest
log of today (default: checked after the newest session, if any) or named with --session, it gets a reduced
check: valve open time vs the expected valve of the day ('day<N>' in the test name, else the day of today's
sessions, else the day after the last filed day), Python reward_cmd vs Saleae openings, Saleae settings, video
frames vs camera pulses (when there is a video), camera trigger stopped at its end. Its .sal / CSV: a capture
with a camera pulse train where the .sal says the test started (also a test recorded into a session's capture).

Usage (tools .venv Python, or SESSION_CHECK.bat with the same options):
  session_check.py                               newest finished session of today (+ the newest log if it is a
                                                 test run)
  session_check.py --all                         every session of today
  session_check.py --session 9669mblack_prday3   one session (or test run) by its Python log name
  session_check.py --date 2026-11-13 --all       another day (works on the filed copies; also with --session)
  session_check.py --project FOLDER ...          another project folder (tests): only that folder is searched

Outputs (only under <project>/analysis/session_checks/):
  day<N>/<session>_check.txt   the report (day_unknown/ when the day cannot be told); <session>_<date>_check.txt
                               when that name already has a report from another date
  day<N>/<session>/            the sync outputs (_events.csv, _frames.csv, _report.txt); <session>_<date>/ as above
  test_runs/<name>_<date>_<time>_check.txt   the report of a test run
  session_checks.csv           one row per session with the key numbers (a re-check of the same session - same
                               day, name and start - replaces its row; test runs are not added)
Raw files and the real sync folder are only read. If a Python log or one of the session's files is still
being written, the tool says which file and does not check that session.
Exit code: 0 all PASS, 1 something to check, 2 nothing to check, 3 a file is still being written (and nothing
else needs checking), 4 all PASS but session_checks.csv could not be updated (open in Excel? close it and run
again).
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import traceback
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
from types import SimpleNamespace

sys.dont_write_bytecode = True        # importing the tools below must not write tools/__pycache__

import numpy as np                    # noqa: E402
import pandas as pd                   # noqa: E402

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
import daily_summary as ds            # noqa: E402  file search + pairing, in_use, NAME_RE, known_mice
import sync_sessions as ss            # noqa: E402  the evening sync and its readers
import quick_analysis as qa           # noqa: E402  GLITCH_MS (valve pulses shorter than this are glitches)

PROJECT = TOOLS.parent
CAM_HZ = ds.CAM_HZ                    # 30.26 Hz camera trigger
PREV_VALVE_MS = 179.8                 # previous cohort PR valve (context only)
PREV_SALEAE = "10 MHz, 1.2 V"         # previous cohort Saleae settings (context only)
CHANNELS = {4: "reward", 6: "camera", 7: "lick"}
FIRST_PULSE_TOL_MS = 50.0             # first camera pulse vs Python start (the sync FAILs beyond this too)
PAD_BEFORE, PAD_AFTER = 0.05, 0.1     # s before Python start / after its end that still count (as in the sync)
CLOCK_TOL_S = 0.1                     # camera train vs Python on the Saleae clock: first pulse ('S') vs Python start
                                      # (normal +-1 ms), last pulse ('T') vs experiment_end (normal -35..0 ms), and
                                      # train length vs Python duration
SAL_LAG_TOL_S = 0.6                   # train start vs where the .sal says Python started (that stamp is ~0.2 s late)
STEP_MAX_S = 5.0                      # no .sal: a train starting up to this much before the clock line's Python
                                      # start is a PC clock change; earlier = the trigger was already running
STEP_MARGIN_S = 1.0                   # matching judged after a clock change: skip this much more (lick hand-over)
LATE_VIDEO_S = 3.0                    # pair_day pairs videos created up to this long after the Python start
STEP_MIN_MS = 100                     # PC clock changes in the Windows event log at least this large are reported
FRAME_EST_TOL = 30                    # no Saleae CSV: video frames vs Python duration x 30.26 Hz, PASS within this
RUN_MATCH_S = 2.5                     # another camera pulse train = a Python run if start and length agree this well
MATCH_FAIL_FRAC = 0.25                # after a PC clock change: more unmatched Python licks than this = the sync's
                                      # lick matching failed (real re-triggers: ~1 %, at most 6.5 % with the 500 ms valve)
REPLAY_TIME_TOL_S = 0.5               # passive replay: reward times vs the replayed log's (previous cohort, 60 complete
                                      # replays: 0.02-0.24 s late at the end, about 2 ms per reward)
SYNC_WARN_KO = [("licks:", "lick 수 불일치"), ("lick delay", "lick 지연"), ("video has", "비디오 프레임"),
                ("video file created", "비디오 짝"), ("video:", "비디오 읽기"), ("video offset", "비디오 offset"),
                ("camera pulses stop", "Saleae 조기 종료"), ("first camera pulse", "첫 펄스"),
                ("gap(s) in camera", "카메라 펄스 빠짐"), ("rewards:", "보상 수"), ("experiment_end", "experiment_end 없음"),
                ("drift", "시계 drift"), ("HIGH when", "채널 HIGH 시작"), ("column", "CSV 채널"),
                ("python start is", "Saleae 짝")]
CSV_COLS = ["checked_at", "day", "session", "cage", "mouse", "start", "overall", "n_check", "check_items",
            "paradigm", "req_start", "req_step", "trials", "lockout_median_s", "injector_run", "session_s",
            "experiment_end", "reward_cmd", "valve_openings", "valve_median_ms", "valve_min_ms", "valve_max_ms",
            "valve_glitches", "rewards_matched", "licks_saleae", "licks_python", "licks_python_valid",
            "python_only_licks", "python_only_valid", "camera_pulses", "video_frames", "frames_minus_pulses",
            "pulses_before_start", "pulses_after_end", "first_pulse_ms", "clock_end_s", "clock_start_s",
            "pc_clock_steps",
            "video_closed_minus_end_s", "other_runs_in_video", "sample_rate_mhz", "threshold", "glitch_1us",
            "sync_status", "python", "sal", "saleae_csv", "video", "report"]
EVENT_NS = "{http://schemas.microsoft.com/win/2004/08/events/event}"

_read_python = ds.read_python          # the original reader, kept before main() installs the wrapper below


def read_python_closed_is_finished(p):
    """daily_summary.read_python, but a log without experiment_end that no program has open any more counts as
    finished: a killed task is paired and checked at once (C reports the missing experiment_end) instead of
    looking 'running' for 5 minutes."""
    x = _read_python(p)
    if x and not x["ended"] and not ds.in_use(p):
        x = dict(x, ended=True)
    return x


# ----------------------------------------------------------------------------- small helpers
def wall(e):
    """Wall-clock time (epoch s) of a log record, None if it has no readable 't'."""
    try:
        return datetime.fromisoformat(e["t"]).timestamp()
    except Exception:
        return None


def hms(ts):
    return datetime.fromtimestamp(ts).strftime("%H:%M:%S") if ts is not None else "-"


def day_of(ts):
    return datetime.fromtimestamp(ts).date()


def start_of(s):
    return s["py"]["start"]


def short_id(s):
    d = s["day"] if s.get("day") is not None else "?"
    sid = f"d{d} {s['cage']} {s['mouse']}"
    return f"{sid} ({s['py']['name']})" if s.get("provisional") else sid


def trains(cam):
    """(first index, last index) of every camera pulse train (a gap > 1 s starts a new train, as in the sync)."""
    starts = ss.pulse_trains(cam)
    return list(zip(starts, [i - 1 for i in starts[1:]] + [len(cam) - 1]))


def log_start(p: Path):
    """Start time (epoch s) from the meta line of a Python log; None if the first line is not a meta line."""
    try:
        with open(p, encoding="utf-8-sig", errors="replace") as f:
            first = json.loads(f.readline())
        return datetime.fromisoformat(first["t"]).timestamp() if first.get("tag") == "meta" else None
    except Exception:
        return None


def filed_logs(project: Path):
    """{day number: [(start date, log name, start epoch s)]} of the Python logs already filed in python/day<N>/."""
    out = {}
    for d in (project / "python").glob("day*"):
        if d.is_dir() and d.name[3:].isdigit():
            for p in d.rglob("*.jsonl"):
                t = log_start(p)
                if t is not None:
                    out.setdefault(int(d.name[3:]), []).append((day_of(t), p.stem, t))
    return out


def csv_first_time(path: Path):
    """Time of the first data row of a Saleae CSV (the capture start, normally 0)."""
    try:
        with open(path, encoding="utf-8-sig", errors="replace") as f:
            f.readline()
            return float(f.readline().split(",")[0])
    except (OSError, ValueError, IndexError):
        return 0.0


def report_start(path: Path):
    """'YYYY-MM-DD HH:MM:SS' Python start written in the first line of an existing report, None if none."""
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            m = re.search(r"\(Python start (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})", f.readline())
        return m.group(1) if m else None
    except OSError:
        return None


def logs_named(ctx, stem):
    """Python logs called <stem> (today's, filed, or under _tests_and_aborted): [(start epoch s or None, path)]."""
    key, out = stem.lower(), []
    paths = list(ctx.found["jsonl"]) + sorted((ctx.project / "_tests_and_aborted" / "python").glob("*.jsonl"))
    paths += sorted((ctx.project / "python").glob("day*/*/*/*.jsonl"))
    for p in paths:
        if p.stem.lower() == key:
            out.append((log_start(p), p))
    return out


def read_protocol(project: Path):
    p = project / "tools" / "protocol_expected.csv"
    if not p.exists():
        p = TOOLS / "protocol_expected.csv"
    rows = []
    if p.exists():
        with open(p, newline="", encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                rows.append({(k or "").strip(): (v or "").strip() for k, v in r.items()})
    return rows, p


def protocol_for(rows, day_n, paradigm=None):
    """The protocol_expected.csv row of this day. Several rows can cover one day (days 6-10: the active PR row and
    the passive replay row): the row whose paradigm list holds the log's paradigm is used, else the first one."""
    if day_n is None:
        return None
    hits = []
    for r in rows:
        try:
            if int(r["from_day"]) <= day_n <= int(r["to_day"]):
                hits.append(r)
        except (KeyError, ValueError):
            continue
    if paradigm:
        for r in hits:
            if paradigm in [p for p in r.get("paradigm", "").split("|") if p]:
                return r
    return hits[0] if hits else None


def num(row, key):
    """Number in a protocol cell; None when the cell is empty (= no expectation)."""
    try:
        return float(row[key]) if row and row.get(key) else None
    except ValueError:
        return None


def read_sal_settings(path: Path):
    """Sample rate, threshold, digital channel names, enabled channels and glitch filters of a Logic 2 .sal."""
    with zipfile.ZipFile(path) as z:
        m = json.loads(z.read("meta.json"))["data"]
    ls = m.get("legacySettings", {})
    names = {r["channel"]["deviceChannel"]: r.get("name", "") for r in m.get("rowsSettings", [])
             if r.get("type") == "channel" and r.get("channel", {}).get("type") == "Digital"}
    enabled = [c["index"] for c in ls.get("enabledChannels", []) if c.get("type") == "Digital"]
    gf = ls.get("glitchFilter", {})
    glitch = {}
    if gf.get("enabled"):
        for c in gf.get("channels", []):
            if c.get("filter", {}).get("enabled"):
                glitch[c["channel"]["deviceChannel"]] = c["filter"].get("widthSec")
    return dict(rate=ls.get("sampleRate", {}).get("digital"),
                threshold=ls.get("digitalThreshold", {}).get("description", ""),
                names=names, enabled=enabled, glitch=glitch)


def read_events(path: Path):
    """The sync's <session>_events.csv with the columns used here as numbers."""
    ev = pd.read_csv(path)
    for c in ("t_saleae_s", "duration_ms", "python_elapsed_s", "valid_or_hit"):
        ev[c] = pd.to_numeric(ev[c], errors="coerce")
    return ev


def utc_epoch(text):
    """'2026-11-12T22:47:14.735953000Z' (event log time, UTC, up to 9 decimals) -> epoch s."""
    main_part, _, frac = text.strip().rstrip("Z").partition(".")
    t = datetime.strptime(main_part, "%Y-%m-%dT%H:%M:%S").replace(tzinfo=timezone.utc).timestamp()
    return t + (float("0." + frac) if frac else 0.0)


def pc_clock_steps(t0, t1):
    """Changes of the PC clock of at least STEP_MIN_MS whose old or new time lies between t0 and t1 (epoch s), from
    the Windows System log (Kernel-General event 1; read-only query with wevtutil). Each: old, new (epoch s, the
    clock before / after), delta_ms, by (program). None when the event log cannot be read."""
    if sys.platform != "win32":
        return None
    utc = lambda t: datetime.fromtimestamp(t, timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.000Z")   # noqa: E731
    q = ("*[System[Provider[@Name='Microsoft-Windows-Kernel-General'] and (EventID=1) and "
         f"TimeCreated[@SystemTime>='{utc(t0 - 3600)}' and @SystemTime<='{utc(t1 + 3600)}']]]")
    try:
        r = subprocess.run(["wevtutil", "qe", "System", f"/q:{q}", "/f:xml"], capture_output=True, text=True,
                           errors="replace", timeout=30)
        if r.returncode != 0:
            return None
        root = ET.fromstring("<events>" + r.stdout + "</events>")
    except Exception:
        return None
    out = []
    for e in root.iter(EVENT_NS + "Event"):
        data = {d.get("Name"): (d.text or "") for d in e.iter(EVENT_NS + "Data")}
        try:
            old, new = utc_epoch(data["OldTime"]), utc_epoch(data["NewTime"])
        except (KeyError, ValueError):
            continue
        if abs(new - old) * 1000 < STEP_MIN_MS or max(old, new) < t0 or min(old, new) > t1:
            continue
        by = Path(data.get("ProcessName", "").replace("\\", "/")).name
        out.append(dict(old=old, new=new, delta_ms=(new - old) * 1000,
                        by="Windows time sync" if by.lower() == "svchost.exe" else (by or "unknown program")))
    return sorted(out, key=lambda x: x["new"])


def step_text(x, start=None):
    when = f", about {min(x['old'], x['new']) - start:.0f} s into the session" if start is not None else ""
    return f"PC clock changed by {x['delta_ms']:+.0f} ms at {hms(x['new'])} ({x['by']}{when})"


class Report:
    """The PASS / CHECK / INFO lines of one session."""
    def __init__(self):
        self.lines, self.items, self.items_ko = [], [], []

    def add(self, level, sec, label, text, ko=None):
        self.lines.append(f"  {level:5s} {sec} {label:14s} {text}")
        if level == "CHECK":
            self.items.append(f"{sec} {label}")
            self.items_ko.append(f"{sec} {ko or label}")

    def note(self, text):
        self.lines.append(f"{'':24s}{text}")


# ----------------------------------------------------------------------------- finding files
def other_runs(project: Path, found, day: date):
    """Every Python run that started that day (sessions, tests, aborted runs), one per start time."""
    paths = list(found["jsonl"]) + sorted((project / "_tests_and_aborted" / "python").glob("*.jsonl"))
    runs = {}
    for p in paths:
        x = _read_python(p)
        if x and day_of(x["start"]) == day:
            runs.setdefault(round(x["start"]), dict(name=x["name"], start=x["start"], end=x["start"] + x["dur"]))
    return sorted(runs.values(), key=lambda r: r["start"])


def paired_paths(sessions):
    return {str(x[k]["path"]).lower() for x in sessions for k in ("sal", "csv", "video") if x.get(k)}


def busy_files(s, sessions, found):
    """Files of this session that a program still has open, plus unpaired recordings saved after its start that
    are still open (an IC Capture video that is still recording cannot be paired yet)."""
    own = [s["py"]["path"]] + [s[k]["path"] for k in ("sal", "csv", "video") if s.get(k)]
    busy = [p for p in own if ds.in_use(p)]
    paired = paired_paths(sessions)
    for kind in ("video", "sal", "csv"):
        if s.get(kind):
            continue                       # this session already has one: another open file is not its file
        for p in found[kind]:
            if str(p).lower() in paired or (kind == "csv" and not p.name.lower().startswith("digital")):
                continue
            try:
                recent = p.stat().st_mtime >= s["py"]["start"] - 60
            except OSError:
                continue
            if recent and ds.in_use(p):
                busy.append(p)
    return busy


def csv_next_to_sal(s, sessions):
    """pair_day found a .sal but no CSV: use digital.csv next to the .sal when no other session has it and its
    first camera pulse sits where the .sal says Python started (+-1 s). Returns (path, note) or (None, None)."""
    dig = s["sal"]["path"].parent / "digital.csv"
    if not dig.is_file() or str(dig).lower() in paired_paths(sessions):
        return None, None
    x = ds.read_csv(dig)
    if not x or x["first_cam"] is None or abs(x["first_cam"] - s.get("sal_lag", 1e9)) > 1.0:
        return None, None
    dur, n_rew = s["py"]["dur"], s["py"]["rewards"]
    in_sess = [r for r in x["rewards"] if x["first_cam"] - 0.1 <= r <= x["first_cam"] + dur + 2]
    note = (f"pair_day did not match this CSV ({x['n_cam']} camera pulses vs {dur * CAM_HZ:.0f} expected for "
            f"{dur:.0f} s, {len(in_sess)} valve pulses vs {n_rew} Python rewards), so the evening run will not file "
            f"it by itself; used here because it sits next to the .sal and its first camera pulse "
            f"({x['first_cam']:.2f} s into the capture) is at the Python start ({s['sal_lag']:.2f} s)")
    return dig, note


def unpaired_video_covering(s, sessions, found, ff):
    """An unpaired video that was being written while this session ran (e.g. a test run recorded first)."""
    paired = paired_paths(sessions)
    t0 = s["py"]["start"]
    for p in found["video"]:
        if str(p).lower() in paired:
            continue
        try:
            if p.stat().st_size <= 1_000_000 or p.stat().st_mtime < t0:
                continue
        except OSError:
            continue
        v = ds.read_video(p, ff)
        if v.get("created") is not None and v["created"] <= t0 <= p.stat().st_mtime:
            return p, t0 - v["created"]
    return None, None


def late_videos(s, sessions, found, ff):
    """Unpaired videos created more than LATE_VIDEO_S after the Python start and before its end (IC Capture
    started late: pair_day does not pair them, so the evening run will not file them): [(path, seconds late)]."""
    paired = paired_paths(sessions)
    t0, t1 = s["py"]["start"], s["py"]["start"] + s["py"]["dur"]
    out = []
    for p in found["video"]:
        if str(p).lower() in paired:
            continue
        try:
            if p.stat().st_size <= 1_000_000 or p.stat().st_mtime < t0:
                continue
        except OSError:
            continue
        v = ds.read_video(p, ff)
        if v.get("created") is not None and t0 + LATE_VIDEO_S < v["created"] < t1:
            out.append((p, v["created"] - t0))
    return sorted(out, key=lambda x: x[1])


def guess_day(ctx, this_date):
    """Day number for a run whose name has none: (number or None, where it comes from)."""
    same = [s["day"] for s in ctx.sessions if s["status"] == "session" and day_of(s["py"]["start"]) == this_date]
    if same:
        n, k = Counter(same).most_common(1)[0]
        return n, f"day number of {k} other session(s) of {this_date}"
    filed_today = sorted(n for n, v in ctx.filed.items() if any(x[0] == this_date for x in v))
    if filed_today:
        return filed_today[-1], f"python/day{filed_today[-1]}/ holds logs of {this_date}"
    prev = [n for n, v in ctx.filed.items() if any(x[0] < this_date for x in v)]
    if prev:
        k = max(prev)
        return k + 1, f"the day after the last filed day {k} ({max(x[0] for x in ctx.filed[k] if x[0] < this_date)})"
    return None, "no other session or filed day to compare with"


def provisional_sessions(ctx):
    """Full runs that pair_day calls 'name not understood' (the evening run will not file them): paired again the
    same way, with a provisional name made of the cage/mouse at the start of the file name and a guessed day,
    so A and C-H can still be checked. Files that belong to a real session are not used."""
    todo = {round(s["py"]["start"]): s for s in ctx.sessions if s["status"] == "name not understood"
            and "test" not in s["py"]["name"].lower()}                # a long test run stays a test run
    if not todo:
        return []
    names = {}                                     # every copy of such a run (same start time) gets the name
    for p in ctx.found["jsonl"]:
        x = _read_python(p)
        s = todo.get(round(x["start"])) if x else None
        if s is None:
            continue
        m = re.match(r"^(\d{4}[mfMF]?)_*([A-Za-z]+\d*)", s["py"]["name"])
        cage, mouse = (m.group(1).lower(), m.group(2).lower()) if m else ("0000", "unknown")
        n, src = guess_day(ctx, day_of(s["py"]["start"]))
        names[str(p).lower()] = (f"{cage}_{mouse}_day{n if n is not None else 0}", n, src)

    def provisional_reader(p):
        x = read_python_closed_is_finished(p)
        hit = names.get(str(p).lower())
        return dict(x, name=hit[0]) if (x and hit) else x

    ds.read_python = provisional_reader
    try:
        again = ds.pair_day(ctx.found, ctx.ffmpeg, [])
    except Exception as exc:                       # such runs are then listed as 'not checked' (CHECK)
        print(f"  (could not pair the run(s) with a name the evening run cannot read: {exc})")
        return []
    finally:
        ds.read_python = read_python_closed_is_finished
    taken = paired_paths([s for s in ctx.sessions if s["status"] == "session"])
    out = []
    for s in again:
        hit = names.get(str(s["py"]["path"]).lower())
        if hit is None or s["status"] != "session":
            continue
        s["py"] = dict(s["py"], name=Path(s["py"]["path"]).stem)
        s["day"], s["provisional"] = hit[1], dict(day_src=hit[2])
        for k in ("sal", "csv", "video"):
            if s.get(k) and str(s[k]["path"]).lower() in taken:
                s["problems"].append(f"{Path(s[k]['path']).name} fits too, but it is paired with another session "
                                     f"- not used")
                s[k] = None
        out.append(s)
    return out


def other_logs(ctx):
    """Python logs of the day that are not checked: test runs, aborted runs, logs still being written, empty or
    unreadable logs. Each: name, path, start, kind (English), ko (Korean), check (counts as something to
    check), busy (still open)."""
    out = []
    checked = {str(s["py"]["path"]).lower() for s in ctx.provisional}
    for s in ctx.sessions:
        if s["status"] == "session" or str(s["py"]["path"]).lower() in checked:
            continue
        p = s["py"]
        busy = s["status"].startswith("RUNNING") or ds.in_use(p["path"])
        test = False
        if busy:
            kind, ko, check = "STILL BEING WRITTEN (a task is running, or the log is open in another program)", \
                "기록 중", False
        elif "test" in p["name"].lower():
            kind, ko, check = f"test run ({p['dur']:.0f} s)", f"테스트 {p['dur']:.0f} s", False
            test = True
        elif s["status"] == "name not understood":
            kind, ko, check = (f"name not understood ({p['dur']:.0f} s): the evening run will NOT file it, and it "
                               f"could not be checked here", "파일 이름 인식 불가", True)
        elif ds.parse_name(p["name"], ds.known_mice()):     # 2026-11-19: the flexible name rule
            kind, ko, check = (f"aborted run with a session name ({p['dur']:.0f} s, shorter than "
                               f"{ds.MIN_SESSION_S} s; the evening run files it under _tests_and_aborted)",
                               f"중단된 실행 {p['dur']:.0f} s", True)
            test = True    # may be the morning valve test run under the session name: valve check too (2026-11-16)
        else:
            kind, ko, check = f"short run ({p['dur']:.0f} s): test or aborted", f"짧은 실행 {p['dur']:.0f} s", False
            test = True    # e.g. a valve test saved as '11_1.jsonl': valve check too (2026-11-16)
        out.append(dict(name=p["name"], path=p["path"], start=p["start"], kind=kind, ko=ko, check=check, busy=busy,
                        test=test, py=p))
    in_pairing = {str(s["py"]["path"]).lower() for s in ctx.sessions}
    for path in ctx.found["jsonl"]:
        if str(path).lower() in in_pairing or _read_python(path):
            continue                    # paired, or another copy of a run that pair_day already has
        try:
            st = path.stat()
        except OSError:
            continue
        busy = ds.in_use(path)
        if busy:
            kind, ko, check = "STILL BEING WRITTEN (the task is starting?)", "기록 중", False
        elif st.st_size == 0:
            kind, ko, check = ("EMPTY (0 bytes): the task stopped before its first log line (n_trials too large, "
                               "COM port or Arduino handshake?)", "빈 로그(0 bytes)", True)
        else:
            kind, ko, check = "not a task log (its first line is not a meta line)", "읽을 수 없는 로그", False
        out.append(dict(name=path.stem, path=path, start=st.st_mtime, kind=kind, ko=ko, check=check, busy=busy,
                        test=False, py=None))
    return sorted(out, key=lambda o: o["start"])


# ----------------------------------------------------------------------------- the checks
def check_files(rep, s, ctx):
    """A: the four files of the session. Returns the .sal, Saleae CSV and video paths (None = not found)."""
    py = s["py"]
    sal_p = s["sal"]["path"] if s.get("sal") else None
    csv_p = s["csv"]["path"] if s.get("csv") else None
    vid_p = s["video"]["path"] if s.get("video") else None
    csv_note = None
    if csv_p is None and sal_p is not None:
        csv_p, csv_note = csv_next_to_sal(s, ctx.all)
    hints = {".sal": "save the Logic capture", "Saleae CSV": "export it from Logic (CSV, time in seconds)",
             "video": "no IC Capture video that starts at the Python start"}
    late = late_videos(s, ctx.all, ctx.found, ctx.ffmpeg) if vid_p is None else []
    for label, p in (("python log", py["path"]), (".sal", sal_p), ("Saleae CSV", csv_p), ("video", vid_p)):
        if p is None and label == "video" and late:
            q, dt = late[0]
            k = int(round(dt * CAM_HZ))
            rep.add("CHECK", "A", label, f"NOT PAIRED - {q}: recording started {dt:.0f} s late (about {k} frames "
                    f"missing); the evening run will not file it", ko=f"비디오 녹화 {dt:.0f} s 늦게 시작(약 {k}프레임 없음)")
            rep.note("IC Capture stamps whole seconds; the frames are not checked below because the video is not "
                     "paired")
            for q2, dt2 in late[1:]:
                rep.note(f"{q2}: also created {dt2:.0f} s after the Python start")
        elif p is None:
            rep.add("CHECK", "A", label, f"NOT FOUND - {hints[label]}", ko=f"{label} 없음")
        elif label == "Saleae CSV" and csv_note:
            rep.add("CHECK", "A", label, str(p), ko="Saleae CSV 짝 확인")
            rep.note(csv_note)
        else:
            rep.add("PASS", "A", label, str(p))
    if vid_p is None:
        p, lead = unpaired_video_covering(s, ctx.all, ctx.found, ctx.ffmpeg)
        if p:
            rep.note(f"{p} was being written during this session but starts {lead:.0f} s before the Python start "
                     f"(another run recorded into its start?) - not used")
    for msg in s["problems"]:
        if re.fullmatch(r"no \w+ found", msg) or "mouse name" in msg or "session_aliases" in msg:
            continue                       # missing files are listed above, names in B
        if "unfinished" in msg:
            rep.add("CHECK", "A", "pairing", msg, ko="비디오 미완료")
        else:
            rep.add("INFO", "A", "pairing", msg)
    return sal_p, csv_p, vid_p


def check_name(rep, s, ctx, row):
    """B: cage / mouse / day number from the file name."""
    name = s["py"]["name"]
    if s.get("provisional"):
        cage, mouse, task = s["cage"], s["mouse"], ""
        ex = f"{cage}{mouse}_{(row or {}).get('phase', 'pr').lower()}day{s['day'] if s['day'] is not None else 'N'}"
        rep.add("CHECK", "B", "name", f"'{name}' cannot be read as <cage><mouse>_<task>day<N> (e.g. {ex}): the "
                f"evening run will NOT file this session - rename the Python log before the evening run",
                ko="파일 이름 인식 불가")
        rep.note(f"checked here with a provisional cage {cage}, mouse {mouse}, day "
                 f"{s['day'] if s['day'] is not None else '?'} ({s['provisional']['day_src']}); its .sal, CSV and "
                 f"video were paired by time like any session")
    else:
        m = ds.parse_name(name, ds.known_mice())      # 2026-11-19: the flexible name rule of daily_summary
        cage, mouse, task = m.group("cage").lower(), m.group("mouse").lower(), (m.group("task") or "").lower()
        for nt in m.notes:
            if nt.startswith("CHECK"):
                rep.add("CHECK", "B", "name", nt, ko="이름에 다른 쥐 이름도 있음")
            else:
                rep.add("INFO", "B", "name", nt)
    known = ds.known_mice()
    if not known:
        rep.add("INFO", "B", "cage/mouse", f"cage {cage}, mouse {mouse}, task '{task or '-'}' "
                f"(no filed mice yet to compare with)")
    elif cage not in known:
        rep.add("CHECK", "B", "cage/mouse", f"cage '{cage}' is not a known cage (known: {', '.join(sorted(known))})",
                ko=f"케이지 이름 '{cage}'")
    elif mouse not in known[cage]:
        rep.add("CHECK", "B", "cage/mouse", f"mouse '{mouse}' is not a known {cage} mouse "
                f"(known: {', '.join(sorted(known[cage]))})", ko=f"마우스 이름 '{mouse}'")
        for p in s["problems"]:
            if "mouse name" in p:
                rep.note(f"pair_day: {p}")
    else:
        rep.add("PASS", "B", "cage/mouse", f"cage {cage}, mouse {mouse}, task '{task or '-'}'")
    for p in s["problems"]:
        if "session_aliases.csv ->" in p:
            rep.note(p)

    # the same cage/mouse twice on one date: usually the output name was not changed for the next mouse
    this_date = day_of(s["py"]["start"])
    if (s["cage"], s["mouse"]) != ("0000", "unknown"):
        twins = [x for x in ctx.all if x is not s and day_of(x["py"]["start"]) == this_date
                 and x["cage"] == s["cage"] and x["mouse"] == s["mouse"]]
        if twins:
            who = ", ".join(f"{x['py']['name']} (Python start {hms(x['py']['start'])})" for x in twins)
            rep.add("CHECK", "B", "same mouse", f"{who} is also {s['cage']}/{s['mouse']} on {this_date} - one of them "
                    f"has the wrong name (output name not changed for the next mouse?)",
                    ko=f"같은 날 {s['cage']} {s['mouse']} 세션 또 있음({twins[0]['py']['name']} "
                       f"{hms(twins[0]['py']['start'])[:5]})")
    # '_N' at the end: the task's logger adds it when the output file already exists
    m_n = re.fullmatch(r"(?P<base>.+?)_(?P<n>\d{1,3})", name)
    if m_n:
        base = m_n.group("base")
        hits = [(t, p) for t, p in logs_named(ctx, base) if str(p).lower() != str(s["py"]["path"]).lower()]
        where = ""
        if hits:
            t, p = hits[0]
            where = f"; {p}" + (f" started {datetime.fromtimestamp(t):%Y-%m-%d %H:%M:%S}" if t is not None else
                                " (empty or unreadable)")
        rep.add("INFO", "B", "name", f"'_{m_n.group('n')}' at the end: the task adds _N when '{base}.jsonl' already "
                f"exists in its output folder - a re-run after an aborted try, or the output name was not changed "
                f"(or the task script not saved with Ctrl+S) before this run{where}")
    if s.get("provisional"):
        return                             # the day number is only a guess: nothing more to compare

    # day number: already used by logs of an earlier date? the same as the other sessions of this date? the next day?
    day_n = s["day"]
    earlier = sorted((d, n) for d, n, _ in ctx.filed.get(day_n, []) if d < this_date)
    same = [x for x in ctx.sessions if x["status"] == "session" and x is not s
            and day_of(x["py"]["start"]) == this_date]
    differ = [x for x in same if x["day"] != day_n]
    counts = Counter([day_n] + [x["day"] for x in same]).most_common()
    majority = len(counts) == 1 or (counts[0][0] == day_n and counts[0][1] > counts[1][1])
    others = ", ".join(f"day {x['day']} ({x['py']['name']})" for x in differ[:4]) + (" ..." if len(differ) > 4 else "")
    prev = {k: max(x[0] for x in v if x[0] < this_date) for k, v in ctx.filed.items()
            if k != day_n and any(x[0] < this_date for x in v)}
    k = max(prev) if prev else None
    # sessions of this date (or logs of this date already filed in python/day<N>/) that use the same number
    same_num = [x for x in same if x["day"] == day_n] + \
        [n for d, n, t in ctx.filed.get(day_n, []) if d == this_date and abs(t - s["py"]["start"]) > 1.0]
    skip = k is not None and not earlier and day_n != k + 1 and not same_num and not (differ and not majority)
    notes = []
    if earlier:
        dates = sorted({d.isoformat() for d, _ in earlier})
        rep.add("CHECK", "B", "day number", f"name says day {day_n}, but python/day{day_n}/ already has {len(earlier)} "
                f"log(s) from {', '.join(dates)} (e.g. {earlier[0][1]}) - day number in the name not updated?",
                ko=f"day{day_n} 번호가 이전 날짜에 이미 사용됨")
    if differ and not majority:
        rep.add("CHECK", "B", "day number", f"name says day {day_n}, other session(s) of {this_date} say {others}",
                ko="같은 날 다른 세션과 day 번호 다름")
    elif differ:
        notes.append(f"other session(s) of {this_date} say {others} - flagged in their own check")
    if skip:
        rep.add("CHECK", "B", "day number", f"name says day {day_n}, but the last filed day before {this_date} is day "
                f"{k} ({prev[k]}) and no other session of {this_date} says day {day_n} - the next day would be day "
                f"{k + 1} (day number typed wrong?)", ko=f"day{day_n} 번호 건너뜀(다음은 day{k + 1})")
    if not earlier and (majority or not differ) and not skip:
        rep.add("PASS", "B", "day number", f"day {day_n}: not used on an earlier date"
                + (f", same as {len(same) - len(differ)} other session(s) of {this_date}" if same else ""))
    if k is not None and not skip:
        txt = f"last filed day before {this_date}: day {k} ({prev[k]})"
        if not earlier and day_n != k + 1:
            txt += f" - the name says day {day_n}, the next day would be {k + 1}"
        notes.append(txt)
    for n in notes:
        rep.note(n)


def check_protocol(rep, log, row, replay):
    """C: protocol numbers from the Python log. Returns the measured values."""
    ev, meta = log["events"], log["meta"] or {}
    tag = lambda t: [e for e in ev if e.get("tag") == t]          # noqa: E731
    out = {}
    paradigm = meta.get("paradigm", "")
    out["paradigm"] = paradigm
    want = [p for p in (row or {}).get("paradigm", "").split("|") if p]
    hdr = ", ".join(f"{k} {meta[k]}" for k in ("post_reward_delay", "reward_duration", "total_lockout_duration",
                                                "increment", "ratio") if k in meta)
    if not want:
        rep.add("INFO", "C", "paradigm", f"{paradigm or '-'} ({hdr})")
    elif paradigm in want:
        rep.add("PASS", "C", "paradigm", f"{paradigm} ({hdr})")
    else:
        rep.add("CHECK", "C", "paradigm", f"{paradigm or '-'}, expected {' or '.join(want)} ({hdr})",
                ko=f"패러다임 {paradigm or '-'}")

    # requirement sequence: req of trial k = start + int(k * step), as the task builds it
    ts = [(e.get("trial"), e.get("req")) for e in tag("trial_start") if e.get("req") is not None]
    out["trials"] = len(ts)
    r0, step = num(row, "req_start"), num(row, "req_step")
    exp = ", ".join(x for x in (f"start {r0:g}" if r0 is not None else "",
                                f"step {step:g}" if step is not None else "") if x)
    if not ts:
        n_ts = len(tag("trial_start"))
        text = f"no req in the {n_ts} trial_start record(s)" + (" (replay / passive log)" if replay else "")
        if r0 is None and step is None:
            rep.add("INFO", "C", "requirement", text)
        else:
            rep.add("CHECK", "C", "requirement", f"{text}; expected {exp}", ko="req 없음")
    else:
        reqs = [r for _, r in ts]
        diffs = sorted({round(b - a, 6) for a, b in zip(reqs, reqs[1:])})
        out["req_start"] = reqs[0]
        out["req_step"] = diffs[0] if len(diffs) == 1 else "/".join(f"{d:g}" for d in diffs)
        shown = ", ".join(f"{r:g}" for r in reqs[:4]) + (f" ... {reqs[-1]:g}" if len(reqs) > 4 else "")
        text = f"req {shown} ({len(ts)} trials)"
        if r0 is None and step is None:
            rep.add("INFO", "C", "requirement", text)
        else:
            start = r0 if r0 is not None else reqs[0]
            if step is None:                       # only the start is given
                wrong = [] if reqs[0] == r0 else [ts[0]]
            else:
                wrong = [(k, r) for k, r in ts if k is not None and r != start + int(k * step)]
            if wrong:
                k, r = wrong[0]
                rep.add("CHECK", "C", "requirement", f"{text}; expected {exp}: {len(wrong)} trial(s) differ "
                        f"(first: trial {k} req {r:g})", ko="요구 횟수 순서")
            else:
                rep.add("PASS", "C", "requirement", f"{text} = {exp}")

    # lock-out = reward_cmd -> next trial_start
    starts = [e["elapsed"] for e in tag("trial_start")]
    gaps = []
    for e in tag("reward_cmd"):
        nxt = next((t for t in starts if t > e["elapsed"]), None)
        if nxt is not None:
            gaps.append(nxt - e["elapsed"])
    want_l, tol_l = num(row, "lockout_s"), num(row, "lockout_tol_s") or 0.0
    log_lock = meta.get("total_lockout_duration")
    if gaps:
        med = float(np.median(gaps))
        out["lockout_median_s"] = round(med, 4)
        text = f"reward_cmd -> next trial_start median {med:.3f} s (min {min(gaps):.3f}, max {max(gaps):.3f}, n {len(gaps)})"
        hdr_txt = f" (log header: total_lockout_duration {log_lock:g} s)" if isinstance(log_lock, (int, float)) else ""
        if want_l is None:
            rep.add("INFO", "C", "lock-out", text + hdr_txt)
        elif abs(med - want_l) <= tol_l + 1e-9:
            rep.add("PASS", "C", "lock-out", f"{text} vs {want_l:g} +- {tol_l:g} s{hdr_txt}")
        else:
            rep.add("CHECK", "C", "lock-out", f"{text} vs {want_l:g} +- {tol_l:g} s{hdr_txt}", ko=f"잠금 {med:.2f} s")
    else:
        rep.add("INFO", "C", "lock-out", "not measured (no reward followed by another trial)")

    # pump
    n_inj, n_rew = len(tag("injector_run")), len(tag("reward_cmd"))
    out["injector_run"] = n_inj
    pump = (row or {}).get("pump", "").lower()
    text = f"injector_run {n_inj}, reward_cmd {n_rew}"
    if pump == "no":
        if n_inj == 0:
            rep.add("PASS", "C", "pump", f"{text} (valve only, as expected)")
        else:
            rep.add("CHECK", "C", "pump", f"{text} - expected 0 (valve only, no pump)", ko=f"펌프 injector_run {n_inj}회")
    elif pump == "yes":
        if n_inj == n_rew:
            rep.add("PASS", "C", "pump", f"{text} (pump with every reward)")
        else:
            rep.add("CHECK", "C", "pump", f"{text} - expected one injector_run per reward", ko="펌프 횟수")
    else:
        rep.add("INFO", "C", "pump", text)

    # session length and experiment_end
    end = next((e for e in reversed(ev) if e.get("tag") == "experiment_end"), None)
    t_meta = wall(meta) if meta else None
    t_end = wall(end) if end else (wall(ev[-1]) if ev else None)
    length = (t_end - t_meta) if (t_meta is not None and t_end is not None) else None
    out["session_s"] = round(length, 3) if length is not None else None
    out["experiment_end"] = bool(end)
    out["py_start_wall"], out["py_end_wall"] = t_meta, t_end
    want_s, tol_s = num(row, "session_s"), num(row, "session_tol_s") or 0.0
    if length is None:
        rep.add("INFO", "C", "session length", "not measured (no readable times in the log)")
    else:
        text = f"{length:.1f} s ({length / 60:.2f} min) from meta to {'experiment_end' if end else 'the LAST log line'}"
        if want_s is None:
            rep.add("INFO", "C", "session length", text)
        elif abs(length - want_s) <= tol_s:
            rep.add("PASS", "C", "session length", f"{text} vs {want_s:g} +- {tol_s:g} s")
        else:
            rep.add("CHECK", "C", "session length", f"{text} vs {want_s:g} +- {tol_s:g} s", ko=f"세션 길이 {length:.0f} s")
    if end:
        rep.add("PASS", "C", "experiment_end", f"present (completed_trials {end.get('completed_trials')}, "
                f"elapsed {end.get('elapsed')} s)")
    else:
        rep.add("CHECK", "C", "experiment_end", "MISSING - the task was killed, so 'T' may not have been sent and the "
                "camera trigger may still be running (see F after end)", ko="experiment_end 없음(강제 종료)")
    return out


def rewards_from_first_trial(events):
    """reward_cmd times (s) from the first trial_start: the clock the replay schedule uses."""
    t0 = next((e["elapsed"] for e in events if e.get("tag") == "trial_start" and "elapsed" in e), None)
    if t0 is None:
        return []
    return [e["elapsed"] - t0 for e in events if e.get("tag") == "reward_cmd" and "elapsed" in e]


def manifest_rows(project: Path):
    """{(cage_folder, mouse): row} of cohort_manifest.csv ({} when missing or unreadable)."""
    try:
        with open(project / "cohort_manifest.csv", newline="", encoding="utf-8-sig") as f:
            return {((r.get("cage_folder") or "").strip().lower(), (r.get("mouse") or "").strip().lower()):
                    {k: (v or "").strip() for k, v in r.items() if k} for r in csv.DictReader(f)}
    except (OSError, csv.Error):
        return {}


def check_replay(rep, log, s, ctx):
    """C, passive (replay) log: the log it replayed (replay_source in its meta) - found, another mouse, started on
    the same date and finished before this session, not itself a replay, with experiment_end - and the replayed
    rewards vs the source's: the same count, at the same times from the first trial_start."""
    meta = log["meta"] or {}
    src = str(meta.get("replay_source") or "").strip()
    if not src:
        rep.add("CHECK", "C", "replay source", "the log meta has no replay_source - the replayed log is unknown",
                ko="재생 원본 이름 없음")
        return
    me, t_me = Path(s["py"]["path"]), s["py"]["start"]
    name = Path(src.replace("\\", "/")).name
    found = [(log_start(p), p) for p in ([Path(src)] if Path(src).is_absolute() else []) + [me.parent / name]
             if p.is_file()]
    found += logs_named(ctx, Path(name).stem)
    seen, cands = {str(me.resolve()).lower()}, []
    for t, p in found:
        k = str(p.resolve()).lower()
        if k not in seen:
            seen.add(k)
            cands.append((t, p))
    if not cands:
        rep.add("CHECK", "C", "replay source", f"'{src}' (replay_source) not found in the task folder, python/day*/ or "
                "_tests_and_aborted/ - which log was replayed cannot be checked", ko="재생 원본 로그 못 찾음")
        return
    before = [(t, p) for t, p in cands if t is not None and t < t_me]
    t_src, p_src = max(before, key=lambda x: x[0]) if before else cands[0]   # the newest copy older than this run
    L = ss.read_python_log(p_src)
    sev, smeta = L["events"], (L["meta"] or {})
    s_par = str(smeta.get("paradigm", ""))
    s_end = next((e for e in reversed(sev) if e.get("tag") == "experiment_end"), None)
    s_len = s_end["elapsed"] if s_end else (sev[-1].get("elapsed") if sev else None)
    try:
        where = p_src.relative_to(ctx.project).as_posix()
    except ValueError:
        where = str(p_src)
    head = f"{name} = {where}" + (f" (started {datetime.fromtimestamp(t_src):%Y-%m-%d %H:%M:%S}, {s_par or '-'})"
                                   if t_src is not None else f" ({s_par or '-'})")
    bad = []                                   # (English, Korean)
    if "replay" in s_par.lower():
        bad.append(("it is itself a replay (passive) log", "원본도 재생 로그"))
    if t_src is None:
        bad.append(("its first line is not a meta line", "원본 첫 줄 이상"))
    else:
        if day_of(t_src) != day_of(t_me):
            bad.append((f"it started on {day_of(t_src)}, this session on {day_of(t_me)}",
                        f"원본 날짜 {day_of(t_src)}"))
        if t_src >= t_me:
            bad.append(("it started after this session", "원본이 더 늦게 시작"))
        elif isinstance(s_len, (int, float)) and t_src + s_len > t_me + 1.0:
            bad.append((f"it was still running when this session started ({hms(t_src + s_len)} vs {hms(t_me)})",
                        "원본이 끝나기 전에 시작"))
    if not s_end:
        bad.append(("it has no experiment_end (the active session crashed or was killed)", "원본 experiment_end 없음"))
    m_src = ds.parse_name(p_src.stem, ds.known_mice())     # 2026-11-19: the flexible name rule
    src_cm = (m_src.group("cage").lower(), m_src.group("mouse").lower()) if m_src else None
    if src_cm and src_cm == (str(s["cage"]).lower(), str(s["mouse"]).lower()):
        bad.append(("it is a log of this same mouse", "같은 마우스의 로그"))
    if bad:
        rep.add("CHECK", "C", "replay source", f"{head}: " + "; ".join(en for en, _ in bad),
                ko="재생 원본: " + ", ".join(k for _, k in bad))
    else:
        rep.add("PASS", "C", "replay source", f"{head}: another mouse, same date, finished before this session")

    # active/passive pairs of cohort_manifest.csv (checked only when they are filled in)
    man = manifest_rows(ctx.project)
    mine = man.get((str(s["cage"]).lower(), str(s["mouse"]).lower()), {})
    theirs = man.get(src_cm, {}) if src_cm else {}
    notes = []
    if mine.get("group") and mine["group"].lower() != "passive":
        notes.append(f"cohort_manifest.csv says this mouse is {mine['group']}")
    if theirs.get("group") and theirs["group"].lower() != "active":
        notes.append(f"cohort_manifest.csv says the source mouse is {theirs['group']}")
    if mine.get("pair_id") and theirs.get("pair_id") and mine["pair_id"] != theirs["pair_id"]:
        notes.append(f"pair_id {mine['pair_id']} (this mouse) vs {theirs['pair_id']} (source mouse)")
    if notes:
        rep.add("CHECK", "C", "replay pair", "; ".join(notes), ko="재생 짝(cohort_manifest) 불일치")
    elif mine.get("pair_id") and theirs.get("pair_id"):
        rep.add("PASS", "C", "replay pair", f"pair_id {mine['pair_id']} for both mice (cohort_manifest.csv)")
    else:
        rep.add("INFO", "C", "replay pair", "pairs not checked (group / pair_id not filled in cohort_manifest.csv)")

    # the replayed rewards: same number, same times from the first trial_start
    a, b = rewards_from_first_trial(log["events"]), rewards_from_first_trial(sev)
    if len(a) == len(b):
        rep.add("PASS", "C", "replay rewards", f"{len(a)} reward_cmd replayed = {len(b)} in the source")
    else:
        why = " (this session stopped before the end of the source?)" if len(a) < len(b) else ""
        rep.add("CHECK", "C", "replay rewards", f"{len(a)} reward_cmd replayed, the source has {len(b)}{why}",
                ko=f"재생 보상 {len(a)}회, 원본 {len(b)}회")
    d = [x - y for x, y in zip(a, b)]
    if d:
        big = max(d, key=abs)
        text = (f"replayed reward times - source times (from the first trial_start): max {big:+.3f} s, last "
                f"{d[-1]:+.3f} s (previous cohort: up to +0.24 s, about 2 ms per reward)")
        if abs(big) <= REPLAY_TIME_TOL_S:
            rep.add("PASS", "C", "replay timing", f"{text}; within {REPLAY_TIME_TOL_S:g} s")
        else:
            rep.add("CHECK", "C", "replay timing", f"{text}; more than {REPLAY_TIME_TOL_S:g} s",
                    ko=f"재생 보상 시각 차이 {big:+.2f} s")
    me_end = next((e for e in reversed(log["events"]) if e.get("tag") == "experiment_end"), None)
    if me_end and isinstance(s_len, (int, float)):
        rep.add("INFO", "C", "replay length", f"this session {me_end['elapsed']:.1f} s, the source {s_len:.1f} s "
                f"({me_end['elapsed'] - s_len:+.2f} s)")


def session_timeline(cam, cap_start, cap_end, clock, py_dur, has_end, sal_lag=None, during=()):
    """Where the session sits in the Saleae capture. t0/t1 = Python start / end put on the Saleae clock with the
    sync's clock line. The camera trigger runs from 'S' (Python start) to 'T' (sent right before experiment_end),
    so the pulse train that starts at the Python start shows the session on the Saleae clock whatever the PC
    clock did. Its start is compared with t0, its end with t1 and its length with the Python duration; when one
    of them is off by more than CLOCK_TOL_S (clock_off: the PC clock changed during the session), that train is
    the session window. A train that began before t0 is this session's 'S' only when it starts where the .sal
    says Python started (no .sal: at most STEP_MAX_S before t0, or by as much as a logged PC clock change);
    otherwise the trigger was already running before Python started (see F before start)."""
    a, b = clock
    T = SimpleNamespace(a=a, b=b, t0=a, t1=a + b * py_dur, train=None, first=None, last=None, n=0, d_end=None,
                        d_start=None, d_len=None, to_cap_end=False, end_known=False, start_is_s=False,
                        start_off=False, end_off=False, len_off=False, clock_off=False, sal_lag=sal_lag,
                        side="all", t_step=None, step_src="", step_between=None)
    T.win = (T.t0 - PAD_BEFORE, T.t1 + PAD_AFTER)
    tr = trains(cam)
    near = [(abs(cam[i] - T.t0), i, j) for i, j in tr if abs(cam[i] - T.t0) <= 1.0]
    holds = [(0.0, i, j) for i, j in tr if cam[i] <= T.t0 <= cam[j] + 1.0 / CAM_HZ]
    pick = min(near) if near else (holds[0] if holds else None)
    if pick is None:
        return T
    _, i, j = pick
    T.train, T.first, T.last, T.n = (i, j), float(cam[i]), float(cam[j]), j - i + 1
    T.d_end = T.last - T.t1
    T.d_start = T.first - T.t0
    T.to_cap_end = cap_end - T.last <= 3.0 / CAM_HZ
    T.end_known = bool(has_end and not T.to_cap_end)                 # the train stopped at 'T' inside the capture
    T.d_len = (T.last - T.first) - b * py_dur if T.end_known else None
    seen = T.first - cap_start > 2.0 / CAM_HZ                         # the train began inside the capture
    if seen and T.d_start < -CLOCK_TOL_S:                             # ... before the clock line's Python start
        if sal_lag is not None:
            seen = abs(T.first - sal_lag) <= SAL_LAG_TOL_S
        else:
            shift = sum(x["delta_ms"] for x in during) / 1000.0
            seen = bool(during and abs(T.d_start - shift) <= 1.5 / CAM_HZ) or -T.d_start <= STEP_MAX_S
    T.start_is_s = seen                                               # its first pulse is this session's 'S'
    T.end_off = bool(T.end_known and abs(T.d_end) > CLOCK_TOL_S)
    T.start_off = bool(seen and abs(T.d_start) > CLOCK_TOL_S)
    T.len_off = bool(seen and T.d_len is not None and abs(T.d_len) > CLOCK_TOL_S)
    T.clock_off = T.end_off or T.start_off or T.len_off
    if T.clock_off:
        T.win = (T.first - PAD_BEFORE, (T.last if T.end_known else T.t1) + PAD_AFTER)
    return T


CLOCK_NOTE = ("the PC clock changed during the session and the sync's single clock line fits neither its start nor "
              "its end (see F clock): Python-vs-Saleae matching is not reliable; counts from the Saleae are still valid")


def step_from_rewards(T, ev, rew_e):
    """No Windows event log entry: find the PC clock change from the rewards. Each Python reward_cmd (elapsed s) is
    put on the Saleae clock twice - with the sync's line, and with the clock of the other part of the session (the
    train start 'S' when the line fits the end; the line moved by the end difference when it fits the start) - and
    counts as before / after the change when a valve opening sits at only one of the two places.
    Returns (last reward before, first reward after the change), elapsed s; None where not found."""
    v = np.sort(ev.loc[(ev["event"] == "reward_valve") & ~(ev["duration_ms"] < qa.GLITCH_MS),
                       "t_saleae_s"].dropna().to_numpy())

    def hit(t, lo, hi):
        k = int(np.searchsorted(v, t + lo))
        return k < len(v) and v[k] <= t + hi

    side = []                              # in log order (Python times can repeat after a clock set back)
    for e in rew_e:
        on_line = hit(T.a + T.b * e, -0.005, 0.005)
        if T.side == "end":                # the other part's clock starts at 'S' (first pulse ~0.5 ms after it)
            on_other = hit(T.first + T.b * e, -0.005, 0.005)
        else:                              # the last pulse comes 0-35 ms before 'T'
            on_other = hit(T.a + T.b * e + T.d_end, -0.005, 0.040)
        side.append(None if on_line == on_other else ("after" if (on_line if T.side == "end" else on_other)
                                                      else "before"))
    bef = [k for k, x in enumerate(side) if x == "before"]
    aft = [k for k, x in enumerate(side) if x == "after"]
    if T.side == "end":                    # the first reward after the change that follows the last one before it
        nxt = [k for k in aft if not bef or k > bef[-1]]
        return (rew_e[bef[-1]] if bef else None), (rew_e[nxt[0]] if nxt else None)
    prv = [k for k in bef if not aft or k < aft[0]]
    return (rew_e[prv[-1]] if prv else None), (rew_e[aft[0]] if aft else None)


def clock_side(T, during, t_start, ev, rew_e=()):
    """After a PC clock change the sync's single clock line fits only one part of the session. Sets T.side: 'all'
    (no change seen), 'start' (the line fits before the change), 'end' (it fits after it) or '' (neither: reward
    and lick matching are not judged), and T.t_step / T.step_src: the Python elapsed time up to which ('start') or
    from which ('end') the matching is judged (from the Windows event log, else from the rewards)."""
    if not (T.clock_off or during):
        return
    start_ok = abs(T.d_start) <= CLOCK_TOL_S if (T.train is not None and T.start_is_s) else None
    end_ok = abs(T.d_end) <= CLOCK_TOL_S if (T.train is not None and T.end_known) else None
    if start_ok is False and end_ok is False:
        T.side = ""
        return
    if start_ok is not False and end_ok is not False and not during:
        return                             # only the length is off (both ends within CLOCK_TOL_S): judge it all
    T.side = "end" if start_ok is False else "start"
    T.step_src = "the PC clock change"
    if during and t_start is not None:
        if T.side == "start":
            T.t_step = min(min(x["old"], x["new"]) for x in during) - t_start
        else:
            T.t_step = max(max(x["old"], x["new"]) for x in during) - t_start + STEP_MARGIN_S
    elif T.train is not None:
        T.step_between = step_from_rewards(T, ev, rew_e)
        last_before, first_after = T.step_between
        if T.side == "start" and last_before is not None:
            T.t_step = last_before - max(0.0, T.d_end)        # later Python times can repeat earlier ones
        elif T.side == "end" and first_after is not None:
            T.t_step = first_after + max(0.0, -T.d_start) + STEP_MARGIN_S
    if T.t_step is None and T.clock_off:
        T.side = ""


def check_valve(rep, ev, n_cmd, sync, row, T, residual):
    """D: valve openings on the Saleae reward channel inside the session, and the sync's reward matching."""
    out = {}
    v = ev[(ev["event"] == "reward_valve") & (ev["t_saleae_s"] >= T.win[0]) & (ev["t_saleae_s"] <= T.win[1])]
    short = v["duration_ms"] < qa.GLITCH_MS
    glitch, opens = v[short], v[~short]                      # NaN (no falling edge) counts as an opening
    d = opens["duration_ms"].dropna()
    out.update(valve_openings=len(opens), valve_glitches=len(glitch), rewards_matched=sync["rewards_matched"])
    want, tol = num(row, "valve_ms"), num(row, "valve_tol_ms") or 0.0
    prev = f"(previous cohort {PREV_VALVE_MS:.0f} ms)"
    if len(d):
        med, lo, hi = float(d.median()), float(d.min()), float(d.max())
        out.update(valve_median_ms=round(med, 2), valve_min_ms=round(lo, 2), valve_max_ms=round(hi, 2))
        text = f"open time median {med:.1f} ms (min {lo:.1f}, max {hi:.1f}), {len(d)} openings"
        if want is None:
            rep.add("INFO", "D", "valve time", f"{text}; no expectation for this day {prev}")
        else:
            outside = int(((d < want - tol) | (d > want + tol)).sum())
            if outside == 0:
                rep.add("PASS", "D", "valve time", f"{text} vs {want:g} +- {tol:g} ms {prev}")
            else:
                rep.add("CHECK", "D", "valve time", f"{text} vs {want:g} +- {tol:g} ms: {outside} opening(s) outside "
                        f"{prev}", ko=f"밸브 {med:.0f} ms")
    else:
        rep.add("CHECK" if n_cmd else "INFO", "D", "valve time", "no valve opening on the Saleae reward channel",
                ko="밸브 펄스 없음")
    n_open = len(opens)
    glitch_txt = f" ({len(glitch)} pulse(s) < {qa.GLITCH_MS:g} ms ignored as glitches)" if len(glitch) else ""
    if n_open == n_cmd:
        rep.add("PASS", "D", "reward count", f"Python reward_cmd {n_cmd} = Saleae openings {n_open}{glitch_txt}")
    else:
        rep.add("CHECK", "D", "reward count", f"Python reward_cmd {n_cmd} vs Saleae openings {n_open}{glitch_txt}",
                ko=f"보상 수 불일치(Python {n_cmd}, Saleae {n_open})")

    # did the sync pair every reward_cmd with an opening, and every opening with a reward_cmd?
    lost = ev[(ev["event"] == "reward_cmd") & (ev["source"] == "python")]       # reward_cmd without an opening
    g_t = glitch.loc[glitch["python_elapsed_s"].notna(), "t_saleae_s"].to_numpy()
    no_cmd = [t for t in opens.loc[opens["python_elapsed_s"].isna(), "t_saleae_s"]
              if not any(0 <= t - g < 0.05 for g in g_t)]     # an opening right after a matched glitch is matched
    text = f"the sync matched {sync['rewards_matched']} of {n_cmd} reward_cmd to an opening{residual}"
    if not len(lost) and not no_cmd:
        rep.add("PASS", "D", "reward match", text)
        return out
    text += f"; {len(lost)} reward_cmd without an opening, {len(no_cmd)} opening(s) without a reward_cmd"
    if T.side == "":
        rep.add("INFO", "D", "reward match", f"{text} - not judged: {CLOCK_NOTE}")
        return out
    if T.t_step is not None:               # judge only where the sync's clock line fits (before / after the change)
        t_lim = T.a + T.b * T.t_step
        if T.side == "end":
            n_in = int((lost["python_elapsed_s"] >= T.t_step).sum()) + sum(1 for t in no_cmd if t >= t_lim)
        else:
            n_in = int((lost["python_elapsed_s"] <= T.t_step).sum()) + sum(1 for t in no_cmd if t <= t_lim)
        if n_in == 0 and T.side == "end":
            rep.add("INFO", "D", "reward match", f"{text} - all before {T.step_src} (elapsed {T.t_step:.0f} s), where "
                    f"the sync's single clock line does not fit: matching there is not reliable; counts from the "
                    f"Saleae are still valid (see F clock)")
            return out
        if n_in == 0:
            rep.add("INFO", "D", "reward match", f"{text} - all after {T.step_src} (elapsed {T.t_step:.0f} s), where "
                    f"the sync's single clock line no longer fits (see F clock)")
            return out
        text += f" ({n_in} of them {'after' if T.side == 'end' else 'before'} {T.step_src})"
    rep.add("CHECK", "D", "reward match", text, ko=f"보상 매칭 안 됨({len(lost) + len(no_cmd)}개)")
    return out


def check_licks(rep, ev, log, sync, is_pr, T, replay):
    """E: Saleae lick onsets vs Python licks; Python-only licks that the task counted (valid=True)."""
    lk = ev[(ev["event"] == "lick") & (ev["source"] == "saleae")]
    n_sal = int(((lk["t_saleae_s"] >= T.win[0]) & (lk["t_saleae_s"] <= T.win[1])).sum())
    if replay:
        n_obs = sum(1 for e in log["events"] if e.get("tag") == "observed_lick")
        rep.add("INFO", "E", "licks", f"Saleae onsets {n_sal} in the session, Python observed_lick {n_obs} (replay "
                f"log: the sync matches only 'lick' records, so only the counts are shown)")
        rep.add("INFO", "E", "Python-only", "not checked (replay log: observed licks do not move the task on)")
        return dict(licks_saleae=n_sal, licks_python=n_obs)
    licks = [e for e in log["events"] if e.get("tag") == "lick"]
    n_valid = sum(1 for e in licks if e.get("valid") is True)
    po = ev[(ev["event"] == "lick") & (ev["source"] == "python")]          # Python licks the sync did not match
    po_valid = int((po["valid_or_hit"] == 1).sum())
    out = dict(licks_saleae=n_sal, licks_python=len(licks), licks_python_valid=n_valid,
               python_only_licks=len(po), python_only_valid=po_valid)
    rep.add("INFO", "E", "licks", f"Saleae onsets {n_sal} in the session, Python {len(licks)} "
            f"(valid {n_valid}, valid=False {len(licks) - n_valid}), matched {sync['licks_matched']}; "
            f"Saleae-only {n_sal - sync['licks_matched']}, Python-only {len(po)}")
    if T.side == "":                       # the sync's clock line fits neither end: the matching cannot be judged
        rep.add("INFO", "E", "Python-only", f"not judged: {len(po)} Python-only licks, {po_valid} of them "
                f"valid=True; {CLOCK_NOTE}")
        return out
    if T.t_step is not None:
        # judge only the side of the clock change where the sync's single clock line fits
        after = T.side == "end"
        sel = po[po["python_elapsed_s"] >= T.t_step] if after else po[po["python_elapsed_s"] <= T.t_step]
        n_bad = int((sel["valid_or_hit"] == 1).sum())
        rest = len(po) - len(sel)
        word = "after" if after else "before"
        text = (f"{n_bad} of the {len(sel)} Python-only licks {word} {T.step_src} (elapsed {T.t_step:.0f} s) were "
                f"valid=True" + (f"; the {rest} {'before' if after else 'after'} it are not judged (the sync's "
                                 f"matching is off there, see F clock)" if rest else ""))
        # an early clock change can also push the sync's lick matcher off on the judged side; then nearly every
        # Python lick there is 'Python-only'. Real Arduino re-triggers are about 1 % of licks, so a large unmatched
        # share is a matching failure, not licks the task counted (added 2026-11-16)
        el = [float(e.get("elapsed", float("nan"))) for e in licks]
        n_side = sum(1 for x in el if (x >= T.t_step if after else x <= T.t_step))
        if n_side and len(sel) > MATCH_FAIL_FRAC * n_side:
            rep.add("INFO", "E", "Python-only", f"not judged: {len(sel)} of the {n_side} Python licks {word} "
                    f"{T.step_src} are unmatched (real re-triggers are about 1 %), so the sync's lick matching failed "
                    f"around the PC clock change; {n_bad} of them valid=True. Counts are still valid: Saleae onsets "
                    f"{n_sal}, Python licks {len(licks)} (valid {n_valid}). The evening sync report shows the same.")
            return out
    else:
        n_bad = po_valid
        text = f"{po_valid} of the {len(po)} Python-only licks were valid=True"
    text += " (counted by the task, not seen on the Saleae lick line)"
    if not is_pr:
        rep.add("INFO", "E", "Python-only", text + " - checked on PR days only")
    elif n_bad == 0:
        rep.add("PASS", "E", "Python-only", text)
    else:
        rep.add("CHECK", "E", "Python-only", text + " - on a PR day they move the ratio on (Arduino re-trigger while "
                "the valve is open?)", ko=f"Saleae에 없는 유효 lick {n_bad}개")
    return out


def check_clock(rep, T, py_dur, has_end, during, log_ok, py_start_wall):
    """F clock: does the session's camera train ('S' to 'T') start at the Python start, end at experiment_end and
    last as long as Python says (PC clock)?"""
    shift = -sum(x["delta_ms"] for x in during) / 1000.0
    explain = lambda diff: ("the PC clock change below explains it" if during and abs(diff - shift) <= 1.5 / CAM_HZ  # noqa: E731
                            else "a PC clock change is logged below, but it does not explain all of it" if during
                            else "PC clock changed during the session? none found in the Windows event log")
    if T.train is None:
        if during:
            rep.add("CHECK", "F", "clock", "the Windows event log shows a PC clock change during the session (no "
                    "camera pulse train starts at the Python start to compare with)", ko="PC 시계 변경")
        else:
            rep.add("INFO", "F", "clock", "not checked: no camera pulse train starts at the Python start")
    elif T.clock_off and T.end_off and not T.start_off:          # the sync's clock line fits the start
        rep.add("CHECK", "F", "clock", f"Python and Saleae disagree by {T.d_end:+.2f} s at the end: the camera "
                f"train ('S' to 'T') lasts {T.last - T.first:.2f} s ({T.n} pulses), Python start -> experiment_end "
                f"{py_dur:.2f} s by the PC clock ({explain(T.d_end)})", ko=f"PC 시계 불일치 {T.d_end:+.2f} s")
    elif T.clock_off:
        diff = T.d_len if T.d_len is not None else -T.d_start
        if T.d_len is not None:
            body = (f"the camera train ('S' to 'T') lasts {T.last - T.first:.2f} s ({T.n} pulses), Python start -> "
                    f"experiment_end {py_dur:.2f} s by the PC clock; on the sync's clock line the train starts "
                    f"{T.d_start * 1000:+.0f} ms from the Python start and ends {T.d_end * 1000:+.0f} ms from "
                    f"experiment_end")
        else:
            body = (f"on the sync's clock line the session's camera train starts {T.d_start * 1000:+.0f} ms from the "
                    f"Python start (its length cannot be compared: "
                    f"{'no experiment_end' if not has_end else 'it runs until the Saleae capture stopped'})")
        rep.add("CHECK", "F", "clock", f"Python and Saleae disagree by {diff:+.2f} s: {body} ({explain(diff)})",
                ko=f"PC 시계 불일치 {diff:+.2f} s")
    elif during:
        rep.add("CHECK", "F", "clock", "the Windows event log shows a PC clock change during the session, but the "
                f"camera train ends {T.d_end * 1000:+.0f} ms from experiment_end", ko="PC 시계 변경")
    elif not has_end:
        rep.add("INFO", "F", "clock", "not checked: no experiment_end, so 'T' was probably not sent")
    elif abs(T.d_end) > CLOCK_TOL_S:
        rep.add("INFO", "F", "clock", "not checked: the camera train runs until the Saleae capture stopped "
                "(see below)")
    elif not T.start_is_s and abs(T.d_start) > CLOCK_TOL_S:
        rep.add("INFO", "F", "clock", f"only the end compared: the camera train ends {T.d_end * 1000:+.0f} ms from "
                f"experiment_end, but it does not start at 'S' ({T.d_start:+.2f} s from the Python start: trigger "
                f"already running, or the capture started late - see F before start / first pulse), so its start "
                f"and length are not compared")
    else:
        rep.add("PASS", "F", "clock", f"camera train ('S' to 'T') {T.last - T.first:.2f} s ({T.n} pulses) ends "
                f"{T.d_end * 1000:+.0f} ms from experiment_end - PC and Saleae clocks agree")
    for x in during:
        rep.note("Windows event log: " + step_text(x, py_start_wall))
    if T.clock_off or during:
        if not log_ok:
            rep.note("Windows event log: not readable here")
        if T.clock_off:
            rep.note("frames, valve and lick counts below use this camera train as the session window")
        if T.step_between and not during:
            lb, fa = T.step_between
            where = (f"between the rewards at elapsed {lb:.0f} s and {fa:.0f} s" if lb is not None and fa is not None
                     else f"after the reward at elapsed {lb:.0f} s" if lb is not None
                     else f"before the reward at elapsed {fa:.0f} s" if fa is not None
                     else "nowhere (no reward fits one clock only)")
            rep.note(f"no Windows event log entry: matching the rewards with both clocks puts the PC clock change "
                     f"{where}")
        if T.side == "end" and T.t_step is not None:
            rep.note(f"the sync puts all Python events on one clock line, and it fits the session after {T.step_src}, "
                     f"so before it its reward/lick matching is off - the evening sync too; D reward match and "
                     f"E Python-only below are judged from elapsed {T.t_step:.0f} s on")
        elif T.side == "":
            rep.note("the sync's single clock line fits neither the start nor the end of the session: D reward match "
                     "and E Python-only below are not judged - the evening sync's matching is off too")
        elif T.t_step is not None and during:
            rep.note(f"the sync puts all Python events on one clock line, so after the PC clock change (elapsed "
                     f"{T.t_step:.0f} s) its reward/lick matching is off - the evening sync too")
        elif T.t_step is not None:
            rep.note(f"the clocks agree at least until {T.step_src} (elapsed {T.t_step:.0f} s); after that the sync's "
                     f"reward/lick matching may be off (one clock line for the whole session) - the evening sync too")
        rep.note("the task's times come from the PC clock (time.time()), which Windows time sync can step; whether "
                 "to change that (task code or PC time settings) is outside this check")


def describe_train(cam, i, j, T, ctx, s, py_start_wall):
    """Another camera pulse train in the capture: (text, Python run that explains it or None)."""
    n, t_a, t_b = j - i + 1, float(cam[i]), float(cam[j])
    t_end = T.last if (T.clock_off and T.end_known) else T.t1
    t_beg = T.first if T.start_off else T.t0
    where = f"starts {t_a - t_end:.1f} s after the end" if t_a > t_end else f"ends {t_beg - t_b:.1f} s before the start"
    for r in ctx.runs:
        if abs(r["start"] - s["py"]["start"]) <= 1.0 or py_start_wall is None:
            continue
        expect = T.a + T.b * (r["start"] - py_start_wall)
        if abs(expect - t_a) <= RUN_MATCH_S and abs((r["end"] - r["start"]) - (t_b - t_a)) <= RUN_MATCH_S:
            return (f"another camera pulse train of {n} pulses {where}: {r['name']} ran {hms(r['start'])}-"
                    f"{hms(r['end'])} into the same Saleae capture"), r
    return f"another camera pulse train of {n} pulses ({t_b - t_a:.1f} s) {where}", None


def check_camera(rep, s, ctx, sync, cam, cap_end, T, has_end, py_start_wall):
    """F (with the Saleae capture): camera pulses vs video frames, first pulse, pulses before / after the session."""
    frames, pulses = sync["video_frames"], sync["camera_pulses"]
    before = int((cam < T.win[0]).sum())
    after = int((cam > T.win[1]).sum())
    # with a PC clock change the session's camera train is the reference ('S' to 'T', or to the Python end when
    # 'T' is not in the capture)
    train_n = (T.n if T.end_known else int(((cam >= T.win[0]) & (cam <= T.win[1])).sum())) if T.clock_off else None
    out = dict(camera_pulses=train_n if T.clock_off else pulses, video_frames=frames,
               pulses_before_start=before, pulses_after_end=after)
    capture_short = T.train is not None and T.to_cap_end and T.d_end < -CLOCK_TOL_S
    if not s.get("video"):
        rep.add("INFO", "F", "frames", f"no video (see A); {train_n if T.clock_off else pulses} camera pulses in the "
                f"session")
    elif frames is None:
        rep.add("CHECK", "F", "frames", f"video frames could not be counted ({s['video']['path']}); {pulses} camera "
                f"pulses in the session", ko="비디오 프레임 수 모름")
    elif capture_short:
        miss = int(round((T.t1 - T.last) * CAM_HZ))
        diff = frames - (pulses + miss)
        text = (f"video {frames} frames vs {pulses} pulses recorded + about {miss} after the Saleae capture stopped "
                f"-> {diff:+d} (estimate)")
        if abs(diff) <= FRAME_EST_TOL:
            rep.add("INFO", "F", "frames", text)
        else:
            rep.add("CHECK", "F", "frames", text, ko=f"비디오 프레임 약 {abs(diff)}개 {'부족' if diff < 0 else '초과'}")
    else:
        ref, what = (train_n, "pulses in the session's camera train ('S' to 'T')" if T.end_known else
                     "pulses in the session's camera train (from 'S' to the Python end)") if T.clock_off else \
            (pulses, "camera pulses in the session")
        diff = frames - ref
        out["frames_minus_pulses"] = diff
        if diff == 0:
            rep.add("PASS", "F", "frames", f"video {frames} frames = {ref} {what}")
        else:
            how = f"{-diff} frames SHORT" if diff < 0 else f"{diff} EXTRA frames"
            rep.add("CHECK", "F", "frames", f"video {frames} frames vs {ref} {what} -> {how}",
                    ko=f"비디오 프레임 {abs(diff)}개 {'부족' if diff < 0 else '초과'}")
        if T.clock_off and pulses != ref:
            rep.note(f"the sync's clock-mapped window holds {pulses} pulses (video {frames - pulses:+d}); the "
                     f"train is the right reference here (see F clock)")
        elif diff and (before or after):
            rep.note(f"whole capture: {len(cam)} pulses ({before} before Python start, {after} after its end) "
                     f"-> video {frames - len(cam):+d} frames")

    in_sess = cam[(cam >= T.win[0]) & (cam <= T.win[1])]
    if len(in_sess) and T.start_off:
        rep.add("INFO", "F", "first pulse", f"not compared: the sync's clock line does not fit the start of the "
                f"session (see F clock); the session's camera train starts {T.first:.3f} s into the capture"
                + (f", where the .sal says Python started ({T.sal_lag:.2f} s; that stamp is good to about 0.2 s)"
                   if T.sal_lag is not None else ""))
    elif len(in_sess):
        first_ms = (in_sess[0] - T.t0) * 1000
        out["first_pulse_ms"] = round(first_ms, 2)
        first_ms = 0.0 if abs(first_ms) < 0.05 else first_ms           # no '-0.0'
        text = f"first camera pulse {first_ms:+.1f} ms from Python start ({T.t0:.3f} s into the capture)"
        if abs(first_ms) <= FIRST_PULSE_TOL_MS:
            rep.add("PASS", "F", "first pulse", text)
        else:
            rep.add("CHECK", "F", "first pulse", f"{text} - more than {FIRST_PULSE_TOL_MS:g} ms (capture started late?)",
                    ko="첫 카메라 펄스 시각")

    # pulses before the session
    tr = trains(cam)
    if T.train is not None and T.first < T.t0 - 1.0 and not T.start_off:
        rep.add("CHECK", "F", "before start", f"the camera trigger was already running {T.t0 - T.first:.1f} s before "
                f"Python started ('T' of an earlier run not received?)", ko="시작 전부터 트리거 동작")
    if T.train is None:
        if before:
            rep.add("INFO", "F", "before start", f"{before} camera pulses before Python start in this capture "
                    f"(another run in the same capture, or the trigger was already running)")
    else:
        for i, j in tr:
            if j < T.train[0]:
                text, _ = describe_train(cam, i, j, T, ctx, s, py_start_wall)
                rep.add("INFO", "F", "before start", text)

    # pulses after the session: the session's own train running on, or other trains
    if T.train is None:
        if after:
            rep.add("CHECK", "F", "after end", f"{after} camera pulses after Python ended", ko=f"종료 후 트리거 {after}펄스")
        else:
            rep.add("PASS", "F", "after end", "no camera pulse after Python ended")
        return out
    i0, j0 = T.train
    own_after = int((cam[i0:j0 + 1] > T.t1 + PAD_AFTER).sum())      # the session's own train after the end
    if T.clock_off and T.end_known:
        rep.add("PASS", "F", "after end", "the camera train stopped at 'T' (experiment_end)" +
                (f"; its last {own_after} pulses lie after the clock-mapped end only because of the clock difference"
                 if own_after else ""))
    elif capture_short:
        rep.add("CHECK", "F", "capture end", f"the Saleae capture stopped {T.t1 - T.last:.1f} s before Python ended "
                f"- the end of the session is not recorded", ko=f"Saleae가 과제 종료 {T.t1 - T.last:.0f} s 전에 멈춤")
    elif own_after:
        last, tail = T.last - T.t1, cap_end - T.last
        why = "'T' was sent with experiment_end, but the trigger ran on" if has_end else \
            "the task was killed (no experiment_end), so 'T' was probably not sent"
        if T.to_cap_end:
            how = (f"the trigger was still running when the Saleae capture stopped ({last:.1f} s after the end) - "
                   f"it may still be running: stop it (reset the Arduino) before the next recording")
        else:
            how = f"the trigger stopped {last:.2f} s after the end ({tail:.1f} s before the capture stopped)"
        rep.add("CHECK", "F", "after end", f"{own_after} camera pulses after Python ended - {how}; {why}",
                ko=f"종료 후 트리거 {own_after}펄스")
    else:
        rep.add("PASS", "F", "after end", "no camera pulse after Python ended")
    for i, j in tr:
        if i > j0:
            text, run = describe_train(cam, i, j, T, ctx, s, py_start_wall)
            if run:
                rep.add("INFO", "F", "after end", text)
            else:
                rep.add("CHECK", "F", "after end", text + " and no Python run of this day explains it (camera trigger "
                        "started without the task?)", ko=f"종료 후 다른 펄스 {j - i + 1}개")
    return out


def check_frames_estimate(rep, ctx, vid_p, py_dur, during):
    """F without the Saleae capture: video frames vs Python duration x 30.26 Hz (an estimate)."""
    v = ss.probe_video(vid_p, ctx.ffmpeg)
    if v.get("frames") is None:
        rep.add("CHECK", "F", "frames", f"video frames could not be counted ({v.get('error')})", ko="비디오 프레임 수 모름")
        return {}
    shift = -sum(x["delta_ms"] for x in during) / 1000.0           # a PC clock set back makes Python's time short
    real = py_dur + shift
    exp = real * CAM_HZ
    diff = v["frames"] - int(round(exp))
    text = (f"video {v['frames']} frames vs about {exp:.0f} expected ({real:.1f} s x {CAM_HZ:g} Hz"
            f"{', PC clock change included' if shift else ''}) -> {diff:+d} (estimate without the Saleae CSV)")
    if abs(diff) <= FRAME_EST_TOL:
        rep.add("PASS", "F", "frames", text)
    else:
        rep.add("CHECK", "F", "frames", text, ko=f"비디오 프레임 약 {abs(diff)}개 {'부족' if diff < 0 else '초과'}(추정)")
    return dict(video_frames=v["frames"])


def check_video_file(rep, s, ctx, vid_p, py_end_wall, steps_after):
    """F from the files only: video closed before the task ended; other Python runs while it was recorded."""
    out = {}
    v_end = vid_p.stat().st_mtime
    created = (s.get("video") or {}).get("created")
    if py_end_wall is not None:
        dv = v_end - py_end_wall
        out["video_closed_minus_end_s"] = round(dv, 2)
        text = f"video last written {hms(v_end)}, Python ended {hms(py_end_wall)} ({dv:+.1f} s)"
        if dv < 0:
            rep.add("CHECK", "F", "video closed", f"{text} - camera stopped {-dv:.1f} s (about {-dv * CAM_HZ:.0f} "
                    f"frames) before the task ended", ko=f"카메라가 과제 종료 {-dv:.1f} s 전에 멈춤")
        else:
            rep.add("PASS", "F", "video closed", text)
        for x in steps_after:
            rep.note(f"Windows event log: {step_text(x)}, after the task ended - this comparison is off by that much")
    if created is not None:
        rep.note(f"video file created {created - s['py']['start']:+.1f} s from Python start "
                 f"(IC Capture stamps whole seconds)")
    v_start = min(s["py"]["start"], created) if created is not None else s["py"]["start"]
    hits = [r for r in ctx.runs if abs(r["start"] - s["py"]["start"]) > 1.0
            and r["start"] < v_end and r["end"] > v_start]
    out["other_runs_in_video"] = "; ".join(r["name"] for r in hits)
    for r in hits:
        d = r["end"] - r["start"]
        rep.add("CHECK", "F", "other runs", f"{r['name']} ran {hms(r['start'])}-{hms(r['end'])} ({d:.1f} s, about "
                f"{d * CAM_HZ:.0f} pulses) while this video was still being written (until {hms(v_end)})",
                ko=f"다른 실행이 같은 비디오에 녹화됨({r['name']})")
    if not hits:
        rep.add("PASS", "F", "other runs", f"no other Python run between {hms(v_start)} and {hms(v_end)}")
    return out


def check_saleae(rep, sal_p, row):
    """G: Saleae settings from the .sal meta.json."""
    out = {}
    try:
        st = read_sal_settings(sal_p)
    except Exception as exc:
        rep.add("CHECK", "G", "settings", f"cannot read meta.json of {sal_p.name}: {exc}", ko=".sal 설정 읽기 실패")
        return out
    rate = st["rate"]
    out.update(sample_rate_mhz=rate / 1e6 if rate else None, threshold=st["threshold"])
    rep.add("INFO", "G", "rate/threshold", f"{rate / 1e6 if rate else 0:g} MHz, {st['threshold'] or '-'} "
            f"(previous cohort: {PREV_SALEAE})")
    wrong = [f"Ch{c} '{st['names'].get(c, '-')}' (expected '{n}')" for c, n in CHANNELS.items()
             if st["names"].get(c, "").lower() != n]
    off = [f"Ch{c} not enabled" for c in CHANNELS if c not in st["enabled"]]
    if wrong or off:
        rep.add("CHECK", "G", "channels", "; ".join(wrong + off), ko="Saleae 채널 이름/켜짐")
    else:
        rep.add("PASS", "G", "channels", "Ch4 reward, Ch6 camera, Ch7 lick (enabled)")
    is_1us = lambda w: bool(w) and abs(w - 1e-6) < 1e-9                       # noqa: E731
    have = ", ".join(f"Ch{c} {w * 1e6:g} us" for c, w in sorted(st["glitch"].items()) if w) or "none"
    out["glitch_1us"] = ";".join(str(c) for c, w in sorted(st["glitch"].items()) if is_1us(w))
    need = [int(c) for c in ((row or {}).get("glitch_1us_channels") or "").replace(",", ";").split(";") if c.strip()]
    if not need:
        rep.add("INFO", "G", "glitch filter", f"{have}; no expectation for this day")
    else:
        miss = [c for c in need if not is_1us(st["glitch"].get(c))]
        want = ", ".join(f"Ch{c}" for c in need)
        if miss:
            rep.add("CHECK", "G", "glitch filter", f"{have}; 1 us expected on {want} - missing on "
                    + ", ".join(f"Ch{c}" for c in miss), ko="글리치 필터 없음(" + ",".join(f"Ch{c}" for c in miss) + ")")
        else:
            rep.add("PASS", "G", "glitch filter", f"{have} (1 us on {want})")
    return out


def check_sync(rep, sync, sync_err, replay, clock_bad):
    """H: status and warnings of the evening sync."""
    if not sync:
        if sync_err:
            rep.add("CHECK", "H", "sync", f"FAIL - {sync_err}", ko="sync 실패")
        else:
            rep.add("INFO", "H", "sync", "not run (no Saleae CSV)")
        return
    warns = [w for w in (sync.get("problems") or "").split(" | ") if w]
    expected = [w for w in warns if replay and ("licks:" in w or "lick delay" in w)]
    rest = [w for w in warns if w not in expected]
    if sync["status"] == "OK":
        rep.add("PASS", "H", "sync", "OK (no warnings)")
    elif not rest:
        rep.add("INFO", "H", "sync", f"{sync['status']} only because a replay log has no 'lick' records (the sync "
                f"matches 'lick' only)")
    else:
        why = sorted({ko for key, ko in SYNC_WARN_KO if any(key in w for w in rest)})
        if clock_bad:
            why.append("PC 시계 영향")
        rep.add("CHECK", "H", "sync", f"{sync['status']} ({len(warns)} warning(s))",
                ko=f"sync {sync['status']}" + (f"({', '.join(why)})" if why else ""))
    for w in warns:
        rep.note(w)
    if clock_bad and rest:
        rep.note("these numbers come from the sync's single clock line, which does not fit the whole session "
                 "(see F clock) - the evening sync will show the same")


# ----------------------------------------------------------------------------- one session
def output_paths(day_dir: Path, name, start_txt, csv_rows):
    """Report file and sync folder of a session: <name>_check.txt and <name>/, or <name>_<date> (then
    <name>_<date>_<time>) when that name already holds the report of a run with another start (a log name used
    again on another date must not replace the earlier date's report)."""
    for tag in ("", f"_{start_txt[:10]}", f"_{start_txt[:10]}_{start_txt[11:].replace(':', '')}"):
        rp = day_dir / f"{name}{tag}_check.txt"
        owner = report_start(rp)
        if owner is None:
            owner = next((r.get("start") for r in csv_rows
                          if str(r.get("report") or "").lower() == str(rp).lower()), None)
        if owner is None or owner == start_txt:
            break
    return rp, day_dir / f"{name}{tag}"


def check_session(s, ctx):
    """Check one paired session; writes the report. Returns (overall, Korean items, csv row)."""
    rep = Report()
    py = s["py"]
    name, day_n = py["name"], s.get("day")
    out_root = ctx.project / "analysis" / "session_checks"
    day_dir = out_root / (f"day{day_n}" if day_n is not None else "day_unknown")       # provisional, day unknown
    start = datetime.fromtimestamp(py["start"])
    report_path, sync_dir = output_paths(day_dir, name, f"{start:%Y-%m-%d %H:%M:%S}", ctx.csv_rows)
    log = ss.read_python_log(py["path"])
    ev = log["events"]
    paradigm = (log["meta"] or {}).get("paradigm", "")
    row = protocol_for(ctx.protocol, day_n, paradigm)          # days 6-10: the active or the passive row
    tags = {e.get("tag") for e in ev}
    replay = "replay" in paradigm.lower() or ("observed_lick" in tags and "lick" not in tags)
    is_pr = bool(row and row.get("phase", "").upper().startswith("PR")) or "progressive" in paradigm.lower()

    sal_p, csv_p, vid_p = check_files(rep, s, ctx)                         # A
    check_name(rep, s, ctx, row)                                            # B
    vals = check_protocol(rep, log, row, replay)                            # C
    if replay:
        check_replay(rep, log, s, ctx)                                      # C: the log this passive session replayed
    n_cmd = sum(1 for e in ev if e.get("tag") == "reward_cmd")
    vals["reward_cmd"] = n_cmd
    end = next((e for e in reversed(ev) if e.get("tag") == "experiment_end"), None)
    py_dur = end["elapsed"] if end else (ev[-1]["elapsed"] if ev else 0.0)
    t_start, t_end = vals.get("py_start_wall"), vals.get("py_end_wall")

    # PC clock changes (Windows event log) during the session, or after it and before the video was closed
    steps, during, after = None, [], []
    if t_start is not None and t_end is not None:
        v_end = vid_p.stat().st_mtime if vid_p else t_end
        steps = pc_clock_steps(t_start - 5, max(t_end, v_end) + 5)
        during = [x for x in steps or [] if max(x["old"], x["new"]) >= t_start and min(x["old"], x["new"]) <= t_end]
        after = [x for x in steps or [] if x not in during and min(x["old"], x["new"]) > t_end]
    vals["pc_clock_steps"] = "; ".join(f"{x['delta_ms']:+.0f} ms at {hms(x['new'])}" for x in during + after)

    # the evening sync into analysis/session_checks (never the real sync folder), then D E F from its outputs
    sync, sync_err, T, why = None, None, None, "no Saleae CSV"
    if csv_p:
        sdict = {"day": f"{day_dir.name}/{s['cage']}/{s['mouse']}", "name": name, "python": py["path"],
                 "saleae_csv": csv_p, "sal": sal_p, "video": vid_p, "notes": [], "found": []}
        try:
            sync = ss.process_session(sdict, sync_dir, ctx.ffmpeg, {})
        except Exception as exc:
            sync_err, why = str(exc), f"sync failed: {exc}"
    if sync:
        text = (sync_dir / f"{name}_report.txt").read_text(encoding="utf-8")
        m = re.search(r"t = (-?[\d.]+) \+ (-?[\d.]+) \* elapsed", text)
        r = re.search(r"residual max ([\d.]+) ms", text)
        residual = f" (residual max {r.group(1)} ms)" if r else ""
        if m:
            sal = ss.read_saleae_csv(csv_p)
            cam, cap_end = sal["channels"]["camera"]["rises"], sal["end"]
            evdf = read_events(sync_dir / f"{name}_events.csv")
            T = session_timeline(cam, csv_first_time(csv_p), cap_end, (float(m.group(1)), float(m.group(2))),
                                 py_dur, bool(end), s.get("sal_lag") if sal_p else None, during)
            clock_side(T, during, t_start, evdf, [e["elapsed"] for e in ev if e.get("tag") == "reward_cmd"])
            vals["clock_end_s"] = round(T.d_end, 3) if T.d_end is not None else None
            vals["clock_start_s"] = round(T.d_start, 3) + 0.0 if T.d_start is not None else None   # no '-0.0'
            vals.update(check_valve(rep, evdf, n_cmd, sync, row, T, residual))                  # D
            vals.update(check_licks(rep, evdf, log, sync, is_pr, T, replay))                     # E
            check_clock(rep, T, py_dur, bool(end), during, steps is not None, t_start)          # F
            vals.update(check_camera(rep, s, ctx, sync, cam, cap_end, T, bool(end), t_start))
        else:
            why = "the sync report has no clock line"
    if T is None:
        for sec in ("D", "E"):
            rep.add("INFO", sec, "not checked", why)
        if during:
            rep.add("CHECK", "F", "clock", "; ".join(step_text(x, t_start) for x in during) + " - the task's times "
                    "after it are shifted, and the evening sync's matching after it is off", ko="PC 시계 변경")
        rep.add("INFO", "F", "Saleae part", f"not checked ({why}): clock, pulses vs frames, first pulse, pulses "
                f"after the end")
        if vid_p:
            vals.update(check_frames_estimate(rep, ctx, vid_p, py_dur, during))
    if vid_p:
        vals.update(check_video_file(rep, s, ctx, vid_p, t_end, after))                          # F (files)
    if sal_p:                                                                                    # G
        vals.update(check_saleae(rep, sal_p, row))
    else:
        rep.add("INFO", "G", "not checked", "no .sal file")
    check_sync(rep, sync, sync_err, replay, bool(T and (T.clock_off or during)))               # H

    overall = "CHECK" if rep.items else "PASS"
    day_txt = f"day{day_n if day_n is not None else '?'}" + (" (provisional)" if s.get("provisional") else "")
    head = [f"=== SESSION CHECK  {day_txt}  {s['cage']}/{s['mouse']}  {name}  "
            f"(Python start {start:%Y-%m-%d %H:%M:%S}, {py['dur'] / 60:.1f} min) ===",
            "protocol row: " + (f"days {row['from_day']}-{row['to_day']} {row.get('phase', '')} "
                                f"(approved: {row.get('approved') or '-'})" if row else
                                f"none for day {day_n if day_n is not None else '?'} - values without an "
                                f"expectation are INFO")]
    tail = f"OVERALL: {overall}" + (f" ({len(rep.items)}): " + "; ".join(rep.items) if rep.items else "")
    ko = "모두 정상" if not rep.items else "확인 필요: " + ", ".join(rep.items_ko)
    day_dir.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(head + rep.lines + [tail, f"{short_id(s)}: {ko}"]) + "\n", encoding="utf-8")
    print("\n".join(head + rep.lines + [tail]))
    print(f"  -> {report_path}")
    row_out = {c: "" for c in CSV_COLS}
    row_out.update({k: v for k, v in vals.items() if k in row_out and v is not None})
    row_out.update(checked_at=datetime.now().isoformat(timespec="seconds"),
                   day=day_n if day_n is not None else "", session=name,
                   cage=s["cage"], mouse=s["mouse"], start=f"{start:%Y-%m-%d %H:%M:%S}", overall=overall,
                   n_check=len(rep.items), check_items="; ".join(rep.items),
                   sync_status=sync["status"] if sync else ("FAIL" if sync_err else ""),
                   python=py["path"], sal=sal_p or "", saleae_csv=csv_p or "", video=vid_p or "", report=report_path)
    return overall, rep.items_ko, row_out


def read_csv_rows(path: Path):
    """Rows of session_checks.csv ([] when missing or unreadable)."""
    try:
        with open(path, newline="", encoding="utf-8-sig") as f:
            return list(csv.DictReader(f))
    except (OSError, csv.Error):
        return []


def update_csv(path: Path, new_rows):
    """Add the rows to session_checks.csv; the same session checked again (same day, name and start) replaces its
    old row. A log name used again on another date gets its own row."""
    try:
        rows = []
        if path.exists():
            with open(path, newline="", encoding="utf-8-sig") as f:
                rows = list(csv.DictReader(f))
        key = lambda r: (str(r.get("day")), r.get("session"), str(r.get("start")))      # noqa: E731
        keys = {key(r) for r in new_rows}
        rows = [r for r in rows if key(r) not in keys] + new_rows
        rows.sort(key=lambda r: (int(r.get("day") or 0), str(r.get("start") or "")))
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=CSV_COLS, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        return True
    except OSError as exc:                 # PermissionError included: open in Excel or another program
        print(f"!! could not update {path} ({exc.__class__.__name__}) - close it (Excel?) and run again; "
              f"the report files above were written")
        return False


# ----------------------------------------------------------------------------- test runs
def test_run_day(ctx, name, this_date):
    """Day whose expectations apply to a test run: (number or None, where it comes from). 'day<N>' anywhere in the
    test name, else the day of today's sessions / filed logs, else the day after the last filed day."""
    m = re.search(r"day(\d+)", name, re.I)
    if m:
        return int(m.group(1)), f"'day{m.group(1)}' in the test name"
    return guess_day(ctx, this_date)


def camera_trains(p: Path):
    """[(first, last)] camera pulse trains (gap > 1 s) of a Saleae CSV, read line by line; [] if it has no camera
    column or cannot be read."""
    rises = []
    try:
        with open(p, encoding="utf-8-sig", errors="replace") as f:
            hdr = [h.strip().lower() for h in f.readline().split(",")]
            if not hdr or not hdr[0].startswith("time") or "camera" not in hdr:
                return []
            ci, prev = hdr.index("camera"), "0"
            for ln in f:
                parts = ln.split(",")
                if len(parts) > ci:
                    v = parts[ci].strip()
                    if v == "1" and prev != "1":
                        rises.append(float(parts[0]))
                    prev = v
    except (OSError, ValueError):
        return []
    cam = np.array(rises)
    return [(float(cam[i]), float(cam[j])) for i, j in trains(cam)] if len(cam) else []


def pair_test_run(t, ctx):
    """.sal, Saleae CSV and video of a test run: (1) a .sal whose capture started up to 3 h before it, with a CSV
    holding a camera pulse train that starts where that .sal says the run started (+-1.5 s; also a test recorded
    into a session's capture); else (2) the .sal started 0-10 min before it and the CSV next to that .sal (no camera
    train: placed by the .sal start); else (3) no .sal: a CSV saved within an hour after the run with exactly one
    camera train as long as the run. Video: created at its start, not paired with a session.
    Returns (.sal path, run start in the capture [s], CSV path, its camera train (first, last), video path)."""
    t0, dur = t["py"]["start"], t["py"]["dur"]

    def saved_after_run(p):
        # a .sal or CSV export written before this run ended cannot hold it (e.g. an earlier test's CSV when the
        # re-test after a firmware upload was not exported) - never borrow it (added 2026-11-16)
        try:
            return p.stat().st_mtime >= t0 + dur - 1.0
        except OSError:
            return False

    sals = sorted(((t0 - x["start"], str(x["path"]), x) for x in (ds.read_sal(p) for p in ctx.found["sal"]) if x
                   and -5 <= t0 - x["start"] <= 3 * 3600 and saved_after_run(x["path"])), key=lambda z: z[:2])
    csvs = sorted(p for p in ctx.found["csv"] if saved_after_run(p) and ds.read_csv(p))
    cache = {}

    def tr_of(p):
        if p not in cache:
            cache[p] = camera_trains(p)
        return cache[p]

    pick = None
    for lag, _, sal in sals:                                   # (1) nearest capture start first
        for p in sorted(csvs, key=lambda q: q.parent != sal["path"].parent):
            hit = [(a, b) for a, b in tr_of(p) if abs(a - lag) <= 1.5]
            if hit:
                pick = (sal["path"], lag, p, hit[0])
                break
        if pick:
            break
    if pick is None and sals and sals[0][0] <= 600:            # (2)
        lag, _, sal = sals[0]
        same = [p for p in csvs if p.parent == sal["path"].parent]
        pick = (sal["path"], lag, same[0] if same else None, None)
    if pick is None and not sals:                              # (3)
        for p in sorted(csvs, key=lambda q: q.stat().st_mtime):
            if 0 <= p.stat().st_mtime - (t0 + dur) <= 3600:
                hit = [(a, b) for a, b in tr_of(p) if abs((b - a) - dur) <= RUN_MATCH_S]
                if len(hit) == 1:
                    pick = (None, None, p, hit[0])
                    break
    taken = paired_paths(ctx.all)
    vids = []
    for p in ctx.found["video"]:
        try:
            if str(p).lower() in taken or p.stat().st_size <= 1_000_000:
                continue
        except OSError:
            continue
        v = ds.read_video(p, ctx.ffmpeg)
        if v.get("created") is not None and -1.5 <= v["created"] - t0 <= LATE_VIDEO_S:
            vids.append((abs(v["created"] - t0), str(p), p))
    return (*(pick or (None, None, None, None)), min(vids)[2] if vids else None)


def check_test_run(t, ctx, files):
    """Reduced check of a test run (e.g. the morning valve test before the first mouse): valve open time vs the
    expected valve of the day, Python reward_cmd vs Saleae openings, Saleae settings, video frames vs camera
    pulses, camera trigger stopped. Writes test_runs/<name>_<date>_<time>_check.txt. Returns (overall, Korean)."""
    rep = Report()
    py = t["py"]
    name, t0 = py["name"], py["start"]
    start = datetime.fromtimestamp(t0)
    day_n, day_src = test_run_day(ctx, name, day_of(t0))
    row = protocol_for(ctx.protocol, day_n)
    log = ss.read_python_log(py["path"])
    ev, meta = log["events"], log["meta"] or {}
    end = next((e for e in reversed(ev) if e.get("tag") == "experiment_end"), None)
    dur = end["elapsed"] if end else (ev[-1]["elapsed"] if ev else 0.0)
    n_cmd = sum(1 for e in ev if e.get("tag") == "reward_cmd")
    sal_p, lag, csv_p, tr_ab, vid_p = files

    rep.add("PASS", "A", "python log", str(py["path"]))
    if sal_p:
        rep.add("PASS", "A", ".sal", str(sal_p))
    else:
        rep.add("CHECK", "A", ".sal", "NOT FOUND - save the Logic capture (the Saleae settings are read from it)",
                ko=".sal 없음")
    if csv_p:
        rep.add("PASS", "A", "Saleae CSV", str(csv_p))
    else:
        rep.add("CHECK", "A", "Saleae CSV", "NOT FOUND - export it from Logic (CSV, time in seconds): the valve open "
                "time is read from it", ko="Saleae CSV 없음(밸브 시간 확인 불가)")
    if vid_p:
        rep.add("PASS", "A", "video", str(vid_p))
    else:
        rep.add("INFO", "A", "video", "none created at the test start (optional for a test run)")
    hdr = ", ".join(f"{k} {meta[k]}" for k in ("reward_duration", "total_lockout_duration", "increment") if k in meta)
    rep.add("INFO", "C", "run", f"{meta.get('paradigm') or '-'}{f' ({hdr})' if hdr else ''}; {n_cmd} reward_cmd in "
            f"{dur:.1f} s" + ("" if end else "; no experiment_end (killed?)"))

    want, tol = num(row, "valve_ms"), num(row, "valve_tol_ms") or 0.0
    valve = None                           # (median ms, verdict) for the Korean summary
    if csv_p:
        sald = ss.read_saleae_csv(csv_p)
        ch = sald["channels"]
        cam = ch["camera"]["rises"] if "camera" in ch else np.array([])
        # where the test run sits in the capture: its camera train ('S'), else the .sal's Python start
        train = next(((i, j) for i, j in trains(cam) if abs(cam[i] - tr_ab[0]) < 1e-6), None) \
            if (tr_ab is not None and len(cam)) else None
        if train:
            w0, pad, where = float(cam[train[0]]), 0.0, f"its camera pulse train ({cam[train[0]]:.3f} s into the capture)"
        elif lag is not None:
            w0, pad, where = lag, 0.5, f"the .sal's Python start ({lag:.2f} s into the capture; no camera train there)"
        else:
            w0, pad, where = None, 0.0, "the whole capture (no .sal to place the run)"
        win = (w0 - PAD_BEFORE - pad, w0 + dur + PAD_AFTER + pad) if w0 is not None else (-np.inf, np.inf)
        if "reward" not in ch:
            rep.add("CHECK", "D", "valve time", f"the Saleae CSV has no 'reward' column (columns: "
                    f"{', '.join(sald['names'])})", ko="Saleae reward 채널 없음")
        else:
            iv = ch["reward"]["intervals"]
            on = np.array([a for a, _ in iv], dtype=float)
            d_ms = np.array([(b - a) * 1000 for a, b in iv], dtype=float)          # NaN: no falling edge
            inside = (on >= win[0]) & (on <= win[1])
            short = d_ms < qa.GLITCH_MS
            opens, glitch = inside & ~short, inside & short
            d = d_ms[opens & ~np.isnan(d_ms)]
            rep.note(f"the test run is placed by {where}")
            if len(d):
                med, lo, hi = float(np.median(d)), float(d.min()), float(d.max())
                text = f"open time median {med:.1f} ms (min {lo:.1f}, max {hi:.1f}), {len(d)} openings"
                if want is None:
                    rep.add("INFO", "D", "valve time", f"{text}; no expectation for day {day_n if day_n else '?'}")
                    valve = (med, "info")
                else:
                    outside = int(((d < want - tol) | (d > want + tol)).sum())
                    if outside == 0:
                        rep.add("PASS", "D", "valve time", f"{text} vs {want:g} +- {tol:g} ms (day {day_n})")
                        valve = (med, "ok")
                    elif abs(med - want) > tol:
                        rep.add("CHECK", "D", "valve time", f"{text} vs {want:g} +- {tol:g} ms (day {day_n}) - the "
                                f"board does not run the {want:g} ms valve: upload the {want:g} ms firmware before the "
                                f"first mouse", ko=f"밸브 {med:.0f} ms - {want:.0f} ms 펌웨어를 올려야 함")
                        valve = (med, "firmware")
                    else:
                        rep.add("CHECK", "D", "valve time", f"{text} vs {want:g} +- {tol:g} ms (day {day_n}): "
                                f"{outside} opening(s) outside", ko=f"밸브 {outside}회 {want:.0f}±{tol:.0f} ms 밖")
                        valve = (med, "some")
            elif n_cmd:
                rep.add("CHECK", "D", "valve time", "no valve opening on the Saleae reward channel during the test run",
                        ko="밸브 펄스 없음")
            else:
                rep.add("CHECK", "D", "valve time", "no reward in the test run - touch the spout until a few rewards "
                        "are given, then check again", ko="보상 없음(밸브 시간 확인 불가)")
            n_open = int(opens.sum())
            g_txt = f" ({int(glitch.sum())} pulse(s) < {qa.GLITCH_MS:g} ms ignored as glitches)" if glitch.any() else ""
            if n_open == n_cmd:
                rep.add("PASS", "D", "reward count", f"Python reward_cmd {n_cmd} = Saleae openings {n_open}{g_txt}")
            else:
                rep.add("CHECK", "D", "reward count", f"Python reward_cmd {n_cmd} vs Saleae openings {n_open}{g_txt}",
                        ko=f"보상 수 불일치(Python {n_cmd}, Saleae {n_open})")
            rest = ~inside & ~short
            if rest.any():
                dr = d_ms[rest & ~np.isnan(d_ms)]
                rep.add("INFO", "D", "other pulses", f"{int(rest.sum())} valve opening(s) outside the test run in this "
                        f"capture (another run, or flushing)" + (f", median {float(np.median(dr)):.1f} ms" if len(dr)
                                                                 else "") + " - not counted")
        # camera: frames vs pulses, trigger stopped at the end
        if train:
            i, j = train
            n_pulses = int(((cam[i:j + 1] >= win[0]) & (cam[i:j + 1] <= w0 + dur + PAD_AFTER)).sum())
            own_after = int((cam[i:j + 1] > w0 + dur + PAD_AFTER).sum())
            if vid_p:
                frames = ss.probe_video(vid_p, ctx.ffmpeg).get("frames")
                if frames is None:
                    rep.add("CHECK", "F", "frames", f"video frames could not be counted ({vid_p})",
                            ko="비디오 프레임 수 모름")
                elif frames == n_pulses:
                    rep.add("PASS", "F", "frames", f"video {frames} frames = {n_pulses} camera pulses of the test run")
                else:
                    diff = frames - n_pulses
                    rep.add("CHECK", "F", "frames", f"video {frames} frames vs {n_pulses} camera pulses of the test run "
                            f"-> {abs(diff)} {'EXTRA' if diff > 0 else 'SHORT'}",
                            ko=f"비디오 프레임 {abs(diff)}개 {'초과' if diff > 0 else '부족'}")
            if own_after:
                running = sald["end"] - float(cam[j]) <= 3.0 / CAM_HZ
                rep.add("CHECK", "F", "after end", f"{own_after} camera pulses after the test run ended"
                        + (" - the trigger was still running when the capture stopped: reset the Arduino before the "
                           "first mouse" if running else " - 'T' came late"), ko=f"종료 후 트리거 {own_after}펄스")
            else:
                rep.add("PASS", "F", "after end", "the camera trigger stopped at the end of the test run")
        else:
            rep.add("INFO", "F", "camera", "no camera pulse train at the test run's start (camera trigger not used?)"
                    + ("; video frames not compared" if vid_p else ""))
    else:
        rep.add("INFO", "D", "not checked", "no Saleae CSV")
    if sal_p:
        check_saleae(rep, sal_p, row)
    else:
        rep.add("INFO", "G", "not checked", "no .sal file")

    overall = "CHECK" if rep.items else "PASS"
    head = [f"=== TEST RUN CHECK  {name}  (Python start {start:%Y-%m-%d %H:%M:%S}, {dur:.0f} s) ===",
            f"expectations of day {day_n if day_n is not None else '?'} (assumed: {day_src}): "
            + (f"days {row['from_day']}-{row['to_day']} {row.get('phase', '')}, valve "
               + (f"{want:g} +- {tol:g} ms" if want is not None else "no expectation") if row else
               "no protocol row - values are INFO")]
    tail = f"OVERALL: {overall}" + (f" ({len(rep.items)}): " + "; ".join(rep.items) if rep.items else "")
    if valve and valve[1] == "ok":
        ko = f"시험 실행: 밸브 {valve[0]:.0f} ms 정상"
    elif valve and valve[1] == "info":
        ko = f"시험 실행: 밸브 {valve[0]:.0f} ms (day{day_n if day_n is not None else '?'} 기대값 없음)"
    elif valve and valve[1] == "firmware":
        ko = f"시험 실행: 밸브 {valve[0]:.0f} ms - {want:.0f} ms 펌웨어를 올려야 함"
    else:
        ko = "시험 실행"
    others = [k for k in rep.items_ko if not (valve and valve[1] == "firmware" and "펌웨어" in k)]
    if others:
        ko += (", " if ko != "시험 실행" else ": ") + "확인 필요: " + ", ".join(others)
    out_dir = ctx.project / "analysis" / "session_checks" / "test_runs"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / f"{name}_{start:%Y-%m-%d_%H%M%S}_check.txt"
    report_path.write_text("\n".join(head + rep.lines + [tail, ko]) + "\n", encoding="utf-8")
    print("\n".join(head + rep.lines + [tail]))
    print(f"  -> {report_path}")
    return overall, ko


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true", help="every session of the day")
    ap.add_argument("--session", help="one session by its Python log name (without .jsonl)")
    ap.add_argument("--date", help="YYYY-MM-DD (default: today)")
    ap.add_argument("--project", help="another project folder (tests): only this folder is searched")
    args = ap.parse_args(argv)
    try:
        sys.stdout.reconfigure(errors="replace")
    except Exception:
        pass
    try:
        day = date.fromisoformat(args.date) if args.date else date.today()
    except ValueError:
        print(f"--date must look like 2026-11-16, not '{args.date}'")
        print("확인할 세션 없음")
        return 2

    project = Path(args.project).resolve() if args.project else PROJECT
    if args.project:                       # tests: search only the given folder
        ds.PROJECT, ds.SEARCH, ds.ALIAS_FILE = project, [(project, 6)], project / "session_aliases.csv"
    ds.read_python = read_python_closed_is_finished

    ctx = SimpleNamespace(project=project, ffmpeg=ds.ffmpeg_exe())
    ctx.protocol, proto_path = read_protocol(project)
    ctx.found = ds.scan(day)
    ctx.sessions = ds.pair_day(ctx.found, ctx.ffmpeg, [])
    ctx.filed = filed_logs(project)
    real = [s for s in ctx.sessions if s["status"] == "session"]
    ctx.provisional = provisional_sessions(ctx)
    ctx.all = sorted(real + ctx.provisional, key=start_of)
    ctx.runs = other_runs(project, ctx.found, day)
    others = other_logs(ctx)
    print(f"session_check  {day}  project {project}")
    print(f"found: {len(real)} session(s)"
          + (f", {len(ctx.provisional)} full run(s) whose name the evening run cannot read (checked with a "
             f"provisional name)" if ctx.provisional else "")
          + f", {len(others)} other Python log(s); expectations: {proto_path}")
    if not ctx.protocol:
        print("  (protocol_expected.csv missing - every expected value is INFO)")

    ctx.csv_rows = read_csv_rows(project / "analysis" / "session_checks" / "session_checks.csv")
    ko_extra = []                          # Korean notes about logs that were not checked
    tests = []                             # a test run to check (reduced check)
    if args.session:
        key = args.session.lower().removesuffix(".jsonl")
        hit = [s for s in ctx.all if s["py"]["name"].lower() == key]
        if not hit:
            o = [o for o in others if o["name"].lower() == key]
            if o and o[0]["test"] and not o[0]["busy"]:
                tests = o[:1]
            else:
                why = f"{o[0]['kind']} - not a session that can be checked" if o else \
                    "no Python log with that name. Logs of that day: " + \
                    ", ".join(sorted({s["py"]["name"] for s in ctx.sessions} | {o["name"] for o in others}))
                print(f"{args.session}: {why}")
                if o and o[0]["busy"]:
                    print("아직 기록 중인 파일이 있어 점검하지 않았습니다: " + o[0]["name"])
                    return 3
                print("확인할 세션 없음")
                return 2
        targets, listed = hit, []
    elif args.all:
        targets, listed = list(ctx.all), others
    else:
        open_logs = [p for p in ctx.found["jsonl"] if ds.in_use(p)]
        if open_logs:
            for p in open_logs:
                print(f"STILL BEING WRITTEN (a task is running, or the log is open in another program): {p}")
            print("Run this again after the session has ended (--all checks the finished sessions of today now).")
            print("아직 기록 중인 파일이 있어 점검하지 않았습니다: " + ", ".join(p.name for p in open_logs))
            return 3
        targets = ctx.all[-1:]
        t_ref = targets[0]["py"]["start"] if targets else None
        listed = [o for o in others if t_ref is None or o["start"] > t_ref]
        if listed and listed[-1]["test"] and not listed[-1]["busy"]:     # the newest log is a test run: check it too
            # an aborted run with a session name stays listed as CHECK as well (it is also a mistake to look at)
            tests, listed = listed[-1:], (listed if listed[-1]["check"] else listed[:-1])
        for o in listed:
            if not o["check"]:
                ko_extra.append(f"최신 로그 {o['name']}({o['ko']}) 점검 안 함")
    if listed:
        print("\nPython logs not checked" + (" (newer than the checked session)" if not args.all else "") + ":")
        for o in listed:
            print(f"  {'CHECK' if o['check'] else 'INFO':5s} {o['name']:30s} {hms(o['start'])}  {o['kind']}")
    if targets and not args.all and not args.session and listed:
        print(f"checking the newest session: {targets[0]['py']['name']}")
    if tests and not args.session:
        print(f"the newest log is a test run: {tests[0]['name']} - reduced check (valve, Saleae, frames)"
              + (f" after the newest session {targets[0]['py']['name']}" if targets else ""))

    rows, results, skipped = [], [], []
    for s in targets:
        busy = busy_files(s, ctx.all, ctx.found)
        if busy:
            print(f"\n{s['py']['name']}: STILL BEING WRITTEN (or open in another program) - not checked:")
            for p in busy:
                print(f"    {p}")
            skipped.append(s)
            if not args.all:
                print("Stop the recording / close the file (IC Capture, Logic, a video player), then run this again.")
                print("아직 기록 중인 파일이 있어 점검하지 않았습니다: " + ", ".join(p.name for p in busy))
                return 3
            continue
        print()
        try:
            overall, items_ko, row = check_session(s, ctx)
        except Exception as exc:
            traceback.print_exc()
            overall, items_ko, row = "CHECK", [f"점검 오류({exc})"], None
        if row:
            rows.append(row)
        results.append((s, overall, items_ko))
    csv_failed = False
    if rows:
        out_csv = project / "analysis" / "session_checks" / "session_checks.csv"
        if update_csv(out_csv, rows):
            print(f"\nsummary table: {out_csv}")
        else:
            csv_failed = True

    test_results, test_busy = [], []       # the test run: (log entry, PASS/CHECK, Korean summary)
    for o in tests:
        files = pair_test_run(o, ctx)
        own = [o["path"]] + [p for p in (files[0], files[2], files[4]) if p]
        busy = [p for p in own if ds.in_use(p)]
        for p in ctx.found["video"]:                      # IC Capture still recording (not readable yet)
            try:
                if p not in own and p.stat().st_mtime >= o["start"] - 60 and \
                        str(p).lower() not in paired_paths(ctx.all) and ds.in_use(p):
                    busy.append(p)
            except OSError:
                continue
        if busy:
            print(f"\n{o['name']}: STILL BEING WRITTEN (or open in another program) - not checked:")
            for p in busy:
                print(f"    {p}")
            print("Stop the recording / close the file (IC Capture, Logic, a video player), then run this again.")
            test_busy.append(o)
            continue
        print()
        try:
            test_results.append((o, *check_test_run(o, ctx, files)))
        except Exception as exc:
            traceback.print_exc()
            test_results.append((o, "CHECK", f"시험 실행 점검 오류({exc})"))

    bad = [(s, items) for s, overall, items in results if overall != "PASS"]
    flagged = [o for o in listed if o["check"]]
    busy_logs = [o for o in listed if o["busy"]] + test_busy
    test_bad = [k for _, overall, k in test_results if overall != "PASS"]
    if results or not (test_results or test_busy):
        print(f"\nchecked {len(results)} session(s): {len(results) - len(bad)} PASS, {len(bad)} CHECK"
              + (f", {len(skipped)} skipped (still being written)" if skipped else "")
              + (f"; {len(flagged)} other log(s) to check" if flagged else ""))
    for o, overall, _ in test_results:
        print(f"{'' if results else chr(10)}checked test run {o['name']}: {overall}"
              + (f"; {len(flagged)} other log(s) to check" if flagged and not results else ""))
    ko = []
    if bad or flagged:
        ko.append("확인 필요: " + "; ".join([f"{short_id(s)} - {', '.join(items)}" for s, items in bad]
                                         + [f"{o['ko']} {o['name']}" for o in flagged]))
    elif results:
        ko.append(f"모두 정상 ({len(results)}개 세션)")
    elif not (test_results or test_busy):
        ko.append("확인할 세션 없음")
    ko += [k for _, _, k in test_results]
    if skipped or busy_logs:
        ko.append("기록 중이라 건너뜀: " + ", ".join([short_id(s) for s in skipped] + [o["name"] for o in busy_logs]))
    if csv_failed:
        ko.append("요약표(session_checks.csv) 갱신 실패 - 파일을 닫고 다시 실행")
    print(" / ".join(ko + ko_extra))
    if bad or flagged or test_bad:
        return 1
    if skipped or busy_logs:
        return 3
    if csv_failed:
        return 4
    return 0 if (results or test_results) else 2


if __name__ == "__main__":
    raise SystemExit(main())
