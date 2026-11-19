#!/usr/bin/env python3
"""
pr_daily_checks.py - report-only check table per mouse for each PR day, run in the evening after quick_analysis.py,
so a valve, lock-out or stop-rule slip on a PR day is seen the same evening and the next day can still be corrected.
It only reads the evening outputs and writes its own reports; nothing else is moved, renamed or changed.

One row per mouse/session (times in s from Python start; Saleae times are put on the Python clock with the sync fit):
  valve     open time of every commanded valve opening on the Saleae reward line (analysis/rewards_long.csv, checked
            against events.csv): median, min, max vs the expected value; CHECK when an opening is outside the
            tolerance. Session medians of one day differing by more than 20 ms = one CHECK line for the day; only the
            sessions outside the main group (the group nearest the expected value, else the largest) are shaded.
  lock-out  median reward_cmd -> next trial_start (Python log), next to the log meta total_lockout_duration and the
            expected value. The task (paradigm), requirement start/step and pump runs (injector_run) are checked too.
  session   length (experiment_end) and how it stopped: a Timed task that ends more than 0.5 s before its auto-stop
            (the expected session length, else 910 s), or whose last trial has no trial_end (= Ctrl+C), is a CHECK
            (INFO for an untimed task, which can only be stopped by hand); last reward, final gap (session end - last
            reward); rewards Python vs Saleae valve openings >= 5 ms (CHECK when they differ); video frames vs camera
            pulses (INFO, with the pulses of the session and of the whole capture).
  capture   the camera pulses (frames.csv) run from Python start to Python end, so pulses that start more than 1 s
            late or stop more than 1 s early mean the Saleae capture did not cover the session: CHECK. Python licks
            and rewards outside the capture are not compared with Saleae; the last-60-s Saleae licks are then n/a.
  Py-only   Python licks with no Saleae lick edge inside the capture (events.csv): all, valid=True, and the rewarded
            trials that counted at least one valid Python-only lick (CHECK on PR days, INFO on FR days). Python licks
            outside the capture are shown apart ('+N out').
  sync      sync FAIL = CHECK, and the Saleae-based cells of that session get '?' and are not checked. A session that
            is missing from sync_summary.csv, or a sync_summary.csv that cannot be read, is a CHECK (status unknown).
  PR        rewards and RL (RequirementLast) for the whole session and for the first 910 s; time to reward
            (reward_cmd - trial_start of the same trial) vs requirement: Theil-Sen slope (s per extra required lick)
            and median s per required lick; previous cohort RL on the same day index (days 3-5), labelled
            'reference only (training and valve differ)'. FR days get every column except these, which say n/a.
Day checks: every mouse of <project>/cohort_manifest.csv (cage_folder, mouse) has a session that day (CHECK if not; a
  log of that mouse filed under another day is named when its start date fits); all sessions of the day started on
  the same date (meta 't'): a session dated unlike the others is a CHECK (a session of another day filed here, e.g.
  because the output name was not changed and the task added '_1'); several sessions of one mouse = INFO, and the
  rows, CHECK lines and figure call them '<mouse> #1', '#2' (by start time).

RL = the last 'req'/'requirement' value in the log = the requirement of the last trial started (the old Step 02
RequirementLast). quick_analysis.pr_scores' pr_breakpoint is the highest 0-based trial index + 1, which for start 1 +1
is the same number: rewards + 1 when the last trial was still running at the stop, rewards when the stop came during
the lock-out right after a reward. RL is compared with quick_analysis' day<N>_metrics.csv (CHECK if they differ).
First 910 s = Python elapsed < 910 s (the clock of the auto-stop); quick_analysis' 15:10 window starts at the first
camera pulse, a few ms later.
tools/protocol_expected.csv: an empty cell = no expectation (value reported as INFO, never CHECK); a row whose
'approved' cell is empty or starts with 'no' is shown but not checked; a day without a row is not checked against
expected values (the data checks above still run) and gets PR columns only when its log is a progressive-ratio task;
a missing or unreadable file is a CHECK. Several rows of one day (days 6-10: active PR and passive replay): each
session is checked against the row of its task (paradigm), else the first row. A passive (replay) session has no PR
columns; its reward count is compared with the log it replayed (replay_source, found in python/day*/; another start
date is a CHECK) - the full replay check (pair, reward times) is session_check.py's. The glitch-filter column is
checked right after each session by session_check.py, not here.

Usage:  python pr_daily_checks.py --day 3          (also --day 3 4 5, or --day 3-5)
  --project DIR         project folder (default: the folder above tools/)
  --out-root DIR        where the outputs go (default: <project>/analysis; for tests)
  --prev-cohort DIR     previous-cohort archive (default: <project>/../PrevCohort_Dec2025_0910_0911_6099_6100)
  --rebuild-reference   compute the previous-cohort reference again
Inputs (read only): sync/sync_summary.csv, sync/day<N>/<cage>/<mouse>/<session>_events.csv, _frames.csv,
  _report.txt, _inputs.json, the session's Python jsonl (python/day<N>/<cage>/<mouse>/<session>.jsonl, else the file
  named in _inputs.json while it still has the size and time the sync read), analysis/rewards_long.csv,
  analysis/days/day<N>/day<N>_metrics.csv, tools/protocol_expected.csv, cohort_manifest.csv, the start line of every
  python/day*/ log (dates); the previous-cohort archive only to build the reference cache. The evening CSV outputs are
  read directly (a file only open in Excel is still read; one locked by another program is named); a Python log that
  another program still has open (daily_summary.in_use: the task still writing it) is skipped and named.
Outputs (nothing is written for a day without any session):
  analysis/days/day<N>/pr_checks_day<N>.csv      one row per session: all numbers, its own CHECK and INFO items
                                                 (n_check = the row's CHECK items; day-level CHECK lines such as the
                                                 valve spread between mice or a missing mouse are in the .txt only)
  analysis/days/day<N>/pr_checks_day<N>.txt      CHECK lines first, then INFO lines, the tables ('!' = CHECK cell,
                                                 '?' = Saleae value of a sync FAIL session, not checked)
  analysis/days/day<N>/figs/<NN>_pr_checks.png   the tables (CHECK cells red, '?' cells yellow), requirement completed
                                                 vs time with the previous-cohort band (PR days), valve open time,
                                                 lock-out per reward
  analysis/reference/prev_cohort_pr_days3to5.csv previous-cohort reference: built on the first run (how: its # lines),
                                                 then reused
Exit code: 0 all fine, 1 CHECK items or a day without sessions, 2 an output could not be written (open in Excel?).
"""
from __future__ import annotations

import sys
sys.dont_write_bytecode = True        # never leave .pyc files in tools/ (the tools below are only imported)

import argparse  # noqa: E402
import io  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402
import re  # noqa: E402
from collections import Counter  # noqa: E402
from datetime import datetime  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
import daily_summary as DS    # noqa: E402  in_use(), CAM_HZ
import sync_sessions as SS    # noqa: E402  read_python_log(), parse_iso()
import quick_analysis as QA   # noqa: E402  load_session(), parse_name(), GLITCH_MS, COMMON_WINDOW_S, mcol()

WINDOW_S = QA.COMMON_WINDOW_S  # 910 s: the PR window and the auto-stop (s from Python start)
LAST_S = 60.0                  # Saleae licks counted in the last 60 s of the session
VALVE_SPREAD_MS = 20.0         # CHECK when the valve medians of one day differ by more than this
COVER_TOL_S = 1.0              # camera pulses starting later / stopping earlier than this: Saleae capture incomplete
EARLY_STOP_S = 0.5             # a Timed task ending more than this before its auto-stop was stopped early
DEFAULT_TOL = {"valve_ms": 20.0, "lockout_s": 0.05, "session_s": 10.0}
REF_DAYS = (3, 4, 5)
REF_NOTE = "reference only (training and valve differ)"
PREV_DIR_NAME = "PrevCohort_Dec2025_0910_0911_6099_6100"
REF_FILE = "prev_cohort_pr_days3to5.csv"
MANIFEST = "cohort_manifest.csv"
QA_LAST_FIG = 11               # quick_analysis numbers its day figures 00-11
CHECK_COLOR, UNREL_COLOR, HEAD_COLOR = "#f9d2cd", "#fbf0c4", "#e8e8e8"
NA = "n/a"
SALEAE_COLS = {"valve", "rewards", "capture", "pyonly", "pytrials", "last60"}   # '?' in a sync FAIL session
KOREAN = {"valve": "밸브 시간", "valve between mice": "밸브 시간 마우스 간 차이", "lock-out": "락아웃",
          "session length": "세션 길이", "stopped early": "자동 종료 전에 끝남(손으로 중지?)",
          "stopped by hand": "마지막 시행 trial_end 없음(손으로 중지?)", "experiment_end": "로그 끝(experiment_end) 없음",
          "rewards Python vs Saleae": "보상 수 Python/Saleae 불일치",
          "Python-only lick counted": "Saleae에 없는 유효 Python 릭이 들어간 보상 시행",
          "Saleae coverage": "Saleae 캡처가 세션 전체를 덮지 않음", "paradigm": "과제(paradigm)",
          "requirement": "요구 릭 수(시작/증가)", "pump": "펌프", "rewards_long.csv": "rewards_long.csv 재생성 필요",
          "quick_analysis PR": "quick_analysis PR 값 차이", "sync FAIL": "싱크 FAIL",
          "sync status unknown": "싱크 상태 알 수 없음(sync_summary.csv)",
          "not in sync_summary": "sync_summary.csv에 없는 세션", "no sync": "싱크 결과 없는 세션",
          "protocol_expected.csv": "protocol_expected.csv 읽기 실패", "missing mouse": "세션 없는 마우스",
          "log date": "다른 날짜에 시작한 세션", "log changed": "싱크 뒤 바뀐 Python 로그",
          "file in use": "사용 중인 파일", "no Python log": "Python 로그 없음", "error": "처리 오류",
          "no sessions": "세션 없음", "replay source": "재생 원본 로그 없음", "replay date": "재생 원본이 다른 날짜",
          "replay rewards": "재생 보상 수가 원본과 다름", "replay": "재생"}


# ----------------------------------------------------------------------------- small helpers
def finite(v) -> bool:
    try:
        return v is not None and not isinstance(v, (str, bool)) and math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def num(v):
    """number from a CSV cell / log field, or None (empty, text, NaN)."""
    if isinstance(v, bool):
        return None
    try:
        x = float(str(v).strip())
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def fnum(v, nd=1):
    return f"{float(v):.{nd}f}" if finite(v) else "-"


def fint(v):
    return f"{int(round(float(v)))}" if finite(v) else "-"


def fg(v):
    """compact number: 1, 0.2, 910."""
    return f"{float(v):g}" if finite(v) else "-"


def med(x):
    x = [float(v) for v in x if finite(v)]
    return float(np.median(x)) if x else math.nan


def read_table(p: Path, what: str, notes: list, **kw):
    """CSV -> DataFrame, or None with a clear message (missing, locked, unreadable). The evening outputs are read
    directly, so a file that is only open in Excel is still read."""
    if not p.is_file():
        notes.append(f"{what} not found: {p}")
        return None
    try:
        return pd.read_csv(p, **kw)
    except PermissionError:
        notes.append(f"{what} could not be read: locked by another program (close it and run again): {p}")
    except Exception as exc:
        notes.append(f"{what} could not be read ({' '.join(str(exc).split())}): {p}")
    return None


def locked(p: Path) -> bool:
    """True when another program holds the file with an exclusive lock (it cannot even be opened to read)."""
    try:
        with open(p, "rb"):
            return False
    except PermissionError:
        return True
    except OSError:
        return False


def being_written(p: Path) -> bool:
    """a raw Python log that another program still has open (the task still writing it): daily_summary.in_use().
    A read-only file is never 'being written' (in_use() asks for write access)."""
    try:
        return p.is_file() and os.access(p, os.W_OK) and DS.in_use(p)
    except Exception:
        return False


def log_head(p: Path):
    """start time (meta 't', with its UTC offset) of a Python log from its first records, or None."""
    try:
        with open(p, encoding="utf-8-sig", errors="replace") as f:
            for _, line in zip(range(5), f):
                if line.strip():
                    rec = json.loads(line)
                    if rec.get("tag") == "meta":
                        return SS.parse_iso(str(rec.get("t", "")))
    except Exception:
        pass
    return None


def day_of(start) -> str:
    """local start date (the log's own UTC offset) as YYYY-MM-DD, or ''."""
    return start.date().isoformat() if start is not None else ""


def theil_sen(x, y):
    """median of all pairwise slopes (pairs with the same x are skipped); NaN with fewer than 3 points."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3:
        return math.nan
    i, j = np.triu_indices(len(x), k=1)
    dx = x[j] - x[i]
    keep = dx != 0
    if keep.sum() < 2:
        return math.nan
    return float(np.median((y[j] - y[i])[keep] / dx[keep]))


def requirement_last(events, t_max=None):
    """old Step 02 RequirementLast: the last 'req' / 'requirement' value in the log (records before t_max only)."""
    last = math.nan
    for e in events:
        if t_max is not None and float(e.get("elapsed", math.inf)) >= t_max:
            continue
        for k in ("requirement", "req"):
            v = num(e.get(k))
            if v is not None:
                last = v
    return last


def short_paradigm(p) -> str:
    s = str(p or "")
    for a, b in (("Experiment", ""), ("Progressive", "P"), ("Fixed", "F"), ("Ratio", "R")):
        s = s.replace(a, b)
    return s or "-"


def parse_days(tokens):
    days = []
    for tok in tokens:
        for part in str(tok).replace(",", " ").split():
            m = re.fullmatch(r"(?:day)?(\d+)(?:-(\d+))?", part.strip().lower())
            if not m:
                raise SystemExit(f"--day: '{part}' is not a day number (use e.g. 3, '3 4 5' or 3-5)")
            a, b = int(m.group(1)), int(m.group(2) or m.group(1))
            days += list(range(min(a, b), max(a, b) + 1))
    return sorted(set(days))


# ----------------------------------------------------------------------------- expected values, roster, log dates
def load_expected(path: Path, day: int, notes: list, problems: list):
    """(the protocol_expected.csv rows that cover this day [list, empty when none], 'row' | 'no row' |
    'unreadable'); a missing or unreadable file goes to problems (= CHECK lines). Several rows of one day (days 6-10:
    active PR and passive replay): each session is checked against the row of its paradigm (pick_row)."""
    P = read_table(path, "tools/protocol_expected.csv", problems, dtype=str, encoding="utf-8-sig")
    if P is None:
        return [], "unreadable"
    try:
        P = P.fillna("")
        P.columns = [str(c).strip().lower() for c in P.columns]
        rows = []
        for _, r in P.iterrows():
            a, b = num(r.get("from_day")), num(r.get("to_day"))
            if a is not None and a <= day <= (b if b is not None else a):
                rows.append(r)
        if not rows:
            return [], "no row"
        Es = []
        for r in rows:
            E = {"from_day": int(num(r["from_day"])), "to_day": int(num(r.get("to_day")) or num(r["from_day"])),
                 "phase": str(r.get("phase", "")).strip(),
                 "paradigms": [p.strip() for p in str(r.get("paradigm", "")).split("|") if p.strip()],
                 "req_start": num(r.get("req_start")), "req_step": num(r.get("req_step")),
                 "pump": str(r.get("pump", "")).strip().lower() or None,
                 "approved": str(r.get("approved", "")).strip(), "notes": str(r.get("notes", "")).strip()}
            for key, tol_key in (("valve_ms", "valve_tol_ms"), ("lockout_s", "lockout_tol_s"),
                                 ("session_s", "session_tol_s")):
                v, t = num(r.get(key)), num(r.get(tol_key))
                if v is not None and t is None:
                    t = DEFAULT_TOL[key]
                    notes.append(f"protocol_expected.csv: {key} {v:g} has no tolerance - {t:g} used")
                E[key], E[tol_key] = v, t
            E["checked"] = bool(E["approved"]) and not E["approved"].lower().startswith("no")
            Es.append(E)
        if len(Es) > 1:
            notes.append(f"protocol_expected.csv: {len(Es)} rows cover day {day} ("
                         + ", ".join(E["phase"] or "?" for E in Es) + ") - each session is checked against the row "
                         "of its task (paradigm), else the first one")
    except Exception as exc:                 # e.g. a renamed column
        problems.append(f"tools/protocol_expected.csv could not be read ({' '.join(str(exc).split())}): {path}")
        return [], "unreadable"
    return Es, "row"


def pick_row(Es: list, paradigm: str):
    """the row whose paradigm list holds this paradigm, else the first row (its paradigm check then flags it)."""
    for E in Es:
        if paradigm and paradigm in E["paradigms"]:
            return E
    return Es[0] if Es else None


def expected_text(E, state="row") -> str:
    if E is None:
        why = ("tools/protocol_expected.csv could not be read" if state == "unreadable" else
               "no row in tools/protocol_expected.csv for this day")
        return (f"{why} - expected values not checked; the data checks (valve spread between mice, rewards Python vs "
                "Saleae, Saleae capture, Python-only licks, early stop, roster, start dates) still apply")
    parts = [f"paradigm {' or '.join(E['paradigms']) or '(none)'}",
             f"requirement start {fg(E['req_start'])} step +{fg(E['req_step'])}"
             if (E["req_start"] is not None or E["req_step"] is not None) else "requirement (none)"]
    for key, tol, unit, lab in (("valve_ms", "valve_tol_ms", "ms", "valve"), ("lockout_s", "lockout_tol_s", "s", "lock-out"),
                                ("session_s", "session_tol_s", "s", "session")):
        parts.append(f"{lab} {fg(E[key])}±{fg(E[tol])} {unit}" if E[key] is not None else f"{lab} (none)")
    parts.append(f"pump {E['pump'] or '(none)'}")
    head = f"days {E['from_day']}-{E['to_day']}, {E['phase'] or '?'}"
    appr = f"approved: {E['approved']}" if E["checked"] else "NOT APPROVED (shown, not checked)"
    return f"{head}: " + " | ".join(parts) + f"   [{appr}]"


def expected_all(D: dict, sep: str) -> str:
    """the expected values of the day: one text per protocol_expected.csv row (days 6-10: active and passive)."""
    Es = D.get("Es") or []
    if len(Es) <= 1:
        return expected_text(D["E"], D["state"])
    return sep.join(expected_text(E) for E in Es)


def is_pr_day(E) -> bool:
    if E is None:
        return False
    return E["phase"].upper().startswith("PR") or any("progressive" in p.lower() for p in E["paradigms"])


def read_roster(project: Path, notes: list):
    """sorted '<cage_folder>_<mouse>' keys of cohort_manifest.csv, or None (then the roster is not checked)."""
    msg = []
    T = read_table(project / MANIFEST, MANIFEST, msg, dtype=str, encoding="utf-8-sig")
    if T is not None:
        T.columns = [str(c).strip().lower() for c in T.columns]
        if {"cage_folder", "mouse"} <= set(T.columns):
            keys = {f"{c}_{m}".lower() for c, m in zip(T["cage_folder"].fillna("").str.strip(),
                                                      T["mouse"].fillna("").str.strip()) if c and m}
            if keys:
                return sorted(keys)
        msg.append(f"{MANIFEST} has no cage_folder / mouse rows: {project / MANIFEST}")
    notes.append(f"{msg[0]} - mouse roster not checked")
    return None


def index_logs(project: Path):
    """every Python log filed under python/day<N>/<cage>/<mouse>/ with its start date, and the most common start
    date of each day folder (None on a tie)."""
    idx = []
    for p in sorted((project / "python").glob("day*/*/*/*.jsonl")):
        m = re.fullmatch(r"day(\d+)", p.parent.parent.parent.name, re.I)
        if m:
            idx.append({"day": int(m.group(1)), "mouse": f"{p.parent.parent.name}_{p.parent.name}".lower(),
                        "path": p, "date": day_of(log_head(p))})
    dates = {}
    for d in sorted({x["day"] for x in idx}):
        cnt = Counter(x["date"] for x in idx if x["day"] == d and x["date"]).most_common()
        dates[d] = cnt[0][0] if cnt and (len(cnt) == 1 or cnt[0][1] > cnt[1][1]) else None
    return idx, dates


# ----------------------------------------------------------------------------- previous cohort reference
def build_reference(prev_root: Path, notes: list):
    """Previous-cohort PR on days 3-5: RL from the session summary, first-910-s values from the jsonl the old
    analysis used (the frame-level ALL_mice_longitudinal.csv is not needed: the jsonl carries the Python clock)."""
    base = prev_root / "01_ALL_mice_longitudinal"
    summ_p = base / "session_summary_from_ALL_mice_longitudinal.csv"
    map_p = base / "day_map_and_session_check.csv"
    jroot = prev_root / "02_addiction_for_model_drive" / "addiction_for_model" / "addiction_concate_Dec_2025_json_only"
    if not summ_p.is_file():
        notes.append(f"previous-cohort reference not available: {summ_p} not found (use --prev-cohort)")
        return None, []
    S = pd.read_csv(summ_p)
    S = S[S["day_index"].isin(REF_DAYS) & S["Session_Paradigm"].astype(str).str.contains("Progressive", case=False)
          & (pd.to_numeric(S["passive"], errors="coerce").fillna(0) == 0)]
    n_all = len(S)
    S = S.sort_values(["mouse_key", "day_index", "session_idx"]).drop_duplicates(["mouse_key", "day_index"], keep="first")
    M = pd.read_csv(map_p) if map_p.is_file() else None
    rows, n_agree, n_json = [], 0, 0
    for s in S.itertuples(index=False):
        cage, color = str(s.mouse_key).split("_", 1)
        jname = ""
        if M is not None:
            m = M[(M["mouse_key"] == s.mouse_key) & (M["day_index"] == s.day_index) & (M["session_idx"] == s.session_idx)]
            if len(m):
                jname = str(m["jsonl"].iloc[0])
        jdir = jroot / f"day{int(s.day_index)}" / cage / color
        jp = jdir / jname if jname else None
        if jp is None or not jp.is_file():
            cand = sorted(jdir.glob("*.jsonl")) if jdir.is_dir() else []
            jp = cand[0] if len(cand) == 1 else None
        row = {"mouse_key": s.mouse_key, "cage": cage, "mouse": color, "day": int(s.day_index),
               "session_idx": int(s.session_idx), "jsonl": jp.name if jp else jname,
               "session_min_camera": round(float(s.session_min_camera), 3),
               "rewards_summary": num(s.inj_on), "req_last_summary": num(s.req_last)}
        if jp is not None:
            ev = SS.read_python_log(jp)["events"]
            rew = [e for e in ev if e.get("tag") == "reward_cmd"]
            req_of = {int(e["trial"]): num(e.get("req")) for e in ev if e.get("tag") == "trial_start" and "trial" in e}
            end = next((e for e in reversed(ev) if e.get("tag") == "experiment_end"), ev[-1] if ev else {})
            row.update(log_end_s=round(float(end.get("elapsed", math.nan)), 3), rewards_jsonl=len(rew),
                       req_last_jsonl=requirement_last(ev),
                       rewards_910s=sum(float(e["elapsed"]) < WINDOW_S for e in rew),
                       req_last_910s=requirement_last(ev, WINDOW_S),
                       reward_times_s=";".join(f"{float(e['elapsed']):.3f}" for e in rew),
                       reward_reqs=";".join(fg(req_of.get(int(e.get("trial", -1)))) for e in rew))
            n_json += 1
            agree = row["rewards_jsonl"] == row["rewards_summary"] and row["req_last_jsonl"] == row["req_last_summary"]
            n_agree += agree
            row["note"] = "" if agree else (f"jsonl differs from the summary: rewards {row['rewards_jsonl']} vs "
                                            f"{fg(row['rewards_summary'])}, RL {fg(row['req_last_jsonl'])} vs "
                                            f"{fg(row['req_last_summary'])}")
        else:
            row["note"] = "jsonl not found - first-910-s values missing"
        rows.append(row)
    how = [
        f"{REF_FILE} - previous cohort (0910/0911/6099/6100, PR days 3-5 = 2026-01-01..03): {REF_NOTE}.",
        f"made {datetime.now():%Y-%m-%d %H:%M} by pr_daily_checks.py; delete this file or use --rebuild-reference to "
        f"compute it again.",
        f"rows: {summ_p} with day_index 3-5, Session_Paradigm 'Progressive Ratio', passive 0 ({n_all} rows; "
        f"{n_all - len(S)} duplicate rows of the same mouse/day dropped, the first session_idx kept) -> {len(S)} sessions.",
        "req_last_summary = its req_last (= RequirementLast, the old pipeline's PR); rewards_summary = its inj_on "
        "(Injector_TTL onsets).",
        f"jsonl: {jroot}\\day<N>\\<cage>\\<mouse>\\<the file named in {map_p.name}> (the log the old analysis used); "
        f"found for {n_json} of {len(S)} sessions.",
        "req_last_jsonl = last 'req'/'requirement' in the jsonl (Step 02 definition); req_last_910s / rewards_910s = the "
        "same over records with elapsed < 910 s (Python start = 0; the frame-level ALL_mice_longitudinal.csv was not "
        "needed); reward_times_s / reward_reqs = elapsed and requirement of every reward_cmd (band in the figure).",
        f"check: jsonl reward count and RequirementLast equal the summary in {n_agree} of {n_json} sessions.",
        "protocol then: start 1 +1, valve only (179.8 ms measured on day 6), lock-out 0.2 s, stopped by hand "
        "(session_min_camera / log_end_s).",
    ]
    return pd.DataFrame(rows), how


def load_reference(out_root: Path, prev_root: Path, rebuild: bool, notes: list):
    """(reference rows, cache path, 'cache' | 'built' | 'missing')."""
    cache = out_root / "reference" / REF_FILE
    if cache.is_file() and not rebuild:
        try:
            lines = [ln for ln in cache.read_text(encoding="utf-8-sig").splitlines() if not ln.startswith("#")]
            return pd.read_csv(io.StringIO("\n".join(lines))), cache, "cache"
        except Exception as exc:
            notes.append(f"reference cache could not be read ({exc}) - computed again")
    ref, how = build_reference(prev_root, notes)
    if ref is None:
        return None, cache, "missing"
    try:
        cache.parent.mkdir(parents=True, exist_ok=True)
        with open(cache, "w", newline="", encoding="utf-8") as f:
            f.write("".join(f"# {h}\n" for h in how))
            ref.to_csv(f, index=False)
    except OSError as exc:
        notes.append(f"reference cache could not be written ({exc}): {cache}")
    return ref, cache, "built"


def ref_stats(ref, day):
    """median / range of RequirementLast (whole session and first 910 s) of one previous-cohort day."""
    if ref is None or day not in REF_DAYS:
        return None
    d = ref[pd.to_numeric(ref["day"], errors="coerce") == day]
    if d.empty:
        return None
    full = pd.to_numeric(d["req_last_summary"], errors="coerce")
    if "req_last_jsonl" in d:
        full = full.fillna(pd.to_numeric(d["req_last_jsonl"], errors="coerce"))
    w = pd.to_numeric(d["req_last_910s"], errors="coerce") if "req_last_910s" in d else pd.Series(dtype=float)
    full, w = full.dropna(), w.dropna()
    out = {"day": day, "n": len(full), "med": full.median(), "lo": full.min(), "hi": full.max(),
           "n910": len(w), "med910": w.median() if len(w) else math.nan,
           "lo910": w.min() if len(w) else math.nan, "hi910": w.max() if len(w) else math.nan}
    return out if out["n"] else None


def ref_text(st) -> str:
    if st is None:
        return NA
    t = f"day {st['day']}: RL median {fg(st['med'])}, range {fg(st['lo'])}-{fg(st['hi'])} (n={st['n']})"
    if st["n910"]:
        t += f"; first 910 s: median {fg(st['med910'])}, range {fg(st['lo910'])}-{fg(st['hi910'])} (n={st['n910']})"
    return t


def ref_band(ref, day, grid):
    """min / 25 % / median / 75 % / max over the previous mice of the requirement completed by each time."""
    if ref is None or "reward_times_s" not in ref:
        return None
    d = ref[pd.to_numeric(ref["day"], errors="coerce") == day]
    curves = []
    for _, r in d.iterrows():
        ts = [num(x) for x in str(r["reward_times_s"]).split(";") if num(x) is not None]
        qs = [num(x) for x in str(r.get("reward_reqs", "")).split(";")]
        if len(qs) != len(ts) or any(q is None for q in qs):
            qs = list(range(1, len(ts) + 1))          # start 1 +1
        ts, qs = np.asarray(ts, float), np.asarray(qs, float)
        curves.append([qs[ts <= g].max() if np.any(ts <= g) else 0.0 for g in grid])
    if not curves:
        return None
    return np.percentile(np.asarray(curves), [0, 25, 50, 75, 100], axis=0), len(curves)


# ----------------------------------------------------------------------------- one session
def add(items, level, col, key, text):
    """one CHECK / INFO item (MARK = shade the cell only); col = the table column it belongs to."""
    items.append({"level": level, "col": col, "key": key, "text": text})


def clock_fit(ev: pd.DataFrame):
    """saleae = a + b * python_elapsed, from the Python-only rows of events.csv (sync_sessions writes them so)."""
    p = ev[(ev["source"] == "python") & (ev["event"] != "lick")][["python_elapsed_s", "t_saleae_s"]].dropna()
    if len(p) >= 2 and p["python_elapsed_s"].nunique() >= 2:
        b, a = np.polyfit(p["python_elapsed_s"].to_numpy(float), p["t_saleae_s"].to_numpy(float), 1)
        return float(a), float(b)
    if len(p):
        return float((p["t_saleae_s"] - p["python_elapsed_s"]).median()), 1.0
    return None


def inside_project(p: Path, project: Path) -> bool:
    a, b = os.path.normcase(str(p.resolve())), os.path.normcase(str(project))
    return a == b or a.startswith(b.rstrip("\\/") + os.sep)


def find_log(ev_path: Path, session: str, project: Path, day_folder: str, items: list):
    """the session's Python jsonl: python/<day folder>/<session>.jsonl first; else the file named in
    <session>_inputs.json, only when it still has the size and time the sync recorded. A log whose size or time
    differs from what the sync read is a CHECK (the sync outputs are older than the log); no log = None."""
    fp = None
    try:
        fp = json.loads(ev_path.with_name(f"{session}_inputs.json").read_text(encoding="utf-8")).get("python")
    except Exception:
        pass

    def as_synced(p: Path):
        if not fp or len(fp) < 3:
            return None                          # unknown: no fingerprint
        st = p.stat()
        return st.st_size == fp[1] and int(st.st_mtime) == int(fp[2])

    local = project / "python" / day_folder / f"{session}.jsonl"
    named = Path(fp[0]) if fp else None
    if local.is_file():
        p = local
    elif named is not None and named.is_file() and as_synced(named):
        p = named
        if named.stem != session:
            add(items, "INFO", "mouse", "log name", f"Python log {named.name} (named in {session}_inputs.json) has "
                "another name than the session")
    else:
        return None
    if as_synced(p) is False:
        st = p.stat()
        add(items, "CHECK", "mouse", "log changed", f"{p.name} is not the file the sync read ({st.st_size} bytes, "
            f"modified {datetime.fromtimestamp(st.st_mtime):%Y-%m-%d %H:%M:%S}; the sync read {fp[1]} bytes, "
            f"{datetime.fromtimestamp(int(fp[2])):%Y-%m-%d %H:%M:%S}) - run the evening sync again")
    if not inside_project(p, project):
        add(items, "INFO", "mouse", "log outside project", f"Python log read from outside the project: {p}")
    return p


def log_numbers(log_p: Path):
    """task, schedule, lock-out, rewards, session length and stop from the Python jsonl -> (row fields, raw lists)."""
    L = SS.read_python_log(log_p)
    ev, meta = L["events"], (L["meta"] or {})
    starts = [e for e in ev if e.get("tag") == "trial_start"]
    rew = [e for e in ev if e.get("tag") == "reward_cmd"]
    end = next((e for e in reversed(ev) if e.get("tag") == "experiment_end"), None)
    session_s = float(end["elapsed"]) if end else (float(ev[-1]["elapsed"]) if ev else math.nan)
    reqs = [num(e.get("req")) for e in starts if num(e.get("req")) is not None]
    steps = np.diff(reqs) if len(reqs) > 1 else np.array([])
    st = np.array([float(e["elapsed"]) for e in starts])
    lock_t, lock = [], []
    for e in rew:                                    # lock-out = reward_cmd -> next trial_start
        k = int(np.searchsorted(st, float(e["elapsed"]), side="right"))
        if k < len(st):
            lock_t.append(float(e["elapsed"]))
            lock.append(float(st[k]) - float(e["elapsed"]))
    last_trial, last_end = "-", "-"
    if starts:                                       # Ctrl+C inside a trial leaves no trial_end for it
        last = starts[-1]
        k = num(last.get("trial"))
        ended = any(e.get("tag") == "trial_end" and float(e["elapsed"]) >= float(last["elapsed"])
                    and (k is None or num(e.get("trial")) == k) for e in ev)
        last_trial, last_end = f"trial {fg(k)}, req {fg(num(last.get('req')))}", ("yes" if ended else "no")
    start = L["start_wall"]
    r = {"paradigm": str(meta.get("paradigm", "")), "replay_source": str(meta.get("replay_source") or ""),
         "log_start": start.isoformat(timespec="seconds") if start else "",
         "log_date": day_of(start), "req_first": reqs[0] if reqs else math.nan,
         "req_step": float(np.median(steps)) if len(steps) else math.nan,
         "injector_runs": sum(e.get("tag") == "injector_run" for e in ev),
         "lockout_n": len(lock), "lockout_median_s": round(med(lock), 4),
         "lockout_min_s": round(min(lock), 4) if lock else math.nan,
         "lockout_max_s": round(max(lock), 4) if lock else math.nan,
         "lockout_meta_s": num(meta.get("total_lockout_duration")),
         "rewards_python": len(rew), "licks_python": sum(e.get("tag") == "lick" for e in ev),
         "experiment_end": "yes" if end else "no", "session_s": round(session_s, 3),
         "last_trial": last_trial, "last_trial_end": last_end,
         "last_reward_s": round(float(rew[-1]["elapsed"]), 3) if rew else math.nan,
         "final_gap_s": round(session_s - float(rew[-1]["elapsed"]), 3) if rew else math.nan}
    raw = {"ev": ev, "starts": starts, "rew": rew, "steps": steps, "lock_t": np.array(lock_t), "lock": np.array(lock)}
    return r, raw


def short_problems(problems: list, session: str, limit: int = 320) -> str:
    """the sync problems joined with ' | ', cut at an item boundary (the full text is in the sync report)."""
    out = ""
    for k, x in enumerate(problems):
        nxt = f"{out} | {x}" if out else x
        if len(nxt) > limit and out:
            return f"{out} | ... ({len(problems) - k} more in {session}_report.txt)"
        out = nxt
    return out


def capture_pulses(ev_path: Path, session: str, cp):
    """camera pulses of the whole Saleae capture = the session's + 'N camera pulses outside this python session'
    from the sync report (NaN when unknown)."""
    if cp is None:
        return math.nan
    try:
        txt = ev_path.with_name(f"{session}_report.txt").read_text(encoding="utf-8", errors="replace")
    except OSError:
        return math.nan
    m = re.search(r"(\d+) camera pulses outside this python session", txt)
    return cp + (int(m.group(1)) if m else 0)


def saleae_numbers(ev_path: Path, C: dict, r: dict, raw: dict, items: list, plot: dict):
    """licks, Python-only licks, valve openings, Saleae coverage and sync status from the sync outputs
    (quick_analysis' loader)."""
    session, day = r["session"], C["day"]
    S = QA.load_session(ev_path, C["sync_root"], C["summary"])
    status = S["status"]
    reliable = status != "FAIL"               # Saleae-based checks of a sync FAIL session are not trusted
    r["sync_status"] = status
    r["saleae_reliable"] = "yes" if reliable else "no (sync FAIL)"
    fit = clock_fit(S["ev"])
    if fit is None:
        raise ValueError("events.csv has no Python rows - the Saleae times cannot be put on the Python clock")
    a, b = fit

    def to_py(t):
        return (pd.to_numeric(pd.Series(t), errors="coerce").to_numpy(float) - a) / b

    session_s = r["session_s"]
    # Saleae coverage: the camera pulses (frames.csv = this session's pulses) run from Python start to Python end
    cov_lo, cov_hi, gaps = -math.inf, math.inf, []
    r["camera_first_s"] = r["camera_last_s"] = r["camera_s"] = math.nan
    r["capture"] = "unknown"
    frames_p = ev_path.with_name(f"{session}_frames.csv")
    if frames_p.is_file():
        cam = to_py(pd.read_csv(frames_p, usecols=["t_saleae_s"])["t_saleae_s"])
        cam = cam[np.isfinite(cam)]
        if len(cam):
            period = float(np.median(np.diff(cam))) if len(cam) > 1 else 1 / DS.CAM_HZ
            r["camera_first_s"], r["camera_last_s"] = round(float(cam[0]), 3), round(float(cam[-1]), 3)
            r["camera_s"] = round(float(cam[-1] - cam[0]), 3)
            c0, c1 = 0.0, session_s
            if cam[0] > COVER_TOL_S:
                cov_lo, c0 = float(cam[0]) - period, float(cam[0])
                gaps.append(f"camera pulses start {cam[0]:.1f} s after Python start")
            if finite(session_s) and session_s - cam[-1] > COVER_TOL_S:
                cov_hi, c1 = float(cam[-1]) + period, float(cam[-1])
                gaps.append(f"camera pulses stop {session_s - cam[-1]:.1f} s before the session end (last pulse at "
                            f"{cam[-1]:.1f} s of {session_s:.1f} s)")
            r["capture"] = f"{c0:.1f}-{c1:.1f}" if gaps else "full"
    else:
        add(items, "INFO", "capture", "no frames.csv", f"{frames_p.name} not found - Saleae coverage not checked")

    def in_capture(t):
        return ~((t < cov_lo) | (t > cov_hi))            # NaN counts as inside

    lk_t = to_py(S["licks"]["t_saleae_s"])
    r["licks_saleae"] = len(lk_t)
    r["saleae_licks_last60s"] = int(np.sum(lk_t > session_s - LAST_S)) \
        if (finite(session_s) and cov_hi == math.inf) else math.nan
    po = S["py_only_licks"]
    keep = in_capture(pd.to_numeric(po["python_elapsed_s"], errors="coerce").to_numpy(float))
    po_in, po_out = po[keep], po[~keep]
    pv = po_in[pd.to_numeric(po_in["valid_or_hit"], errors="coerce") == 1]
    rewarded = {int(num(e.get("trial"))) for e in raw["rew"] if num(e.get("trial")) is not None}
    vtrials = set(pd.to_numeric(pv["trial"], errors="coerce").dropna().astype(int)) & rewarded
    r["py_only_licks"], r["py_only_valid"], r["trials_with_valid_py_only"] = len(po_in), len(pv), len(vtrials)
    r["py_licks_outside_capture"] = len(po_out)
    rew_t = np.array([float(e["elapsed"]) for e in raw["rew"]], float)
    r["rewards_python_in_capture"] = int(in_capture(rew_t).sum())

    # valve openings: rewards_long.csv (quick_analysis), cross-checked with the same rules on events.csv
    V, src = S["valves"], "events.csv"
    if C["rl"] is not None:
        R = C["rl"][(pd.to_numeric(C["rl"]["day"], errors="coerce") == day) & (C["rl"]["session"].astype(str) == session)]
        if R.empty:
            if status == "FAIL":
                why = "quick_analysis skips sync FAIL sessions"
            elif not (V["kind"] == "reward").any():
                why = "no valve openings in this session"
            else:
                why = "quick_analysis.py not run after the sync?"
            add(items, "INFO", "valve", "rewards_long.csv", f"rewards_long.csv has no rows for this session ({why}) - "
                "valve times taken from events.csv with the same rules")
        else:
            n1, n2 = int((R["kind"] == "reward").sum()), int((V["kind"] == "reward").sum())
            m1 = med(R.loc[R["kind"] == "reward", "duration_ms"])
            m2 = med(V.loc[V["kind"] == "reward", "duration_ms"])
            if n1 == n2 and len(R) == len(V) and (not finite(m1) or abs(m1 - m2) < 0.5):
                V, src = R, "rewards_long.csv"
            else:
                add(items, "CHECK", "valve", "rewards_long.csv", f"rewards_long.csv ({len(R)} openings, {n1} commanded, "
                    f"median {fnum(m1, 2)} ms) differs from events.csv ({len(V)}, {n2}, {fnum(m2, 2)} ms) - the table "
                    "uses events.csv; run quick_analysis.py again")
    dur = pd.to_numeric(V["duration_ms"], errors="coerce")
    cmd = dur[(V["kind"] == "reward").to_numpy()].dropna()
    raw["valve_cmd_ms"] = cmd.to_numpy(float)
    r["valve_source"] = src
    r["valve_n"] = len(cmd)
    r["valve_median_ms"] = round(float(cmd.median()), 2) if len(cmd) else math.nan
    r["valve_min_ms"] = round(float(cmd.min()), 2) if len(cmd) else math.nan
    r["valve_max_ms"] = round(float(cmd.max()), 2) if len(cmd) else math.nan
    r["valve_extra_openings"] = int((V["kind"] == "extra_opening").sum())
    r["valve_short_pulses"] = int((V["kind"] == "short_pulse").sum())
    r["valve_glitches"] = int((dur < QA.GLITCH_MS).sum())
    r["rewards_saleae"] = int((dur >= QA.GLITCH_MS).sum())
    ok = dur.notna().to_numpy()
    plot["valve_t"], plot["valve_ms"] = to_py(V["t_saleae_s"])[ok], dur.to_numpy(float)[ok]

    # checks on the data themselves (no expected value needed)
    if gaps:
        out_valid = int((pd.to_numeric(po_out["valid_or_hit"], errors="coerce") == 1).sum())
        add(items, "CHECK" if reliable else "INFO", "capture", "Saleae coverage",
            f"Saleae capture does not cover the whole session: {'; '.join(gaps)} - {len(po_out)} Python lick(s) "
            f"({out_valid} valid) and {len(rew_t) - r['rewards_python_in_capture']} reward(s) fall outside it and are "
            "not compared with Saleae" + ("; Saleae licks of the last 60 s n/a" if cov_hi < math.inf else ""))
    if reliable and r["rewards_saleae"] != r["rewards_python_in_capture"]:
        where = f" inside the Saleae capture ({r['rewards_python']} in the log)" if gaps else ""
        add(items, "CHECK", "rewards", "rewards Python vs Saleae", f"rewards Python {r['rewards_python_in_capture']}"
            f"{where} vs Saleae {r['rewards_saleae']} valve openings >= {QA.GLITCH_MS:g} ms (commanded {len(cmd)}, "
            f"extra {r['valve_extra_openings']}, short {r['valve_short_pulses']})")
    if r["valve_glitches"]:
        add(items, "INFO", "valve", "valve glitches", f"{r['valve_glitches']} valve pulse(s) < {QA.GLITCH_MS:g} ms "
            "(not counted)")
    sm = C["summary"]
    hit = sm.loc[(sm["day"].astype(str) == r["day_folder"]) & (sm["session"].astype(str) == session), "problems"]
    problems = [x.strip() for x in str(hit.iloc[0]).split(" | ") if x.strip()] if len(hit) else []
    if status == "FAIL":
        add(items, "CHECK", "sync", "sync FAIL", f"{short_problems(problems, session)} - the Saleae-based cells of "
            "this session (valve, rewards Saleae, capture, Py-only, last 60 s) are not reliable ('?') and not checked")
    elif status == "?" and C["summary_ok"]:
        add(items, "CHECK", "sync", "not in sync_summary", "not in sync/sync_summary.csv - stale sync outputs (log moved "
            "or renamed?) or an unfinished sync; run the evening sync again")
    elif status == "WARN":
        other = [x for x in problems if "video has" not in x]       # the video counts get their own line below
        if other:
            add(items, "INFO", "sync", "sync WARN", short_problems(other, session))
    vf, cp = num(S["video_frames"]), num(S["camera_pulses"])
    cap = capture_pulses(ev_path, session, cp)
    r["video_frames"], r["camera_pulses"], r["camera_pulses_capture"] = \
        (vf if vf is not None else math.nan), (cp if cp is not None else math.nan), cap
    if vf is not None and cp is not None and vf != cp:
        d = int(vf - cp)
        txt = f"video {int(vf)} frames vs {int(cp)} camera pulses in the session"
        if finite(cap) and cap != cp:
            txt += f" ({int(cap)} in the whole capture, {int(cap - cp)} outside the Python session)"
        txt += f": {abs(d)} {'more' if d > 0 else 'fewer'} than the session pulses"
        if finite(cap) and cap != cp:
            dc = int(vf - cap)
            txt += f" ({abs(dc)} {'more' if dc > 0 else 'fewer'} than the whole capture)"
        add(items, "INFO", "video", "video frames", f"{txt} (sync {status})")


def pr_numbers(C: dict, r: dict, raw: dict, items: list, plot: dict):
    """PR columns: rewards and RL (whole session / first 910 s), time to reward vs requirement, reference."""
    ev, starts, rew = raw["ev"], raw["starts"], raw["rew"]
    start_of = {int(num(e["trial"])): float(e["elapsed"]) for e in starts if num(e.get("trial")) is not None}
    req_of = {int(num(e["trial"])): num(e.get("req")) for e in starts if num(e.get("trial")) is not None}
    xs, ys, rt = [], [], []
    for e in rew:
        k = int(num(e.get("trial"))) if num(e.get("trial")) is not None else -1
        if k in start_of and req_of.get(k):
            xs.append(req_of[k])
            ys.append(float(e["elapsed"]) - start_of[k])
            rt.append(float(e["elapsed"]))
    r.update(pr_rewards_full=len(rew), pr_RL_full=requirement_last(ev),
             pr_rewards_910s=sum(float(e["elapsed"]) < WINDOW_S for e in rew),
             pr_RL_910s=requirement_last(ev, WINDOW_S), ttr_n=len(xs),
             ttr_slope_s_per_req=round(theil_sen(xs, ys), 4),
             ttr_median_s_per_req=round(med(np.asarray(ys) / np.asarray(xs)), 4) if xs else math.nan)
    plot["rew_t"], plot["rew_req"] = np.asarray(rt), np.asarray(xs)
    st = C["ref_stats"]
    r.update(prev_ref_day=st["day"] if st else NA, prev_RL_median=st["med"] if st else NA,
             prev_RL_min=st["lo"] if st else NA, prev_RL_max=st["hi"] if st else NA, prev_RL_n=st["n"] if st else NA,
             prev_RL910_median=st["med910"] if st else NA, prev_RL910_min=st["lo910"] if st else NA,
             prev_RL910_max=st["hi910"] if st else NA, prev_ref_note=REF_NOTE if st else NA)
    q = C["qa"]
    edge = [float(e["elapsed"]) for e in starts if WINDOW_S <= float(e["elapsed"]) < WINDOW_S + 0.1]
    for win, col, tag in (("full_session", "pr_RL_full", "full"), ("first_15m10s", "pr_RL_910s", "910s")):
        qv = None
        if q is not None and "pr_breakpoint" in q:
            qq = q[(q["session"].astype(str) == r["session"]) & (q["window"].astype(str) == win)]
            qv = num(qq["pr_breakpoint"].iloc[0]) if len(qq) else None
        r["qa_pr_breakpoint_" + tag] = qv if qv is not None else math.nan
        if qv is not None and finite(r[col]) and qv != r[col]:
            txt = f"quick_analysis pr_breakpoint ({win}) {fg(qv)} vs RL {fg(r[col])} here"
            if tag == "910s" and edge:     # a trial started just after 910 s Python time: the two clocks differ
                add(items, "INFO", "pr_910", "quick_analysis PR", f"{txt}: a trial started at {edge[0]:.3f} s, at "
                    "the edge of the window (quick_analysis counts from the first camera pulse, a few ms later)")
            else:
                add(items, "CHECK", "pr_full" if tag == "full" else "pr_910", "quick_analysis PR", txt)
    if finite(r["pr_RL_full"]) and r["req_first"] == 1 and r["req_step"] == 1 and r["pr_RL_full"] != len(rew) + 1:
        why = "the session stopped during the lock-out right after a reward" if r["pr_RL_full"] == len(rew) else \
            "check the log"
        add(items, "INFO", "pr_full", "RL vs rewards", f"RL {fg(r['pr_RL_full'])} is not rewards + 1 "
            f"({len(rew)} rewards): {why}")


def expected_checks(C: dict, r: dict, raw: dict, items: list):
    """compare with the protocol_expected.csv row (CHECK only when the row is approved, else INFO), plus the stop
    rule and the Python-only licks."""
    E, lvl = C["E"], ("CHECK" if C["strict"] else "INFO")
    reliable = r.get("saleae_reliable") == "yes"
    cmd = raw.get("valve_cmd_ms", np.array([]))
    if E is not None and E["valve_ms"] is not None and len(cmd):
        lo, hi = E["valve_ms"] - E["valve_tol_ms"], E["valve_ms"] + E["valve_tol_ms"]
        n_out = int(((cmd < lo) | (cmd > hi)).sum())
        r["valve_outside_n"] = n_out
        if n_out and reliable:
            add(items, lvl, "valve", "valve", f"valve median {fnum(r['valve_median_ms'])} ms (min "
                f"{fnum(r['valve_min_ms'])}, max {fnum(r['valve_max_ms'])}); {n_out} of {len(cmd)} openings outside "
                f"{fg(E['valve_ms'])}±{fg(E['valve_tol_ms'])} ms")
    if E is not None and E["lockout_s"] is not None and finite(r["lockout_median_s"]):
        if abs(r["lockout_median_s"] - E["lockout_s"]) > E["lockout_tol_s"] + 1e-9:
            add(items, lvl, "lockout", "lock-out", f"lock-out {r['lockout_median_s']:.3f} s (reward_cmd -> next "
                f"trial_start, median of {r['lockout_n']}) outside {fg(E['lockout_s'])}±{fg(E['lockout_tol_s'])} s")
    tol_meta = E["lockout_tol_s"] if (E is not None and E["lockout_tol_s"]) else DEFAULT_TOL["lockout_s"]
    if r["lockout_meta_s"] is not None and finite(r["lockout_median_s"]) and \
            abs(r["lockout_meta_s"] - r["lockout_median_s"]) > tol_meta:
        pump = "; the pump run (injector_run) is not in the meta" if r["injector_runs"] else ""
        add(items, "INFO", "lockout", "lock-out meta vs measured", f"log meta total_lockout_duration "
            f"{fg(r['lockout_meta_s'])} s differs from the measured lock-out (table) by more than {fg(tol_meta)} s{pump}")
    if r["experiment_end"] != "yes":
        add(items, "CHECK", "session", "experiment_end", f"no experiment_end in the log (crash, or still running?) - "
            f"last record at {fnum(r['session_s'])} s")
    # how the session stopped: a Timed task must stop by itself at its auto-stop (the expected length, else 910 s)
    timed = r["paradigm"].lower().startswith("timed")
    auto_s = E["session_s"] if (E is not None and E["session_s"] is not None) else WINDOW_S
    r["auto_stop_s"] = auto_s if timed else math.nan
    early = timed and finite(r["session_s"]) and r["session_s"] < auto_s - EARLY_STOP_S
    no_end = r["last_trial_end"] == "no"
    how = ([f"a Timed task ending before its {fg(auto_s)} s auto-stop"] if early else []) + \
        ([f"the last trial ({r['last_trial']}) has no trial_end"] if no_end else [])
    r["stop"] = ("Ctrl+C" if no_end else "early" if early else "auto" if timed else "untimed")
    length_bad = E is not None and E["session_s"] is not None and finite(r["session_s"]) and \
        abs(r["session_s"] - E["session_s"]) > E["session_tol_s"]
    if length_bad:
        add(items, lvl, "session", "session length", f"session {r['session_s']:.1f} s (Python experiment_end) outside "
            f"{fg(E['session_s'])}±{fg(E['session_tol_s'])} s" + (f"; {' and '.join(how)}" if how else ""))
    if timed and (early or no_end) and not (length_bad and lvl == "CHECK"):
        add(items, "CHECK", "session", "stopped early" if early else "stopped by hand",
            f"session {fnum(r['session_s'])} s: {' and '.join(how)} - stopped by hand (Ctrl+C) or crashed")
    elif no_end and not timed and not length_bad:
        add(items, "INFO", "session", "stopped by hand", f"untimed task stopped by hand (Ctrl+C): the last trial "
            f"({r['last_trial']}) has no trial_end")
    if E is not None and E["paradigms"] and r["paradigm"] not in E["paradigms"]:
        add(items, lvl, "task", "paradigm", f"paradigm {r['paradigm'] or '(none)'} not {' or '.join(E['paradigms'])}")
    if E is not None and E["req_start"] is not None and finite(r["req_first"]) and r["req_first"] != E["req_start"]:
        add(items, lvl, "req", "requirement", f"first requirement {fg(r['req_first'])}, expected {fg(E['req_start'])}")
    steps = raw["steps"]
    if E is not None and E["req_step"] is not None and len(steps):
        bad = int(np.sum(steps != E["req_step"]))
        if bad:
            add(items, lvl, "req", "requirement", f"requirement step +{fg(r['req_step'])} (median; {bad} of "
                f"{len(steps)} steps differ), expected +{fg(E['req_step'])}")
    if E is not None and E["pump"] == "no" and r["injector_runs"]:
        add(items, lvl, "pump", "pump", f"{r['injector_runs']} injector_run in the log, expected no pump (valve only)")
    if E is not None and E["pump"] == "yes" and r["rewards_python"] and not r["injector_runs"]:
        add(items, lvl, "pump", "pump", "no injector_run in the log, expected the pump")
    if r["trials_with_valid_py_only"] and reliable:
        add(items, "CHECK" if r.get("pr") else "INFO", "pytrials", "Python-only lick counted",
            f"{r['trials_with_valid_py_only']} rewarded trial(s) counted a valid Python lick that has no Saleae lick "
            f"edge inside the capture ({r['py_only_valid']} valid of {r['py_only_licks']} Python-only licks) - "
            "the Python count and the Saleae lick line disagree in these trials")


def base_row(ev_path: Path, C: dict) -> dict:
    session = ev_path.name[: -len("_events.csv")]
    day_folder = ev_path.parent.relative_to(C["sync_root"]).as_posix()
    mouse, _, _ = QA.parse_name(session, day_folder)
    parts = day_folder.split("/")
    return {"day": C["day"], "cage": parts[1] if len(parts) >= 3 else "", "mouse": mouse, "label": mouse,
            "session": session, "day_folder": day_folder}


def check_session(ev_path: Path, C: dict, r: dict):
    """numbers (into r) + CHECK/INFO items + plot data of one session. C = the day's context."""
    session, day_folder = r["session"], r["day_folder"]
    items, plot = [], {"mouse": r["mouse"]}
    log_p = find_log(ev_path, session, C["project"], day_folder, items)
    busy = [p.name for p in (ev_path, ev_path.with_name(f"{session}_frames.csv")) if p.is_file() and locked(p)]
    if log_p is not None and (being_written(log_p) or locked(log_p)):
        busy.append(log_p.name)
    if busy:
        add(items, "CHECK", "mouse", "file in use", f"{', '.join(busy)} still open in another program (being written "
            "or locked) - session skipped; close it and run again")
        return items, None
    if log_p is None:
        add(items, "CHECK", "mouse", "no Python log", f"no Python log for these sync outputs (python/{day_folder}/"
            f"{session}.jsonl not found, and no matching file named in {session}_inputs.json) - session skipped; if "
            "the log was moved away (test / aborted run), these sync outputs are stale")
        return items, None
    r["python_log"] = str(log_p)
    fields, raw = log_numbers(log_p)
    r.update(fields)
    plot["lock_t"], plot["lock_s"] = raw["lock_t"], raw["lock"]
    saleae_numbers(ev_path, C, r, raw, items, plot)
    Es = C.get("Es") or ([C["E"]] if C["E"] is not None else [])
    E = pick_row(Es, r["paradigm"])            # days 6-10: the active PR row or the passive replay row
    r["_row"] = Es.index(E) if E is not None else None
    C = dict(C, E=E, strict=bool(E and E["checked"]), pr_day=is_pr_day(E))
    if C["E"] is not None:
        r["pr"] = bool(C["pr_day"])
    else:                     # a day without a protocol row: PR columns when the log is a progressive-ratio task
        r["pr"] = "progressive" in r["paradigm"].lower() or bool(len(raw["steps"]) and np.any(raw["steps"] != 0))
    if "replay" in r["paradigm"].lower():
        r["pr"] = False       # a passive replay has no requirement: its trials are the partner's
        replay_check(r, log_p, C["project"], items)
    if r["pr"]:
        pr_numbers(C, r, raw, items, plot)
    else:
        for c in PR_COLS:
            r[c] = NA
    expected_checks(C, r, raw, items)
    r["done"] = True
    return items, plot


def replay_check(r: dict, log_p: Path, project: Path, items: list):
    """passive (replay) log: the replayed log (replay_source) - found in python/day*/, same start date - and its
    reward count vs this session's. The full check (pair, times) is session_check.py's."""
    src = str(r.get("replay_source") or "").strip()
    if not src:
        add(items, "CHECK", "task", "replay source", "replay log without replay_source in its meta - the replayed log "
            "is unknown")
        return
    stem = Path(src.replace("\\", "/")).stem
    me = str(log_p.resolve()).lower()
    cands = [p for p in [log_p.parent / f"{stem}.jsonl"] + sorted((project / "python").glob(f"day*/*/*/{stem}.jsonl"))
             if p.is_file() and str(p.resolve()).lower() != me]
    if not cands:
        add(items, "CHECK", "task", "replay source", f"replayed log '{src}' not found in python/day*/ - its rewards are "
            "not compared")
        return
    logs = [(SS.read_python_log(p), p) for p in dict.fromkeys(cands)]
    mine = SS.read_python_log(log_p)["start_wall"]
    older = [(L, p) for L, p in logs if L["start_wall"] and mine and L["start_wall"] < mine]
    L, p = max(older, key=lambda x: x[0]["start_wall"]) if older else logs[0]
    n_src = sum(e.get("tag") == "reward_cmd" for e in L["events"])
    r["replay_source_rewards"] = n_src
    start = L["start_wall"]
    txt = f"replay of {p.relative_to(project).as_posix() if p.is_relative_to(project) else p}" + \
        (f" (started {start:%Y-%m-%d %H:%M})" if start else "")
    if start and r.get("log_date") and day_of(start) != r["log_date"]:
        add(items, "CHECK", "task", "replay date", f"{txt}: the replayed log started on {day_of(start)}, this session "
            f"on {r['log_date']}")
    if n_src != r["rewards_python"]:
        add(items, "CHECK", "rewards", "replay rewards", f"{txt}: {r['rewards_python']} rewards replayed, the source "
            f"has {n_src}")
    else:
        add(items, "INFO", "rewards", "replay", f"{txt}: {n_src} rewards, as in the source")


PR_COLS = ["pr_rewards_full", "pr_RL_full", "pr_rewards_910s", "pr_RL_910s", "qa_pr_breakpoint_full",
           "qa_pr_breakpoint_910s", "ttr_n", "ttr_slope_s_per_req", "ttr_median_s_per_req", "prev_ref_day",
           "prev_RL_median", "prev_RL_min", "prev_RL_max", "prev_RL_n", "prev_RL910_median", "prev_RL910_min",
           "prev_RL910_max", "prev_ref_note"]
CSV_COLS = ["day", "cage", "mouse", "label", "session", "log_date", "log_start", "sync_status", "saleae_reliable",
            "paradigm", "req_first", "req_step", "injector_runs",
            "valve_source", "valve_n", "valve_median_ms", "valve_min_ms", "valve_max_ms", "valve_expected_ms",
            "valve_tol_ms", "valve_outside_n", "valve_vs_main_group_ms", "valve_extra_openings", "valve_short_pulses",
            "valve_glitches", "lockout_n", "lockout_median_s", "lockout_min_s", "lockout_max_s", "lockout_meta_s",
            "lockout_expected_s", "lockout_tol_s", "licks_saleae", "licks_python", "py_only_licks", "py_only_valid",
            "trials_with_valid_py_only", "py_licks_outside_capture", "rewards_python", "rewards_python_in_capture",
            "rewards_saleae", "experiment_end", "session_s", "session_expected_s", "session_tol_s", "auto_stop_s",
            "stop", "last_trial", "last_trial_end", "capture", "camera_first_s", "camera_last_s", "camera_s",
            "video_frames", "camera_pulses", "camera_pulses_capture", "last_reward_s", "final_gap_s",
            "saleae_licks_last60s"] + PR_COLS + ["n_check", "checks", "infos", "python_log",
                                                 # days 6-10 (2026-11-18): the passive replay's source and its rewards,
                                                 # so the PR table can show "passive rewards = the partner's"
                                                 "replay_source", "replay_source_rewards"]

# table cells: (key, header); the .txt and the figure show two tables
SHOW = [("mouse", "mouse"), ("task", "task"), ("req", "req\nstart+step"), ("pump", "pump\nruns"),
        ("valve", "valve ms\nmedian (min-max)"), ("lockout", "lock-out s\nmeasured (meta)"),
        ("session", "session s\n(stop)"), ("rewards", "rewards\nPy / Saleae"), ("sync", "sync"),
        ("capture", "Saleae capture\ns covered"), ("video", "video frames\n(vs pulses)"),
        ("last", "last reward s\n(final gap s)"), ("last60", "Saleae licks\nlast 60 s"),
        ("pyonly", "Py-only licks\nall (valid)"), ("pytrials", "trials done with\nvalid Py-only lick"),
        ("pr_full", "PR whole session\nrewards / RL"), ("pr_910", "PR first 910 s\nrewards / RL"),
        ("slope", "s per extra req.\nlick (Theil-Sen)"), ("s_per", "median s per\nreq. lick")]
TABLES = [("protocol and session", ["mouse", "task", "req", "pump", "valve", "lockout", "session", "rewards", "sync",
                                    "capture", "video"]),
          ("session end, Python-only licks, PR",
           ["mouse", "last", "last60", "pyonly", "pytrials", "pr_full", "pr_910", "slope", "s_per"])]


def show_cells(r: dict) -> dict:
    """formatted cells of one row for the .txt tables and the figure ('?' = Saleae value of a sync FAIL session)."""
    if not r.get("done"):                         # session skipped (file in use, no log, error)
        return {k: (r["label"] if k == "mouse" else "-") for k, _ in SHOW}
    rew_py = f"{r['rewards_python_in_capture']} of {r['rewards_python']}" \
        if r.get("rewards_python_in_capture", r["rewards_python"]) != r["rewards_python"] else f"{r['rewards_python']}"
    vf, cp = r.get("video_frames"), r.get("camera_pulses")
    video = "none" if not finite(vf) else (f"{int(vf)} (=)" if vf == cp else
                                           f"{int(vf)} ({int(vf - cp):+d})" if finite(cp) else f"{int(vf)}")
    out = r.get("py_licks_outside_capture") or 0
    c = {"mouse": r["label"], "task": short_paradigm(r["paradigm"]),
         "req": f"{fg(r['req_first'])}+{fg(r['req_step'])}", "pump": fint(r["injector_runs"]),
         "valve": f"{fnum(r['valve_median_ms'])} ({fnum(r['valve_min_ms'])}-{fnum(r['valve_max_ms'])})"
                  if r.get("valve_n") else "-",
         "lockout": f"{fnum(r['lockout_median_s'], 3)} ({fg(r['lockout_meta_s'])})",
         "session": f"{fnum(r['session_s'])} ({r.get('stop', '-')})", "rewards": f"{rew_py} / {r['rewards_saleae']}",
         "sync": str(r.get("sync_status", "-")), "capture": str(r.get("capture", "-")), "video": video,
         "last": f"{fnum(r['last_reward_s'])} ({fnum(r['final_gap_s'])})",
         "last60": fint(r["saleae_licks_last60s"]) if finite(r["saleae_licks_last60s"]) else NA,
         "pyonly": f"{r['py_only_licks']} ({r['py_only_valid']})" + (f" +{out} out" if out else ""),
         "pytrials": fint(r["trials_with_valid_py_only"])}
    if r.get("pr"):
        c.update(pr_full=f"{fint(r['pr_rewards_full'])} / {fg(r['pr_RL_full'])}",
                 pr_910=f"{fint(r['pr_rewards_910s'])} / {fg(r['pr_RL_910s'])}",
                 slope=fnum(r["ttr_slope_s_per_req"], 2), s_per=fnum(r["ttr_median_s_per_req"], 2))
    else:
        c.update(pr_full=NA, pr_910=NA, slope=NA, s_per=NA)
    if r.get("saleae_reliable", "yes") != "yes":
        for k in SALEAE_COLS:
            c[k] += " ?"
    return c


# ----------------------------------------------------------------------------- one day
def check_day(day: int, G: dict):
    """all sessions of one day -> rows, CHECK/INFO items, plot data. G = shared tables and paths."""
    project, sync_root = G["project"], G["project"] / "sync"
    notes, problems = [], []
    Es, state = load_expected(project / "tools" / "protocol_expected.csv", day, notes, problems)
    E = Es[0] if Es else None                  # day-level uses (valve between mice); sessions use pick_row()
    pr_day = any(is_pr_day(e) for e in Es)
    st_ref = ref_stats(G["ref"], day) if (pr_day or E is None) else None
    qa = None
    if pr_day or E is None:
        qa_p = project / "analysis" / "days" / f"day{day}" / f"day{day}_metrics.csv"
        qa = read_table(qa_p, f"quick_analysis day{day}_metrics.csv", notes)
        if qa is None:
            notes[-1] += " (quick_analysis.py not run, or it stopped with an error) - RL not compared with quick_analysis"
    C = {"day": day, "project": project, "sync_root": sync_root, "summary": G["summary"],
         "summary_ok": G["summary_ok"], "rl": G["rl"], "qa": qa, "E": E, "Es": Es, "strict": bool(E and E["checked"]),
         "pr_day": pr_day, "ref_stats": st_ref}
    day_dir = sync_root / f"day{day}"
    ev_paths = sorted(day_dir.rglob("*_events.csv")) if day_dir.is_dir() else []
    rows, items_by_row, plots, day_items = [], [], [], []
    for p in ev_paths:
        r = base_row(p, C)
        try:
            items, plot = check_session(p, C, r)
        except Exception as exc:          # one broken session must not stop the others
            items, plot = [], None
            add(items, "CHECK", "mouse", "error", f"could not be processed: {' '.join(str(exc).split())}")
        rows.append(r)
        items_by_row.append(items)
        plots.append(plot)

    # Python logs of this day without sync outputs (sync skipped or FAIL)
    synced = {r["session"] for r in rows}
    py_dir = project / "python" / f"day{day}"
    unsynced, n_logs = [], 0
    for log in sorted(py_dir.rglob("*.jsonl")) if py_dir.is_dir() else []:
        n_logs += 1
        if log.stem in synced:
            continue
        rel = log.relative_to(py_dir).parts
        unsynced.append({"mouse": f"{rel[0]}_{rel[1]}".lower() if len(rel) >= 3 else log.stem, "path": log,
                         "date": day_of(log_head(log))})
        sm = G["summary"]
        hit = sm[sm["session"].astype(str) == log.stem]
        why = f"sync {hit['status'].iloc[0]}: {str(hit['problems'].iloc[0])[:200]}" if len(hit) else \
            "not in sync_summary.csv (no Saleae CSV found, or the sync was not run)"
        add(day_items, "CHECK", "", "no sync", f"{log.relative_to(project)} has no sync outputs - not in the table ({why})")
    empty = not ev_paths and not n_logs
    if empty:
        add(day_items, "CHECK", "", "no sessions", f"no sync outputs in {day_dir} and no Python logs in {py_dir} "
                                                   "(run the evening sync and quick_analysis first) - nothing written")
    elif not G["summary_ok"]:
        add(day_items, "CHECK", "sync", "sync status unknown", f"{G['summary_note']} - the sync status of every session "
            "is unknown ('?'): sync FAIL sessions cannot be told apart; run the evening sync again")

    # several sessions of one mouse: rows #1, #2 by start time
    by_mouse = {}
    for r in rows:
        by_mouse.setdefault(r["mouse"], []).append(r)
    for m, rs in by_mouse.items():
        rs.sort(key=lambda x: (x.get("log_start") or "~", x["session"]))
        for k, x in enumerate(rs, 1):
            x["label"] = m if len(rs) == 1 else f"{m} #{k}"
        if len(rs) > 1:
            add(day_items, "INFO", "", "several sessions", f"{len(rs)} sessions of {m} on day {day}: " + ", ".join(
                f"#{k} {x['session']} ({fnum(x.get('session_s'))} s)" for k, x in enumerate(rs, 1)))

    # start dates: all sessions of a day should have started on the same date
    dated = [r.get("log_date") for r in rows if r.get("log_date")] + [x["date"] for x in unsynced if x["date"]]
    cnt = Counter(dated).most_common()
    day_date = cnt[0][0] if cnt and (len(cnt) == 1 or cnt[0][1] > cnt[1][1]) else None
    if len(cnt) > 1:
        def also(d):
            m = [str(k) for k, v in G["folder_dates"].items() if v == d and k != day]
            return f" (= the start date of day {', '.join(m)})" if m else ""
        why = "a session of another day filed under day {} (output name not changed, so the task added '_1'?)"
        if day_date is None:
            add(day_items, "CHECK", "", "log date", f"the sessions of day {day} started on {len(cnt)} different dates ("
                + ", ".join(f"{d}: {n}" for d, n in sorted(cnt)) + ") - sessions of another day filed here?")
        else:
            for r, items in zip(rows, items_by_row):
                d = r.get("log_date")
                if d and d != day_date:
                    add(items, "CHECK", "mouse", "log date", f"started {d}{also(d)}, the other sessions of day {day} on "
                        f"{day_date} - " + why.format(day))
            for x in unsynced:
                if x["date"] and x["date"] != day_date:
                    add(day_items, "CHECK", "", "log date", f"{x['path'].relative_to(project)} started {x['date']}"
                        f"{also(x['date'])}, the other sessions of day {day} on {day_date} - " + why.format(day))

    # roster: every mouse of cohort_manifest.csv has a session this day
    present = {r["mouse"] for r in rows} | {x["mouse"] for x in unsynced}
    if G["roster"] is not None and not empty:
        for m in G["roster"]:
            if m in present:
                continue
            cage, mouse = m.split("_", 1)
            txt = (f"no session of {m} on day {day} (no sync outputs in sync/day{day}/{cage}/{mouse}/, no Python log in "
                   f"python/day{day}/{cage}/{mouse}/)")
            # a log of this mouse that looks out of place in its own day folder (dated unlike that folder's other
            # logs) and fits this day's date: probably this day's session saved under another day's name
            hint = [x for x in G["log_index"] if x["mouse"] == m and x["day"] != day and x["date"]
                    and G["folder_dates"].get(x["day"]) not in (None, x["date"])
                    and (day_date is None or x["date"] == day_date)]
            if hint:
                txt += " - but " + "; ".join(f"{x['path'].relative_to(project).as_posix()} started {x['date']}"
                                             for x in hint) + " (filed under another day: output name not changed?)"
            add(day_items, "CHECK", "", "missing mouse", txt)
        extra = sorted(present - set(G["roster"]))
        if extra:
            add(day_items, "INFO", "", "not in roster", f"not in {MANIFEST}: {', '.join(extra)}")

    # valve medians of the mice of this day (sessions with reliable Saleae data)
    vals = sorted((r["valve_median_ms"], r["label"], k) for k, r in enumerate(rows)
                  if r.get("done") and finite(r.get("valve_median_ms")) and r.get("saleae_reliable") == "yes")
    if len(vals) >= 2 and vals[-1][0] - vals[0][0] > VALVE_SPREAD_MS:
        groups = [[vals[0]]]                 # medians within 20 ms of the previous one form a group
        for v in vals[1:]:
            if v[0] - groups[-1][-1][0] <= VALVE_SPREAD_MS:
                groups[-1].append(v)
            else:
                groups.append([v])
        centers = [med([v for v, _, _ in g]) for g in groups]
        sizes = [len(g) for g in groups]
        if E is not None and E["valve_ms"] is not None:
            main = int(np.argmin([abs(c - E["valve_ms"]) for c in centers]))
            why = f"the group nearest the expected {fg(E['valve_ms'])} ms"
        elif sizes.count(max(sizes)) == 1:
            main, why = int(np.argmax(sizes)), "the largest group"
        else:
            main, why = None, "no group is larger than the others"
        txt = "; ".join(f"~{c:.0f} ms: {', '.join(lab for _, lab, _ in g)}" for c, g in zip(centers, groups))
        shade = "shaded: all (" + why + ")" if main is None else f"shaded: the sessions outside ~{centers[main]:.0f} ms " \
                                                                f"({why})"
        add(day_items, "CHECK", "valve", "valve between mice", f"valve medians differ between mice by "
            f"{vals[-1][0] - vals[0][0]:.1f} ms ({txt}); {shade}")
        for gi, g in enumerate(groups):
            if gi == main:
                continue
            for _, _, k in g:
                rows[k]["valve_vs_main_group_ms"] = round(rows[k]["valve_median_ms"] - centers[main], 1) \
                    if main is not None else math.nan
                add(items_by_row[k], "MARK", "valve", "valve between mice", "")

    for n in problems:
        add(day_items, "CHECK", "", "protocol_expected.csv", f"{n} - nothing is checked against the expected values")
    for n in notes:
        add(day_items, "INFO", "", "input", n)
    for Ei in Es:
        empty_cells = [lab for key, lab in (("valve_ms", "valve"), ("lockout_s", "lock-out"),
                                            ("session_s", "session length"), ("req_start", "requirement start"),
                                            ("req_step", "requirement step")) if Ei[key] is None]
        if not Ei["paradigms"]:
            empty_cells.append("paradigm")
        if Ei["pump"] not in ("yes", "no"):
            empty_cells.append("pump")
        if empty_cells:
            which = f"{Ei['phase'] or '?'} row: " if len(Es) > 1 else ""
            add(day_items, "INFO", "", "no expectation", f"{which}no expected value for: {', '.join(empty_cells)} - "
                                                         "values reported only (protocol_expected.csv cells empty)")
    for r in rows:   # expected values into the CSV rows (the row each session was checked against)
        k = r.pop("_row", None)
        Er = Es[k] if (k is not None and k < len(Es)) else E
        for key, col in (("valve_ms", "valve_expected_ms"), ("valve_tol_ms", "valve_tol_ms"),
                         ("lockout_s", "lockout_expected_s"), ("lockout_tol_s", "lockout_tol_s"),
                         ("session_s", "session_expected_s"), ("session_tol_s", "session_tol_s")):
            r[col] = Er[key] if (Er is not None and Er[key] is not None) else math.nan
    for r, items in zip(rows, items_by_row):
        r["n_check"] = sum(i["level"] == "CHECK" for i in items)
        r["checks"] = "; ".join(f"{i['key']}: {i['text']}" for i in items if i["level"] == "CHECK")
        r["infos"] = "; ".join(f"{i['key']}: {i['text']}" for i in items if i["level"] == "INFO")
    order = sorted(range(len(rows)), key=lambda k: (rows[k].get("cage", ""), rows[k]["mouse"],
                                                    rows[k].get("log_start") or "~", rows[k]["session"]))
    return {"day": day, "E": E, "Es": Es, "state": state, "pr_day": pr_day, "ref_stats": st_ref, "day_date": day_date,
            "dates": cnt, "present": present, "rows": [rows[k] for k in order],
            "items": [items_by_row[k] for k in order], "plots": [plots[k] for k in order], "day_items": day_items,
            "empty": empty}


# ----------------------------------------------------------------------------- outputs
def who(r: dict) -> str:
    return f"{r['label']} ({r['session']})" if "#" in r["label"] else r["label"]


def item_lines(D: dict, level: str):
    """CHECK lines (one per item) or INFO lines (the same text for several sessions on one line) of one day."""
    lines = [f"  {level:5s} day {D['day']:<13} {i['key']}: {i['text']}" for i in D["day_items"] if i["level"] == level]
    if level == "CHECK":
        for r, items in zip(D["rows"], D["items"]):
            for i in items:
                if i["level"] == "CHECK":
                    lines.append(f"  CHECK {who(r):<17} {i['key']}: {i['text']}")
        return lines
    groups = {}
    for r, items in zip(D["rows"], D["items"]):
        for i in items:
            if i["level"] == "INFO":
                groups.setdefault((i["key"], i["text"]), []).append(r["label"])
    for (key, text), labels in groups.items():
        if len(labels) == 1:
            w = labels[0]
        elif len(labels) == len(D["rows"]) and len(labels) > 2:
            w = f"all {len(labels)} sessions"
        else:
            w = f"{len(labels)} sessions: {', '.join(labels)}"
        lines.append(f"  INFO  {key}: {text}   [{w}]")
    return lines


def korean_summary(results, unwritten=()) -> str:
    parts = [f"출력 파일 저장 실패 {len(unwritten)}건 ({', '.join(unwritten)}: 열려 있으면 닫고 다시 실행)"] if unwritten else []
    for D in results:
        cnt = {}
        for i in D["day_items"] + [i for items in D["items"] for i in items]:
            if i["level"] == "CHECK":
                cnt[i["key"]] = cnt.get(i["key"], 0) + 1
        if cnt:
            parts.append(f"day{D['day']} " + ", ".join(f"{KOREAN.get(k, k)} {n}건" for k, n in cnt.items()))
    if not parts:
        return "모두 정상 (" + ", ".join(f"day{D['day']}" for D in results) + ")"
    return "확인 필요: " + "; ".join(parts)


def shaded(D: dict):
    """per row: the table columns that have a CHECK item (or a MARK)."""
    return [{i["col"] for i in items if i["level"] in ("CHECK", "MARK")} for items in D["items"]]


def unreliable(D: dict):
    """per row: the Saleae-based columns of a sync FAIL session ('?')."""
    return [SALEAE_COLS if r.get("saleae_reliable", "yes") != "yes" else set() for r in D["rows"]]


def day_kind(D: dict) -> str:
    if D["pr_day"]:
        return "PR day"
    base = {"unreadable": "protocol_expected.csv unreadable", "no row": "no protocol row"}.get(D["state"])
    if any(r.get("pr") for r in D["rows"]):
        return f"{base or 'not a PR day'}: PR columns from the progressive-ratio logs"
    return f"{base or 'not a PR day'}: PR columns n/a"


def roster_text(D: dict, G: dict) -> str:
    if G["roster"] is None:
        return f"roster: not checked ({MANIFEST} not read)"
    have = len(set(G["roster"]) & D["present"])
    return f"roster ({MANIFEST}): {len(G['roster'])} mice, {have} with a session"


def write_txt(D: dict, path: Path, G: dict):
    rows, E = D["rows"], D["E"]
    checks, infos = item_lines(D, "CHECK"), item_lines(D, "INFO")
    dates = ", ".join(f"{d} ({n})" for d, n in D["dates"]) or "-"
    lines = [f"PR DAILY CHECKS - day {D['day']}   (made {datetime.now():%Y-%m-%d %H:%M} by pr_daily_checks.py; "
             "report only, nothing was changed)",
             f"project: {G['project']}",
             f"expected ({day_kind(D)}): {expected_all(D, chr(10) + '    ')}",
             f"sessions in the table: {len(rows)} (sync/day{D['day']}); start date(s): {dates}; {roster_text(D, G)}", "",
             f"CHECK ({len(checks)}):"] + (checks or ["  none"]) + ["", f"INFO ({len(infos)}):"] + (infos or ["  none"])
    cells, marks, unrel = [show_cells(r) for r in rows], shaded(D), unreliable(D)
    head = dict(SHOW)
    for title, keys in TABLES:
        body = [[c[k] + (" !" if k in m else "") for k in keys] for c, m in zip(cells, marks)]
        T = pd.DataFrame(body, columns=[head[k].replace("\n", " ") for k in keys])
        lines += ["", f"TABLE - {title}   ('!' = CHECK, '?' = Saleae value of a sync FAIL session; times in s from "
                      "Python start; RL = RequirementLast)"]
        with pd.option_context("display.width", 400, "display.max_columns", 50, "display.max_colwidth", 40):
            lines += (T.to_string(index=False).splitlines() if len(T) else ["  (no sessions)"])
    if D["pr_day"] or any(r.get("pr") for r in rows):
        lines += ["", f"Previous cohort (0910/0911/6099/6100) RequirementLast, {REF_NOTE}:"]
        ref = G.get("ref")
        if ref is None:
            lines.append("  not available (" + "; ".join(n for n in G["notes"] if "reference" in n) + ")")
        else:
            for d in REF_DAYS:
                lines.append("  " + ref_text(ref_stats(ref, d)) + ("   <- same day index" if d == D["day"] else ""))
            lines.append(f"  source: {G.get('ref_path')} ({G.get('ref_how')})")
    lines += ["", "Definitions: valve = Saleae reward-line pulse width of each commanded opening (quick_analysis kind "
              "'reward'); lock-out = reward_cmd -> next trial_start in the Python log (meta = log meta "
              "total_lockout_duration); session (stop): auto = Timed task stopped by itself, Ctrl+C = the last trial "
              "has no trial_end, early = Timed task ended before its auto-stop, untimed = not a Timed task; rewards "
              "Py = reward_cmd ('a of b' = a inside the Saleae capture); Saleae capture = the part of the session "
              "covered by the camera pulses ('full' = all); video = frames (frames - camera pulses of the session); "
              "Py-only = Python lick with no Saleae lick edge inside the capture ('+N out' = Python licks outside the "
              "capture, not checked); final gap = session end (experiment_end) - last reward_cmd; RL = last 'req' in "
              "the log (= quick_analysis pr_breakpoint for start 1 +1); first 910 s = Python elapsed < 910 s; "
              "Theil-Sen slope = time to reward (reward_cmd - trial_start) per extra required lick; a sync FAIL "
              "session's Saleae cells ('?') are not reliable and not checked.",
              "", "Summary: " + korean_summary([D])]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def fig_number(figdir: Path) -> Path:
    """reuse this tool's own figure number; otherwise the next number after quick_analysis' 00-11 and any other."""
    own = sorted(figdir.glob("*_pr_checks.png"))
    if own:
        return own[0]
    nums = [int(m.group(1)) for p in figdir.glob("*.png") if (m := re.match(r"^(\d+)_", p.name))]
    return figdir / f"{max(nums + [QA_LAST_FIG]) + 1:02d}_pr_checks.png"


def table_panel(ax, keys, cells, marks, unrel, title):
    """one table of the figure: CHECK cells red, '?' cells (sync FAIL) yellow, n/a grey."""
    ax.axis("off")
    head = dict(SHOW)
    labels = [head[k] for k in keys]
    body = [[c[k] for k in keys] for c in cells] or [["(no sessions)"] + [""] * (len(keys) - 1)]
    w = np.array([max([len(s) for s in labels[j].split("\n")] + [len(str(row[j])) for row in body]) + 2
                  for j in range(len(keys))], float)
    tb = ax.table(cellText=body, colLabels=labels, cellLoc="center", bbox=[0, 0, 1, 1], colWidths=list(w / w.sum()))
    tb.auto_set_font_size(False)
    tb.set_fontsize(8.5)
    for (i, j), cell in tb.get_celld().items():
        if i == 0:
            cell.set_facecolor(HEAD_COLOR)
            cell.set_text_props(weight="bold")
            cell.set_height(cell.get_height() * 2)
            continue
        if i - 1 >= len(marks):
            continue
        if keys[j] in marks[i - 1]:
            cell.set_facecolor(CHECK_COLOR)
        elif keys[j] in unrel[i - 1]:
            cell.set_facecolor(UNREL_COLOR)
        if cell.get_text().get_text() == NA:
            cell.get_text().set_color("0.55")
    ax.set_title(title, fontsize=9.5, loc="left")


def later(r: dict):
    """the session number (2, 3, ..) of a mouse's later session of the day, None for its only / first session."""
    m = re.search(r"#(\d+)$", r["label"])
    return int(m.group(1)) if m and m.group(1) != "1" else None


def scatter_panel(ax, D, key_t, key_y, expected, tol, ylabel, title):
    """one dot per valve opening / reward per session, with the expected band."""
    if expected is not None:
        ax.axhspan(expected - tol, expected + tol, color="g", alpha=0.12, label=f"expected {fg(expected)}±{fg(tol)}")
    for r, p in zip(D["rows"], D["plots"]):
        if p is not None and len(p.get(key_y, [])):
            ax.plot(np.asarray(p[key_t]) / 60, p[key_y], "." if later(r) is None else "x", ms=3,
                    color=QA.mcol(r["mouse"]), label=r["label"])
    ax.set_xlabel("time from Python start (min)")
    ax.set_ylabel(ylabel)
    ax.set_title(title + ("" if expected is not None else " (no expected value)"), fontsize=9, loc="left")
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3)


def write_fig(D: dict, path: Path, G: dict):
    rows, E = D["rows"], D["E"]
    cells, marks, unrel = [show_cells(r) for r in rows], shaded(D), unreliable(D)
    n = max(len(rows), 1)
    table_h, plot_h = 0.28 * (n + 2), 4.4            # inches: header (2 units) + one unit per row; plot row
    pr = D["pr_day"] or any(r.get("pr") for r in rows)
    widths = [1.5, 1, 1] if pr else [1, 1]
    fig = plt.figure(figsize=(22, 2 * table_h + plot_h + 2.2), layout="constrained")
    gs = fig.add_gridspec(3, len(widths), height_ratios=[table_h, table_h, plot_h], width_ratios=widths)
    n_check = len(item_lines(D, "CHECK"))
    for gi, (title, keys) in enumerate(TABLES):
        if gi == 0:
            title = (f"day {D['day']} - PR daily checks: {len(rows)} sessions, {n_check} CHECK item(s) (red cells; "
                     f"yellow '?' = Saleae value of a sync FAIL session; see pr_checks_day{D['day']}.txt)\n"
                     f"expected: {expected_all(D, chr(10) + '          ')}\n{title}")
        table_panel(fig.add_subplot(gs[gi, :]), keys, cells, marks, unrel, title)
    col = 0
    if pr:                      # requirement completed vs time, with the previous-cohort band
        a1 = fig.add_subplot(gs[2, 0])
        col = 1
        grid = np.arange(0, 1001, 5.0)
        band = ref_band(G.get("ref"), D["day"], grid) if D["day"] in REF_DAYS else None
        if band is not None:
            (lo, q1, md, q3, hi), nb = band
            a1.fill_between(grid / 60, lo, hi, color="0.6", alpha=0.18, lw=0,
                            label=f"previous cohort day {D['day']} range (n={nb})")
            a1.fill_between(grid / 60, q1, q3, color="0.45", alpha=0.25, lw=0, label="previous cohort 25-75 %")
            a1.plot(grid / 60, md, color="0.3", ls="--", lw=1.2, label=f"previous cohort median - {REF_NOTE}")
        for r, p in zip(rows, D["plots"]):
            if p is None or "rew_t" not in p:
                continue
            a1.step(np.r_[0, p["rew_t"] / 60], np.r_[0, p["rew_req"]], where="post", color=QA.mcol(r["mouse"]), lw=1.4,
                    ls="-" if later(r) is None else "--", label=f"{r['label']} (RL {fg(r.get('pr_RL_full'))})")
        a1.axvline(WINDOW_S / 60, color="r", ls="--", lw=0.8)
        a1.set_xlabel("time from Python start (min)   red dashed line = 910 s")
        a1.set_ylabel("requirement completed (licks)")
        st = D["ref_stats"]
        a1.set_title("requirement completed vs time" + (f"   (previous cohort {ref_text(st)})" if st else
                     "   (no previous-cohort reference for this day)"), fontsize=9, loc="left")
        if a1.get_legend_handles_labels()[0]:
            a1.legend(fontsize=7, loc="upper left")
        a1.grid(alpha=0.3)
    ev = E["valve_ms"] if E is not None else None
    scatter_panel(fig.add_subplot(gs[2, col]), D, "valve_t", "valve_ms", ev, E["valve_tol_ms"] if ev is not None else 0,
                  "valve open time (ms, Saleae)", "valve open time of every opening")
    el = E["lockout_s"] if E is not None else None
    scatter_panel(fig.add_subplot(gs[2, col + 1]), D, "lock_t", "lock_s", el, E["lockout_tol_s"] if el is not None else 0,
                  "lock-out (s)", "lock-out of every reward (reward_cmd -> next trial_start)")
    fig.savefig(path, dpi=110)
    plt.close(fig)


def write_csv(D: dict, path: Path, G: dict):
    T = pd.DataFrame([{k: v for k, v in r.items() if k != "done"} for r in D["rows"]])
    for c in CSV_COLS:
        if c not in T:
            T[c] = math.nan
    T[CSV_COLS].to_csv(path, index=False, encoding="utf-8-sig")


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--day", nargs="+", required=True, help="day number(s): 3, or 3 4 5, or 3-5")
    ap.add_argument("--project", type=Path, default=TOOLS.parent, help="project folder (default: the folder above tools/)")
    ap.add_argument("--out-root", type=Path, help="output root (default: <project>/analysis)")
    ap.add_argument("--prev-cohort", type=Path, help=f"previous-cohort archive (default: <project>/../{PREV_DIR_NAME})")
    ap.add_argument("--rebuild-reference", action="store_true", help="compute the previous-cohort reference again")
    a = ap.parse_args(argv)
    try:
        sys.stdout.reconfigure(errors="replace")
    except Exception:
        pass
    days = parse_days(a.day)
    project = a.project.resolve()
    G = {"project": project, "out_root": (a.out_root or project / "analysis").resolve(), "notes": []}
    print(f"PR daily checks - project {project}, day(s) {', '.join(map(str, days))}")
    if not project.is_dir():
        print(f"  project folder not found: {project}\n\n확인 필요: 프로젝트 폴더 없음")
        return 1
    prev_root = a.prev_cohort or project.parent / PREV_DIR_NAME
    G["ref"], G["ref_path"], G["ref_how"] = load_reference(G["out_root"], prev_root, a.rebuild_reference, G["notes"])
    print(f"  previous-cohort reference: {G['ref_how']} - {G['ref_path']}")
    snotes = []
    summ = read_table(project / "sync" / "sync_summary.csv", "sync/sync_summary.csv", snotes, dtype=str,
                      encoding="utf-8-sig")
    if summ is not None and not {"day", "session", "status"} <= set(summ.columns):
        snotes.append(f"sync/sync_summary.csv has no day / session / status columns: {project / 'sync'}")
        summ = None
    G["summary_ok"], G["summary_note"] = summ is not None, (snotes[-1] if snotes else "")
    if summ is None:
        summ = pd.DataFrame(columns=["day", "session", "status", "video_frames", "camera_pulses", "problems"])
    for c in ("video_frames", "camera_pulses", "problems"):
        if c not in summ:
            summ[c] = ""
    G["summary"] = summ.fillna("")
    G["rl"] = read_table(project / "analysis" / "rewards_long.csv", "analysis/rewards_long.csv", G["notes"])
    if G["rl"] is None:
        G["notes"].append("valve times taken from events.csv (same rules as quick_analysis)")
    G["roster"] = read_roster(project, G["notes"])
    G["log_index"], G["folder_dates"] = index_logs(project)
    rc, results, unwritten = 0, [], []
    for day in days:
        D = check_day(day, G)
        for n in G["notes"]:                    # input notes about shared files go into every day's report
            if all(n not in i["text"] for i in D["day_items"]):
                add(D["day_items"], "INFO", "", "input", n)
        results.append(D)
        checks = item_lines(D, "CHECK")
        print(f"\n=== day {day}: {len(D['rows'])} session(s), {day_kind(D)}, {len(checks)} CHECK, "
              f"{len(item_lines(D, 'INFO'))} INFO line(s)")
        for ln in checks:
            print(ln)
        if checks and rc == 0:
            rc = 1
        if D["empty"]:
            continue
        out_dir = G["out_root"] / "days" / f"day{day}"
        figdir = out_dir / "figs"
        figdir.mkdir(parents=True, exist_ok=True)
        csv_p, txt_p, fig_p = out_dir / f"pr_checks_day{day}.csv", out_dir / f"pr_checks_day{day}.txt", fig_number(figdir)
        for fn, p in ((write_csv, csv_p), (write_txt, txt_p), (write_fig, fig_p)):
            try:
                fn(D, p, G)
                print(f"  -> {p}")
            except PermissionError:
                print(f"  !! could not write {p} - close it (Excel / viewer?) and run again")
                unwritten.append(p.name)
                rc = 2
    print("\n" + korean_summary(results, unwritten))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
