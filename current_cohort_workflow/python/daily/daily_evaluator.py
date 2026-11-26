#!/usr/bin/env python3
"""
daily_evaluator.py - an INDEPENDENT check of the day's routine (added 2026-11-22; user: "I sometimes find you read
wrongly file so please alsway every day please make evaluator whether you read file correctly, and QC correctly and
analysiscorrectly").

The routine's own readers are NOT reused: the raw files are read again with this file's small readers (python log
lines, digital.csv edges, video frame count + creation time, .sal capture start, pupil frames) and compared with what
the routine wrote. Report only - nothing is changed. FAIL = stop and look before using the day; WARN = look;
INFO = for the record. The rules it compares against (which tree, which option lists) come from concat_choice.py.

A  FILES - was the right file read?
 A1 present    one python log, one .sal, one digital.csv and one video per mouse in the day folders (a session
               without video must be in behavioronly_sessions.csv)
 A2 names      every file's name shows its folder's mouse (cage + colour anywhere) and day N
 A3 log        the log's date = the day's date (one date per day, no other day has it); its paradigm = the expected
               one (tools/protocol_expected.csv, group from cohort_manifest.csv); a passive log replays the manifest
               partner's log of the same day; the log ends with experiment_end
 A4 video      the video file was created within 2 s of the python start (mp4 creation_time)
 A5 saleae     the capture started before the python start and holds the whole session; every python reward has its
               valve opening (offset = the most common time difference), no extra openings inside the session
 A6 copy       the project copy of the log has the same bytes as the original in Downloads\\behavior_task_windows
 A7 unique     no file content is used for two mouse-days (size + first + last MB)
 A8 assays     TI file value (B1) = tail_immersion_all.csv; one Straub / TST / HOT video per mouse and day in the raw
               folders (two = FAIL; straub_video/name_fixes.csv corrections applied); recordings of that date without
               a mouse name are listed
B  QC
 B1 sync       sync_summary.csv has the session; its status and problems
 B2 alignment  the routine's video frames / camera pulses = this count; unequal -> the session is in
               alignment_fixes.csv and its table is the sync-aligned copy
 B3 pupil      out % and median size recomputed from the pupil file the check used = pupil_check table; that file
               has one row per video frame
 B4 verdict    a verdict is written; out >= 10 % with 'ok' needs a decision note
 B5 table      Step 02 read the same pupil table as the check (session_inventory session_dir = table_run)
C  ANALYSIS
 C1 PR         last requirement and rewards from the log = pr_checks = PR table = Step 02 (RequirementLast, rewards)
 C2 passive    drops = the manifest partner's rewards of that day
 C3 rewards    Saleae rewards in pr_checks = python rewards
 C4 TI         assays file = TI table = Step 02 Immersion_Latency_s
 C5 straub     straub_dayN.csv auto % = recomputed from the per-frame CSV; one per mouse
 C6 step02     Step 02 rows = video frames (camera pulses without video); isPassive = the log's paradigm; a day in
               tools/step02_window_days.csv (day 13 = 20 min, user 2026-11-26): rows = camera pulses of the first N s
               and the window marker names the whole table (= all video frames)
Outputs: analysis\\days\\dayN\\evaluator_dayN.csv / .txt / .png. Exit code 1 if any FAIL.
usage: python daily_evaluator.py --day N        python daily_evaluator.py --days 1 2 3 ...
"""
from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import re
import subprocess
import sys
import textwrap
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

T = Path(__file__).resolve().parent
P = T.parent
sys.path.insert(0, str(T))
import concat_choice as CC  # noqa: E402

try:                                       # Windows idle priority: never slow the rig
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x40)
except Exception:
    pass

HOME = Path.home()
MICE = [(c, m) for c in ("9669m", "9670f") for m in ("black", "green1", "green2", "white")]
ORDER = ["A1 present", "A2 names", "A3 log", "A4 video", "A5 saleae", "A6 copy", "A7 unique", "A8 assays",
         "B1 sync", "B2 alignment", "B3 pupil", "B4 verdict", "B5 table",
         "C1 PR", "C2 passive", "C3 rewards", "C4 TI", "C5 straub", "C6 step02"]
RANK = {"FAIL": 3, "WARN": 2, "PASS": 1, "INFO": 0}
# the pupil check's rule (pupil_check_examples.py); checked against its source at start
JUMP_POS, JUMP_SIZE, WIN, ODD_HI, ODD_LO = 0.5, 0.30, 61, 1.5, 0.5
OFF_EYE_PX, OFF_WIN = 20, 1800
RULE_LINES = ("JUMP_POS, JUMP_SIZE, WIN = 0.5, 0.30, 61", "ODD_HI, ODD_LO = 1.5, 0.5", "OFF_EYE_PX, OFF_WIN = 20, 1800")
VALVE_MIN_S, VALVE_MAX_S, MATCH_TOL_S = 0.05, 5.0, 0.15
VIDEO_PASS_S, VIDEO_WARN_S = 2.0, 60.0
OUT_TOL_PP, MED_TOL = 0.3, 0.01
NOISY_PCT = 10.0
ASSAY_MIN_S = {"straub": 300, "tst": 180, "hot": 30}
UNNAMED = re.compile(r"^(\d{4}-\d{2}-\d{2})_\d{2}-\d{2}-\d{2}\.(mp4|avi|mov)$", re.I)


# ----------------------------------------------------------------------------------------------- small readers
def who(name: str):
    """(cage, mouse, day) shown in a file name: cage + colour anywhere ('9670fgreen' = green1), 'dayN' anywhere"""
    s = name.lower()
    m = re.search(r"(?<!\d)(\d{4})([mf])?.*?(green1|green2|green|white|black)", s)
    d = re.search(r"day\s*_?(\d+)", s)
    cages = {c[:4]: c for c, _ in MICE}                  # '9669' without the sex letter = 9669m (2026-11-22 names)
    cage = (m.group(1) + m.group(2) if m.group(2) else cages.get(m.group(1))) if m else None
    mouse = ("green1" if m.group(3) == "green" else m.group(3)) if m else None
    return cage, mouse, (int(d.group(1)) if d else None)


def read_log(p: Path) -> dict:
    recs, bad = [], 0
    for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
        ln = ln.strip()
        if not ln:
            continue
        try:
            recs.append(json.loads(ln))
        except ValueError:
            bad += 1
    start = next((r["t"] for r in recs if isinstance(r, dict) and "t" in r), None)
    meta = next((r for r in recs if r.get("tag") == "meta"), {})
    rew = [float(r.get("elapsed", np.nan)) for r in recs if r.get("tag") == "reward_cmd"]
    req = np.nan
    for r in recs:                                   # the last 'requirement' / 'req' in the file
        for k in ("requirement", "req"):
            if k in r:
                try:
                    req = float(r[k])
                except (TypeError, ValueError):
                    pass
    end = [r for r in recs if r.get("tag") == "experiment_end"]
    el = [float(r["elapsed"]) for r in recs if "elapsed" in r]
    try:
        st = datetime.fromisoformat(start) if start else None
    except ValueError:
        st = None
    return dict(start=st, paradigm=str(meta.get("paradigm", "")), replay=str(meta.get("replay_source") or ""),
                rewards=np.array(rew, float), last_req=req, ended=bool(end), bad=bad,
                session_s=float(end[-1]["elapsed"]) if end else (max(el) if el else np.nan))


def read_digital(p: Path) -> dict:
    D = pd.read_csv(p)
    t = D.iloc[:, 0].to_numpy(float)
    out = {"end": float(t[-1]) if len(t) else np.nan}
    for key, names in (("valve", ("reward", "valve", "injector")), ("camera", ("camera",))):
        col = next((c for c in D.columns[1:] if any(n in str(c).lower() for n in names)), None)
        if col is None:
            out[key] = None
            continue
        v = D[col].to_numpy(int)
        up = np.flatnonzero((v[1:] == 1) & (v[:-1] == 0)) + 1
        dn = np.flatnonzero((v[1:] == 0) & (v[:-1] == 1)) + 1
        rises, falls = t[up], t[dn]
        k = np.searchsorted(falls, rises)
        dur = np.full(len(rises), np.nan)
        ok = k < len(falls)
        dur[ok] = falls[k[ok]] - rises[ok]
        out[key] = (rises, dur)
    return out


def ffmpeg():
    c = sorted((T / ".venv" / "Lib" / "site-packages" / "imageio_ffmpeg" / "binaries").glob("ffmpeg*.exe"))
    return str(c[-1]) if c else None


def video_info(p: Path, ff) -> dict:
    import cv2
    cap = cv2.VideoCapture(str(p))
    n, fps = int(cap.get(cv2.CAP_PROP_FRAME_COUNT)), float(cap.get(cv2.CAP_PROP_FPS) or 0)
    cap.release()
    created = None
    if ff:
        r = subprocess.run([ff, "-hide_banner", "-i", str(p)], capture_output=True, text=True, errors="replace")
        m = re.search(r"creation_time\s*:\s*(\S+)", r.stderr)
        if m:
            try:
                created = datetime.fromisoformat(m.group(1).replace("Z", "+00:00"))
            except ValueError:
                created = None
    return dict(frames=n, fps=fps, created=created)


def sal_start(p: Path):
    try:
        with zipfile.ZipFile(p) as z:
            st = json.loads(z.read("meta.json"))["data"]["captureStartTime"]
        ms = float(st["unixTimeMilliseconds"]) + float(st.get("fractionalMilliseconds", 0))
        return datetime.fromtimestamp(ms / 1000.0, tz=timezone.utc)
    except Exception:
        return None


def quick_hash(p: Path, mb: int = 1) -> str:
    h, n, k = hashlib.sha1(), p.stat().st_size, mb << 20
    h.update(str(n).encode())
    with open(p, "rb") as fh:
        h.update(fh.read(k))
        if n > 2 * k:
            fh.seek(-k, 2)
            h.update(fh.read())
    return h.hexdigest()


def full_hash(p: Path) -> str:
    return hashlib.sha1(p.read_bytes()).hexdigest()


def match_valve(e: np.ndarray, r: np.ndarray, session_s: float) -> dict:
    """python reward times e (s from python start) vs valve openings r (s in the capture): offset = the most common
    difference r - e (20-ms bins, refined by the median near it); matched = rewards with an opening within 0.15 s"""
    if len(e) == 0 or len(r) == 0:
        return dict(offset=np.nan, matched=0, in_window=int(len(r)))
    dif = (r[None, :] - e[:, None]).ravel()
    lo, hi = float(dif.min()), float(dif.max())
    h, edges = np.histogram(dif, bins=max(1, int((hi - lo) / 0.02) + 1), range=(lo, hi + 0.02))
    k = int(np.argmax(h))
    off = 0.5 * (edges[k] + edges[k + 1])
    near = dif[np.abs(dif - off) < 0.1]
    off = float(np.median(near)) if len(near) else off
    matched = int(sum(float(np.min(np.abs(r - (x + off)))) < MATCH_TOL_S for x in e))
    inw = int(((r >= off - 0.5) & (r <= off + (session_s if np.isfinite(session_s) else 0) + 0.5)).sum())
    return dict(offset=off, matched=matched, in_window=inw)


def roll_med(x, w, mp):
    return pd.Series(x).rolling(w, center=True, min_periods=mp).median().to_numpy()


def package_out(clean: Path):
    C = pd.read_csv(clean, usecols=["Diameter_pixels", "SmoothedDiameter_px", "CenterX_px", "CenterY_px"])
    dr = pd.to_numeric(C.Diameter_pixels, errors="coerce").to_numpy(float)
    d = pd.to_numeric(C.SmoothedDiameter_px, errors="coerce").to_numpy(float)
    cx = pd.to_numeric(C.CenterX_px, errors="coerce").to_numpy(float)
    cy = pd.to_numeric(C.CenterY_px, errors="coerce").to_numpy(float)
    held = np.r_[False, dr[1:] == dr[:-1]] | ~np.isfinite(dr)
    med = float(np.nanmedian(d))
    mp = max(5, WIN // 6)
    mx, my, md = (roll_med(np.where(held, np.nan, v), WIN, mp) for v in (cx, cy, dr))
    with np.errstate(invalid="ignore", divide="ignore"):
        jump = ~held & ((np.hypot(cx - mx, cy - my) > JUMP_POS * med) | (np.abs(dr / md - 1) > JUMP_SIZE))
        odd = ~held & ~jump & ((d > ODD_HI * med) | (d < ODD_LO * med))
    return len(C), 100.0 * float((held | jump | odd).mean()), med


def eye_out(frames: Path, clean: Path):
    W = pd.read_csv(frames, usecols=["found", "cx", "cy"])
    fnd = W.found.to_numpy() == 1
    cx, cy = np.where(fnd, W.cx, np.nan), np.where(fnd, W.cy, np.nan)
    with np.errstate(invalid="ignore"):
        off = fnd & (np.hypot(cx - roll_med(cx, OFF_WIN, 100), cy - roll_med(cy, OFF_WIN, 100)) > OFF_EYE_PX)
    C = pd.read_csv(clean, usecols=["SmoothedDiameter_px"])
    med = float(np.nanmedian(pd.to_numeric(C.SmoothedDiameter_px, errors="coerce")))
    return len(W), 100.0 * float((~fnd | off).mean()), med


def num(v):
    try:
        x = float(v)
        return x if np.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def newest(pattern: str, folder: Path, day: int, first_last=re.compile(r"days?(\d+)-(\d+)\.csv$")):
    """newest table of a 'days<a>-<b>' series that includes day (files with '_before_' are backups)"""
    best = None
    for f in folder.glob(pattern):
        if "_before_" in f.name:
            continue
        m = first_last.search(f.name)
        if m and int(m.group(1)) <= day <= int(m.group(2)):
            if best is None or f.stat().st_mtime > best.stat().st_mtime:
                best = f
    return best


# ----------------------------------------------------------------------------------------------- context
class Ctx:
    def __init__(self, days):
        self.ff = ffmpeg()
        self.root = CC.active_pupil_root(P)
        self.lists = {k: CC.approved_sessions(self.root, n, warn=lambda *a: None) for k, n in CC.LISTS.items()}
        M = pd.read_csv(P / "cohort_manifest.csv", dtype=str).fillna("")
        self.group = {(r.cage_folder, r.mouse): r.group for r in M.itertuples()}
        self.partner = {}
        for r in M.itertuples():
            other = [q for q in M.itertuples() if q.pair_id == r.pair_id and q.mouse != r.mouse]
            self.partner[(r.cage_folder, r.mouse)] = (other[0].cage_folder, other[0].mouse) if other else None
        self.protocol = pd.read_csv(T / "protocol_expected.csv", dtype=str).fillna("")
        sp = P / "sync" / "sync_summary.csv"
        S = pd.read_csv(sp, dtype=str, encoding="utf-8-sig").fillna("") if sp.exists() else pd.DataFrame()
        self.sync = {r.day: r for r in S.itertuples()} if len(S) else {}
        vf = P / "analysis" / "pupil_check" / "pupil_check_verdicts.csv"
        V = pd.read_csv(vf, dtype=str, encoding="utf-8-sig").fillna("") if vf.exists() else pd.DataFrame()
        self.verdict = {(int(r.day), r.cage.strip().lower(), r.mouse.strip().lower()): (r.verdict.strip(), r.note.strip())
                        for r in V.itertuples()} if len(V) else {}
        tf = P / "assays" / "tail_immersion_all.csv"
        TI = pd.read_csv(tf, dtype=str, encoding="utf-8-sig").fillna("") if tf.exists() else pd.DataFrame()
        self.ti_all = {(int(r.day), r.cage, r.mouse): num(r.latency_s) for r in TI.itertuples()} if len(TI) else {}
        runs = sorted((P / "step02_base" / "longitudinal_outputs").glob("run_*/session_inventory.csv"),
                      key=lambda f: f.stat().st_mtime)
        self.step02_run = runs[-1].parent.name if runs else ""
        I = pd.read_csv(runs[-1], dtype=str).fillna("") if runs else pd.DataFrame()
        self.inv = {(int(r.day_index), r.cage_folder, r.mouse_id): r for r in I.itertuples()} if len(I) else {}
        # dates of every day folder's logs (A3) and the content index of all day folders (A7)
        self.dates, self.logs, self.hash = {}, {}, {}
        for dd in sorted((P / "python").glob("day*")):
            if not dd.name[3:].isdigit():
                continue
            n = int(dd.name[3:])
            for f in dd.rglob("*.jsonl"):
                L = read_log(f)
                self.logs[f] = L
                if L["start"] is not None:
                    self.dates.setdefault(n, set()).add(L["start"].date())
        for kind, pats in (("python", ("*.jsonl",)), ("salae", ("*.sal", "*.csv")), ("video", ("*.mp4", "*.avi"))):
            for pat in pats:
                for f in (P / kind).glob(f"day*/*/*/{pat}"):
                    if f.name.startswith("_"):
                        continue
                    self.hash.setdefault(full_hash(f) if f.suffix == ".jsonl" else quick_hash(f), []).append(f)
        # days whose Step 02 table keeps only the first N s (tools/step02_window.py, added 2026-11-26 for day 13)
        wf = T / "step02_window_days.csv"
        W = pd.read_csv(wf, dtype=str, encoding="utf-8-sig").fillna("") if wf.exists() else pd.DataFrame()
        self.window = {int(r.day): float(r.seconds) for r in W.itertuples() if str(r.approved).strip()} \
            if len(W) else {}
        lp = T / "local_paths.json"
        L = json.loads(lp.read_text(encoding="utf-8")) if lp.exists() else {}
        self.raw_dirs = [Path(p) for p in L.get("raw_dirs", [])]
        self.assay_dirs = list(dict.fromkeys([Path(p) for p in L.get("assay_video_dirs", [])] + self.raw_dirs +
                                             [HOME / "Videos" / "Logitech" / "LogiCapture"]))
        self.orig_dirs = [HOME / "Downloads" / "behavior_task_windows",
                          HOME / "Downloads" / "behavior_task_windows" / "behavior_task_windows"]
        self.rule_ok = all(s in (T / "pupil_check_examples.py").read_text(encoding="utf-8") for s in RULE_LINES)
        nf = P / "straub_video" / "name_fixes.csv"          # recordings saved under a wrong name (user-confirmed)
        F = pd.read_csv(nf, dtype=str, encoding="utf-8-sig").fillna("") if nf.exists() else pd.DataFrame()
        self.name_fix = {r.file.strip().lower(): (r.cage.strip(), r.mouse.strip(), int(r.day))
                         for r in F.itertuples() if r.approved.strip()} if len(F) else {}

    def expected_paradigms(self, day, cage, mouse):
        rows = [r for r in self.protocol.itertuples() if num(r.from_day) and num(r.to_day)
                and int(num(r.from_day)) <= day <= int(num(r.to_day))]
        if not rows:
            return None
        if len(rows) > 1:
            passive = self.group.get((cage, mouse), "") == "Passive"
            rows = [r for r in rows if (r.phase.strip().lower() == "passive") == passive] or rows
        return set(rows[0].paradigm.split("|"))


# ----------------------------------------------------------------------------------------------- the checks
def evaluate_day(day: int, X: Ctx) -> list:
    rows = []

    def add(check, cage, mouse, status, detail):
        rows.append(dict(day=day, check=check, cage=cage, mouse=mouse, status=status, detail=detail))

    dset = X.dates.get(day, set())
    others = {d for n, s in X.dates.items() if n != day for d in s}
    if len(dset) == 1 and not (dset & others):
        add("A3 log", "-", "-", "PASS", f"all logs of day {day} are from {next(iter(dset))}")
    elif dset:
        add("A3 log", "-", "-", "FAIL", f"day {day} logs have dates {sorted(map(str, dset))}"
            + (f"; shared with other days: {sorted(map(str, dset & others))}" if dset & others else ""))
    day_date = next(iter(dset)) if len(dset) == 1 else None
    cands = [f for f in (newest("pupil_check_days*.csv", P / "analysis" / "pupil_check", day),
                         P / "analysis" / "pupil_check" / f"pupil_check_day{day}.csv") if f is not None and f.exists()]
    pc_file = max(cands, key=lambda f: f.stat().st_mtime) if cands else P / "analysis" / "pupil_check" / "none.csv"
    PC = pd.read_csv(pc_file, dtype=str, encoding="utf-8-sig").fillna("") if pc_file.exists() else pd.DataFrame()
    pcheck = {(int(r.day), r.cage, r.mouse): r for r in PC.itertuples()} if len(PC) else {}
    prf = P / "analysis" / "days" / f"day{day}" / f"pr_checks_day{day}.csv"
    PR = pd.read_csv(prf, dtype=str, encoding="utf-8-sig").fillna("") if prf.exists() else pd.DataFrame()
    ptf = newest("pr_table_days*.csv", P / "analysis" / "days", day)
    PT = pd.read_csv(ptf, dtype=str, encoding="utf-8-sig").fillna("") if ptf else pd.DataFrame()
    tif = newest("ti_table_days*.csv", P / "analysis" / "days", day)
    TT = pd.read_csv(tif, dtype=str, encoding="utf-8-sig").fillna("") if tif else pd.DataFrame()
    stf = P / "analysis" / "days" / f"day{day}" / f"straub_day{day}.csv"
    ST = pd.read_csv(stf, dtype=str, encoding="utf-8-sig").fillna("") if stf.exists() else pd.DataFrame()
    in_step02 = any(k[0] == day for k in X.inv)

    # assay videos of the day in the raw folders (A8, per mouse + unnamed list)
    assay = {}
    unnamed = []
    for root in X.assay_dirs:
        if not root.is_dir():
            continue
        for f in [*root.glob("*.*"), *root.glob("*/*.*")]:
            if f.suffix.lower() not in (".mp4", ".avi", ".mov"):
                continue
            mu = UNNAMED.match(f.name)
            if mu:
                if day_date is not None and mu.group(1) == str(day_date):
                    unnamed.append(f)
                continue
            low = f.name.lower()
            kind = next((k for k in ("straub", "tst", "hot") if k in low), None)
            c, m, d = X.name_fix.get(low) or who(f.name)
            if kind and c and d == day:
                assay.setdefault((kind, c, m), []).append(f)
    if unnamed:
        add("A8 assays", "-", "-", "INFO", f"{len(unnamed)} recording(s) of {day_date} without a mouse name (not used): "
            + ", ".join(sorted(u.name for u in unnamed)))

    for cage, mouse in MICE:
        key = CC.norm_key(day, cage, mouse)
        pyd, sad, vid = (P / k / f"day{day}" / cage / mouse for k in ("python", "salae", "video"))
        logs = sorted(pyd.glob("*.jsonl")) if pyd.is_dir() else []
        sals = sorted(sad.glob("*.sal")) if sad.is_dir() else []
        digs = sorted(sad.glob("*.csv")) if sad.is_dir() else []
        vids = sorted([*vid.glob("*.mp4"), *vid.glob("*.avi")]) if vid.is_dir() else []
        novid_ok = key in X.lists.get(CC.NOVID, set())
        # A1
        prob = [f"{len(x)} {w}" for x, w in ((logs, "python logs"), (sals, ".sal"), (digs, "Saleae csv")) if len(x) != 1]
        if len(vids) != 1 and not (len(vids) == 0 and novid_ok):
            prob.append(f"{len(vids)} videos")
        if not (logs or sals or vids):
            add("A1 present", cage, mouse, "FAIL", "no files of this mouse-day in python/salae/video")
            continue
        add("A1 present", cage, mouse, "FAIL" if prob else "PASS",
            "; ".join(prob) if prob else ("log + .sal + csv + video" if vids else "log + .sal + csv, no video (listed)"))
        # A2
        bad = []
        for f in [*logs, *sals, *vids]:
            c, m, d = who(f.name)
            if (c, m) != (cage, mouse) or d != day:
                bad.append(f"{f.name} -> {c} {m} day {d}")
        add("A2 names", cage, mouse, "FAIL" if bad else "PASS", "; ".join(bad) if bad else "names = folder")
        if not logs:
            continue
        L = X.logs.get(logs[0]) or read_log(logs[0])
        replay = "replay" in L["paradigm"].lower()
        # A3
        msg, st = [], "PASS"
        if L["start"] is None:
            msg.append("no start time in the log")
            st = "FAIL"
        exp = X.expected_paradigms(day, cage, mouse)
        if exp is not None and L["paradigm"] not in exp:
            msg.append(f"paradigm {L['paradigm'] or '?'} (expected {' or '.join(sorted(exp))})")
            st = "FAIL"
        if replay:
            pc_, pm_, pd_ = who(L["replay"])
            want = X.partner.get((cage, mouse))
            if want is None or (pc_, pm_) != want or pd_ != day:
                msg.append(f"replays {L['replay'] or '?'} (manifest partner {want[0] + ' ' + want[1] if want else '?'}, "
                           f"day {day})")
                st = "FAIL"
        if not L["ended"]:
            msg.append("no experiment_end (stopped by hand?)")
            st = "FAIL" if st == "FAIL" else "WARN"
        if L["bad"]:
            msg.append(f"{L['bad']} unreadable line(s)")
            st = "FAIL" if st == "FAIL" else "WARN"
        add("A3 log", cage, mouse, st, "; ".join(msg) if msg else
            f"{L['paradigm']}" + (f", replays {Path(L['replay']).stem}" if replay else "") + f", {L['session_s']:.0f} s")
        # A4 / video counts
        vi = video_info(vids[0], X.ff) if vids else None
        if vi is not None:
            if vi["created"] is None or L["start"] is None:
                add("A4 video", cage, mouse, "WARN", "no creation time in the video (or no log start)")
            else:
                dv = (vi["created"] - L["start"]).total_seconds()
                st = "PASS" if abs(dv) <= VIDEO_PASS_S else ("WARN" if abs(dv) <= VIDEO_WARN_S else "FAIL")
                add("A4 video", cage, mouse, st, f"video created {dv:+.1f} s from the python start"
                    + ("" if st == "PASS" else " - wrong file pair?" if st == "FAIL" else ""))
        # A5
        dg = read_digital(digs[0]) if digs else None
        cam = (dg["camera"][0] if dg and dg["camera"] is not None else np.array([]))
        off = np.nan
        if dg is not None:
            msg, st = [], "PASS"
            ss = sal_start(sals[0]) if sals else None
            if ss is not None and L["start"] is not None:
                rel = (L["start"] - ss).total_seconds()
                if rel < 0:
                    msg.append(f"Saleae started {-rel:.1f} s AFTER the python start")
                    st = "FAIL"
                elif np.isfinite(L["session_s"]) and rel + L["session_s"] > dg["end"] + 1:
                    msg.append(f"capture ends {rel + L['session_s'] - dg['end']:.0f} s before the python end")
                    st = "WARN"
                else:
                    msg.append(f"python start {rel:.1f} s into the capture")
            else:
                msg.append("no .sal start time")
                st = "WARN"
            if dg["valve"] is None:
                msg.append("no valve channel in the csv")
                st = "FAIL"
            else:
                r, dur = dg["valve"]
                rv = r[(dur >= VALVE_MIN_S) & (dur <= VALVE_MAX_S)]
                mv = match_valve(L["rewards"], rv, L["session_s"])
                off = mv["offset"]
                n = len(L["rewards"])
                gone = 0
                if ss is not None and L["start"] is not None and np.isfinite(mv["offset"]):
                    gone = int((L["rewards"] + mv["offset"] > dg["end"]).sum())   # after the capture end
                if mv["matched"] == n and mv["in_window"] == n:
                    msg.append(f"valve openings = python rewards ({n}, all matched)")
                elif mv["matched"] + gone == n and mv["in_window"] == mv["matched"]:
                    msg.append(f"{mv['matched']} of {n} rewards matched; {gone} after the capture end")
                    st = "FAIL" if st == "FAIL" else "WARN"
                elif mv["matched"] == n:
                    msg.append(f"all {n} rewards matched, but {mv['in_window'] - n} extra valve opening(s) inside the "
                               "session")
                    st = "FAIL" if st == "FAIL" else "WARN"
                else:
                    msg.append(f"only {mv['matched']} of {n} python rewards have a valve opening "
                               f"({mv['in_window']} openings in the session)")
                    st = "FAIL"
            add("A5 saleae", cage, mouse, st, "; ".join(msg))
        # A6
        origs = [d / logs[0].name for d in X.orig_dirs if (d / logs[0].name).is_file()]
        if not origs:
            add("A6 copy", cage, mouse, "INFO", "original log not found in Downloads (moved / renamed)")
        else:
            O = None if full_hash(origs[0]) == full_hash(logs[0]) else read_log(origs[0])
            if O is None:
                add("A6 copy", cage, mouse, "PASS", "copy = original")
            elif O["start"] is not None and L["start"] is not None and O["start"] != L["start"]:
                add("A6 copy", cage, mouse, "WARN",
                    f"{origs[0].name} in Downloads is ANOTHER run that reused the name (start {O['start']:%H:%M:%S}, "
                    f"{len(O['rewards'])} rewards); the project copy is the session of {L['start']:%H:%M:%S} "
                    f"({len(L['rewards'])} rewards)")
            else:
                add("A6 copy", cage, mouse, "FAIL", f"copy differs from {origs[0]} with the same start (cut or changed)")
        # A7
        dup = []
        for f in [*logs, *sals, *digs, *vids]:
            hk = full_hash(f) if f.suffix == ".jsonl" else quick_hash(f)
            dup += [str(g.relative_to(P)) for g in X.hash.get(hk, []) if g != f]
        add("A7 unique", cage, mouse, "FAIL" if dup else "PASS",
            ("same content as " + "; ".join(sorted(set(dup)))) if dup else "no copy of these files elsewhere")
        # A8 (TI + assay videos of this mouse)
        msg, st = [], None
        ad = P / "assays" / f"day{day}" / cage / mouse
        tif_ = sorted(ad.glob("*immer*.csv")) if ad.is_dir() else []
        tv = None
        if tif_:
            first = tif_[-1].read_text(encoding="utf-8-sig").splitlines()[0].split(",")
            tv = num(first[1]) if len(first) > 1 else None
        ta = X.ti_all.get((day, cage, mouse))
        if tv is None and ta is None:
            pass
        elif tv is not None and ta is not None and abs(tv - ta) < 1e-9:
            msg.append(f"TI {tv:g} s (file = table)")
            st = "PASS"
        else:
            msg.append(f"TI file {tv} vs tail_immersion_all {ta}")
            st = "FAIL"
        for kind in ("straub", "tst", "hot"):
            fs = assay.get((kind, cage, mouse), [])
            if len(fs) > 1:
                msg.append(f"{len(fs)} {kind} videos: " + ", ".join(f.name for f in fs))
                st = "FAIL"
            elif fs:
                import cv2
                c2 = cv2.VideoCapture(str(fs[0]))
                secs = c2.get(cv2.CAP_PROP_FRAME_COUNT) / (c2.get(cv2.CAP_PROP_FPS) or 30)
                c2.release()
                short = secs < ASSAY_MIN_S[kind] - 1
                msg.append(f"{kind} {fs[0].name} {secs:.0f} s" + (f" (< {ASSAY_MIN_S[kind]} s)" if short else ""))
                st = st if st == "FAIL" else ("WARN" if short else (st or "PASS"))
        if msg:
            add("A8 assays", cage, mouse, st or "INFO", "; ".join(msg))
        # B1
        sr = X.sync.get(f"day{day}/{cage}/{mouse}")
        if sr is None:
            add("B1 sync", cage, mouse, "FAIL", "session not in sync_summary.csv")
        else:
            probs = [q.strip() for q in str(sr.problems).split(" | ") if q.strip()]
            listed = key in X.lists.get(CC.SYNC, set()) or key in X.lists.get(CC.NOVID, set())
            known = [q for q in probs if (replay and "lick" in q.lower())
                     or (listed and re.search(r"video has \d+ (fewer|MORE) frames|video lost \d+ frames|no video file", q))]
            rest = [q for q in probs if q not in known]
            if sr.status == "OK":
                add("B1 sync", cage, mouse, "PASS", "OK")
            elif not rest:
                add("B1 sync", cage, mouse, "INFO", "known, handled: " + " | ".join(
                    "passive lick line (intraoral spout, not licks)" if "lick" in q.lower() else q[:90] for q in known))
            else:
                add("B1 sync", cage, mouse, "WARN", sr.status + ": " + " | ".join(rest)[:200])
        # B2 (camera pulses inside the python session: python start - 0.05 s .. end + 0.1 s, as sync_sessions.py)
        if np.isfinite(off) and np.isfinite(L["session_s"]) and len(cam):
            cam = cam[(cam >= off - 0.05) & (cam <= off + L["session_s"] + 0.1)]
        pcr = pcheck.get((day, cage, mouse))
        table_run = pcr.table_run if pcr is not None else ""
        if vi is None:
            add("B2 alignment", cage, mouse, "PASS" if novid_ok else "FAIL",
                "no video; behaviour-only (listed)" if novid_ok else "no video and not in behavioronly_sessions.csv")
        else:
            msg, st = [], "PASS"
            sp = num(sr.camera_pulses) if sr is not None else None
            if sr is not None and (str(vi["frames"]) != sr.video_frames or sp is None or abs(sp - len(cam)) > 1):
                msg.append(f"routine counted video {sr.video_frames} / pulses {sr.camera_pulses}, evaluator "
                           f"{vi['frames']} / {len(cam)}")
                st = "FAIL"
            npul = int(sp) if sp is not None and abs(sp - len(cam)) <= 1 else len(cam)   # +-1 pulse at the window edge
            if vi["frames"] == npul:
                msg.append(f"video frames = camera pulses ({npul})")
            elif key in X.lists.get(CC.SYNC, set()):
                ok = table_run.endswith("_syncaligned")
                msg.append(f"video {vi['frames']} vs pulses {npul}: in alignment_fixes.csv, table "
                           f"{table_run or '?'}" + ("" if ok else " is NOT the sync-aligned copy"))
                st = st if ok else "FAIL"
            else:
                msg.append(f"video {vi['frames']} vs pulses {npul} and NOT in alignment_fixes.csv")
                st = "FAIL"
            add("B2 alignment", cage, mouse, st, "; ".join(msg))
        # B3 / B4
        if pcr is None:
            if vi is not None:
                add("B3 pupil", cage, mouse, "FAIL", f"no row in {pc_file.name}")
        elif vi is None:
            add("B3 pupil", cage, mouse, "INFO", "no video: no pupil (behaviour-only session)")
        else:
            sd = CC.session_dir(X.root, day, cage, mouse)
            src = pcr.source.lower()
            got = None
            try:
                if src.startswith("package") and sd is not None:
                    got = package_out(sd / "pupil_outputs_general" / "pupil_diameter_FINAL_U-Net_results_CLEAN.csv")
                elif src.startswith("whole eye") and sd is not None:
                    got = eye_out(P / "analysis" / "wholeeye" / f"day{day}" / cage / mouse / "wholeeye_frames.csv",
                                  CC.wholeeye_pupil_csv(sd, X.root))
            except Exception as exc:                     # noqa: BLE001
                add("B3 pupil", cage, mouse, "FAIL", f"pupil file could not be read again: {exc}")
                got = False
            if got is None:
                add("B3 pupil", cage, mouse, "INFO", f"source '{pcr.source}': not recomputed here")
            elif got:
                n, out, med = got
                po, pm, pn = num(pcr.out_pct), num(pcr.median_px), num(pcr.n_frames)
                msg, st = [], "PASS"
                if po is None or abs(out - po) > OUT_TOL_PP:
                    msg.append(f"out {out:.2f} % here vs {pcr.out_pct} % in the check")
                    st = "FAIL"
                if pm is None or abs(med - pm) > MED_TOL * max(pm, 1e-6):
                    msg.append(f"median {med:.2f} px here vs {pcr.median_px} px")
                    st = "FAIL"
                if pn is None or int(pn) != n:
                    msg.append(f"{n} frames here vs {pcr.n_frames} in the check")
                    st = "FAIL"
                if vi is not None and n != vi["frames"]:
                    msg.append(f"pupil file has {n} rows but the video {vi['frames']} frames")
                    st = "FAIL"
                add("B3 pupil", cage, mouse, st, "; ".join(msg) if msg else
                    f"{pcr.source.split(' (')[0]}: out {out:.2f} %, median {med:.1f} px, {n} frames = check")
            vd, note = X.verdict.get((day, cage, mouse), ("", ""))
            po, wtxt = num(pcr.out_pct) or 0.0, ""
            p910 = num(getattr(pcr, "out_first910s_pct", None))
            if X.window.get(day) == 910.0 and p910 is not None:
                # 2026-11-26 (day 13 = 20 min): the verdict is judged on the analysed part, the first 15 min 10 s
                # (pupil_check_examples' out_first910s_pct), as every other measure of a window day
                po, wtxt = p910, " in the first 15:10"
            if not vd or vd.upper() == "LOOK":
                add("B4 verdict", cage, mouse, "WARN", "no verdict written yet (pupil_check_verdicts.csv)")
            elif po >= NOISY_PCT and vd.lower() == "ok" and not re.search(r"decision|kept|user|approved", note, re.I):
                add("B4 verdict", cage, mouse, "WARN", f"out {po:.1f} %{wtxt} but verdict 'ok' without a decision note")
            else:
                add("B4 verdict", cage, mouse, "PASS" if vd.lower() == "ok" else "INFO", f"{vd} (out {po:.2f} %{wtxt})")
        # B5 / C6 (Step 02)
        iv = X.inv.get((day, cage, mouse))
        if not in_step02:
            add("B5 table", cage, mouse, "WARN", f"Step 02 {X.step02_run or '-'} does not have day {day} yet")
        elif iv is None:
            add("B5 table", cage, mouse, "FAIL", f"session missing in Step 02 {X.step02_run}")
        else:
            sdn = Path(iv.session_dir).name
            ok = (not table_run) or sdn == table_run
            add("B5 table", cage, mouse, "PASS" if ok else "FAIL",
                f"Step 02 {X.step02_run} read {sdn}" + ("" if ok else f", the check used {table_run}"))
            msg, st = [], "PASS"
            want_rows = vi["frames"] if vi is not None else len(cam)
            what = "video frames" if vi is not None else "camera pulses"
            win = X.window.get(day)
            if win is not None:
                # 2026-11-26 (day 13 = 20-min sessions; user: PR = the whole session, every other measure = the first
                # 15 min 10 s): step02_window.py keeps the rows with CamTime_rel_s < win in a listed day's Step 02
                # table -> expected rows = camera pulses within win s of the first pulse (pulse k = row k); the marker
                # next to the table must name the whole table = all video frames
                full_rows = want_rows
                want_rows = int(np.sum((cam - cam[0]) < win)) if len(cam) else want_rows
                what = f"camera pulses of the first {win:g} s (whole table {full_rows} {what})"
                try:
                    M = json.loads((Path(iv.session_dir) / "combined_pupil_digital.window.json").read_text(encoding="utf-8"))
                except (OSError, ValueError):
                    M = None
                if M is None:
                    msg.append("no window marker next to the Step 02 table")
                    st = "FAIL"
                elif int(num(M.get("rows_whole")) or -1) != full_rows:
                    msg.append(f"window marker: whole table {M.get('rows_whole')} rows vs {full_rows} {what}")
                    st = "FAIL"
            if num(iv.rows) is None or abs(int(num(iv.rows)) - want_rows) > (1 if win is not None else 0):
                msg.append(f"{iv.rows} rows vs {want_rows} {what}")
                st = "FAIL"
            elif win is not None:
                msg.append(f"{iv.rows} rows = {what}")
            if str(iv.isPassive).strip() not in ("1" if replay else "0",):
                msg.append(f"isPassive {iv.isPassive} but the log is {'a replay' if replay else 'not a replay'}")
                st = "FAIL"
            if Path(iv.jsonl_file).name != logs[0].name:
                msg.append(f"Step 02 log {Path(iv.jsonl_file).name}")
                st = "FAIL"
            add("C6 step02", cage, mouse, st, "; ".join(msg) if msg else f"rows {iv.rows}, log + paradigm match")
        # C1
        msg, st = [], "PASS"
        n, req = len(L["rewards"]), L["last_req"]
        pr = PR[(PR.cage == cage) & PR.mouse.str.endswith("_" + mouse)] if len(PR) else PR
        if len(PR) and pr.empty:
            msg.append("no pr_checks row")
            st = "FAIL"
        elif len(pr):
            r = pr.iloc[-1]
            if num(r.rewards_python) != n:
                msg.append(f"pr_checks rewards {r.rewards_python} vs log {n}")
                st = "FAIL"
            if not replay and num(r.pr_RL_full) is not None and num(r.pr_RL_full) != req:
                msg.append(f"pr_checks PR {r.pr_RL_full} vs log {req:g}")
                st = "FAIL"
        pt = PT[(PT.day == str(day)) & (PT.cage == cage) & (PT.mouse == mouse)] if len(PT) else PT
        if len(pt):
            r = pt.iloc[-1]
            if num(r.rewards_python) != n:
                msg.append(f"PR table rewards {r.rewards_python} vs log {n}")
                st = "FAIL"
            if not replay and num(r.pr) is not None and num(r.pr) != req:
                msg.append(f"PR table {r.pr} vs log {req:g}")
                st = "FAIL"
        if iv is not None:
            if not replay and num(iv.RequirementLast) is not None and num(iv.RequirementLast) != req:
                msg.append(f"Step 02 RequirementLast {iv.RequirementLast} vs log {req:g}")
                st = "FAIL"
            k = len([x for x in str(iv.RewardTimes_s).split(",") if x.strip()])
            if str(iv.RewardTimes_s).strip() and k != n:
                msg.append(f"Step 02 {k} reward times vs log {n}")
                st = "FAIL"
        add("C1 PR", cage, mouse, st, "; ".join(msg) if msg else
            (f"rewards {n}" + ("" if replay else f", last requirement {req:g}") + " = pr_checks = PR table"
             + (" = Step 02" if iv is not None else "")))
        # C2
        if replay:
            pa = X.partner.get((cage, mouse))
            plogs = sorted((P / "python" / f"day{day}" / pa[0] / pa[1]).glob("*.jsonl")) if pa else []
            if len(plogs) != 1:
                add("C2 passive", cage, mouse, "FAIL", f"partner {pa} log not found once ({len(plogs)})")
            else:
                pn = len((X.logs.get(plogs[0]) or read_log(plogs[0]))["rewards"])
                add("C2 passive", cage, mouse, "PASS" if pn == n else "FAIL",
                    f"drops {n} {'=' if pn == n else '!='} partner {pa[0]} {pa[1]} rewards {pn}")
        # C3
        if len(pr):
            r = pr.iloc[-1]
            a_, b_ = num(r.rewards_saleae), num(r.rewards_python)
            add("C3 rewards", cage, mouse, "PASS" if a_ is not None and a_ == b_ else "FAIL",
                f"pr_checks Saleae {r.rewards_saleae} vs python {r.rewards_python}")
        # C4
        if tv is not None or ta is not None:
            msg, st = [], "PASS"
            tt = TT[(TT.day == str(day)) & (TT.cage == cage) & (TT.mouse == mouse)] if len(TT) else TT
            tval = num(tt.iloc[-1].latency_s) if len(tt) else None
            if tval is None or tv is None or abs(tval - tv) > 1e-9:
                msg.append(f"TI table {tval} vs assays file {tv}")
                st = "FAIL"
            if iv is not None:
                sv = num(iv.Immersion_Latency_s)
                if sv is None or tv is None or abs(sv - tv) > 1e-9:
                    msg.append(f"Step 02 Immersion {iv.Immersion_Latency_s} vs file {tv}")
                    st = "FAIL"
            elif in_step02:
                pass
            add("C4 TI", cage, mouse, st, "; ".join(msg) if msg else
                f"{tv:g} s: file = TI table" + (" = Step 02" if iv is not None else ""))
        # C5
        sdir = P / "analysis" / "straub_auto" / f"day{day}"
        fr = [f for f in sdir.glob("*_frames.csv") if who(f.name)[:2] == (cage, mouse)] if sdir.is_dir() else []
        if len(fr) > 1:
            add("C5 straub", cage, mouse, "FAIL", f"{len(fr)} automatic Straub results: " + ", ".join(f.name for f in fr))
        elif fr:
            D = pd.read_csv(fr[0], usecols=["straub"])
            s = pd.to_numeric(D.straub, errors="coerce")
            mine = round(100 * float(s[s.notna()].mean()), 1) if s.notna().any() else None
            sr_ = ST[(ST.cage == cage) & (ST.mouse == mouse)] if len(ST) else ST
            tab = num(sr_.iloc[-1].get("auto_straub_pct")) if len(sr_) else None
            ok = mine is not None and tab is not None and abs(mine - tab) <= 0.1
            add("C5 straub", cage, mouse, "PASS" if ok else "FAIL",
                f"auto Straub {mine} % from {fr[0].name} vs table {tab}")
    return rows


# ----------------------------------------------------------------------------------------------- outputs
def summary_lines(day, R, X):
    c = R.status.value_counts().to_dict() if len(R) else {}
    head = (f"EVALUATOR day {day}: {c.get('FAIL', 0)} FAIL, {c.get('WARN', 0)} WARN, {c.get('PASS', 0)} PASS, "
            f"{c.get('INFO', 0)} INFO  (independent re-read of the raw files; Step 02 {X.step02_run or '-'})")
    out = [head]
    if not X.rule_ok:
        out.append("  !! the pupil check's rule in pupil_check_examples.py changed - update this file's constants")
    for s in ("FAIL", "WARN"):
        for r in R[R.status == s].itertuples():
            who_ = "day" if r.cage == "-" else f"{r.cage} {r.mouse}"
            out.append(f"  {s} {r.check:13s} {who_:13s} {r.detail}")
    return out


def figure(day, R, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = [f"{c} {m}" for c, m in MICE] + ["day"]
    col_of = {f"{c} {m}": i for i, (c, m) in enumerate(MICE)}
    face = {"FAIL": "#f4b6b6", "WARN": "#ffe8a3", "PASS": "#cfe8cf", "INFO": "#e3e7eb"}
    word = {"FAIL": "FAIL", "WARN": "look", "PASS": "ok", "INFO": "i"}
    bad = R[R.status.isin(["FAIL", "WARN"])]
    lines = []
    for r in bad.itertuples():
        txt = f"{r.status} {r.check} {'day' if r.cage == '-' else r.cage + ' ' + r.mouse}: {r.detail}"
        lines += textwrap.wrap(txt, 175)[:2]
    lines = lines[:26]
    h = 0.34 * (len(ORDER) + 2) + 0.2 * len(lines) + 1.2
    fig = plt.figure(figsize=(15, h), dpi=110)
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
    ax.axis("off")
    ny = len(ORDER) + 1
    top = h - 0.9
    ax.set_xlim(0, 15)
    ax.set_ylim(0, h)
    c = R.status.value_counts().to_dict()
    ax.text(0.1, h - 0.4, f"Evaluator day {day}  -  {c.get('FAIL', 0)} FAIL · {c.get('WARN', 0)} look · "
            f"{c.get('PASS', 0)} ok  (independent re-read of the raw files: files read right? QC right? analysis right?)",
            fontsize=12, fontweight="bold", color="#1f3a5f", va="center")
    x0, cw, rh = 2.3, (15 - 2.4) / len(cols), 0.34
    for j, cn in enumerate(cols):
        ax.add_patch(plt.Rectangle((x0 + j * cw, top - rh), cw, rh, facecolor="#1f3a5f", edgecolor="white"))
        ax.text(x0 + (j + 0.5) * cw, top - rh / 2, cn, ha="center", va="center", fontsize=8.5, color="white")
    for i, ck in enumerate(ORDER, start=1):
        y = top - (i + 1) * rh
        ax.text(0.1, y + rh / 2, ck, va="center", fontsize=9, color="#1f2933")
        g = R[R.check == ck]
        for j, cn in enumerate(cols):
            if cn == "day":
                gg = g[g.cage == "-"]
            else:
                cg, ms = cn.split(" ")
                gg = g[(g.cage == cg) & (g.mouse == ms)]
            if gg.empty:
                ax.add_patch(plt.Rectangle((x0 + j * cw, y), cw, rh, facecolor="white", edgecolor="#d5d9de"))
                continue
            s = max(gg.status, key=lambda v: RANK[v])
            ax.add_patch(plt.Rectangle((x0 + j * cw, y), cw, rh, facecolor=face[s], edgecolor="#d5d9de"))
            ax.text(x0 + (j + 0.5) * cw, y + rh / 2, word[s], ha="center", va="center", fontsize=8.5,
                    fontweight="bold" if s == "FAIL" else "normal")
    y = top - (ny + 1) * rh - 0.15
    for ln in lines or ["no FAIL or look items"]:
        ax.text(0.1, y, ln, fontsize=7.6, color="#7a1f1f" if ln.startswith("FAIL") else "#5a4a00", va="top",
                family="monospace")
        y -= 0.2
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--day", type=int)
    ap.add_argument("--days", type=int, nargs="+")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    days = a.days or ([a.day] if a.day else None)
    if not days:
        ap.error("give --day N or --days ...")
    X = Ctx(days)
    worst = 0
    for d in days:
        R = pd.DataFrame(evaluate_day(d, X), columns=["day", "check", "cage", "mouse", "status", "detail"])
        R["order"] = R.check.map({k: i for i, k in enumerate(ORDER)})
        R = R.sort_values(["order", "cage", "mouse"]).drop(columns="order")
        od = P / "analysis" / "days" / f"day{d}"
        od.mkdir(parents=True, exist_ok=True)
        R.to_csv(od / f"evaluator_day{d}.csv", index=False, encoding="utf-8-sig")
        lines = summary_lines(d, R, X)
        (od / f"evaluator_day{d}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
        figure(d, R, od / f"evaluator_day{d}.png")
        print("\n".join(lines if not a.quiet else lines[:1]))
        print(f"  -> {od / f'evaluator_day{d}.csv'} / .txt / .png")
        worst = max(worst, int((R.status == "FAIL").any()))
    return worst


if __name__ == "__main__":
    raise SystemExit(main())
