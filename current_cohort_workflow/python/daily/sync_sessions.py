#!/usr/bin/env python3
"""
sync_sessions.py - align Python behaviour logs, Saleae captures and camera videos.

For every session it finds (same file name in the same day folder of the
<batch>_python, <batch>_salae and <batch>_video folders) it writes, into
<batch>_sync/<day>/:

    <session>_events.csv   every lick / reward / python event with its time and video frame
    <session>_frames.csv   every camera pulse with its real time and its video frame
    <session>_report.txt   the sync checks for that session (status OK / WARN / FAIL)

and <batch>_sync/sync_summary.csv with one line per session.

Master clock = the Saleae capture. Video frame = camera pulse index (minus a video offset
if the recording started late).
Usage:  python sync_sessions.py [--root FOLDER] [--force]
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import zipfile
from bisect import bisect_left, bisect_right
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

TOOL_VERSION = "1.1"

# ----------------------------------------------------------------------------- settings
CHANNEL_NAMES = {"camera": "camera", "lick": "lick", "reward": "reward"}  # Saleae channel names
VIDEO_EXTS = (".mp4", ".avi", ".mkv", ".mov")
REWARD_TOL_MS = 50.0          # python reward_cmd vs Saleae valve edge (before the fit)
DRIFT_PPM = 300.0             # max believable clock drift PC vs Saleae
LICK_DELAY0 = 0.016           # typical python lick delay (s) behind the Saleae edge
LICK_BAND = 0.030             # +- band around the running delay estimate when matching licks
FAIL_RESID_MS = 5.0           # rewards further than this from the fit -> FAIL
TRAIN_GAP_S = 1.0             # camera pulses separated by more than this start a new pulse train
VIDEO_START_TOL_S = 120.0     # video file creation time vs python start (pairing check)
VALVE_LICK_PAD_S = 0.03       # a lick overlapping valve-open (+30 ms) is flagged
BOUNCE_S = 0.001              # reward-channel pulses closer than this are one valve opening (contact bounce)
PY_ONLY_SKIP = {"meta", "lick", "reward_cmd"}
OFFSETS_FILE = "video_offsets.csv"   # optional, in <batch>_sync: day,session,first_pulse_in_video


# ----------------------------------------------------------------------------- helpers
def parse_iso(s: str):
    try:
        return datetime.fromisoformat(s.strip().replace("Z", "+00:00"))
    except Exception:
        return None


def find_ffmpeg():
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:
        from shutil import which
        return which("ffmpeg")


def fmt(x, nd=6):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    if isinstance(x, (bool, np.bool_)):
        return "1" if x else "0"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, (float, np.floating)):
        return f"{float(x):.{nd}f}"
    return str(x)


# ----------------------------------------------------------------------------- readers
def read_python_log(path: Path) -> dict:
    events = []
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                events.append(json.loads(line))
            except json.JSONDecodeError:
                continue  # e.g. a half-written last line
    meta = next((e for e in events if e.get("tag") == "meta"), None)
    start_wall = parse_iso(meta["t"]) if meta and "t" in meta else None
    for e in events:  # older logs may lack 'elapsed'
        if "elapsed" not in e and start_wall and "t" in e:
            t = parse_iso(e["t"])
            if t:
                e["elapsed"] = (t - start_wall).total_seconds()
    events = [e for e in events if isinstance(e.get("elapsed"), (int, float))]
    return {"events": events, "meta": meta, "start_wall": start_wall}


def read_saleae_csv(path: Path) -> dict:
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        header = next(csv.reader(f))
    names = [h.strip().lower() for h in header]
    try:
        data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2, encoding="utf-8-sig")
    except ValueError as exc:
        raise ValueError(f"cannot read {path.name} as numbers - export the CSV with time in seconds "
                         f"(not ISO8601 timestamps): {exc}")
    t = data[:, 0]
    chans = {}
    for i, n in enumerate(names[1:], start=1):
        v = data[:, i]
        d = np.diff(v)
        rises = t[1:][d > 0]
        falls = t[1:][d < 0]
        starts_high = bool(v[0] > 0)
        intervals = []
        fi = 0
        for r in rises:
            while fi < len(falls) and falls[fi] < r:
                fi += 1
            intervals.append((float(r), float(falls[fi]) if fi < len(falls) else float("nan")))
        if n == CHANNEL_NAMES["reward"] and intervals:
            # valve contact bounce: a few-ns pulse right before the real opening -> merge pulses < 1 ms apart
            merged = [intervals[0]]
            for on, off in intervals[1:]:
                p_on, p_off = merged[-1]
                if not np.isnan(p_off) and on - p_off < BOUNCE_S:
                    merged[-1] = (p_on, off)
                else:
                    merged.append((on, off))
            intervals = merged
            rises = np.array([iv[0] for iv in intervals])
        chans[n] = {"rises": rises, "intervals": intervals, "starts_high": starts_high}
    return {"channels": chans, "names": names[1:], "end": float(t[-1])}


def read_sal_start(path: Path):
    """Wall-clock capture start stored inside a Logic 2 .sal file (software stamp, ~0.2 s accuracy)."""
    try:
        with zipfile.ZipFile(path) as z:
            meta = json.loads(z.read("meta.json"))
        st = meta["data"]["captureStartTime"]
        ms = float(st["unixTimeMilliseconds"]) + float(st.get("fractionalMilliseconds", 0))
        return datetime.fromtimestamp(ms / 1000.0, tz=timezone.utc)
    except Exception:
        return None


def probe_video(path: Path, ffmpeg) -> dict:
    out = {"frames": None, "fps_tag": None, "created": None, "error": None}
    if not ffmpeg:
        out["error"] = "ffmpeg not found (pip install imageio-ffmpeg)"
        return out
    try:
        head = subprocess.run([ffmpeg, "-hide_banner", "-i", str(path)],
                              capture_output=True, text=True, errors="replace").stderr
        m = re.search(r"Video:.*?(\d+(?:\.\d+)?) fps", head)
        out["fps_tag"] = float(m.group(1)) if m else None
        m = re.search(r"creation_time\s*:\s*(\S+)", head)
        out["created"] = parse_iso(m.group(1)) if m else None
        cnt = subprocess.run([ffmpeg, "-hide_banner", "-i", str(path), "-map", "0:v:0",
                              "-c", "copy", "-f", "null", "-"],
                             capture_output=True, text=True, errors="replace").stderr
        frames = re.findall(r"frame=\s*(\d+)", cnt)
        out["frames"] = int(frames[-1]) if frames else None
        if out["frames"] is None:
            out["error"] = "could not count frames"
    except Exception as exc:
        out["error"] = str(exc)
    return out


def read_offsets(out_root: Path) -> dict:
    p = out_root / OFFSETS_FILE
    res = {}
    if p.exists():
        try:
            with open(p, newline="", encoding="utf-8-sig") as f:
                for r in csv.DictReader(f):
                    try:
                        res[(r.get("day", "").strip(), r.get("session", "").strip())] = \
                            int(r["first_pulse_in_video"])
                    except (KeyError, ValueError):
                        continue
        except OSError:
            pass
    return res


# ----------------------------------------------------------------------------- matching
def match_rewards(py_rew, valve_rises, a0):
    """Sequentially match python reward_cmd times to Saleae valve rises (monotonic),
    re-fitting the clock as it goes. Returns {py_index: valve_index}."""
    matches = {}
    pairs = []
    a, b = a0, 1.0
    j0 = 0
    for i, e in enumerate(py_rew):
        pred = a + b * e
        tol = REWARD_TOL_MS / 1000.0 + DRIFT_PPM * 1e-6 * abs(e)
        lo = bisect_left(valve_rises, pred - tol, lo=j0)
        best = None
        k = lo
        while k < len(valve_rises) and valve_rises[k] <= pred + tol:
            if best is None or abs(valve_rises[k] - pred) < abs(valve_rises[best] - pred):
                best = k
            k += 1
        if best is not None:
            matches[i] = best
            j0 = best + 1
            pairs.append((e, float(valve_rises[best])))
            if len(pairs) >= 2 and pairs[-1][0] - pairs[0][0] > 30.0:
                ee = np.array([p[0] for p in pairs]); ss = np.array([p[1] for p in pairs])
                b, a = np.polyfit(ee, ss, 1)
                if abs(b - 1) > DRIFT_PPM * 1e-6:
                    b = 1.0
                    a = float(np.median(ss - ee))
            else:
                a = float(np.median([p[1] - p[0] for p in pairs]))
    return matches


def match_licks(py_t, rises, delay0=LICK_DELAY0):
    """Monotonic matching of python lick times (on the Saleae clock) to Saleae lick rises.
    The expected delay is tracked with a running median, so slow clock drift is absorbed.
    Returns {rise_index: py_index}."""
    matches = {}
    j0 = 0
    recent = []
    m = delay0
    for i, t in enumerate(py_t):
        lo_t, hi_t = t - (m + LICK_BAND), t - (m - LICK_BAND)
        j = bisect_left(rises, lo_t, lo=j0)
        best = None
        while j < len(rises) and rises[j] <= hi_t:
            if best is None or abs((t - rises[j]) - m) < abs((t - rises[best]) - m):
                best = j
            j += 1
        if best is not None:
            matches[best] = i
            j0 = best + 1
            recent.append(t - rises[best])
            if len(recent) > 25:
                recent.pop(0)
            m = float(np.median(recent))
    return matches


def pulse_trains(cam):
    """Indices where a new camera pulse train starts (gap > TRAIN_GAP_S)."""
    if len(cam) == 0:
        return []
    starts = [0] + [int(k) + 1 for k in np.where(np.diff(cam) > TRAIN_GAP_S)[0]]
    return starts


def choose_anchor(cam, valve_r, lick_r, py_rew, py_lick, py_dur):
    """Find where Python's elapsed=0 sits on the Saleae clock.
    Candidates: start of every camera pulse train, and the most common valve-reward_cmd offset."""
    cands = []
    starts = pulse_trains(cam)
    ends = [s - 1 for s in starts[1:]] + [len(cam) - 1]
    for s, e in zip(starts, ends):
        cands.append(("camera train", float(cam[s]), float(cam[e] - cam[s])))
    if len(py_rew) and len(valve_r):
        diffs = (valve_r[None, :] - py_rew[:, None]).ravel()
        diffs = diffs[np.isfinite(diffs)]
        if len(diffs):
            rb = np.round(diffs / 0.002)
            vals, counts = np.unique(rb, return_counts=True)
            mode = float(vals[np.argmax(counts)] * 0.002)
            close = diffs[np.abs(diffs - mode) <= 0.004]
            cands.append(("reward offset", float(np.median(close)), None))
    best = None
    for kind, a0, train_dur in cands:
        rm = match_rewards(py_rew, valve_r, a0) if len(py_rew) else {}
        lm = match_licks(a0 + py_lick, lick_r) if len(py_lick) else {}
        dur_err = abs(train_dur - py_dur) if (train_dur is not None and py_dur is not None) else 1e9
        score = (len(rm), len(lm), -dur_err)
        if best is None or score > best[0]:
            best = (score, kind, a0)
    return best[1], best[2]


def fit_clock(py_rew, valve_r, rmatch, py_lick, lick_r, lmatch, a_anchor):
    """Final clock: saleae = a + b * elapsed.
    Slope from licks (many points) or rewards; intercept from rewards (exact valve edge)
    or the camera anchor."""
    b, b_src = 1.0, "not measured (assumed 1)"
    lp = [(py_lick[i], lick_r[j]) for j, i in lmatch.items()]
    rp = [(py_rew[i], valve_r[j]) for i, j in rmatch.items()]
    if len(lp) >= 50 and max(p[0] for p in lp) - min(p[0] for p in lp) > 60:
        e = np.array([p[0] for p in lp]); s = np.array([p[1] for p in lp])
        bl, _ = np.polyfit(e, s, 1)
        if abs(bl - 1) <= DRIFT_PPM * 1e-6:
            b, b_src = float(bl), f"{len(lp)} licks"
    if b_src.startswith("not") and len(rp) >= 2 and max(p[0] for p in rp) - min(p[0] for p in rp) > 30:
        e = np.array([p[0] for p in rp]); s = np.array([p[1] for p in rp])
        br, _ = np.polyfit(e, s, 1)
        if abs(br - 1) <= DRIFT_PPM * 1e-6:
            b, b_src = float(br), f"{len(rp)} rewards"
    if rp:
        a = float(np.median([s - b * e for e, s in rp]))
        a_src = f"{len(rp)} reward valve edges"
    else:
        a = a_anchor
        a_src = "camera / anchor"
    return a, b, a_src, b_src


# ----------------------------------------------------------------------------- one session
def process_session(s: dict, out_dir: Path, ffmpeg, offsets: dict) -> dict:
    name, day = s["name"], s["day"]
    warn, info = [], list(s.get("found", []))
    status = "OK"
    for n in s.get("notes", []):
        warn.append(f"[WARN] {n}")
        status = "WARN"

    def flag(level, msg):
        nonlocal status
        warn.append(f"[{level}] {msg}")
        if level == "FAIL":
            status = "FAIL"
        elif status == "OK":
            status = "WARN"

    py = read_python_log(s["python"])
    sal = read_saleae_csv(s["saleae_csv"])
    ch = sal["channels"]
    if CHANNEL_NAMES["camera"] not in ch:
        raise ValueError(f"Saleae CSV has no '{CHANNEL_NAMES['camera']}' column (columns: {sal['names']})")
    for k in ("lick", "reward"):
        if CHANNEL_NAMES[k] not in ch:
            flag("WARN", f"Saleae CSV has no '{CHANNEL_NAMES[k]}' column (columns: {sal['names']})")
    for k, v in CHANNEL_NAMES.items():
        if v in ch and ch[v]["starts_high"]:
            flag("WARN", f"Saleae channel '{v}' was already HIGH when the capture started "
                         "(start Saleae before Python); that first partial pulse is ignored")

    cam_all = ch[CHANNEL_NAMES["camera"]]["rises"]
    cam_iv_all = ch[CHANNEL_NAMES["camera"]]["intervals"]
    lick_iv = ch.get(CHANNEL_NAMES["lick"], {"intervals": []})["intervals"]
    valve_iv = ch.get(CHANNEL_NAMES["reward"], {"intervals": []})["intervals"]
    lick_r = np.array([iv[0] for iv in lick_iv])
    valve_r = np.array([iv[0] for iv in valve_iv])
    if len(cam_all) == 0:
        raise ValueError("no camera pulses in the Saleae capture")

    ev = py["events"]
    py_rew_ev = [e for e in ev if e.get("tag") == "reward_cmd"]
    py_lick_ev = [e for e in ev if e.get("tag") == "lick"]
    py_rew = np.array([e["elapsed"] for e in py_rew_ev], dtype=float)
    py_lick = np.array([e["elapsed"] for e in py_lick_ev], dtype=float)
    end_ev = next((e for e in reversed(ev) if e.get("tag") == "experiment_end"), None)
    py_dur = end_ev["elapsed"] if end_ev else (ev[-1]["elapsed"] if ev else None)
    if not end_ev:
        flag("WARN", "python log has no experiment_end (session crashed or log cut short?)")

    # ---- where is python's elapsed=0 on the Saleae clock?
    kind, a0 = choose_anchor(cam_all, valve_r, lick_r, py_rew, py_lick, py_dur)
    n_trains = len(pulse_trains(cam_all))
    if n_trains > 1:
        info.append(f"Saleae capture has {n_trains} separate camera pulse trains; "
                    f"this python log was placed by {kind}")
    rmatch = match_rewards(py_rew, valve_r, a0)
    lmatch = match_licks(a0 + py_lick, lick_r)
    a, b, a_src, b_src = fit_clock(py_rew, valve_r, rmatch, py_lick, lick_r, lmatch, a0)
    # refine once with the final clock
    rmatch = match_rewards(py_rew, valve_r, a)
    lmatch = match_licks(a + b * py_lick, lick_r)
    a, b, a_src, b_src = fit_clock(py_rew, valve_r, rmatch, py_lick, lick_r, lmatch, a)

    def to_sal(e):
        return a + b * e

    info.append(f"python -> saleae: t = {a:.6f} + {b:.9f} * elapsed  (offset from {a_src}; "
                f"drift from {b_src})")
    if b_src.startswith("not") and (py_dur or 0) > 120:
        flag("WARN", "clock drift could not be measured (too few licks/rewards) - "
                     "python-only event times may be off by a few ms late in the session")

    # ---- this session's camera pulses
    t_start = to_sal(0.0)
    t_end = to_sal(py_dur) if py_dur is not None else float(cam_all[-1])
    sel = (cam_all >= t_start - 0.05) & (cam_all <= t_end + 0.1)
    cam = cam_all[sel]
    cam_iv = [iv for iv, keep in zip(cam_iv_all, sel) if keep]
    if len(cam) == 0:
        raise ValueError("no camera pulses inside the python session - wrong file pair, "
                         "or the Saleae capture did not cover the session")
    period = float(np.median(np.diff(cam))) if len(cam) > 1 else 0.033
    if len(cam_all) != len(cam):
        info.append(f"{len(cam_all) - len(cam)} camera pulses outside this python session were ignored")
    if abs(cam[0] - t_start) > 0.05:
        flag("FAIL", f"first camera pulse of the session is {1000*(cam[0] - t_start):.0f} ms from python start "
                     "- Saleae capture probably started after Python (camera start not recorded)")
    if py_dur is not None and t_end - cam[-1] > 3 * period:
        flag("WARN", f"camera pulses stop {t_end - cam[-1]:.2f} s before python ended - "
                     "Saleae capture stopped too early?")
    if len(cam) > 1:
        d = np.diff(cam) * 1000
        info.append(f"camera pulses {len(cam)}: mean {1000/d.mean():.3f} Hz, interval {d.min():.2f}-{d.max():.2f} ms")
        big = int((d > 1.5 * np.median(d)).sum())
        if big:
            flag("WARN", f"{big} gap(s) in camera pulses longer than 1.5x the normal interval")

    # ---- reward / lick quality
    in_sess = lambda t: t_start - 0.05 <= t <= t_end + 0.1   # noqa: E731
    valve_in = [j for j, t in enumerate(valve_r) if in_sess(t)]
    resid = [(valve_r[j] - to_sal(py_rew[i])) * 1000 for i, j in rmatch.items()]
    if len(py_rew_ev) or len(valve_in):
        msg = (f"rewards: python {len(py_rew_ev)}, saleae {len(valve_in)} in session, matched {len(rmatch)}")
        if resid:
            msg += f", residual max {max(abs(r) for r in resid):.2f} ms"
        if resid and max(abs(r) for r in resid) > FAIL_RESID_MS:
            flag("FAIL", msg + " - rewards do not line up: probably the wrong Saleae/Python pair")
        elif len(py_rew_ev) >= 3 and len(rmatch) < 0.8 * len(py_rew_ev):
            if len(valve_in) == 0:
                flag("WARN", msg + " - no valve pulses on the Saleae 'reward' channel (wiring?)")
            else:
                flag("FAIL", msg + " - most rewards do not line up: probably the wrong Saleae/Python pair")
        elif len(rmatch) == len(py_rew_ev) == len(valve_in):
            info.append(msg)
        else:
            flag("WARN", msg)
    outside_v = len(valve_r) - len(valve_in)
    if outside_v:
        info.append(f"{outside_v} valve pulse(s) outside the python session (e.g. flushing) - listed, not matched")

    lick_in = [j for j, t in enumerate(lick_r) if in_sess(t)]
    msg = f"licks: saleae {len(lick_in)} in session, python {len(py_lick_ev)}, matched {len(lmatch)}"
    delays = [(to_sal(py_lick[i]) - lick_r[j]) * 1000 for j, i in lmatch.items()]
    med_delay = float(np.median(delays)) / 1000 if delays else LICK_DELAY0
    if delays:
        msg += (f", python late by {np.median(delays):.1f} ms "
                f"(range {min(delays):.1f}-{max(delays):.1f})")
    n_sal_only = len(lick_in) - len(lmatch)
    n_py_only = len(py_lick_ev) - len(lmatch)
    # a few python-only lockout licks while the wet spout keeps the line high are normal
    if n_sal_only <= max(2, 0.005 * len(lick_in)) and n_py_only <= max(5, 0.05 * len(py_lick_ev)):
        info.append(msg + (f" ({n_py_only} python-only licks, normal while the spout is wet)" if n_py_only else ""))
    else:
        flag("WARN", msg)
    if delays and (max(delays) - min(delays)) > 40:
        flag("WARN", "python lick delay varies by more than 40 ms - check the lick matching")
    outside_l = len(lick_r) - len(lick_in)
    if outside_l:
        info.append(f"{outside_l} lick(s) before python start / after python end - listed, not matched")
    if end_ev:
        info.append(f"python experiment_end at {t_end:.3f} s, last camera pulse {cam[-1]:.3f} s")

    # ---- video
    vid = probe_video(s["video"], ffmpeg) if s.get("video") else {"frames": None, "error": "no video file",
                                                                    "created": None, "fps_tag": None}
    n_frames = vid.get("frames")
    pstart = py["start_wall"]
    offset = 0
    if vid.get("error"):
        flag("WARN", f"video: {vid['error']} - frame numbers are camera pulse numbers")
    else:
        deficit = len(cam) - n_frames
        key = (day, name)
        if key in offsets:
            offset = offsets[key]
            flag("WARN", f"video offset set by hand in {OFFSETS_FILE}: video frame 0 = camera pulse {offset} "
                         f"(video {n_frames} frames, pulses {len(cam)})")
        elif deficit == 0:
            info.append(f"video frames {n_frames} = camera pulses {len(cam)}  -> video frame k = pulse k")
        elif deficit > 0:
            # mp4 creation_time is truncated to whole seconds: a file opened L s after python start
            # shows L-1 .. L s; one opened at/before python start shows about -1 .. +0.1 s.
            dv = (vid["created"] - pstart).total_seconds() if (pstart and vid.get("created")) else None
            late_s = deficit * period
            if dv is not None and dv > 0.2 and late_s - 1.2 <= dv <= late_s + 0.5:
                offset = deficit
                flag("WARN", f"video has {deficit} fewer frames than camera pulses and the file was created "
                             f"{dv:.1f} s after python start (~{late_s:.1f} s expected) -> recording STARTED LATE: "
                             f"video frame 0 = pulse {offset}. To override, edit {OFFSETS_FILE}")
            elif dv is not None and dv <= 0.2:
                flag("WARN", f"video has {deficit} fewer frames than camera pulses; the file was created at "
                             f"python start -> recording STOPPED EARLY: video frame 0 = pulse 0 "
                             f"(the last {deficit} pulses have no video). To override, edit {OFFSETS_FILE}")
            else:
                why = "the video has no creation time" if dv is None else \
                      f"the file creation time ({dv:+.1f} s from python start) does not fit either case"
                flag("WARN", f"video has {deficit} fewer frames than camera pulses and {why} -> ASSUMED the "
                             f"recording stopped early (video frame 0 = pulse 0). If it started late instead, "
                             f"set first_pulse_in_video={deficit} in {OFFSETS_FILE}")
        else:
            flag("WARN", f"video has {-deficit} MORE frames than camera pulses (camera not in trigger mode, or "
                         "the recording includes another run?) - frame numbers assume video frame 0 = pulse 0")
        if vid.get("fps_tag"):
            info.append(f"video file says {vid['fps_tag']:g} fps - ignore it, real frame times come from Saleae")
    if pstart and vid.get("created"):
        dv = (vid["created"] - pstart).total_seconds()
        if abs(dv) <= VIDEO_START_TOL_S:
            info.append(f"video file created {dv:+.1f} s from python start")
        else:
            flag("WARN", f"video file created {dv:+.0f} s from python start - wrong file pair?")
    sal_start = read_sal_start(s["sal"]) if s.get("sal") else None
    if pstart and sal_start:
        rel = (pstart - sal_start).total_seconds()
        if -1.0 <= rel <= sal["end"] + 1.0:
            info.append(f"python started {rel:.1f} s into the Saleae capture (pulse-based: {t_start:.1f} s)")
        else:
            flag("WARN", f"python start is {rel:.1f} s from the Saleae capture start, outside the capture - "
                         "wrong file pair?")
    elif not s.get("sal"):
        info.append("no .sal file next to the CSV - skipped the Saleae start-time check")

    have_video = n_frames is not None and not vid.get("error")

    def pulse_of(t):
        if t is None or np.isnan(t):
            return None
        k = bisect_right(cam, t) - 1
        if k < 0 or t - cam[k] > 1.5 * period:
            return None
        return k

    def frame_of(t):
        k = pulse_of(t)
        if k is None:
            return None
        f = k - offset
        if f < 0 or (have_video and f >= n_frames):
            return None
        return f

    def in_valve(on, off):
        for v_on, v_off in valve_iv:
            v_off = v_off if not np.isnan(v_off) else v_on
            if on <= v_off + VALVE_LICK_PAD_S and (np.isnan(off) or off >= v_on):
                return True
        return False

    def where(t):
        if t < t_start - 0.05:
            return "before python start"
        if t > t_end + 0.1:
            return "after python end"
        return ""

    rows = []
    lick_by_rise = {j: py_lick_ev[i] for j, i in lmatch.items()}
    for j, (on, off) in enumerate(lick_iv):
        pe = lick_by_rise.get(j)
        note = where(on) or ("" if pe else "no matching python lick")
        rows.append({"event": "lick", "source": "saleae", "trial": pe.get("trial") if pe else None,
                     "valid": pe.get("valid") if pe else None, "t_on": on, "t_off": off,
                     "during_valve": in_valve(on, off),
                     "py_elapsed": pe["elapsed"] if pe else None,
                     "py_delay_ms": (to_sal(pe["elapsed"]) - on) * 1000 if pe else None, "note": note})
    matched_py = set(lmatch.values())
    for i, e in enumerate(py_lick_ev):
        if i not in matched_py:
            rows.append({"event": "lick", "source": "python", "trial": e.get("trial"), "valid": e.get("valid"),
                         "t_on": to_sal(e["elapsed"]) - med_delay, "t_off": None, "during_valve": None,
                         "py_elapsed": e["elapsed"], "py_delay_ms": None,
                         "note": f"python lick with no Saleae edge (python time - {1000*med_delay:.0f} ms)"})
    rew_by_valve = {j: py_rew_ev[i] for i, j in rmatch.items()}
    for j, (on, off) in enumerate(valve_iv):
        pe = rew_by_valve.get(j)
        note = where(on) or ("" if pe else "no matching python reward_cmd")
        rows.append({"event": "reward_valve", "source": "saleae", "trial": pe.get("trial") if pe else None,
                     "valid": None, "t_on": on, "t_off": off, "during_valve": None,
                     "py_elapsed": pe["elapsed"] if pe else None,
                     "py_delay_ms": (to_sal(pe["elapsed"]) - on) * 1000 if pe else None, "note": note})
    for i, e in enumerate(py_rew_ev):
        if i not in rmatch:
            rows.append({"event": "reward_cmd", "source": "python", "trial": e.get("trial"), "valid": None,
                         "t_on": to_sal(e["elapsed"]), "t_off": None, "during_valve": None,
                         "py_elapsed": e["elapsed"], "py_delay_ms": None,
                         "note": "python reward with no Saleae valve edge"})
    for e in ev:
        tag = e.get("tag")
        if tag in PY_ONLY_SKIP:
            continue
        if tag == "trial_end" and "requirement" in e:   # duplicate summary line written by the experiment loop
            continue
        rows.append({"event": tag, "source": "python", "trial": e.get("trial"),
                     "valid": e.get("hit") if tag == "trial_end" else None,
                     "t_on": to_sal(e["elapsed"]), "t_off": None, "during_valve": None,
                     "py_elapsed": e["elapsed"], "py_delay_ms": None,
                     "req": e.get("req"), "note": "time from python log (~1 ms)"})
    rows.sort(key=lambda r: (r["t_on"], r["event"]))

    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = float(cam[0])
    with open(out_dir / f"{name}_events.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["event", "source", "trial", "valid_or_hit", "t_saleae_s", "t_from_first_pulse_s",
                    "frame", "t_off_saleae_s", "frame_off", "duration_ms", "during_valve",
                    "python_elapsed_s", "python_delay_ms", "note", "req"])
        for r in rows:
            off = r["t_off"]
            has_off = off is not None and not np.isnan(off)
            w.writerow([r["event"], r["source"], fmt(r["trial"]), fmt(r["valid"]), fmt(r["t_on"]),
                        fmt(r["t_on"] - t0), fmt(frame_of(r["t_on"])), fmt(off if has_off else None),
                        fmt(frame_of(off) if has_off else None),
                        fmt((off - r["t_on"]) * 1000 if has_off else None, 2),
                        fmt(r["during_valve"]), fmt(r["py_elapsed"], 3), fmt(r["py_delay_ms"], 2), r["note"],
                        fmt(r.get("req"), 4)])

    with open(out_dir / f"{name}_frames.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pulse", "frame", "t_saleae_s", "t_from_first_pulse_s", "interval_ms", "pulse_high_ms"])
        for k, (on, off) in enumerate(cam_iv):
            fr = k - offset
            fr_ok = fr >= 0 and (not have_video or fr < n_frames)
            w.writerow([k, fmt(fr if fr_ok else None), fmt(on), fmt(on - t0),
                        fmt((on - cam_iv[k - 1][0]) * 1000 if k else None, 3),
                        fmt((off - on) * 1000 if not np.isnan(off) else None, 3)])

    write_report(out_dir, s, status, warn, info)
    return {"day": day, "session": name, "status": status,
            "video_frames": n_frames, "camera_pulses": len(cam), "video_offset": offset,
            "camera_hz": round(1 / np.mean(np.diff(cam)), 3) if len(cam) > 1 else None,
            "licks_saleae": len(lick_in), "licks_python": len(py_lick_ev), "licks_matched": len(lmatch),
            "rewards_saleae": len(valve_in), "rewards_python": len(py_rew_ev), "rewards_matched": len(rmatch),
            "duration_s": round(float(cam[-1] - cam[0]), 1), "problems": " | ".join(warn)}


def write_report(out_dir: Path, s: dict, status: str, warn: list, info: list):
    out_dir.mkdir(parents=True, exist_ok=True)
    day, name = s["day"], s["name"]
    with open(out_dir / f"{name}_report.txt", "w", encoding="utf-8") as f:
        f.write(f"SYNC REPORT  {day + '/' if day else ''}{name}   status: {status}   (tool {TOOL_VERSION})\n")
        f.write(f"python : {s['python']}\nsaleae : {s.get('saleae_csv') or '-'}\n")
        f.write(f"sal    : {s.get('sal') or '-'}\nvideo  : {s.get('video') or '-'}\n\n")
        for line in warn:
            f.write(line + "\n")
        for line in info:
            f.write("[ok] " + line + "\n")
        f.write("\nNotes: video frame numbers start at 0. Frame time = rising edge of the camera pulse; if the "
                "camera is set to trigger on the falling edge add ~16 ms. Lick on/off times that overlap a valve "
                "opening (during_valve=1) are not reliable contact times.\n")


# ----------------------------------------------------------------------------- discovery
def find_one(folder: Path, name: str, exts, prefer=None):
    for ext in exts:
        p = folder / f"{name}{ext}"
        if p.is_file():
            return p, None
    sub = folder / name
    if sub.is_dir():
        hits = [p for ext in exts for p in sorted(sub.glob(f"*{ext}"))]
        if prefer:
            pref = [p for p in hits if p.name.lower() == prefer]
            if pref:
                return pref[0], None
        hits = [p for p in hits if not p.name.lower().startswith("analog")]
        if len(hits) == 1:
            return hits[0], None
        if len(hits) > 1:
            return None, f"{len(hits)} candidate files in {sub.name}/ - name it {name}{exts[0]}"
    return None, None


def discover(root: Path):
    batches = []
    if (root / "python").is_dir():          # project layout: <project>/python, salae, video -> sync
        sal_dir = next((root / s for s in ("salae", "saleae") if (root / s).is_dir()), None)
        batches.append({"prefix": root.name, "python": root / "python", "saleae": sal_dir,
                        "video": root / "video" if (root / "video").is_dir() else None, "out": root / "sync"})
    for py_dir in sorted(root.glob("*_python")):
        if not py_dir.is_dir():
            continue
        prefix = py_dir.name[: -len("_python")]
        sal_dir = next((root / f"{prefix}{s}" for s in ("_salae", "_saleae") if (root / f"{prefix}{s}").is_dir()), None)
        vid_dir = root / f"{prefix}_video"
        batches.append({"prefix": prefix, "python": py_dir, "saleae": sal_dir,
                        "video": vid_dir if vid_dir.is_dir() else None, "out": root / f"{prefix}_sync"})
    return batches


def locate(batch_dir, rel_dir, name, exts, prefer=None):
    """Find <name>.<ext>: first in the mirrored folder, then anywhere under the batch folder
    (only if exactly one file has that name). Returns (path, problem_note, found_note)."""
    if batch_dir is None:
        return None, None, None
    p, note = find_one(batch_dir / rel_dir, name, exts, prefer)
    if p or note:
        return p, note, None
    wanted = {f"{name}{e}".lower() for e in exts}
    hits = [q for q in batch_dir.rglob("*") if q.is_file() and q.name.lower() in wanted]
    if len(hits) == 1:
        return hits[0], None, f"{hits[0].name} found in {hits[0].parent.relative_to(batch_dir).as_posix()}/ " \
                              f"(not in the matching folder {rel_dir.as_posix()}/)"
    if len(hits) > 1:
        return None, f"{len(hits)} files named {name}{exts[0]} in {batch_dir.name} - keep one, or set it in " \
                     f"{MANIFEST_FILE}", None
    return None, None, None


MANIFEST_FILE = "session_files.csv"   # optional, in <batch>_sync: fix pairings without renaming data
MANIFEST_COLS = ("python", "saleae_csv", "sal", "video")


def read_manifest(batch, root: Path) -> dict:
    """<batch>_sync/session_files.csv: day,session,python,saleae_csv,sal,video,skip.
    Filled cells override what the tool found; paths can be absolute or relative to the root folder."""
    p = batch["out"] / MANIFEST_FILE
    res = {}
    if not p.exists():
        return res
    try:
        with open(p, newline="", encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                day = (r.get("day") or "").strip().replace("\\", "/").strip("/")
                ses = (r.get("session") or "").strip()
                if not ses:
                    continue
                o = {"skip": (r.get("skip") or "").strip().lower() in ("1", "yes", "y", "true", "x")}
                for c in MANIFEST_COLS:
                    v = (r.get(c) or "").strip().strip('"')
                    if v:
                        q = Path(v)
                        o[c] = q if q.is_absolute() else (root / q)
                res[(day, ses)] = o
    except OSError as exc:
        print(f"  (could not read {p.name}: {exc})")
    return res


def sessions_in(batch, root: Path):
    manifest = read_manifest(batch, root)
    out = []
    seen = set()
    for log in sorted(batch["python"].rglob("*.jsonl")):
        rel_dir = log.parent.relative_to(batch["python"])
        day = "" if str(rel_dir) == "." else rel_dir.as_posix()
        out.append(_session(batch, rel_dir, day, log.stem, log, manifest.get((day, log.stem), {})))
        seen.add((day, log.stem))
    for (day, ses), o in manifest.items():          # sessions given only in the manifest
        if (day, ses) not in seen and o.get("python"):
            out.append(_session(batch, Path(day) if day else Path("."), day, ses, o["python"], o))
    return out


def _session(batch, rel_dir, day, name, log, o):
    notes, found = [], []
    sal_p, n1, f1 = locate(batch["saleae"], rel_dir, name, (".sal",))
    csv_p, n2, f2 = locate(batch["saleae"], rel_dir, name, (".csv",), prefer="digital.csv")
    if csv_p is None and n2 is None and sal_p is not None:
        # Logic 2 export saved as digital.csv next to <name>.sal
        dig = sal_p.parent / "digital.csv"
        if dig.is_file() and len(list(sal_p.parent.glob("*.sal"))) == 1:
            csv_p = dig
            found.append(f"Saleae CSV = digital.csv next to {sal_p.name}")
    vid_p, n3, f3 = locate(batch["video"], rel_dir, name, VIDEO_EXTS)
    s = {"day": day, "name": name, "python": Path(log), "saleae_csv": csv_p, "sal": sal_p, "video": vid_p,
         "skip": bool(o.get("skip"))}
    for c in MANIFEST_COLS:
        if o.get(c):
            s[c] = o[c]
            found.append(f"{c} set in {MANIFEST_FILE}: {o[c]}")
    for k in ("python", "saleae_csv", "sal", "video"):
        if s.get(k) is not None and not Path(s[k]).is_file():
            notes.append(f"{k} file does not exist: {s[k]}")
            s[k] = None if k != "python" else s[k]
    s["notes"] = [n for n in (n1, n2, n3) if n] + notes
    s["found"] = [f for f in (f1, f2, f3) if f] + found
    return s


def _tool_hash() -> str:
    import hashlib
    try:
        return hashlib.md5(Path(__file__).read_bytes()).hexdigest()
    except OSError:
        return TOOL_VERSION


def fingerprint(s) -> dict:
    fp = {"tool": TOOL_VERSION, "tool_hash": _tool_hash()}
    for k in ("python", "saleae_csv", "sal", "video"):
        p = s.get(k)
        if p:
            st = p.stat()
            fp[k] = [str(p), st.st_size, int(st.st_mtime)]
        else:
            fp[k] = None
    return fp


def up_to_date(s, out_dir: Path) -> bool:
    fpf = out_dir / f"{s['name']}_inputs.json"
    try:
        return json.loads(fpf.read_text(encoding="utf-8")) == fingerprint(s)
    except Exception:
        return False


# ----------------------------------------------------------------------------- main
SUMMARY_COLS = ["day", "session", "status", "video_frames", "camera_pulses", "video_offset", "camera_hz",
                "licks_saleae", "licks_python", "licks_matched", "rewards_saleae", "rewards_python",
                "rewards_matched", "duration_s", "problems"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent,
                    help="folder that contains the <batch>_python/_salae/_video folders "
                         "(default: the folder this tool sits in)")
    ap.add_argument("--force", action="store_true", help="redo sessions that are already done")
    args = ap.parse_args(argv)

    ffmpeg = find_ffmpeg()
    batches = discover(args.root)
    if not batches:
        print(f"No '<batch>_python' folders found in {args.root}")
        return 1
    rc = 0
    for batch in batches:
        print(f"\n=== {batch['prefix']} ===")
        if not batch["saleae"]:
            print(f"  no {batch['prefix']}_salae folder - skipped")
            continue
        summary_path = batch["out"] / "sync_summary.csv"
        old = {}
        if summary_path.exists():
            try:
                with open(summary_path, newline="", encoding="utf-8-sig") as f:
                    for r in csv.DictReader(f):
                        if r.get("session") is not None:
                            old[(r.get("day", ""), r["session"])] = r
            except Exception as exc:
                print(f"  (could not read old {summary_path.name}: {exc} - it will be rebuilt)")
        offsets = read_offsets(batch["out"])
        results = {}
        for s in sessions_in(batch, args.root):
            label = f"{s['day'] + '/' if s['day'] else ''}{s['name']}"
            key = (s["day"], s["name"])
            out_dir = batch["out"] / s["day"] if s["day"] else batch["out"]
            if s.get("skip"):
                print(f"  {label:45s} SKIPPED - marked skip in {MANIFEST_FILE}")
                continue
            if not s["saleae_csv"]:
                why = "; ".join(s["notes"]) or f"no Saleae CSV named {s['name']}.csv in {batch['saleae'].name}/{s['day']}"
                print(f"  {label:45s} SKIPPED - {why}")
                continue
            if not args.force and key in old and up_to_date(s, out_dir) and \
                    not (key in offsets and str(offsets[key]) != str(old[key].get("video_offset"))):
                results[key] = old[key]
                print(f"  {label:45s} already done ({old[key].get('status')})")
                continue
            try:
                r = process_session(s, out_dir, ffmpeg, offsets)
            except Exception as exc:
                r = {"day": s["day"], "session": s["name"], "status": "FAIL", "problems": f"error: {exc}"}
                for suffix in ("_events.csv", "_frames.csv"):
                    try:
                        (out_dir / f"{s['name']}{suffix}").unlink()
                    except OSError:
                        pass
                write_report(out_dir, s, "FAIL", [f"[FAIL] error: {exc}"], [])
            try:
                (out_dir / f"{s['name']}_inputs.json").write_text(json.dumps(fingerprint(s)), encoding="utf-8")
            except OSError:
                pass
            results[key] = r
            if r["status"] == "FAIL":
                rc = max(rc, 1)
            print(f"  {label:45s} {r['status']}" + (f"  - {r['problems']}" if r.get("problems") else ""))
        if results:
            batch["out"].mkdir(parents=True, exist_ok=True)
            try:
                with open(summary_path, "w", newline="", encoding="utf-8-sig") as f:
                    w = csv.DictWriter(f, fieldnames=SUMMARY_COLS, extrasaction="ignore")
                    w.writeheader()
                    for key in sorted(results):
                        w.writerow({c: results[key].get(c, "") for c in SUMMARY_COLS})
                print(f"  -> results in {batch['out']}")
            except PermissionError:
                print(f"  !! could not write {summary_path} - close it in Excel and run again")
                rc = 2
    return rc


if __name__ == "__main__":
    sys.exit(main())
