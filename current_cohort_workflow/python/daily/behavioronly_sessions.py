#!/usr/bin/env python3
"""
behavioronly_sessions.py - a Step 02 table for a session WITHOUT video, as an OPTION for listed sessions (2026-11-16).

Why: the original Step 02 reads one concat_out per session (rows = video frames with the pupil, the camera time and
the lick / valve states). A session whose video was not saved has no concat_out, so Step 02 drops the WHOLE session,
and its behavior (licks, rewards, PR) is missing from the longitudinal table although the Saleae file and the task
log are fine. The user chose on 2026-11-16 to keep the behavior and leave only the pupil absent (day 3 9670f black:
"행동만 넣기").

What, for every session approved in <pupil_root>/behavioronly_sessions.csv (day,cage,mouse,approved,reason):
  1) a pupil file with one row per Saleae camera pulse of the session (sync_sessions.py's <session>_frames.csv) and
     EMPTY pupil columns (Timestamp_seconds = row / 60 as in the package; diameter and centre NaN);
  2) the original alignment code (current_cohort_workflow/matlab/pupil_align_ttl_sync.m = the package's
     pupil_align_ttl with the frame <-> pulse match from frames.csv) builds combined_pupil_digital.csv,
     lick_onsets.csv, injector_onsets.csv and summary.txt from the session's digital.csv, so the lick / valve columns
     are made exactly as for every other session;
  3) checks: rows = camera pulses; pupil columns all empty; camera times = frames.csv; Camera_TTL = 1 at every row;
     lick and valve states at every pulse = the Saleae intervals of <session>_events.csv (at most 3 rows off); reward
     onsets = the Saleae valve openings during the camera pulses;
  4) only then the folder is renamed concat_out_<stamp>_behavioronly and fix_source.json is written. concat_choice.py
     gives it to prepare_step02_base.py while the session stays approved. Log: <pupil_root>/behavioronly_log.csv.
A session that has a package table (a video) is never touched. Nothing of the package or of Step 02 is changed.

usage: python behavioronly_sessions.py --project <project> [--dry-run]
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import concat_choice as CC  # noqa: E402
import fix_concat_alignment as FCA  # noqa: E402  (MATLAB helpers, state_at, onsets, fingerprint)

TMP_PREFIX = "_behavioronly_tmp_"
EMPTY_DIR = "pupil_outputs_behavioronly"
PUPIL_EMPTY = ["Diameter_px", "CenterX_px", "CenterY_px"]
LOG_COLS = ["when", "day", "cage", "mouse", "status", "run", "rows", "lick_state_mismatch", "valve_state_mismatch",
            "lick_onsets", "reward_onsets", "saleae_valve_openings", "note"]


def find_dir(base: Path, day, cage, mouse):
    d = Path(base) / f"day{int(day)}"
    if not d.is_dir():
        return None
    c = next((p for p in d.iterdir() if p.is_dir() and p.name.lower() == cage), None)
    return next((p for p in c.iterdir() if p.is_dir() and p.name.lower() == mouse), None) if c else None


def digital_csv(project: Path, sync_dir: Path, salae_dir):
    """the digital.csv sync_sessions.py used (its inputs.json; a pre-move path is re-rooted at the project)"""
    for j in sorted(sync_dir.glob("*_inputs.json")):
        try:
            p = Path(json.loads(j.read_text(encoding="utf-8"))["saleae_csv"][0])
        except (OSError, ValueError, KeyError, TypeError, IndexError):
            continue
        if not p.exists() and "salae" in p.parts:
            p = project.joinpath(*p.parts[p.parts.index("salae"):])
        if p.exists():
            return p
    hits = sorted(Path(salae_dir).rglob("digital.csv")) if salae_dir else []
    return hits[0] if len(hits) == 1 else None


def write_empty_pupil(path: Path, n: int):
    k = np.arange(n)
    pd.DataFrame({"Timestamp_seconds": k / 60.0, "Diameter_pixels": np.nan, "CenterX_px": np.nan,
                  "CenterY_px": np.nan, "UsedScoreFallback": 0, "SmoothedDiameter_px": np.nan}) \
        .to_csv(path, index=False, float_format="%.7f", na_rep="NaN")


def verify(run: Path, frames_csv: Path, events_csv, n: int) -> dict:
    N = pd.read_csv(run / "combined_pupil_digital.csv")
    want = set(FCA.PUPIL_COLS) | set(FCA.DIGITAL_COLS)
    if set(N.columns) != want:
        raise ValueError(f"columns {list(N.columns)} are not the package's columns")
    if len(N) != n:
        raise ValueError(f"{len(N)} rows, the session has {n} camera pulses")
    if N[PUPIL_EMPTY].notna().any().any():
        raise ValueError("pupil columns are not empty")
    F = pd.read_csv(frames_csv, usecols=["frame", "t_saleae_s"]).dropna()
    F["frame"] = F["frame"].astype(int)
    m = N.assign(frame=N["Frame"] - 1).merge(F, on="frame", how="left")
    if m["t_saleae_s"].isna().any():
        raise ValueError(f"{int(m['t_saleae_s'].isna().sum())} rows without a camera pulse")
    dt = float(np.max(np.abs(m["CamTime_s"] - m["t_saleae_s"])))
    if not np.isfinite(dt) or dt > 1e-6:
        raise ValueError(f"camera times differ from frames.csv by up to {dt:.6f} s")
    if not (m["Camera_TTL"] == 1).all():
        raise ValueError("Camera_TTL is not 1 at every camera pulse")
    info = dict(rows=len(N), max_time_diff_s=dt, lick_state_mismatch=np.nan, valve_state_mismatch=np.nan,
                lick_onsets=FCA.onsets(N["Lick_TTL"]), reward_onsets=FCA.onsets(N["Injector_TTL"]),
                saleae_valve_openings=np.nan)
    if events_csv is not None:
        E = pd.read_csv(events_csv, usecols=["event", "source", "t_saleae_s", "t_off_saleae_s", "frame"])
        E = E[E["source"] == "saleae"]
        t = m["CamTime_s"].to_numpy(float)
        for ev, col, key in (("lick", "Lick_TTL", "lick_state_mismatch"),
                             ("reward_valve", "Injector_TTL", "valve_state_mismatch")):
            iv = E.loc[E["event"] == ev, ["t_saleae_s", "t_off_saleae_s"]].to_numpy(float)
            got = m[col].to_numpy(float) > 0.5
            info[key] = int(np.sum(FCA.state_at(t, [tuple(r) for r in iv]) != got))
            if info[key] > FCA.MAX_STATE_MISMATCH:
                raise ValueError(f"{col} differs from the Saleae {ev} intervals at {info[key]} camera pulses")
        info["saleae_valve_openings"] = int(E.loc[E["event"] == "reward_valve", "frame"].notna().sum())
        if info["reward_onsets"] != info["saleae_valve_openings"]:
            raise ValueError(f"{info['reward_onsets']} reward onsets, but {info['saleae_valve_openings']} Saleae "
                             "valve openings during the camera pulses")
    return info


def log_row(logp: Path, row: dict):
    new = not logp.exists()
    with open(logp, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=LOG_COLS, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow(row)


def one(project: Path, root: Path, key, row: dict, dry: bool) -> str:
    day, cage, mouse = key
    sync_dir = find_dir(project / "sync", day, cage, mouse)
    salae_dir = find_dir(project / "salae", day, cage, mouse)
    if sync_dir is None:
        raise RuntimeError("no sync folder (run sync_sessions.py first)")
    sd = CC.session_dir(root, day, cage, mouse) or (root / f"day{day}" / sync_dir.parent.name / sync_dir.name)
    if sd.exists() and CC.package_runs(sd):
        return "has a package table (a video) -> not touched"
    fcs = sorted(sync_dir.glob("*_frames.csv"))
    evs = sorted(sync_dir.glob("*_events.csv"))
    if len(fcs) != 1:
        raise RuntimeError(f"{len(fcs)} frames.csv in the sync folder")
    dig = digital_csv(project, sync_dir, salae_dir)
    if dig is None:
        raise RuntimeError("the session's digital.csv was not found")
    F = pd.read_csv(fcs[0], usecols=["frame"]).dropna()
    n = int(F["frame"].max()) + 1
    inputs = {"digital": FCA.fingerprint(dig), "frames": FCA.fingerprint(fcs[0])}
    have = sorted(p for p in sd.glob(f"concat_out_*{CC.NOVID}") if (p / "fix_source.json").exists()) \
        if sd.exists() else []
    if have:
        old = json.loads((have[-1] / "fix_source.json").read_text(encoding="utf-8")).get("inputs", {})
        if all(old.get(k, {}).get("sha1") == v["sha1"] for k, v in inputs.items()):
            return f"up to date ({have[-1].name})"
    if dry:
        return f"would write a behavior-only table ({n} camera pulses) from {dig.name}"
    mexe, mdir = FCA.matlab_exe(), FCA.matlab_dir()
    if not mexe or not mdir:
        raise RuntimeError("MATLAB or pupil_align_ttl_sync.m not found")
    (sd / EMPTY_DIR).mkdir(parents=True, exist_ok=True)
    pup = sd / EMPTY_DIR / "empty_pupil_rows.csv"
    write_empty_pupil(pup, n)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    tmp, final = sd / f"{TMP_PREFIX}{stamp}", sd / f"concat_out_{stamp}{CC.NOVID}"
    q = lambda s: str(s).replace("'", "''")  # noqa: E731
    cmd = (f"addpath('{q(mdir)}'); pupil_align_ttl_sync('{q(dig)}', '{q(pup)}', '{q(sd)}', '{q(fcs[0])}', "
           f"'{q(tmp.name)}');")
    try:
        res = subprocess.run([mexe, "-batch", cmd], capture_output=True, text=True, errors="replace",
                             timeout=FCA.MATLAB_TIMEOUT_S)
        if res.returncode != 0 or not (tmp / "combined_pupil_digital.csv").exists():
            raise ValueError(f"MATLAB exit {res.returncode}: " + FCA.clean((res.stderr or res.stdout)[-600:]))
        info = verify(tmp, fcs[0], evs[0] if len(evs) == 1 else None, n)
        tmp.rename(final)
        for fn in ("summary.txt", "log.txt"):
            fp = final / fn
            if fp.exists():
                fp.write_bytes(fp.read_bytes().replace(tmp.name.encode(), final.name.encode()))
        (final / "fix_source.json").write_text(json.dumps(
            {"made_by": "behavioronly_sessions.py + pupil_align_ttl_sync.m", "inputs": inputs,
             "pupil_rows_file": str(pup), "approved": row.get("approved", ""), "reason": row.get("reason", ""),
             "note": "no video: rows = the Saleae camera pulses of the session; pupil columns empty (NaN)",
             "checks": info}, indent=1, default=float), encoding="utf-8")
    except Exception:
        if tmp.exists():
            try:
                tmp.rename(sd / f"FAILED_behavioronly_{stamp}")
            except OSError:
                pass
        raise
    row.update(_run=final.name, **info)
    return (f"wrote {final.name} | rows {info['rows']} | lick / valve state mismatches "
            f"{info['lick_state_mismatch']} / {info['valve_state_mismatch']} | lick onsets {info['lick_onsets']}, "
            f"reward onsets {info['reward_onsets']} (Saleae valve openings {info['saleae_valve_openings']})")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    project = Path(a.project).resolve()
    root = CC.active_pupil_root(project)      # 2026-11-17: the tree the tables use (pupil_active_root.txt)
    rows = {CC.norm_key(r["day"], r["cage"], r["mouse"]): r
            for r in CC.read_rows(root / CC.LISTS[CC.NOVID]) if CC.is_approved(r.get("approved"))}
    if not rows:
        print("  (no session approved in behavioronly_sessions.csv)")
        return 0
    bad = 0
    for key in sorted(rows):
        r = dict(rows[key])
        try:
            msg = one(project, root, key, r, a.dry_run)
            status = "ok"
        except Exception as exc:  # noqa: BLE001
            msg, status, bad = f"FAILED ({FCA.clean(str(exc))})", "FAILED", bad + 1
        print(f"  behavior only day{key[0]} {key[1]}_{key[2]}: {msg}")
        if not a.dry_run and not msg.startswith("up to date"):
            log_row(root / "behavioronly_log.csv",
                    dict(when=datetime.now().isoformat(timespec="seconds"), day=key[0], cage=key[1], mouse=key[2],
                         status=status if status != "ok" else msg.split(" |")[0], run=r.get("_run", ""),
                         note=r.get("reason", ""), **{k: r.get(k, "") for k in LOG_COLS[6:12]}))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
