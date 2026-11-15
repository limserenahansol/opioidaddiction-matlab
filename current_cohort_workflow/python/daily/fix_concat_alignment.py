"""Sync-aligned copy of the package's concat_out (the Step 02 input) for sessions the user approved.

Why: when a video has more or fewer frames than camera TTL edges, the package's pupil_align_ttl (the original Step 01a
logic, not changed) matches them at the END. On the current rig the video starts with the first camera pulse (the
recording stopped early, or extra frames were written at the end), so the package's concat_out then has licks and
rewards at the wrong video frames. pupil_session_quality.py lists such sessions (align_check MISMATCH).

What: for every session in <pupil_out_root>/alignment_fixes.csv (one row per session the user approved: day, cage,
mouse, approved, reason; a row with an empty 'approved' cell, or one starting with 'no', is skipped) this runs
current_cohort_workflow/matlab/pupil_align_ttl_sync.m - the original function with one change: video frame k is
matched to the camera pulse that sync_sessions.py gives for frame k (<session>_frames.csv). Nothing of the package's
output is changed or removed.
  1) MATLAB writes into a temporary folder _syncaligned_tmp_<stamp> (no tool reads folders with that name);
  2) the result is checked: same rows and pupil columns as the package's table; every row's camera time equals the
     time sync_sessions.py gives that video frame (rows without a pulse: empty); CamTime_rel_s; the lick and valve
     line states at every camera pulse equal the Saleae lick / valve intervals of <session>_events.csv; one reward
     onset per Saleae valve opening inside the video;
  3) only then it is renamed concat_out_<stamp>_syncaligned (paths in its summary.txt / log.txt are updated) and a
     fix_source.json is written. A folder that fails a check is renamed FAILED_syncaligned_<stamp>.
The newest valid concat_out is the one prepare_step02_base.py copies to step02_base, so Step 02 then reads the
corrected table. Every fix attempt is logged in alignment_fixes_log.csv.

Safe to run every day (daily_summary.py runs it for all approved sessions): a session is redone only when its newest
concat_out is not an up-to-date sync-aligned copy of the package's newest concat_out (for example after a pupil redo,
which moves all concat_out folders to _redo_backup), or when the inputs' contents changed.

    python fix_concat_alignment.py --project <project> [--days 1 2] [--dry-run]
"""
import argparse, csv, hashlib, json, re, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import concat_choice as CC  # noqa: E402
SUFFIX = "_syncaligned"
DERIVED = ("_syncaligned", "_fullpupil")   # folders made by this tool / by fullpupil_sessions.py
TMP_PREFIX = "_syncaligned_tmp_"
MATLAB_TIMEOUT_S = 1800
PUPIL_COLS = ["Frame", "PupilTimestamp_s", "Diameter_px", "CenterX_px", "CenterY_px", "UsedScoreFallback",
              "PlotTime_s_30fps"]
DIGITAL_COLS = ["CamTime_s", "CamTime_rel_s", "Camera_TTL", "Lick_TTL", "Injector_TTL"]
MAX_STATE_MISMATCH = 3   # rows; a wrong frame<->pulse match gives hundreds
LOG_COLS = ["when", "day", "cage", "mouse", "status", "package_run", "package_mode", "package_offset", "fixed_run",
            "rows", "rows_without_pulse", "first_row_with_pulse", "last_row_with_pulse", "max_time_diff_s",
            "lick_state_mismatch", "valve_state_mismatch", "lick_onsets_package", "lick_onsets_fixed",
            "injector_onsets_package", "injector_onsets_fixed", "rewards_in_video_sync", "note"]


def matlab_exe():
    return shutil.which("matlab") or next((str(p) for p in sorted(Path("C:/Program Files/MATLAB").glob("R*/bin/matlab.exe"))), None)


def matlab_dir():
    """folder with pupil_align_ttl_sync.m: tools/opioidaddiction-matlab/current_cohort_workflow/matlab (project copy)
    or ../../matlab (this file inside the repo's current_cohort_workflow/python/daily)."""
    for c in (HERE / "opioidaddiction-matlab" / "current_cohort_workflow" / "matlab", HERE.parent.parent / "matlab"):
        if (c / "pupil_align_ttl_sync.m").exists():
            return c
    return None


def summary_field(run: Path, key: str) -> str:
    p = run / "summary.txt"
    if not p.exists():
        return ""
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(key + ":"):
            return line.split(":", 1)[1].strip()
    return ""


def fingerprint(p: Path) -> dict:
    h = hashlib.sha1()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return {"path": str(p), "size": p.stat().st_size, "sha1": h.hexdigest()}


def onsets(x) -> int:
    b = np.nan_to_num(np.asarray(x, float)) > 0.5
    return int(np.sum(b & ~np.r_[False, b[:-1]]))


def valid_fixed(run: Path) -> bool:
    return run.name.endswith(SUFFIX) and (run / "fix_source.json").exists() and (run / "combined_pupil_digital.csv").exists()


def runs_of(sdir: Path):
    """(newest checked sync-aligned copy or None, package folders oldest first). Folders made by fullpupil_sessions.py
    (_fullpupil) are neither: they are made FROM the newest of these."""
    allr = sorted(p for p in sdir.glob("concat_out_*") if p.is_dir())
    pkg = [p for p in allr if not p.name.endswith(DERIVED)]
    fixed = [p for p in allr if valid_fixed(p)]
    return (fixed[-1] if fixed else None), pkg


def state_at(times, intervals):
    """1 where a time lies inside one of the (on, off) intervals (on <= t < off; open interval if off is NaN)."""
    out = np.zeros(len(times), bool)
    if len(intervals) == 0:
        return out
    on = np.asarray([iv[0] for iv in intervals], float)
    off = np.asarray([iv[1] if np.isfinite(iv[1]) else np.inf for iv in intervals], float)
    o = np.argsort(on); on, off = on[o], off[o]
    k = np.searchsorted(on, times, side="right") - 1       # last interval that started at or before t
    ok = k >= 0
    out[ok] = times[ok] < off[k[ok]]
    return out


def verify(new: Path, pkg: Path, frames_csv: Path, events_csv) -> dict:
    """checks listed in the module docstring; raises ValueError on the first failure."""
    N = pd.read_csv(new / "combined_pupil_digital.csv")
    P = pd.read_csv(pkg / "combined_pupil_digital.csv")
    if list(N.columns) != list(P.columns):
        raise ValueError(f"columns differ from the package's table: {list(N.columns)}")
    if len(N) != len(P):
        raise ValueError(f"{len(N)} rows, the package's table has {len(P)}")
    for c in PUPIL_COLS:
        if not np.array_equal(N[c].to_numpy(float), P[c].to_numpy(float), equal_nan=True):
            raise ValueError(f"pupil column {c} differs from the package's table")
    F = pd.read_csv(frames_csv, usecols=["frame", "t_saleae_s"]).dropna()
    F["frame"] = F["frame"].astype(int)
    m = N[["Frame"] + DIGITAL_COLS].assign(frame=N["Frame"] - 1).merge(F, on="frame", how="left")
    has = m["t_saleae_s"].notna().to_numpy()
    if not has.any():
        raise ValueError("no row has a camera pulse in frames.csv")
    dt = float(np.max(np.abs(m.loc[has, "CamTime_s"] - m.loc[has, "t_saleae_s"])))
    if not np.isfinite(dt) or dt > 1e-6:
        raise ValueError(f"camera times differ from frames.csv by up to {dt:.6f} s")
    if m.loc[~has, DIGITAL_COLS].notna().any().any():
        raise ValueError("rows without a camera pulse have digital values")
    if m.loc[has, ["CamTime_s", "Camera_TTL", "Lick_TTL"]].isna().any().any():
        raise ValueError("rows with a camera pulse have empty digital values")
    if not (m.loc[has, "Camera_TTL"] == 1).all():
        raise ValueError("Camera_TTL is not 1 at every camera pulse")
    cam = m["CamTime_s"].to_numpy(float)
    rel = cam - cam[np.flatnonzero(has)[0]]
    drel = float(np.nanmax(np.abs(m["CamTime_rel_s"].to_numpy(float) - rel)))
    if drel > 1e-6:
        raise ValueError(f"CamTime_rel_s is off by up to {drel:.6f} s")
    rows = np.flatnonzero(has) + 1
    info = dict(rows=len(N), rows_without_pulse=int((~has).sum()), first_row_with_pulse=int(rows[0]),
                last_row_with_pulse=int(rows[-1]), max_time_diff_s=dt, lick_state_mismatch=np.nan,
                valve_state_mismatch=np.nan, lick_onsets_package=onsets(P["Lick_TTL"]),
                lick_onsets_fixed=onsets(N["Lick_TTL"]), injector_onsets_package=onsets(P["Injector_TTL"]),
                injector_onsets_fixed=onsets(N["Injector_TTL"]), rewards_in_video_sync=np.nan)
    if events_csv is not None:
        E = pd.read_csv(events_csv, usecols=["event", "source", "t_saleae_s", "t_off_saleae_s", "frame"])
        E = E[E["source"] == "saleae"]
        t = cam[has]
        for ev, col, key in (("lick", "Lick_TTL", "lick_state_mismatch"), ("reward_valve", "Injector_TTL", "valve_state_mismatch")):
            if N[col].isna().all():
                continue                                    # no such line in this digital.csv
            iv = E.loc[E["event"] == ev, ["t_saleae_s", "t_off_saleae_s"]].to_numpy(float)
            want = state_at(t, [tuple(r) for r in iv])
            got = m.loc[has, col].to_numpy(float) > 0.5
            info[key] = int(np.sum(want != got))
            if info[key] > MAX_STATE_MISMATCH:
                raise ValueError(f"{col} differs from the Saleae {ev} intervals at {info[key]} camera pulses")
        V = E[E["event"] == "reward_valve"]
        info["rewards_in_video_sync"] = int(V["frame"].notna().sum())
        if not N["Injector_TTL"].isna().all() and info["injector_onsets_fixed"] != info["rewards_in_video_sync"]:
            raise ValueError(f"{info['injector_onsets_fixed']} reward onsets, but {info['rewards_in_video_sync']} "
                             "Saleae valve openings inside the video")
    return info


def log_row(logp: Path, row: dict):
    try:
        write_header = not logp.exists()
        with open(logp, "a", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=LOG_COLS, extrasaction="ignore")
            if write_header:
                w.writeheader()
            w.writerow(row)
    except OSError as exc:   # e.g. the log is open in Excel: report, do not lose the result
        print(f"  !! could not write {logp.name} ({exc}); row: {row}")


def clean(s: str) -> str:
    return re.sub(r"[\x00-\x08\x0b-\x1f\x7f]", "", s).replace("\n", " ").strip()


def is_approved(v: str) -> bool:
    return CC.is_approved(v)          # 'yes...' or a date YYYY-MM-DD (concat_choice.py)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True, type=Path)
    ap.add_argument("--pupil-root", type=Path, help="default <project>/pupil_tst_final")
    ap.add_argument("--sync-root", type=Path, help="default <project>/sync")
    ap.add_argument("--days", type=int, nargs="*", default=[])
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    root = a.pupil_root or a.project / "pupil_tst_final"
    sync_root = a.sync_root or a.project / "sync"
    reg = root / "alignment_fixes.csv"
    if not reg.exists():
        print(f"alignment fixes: none approved ({reg.name} not found)")
        return 0
    try:
        R = pd.read_csv(reg, dtype=str, encoding="utf-8-sig").fillna("")
        R.columns = [str(c).strip().lower() for c in R.columns]
        if not {"day", "cage", "mouse"} <= set(R.columns):
            raise ValueError(f"columns {list(R.columns)}")
    except Exception as exc:   # a broken list must not stop the daily run (nothing is fixed then)
        print(f"  !! {reg.name} could not be read ({exc}); no alignment fix applied")
        return 1
    R = R[R["day"].str.strip() != ""]
    skipped = R[~R["approved"].map(is_approved)] if "approved" in R else R.iloc[:0]
    for _, r in skipped.iterrows():
        print(f"  day{r['day']} {r['cage']}_{r['mouse']}: 'approved' is empty or 'no' -> not applied")
    R = R[R["approved"].map(is_approved)] if "approved" in R else R.iloc[:0]
    if a.days:
        R = R[R["day"].astype(int).isin(a.days)]
    if R.empty:
        print("alignment fixes: no approved session for these days")
        return 0
    mexe, mdir = matlab_exe(), matlab_dir()
    logp = root / "alignment_fixes_log.csv"
    bad = 0
    for _, r in R.iterrows():
        day, cage, mouse = int(r["day"]), r["cage"].strip(), r["mouse"].strip()
        key = f"day{day} {cage}_{mouse}"
        sdir = root / f"day{day}" / cage / mouse
        newest, pkg = runs_of(sdir) if sdir.exists() else (None, [])   # newest = newest checked sync-aligned copy
        fcs = sorted((sync_root / f"day{day}" / cage / mouse).glob("*_frames.csv"))
        evs = sorted((sync_root / f"day{day}" / cage / mouse).glob("*_events.csv"))
        if not pkg or len(fcs) != 1:
            print(f"  !! {key}: {'no package concat_out' if not pkg else f'{len(fcs)} sync frames.csv'} -> not fixed")
            bad += 1
            continue
        P = pkg[-1]
        dig, pup = summary_field(P, "Digital file"), summary_field(P, "Pupil file")
        have_inputs = bool(dig and pup and Path(dig).exists() and Path(pup).exists())
        inputs = {"digital": fingerprint(Path(dig)), "pupil": fingerprint(Path(pup)),
                  "frames": fingerprint(fcs[0])} if have_inputs else None
        if newest is not None and valid_fixed(newest):
            src = json.loads((newest / "fix_source.json").read_text(encoding="utf-8"))
            if src.get("package_run") == P.name:
                old = src.get("inputs", {})
                same = inputs is not None and all(old.get(k, {}).get("sha1") == v["sha1"] for k, v in inputs.items())
                if same:
                    print(f"  {key}: up to date ({newest.name})")
                    continue
                if inputs is None:
                    print(f"  {key}: up to date ({newest.name}); the inputs are no longer at the paths in "
                          f"{P.name}\\summary.txt, so they could not be re-checked")
                    continue
        if not have_inputs:
            print(f"  !! {key}: inputs of {P.name} not found ({dig} | {pup}) -> not fixed")
            bad += 1
            continue
        if a.dry_run:
            print(f"  {key}: would write a sync-aligned copy of {P.name} "
                  f"(package: {summary_field(P, 'Alignment mode')}, offset {summary_field(P, 'Offset')})")
            continue
        if not mexe or not mdir:
            print(f"  !! {key}: {'MATLAB' if not mexe else 'pupil_align_ttl_sync.m'} not found -> not fixed")
            bad += 1
            continue
        stamp_s = datetime.now().strftime("%Y%m%d_%H%M%S")
        tmp, final = sdir / f"{TMP_PREFIX}{stamp_s}", sdir / f"concat_out_{stamp_s}{SUFFIX}"
        q = lambda s: str(s).replace("'", "''")  # noqa: E731
        cmd = (f"addpath('{q(mdir)}'); pupil_align_ttl_sync('{q(dig)}', '{q(pup)}', '{q(sdir)}', '{q(fcs[0])}', "
               f"'{q(tmp.name)}');")
        print(f"  {key}: writing a sync-aligned copy of {P.name} (MATLAB) ...", flush=True)
        row = dict(when=datetime.now().strftime("%Y-%m-%d %H:%M:%S"), day=day, cage=cage, mouse=mouse,
                   package_run=P.name, package_mode=summary_field(P, "Alignment mode"),
                   package_offset=summary_field(P, "Offset"), fixed_run=final.name, note=r.get("reason", ""))
        try:
            try:
                res = subprocess.run([mexe, "-batch", cmd], capture_output=True, text=True, errors="replace",
                                     timeout=MATLAB_TIMEOUT_S)
            except subprocess.TimeoutExpired:
                raise ValueError(f"MATLAB did not finish within {MATLAB_TIMEOUT_S} s")
            if res.returncode != 0 or not (tmp / "combined_pupil_digital.csv").exists():
                raise ValueError(f"MATLAB exit {res.returncode}: " + clean((res.stderr or res.stdout)[-600:]))
            info = verify(tmp, P, fcs[0], evs[0] if len(evs) == 1 else None)
            if final.exists():
                raise ValueError(f"{final.name} already exists")
            tmp.rename(final)                                   # only now does a reader see it
            for fn in ("summary.txt", "log.txt"):              # paths inside the files: temporary -> final name
                fp = final / fn
                if fp.exists():
                    fp.write_bytes(fp.read_bytes().replace(tmp.name.encode(), final.name.encode()))
            (final / "fix_source.json").write_text(json.dumps(
                {"package_run": P.name, "inputs": inputs, "approved": r.get("approved", ""),
                 "reason": r.get("reason", ""), "made_by": "fix_concat_alignment.py + pupil_align_ttl_sync.m",
                 "checks": info}, indent=1, default=float), encoding="utf-8")
            row.update(status="ok", **info)
            print(f"     ok: {final.name} | rows {info['rows']}, without pulse {info['rows_without_pulse']} "
                  f"(rows with pulse {info['first_row_with_pulse']}-{info['last_row_with_pulse']}) | "
                  f"lick/valve state mismatches {info['lick_state_mismatch']}/{info['valve_state_mismatch']} | "
                  f"lick onsets {info['lick_onsets_package']} -> {info['lick_onsets_fixed']}, "
                  f"reward onsets {info['injector_onsets_package']} -> {info['injector_onsets_fixed']} "
                  f"(Saleae valve openings inside the video: {info['rewards_in_video_sync']})")
        except Exception as exc:   # keep what was written for a look, under a name no tool reads
            bad += 1
            row.update(status=f"FAILED: {clean(str(exc))}", fixed_run="")
            where = tmp if tmp.exists() else None
            if where is not None:
                try:
                    where.rename(sdir / f"FAILED_syncaligned_{stamp_s}")
                except OSError:
                    pass                                        # the temporary name is not read by any tool either
            print(f"  !! {key}: FAILED ({clean(str(exc))})")
        log_row(logp, row)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
