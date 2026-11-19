#!/usr/bin/env python3
"""
darkcontrast_sessions.py - the repaired pupil (pupil_dark_contrast.py) as the pupil columns of the Step 02 table, for
the sessions the user approved (added 2026-11-18; the user: "제 추천은 9670f green2 5일에만 넣는 것 okay" = only
9670f green2 day 5). Same pattern as wholeeye_sessions.py / besteffort_sessions.py.

Approved sessions: rows of <pupil tree>/darkcontrast_sessions.csv (day,cage,mouse,approved,reason); the pupil tree =
concat_choice.active_pupil_root (pupil_active_root.txt). For every approved session:
  * the repair is (re)made when pupil_outputs_darkcontrast/pupil_dark_CLEAN.csv is missing or older than the
    package's results (pupil_dark_contrast.py --all); its source.json must say usable (else nothing is written);
  * base table = the session's table without pupil options (concat_choice.py: the package's concat_out or its
    approved sync-aligned copy);
  * a copy concat_out_<stamp>_darkcontrast in which ONLY Diameter_px (<- SmoothedDiameter_px), CenterX_px and
    CenterY_px are replaced, by the row's video Frame; every other column is copied as text, unchanged; the copy is
    checked row by row after writing, including that rows of frames the repair did not touch keep the base values;
  * fix_source.json records the base table, the SHA-1 of the repaired file and how many rows changed; an up-to-date
    copy is not made again.
The package's own output is never changed. Log: <pupil tree>/darkcontrast_log.csv.
usage: python darkcontrast_sessions.py --project <project> [--dry-run]
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import concat_choice as CC  # noqa: E402

PUPIL_COLS = {"Diameter_px": "SmoothedDiameter_px", "CenterX_px": "CenterX_px", "CenterY_px": "CenterY_px"}
SKIP_FILES = (".fig", ".png")          # the figure of the base table shows the old pupil: not copied
TOOL = Path(__file__).resolve().parent / "pupil_dark_contrast.py"


def read_text_csv(p: Path):
    raw = p.read_bytes()
    nl = "\r\n" if b"\r\n" in raw[:4096] else "\n"
    rows = list(csv.reader(raw.decode("utf-8").splitlines()))
    return rows[0], rows[1:], nl


def fmt(v):
    return "NaN" if not np.isfinite(v) else f"{v:.10g}"


def num(x):
    try:
        return float(x)
    except ValueError:
        return np.nan


def base_table(sd: Path, root: Path):
    runs = [p for p in CC.usable_runs(sd, root) if not p.name.endswith((CC.DARK, CC.WHOLE, CC.BEST, CC.FULL, CC.NOVID))]
    return runs[-1] if runs else None


def ensure_repair(sd: Path, day, cage, mouse, python_exe: str) -> Path:
    dc = sd / CC.DARK_DIR / CC.DARK_CSV
    pk = sd / "pupil_outputs_general" / "pupil_diameter_FINAL_U-Net_results2.csv"
    if not dc.exists() or (pk.exists() and pk.stat().st_mtime > dc.stat().st_mtime):
        print(f"  dark contrast day{day} {cage}_{mouse}: making the repair (pupil_dark_contrast.py --all) ...")
        rc = subprocess.run([python_exe, "-B", str(TOOL), "--day", str(day), "--cage", cage, "--mouse", mouse,
                             "--all", "--quiet"]).returncode
        if rc or not dc.exists():
            raise RuntimeError("pupil_dark_contrast.py failed")
    return dc


def one(sd: Path, root: Path, dry: bool, python_exe: str) -> str:
    key = CC.session_key(sd)
    dc = ensure_repair(sd, key[0], sd.parent.name, sd.name, python_exe) if not dry else sd / CC.DARK_DIR / CC.DARK_CSV
    sj = sd / CC.DARK_DIR / "source.json"
    info = json.loads(sj.read_text(encoding="utf-8")) if sj.exists() else {}
    if not info.get("usable"):
        raise RuntimeError(f"the repair is not usable for this session ({info.get('not_usable_because')}) - "
                           "nothing written; look at its check_sheet.png first")
    B = base_table(sd, root)
    if B is None:
        raise RuntimeError("no concat_out to start from")
    want = {"base_run": B.name, "dark_clean": CC.sha1(dc)}
    have = sorted(p for p in sd.glob(f"concat_out_*{CC.DARK}") if (p / "fix_source.json").exists())
    if have:
        s = json.loads((have[-1] / "fix_source.json").read_text(encoding="utf-8"))
        if s.get("base_run") == want["base_run"] and s.get("dark_clean") == want["dark_clean"]:
            return f"up to date ({have[-1].name})"
    if dry:
        return f"would write a dark-contrast copy of {B.name}"
    hdr, rows, nl = read_text_csv(B / "combined_pupil_digital.csv")
    T = pd.read_csv(dc)
    F = pd.read_csv(sd / CC.DARK_DIR / "pupil_dark_frames.csv", usecols=["frame", "trace_changed"])
    ci = {c: i for i, c in enumerate(hdr)}
    for c in list(PUPIL_COLS) + ["Frame"]:
        if c not in ci:
            raise RuntimeError(f"column {c} missing in {B.name}")
    src = {c: T[f].to_numpy(float) for c, f in PUPIL_COLS.items()}
    frames = np.array([num(r[ci["Frame"]]) for r in rows])
    ok = np.isfinite(frames) & (frames >= 1) & (frames <= len(T))
    if ok.mean() < 0.99:
        raise RuntimeError(f"only {100 * ok.mean():.1f} % of {B.name}'s rows have a video frame in the repaired trace")
    idx = np.where(ok, frames, 1).astype(int) - 1
    vals = {c: np.where(ok, v[idx], np.nan) for c, v in src.items()}
    # smoothing reaches 5 frames around a changed frame: rows whose frame is farther than that must keep the base value
    chg = F.trace_changed.to_numpy() > 0
    near = np.convolve(chg.astype(float), np.ones(11), mode="same") > 0
    touched = np.where(ok, near[idx], True)
    new = []
    for k, r in enumerate(rows):
        r2 = list(r)
        for c in PUPIL_COLS:
            r2[ci[c]] = fmt(vals[c][k])
        new.append(r2)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    tmp, final = sd / f"_darkcontrast_tmp_{stamp}", sd / f"concat_out_{stamp}{CC.DARK}"
    tmp.mkdir()
    try:
        for f in B.iterdir():
            if f.is_file() and f.name != "combined_pupil_digital.csv" and f.suffix.lower() not in SKIP_FILES:
                shutil.copy2(f, tmp / f.name)
        with open(tmp / "combined_pupil_digital.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, lineterminator=nl)
            w.writerow(hdr)
            w.writerows(new)
        h2, r2, _ = read_text_csv(tmp / "combined_pupil_digital.csv")
        keep = [i for c, i in ci.items() if c not in PUPIL_COLS]
        if h2 != hdr or len(r2) != len(rows) or any(a[i] != b[i] for a, b in zip(r2, rows) for i in keep):
            raise RuntimeError("non-pupil columns changed")
        diff_rows = big_rows = 0
        for c in PUPIL_COLS:
            got = np.array([num(x[ci[c]]) for x in r2])
            same = np.isclose(got, vals[c], equal_nan=True, rtol=1e-8, atol=1e-8)
            if not same.all():
                raise RuntimeError(f"{c} is not the repaired value in {int((~same).sum())} rows")
            base = np.array([num(x[ci[c]]) for x in rows])
            off = ~np.isclose(got, base, equal_nan=True, rtol=1e-7, atol=1e-6)
            if (off & ~touched).any():
                raise RuntimeError(f"{c} changed in {int((off & ~touched).sum())} rows the repair did not touch")
            diff_rows = max(diff_rows, int(off.sum()))
            if c == "Diameter_px":
                big_rows = int(np.sum(np.abs(got - base) > 1.0))
        srcj = dict(tool="darkcontrast_sessions.py", made=datetime.now().isoformat(timespec="seconds"), **want,
                    columns_replaced=list(PUPIL_COLS),
                    diameter_column="Diameter_px <- SmoothedDiameter_px of pupil_dark_CLEAN.csv (the package's own "
                                    "values on every frame the repair does not touch)",
                    rows_with_frame_pct=round(100 * float(ok.mean()), 2), rows_changed=diff_rows,
                    rows_changed_pct=round(100 * diff_rows / max(1, len(rows)), 2),
                    rows_diameter_changed_gt_1px=big_rows,
                    rows_diameter_changed_gt_1px_pct=round(100 * big_rows / max(1, len(rows)), 2),
                    repair=dict(replaced_pct=info.get("actions", {}).get("replaced_pct"),
                                size_factor=info.get("size_factor_package_over_dark"),
                                median_diam_px=info.get("median_diam_px"), made=info.get("made")),
                    note="pupil repaired where the package sat on another dark blob (user approved 2026-11-18); "
                         "all other columns copied unchanged from the base table; figures of the base table not copied")
        (tmp / "fix_source.json").write_text(json.dumps(srcj, indent=1), encoding="utf-8")
        tmp.rename(final)
    except Exception:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    return f"wrote {final.name} (from {B.name}; {diff_rows} of {len(rows)} rows have a changed pupil)"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--python", default=sys.executable, help="python for pupil_dark_contrast.py (needs cv2, scipy)")
    a = ap.parse_args(argv)
    root = CC.active_pupil_root(a.project)
    keys = CC.approved_sessions(root, CC.LISTS[CC.DARK])
    if not keys:
        print(f"  (no session approved for the dark-contrast repair in {root.name})")
        return 0
    log, bad = [], 0
    for day, cage, mouse in sorted(keys):
        sd = CC.session_dir(root, day, cage, mouse)
        if sd is None or not CC.package_runs(sd):
            print(f"  dark contrast day{day} {cage}_{mouse}: no package table (no video?) - nothing to do")
            continue
        try:
            msg = one(sd, root, a.dry_run, a.python)
        except Exception as exc:  # noqa: BLE001
            msg, bad = f"FAILED ({exc})", bad + 1
        print(f"  dark contrast day{day} {cage}_{mouse}: {msg}")
        log.append([datetime.now().isoformat(timespec="seconds"), day, cage, mouse, msg])
    if not a.dry_run and log:
        lp = root / "darkcontrast_log.csv"
        new = not lp.exists()
        with open(lp, "a", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            if new:
                w.writerow(["when", "day", "cage", "mouse", "result"])
            w.writerows(log)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
