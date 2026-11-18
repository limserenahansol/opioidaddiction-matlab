#!/usr/bin/env python3
"""
besteffort_sessions.py - the best-effort pupil of hard eyes as an OPTION for listed sessions (added 2026-11-16).

The user (2026-11-16) approved putting the validated best-effort traces (tools/pupil_besteffort.py) into the analysis
tables and Step 02. Only sessions approved in <pupil_root>/besteffort_sessions.csv (day,cage,mouse,approved,reason)
are touched; approved = 'yes...' or a date. Removing a row (or 'no ...') takes the session back to the package's table.

For every approved session:
  * base table = the session's current table without this option (concat_choice.py: the package's concat_out, or its
    approved sync-aligned copy);
  * a copy concat_out_<stamp>_besteffort in which ONLY the pupil columns (Diameter_px, CenterX_px, CenterY_px,
    UsedScoreFallback) are replaced row by row from pupil_outputs_besteffort/besteffort_CLEAN.csv (frames where the
    pupil was not visible -> NaN). Every other column (Frame, PupilTimestamp_s, CamTime_s, Camera_TTL, Lick_TTL,
    Injector_TTL, ...) is copied as text, unchanged; the copy is checked row by row after writing.
  * fix_source.json records the base table and the SHA-1 of the trace; an up-to-date copy is not made again.
The package's own output is never changed. Log: <pupil_root>/besteffort_log.csv.

usage: python besteffort_sessions.py --project <project> [--dry-run]
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import concat_choice as CC  # noqa: E402

PUPIL_COLS = {"Diameter_px": "Diameter_pixels", "CenterX_px": "CenterX_px", "CenterY_px": "CenterY_px"}
SKIP_FILES = (".fig", ".png")          # the figure of the base table shows the old pupil: not copied


def read_text_csv(p: Path):
    raw = p.read_bytes()
    nl = "\r\n" if b"\r\n" in raw[:4096] else "\n"
    rows = list(csv.reader(raw.decode("utf-8").splitlines()))
    return rows[0], rows[1:], nl


def fmt(v):
    return "NaN" if not np.isfinite(v) else f"{v:.10g}"


def base_table(sd: Path, root: Path):
    runs = [p for p in CC.usable_runs(sd, root) if not p.name.endswith(CC.BEST)]
    return runs[-1] if runs else None


def one(sd: Path, root: Path, dry: bool) -> str:
    be = sd / "pupil_outputs_besteffort" / "besteffort_CLEAN.csv"
    chk = sd / "pupil_outputs_besteffort" / "besteffort_check.json"
    if not (be.exists() and chk.exists()):
        raise RuntimeError("no best-effort trace (run pupil_besteffort.py first)")
    B = base_table(sd, root)
    if B is None:
        raise RuntimeError("no concat_out to start from")
    want = {"base_run": B.name, "besteffort_clean": CC.sha1(be)}
    have = sorted(p for p in sd.glob(f"concat_out_*{CC.BEST}") if (p / "fix_source.json").exists())
    if have:
        s = json.loads((have[-1] / "fix_source.json").read_text(encoding="utf-8"))
        if s.get("base_run") == want["base_run"] and s.get("besteffort_clean") == want["besteffort_clean"]:
            return f"up to date ({have[-1].name})"
    if dry:
        return f"would write a best-effort copy of {B.name}"
    hdr, rows, nl = read_text_csv(B / "combined_pupil_digital.csv")
    T = pd.read_csv(be)
    if len(T) != len(rows):
        raise RuntimeError(f"the best-effort trace has {len(T)} rows, {B.name} has {len(rows)}")
    ci = {c: i for i, c in enumerate(hdr)}
    for c in list(PUPIL_COLS) + ["UsedScoreFallback"]:
        if c not in ci:
            raise RuntimeError(f"column {c} missing in {B.name}")
    vals = {c: T[f].to_numpy(float) for c, f in PUPIL_COLS.items()}
    new = []
    for k, r in enumerate(rows):
        r2 = list(r)
        for c in PUPIL_COLS:
            r2[ci[c]] = fmt(vals[c][k])
        r2[ci["UsedScoreFallback"]] = "0"
        new.append(r2)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    tmp, final = sd / f"_besteffort_tmp_{stamp}", sd / f"concat_out_{stamp}{CC.BEST}"
    tmp.mkdir()
    try:
        for f in B.iterdir():
            if f.is_file() and f.name != "combined_pupil_digital.csv" and f.suffix.lower() not in SKIP_FILES:
                shutil.copy2(f, tmp / f.name)
        with open(tmp / "combined_pupil_digital.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, lineterminator=nl)
            w.writerow(hdr)
            w.writerows(new)
        # check: same header and rows; every non-pupil column identical as text; pupil = the trace
        h2, r2, _ = read_text_csv(tmp / "combined_pupil_digital.csv")
        keep = [i for c, i in ci.items() if c not in PUPIL_COLS and c != "UsedScoreFallback"]
        if h2 != hdr or len(r2) != len(rows) or any(a[i] != b[i] for a, b in zip(r2, rows) for i in keep):
            raise RuntimeError("non-pupil columns changed")
        for c in PUPIL_COLS:
            got = np.array([float(x[ci[c]]) for x in r2])
            ok = np.isclose(got, vals[c], equal_nan=True, rtol=1e-8, atol=1e-8)
            if not ok.all():
                raise RuntimeError(f"{c} is not the best-effort value in {int((~ok).sum())} rows")
        src = dict(tool="besteffort_sessions.py", made=datetime.now().isoformat(timespec="seconds"), **want,
                   besteffort_check=json.loads(chk.read_text(encoding="utf-8")).get("besteffort_kfold"),
                   columns_replaced=list(PUPIL_COLS) + ["UsedScoreFallback"],
                   visible_frames_pct=round(100 * float(np.isfinite(vals["Diameter_px"]).mean()), 1),
                   note="all other columns copied unchanged from the base table; figures of the base table not copied")
        (tmp / "fix_source.json").write_text(json.dumps(src, indent=1), encoding="utf-8")
        tmp.rename(final)
    except Exception:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    return f"wrote {final.name} (from {B.name}; pupil visible in {src['visible_frames_pct']} % of rows)"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    root = CC.active_pupil_root(a.project)    # 2026-11-17: the tree the tables use (pupil_active_root.txt)
    keys = CC.approved_sessions(root, CC.LISTS[CC.BEST])
    if not keys:
        print("  (no session approved in besteffort_sessions.csv)")
        return 0
    log, bad = [], 0
    for day, cage, mouse in sorted(keys):
        sd = CC.session_dir(root, day, cage, mouse)
        if sd is None:
            print(f"  !! day{day} {cage}_{mouse}: no pupil folder")
            bad += 1
            continue
        try:
            msg = one(sd, root, a.dry_run)
        except Exception as exc:  # noqa: BLE001
            msg, bad = f"FAILED ({exc})", bad + 1
        print(f"  best effort day{day} {cage}_{mouse}: {msg}")
        log.append([datetime.now().isoformat(timespec="seconds"), day, cage, mouse, msg])
    if not a.dry_run and log:
        lp = root / "besteffort_log.csv"
        new = not lp.exists()
        with open(lp, "a", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            if new:
                w.writerow(["when", "day", "cage", "mouse", "result"])
            w.writerows(log)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
