#!/usr/bin/env python3
"""
wholeeye_sessions.py - the whole eye as the pupil columns of the Step 02 table, for mice with eye damage (added
2026-11-17; the user: "9670 black has eye damage ... use whole eye as alternate", "9669 black also ... both",
"네, 넣어 주세요"). Same pattern as besteffort_sessions.py.

Approved sessions: every session of the mice in <pupil tree>/wholeeye_mice.csv (cage,mouse,approved,reason) plus rows
of wholeeye_sessions.csv (day,cage,mouse,approved,reason); a 'no' row there takes that session back to the package's
table. The pupil tree = concat_choice.active_pupil_root (pupil_active_root.txt).

For every approved session that has pupil_outputs_wholeeye/wholeeye_CLEAN.csv (wholeeye_clean.py):
  * base table = the session's table without options of the pupil columns (concat_choice.py: the package's concat_out
    or its approved sync-aligned copy);
  * a copy concat_out_<stamp>_wholeeye in which ONLY Diameter_px (<- SmoothedDiameter_px, gap-filled like the
    package's), CenterX_px, CenterY_px and UsedScoreFallback (0) are replaced, by the row's video Frame; every other
    column is copied as text, unchanged; the copy is checked row by row after writing;
  * fix_source.json records the base table and the SHA-1 of the whole-eye file; an up-to-date copy is not made again.
The package's own output is never changed. Log: <pupil tree>/wholeeye_log.csv.
usage: python wholeeye_sessions.py --project <project> [--dry-run]
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

PUPIL_COLS = {"Diameter_px": "SmoothedDiameter_px", "CenterX_px": "CenterX_px", "CenterY_px": "CenterY_px"}
SKIP_FILES = (".fig", ".png")          # the figure of the base table shows the old pupil: not copied


def read_text_csv(p: Path):
    raw = p.read_bytes()
    nl = "\r\n" if b"\r\n" in raw[:4096] else "\n"
    rows = list(csv.reader(raw.decode("utf-8").splitlines()))
    return rows[0], rows[1:], nl


def fmt(v):
    return "NaN" if not np.isfinite(v) else f"{v:.10g}"


def base_table(sd: Path, root: Path):
    runs = [p for p in CC.usable_runs(sd, root) if not p.name.endswith((CC.WHOLE, CC.BEST, CC.FULL, CC.NOVID))]
    return runs[-1] if runs else None


def one(sd: Path, root: Path, dry: bool) -> str:
    we = sd / CC.WHOLE_DIR / CC.WHOLE_CSV
    if not we.exists():
        raise RuntimeError("no whole-eye trace (run eye_track_all.m and wholeeye_clean.py first)")
    B = base_table(sd, root)
    if B is None:
        raise RuntimeError("no concat_out to start from")
    want = {"base_run": B.name, "wholeeye_clean": CC.sha1(we)}
    have = sorted(p for p in sd.glob(f"concat_out_*{CC.WHOLE}") if (p / "fix_source.json").exists())
    if have:
        s = json.loads((have[-1] / "fix_source.json").read_text(encoding="utf-8"))
        if s.get("base_run") == want["base_run"] and s.get("wholeeye_clean") == want["wholeeye_clean"]:
            return f"up to date ({have[-1].name})"
    if dry:
        return f"would write a whole-eye copy of {B.name}"
    hdr, rows, nl = read_text_csv(B / "combined_pupil_digital.csv")
    T = pd.read_csv(we)
    ci = {c: i for i, c in enumerate(hdr)}
    for c in list(PUPIL_COLS) + ["UsedScoreFallback", "Frame"]:
        if c not in ci:
            raise RuntimeError(f"column {c} missing in {B.name}")
    src = {c: T[f].to_numpy(float) for c, f in PUPIL_COLS.items()}
    frames = np.array([float(r[ci["Frame"]]) if r[ci["Frame"]] not in ("", "NaN") else np.nan for r in rows])
    ok = np.isfinite(frames) & (frames >= 1) & (frames <= len(T))
    if ok.mean() < 0.99:
        raise RuntimeError(f"only {100 * ok.mean():.1f} % of {B.name}'s rows have a video frame in the whole-eye trace")
    idx = np.where(ok, frames, 1).astype(int) - 1
    vals = {c: np.where(ok, v[idx], np.nan) for c, v in src.items()}
    new = []
    for k, r in enumerate(rows):
        r2 = list(r)
        for c in PUPIL_COLS:
            r2[ci[c]] = fmt(vals[c][k])
        r2[ci["UsedScoreFallback"]] = "0"
        new.append(r2)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    tmp, final = sd / f"_wholeeye_tmp_{stamp}", sd / f"concat_out_{stamp}{CC.WHOLE}"
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
        keep = [i for c, i in ci.items() if c not in PUPIL_COLS and c != "UsedScoreFallback"]
        if h2 != hdr or len(r2) != len(rows) or any(a[i] != b[i] for a, b in zip(r2, rows) for i in keep):
            raise RuntimeError("non-pupil columns changed")
        for c in PUPIL_COLS:
            got = np.array([float(x[ci[c]]) for x in r2])
            same = np.isclose(got, vals[c], equal_nan=True, rtol=1e-8, atol=1e-8)
            if not same.all():
                raise RuntimeError(f"{c} is not the whole-eye value in {int((~same).sum())} rows")
        js = sd / CC.WHOLE_DIR / "wholeeye_source.json"
        info = json.loads(js.read_text(encoding="utf-8")) if js.exists() else {}
        srcj = dict(tool="wholeeye_sessions.py", made=datetime.now().isoformat(timespec="seconds"), **want,
                    columns_replaced=list(PUPIL_COLS) + ["UsedScoreFallback"],
                    diameter_column="Diameter_px <- SmoothedDiameter_px of wholeeye_CLEAN.csv (equivalent diameter "
                                    "of the eye opening, gap-filled)",
                    rows_with_frame_pct=round(100 * float(ok.mean()), 2), eye_found_pct=info.get("found_pct"),
                    model=info.get("model"),
                    note="whole eye instead of the pupil (eye damage, user 2026-11-17); all other columns copied "
                         "unchanged from the base table; figures of the base table not copied")
        (tmp / "fix_source.json").write_text(json.dumps(srcj, indent=1), encoding="utf-8")
        tmp.rename(final)
    except Exception:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    return f"wrote {final.name} (from {B.name}; eye found in {info.get('found_pct')} % of frames)"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    root = CC.active_pupil_root(a.project)
    keys = CC.approved_sessions(root, CC.LISTS[CC.WHOLE])
    if not keys:
        print(f"  (no session approved for the whole eye in {root.name})")
        return 0
    log, bad = [], 0
    for day, cage, mouse in sorted(keys):
        sd = CC.session_dir(root, day, cage, mouse)
        if sd is None:
            continue
        if not CC.package_runs(sd):
            print(f"  whole eye day{day} {cage}_{mouse}: no package table (no video?) - nothing to do")
            continue
        try:
            msg = one(sd, root, a.dry_run)
        except Exception as exc:  # noqa: BLE001
            msg, bad = f"FAILED ({exc})", bad + 1
        print(f"  whole eye day{day} {cage}_{mouse}: {msg}")
        log.append([datetime.now().isoformat(timespec="seconds"), day, cage, mouse, msg])
    if not a.dry_run and log:
        lp = root / "wholeeye_log.csv"
        new = not lp.exists()
        with open(lp, "a", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            if new:
                w.writerow(["when", "day", "cage", "mouse", "result"])
            w.writerows(log)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
