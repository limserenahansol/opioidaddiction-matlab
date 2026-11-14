"""Build the folder that the ORIGINAL Step 02 reads (BASE\\dayN\\<cage>\\<mouse>\\concat_out_*\\ + the session's .jsonl
+ assay files), from the outputs of the pupil/TST package and the arranged task logs. Nothing is changed in the source
folders: files are copied (small CSV/JSONL/MAT files only).

    python prepare_step02_base.py --pupil-out <pupil_out_root> --python-root <python> --base <step02_base>
        [--days 1 2] [--assays <assay_root>] [--apply]

Rules (same meaning as in the previous cohort):
  * one session per mouse-day: only the NEWEST concat_out_* of the package output is copied. Step 02 counts every
    concat_out folder as a separate session, so copying several would duplicate the session (as happened in the
    previous cohort: 6099_white day 3 x3, 0911_orange day 7 x2).
  * one .jsonl per mouse-day: the task log arranged in <python-root>\\dayN\\<cage>\\<mouse>\\. Step 02 picks the newest
    .jsonl whose name contains the mouse or cage name, so a second log in the folder would be ambiguous: it stops.
  * assays (optional): files already named the way Step 02 looks for them are copied from
    <assays>\\dayN\\<cage>\\<mouse>\\:  ManualScoringResults*tst*.mat / *hot*.mat / *straub*.mat (or the TST package's
    AutoScoringResults*tst*.mat, see tst_step02_copy.py) and the tail-immersion file.
Dry run unless --apply. Writes step02_base_manifest.csv into --base.
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
from pathlib import Path

ASSAY_PATTERNS = ("ManualScoringResults*tst*.mat", "ManualScoringResults*hot*.mat", "ManualScoringResults*straub*.mat",
                  "AutoScoringResults*tst*.mat", "*immersion*.xlsx", "*immersion*.csv", "*Immersion*.xlsx", "*Immersion*.csv")


def day_dirs(root: Path, days):
    out = []
    for d in sorted(root.glob("day*")):
        m = re.fullmatch(r"day(\d+)", d.name)
        if d.is_dir() and m and (not days or int(m.group(1)) in days):
            out.append((int(m.group(1)), d))
    return sorted(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pupil-out", required=True, type=Path, help="pupil_out_root of the pupil/TST package")
    ap.add_argument("--python-root", required=True, type=Path, help="arranged task logs: <root>\\dayN\\<cage>\\<mouse>\\*.jsonl")
    ap.add_argument("--base", required=True, type=Path, help="Step 02 BASE folder to (re)fill")
    ap.add_argument("--assays", type=Path, default=None, help="optional: <root>\\dayN\\<cage>\\<mouse>\\ assay files")
    ap.add_argument("--days", type=int, nargs="*", default=[])
    ap.add_argument("--apply", action="store_true", help="copy (default: only show what would be copied)")
    a = ap.parse_args(argv)

    rows, problems = [], []
    for di, ddir in day_dirs(a.pupil_out, a.days):
        for cage in sorted(p for p in ddir.iterdir() if p.is_dir()):
            for mouse in sorted(p for p in cage.iterdir() if p.is_dir()):
                key = f"day{di} {cage.name}_{mouse.name}"
                cos = sorted(mouse.glob("concat_out_*"))
                cos = [c for c in cos if (c / "combined_pupil_digital.csv").exists()]
                if not cos:
                    problems.append(f"{key}: no concat_out_* with combined_pupil_digital.csv (pupil step not done?)")
                    continue
                newest = cos[-1]
                if len(cos) > 1:
                    problems.append(f"{key}: {len(cos)} concat_out folders; only the newest {newest.name} is used")
                logs = sorted((a.python_root / ddir.name / cage.name / mouse.name).glob("*.jsonl"))
                if len(logs) != 1:
                    problems.append(f"{key}: {len(logs)} .jsonl files in python folder (need exactly 1) -> skipped")
                    continue
                dst = a.base / ddir.name / cage.name / mouse.name
                todo = [(newest, dst / newest.name, "concat_out"), (logs[0], dst / logs[0].name, "jsonl")]
                if a.assays:
                    src_a = a.assays / ddir.name / cage.name / mouse.name
                    for pat in ASSAY_PATTERNS:
                        for f in sorted(src_a.glob(pat)):
                            todo.append((f, dst / f.name, "assay"))
                for src, tgt, kind in todo:
                    exists = tgt.exists()
                    action = "exists" if exists else ("copied" if a.apply else "would copy")
                    if a.apply and not exists:
                        tgt.parent.mkdir(parents=True, exist_ok=True)
                        if src.is_dir():
                            shutil.copytree(src, tgt)
                        else:
                            shutil.copy2(src, tgt)
                    rows.append(dict(session=key, kind=kind, source=str(src), target=str(tgt), action=action))
                old = [c for c in dst.glob("concat_out_*") if c.name != newest.name] if dst.exists() else []
                for c in old:
                    problems.append(f"{key}: BASE already has another {c.name}: Step 02 would count it as a 2nd session - remove it")
    if a.apply:
        a.base.mkdir(parents=True, exist_ok=True)
        with open(a.base / "step02_base_manifest.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["session", "kind", "source", "target", "action"])
            w.writeheader(); w.writerows(rows)
    n_ses = len({r["session"] for r in rows})
    print(f"{n_ses} sessions, {len(rows)} files ({'applied' if a.apply else 'dry run'})")
    for p in problems:
        print("NOTE:", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
