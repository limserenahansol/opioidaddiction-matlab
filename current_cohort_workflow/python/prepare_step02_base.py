"""Build the folder that the ORIGINAL Step 02 reads (BASE\\dayN\\<cage>\\<mouse>\\concat_out_*\\ + the session's .jsonl
+ assay files), from the outputs of the pupil/TST package and the arranged task logs. Nothing is changed in the source
folders: files are copied (small CSV/JSONL/MAT files only).

    python prepare_step02_base.py --pupil-out <pupil_out_root> --python-root <python> --base <step02_base>
        [--days 1 2] [--assays <assay_root>] [--apply]

Rules (same meaning as in the previous cohort):
  * one session per mouse-day: only the NEWEST concat_out_* of the package output is copied. Step 02 counts every
    concat_out folder as a separate session, so copying several would duplicate the session (as happened in the
    previous cohort: 6099_white day 3 x3, 0911_orange day 7 x2).
  * the session's table is chosen by daily/concat_choice.py (the same rule as every other tool): the package's
    newest concat_out, or a corrected copy - concat_out_<stamp>_syncaligned (fix_concat_alignment.py) or
    concat_out_<stamp>_fullpupil (fullpupil_sessions.py) - that passed its checks (fix_source.json) and whose session is
    approved NOW in <pupil-out>/alignment_fixes.csv or fullpupil_sessions.csv. Corrected copies are options for rare
    sessions; removing a session from its list takes it back to the package's table.
    A listed session whose corrected copy is missing (e.g. after a pupil redo) is NOT updated until its tool has run;
    a session in both lists needs a full-pupil copy built on its sync-aligned table.
  * an older copy already in BASE (replaced by a corrected copy, a pupil redo, or a withdrawn approval) is moved to
    BASE/_superseded_<date>/dayN/<cage>/<mouse>/ (kept, outside the dayN folders that Step 02 reads); every move
    is also appended to BASE/step02_base_moves.csv.
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
from datetime import date, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "daily"))
import concat_choice as CC  # noqa: E402  (which concat_out is the session's table; approvals)

SYNC_SUFFIX, FULL_SUFFIX, DERIVED = CC.SYNC, CC.FULL, CC.DERIVED
ASSAY_PATTERNS = ("ManualScoringResults*tst*.mat", "ManualScoringResults*hot*.mat", "ManualScoringResults*straub*.mat",
                  "AutoScoringResults*tst*.mat", "*immersion*.xlsx", "*immersion*.csv", "*Immersion*.xlsx", "*Immersion*.csv")


def unique_dir(p: Path) -> Path:
    """p, or p_2, p_3, ... if it exists (a second move of a folder with the same name must not nest into the first)."""
    if not p.exists():
        return p
    k = 2
    while (p.parent / f"{p.name}_{k}").exists():
        k += 1
    return p.parent / f"{p.name}_{k}"


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

    rows, problems, moves = [], [], []
    approvals = {suf: CC.approved_sessions(a.pupil_out, name) for suf, name in CC.LISTS.items()}   # never raises
    for di, ddir in day_dirs(a.pupil_out, a.days):
        for cage in sorted(p for p in ddir.iterdir() if p.is_dir()):
            for mouse in sorted(p for p in cage.iterdir() if p.is_dir()):
                key = f"day{di} {cage.name}_{mouse.name}"
                sk = CC.norm_key(di, cage.name, mouse.name)
                allc = [c for c in sorted(mouse.glob("concat_out_*")) if (c / "combined_pupil_digital.csv").exists()]
                usable = CC.usable_runs(mouse, a.pupil_out, approvals)
                for c in allc:
                    if c not in usable:   # a corrected copy that is not checked or not approved now is not used
                        why = "has no fix_source.json (not checked)" if not (c / "fix_source.json").exists() else                             "is not approved now (fix list)"
                        problems.append(f"{key}: {c.name} {why} -> not used")
                if not usable:
                    problems.append(f"{key}: no concat_out_* with combined_pupil_digital.csv (pupil step not done?)")
                    continue
                newest = usable[-1]
                if len(usable) > 1:
                    problems.append(f"{key}: {len(usable)} usable concat_out folders; the session's table is {newest.name}")
                need_full, need_sync = sk in approvals[FULL_SUFFIX], sk in approvals[SYNC_SUFFIX]
                need = FULL_SUFFIX if need_full else (SYNC_SUFFIX if need_sync else "")
                if need and not newest.name.endswith(need):
                    # e.g. the package redid the session after the fix: its new table is END-aligned / U-Net-only again
                    tool = "fullpupil_sessions.py" if need == FULL_SUFFIX else "fix_concat_alignment.py"
                    problems.append(f"!! {key}: listed {need[1:]} copy missing for the newest {newest.name} -> "
                                    f"this session is NOT updated; run {tool} first")
                    continue
                if need_full and need_sync and not str(CC._source(newest).get("base_run", "")).endswith(SYNC_SUFFIX):
                    problems.append(f"!! {key}: in both lists, but {newest.name} is not built on the sync-aligned table -> "
                                    "NOT updated; run fix_concat_alignment.py, then fullpupil_sessions.py")
                    continue
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
                    # BASE has an older copy than the newest source table (a sync-aligned copy replaced the package's
                    # table, or the package redid the session): keep it, but outside the dayN folders so Step 02 does
                    # not count it as a 2nd session
                    why = ("the sync-aligned " if newest.name.endswith(SYNC_SUFFIX) else
                           "the full-pupil " if newest.name.endswith(FULL_SUFFIX) else "the package's ")
                    if a.apply:
                        tgt = unique_dir(a.base / f"_superseded_{date.today():%Y%m%d}" / ddir.name / cage.name /
                                         mouse.name / c.name)
                        tgt.parent.mkdir(parents=True, exist_ok=True)
                        shutil.move(str(c), str(tgt))
                        mv = dict(when=f"{datetime.now():%Y-%m-%d %H:%M:%S}", session=key, kind="superseded",
                                  source=str(c), target=str(tgt), action="moved", replaced_by=newest.name)
                        rows.append(mv); moves.append(mv)
                        problems.append(f"{key}: {c.name} replaced by {why}{newest.name}; old copy moved to {tgt}")
                    else:
                        problems.append(f"{key}: {c.name} would be moved to _superseded_<date> (replaced by {why}{newest.name})")
    if a.apply:
        a.base.mkdir(parents=True, exist_ok=True)
        with open(a.base / "step02_base_manifest.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["session", "kind", "source", "target", "action"], extrasaction="ignore")
            w.writeheader(); w.writerows(rows)
        if moves:   # the manifest is rewritten on every run; moves are also kept in an append-only log
            mp = a.base / "step02_base_moves.csv"
            new = not mp.exists()
            with open(mp, "a", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=["when", "session", "source", "target", "replaced_by"], extrasaction="ignore")
                if new:
                    w.writeheader()
                w.writerows(moves)
    n_ses = len({r["session"] for r in rows})
    print(f"{n_ses} sessions, {len(rows)} files ({'applied' if a.apply else 'dry run'})")
    for p in problems:
        print("NOTE:", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
