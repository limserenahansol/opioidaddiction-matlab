"""step02_window.py - the Step 02 input of LISTED days cut to a time window (added 2026-11-26 for day 13).

Day 13 (2026-11-26) ran 20-min sessions (main_hansol.py SESSION_SECONDS = 20 * 60 + 10 = 1210 s). The user's rule,
same day: "PR you can use whole session whole 20min or whole 15min but for the others you can use only first 15min 10s".

The ORIGINAL Step 02 (unchanged) takes
  * PR from the task log: RequirementLast = the last requirement in the WHOLE .jsonl -> the .jsonl is not touched;
  * every other measure (licks, rewards via Injector_TTL, trials per frame, pupil) from
    BASE/dayN/<cage>/<mouse>/concat_out_*/combined_pupil_digital.csv -> this table is cut here.
For a day in tools/step02_window_days.csv (day, seconds, time_column, approved, reason) each such table keeps its rows
BEFORE the first row whose <time_column> is >= <seconds> (default CamTime_rel_s = camera time from the first camera
pulse, the same clock as quick_analysis' first_15m10s window). The kept rows are written byte for byte.

Nothing is lost: the whole table is first copied to BASE/_full_tables/dayN/<cage>/<mouse>/<concat_out>/ (outside the
dayN folders; Step 02 reads only folders named ^day\\d+$), and a marker combined_pupil_digital.window.json next to the
cut table records the window, the rows kept / dropped, the last kept time and the size + SHA-1 of the whole and of the
cut table. A table whose marker matches is left alone, so this can run after every prepare_step02_base.py
(daily_summary.py step "STEP 02 WINDOW"). A NEW concat_out of a session (package re-run) is copied whole by
prepare_step02_base.py and cut again here. A table whose marker does NOT match is not touched (reported).

Usage: python step02_window.py [--project P] [--apply]      (default: dry run, shows what would be cut)
Exit code 0 = fine, 1 = a problem was reported (nothing is cut for that table).
"""
import argparse
import csv
import hashlib
import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

P = Path(__file__).resolve().parent.parent
LIST = Path("tools") / "step02_window_days.csv"
MARK = "combined_pupil_digital.window.json"
LOG = "step02_window_log.csv"


def sha1(p: Path) -> str:
    h = hashlib.sha1()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_list(project: Path, probs: list) -> list:
    f = project / LIST
    if not f.exists():
        print(f"  ({LIST.as_posix()} not found - nothing to cut)")
        return []
    rows = []
    with open(f, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            try:
                d, s = int(str(r.get("day", "")).strip()), float(str(r.get("seconds", "")).strip())
            except ValueError:
                probs.append(f"{LIST.as_posix()}: unreadable row {dict(r)}")
                continue
            if not (r.get("approved") or "").strip():
                print(f"  day {d}: not approved in {LIST.as_posix()} - not cut")
                continue
            rows.append(dict(day=d, seconds=s, col=(r.get("time_column") or "").strip() or "CamTime_rel_s",
                             reason=(r.get("reason") or "").strip(), approved=r["approved"].strip()))
    return rows


def cut_lines(path: Path, col: str, seconds: float):
    """(header line, kept data lines, number of data lines, last kept time, first dropped time); lines as bytes, kept
    byte for byte. The table is in frame order, so everything from the first row at >= seconds on is dropped."""
    lines = path.read_bytes().splitlines(keepends=True)
    if not lines:
        raise ValueError("empty file")
    head = lines[0]
    names = next(csv.reader([head.decode("utf-8-sig").rstrip("\r\n")]))
    if col not in names:
        raise ValueError(f"no column {col} in the header")
    k = names.index(col)
    data = [ln for ln in lines[1:] if ln.strip()]
    kept, last, first_drop = [], None, None
    for ln in data:
        try:
            v = float(next(csv.reader([ln.decode("utf-8").rstrip("\r\n")]))[k])
        except (ValueError, IndexError):
            v = None                      # an empty time cell before the cut stays (frame order)
        if v is not None and v >= seconds:
            first_drop = v
            break
        kept.append(ln)
        if v is not None:
            last = v
    return head, kept, len(data), last, first_drop


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", default=str(P))
    ap.add_argument("--apply", action="store_true", help="cut (default: only show what would be cut)")
    a = ap.parse_args(argv)
    project = Path(a.project)
    base = project / "step02_base"
    probs, log_rows = [], []
    for D in read_list(project, probs):
        tables = sorted((base / f"day{D['day']}").glob("*/*/concat_out_*/combined_pupil_digital.csv"))
        print(f"  day {D['day']}: {len(tables)} Step 02 table(s); keep {D['col']} < {D['seconds']:g} s "
              f"(approved: {D['approved']})")
        if not tables:
            print("    (no table yet - the pupil step makes them)")
        for t in tables:
            cage, mouse, conc = t.parts[-4], t.parts[-3], t.parts[-2]
            key = f"day{D['day']} {cage}/{mouse} {conc}"
            mk = t.with_name(MARK)
            if mk.exists():
                try:
                    M = json.loads(mk.read_text(encoding="utf-8"))
                except ValueError:
                    M = {}
                if (M.get("seconds") == D["seconds"] and M.get("time_column") == D["col"]
                        and M.get("cut_sha1") == sha1(t)):
                    print(f"    {key}: already cut ({M.get('rows_kept')} of {M.get('rows_whole')} rows)")
                else:
                    probs.append(f"{key}: {MARK} does not match the table or the list - NOT touched, look by hand")
                continue
            try:
                head, kept, n_all, last, first_drop = cut_lines(t, D["col"], D["seconds"])
            except Exception as e:  # noqa: BLE001 - report and go on with the other tables
                probs.append(f"{key}: not cut ({' '.join(str(e).split())})")
                continue
            if first_drop is None:
                # 2026-11-26 (review): a table already shorter than the window (e.g. an early stop) gets a marker
                # too, so the evaluator's C6 sees that it was checked (rows kept = rows whole; the table is untouched)
                print(f"    {key}: all {n_all} rows are before {D['seconds']:g} s - nothing to cut"
                      + ("" if a.apply else " (would write the marker)"))
                if a.apply:
                    try:
                        sh = sha1(t)
                        M = dict(day=D["day"], cage=cage, mouse=mouse, concat_out=conc, time_column=D["col"],
                                 seconds=D["seconds"], rows_whole=n_all, rows_kept=n_all, last_kept_s=last,
                                 first_dropped_s=None, whole_table=str(t), whole_size=t.stat().st_size, whole_sha1=sh,
                                 cut_size=t.stat().st_size, cut_sha1=sh, approved=D["approved"], reason=D["reason"],
                                 note="nothing to cut: the whole table is shorter than the window",
                                 when=datetime.now().isoformat(timespec="seconds"), tool="tools/step02_window.py")
                        mk.write_text(json.dumps(M, indent=1), encoding="utf-8")
                    except OSError as e:
                        probs.append(f"{key}: marker not written ({e})")
                continue
            full = base / "_full_tables" / f"day{D['day']}" / cage / mouse / conc / t.name
            msg = (f"{n_all} -> {len(kept)} rows (last kept {D['col']} {last:.3f} s, first dropped {first_drop:.3f} s)"
                   if last is not None else f"{n_all} -> {len(kept)} rows")
            if not a.apply:
                print(f"    {key}: would cut {msg}; the whole table -> {full.relative_to(project).as_posix()}")
                continue
            # 2026-11-26 (review): a file held open elsewhere (viewer, Excel, MATLAB, virus scan) must not stop the
            # other tables -> report it and go on (the table stays whole; the next run cuts it)
            try:
                whole_sha, whole_size = sha1(t), t.stat().st_size
                if full.exists():
                    if sha1(full) != whole_sha:
                        probs.append(f"{key}: {full} exists and differs from this table - NOT cut, look by hand")
                        continue
                else:
                    full.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(t, full)
                    if sha1(full) != whole_sha:
                        probs.append(f"{key}: the copy of the whole table did not match - NOT cut")
                        continue
                tmp = t.with_name(t.name + ".cut_tmp")
                tmp.write_bytes(head + b"".join(kept))
                tmp.replace(t)
                M = dict(day=D["day"], cage=cage, mouse=mouse, concat_out=conc, time_column=D["col"],
                         seconds=D["seconds"], rows_whole=n_all, rows_kept=len(kept), last_kept_s=last,
                         first_dropped_s=first_drop, whole_table=str(full), whole_size=whole_size, whole_sha1=whole_sha,
                         cut_size=t.stat().st_size, cut_sha1=sha1(t), approved=D["approved"], reason=D["reason"],
                         when=datetime.now().isoformat(timespec="seconds"), tool="tools/step02_window.py")
                mk.write_text(json.dumps(M, indent=1), encoding="utf-8")
            except OSError as e:
                probs.append(f"{key}: NOT cut ({' '.join(str(e).split())}) - Step 02 would read the whole table; "
                             "close the program that holds it and run again")
                continue
            log_rows.append([M["when"], D["day"], cage, mouse, conc, D["col"], D["seconds"], n_all, len(kept),
                             last, first_drop, str(full)])
            print(f"    {key}: cut {msg}; the whole table is kept in {full.relative_to(project).as_posix()}")
    if log_rows:
        lp = base / LOG
        new = not lp.exists()
        with open(lp, "a", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            if new:
                w.writerow(["when", "day", "cage", "mouse", "concat_out", "time_column", "seconds", "rows_whole",
                            "rows_kept", "last_kept_s", "first_dropped_s", "whole_table"])
            w.writerows(log_rows)
    for p in probs:
        print("  !! " + p)
    return 1 if probs else 0


if __name__ == "__main__":
    sys.exit(main())
