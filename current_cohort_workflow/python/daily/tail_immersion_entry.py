#!/usr/bin/env python3
"""
tail_immersion_entry.py - the tail-immersion latencies the user reports each day (from day 5, 2026-11-18; user: "from
today day 5 I will have every day tail immersion task so ill let you know every day the latency of tail (As what we
did also for previous cohort)") written into the files Step 02 reads, as the previous cohort had them.

Per mouse and day:  assays/day<N>/<cage>/<mouse>/<cage><mouse>_immersion_day<N>.csv
    row 1  Immersion_Latency_s,<latency>   Step 02 (findImmersionLatency) reads ONLY cell B1 of the newest *immer* file
    row 2  trials_s,<each value>           only when the user gave single trials too
    row 3  entered,<when>   row 4  source,<text>
prepare_step02_base.py (--assays <project>/assays, every daily run) copies it next to the task log in step02_base, so
the Step 02 table gets Immersion_Latency_s for that mouse-day (Step 02 flags values outside 0-60 s).
All days in one table: assays/tail_immersion_all.csv (day, cage, mouse, latency_s, trials_s, note, entered, source).

Rules: the value that goes into B1 is the one the user gives (--value). Single trials alone (--trials) are not turned
into a latency unless the user said how (--summary mean|median). A mouse-day that already has a DIFFERENT value is not
changed: both are printed, and --replace is needed after the user says which one is right. The mouse must have a task
log folder python/day<N>/<cage>/<mouse> - or, while day N is not filed yet (the latencies come during the day, the
evening routine files the logs), one on an earlier day.
usage:
  python tail_immersion_entry.py --day 5 --value 9670f black 3.2 --value 9670f green1 2.8 [--source "chat 11-18"] [--apply]
  python tail_immersion_entry.py --day 5 --trials 9670f black 3.1 3.3 --summary mean [--apply]
Dry run unless --apply.
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from datetime import datetime
from pathlib import Path

P = Path(__file__).resolve().parent.parent
MASTER_COLS = ["day", "cage", "mouse", "latency_s", "trials_s", "note", "entered", "source"]


def num(s: str) -> float:
    v = float(str(s).strip().rstrip("sS").replace(",", "."))
    if not (0 <= v <= 60):
        raise ValueError(f"{s}: outside 0-60 s (Step 02 would flag it)")
    return v


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--day", type=int, required=True)
    ap.add_argument("--value", nargs=3, action="append", default=[], metavar=("CAGE", "MOUSE", "LATENCY"))
    ap.add_argument("--trials", nargs="+", action="append", default=[], metavar="CAGE MOUSE V1 [V2 ...]")
    ap.add_argument("--summary", choices=("mean", "median"), help="only when the user said how trials become one value")
    ap.add_argument("--note", default="")
    ap.add_argument("--source", default="user message")
    ap.add_argument("--replace", action="store_true", help="only after the user said the new value is the right one")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--project", default=str(P))
    a = ap.parse_args(argv)
    project = Path(a.project)
    day = a.day
    entries, bad = {}, 0
    for cage, mouse, v in a.value:
        try:
            entries[(cage.lower(), mouse.lower())] = {"latency": num(v), "trials": []}
        except ValueError as exc:
            print(f"  !! day{day} {cage} {mouse}: {exc} - not written")
            bad += 1
    for t in a.trials:
        if len(t) < 3:
            ap.error("--trials needs CAGE MOUSE and at least one value")
        k = (t[0].lower(), t[1].lower())
        try:
            vals = [num(x) for x in t[2:]]
        except ValueError as exc:
            print(f"  !! day{day} {t[0]} {t[1]}: trial {exc} - not written")
            bad += 1
            continue
        e = entries.setdefault(k, {"latency": None, "trials": []})
        e["trials"] = vals
    day_arranged = (project / "python" / f"day{day}").is_dir()
    for k, e in entries.items():
        # the mouse must exist: in that day's task-log folders, or - before the evening routine has filed that day
        # (the latencies come during the day) - in the folders of an earlier day
        here = (project / "python" / f"day{day}" / k[0] / k[1]).is_dir()
        before = sorted(int(p.parent.parent.name[3:]) for p in (project / "python").glob(f"day*/{k[0]}/{k[1]}")
                        if p.is_dir() and p.parent.parent.name[3:].isdigit() and int(p.parent.parent.name[3:]) < day)
        if not here and (day_arranged or not before):
            print(f"  !! day{day} {k[0]} {k[1]}: no task-log folder python/day{day}/{k[0]}/{k[1]}"
                  f"{'' if day_arranged else ' or on an earlier day'} - check the name (not written)")
            bad += 1
            e["skip"] = True
            continue
        if not here:
            e["known"] = f"day{day} not filed yet; mouse known from day{before[-1]}"
        if e["latency"] is None:
            if not a.summary:
                print(f"  !! day{day} {k[0]} {k[1]}: only single trials {e['trials']} - which value goes into the table "
                      f"(the user decides: --value, or --summary mean|median)?")
                bad += 1
                continue
            e["latency"] = (statistics.mean if a.summary == "mean" else statistics.median)(e["trials"])
            e["how"] = f"{a.summary} of {len(e['trials'])} trials"
    master = project / "assays" / "tail_immersion_all.csv"
    rows = []
    if master.exists():
        with open(master, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
    have = {(int(r["day"]), r["cage"], r["mouse"]): r for r in rows}
    now = datetime.now().isoformat(timespec="seconds")
    todo = []
    for (cage, mouse), e in sorted(entries.items()):
        if e["latency"] is None or e.get("skip"):
            continue
        old = have.get((day, cage, mouse))
        if old is not None and abs(float(old["latency_s"]) - e["latency"]) > 1e-9 and not a.replace:
            print(f"  !! day{day} {cage} {mouse}: the table already has {old['latency_s']} s ({old['entered']}, "
                  f"{old['source']}); new value {e['latency']} s NOT written - which one is right? (--replace)")
            bad += 1
            continue
        todo.append((cage, mouse, e))
    for cage, mouse, e in todo:
        f = project / "assays" / f"day{day}" / cage / mouse / f"{cage}{mouse}_immersion_day{day}.csv"
        trials = " ".join(f"{x:g}" for x in e["trials"])
        print(f"  {'write' if a.apply else 'would write'} day{day} {cage} {mouse}: {e['latency']:g} s"
              f"{' (trials ' + trials + ')' if trials else ''}{' = ' + e['how'] if e.get('how') else ''}"
              f"{' [' + e['known'] + ']' if e.get('known') else ''} -> "
              f"{f.relative_to(project).as_posix()}")
        if not a.apply:
            continue
        f.parent.mkdir(parents=True, exist_ok=True)
        with open(f, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["Immersion_Latency_s", f"{e['latency']:g}"])
            if e["trials"]:
                w.writerow(["trials_s", *[f"{x:g}" for x in e["trials"]]])
            w.writerow(["entered", now])
            w.writerow(["source", a.source + (f"; {e['how']}" if e.get("how") else "")])
        have[(day, cage, mouse)] = dict(day=day, cage=cage, mouse=mouse, latency_s=f"{e['latency']:g}", trials_s=trials,
                                        note=a.note, entered=now,
                                        source=a.source + (f"; {e['how']}" if e.get("how") else ""))
    if a.apply and todo:
        master.parent.mkdir(parents=True, exist_ok=True)
        tmp = master.with_name(master.name + ".tmp")
        with open(tmp, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=MASTER_COLS)
            w.writeheader()
            for k in sorted(have, key=lambda x: (int(x[0]), x[1], x[2])):
                w.writerow({c: have[k].get(c, "") for c in MASTER_COLS})
        tmp.replace(master)
        print(f"  table: {master.relative_to(project).as_posix()} ({len(have)} mouse-days)")
    if not a.apply:
        print("  dry run - nothing written")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
