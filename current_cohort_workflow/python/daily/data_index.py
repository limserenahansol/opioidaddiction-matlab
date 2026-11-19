#!/usr/bin/env python3
"""
data_index.py - ONE list of the data files the daily routine has saved, with paths relative to the project (added
2026-11-18; user: "not just figure but we are making csv or output data per day so in order to use and recall this
for 18days later all combined analysis as we do not want to run again this daily analysis for 18days").
The analysis after day 18 finds its inputs through this list and never re-runs sync, pupil tracking or the session
checks. Reads the project; writes only the index files.

Outputs
  analysis/data_index.csv                      all days (day = "all" for tables that hold every day)
  analysis/days/day<N>/data_index_day<N>.csv   the rows of day N
Columns: day, cage, mouse, session, kind, path (relative to the project, '/'), rows (data rows of a CSV), bytes,
  modified, final, note
  final = yes : the one to use where there is a choice (the Step 02 input that concat_choice.py picks = the newest
                usable concat_out with the approved options; the newest Step 02 run; the newest QC table)
          no  : kept, but not the one to use (an older concat_out, an older Step 02 run, ...)
          ""  : nothing to choose between (one file per session)
Kinds, per session: task_log, saleae_sal, saleae_csv, video, eye_box, sync_events, sync_frames, sync_inputs,
  sync_report, pupil_clean, pupil_qc, wholeeye_clean, step02_input, wholeeye_frames, wholeeye_qc, pupil_aligned,
  session_check, step02_base_input, step02_base_log, tst_score_step02 / hot_score_step02 / tail_immersion_step02 /
  assay (assays/dayN), tst_video, hot_video, hot_scoring_copy, hot_hand_score, tst_auto_result, tst_auto_info,
  tst_mouse_qc
per day: day_metrics, pr_checks, day_summary, assays_day, day_other
all days: master_<name> (analysis/*.csv), step02_table, step02_other, cohort_table, qc_summary, qc_issues,
  qc_decisions, qc_open_issues, next_steps, routine_log, session_checks_table, cohort_manifest, session_log, tail_immersion_table,
  arrangement_log
Not listed: figures, folders named _superseded* or _tests_and_aborted, the staging folder of a model retrain, the old
pupil tree (only the active tree of pupil_active_root.txt).
A session that has a task log but misses a usual file is printed as MISSING (report only; a behaviour-only session in
behavioronly_sessions.csv is not expected to have a video, pupil or eye box).

Loader for the later analysis:
  import data_index as DI
  DI.paths("sync_events", day=4)         -> [Path, ...]   (rows with final "yes" or "", unless final_only=False)
  DI.load_index()                        -> pandas DataFrame of analysis/data_index.csv (all text)
usage: python data_index.py [--project P] [--quiet]
"""
from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

T = Path(__file__).resolve().parent
sys.path.insert(0, str(T))
import concat_choice as CC  # noqa: E402

P = T.parent
COLS = ["day", "cage", "mouse", "session", "kind", "path", "rows", "bytes", "modified", "final", "note"]
DAY_RE = re.compile(r"^day(\d+)$")
NAME_RE = re.compile(r"^(\d+[mf])(black|white|green\d)", re.I)
SKIP = ("_superseded", "_tests_and_aborted")
MASTER = ("master_day_level", "session_metrics", "trials", "licks_long", "rewards_long", "per_minute",
          "peri_reward_psth")
EXPECT = ("task_log", "saleae_csv", "sync_events", "sync_frames", "session_check", "step02_input",
          "step02_base_input")
EXPECT_VIDEO = ("video", "eye_box", "pupil_clean")


def csv_rows(f: Path) -> int:
    n, last = 0, b"\n"
    with open(f, "rb") as fh:
        while True:
            b = fh.read(1 << 20)
            if not b:
                break
            n += b.count(b"\n")
            last = b[-1:]
    if last != b"\n":
        n += 1
    return max(0, n - 1)


class Index:
    def __init__(self, project: Path):
        self.p = project
        self.rows = []

    def add(self, f: Path, kind, day="", cage="", mouse="", session="", final="", note=""):
        if any(s in part for part in f.relative_to(self.p).parts for s in SKIP):
            return
        st = f.stat()
        self.rows.append(dict(day=day, cage=cage, mouse=mouse, session=session, kind=kind,
                              path=f.relative_to(self.p).as_posix(),
                              rows=csv_rows(f) if f.suffix.lower() == ".csv" else "", bytes=st.st_size,
                              modified=datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds"),
                              final=final, note=note))


def sessions(area: Path):
    """(day, cage, mouse, folder) of <area>/day<N>/<cage>/<mouse>/"""
    if not area.is_dir():
        return
    for d in sorted(area.glob("day*")):
        m = DAY_RE.match(d.name)
        if not m or not d.is_dir():
            continue
        for c in sorted(x for x in d.iterdir() if x.is_dir()):
            for mo in sorted(x for x in c.iterdir() if x.is_dir()):
                yield int(m.group(1)), c.name, mo.name, mo


def from_name(stem: str):
    """(cage, mouse, day) from a session name such as 9670fblack_prday4 or 9669mgreen1__fr3day1"""
    m, d = NAME_RE.match(stem), re.search(r"day(\d+)", stem, re.I)
    if not m:   # 2026-11-19 (user): cage and colour anywhere in the name, e.g. passive_9670fblack_day6
        m = re.search(r"(?<!\d)(\d{4}[mf]).*?(black|white|green\d)", stem, re.I)
    return (m.group(1).lower(), m.group(2).lower()) if m else ("", ""), (int(d.group(1)) if d else "")


def option_of(name: str) -> str:
    for s, txt in ((CC.WHOLE, "whole-eye copy"), (CC.DARK, "dark-contrast repair copy"), (CC.BEST, "best-effort copy"),
                   (CC.FULL, "full-pupil copy"),
                   (CC.NOVID, "behaviour-only copy"), (CC.SYNC, "sync-aligned copy")):
        if name.endswith(s):
            return txt
    return "package run"


def newest(paths):
    return max(paths, key=lambda x: x.name) if paths else None


def build(project: Path) -> pd.DataFrame:
    X = Index(project)
    ses = {}
    for day, cage, mouse, f in sessions(project / "python"):
        for g in sorted(f.glob("*.jsonl")):
            ses.setdefault((day, cage, mouse), []).append(g.stem)
            X.add(g, "task_log", day, cage, mouse, g.stem)
    sname = lambda d, c, m: ";".join(ses.get((d, c, m), []))   # noqa: E731
    stem2 = {st: k for k, v in ses.items() for st in v}   # session name -> (day, cage, mouse) of its python folder

    def by_name(stem, day=""):
        """(day, cage, mouse, final, note) of a file named after a session: the python folder of the task log with
        that name wins (names typed wrong at the rig, e.g. day 2 '9670fgreen_fr5day2_1' = 9670f green1)"""
        if stem in stem2:
            return (*stem2[stem], "", "")
        (c, m), d = from_name(stem)
        return day or d, c, m, "no", "no task log of this cohort with this name"
    for day, cage, mouse, f in sessions(project / "salae"):
        for g in sorted(f.glob("*.sal")):
            X.add(g, "saleae_sal", day, cage, mouse, sname(day, cage, mouse))
        for g in sorted(f.glob("*.csv")):
            X.add(g, "saleae_csv", day, cage, mouse, sname(day, cage, mouse))
    for day, cage, mouse, f in sessions(project / "video"):
        for g in sorted(f.glob("*.mp4")):
            X.add(g, "video", day, cage, mouse, sname(day, cage, mouse))
        if (f / "roi.mat").exists():
            X.add(f / "roi.mat", "eye_box", day, cage, mouse, sname(day, cage, mouse))
    sync_kinds = (("_events.csv", "sync_events"), ("_frames.csv", "sync_frames"), ("_inputs.json", "sync_inputs"),
                  ("_report.txt", "sync_report"))
    for day, cage, mouse, f in sessions(project / "sync"):
        for g in sorted(f.iterdir()):
            for suf, kind in sync_kinds:
                if g.is_file() and g.name.endswith(suf):
                    X.add(g, kind, day, cage, mouse, g.name[:-len(suf)])
    root = CC.active_pupil_root(project)
    approvals = {s: CC.approved_sessions(root, n, warn=None) for s, n in CC.LISTS.items()}
    for day, cage, mouse, sd in sessions(root):
        s = sname(day, cage, mouse)
        gen = sd / "pupil_outputs_general"
        if (gen / "pupil_diameter_FINAL_U-Net_results_CLEAN.csv").exists():
            X.add(gen / "pupil_diameter_FINAL_U-Net_results_CLEAN.csv", "pupil_clean", day, cage, mouse, s,
                  note=f"package pupil trace ({root.name})")
        if (gen / "qc.json").exists():
            X.add(gen / "qc.json", "pupil_qc", day, cage, mouse, s, note="package eye box + quality numbers")
        we = sd / CC.WHOLE_DIR / CC.WHOLE_CSV
        if we.exists():
            X.add(we, "wholeeye_clean", day, cage, mouse, s, note="whole eye in the package CLEAN format")
        chosen = CC.newest_usable(sd, root, approvals)
        for run in sorted(sd.glob("concat_out_*")):
            t = run / "combined_pupil_digital.csv"
            if not t.exists():
                continue
            ok = chosen is not None and run.name == chosen.name
            X.add(t, "step02_input", day, cage, mouse, s, "yes" if ok else "no",
                  option_of(run.name) + ("" if ok else "; not the one Step 02 uses (older, or option not approved)"))
    for day, cage, mouse, f in sessions(project / "analysis" / "wholeeye"):
        if (f / "wholeeye_frames.csv").exists():
            X.add(f / "wholeeye_frames.csv", "wholeeye_frames", day, cage, mouse, sname(day, cage, mouse),
                  note="eye model output per video frame")
        if (f / "wholeeye_qc.json").exists():
            X.add(f / "wholeeye_qc.json", "wholeeye_qc", day, cage, mouse, sname(day, cage, mouse),
                  note="which eye model, crop, % found")
    for g in sorted((project / "analysis" / "pupil_aligned").glob("*_pupil_aligned.csv")):
        stem = g.name[:-len("_pupil_aligned.csv")]
        day, cage, mouse, fin, note = by_name(stem)
        X.add(g, "pupil_aligned", day, cage, mouse, stem, fin, note)
    sc = project / "analysis" / "session_checks"
    for d in sorted(sc.glob("day*")):
        m = DAY_RE.match(d.name)
        if not m:
            continue
        for g in sorted(d.glob("*_check.txt")):
            stem = g.stem[:-len("_check")]
            day, cage, mouse, fin, note = by_name(stem, int(m.group(1)))
            X.add(g, "session_check", day, cage, mouse, stem, fin, note)
    base = project / "step02_base"
    for day, cage, mouse, f in sessions(base):
        for t in sorted(f.glob("concat_out_*/combined_pupil_digital.csv")):
            X.add(t, "step02_base_input", day, cage, mouse, sname(day, cage, mouse), "yes",
                  "copy that the original Step 02 reads (" + option_of(t.parent.name) + ")")
        for g in sorted(f.glob("*.jsonl")):
            X.add(g, "step02_base_log", day, cage, mouse, g.stem, "yes", "task log copy that Step 02 reads")
    runs = sorted((base / "longitudinal_outputs").glob("run_*"))
    for r in runs:
        last = r == runs[-1]
        for g in sorted(r.rglob("*.csv")):
            X.add(g, "step02_table" if g.name == "ALL_mice_longitudinal.csv" else "step02_other", "all",
                  final="yes" if last else "no", note=f"Step 02 {r.name}" + ("" if last else " (older run)"))
    cruns = sorted((project / "analysis" / "cohort_figs").glob("run_*"))
    for r in cruns:
        last = r == cruns[-1]
        for g in sorted(r.rglob("*.csv")):
            X.add(g, "cohort_table", "all", final="yes" if last else "no",
                  note=f"previous-cohort figures {r.name}" + ("" if last else " (older run)"))
    days_dir = project / "analysis" / "days"
    for d in sorted(days_dir.glob("day*")):
        m = DAY_RE.match(d.name)
        if not m or not d.is_dir():
            continue
        n = int(m.group(1))
        known = {f"day{n}_metrics.csv": "day_metrics", f"pr_checks_day{n}.csv": "pr_checks",
                 f"day{n}_summary.txt": "day_summary", f"assays_day{n}.csv": "assays_day"}
        for g in sorted(d.rglob("*")):
            if not g.is_file() or "figs" in g.relative_to(d).parts or g.name.startswith("data_index"):
                continue
            if g.name in known:
                X.add(g, known[g.name], n)
            elif g.suffix.lower() in (".csv", ".json"):
                (cage, mouse), _ = from_name(g.parent.name.replace("frame_drop_", ""))
                X.add(g, "day_other", n, cage, mouse, note=g.parent.name if g.parent != d else "")
    for stem in MASTER:
        g = project / "analysis" / f"{stem}.csv"
        if g.exists():
            X.add(g, f"master_{stem}", "all", final="yes", note="rebuilt from the per-session files every day")
    for pat, kind in (("qc_summary_days*.csv", "qc_summary"), ("qc_issues_days*.csv", "qc_issues"),
                      ("pr_table_days*.csv", "pr_table")):          # pr_table: pr_table_days.py (2026-11-18)
        fs = sorted(days_dir.glob(pat), key=lambda x: int(re.findall(r"(\d+)\.csv$", x.name)[0]))
        for g in fs:
            X.add(g, kind, "all", final="yes" if g == fs[-1] else "no")
    # pupil check (2026-11-18, pupil_check_examples.py): frames out of the pupil the tables use + Claude's verdicts
    pc = project / "analysis" / "pupil_check"
    fs = sorted(pc.glob("pupil_check_days*.csv"), key=lambda x: int(re.findall(r"(\d+)\.csv$", x.name)[0])) \
        if pc.is_dir() else []
    for g in fs:
        X.add(g, "pupil_check", "all", final="yes" if g == fs[-1] else "no",
              note="frames out per session (held, jumps, odd size / no eye, off eye) + Claude's check")
    if (pc / "pupil_check_verdicts.csv").exists():
        X.add(pc / "pupil_check_verdicts.csv", "pupil_check_verdicts", "all", note="Claude's check per session")
    for name, kind in (("qc_decisions.csv", "qc_decisions"), ("qc_open_issues.csv", "qc_open_issues"),
                       ("next_steps.csv", "next_steps"), ("routine_log.csv", "routine_log")):
        if (days_dir / name).exists():
            X.add(days_dir / name, kind, "all")
    for g, kind in ((sc / "session_checks.csv", "session_checks_table"), (project / "cohort_manifest.csv",
                     "cohort_manifest"), (project / "session_log.csv", "session_log"),
                    (project / "_arrangement_log.csv", "arrangement_log"),
                    (project / "assays" / "tail_immersion_all.csv", "tail_immersion_table")):
        if g.exists():
            X.add(g, kind, "all")
    for day, cage, mouse, f in sessions(project / "assays"):
        for g in sorted(x for x in f.rglob("*") if x.is_file()):
            kind = "assay"
            low = g.name.lower()
            if low.startswith("manualscoringresults") and "tst" in low:
                kind = "tst_score_step02"          # _AUTO = automatic (assays_day.py), else a hand score
            elif low.startswith("manualscoringresults") and "hot" in low:
                kind = "hot_score_step02"
            elif "immer" in low:
                kind = "tail_immersion_step02"
            X.add(g, kind, day, cage, mouse, sname(day, cage, mouse),
                  note="automatic TST, first 180 s" if "_auto" in low else "")
    # assay videos (2026-11-18): originals copied by assays_day.py, HOT first-30-s scoring copies + hand scores, TST
    # automatic results (the chosen run is the AutoScoringResults file in the video's folder root)
    for root_name, kind in (("tst_video", "tst_video"), ("hp_video", "hot_video")):
        for day, cage, mouse, f in sessions(project / root_name):
            for g in sorted(f.glob("*")):
                if g.is_file() and g.suffix.lower() in (".mp4", ".avi", ".mov", ".mkv"):
                    X.add(g, kind, day, cage, mouse, sname(day, cage, mouse))
    for d in sorted((project / "hp_video").glob("day*_first30s")) if (project / "hp_video").exists() else []:
        m = re.match(r"day(\d+)_first30s$", d.name)
        if not m:
            continue
        n = int(m.group(1))
        for g in sorted(d.glob("*.mp4")) + sorted((d / "output").glob("ManualScoringResults_*.mat")):
            (cage, mouse), _ = from_name(g.stem.replace("ManualScoringResults_", ""))
            X.add(g, "hot_scoring_copy" if g.suffix.lower() == ".mp4" else "hot_hand_score", n, cage, mouse,
                  sname(n, cage, mouse), note="first 30 s")
    for d in sorted((project / "analysis" / "tst").glob("day*")) if (project / "analysis" / "tst").exists() else []:
        m = DAY_RE.match(d.name)
        if not m:
            continue
        n = int(m.group(1))
        for vd in sorted(x for x in d.iterdir() if x.is_dir() and x.name != "examples"):
            (cage, mouse), _ = from_name(vd.name)
            for g, kind in ((vd / f"AutoScoringResults_{vd.name}.mat", "tst_auto_result"),
                            (vd / f"{vd.name}_lock.json", "tst_auto_info"),
                            (vd / f"{vd.name}_mouse_qc.csv", "tst_mouse_qc")):
                if g.exists():
                    X.add(g, kind, n, cage, mouse, sname(n, cage, mouse),
                          "yes" if kind == "tst_auto_result" else "", "first 180 s, chosen run")
    D = pd.DataFrame(X.rows, columns=COLS)
    return D, approvals


def missing(D: pd.DataFrame, approvals) -> list:
    out = []
    novid = approvals.get(CC.NOVID, set())
    S = D[D.day.astype(str) != "all"]
    for (day, cage, mouse), g in S.groupby(["day", "cage", "mouse"]):
        if "task_log" not in set(g.kind):
            continue
        need = list(EXPECT)
        if CC.norm_key(day, cage, mouse) not in novid:
            need += list(EXPECT_VIDEO)
        have = set(g.kind[g.final != "no"])
        lost = [k for k in need if k not in have]
        if lost:
            out.append(f"  MISSING day{day} {cage}_{mouse}: {', '.join(lost)}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", default=str(P))
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    project = Path(a.project)
    D, approvals = build(project)
    out = project / "analysis" / "data_index.csv"
    tmp = out.with_name(out.name + ".tmp")
    D.to_csv(tmp, index=False)
    tmp.replace(out)
    for n in sorted({d for d in D.day if str(d) != "all"}, key=int):
        dd = project / "analysis" / "days" / f"day{n}"
        if dd.is_dir():
            D[D.day.astype(str) == str(n)].to_csv(dd / f"data_index_day{n}.csv", index=False)
    print(f"  data index: {len(D)} files -> {out.relative_to(project)} (+ analysis/days/dayN/data_index_dayN.csv)")
    if not a.quiet:
        k = D.groupby("kind").size()
        print("  " + ", ".join(f"{i} {v}" for i, v in k.items()))
    for line in missing(D, approvals):
        print(line)
    return 0


def load_index(project=None) -> pd.DataFrame:
    p = Path(project) if project else P
    return pd.read_csv(p / "analysis" / "data_index.csv", dtype=str, keep_default_na=False)


def paths(kind, day=None, cage=None, mouse=None, final_only=True, project=None) -> list:
    p = Path(project) if project else P
    D = load_index(p)
    D = D[D.kind == kind]
    for col, v in (("day", day), ("cage", cage), ("mouse", mouse)):
        if v is not None:
            D = D[D[col] == str(v)]
    if final_only:
        D = D[D.final != "no"]
    return [p / x for x in D.path]


if __name__ == "__main__":
    raise SystemExit(main())
