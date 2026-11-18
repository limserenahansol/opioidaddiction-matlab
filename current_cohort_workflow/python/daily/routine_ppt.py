#!/usr/bin/env python3
"""routine_ppt.py - the daily routine PowerPoint (user 2026-11-17: "PPT as routine summary ... QC, problem fixed,
folder address where I can find all plots and representative plots as fundamental, what to do better next";
complicated analysis waits until after day 18). Report only: reads what the routine made, changes nothing.

  1  overview: counts, today's PR, decisions needed       2  QC per session (green / yellow / red)
  3  data problems found and what was done                 4  routine problems fixed (rig, files, software)
  5  where to find all plots                                6-10  basic plots (previous-cohort figures + daily figures)
  11 what to do better next
  (+ after slide 2, on assay days: TI / TST / HOT of the day from analysis\\days\\dayN\\assays_dayN.csv; added 2026-11-18)

Sources (inside the project): analysis\\days\\qc_summary_days*.csv and qc_issues_days*.csv (tools\\qc_summary.py),
analysis\\days\\routine_log.csv and next_steps.csv (kept by Claude), analysis\\days\\dayN\\pr_checks_dayN.csv,
analysis\\reference\\prev_cohort_pr_days3to5.csv, analysis\\figs_across_days\\, analysis\\cohort_figs\\run_###\\,
analysis\\days\\dayN\\figs\\.
usage:  python routine_ppt.py [--days 1 2 3 4] [--out FILE.pptx]
"""
from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
import make_summary_ppt as MS  # noqa: E402  (look and helpers of the first PPT tool)
import qc_summary as QS  # noqa: E402

P = TOOLS.parent
DAYS_DIR = P / "analysis" / "days"
NAVY, TEXT, MUTED, LINE = MS.NAVY, MS.TEXT, MS.MUTED, MS.LINE
MICE = QS.MICE
QC_FILL = {"green": "E6F4EA", "yellow": "FFF6DB", "red": "FDE2E1"}
QC_TEXT = {"green": "1E6B34", "yellow": "7A5A00", "red": "B42318"}
STATUS = {  # fill, text colour, label
    "open": ("FDE2E1", "B42318", "needs you / a look"),
    "fixed": ("E6F4EA", "1E6B34", "fixed"),
    "approved": ("ECE8F8", "4B3A8C", "approved option"),
    "noted": ("FFF6DB", "7A5A00", "known, no change"),
    "done": ("E6F4EA", "1E6B34", "done"),
}
ORDER = ["open", "fixed", "approved", "noted", "done"]
# qc_issues category -> (problem, what was done, status). None as problem = use the issue's own title.
CAT = {
    "off_eye": ("Whole eye: when the eye closes, the mask lands on the spout / whiskers / fur",
                "Needs your decision: count those frames as 'no eye'?", "open"),
    "off_eye_kept": ("Whole eye: when the eye closes, the mask lands on the spout / whiskers / fur",
                     "Your decision 11-18: keep those frames as they are (qc_decisions.csv)", "approved"),
    "open_issue": (None, "Needs your decision: is the white part eye?", "open"),
    "rewards": ("Rewards: Python and Saleae differ", "Needs a look", "open"),
    "frames_open": ("Video frames differ from camera pulses (not fixed)", "Needs a look", "open"),
    "no_video_open": ("No video found", "Needs a look", "open"),
    "lick_unknown": ("Task counted a lick with no Saleae lick near it", "Needs a look", "open"),
    "check_open": ("Session-check item not explained", "Needs a look", "open"),
    "pupil_open": ("Pupil grade missing or bad", "Needs a look", "open"),
    "frames": ("Video frames differ from camera pulses", "Sync-aligned table used in Step 02 (you approved)", "fixed"),
    "clock": ("PC clock jumped during the session (Windows time sync)", "The sync uses the camera pulse train", "fixed"),
    "lick_valve_close": ("Lick counted at valve close (Arduino re-trigger, lock-out 0.20 s)",
                         "Lock-out 0.22 s since then; none counted at 0.22 s", "fixed"),
    "no_video": ("No video saved", "Behaviour-only table (you approved); pupil left empty", "approved"),
    "whole_eye": ("Eye damage: the pupil is not visible", "Whole eye instead of the pupil (you approved 11-17)",
                  "approved"),
    "length": ("Session longer than 15 min 10 s", "All daily metrics use the first 15 min 10 s", "noted"),
    "lick_start": ("A lick from before the start was counted at t = 0 (serial buffer)",
                   "Known; the previous cohort had it too (7 of 42 PR sessions); task code kept for comparability",
                   "noted"),
    "lick_late": ("A real lick (on the Saleae) that Python read late", "Nothing to change", "noted"),
    "licks_extra": ("Python logged extra licks while the 500-ms valve was open (not counted)",
                    "Tables use the Saleae lick line; 200-ms valve since day 3", "noted"),
    "filename": ("File name without the mouse", "Filed by start time (the only free mouse)", "noted"),
    "pupil_noisy": ("Pupil noisy (spout / whiskers cross the eye)", "Graded by eye; used, marked noisy", "noted"),
}
SHOW_SHORT = {"frames", "frames_open", "length", "off_eye", "off_eye_kept", "licks_extra", "lick_start", "lick_valve_close",
              "lick_late", "lick_unknown", "pupil_noisy", "rewards", "filename", "check_open"}


# ---------------------------------------------------------------- helpers
def picture(slide, path: Path, x, y, w, h, name):
    from PIL import Image
    with Image.open(path) as im:
        pw, ph = im.size
    sc = min(w / pw, h / ph)
    W, H = pw * sc, ph * sc
    pic = slide.shapes.add_picture(str(path), Inches(x + (w - W) / 2), Inches(y + (h - H) / 2), Inches(W), Inches(H))
    pic.name = name
    return pic


def table(slide, x, y, widths, header, rows, name, size=10, head_size=11):
    """rows: list of lists of (text, fill, colour, bold); returns the estimated height (in)"""
    heights = []
    for r in rows:
        n = max(MS.text_lines(c[0], w, size) for c, w in zip(r, widths))
        heights.append(max(0.27, 0.07 + 0.165 * n * size / 10))
    head_h = 0.32
    shp = slide.shapes.add_table(len(rows) + 1, len(header), Inches(x), Inches(y), Inches(sum(widths)),
                                 Inches(head_h + sum(heights)))
    shp.name = name
    t = shp.table
    for j, w in enumerate(widths):
        t.columns[j].width = Inches(w)
    t.rows[0].height = Inches(head_h)
    for j, h in enumerate(header):
        MS.style_cell(t.cell(0, j), [[(h, {"bold": True, "color": "FFFFFF"})]], NAVY, size=head_size)
    for i, (r, h) in enumerate(zip(rows, heights), start=1):
        t.rows[i].height = Inches(h)
        for j, (txt, fill, col, bold) in enumerate(r):
            MS.style_cell(t.cell(i, j), [[(txt, {"color": col, "bold": bold})]], fill, size=size)
    return head_h + sum(heights)


def legend(slide, x, y, items):
    for fill, col, label in items:
        sw = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y + 0.04), Inches(0.24), Inches(0.18))
        sw.name = f"Legend {label}"
        sw.fill.solid()
        sw.fill.fore_color.rgb = MS.rgb(fill)
        sw.line.color.rgb = MS.rgb(LINE)
        sw.shadow.inherit = False
        w = 0.075 * len(label) + 0.15
        MS.add_text(slide, x + 0.3, y, w, 0.26, label, size=10, color=MUTED, name=f"Legend text {label}")
        x += 0.3 + w + 0.2


def box(slide, x, y, w, h, title, paragraphs, accent, name, size=12):
    card = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    card.name = name
    card.adjustments[0] = 0.05
    card.fill.solid()
    card.fill.fore_color.rgb = MS.rgb("F6F7F9")
    card.line.fill.background()
    card.shadow.inherit = False
    MS.add_text(slide, x + 0.2, y + 0.12, w - 0.4, 0.34, title, size=14, color=accent, bold=True, name=name + " title")
    MS.add_text(slide, x + 0.2, y + 0.5, w - 0.4, h - 0.6, paragraphs, size=size, color=TEXT, name=name + " text",
                spacing_after=4)


def two_up(prs, a: Path, b: Path, title, caption, stacked=False):
    """two figures on one slide: side by side, or one above the other (stacked); missing files are skipped"""
    have = [f for f in (a, b) if f.exists()]
    if not have:
        return
    s = prs.slides.add_slide(prs.slide_layouts[6])
    MS.add_text(s, 0.5, 0.25, 12.33, 0.5, title, size=22, color=NAVY, bold=True, name="Title")
    MS.add_text(s, 0.5, 0.78, 12.33, 0.34, caption, size=12, color=MUTED, name="Caption")
    for k, f in enumerate(have):
        if stacked:
            picture(s, f, 0.5, 1.2 + k * 2.98, 12.33, 2.9, f.stem)
        else:
            picture(s, f, 0.5 + k * 6.25, 1.25, 6.08, 5.8, f.stem)
    MS.add_text(s, 0.5, 7.18, 12.33, 0.25, str(have[0].parent), size=8, color=MUTED, name="Source")


def sess_text(g, cat):
    """sessions of one problem, grouped by mouse: '9670f white d1 (-212 frames); 9669m green2 d4 (-33 frames)'"""
    parts = []
    for (c, m), q in g.groupby(["cage", "mouse"], sort=False):
        ds = [f"d{int(r.day)}" + (f" ({r.short})" if cat in SHOW_SHORT and str(r.short) not in ("", "nan") else "")
              for r in q.itertuples()]
        parts.append(f"{c} {m}: " + ", ".join(ds))
    return "; ".join(parts)


def pr_table(days):
    """{day: DataFrame(cage, color, rewards, pr)} from pr_checks_dayN.csv (PR = pr_RL_full = last requirement)"""
    out = {}
    for d in days:
        f = DAYS_DIR / f"day{d}" / f"pr_checks_day{d}.csv"
        if f.exists():
            X = pd.read_csv(f, encoding="utf-8-sig", dtype={"cage": str})
            X["color"] = X["mouse"].astype(str).str.split("_").str[-1]
            out[d] = X
    return out


def newest_run():
    runs = sorted((P / "step02_base" / "longitudinal_outputs").glob("run_*"))
    return runs[-1] if runs else None


def cohort_dir(run):
    d = P / "analysis" / "cohort_figs" / run.name / "bout1s_excl4-6-11-14-17" if run else None
    return d if d and d.exists() else None


def count_files(folder: Path, pattern="*.png"):
    return len(list(folder.rglob(pattern))) if folder.exists() else 0


# ---------------------------------------------------------------- slides
def assays_slide(prs, blank, N, foot):
    """Assays of day N (added 2026-11-18, TST/HP days 5, 10, 13, 16, 18; tail immersion every day from day 5):
    analysis\\days\\dayN\\assays_dayN.csv from assays_day.py. Green = present (TST: the tracked blob is the mouse),
    yellow = TST scored after the camera-view fix, red = missing."""
    f = DAYS_DIR / f"day{N}" / f"assays_day{N}.csv"
    if not f.exists():
        return False
    A = pd.read_csv(f, dtype={"cage": str})
    s = prs.slides.add_slide(blank)
    MS.add_title(s, f"Assays day {N}: tail immersion · TST · hot plate",
                 "TI = the latency you sent (s) · TST = automatic, first 3 min, with a per-frame mouse check · "
                 "HOT = your hand score, first 30 s")
    hot_cols = [c for c in A.columns if c.startswith("hot_pct_")]
    names = [c.replace("hot_pct_", "").replace("_", " ") for c in hot_cols]
    header = ["Mouse", "TI (s)", "TST non-moving", "TST check"] + [f"HOT {n}" for n in names]
    rest = 12.33 - 1.45 - 0.8 - 1.2 - 2.1
    widths = [1.45, 0.8, 1.2, 2.1] + [rest / max(1, len(hot_cols))] * len(hot_cols)

    def cell(txt, ok=True, warn=False):
        st = "red" if not ok else ("yellow" if warn else "green")
        return (txt, QC_FILL[st], QC_TEXT[st], False)

    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return float("nan")

    rows = []
    for c, m in MICE:
        g = A[(A.cage == c) & (A.mouse == m)]
        r = g.iloc[0].to_dict() if len(g) else {}
        ti, tst, nm = num(r.get("tail_immersion_s")), num(r.get("tst_pct_nonmoving")), num(r.get("tst_not_mouse_pct"))
        run, src = str(r.get("tst_run", "")), str(r.get("tst_source", ""))
        fix = run == "viewfix"
        chk = ("hand score" if src == "hand" else
               f"{'camera-view fix' if fix else 'untouched'}; not mouse {nm:g} %" if nm == nm else run)
        row = [(f"{c} {m}", "FFFFFF", TEXT, True), cell(f"{ti:.2f}" if ti == ti else "missing", ti == ti),
               cell(f"{tst:.1f} %" if tst == tst else "missing", tst == tst, fix),
               cell(chk if tst == tst else "-", tst == tst, fix)]
        hot_ok = str(r.get("hot_file", "")) not in ("", "nan")
        for hc in hot_cols:
            v = num(r.get(hc))
            row.append(cell(f"{v:.1f} %" if hot_ok and v == v else ("missing" if not hot_ok else "-"), hot_ok))
        rows.append(row)
    table(s, 0.5, 1.45, widths, header, rows, "Assays table", size=10, head_size=10)
    MS.add_text(s, 0.5, 6.2, 12.3, 0.75,
                f"The original Step 02 reads them from assays\\day{N}\\<cage>\\<mouse>\\ (TI: *_immersion_day{N}.csv, "
                f"cell B1 · TST: ManualScoringResults_*_TST_day{N}_AUTO.mat · HOT: ManualScoringResults_*_HOT_day{N}_"
                f"first30s.mat) · TST examples: analysis\\tst\\day{N}\\examples\\ · table: analysis\\days\\day{N}\\"
                f"assays_day{N}.csv · yellow = the camera moved during that TST video; the moved frames were aligned "
                "before scoring", size=10, color=MUTED, name="Assays note")
    MS.add_footer(s, foot)
    return True


def build(days, out: Path):
    tag = f"days{days[0]}-{days[-1]}" if len(days) > 1 else f"day{days[0]}"
    qf, iff = DAYS_DIR / f"qc_summary_{tag}.csv", DAYS_DIR / f"qc_issues_{tag}.csv"
    QS.main(["--days", *map(str, days)])        # refresh the QC picture / tables first (about 10 s)
    Q = pd.read_csv(qf, dtype={"cage": str})
    I = pd.read_csv(iff, dtype={"cage": str}).fillna("")
    S = pd.read_csv(P / "analysis" / "session_checks" / "session_checks.csv", dtype={"cage": str}, encoding="utf-8-sig")
    S = S[S.cage.isin(QS.CAGES) & S.day.isin(days)].sort_values("checked_at").drop_duplicates(
        ["day", "cage", "mouse"], keep="last")
    PRT = pr_table(days)
    run = newest_run()
    cfd = cohort_dir(run)
    N, made = days[-1], datetime.now().strftime("%Y-%m-%d %H:%M")
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    blank = prs.slide_layouts[6]
    foot = f"Routine summary day {N} · made {made} · project: {P}"

    # ---- 1 overview
    s = prs.slides.add_slide(blank)
    MS.add_title(s, f"Day {N} routine summary",
                 f"Days {days[0]}–{N} · {len(MICE)} mice (9669m M, 9670f F) · QC first; the full analysis comes after day 18")
    n = len(Q)
    nr, ny = int((Q.status == "red").sum()), int((Q.status == "yellow").sum())
    rew_ok = int((S.reward_cmd == S.valve_openings).sum())
    has_vid = S.video.astype(str).str.strip().ne("") & S.video.notna()
    feq = int(((S.frames_minus_pulses == 0) & has_vid).sum())
    fneq = int(((S.frames_minus_pulses != 0) & S.frames_minus_pulses.notna() & has_vid).sum())
    cards = [(f"{n}", "sessions checked", f"days {days[0]}–{N}, 8 mice a day", NAVY),
             (f"{nr}", "red = needs a look", f"green {n - nr - ny} · yellow {ny} (known / fixed)", "B42318"),
             (f"{rew_ok}/{n}", "rewards Python = Saleae", "valve openings vs reward commands", "1E6B34"),
             (f"{feq}+{fneq}", "video frames = pulses", f"+{fneq} fixed · {int((~has_vid).sum())} without video", "557A9E")]
    for k, (big, lab, sub, col) in enumerate(cards):
        MS.stat_card(s, 0.5 + k * 3.13, 1.45, 3.0, 1.1, big, lab, sub, col, f"Card {k + 1}")
    # today
    lines = []
    prd = [d for d in days if d in PRT and PRT[d].pr_RL_full.notna().any()]
    if N in prd:
        X = PRT[N]
        lines.append([("PR (last requirement): ", {"bold": True}),
                      (f"median {X.pr_RL_full.median():g}, range {X.pr_RL_full.min():.0f}–{X.pr_RL_full.max():.0f}", {})])
        lines.append([("Rewards: ", {"bold": True}), (f"median {X.rewards_python.median():g}", {})])
        if len(prd) > 1:
            Y = PRT[prd[-2]]
            m = X.merge(Y[["cage", "color", "pr_RL_full"]], on=["cage", "color"], suffixes=("", "_prev"))
            m["dpr"] = m.pr_RL_full - m.pr_RL_full_prev
            big = m.loc[m.dpr.abs().idxmax()]
            lines.append([("vs day %d: " % prd[-2], {"bold": True}),
                          (f"PR median {Y.pr_RL_full.median():g} → {X.pr_RL_full.median():g}; biggest change "
                           f"{big.cage} {big.color} {big.pr_RL_full_prev:.0f} → {big.pr_RL_full:.0f}", {})])
        ref = MS.read_csv(P / "analysis" / "reference" / "prev_cohort_pr_days3to5.csv", comment="#")
        rr = ref[ref["day"] == N]["req_last_910s"].dropna() if len(ref) else pd.Series(dtype=float)
        if len(rr):
            lines.append([("Previous cohort day %d: " % N, {"bold": True}),
                          (f"median {rr.median():g}, range {rr.min():.0f}–{rr.max():.0f} (reference only)", {})])
    else:
        sub = S[S.day == N]
        lines.append([("Rewards: ", {"bold": True}), (f"median {sub.reward_cmd.median():g}", {})])
    qn = Q[Q.day == N]
    lines.append([("QC today: ", {"bold": True}),
                  (f"{len(qn)} sessions · green {int((qn.status == 'green').sum())}, yellow "
                   f"{int((qn.status == 'yellow').sum())}, red {int((qn.status == 'red').sum())}", {})])
    box(s, 0.5, 2.8, 6.1, 3.0, f"Today (day {N})", lines, NAVY, "Today box")
    dec = []
    for cat in ("off_eye", "open_issue", "rewards", "frames_open", "no_video_open", "lick_unknown", "check_open",
                "pupil_open"):
        g = I[(I.cat == cat) & (I.level == "red")]
        if len(g):
            prob = CAT[cat][0] or g.short.iloc[0]
            dec.append([(prob + ": ", {"bold": True}), (sess_text(g, cat), {})])
    NS = MS.read_csv(DAYS_DIR / "next_steps.csv")
    for r in NS[(NS.get("who", pd.Series(dtype=str)) == "you") & (NS.get("status", "") == "open")].itertuples() \
            if len(NS) else []:
        if r.what.startswith("Before day"):
            dec.append([("Soon: ", {"bold": True}), (r.what, {})])
    box(s, 6.75, 2.8, 6.08, 3.0, "Needs your decision", dec or ["Nothing open."], "B42318", "Decision box", size=11)
    has_assays = (DAYS_DIR / f"day{N}" / f"assays_day{N}.csv").exists()
    kA = 1 if has_assays else 0
    MS.add_text(s, 0.5, 6.0, 12.3, 0.4, "Slides: 2 QC per session · " + ("3 assays (TI, TST, HOT) · " if has_assays else "")
                + f"{3 + kA} data problems and fixes · {4 + kA} routine problems fixed · {5 + kA} where all plots are · "
                f"{6 + kA}–{11 + kA} basic plots · {12 + kA} what to do better next", size=12, color=MUTED,
                name="Contents")
    MS.add_footer(s, foot)

    # ---- 2 QC per session
    s = prs.slides.add_slide(blank)
    MS.add_title(s, "QC per session", "Each session vs the protocol, the Saleae and the camera")
    legend(s, 6.55, 1.0, [(QC_FILL["green"], QC_TEXT["green"], "all fine"),
                          (QC_FILL["yellow"], QC_TEXT["yellow"], "known / fixed / approved"),
                          (QC_FILL["red"], QC_TEXT["red"], "needs a look")])
    w0 = 1.55
    wd = (12.33 - w0) / len(days)
    big_pt = 11 if len(days) <= 6 else 8
    shp = s.shapes.add_table(len(MICE) + 1, len(days) + 1, Inches(0.5), Inches(1.45), Inches(12.33),
                             Inches(0.32 + 0.56 * len(MICE)))
    shp.name = "QC table"
    t = shp.table
    t.columns[0].width = Inches(w0)
    for j in range(len(days)):
        t.columns[j + 1].width = Inches(wd)
    MS.style_cell(t.cell(0, 0), [[("Mouse", {"bold": True, "color": "FFFFFF"})]], NAVY, size=11)
    for j, d in enumerate(days):
        task = "PR" if d in prd else "FR"
        MS.style_cell(t.cell(0, j + 1), [[(f"Day {d} · {task}", {"bold": True, "color": "FFFFFF"})]], NAVY, size=11,
                      align=PP_ALIGN.CENTER)
    for i, (c, m) in enumerate(MICE, start=1):
        t.rows[i].height = Inches(0.56)
        MS.style_cell(t.cell(i, 0), [[(f"{c} ", {"color": MUTED}), (m, {"bold": True})]], "FFFFFF", size=11)
        for j, d in enumerate(days):
            r = Q[(Q.day == d) & (Q.cage == c) & (Q.mouse == m)]
            if r.empty:
                MS.style_cell(t.cell(i, j + 1), "-", "FFFFFF", size=big_pt)
                continue
            r = r.iloc[0]
            st, fl = r["status"], str(r["flags"])        # r.flags would be the pandas Series' own flags
            l1 = ((f"R {int(r['rewards'])}" if pd.notna(r["rewards"]) else "R -")
                  + (f" · PR {int(r['pr'])}" if pd.notna(r["pr"]) else ""))
            l2 = fl if fl not in ("", "nan") else "all fine"
            MS.style_cell(t.cell(i, j + 1), [[(l1, {"bold": True, "color": QC_TEXT[st]})],
                                             [(l2, {"size": big_pt - 2, "color": TEXT})]],
                          QC_FILL[st], size=big_pt, align=PP_ALIGN.CENTER)
    MS.add_text(s, 0.5, 6.55, 12.3, 0.45,
                "R = rewards (Python = Saleae in every session) · PR = last ratio requirement (= rewards + 1 in a full "
                f"910-s session, the previous cohort's PR) · flags are explained on slide {3 + kA} · the full picture with "
                f"every note: analysis\\days\\qc_summary_{tag}.png", size=10, color=MUTED, name="QC note")
    MS.add_footer(s, foot)

    # ---- 2b assays of the day (TI, TST, HOT), when there are any
    if has_assays:
        assays_slide(prs, blank, N, foot)

    # ---- 3 data problems and what was done
    s = prs.slides.add_slide(blank)
    MS.add_title(s, "Data problems found and what was done",
                 f"From the QC of days {days[0]}–{N}; one row per kind of problem (sessions: day and detail)")
    rows = []
    groups = []
    for cat, g in I[I.cat.isin(CAT)].groupby("cat", sort=False):
        prob, done, st = CAT[cat]
        if cat == "open_issue":
            for title, gg in g.groupby("short", sort=False):
                groups.append((st, title, sess_text(gg, cat), done))
            continue
        groups.append((st, prob, sess_text(g, cat), done))
    groups.sort(key=lambda x: ORDER.index(x[0]))
    for st, prob, sess, done in groups:
        fill, col, lab = STATUS[st]
        rows.append([(prob, "FFFFFF", TEXT, True), (sess, "FFFFFF", TEXT, False), (done, "FFFFFF", TEXT, False),
                     (lab, fill, col, True)])
    size = 10 if len(rows) <= 12 else 9
    h = table(s, 0.5, 1.42, [3.55, 4.3, 3.28, 1.2], ["Problem", "Sessions", "What was done", "Status"], rows,
              "Problems table", size=size)
    if 1.42 + h > 6.9:
        print(f"  !! slide 3 table may run off the slide (estimated {h:.1f} in)")
    MS.add_footer(s, f"Source: analysis\\days\\qc_issues_{tag}.csv (qc_summary.py) · {foot}")

    # ---- 4 routine problems fixed
    s = prs.slides.add_slide(blank)
    MS.add_title(s, "Routine problems fixed (rig, files, software)",
                 "Problems in the way the days ran or were analysed, and what was done (newest last)")
    RL = MS.read_csv(DAYS_DIR / "routine_log.csv")
    if len(RL):
        RL = RL[pd.to_numeric(RL.day, errors="coerce").isin(days)].tail(14)
        rows = []
        for r in RL.itertuples():
            fill, col, lab = STATUS.get(str(r.status), STATUS["noted"])
            rows.append([(str(r.date)[5:], "FFFFFF", MUTED, False), (str(r.day), "FFFFFF", MUTED, False),
                         (str(r.problem), "FFFFFF", TEXT, False), (str(r.done), "FFFFFF", TEXT, False),
                         (lab, fill, col, True)])
        table(s, 0.5, 1.42, [0.7, 0.45, 5.0, 4.98, 1.2], ["Date", "Day", "Problem", "What was done", "Status"], rows,
              "Routine table", size=9.5)
    MS.add_footer(s, f"Source: analysis\\days\\routine_log.csv · {foot}")

    # ---- 5 where all plots are
    s = prs.slides.add_slide(blank)
    MS.add_title(s, "Where to find all plots and tables", f"All paths are inside {P}")
    rn = run.name if run else "run_###"
    places = [
        ("QC picture (every session, every note)", f"analysis\\days\\qc_summary_{tag}.png  (+ .csv, qc_issues_{tag}.csv)",
         DAYS_DIR / f"qc_summary_{tag}.png"),
        ("Day figures: dashboard, rasters, licks, PR, pupil, eye-box check", "analysis\\days\\dayN\\figs\\",
         DAYS_DIR),
        ("PR table per day (valve, lock-out, licks, PR vs previous cohort)", "analysis\\days\\dayN\\pr_checks_dayN.txt",
         DAYS_DIR),
        ("Session checks (one report per session)", "analysis\\session_checks\\dayN\\*_check.txt",
         P / "analysis" / "session_checks"),
        ("Across days, per mouse", "analysis\\figs_across_days\\  ·  analysis\\master_day_level.csv",
         P / "analysis" / "figs_across_days"),
        ("Previous-cohort figures (PR01, TJ02, PU01/03/04/06, QC03)",
         f"analysis\\cohort_figs\\{rn}\\bout1s_excl4-6-11-14-17\\  (3 more day-rule sets next to it)",
         P / "analysis" / "cohort_figs" / rn),
        ("Step 02 table (input of every later analysis)",
         f"step02_base\\longitudinal_outputs\\{rn}\\ALL_mice_longitudinal.csv", None),
        ("Data index: every saved data file (which one is final)",
         "analysis\\data_index.csv  ·  analysis\\days\\dayN\\data_index_dayN.csv  (loader: tools\\data_index.py)", None),
        ("Pupil QC sheets (per session)",
         "pupil_tst_teach\\dayN\\<cage>\\<mouse>\\pupil_outputs_general\\qc_sheet.png  ·  pupil_tst_teach\\QC_index.html",
         P / "pupil_tst_teach"),
        ("Whole eye (9669m black, 9670f black)", "analysis\\wholeeye\\  (review sheets: qc\\review_20261117\\ · "
                                                 "old vs new eye model: qc\\compare_black2_black3\\)",
         P / "analysis" / "wholeeye"),
        ("One example frame per session", "analysis\\days\\summary_day1234_examples.png", None),
        ("Problem / next-step lists used in this deck", "analysis\\days\\routine_log.csv · next_steps.csv · "
                                                      "qc_open_issues.csv", None),
        ("These decks", "analysis\\ppt\\", P / "analysis" / "ppt"),
    ]
    patterns = {"Day figures": (DAYS_DIR, "day*/figs/*.png"), "PR table": (DAYS_DIR, "day*/pr_checks_day*.txt"),
                "Session checks": (P / "analysis" / "session_checks", "day*/*_check.txt"),
                "Across days": (P / "analysis" / "figs_across_days", "*.png"),
                "Previous-cohort": (P / "analysis" / "cohort_figs" / rn, "*/*.png"),
                "Pupil QC": (P / "pupil_tst_teach", "day*/*/*/pupil_outputs_general/qc_sheet.png"),
                "Whole eye": (P / "analysis" / "wholeeye", "**/*.png"),
                "These decks": (P / "analysis" / "ppt", "*.pptx")}
    rows = []
    for what, where, _ in places:
        spec = next((v for k, v in patterns.items() if what.startswith(k)), None)
        nfig = ""
        if spec is not None and spec[0].is_dir():
            nfig = f"{len(list(spec[0].glob(spec[1])))} {spec[1].rsplit('.', 1)[-1]}"
        rows.append([(what, "FFFFFF", TEXT, True), (where, "FFFFFF", TEXT, False), (nfig, "FFFFFF", MUTED, False)])
    table(s, 0.5, 1.42, [3.9, 7.33, 1.1], ["What", "Where (inside the project)", "Files"], rows, "Folders table",
          size=10)
    MS.add_footer(s, foot)

    # ---- 6-10 basic plots
    across = P / "analysis" / "figs_across_days" / "across_days_metrics.png"
    if across.exists():
        MS.image_slide(prs, across, "Basic plots 1 · across days, per mouse",
                       "First 15 min 10 s of each session · solid = female, dashed = male · PR breakpoint = last "
                       "requirement · pupil panel: 9669m black and 9670f black are the WHOLE EYE (eye damage), not the pupil")
    if cfd is not None and (cfd / "PR01_daily_PR_spaghetti.png").exists():
        note = ""
        if N in prd and len(prd) > 1:
            note = (f"day {N} median {PRT[N].pr_RL_full.median():g} (day {prd[-2]}: "
                    f"{PRT[prd[-2]].pr_RL_full.median():g}) · ")
        MS.image_slide(prs, cfd / "PR01_daily_PR_spaghetti.png", "Basic plots 2 · PR per day (previous cohort's PR01)",
                       note + "groups are not set yet, so every mouse shows as 'Active' (pairs before day 6) · "
                              f"from Step 02 {rn}")
    tj, pu6 = (cfd / "TJ02_lick_density_vs_PR_level.png", cfd / "PU06_reward_locked_pupil.png") if cfd else (None, None)
    if tj is not None and (tj.exists() or pu6.exists()):
        s = prs.slides.add_slide(blank)
        MS.add_text(s, 0.5, 0.25, 12.33, 0.5, "Basic plots 3 · licking vs PR level (TJ02) and pupil at reward (PU06)",
                    size=22, color=NAVY, bold=True, name="Title")
        MS.add_text(s, 0.5, 0.78, 12.33, 0.34, "Previous cohort's figures, PR days so far (Pre phase) · mean ± SEM over "
                    "mice · PU06 leaves out the 2 whole-eye mice", size=12, color=MUTED, name="Caption")
        if tj.exists():
            picture(s, tj, 0.5, 1.25, 6.0, 5.85, "TJ02")
        if pu6.exists():
            picture(s, pu6, 6.83, 1.25, 6.0, 5.85, "PU06")
        MS.add_text(s, 0.5, 7.18, 12.33, 0.25, str(cfd), size=8, color=MUTED, name="Source")
    if cfd is not None and (cfd / "PU01_PU03_pupil_trajectory.png").exists():
        MS.image_slide(prs, cfd / "PU01_PU03_pupil_trajectory.png", "Basic plots 4 · pupil across days (PU01 / PU03)",
                       "Sessions with > 50 % tracked frames, PR days · 6 mice with a visible pupil (9669m black and "
                       "9670f black = whole eye, left out)")
    figs = DAYS_DIR / f"day{N}" / "figs"
    two_up(prs, figs / "09_pr_breakpoint.png", figs / "02_cumulative_licks_rewards.png",
           f"Basic plots 5 · today (day {N}): PR over the session and cumulative licks / rewards",
           "Requirement completed over time and the breakpoint per mouse (top); cumulative licks and rewards (bottom) · "
           "first 15 min 10 s", stacked=True)
    two_up(prs, figs / "05_peri_reward_licking.png", figs / "11_peri_reward_pupil.png",
           f"Basic plots 6 · today (day {N}): licking and pupil around reward",
           "Mean ± SEM over rewards · pupil: 9669m black and 9670f black are the WHOLE EYE (eye damage), not the pupil")

    # ---- 11 what to do better next
    s = prs.slides.add_slide(blank)
    MS.add_title(s, "What to do better next", "From the problems above · open items only (done items stay in "
                                              "analysis\\days\\next_steps.csv)")
    for k, (who, title, col) in enumerate((("you", "At the rig (you)", NAVY), ("Claude", "In the routine (Claude)",
                                                                                 "557A9E"))):
        items = NS[(NS.who == who) & (NS.status != "done")] if len(NS) else NS
        paras = []
        for r in items.itertuples():
            paras.append([(f"{r.what}", {"bold": True, "size": 12})])
            paras.append([(f"why: {r.why}", {"size": 10, "color": MUTED})])
        box(s, 0.5 + k * 6.25, 1.42, 6.08, 5.45, title, paras or ["Nothing open."], col, f"Next {who}", size=12)
    MS.add_footer(s, f"Source: analysis\\days\\next_steps.csv · {foot}")

    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        prs.save(out)
    except PermissionError:                     # the deck is open in PowerPoint: save next to it
        out = out.with_name(out.stem + datetime.now().strftime("_%H%M") + out.suffix)
        prs.save(out)
    print(f"{out}\n  {len(prs.slides)} slides")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--days", type=int, nargs="+", help="default: every day in session_checks.csv")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    if not a.days:
        S = pd.read_csv(P / "analysis" / "session_checks" / "session_checks.csv", dtype={"cage": str},
                        encoding="utf-8-sig")
        a.days = sorted(S[S.cage.isin(QS.CAGES)].day.astype(int).unique())
    out = a.out or P / "analysis" / "ppt" / f"routine_day{a.days[-1]}_{datetime.now():%Y%m%d}.pptx"
    build(a.days, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
