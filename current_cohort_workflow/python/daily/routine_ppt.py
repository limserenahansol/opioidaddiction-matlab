#!/usr/bin/env python3
"""routine_ppt.py - the daily routine PowerPoint (user 2026-11-17: "PPT as routine summary ... QC, problem fixed,
folder address where I can find all plots and representative plots as fundamental, what to do better next";
complicated analysis waits until after day 18). Report only: reads what the routine made, changes nothing.

  1  overview: counts, today's PR and TI, decisions needed 2  PR score per mouse and day, all days (2026-11-18)
  2b tail immersion per mouse and day, all days + plot (2026-11-20; user: "in summary we should also show TI data
     aswell ofc in all days"; tools\\ti_table_days.py)
  3  QC per session (green / yellow / red)                 (+ on assay days: TI / TST / HOT of the day, 2026-11-18)
  4  pupil check per session: Claude's check and the frames out, + the overview picture (2026-11-18)
  5  data problems found and what was done                 6  routine problems fixed (rig, files, software)
  7  where to find all plots                                8-13  basic plots (previous-cohort figures + daily figures)
  14 what to do better next
  appendix: the pupil-check picture of every session of --pupil-days (default: the last day), one per slide
  (slide numbers move with the assay slide and the number of overview pictures; slide 1 lists them)

Sources (inside the project): analysis\\days\\qc_summary_days*.csv and qc_issues_days*.csv (tools\\qc_summary.py),
analysis\\days\\routine_log.csv and next_steps.csv (kept by Claude), analysis\\days\\dayN\\pr_checks_dayN.csv,
analysis\\reference\\prev_cohort_pr_days3to5.csv, analysis\\figs_across_days\\, analysis\\cohort_figs\\run_###\\,
analysis\\days\\dayN\\figs\\, analysis\\days\\pr_table_days*.csv (tools\\pr_table_days.py), analysis\\pupil_check\\
(tools\\pupil_check_examples.py; both are refreshed by this script, like the QC picture).
usage:  python routine_ppt.py [--days 1 2 3 4] [--pupil-days 1 2 3 4 5] [--out FILE.pptx]
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
import pr_table_days as PT  # noqa: E402  (PR score per mouse and day, 2026-11-18)
import pupil_check_examples as PCE  # noqa: E402  (pupil check pictures + Claude's verdicts, 2026-11-18)
import ti_table_days as TI  # noqa: E402  (tail immersion per mouse and day, 2026-11-20)

P = TOOLS.parent
DAYS_DIR = P / "analysis" / "days"
NAVY, TEXT, MUTED, LINE = MS.NAVY, MS.TEXT, MS.MUTED, MS.LINE
MICE = QS.MICE
TST_HP_DAYS = {5, 10, 13, 16, 18}  # TST + hot plate days (tail immersion every day from day 5); see assays_day.py
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
    "pupil_noisy": ("Pupil noisy: ≥ 10 % of frames out (spout / whiskers / reflections) or wrong spikes",
                    "Checked by eye on the pupil-check pictures; used in the tables, marked noisy", "noted"),
    "pupil_pending": ("Pupil step not run yet for this day (light run right after the experiment)",
                      "Runs next: eye boxes, the pupil package and the whole eye, then the pupil check", "noted"),
}
SHOW_SHORT = {"frames", "frames_open", "length", "off_eye", "off_eye_kept", "licks_extra", "lick_start", "lick_valve_close",
              "lick_late", "lick_unknown", "pupil_noisy", "rewards", "filename", "check_open"}


def off_eye_decisions_text():
    """'What was done' for kept off-eye frames, read from analysis/days/qc_decisions.csv (qc_summary.decisions()):
    every mouse, its day (or every day) and the date of the user's answer. Added 2026-11-19: the fixed text said
    '11-18' for every session, which was wrong once 9670f black day 6 was decided on 11-19. None = no decision."""
    parts = []
    for (c, m, cat), L in sorted(QS.decisions().items()):
        if cat != "off_eye":
            continue
        for day, _dec, when in sorted(L, key=lambda x: (x[0] is not None, x[0] or 0)):
            w = when[5:] if len(when) == 10 else when        # 2026-11-18 -> 11-18
            parts.append(f"{c} {m} {'every day' if day is None else f'day {day}'} ({w})")
    if not parts:
        return None
    return "Your decision: keep those frames as they are - " + "; ".join(parts) + " (qc_decisions.csv)"


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
                 "HOT = your hand score, first 30 s" +
                 ("" if N in TST_HP_DAYS else " · TST and hot plate only on days 5, 10, 13, 16, 18"))
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

    # 2026-11-19: TST is only on the TST/HP days; on the other days the TST cells say so (grey), not red "missing"
    tst_day = N in TST_HP_DAYS
    rows = []
    for c, m in MICE:
        g = A[(A.cage == c) & (A.mouse == m)]
        r = g.iloc[0].to_dict() if len(g) else {}
        ti, tst, nm = num(r.get("tail_immersion_s")), num(r.get("tst_pct_nonmoving")), num(r.get("tst_not_mouse_pct"))
        run, src = str(r.get("tst_run", "")), str(r.get("tst_source", ""))
        fix = run == "viewfix"
        chk = ("hand score" if src == "hand" else
               f"{'camera-view fix' if fix else 'untouched'}; not mouse {nm:g} %" if nm == nm else run)
        if tst_day or tst == tst:
            tst_cells = [cell(f"{tst:.1f} %" if tst == tst else "missing", tst == tst, fix),
                         cell(chk if tst == tst else "-", tst == tst, fix)]
        else:
            tst_cells = [("no TST today", "F2F4F7", MUTED, False), ("-", "F2F4F7", MUTED, False)]
        row = [(f"{c} {m}", "FFFFFF", TEXT, True), cell(f"{ti:.2f}" if ti == ti else "missing", ti == ti),
               *tst_cells]
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


def _blue(v):
    """fill for a PR cell: v = 0 (lowest PR of the table) .. 1 (highest)"""
    a, b = (0xEE, 0xF4, 0xFB), (0x8F, 0xBB, 0xE3)
    return "".join(f"{int(round(x + (y - x) * max(0.0, min(1.0, v)))):02X}" for x, y in zip(a, b))


def pr_slide(prs, blank, days, foot, Tb):
    """PR score of every mouse on every day so far (added 2026-11-18; user: "all days' PR score table show me them";
    the user's order after an experiment: PR and QC first). Tb = pr_table_days.table(days): PR = pr_RL_full of
    pr_checks_dayN.csv (last requirement = rewards + 1 in a full 910-s session); FR days = rewards in the first
    15 min 10 s; days 6-10: a passive mouse replays its active partner's log (no PR of its own; its rewards are
    compared with the partner's); rewards Python = Saleae per session; the previous cohort's PR (first 910 s)."""
    s = prs.slides.add_slide(blank)
    MS.add_title(s, "PR score per mouse and day",
                 "PR = last ratio requirement = rewards + 1 (full 910-s session; the previous cohort's PR) · FR days: "
                 "rewards in the first 15 min 10 s · P = passive: R = drops delivered; small line = its partner's R "
                 "and PR")
    tasks = {d: PT.day_type(Tb, d) for d in days}
    two = len(days) <= 8
    w0, w1 = 1.6, 0.95
    wd = (12.33 - w0 - w1) / len(days)
    nrow = len(MICE) + 3
    rh = 0.44 if two else 0.36
    shp = s.shapes.add_table(nrow, len(days) + 2, Inches(0.5), Inches(1.45), Inches(12.33), Inches(0.34 + rh * (nrow - 1)))
    shp.name = "PR table"
    t = shp.table
    t.columns[0].width, t.columns[1].width = Inches(w0), Inches(w1)
    for j in range(len(days)):
        t.columns[j + 2].width = Inches(wd)
    t.rows[0].height = Inches(0.34)
    for j, h in enumerate(["Mouse", "Group"] + [f"Day {d} · {tasks[d]}" for d in days]):
        MS.style_cell(t.cell(0, j), [[(h, {"bold": True, "color": "FFFFFF"})]], NAVY, size=10 if two else 8,
                      align=PP_ALIGN.CENTER if j >= 2 else PP_ALIGN.LEFT)
    prv = Tb.pr.dropna()
    lo, hi = (prv.min(), prv.max()) if len(prv) else (0, 1)
    big, small = (12, 8) if two else (9, 7)
    for i, (c, m) in enumerate(MICE, start=1):
        t.rows[i].height = Inches(rh)
        g = Tb[(Tb.cage == c) & (Tb.mouse == m)]
        MS.style_cell(t.cell(i, 0), [[(f"{c} ", {"color": MUTED}), (m, {"bold": True})]], "FFFFFF", size=10)
        grp = str(g.group.iloc[0]) if len(g) and str(g.group.iloc[0]) not in ("", "nan") else "-"
        MS.style_cell(t.cell(i, 1), grp, "FFFFFF", size=9, color=MUTED)
        for j, d in enumerate(days):
            r = g[g.day == d]
            cell = t.cell(i, j + 2)
            if r.empty or r.task.iloc[0] == "missing":
                MS.style_cell(cell, "-", "F2F2F2", size=big, align=PP_ALIGN.CENTER)
                continue
            r = r.iloc[0]
            flag = "" if r.rewards_match else "  Python ≠ Saleae!"
            if r.task == "PR":
                lines = [[(f"{r.pr:.0f}", {"bold": True})]]
                if two:
                    lines.append([(f"rewards {r.rewards_python:.0f}{flag}", {"size": small, "color": MUTED})])
                MS.style_cell(cell, lines, _blue((r.pr - lo) / max(1e-6, hi - lo)), size=big, align=PP_ALIGN.CENTER)
            elif r.task == "PASSIVE":
                same = r.replay_match is True or r.replay_match == 1
                lines = [[(f"P · R {r.rewards_python:.0f}", {})]]
                if two:
                    # 2026-11-19 (user: "네, 같이 적기"): the partner's R and PR, so P's R is not read as "PR - 1";
                    # the longest text that fits the column (Calibri 8 pt: about 0.42 em per character, measured on
                    # the LibreOffice render; 0.2 in for the cell margins)
                    name = str(getattr(r, "partner", "") or PT.partner_of(r.replay_source))
                    pc, _, pcol = name.partition(" ")
                    short = pcol if pc == c and pcol else name
                    n_src, ppr = r.replay_source_rewards, getattr(r, "partner_pr", float("nan"))
                    nums = (f"R {n_src:.0f}" if pd.notna(n_src) else "R ?") + (f" · PR {ppr:.0f}" if pd.notna(ppr)
                                                                                else "")
                    tail = ("" if same else " (≠ partner!)") + flag
                    opts = [f"partner {short}: {nums}{tail}", f"{short}: {nums}{tail}", f"{nums}{tail}"]
                    fit = [x for x in opts if len(x) * small * 0.42 / 72 <= wd - 0.2]
                    lines.append([(fit[0] if fit else opts[-1],
                                   {"size": small, "color": MUTED if same else "B42318"})])
                MS.style_cell(cell, lines, "ECE8F8" if same else QC_FILL["red"], size=big, align=PP_ALIGN.CENTER)
            else:
                lines = [[(f"R {r.rewards_910s:.0f}", {})]]
                if two:
                    sub = (f"untimed {r.session_s / 60:.0f} min: {r.rewards_python:.0f} in all" if r.session_s > 930
                           else r.task) + flag
                    lines.append([(sub, {"size": small, "color": "9A5B00" if r.session_s > 930 else MUTED})])
                MS.style_cell(cell, lines, "EEF1F5", size=big, align=PP_ALIGN.CENTER)
    for k, (lab, fill) in enumerate((("Median (range)", "FFF6DB"), ("Previous cohort (first 15:10)", "F6F7F9"))):
        i = len(MICE) + 1 + k
        t.rows[i].height = Inches(rh)
        MS.style_cell(t.cell(i, 0), [[(lab, {"bold": True})]], fill, size=9)
        MS.style_cell(t.cell(i, 1), "", fill, size=9)
        for j, d in enumerate(days):
            g = Tb[Tb.day == d]
            if k == 0:
                is_pr = tasks[d].startswith("PR")
                x = g.pr.dropna() if is_pr else g.rewards_910s.dropna()
                pre = "" if is_pr else "R "
                act = " (active)" if tasks[d] == "PR + passive" else ""
                txt = [[(f"{pre}{x.median():g}", {"bold": True})], [(f"{x.min():.0f}–{x.max():.0f}{act}",
                                                                      {"size": small, "color": MUTED})]] if len(x) else "-"
            else:
                pm = g.prev_pr_median.dropna()
                txt = ([[(f"{float(pm.iloc[0]):g}", {})], [(f"{float(g.prev_pr_min.dropna().iloc[0]):.0f}–"
                                                            f"{float(g.prev_pr_max.dropna().iloc[0]):.0f}",
                                                            {"size": small, "color": MUTED})]]
                       if len(pm) else "-")
            MS.style_cell(t.cell(i, j + 2), txt if two or txt == "-" else [txt[0]], fill, size=big,
                          align=PP_ALIGN.CENTER)
    done = Tb[Tb.task != "missing"]
    n_ok = int(done.rewards_match.astype(bool).sum())
    pas = done[done.task == "PASSIVE"]
    p_txt = (f" · passive rewards = the partner's in {int((pas.replay_match == True).sum())}/{len(pas)}"  # noqa: E712
             if len(pas) else "")
    MS.add_text(s, 0.5, 6.42, 12.3, 0.5,
                f"Rewards Python = Saleae in {n_ok}/{len(done)} sessions{p_txt} · previous cohort = reference only "
                "(training and valve differ) · table and picture: analysis\\days\\pr_table_" +
                (f"days{days[0]}-{days[-1]}" if len(days) > 1 else f"day{days[0]}") + ".csv / .png",
                size=10, color=MUTED, name="PR note")
    MS.add_footer(s, foot)


def _orange(v):
    """fill for a TI cell: v = 0 (shortest latency of the table) .. 1 (longest)"""
    a, b = (0xFF, 0xF4, 0xE8), (0xF6, 0xB2, 0x6B)
    return "".join(f"{int(round(x + (y - x) * max(0.0, min(1.0, v)))):02X}" for x, y in zip(a, b))


def ti_slide(prs, blank, TIt, foot):
    """Tail-immersion latency of every mouse on every TI day so far (added 2026-11-20; user: "in summary we should
    also show TI data aswell ofc in all days"). TIt = ti_table_days.table(): the latency the user sent (s) and its
    change from day 5 (the last natural-reward day before morphine). A native table, and the plot
    analysis\\days\\ti_plot_days*.png (ti_table_days.py) next to it, or on a second slide after 8 TI days."""
    days = sorted(int(d) for d in TIt.day.unique())
    tag = TI.tag_of(days)
    s = prs.slides.add_slide(blank)
    MS.add_title(s, "Tail immersion per mouse and day",
                 "Latency you sent (s) · small line = change from day 5 (before morphine) · Baseline 1–5 natural reward"
                 " · Yoked 6–10 morphine · Equal Access 11–13 · Withdrawal 14–16 · Re-exposure 17–18")
    side = len(days) <= 8                       # the plot fits next to the table
    w0, w1 = 1.45, 0.8
    wd = 0.92 if side else (12.33 - w0 - w1) / len(days)
    W = w0 + w1 + wd * len(days)
    nrow = len(MICE) + 3                        # header + mice + Active / Passive median
    rh = 0.42
    shp = s.shapes.add_table(nrow, len(days) + 2, Inches(0.5), Inches(1.45), Inches(W), Inches(0.34 + rh * (nrow - 1)))
    shp.name = "TI table"
    t = shp.table
    t.columns[0].width, t.columns[1].width = Inches(w0), Inches(w1)
    for j in range(len(days)):
        t.columns[j + 2].width = Inches(wd)
    t.rows[0].height = Inches(0.34)
    for j, h in enumerate(["Mouse", "Group"] + [f"Day {d}" for d in days]):
        MS.style_cell(t.cell(0, j), [[(h, {"bold": True, "color": "FFFFFF"})]], NAVY, size=10,
                      align=PP_ALIGN.CENTER if j >= 2 else PP_ALIGN.LEFT)
    v = TIt.latency_s.dropna()
    lo, hi = (v.min(), v.max()) if len(v) else (0, 1)
    small = 7 if wd < 1.0 else 8
    for i, (c, m) in enumerate(MICE, start=1):
        t.rows[i].height = Inches(rh)
        g = TIt[(TIt.cage == c) & (TIt.mouse == m)]
        MS.style_cell(t.cell(i, 0), [[(f"{c} ", {"color": MUTED}), (m, {"bold": True})]], "FFFFFF", size=10)
        grp = str(g.group.iloc[0]) if len(g) and str(g.group.iloc[0]) not in ("", "nan") else "-"
        MS.style_cell(t.cell(i, 1), grp, "FFFFFF", size=9, color=MUTED)
        for j, d in enumerate(days):
            r = g[g.day == d]
            cell = t.cell(i, j + 2)
            if r.empty or pd.isna(r.latency_s.iloc[0]):
                sent = bool(len(r)) and bool(r.session.iloc[0])     # a session but no TI value yet
                MS.style_cell(cell, "not sent" if sent else "-", QC_FILL["yellow"] if sent else "F2F2F2",
                              size=8 if sent else 11, color=QC_TEXT["yellow"] if sent else TEXT, align=PP_ALIGN.CENTER)
                continue
            r = r.iloc[0]
            lines = [[(f"{r.latency_s:.2f}", {"bold": True})]]
            if d > TI.FIRST_TI_DAY and pd.notna(r.change_vs_day5_s):
                lines.append([(f"{r.change_vs_day5_s:+.2f} vs d{TI.FIRST_TI_DAY}", {"size": small, "color": MUTED})])
            MS.style_cell(cell, lines, _orange((r.latency_s - lo) / max(1e-6, hi - lo)), size=11, align=PP_ALIGN.CENTER)
    for k, grp in enumerate(("Active", "Passive")):
        i = len(MICE) + 1 + k
        t.rows[i].height = Inches(rh)
        MS.style_cell(t.cell(i, 0), [[(f"{grp} median", {"bold": True})]], "FFF6DB", size=9)
        MS.style_cell(t.cell(i, 1), "", "FFF6DB", size=9)
        for j, d in enumerate(days):
            x = TIt[(TIt.day == d) & (TIt.group == grp)].latency_s.dropna()
            txt = ([[(f"{x.median():.2f}", {"bold": True})], [(f"n {len(x)}", {"size": small, "color": MUTED})]]
                   if len(x) else "-")
            MS.style_cell(t.cell(i, j + 2), txt, "FFF6DB", size=11, align=PP_ALIGN.CENTER)
    plot = DAYS_DIR / f"ti_plot_{tag}.png"
    if side and plot.exists():
        x0 = 0.5 + W + 0.25
        picture(s, plot, x0, 1.4, 12.83 - x0, 4.85, "TI plot")
    MS.add_text(s, 0.5, 6.42, 12.3, 0.5,
                "Values as you sent them in chat (assays\\tail_immersion_all.csv); Step 02 reads the same value (B1 of "
                "assays\\dayN\\<cage>\\<mouse>\\*_immersion_dayN.csv) · '-' = no session that day yet · 'not sent' = a "
                f"session but no TI value · table and picture: analysis\\days\\ti_table_{tag}.csv / .png",
                size=10, color=MUTED, name="TI note")
    MS.add_footer(s, foot)
    if not side and plot.exists():
        s2 = prs.slides.add_slide(blank)
        MS.add_title(s2, "Tail immersion per day: plot", "One thin line per mouse, thick = group median · shaded = "
                     "morphine days · Active red, Passive blue")
        picture(s2, plot, 0.5, 1.45, 12.33, 4.9, "TI plot")
        MS.add_footer(s2, foot)


def pupil_check_slides(prs, blank, days, foot):
    """Claude's check of the pupil the tables use, every session so far (added 2026-11-18; user: "show me png files
    as example of 8 mice x 5 days ... your confirmation this session is okay if it is noisy how many frames are
    out"): tools\\pupil_check_examples.py -> analysis\\pupil_check\\ (one picture per session, the verdicts Claude
    writes in pupil_check_verdicts.csv after looking, overview pictures). A session without a verdict is red."""
    PCE.main(["--days", *map(str, days)])
    tag = f"days{days[0]}-{days[-1]}" if len(days) > 1 else f"day{days[0]}"
    C = pd.read_csv(PCE.OUT / f"pupil_check_{tag}.csv", dtype={"cage": str}, encoding="utf-8-sig")
    C["verdict"] = C.verdict.fillna("").astype(str)
    s = prs.slides.add_slide(blank)
    MS.add_title(s, "Pupil check per session (Claude looked at every picture)",
                 "ok = on the pupil and < 10 % out · noisy = ≥ 10 % out or wrong spikes")
    legend(s, 7.3, 1.0, [(QC_FILL["green"], QC_TEXT["green"], "ok"), (QC_FILL["yellow"], QC_TEXT["yellow"], "noisy"),
                         ("ECECEC", "555555", "no video"), (QC_FILL["red"], QC_TEXT["red"], "not looked at yet")])
    w0 = 1.55
    wd = (12.33 - w0) / len(days)
    pt = 11 if len(days) <= 6 else 8
    rh = 0.5
    shp = s.shapes.add_table(len(MICE) + 1, len(days) + 1, Inches(0.5), Inches(1.45), Inches(12.33),
                             Inches(0.32 + rh * len(MICE)))
    shp.name = "Pupil check table"
    t = shp.table
    t.columns[0].width = Inches(w0)
    for j in range(len(days)):
        t.columns[j + 1].width = Inches(wd)
    MS.style_cell(t.cell(0, 0), [[("Mouse", {"bold": True, "color": "FFFFFF"})]], NAVY, size=11)
    for j, d in enumerate(days):
        MS.style_cell(t.cell(0, j + 1), [[(f"Day {d}", {"bold": True, "color": "FFFFFF"})]], NAVY, size=11,
                      align=PP_ALIGN.CENTER)
    for i, (c, m) in enumerate(MICE, start=1):
        t.rows[i].height = Inches(rh)
        MS.style_cell(t.cell(i, 0), [[(f"{c} ", {"color": MUTED}), (m, {"bold": True})]], "FFFFFF", size=11)
        for j, d in enumerate(days):
            r = C[(C.day == d) & (C.cage == c) & (C.mouse == m)]
            cell = t.cell(i, j + 1)
            if r.empty:
                MS.style_cell(cell, "-", "FFFFFF", size=pt, align=PP_ALIGN.CENTER)
                continue
            r = r.iloc[0]
            if pd.isna(r.get("out_pct")):
                novid = "no video" in str(r.get("source", ""))
                MS.style_cell(cell, [[("no video" if novid else "not run yet", {"bold": True, "color": "555555"})],
                                     [("behaviour only" if novid else "pupil step", {"size": pt - 2, "color": MUTED})]],
                              "ECECEC", size=pt, align=PP_ALIGN.CENTER)
                continue
            v = r.verdict.lower()
            st = {"ok": "green", "noisy": "yellow"}.get(v, "red")
            eye = " · eye" if str(r.grade) == "whole_eye" else ""
            MS.style_cell(cell, [[((v.upper() if v else "LOOK"), {"bold": True, "color": QC_TEXT[st]})],
                                 [(f"out {r.out_pct:.1f} %{eye}", {"size": pt - 2, "color": TEXT})]],
                          QC_FILL[st], size=pt, align=PP_ALIGN.CENTER)
    vc = C[C.out_pct.notna()].verdict.str.lower().value_counts()
    nv = int(C.out_pct.isna().sum())
    MS.add_text(s, 0.5, 6.28, 12.3, 0.62,
                [[(f"ok {int(vc.get('ok', 0))} · noisy {int(vc.get('noisy', 0))} · no video {nv}"
                   + (f" · NOT LOOKED AT {int((C.out_pct.notna() & (C.verdict == '')).sum())}"
                      if (C.out_pct.notna() & (C.verdict == "")).any() else ""), {"bold": True})],
                 "out = frames the table fills from the neighbours: held (no pupil / refused jump) + jumps + odd size "
                 "· whole eye (eye): no eye + off eye · pictures: analysis\\pupil_check\\dayN\\ · numbers: "
                 f"pupil_check_{tag}.csv"], size=10, color=MUTED, name="Pupil check note")
    MS.add_footer(s, foot)
    for f in PCE.overview_files(days):
        if f.exists():
            MS.image_slide(prs, f, "Pupil check · one frame per session",
                           "green circle = the pupil in the tables (9670f green2 day 5 = dark-contrast repair) · pink = "
                           "whole eye (eye damage) · frame = the spread frame nearest the session's median size")
    return C


def pupil_appendix(prs, days, C):
    """one slide per session picture of the given days (the user's 40 examples of 2026-11-18; default the last day)"""
    for d in days:
        for c, m in MICE:
            f = PCE.OUT / f"day{d}" / f"day{d}_{c}_{m}.png"
            if not f.exists():
                continue
            r = C[(C.day == d) & (C.cage == c) & (C.mouse == m)] if C is not None else pd.DataFrame()
            v = str(r.verdict.iloc[0]).upper() if len(r) and str(r.verdict.iloc[0]) else ""
            note = str(r.note.iloc[0]) if len(r) and str(r.note.iloc[0]) not in ("", "nan") else ""
            title = f"Pupil check · day {d} · {c} {m}" + (f" · {v}" if v else "")
            if len(r) and pd.notna(r.out_pct.iloc[0]):
                note = f"out {r.out_pct.iloc[0]:.1f} % of frames · " + note
            MS.image_slide(prs, f, title, note or "no pupil (see the picture)")


def build(days, out: Path, pupil_days=None):
    tag = f"days{days[0]}-{days[-1]}" if len(days) > 1 else f"day{days[0]}"
    qf, iff = DAYS_DIR / f"qc_summary_{tag}.csv", DAYS_DIR / f"qc_issues_{tag}.csv"
    QS.main(["--days", *map(str, days)])        # refresh the QC picture / tables first (about 10 s)
    Q = pd.read_csv(qf, dtype={"cage": str})
    I = pd.read_csv(iff, dtype={"cage": str}).fillna("")
    S = pd.read_csv(P / "analysis" / "session_checks" / "session_checks.csv", dtype={"cage": str}, encoding="utf-8-sig")
    S = S[S.cage.isin(QS.CAGES) & S.day.isin(days)].sort_values("checked_at").drop_duplicates(
        ["day", "cage", "mouse"], keep="last")
    PRT = pr_table(days)
    PT.main(["--days", *map(str, days)])        # refresh analysis\\days\\pr_table_days*.csv/.png (a few s)
    PTb = PT.table(days)                        # one row per mouse and day: PR / passive / FR (2026-11-18)
    ti_days = [d for d in days if d >= TI.FIRST_TI_DAY]    # tail immersion from day 5 (2026-11-20)
    TIt = TI.table(ti_days) if ti_days else pd.DataFrame()
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
        pas_n = PTb[(PTb.day == N) & (PTb.task == "PASSIVE")]
        act = " (active mice)" if len(pas_n) else ""
        lines.append([("PR (last requirement): ", {"bold": True}),
                      (f"median {X.pr_RL_full.median():g}, range {X.pr_RL_full.min():.0f}–{X.pr_RL_full.max():.0f}"
                       + act, {})])
        lines.append([("Rewards: ", {"bold": True}),
                      (f"median {X[X.pr_RL_full.notna()].rewards_python.median():g}{act}", {})])
        if len(pas_n):                          # days 6-10: the passive replays
            nm = int((pas_n.replay_match == True).sum())  # noqa: E712
            lines.append([("Passive: ", {"bold": True}),
                          (f"{len(pas_n)} replays; rewards = the partner's in {nm}/{len(pas_n)}", {})])
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
    tin = TIt[TIt.day == N] if len(TIt) else TIt
    if len(tin) and tin.latency_s.notna().any():   # 2026-11-20: today's tail immersion in the overview too
        parts = [f"{grp} {tin[tin.group == grp].latency_s.median():.2f}" for grp in ("Active", "Passive")
                 if tin[tin.group == grp].latency_s.notna().any()]
        b5 = TIt[TIt.day == TI.FIRST_TI_DAY].latency_s.dropna()
        lines.append([("Tail immersion: ", {"bold": True}),
                      (f"median {' · '.join(parts)} s ({int(tin.latency_s.notna().sum())} mice"
                       + (f"; day {TI.FIRST_TI_DAY} {b5.median():.2f} s" if len(b5) and N > TI.FIRST_TI_DAY else "")
                       + ")", {})])
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
    # the contents line is filled in at the end, from the slides actually made (2026-11-18: more slides, which move
    # with the assay slide, the overview pictures and the appendix)
    contents = MS.add_text(s, 0.5, 5.95, 12.3, 0.8, "Slides:", size=12, color=MUTED, name="Contents")
    MS.add_footer(s, foot)
    toc = []                                    # (first slide number, label) of every block of slides

    def block(label):
        toc.append((len(prs.slides) + 1, label))

    # ---- 2 PR score per mouse and day (all days so far; added 2026-11-18)
    block("PR per day")
    pr_slide(prs, blank, days, foot, PTb)

    # ---- 2b tail immersion per mouse and day (all TI days so far; added 2026-11-20, user: "in summary we should also
    # show TI data aswell ofc in all days")
    if len(TIt) and TIt.latency_s.notna().any():
        TI.main(["--days", *map(str, ti_days)])  # refresh analysis\\days\\ti_table_days*.csv/.png (a few s)
        block("Tail immersion per day")
        ti_slide(prs, blank, TIt, foot)

    # ---- 3 QC per session
    block("QC per session")
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
        task = PT.day_type(PTb, d, short=True) or ("PR" if d in prd else "FR")
        task = "FR" if task.startswith("FR") else task
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
            is_p = ((PTb.day == d) & (PTb.cage == c) & (PTb.mouse == m) & (PTb.task == "PASSIVE")).any()
            l1 = ((f"R {int(r['rewards'])}" if pd.notna(r["rewards"]) else "R -")
                  + (f" · PR {int(r['pr'])}" if pd.notna(r["pr"]) else "") + (" · passive" if is_p else ""))
            l2 = fl if fl not in ("", "nan") else "all fine"
            MS.style_cell(t.cell(i, j + 1), [[(l1, {"bold": True, "color": QC_TEXT[st]})],
                                             [(l2, {"size": big_pt - 2, "color": TEXT})]],
                          QC_FILL[st], size=big_pt, align=PP_ALIGN.CENTER)
    MS.add_text(s, 0.5, 6.55, 12.3, 0.45,
                "R = rewards (Python = Saleae in every session) · PR = last ratio requirement (= rewards + 1 in a full "
                "910-s session, the previous cohort's PR) · flags are explained on the 'Data problems' slide · the full "
                f"picture with every note: analysis\\days\\qc_summary_{tag}.png", size=10, color=MUTED, name="QC note")
    MS.add_footer(s, foot)

    # ---- 3b assays of the day (TI, TST, HOT), when there are any
    if has_assays:
        block("assays (TI, TST, HOT)")
        assays_slide(prs, blank, N, foot)

    # ---- 4 pupil check per session + overview picture(s) (added 2026-11-18)
    block("pupil check")
    PC = pupil_check_slides(prs, blank, days, foot)

    # ---- 5 data problems and what was done
    block("data problems and fixes")
    s = prs.slides.add_slide(blank)
    MS.add_title(s, "Data problems found and what was done",
                 f"From the QC of days {days[0]}–{N}; one row per kind of problem (sessions: day and detail)")
    rows = []
    groups = []
    for cat, g in I[I.cat.isin(CAT)].groupby("cat", sort=False):
        prob, done, st = CAT[cat]
        if cat == "off_eye_kept":          # 2026-11-19: who and when from qc_decisions.csv, not a fixed date
            done = off_eye_decisions_text() or done
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

    # ---- 6 routine problems fixed
    block("routine problems fixed")
    s = prs.slides.add_slide(blank)
    MS.add_title(s, "Routine problems fixed (rig, files, software)",
                 "Problems in the way the days ran or were analysed, and what was done (newest last)")
    RL = MS.read_csv(DAYS_DIR / "routine_log.csv")
    if len(RL):
        RL = RL[pd.to_numeric(RL.day, errors="coerce").isin(days)]
        # newest rows that fit above the footer (2026-11-18: long rows ran off the slide); the same height estimate
        # as table()
        wid, n_fit = [0.7, 0.45, 5.0, 4.98, 1.2], 0
        for n in range(min(14, len(RL)), 0, -1):
            hs = [max(0.27, 0.07 + 0.165 * max(MS.text_lines(str(v), w, 9.5) for v, w in
                                                 zip((r.date, r.day, r.problem, r.done, r.status), wid)) * 0.95)
                  for r in RL.tail(n).itertuples()]
            if 0.32 + sum(hs) <= 5.4:
                n_fit = n
                break
        RL = RL.tail(max(1, n_fit))
        rows = []
        for r in RL.itertuples():
            fill, col, lab = STATUS.get(str(r.status), STATUS["noted"])
            rows.append([(str(r.date)[5:], "FFFFFF", MUTED, False), (str(r.day), "FFFFFF", MUTED, False),
                         (str(r.problem), "FFFFFF", TEXT, False), (str(r.done), "FFFFFF", TEXT, False),
                         (lab, fill, col, True)])
        table(s, 0.5, 1.42, [0.7, 0.45, 5.0, 4.98, 1.2], ["Date", "Day", "Problem", "What was done", "Status"], rows,
              "Routine table", size=9.5)
    MS.add_footer(s, f"Source: analysis\\days\\routine_log.csv · {foot}")

    # ---- 7 where all plots are
    block("where all plots are")
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
        ("PR score of every mouse on every day (one table)", f"analysis\\days\\pr_table_{tag}.csv / .png", None),
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
        ("Pupil check per session (Claude's check, frames out)", "analysis\\pupil_check\\dayN\\dayN_<cage>_<mouse>.png"
                                                                 f"  ·  pupil_check_{tag}.csv  ·  overview pictures",
         P / "analysis" / "pupil_check"),
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
                "Pupil check": (P / "analysis" / "pupil_check", "day*/*.png"),
                "These decks": (P / "analysis" / "ppt", "*.pptx")}
    rows = []
    for what, where, _ in places:
        spec = next((v for k, v in patterns.items() if what.startswith(k)), None)
        nfig = ""
        if spec is not None and spec[0].is_dir():
            nfig = f"{len(list(spec[0].glob(spec[1])))} {spec[1].rsplit('.', 1)[-1]}"
        rows.append([(what, "FFFFFF", TEXT, True), (where, "FFFFFF", TEXT, False), (nfig, "FFFFFF", MUTED, False)])
    table(s, 0.5, 1.42, [3.9, 7.33, 1.1], ["What", "Where (inside the project)", "Files"], rows, "Folders table",
          size=10 if len(rows) <= 13 else 9.5)
    MS.add_footer(s, foot)

    # ---- 8-13 basic plots
    block("basic plots")
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
        man = MS.read_csv(P / "cohort_manifest.csv")
        grp_set = len(man) and "group" in man.columns and man.group.astype(str).str.strip().isin(
            ["Active", "Passive"]).any()
        gtxt = ("Active / Passive = the groups in cohort_manifest.csv (your pairs 11-16, groups 11-18) · " if grp_set
                else "groups are not set yet, so every mouse shows as 'Active' · ")
        MS.image_slide(prs, cfd / "PR01_daily_PR_spaghetti.png", "Basic plots 2 · PR per day (previous cohort's PR01)",
                       note + gtxt + f"from Step 02 {rn}")
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

    # ---- 14 what to do better next
    block("what to do better next")
    s = prs.slides.add_slide(blank)
    MS.add_title(s, "What to do better next", "From the problems above · open items only (done items stay in "
                                              "analysis\\days\\next_steps.csv)")
    for k, (who, title, col) in enumerate((("you", "At the rig (you)", NAVY), ("Claude", "In the routine (Claude)",
                                                                                 "557A9E"))):
        items = NS[(NS.who == who) & (NS.status != "done")] if len(NS) else NS
        paras = []
        # smaller text when there are many open items, so the box holds them (2026-11-18: 9 items ran off the slide)
        chars = sum(len(str(r.what)) + len(str(r.why)) for r in items.itertuples())
        big = 12 if chars <= 900 else 10.5 if chars <= 1300 else 9.5 if chars <= 1700 else 8.5 if chars <= 2200 else 8
        for r in items.itertuples():
            paras.append([(f"{r.what}", {"bold": True, "size": big})])
            paras.append([(f"why: {r.why}", {"size": big - 2, "color": MUTED})])
        box(s, 0.5 + k * 6.25, 1.42, 6.08, 5.45, title, paras or ["Nothing open."], col, f"Next {who}", size=12)
    MS.add_footer(s, f"Source: analysis\\days\\next_steps.csv · {foot}")

    # ---- appendix: the pupil-check picture of every session of pupil_days (default: the last day)
    pdays = [d for d in (pupil_days or [N]) if d in days]
    if pdays:
        block(f"pupil pictures day{'s' if len(pdays) > 1 else ''} {pdays[0]}" +
              (f"–{pdays[-1]}" if len(pdays) > 1 else "") + " (appendix)")
        pupil_appendix(prs, pdays, PC)

    # the contents line of slide 1
    parts = []
    for k, (start, label) in enumerate(toc):
        end = (toc[k + 1][0] - 1) if k + 1 < len(toc) else len(prs.slides)
        if end >= start:
            parts.append((f"{start}" if end == start else f"{start}–{end}") + " " + label)
    contents.text_frame.paragraphs[0].runs[0].text = "Slides: " + " · ".join(parts)

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
    ap.add_argument("--pupil-days", type=int, nargs="+",
                    help="days whose per-session pupil-check pictures go in the appendix (default: the last day)")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    if not a.days:
        S = pd.read_csv(P / "analysis" / "session_checks" / "session_checks.csv", dtype={"cage": str},
                        encoding="utf-8-sig")
        a.days = sorted(S[S.cage.isin(QS.CAGES)].day.astype(int).unique())
    out = a.out or P / "analysis" / "ppt" / f"routine_day{a.days[-1]}_{datetime.now():%Y%m%d}.pptx"
    build(a.days, out, a.pupil_days)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
