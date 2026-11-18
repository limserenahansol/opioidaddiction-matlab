#!/usr/bin/env python3
"""
make_summary_ppt.py - a 3-slide PowerPoint summary of chosen days (report only: it reads the project's tables and
changes nothing). Added 2026-11-16.

  slide 1  rewards, licks per minute and PR breakpoint per mouse and day
  slide 2  lick microstructure (bouts) and pupil (grade, median diameter, reward response) per mouse and day
  slide 3  the daily routine checklist: every step of the end-of-day routine and what it found on each day

All behaviour numbers are quick_analysis's first 15 min 10 s of each session, counted only while the Python task ran
(licks / valve openings before Python started or after it ended are dropped; the camera is triggered only during the
run, so video frames = camera pulses of the run).

Sources (all inside the project): analysis\\session_metrics.csv, analysis\\session_checks\\session_checks.csv,
sync\\sync_summary.csv, pupil_tst_final\\pupil_session_quality.csv, pupil_tst_final\\alignment_fixes.csv,
pupil_tst_final\\fullpupil_sessions.csv, video\\roi_log.csv, analysis\\days\\dayN\\pr_checks_dayN.txt,
analysis\\reference\\prev_cohort_pr_days3to5.csv, step02_base\\dayN, tools\\protocol_expected.csv.

usage:  python make_summary_ppt.py --days 1 2 3 [--project FOLDER] [--out FILE.pptx]
needs:  python-pptx (pip install python-pptx)
"""
from __future__ import annotations

import argparse
import json
import math
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from pptx import Presentation
from pptx.chart.data import CategoryChartData
from pptx.dml.color import RGBColor
from pptx.enum.chart import XL_CHART_TYPE, XL_LABEL_POSITION, XL_LEGEND_POSITION
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt

TOOLS = Path(__file__).resolve().parent
PROJECT_DEFAULT = TOOLS.parent
LOG_FOLDER = Path.home() / "Downloads" / "behavior_task_windows"     # where the task writes its logs (empty-log count)

# ---------------------------------------------------------------- look
FONT = "Calibri"
NAVY, TEXT, MUTED, GRID, LINE = "16263A", "2B3440", "6B7785", "E3E7EC", "C9D1DA"
FR_COLORS = ["A9BCD0", "557A9E", "2F4B6E"]                  # training days (FR): blue-greys
PR_COLORS = ["E4572E", "F0883E", "B23A48", "F2B134", "8C2F39"]  # PR days: warm
STATUS = {  # fill, text, label
    "pass": ("E6F4EA", "1E6B34", "PASS"),
    "check": ("FFF0D9", "8A5200", "CHECK (reported, look at it)"),
    "fixed": ("ECE8F8", "4B3A8C", "approved fix / option applied"),
    "na": ("F2F4F6", "6B7785", "not applicable / not run"),
    "info": ("FFFFFF", TEXT, "info"),
}
GRADE_STYLE = {  # fill, text
    "ok": ("E6F4EA", "1E6B34"), "noisy": ("FFF6DB", "7A5A00"), "guess": ("EEF0F3", "5F6B78"),
    "full_pupil": ("E3EEFA", "1F4E79"), "no_pupil": ("EEF0F3", "5F6B78"), "check": ("FFF0D9", "8A5200"),
    "best_effort": ("ECE8F8", "4B3A8C"),
    "no video": ("F7F7F8", "8A939E"), "-": ("FFFFFF", MUTED),
}
CAGE_SEX = {"9669m": "M", "9670f": "F"}
CHECK_GROUPS = "A files · B name · C protocol · D valve/rewards · E licks · F camera/video · G Saleae settings · H sync"


def rgb(h):
    return RGBColor.from_string(h)


# ---------------------------------------------------------------- data
def read_csv(p: Path, **kw) -> pd.DataFrame:
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p, encoding="utf-8-sig", **kw)


def log_lockout(project: Path, py) -> float:
    """total_lockout_duration written in the task log (session_checks may hold the pre-OneDrive path)."""
    if not isinstance(py, str) or not py:
        return np.nan
    p = Path(py)
    if not p.exists() and "python" in p.parts:
        p = project.joinpath(*p.parts[p.parts.index("python"):])
    if not p.exists():
        return np.nan
    m = re.search(r'"total_lockout_duration":\s*([0-9.]+)', p.read_text(encoding="utf-8", errors="replace"))
    return float(m.group(1)) if m else np.nan


def roster(project: Path):
    man = read_csv(project / "cohort_manifest.csv")
    if len(man):
        cages = sorted(man["cage_folder"].astype(str).unique())
        order = [(r.cage_folder, r.mouse) for r in man.sort_values(["cage_folder", "mouse"]).itertuples()]
    else:
        order = [(c, m) for c in ("9669m", "9670f") for m in ("black", "green1", "green2", "white")]
        cages = ["9669m", "9670f"]
    return cages, order


def load(project: Path, days):
    M = read_csv(project / "analysis" / "session_metrics.csv")
    M = M[(M["window"] == "first_15m10s") & (M["day"].isin(days))].copy()
    M["cage"] = M["mouse"].str.split("_").str[0]
    M["color"] = M["mouse"].str.split("_").str[1]
    # session_metrics.csv keeps 4 significant digits (10.549 -> 10.550 shows as 10.6); recompute from the counts
    ok = M["bout_count"].notna() & (M["duration_min"] > 0)
    M.loc[ok, "bouts_per_min"] = M.loc[ok, "bout_count"] / M.loc[ok, "duration_min"]
    SC = read_csv(project / "analysis" / "session_checks" / "session_checks.csv")
    if len(SC):
        SC = SC[SC["day"].isin(days)].sort_values("checked_at").drop_duplicates(["day", "session"], keep="last")
        SC["cage"] = SC["cage"].astype(str)
    SY = read_csv(project / "sync" / "sync_summary.csv")
    if len(SY):
        SY["day_n"] = SY["day"].str.extract(r"day(\d+)")[0].astype(int)
        SY = SY[SY["day_n"].isin(days)]
    PQ = read_csv(project / "pupil_tst_final" / "pupil_session_quality.csv")
    if len(PQ):
        PQ = PQ[PQ["day"].isin(days)]
        PQ["cage"] = PQ["cage"].astype(str)
    return M, SC, SY, PQ


def series_by_day(M, order, days, col):
    out = []
    for d in days:
        g = M[M["day"] == d].set_index(["cage", "color"])
        vals = [float(g.loc[k, col]) if k in g.index and pd.notna(g.loc[k, col]) else None for k in order]
        task = M.loc[M["day"] == d, "task"].dropna()
        out.append((d, task.iloc[0] if len(task) else "", vals))
    return out


def day_colors(series):
    fr = iter(FR_COLORS)
    pr = iter(PR_COLORS)
    return [next(pr, "E4572E") if str(task).upper().startswith("PR") else next(fr, "557A9E") for _, task, _ in series]


def is_pr_day(M, d):
    t = M.loc[M["day"] == d, "task"].dropna()
    return bool(len(t)) and str(t.iloc[0]).upper().startswith("PR")


# ---------------------------------------------------------------- pptx helpers
def add_text(slide, x, y, w, h, paragraphs, size=14, color=TEXT, bold=False, align=PP_ALIGN.LEFT,
             anchor=MSO_ANCHOR.TOP, name=None, margin=0.0, spacing_after=0):
    """paragraphs: str, or list of paragraphs; a paragraph is str or list of (text, dict(bold/color/size/italic))."""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    if name:
        tb.name = name
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    for side in ("margin_left", "margin_right", "margin_top", "margin_bottom"):
        setattr(tf, side, Inches(margin))
    if isinstance(paragraphs, str):
        paragraphs = [paragraphs]
    for i, para in enumerate(paragraphs):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(spacing_after)
        runs = [(para, {})] if isinstance(para, str) else para
        for text, o in runs:
            r = p.add_run()
            r.text = text
            f = r.font
            f.name = FONT
            f.size = Pt(o.get("size", size))
            f.bold = o.get("bold", bold)
            f.italic = o.get("italic", False)
            f.color.rgb = rgb(o.get("color", color))
    return tb


def add_title(slide, title, subtitle):
    add_text(slide, 0.5, 0.32, 12.33, 0.62, title, size=30, color=NAVY, bold=True, name="Title")
    add_text(slide, 0.5, 0.94, 12.33, 0.36, subtitle, size=13, color=MUTED, name="Subtitle")


def add_footer(slide, text, y=7.0, h=0.4):
    add_text(slide, 0.5, y, 12.33, h, text, size=10, color=MUTED, name="Footer")


def column_chart(slide, x, y, w, h, title, order, series, colors, number_format="0", label_size=8,
                 vmax=None, legend=True, name=None):
    cd = CategoryChartData(number_format=number_format)
    cages = []
    for cage, _ in order:
        if cage not in cages:
            cages.append(cage)
    for cage in cages:
        c = cd.add_category(f"{cage} ({CAGE_SEX.get(cage, '?')})")
        for cg, color in order:
            if cg == cage:
                c.add_sub_category(color)
    for d, task, vals in series:
        cd.add_series(f"Day {d} · {task}", vals)
    gf = slide.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(x), Inches(y), Inches(w), Inches(h), cd)
    if name:
        gf.name = name
    ch = gf.chart
    ch.font.name = FONT
    ch.font.size = Pt(9)
    ch.font.color.rgb = rgb(TEXT)
    ch.has_title = True
    ch.chart_title.text_frame.text = title
    tr = ch.chart_title.text_frame.paragraphs[0].runs[0].font
    tr.size, tr.bold, tr.name = Pt(12), True, FONT
    tr.color.rgb = rgb(TEXT)
    ch.has_legend = legend and len(series) > 1
    if ch.has_legend:
        ch.legend.position = XL_LEGEND_POSITION.TOP
        ch.legend.include_in_layout = False
        ch.legend.font.size = Pt(9)
    plot = ch.plots[0]
    plot.gap_width = 55
    plot.overlap = -8
    plot.has_data_labels = True
    dl = plot.data_labels
    dl.font.size = Pt(label_size)
    dl.font.color.rgb = rgb(TEXT)
    dl.number_format = number_format
    dl.number_format_is_linked = False
    dl.position = XL_LABEL_POSITION.OUTSIDE_END
    for s, col in zip(plot.series, colors):
        s.format.fill.solid()
        s.format.fill.fore_color.rgb = rgb(col)
        s.format.line.fill.background()
        s.invert_if_negative = False
    va = ch.value_axis
    va.has_major_gridlines = True
    va.major_gridlines.format.line.color.rgb = rgb(GRID)
    va.major_gridlines.format.line.width = Pt(0.75)
    va.format.line.fill.background()
    va.tick_labels.font.size = Pt(9)
    va.tick_labels.font.color.rgb = rgb(MUTED)
    va.minimum_scale = 0
    if vmax:
        va.maximum_scale = vmax
    ca = ch.category_axis
    ca.has_major_gridlines = False
    ca.format.line.color.rgb = rgb(LINE)
    ca.tick_labels.font.size = Pt(9)
    ca.tick_labels.font.color.rgb = rgb(TEXT)
    return gf


def nice_max(vals, pad=1.18):
    v = [x for s in vals for x in s if x is not None]
    if not v:
        return None
    m = max(v) * pad
    step = 10 ** np.floor(np.log10(m)) / 2 if m > 0 else 1
    return float(np.ceil(m / step) * step)


def stat_card(slide, x, y, w, h, big, label, sub, accent, name):
    card = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    card.name = name
    card.adjustments[0] = 0.08
    card.fill.solid()
    card.fill.fore_color.rgb = rgb("F6F7F9")
    card.line.fill.background()
    card.shadow.inherit = False
    add_text(slide, x + 0.18, y + 0.08, 1.5, h - 0.16, big, size=36, color=accent, bold=True,
             anchor=MSO_ANCHOR.MIDDLE, name=name + " value")
    add_text(slide, x + 1.62, y + 0.1, w - 1.75, h - 0.2,
             [[(label, {"bold": True, "size": 13, "color": TEXT})], [(sub, {"size": 11, "color": MUTED})]],
             anchor=MSO_ANCHOR.MIDDLE, name=name + " label", spacing_after=2)


def style_cell(cell, lines, fill, size=10, color=TEXT, bold=False, align=PP_ALIGN.LEFT):
    cell.fill.solid()
    cell.fill.fore_color.rgb = rgb(fill)
    cell.margin_left = cell.margin_right = Inches(0.06)
    cell.margin_top = cell.margin_bottom = Inches(0.02)
    cell.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf = cell.text_frame
    tf.word_wrap = True
    if isinstance(lines, str):
        lines = [[(lines, {})]]
    for i, line in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        if isinstance(line, str):
            line = [(line, {})]
        for t, o in line:
            r = p.add_run()
            r.text = t
            r.font.name = FONT
            r.font.size = Pt(o.get("size", size))
            r.font.bold = o.get("bold", bold)
            r.font.color.rgb = rgb(o.get("color", color))


def text_lines(text, width_in, size_pt):
    """Rough number of wrapped lines for Calibri (average glyph ~0.5 em)."""
    chars_per_line = max(8, int((width_in - 0.14) * 72 / (size_pt * 0.45)))
    return max(1, math.ceil(len(text) / chars_per_line))


# ---------------------------------------------------------------- routine checklist
def rng(vals, fmt="{:.0f}"):
    v = [x for x in vals if pd.notna(x)]
    if not v:
        return "-"
    lo, hi = min(v), max(v)
    return fmt.format(lo) if fmt.format(lo) == fmt.format(hi) else f"{fmt.format(lo)}–{fmt.format(hi)}"


def checklist(project: Path, days, M, SC, SY, PQ, order):
    cages = {c for c, _ in order}
    sc_all = SC.copy()
    SC = SC[SC["cage"].isin(cages)] if len(SC) else SC
    expected = read_csv(TOOLS / "protocol_expected.csv")
    fixes = read_csv(project / "pupil_tst_final" / "alignment_fixes.csv")
    fullp = read_csv(project / "pupil_tst_final" / "fullpupil_sessions.csv")
    roi_log = read_csv(project / "video" / "roi_log.csv")
    ref = read_csv(project / "analysis" / "reference" / "prev_cohort_pr_days3to5.csv", comment="#")
    dates = {}
    for d in days:
        s = SC[SC["day"] == d]["start"].dropna().astype(str)
        dates[d] = s.str[:10].mode().iloc[0] if len(s) else ""

    def exp_row(d):
        if not len(expected):
            return None
        r = expected[(expected["from_day"] <= d) & (expected["to_day"] >= d)]
        r = r[r["phase"].astype(str).str.lower() != "passive"]
        return r.iloc[0] if len(r) else None

    rows = []

    def add(step, cells):
        rows.append((step, cells))

    # 1 files
    cells = {}
    for d in days:
        s = SC[SC["day"] == d]
        n = len(s)
        nvid = int(s["video"].notna().sum())
        missing_v = [f"{r.cage} {r.mouse}" for r in s.itertuples() if pd.isna(r.video)]
        date = dates[d]
        aside = sum(1 for p in (project / "_tests_and_aborted" / "python").glob("*.jsonl")
                    if datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d") == date)
        empty = sum(1 for p in LOG_FOLDER.glob("*.jsonl")
                    if p.stat().st_size == 0 and datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d") == date)
        txt = f"log + Saleae {n}/{len(order)} · video {nvid}/{n}"
        if missing_v:
            txt += f" ({', '.join(missing_v)}: none)"
        extra = [f"{aside} short/test runs" if aside else "", f"{empty} empty logs" if empty else ""]
        extra = [e for e in extra if e]
        if extra:
            txt += " · aside: " + ", ".join(extra)
        cells[d] = (txt, "check" if missing_v or n < len(order) else "pass")
    add("Files filed (log · Saleae · video)", cells)

    # 2 protocol
    cells = {}
    for d in days:
        s = SC[SC["day"] == d]
        task = M.loc[M["day"] == d, "task"].dropna()
        task = task.iloc[0] if len(task) else ""
        byl = s.groupby(s["lockout_median_s"].round(2))
        parts = []
        for lk, g in byl:
            parts.append(f"{lk:.2f} s" + (f" ({g.cage.iloc[0]} {g.mouse.iloc[0]})" if len(g) == 1 and byl.ngroups > 1
                                          else (f" ({len(g)})" if byl.ngroups > 1 else "")))
        long = [f"{r.cage} {r.mouse} {r.session_s:.0f} s" for r in s.itertuples()
                if pd.notna(r.session_s) and r.session_s > 920]
        req = " 1+1" if str(task).upper().startswith("PR") else ""
        txt = f"{task}{req} · lock-out {' / '.join(parts)} (measured, reward -> next trial) · 910 s"
        logl = sorted({round(v, 2) for v in (log_lockout(project, r.python) for r in s.itertuples()) if pd.notna(v)})
        differs = any(abs(v - lk) > 0.05 for v in logl for lk in byl.groups)
        if differs:
            txt += f" · the log says {' / '.join(f'{v:.2f}' for v in logl)} s"
        if long:
            txt += f" · longer: {', '.join(long)}"
        bad = s["check_items"].fillna("").str.contains(r"(?:^|; )C ")
        cells[d] = (txt, "check" if bad.any() else ("info" if differs else "pass"))
    add("Protocol in the task log", cells)

    # 3 valve
    cells = {}
    for d in days:
        s = SC[SC["day"] == d]
        v = s["valve_median_ms"].dropna()
        groups = v.round(-1).value_counts().sort_index()
        e = exp_row(d)
        exp_ms = e["valve_ms"] if e is not None else np.nan
        if len(groups) > 1:
            txt = " · ".join(f"{int(k)} ms ×{n}" for k, n in groups.items()) + " (changed mid-day)"
            st = "check"
        else:
            txt = f"{rng(v, '{:.1f}')} ms ×{len(v)}"
            st = "pass" if pd.notna(exp_ms) and abs(v.median() - exp_ms) <= 20 else "info"
        if pd.notna(exp_ms):
            txt += f" (expected {exp_ms:.0f})"
        cells[d] = (txt, st)
    add("Valve open time (Saleae)", cells)

    # 4 rewards
    cells = {}
    for d in days:
        s = SC[SC["day"] == d]
        eq = int((s["reward_cmd"] == s["valve_openings"]).sum())
        cells[d] = (f"{eq}/{len(s)} sessions equal · {rng(s['reward_cmd'])} rewards per session",
                    "pass" if eq == len(s) else "check")
    add("Rewards: Python = Saleae", cells)

    # 5 licks
    cells = {}
    for d in days:
        s = SC[SC["day"] == d]
        y = SY[SY["day_n"] == d] if len(SY) else SY
        sal_only = (y["licks_saleae"] - y["licks_matched"]).clip(lower=0) if len(y) else pd.Series(dtype=float)
        val = s[s["python_only_valid"].fillna(0) > 0]
        flagged = s["check_items"].fillna("").str.contains("E Python-only")
        txt = f"Python-only {rng(s['python_only_licks'])}/session · Saleae-only {rng(sal_only)}"
        if len(val):
            txt += " · valid, counted: " + ", ".join(f"{r.cage} {r.mouse} {int(r.python_only_valid)}"
                                                    for r in val.itertuples())
        cells[d] = (txt, "check" if flagged.any() else "pass")
    add("Licks: Python vs Saleae", cells)

    # 6 frames
    cells = {}
    for d in days:
        s = SC[SC["day"] == d]
        v = s[s["video"].notna()]
        eq = int((v["frames_minus_pulses"] == 0).sum())
        diff = [f"{r.cage} {r.mouse} {int(r.frames_minus_pulses):+d}" for r in v.itertuples()
                if pd.notna(r.frames_minus_pulses) and r.frames_minus_pulses != 0]
        f = fixes[fixes["day"] == d] if len(fixes) else fixes
        txt = f"{eq}/{len(v)} equal"
        st = "pass"
        if diff:
            txt += f" · {', '.join(diff)}"
            st = "fixed" if len(f) >= len(diff) else "check"
        if len(v) < len(s):
            txt += f" · {len(s) - len(v)} session without video"
            st = "check" if st == "pass" else st
        cells[d] = (txt, st)
    add("Video frames = camera pulses", cells)

    # 7 sync
    cells = {}
    for d in days:
        y = SY[SY["day_n"] == d] if len(SY) else SY
        n_ok = int((y["status"] == "OK").sum())
        warn = y[y["status"] != "OK"]
        why = []
        for r in warn.itertuples():
            pr = str(r.problems)
            why.append("no video file" if "no video file" in pr else "video length" if "frames" in pr else
                       "lick matching" if "licks" in pr else "other")
        txt = f"{n_ok} OK"
        if len(warn):
            txt += f" · {len(warn)} WARN ({', '.join(sorted(set(why)))})"
        cells[d] = (txt, "pass" if not len(warn) else "check")
    add("Sync: one clock (Python · Saleae · video)", cells)

    # 8 session check
    cells = {}
    for d in days:
        s = SC[SC["day"] == d]
        n_pass = int((s["overall"] == "PASS").sum())
        bad = []
        for r in s[s["overall"] != "PASS"].itertuples():
            letters = sorted({it.strip()[0] for it in str(r.check_items).split(";") if it.strip()})
            bad.append(f"{r.cage} {r.mouse} ({' '.join(letters)})")
        txt = f"{n_pass} PASS · {len(s) - n_pass} CHECK" + (": " + ", ".join(bad) if bad else "")
        cells[d] = (txt, "pass" if n_pass == len(s) else "check")
    add("Session check A–H", cells)

    # 9 PR checks
    cells = {}
    for d in days:
        if not is_pr_day(M, d):
            cells[d] = ("— (FR training day)", "na")
            continue
        p = project / "analysis" / "days" / f"day{d}" / f"pr_checks_day{d}.txt"
        n_check = None
        if p.exists():
            m = re.search(r"^CHECK \((\d+)\)", p.read_text(encoding="utf-8", errors="replace"), re.M)
            n_check = int(m.group(1)) if m else 0
        pr = M[M["day"] == d]["pr_breakpoint"]
        txt = f"PR {rng(pr)} (median {pr.median():g})"
        if len(ref) and "req_last_910s" in ref:
            rr = ref[ref["day"] == d]["req_last_910s"].dropna()
            if len(rr):
                txt += f" · previous cohort day {d}: median {rr.median():g}"
        if n_check is not None:
            txt += f" · {n_check} CHECK" + (" (the valid Python-only licks)" if n_check else "")
        cells[d] = (txt, "check" if n_check else "pass")
    add("PR checks (per mouse)", cells)

    # 10 analysis + figures
    cells = {}
    for d in days:
        figs = sorted((project / "analysis" / "days" / f"day{d}" / "figs").glob("*.png"))
        summ = (project / "analysis" / "days" / f"day{d}" / f"day{d}_summary.txt").exists()
        cells[d] = (f"{len(figs)} figures + day{d}_summary.txt" if summ else f"{len(figs)} figures",
                    "pass" if figs and summ else "check")
    add("Daily analysis + figures", cells)

    # 11 ROI
    cells = {}
    for d in days:
        s = SC[SC["day"] == d]
        nvid = int(s["video"].notna().sum())
        r = roi_log[roi_log["day"] == f"day{d}"] if len(roi_log) else roi_log
        if not len(r):
            cells[d] = ("package ROI, no change logged", "info")
            continue
        r = r.assign(date=r["when"].str[:10])
        parts = []
        for i, (date, g) in enumerate(sorted(r.groupby("date"), key=lambda z: z[0])):
            n = g[["cage", "mouse"]].drop_duplicates().shape[0]
            if i == 0 and n >= nvid:
                parts.append(f"set {date[5:]}")
            else:
                parts.append(f"re-checked {date[5:]}: {n} re-centred" + (f", {nvid - n} kept" if i == 0 else ""))
        txt = " · ".join(parts)
        if r["how"].astype(str).str.contains("outside or at the edge").any():
            txt += " (the day-2 boxes missed the eye)"
        cells[d] = (txt, "pass")
    add("Eye box (ROI) check", cells)

    # 12 pupil tracked + grade
    cells = {}
    for d in days:
        q = PQ[PQ["day"] == d] if len(PQ) else PQ
        s = SC[SC["day"] == d]
        nvid = int(s["video"].notna().sum())
        if not len(q):
            cells[d] = (f"not run yet (0/{nvid})", "check")
            continue
        g = q["grade"].value_counts()
        order_g = ["ok", "noisy", "best_effort", "full_pupil", "guess", "no_pupil", "check"]
        txt = f"{len(q)}/{nvid} tracked · " + " · ".join(f"{g[k]} {k.replace('_', ' ')}" for k in order_g if k in g)
        cells[d] = (txt, "check" if ("check" in g or len(q) < nvid) else "pass")
    add("Pupil: tracked · grade", cells)

    # 13 alignment
    cells = {}
    for d in days:
        q = PQ[PQ["day"] == d] if len(PQ) else PQ
        if not len(q):
            cells[d] = ("—", "na")
            continue
        a = q["align_check"].fillna("").astype(str)
        fixed = q[a.str.contains("sync-aligned")]
        mism = q[a.str.upper().str.startswith("MISMATCH")]
        txt = f"{int((a.str.startswith('ok') & ~a.str.contains('sync-aligned')).sum())} ok"
        st = "pass"
        if len(fixed):
            fc = fixed["cage"].iloc[0]
            txt += f" · {len(fixed)} fixed: {fc} " + ", ".join(fixed["mouse"].astype(str)) + " (sync-aligned)"
            st = "fixed"
        if len(mism):
            txt += f" · {len(mism)} MISMATCH (ask first)"
            st = "check"
        cells[d] = (txt, st)
    add("Alignment: package table vs sync", cells)

    # 14 full pupil
    cells = {}
    for d in days:
        f = fullp[fullp["day"] == d] if len(fullp) else fullp
        q = PQ[(PQ["day"] == d) & (PQ["grade"] == "full_pupil")] if len(PQ) else PQ
        appr = f[f["approved"].astype(str).str.lower().str.startswith(("yes", "20"))] if len(f) else f
        withdrawn = f[f["approved"].astype(str).str.lower().str.startswith("no")] if len(f) else f
        names = sorted(set(q["mouse"].astype(str)) | (set(appr["mouse"].astype(str)) if len(appr) else set()))
        if names:
            cage = q["cage"].iloc[0] if len(q) else str(appr["cage"].iloc[0])
            cells[d] = (f"{len(names)}: {cage} {'/'.join(names)} (whole dark eye)", "fixed")
        elif len(withdrawn):
            cells[d] = (f"withdrawn 11-16 ({withdrawn['cage'].iloc[0]} {'/'.join(withdrawn['mouse'].astype(str))}): "
                        "the black pupil inside is the pupil -> package pupil", "check")
        else:
            cells[d] = ("none", "na")
    add("Full-pupil option", cells)

    # 14b best effort (2026-11-16)
    cells = {}
    be_list = read_csv(project / "pupil_tst_final" / "besteffort_sessions.csv")
    for d in days:
        fs = sorted((project / "pupil_tst_final" / f"day{d}").glob("*/*/pupil_outputs_besteffort/besteffort_check.json"))
        used = be_list[(be_list["day"] == d) & be_list["approved"].astype(str).str.lower().str.startswith(("yes", "20"))] \
            if len(be_list) else be_list
        if len(used):
            cells[d] = (f"{len(used)} in use (approved 11-16): " +
                        ", ".join(f"{r.cage} {r.mouse}" for r in used.itertuples()), "fixed")
        elif fs:
            cells[d] = (f"{len(fs)} made, not in the tables: " + ", ".join(f"{f.parts[-4]} {f.parts[-3]}" for f in fs),
                        "info")
        else:
            cells[d] = ("none", "na")
    add("Best-effort pupil (hard eyes)", cells)

    # 15 step 02 input
    cells = {}
    runs = sorted((project / "step02_base" / "longitudinal_outputs").glob("run_*"))
    run_t = runs[-1].stat().st_mtime if runs else None
    for d in days:
        base = project / "step02_base" / f"day{d}"
        dirs = [p for p in base.glob("*/*") if p.is_dir() and any(p.glob("concat_out_*"))] if base.exists() else []
        if not dirs:
            cells[d] = ("not copied yet", "check")
            continue
        newest = max(c.stat().st_mtime for p in dirs for c in p.glob("concat_out_*"))
        nvid = int(SC[SC["day"] == d]["video"].notna().sum()) if len(SC) else len(dirs)
        nses = int((SC["day"] == d).sum()) if len(SC) else len(dirs)
        nbo = sum(1 for p in dirs if any(c.name.endswith("_behavioronly") for c in p.glob("concat_out_*")))
        txt = f"{len(dirs)}/{nses} sessions copied" + (f" ({nbo} behavior only: no video)" if nbo else "")
        rebuilt = bool(runs and run_t and run_t >= newest)
        txt += f" · in {runs[-1].name}" if rebuilt else " · Step 02 table not rebuilt yet"
        cells[d] = (txt, "pass" if (len(dirs) >= nvid and rebuilt) else "check")
    add("Step 02 input (step02_base)", cells)
    return rows, dates


# ---------------------------------------------------------------- slides
def build(project: Path, days, out: Path, with_appendix: bool = True):
    cages, order = roster(project)
    M, SC, SY, PQ = load(project, days)
    made = datetime.now().strftime("%Y-%m-%d %H:%M")
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    blank = prs.slide_layouts[6]
    dtxt = f"Days {days[0]}–{days[-1]}" if len(days) > 1 else f"Day {days[0]}"
    window = "first 15 min 10 s of each session · counted only while the Python task ran (camera on)"

    # ---- slide 1
    s1 = prs.slides.add_slide(blank)
    add_title(s1, f"{dtxt} · rewards, licking and progressive ratio",
              f"{len(order)} mice ({', '.join(f'{c} {CAGE_SEX.get(c, chr(63))}' for c in cages)}) · {window}")
    ser = series_by_day(M, order, days, "rewards")
    cols = day_colors(ser)
    column_chart(s1, 0.45, 1.4, 7.95, 2.8, "Rewards per session", order, ser, cols,
                 vmax=nice_max([v for _, _, v in ser]), name="Rewards chart")
    ser = series_by_day(M, order, days, "licks_per_min")
    column_chart(s1, 0.45, 4.22, 7.95, 2.7, "Licks per minute", order, ser, cols,
                 vmax=nice_max([v for _, _, v in ser]), legend=False, name="Lick rate chart")
    pr_days = [d for d in days if is_pr_day(M, d)]
    if pr_days:
        ser_pr = series_by_day(M, order, pr_days, "pr_breakpoint")
        cols_pr = [c for (d, _, _), c in zip(ser, cols) if d in pr_days]
        column_chart(s1, 8.6, 1.4, 4.28, 3.15, "PR breakpoint (last requirement reached)", order, ser_pr, cols_pr,
                     vmax=nice_max([v for _, _, v in ser_pr], 1.22), name="PR chart")
        last = pr_days[-1]
        v = M[M["day"] == last]["pr_breakpoint"].dropna()
        stat_card(s1, 8.6, 4.7, 4.28, 1.05, f"{v.median():g}", f"Day {last} PR, median",
                  f"range {v.min():.0f}–{v.max():.0f} · n = {len(v)}", PR_COLORS[0], "Card new cohort")
        ref = read_csv(project / "analysis" / "reference" / "prev_cohort_pr_days3to5.csv", comment="#")
        rr = ref[ref["day"] == last]["req_last_910s"].dropna() if len(ref) else pd.Series(dtype=float)
        if len(rr):
            stat_card(s1, 8.6, 5.87, 4.28, 1.05, f"{rr.median():g}", f"Previous cohort, day {last}",
                      f"range {rr.min():.0f}–{rr.max():.0f} · n = {len(rr)} · reference only", "557A9E",
                      "Card previous cohort")
    add_footer(s1, "Licks = Saleae lick onsets (Python logs a few more, see slide 4); rewards = valve openings "
                   "(= Python reward commands); PR = last requirement reached in 910 s (start 1, +1). Previous cohort "
                   f"(0910/0911/6099/6100): first 910 s; training and valve (180 ms) differ. "
                   f"Source: analysis\\session_metrics.csv · made {made}")

    # ---- slide 2
    s2 = prs.slides.add_slide(blank)
    add_title(s2, f"{dtxt} · lick microstructure and pupil",
              "Bout = licks ≤ 0.25 s apart (a gap ≥ 0.5 s ends it) · pupil from the camera-on frames only "
              "(video frame = camera pulse)")
    ser = series_by_day(M, order, days, "bouts_per_min")
    column_chart(s2, 0.45, 1.4, 6.45, 2.8, "Lick bouts per minute", order, ser, cols, number_format="0.0",
                 label_size=7, vmax=nice_max([v for _, _, v in ser]), name="Bout rate chart")
    ser = series_by_day(M, order, days, "mean_licks_per_bout")
    column_chart(s2, 0.45, 4.22, 6.45, 2.7, "Licks per bout (mean)", order, ser, cols, number_format="0.0",
                 label_size=7, vmax=nice_max([v for _, _, v in ser]), legend=False, name="Bout size chart")
    x0, y0, wt = 7.15, 1.42, 5.73
    add_text(s2, x0, y0, wt, 0.3, "Pupil per session: grade · median diameter · reward response",
             size=12, color=TEXT, bold=True, name="Pupil table title")
    tbl_shape = s2.shapes.add_table(len(order) + 1, len(days) + 1, Inches(x0), Inches(y0 + 0.36), Inches(wt),
                                    Inches(0.36 + 0.5 * len(order)))
    tbl_shape.name = "Pupil table"
    tbl = tbl_shape.table
    tbl.columns[0].width = Inches(1.25)
    for j in range(len(days)):
        tbl.columns[j + 1].width = Inches((wt - 1.25) / len(days))
    tbl.rows[0].height = Inches(0.36)
    style_cell(tbl.cell(0, 0), [[("Mouse", {"bold": True, "color": "FFFFFF"})]], NAVY, size=11)
    for j, d in enumerate(days):
        task = M.loc[M["day"] == d, "task"].dropna()
        style_cell(tbl.cell(0, j + 1), [[(f"Day {d} · {task.iloc[0] if len(task) else ''}",
                                          {"bold": True, "color": "FFFFFF"})]], NAVY, size=11, align=PP_ALIGN.CENTER)
    pq = PQ.set_index(["day", "cage", "mouse"]) if len(PQ) else None
    for i, (cage, color) in enumerate(order, start=1):
        tbl.rows[i].height = Inches(0.5)
        style_cell(tbl.cell(i, 0), [[(f"{cage} ", {"color": MUTED}), (color, {"bold": True})]], "FFFFFF", size=11)
        for j, d in enumerate(days):
            mrow = M[(M["day"] == d) & (M["cage"] == cage) & (M["color"] == color)]
            has_session = len(mrow) > 0
            s = SC[(SC["day"] == d) & (SC["cage"] == cage) & (SC["mouse"] == color)] if len(SC) else SC
            no_video = has_session and len(s) > 0 and s["video"].isna().all()
            key = (d, cage, color)
            if no_video:
                grade, l2 = "no video", "pupil absent"
            elif pq is not None and key in pq.index:
                q = pq.loc[key]
                q = q.iloc[-1] if isinstance(q, pd.DataFrame) else q
                grade = str(q["grade"])
                px = mrow["pupil_median_px"].iloc[0] if has_session else np.nan
                resp = mrow["pupil_reward_response_pct"].iloc[0] if has_session else np.nan
                be = project / "pupil_tst_final" / f"day{d}" / cage / color / "pupil_outputs_besteffort" / \
                    "besteffort_check.json"
                if grade in ("guess", "no_pupil", "check"):
                    l2 = "not used"
                    if be.exists():
                        l2 = f"best effort {json.loads(be.read_text(encoding='utf-8'))['median_diam_px']:.0f} px (not approved)"
                elif pd.notna(px):
                    l2 = f"{px:.0f} px · {resp:+.1f} %" if pd.notna(resp) else f"{px:.0f} px"
                else:
                    l2 = "-"
            elif has_session:
                grade, l2 = "-", "not tracked yet"
            else:
                grade, l2 = "-", "no session"
            fill, tcol = GRADE_STYLE.get(grade, ("FFFFFF", TEXT))
            style_cell(tbl.cell(i, j + 1), [[(grade.replace("_", " "), {"bold": True, "color": tcol})],
                                            [(l2, {"size": 10, "color": TEXT if grade != "no video" else MUTED})]],
                       fill, size=11, align=PP_ALIGN.CENTER)
    add_text(s2, x0, y0 + 0.36 + 0.36 + 0.5 * len(order) + 0.08, wt, 0.62,
             "Reward response = mean % change 0.5–2.5 s after each reward vs the 1 s before. "
             "best effort = hard eye tracked with the checked best-effort tracker (approved 11-16); guess = not used. "
             "Full pupil (whole dark eye) was withdrawn on 11-16: the black pupil inside is tracked.",
             size=10, color=MUTED, name="Pupil note")
    add_footer(s2, f"Source: analysis\\session_metrics.csv, pupil_tst_final\\pupil_session_quality.csv · made {made}")

    # ---- slide 3
    s3 = prs.slides.add_slide(blank)
    rows, dates = checklist(project, days, M, SC, SY, PQ, order)
    add_title(s3, f"Daily routine checklist · {dtxt.lower()}",
              "Every end-of-day step and what it found · counts only inside the Python run (camera on) · "
              "checks report only")
    x_t, y_t, w_t = 0.48, 1.4, 12.38
    wcol0 = 2.45
    wday = (w_t - wcol0) / len(days)
    cell_pt, head_h = 9.5, 0.3
    heights = []
    for step, cells in rows:
        n = max([text_lines(step, wcol0, 10)] + [text_lines(cells.get(d, ("-", ""))[0], wday, cell_pt) for d in days])
        heights.append(max(0.26, 0.05 + 0.15 * n))
    th = head_h + sum(heights)
    t = s3.shapes.add_table(len(rows) + 1, len(days) + 1, Inches(x_t), Inches(y_t), Inches(w_t), Inches(th))
    t.name = "Routine table"
    tb = t.table
    tb.columns[0].width = Inches(wcol0)
    for j in range(len(days)):
        tb.columns[j + 1].width = Inches(wday)
    tb.rows[0].height = Inches(head_h)
    style_cell(tb.cell(0, 0), [[("Routine step", {"bold": True, "color": "FFFFFF"})]], NAVY, size=11)
    for j, d in enumerate(days):
        task = M.loc[M["day"] == d, "task"].dropna()
        style_cell(tb.cell(0, j + 1), [[(f"Day {d} · {task.iloc[0] if len(task) else ''} · {dates.get(d, '')}",
                                         {"bold": True, "color": "FFFFFF"})]], NAVY, size=11)
    for i, ((step, cells), h) in enumerate(zip(rows, heights), start=1):
        tb.rows[i].height = Inches(h)
        style_cell(tb.cell(i, 0), [[(step, {"bold": True})]], "FFFFFF", size=10)
        for j, d in enumerate(days):
            txt, st = cells.get(d, ("-", "na"))
            fill, tcol, _ = STATUS[st]
            style_cell(tb.cell(i, j + 1), [[(txt, {"color": tcol if st in ("check", "na") else TEXT})]], fill,
                       size=cell_pt)
    lx, ly = 0.5, y_t + th + 0.08
    for st in ("pass", "check", "fixed", "na"):
        fill, tcol, label = STATUS[st]
        sw = s3.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(lx), Inches(ly + 0.04), Inches(0.24),
                                 Inches(0.18))
        sw.name = f"Legend {label}"
        sw.fill.solid()
        sw.fill.fore_color.rgb = rgb(fill)
        sw.line.color.rgb = rgb(LINE)
        sw.shadow.inherit = False
        wlab = 0.068 * len(label) + 0.1
        add_text(s3, lx + 0.3, ly, wlab, 0.26, label, size=10, color=MUTED, name=f"Legend text {label}")
        lx += 0.3 + wlab + 0.25
    add_footer(s3, f"Session-check groups: {CHECK_GROUPS}. Sources: session_checks.csv, sync_summary.csv, "
                   "pr_checks_dayN.txt, pupil_session_quality.csv, alignment_fixes.csv, fullpupil_sessions.csv, "
                   f"video\\roi_log.csv, step02_base · made {made}", y=6.98, h=0.42)

    # ---- today's extra results (2026-11-16) and the appendix with every figure
    slide_licks(prs, project, days, made)
    slide_dayrules(prs, project, made)
    slide_besteffort(prs, project, made)
    n_fig = appendix(prs, project, days) if with_appendix else 0
    out.parent.mkdir(parents=True, exist_ok=True)
    prs.save(out)
    print(f"slides: {len(prs.slides)} ({n_fig} figure slides)")
    return out, y_t + th


# ---------------------------------------------------------------- extra result slides + appendix (2026-11-16)
CAPTIONS = {
    "00_day_dashboard": "Day dashboard: summary bars, rasters, cumulative and per-minute licks / rewards, peri-reward licking, ILI, PR",
    "00_day_summary_table": "Day summary table per mouse (first 15 min 10 s)",
    "01_lick_reward_raster": "Lick and reward raster per mouse (counted vs lock-out licks)",
    "01_pupil_roi_check": "Eye-box (ROI) check: 3 frames per mouse with the box the pupil tracker used",
    "02_cumulative_licks_rewards": "Cumulative licks and rewards over the session",
    "03_rate_per_minute": "Licks and rewards per minute (1-min bins; the last bin is the partial 10 s)",
    "04_inter_lick_interval": "Inter-lick intervals (dotted 0.25 s = in-bout, dashed 0.5 s = bout end; < 80 ms = likely double contacts)",
    "05_peri_reward_licking": "Licking around reward (mean +- SEM over rewards)",
    "06_summary_first_15m10s": "Per-mouse summary metrics (first 15 min 10 s, the same window for every mouse)",
    "07_time_to_reward": "Time to reward per trial",
    "09_pr_breakpoint": "PR: requirement completed over time and the breakpoint per mouse",
    "10_pupil_timeline": "Pupil diameter over the session (sessions whose grade is used)",
    "11_peri_reward_pupil": "Pupil around reward (% change; sessions whose grade is used)",
    "12_pr_checks": "PR daily checks: valve time, lock-out, Python-only licks, session end, PR vs the previous cohort",
    "across_days_metrics": "Across days: per-mouse metrics (solid = female, dashed = male)",
    "day3_roi_day2_boxes_on_day3": "Day-2 eye boxes (yellow) drawn on day-3 frames: every box missed part of the eye",
    "day3_roi_new_boxes": "Day-3 eye boxes set by eye (magenta) vs day 2 (yellow dashed): 5 time points + a 25-frame median",
    "qc_sheet": "Pupil QC sheet of the package (yellow = tracked pupil, green = U-Net mask)",
    "besteffort_trace": "Best-effort pupil (red) vs the package (grey); lower panel = % of frames with a visible pupil",
    "pupil_examples": "The pupil the tables use, 6 frames per session (green circle = tracked pupil); for your corrections",
    "fullpupil_sheet": "NOT USED - withdrawn 2026-11-16: the red outline (whole dark area) is not the pupil; you "
                       "showed the black pupil inside it, which the package tracks",
}


def image_slide(prs, path: Path, title: str, caption: str, section_note: str = ""):
    """One figure per slide, aspect kept. A tall figure (height > 1.3 x width) is shown in two overlapping halves
    (top 0-55 %, bottom 45-100 %) on two slides so it stays readable."""
    import io
    from PIL import Image
    with Image.open(path) as im:
        im.load()
        w, h = im.size
        parts = [(0, h, "")]
        if h > 1.3 * w:
            parts = [(0, int(0.55 * h), " (top)"), (int(0.45 * h), h, " (bottom)")]
        crops = []
        for y0, y1, tag in parts:
            if (y0, y1) == (0, h):
                crops.append((None, w, h, tag))
            else:
                buf = io.BytesIO()
                im.crop((0, y0, w, y1)).save(buf, format="PNG")
                buf.seek(0)
                crops.append((buf, w, y1 - y0, tag))
    for buf, cw, ch, tag in crops:
        s = prs.slides.add_slide(prs.slide_layouts[6])
        add_text(s, 0.5, 0.25, 12.33, 0.5, title + tag, size=22, color=NAVY, bold=True, name="Title")
        add_text(s, 0.5, 0.78, 12.33, 0.34, caption + (f"  ·  {section_note}" if section_note else ""), size=12,
                 color=MUTED, name="Caption")
        bx, by, bw, bh = 0.5, 1.2, 12.33, 6.0
        sc = min(bw / cw, bh / ch)
        pw, ph = cw * sc, ch * sc
        pic = s.shapes.add_picture(buf if buf is not None else str(path), Inches(bx + (bw - pw) / 2),
                                   Inches(by + (bh - ph) / 2), Inches(pw), Inches(ph))
        pic.name = path.stem + tag.replace(" ", "_")
        add_text(s, 0.5, 7.18, 12.33, 0.25, str(path), size=8, color=MUTED, name="Source")
    return len(crops)


def divider_slide(prs, title, sub):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    bg = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, prs.slide_width, prs.slide_height)
    bg.fill.solid()
    bg.fill.fore_color.rgb = rgb(NAVY)
    bg.line.fill.background()
    bg.name = "Background"
    add_text(s, 0.8, 2.7, 11.7, 1.0, title, size=40, color="FFFFFF", bold=True, name="Title")
    add_text(s, 0.8, 3.75, 11.7, 1.2, sub, size=16, color="C9D1DA", name="Subtitle")


def slide_licks(prs, project: Path, days, made):
    p = project / "analysis" / "checks" / "lick_mismatch_per_session.csv"
    if not p.exists():
        return
    T = pd.read_csv(p)
    T = T[T["day"].isin(days)].copy()
    T["cage"] = T["mouse"].str.split("_").str[0]
    T["color"] = T["mouse"].str.split("_").str[1]
    T["valve"] = T["py_only_during_valve_pad30"]
    T["start"] = T["py_only_at_start"]
    T["other"] = (T["python_only"] - T["valve"] - T["start"]).clip(lower=0)
    s = prs.slides.add_slide(prs.slide_layouts[6])
    add_title(s, "Why Python logs a few more licks than the Saleae",
              "Python-only licks per session, by cause (counted inside the Python run); the tables use the Saleae lick line")
    cd = CategoryChartData(number_format="0")
    for d in days:
        c = cd.add_category(f"Day {d}")
        for r in T[T["day"] == d].sort_values(["cage", "color"]).itertuples():
            c.add_sub_category(f"{r.cage[-1].upper()} {r.color}")
    rows = [r for d in days for r in T[T["day"] == d].sort_values(["cage", "color"]).itertuples()]
    cd.add_series("valve open (Arduino re-trigger)", [float(r.valve) for r in rows])
    cd.add_series("before Python start (serial queue)", [float(r.start) for r in rows])
    cd.add_series("other", [float(r.other) for r in rows])
    gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_STACKED, Inches(0.45), Inches(1.4), Inches(8.0), Inches(5.5), cd)
    gf.name = "Python-only licks chart"
    ch = gf.chart
    ch.font.name, ch.font.size = FONT, Pt(9)
    ch.has_title = True
    ch.chart_title.text_frame.text = "Python-only licks per session (m = 9669m, f = 9670f)"
    f = ch.chart_title.text_frame.paragraphs[0].runs[0].font
    f.size, f.bold, f.name = Pt(12), True, FONT
    ch.has_legend = True
    ch.legend.position = XL_LEGEND_POSITION.TOP
    ch.legend.include_in_layout = False
    ch.legend.font.size = Pt(9)
    plot = ch.plots[0]
    plot.gap_width = 40
    for sr, col in zip(plot.series, ["557A9E", "E4572E", "C9D1DA"]):
        sr.format.fill.solid()
        sr.format.fill.fore_color.rgb = rgb(col)
    va = ch.value_axis
    va.has_major_gridlines = True
    va.major_gridlines.format.line.color.rgb = rgb(GRID)
    va.format.line.fill.background()
    va.tick_labels.font.size = Pt(9)
    va.tick_labels.font.color.rgb = rgb(MUTED)
    ca = ch.category_axis
    ca.tick_labels.font.size = Pt(8)
    ca.format.line.color.rgb = rgb(LINE)
    tot_py, tot_valve, tot_start = int(T["python_only"].sum()), int(T["valve"].sum()), int(T["start"].sum())
    tot_sal_only, tot_py_all = int(T["saleae_only"].sum()), int(T["python_licks"].sum())
    valid = int(T["py_only_valid_true"].sum())
    stat_card(s, 8.7, 1.4, 4.18, 1.05, f"{100 * tot_valve / max(1, tot_py):.0f} %", "while the valve was open",
              f"{tot_valve} of {tot_py} Python-only licks (incl. 2 within 1 ms of closing; "
              f"{100 * tot_py / max(1, tot_py_all):.2f} % of all)", "557A9E", "Card valve")
    stat_card(s, 8.7, 2.57, 4.18, 1.05, f"{tot_start}", "licks from before the start",
              "queued on the serial port, read at t = 0", PR_COLORS[0], "Card start")
    stat_card(s, 8.7, 3.74, 4.18, 1.05, f"{tot_sal_only}", "Saleae-only, none missed by the task",
              "an 8-ns spike; 2 contacts at the very end of Python", "6B7785", "Card saleae only")
    add_text(s, 8.7, 4.95, 4.18, 2.0,
             [[("Why: ", {"bold": True}), ("the Saleae line stays HIGH (one contact) while the valve is open, but the "
                                          "Arduino reads one LOW sample and sends a new LICK 10 ms later. ", {}),
               (f"{valid} of {tot_py} were valid and counted by the task", {"bold": True}),
               (" (2 at valve close with the 0.20-s lock-out in 9670f green1; 3 start licks). ", {}),
               ("Same mechanism in the previous cohort (day 1: 1.28 %; its PR days were lower: day 6 0.52 %, "
                "day 13 0 %); new day 3: 1.73 %.", {})]],
             size=11, color=TEXT, name="Lick explanation")
    add_footer(s, "Sources: analysis\\checks\\lick_mismatch_per_session.csv (24 sessions; recounted independently, "
                  "all numbers match); firmware pr_fr_task.ino; core_experiment.py:65-72 (serial reader starts when the "
                  f"port opens) · made {made}")


def slide_dayrules(prs, project: Path, made):
    p = project / "analysis" / "reference" / "prev_cohort_day_rules.csv"
    if not p.exists():
        return
    R = pd.read_csv(p).fillna("")
    s = prs.slides.add_slide(prs.slide_layouts[6])
    add_title(s, "Which days the previous cohort used, and our day sets",
              "Checked in the original MATLAB code and matched to the lab-meeting figures (pixel / MD5); "
              "no code says 'first day of phase'")
    cols = [("figure (deck slides)", 4.4), ("excluded days", 1.25), ("days used", 4.0)]
    t = s.shapes.add_table(len(R) + 1, len(cols), Inches(0.48), Inches(1.4), Inches(sum(w for _, w in cols)),
                           Inches(0.3 + 0.37 * len(R)))
    t.name = "Previous cohort day rules"
    tb = t.table
    for j, (c, w) in enumerate(cols):
        tb.columns[j].width = Inches(w)
        style_cell(tb.cell(0, j), [[(c, {"bold": True, "color": "FFFFFF"})]], NAVY, size=10)
    for i, r in enumerate(R.itertuples(), start=1):
        ex = str(r.excluded_days).replace(";", " ")
        fill = "FFF0D9" if ex.strip() not in ("", "none") else "FFFFFF"
        style_cell(tb.cell(i, 0), [[(f"{r.family}: ", {"bold": True, "size": 9}), (r._2, {"size": 9})]], "FFFFFF", size=9)
        style_cell(tb.cell(i, 1), [[(ex, {"bold": True, "size": 9})]], fill, size=9)
        style_cell(tb.cell(i, 2), [[(r.days_used, {"size": 9})]], "FFFFFF", size=9)
    add_text(s, 10.25, 1.4, 2.65, 5.6,
             [[("Our day sets", {"bold": True, "size": 14, "color": NAVY})],
              [("1  all days 3-18", {"size": 11})],
              [("2  first day of each phase out: 4 6 11 14 (option +17)", {"size": 11})],
              [("3  every list the previous cohort used: 4 6 11 14 17 · 3 6 11 14 · 4 6 11 14 · day ≥ 4 "
                "(+ the calendar first days 3 6 11 14 17)", {"size": 11})],
              [("tools\\phase_daysets.py, tools\\day_sets.csv -> analysis\\day_sets\\", {"size": 9, "color": MUTED})],
              [("Only written reason: day 3 is the first PR day and unstable (Active PR median 42.5 on day 3 vs "
                "28.5 on day 4).", {"size": 10})],
              [("Decided 11-16: day 4 is right for set 2 = [4 6 11 14]. Every set (8 in tools\\day_sets.csv) "
                "is computed for every analysis and every metric, not only for the figure it came from.",
                {"size": 10, "bold": True, "color": "1E6B34"})]],
             size=11, name="Day sets", spacing_after=6)
    add_footer(s, "Sources: previous-cohort code (step05/06/07/08/10), run_009 outputs, lab-meeting deck; "
                  "analysis\\reference\\prev_cohort_day_rules.csv and prev_cohort_day_rules_research_20261116.json "
                  f"(50 claims checked: 44 confirmed, 3 partly, 3 refuted) · made {made}")


def slide_besteffort(prs, project: Path, made):
    files = sorted((project / "pupil_tst_final").glob("day*/*/*/pupil_outputs_besteffort/besteffort_check.json"))
    if not files:
        return
    s = prs.slides.add_slide(prs.slide_layouts[6])
    add_title(s, "Best-effort pupil for hard eyes (narrowed / inflamed / cloudy)",
              "Every candidate of every frame, a prior from 2 independent labels per frame, and a whole-session "
              "path; frames without a visible pupil become gaps")
    head = ["session", "labelled frames", "package correct / wrong / false", "best effort correct / wrong / false",
            "not visible, said so", "pupil visible", "median diam / width"]
    widths = [1.9, 2.4, 1.75, 1.75, 1.3, 1.0, 1.4]
    tot = dict(pc=0, pw=0, pf=0, pm=0, bc=0, bw=0, bf=0, bm=0, bn=0, pn=0, n=0)
    be_list = read_csv(project / "pupil_tst_final" / "besteffort_sessions.csv")
    approved = set()
    if len(be_list):
        ok = be_list["approved"].astype(str).str.lower().str.startswith(("yes", "20"))
        approved = {f"day{r.day}_{r.cage}_{r.mouse}" for r in be_list[ok].itertuples()}
    rows, waiting = [], []
    for f in files:
        d = json.loads(f.read_text(encoding="utf-8"))
        b, p = d["besteffort_kfold"], d["package"]
        if d["session"] not in approved:
            waiting.append(d["session"].replace("_", " "))
        rows.append([d["session"].replace("_", " ") + ("" if d["session"] in approved else " (not approved)"),
                     f"{d['labelled_frames']}: {d['labels']['pupil']} outlined, "
                     f"{d['labels'].get('visible_no_contour', 0)} seen not outlined, "
                     f"{d['labels']['not_visible']} not visible",
                     f"{p['correct']} / {p['wrong']} / {p['false_pupil']}",
                     f"{b['correct']} / {b['wrong']} / {b['false_pupil']}",
                     f"{b['correct_none']} (package {p['correct_none']})",
                     f"{d['visible_pct']:.0f} %",
                     f"{d['median_diam_px']:.0f} / {d['median_width_px']:.0f} px"])
        for k, src, key in (("pc", p, "correct"), ("pw", p, "wrong"), ("pf", p, "false_pupil"), ("pm", p, "missed"),
                            ("bc", b, "correct"), ("bw", b, "wrong"), ("bf", b, "false_pupil"), ("bm", b, "missed"),
                            ("bn", b, "correct_none"), ("pn", p, "correct_none")):
            tot[k] += src[key]
    t = s.shapes.add_table(len(rows) + 1, len(head), Inches(0.48), Inches(1.45), Inches(sum(widths)),
                           Inches(0.45 + 0.42 * len(rows)))
    t.name = "Best effort validation"
    tb = t.table
    for j, (h_, w) in enumerate(zip(head, widths)):
        tb.columns[j].width = Inches(w)
        style_cell(tb.cell(0, j), [[(h_, {"bold": True, "color": "FFFFFF"})]], NAVY, size=9)
    for i, r in enumerate(rows, start=1):
        for j, v in enumerate(r):
            style_cell(tb.cell(i, j), [[(v, {"bold": j in (0, 3)})]], "E6F4EA" if j == 3 else "FFFFFF", size=10)
    pk_acc = (tot["pc"] + tot["pn"]) / max(1, tot["pc"] + tot["pw"] + tot["pf"] + tot["pm"] + tot["pn"])
    be_acc = (tot["bc"] + tot["bn"]) / max(1, tot["bc"] + tot["bw"] + tot["bf"] + tot["bm"] + tot["bn"])
    yc = max(3.95, 1.45 + 0.45 + 0.43 * len(rows) + 0.2)
    stat_card(s, 0.5, yc, 4.0, 1.05, f"{100 * be_acc:.0f} %", "best effort, right on labelled frames",
              f"package {100 * pk_acc:.0f} % (4-fold over each session's labels)", "1E6B34", "Card BE accuracy")
    stat_card(s, 4.7, yc, 4.0, 1.05, f"{tot['bf'] + tot['bw']}", "wrong picks, best effort",
              f"package {tot['pf'] + tot['pw']} (of the same frames)", PR_COLORS[0], "Card wrong")
    status = (f"In the tables and Step 02 since 11-16 (you approved): {len(rows) - len(waiting)} sessions in "
              "besteffort_sessions.csv. Day 2 9669m black: 2 fewer correct frames than the package (14 vs 16 of 24).")
    if waiting:
        status += f" Waiting for your OK (not in any table): {', '.join(waiting)}."
    add_text(s, 9.0, yc, 3.85, 7.0 - yc,
             [[("How", {"bold": True, "size": 13, "color": NAVY})],
              [("1  pupil_all_candidates.m: all U-Net blobs + dark regions of every frame (package code)", {})],
              [("2  prior per session from its labelled frames: position, size, darkness", {})],
              [("3  Viterbi path: cost of each candidate + jumps; 'not visible' is a state", {})],
              [("4  checked: 4-fold over the labels; settings chosen on all 4 sessions together", {})],
              [(status, {"bold": True, "color": "4B3A8C"})]],
             size=10, name="BE method", spacing_after=4)
    add_text(s, 0.5, yc + 1.25, 8.3, 7.0 - (yc + 1.25),
             [[("Reading the table: ", {"bold": True}),
               ("correct = the picked blob is the labelled pupil; wrong = another blob; false = a pupil picked where the "
                "labellers saw none. Width = horizontal size of the blob; in a narrowed eye the lids cut the top and "
                "bottom of the pupil, so the width is the better size there (9669m black).", {})]],
             size=10, color=TEXT, name="BE note")
    add_footer(s, "Sources: pupil_tst_final\\dayN\\cage\\mouse\\pupil_outputs_besteffort\\besteffort_check.json; "
                  f"tools\\pupil_besteffort.py; labels: 2 labellers + adjudication per frame · made {made}")


def appendix(prs, project: Path, days):
    divider_slide(prs, "Appendix · all figures",
                  "Every figure of the daily routine (days " + ", ".join(str(d) for d in days) +
                  "), the across-day figure, today's eye-box and pupil checks")
    n = 0
    for d in days:
        figs = sorted((project / "analysis" / "days" / f"day{d}" / "figs").glob("*.png"))
        extra = sorted((project / "analysis" / "days" / f"day{d}" / "figs_extra").glob("*.png"))
        for p in figs + extra:
            n += image_slide(prs, p, f"Day {d} · {p.stem.split('_', 1)[-1].replace('_', ' ')}",
                             CAPTIONS.get(p.stem) or CAPTIONS.get(p.stem.split("_", 1)[-1], p.stem.replace("_", " ")))
        for kind, sub in (("qc_sheet", "pupil_outputs_general"), ("besteffort_trace", "pupil_outputs_besteffort"),
                          ("fullpupil_sheet", "pupil_outputs_fullpupil")):
            for p in sorted((project / "pupil_tst_final" / f"day{d}").glob(f"*/*/{sub}/{kind}.png")):
                cage, mouse = p.parts[-4], p.parts[-3]
                n += image_slide(prs, p, f"Day {d} · {cage} {mouse} · {kind.replace('_', ' ')}", CAPTIONS[kind])
    for p in sorted((project / "analysis" / "figs_across_days").glob("*.png")):
        n += image_slide(prs, p, "Across days", CAPTIONS.get(p.stem, p.stem))
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", default=str(PROJECT_DEFAULT))
    ap.add_argument("--days", type=int, nargs="+", required=True)
    ap.add_argument("--out", help="default: <project>\\analysis\\ppt\\days<first>-<last>_summary_<date>.pptx")
    ap.add_argument("--no-appendix", action="store_true", help="leave out the figure slides")
    a = ap.parse_args(argv)
    project = Path(a.project)
    days = sorted(set(a.days))
    out = Path(a.out) if a.out else project / "analysis" / "ppt" / \
        f"days{days[0]}-{days[-1]}_summary_{datetime.now():%Y%m%d}.pptx"
    p, table_bottom = build(project, days, out, with_appendix=not a.no_appendix)
    print(f"written {p}  (slide-3 table ends at about {table_bottom:.2f} in of 7.5)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
