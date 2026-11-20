"""qc_summary.py - ONE QC picture for all days so far = the daily report (user 2026-11-17: "daily routine need to focus
on QC ... later I will do analysis after 18 days"). Reads only what the routine already made; changes nothing.
Per mouse x day:
  rewards (Python = Saleae?) and PR (last ratio requirement = pr_RL_full of pr_checks, as the previous cohort's PR)
  valve time and session length
  licks Saleae / Python, and how many Python-only licks the task COUNTED (with the reason of each)
  video frames vs camera pulses (and whether the approved fix is in the table)
  pupil grade (active pupil tree; whole eye for the mice in wholeeye_mice.csv)
Every CHECK item of the session check is either explained (yellow, with a note) or shown red with its text.
A red flag the user has decided on (analysis/days/qc_decisions.csv: day, cage, mouse, cat, decision, decided; empty
day = every day of that mouse) is shown yellow with the decision (added 2026-11-18).
Colour: green = all fine; yellow = deviation that is known, fixed or approved (note below); red = needs a look.
More than 6 days: short cells (colour + rewards/PR + flags); the notes list stays complete.
Outputs: analysis/days/qc_summary_days<first>-<last>.png and .csv
usage: python qc_summary.py [--days 1 2 3 4]"""
import argparse
import re
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

T = Path(__file__).resolve().parent
sys.path.insert(0, str(T))
import concat_choice as CC  # noqa: E402

P = T.parent
CAGES = ("9669m", "9670f")
MICE = [(c, m) for c in CAGES for m in ("black", "green1", "green2", "white")]
COL = {"green": "#d9f2d9", "yellow": "#fff1c2", "red": "#f8c9c9"}
STAT_RE = re.compile(r"^\s{2}(PASS|INFO|CHECK)\s+(.{16}) (.*)$")


def check_texts(report):
    """{item: text} of the CHECK lines (with their indented continuation lines) of a session check report"""
    out = {}
    try:
        L = Path(str(report)).read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return out
    for i, s in enumerate(L):
        m = STAT_RE.match(s)
        if not m or m.group(1) != "CHECK":
            continue
        txt, j = [m.group(3).strip()], i + 1
        while j < len(L) and L[j].startswith(" " * 10) and not STAT_RE.match(L[j]):
            txt.append(L[j].strip())
            j += 1
        out[m.group(2).strip()] = " ".join(txt)
    return out


def py_only_counted(day, cage, mouse, session):
    """[(kind, text)] for every Python lick the task counted (valid) that has no Saleae lick edge"""
    f = P / "sync" / f"day{day}" / cage / mouse / f"{session}_events.csv"
    if not f.exists():
        return [("unknown", "events file missing")]
    E = pd.read_csv(f)
    lick = E[E.event == "lick"]
    py = lick[(lick.source == "python") & lick.python_delay_ms.isna()
              & (pd.to_numeric(lick.valid_or_hit, errors="coerce") == 1)]
    sal_un = lick[(lick.source == "saleae")
                  & lick.note.astype(str).str.contains("no matching python lick")].t_saleae_s.to_numpy(float)
    valve = E[E.event == "reward_valve"].t_saleae_s.to_numpy(float)
    out = []
    for _, r in py.iterrows():
        el, t = float(r.python_elapsed_s), float(r.t_saleae_s)
        prev = valve[valve <= t]
        dv = t - prev.max() if len(prev) else np.inf
        if el < 0.1:
            out.append(("start", f"{el:.3f} s = a lick from before the start (serial buffer)"))
        elif 0.15 <= dv <= 0.26:
            out.append(("valve_close", f"{el:.0f} s = at valve close ({1000 * dv:.0f} ms after the reward)"))
        elif len(sal_un) and np.min(np.abs(sal_un - t)) < 0.15:
            out.append(("late", f"{el:.0f} s = a real lick (on the Saleae) that Python read late"))
        else:
            out.append(("unknown", f"{el:.1f} s = no Saleae lick near it"))
    return out


def wholeeye_offeye(d, c, m):
    """(% frames whose eye mask is > 20 px from the eye centre of the surrounding minute, the same % in the worst tenth
    of the session) from the eye model's per-frame output; None without that file. Off-eye masks are the spout,
    whiskers or fur picked when the damaged eye is closed (review 2026-11-17)."""
    f = P / "analysis" / "wholeeye" / f"day{d}" / c / m / "wholeeye_frames.csv"
    if not f.exists():
        return None
    W = pd.read_csv(f, usecols=["found", "cx", "cy"])
    fnd = W.found.to_numpy() == 1
    cx, cy = np.where(fnd, W.cx, np.nan), np.where(fnd, W.cy, np.nan)
    med = lambda x: pd.Series(x).rolling(1800, center=True, min_periods=100).median().to_numpy()
    with np.errstate(invalid="ignore"):
        off = fnd & (np.hypot(cx - med(cx), cy - med(cy)) > 20)
    return 100 * off.mean(), 100 * max(x.mean() for x in np.array_split(off, 10))


def open_issues():
    """{(day, cage, mouse): [(title, issue)]} from analysis/days/qc_open_issues.csv (status not 'resolved')"""
    f = P / "analysis" / "days" / "qc_open_issues.csv"
    out = {}
    if f.exists():
        for _, r in pd.read_csv(f, dtype=str, encoding="utf-8-sig").fillna("").iterrows():
            if r.status.strip().lower() != "resolved":
                out.setdefault(CC.norm_key(r.day, r.cage, r.mouse), []).append((r.get("title", "") or r.issue[:60], r.issue))
    return out


def decisions():
    """{(cage, mouse, cat): [(day or None, decision, date)]} from analysis/days/qc_decisions.csv: what the user chose
    for a red QC flag (empty day = every day of that mouse). Missing file = no decisions."""
    out = {}
    for r in CC.read_rows(P / "analysis" / "days" / "qc_decisions.csv", need=("cage", "mouse", "cat", "decision")):
        day = int(r["day"]) if str(r.get("day", "")).strip() else None
        out.setdefault((r["cage"].strip().lower(), r["mouse"].strip().lower(), r["cat"].strip().lower()), []).append(
            (day, r["decision"].strip(), str(r.get("decided", "")).strip()))
    return out


def pupil_check_out():
    """{(day, cage, mouse): % of frames out} from analysis/pupil_check/pupil_check_day*.csv (pupil_check_examples.py,
    added 2026-11-18: held + jumps + odd size of the pupil the tables use; the newest file wins). Empty if none."""
    out = {}
    files = sorted((P / "analysis" / "pupil_check").glob("pupil_check_day*.csv"), key=lambda p: p.stat().st_mtime)
    for f in files:
        try:
            X = pd.read_csv(f, dtype={"cage": str}, encoding="utf-8-sig")
        except (OSError, ValueError):
            continue
        for r in X.itertuples():
            if pd.notna(getattr(r, "out_pct", np.nan)):
                out[(int(r.day), str(r.cage), str(r.mouse))] = float(r.out_pct)
    return out


def decided(keys, d, c, m, cat):
    """'decision (date)' if the user decided this flag for this session, else None"""
    for day, dec, when in keys.get("decisions", {}).get((str(c).lower(), str(m).lower(), cat), []):
        if day is None or day == int(d):
            return f"{dec} ({when})" if when else dec
    return None


def one_session(s, d, c, m, q, pr, keys):
    """(level 0/1/2, cell lines, short flags, notes) for one session"""
    F, flags = [], []                         # F: (level, category, short, note)
    items = [x.strip() for x in str(s.check_items).split(";") if x.strip()] if pd.notna(s.check_items) else []
    # the stored path can be from before the project moved into OneDrive: look in this project's folder first
    rep = P / "analysis" / "session_checks" / f"day{d}" / Path(str(s.report).replace("\\", "/")).name
    texts = check_texts(rep if rep.exists() else s.report) if items else {}
    handled = set()
    key = CC.norm_key(d, c, m)
    num = lambda v: int(v) if pd.notna(v) else None
    # rewards
    rc, vo, rm = num(s.reward_cmd), num(s.valve_openings), num(s.rewards_matched)
    rew_ok = rc is not None and rc == vo and (rm is None or rm == rc)
    if not rew_ok:
        F.append((2, "rewards", f"Py {rc} / Sal {vo}", f"rewards differ: Python {rc}, Saleae valve openings {vo}, matched {rm}"))
        flags.append("REWARDS")
    handled |= {"D reward count", "D reward match"} if not rew_ok else set()
    # session length
    if "C session length" in items:
        longer = float(s.session_s) > 910
        F.append((1 if longer else 2, "length", f"{float(s.session_s) / 60:.1f} min",
                  f"session ran {float(s.session_s) / 60:.1f} min (protocol 15.2 min)"))
        flags.append("length")
        handled.add("C session length")
    # video frames vs camera pulses
    has_video = isinstance(s.video, str) and bool(s.video.strip())
    fmp = num(s.frames_minus_pulses)
    if not has_video:
        ok = key in keys["behavior_only"]
        frames_txt = "no video" + (" (approved)" if ok else "")
        F.append((1, "no_video", "no video", "no video; the approved behaviour-only table is used (Saleae + task log "
                  "complete)") if ok else (2, "no_video_open", "no video", "no video found"))
        flags.append("no video")
        handled.add("A video")
    elif fmp:
        ok = key in keys["sync_fix"]
        frames_txt = f"frames {fmp:+d} vs pulses ({'fixed' if ok else 'NOT fixed'})"
        F.append((1 if ok else 2, "frames" if ok else "frames_open", f"{fmp:+d} frames",
                  f"video has {abs(fmp)} {'more' if fmp > 0 else 'fewer'} frames than camera pulses - "
                  + ("the approved alignment fix is in the table" if ok else "NOT fixed, needs a look")))
        flags.append("frames fixed" if ok else "FRAMES")
        handled.add("F frames")
        for it in ("F other runs", "F video closed"):
            if it in items and ok:
                F.append((1, "frames_cause", it, f"  (cause) {texts.get(it, it)[:130]}"))
                handled.add(it)
    else:
        frames_txt = "frames = pulses"
    if "F clock" in items and "PC clock changed" in texts.get("F clock", ""):
        F.append((1, "clock", "PC clock jump", "PC clock jumped during the session (Windows time sync); the sync uses "
                  "the camera train"))
        handled.add("F clock")
    # licks
    n_cnt = num(s.python_only_valid) or 0
    lick_txt = f"licks Sal {num(s.licks_saleae)}, Py {num(s.licks_python)}" + (f" ({n_cnt} Py-only)" if n_cnt else "")
    if "replay" in str(getattr(s, "paradigm", "")).lower():
        # 2026-11-19 (user): the Passive spout sits inside the mouth - no licking is possible, so the lick line of a
        # passive replay session shows contacts, not licks
        lick_txt = f"contacts Sal {num(s.licks_saleae)}, Py {num(s.licks_python)} (not licks)"
    if n_cnt:
        po = py_only_counted(d, c, m, s.session)
        lock = float(s.lockout_median_s) if pd.notna(s.lockout_median_s) else np.nan
        lv = 1
        for kind, _ in po:
            if kind == "unknown" or (kind == "valve_close" and not lock < 0.21):
                lv = 2
        why = "; ".join(t for _, t in po)
        extra = (" - known at lock-out 0.20 s (the 0.22-s lock-out stops these)"
                 if any(k == "valve_close" for k, _ in po) and lock < 0.21 else "")
        for kind in ("start", "valve_close", "late", "unknown"):
            n_k = sum(1 for k, _ in po if k == kind)
            if n_k:
                bad = kind == "unknown" or (kind == "valve_close" and not lock < 0.21)
                F.append((2 if bad else 1, "lick_" + kind, f"{n_k} lick{'s' if n_k > 1 else ''}",
                          f"the task counted {n_k} lick(s) the Saleae did not see: "
                          + "; ".join(tx for k, tx in po if k == kind) + (extra if kind == "valve_close" else "")))
        flags.append("licks" if lv == 1 else "LICKS")
        handled.add("E Python-only")
    if "H sync" in items:
        t = texts.get("H sync", "")
        if "[WARN] licks:" in t:
            extra_py = (num(s.licks_python) or 0) - (num(s.licks_saleae) or 0)
            F.append((1 if not n_cnt or "E Python-only" in handled else 2, "licks_extra", f"+{extra_py} licks",
                      f"Python logged {extra_py} more licks than the Saleae; {n_cnt} of them counted by the task"))
            handled.add("H sync")
        elif ("frames than camera pulses" in t and "F frames" in handled) or ("no video file" in t and "A video" in handled):
            handled.add("H sync")
    if "B cage/mouse" in items and "filed as" in texts.get("B cage/mouse", ""):
        F.append((1, "filename", f"'{s.session}'",
                  f"file name '{s.session}' does not name the mouse; filed as {c} {m} by start time"))
        handled.add("B cage/mouse")
    for it in items:
        if it not in handled:
            F.append((2, "check_open", it, f"CHECK {it}: {texts.get(it, '(see the session check report)')[:150]}"))
            flags.append(it.upper())
    # pupil
    g = str(q.iloc[0].grade) if len(q) else ("-" if not has_video else "missing")
    # 2026-11-18: right after an experiment the light run makes the QC before the pupil step; a day with no pupil
    # grade at all is "not run yet" (yellow), not 8 red sessions (day 5 showed 8 reds). Missing in a day that has
    # other grades stays red.
    if g == "missing" and int(d) not in keys.get("pupil_days", set()):
        g = "pending"
    nb = float(q.iloc[0].no_blob) if len(q) and pd.notna(q.iloc[0].no_blob) else np.nan
    if g == "pending":
        F.append((1, "pupil_pending", "pupil not run yet",
                  "pupil: the package / whole-eye step has not run for this day yet (it runs after the PR + QC report)"))
        flags.append("pupil not run yet")
        pupil_txt = "pupil not run yet"
    elif g == "whole_eye":
        F.append((1, "whole_eye", "whole eye", "whole eye instead of the pupil (eye damage; approved in wholeeye_mice.csv)"))
        flags.append("whole eye")
        pupil_txt = "whole eye (approved)"
        oe = wholeeye_offeye(d, c, m)
        if oe is not None:
            pupil_txt = f"whole eye, {oe[0]:.0f}% off eye"
            dec = decided(keys, d, c, m, "off_eye") if oe[1] > 25 else None
            if dec:
                F.append((1, "off_eye_kept", f"{oe[0]:.0f}% off, worst tenth {oe[1]:.0f}%",
                          f"whole eye: {oe[0]:.0f}% of frames have the mask off the eye (spout/whiskers/fur), "
                          f"{oe[1]:.0f}% in the worst tenth of the session - your decision: {dec}"))
                flags.append("off eye kept")
            elif oe[1] > 25:
                F.append((2, "off_eye", f"{oe[0]:.0f}% off, worst tenth {oe[1]:.0f}%",
                          f"whole eye: {oe[0]:.0f}% of frames have the mask off the eye (spout/whiskers/fur), "
                             f"{oe[1]:.0f}% in the worst tenth of the session - needs a look"))
                flags.append("OFF EYE")
    elif g == "noisy":
        why = str(q.iloc[0].reason).replace("new package: ", "") if pd.notna(q.iloc[0].reason) else ""
        po = keys.get("pupil_out", {}).get((int(d), str(c), str(m)))   # 2026-11-18: frames out of the table's pupil
        short = f"{po:.0f}% out" if po is not None else (f"{100 * nb:.0f}% no pupil" if np.isfinite(nb) else "noisy")
        F.append((1, "pupil_noisy", short,
                  f"pupil noisy (graded by eye): {textwrap.shorten(why, 120, placeholder=' ...')}"))
        flags.append("pupil noisy")
        pupil_txt = (f"pupil noisy, {po:.0f}% out" if po is not None else
                     f"pupil noisy, {100 * nb:.0f}% none" if np.isfinite(nb) else "pupil noisy")
    elif g == "ok":
        po = keys.get("pupil_out", {}).get((int(d), str(c), str(m)))
        pupil_txt = (f"pupil ok, {po:.0f}% out" if po is not None else
                     f"pupil ok, {100 * nb:.0f}% none" if np.isfinite(nb) else "pupil ok")
    elif g == "-":
        pupil_txt = "pupil -"
    else:
        F.append((2, "pupil_open", g, f"pupil grade '{g}'"))
        flags.append("PUPIL")
        pupil_txt = f"pupil {g}"
    for title, issue in keys["open"].get(key, []):
        F.append((2, "open_issue", title, f"open QC issue: {issue}"))
        flags.append("OPEN")
    level = max([f[0] for f in F], default=0)
    vm = float(s.valve_median_ms) if pd.notna(s.valve_median_ms) else np.nan
    pr_txt = f"  PR {pr}" if pr is not None else ""
    lines = [f"rewards {rc} {'= Sal' if rew_ok else '!= Sal'}{pr_txt}",
             f"valve {vm:.0f} ms, {float(s.session_s) / 60:.1f} min", lick_txt, frames_txt, pupil_txt]
    return level, lines, flags, [f[3] for f in F], F


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--days", type=int, nargs="+")
    a = ap.parse_args(argv)
    S = pd.read_csv(P / "analysis" / "session_checks" / "session_checks.csv", dtype={"cage": str}, encoding="utf-8-sig")
    S = S[S.cage.isin(CAGES)].sort_values("checked_at").drop_duplicates(["day", "cage", "mouse"], keep="last")
    days = a.days or sorted(S.day.astype(int).unique())
    root = CC.active_pupil_root(P)
    qf = root / "pupil_session_quality.csv"
    Q = pd.read_csv(qf, dtype={"cage": str}) if qf.exists() else pd.DataFrame(columns=["day", "cage", "mouse"])
    keys = {"sync_fix": CC.approved_sessions(root, CC.LISTS[CC.SYNC]),
            "behavior_only": CC.approved_sessions(root, "behavioronly_sessions.csv"), "open": open_issues(),
            "decisions": decisions(), "pupil_out": pupil_check_out(),
            "pupil_days": set(pd.to_numeric(Q["day"], errors="coerce").dropna().astype(int)) if len(Q) else set()}
    PRS = {}
    for d in days:
        f = P / "analysis" / "days" / f"day{d}" / f"pr_checks_day{d}.csv"
        if f.exists():
            X = pd.read_csv(f, encoding="utf-8-sig", dtype={"cage": str})
            for _, r in X.iterrows():
                if pd.notna(r.get("pr_RL_full")):
                    PRS[(d, str(r.cage), str(r.mouse).split("_")[-1].lower())] = int(r.pr_RL_full)
    rows, issues = [], []
    for d in days:
        for c, m in MICE:
            s = S[(S.day == d) & (S.cage == c) & (S.mouse == m)]
            q = Q[(Q.day == d) & (Q.cage == c) & (Q.mouse.astype(str).str.lower() == m)]
            if s.empty:
                rows.append(dict(day=d, cage=c, mouse=m, status="red", cell="no session check\n(session missing?)",
                                 flags="MISSING", notes="no session check row (session missing?)"))
                continue
            lv, lines, flags, notes, FF = one_session(s.iloc[0], d, c, m, q, PRS.get((d, c, m)), keys)
            issues += [dict(day=d, cage=c, mouse=m, level=["green", "yellow", "red"][f[0]], cat=f[1], short=f[2],
                            note=f[3]) for f in FF]
            rows.append(dict(day=d, cage=c, mouse=m, status=["green", "yellow", "red"][lv], cell="\n".join(lines),
                             flags=", ".join(flags), notes=" | ".join(notes),
                             rewards=s.iloc[0].reward_cmd, pr=PRS.get((d, c, m)), session=s.iloc[0].session))
    R = pd.DataFrame(rows)
    tag = f"days{days[0]}-{days[-1]}" if len(days) > 1 else f"day{days[0]}"
    out = P / "analysis" / "days" / f"qc_summary_{tag}.png"
    R.to_csv(out.with_suffix(".csv"), index=False)
    pd.DataFrame(issues, columns=["day", "cage", "mouse", "level", "cat", "short", "note"]).to_csv(
        P / "analysis" / "days" / f"qc_issues_{tag}.csv", index=False)
    # general notes: what PR means, valve per day
    V = S[S.day.isin(days)].copy()
    V["v"] = (V.valve_median_ms / 100).round() * 100
    vtxt = "; ".join(f"day {d}: " + ", ".join(f"{int(v)} ms x{n}" for v, n in g.v.value_counts().sort_index().items())
                     for d, g in V.groupby("day"))
    general = ["PR = last ratio requirement (pr_RL_full = rewards + 1 in a full 910-s session), as the previous cohort's PR. "
               "'(n Py-only)' = licks the task counted that the Saleae lick line did not see.",
               f"valve time per day: {vtxt}"]
    notes = [f"day{r.day} {r.cage} {r.mouse}: {r.notes}" for r in R.itertuples() if str(r.notes) not in ("", "nan")]
    wrapped = []
    for n in general + [""] + notes:
        wrapped += textwrap.wrap(n, 190, subsequent_indent="      ") or [""]
    compact = len(days) > 6
    cw = 1.25 if compact else 3.2
    fig_w = 2.0 + cw * len(days)
    grid_h, notes_h = (0.55 if compact else 1.25) * len(MICE), 0.16 * len(wrapped) + 0.3
    fig = plt.figure(figsize=(max(fig_w, 14), grid_h + notes_h + 0.9))
    gs = fig.add_gridspec(2, 1, height_ratios=[grid_h, notes_h], hspace=0.04)
    ax = fig.add_subplot(gs[0])
    ax.set_xlim(0, len(days)); ax.set_ylim(len(MICE), -0.35); ax.axis("off")
    for i, (c, m) in enumerate(MICE):
        ax.text(-0.03, i + 0.5, f"{c} {m}", ha="right", va="center", fontsize=11, fontweight="bold")
        for j, d in enumerate(days):
            r = R[(R.day == d) & (R.cage == c) & (R.mouse == m)].iloc[0]
            ax.add_patch(plt.Rectangle((j + 0.02, i + 0.04), 0.96, 0.92, color=COL[r.status]))
            if compact:
                # r["flags"], not r.flags: on a pandas row, .flags is pandas' own Flags object (fixed 2026-11-20,
                # the first day with more than 6 days showed "<Flags(allows_duplicate_labels=True)>" in every cell)
                fl = r.get("flags")
                fl = "" if pd.isna(fl) else str(fl)
                top = (f"R{int(r.rewards)}" if pd.notna(r.get("rewards")) else "") + \
                      (f" PR{int(r.pr)}" if pd.notna(r.get("pr")) else "")
                txt = "\n".join([top] + (textwrap.wrap(fl, 21)[:2] or [""]))
                ax.text(j + 0.06, i + 0.5, txt, va="center", fontsize=7, family="monospace")
            else:
                ax.text(j + 0.05, i + 0.5, r.cell, va="center", fontsize=8.3, family="monospace")
    for j, d in enumerate(days):
        ax.text(j + 0.5, -0.12, f"day {d}", ha="center", fontsize=12, fontweight="bold")
    n_red, n_y = int((R.status == "red").sum()), int((R.status == "yellow").sum())
    ax.set_title(f"QC summary, {tag.replace('days', 'days ')}: {len(R)} sessions - green {len(R) - n_red - n_y} (all fine), "
                 f"yellow {n_y} (known / fixed / approved, see notes), red {n_red} (needs a look)",
                 fontsize=12, loc="left")
    axn = fig.add_subplot(gs[1]); axn.axis("off")
    axn.text(0, 1, "\n".join(wrapped), va="top", ha="left", fontsize=7.6, family="monospace", transform=axn.transAxes)
    fig.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(out)
    print(R.groupby(["day", "status"]).size().unstack(fill_value=0).to_string())
    for r in R[R.status == "red"].itertuples():
        print(f"  RED day{r.day} {r.cage} {r.mouse}: {r.notes}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
