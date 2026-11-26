#!/usr/bin/env python3
"""
pr_table_days.py - one table of every mouse's PR score on every day so far (added 2026-11-18; user: "all days' PR
score table show me them"). Report only: reads the routine's tables, changes nothing.

  PR (PR days)  pr_RL_full of analysis/days/dayN/pr_checks_dayN.csv = the last ratio requirement = rewards + 1 over
                the WHOLE session (910 s; day 13 = 1210 s, user 2026-11-26: "PR you can use whole session")
                (the previous cohort's Step 02 RequirementLast; Claude chose it on 2026-11-14 to
                match the old pipeline, and the user kept it on 2026-11-19: "네, 그대로". The old Step 10 used max
                Trial = rewards instead)
  FR days       no PR: the cell shows the ratio (FR3 day 1, FR5 day 2) and the rewards in the first 15 min 10 s
                (analysis/session_metrics.csv, window first_15m10s; two day-1 sessions ran untimed and longer)
  rewards check Python reward commands = Saleae valve openings (pr_checks_dayN.csv), per session
  reference     the previous cohort's PR of that day in its first 910 s (median and range; pr_checks_dayN.csv
                prev_RL910_*; reference only - training and valve differ)
  passive       (days 6-10) the drops delivered (R) and, since 2026-11-19 (user: "네, 같이 적기"), the replayed
                active partner's R and PR in the same cell (columns partner, partner_pr). A passive mouse has no PR of
                its own; its R equals the partner's R, which is the partner's PR - 1, so R alone looked like "PR - 1".
Outputs: analysis/days/pr_table_days<a>-<b>.csv (one row per mouse and day) and pr_table_days<a>-<b>.png

usage: python pr_table_days.py [--days 1 2 3 4 5]      (default: every day with a pr_checks file)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

T = Path(__file__).resolve().parent
P = T.parent
DAYS_DIR = P / "analysis" / "days"
MICE = [(c, m) for c in ("9669m", "9670f") for m in ("black", "green1", "green2", "white")]


def groups():
    f = P / "cohort_manifest.csv"
    if not f.exists():
        return {}
    G = pd.read_csv(f, dtype=str).fillna("")
    return {(r.cage_folder, r.mouse): (r.group, r.pair_id) for r in G.itertuples()}


def all_days():
    return sorted(int(p.name[3:]) for p in DAYS_DIR.glob("day*") if p.name[3:].isdigit()
                  and (p / f"pr_checks_{p.name}.csv").exists())


def table(days):
    """one row per mouse and day"""
    SM = pd.read_csv(P / "analysis" / "session_metrics.csv", dtype={"cage": str})
    SM = SM[SM.window == "first_15m10s"]
    G = groups()
    rows = []
    for d in days:
        f = DAYS_DIR / f"day{d}" / f"pr_checks_day{d}.csv"
        X = pd.read_csv(f, encoding="utf-8-sig", dtype={"cage": str}) if f.exists() else pd.DataFrame()
        for c, m in MICE:
            r = X[(X.cage == c) & (X.mouse.astype(str).str.endswith("_" + m))] if len(X) else X
            grp, pair = G.get((c, m), ("", ""))
            base = dict(day=d, cage=c, mouse=m, group=grp, pair_id=pair)
            if r.empty:
                rows.append(dict(base, task="missing"))
                continue
            r = r.iloc[-1]
            s = SM[(SM.day == d) & (SM.mouse == f"{c}_{m}")]
            task = str(s.task.iloc[-1]) if len(s) else ""
            is_pr = pd.notna(r.get("pr_RL_full"))
            # days 6-10 (2026-11-18): a passive mouse replays its active partner's log of the day (no PR of its own)
            is_passive = "replay" in str(r.get("paradigm", "")).lower()
            src = str(r.get("replay_source", "") or "") if is_passive else ""
            n_src = r.get("replay_source_rewards", np.nan) if is_passive else np.nan
            rows.append(dict(
                base, task="PR" if is_pr else ("PASSIVE" if is_passive else (task or str(r.get("paradigm", "")))),
                session=r.get("session", ""),
                replay_source="" if src in ("", "nan") else src,
                replay_source_rewards=float(n_src) if pd.notna(n_src) else np.nan,
                replay_match=bool(pd.notna(n_src) and int(r.rewards_python) == int(float(n_src))) if is_passive
                else np.nan,
                pr=float(r.pr_RL_full) if is_pr else np.nan,
                pr_910s=float(r.pr_RL_910s) if pd.notna(r.get("pr_RL_910s")) else np.nan,
                rewards_910s=float(s.rewards.iloc[-1]) if len(s) else np.nan,
                rewards_python=int(r.rewards_python), rewards_saleae=int(r.rewards_saleae) if pd.notna(
                    r.rewards_saleae) else np.nan,
                rewards_match=bool(pd.notna(r.rewards_saleae) and int(r.rewards_python) == int(r.rewards_saleae)),
                session_s=float(r.session_s) if pd.notna(r.get("session_s")) else np.nan,
                # the previous cohort's PR in ITS first 910 s (it was hand-stopped at 13-17 min), as this cohort's
                # 910-s sessions; the whole-session value differs only on day 3 (39 vs 38.5)
                prev_pr_median=r.get("prev_RL910_median", r.get("prev_RL_median", np.nan)),
                prev_pr_min=r.get("prev_RL910_min", r.get("prev_RL_min", np.nan)),
                prev_pr_max=r.get("prev_RL910_max", r.get("prev_RL_max", np.nan)), prev_pr_n=r.get("prev_RL_n", np.nan)))
    Tb = pd.DataFrame(rows)
    if Tb.empty:
        return Tb
    # 2026-11-19 (user: "네, 같이 적기"): the replayed active partner and its PR next to a passive mouse's R
    Tb["partner"], Tb["partner_pr"] = "", np.nan
    for k, r in Tb[Tb.task == "PASSIVE"].iterrows():
        src = str(r.get("replay_source", "") or "")
        if src in ("", "nan"):
            continue
        name = partner_of(src)
        Tb.at[k, "partner"] = name
        pc, _, pm = name.partition(" ")
        q = Tb[(Tb.day == r.day) & (Tb.cage == pc) & (Tb.mouse == pm) & (Tb.task == "PR")]
        if len(q) and pd.notna(q.pr.iloc[0]):
            Tb.at[k, "partner_pr"] = float(q.pr.iloc[0])
    return Tb


def day_type(Tb: pd.DataFrame, d, short=False):
    """'PR', 'PR + passive' (days 6-10), 'passive', or the FR task of the day ('FR3', 'FR5')"""
    t = set(Tb[(Tb.day == d) & (Tb.task != "missing")].task)
    if "PR" in t and "PASSIVE" in t:
        return "PR+P" if short else "PR + passive"
    if "PR" in t:
        return "PR"
    if "PASSIVE" in t:
        return "P" if short else "passive"
    fr = sorted(x for x in t if x)
    return fr[0] if fr else ""


def partner_of(src: str):
    """'9670fgreen2_prday6.jsonl' -> '9670f green2' (the active partner a passive log replayed)"""
    import re
    m = re.match(r"(\d{4}[mf])_?(green1|green2|green|white|black)", Path(str(src)).name.lower())
    if not m:   # 2026-11-19 (user): cage and colour anywhere in the name
        m = re.search(r"(?<!\d)(\d{4}[mf]).*?(green1|green2|green|white|black)", Path(str(src)).name.lower())
    if not m:
        return Path(str(src)).stem
    mouse = "green1" if m.group(2) == "green" else m.group(2)
    return f"{m.group(1)} {mouse}"


def figure(Tb: pd.DataFrame, days, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nr, nc = len(MICE) + 3, len(days) + 1                 # header + mice + median + previous cohort
    fig = plt.figure(figsize=(1.95 * nc + 0.3, 0.72 * (nr + 1.6)), dpi=130)
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
    ax.set_xlim(-0.05, nc + 0.05)
    ax.set_ylim(nr + 0.75, -0.85)
    ax.axis("off")
    pr_vals = Tb.pr.dropna()
    lo, hi = (pr_vals.min(), pr_vals.max()) if len(pr_vals) else (0, 1)
    cmap = plt.get_cmap("Blues")

    def cell(x, y, txt, sub="", face="white", bold=False, col="#1f2933", sub2=""):
        ax.add_patch(plt.Rectangle((x, y), 1, 1, facecolor=face, edgecolor="#c9ced6", lw=0.8))
        ty = 0.52 if not sub else (0.36 if sub2 else 0.42)
        ax.text(x + 0.5, y + ty, txt, ha="center", va="center", fontsize=14 if sub else 12,
                fontweight="bold" if bold else "normal", color=col)
        if sub:
            ax.text(x + 0.5, y + (0.66 if sub2 else 0.78), sub, ha="center", va="center", fontsize=8, color="#52606d")
        if sub2:
            ax.text(x + 0.5, y + 0.86, sub2, ha="center", va="center", fontsize=7.5, color="#9a5b00")

    ax.text(0.0, -0.45, "PR score per mouse and day  ·  PR = last ratio requirement = rewards + 1  ·  FR days: rewards "
            "in the first 15 min 10 s", fontsize=12, fontweight="bold", color="#1f3a5f", va="center")
    for j, d in enumerate(days):
        cell(j + 1, 0, f"day {d} · {day_type(Tb, d)}", face="#1f3a5f", bold=True, col="white")
    cell(0, 0, "mouse (group)", face="#1f3a5f", bold=True, col="white")
    for i, (c, m) in enumerate(MICE, start=1):
        g = Tb[(Tb.cage == c) & (Tb.mouse == m)]
        grp = g.group.iloc[0] if len(g) and str(g.group.iloc[0]) else ""
        cell(0, i, f"{c} {m}", sub=grp, bold=True)
        for j, d in enumerate(days):
            r = g[g.day == d]
            if r.empty or r.task.iloc[0] == "missing":
                cell(j + 1, i, "-", face="#f2f2f2")
                continue
            r = r.iloc[0]
            ok = "Python = Saleae" if r.rewards_match else "Python ≠ Saleae!"
            if r.task == "PR":
                v = (r.pr - lo) / max(1e-6, hi - lo)
                cell(j + 1, i, f"{r.pr:.0f}", sub=f"rewards {r.rewards_python:.0f} · {ok}",
                     face=matplotlib.colors.to_hex(cmap(0.12 + 0.5 * v)), bold=True)
            elif r.task == "PASSIVE":
                # 2026-11-19 (user: "네, 같이 적기"): the partner's R and PR in the cell, so R is not read as "PR - 1"
                n_src = r.replay_source_rewards
                same = r.replay_match is True or r.replay_match == 1
                if r.partner:
                    pc, _, pcol = r.partner.partition(" ")
                    short = pcol if pc == c and pcol else r.partner        # same cage: the colour is enough
                    part = f"partner {short}: " + (f"R {n_src:.0f}" if pd.notna(n_src) else "R ?") + \
                           (f" · PR {r.partner_pr:.0f}" if pd.notna(r.partner_pr) else "")
                else:
                    part = "replay source not found"
                warn = [w for w in (None if same else ("partner ?" if pd.isna(n_src) else "≠ partner!"),
                                    None if r.rewards_match else ok) if w]
                cell(j + 1, i, f"passive · R {r.rewards_python:.0f}", sub=part, face="#ece8f8", sub2=" · ".join(warn))
            else:
                full = "" if not (r.session_s > 930) else f"untimed {r.session_s / 60:.0f}-min run: {r.rewards_python:.0f} in all"
                cell(j + 1, i, f"R {r.rewards_910s:.0f}", sub=f"{r.task} · {ok}", face="#eef1f5", sub2=full)
    # median + previous cohort
    i = len(MICE) + 1
    cell(0, i, "median (range)", bold=True, face="#fff6db")
    cell(0, i + 1, "previous cohort", sub="same day, first 15:10, reference", bold=True, face="#f6f7f9")
    for j, d in enumerate(days):
        g = Tb[Tb.day == d]
        if (g.task == "PR").any():
            x = g.pr.dropna()
            act = " (active)" if (g.task == "PASSIVE").any() else ""
            cell(j + 1, i, f"{x.median():g}", sub=f"{x.min():.0f}-{x.max():.0f}{act}", face="#fff6db", bold=True)
            pm = g.prev_pr_median.dropna()
            if len(pm):
                cell(j + 1, i + 1, f"{float(pm.iloc[0]):g}", sub=f"{float(g.prev_pr_min.dropna().iloc[0]):.0f}-"
                     f"{float(g.prev_pr_max.dropna().iloc[0]):.0f}", face="#f6f7f9")
            else:
                cell(j + 1, i + 1, "-", face="#f6f7f9")
        else:
            x = g.rewards_910s.dropna()
            cell(j + 1, i, f"R {x.median():g}", sub=f"{x.min():.0f}-{x.max():.0f}", face="#fff6db")
            cell(j + 1, i + 1, "-", sub="no PR on FR days", face="#f6f7f9")
    ax.text(0.0, nr + 0.4, "Sources: analysis\\days\\dayN\\pr_checks_dayN.csv (PR, rewards Python vs Saleae, previous "
            "cohort) · analysis\\session_metrics.csv (FR rewards, first 15 min 10 s) · tools\\pr_table_days.py",
            fontsize=8, color="#7b8794", va="center")
    fig.savefig(out, dpi=130)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--days", type=int, nargs="+")
    a = ap.parse_args(argv)
    days = a.days or all_days()
    Tb = table(days)
    tag = f"days{days[0]}-{days[-1]}" if len(days) > 1 else f"day{days[0]}"
    out_csv, out_png = DAYS_DIR / f"pr_table_{tag}.csv", DAYS_DIR / f"pr_table_{tag}.png"
    Tb.to_csv(out_csv, index=False, encoding="utf-8-sig")
    figure(Tb, days, out_png)
    bad = Tb[(Tb.task != "missing") & ~Tb.rewards_match.fillna(False).astype(bool)]
    pas = Tb[Tb.task == "PASSIVE"]
    bad_rep = pas[pas.replay_match != True]                    # noqa: E712  (NaN = source not found)
    print(f"PR table -> {out_csv}\n            {out_png}")
    for d in days:
        g = Tb[(Tb.day == d) & (Tb.task == "PR")]
        if len(g):
            print(f"  day {d}: PR median {g.pr.median():g} (range {g.pr.min():.0f}-{g.pr.max():.0f})"
                  + (" - active mice" if ((Tb.day == d) & (Tb.task == "PASSIVE")).any() else ""))
    if len(bad):
        print(f"  !! rewards Python != Saleae in {len(bad)} session(s): "
              + ", ".join(f"d{r.day} {r.cage} {r.mouse}" for r in bad.itertuples()))
    if len(bad_rep):
        print(f"  !! passive rewards != the replayed partner's log (or the log not found) in {len(bad_rep)} session(s): "
              + ", ".join(f"d{r.day} {r.cage} {r.mouse}" for r in bad_rep.itertuples()))
    return 1 if len(bad) or len(bad_rep) else 0


if __name__ == "__main__":
    raise SystemExit(main())
