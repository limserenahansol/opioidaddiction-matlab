#!/usr/bin/env python3
"""
quick_analysis.py - quick & dirty first-pass analysis of the synced sessions.

Reads <batch>_sync/**/<session>_events.csv (made by sync_sessions.py) and writes
<batch>_analysis/:
    session_metrics.csv   one row per session x window (full session, first 15:10)
    licks_long.csv        every Saleae lick (time, frame, trial, valid, ILI, ...)
    rewards_long.csv      every valve opening (commanded reward / extra opening / glitch)
    trials.csv            one row per trial (time to reward, licks in trial)
    per_minute.csv        licks and rewards per 1-min bin
    figs/*.png            overview plots
Lick metrics follow the MATLAB step04_lick pipeline: bout = ILI < 0.25 s, a gap >= 0.5 s ends a bout.
Usage:  python quick_analysis.py [--root FOLDER]
"""
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BOUT_WITHIN_IEI = 0.25     # s (same as run_lick_mega_pipeline.m)
BOUT_GAP_IEI = 0.5         # s
COMMON_WINDOW_S = 910.0    # 15 min 10 s - compare all sessions over the same time
GLITCH_MS = 5.0            # valve pulses shorter than this are electrical glitches
FULL_OPEN_MS = 100.0       # valve pulses at least this long are real openings
PSTH_WIN = (-5.0, 10.0)    # s around reward
PSTH_BIN = 0.1             # s


# ----------------------------------------------------------------------------- helpers
def parse_name(session: str, day_folder: str):
    """mouse / task / day. With the project layout day<N>/<cohort>/<animal>/ the mouse key comes from the
    folders (e.g. 9670f_green1), so the same mouse gets the same key every day."""
    m = re.match(r"^(?P<mouse>.+?)_*(?P<task>fr\d+|pr\d*)?_*day(?P<day>\d+)", session, re.I)
    mouse = m.group("mouse").strip("_") if m else session
    task = (m.group("task") or "").upper() if m else ""
    parts = [p for p in day_folder.split("/") if p]
    d = re.search(r"day(\d+)", parts[0] if parts else "", re.I)
    if d:
        day = int(d.group(1))
    elif m:
        day = int(m.group("day"))
    else:
        day = np.nan
    if len(parts) >= 3:
        mouse = f"{parts[1]}_{parts[2]}".lower()
    return mouse, task, day


def lick_bouts(on: np.ndarray, dur_s: float):
    """Port of compute_session_features_full (MATLAB step04_lick)."""
    out = dict(bout_count=0, bouts_per_min=np.nan, mean_licks_per_bout=np.nan,
               within_bout_ili_ms=np.nan, between_bout_gap_s=np.nan, burst_fraction=np.nan,
               duty_cycle=np.nan)
    if len(on) < 2:
        return out
    iei = np.diff(on)
    gap = np.r_[np.inf, iei]
    in_bout, bc, c, lick_in_bout = False, 0, 0, 0
    wb, bb = [], []
    for k in range(1, len(gap)):
        g = gap[k]
        if not in_bout and g < BOUT_WITHIN_IEI:
            in_bout, bc, c = True, bc + 1, 1
            wb.append(g)
        elif in_bout:
            if g < BOUT_WITHIN_IEI:
                c += 1
                wb.append(g)
            elif g >= BOUT_GAP_IEI:
                in_bout = False
                bb.append(g)
        if not in_bout and c > 0:
            lick_in_bout += c + 1
            c = 0
    mins = max(dur_s / 60, 1e-9)
    out.update(bout_count=bc, bouts_per_min=bc / mins,
               mean_licks_per_bout=lick_in_bout / max(bc, 1) if bc else np.nan,
               within_bout_ili_ms=float(np.median(wb)) * 1000 if wb else np.nan,
               between_bout_gap_s=float(np.median(bb)) if bb else np.nan,
               burst_fraction=float(np.mean(iei < BOUT_WITHIN_IEI)),
               duty_cycle=float(np.sum(wb)) / max(dur_s, 1e-9))
    return out


def load_session(ev_path: Path, sync_root: Path, summary: pd.DataFrame) -> dict:
    session = ev_path.name[: -len("_events.csv")]
    day_folder = ev_path.parent.relative_to(sync_root).as_posix()
    mouse, task, day = parse_name(session, day_folder)
    ev = pd.read_csv(ev_path)
    for c in ("t_from_first_pulse_s", "t_saleae_s", "t_off_saleae_s", "frame", "frame_off", "duration_ms",
              "python_elapsed_s", "python_delay_ms", "trial", "valid_or_hit", "during_valve"):
        ev[c] = pd.to_numeric(ev[c], errors="coerce")
    ev["note"] = ev["note"].fillna("").astype(str)
    ev["t"] = ev["t_from_first_pulse_s"]
    ev["t_off"] = ev["t"] + ev["duration_ms"] / 1000.0
    end_rows = ev[ev["event"] == "experiment_end"]
    t_end = float(end_rows["t"].iloc[-1]) if len(end_rows) else float(ev["t"].max())
    srow = summary[(summary["day"].astype(str) == day_folder) & (summary["session"].astype(str) == session)]
    status = srow["status"].iloc[0] if len(srow) else "?"
    vf = pd.to_numeric(srow["video_frames"], errors="coerce").iloc[0] if len(srow) else np.nan
    cp = pd.to_numeric(srow["camera_pulses"], errors="coerce").iloc[0] if len(srow) else np.nan
    frames_exact = bool(len(srow) and vf == cp)
    in_sess = ~ev["note"].str.contains("before python start|after python end", regex=True)
    licks = ev[(ev["event"] == "lick") & (ev["source"] == "saleae") & in_sess].copy()
    py_only_licks = ev[(ev["event"] == "lick") & (ev["source"] == "python")].copy()
    valves = ev[(ev["event"] == "reward_valve") & in_sess].copy()
    valves["kind"] = np.where(valves["python_elapsed_s"].notna() & (valves["duration_ms"] >= GLITCH_MS), "reward",
                     np.where(valves["duration_ms"] < GLITCH_MS, "glitch",
                     np.where(valves["duration_ms"] >= FULL_OPEN_MS, "extra_opening", "short_pulse")))
    # a reward whose matched edge was a 0-ms glitch right before the real opening: keep one "reward" per command
    cmd = valves[valves["python_elapsed_s"].notna()]
    glitch_cmd = cmd[cmd["duration_ms"] < GLITCH_MS]
    for idx, r in glitch_cmd.iterrows():
        nxt = valves[(valves.index != idx) & (valves["t"] >= r["t"]) & (valves["t"] < r["t"] + 0.05)
                     & (valves["duration_ms"] >= FULL_OPEN_MS)]
        if len(nxt):
            valves.loc[nxt.index[0], ["kind", "trial", "python_elapsed_s"]] = ["reward", r["trial"], r["python_elapsed_s"]]
            valves.loc[idx, "kind"] = "glitch"
        else:
            valves.loc[idx, "kind"] = "reward"          # keep the command even if the edge was odd
    if "req" not in ev.columns:
        ev["req"] = np.nan
    ev["req"] = pd.to_numeric(ev["req"], errors="coerce")
    trials = ev[(ev["event"] == "trial_start")][["trial", "t", "req"]].rename(columns={"t": "t_start"})
    reqs = trials["req"].dropna().unique()
    is_pr = len(reqs) > 1 or task.startswith("PR")
    if not task and len(reqs) == 1:
        task = f"FR{reqs[0]:g}"
    elif not task and is_pr:
        task = "PR"
    return dict(session=session, day_folder=day_folder, mouse=mouse, task=task, day=day, ev=ev,
                licks=licks, py_only_licks=py_only_licks, valves=valves, trials=trials, is_pr=is_pr,
                t_end=t_end, status=status, frames_exact=frames_exact, video_frames=vf, camera_pulses=cp)


def metrics(S: dict, t_max: float) -> dict:
    lk = S["licks"][S["licks"]["t"] < t_max]
    vv = S["valves"][S["valves"]["t"] < t_max]
    dur = min(S["t_end"], t_max)
    rew = vv[vv["kind"] == "reward"]
    on = np.sort(lk["t"].to_numpy())
    ili = np.diff(on)
    m = dict(mouse=S["mouse"], day=S["day"], task=S["task"], session=S["session"],
             window="first_15m10s" if t_max < 1e8 else "full_session", duration_min=dur / 60,
             sync_status=S["status"], video_frames_exact=S["frames_exact"],
             licks=len(on), licks_per_min=len(on) / max(dur / 60, 1e-9),
             licks_valid=int((lk["valid_or_hit"] == 1).sum()), licks_lockout=int((lk["valid_or_hit"] == 0).sum()),
             licks_unlabeled=int(lk["valid_or_hit"].isna().sum()),
             rewards=len(rew), rewards_per_min=len(rew) / max(dur / 60, 1e-9),
             licks_per_reward=len(on) / len(rew) if len(rew) else np.nan,
             extra_valve_openings=int((vv["kind"] == "extra_opening").sum()),
             valve_glitches=int((vv["kind"] == "glitch").sum()),
             python_only_licks=int((S["py_only_licks"]["t"] < t_max).sum()),
             first_lick_s=float(on[0]) if len(on) else np.nan,
             ili_median_ms=float(np.median(ili)) * 1000 if len(ili) else np.nan,
             ili_q10_ms=float(np.quantile(ili, 0.1)) * 1000 if len(ili) else np.nan,
             ili_q90_s=float(np.quantile(ili, 0.9)) if len(ili) else np.nan,
             ili_under_80ms_pct=float(np.mean(ili < 0.08)) * 100 if len(ili) else np.nan,   # likely double contacts
             lick_contact_median_ms=float(lk.loc[lk["during_valve"] == 0, "duration_ms"].median()))
    m.update(lick_bouts(on, dur))
    m["lick_rate_in_bout_hz"] = 1000.0 / m["within_bout_ili_ms"] if m["within_bout_ili_ms"] == m["within_bout_ili_ms"] else np.nan
    # time to reward per trial (trial start -> commanded valve opening)
    tr = S["trials"].merge(rew[["trial", "t"]].rename(columns={"t": "t_reward"}), on="trial", how="inner")
    tr = tr[tr["t_start"] < t_max]
    ttr = (tr["t_reward"] - tr["t_start"]).to_numpy()
    m["time_to_reward_median_s"] = float(np.median(ttr)) if len(ttr) else np.nan
    m["time_to_reward_iqr_s"] = float(np.subtract(*np.quantile(ttr, [0.75, 0.25]))) if len(ttr) else np.nan
    m.update(pr_scores(S, t_max))
    return m


def pr_scores(S: dict, t_max: float) -> dict:
    """PR scores. pr_breakpoint = highest trial (ratio step) reached, 1-based - the lab's official
    MATLAB Step 10 definition ('max Trial'). Also the last completed requirement and RequirementLast."""
    out = dict(is_pr=S["is_pr"], pr_breakpoint=np.nan, pr_last_completed_req=np.nan, pr_final_req=np.nan,
               pr_trials_completed=np.nan, pr_last_reward_min=np.nan)
    if not S["is_pr"]:
        return out
    tr = S["trials"][S["trials"]["t_start"] < t_max]
    rew = S["valves"][(S["valves"]["kind"] == "reward") & (S["valves"]["t"] < t_max)]
    done = tr[tr["trial"].isin(rew["trial"])]
    if len(tr):
        out["pr_breakpoint"] = int(tr["trial"].max()) + 1
        out["pr_final_req"] = float(tr.sort_values("t_start")["req"].iloc[-1])
    out["pr_trials_completed"] = len(done)
    if len(done):
        out["pr_last_completed_req"] = float(done["req"].max())
        out["pr_last_reward_min"] = float(rew["t"].max()) / 60
    return out


# ----------------------------------------------------------------------------- plots
MOUSE_COLORS = {}          # mouse -> colour, the same in every figure and every day


def mcol(mouse):
    if mouse not in MOUSE_COLORS:
        cmap = plt.get_cmap("tab10" if len(MOUSE_COLORS) < 10 else "tab20")
        MOUSE_COLORS[mouse] = cmap(len(MOUSE_COLORS) % 20)
    return MOUSE_COLORS[mouse]


def colors(n):
    cmap = plt.get_cmap("tab10")
    return [cmap(i % 10) for i in range(n)]


def label(S):
    return f"{S['mouse']} d{S['day']}" if S["day"] == S["day"] else S["session"]


def fig_rasters(SS, out):
    fig, axes = plt.subplots(len(SS), 1, figsize=(14, 1.6 * len(SS) + 0.8), sharex=True, squeeze=False)
    for ax, S in zip(axes[:, 0], SS):
        col = mcol(S["mouse"])
        lk = S["licks"]
        ax.eventplot(lk.loc[lk["valid_or_hit"] == 1, "t"] / 60, lineoffsets=1, linelengths=0.8, colors=[col], lw=0.5)
        ax.eventplot(lk.loc[lk["valid_or_hit"] != 1, "t"] / 60, lineoffsets=0, linelengths=0.8, colors=["0.6"], lw=0.5)
        rw = S["valves"][S["valves"]["kind"] == "reward"]["t"] / 60
        ax.eventplot(rw, lineoffsets=2, linelengths=0.8, colors=["k"], lw=0.8)
        ex = S["valves"][S["valves"]["kind"] == "extra_opening"]["t"] / 60
        if len(ex):
            ax.plot(ex, np.full(len(ex), 2.6), "rv", ms=5, label="extra valve opening")
        ax.axvline(COMMON_WINDOW_S / 60, color="r", ls="--", lw=0.8)
        ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["lockout lick", "counted lick", "reward"], fontsize=8)
        ax.set_title(f"{label(S)}   ({len(lk)} licks, {int((S['valves']['kind'] == 'reward').sum())} rewards, "
                     f"{S['t_end']/60:.1f} min)", fontsize=9, loc="left")
        ax.set_ylim(-0.6, 3.0)
    axes[-1, 0].set_xlabel("time from session start (min)   red dashed line = 15:10")
    fig.tight_layout(); fig.savefig(out / "01_lick_reward_raster.png", dpi=130); plt.close(fig)


def fig_cumulative(SS, out):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for S in SS:
        col = mcol(S["mouse"])
        t = np.sort(S["licks"]["t"].to_numpy()) / 60
        axes[0].step(np.r_[0, t], np.arange(len(t) + 1), where="post", color=col, label=label(S))
        r = np.sort(S["valves"][S["valves"]["kind"] == "reward"]["t"].to_numpy()) / 60
        axes[1].step(np.r_[0, r], np.arange(len(r) + 1), where="post", color=col, label=label(S))
    for ax, tt in zip(axes, ("cumulative licks", "cumulative rewards")):
        ax.axvline(COMMON_WINDOW_S / 60, color="r", ls="--", lw=0.8)
        ax.set_xlabel("time (min)"); ax.set_ylabel(tt); ax.set_title(tt); ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "02_cumulative_licks_rewards.png", dpi=130); plt.close(fig)


def per_minute(S):
    nb = int(np.ceil(S["t_end"] / 60))
    edges = np.arange(nb + 1) * 60.0
    lk, _ = np.histogram(S["licks"]["t"], edges)
    rw, _ = np.histogram(S["valves"][S["valves"]["kind"] == "reward"]["t"], edges)
    width = np.minimum(edges[1:], S["t_end"]) - edges[:-1]
    keep = width > 0    # last partial minute is kept (rates are normalised by its real length)
    return pd.DataFrame(dict(mouse=S["mouse"], day=S["day"], session=S["session"], minute=np.arange(nb)[keep] + 1,
                             licks=lk[keep], rewards=rw[keep],
                             licks_per_min=lk[keep] / (width[keep] / 60), rewards_per_min=rw[keep] / (width[keep] / 60)))


def fig_per_minute(SS, PM, out):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for S in SS:
        col = mcol(S["mouse"])
        d = PM[PM["session"] == S["session"]]
        axes[0].plot(d["minute"], d["licks_per_min"], "-o", ms=3, color=col, label=label(S))
        axes[1].plot(d["minute"], d["rewards_per_min"], "-o", ms=3, color=col, label=label(S))
    for ax, tt in zip(axes, ("licks / min", "rewards / min")):
        ax.axvline(COMMON_WINDOW_S / 60, color="r", ls="--", lw=0.8)
        ax.set_xlabel("minute of session"); ax.set_ylabel(tt); ax.set_title(tt + " (1-min bins)"); ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "03_rate_per_minute.png", dpi=130); plt.close(fig)


def fig_ili(SS, out):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    bins = np.logspace(np.log10(0.02), np.log10(60), 70)
    for S in SS:
        col = mcol(S["mouse"])
        ili = np.diff(np.sort(S["licks"]["t"].to_numpy()))
        if len(ili):
            ax.hist(ili, bins=bins, histtype="step", density=True, color=col, lw=1.5, label=f"{label(S)} (n={len(ili)})")
    ax.axvline(BOUT_WITHIN_IEI, color="k", ls=":", lw=1); ax.axvline(BOUT_GAP_IEI, color="k", ls="--", lw=1)
    ax.axvspan(0.02, 0.08, color="r", alpha=0.07)
    ax.text(0.022, ax.get_ylim()[1] * 0.92, "< 80 ms: likely\ndouble contacts", fontsize=8, color="r", va="top")
    ax.set_xscale("log"); ax.set_xlabel("inter-lick interval (s)   dotted 0.25 s = in-bout, dashed 0.5 s = bout end")
    ax.set_ylabel("density"); ax.set_title("inter-lick interval distribution"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out / "04_inter_lick_interval.png", dpi=130); plt.close(fig)


def fig_psth(SS, out):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    edges = np.arange(PSTH_WIN[0], PSTH_WIN[1] + PSTH_BIN / 2, PSTH_BIN)
    centers = edges[:-1] + PSTH_BIN / 2
    rows = []
    for S in SS:
        col = mcol(S["mouse"])
        lt = np.sort(S["licks"]["t"].to_numpy())
        rt = S["valves"][S["valves"]["kind"] == "reward"]["t"].to_numpy()
        if not len(rt):
            continue
        H = np.array([np.histogram(lt - r, edges)[0] for r in rt]) / PSTH_BIN
        mean = H.mean(0); sem = H.std(0, ddof=1) / np.sqrt(len(rt)) if len(rt) > 1 else np.zeros_like(mean)
        ax.plot(centers, mean, color=col, lw=1.4, label=f"{label(S)} (n={len(rt)} rewards)")
        ax.fill_between(centers, mean - sem, mean + sem, color=col, alpha=0.2, lw=0)
        rows.append(pd.DataFrame(dict(mouse=S["mouse"], day=S["day"], session=S["session"], t_rel_s=centers,
                                      lick_rate_hz=mean, sem=sem)))
    durs = pd.concat([S["valves"].loc[S["valves"]["kind"] == "reward", "duration_ms"] for S in SS])
    vopen = float(np.nanmedian(durs)) / 1000 if len(durs) else 0.2      # measured valve open time
    ax.axvline(0, color="k", lw=1)
    ax.axvspan(0, vopen, color="b", alpha=0.08, label=f"valve open ({vopen*1000:.0f} ms, measured)")
    ax.set_xlabel("time from reward valve opening (s)"); ax.set_ylabel("lick rate (licks/s)")
    ax.set_title("licking around reward (mean +- SEM over rewards)"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out / "05_peri_reward_licking.png", dpi=130); plt.close(fig)
    return pd.concat(rows) if rows else pd.DataFrame()


def fig_summary(M, out):
    m = M[M["window"] == "first_15m10s"].reset_index(drop=True)
    keys = [("licks_per_min", "licks / min"), ("rewards", "rewards (first 15:10)"),
            ("licks_per_reward", "licks / reward"), ("bouts_per_min", "lick bouts / min"),
            ("mean_licks_per_bout", "licks / bout"), ("lick_rate_in_bout_hz", "lick rate in bout (Hz)"),
            ("time_to_reward_median_s", "time to reward, median (s)"), ("licks_lockout", "lockout licks")]
    fig, axes = plt.subplots(2, 4, figsize=(max(15, 2.1 * len(m) * 4 / 3), 7))
    multi_day = m["day"].nunique() > 1
    labs = [str(r.mouse).replace("_", "\n") + (f"\nd{r.day}" if multi_day else "") for r in m.itertuples()]
    cols = [mcol(x) for x in m["mouse"]]
    for ax, (k, t) in zip(axes.ravel(), keys):
        ax.bar(range(len(m)), m[k], color=cols)
        for i, v in enumerate(m[k]):
            if v == v:
                ax.text(i, v, f"{v:.3g}", ha="center", va="bottom", fontsize=8)
        ax.set_xticks(range(len(m))); ax.set_xticklabels(labs, fontsize=8); ax.set_title(t, fontsize=10)
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("session summary - first 15 min 10 s of each session (same time window for every mouse)")
    fig.tight_layout(); fig.savefig(out / "06_summary_first_15m10s.png", dpi=130); plt.close(fig)


def fig_trials(SS, out):
    fig, ax = plt.subplots(figsize=(12, 4.2))
    for S in SS:
        col = mcol(S["mouse"])
        rew = S["valves"][S["valves"]["kind"] == "reward"][["trial", "t"]].rename(columns={"t": "t_reward"})
        tr = S["trials"].merge(rew, on="trial")
        ax.semilogy(tr["t_start"] / 60, tr["t_reward"] - tr["t_start"], ".", ms=4, color=col, label=label(S))
    ax.axvline(COMMON_WINDOW_S / 60, color="r", ls="--", lw=0.8)
    ax.set_xlabel("trial start (min)"); ax.set_ylabel("time to reward (s, log)")
    ax.set_title("time from trial start to reward (FR3: 3 licks)"); ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(out / "07_time_to_reward.png", dpi=130); plt.close(fig)


def fig_by_day(M, out):
    m = M[M["window"] == "first_15m10s"]
    if m["day"].nunique() < 2:
        return False
    keys = [("licks_per_min", "licks / min"), ("rewards", "rewards (first 15:10)"),
            ("licks_per_reward", "licks / reward")]
    if m["is_pr"].any():
        keys.append(("pr_breakpoint", "PR breakpoint (first 15:10, PR days)"))
    fig, axes = plt.subplots(1, len(keys), figsize=(5 * len(keys), 4.2))
    for ax, (k, t) in zip(axes, keys):
        for mouse, d in m.groupby("mouse"):
            col = mcol(mouse)
            d = d.sort_values("day")
            ax.plot(d["day"], d[k], "-o", color=col, label=mouse)
        ax.set_xlabel("day"); ax.set_title(t); ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "08_by_day.png", dpi=130); plt.close(fig)
    return True


def fig_pr(SS, M, out):
    """PR: requirement completed over time (staircase) + breakpoint per session."""
    pr = [S for S in SS if S["is_pr"]]
    if not pr:
        return False
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.6), gridspec_kw=dict(width_ratios=[2, 1]))
    for S in pr:
        col = mcol(S["mouse"])
        rew = S["valves"][S["valves"]["kind"] == "reward"][["trial", "t"]]
        d = S["trials"].merge(rew, on="trial").sort_values("t")
        axes[0].step(np.r_[0, d["t"] / 60], np.r_[0, d["req"]], where="post", color=col, label=label(S))
    axes[0].axvline(COMMON_WINDOW_S / 60, color="r", ls="--", lw=0.8)
    axes[0].set_xlabel("time (min)"); axes[0].set_ylabel("ratio requirement completed")
    axes[0].set_title("PR: requirement completed over time"); axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3)
    m = M[(M["window"] == "full_session") & (M["is_pr"])].reset_index(drop=True)
    axes[1].bar(range(len(m)), m["pr_breakpoint"], color=[mcol(x) for x in m["mouse"]], label="breakpoint (max trial)")
    axes[1].plot(range(len(m)), m["pr_last_completed_req"], "kD", label="last completed requirement")
    for i, v in enumerate(m["pr_breakpoint"]):
        axes[1].text(i, v, f"{v:g}", ha="center", va="bottom", fontsize=9)
    axes[1].set_xticks(range(len(m))); axes[1].set_xticklabels([f"{r.mouse}\nd{r.day}" for r in m.itertuples()], fontsize=8)
    axes[1].set_title("PR breakpoint (full session)"); axes[1].legend(fontsize=7); axes[1].grid(axis="y", alpha=0.3)
    fig.tight_layout(); fig.savefig(out / "09_pr_breakpoint.png", dpi=130); plt.close(fig)
    return True


# ----------------------------------------------------------------------------- pupil
PUPIL_DIAM_COLS = ("smootheddiameter_px", "diameter_pixels", "diameter_px", "diameter", "pupil_diameter")


ONE_CLICK_CSV = "pupil_diameter_FINAL_U-Net_results_CLEAN.csv"


def find_pupil(pupil_root: Path, day_folder: str, session: str):
    """1) <project>/pupil/.../<session>_pupil*.csv
       2) the MATLAB one-click output: <project>/pupil/<dayN>/<cage>/<mouse>/pupil_outputs_general/
          pupil_diameter_FINAL_U-Net_results_CLEAN.csv  (one session per mouse-day folder)"""
    if pupil_root is None or not pupil_root.is_dir():
        return None
    hits = [p for p in pupil_root.rglob("*.csv")
            if p.stem.lower().startswith(session.lower() + "_pupil") or p.parent.name.lower() == session.lower()]
    if hits:
        return sorted(hits)[0]
    d = pupil_root / day_folder
    if d.is_dir():
        oc = sorted(d.rglob(ONE_CLICK_CSV), key=lambda p: p.stat().st_mtime)
        if oc:
            return oc[-1]
    return None


def pupil_grade(pupil_root: Path, day_folder: str):
    """(grade, use_for_analysis, reason) from pupil/<dayN>/pupil_report_<dayN>.csv written by tools/pupil_py/
    run_pupil_day.py; (None, True, '') when there is no report (e.g. MATLAB one-click output)."""
    parts = [p for p in day_folder.split("/") if p]
    if pupil_root is None or len(parts) < 3:
        return None, True, ""
    # pupil/TST package (RUN_ALL_18DAYS): its qc.json says how often no pupil was found (2026-11-14)
    qj = pupil_root / parts[0] / parts[1] / parts[2] / "pupil_outputs_general" / "qc.json"
    if qj.exists():
        try:
            q = json.loads(qj.read_text(encoding="utf-8"))
            nb = float(q.get("fracNoBlob", 0) or 0)
            if nb > 0.5:
                return "no_pupil", False, f"package QC: no pupil found in {100 * nb:.0f} % of frames"
            return ("ok" if nb <= 0.1 else "partial"), True, f"package QC: no pupil in {100 * nb:.0f} % of frames"
        except Exception:
            pass
    rep = pupil_root / parts[0] / f"pupil_report_{parts[0]}.csv"
    if not rep.exists():
        return None, True, ""
    try:
        R = pd.read_csv(rep, dtype={"cage": str, "mouse": str}, encoding="utf-8-sig")
    except Exception:
        return None, True, ""
    r = R[(R["cage"].astype(str) == parts[1]) & (R["mouse"].astype(str).str.lower() == parts[2].lower())]
    if r.empty:
        return None, True, ""
    r = r.iloc[-1]
    use = str(r.get("use_for_analysis", "True")).strip().lower() in ("true", "1", "yes")
    reason = str(r.get("grade_reason", "") if pd.notna(r.get("grade_reason", "")) else "")
    return str(r.get("grade", "")), use, reason


def pupil_metrics(S: dict, P: pd.DataFrame, t_max: float) -> dict:
    """median diameter and the reward response (% change 0.5-2.5 s after reward vs -1..0 s before)."""
    out = dict(pupil_median_px=np.nan, pupil_reward_response_pct=np.nan, pupil_rewards_used=0)
    if P is None or not len(P):
        return out
    tt, dd = P["t_s"].to_numpy(), P["diameter"].to_numpy()
    ok = np.isfinite(tt) & np.isfinite(dd) & (tt < t_max)
    if ok.sum() < 100:
        return out
    tt, dd = tt[ok], dd[ok]
    out["pupil_median_px"] = float(np.median(dd))
    resp = []
    for r in S["valves"][(S["valves"]["kind"] == "reward") & (S["valves"]["t"] < t_max - 3)]["t"]:
        base = dd[(tt >= r - 1) & (tt < r)]
        post = dd[(tt >= r + 0.5) & (tt < r + 2.5)]
        if len(base) > 5 and len(post) > 5 and np.mean(base) > 0:
            resp.append((np.mean(post) - np.mean(base)) / np.mean(base) * 100)
    if resp:
        out["pupil_reward_response_pct"] = float(np.mean(resp))
        out["pupil_rewards_used"] = len(resp)
    return out


def load_pupil(S: dict, pupil_csv: Path, frames_csv: Path):
    """Pupil rows = video frames in order (row 0 = video frame 0 = first camera pulse after the video offset).
    Aligned from the START: today's checks showed extra/missing frames are at the END of the videos,
    so the old MATLAB end-alignment would be wrong."""
    P = pd.read_csv(pupil_csv)
    low = {c.lower(): c for c in P.columns}
    dcol = next((low[c] for c in PUPIL_DIAM_COLS if c in low), None) or \
        next((c for c in P.columns if "diam" in c.lower()), None)
    if dcol is None:
        return None, f"no diameter column in {pupil_csv.name} (columns: {list(P.columns)[:8]})"
    fcol = low.get("frame")
    frame = pd.to_numeric(P[fcol], errors="coerce").to_numpy() if fcol else np.arange(len(P))
    if fcol and np.nanmin(frame) == 1:
        frame = frame - 1                     # MATLAB-style 1-based frame numbers
    F = pd.read_csv(frames_csv)
    F = F[pd.to_numeric(F["frame"], errors="coerce").notna()]
    t_of = dict(zip(F["frame"].astype(int), F["t_from_first_pulse_s"].astype(float)))
    t = np.array([t_of.get(int(f), np.nan) if f == f else np.nan for f in frame])
    diam = pd.to_numeric(P[dcol], errors="coerce").to_numpy(float)
    raw = low.get("diameter_pixels")
    if raw is not None and raw != dcol:
        # MATLAB/port CLEAN files fill the frames before the first detection by linear EXTRAPOLATION
        # (fillmissing 'linear'), i.e. made-up values -> drop them here
        diam = np.where(np.isfinite(pd.to_numeric(P[raw], errors="coerce").to_numpy(float)), diam, np.nan)
    out = pd.DataFrame(dict(frame=frame, t_s=t, diameter=diam))
    n_vid = len(F)
    note = f"{len(P)} pupil rows, {n_vid} video frames with a time"
    if len(P) != n_vid:
        note += f" - MISMATCH ({len(P) - n_vid:+d}); aligned from the first frame"
    return out, note


def fig_pupil(SP, out):
    if not SP:
        return
    fig, axes = plt.subplots(len(SP), 1, figsize=(14, 2.2 * len(SP) + 0.6), sharex=True, squeeze=False)
    for ax, (S, P) in zip(axes[:, 0], SP):
        col = mcol(S["mouse"])
        ax.plot(P["t_s"] / 60, P["diameter"], color=col, lw=0.5)
        r = S["valves"][S["valves"]["kind"] == "reward"]["t"] / 60
        ax.plot(r, np.full(len(r), np.nanpercentile(P["diameter"], 99)), "k|", ms=8, label="reward")
        ax.set_ylabel("pupil (px)"); ax.set_title(label(S), fontsize=9, loc="left")
    axes[-1, 0].set_xlabel("time from session start (min)")
    fig.tight_layout(); fig.savefig(out / "10_pupil_timeline.png", dpi=130); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    edges = np.arange(PSTH_WIN[0], PSTH_WIN[1] + 1e-9, 0.1)
    for S, P in SP:
        col = mcol(S["mouse"])
        tt, dd = P["t_s"].to_numpy(), P["diameter"].to_numpy()
        ok = np.isfinite(tt) & np.isfinite(dd)
        tt, dd = tt[ok], dd[ok]
        rows = []
        for r in S["valves"][S["valves"]["kind"] == "reward"]["t"]:
            y = np.interp(r + edges, tt, dd, left=np.nan, right=np.nan)
            base = np.nanmean(y[(edges >= -1) & (edges < 0)])
            if np.isfinite(base) and base > 0:
                rows.append((y - base) / base * 100)
        if rows:
            H = np.array(rows)
            mean = np.nanmean(H, 0); sem = np.nanstd(H, 0) / np.sqrt(np.sum(np.isfinite(H), 0).clip(1))
            ax.plot(edges, mean, color=col, label=f"{label(S)} (n={len(rows)})")
            ax.fill_between(edges, mean - sem, mean + sem, color=col, alpha=0.2, lw=0)
    ax.axvline(0, color="k", lw=1); ax.axhline(0, color="0.5", lw=0.5)
    ax.set_xlabel("time from reward (s)"); ax.set_ylabel("pupil change from -1..0 s baseline (%)")
    ax.set_title("pupil around reward (mean +- SEM)"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out / "11_peri_reward_pupil.png", dpi=130); plt.close(fig)


def fig_dashboard(out, title="", name="00_dashboard.png"):
    import matplotlib.image as mpimg
    names = ["06_summary_first_15m10s.png", "01_lick_reward_raster.png", "02_cumulative_licks_rewards.png",
             "03_rate_per_minute.png", "05_peri_reward_licking.png", "04_inter_lick_interval.png",
             "09_pr_breakpoint.png", "11_peri_reward_pupil.png"]
    names = [n for n in names if (out / n).exists()]
    rows = int(np.ceil(len(names) / 2))
    fig, axes = plt.subplots(rows, 2, figsize=(20, 6 * rows), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, n in zip(axes.ravel(), names):
        ax.imshow(mpimg.imread(out / n))
    if title:
        fig.suptitle(title, fontsize=20)
    fig.tight_layout(); fig.savefig(out / name, dpi=80); plt.close(fig)


def fig_across_days(MD, out):
    """Across-day spaghetti (first 15 min 10 s of every session): one line per mouse."""
    keys = [("licks_per_min", "licks / min"), ("rewards", "rewards"), ("licks_per_reward", "licks / reward"),
            ("bouts_per_min", "lick bouts / min"), ("mean_licks_per_bout", "licks / bout"),
            ("time_to_reward_median_s", "time to reward (s)")]
    if MD["pr_breakpoint"].notna().any():
        keys.append(("pr_breakpoint", "PR breakpoint"))
    if MD["pupil_reward_response_pct"].notna().any():
        keys.append(("pupil_reward_response_pct", "pupil reward response (%)"))
    cols = 4
    rows = int(np.ceil(len(keys) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), squeeze=False)
    for ax in axes.ravel()[len(keys):]:
        ax.axis("off")
    for ax, (k, t) in zip(axes.ravel(), keys):
        for mouse, d in MD.groupby("mouse"):
            d = d.sort_values("day")
            ls = "-" if str(mouse).split("_")[0].endswith("f") else "--"
            ax.plot(d["day"], d[k], ls, marker="o", color=mcol(mouse), label=mouse)
        ax.set_xlabel("day"); ax.set_title(t); ax.grid(alpha=0.3)
        ax.set_xticks(sorted(MD["day"].dropna().unique()))
    axes[0, 0].legend(fontsize=7, ncol=2)
    fig.suptitle("across days - first 15 min 10 s of each session (solid = female, dashed = male)")
    fig.tight_layout(); fig.savefig(out / "across_days_metrics.png", dpi=110); plt.close(fig)


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    args = ap.parse_args(argv)
    projects = []
    if (args.root / "sync" / "sync_summary.csv").exists():        # project layout
        # pupil: the user's pupil/TST package output (pupil_tst_final) when it exists, else the older pupil folder
        proot = args.root / "pupil_tst_final" if (args.root / "pupil_tst_final").is_dir() else args.root / "pupil"
        projects.append((args.root.name, args.root / "sync", args.root / "analysis", proot, True))
    for sync_root in sorted(args.root.glob("*_sync")):             # older <batch>_sync layout
        prefix = sync_root.name[: -len("_sync")]
        projects.append((prefix, sync_root, args.root / f"{prefix}_analysis", args.root / f"{prefix}_pupil", False))
    for prefix, sync_root, out, pupil_root, is_project in projects:
        sp = sync_root / "sync_summary.csv"
        if not sp.exists():
            continue
        summary = pd.read_csv(sp, dtype=str, encoding="utf-8-sig").fillna("")
        SS = []
        for ev_path in sorted(sync_root.rglob("*_events.csv")):
            S = load_session(ev_path, sync_root, summary)
            if S["status"] == "FAIL":
                print(f"  skip {S['session']} (sync FAIL)")
                continue
            SS.append(S)
        if not SS:
            continue
        SS.sort(key=lambda S: (S["day"], S["mouse"]))
        figs = out / "figs"
        out.mkdir(parents=True, exist_ok=True)
        print(f"=== {prefix}: {len(SS)} session(s) -> {out}")

        # ---- pupil (optional): <project>/pupil/.../<session>_pupil.csv
        SP, pupil_notes, pupil_grades = [], {}, {}
        for S in SS:
            pc = find_pupil(pupil_root, S["day_folder"], S["session"])
            if pc is None:
                continue
            g, use, why = pupil_grade(pupil_root, S["day_folder"])
            if g is not None:
                pupil_grades[S["session"]] = g
            if not use:
                pupil_notes[S["session"]] = f"pupil NOT used: grade '{g}' ({why})"
                print(f"  pupil {S['session']}: {pupil_notes[S['session']]}")
                continue
            fc = sync_root / S["day_folder"] / f"{S['session']}_frames.csv"
            try:
                P, note = load_pupil(S, pc, fc)
            except Exception as exc:
                P, note = None, f"could not read {pc.name}: {exc}"
            if g is not None:
                note = f"grade '{g}'" + (f" ({why})" if why else "") + "; " + note
            pupil_notes[S["session"]] = note
            print(f"  pupil {S['session']}: {note}")
            if P is not None:
                SP.append((S, P))
                (out / "pupil_aligned").mkdir(exist_ok=True)
                P.to_csv(out / "pupil_aligned" / f"{S['session']}_pupil_aligned.csv", index=False, float_format="%.6g")

        PUP = {S["session"]: P for S, P in SP}
        mrows = []
        for tmax in (1e9, COMMON_WINDOW_S):
            for S in SS:
                r = metrics(S, tmax)
                r.update(pupil_metrics(S, PUP.get(S["session"]), tmax))
                r["pupil_grade"] = pupil_grades.get(S["session"], "")
                parts = S["day_folder"].split("/")
                r["cage"] = parts[1] if len(parts) >= 3 else ""
                r["sex"] = {"f": "F", "m": "M"}.get(r["cage"][-1:].lower(), "") if r["cage"] else ""
                mrows.append(r)
        M = pd.DataFrame(mrows)
        M = M.sort_values(["window", "day", "mouse"])
        M.to_csv(out / "session_metrics.csv", index=False, float_format="%.4g")

        L = []
        for S in SS:
            lk = S["licks"].sort_values("t").copy()
            lk["ili_s"] = lk["t"].diff()
            lk.insert(0, "session", S["session"]); lk.insert(0, "day", S["day"]); lk.insert(0, "mouse", S["mouse"])
            lk["frame_reliable"] = S["frames_exact"]
            L.append(lk[["mouse", "day", "session", "t", "t_off", "duration_ms", "ili_s", "trial", "valid_or_hit",
                         "during_valve", "frame", "frame_off", "frame_reliable", "t_saleae_s"]]
                     .rename(columns={"t": "t_s", "t_off": "t_off_s", "valid_or_hit": "counted"}))
        pd.concat(L).to_csv(out / "licks_long.csv", index=False, float_format="%.6g")

        R = []
        for S in SS:
            v = S["valves"].sort_values("t").copy()
            v.insert(0, "session", S["session"]); v.insert(0, "day", S["day"]); v.insert(0, "mouse", S["mouse"])
            v["frame_reliable"] = S["frames_exact"]
            R.append(v[["mouse", "day", "session", "kind", "t", "duration_ms", "trial", "frame", "frame_off",
                        "frame_reliable", "t_saleae_s"]].rename(columns={"t": "t_s"}))
        pd.concat(R).to_csv(out / "rewards_long.csv", index=False, float_format="%.6g")

        T = []
        for S in SS:
            rew = S["valves"][S["valves"]["kind"] == "reward"][["trial", "t", "frame"]].rename(
                columns={"t": "t_reward_s", "frame": "reward_frame"})
            tr = S["trials"].rename(columns={"t_start": "t_start_s"}).merge(rew, on="trial", how="left")
            lk = S["licks"]
            tr["licks_counted"] = [int(((lk["trial"] == k) & (lk["valid_or_hit"] == 1)).sum()) for k in tr["trial"]]
            tr["licks_lockout"] = [int(((lk["trial"] == k) & (lk["valid_or_hit"] == 0)).sum()) for k in tr["trial"]]
            tr["time_to_reward_s"] = tr["t_reward_s"] - tr["t_start_s"]
            tr.insert(0, "session", S["session"]); tr.insert(0, "day", S["day"]); tr.insert(0, "mouse", S["mouse"])
            T.append(tr)
        pd.concat(T).to_csv(out / "trials.csv", index=False, float_format="%.6g")

        PM = pd.concat([per_minute(S) for S in SS])
        PM.to_csv(out / "per_minute.csv", index=False, float_format="%.4g")

        # ---- one summary folder per day: analysis/days/day<N>/
        import shutil as _sh
        if is_project and figs.exists():
            _sh.rmtree(figs, ignore_errors=True)           # old all-session figures (now per day)
        PS = []
        for d in sorted({S["day"] for S in SS}):
            SSd = [S for S in SS if S["day"] == d]
            Md = M[M["day"] == d]
            dd = out / "days" / f"day{d}"
            fd = dd / "figs"
            fd.mkdir(parents=True, exist_ok=True)
            fig_summary(Md, fd)
            fig_rasters(SSd, fd)
            fig_cumulative(SSd, fd)
            fig_per_minute(SSd, PM[PM["day"] == d], fd)
            fig_ili(SSd, fd)
            P = fig_psth(SSd, fd)
            if len(P):
                PS.append(P)
            fig_trials(SSd, fd)
            fig_pr(SSd, Md, fd)
            fig_pupil([(S, P_) for S, P_ in SP if S["day"] == d], fd)
            fig_dashboard(fd, title=f"day {d}  -  {len(SSd)} mice  (first 15 min 10 s unless noted)",
                          name="00_day_dashboard.png")
            Md.to_csv(dd / f"day{d}_metrics.csv", index=False, float_format="%.4g")
            write_day_summary(dd / f"day{d}_summary.txt", d, Md, SSd, pupil_notes)
            print(f"  day {d}: {len(SSd)} sessions -> {dd}")
        if PS:
            pd.concat(PS).to_csv(out / "peri_reward_psth.csv", index=False, float_format="%.5g")
        # ---- across days: one row per mouse x day (first 15 min 10 s)
        MD = M[M["window"] == "first_15m10s"].copy()
        master_cols = ["day", "cage", "sex", "mouse", "task", "session", "duration_min", "licks", "licks_per_min",
                       "rewards", "rewards_per_min", "licks_per_reward", "licks_valid", "licks_lockout",
                       "bout_count", "bouts_per_min", "mean_licks_per_bout", "within_bout_ili_ms",
                       "lick_rate_in_bout_hz", "time_to_reward_median_s", "ili_under_80ms_pct",
                       "pr_breakpoint", "pr_last_completed_req", "pr_final_req",
                       "pupil_median_px", "pupil_reward_response_pct", "pupil_grade", "sync_status",
                       "video_frames_exact"]
        MD[master_cols].sort_values(["day", "cage", "mouse"]).to_csv(out / "master_day_level.csv", index=False,
                                                                     float_format="%.4g")
        fa = out / "figs_across_days"
        fa.mkdir(exist_ok=True)
        if MD["day"].nunique() >= 2:
            fig_across_days(MD, fa)
        else:
            (fa / "README.txt").write_text("Across-day figures appear here once there are sessions from 2 or more days.\n",
                                           encoding="utf-8")
        if is_project:
            write_session_log(args.root, sync_root, summary, M, pupil_notes)
        cols = ["mouse", "day", "task", "duration_min", "licks", "licks_per_min", "rewards", "licks_per_reward",
                "bouts_per_min", "mean_licks_per_bout", "time_to_reward_median_s", "pr_breakpoint"]
        with pd.option_context("display.width", 220, "display.max_columns", 30):
            for w in ("full_session", "first_15m10s"):
                print(f"\n  [{w}]")
                print(M[M["window"] == w][cols].round(2).to_string(index=False))
    return 0


def write_day_summary(path: Path, d, Md: pd.DataFrame, SSd, pupil_notes: dict):
    """Plain-text one-page summary of a day (all mice)."""
    m = Md[Md["window"] == "first_15m10s"].sort_values(["cage", "mouse"])
    cols = ["mouse", "task", "duration_min", "licks", "licks_per_min", "rewards", "licks_per_reward",
            "bouts_per_min", "mean_licks_per_bout", "time_to_reward_median_s", "pr_breakpoint",
            "pupil_reward_response_pct", "pupil_grade", "sync_status"]
    cols = [c for c in cols if c in Md.columns]
    full = Md[Md["window"] == "full_session"].set_index("mouse")
    lines = [f"DAY {d} SUMMARY  ({len(SSd)} mice)   made {datetime.now():%Y-%m-%d %H:%M}",
             "All numbers are for the first 15 min 10 s of each session (same window for every mouse).", ""]
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        lines.append(m[cols].round(2).to_string(index=False))
    lines += ["", "Group means (first 15:10):"]
    for k in ("licks_per_min", "rewards", "licks_per_reward", "pr_breakpoint", "pupil_reward_response_pct"):
        v = m[k].dropna()
        if len(v):
            lines.append(f"  {k:28s} mean {v.mean():8.2f}   SD {v.std(ddof=1) if len(v) > 1 else 0:7.2f}   n={len(v)}")
    for sex, g in m.groupby("sex"):
        if sex:
            lines.append(f"  {sex}: licks/min {g['licks_per_min'].mean():.1f}, rewards {g['rewards'].mean():.1f} (n={len(g)})")
    long = [mo for mo in m["mouse"] if mo in full.index and full.loc[mo, "duration_min"] > 15.5]
    if long:
        lines += ["", "Sessions longer than 15:10 (full-session numbers are in day*_metrics.csv): " + ", ".join(long)]
    notes = [f"  {S['mouse']}: {pupil_notes[S['session']]}" for S in SSd if S["session"] in pupil_notes]
    lines += ["", "Pupil: " + ("" if notes else "no pupil data for this day (put one-click output in pupil\\)")] + notes
    warn = [f"  {r.mouse}: {r.sync_status}" for r in m.itertuples() if r.sync_status != "OK"]
    if warn:
        lines += ["", "Sync status not OK (see sync\\...\\_report.txt):"] + warn
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_session_log(root: Path, sync_root: Path, summary: pd.DataFrame, M: pd.DataFrame, pupil_notes: dict):
    """<project>/session_log.csv - one row per session found in python/ (synced or not), to track everything."""
    try:
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import sync_sessions as ss
        batch = ss.discover(root)[0]
        found = ss.sessions_in(batch, root)
    except Exception as exc:
        print(f"  (session_log skipped: {exc})")
        return
    full = M[M["window"] == "full_session"].set_index("session")
    first = M[M["window"] == "first_15m10s"].set_index("session")
    pupil_root = root / "pupil_tst_final" if (root / "pupil_tst_final").is_dir() else root / "pupil"
    rows = []
    for s in found:
        name, day = s["name"], s["day"]
        mouse, task, dayn = parse_name(name, day)
        start, dur = "", np.nan
        try:
            import json
            with open(s["python"], encoding="utf-8-sig") as f:
                lines = [ln for ln in f if ln.strip()]
            start = json.loads(lines[0]).get("t", "")[:19].replace("T", " ")
            dur = round(float(json.loads(lines[-1]).get("elapsed", np.nan)) / 60, 2)
        except Exception:
            pass
        sm = summary[(summary["day"] == day) & (summary["session"] == name)]
        r = dict(day=dayn, mouse=mouse, task=task, session=name, start=start, python_min=dur,
                 python="yes", saleae_csv="yes" if s.get("saleae_csv") else "MISSING",
                 sal="yes" if s.get("sal") else "missing", video="yes" if s.get("video") else "MISSING",
                 pupil="yes" if find_pupil(pupil_root, day, name) else "",
                 sync_status=sm["status"].iloc[0] if len(sm) else ("SKIPPED" if not s.get("saleae_csv") else "not run"),
                 video_frames=sm["video_frames"].iloc[0] if len(sm) else "",
                 camera_pulses=sm["camera_pulses"].iloc[0] if len(sm) else "")
        if name in full.index:
            f1, f2 = full.loc[name], first.loc[name]
            r.update(licks=int(f1["licks"]), rewards=int(f1["rewards"]),
                     licks_per_min_15m10s=round(float(f2["licks_per_min"]), 2), rewards_15m10s=int(f2["rewards"]),
                     pr_breakpoint=f1["pr_breakpoint"])
        r["notes"] = "; ".join(x for x in [sm["problems"].iloc[0] if len(sm) else "", pupil_notes.get(name, "")] if x)
        rows.append(r)
    log = pd.DataFrame(rows).sort_values(["day", "mouse", "start"])
    try:
        log.to_csv(root / "session_log.csv", index=False, encoding="utf-8-sig")
        print(f"  -> session log: {root / 'session_log.csv'} ({len(log)} sessions)")
    except PermissionError:
        print("  !! could not write session_log.csv - close it in Excel and run again")


if __name__ == "__main__":
    raise SystemExit(main())
