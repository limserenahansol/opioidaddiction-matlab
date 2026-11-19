#!/usr/bin/env python3
r"""
explore_days1to5.py - EXPLORATORY tests of the post-hoc ideas on the saved day 1-5 outputs (9669m / 9670f).

Made 2026-11-19 (Claude) after the user's "test your ideas". The ideas and papers are in
analysis/posthoc_plan/NEW_IDEAS_20261119.md (the user: "you can try this kind of posthoc analysis as exploratory
analysis but you can refer to other papers").
Reads ONLY saved per-day outputs (never re-runs a daily step) and writes only analysis/exploratory/<YYYYMMDD>/.
Every figure and table says EXPLORATORY. Nothing here goes into the daily tables, Step 02, the routine PPT or the
plan's main analyses.

  E1 (idea I03)  PR as a function of licks (+1 schedule: n rewards <-> n(n+1)/2 counted licks); sessions still
                 licking at 910 s (censored); quit time = last lick before a final pause >= 30 s.
  E2 (I02)       clock time of the session start vs PR and pupil, days 3-5 (cage = sex = time block).
  E3 (L04)       lick clusters (ILI < 0.5 s; Dwyer 2012) and bouts (gap >= 1 s): number ("seeking") and size
                 ("liking"), days 1-5.
  E4 (L05)       hot-plate first-response latency from the day-5 hand scores, next to tail immersion.
  E5 (I04)       day-to-day reliability over days 3-5: ICC(3,1) (consistency) and ICC(2,1) (agreement), 95 % CI.
  E6 (L06)       engaged / disengaged states: 2-state Poisson HMM on 1-s lick counts, days 3-5.
  E7 (I06/I08)   pupil GLM, days 3-5: reward kernel with and without lick regressors (does licking explain the
                 reward-locked dilation? Yamada & Toda 2022). Whole-eye mice apart; noisy sessions also left out.

Inputs: analysis/master_day_level.csv, licks_long.csv, rewards_long.csv, days/day<N>/pr_checks_day<N>.csv,
cohort_manifest.csv, assays/day5/<cage>/<mouse>/ManualScoringResults_*HOT*.mat, assays/tail_immersion_all.csv,
the pupil traces the tables use (pupil_check_examples.session = the same frames-out rules as the pupil-check
pictures) and analysis/pupil_check/pupil_check_verdicts.csv.
Runs at Windows idle priority with background I/O and one thread.
usage: python explore_days1to5.py [--out-date YYYYMMDD] [--only E1,E3,...]
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import ctypes  # noqa: E402
import datetime as dt  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import traceback  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TOOLS = Path(__file__).resolve().parent.parent
P = TOOLS.parent
A = P / "analysis"
sys.path.insert(0, str(TOOLS))

WIN_S = 910.0
CLUSTER_ILI_S = 0.5
BOUT_GAP_S = 1.0
CRIT_PAUSE_S = 30.0
LAST_S = 60.0
PR_DAYS = (3, 4, 5)
ALL_DAYS = (1, 2, 3, 4, 5)
CAM_HZ = 30.26
TAG = "EXPLORATORY - days 1-5, n = 8 mice"
TAX = np.arange(-int(round(2 * CAM_HZ)), int(round(6 * CAM_HZ)))     # frames around a reward (-2 .. +6 s)
GCOL = {"Active": "#c0392b", "Passive": "#2471a3"}


def background():
    try:
        k = ctypes.windll.kernel32
        h = k.GetCurrentProcess()
        k.SetPriorityClass(h, 0x40)            # IDLE_PRIORITY_CLASS
        k.SetPriorityClass(h, 0x00100000)      # PROCESS_MODE_BACKGROUND_BEGIN: very low I/O priority
    except Exception:  # noqa: BLE001
        pass


def manifest() -> pd.DataFrame:
    m = pd.read_csv(P / "cohort_manifest.csv")
    m["mouse_id"] = m["cage_folder"].str.lower() + "_" + m["mouse"].str.lower()
    return m.set_index("mouse_id")


def lab(mid: str) -> str:
    c, m = mid.split("_", 1)
    return f"{c} {m}"


def mouse_colors(ids):
    cm = plt.get_cmap("tab10")
    return {m: cm(i % 10) for i, m in enumerate(sorted(ids))}


def spearman(x, y):
    from scipy.stats import spearmanr
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4:
        return np.nan, np.nan, int(ok.sum())
    r, p = spearmanr(x[ok], y[ok])
    return float(r), float(p), int(ok.sum())


def save(fig, path: Path):
    fig.text(0.995, 0.005, TAG, ha="right", va="bottom", fontsize=7, color="0.45")
    fig.savefig(path, dpi=130)
    plt.close(fig)


def readme(d: Path, title: str, lines: list):
    (d / "README.md").write_text(f"# {title}\n\n**{TAG}.** Made by tools/exploratory/explore_days1to5.py "
                                 f"({dt.datetime.now():%Y-%m-%d %H:%M}).\n\n" + "\n".join(lines) + "\n",
                                 encoding="utf-8")


def runs(t: np.ndarray, gap: float) -> np.ndarray:
    """sizes of the runs of event times split where the interval is >= gap"""
    if len(t) == 0:
        return np.array([], int)
    br = np.flatnonzero(np.diff(t) >= gap)
    return np.diff(np.r_[0, br + 1, len(t)])


def tri_n(v: float) -> int:
    """largest n with n(n+1)/2 <= v: the rewards the +1 schedule gives for v counted licks"""
    return int(math.floor((math.sqrt(8 * v + 1) - 1) / 2)) if v >= 0 else -1


# ------------------------------------------------------------------------------------------------ E1
def e1(M, L, R, man, out: Path, res: dict):
    d = out / "E1_pr_licks"
    d.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in M[M.day.isin(PR_DAYS)].itertuples():
        lk = L[(L.session == s.session) & (L.mouse == s.mouse)]
        w = lk[(lk.t_s >= 0) & (lk.t_s <= WIN_S)]
        t = np.sort(w.t_s.to_numpy(float))
        v = int((w.counted == 1).sum())
        rw = np.sort(R[(R.session == s.session) & (R.mouse == s.mouse)].t_s.to_numpy(float))
        rw = rw[(rw >= 0) & (rw <= WIN_S)]
        n = len(rw)
        iri = np.diff(rw)
        last = float(t[-1]) if len(t) else np.nan
        quit_ = bool(len(t)) and (WIN_S - last >= CRIT_PAUSE_S)
        rows.append(dict(
            day=int(s.day), mouse=s.mouse, group=man.loc[s.mouse, "group"], sex=man.loc[s.mouse, "sex"],
            session=s.session, rewards_910=n, PR_RL_910=n + 1, licks_910=len(t), counted_licks_910=v,
            master_rewards=s.rewards, master_licks_valid=s.licks_valid,
            rewards_from_counted=tri_n(v), rewards_from_all_licks=tri_n(len(t)),
            progress_after_last_reward=v - n * (n + 1) // 2,
            last_lick_s=last, final_gap_s=(WIN_S - rw[-1]) if n else np.nan,
            recent_iri_median_s=float(np.median(iri[-5:])) if len(iri) else np.nan,
            licks_last60s=int(np.sum(t >= WIN_S - LAST_S)),
            quit=quit_, quit_min=(last / 60) if quit_ else np.nan, censored_at_910s=not quit_))
    E = pd.DataFrame(rows)
    E["rule_holds"] = (E.rewards_from_counted == E.rewards_910)
    E["matches_master"] = (E.rewards_910 == E.master_rewards) & (E.counted_licks_910 == E.master_licks_valid)
    E.to_csv(d / "E1_pr_licks.csv", index=False)
    r_all = spearman(E.licks_910, E.rewards_910)
    r_cnt = spearman(E.counted_licks_910, E.rewards_910)
    res["E1"] = dict(n_sessions=len(E), rule_holds=int(E.rule_holds.sum()),
                     from_all_licks_within1=int((abs(E.rewards_from_all_licks - E.rewards_910) <= 1).sum()),
                     matches_master=int(E.matches_master.sum()),
                     rho_licks_rewards=r_all[0], rho_counted_rewards=r_cnt[0],
                     censored=int(E.censored_at_910s.sum()),
                     censored_by_day={int(k): int(v) for k, v in E.groupby("day").censored_at_910s.sum().items()},
                     quit_min_median=float(np.nanmedian(E.quit_min)),
                     mismatch_rows=E[~E.rule_holds | ~E.matches_master][
                         ["day", "mouse", "rewards_910", "counted_licks_910", "rewards_from_counted",
                          "master_rewards", "master_licks_valid"]].to_dict("records"))
    col = mouse_colors(E.mouse.unique())
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    xs = np.arange(0, E.counted_licks_910.max() + 30)
    ax[0].plot(xs, [tri_n(x) for x in xs], color="0.6", lw=1, label="rule: n rewards = n(n+1)/2 counted licks")
    for r in E.itertuples():
        ax[0].scatter(r.counted_licks_910, r.rewards_910, color=col[r.mouse], marker="osD"[r.day - 3], s=40,
                      edgecolor="k" if r.censored_at_910s else "none")
    ax[0].set_xlabel("counted licks in the first 910 s")
    ax[0].set_ylabel("rewards in the first 910 s (PR = rewards + 1)")
    ax[0].set_title(f"E1a  PR follows the lick count  (rule holds {int(E.rule_holds.sum())}/{len(E)})", fontsize=10)
    ax[0].legend(fontsize=8, loc="lower right")
    ax[0].text(0.02, 0.97, "marker = day 3 o / 4 s / 5 D; black edge = still licking at 910 s", transform=ax[0].transAxes,
               fontsize=7.5, va="top")
    ids = sorted(E.mouse.unique())
    for i, m in enumerate(ids):
        for r in E[E.mouse == m].itertuples():
            y = r.quit_min if r.quit else WIN_S / 60
            ax[1].scatter(r.day + (i - 3.5) * 0.07, y, color=col[m], marker="o" if r.quit else "^", s=45)
    ax[1].axhline(WIN_S / 60, color="0.5", ls=":", lw=1)
    ax[1].set_xticks(PR_DAYS)
    ax[1].set_xlabel("day")
    ax[1].set_ylabel("quit time (min; triangle = still licking at 15:10)")
    ax[1].set_title(f"E1b  quit = last lick before a final pause >= {CRIT_PAUSE_S:.0f} s", fontsize=10)
    for m in ids:
        ax[1].scatter([], [], color=col[m], label=lab(m))
    ax[1].legend(fontsize=7, ncol=2, loc="lower left")
    fig.tight_layout()
    save(fig, d / "E1_pr_vs_licks.png")
    readme(d, "E1 - PR as a function of licks (idea I03)", [
        "- Schedule start 1, increment 1: trial k needs k counted licks, so n rewards need n(n+1)/2 counted licks.",
        "- `rule_holds`: rewards computed from the counted licks = the real rewards.",
        "- `censored_at_910s`: no final pause of 30 s or more before 910 s, so the mouse was still licking and the PR is "
        "a lower bound.",
        "- `quit_min`: time of the last lick when the session ends with a pause of 30 s or more.",
        "- Counted licks and rewards come from licks_long.csv / rewards_long.csv (Saleae, first 910 s); master_* columns "
        "are master_day_level.csv for comparison.",
        "- References: Richardson & Roberts 1996; PR breakpoint definition and scaling, J Neurosci Methods 2021."])
    return E


# ------------------------------------------------------------------------------------------------ E2
def e2(M, E1, man, out: Path, res: dict):
    d = out / "E2_time_of_day"
    d.mkdir(parents=True, exist_ok=True)
    C = pd.concat([pd.read_csv(A / "days" / f"day{k}" / f"pr_checks_day{k}.csv",
                               usecols=["day", "mouse", "session", "log_start"]) for k in PR_DAYS])
    C["start_h"] = [(lambda x: x.hour + x.minute / 60 + x.second / 3600)(dt.datetime.fromisoformat(s))
                    for s in C.log_start]
    X = E1[["day", "mouse", "session", "group", "sex", "rewards_910"]].merge(C, on=["day", "mouse", "session"],
                                                                             how="left")
    X = X.merge(M[["day", "mouse", "pupil_median_px", "pupil_grade"]], on=["day", "mouse"], how="left")
    X["cage"] = X.mouse.str.split("_").str[0]
    for c in ("start_h", "rewards_910"):
        X[c + "_dev"] = X[c] - X.groupby("mouse")[c].transform("mean")
    pup = X.pupil_grade != "whole_eye"
    X["pupil_dev_pct"] = np.where(pup, 100 * (X.pupil_median_px / X.groupby("mouse").pupil_median_px.transform("mean")
                                              - 1), np.nan)
    X.to_csv(d / "E2_time_of_day.csv", index=False)
    rr = spearman(X.start_h_dev, X.rewards_910_dev)
    rp = spearman(X.start_h_dev, X.pupil_dev_pct)
    cage = X.groupby("cage").agg(start_h=("start_h", "mean"), rewards=("rewards_910", "mean"))
    res["E2"] = dict(n=len(X), within_mouse_rho_time_rewards=rr[0], p_rewards=rr[1],
                     within_mouse_rho_time_pupil=rp[0], p_pupil=rp[1], n_pupil=rp[2],
                     cage_mean_start_h={k: round(float(v), 2) for k, v in cage.start_h.items()},
                     cage_mean_rewards={k: round(float(v), 1) for k, v in cage.rewards.items()},
                     within_mouse_start_range_h=float(X.groupby("mouse").start_h.agg(lambda s: s.max() - s.min())
                                                      .median()))
    col = mouse_colors(X.mouse.unique())
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for m, g in X.groupby("mouse"):
        g = g.sort_values("start_h")
        ax[0].plot(g.start_h, g.rewards_910, "-o", color=col[m], label=lab(m), ms=5)
    ax[0].set_xlabel("session start (clock hour, days 3-5)")
    ax[0].set_ylabel("rewards in 910 s")
    ax[0].set_title("E2a  9670f (female) always runs first, 9669m (male) after", fontsize=10)
    ax[0].legend(fontsize=7, ncol=2)
    ax[1].scatter(X.start_h_dev * 60, X.rewards_910_dev, c=[col[m] for m in X.mouse], s=35)
    ax[1].axhline(0, color="0.6", lw=0.8)
    ax[1].axvline(0, color="0.6", lw=0.8)
    ax[1].set_xlabel("start time minus the mouse's mean (min)")
    ax[1].set_ylabel("rewards minus the mouse's mean")
    ax[1].set_title(f"E2b  within mouse: Spearman {rr[0]:+.2f} (p {rr[1]:.2f}, n {rr[2]})", fontsize=10)
    fig.tight_layout()
    save(fig, d / "E2_time_of_day.png")
    readme(d, "E2 - time of day (idea I02)", [
        "- Start = log_start of pr_checks_day<N>.csv (local clock).",
        "- Between cages, sex and time block are the same thing (9670f always first), so they cannot be separated.",
        "- Within a mouse the start time varies a little from day to day; E2b shows whether those small shifts go with "
        "the rewards. This is weak by design: a few tens of minutes, 3 days.",
        "- Pupil: whole-eye mice left out (a different measure)."])
    return X


# ------------------------------------------------------------------------------------------------ E3
def e3(M, L, man, out: Path, res: dict):
    d = out / "E3_lick_clusters"
    d.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in M[M.day.isin(ALL_DAYS)].itertuples():
        lk = L[(L.session == s.session) & (L.mouse == s.mouse)]
        t = np.sort(lk[(lk.t_s >= 0) & (lk.t_s <= WIN_S)].t_s.to_numpy(float))
        cl, bo = runs(t, CLUSTER_ILI_S), runs(t, BOUT_GAP_S)
        rows.append(dict(day=int(s.day), mouse=s.mouse, group=man.loc[s.mouse, "group"], task=s.task,
                         session=s.session, n_licks=len(t), rewards=s.rewards,
                         n_clusters=len(cl), mean_cluster=float(cl.mean()) if len(cl) else np.nan,
                         median_cluster=float(np.median(cl)) if len(cl) else np.nan,
                         n_clusters_ge3=int((cl >= 3).sum()),
                         mean_cluster_ge3=float(cl[cl >= 3].mean()) if (cl >= 3).any() else np.nan,
                         n_bouts_1s=len(bo), mean_bout_1s=float(bo.mean()) if len(bo) else np.nan,
                         master_bout_count=s.bout_count, master_mean_licks_per_bout=s.mean_licks_per_bout))
    E = pd.DataFrame(rows)
    E.to_csv(d / "E3_lick_clusters.csv", index=False)
    pr = E[E.day.isin(PR_DAYS)]
    rn = spearman(pr.n_clusters, pr.rewards)
    rs = spearman(pr.mean_cluster, pr.rewards)
    dec = []
    for m, g in pr.groupby("mouse"):
        a, b = g[g.day == 3], g[g.day == 5]
        if len(a) and len(b):
            a, b = a.iloc[0], b.iloc[0]
            dec.append(dict(mouse=m, dlog_licks=math.log(b.n_licks / a.n_licks),
                            dlog_clusters=math.log(b.n_clusters / a.n_clusters),
                            dlog_cluster_size=math.log(b.mean_cluster / a.mean_cluster)))
    D = pd.DataFrame(dec)
    D.to_csv(d / "E3_day3_to_day5_decomposition.csv", index=False)
    res["E3"] = dict(n_sessions=len(E), rho_clusters_rewards=rn[0], p_clusters=rn[1],
                     rho_cluster_size_rewards=rs[0], p_size=rs[1],
                     median_cluster_size_by_day={int(k): round(float(v), 2) for k, v in
                                                 E.groupby("day").mean_cluster.median().items()},
                     median_n_clusters_by_day={int(k): int(v) for k, v in E.groupby("day").n_clusters.median().items()},
                     share_of_change_from_number=float(np.median(D.dlog_clusters / D.dlog_licks)) if len(D) else np.nan)
    col = mouse_colors(E.mouse.unique())
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for m, g in E.groupby("mouse"):
        g = g.sort_values("day")
        gc = GCOL.get(man.loc[m, "group"], "k")
        ax[0].plot(g.day, g.n_clusters, "-o", color=gc, alpha=0.75, ms=4)
        ax[1].plot(g.day, g.mean_cluster, "-o", color=gc, alpha=0.75, ms=4)
    for gname, gc in GCOL.items():
        ax[0].plot([], [], color=gc, label=f"{gname} from day 6")
    ax[0].legend(fontsize=8)
    ax[0].set_title("E3a  number of lick clusters (ILI < 0.5 s) = 'seeking'", fontsize=10)
    ax[1].set_title("E3b  licks per cluster = 'liking' (Dwyer 2012)", fontsize=10)
    for a in ax[:2]:
        a.set_xticks(ALL_DAYS)
        a.set_xlabel("day (1-2 FR, 3-5 PR)")
    ax[0].set_ylabel("clusters in 910 s")
    ax[1].set_ylabel("mean licks per cluster")
    ax[2].scatter(pr.n_clusters, pr.rewards, c=[col[m] for m in pr.mouse], s=30, label="clusters")
    ax2 = ax[2].twiny()
    ax2.scatter(pr.mean_cluster, pr.rewards, c=[col[m] for m in pr.mouse], s=30, marker="x")
    ax[2].set_xlabel("clusters (dots)")
    ax2.set_xlabel("licks per cluster (crosses)")
    ax[2].set_ylabel("rewards (days 3-5)")
    ax[2].set_title(f"E3c  rewards vs number rho {rn[0]:+.2f} / size rho {rs[0]:+.2f}", fontsize=10, y=1.12)
    fig.tight_layout()
    save(fig, d / "E3_lick_clusters.png")
    readme(d, "E3 - lick clusters: number vs size (idea L04)", [
        "- Clusters: runs of Saleae licks with inter-lick intervals < 0.5 s (Dwyer 2012); bouts: gap >= 1 s.",
        "- Cluster size is read as a 'liking' (hedonic) index, cluster number as initiation / 'seeking' (Higgs & "
        "Cooper 1998; the mouse nociceptin lick study of 2016).",
        "- Days 1-5 are all before morphine, and the groups are the ones assigned for days 6-10.",
        "- `E3_day3_to_day5_decomposition.csv`: log change of licks = log change of the number of clusters + log change "
        "of the cluster size."])
    return E


# ------------------------------------------------------------------------------------------------ E4
def e4(man, out: Path, res: dict):
    import scipy.io as sio
    d = out / "E4_hp_latency"
    d.mkdir(parents=True, exist_ok=True)
    rows = []
    for f in sorted((P / "assays" / "day5").glob("*/*/ManualScoringResults_*.mat")):
        if "hot" not in f.name.lower():
            continue
        cage, mouse = f.parent.parent.name.lower(), f.parent.name.lower()
        S = sio.loadmat(f, squeeze_me=True)
        sc = np.asarray(S["score"]).ravel().astype(int)
        names = [str(x) for x in np.atleast_1d(S["behaviorNames"])]
        fr = float(np.asarray(S.get("frameRate", 30.0)).ravel()[0])
        low = [n.lower().strip() for n in names]
        strict = [i for i, n in enumerate(low) if ("flinch" in n or "hind" in n or "jump" in n)]
        broad = strict + [i for i, n in enumerate(low) if n == "licking"]

        def first(ids):
            j = np.flatnonzero(np.isin(sc, ids))
            return float(j[0] / fr) if len(j) else np.nan

        # user 2026-11-19: flinching is the REFLEXIVE response (like the tail flick); licking, hindlimb licking and
        # rearing are AFFECTIVE (not reflexive). Jump (escape) is kept apart.
        reflex = [i for i, n in enumerate(low) if "flinch" in n]
        affect = [i for i, n in enumerate(low) if n in ("licking", "hindlimb licking", "rearing")]
        jump = [i for i, n in enumerate(low) if "jump" in n]

        def onsets(ids):
            on = np.isin(sc, ids).astype(int)
            return int(np.sum(np.diff(np.r_[0, on]) == 1))

        r = dict(cage=cage, mouse=mouse, mouse_id=f"{cage}_{mouse}", group=man.loc[f"{cage}_{mouse}", "group"],
                 file=f.name, frames=len(sc), fps=fr, seconds=round(len(sc) / fr, 2), behaviour_names="|".join(names),
                 latency_reflexive_flinch_s=first(reflex), latency_affective_s=first(affect),
                 latency_jump_s=first(jump), n_flinch_onsets=onsets(reflex), n_affective_onsets=onsets(affect),
                 pct_affective=round(100 * float(np.mean(np.isin(sc, affect))), 2),
                 latency_flinch_hindlick_jump_s=first(strict), latency_incl_licking_s=first(broad))
        for i, n in enumerate(names):
            r[f"pct_{n.replace(' ', '_')}"] = round(100 * float(np.mean(sc == i)), 2)
        rows.append(r)
    H = pd.DataFrame(rows)
    TI = pd.read_csv(P / "assays" / "tail_immersion_all.csv")
    TI = TI[TI.day == 5].copy()
    TI["mouse_id"] = TI.cage.str.lower() + "_" + TI.mouse.str.lower()
    H = H.merge(TI[["mouse_id", "latency_s"]].rename(columns={"latency_s": "tail_immersion_s"}), on="mouse_id",
                how="left")
    H.to_csv(d / "E4_hp_latency_day5.csv", index=False)
    rs = spearman(H.latency_flinch_hindlick_jump_s, H.tail_immersion_s)
    rb = spearman(H.latency_incl_licking_s, H.tail_immersion_s)
    rr = spearman(H.latency_reflexive_flinch_s, H.tail_immersion_s)
    ra = spearman(H.latency_affective_s, H.tail_immersion_s)
    rpa = spearman(H.pct_affective, H.tail_immersion_s)
    rra = spearman(H.latency_reflexive_flinch_s, H.latency_affective_s)
    res["E4"] = dict(n=len(H), names=H.behaviour_names.iloc[0] if len(H) else "",
                     rho_reflexive_vs_TI=rr[0], p_reflexive=rr[1], n_reflexive=rr[2],
                     rho_affective_vs_TI=ra[0], p_affective=ra[1], n_affective=ra[2],
                     rho_pct_affective_vs_TI=rpa[0], rho_reflexive_vs_affective=rra[0],
                     no_flinch_in_30s=int(H.latency_reflexive_flinch_s.isna().sum()),
                     latency_reflexive_s={lab(m): (None if pd.isna(v) else round(float(v), 2))
                                          for m, v in zip(H.mouse_id, H.latency_reflexive_flinch_s)},
                     latency_affective_s={lab(m): (None if pd.isna(v) else round(float(v), 2))
                                          for m, v in zip(H.mouse_id, H.latency_affective_s)},
                     tail_immersion_s={lab(m): (None if pd.isna(v) else round(float(v), 2))
                                       for m, v in zip(H.mouse_id, H.tail_immersion_s)},
                     latency_strict_s={lab(m): (None if pd.isna(v) else round(float(v), 2))
                                       for m, v in zip(H.mouse_id, H.latency_flinch_hindlick_jump_s)},
                     latency_broad_s={lab(m): (None if pd.isna(v) else round(float(v), 2))
                                      for m, v in zip(H.mouse_id, H.latency_incl_licking_s)},
                     no_response_in_30s=int(H.latency_flinch_hindlick_jump_s.isna().sum()),
                     rho_strict_vs_TI=rs[0], p_strict=rs[1], rho_broad_vs_TI=rb[0])
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.5))
    for ax, col, rho, ttl in ((axs[0], "latency_reflexive_flinch_s", rr, "E4a  REFLEXIVE: first flinch"),
                              (axs[1], "latency_affective_s", ra,
                               "E4b  AFFECTIVE: first licking / hindlimb licking / rearing")):
        for r in H.itertuples():
            v = getattr(r, col)
            y = v if pd.notna(v) else 30.5
            ax.scatter(r.tail_immersion_s, y, color=GCOL.get(r.group, "k"), s=45, marker="o" if pd.notna(v) else "^")
            ax.annotate(lab(r.mouse_id), (r.tail_immersion_s, y), fontsize=7, xytext=(3, 3),
                        textcoords="offset points")
        ax.set_xlabel("tail immersion latency, day 5 (s)")
        ax.set_ylabel("hot plate latency from placement (s)")
        ax.set_title(f"{ttl}\nSpearman vs TI {rho[0]:+.2f} (p {rho[1]:.2f}, n {rho[2]}); triangle = none in 30 s",
                     fontsize=9.5)
    fig.tight_layout()
    save(fig, d / "E4_hp_latency_vs_ti.png")
    readme(d, "E4 - hot-plate first-response latency (idea L05)", [
        "- From the hand scores (lab scorer, 30 s): the first frame scored as Flinching / Hindlimb licking / Jump "
        "(`latency_flinch_hindlick_jump_s`), and also with forepaw Licking counted (`latency_incl_licking_s`).",
        "- Day 5 = the baseline before morphine. The plan uses % time per behaviour (HOT_Pct); this adds the standard "
        "latency (Balter & Dykstra 2013: hot-plate latency tracks spontaneous morphine withdrawal in mice).",
        "- No response within 30 s = NaN (would be censored at 30 s).",
        "- **User 2026-11-19:** flinching is the REFLEXIVE response (closer to the tail-immersion reflex); licking, "
        "hindlimb licking and rearing are AFFECTIVE responses, not reflexes. So the latencies are now split: "
        "`latency_reflexive_flinch_s` (compare with TI) and `latency_affective_s`; jump is kept apart "
        "(`latency_jump_s`). The old mixed columns are kept for comparison.",
        "- The video starts when the mouse is put on the plate (user 2026-11-19), so frame 0 = placement."])
    return H


# ------------------------------------------------------------------------------------------------ E6
def hmm2_poisson(c, n_iter=300, tol=1e-6):
    """2-state Poisson HMM (Baum-Welch, scaled). State 1 = engaged (higher rate)."""
    from scipy.special import gammaln
    c = np.asarray(c, float)
    T = len(c)
    pos = c[c > 0]
    lam = np.array([max(0.1 * c.mean(), 1e-3), max(pos.mean() if len(pos) else 1.0, 0.5)])
    Atr = np.array([[0.95, 0.05], [0.05, 0.95]])
    pi = np.array([0.5, 0.5])
    lg = gammaln(c + 1)
    ll_old, it = -np.inf, 0
    for it in range(n_iter):
        logB = c[:, None] * np.log(lam[None, :]) - lam[None, :] - lg[:, None]
        sh = logB.max(1)
        B = np.exp(logB - sh[:, None])
        al = np.zeros((T, 2))
        s = np.zeros(T)
        a = pi * B[0]
        s[0] = a.sum()
        al[0] = a / s[0]
        for t in range(1, T):
            a = (al[t - 1] @ Atr) * B[t]
            s[t] = a.sum()
            al[t] = a / s[t]
        ll = float(np.log(s).sum() + sh.sum())
        be = np.ones((T, 2))
        for t in range(T - 2, -1, -1):
            be[t] = (Atr @ (B[t + 1] * be[t + 1])) / s[t + 1]
        g = al * be
        g /= g.sum(1, keepdims=True)
        xi = (al[:-1, :, None] * Atr[None] * (B[1:] * be[1:])[:, None, :]) / s[1:, None, None]
        Atr = xi.sum(0)
        Atr /= Atr.sum(1, keepdims=True)
        pi = g[0]
        lam = np.maximum((g * c[:, None]).sum(0) / g.sum(0), 1e-4)
        if ll - ll_old < tol:
            break
        ll_old = ll
    o = np.argsort(lam)
    return dict(lam=lam[o], A=Atr[np.ix_(o, o)], gamma=g[:, o], ll=ll, iters=it + 1)


def e6(M, L, E1, man, out: Path, res: dict):
    d = out / "E6_engagement_hmm"
    d.mkdir(parents=True, exist_ok=True)
    rows, keep = [], {}
    for s in M[M.day.isin(PR_DAYS)].itertuples():
        lk = L[(L.session == s.session) & (L.mouse == s.mouse)]
        t = lk[(lk.t_s >= 0) & (lk.t_s < WIN_S)].t_s.to_numpy(float)
        cnt = np.bincount(np.floor(t).astype(int), minlength=int(WIN_S))[: int(WIN_S)]
        h = hmm2_poisson(cnt)
        st = (h["gamma"][:, 1] > 0.5).astype(int)
        eng = runs(np.flatnonzero(st == 1).astype(float), 1.5)          # engaged bouts (consecutive seconds)
        ts = np.sort(t)
        gaps = np.diff(np.r_[0.0, ts, WIN_S])
        rows.append(dict(day=int(s.day), mouse=s.mouse, group=man.loc[s.mouse, "group"], session=s.session,
                         rewards=s.rewards, licks=len(t), engaged_frac=float(st.mean()),
                         engaged_min=float(st.sum() / 60), n_engaged_bouts=len(eng),
                         mean_engaged_bout_s=float(eng.mean()) if len(eng) else 0.0,
                         rate_engaged_hz=float(h["lam"][1]), rate_disengaged_hz=float(h["lam"][0]),
                         dwell_engaged_s=float(1 / max(1 - h["A"][1, 1], 1e-9)),
                         dwell_disengaged_s=float(1 / max(1 - h["A"][0, 0], 1e-9)),
                         pause30_frac=float(gaps[gaps >= CRIT_PAUSE_S].sum() / WIN_S), iters=h["iters"]))
        keep[(int(s.day), s.mouse)] = (cnt, h["gamma"][:, 1])
    E = pd.DataFrame(rows)
    E.to_csv(d / "E6_engagement_hmm.csv", index=False)
    re_ = spearman(E.engaged_min, E.rewards)
    rp = spearman(1 - E.pause30_frac, E.rewards)
    res["E6"] = dict(n=len(E), rho_engaged_min_rewards=re_[0], rho_not_paused30_rewards=rp[0],
                     median_engaged_frac=float(E.engaged_frac.median()),
                     median_rate_engaged_hz=float(E.rate_engaged_hz.median()),
                     median_dwell_engaged_s=float(E.dwell_engaged_s.median()),
                     median_dwell_disengaged_s=float(E.dwell_disengaged_s.median()))
    ex = E.sort_values("rewards").iloc[len(E) // 2]
    cnt, pe = keep[(int(ex.day), ex.mouse)]
    fig = plt.figure(figsize=(13, 6.4))
    ax0 = fig.add_axes([0.06, 0.58, 0.9, 0.34])
    ax0.bar(np.arange(len(cnt)) / 60, cnt, width=1 / 60, color="0.35")
    ax0.set_ylabel("licks per s")
    axp = ax0.twinx()
    axp.plot(np.arange(len(pe)) / 60, pe, color="#d35400", lw=1)
    axp.set_ylabel("P(engaged)", color="#d35400")
    axp.set_ylim(-0.05, 1.05)
    ax0.set_title(f"E6a  example: {lab(ex.mouse)} day {int(ex.day)} - engaged {100 * ex.engaged_frac:.0f} % of the "
                  f"session, {int(ex.n_engaged_bouts)} engaged bouts", fontsize=10)
    ax0.set_xlabel("min")
    ax1 = fig.add_axes([0.06, 0.08, 0.4, 0.36])
    ax2 = fig.add_axes([0.56, 0.08, 0.4, 0.36])
    col = mouse_colors(E.mouse.unique())
    ax1.scatter(E.engaged_min, E.rewards, c=[col[m] for m in E.mouse], s=30)
    ax1.set_xlabel("engaged time (min, HMM)")
    ax1.set_ylabel("rewards")
    ax1.set_title(f"E6b  rewards vs engaged time: rho {re_[0]:+.2f}", fontsize=10)
    for m, g in E.groupby("mouse"):
        g = g.sort_values("day")
        ax2.plot(g.day, g.dwell_engaged_s, "-o", color=GCOL.get(man.loc[m, "group"], "k"), alpha=0.75, ms=4)
    ax2.set_xticks(PR_DAYS)
    ax2.set_xlabel("day")
    ax2.set_ylabel("mean engaged dwell (s)")
    ax2.set_title("E6c  how long an engaged state lasts (red = Active from day 6)", fontsize=10)
    save(fig, d / "E6_engagement_hmm.png")
    readme(d, "E6 - engaged / disengaged states (idea L06)", [
        "- 2-state Poisson hidden Markov model on licks per 1-s bin (first 910 s), fitted per session (Baum-Welch).",
        "- Engaged = the higher-rate state; dwell = 1 / (1 - stay probability) in s.",
        "- `pause30_frac` = the time in pauses of 30 s or more (the fixed threshold of plan N04), for comparison.",
        "- Reference: Ashwood et al. 2022 (GLM-HMM; states instead of fixed lapse / pause rules)."])
    return E


# ------------------------------------------------------------------------------------------------ E5
def icc(Y):
    from scipy.stats import f as Fd
    n, k = Y.shape
    gm = Y.mean()
    ssr = k * ((Y.mean(1) - gm) ** 2).sum()
    ssc = n * ((Y.mean(0) - gm) ** 2).sum()
    sse = ((Y - gm) ** 2).sum() - ssr - ssc
    msr, msc, mse = ssr / (n - 1), ssc / (k - 1), sse / ((n - 1) * (k - 1))
    i31 = (msr - mse) / (msr + (k - 1) * mse)
    i21 = (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
    Fv = msr / mse if mse > 0 else np.inf
    d1, d2 = n - 1, (n - 1) * (k - 1)
    FL, FU = Fv / Fd.ppf(0.975, d1, d2), Fv * Fd.ppf(0.975, d2, d1)
    return float(i31), float((FL - 1) / (FL + k - 1)), float((FU - 1) / (FU + k - 1)), float(i21)


def e5(M, E3, E6, out: Path, res: dict):
    d = out / "E5_reliability"
    d.mkdir(parents=True, exist_ok=True)
    X = M[M.day.isin(PR_DAYS)][["day", "mouse", "rewards", "licks", "licks_per_min", "mean_licks_per_bout",
                                "lick_rate_in_bout_hz", "within_bout_ili_ms", "time_to_reward_median_s",
                                "pupil_median_px", "pupil_reward_response_pct", "pupil_grade"]].copy()
    X = X.merge(E3[["day", "mouse", "n_clusters", "mean_cluster"]], on=["day", "mouse"], how="left")
    if E6 is not None:
        X = X.merge(E6[["day", "mouse", "engaged_min", "rate_engaged_hz"]], on=["day", "mouse"], how="left")
    V = pd.read_csv(A / "pupil_check" / "pupil_check_verdicts.csv")
    V["mouse_id"] = V.cage.str.lower() + "_" + V.mouse.str.lower()
    noisy = {(int(r.day), r.mouse_id) for r in V.itertuples() if str(r.verdict).strip() == "noisy"}
    whole = set(X[X.pupil_grade == "whole_eye"].mouse)
    rows = []
    meas = ["rewards", "licks", "licks_per_min", "n_clusters", "mean_cluster", "mean_licks_per_bout",
            "lick_rate_in_bout_hz", "within_bout_ili_ms", "time_to_reward_median_s", "engaged_min", "rate_engaged_hz",
            "pupil_median_px", "pupil_reward_response_pct"]
    for c in meas:
        if c not in X:
            continue
        for version in (("all", None), ("pupil: ok sessions only", "ok")) if c.startswith("pupil") else (("all", None),):
            Z = X[~X.mouse.isin(whole)] if c.startswith("pupil") else X
            if version[1] == "ok":
                Z = Z[[(int(a), b) not in noisy for a, b in zip(Z.day, Z.mouse)]]
            W = Z.pivot_table(index="mouse", columns="day", values=c).dropna()
            if len(W) < 3 or W.shape[1] < 3:
                rows.append(dict(measure=c, version=version[0], n_mice=len(W)))
                continue
            i31, lo, hi, i21 = icc(W.to_numpy(float))
            rho = [spearman(W[a], W[b])[0] for a, b in ((3, 4), (4, 5), (3, 5))]
            rows.append(dict(measure=c, version=version[0], n_mice=len(W), icc31=i31, icc31_lo=lo, icc31_hi=hi,
                             icc21=i21, rho_d3d4=rho[0], rho_d4d5=rho[1], rho_d3d5=rho[2]))
    I = pd.DataFrame(rows)
    I.to_csv(d / "E5_reliability_days3-5.csv", index=False)
    res["E5"] = {f"{r.measure} [{r.version}]": dict(n=int(r.n_mice), icc31=None if pd.isna(r.get("icc31", np.nan))
                                                    else round(float(r.icc31), 2),
                                                    ci=None if pd.isna(r.get("icc31", np.nan)) else
                                                    [round(float(r.icc31_lo), 2), round(float(r.icc31_hi), 2)])
                 for _, r in I.iterrows()}
    J = I.dropna(subset=["icc31"]) if "icc31" in I else I.iloc[0:0]
    fig, ax = plt.subplots(figsize=(8.5, 0.42 * len(J) + 1.6))
    y = np.arange(len(J))
    ax.errorbar(J.icc31, y, xerr=[J.icc31 - J.icc31_lo, J.icc31_hi - J.icc31], fmt="o", color="k", ecolor="0.5",
                capsize=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{m} ({v}, n {n})" for m, v, n in zip(J.measure, J.version, J.n_mice)], fontsize=8)
    ax.axvline(0, color="0.6", lw=0.8)
    ax.axvline(0.75, color="#27ae60", lw=0.8, ls=":")
    ax.axvline(0.5, color="#f39c12", lw=0.8, ls=":")
    ax.set_xlim(-1, 1)
    ax.set_xlabel("ICC(3,1) over days 3-5 with 95 % CI (dotted: 0.5 moderate, 0.75 good)")
    ax.set_title("E5  which measures are stable mouse traits before morphine?", fontsize=10)
    ax.invert_yaxis()
    fig.tight_layout()
    save(fig, d / "E5_reliability.png")
    readme(d, "E5 - day-to-day reliability, days 3-5 (idea I04)", [
        "- ICC(3,1): two-way, consistency, single day; ICC(2,1): absolute agreement. 95 % CI from the F distribution "
        "(Shrout & Fleiss).",
        "- Mice with a missing day are left out of that measure. Pupil: whole-eye mice left out; the 'ok sessions only' "
        "row also leaves out sessions graded noisy (pupil_check_verdicts.csv).",
        "- n = 8 mice (6 for pupil), so the CIs are wide."])
    return I


# ------------------------------------------------------------------------------------------------ E7
def lagged(e: np.ndarray, j0: int, j1: int, wb: int) -> np.ndarray:
    """columns j0..j1-1: events whose lag (t - u) lies in [j*wb, (j+1)*wb) frames (negative j = before the event)"""
    n = len(e)
    Cs = np.r_[0.0, np.cumsum(e)]
    t = np.arange(n)
    cols = []
    for j in range(j0, j1):
        a = np.clip(t - (j + 1) * wb + 1, 0, n)
        b = np.clip(t - j * wb + 1, 0, n)
        cols.append(Cs[b] - Cs[a])
    return np.stack(cols, 1)


def e7(M, L, R, man, out: Path, res: dict):
    import concat_choice as CC
    import pupil_check_examples as PCE
    d = out / "E7_pupil_glm"
    d.mkdir(parents=True, exist_ok=True)
    root = CC.active_pupil_root(P)
    G = PCE.grades(root)
    V = PCE.verdicts()
    wb = 8                                            # 8 frames = 0.264 s
    rj = (-8, 23)                                     # reward lags -2.1 .. +6.1 s
    lj = (-4, 12)                                     # lick lags -1.1 .. +3.2 s
    lag_r = (np.arange(*rj) + 0.5) * wb / CAM_HZ
    lag_l = (np.arange(*lj) + 0.5) * wb / CAM_HZ
    rows, ker = [], {}
    for s in M[M.day.isin(PR_DAYS)].itertuples():
        cage, mouse = s.mouse.split("_", 1)
        try:
            ss = PCE.session(int(s.day), cage, mouse, root, G)
        except Exception as ex:  # noqa: BLE001
            rows.append(dict(day=s.day, mouse=s.mouse, status=f"no trace ({type(ex).__name__})"))
            continue
        if ss.get("kind") not in ("pupil", "eye"):
            rows.append(dict(day=s.day, mouse=s.mouse, status=ss.get("note", "no trace")))
            continue
        dd = np.asarray(ss["d"], float)
        n = min(len(dd), int(WIN_S * CAM_HZ))
        dd, bad = dd[:n], np.asarray(ss["out"], bool)[:n]
        med = float(np.nanmedian(dd[~bad]))
        y = 100 * (dd / med - 1)
        ok = np.isfinite(y) & ~bad
        lk = L[(L.session == s.session) & (L.mouse == s.mouse) & (L.frame_reliable == True)].frame.dropna()  # noqa: E712
        rw = R[(R.session == s.session) & (R.mouse == s.mouse)].frame.dropna()
        el, er = np.zeros(n), np.zeros(n)
        fl, fr_ = lk.astype(int).to_numpy(), rw.astype(int).to_numpy()
        np.add.at(el, fl[(fl >= 0) & (fl < n)], 1)
        np.add.at(er, fr_[(fr_ >= 0) & (fr_ < n)], 1)
        Xr, Xl = lagged(er, *rj, wb), lagged(el, *lj, wb)
        K = int(2 * n / (60 * CAM_HZ))
        tt = np.arange(n)
        Xd = np.stack([np.ones(n)] + [np.cos(np.pi * k * (tt + 0.5) / n) for k in range(1, K + 1)], 1)

        def fit(Xev):
            X = np.hstack([Xev, Xd])[ok]
            yy = y[ok]
            pe = Xev.shape[1]
            lam = (1e-2 * float(np.mean(np.sum(Xev[ok] ** 2, 0))) + 1e-9) if pe else 0.0
            Xa = np.vstack([X, np.hstack([np.sqrt(lam) * np.eye(pe), np.zeros((pe, Xd.shape[1]))])])
            ya = np.r_[yy, np.zeros(pe)]
            b = np.linalg.lstsq(Xa, ya, rcond=None)[0]
            r2 = 1 - np.sum((yy - X @ b) ** 2) / np.sum((yy - yy.mean()) ** 2)
            return b[:pe], float(r2)

        bA, r2A = fit(Xr)
        bB, r2B = fit(np.hstack([Xr, Xl]))
        b0, r20 = fit(np.zeros((n, 0)))
        bBr, bBl = bB[: Xr.shape[1]], bB[Xr.shape[1]:]
        win = (lag_r >= 0.5) & (lag_r <= 2.5)
        pre = (lag_r >= -2.0) & (lag_r < 0)
        # naive reward-locked average, -1..0 s baseline (as quick_analysis)
        tax = TAX
        H = []
        for f0 in fr_:
            idx = f0 + tax
            if idx[0] < 0 or idx[-1] >= n:
                continue
            seg = np.where(ok[idx], y[idx], np.nan)
            base = np.nanmean(seg[(tax >= -CAM_HZ) & (tax < 0)])
            if np.isfinite(base):
                H.append(seg - base)
        naive = np.nanmean(np.array(H), 0) if H else np.full(len(tax), np.nan)
        verdict = V.get((int(s.day), cage, mouse), ("-", ""))[0]
        rows.append(dict(day=int(s.day), mouse=s.mouse, group=man.loc[s.mouse, "group"], kind=ss["kind"],
                         verdict=verdict, status="ok", frames=n, frames_used=int(ok.sum()), rewards=len(fr_),
                         licks=len(fl), r2_drift=r20, r2_reward=r2A, r2_reward_lick=r2B,
                         reward_peak_no_lick_model=float(bA[win].mean()),
                         reward_peak_with_lick_model=float(bBr[win].mean()),
                         reward_pre_no_lick_model=float(bA[pre].mean()),
                         reward_pre_with_lick_model=float(bBr[pre].mean()),
                         lick_kernel_peak=float(bBl[(lag_l >= 0.3) & (lag_l <= 1.5)].mean()),
                         naive_peak=float(np.nanmean(naive[(tax / CAM_HZ >= 0.5) & (tax / CAM_HZ <= 2.5)]))))
        ker[(int(s.day), s.mouse)] = dict(A=bA, Br=bBr, Bl=bBl, naive=naive, kind=ss["kind"], verdict=verdict)
    E = pd.DataFrame(rows)
    E.to_csv(d / "E7_pupil_glm.csv", index=False)
    np.savez(d / "E7_kernels.npz", lag_reward_s=lag_r, lag_lick_s=lag_l, naive_t_s=TAX / CAM_HZ, **{f"d{k[0]}_{k[1]}_{nm}": v[nm] for k, v in ker.items() for nm in ("A", "Br", "Bl", "naive")})

    def mouse_mean(sel, nm):
        by = {}
        for (day, m), v in ker.items():
            if sel(v):
                by.setdefault(m, []).append(v[nm])
        return {m: np.nanmean(np.array(x), 0) for m, x in by.items()}

    out_sets = {}
    for tag, sel in (("pupil mice, all sessions", lambda v: v["kind"] == "pupil"),
                     ("pupil mice, ok sessions", lambda v: v["kind"] == "pupil" and v["verdict"] == "ok")):
        mA, mBr, mBl, mN = (mouse_mean(sel, nm) for nm in ("A", "Br", "Bl", "naive"))
        out_sets[tag] = (mA, mBr, mBl, mN)
    Ep = E[(E.get("kind") == "pupil") & (E.status == "ok")] if "kind" in E else E.iloc[0:0]
    per_mouse = Ep.groupby("mouse")[["reward_peak_no_lick_model", "reward_peak_with_lick_model",
                                     "reward_pre_no_lick_model", "reward_pre_with_lick_model", "naive_peak",
                                     "r2_reward", "r2_reward_lick"]].mean()
    per_mouse.to_csv(d / "E7_per_mouse_days3-5.csv")
    keep_frac = (per_mouse.reward_peak_with_lick_model / per_mouse.reward_peak_no_lick_model).replace(
        [np.inf, -np.inf], np.nan)
    res["E7"] = dict(sessions=int((E.status == "ok").sum()) if "status" in E else 0,
                     pupil_mice=int(per_mouse.shape[0]),
                     median_reward_peak_no_lick=float(per_mouse.reward_peak_no_lick_model.median()),
                     median_reward_peak_with_lick=float(per_mouse.reward_peak_with_lick_model.median()),
                     median_fraction_left=float(keep_frac.median()),
                     median_naive_peak=float(per_mouse.naive_peak.median()),
                     median_r2_reward=float(per_mouse.r2_reward.median()),
                     median_r2_reward_lick=float(per_mouse.r2_reward_lick.median()),
                     not_used=E[E.status != "ok"][["day", "mouse", "status"]].to_dict("records") if "status" in E
                     else [])
    fig, ax = plt.subplots(2, 2, figsize=(13, 8.4))
    mA, mBr, mBl, mN = out_sets["pupil mice, all sessions"]
    t_n = TAX / CAM_HZ

    def band(a, x, Ms, color, label):
        if not Ms:
            return
        Z = np.array(list(Ms.values()))
        mu, se = np.nanmean(Z, 0), np.nanstd(Z, 0) / np.sqrt(max(len(Z), 1))
        a.plot(x, mu, color=color, label=f"{label} (n {len(Z)} mice)")
        a.fill_between(x, mu - se, mu + se, color=color, alpha=0.2, lw=0)

    for m, v in mN.items():
        ax[0, 0].plot(t_n, v, lw=0.8, alpha=0.6, label=lab(m))
    ax[0, 0].axvline(0, color="k", lw=0.8)
    ax[0, 0].set_title("E7a  plain reward-locked pupil (each mouse, days 3-5)", fontsize=10)
    ax[0, 0].set_xlabel("s from reward")
    ax[0, 0].set_ylabel("pupil change from -1..0 s (% of median)")
    ax[0, 0].legend(fontsize=7)
    band(ax[0, 1], lag_r, mA, "#7f8c8d", "reward kernel, reward-only model")
    band(ax[0, 1], lag_r, mBr, "#c0392b", "reward kernel, licks also in the model")
    mA2, mBr2, _, _ = out_sets["pupil mice, ok sessions"]
    band(ax[0, 1], lag_r, mBr2, "#e67e22", "same, ok sessions only")
    ax[0, 1].axvline(0, color="k", lw=0.8)
    ax[0, 1].axhline(0, color="0.6", lw=0.6)
    ax[0, 1].set_title("E7b  GLM reward kernel: how much is left once licking is modelled", fontsize=10)
    ax[0, 1].set_xlabel("lag from reward (s)")
    ax[0, 1].set_ylabel("pupil (% of median) per reward")
    ax[0, 1].legend(fontsize=7)
    band(ax[1, 0], lag_l, mBl, "#2471a3", "lick kernel")
    ax[1, 0].axvline(0, color="k", lw=0.8)
    ax[1, 0].axhline(0, color="0.6", lw=0.6)
    ax[1, 0].set_title("E7c  GLM lick kernel (per lick)", fontsize=10)
    ax[1, 0].set_xlabel("lag from lick (s)")
    ax[1, 0].set_ylabel("pupil (% of median) per lick")
    ax[1, 0].legend(fontsize=7)
    xm = np.arange(len(per_mouse))
    ax[1, 1].bar(xm - 0.2, per_mouse.reward_peak_no_lick_model, 0.4, color="#7f8c8d", label="reward-only model")
    ax[1, 1].bar(xm + 0.2, per_mouse.reward_peak_with_lick_model, 0.4, color="#c0392b", label="licks also in model")
    ax[1, 1].set_xticks(xm)
    ax[1, 1].set_xticklabels([lab(m) for m in per_mouse.index], rotation=30, fontsize=8)
    ax[1, 1].axhline(0, color="0.6", lw=0.6)
    ax[1, 1].set_ylabel("reward kernel 0.5-2.5 s (%)")
    ax[1, 1].set_title(f"E7d  per mouse; median part left {100 * np.nanmedian(keep_frac):.0f} %", fontsize=10)
    ax[1, 1].legend(fontsize=7)
    fig.tight_layout()
    save(fig, d / "E7_pupil_glm.png")
    readme(d, "E7 - pupil GLM: reward vs licking (ideas I06 / I08)", [
        "- Per session (days 3-5, first 910 s): pupil = % of the session median (frames out by the pupil-check rules "
        "held / jump / odd size removed).",
        "- Regressors: rewards (lags -2.1 to +6.1 s) and licks (lags -1.1 to +3.2 s) in 0.264-s bins (FIR), plus slow "
        "drift (cosines, periods >= 60 s). Small ridge on the event columns.",
        "- Reward kernel from a reward-only model vs a reward + lick model. If most of the reward response disappears "
        "when licks are in the model, the plain reward-locked dilation is mostly licking (Yamada & Toda 2022; "
        "deconvolution: Hoeks & Levelt 1993, Wierda et al. 2012).",
        "- Licks and rewards are correlated (rewards come at licks), so the split is approximate.",
        "- Whole-eye mice (9669m black, 9670f black) are in the CSV but not in the means."])
    return E


# ------------------------------------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-date", default=dt.date.today().strftime("%Y%m%d"))
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    background()
    only = {x.strip().upper() for x in a.only.split(",") if x.strip()}
    out = A / "exploratory" / a.out_date
    out.mkdir(parents=True, exist_ok=True)
    man = manifest()
    M = pd.read_csv(A / "master_day_level.csv")
    M = M[M.day.isin(ALL_DAYS)].copy()
    L = pd.read_csv(A / "licks_long.csv", usecols=["mouse", "day", "session", "t_s", "frame", "frame_reliable",
                                                   "counted", "during_valve"])
    R = pd.read_csv(A / "rewards_long.csv", usecols=["mouse", "day", "session", "kind", "t_s", "frame"])
    R = R[R.kind == "reward"]
    res, errs = {}, {}
    E1 = E3 = E6 = None

    def run(key, fn, *args):
        if only and key not in only:
            return None
        try:
            print(f"--- {key}", flush=True)
            return fn(*args)
        except Exception:  # noqa: BLE001
            errs[key] = traceback.format_exc()
            print(errs[key], flush=True)
            return None

    E1 = run("E1", e1, M, L, R, man, out, res)
    if E1 is not None:
        run("E2", e2, M, E1, man, out, res)
    E3 = run("E3", e3, M, L, man, out, res)
    run("E4", e4, man, out, res)
    E6 = run("E6", e6, M, L, E1, man, out, res)
    if E3 is not None:
        run("E5", e5, M, E3, E6, out, res)
    run("E7", e7, M, L, R, man, out, res)
    rj = out / "results.json"
    if only and rj.exists():                       # a partial run keeps the other tests' results
        old = json.loads(rj.read_text(encoding="utf-8"))
        res = {**old.get("results", {}), **res}
        errs = {**{k: v for k, v in old.get("errors", {}).items() if k not in only}, **errs}
    rj.write_text(json.dumps(dict(made=dt.datetime.now().isoformat(timespec="seconds"), tag=TAG, results=res,
                                  errors=errs), indent=2, default=str), encoding="utf-8")
    print(json.dumps(res, indent=1, default=str))
    print("errors:", list(errs) or "none")
    return 1 if errs else 0


if __name__ == "__main__":
    sys.exit(main())
