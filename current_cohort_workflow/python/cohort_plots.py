"""Cohort figures in the style of the previous analysis, from one (or a combined) ALL_mice_longitudinal.csv.

    python cohort_plots.py --longitudinal <run_###>\\ALL_mice_longitudinal.csv --out <figs folder>
                           [--manifest cohort_manifest.csv] [--cohort-col cohort]

Every figure id matches docs/ANALYSIS_METHODS.md (catalog of the previous cohort's deck). Definitions follow the previous
pipeline: lick / reward = 0->1 change of Lick_TTL / Injector_TTL between camera frames, time = CamTime_rel_s,
PR breakpoint = RequirementLast. Analysis rules of the previous analysis:
  phases Pre 3-5, During 6-10, Post 11-13, Withdrawal 14-16, Re-exposure 17-18 (days 1-2 = FR training, not used);
  phase medians leave out the transition days 4, 6, 11, 14, 17 (step05 cfg.transitionDays), sensitivity = all days;
  passive licks on days 6-10 are masked (licks forced invalid), passive rewards are kept (yoked);
  group from the manifest (else: a mouse that is passive on any day = Passive).
Statistics: Mann-Whitney U (Active vs Passive) per phase, Holm-adjusted across phases; ANCOVA outcome ~ Pre + Group.

Writes <out>\\mouse_day_table.csv (one row per mouse-day; the input of every figure) and PNG figures:
  PR01 PR02 PR03 PR06 PR07 PR11 | LK01 LK02 LK05 | TJ02 | PU01 PU03 PU04 PU06 | AS01 AS02 AS03 AS04 | MV15 | QC03
A figure is skipped (and listed) when its data are not there yet (e.g. before day 11 there is no Post phase).
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

PHASES = [("Pre", 3, 5), ("During", 6, 10), ("Post", 11, 13), ("Withdrawal", 14, 16), ("Re-exposure", 17, 18)]
TRANSITION_DAYS = {4, 6, 11, 14, 17}
GCOL = {"Active": "#c0392b", "Passive": "#2c6fbb"}
BOUT_GAP_S = 1.0          # bout = licks separated by < 1 s (analyze_passive_active_dashboard_dec2)
COLS = ["mouse_key", "day_index", "session_idx", "Session_Paradigm", "isPassive", "RequirementLast",
        "Immersion_Latency_s", "CamTime_rel_s", "Lick_TTL", "Injector_TTL", "Trial", "TrialRequirement", "Diameter_px",
        "TST_Pct_Non_moving", "HOT_Pct_Non_moving", "STRAUB_Pct_Non_moving", "cohort"]


def phase_of(d):
    for name, a, b in PHASES:
        if a <= d <= b:
            return name
    return None


# ----------------------------------------------------------------------------- read the frame table once
def read_sessions(path: Path):
    """per session: arrays of time, lick, reward, trial, requirement, pupil (float32) + session constants."""
    hdr = list(pd.read_csv(path, nrows=0).columns)
    use = [c for c in COLS if c in hdr]
    ARR = ["CamTime_rel_s", "Lick_TTL", "Injector_TTL", "Trial", "TrialRequirement", "Diameter_px"]
    ONE = ["isPassive", "RequirementLast", "Immersion_Latency_s", "TST_Pct_Non_moving", "HOT_Pct_Non_moving", "STRAUB_Pct_Non_moving"]
    parts, const = {}, {}
    for ch in pd.read_csv(path, usecols=use, chunksize=1_000_000, low_memory=False):
        for c in COLS:
            if c not in ch:
                ch[c] = np.nan
        ch["cohort"] = ch["cohort"].fillna("cohort").astype(str)
        for c in ARR + ONE:
            ch[c] = pd.to_numeric(ch[c], errors="coerce")
        for k, g in ch.groupby(["cohort", "mouse_key", "day_index", "session_idx"], sort=False):
            parts.setdefault(k, []).append(g[ARR].to_numpy(np.float32))       # only the per-frame numbers are kept
            m = g[ONE].max()
            c0 = const.setdefault(k, dict(paradigm=str(g.Session_Paradigm.iloc[0]), **{c: np.nan for c in ONE}))
            for c in ONE:
                if np.isfinite(m[c]):
                    c0[c] = m[c] if not np.isfinite(c0[c]) else max(c0[c], m[c])
    out = {}
    for k, gl in parts.items():
        A = np.concatenate(gl); c0 = const[k]
        out[k] = dict(t=A[:, 0], lick=np.nan_to_num(A[:, 1]) > 0.5, rew=np.nan_to_num(A[:, 2]) > 0.5, trial=A[:, 3],
                      req=A[:, 4], pupil=A[:, 5], paradigm=c0["paradigm"], passive=c0["isPassive"],
                      req_last=c0["RequirementLast"], ti=c0["Immersion_Latency_s"], tst=c0["TST_Pct_Non_moving"],
                      hot=c0["HOT_Pct_Non_moving"], straub=c0["STRAUB_Pct_Non_moving"])
    return out


def onsets(t, x):
    return t[np.nonzero(np.diff(np.r_[False, x].astype(int)) == 1)[0]]


def bouts(lt, gap=BOUT_GAP_S):
    if len(lt) == 0:
        return 0, 0.0
    br = np.r_[0, np.nonzero(np.diff(lt) >= gap)[0] + 1, len(lt)]
    durs = [lt[b - 1] - lt[a] for a, b in zip(br[:-1], br[1:]) if b - a >= 2]
    return len(durs), float(np.sum(durs))


def mouse_day_table(S, manifest):
    rows, pupil_snips, trials = [], [], []
    for (coh, mk, day, ses), s in sorted(S.items()):
        t = s["t"].astype(float)
        if len(t) < 10:
            continue
        mins = (np.nanmax(t) - np.nanmin(t)) / 60
        lt, rt = onsets(t, s["lick"]), onsets(t, s["rew"])
        nb, bdur = bouts(lt)
        p = s["pupil"].astype(float)
        ok = np.isfinite(p)
        resp = np.nan
        if ok.sum() > 100 and len(rt):
            ss = []
            for r in rt:
                w = (t >= r - 2) & (t <= r + 2) & ok
                if w.sum() < 60:
                    continue
                base = np.nanmean(p[w & (t < r)]); post = np.nanmean(p[w & (t >= r + 0.5) & (t <= r + 2)])
                ss.append(post - base)
                tt = t[w] - r
                pupil_snips.append(dict(cohort=coh, mouse_key=mk, day_index=int(day), t=np.round(tt / 0.2) * 0.2, y=p[w] - base))
            resp = np.nanmean(ss) if ss else np.nan
        # per-trial lick density over the full trial window (licks / s between two reward onsets), by PR level
        if len(rt) >= 2:
            for j in range(len(rt) - 1):
                n = np.sum((lt > rt[j]) & (lt <= rt[j + 1]))
                rq = s["req"][np.searchsorted(t, rt[j])] if np.searchsorted(t, rt[j]) < len(t) else np.nan
                trials.append(dict(cohort=coh, mouse_key=mk, day_index=int(day), trial=j + 1, requirement=float(rq),
                                   lick_density_hz=n / max(1e-3, rt[j + 1] - rt[j])))
        rows.append(dict(cohort=coh, mouse_key=mk, day_index=int(day), session_idx=int(ses), paradigm=s["paradigm"],
                         is_passive_session=s["passive"], RequirementLast=s["req_last"], minutes=mins,
                         lick_n=len(lt), lick_freq_per_min=len(lt) / mins if mins > 0 else np.nan,
                         rew_n=len(rt), rew_freq_per_min=len(rt) / mins if mins > 0 else np.nan,
                         req_speed_per_min=s["req_last"] / mins if mins > 0 else np.nan,
                         bout_n=nb, bout_totalDur_s=bdur, pupil_mean_px=np.nanmean(p) if ok.any() else np.nan,
                         pupil_coverage=ok.mean(), pupil_reward_resp_px=resp, Immersion_Latency_s=s["ti"],
                         TST_Pct_Non_moving=s["tst"], HOT_Pct_Non_moving=s["hot"], STRAUB_Pct_Non_moving=s["straub"]))
    D = pd.DataFrame(rows)
    # one row per mouse-day: duplicate sessions (Step 02 counts every concat_out) collapse to the first
    D = D.sort_values(["cohort", "mouse_key", "day_index", "session_idx"]).drop_duplicates(
        ["cohort", "mouse_key", "day_index", "lick_n", "rew_n", "minutes"])
    D = D.groupby(["cohort", "mouse_key", "day_index"], as_index=False).first()
    grp = D.groupby(["cohort", "mouse_key"]).is_passive_session.max().map({1.0: "Passive", 0.0: "Active"})
    D["group"] = [grp.get((c, m), "Active") for c, m in zip(D.cohort, D.mouse_key)]
    D["sex"], D["pair_id"] = "", ""
    if manifest is not None:
        M = manifest.set_index("mouse_key")
        for col in ("group", "sex", "pair_id"):
            if col in M:
                D[col] = [M[col].get(m, D.loc[i, col]) if pd.notna(M[col].get(m, np.nan)) else D.loc[i, col]
                          for i, m in zip(D.index, D.mouse_key)]
    D["phase"] = D.day_index.map(phase_of)
    # passive licking on days 6-10 is not natural licking (lick line locked out): mask lick-derived metrics
    mask = (D.group == "Passive") & D.day_index.between(6, 10)
    for c in ("lick_n", "lick_freq_per_min", "bout_n", "bout_totalDur_s"):
        D.loc[mask, c] = np.nan
    T = pd.DataFrame(trials)
    if len(T):
        T = T.merge(D[["cohort", "mouse_key", "day_index", "group", "phase"]], on=["cohort", "mouse_key", "day_index"])
        T = T[~((T.group == "Passive") & T.day_index.between(6, 10))]
    P = pupil_snips
    return D, T, P


# ----------------------------------------------------------------------------- helpers
def holm(p):
    p = np.asarray(p, float); n = np.sum(np.isfinite(p)); o = np.argsort(np.where(np.isfinite(p), p, np.inf))
    adj = np.full(p.size, np.nan); run = 0
    for r, i in enumerate(o):
        if not np.isfinite(p[i]):
            continue
        run = max(run, min(1, (n - r) * p[i])); adj[i] = run
    return adj


def mwu(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return stats.mannwhitneyu(a, b, alternative="two-sided").pvalue


def stars(p):
    return "" if not np.isfinite(p) else ("***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s.")


def phase_median(D, col, exclude_transition=True):
    X = D[D.day_index >= 3].copy()
    if exclude_transition:
        X = X[~X.day_index.isin(TRANSITION_DAYS)]
    X = X.dropna(subset=["phase"])
    return X.groupby(["cohort", "mouse_key", "group", "phase"], as_index=False)[col].median()


def save(fig, out, name, made):
    fig.tight_layout(); fig.savefig(out / f"{name}.png", dpi=130); plt.close(fig); made.append(name)


def box_by_phase(ax, M, col, phases, title):
    pos, pv = 0, []
    for ph in phases:
        for g, dx in (("Active", -0.18), ("Passive", 0.18)):
            v = M[(M.phase == ph) & (M.group == g)][col].dropna().to_numpy()
            if len(v):
                ax.boxplot(v, positions=[pos + dx], widths=0.3, showfliers=False, patch_artist=True,
                           boxprops=dict(facecolor=GCOL[g], alpha=0.25), medianprops=dict(color=GCOL[g]))
                ax.plot(np.full(len(v), pos + dx) + np.random.default_rng(1).uniform(-0.06, 0.06, len(v)), v, "o",
                        color=GCOL[g], ms=4)
        pv.append(mwu(M[(M.phase == ph) & (M.group == "Active")][col], M[(M.phase == ph) & (M.group == "Passive")][col]))
        pos += 1
    adj = holm(pv)
    top = np.nanmax(M[col]) if M[col].notna().any() else 1
    for i, (p, q) in enumerate(zip(pv, adj)):
        if np.isfinite(p):
            ax.text(i, top * 1.05 if top > 0 else top + 0.5, f"p={p:.3g}\n{stars(q)} (Holm)", ha="center", fontsize=7)
    ax.set_xticks(range(len(phases))); ax.set_xticklabels(phases); ax.set_title(title, fontsize=10, loc="left")
    ax.axhline(0, color="0.7", lw=0.5)


# ----------------------------------------------------------------------------- figures
def figures(D, T, P, out):
    made, skipped = [], []
    have = lambda ph: (D.phase == ph).any()
    present = [ph for ph, _, _ in PHASES if have(ph)]
    D4 = D[D.day_index >= 3]

    # PR01 daily PR spaghetti by group
    if D4.RequirementLast.notna().any():
        fig, axs = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
        for ax, g in zip(axs, ("Active", "Passive")):
            X = D4[D4.group == g]
            X = X[~((g == "Passive") & X.day_index.between(6, 10))]
            for mk, m in X.groupby("mouse_key"):
                ax.plot(m.day_index, m.RequirementLast, "-o", color=GCOL[g], alpha=0.35, ms=3, lw=0.8)
            mm = X.groupby("day_index").RequirementLast.mean()
            ax.plot(mm.index, mm.values, "-", color="k", lw=2.2, label="group mean")
            for _, a, b in PHASES[1:]:
                ax.axvline(a - 0.5, color="0.8", ls=":")
            ax.set_title(f"PR01  {g}: daily PR (RequirementLast)", loc="left", fontsize=10); ax.set_xlabel("day")
        axs[0].set_ylabel("RequirementLast"); save(fig, out, "PR01_daily_PR_spaghetti", made)
    else:
        skipped.append("PR01 (no PR days yet)")

    # PR02 phase medians (transition days excluded) + PR03 delta from Pre
    M = phase_median(D, "RequirementLast")
    ph = [p for p in present if p != "During" or (M[(M.phase == "During") & (M.group == "Active")].shape[0])]
    if len(M):
        fig, ax = plt.subplots(figsize=(1.9 * max(2, len(ph)) + 2, 4))
        box_by_phase(ax, M[M.phase.isin(ph)], "RequirementLast", ph, "PR02  PR phase medians (transition days 4,6,11,14,17 excluded)")
        ax.set_ylabel("median RequirementLast"); save(fig, out, "PR02_PR_phase_medians", made)
    pre = M[M.phase == "Pre"].set_index(["cohort", "mouse_key"]).RequirementLast
    if len(pre) and len(present) > 1:
        M["delta_from_pre"] = M.RequirementLast - [pre.get((c, m), np.nan) for c, m in zip(M.cohort, M.mouse_key)]
        fig, ax = plt.subplots(figsize=(1.9 * len(ph) + 2, 4))
        box_by_phase(ax, M[M.phase.isin(ph) & (M.phase != "Pre")], "delta_from_pre", [p for p in ph if p != "Pre"],
                     "PR03  delta PR from Pre (phase median - Pre median)")
        ax.set_ylabel("delta RequirementLast"); save(fig, out, "PR03_delta_PR_from_Pre", made)
    else:
        skipped.append("PR03 (needs Pre and a later phase)")

    # PR06 trajectory day 11-18 mean +- SEM
    X = D[D.day_index.between(11, 18)]
    if len(X):
        fig, ax = plt.subplots(figsize=(8, 4)); pv = []
        days = sorted(X.day_index.unique())
        for g in ("Active", "Passive"):
            st = X[X.group == g].groupby("day_index").RequirementLast.agg(["mean", "sem"])
            ax.errorbar(st.index, st["mean"], yerr=st["sem"], color=GCOL[g], marker="o", capsize=3, label=g)
        pv = [mwu(X[(X.day_index == d) & (X.group == "Active")].RequirementLast, X[(X.day_index == d) & (X.group == "Passive")].RequirementLast) for d in days]
        for d, q in zip(days, holm(pv)):
            ax.text(d, ax.get_ylim()[1] * 0.97, stars(q), ha="center", fontsize=8)
        ax.axvline(13.5, color="0.7", ls=":"); ax.axvline(16.5, color="0.7", ls=":")
        ax.set_title("PR06  Post -> Withdrawal -> Re-exposure (mean +- SEM; stars = Holm-adjusted MWU)", loc="left", fontsize=10)
        ax.set_xlabel("day"); ax.set_ylabel("RequirementLast"); ax.legend(); save(fig, out, "PR06_PR_trajectory_post_to_reexposure", made)
    else:
        skipped.append("PR06 (needs day 11-18)")

    # PR07 + PR11 baseline vs morphine effect, ANCOVA outcome ~ Pre + Group
    lev = M.pivot_table(index=["cohort", "mouse_key", "group"], columns="phase", values="RequirementLast").reset_index()
    outs = [c for c in ("Post", "Re-exposure") if c in lev]
    if "Pre" in lev and outs:
        fig, axs = plt.subplots(1, len(outs), figsize=(5.5 * len(outs), 4.2), squeeze=False)
        for ax, oc in zip(axs[0], outs):
            L = lev.dropna(subset=["Pre", oc]).copy(); L["eff"] = L[oc] - L.Pre
            for g in ("Active", "Passive"):
                q = L[L.group == g]; ax.plot(q.Pre, q.eff, "o", color=GCOL[g], label=g)
                for _, r in q.iterrows():
                    ax.annotate(r.mouse_key, (r.Pre, r.eff), fontsize=6, alpha=0.7)
            ax.axhline(0, color="0.6", ls="--", lw=0.8)
            # ANCOVA (common slope): eff ~ Pre + Group
            if L.group.nunique() == 2 and len(L) >= 5:
                X1 = np.c_[np.ones(len(L)), L.Pre, (L.group == "Active").astype(float)]
                b, res, *_ = np.linalg.lstsq(X1, L.eff, rcond=None)
                X0 = X1[:, :2]; b0, *_ = np.linalg.lstsq(X0, L.eff, rcond=None)
                rss1 = np.sum((L.eff - X1 @ b) ** 2); rss0 = np.sum((L.eff - X0 @ b0) ** 2)
                df = len(L) - 3; F = (rss0 - rss1) / (rss1 / df) if rss1 > 0 else np.nan
                pF = 1 - stats.f.cdf(F, 1, df) if np.isfinite(F) else np.nan
                xs = np.linspace(L.Pre.min(), L.Pre.max(), 10)
                for g, gv in (("Active", 1), ("Passive", 0)):
                    ax.plot(xs, b[0] + b[1] * xs + b[2] * gv, color=GCOL[g], lw=1)
                ax.set_title(f"PR07/PR11  {oc} - Pre vs Pre (ANCOVA group p={pF:.3g}; MWU p={mwu(L[L.group=='Active'].eff, L[L.group=='Passive'].eff):.3g})",
                             loc="left", fontsize=9)
            ax.set_xlabel("Pre median PR"); ax.set_ylabel(f"{oc} - Pre")
        axs[0][0].legend(); save(fig, out, "PR07_PR11_baseline_vs_morphine_effect", made)
    else:
        skipped.append("PR07/PR11 (needs Pre and Post or Re-exposure)")

    # LK01 cumulative licking within session by epoch (from the frame table is too heavy: use per-minute rate)
    # LK02 lick frequency delta from Pre, LK05 bout total duration delta from Pre
    for fid, col, lab in (("LK02", "lick_freq_per_min", "licks / min"), ("LK05", "bout_totalDur_s", "bout total duration (s)")):
        Mm = phase_median(D, col)
        pre = Mm[Mm.phase == "Pre"].set_index(["cohort", "mouse_key"])[col]
        later = [p for p in present if p != "Pre"]
        if len(pre) and later:
            Mm["d"] = Mm[col] - [pre.get((c, m), np.nan) for c, m in zip(Mm.cohort, Mm.mouse_key)]
            fig, ax = plt.subplots(figsize=(1.9 * len(later) + 2, 4))
            box_by_phase(ax, Mm[Mm.phase.isin(later)], "d", later, f"{fid}  {lab}: phase median - Pre (passive days 6-10 masked)")
            ax.set_ylabel(f"delta {lab}"); save(fig, out, f"{fid}_{col}_delta_from_Pre", made)
        else:
            skipped.append(f"{fid} (needs Pre and a later phase)")

    # TJ02 lick density over the full trial window vs PR level (requirement), by phase
    if len(T) and T.phase.notna().any():
        phs = [p for p in present if (T.phase == p).any()]
        fig, axs = plt.subplots(1, len(phs), figsize=(3.6 * len(phs), 3.6), sharey=True, squeeze=False)
        for ax, p in zip(axs[0], phs):
            for g in ("Active", "Passive"):
                q = T[(T.phase == p) & (T.group == g)].copy()
                q["lvl"] = q.requirement.clip(upper=40)
                st = q.groupby("lvl").lick_density_hz.agg(["mean", "sem", "count"])
                st = st[st["count"] >= 2]
                if len(st):
                    ax.plot(st.index, st["mean"], color=GCOL[g], label=g); ax.fill_between(st.index, st["mean"] - st["sem"], st["mean"] + st["sem"], color=GCOL[g], alpha=0.2)
            ax.set_title(f"TJ02  {p}", loc="left", fontsize=10); ax.set_xlabel("PR level (TrialRequirement)")
        axs[0][0].set_ylabel("licks / s in trial window"); axs[0][0].legend(); save(fig, out, "TJ02_lick_density_vs_PR_level", made)
    else:
        skipped.append("TJ02 (needs PR trials)")

    # PU01 session pupil trajectory, PU03 % change from Pre, PU04 phase means
    Pp = D[(D.day_index >= 3) & (D.pupil_coverage > 0.5)]
    if Pp.pupil_mean_px.notna().any():
        fig, axs = plt.subplots(1, 2, figsize=(12, 4))
        for g in ("Active", "Passive"):
            q = Pp[Pp.group == g]
            for mk, m in q.groupby("mouse_key"):
                axs[0].plot(m.day_index, m.pupil_mean_px, "-", color=GCOL[g], alpha=0.3)
            st = q.groupby("day_index").pupil_mean_px.agg(["mean", "sem"]); axs[0].errorbar(st.index, st["mean"], st["sem"], color=GCOL[g], lw=2, label=g)
            base = q[q.phase == "Pre"].groupby("mouse_key").pupil_mean_px.mean()
            q = q.assign(pct=100 * (q.pupil_mean_px / q.mouse_key.map(base) - 1))
            st = q.groupby("day_index").pct.agg(["mean", "sem"]); axs[1].errorbar(st.index, st["mean"], st["sem"], color=GCOL[g], lw=2, label=g)
        axs[0].set_title("PU01  session pupil diameter (px), sessions with > 50 % tracked frames", loc="left", fontsize=9)
        axs[1].set_title("PU03  pupil % change from the mouse's Pre mean", loc="left", fontsize=9)
        for ax in axs:
            ax.set_xlabel("day"); ax.legend()
        save(fig, out, "PU01_PU03_pupil_trajectory", made)
        Mp = phase_median(Pp, "pupil_mean_px", exclude_transition=False)
        if len(present):
            fig, ax = plt.subplots(figsize=(1.9 * len(present) + 2, 4))
            box_by_phase(ax, Mp, "pupil_mean_px", present, "PU04  pupil phase median (px)")
            save(fig, out, "PU04_pupil_phase_median", made)
    else:
        skipped.append("PU01/PU03/PU04 (no sessions with > 50 % tracked pupil frames)")

    # PU06 reward-locked pupil (-2..2 s, -2..0 s baseline) by phase and group
    if P:
        S = pd.DataFrame([dict(cohort=s["cohort"], mouse_key=s["mouse_key"], day_index=s["day_index"], t=t, y=y)
                          for s in P for t, y in zip(s["t"], s["y"])])
        S = S.merge(D[["cohort", "mouse_key", "day_index", "group", "phase"]], on=["cohort", "mouse_key", "day_index"])
        phs = [p for p in present if (S.phase == p).any()]
        if phs:
            fig, axs = plt.subplots(1, len(phs), figsize=(3.6 * len(phs), 3.4), sharey=True, squeeze=False)
            for ax, p in zip(axs[0], phs):
                for g in ("Active", "Passive"):
                    q = S[(S.phase == p) & (S.group == g)]
                    if len(q):
                        mm = q.groupby(["mouse_key", "t"]).y.mean().groupby("t").agg(["mean", "sem"])
                        ax.plot(mm.index, mm["mean"], color=GCOL[g], label=g); ax.fill_between(mm.index, mm["mean"] - mm["sem"], mm["mean"] + mm["sem"], color=GCOL[g], alpha=0.2)
                ax.axvline(0, color="k", lw=0.6); ax.set_title(f"PU06  {p}", loc="left", fontsize=10); ax.set_xlabel("s from reward")
            axs[0][0].set_ylabel("pupil change (px, -2..0 s baseline)"); axs[0][0].legend(); save(fig, out, "PU06_reward_locked_pupil", made)
    else:
        skipped.append("PU06 (no tracked pupil around rewards)")

    # AS01 tail immersion per day; AS02 HP, AS03 TST, AS04 Straub (delta from Pre / raw by phase)
    if D.Immersion_Latency_s.notna().any():
        fig, ax = plt.subplots(figsize=(8, 4)); X = D.dropna(subset=["Immersion_Latency_s"]); days = sorted(X.day_index.unique())
        for g in ("Active", "Passive"):
            st = X[X.group == g].groupby("day_index").Immersion_Latency_s.agg(["mean", "sem"]); ax.errorbar(st.index, st["mean"], st["sem"], color=GCOL[g], marker="o", capsize=3, label=g)
        pv = [mwu(X[(X.day_index == d) & (X.group == "Active")].Immersion_Latency_s, X[(X.day_index == d) & (X.group == "Passive")].Immersion_Latency_s) for d in days]
        for d, q in zip(days, holm(pv)):
            ax.text(d, ax.get_ylim()[1] * 0.97, stars(q), ha="center", fontsize=8)
        ax.set_title("AS01  tail immersion latency (s), mean +- SEM, Holm-adjusted MWU", loc="left", fontsize=10); ax.legend(); ax.set_xlabel("day")
        save(fig, out, "AS01_tail_immersion_by_day", made)
    else:
        skipped.append("AS01 (no Immersion_Latency_s yet)")
    for fid, col in (("AS02", "HOT_Pct_Non_moving"), ("AS03", "TST_Pct_Non_moving"), ("AS04", "STRAUB_Pct_Non_moving")):
        if D[col].notna().any():
            Mm = phase_median(D, col, exclude_transition=False)
            phs = [p for p in present if (Mm.phase == p).any()]
            fig, ax = plt.subplots(figsize=(1.9 * len(phs) + 2, 4)); box_by_phase(ax, Mm, col, phs, f"{fid}  {col} by phase")
            save(fig, out, f"{fid}_{col}_by_phase", made)
        else:
            skipped.append(f"{fid} (no {col} yet)")

    # MV15 behavioural-space divergence: cosine / Euclidean distance of the Active and Passive centroids per phase
    feats = ["RequirementLast", "lick_freq_per_min", "bout_n", "bout_totalDur_s", "rew_n", "rew_freq_per_min", "req_speed_per_min"]
    Z = D[(D.day_index >= 3)].dropna(subset=["phase"]).copy()
    if len(Z) and Z.group.nunique() == 2:
        Zs = Z[feats].apply(lambda c: (c - c.mean()) / (c.std() if c.std() > 0 else 1))
        Z[feats] = Zs
        rows = []
        for p in present:
            a = Z[(Z.phase == p) & (Z.group == "Active")][feats].mean(); b = Z[(Z.phase == p) & (Z.group == "Passive")][feats].mean()
            ok = a.notna() & b.notna()
            if ok.sum() >= 2:
                cos = float(np.dot(a[ok], b[ok]) / (np.linalg.norm(a[ok]) * np.linalg.norm(b[ok]) + 1e-12)); dist = float(np.linalg.norm(a[ok] - b[ok]))
                rows.append((p, cos, dist, int(ok.sum())))
        if rows:
            R = pd.DataFrame(rows, columns=["phase", "cosine", "euclidean", "n_features"])
            fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
            axs[0].bar(R.phase, R.cosine, color="0.4"); axs[0].set_title("MV15  cosine similarity of group centroids", loc="left", fontsize=10)
            axs[1].bar(R.phase, R.euclidean, color="0.4"); axs[1].set_title("MV15  Euclidean distance (z-scored features)", loc="left", fontsize=10)
            save(fig, out, "MV15_group_centroid_divergence", made)
            R.to_csv(out / "MV15_group_centroid_divergence.csv", index=False)
    else:
        skipped.append("MV15 (needs both groups)")

    # QC03 coverage: mouse x day
    cov = D.pivot_table(index="mouse_key", columns="day_index", values="rew_n", aggfunc="count").fillna(0)
    fig, ax = plt.subplots(figsize=(0.45 * max(6, cov.shape[1]) + 3, 0.35 * max(4, cov.shape[0]) + 1.5))
    ax.imshow(cov.values > 0, aspect="auto", cmap="Greens", vmin=0, vmax=1)
    ax.set_xticks(range(cov.shape[1])); ax.set_xticklabels(cov.columns); ax.set_yticks(range(cov.shape[0])); ax.set_yticklabels(cov.index, fontsize=7)
    ax.set_title("QC03  sessions present (mouse x day)", loc="left", fontsize=10); save(fig, out, "QC03_coverage", made)
    return made, skipped


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--longitudinal", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--manifest", type=Path, help="CSV mouse_key,group,sex,pair_id (optional)")
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    S = read_sessions(a.longitudinal)
    man = pd.read_csv(a.manifest) if a.manifest else None
    D, T, P = mouse_day_table(S, man)
    D.to_csv(a.out / "mouse_day_table.csv", index=False)
    if len(T):
        T.to_csv(a.out / "trial_table.csv", index=False)
    made, skipped = figures(D, T, P, a.out)
    print(f"{len(D)} mouse-days, {D.mouse_key.nunique()} mice, days {D.day_index.min()}-{D.day_index.max()}")
    print("made:", ", ".join(made))
    print("skipped:", "; ".join(skipped) or "-")


if __name__ == "__main__":
    main()
