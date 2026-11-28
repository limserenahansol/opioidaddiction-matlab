#!/usr/bin/env python3
"""
pr_relationships.py - does a higher PR go with a higher tail immersion, TST, hot plate, Straub tail or pupil?
(added 2026-11-23; user: "analysis for now and for future as well 1) relationship comparison PR vs TI , PR vs TST ,
PR vs HOT , PR vs Straub ,PR vs Pupil to see higher PR session or mouse has higher TI etc... please do this for current
cohort daily and also for previous cohort data as well !"; then "need to label passiive and active and morphine and
natural reward separately"). Report only; exploratory (small n, repeated measures, no correction for multiple tests).
Backups: pr_relationships_before_20261123_split.py (first version, groups mixed), _before_20261123_review.py (before
the fixes of the independent review workflow wf_ca5ec6d6-fa0, 2026-11-23).

Both cohorts go through the SAME mouse-day table (cohort_plots.read_sessions + mouse_day_table on the original Step 02
ALL_mice_longitudinal.csv), so the definitions match:
  reward    by day, as ti_table_days.PHASES: natural reward = Baseline 1-5 + Withdrawal 14-16; morphine = Yoked 6-10 +
            Equal Access 11-13 + Re-exposure 17-18
  x = PR    own PR sessions (Active mice from day 3; Passive mice on days 3-5 and from day 11): RequirementLast
            (= rewards + 1, the deck's PR; on morphine days = drops earned + 1: the drop count is a dose proxy, the dose itself was never measured)
            Passive mice on the Yoked days 6-10 have no PR of their own (the partner's replay delivers the drops):
            x = drops delivered in the analysed session + 1 (= the PR of the replayed partner session, so a pair shares
            x). Current cohort: the routine's PR table (partner_pr; checked = delivered + 1); previous cohort: reward
            onsets in the Step 02 table + 1, with the exceptions in analysis/prev_cohort/prev_yoked_x_notes.csv (6100
            day 10 replayed the partner's day-9 log; 4 days had an aborted replay with extra drops before). Drawn as
            open squares and tested on their own, never pooled with an own PR
  TI        Immersion_Latency_s                      TST / HOT  TST_Pct_Non_moving / HOT_Pct_Non_moving (%)
  Straub    tail-up % = 100 - STRAUB_Pct_Non_moving (hand scores; the user's rule; the old deck's VTTA_PK1 plot used
            Non_moving, so its sign is reversed); current cohort: straub_dayN.csv straub_final_pct (automatic v9.1, or
            the user's hand score over the counting frames), because Step 02 has no automatic Straub
  pupil     session mean pupil (px) as % of the mouse's own days 1-5 mean (px differ between mice: camera geometry);
            pupil change at reward = mean pupil 0.5-2 s after a reward minus the 2 s before (-2..0 s, the cohort_plots
            PU06 window), as % of that days 1-5 mean (NOT the results picture's per-reward 0.5-2.5 s vs -1..0 s %)
            whole-eye mice (eye damage: 9669m black = Active, 9670f black = Passive) are left out of the pupil panels
Session groups, each tested on its own: Active · natural reward, Passive · natural reward, Both groups · natural reward,
Active · morphine, Passive · morphine (own PR, from day 11), Both groups · morphine (own PR), Passive · morphine yoked;
csv only: natural reward split into Baseline 3-5 and Withdrawal 14-16 (they can disagree).
Statistics:
  mouse     do mice with a higher PR have a higher X?  Spearman over mice of the mouse means of DAY-CENTRED values
            (each session minus that day's mean over the mice of the same session group, days with >= 2 mice; the
            lab's main test is mouse-level); exact permutation p for <= 8 mice (4 mice: the smallest possible p is
            0.083; 3 mice: no test), t approximation above. The big dots in the pictures are the plain mouse means
  same day  on the same day, does the session with the higher PR have the higher X?  both centred on each day's mean
            (removes the shared change over days, e.g. the morphine effect); p from a WHOLE-MOUSE permutation (each
            mouse's x series moves to another mouse of the same group; exact for <= 40320 relabelings, else 5000 random),
            so a mouse's repeated sessions are not counted as independent. When every mouse has one session (one-day
            measures) the same-day r uses the same points as the mouse level: then it is left empty ("= mouse level")
  both-groups rows are centred on day x group (each mouse vs its own group that day), so a difference between Active
            and Passive mice cannot make the correlation; the per-group rows are the main ones
  csv only  sameday_r_cage (centred on day x cage as well: the current cohort's cage = sex, 9669m males, 9670f females);
            leave-one-mouse-out / leave-one-day-out for the p < 0.05 results of day set 'all' (loo_*);
            within_r_naive (centred on each mouse; carries the change over days, naive t p), within_dayadj_r (mouse +
            day fixed effects, descriptive), pooled_rho_naive (all sessions as independent, descriptive)
Day sets: every statistic for each of the lab's 8 day sets (tools/day_sets.csv); the pictures show 'all'. No set leaves
out an assay day (5, 10, 13, 16, 18) or a Straub day, so TST / hot plate / Straub results are the same in every set.
Previous cohort calendar: day 8 ran on 01-07 (cages 0911, 6100) or 01-12 (0910, 6099), and 0910 / 6099 ran day 9 the
same day about 8 h later; centring is by day index. 8 Passive mice replayed 6 partners (0910_black and 6100_black twice).
usage:
  python pr_relationships.py --build-prev <prev ALL_mice_longitudinal.csv>   (writes analysis/prev_cohort/
         prev_cohort_mouse_day_table.csv; an older cache is moved to _superseded_<time>; the 3.9 GB file is only read)
  python pr_relationships.py [--current <current ALL_mice_longitudinal.csv>]  (default: the newest Step 02 run)
Outputs (analysis/days): pr_relationships_current_days1-N.png, pr_relationships_previous_cohort.png,
  pr_relationships_daysets_days1-N.png, pr_relationships_stats_days1-N.csv, pr_relationships_data_days1-N.csv,
  pr_relationships_inputs_days1-N.json (what the outputs were made from; the deck re-runs only when it changed)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import sys
from datetime import datetime
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

T = Path(__file__).resolve().parent
P = T.parent
DAYS_DIR = P / "analysis" / "days"
PREV_DIR = P / "analysis" / "prev_cohort"
PREV_CACHE = PREV_DIR / "prev_cohort_mouse_day_table.csv"
PREV_PAIRS = PREV_DIR / "prev_cohort_pairs.csv"
PREV_NOTES = PREV_DIR / "prev_yoked_x_notes.csv"
sys.path.insert(0, str(T))
import ti_table_days as TTD  # noqa: E402  (phase + reward of each day: the routine's one table)

sys.path.insert(1, str(T / "opioidaddiction-matlab" / "current_cohort_workflow" / "python"))
WHOLE_EYE = {"9669_black", "9670_black", "9669m_black", "9670f_black"}
CAGE_SEX = {"current": {"9669": "M", "9670": "F"}}            # current cohort: cage = sex
GCOL = {"Active": "#c0392b", "Passive": "#2c6fbb"}       # as cohort_plots.py / ti_table_days.py
NAT, MOR = "natural reward", "morphine"
EXACT_MAX, N_RAND = 40320, 5000
MEAS = [("ti", "Tail immersion (s)"), ("tst", "TST non-moving (%)"), ("hot", "Hot plate non-moving (%)"),
        ("straub", "Straub tail up (%)"), ("pupil_rel", "Pupil size (% of own days 1-5)"),
        ("pupil_resp_rel", "Pupil change at reward (0.5-2 vs -2-0 s, %)")]
SHORT = {"ti": "TI", "tst": "TST", "hot": "Hot plate", "straub": "Straub", "pupil_rel": "Pupil size",
         "pupil_resp_rel": "Pupil at reward"}
# (session group, reward, who, phase): who = Active / Passive / both = own PR sessions; yoked = Passive on the Yoked
# days; phase None = every day of that reward. The first 7 are the main groups (pictures, deck); the rest csv only
SETS = [("Active · natural reward", NAT, "Active", None),
        ("Passive · natural reward", NAT, "Passive", None),
        ("Both groups · natural reward", NAT, "both", None),
        ("Active · morphine", MOR, "Active", None),
        ("Passive · morphine (own PR)", MOR, "Passive", None),
        ("Both groups · morphine (own PR)", MOR, "both", None),
        ("Passive · morphine yoked (delivered + 1)", MOR, "yoked", None)]
MAIN = [s[0] for s in SETS]
for _who in ("Active", "Passive", "both"):
    for _ph in ("Baseline", "Withdrawal"):
        SETS.append((f"{'Both groups' if _who == 'both' else _who} · natural reward · {_ph}", NAT, _who, _ph))
SET_TINY = {"Active · natural reward": "A·nat", "Passive · natural reward": "P·nat",
            "Both groups · natural reward": "both·nat", "Active · morphine": "A·mor",
            "Passive · morphine (own PR)": "P·mor", "Both groups · morphine (own PR)": "both·mor",
            "Passive · morphine yoked (delivered + 1)": "P yoked"}
SET_SHORT = {"Active · natural reward": "Active·natural", "Passive · natural reward": "Passive·natural",
             "Both groups · natural reward": "both·natural", "Active · morphine": "Active·morphine",
             "Passive · morphine (own PR)": "Passive·morphine", "Both groups · morphine (own PR)": "both·morphine",
             "Passive · morphine yoked (delivered + 1)": "Passive yoked"}


def mouse_day(path: Path, manifest: Path | None = None) -> pd.DataFrame:
    import cohort_plots as CP          # the validated per-session reader (previous cohort deck numbers reproduced)
    S = CP.read_sessions(path)
    # groups from the manifest when given (review 2026-11-24: before the first replay day every mouse read 'Active')
    man = pd.read_csv(manifest) if manifest is not None and Path(manifest).exists() else None
    D, _, _ = CP.mouse_day_table(S, man)
    keep = ["cohort", "mouse_key", "day_index", "paradigm", "is_passive_session", "group", "RequirementLast", "rew_n",
            "minutes", "pupil_mean_px", "pupil_coverage", "pupil_reward_resp_px", "Immersion_Latency_s",
            "TST_Pct_Non_moving", "HOT_Pct_Non_moving", "STRAUB_Pct_Non_moving"]
    return D[[c for c in keep if c in D.columns]].copy()


def newest_current() -> Path | None:
    runs = sorted((P / "step02_base" / "longitudinal_outputs").glob("run_*/ALL_mice_longitudinal.csv"),
                  key=lambda p: p.stat().st_mtime)
    return runs[-1] if runs else None


def day_sets() -> list[tuple[str, set]]:
    f = T / "day_sets.csv"
    if not f.exists():
        return [("all", set())]
    out = []
    for _, r in pd.read_csv(f, dtype=str).fillna("").iterrows():
        out.append((r["set"], {int(x) for x in str(r["exclude_days"]).split(";") if x.strip().isdigit()}))
    return out


def prepare(D: pd.DataFrame, cohort: str) -> pd.DataFrame:
    D = D.copy()
    D["cohort"] = cohort
    D["day"] = D.day_index.astype(int)
    D["cage"] = D.mouse_key.astype(str).str[:4]
    ph = [TTD.phase(int(d)) for d in D.day]
    D["phase"] = [p[0] for p in ph]
    D["reward"] = [p[1] for p in ph]
    D["pr"] = pd.to_numeric(D.RequirementLast, errors="coerce")
    par = D.paradigm.astype(str).str.lower()
    passive_ses = pd.to_numeric(D.is_passive_session, errors="coerce").fillna(0) >= 0.5
    D["own"] = ~passive_ses & ~par.str.contains("fixed") & (D.day >= 3) & D.pr.gt(0)
    D["yoked"] = passive_ses & D.reward.eq(MOR) & D.group.eq("Passive")
    deliv = (pd.to_numeric(D["rew_n"], errors="coerce") + 1) if "rew_n" in D else pd.Series(np.nan, index=D.index)
    D["x"] = np.where(D.own, D.pr, np.nan)
    D["x_source"] = np.where(D.own, "own PR (RequirementLast)", "")
    D.loc[D.yoked, "x"] = deliv[D.yoked]
    D.loc[D.yoked & deliv.notna(), "x_source"] = "drops delivered + 1 (Step 02 reward onsets + 1)"
    D["cat"] = np.select([D.own & D.group.eq("Active"), D.own & D.group.eq("Passive"), D.yoked],
                         ["Active", "Passive", "Passive yoked"], "")
    D["ti"] = pd.to_numeric(D.Immersion_Latency_s, errors="coerce")
    D["tst"] = pd.to_numeric(D.TST_Pct_Non_moving, errors="coerce")
    D["hot"] = pd.to_numeric(D.HOT_Pct_Non_moving, errors="coerce")
    D["straub"] = 100 - pd.to_numeric(D.STRAUB_Pct_Non_moving, errors="coerce")
    cov_ok = pd.to_numeric(D.pupil_coverage, errors="coerce") > 0.5
    px = pd.to_numeric(D.pupil_mean_px, errors="coerce").where(cov_ok)
    base = px.where(D.day <= 5).groupby(D.mouse_key).transform("mean")
    D["pupil_rel"] = 100 * px / base
    D["pupil_resp_rel"] = (100 * pd.to_numeric(D.pupil_reward_resp_px, errors="coerce") / base).where(cov_ok)
    eye = D.mouse_key.astype(str).isin(WHOLE_EYE)
    D.loc[eye, ["pupil_rel", "pupil_resp_rel"]] = np.nan
    return D


def previous_notes(D: pd.DataFrame) -> pd.DataFrame:
    """previous cohort: partner of each Passive mouse (replay logs) and the Yoked-day exceptions (notes csv)"""
    D["partner"] = ""
    if PREV_PAIRS.exists():
        pr = pd.read_csv(PREV_PAIRS, dtype=str).drop_duplicates("passive_mouse_key").set_index("passive_mouse_key")
        D["partner"] = D.mouse_key.map(pr.partner_mouse_key).fillna("")
    if PREV_NOTES.exists():
        for _, r in pd.read_csv(PREV_NOTES, dtype={"mouse_key": str}).iterrows():
            m = (D.mouse_key == r.mouse_key) & (D.day == int(r.day))
            D.loc[m, "x_source"] = D.loc[m, "x_source"].astype(str) + " | note: " + str(r.note)
    return D


def current_straub(D: pd.DataFrame) -> pd.DataFrame:
    """the current cohort's Straub = straub_dayN.csv final value (Step 02 holds only hand scores)"""
    for f in straub_files():
        S = pd.read_csv(f, dtype={"cage": str})
        d = int(re.fullmatch(r"straub_day(\d+)", f.stem).group(1))
        for _, r in S.iterrows():
            v = pd.to_numeric(r.get("straub_final_pct"), errors="coerce")
            if pd.isna(v):
                v = pd.to_numeric(r.get("auto_straub_pct"), errors="coerce")
            mk = f"{str(r.cage)[:4]}_{r.mouse}"
            m = (D.mouse_key.astype(str) == mk) & (D.day == d)
            if m.any() and pd.notna(v):
                D.loc[m, "straub"] = float(v)
    return D


def straub_files() -> list[Path]:
    out = []
    for f in sorted(DAYS_DIR.glob("day*/straub_day*.csv")):
        mm = re.fullmatch(r"straub_day(\d+)", f.stem)       # skip backups like straub_day9_before_....csv
        if mm and f.parent.name == f"day{mm.group(1)}":
            out.append(f)
    return out


def pr_table_file(nd: int) -> Path | None:
    f = DAYS_DIR / f"pr_table_days1-{nd}.csv"
    if f.exists():
        return f
    c = [(int(re.fullmatch(r"pr_table_days1-(\d+)", p.stem).group(1)), p) for p in DAYS_DIR.glob("pr_table_days1-*.csv")
         if re.fullmatch(r"pr_table_days1-(\d+)", p.stem)]
    return max(c)[1] if c else None


def current_vs_pr_table(D: pd.DataFrame, nd: int) -> pd.DataFrame:
    """own PR (Step 02) must equal the routine's PR table; the Yoked days take the PR table's partner_pr (checked = the
    passive session's own delivered drops + 1, and compared with the Step 02 reward onsets + 1). Mismatches are printed."""
    f = pr_table_file(nd)
    D["partner"] = ""
    if f is None:
        print("  !! no PR table found: Yoked x = Step 02 reward onsets + 1")
        return D
    Tb = pd.read_csv(f, dtype={"cage": str})
    Tb["mouse_key"] = Tb.cage.str[:4] + "_" + Tb.mouse.astype(str)
    Tb = Tb.drop_duplicates(["mouse_key", "day"]).set_index(["mouse_key", "day"])
    n_own = bad_own = n_y = bad_dl = bad_s02 = 0
    for i in D.index:
        k = (str(D.at[i, "mouse_key"]), int(D.at[i, "day"]))
        if k not in Tb.index:
            continue
        t = Tb.loc[k]
        if D.at[i, "own"]:
            n_own += 1
            if pd.notna(t["pr"]) and abs(float(t["pr"]) - float(D.at[i, "pr"])) > 0.5:
                bad_own += 1
                print(f"  !! own PR {k}: Step 02 {D.at[i, 'pr']:.0f} vs PR table {t['pr']:.0f}")
        elif D.at[i, "yoked"] and pd.notna(t.get("partner_pr")):
            n_y += 1
            pp = float(t["partner_pr"])
            dl = pd.to_numeric(t.get("rewards_saleae"), errors="coerce")
            if pd.notna(dl) and abs(float(dl) + 1 - pp) > 0.5:
                bad_dl += 1
                print(f"  !! yoked {k}: partner's PR {pp:.0f} but drops delivered {dl:.0f} (+1 = {dl + 1:.0f})")
            s02 = D.at[i, "x"]
            if pd.isna(s02) or abs(float(s02) - pp) > 0.5:
                bad_s02 += 1
                print(f"  .. yoked {k}: partner's PR {pp:.0f}, Step 02 reward onsets + 1 = {s02}")
            D.at[i, "x"] = pp
            D.at[i, "x_source"] = "partner's PR (PR table; = drops delivered + 1)"
            if pd.notna(t.get("partner")):
                D.at[i, "partner"] = str(t["partner"]).replace("m ", "_").replace("f ", "_")
    print(f"  current cohort vs {f.name}: own PR = PR table in {n_own - bad_own}/{n_own} sessions; Yoked x: partner's PR "
          f"= delivered + 1 in {n_y - bad_dl}/{n_y}, = Step 02 reward onsets + 1 in {n_y - bad_s02}/{n_y}")
    return D


_PERM: dict[int, np.ndarray] = {}


def spearman(x, y):
    """Spearman rho; exact two-sided permutation p for n <= 8 (every pairing), scipy's t approximation above"""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n < 4 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan, np.nan, ""
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    r = float(np.corrcoef(rx, ry)[0, 1])
    if n <= 8:
        if n not in _PERM:
            _PERM[n] = np.array(list(permutations(range(n))))
        Y = ry[_PERM[n]]
        xc = rx - rx.mean()
        Yc = Y - ry.mean()
        rp = (Yc @ xc) / (np.sqrt((Yc ** 2).sum(axis=1)) * np.sqrt((xc ** 2).sum()))
        return r, float(np.mean(np.abs(rp) >= abs(r) - 1e-9)), "exact permutation"
    rr, pp = stats.spearmanr(x, y)
    return float(rr), float(pp), "t approximation"


def rm_corr(g: pd.DataFrame, x: str, y: str, by):
    """correlation of x and y after centring both on each group's mean (by = mouse: within-mouse; by = day, or
    [day, group] for the both-groups rows: same day); naive t p (sessions as independent)"""
    h = g.dropna(subset=[x, y])
    h = h[h.groupby(by)[x].transform("size") >= 2]
    k = h.groupby(by).ngroups if len(h) else 0
    n = len(h)
    df = n - k - 1
    if df < 2:
        return np.nan, np.nan, n, k
    xc = h[x] - h.groupby(by)[x].transform("mean")
    yc = h[y] - h.groupby(by)[y].transform("mean")
    if xc.std() == 0 or yc.std() == 0:
        return np.nan, np.nan, n, k
    r = float(np.corrcoef(xc, yc)[0, 1])
    t = r * np.sqrt(df / max(1e-12, 1 - r * r))
    return r, float(2 * stats.t.sf(abs(t), df)), n, k


_RNG = np.random.default_rng(20261123)
_PCACHE: dict = {}


def _relabelings(sizes: list[int]):
    """every relabeling that moves mice only inside their own group (groups are contiguous blocks), or N_RAND random"""
    total = int(np.prod([math.factorial(s) for s in sizes]))
    starts = np.cumsum([0] + sizes[:-1])
    if total <= EXACT_MAX:
        idx = np.zeros((1, 0), dtype=int)
        for s0, s in zip(starts, sizes):
            q = np.array(list(permutations(range(s0, s0 + s))), dtype=int)
            idx = np.concatenate([np.repeat(idx, len(q), axis=0), np.tile(q, (len(idx), 1))], axis=1)
        return idx, True
    blocks = [np.array([_RNG.permutation(np.arange(s0, s0 + s)) for _ in range(N_RAND)]) for s0, s in zip(starts, sizes)]
    return np.concatenate(blocks, axis=1), False


def _perm_r(XC: np.ndarray, YC: np.ndarray, idx: np.ndarray) -> np.ndarray:
    out = []
    yok = np.isfinite(YC)
    for a in range(0, len(idx), 4096):
        Xp = XC[idx[a:a + 4096]]
        M = np.isfinite(Xp) & yok[None]
        n = M.sum((1, 2)).astype(float)
        x = np.where(M, Xp, 0.0)
        y = np.where(M, YC[None], 0.0)
        sx, sy = x.sum((1, 2)), y.sum((1, 2))
        with np.errstate(invalid="ignore", divide="ignore"):
            cov = (x * y).sum((1, 2)) - sx * sy / n
            vx = (x * x).sum((1, 2)) - sx * sx / n
            vy = (y * y).sum((1, 2)) - sy * sy / n
            out.append(cov / np.sqrt(vx * vy))
    return np.concatenate(out)


def sameday_test(g: pd.DataFrame, key: str, cen: list):
    """same-day r (centred on cen) with a whole-mouse permutation p (mice swapped inside their group)"""
    h = g.dropna(subset=["x", key])
    h = h[h.groupby(cen)["x"].transform("size") >= 2].copy()
    n, k = len(h), (h.groupby(cen).ngroups if len(h) else 0)
    if n - k - 1 < 2 or h.mouse_key.nunique() < 3:
        return np.nan, np.nan, n, k
    h["xc"] = h.x - h.groupby(cen).x.transform("mean")
    h["yc"] = h[key] - h.groupby(cen)[key].transform("mean")
    if h.xc.std() == 0 or h.yc.std() == 0:
        return np.nan, np.nan, n, k
    r = float(np.corrcoef(h.xc, h.yc)[0, 1])
    mice = h.drop_duplicates("mouse_key").sort_values(["group", "mouse_key"])[["mouse_key", "group"]]
    order = {m: i for i, m in enumerate(mice.mouse_key)}
    days = sorted(h.day.unique())
    dpos = {d: j for j, d in enumerate(days)}
    XC = np.full((len(order), len(days)), np.nan)
    YC = np.full_like(XC, np.nan)
    for mk, d, xv, yv in zip(h.mouse_key, h.day, h.xc, h.yc):
        XC[order[mk], dpos[d]], YC[order[mk], dpos[d]] = xv, yv
    sizes = mice.groupby("group", sort=True).size().tolist() if "group" in cen else [len(order)]
    ck = (XC.tobytes(), YC.tobytes(), tuple(sizes))
    if ck not in _PCACHE:
        idx, exact = _relabelings(sizes)
        rp = _perm_r(XC, YC, idx)
        hit = np.abs(rp) >= abs(r) - 1e-9
        _PCACHE[ck] = float(np.mean(hit)) if exact else float((1 + hit.sum()) / (1 + len(rp)))
    return r, _PCACHE[ck], n, k


def fe2_r(g: pd.DataFrame, key: str) -> float:
    """r of x and X after removing mouse and day effects (two-way fixed effects); descriptive"""
    h = g.dropna(subset=["x", key])
    if len(h) < 6 or h.mouse_key.nunique() < 2 or h.day.nunique() < 2:
        return np.nan
    Z = np.column_stack([pd.get_dummies(h.mouse_key).to_numpy(float), pd.get_dummies(h.day).to_numpy(float)[:, 1:]])
    res = []
    for v in (h.x.to_numpy(float), h[key].to_numpy(float)):
        b = np.linalg.lstsq(Z, v, rcond=None)[0]
        res.append(v - Z @ b)
    if res[0].std() < 1e-9 or res[1].std() < 1e-9:
        return np.nan
    return float(np.corrcoef(res[0], res[1])[0, 1])


def subset(D: pd.DataFrame, reward: str, who: str, phase=None) -> pd.DataFrame:
    m = D.reward.eq(reward)
    if phase:
        m &= D.phase.eq(phase)
    if who == "yoked":
        return D[m & D.yoked]
    m &= D.own
    if who in ("Active", "Passive"):
        m &= D.group.eq(who)
    return D[m]


def mouse_test(g: pd.DataFrame, key: str, cen: list):
    g2 = g[g.groupby(cen)[key].transform("size") >= 2].copy()
    g2["x_d"] = g2.x - g2.groupby(cen).x.transform("mean")
    g2["v_d"] = g2[key] - g2.groupby(cen)[key].transform("mean")
    M = g2.groupby("mouse_key").agg(x=("x_d", "mean"), v=("v_d", "mean"))
    r1, p1, how = spearman(M.x, M.v) if len(M) else (np.nan, np.nan, "")
    return r1, p1, how, len(M)


def stat_row(g: pd.DataFrame, key: str, cen: list) -> dict:
    g = g.dropna(subset=["x", key])
    one_each = len(g) > 0 and int(g.groupby("mouse_key").size().max()) == 1
    partners = g.partner[g.yoked & g.partner.astype(str).ne("")] if "partner" in g else pd.Series(dtype=str)
    out = dict(n_sessions=len(g), n_mice_any=g.mouse_key.nunique(), centred_on="+".join(cen),
               n_partners=partners.nunique() if len(partners) else np.nan,
               groups_in="+".join(sorted(g.group.unique())), days=",".join(map(str, sorted(g.day.unique()))))
    # mouse level: each session's x and X centred on THAT DAY's mean over the mice of this session group measured that
    # day (days with >= 2 mice; both-groups rows: that day's mean of the mouse's own group, so a difference between
    # Active and Passive cannot make the correlation), then averaged per mouse = the mouse's usual standing vs the
    # others, so the test does not depend on which days a mouse had
    r1, p1, how, nm = mouse_test(g, key, cen)
    if one_each:      # one session per mouse: the same-day r would use the very same points as the mouse level
        r3, p3, n3, k3, note = np.nan, np.nan, len(g), g.groupby(cen).ngroups, "= mouse level (1 session per mouse)"
    else:
        r3, p3, n3, k3 = sameday_test(g, key, cen)
        note = ""
    rc = rm_corr(g, "x", key, cen + ["cage"])[0] if not one_each else np.nan
    r2, p2, n2, k2 = rm_corr(g, "x", key, "mouse_key")
    r0, p0, _ = spearman(g.x, g[key]) if len(g) else (np.nan, np.nan, "")
    out.update(n_mice=nm, mouse_rho=r1, mouse_p=p1, mouse_p_method=how,
               mouse_p_min_possible=(2 / math.factorial(nm) if 4 <= nm <= 8 else np.nan),
               sameday_r=r3, sameday_p=p3, sameday_n=n3, sameday_days=k3, sameday_note=note,
               sameday_p_method="whole-mouse permutation" if np.isfinite(p3) else "",
               sameday_r_cage=rc, within_r_naive=r2, within_p_naive=p2, within_n=n2, within_mice=k2,
               within_dayadj_r=fe2_r(g, key), pooled_rho_naive=r0, pooled_p_naive=p0)
    return out


def loo(g: pd.DataFrame, key: str, cen: list, r: dict) -> dict:
    """leave one mouse / one day out for the tests with p < 0.05; reports the worst p and who"""
    g = g.dropna(subset=["x", key])
    out = {}
    tests = [("mouse", "mouse_p"), ("sameday", "sameday_p")]
    for name, pk in tests:
        if not (np.isfinite(r[pk]) and r[pk] < 0.05):
            continue
        worst = (0.0, "")
        for what, col in (("mouse", "mouse_key"), ("day", "day")):
            if g[col].nunique() < 2:            # one day only: leaving it out leaves nothing to test
                continue
            for v in sorted(g[col].unique()):
                h = g[g[col] != v]
                p = mouse_test(h, key, cen)[1] if name == "mouse" else sameday_test(h, key, cen)[1]
                if np.isfinite(p) and p > worst[0]:      # a drop that leaves too few to test is not counted
                    worst = (p, f"without {what} {v}")
        out[f"loo_{name}_max_p"], out[f"loo_{name}_worst"] = worst
    return out


def relate(D: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for coh, C in D.groupby("cohort"):
        for set_name, excl in day_sets():
            Cs = C[~C.day.isin(excl)]
            for gname, rew, who, ph in SETS:
                g = subset(Cs, rew, who, ph)
                cen = ["day", "group"] if who == "both" else ["day"]
                for key, label in MEAS:
                    st = stat_row(g, key, cen)
                    if set_name == "all" and gname in MAIN:
                        st.update(loo(g, key, cen, st))
                    rows.append(dict(cohort=coh, day_set=set_name, excluded_days=";".join(map(str, sorted(excl))),
                                     session_group=gname, reward=rew, phase=ph or "", measure=key, label=label, **st))
    return pd.DataFrame(rows)


def pf(p) -> str:
    if not np.isfinite(p):
        return "-"
    return "<0.001" if p < 0.001 else (f"{p:.3f}" if p < 0.01 else f"{p:.2f}")


def pf3(p) -> str:
    if not np.isfinite(p):
        return ""
    return "<.001" if p < 0.001 else f"{p:.3f}".lstrip("0")


def shown(R: pd.DataFrame) -> pd.Series:
    """rows worth showing: the main groups; 'Both groups' only when both groups are in it (otherwise it repeats one)"""
    return R.session_group.isin(MAIN) & (~R.session_group.str.startswith("Both") | R.groups_in.eq("Active+Passive"))


def _items(R: pd.DataFrame, coh: str) -> list:
    """one item per session group x measure with a p < 0.05 test (day set 'all'), the lab's order: per-group mouse level,
    then both-groups mouse level, then same-day only"""
    A = R[(R.cohort == coh) & (R.day_set == "all") & shown(R)]
    out = []
    for _, r in A.iterrows():
        hm = np.isfinite(r.mouse_p) and r.mouse_p < 0.05
        hs = np.isfinite(r.sameday_p) and r.sameday_p < 0.05
        if not (hm or hs):
            continue
        rank = (0 if hm and not r.session_group.startswith("Both") else 1 if hm else 2, min(
            r.mouse_p if hm else 1, r.sameday_p if hs else 1))
        out.append((rank, r))
    return [r for _, r in sorted(out, key=lambda t: t[0])]


def _phase_note(R: pd.DataFrame, r) -> str:
    if r.reward != NAT:
        return ""
    bits = []
    for ph in ("Baseline", "Withdrawal"):
        q = R[(R.cohort == r.cohort) & (R.day_set == "all") & (R.session_group == f"{r.session_group} · {ph}")
              & (R.measure == r.measure)]
        if len(q) and q.n_sessions.iloc[0] > 0 and np.isfinite(q.mouse_rho.iloc[0]):
            bits.append(f"{ph} {q.mouse_rho.iloc[0]:+.2f}")
    return f" [{', '.join(bits)}]" if len(bits) == 2 else ""


def found_text(R: pd.DataFrame, coh: str, max_items: int = 3) -> str:
    its = _items(R, coh)
    if not its:
        extra = " (4 mice per group: the smallest possible mouse-level p is 0.083)" if coh == "current" else ""
        return "no relation with p < 0.05" + extra
    txt, used = [], set()
    for r in its:
        k = (r.measure, r.session_group)
        if k in used:
            continue
        parts = []
        if np.isfinite(r.mouse_rho):
            parts.append(f"ρ {r.mouse_rho:+.2f} (p {pf3(r.mouse_p)}, {int(r.n_mice)} mice)")  # mouse level
        if np.isfinite(r.sameday_r):
            parts.append(f"day r {r.sameday_r:+.2f} (p {pf3(r.sameday_p)})")
        s = f"{SHORT[r.measure]} {SET_TINY[r.session_group]} " + " / ".join(parts)
        # dose link: on morphine days the Active and its yoked partner share x by design
        if r.session_group in ("Active · morphine", "Passive · morphine yoked (delivered + 1)"):
            other = ("Passive · morphine yoked (delivered + 1)" if r.session_group == "Active · morphine"
                     else "Active · morphine")
            if any((q.measure, q.session_group) == (r.measure, other) for q in its):
                s += " + same in yoked pair (shared x = the same drop count)"
                used.add((r.measure, other))
        bad = [w for t in ("mouse", "sameday") for p, w in [(r.get(f"loo_{t}_max_p", np.nan),
                                                            r.get(f"loo_{t}_worst", ""))]
               if isinstance(p, float) and np.isfinite(p) and p >= 0.05]
        if bad:
            s += f" (rests on {bad[0].replace('without ', '')})"
        if coh == "current" and np.isfinite(r.sameday_r_cage) and np.isfinite(r.sameday_r):
            s += f" (in-cage r {r.sameday_r_cage:+.2f})"
        s += _phase_note(R, r)
        txt.append(s)
        used.add(k)
    more = len(txt) - max_items
    return "p < 0.05: " + "; ".join(txt[:max_items]) + (f"; +{more} more in the csv" if more > 0 else "")


def robust_text(R: pd.DataFrame, coh: str, max_items: int = 4) -> str:
    """how many p < 0.05 results of day set 'all' stay p < 0.05 in every day set - counted only where some set really
    leaves sessions out (assay / Straub days are in every set). Every set's row counts, also where a 'Both groups' row
    holds one group only or a group has no session in that set (review 2026-11-24: shown() used to drop those sets,
    so 'TI both·mor', which 6 sets change from 28 to 16 Active-only sessions, was called 'not touched')."""
    A = R[(R.cohort == coh) & R.session_group.isin(MAIN)]
    stay, fall, fixed = [], [], []
    for r in _items(R, coh):
        o = A[(A.session_group == r.session_group) & (A.measure == r.measure)]
        lab = f"{SHORT[r.measure]} {SET_TINY[r.session_group]}"
        if (o.days.nunique(dropna=False) == 1 and o.n_sessions.nunique() == 1
                and o.groups_in.nunique(dropna=False) == 1):
            fixed.append(lab)
            continue
        pk = "mouse_p" if np.isfinite(r.mouse_p) and r.mouse_p < 0.05 else "sameday_p"
        ok = o[pk] < 0.05
        if r.session_group.startswith("Both"):           # a set that leaves one group out does not re-test it (11-24)
            ok = ok & o.groups_in.eq("Active+Passive")
        (stay if bool(ok.all()) else fall).append(lab)
    if not (stay or fall or fixed):
        return "no p < 0.05 result in 'all'"
    t = f"{len(stay)} of {len(stay) + len(fall)} testable stay p < 0.05 in all 8 sets"
    if stay:
        t += f" ({', '.join(stay[:max_items])}{', ...' if len(stay) > max_items else ''})"
    if fixed:
        t += f"; {len(fixed)} not touched by any set (assay / Straub days are in every set)"
    return t


STYLE = {"Active": dict(m="o", c=GCOL["Active"], fill=True, z=4),
         "Passive": dict(m="o", c=GCOL["Passive"], fill=True, z=3),
         "Passive yoked": dict(m="s", c=GCOL["Passive"], fill=False, z=3)}
ROW_TITLE = {NAT: "NATURAL REWARD days  (Baseline 3-5 · Withdrawal 14-16)",
             MOR: "MORPHINE days  (Yoked 6-10 · Equal Access 11-13 · Re-exposure 17-18)  ·  PR = drops earned + 1 (drop count = dose proxy; dose not measured)"}
PANEL_SETS = {NAT: [("Active · natural reward", "Active", "Active"), ("Passive · natural reward", "Passive", "Passive"),
                    ("Both groups · natural reward", "Both", None)],
              MOR: [("Active · morphine", "Active", "Active"), ("Passive · morphine (own PR)", "P own", "Passive"),
                    ("Both groups · morphine (own PR)", "Both", None),
                    ("Passive · morphine yoked (delivered + 1)", "P yoked", "Passive yoked")]}


def cell(v, p, n, note="") -> str:
    if note:
        return "= mouse level"
    if not np.isfinite(v):
        return f"   -   ({int(n)})" if n else "   -"
    return f"{v:+.2f} p{pf3(p)}{'*' if p < 0.05 else ' '}({int(n)})"


def panel_header(ax, RA, rew, key):
    """a small table above the panel: one line per session group, mouse level and same day; a coloured mark beside
    the group name carries the identity (text stays dark)"""
    rows = []
    for gname, tag, cat in PANEL_SETS[rew]:
        r = RA[(RA.session_group == gname) & (RA.measure == key)]
        if len(r) and int(r.n_sessions.iloc[0]) > 0:
            r = r.iloc[0]
            rows.append((tag, cat, cell(r.mouse_rho, r.mouse_p, r.n_mice),
                         cell(r.sameday_r, r.sameday_p, r.sameday_n, r.sameday_note if isinstance(r.sameday_note, str)
                              else "")))
    if not rows:
        return
    kw = dict(xycoords="axes fraction", textcoords="offset points", fontsize=7.8, family="DejaVu Sans Mono",
              va="bottom", ha="left", annotation_clip=False)
    lh = 10.2
    y0 = 4 + lh * len(rows)
    ax.annotate(f"  {'':<7} {'mouse ρ (mice)':<17} same day (n)", xy=(0, 1), xytext=(0, y0), color="#7b8794", **kw)
    for k, (tag, cat, a, b) in enumerate(rows):
        y = y0 - lh * (k + 1)
        if cat:
            st = STYLE[cat]
            ax.annotate("●" if st["fill"] else "□", xy=(0, 1), xytext=(0, y), color=st["c"], **kw)
        ax.annotate(f"  {tag:<7} {a:<17} {b}", xy=(0, 1), xytext=(0, y), color="#1f2933", **kw)


def cohort_figure(D: pd.DataFrame, R: pd.DataFrame, coh: str, out_png: Path, title: str, note: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    A = D[D.cohort == coh]
    RA = R[(R.cohort == coh) & (R.day_set == "all") & shown(R)]
    sex = CAGE_SEX.get(coh, {})
    fig, axes = plt.subplots(2, len(MEAS), figsize=(3.95 * len(MEAS), 11.2), squeeze=False)
    fig.subplots_adjust(left=0.045, right=0.975, top=0.75, bottom=0.06, wspace=0.27, hspace=0.84)
    for i, rew in enumerate((NAT, MOR)):
        for j, (key, label) in enumerate(MEAS):
            ax = axes[i, j]
            n_pts = 0
            for cat in ("Passive yoked", "Passive", "Active"):            # Active last = on top
                h = A[(A.cat == cat) & (A.reward == rew)].dropna(subset=["x", key])
                st = STYLE[cat]
                face = st["c"] if st["fill"] else "none"
                for mk, hm in h.groupby("mouse_key"):
                    ax.plot(hm.x, hm[key], st["m"], ls="none", ms=3.8, color=st["c"], mfc=face, alpha=0.5, zorder=2)
                    ax.plot([hm.x.mean()], [hm[key].mean()], st["m"], ls="none", ms=9, color=st["c"], mfc=face,
                            mec=st["c"] if not st["fill"] else "white", mew=1.6, zorder=st["z"])
                    s = sex.get(str(mk)[:4])
                    if s:
                        ax.annotate(s, (hm.x.mean(), hm[key].mean()), xytext=(5, 3), textcoords="offset points",
                                    fontsize=6.5, color="#52606d", zorder=5)
                n_pts += len(h)
            panel_header(ax, RA, rew, key)
            ax.set_title(label, loc="left", fontsize=10.5, color="#1f3a5f", fontweight="bold", pad=64)
            ax.set_xlabel("PR (last requirement)" if rew == NAT else "PR = drops earned + 1  ·  □ drops delivered + 1", fontsize=8.5)
            ax.set_ylabel(label, fontsize=8.5)
            ax.tick_params(labelsize=8)
            ax.grid(color="#e4e7eb", lw=0.6)
            ax.spines[["top", "right"]].set_visible(False)
            if n_pts == 0:
                ax.text(0.5, 0.5, "no sessions yet", transform=ax.transAxes, ha="center", va="center", color="#7b8794",
                        fontsize=9)
        top = axes[i, 0].get_position().y1
        fig.text(0.006, top + (0.112 if i == 0 else 0.108), ROW_TITLE[rew], fontsize=11.5, fontweight="bold",
                 color="#1f3a5f", ha="left", va="bottom")
        fig.add_artist(Line2D([0.006, 0.975], [top + (0.108 if i == 0 else 0.104)] * 2, color="#c9d1da", lw=1.0))
    hd = [Line2D([], [], marker="o", ls="none", color=GCOL["Active"], ms=7, label="Active mouse · own PR"),
          Line2D([], [], marker="o", ls="none", color=GCOL["Passive"], ms=7,
                 label="Passive mouse · own PR (days 3-5, from day 11)"),
          Line2D([], [], marker="s", ls="none", color=GCOL["Passive"], mfc="none", mew=1.6, ms=7,
                 label="Passive mouse · Yoked days 6-10: x = drops delivered + 1 (= PR of the replayed partner session; "
                       "a pair shares x)")]
    fig.legend(handles=hd, loc="upper left", ncol=3, fontsize=9, frameon=False, bbox_to_anchor=(0.0, 0.94))
    fig.suptitle(title, fontsize=14, color="#1f3a5f", x=0.006, ha="left", y=0.996, fontweight="bold")
    fig.text(0.006, 0.966, "small = session · big = mouse mean · mouse ρ = Spearman of day-centred mouse means (exact p "
             "for <= 8 mice; 4 mice: smallest possible p .083; 3 mice: no test) · same day = sessions centred on each "
             "day, p from swapping whole mice (= mouse level when 1 session per mouse) · * p < 0.05 · Both = each mouse "
             "vs its own group that day · exploratory, no correction for multiple tests\n" + note, fontsize=8.6,
             color="#52606d", ha="left", va="top", linespacing=1.35)
    fig.savefig(out_png, dpi=100)
    plt.close(fig)


def dayset_figure(R: pd.DataFrame, out_png: Path, title: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HM = [("Active · natural reward", "A·nat"), ("Passive · natural reward", "P·nat"),
          ("Both groups · natural reward", "both·nat"), ("Active · morphine", "A·mor"),
          ("Passive · morphine (own PR)", "P own·mor"), ("Both groups · morphine (own PR)", "both·mor"),
          ("Passive · morphine yoked (delivered + 1)", "P yoked")]          # both·mor added 2026-11-24 (review)
    sets = day_sets()
    cohorts = [c for c in ("current", "previous") if c in set(R.cohort)]
    st_ = [("mouse_rho", "mouse_p", "mouse ρ"), ("sameday_r", "sameday_p", "same-day r")]
    fig, axes = plt.subplots(1, len(st_) * len(cohorts), figsize=(6.4 * len(st_) * len(cohorts), 12.6),   # same height with 7 rows (wider on the slide)
                             squeeze=False)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#eceff2")
    short_set = {"all": "all", "no_first_day": "no 1st day", "no_first_day_17": "no 1st day +17",
                 "prev_phase_medians": "prev phase medians", "prev_dashboard": "prev dashboard",
                 "prev_pupil": "prev pupil / TI", "prev_step10": "prev step10", "no_phase_start_17": "no phase start +17"}
    im = None
    for j, coh in enumerate(cohorts):
        for i, (vk, pk, stitle) in enumerate(st_):
            ax = axes[0, j * len(st_) + i]
            rows = [(key, g, f"{SHORT[key]} | {gl}") for key, _ in MEAS for g, gl in HM]
            V = np.full((len(rows), len(sets)), np.nan)
            Pm = np.full_like(V, np.nan)
            Nn = np.full_like(V, np.nan)
            one = np.zeros(V.shape, dtype=bool)      # a 'Both groups' cell that holds one group only (review 11-24)
            same = []
            for a, (key, g, _) in enumerate(rows):
                o = R[(R.cohort == coh) & (R.session_group == g) & (R.measure == key)]
                same.append(len(o) > 0 and o.days.nunique() == 1 and o.n_sessions.nunique() == 1)
                for b, (sn, _) in enumerate(sets):
                    r = o[o.day_set == sn]
                    if len(r):
                        if (g.startswith("Both") and str(r.groups_in.iloc[0]) != "Active+Passive"
                                and pd.to_numeric(r.n_sessions.iloc[0], errors="coerce") > 0):
                            one[a, b] = True             # the set leaves one group out: the cell would repeat it
                            continue
                        V[a, b], Pm[a, b] = r[vk].iloc[0], r[pk].iloc[0]
                        Nn[a, b] = r["n_mice" if vk == "mouse_rho" else "sameday_n"].iloc[0]
            keep = ~np.all(np.isnan(V), axis=1)
            V, Pm, Nn, one = V[keep], Pm[keep], Nn[keep], one[keep]
            labels = [lab for (k_, g_, lab), kk in zip(rows, keep) if kk]
            same = [s for s, kk in zip(same, keep) if kk]
            im = ax.imshow(np.ma.masked_invalid(V), cmap=cmap, vmin=-1, vmax=1, aspect="auto")
            for a in range(V.shape[0]):
                for b in range(V.shape[1]):
                    if np.isfinite(V[a, b]):
                        star = "*" if np.isfinite(Pm[a, b]) and Pm[a, b] < 0.05 else ""
                        ax.text(b, a, f"{V[a, b]:+.2f}{star}", ha="center", va="center", fontsize=8.0,
                                color="white" if abs(V[a, b]) > 0.6 else "#1f2933", fontweight="bold" if star else None)
                    elif one[a, b]:
                        ax.text(b, a, "1 group", ha="center", va="center", fontsize=6.6, color="#9aa5b1")
                    else:
                        ax.text(b, a, "-", ha="center", va="center", fontsize=8, color="#9aa5b1")
            ax.set_yticks(range(len(labels)))
            ax.set_yticklabels([f"{lab}  (n {int(np.nanmax(Nn[a])) if np.isfinite(Nn[a]).any() else 0})"
                                + ("  =" if same[a] else "") for a, lab in enumerate(labels)], fontsize=8.4)
            for a, t in enumerate(ax.get_yticklabels()):
                if same[a]:
                    t.set_color("#9aa5b1")
            ax.set_xticks(range(len(sets)))
            ax.set_xticklabels([short_set.get(sn, sn) for sn, _ in sets], fontsize=8.2, rotation=55, ha="left",
                               rotation_mode="anchor")
            ax.xaxis.tick_top()
            for a in range(1, len(labels)):
                if labels[a].split(" | ")[0] != labels[a - 1].split(" | ")[0]:
                    ax.axhline(a - 0.5, color="white", lw=3)
            ax.set_title(f"{'CURRENT' if coh == 'current' else 'PREVIOUS'} cohort · {stitle}", loc="left",
                         fontsize=12, color="#1f3a5f", fontweight="bold", pad=82)
    fig.suptitle(title, fontsize=14, color="#1f3a5f", x=0.005, ha="left", fontweight="bold")
    fig.text(0.005, 0.955, "A = Active, P = Passive, nat = natural reward days, mor = morphine days, P yoked = Passive on "
             "the Yoked days (x = drops delivered + 1) · cell = correlation (* p < 0.05, exploratory; n = mice for ρ, "
             "sessions for same day, in 'all') · grey label with '=' = no day set changes its days (assay / Straub days)\n"
             "columns = the lab's 8 day sets (tools\\day_sets.csv: all; without days 4,6,11,14 (+17); the previous "
             "cohort's phase medians -4,6,11,14,17, dashboard -3,6,11,14, pupil / TI -4,6,11,14, step10 -3; "
             "-3,6,11,14,17) · grey cell = too few · '1 group' = the set leaves one group out of a 'both' row · red = "
             "higher PR with higher value", fontsize=9.5, color="#52606d",
             ha="left", va="top")
    fig.tight_layout(rect=(0, 0, 0.985, 0.935), w_pad=1.2)
    if im is not None:
        cax = fig.add_axes([0.988, 0.25, 0.006, 0.4])
        fig.colorbar(im, cax=cax).ax.tick_params(labelsize=7)
    fig.savefig(out_png, dpi=100)
    plt.close(fig)


def _md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def input_signature(nd: int) -> dict:
    """what the outputs depend on: the Step 02 run (path, size, time), the PR table, the Straub tables, the previous
    cohort's tables and this code (content hashes); the deck re-runs the tool only when this changes"""
    sig = {}
    cur = newest_current()
    if cur is not None and cur.exists():
        s = cur.stat()
        sig["step02"] = [str(cur), s.st_size, round(s.st_mtime, 1)]
    pt = pr_table_file(nd)
    if pt is not None:
        sig["pr_table"] = [pt.name, _md5(pt)]
    sig["straub"] = {f.name: _md5(f) for f in straub_files()}
    for k, f in (("prev_cache", PREV_CACHE), ("prev_pairs", PREV_PAIRS), ("prev_notes", PREV_NOTES),
                 ("day_sets", T / "day_sets.csv"), ("code", Path(__file__).resolve())):
        sig[k] = _md5(f) if f.exists() else ""
    return sig


def signature_file(nd: int) -> Path:
    return DAYS_DIR / f"pr_relationships_inputs_days1-{nd}.json"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-prev", type=Path, help="previous cohort ALL_mice_longitudinal.csv (read once, cached)")
    ap.add_argument("--current", type=Path, help="current cohort ALL_mice_longitudinal.csv (default: newest Step 02 run)")
    a = ap.parse_args(argv)
    if a.build_prev:
        PREV_DIR.mkdir(parents=True, exist_ok=True)
        D = mouse_day(a.build_prev)
        if PREV_CACHE.exists():                            # never overwrite: the older cache moves aside
            sup = PREV_DIR / f"_superseded_{datetime.now():%Y%m%d_%H%M%S}"
            k = 1
            while (sup / PREV_CACHE.name).exists():
                sup = PREV_DIR / f"_superseded_{datetime.now():%Y%m%d_%H%M%S}_{k}"
                k += 1
            sup.mkdir(parents=True, exist_ok=True)
            shutil.move(str(PREV_CACHE), str(sup / PREV_CACHE.name))
        D.to_csv(PREV_CACHE, index=False)
        print(PREV_CACHE, len(D), "mouse-days")
        return 0
    cur = a.current or newest_current()
    parts, nd = [], 0
    if cur is not None and cur.exists():
        C = prepare(mouse_day(cur, P / "cohort_manifest.csv"), "current")
        nd = int(C.day.max())
        C = current_straub(C)
        C = current_vs_pr_table(C, nd)
        parts.append(C)
        print(f"  current cohort: {cur}")
    if PREV_CACHE.exists():
        Pv = pd.read_csv(PREV_CACHE)
        if "rew_n" not in Pv:
            print("  !! the previous-cohort cache has no reward counts (rew_n): Yoked days left out; rebuild it with "
                  "--build-prev")
        parts.append(previous_notes(prepare(Pv, "previous")))
    if not parts:
        print("no data")
        return 1
    D = pd.concat(parts, ignore_index=True)
    D["partner"] = D["partner"].fillna("") if "partner" in D else ""
    R = relate(D)
    tag = f"days1-{nd}"
    DAYS_DIR.mkdir(parents=True, exist_ok=True)
    keep = ["cohort", "mouse_key", "cage", "group", "partner", "day", "phase", "reward", "cat", "own", "yoked", "x",
            "x_source", "pr", "rew_n", "ti", "tst", "hot", "straub", "pupil_rel", "pupil_resp_rel"]
    D[[c for c in keep if c in D.columns]].to_csv(DAYS_DIR / f"pr_relationships_data_{tag}.csv", index=False,
                                                 encoding="utf-8-sig")
    R.to_csv(DAYS_DIR / f"pr_relationships_stats_{tag}.csv", index=False, encoding="utf-8-sig")
    outs = []
    if (D.cohort == "current").any():
        o = DAYS_DIR / f"pr_relationships_current_{tag}.png"
        cohort_figure(D, R, "current", o, f"CURRENT cohort (days 1-{nd}) · PR vs tail immersion, TST, hot plate, Straub "
                      "tail and pupil · Active / Passive and natural reward / morphine apart",
                      "4 mice per group: 2 males (9669m, M) + 2 females (9670f, F) - cage = sex · pupil: the whole-eye "
                      "mice 9669m black (Active) and 9670f black (Passive) are left out, so 3 per group")
        outs.append(o)
    if (D.cohort == "previous").any():
        o = DAYS_DIR / "pr_relationships_previous_cohort.png"
        cohort_figure(D, R, "previous", o, "PREVIOUS cohort (Dec 2025; 14 mice: 6 Active, 8 Passive; days 3-18) · PR vs "
                      "tail immersion, TST, hot plate, Straub tail and pupil · Active / Passive and natural reward / "
                      "morphine apart",
                      "8 Passive mice replayed 6 partners (0910_black and 6100_black twice) · 6100 orange / red day 10 "
                      "replayed the partner's day-9 log; 4 Yoked days had an aborted replay with extra drops "
                      "(prev_yoked_x_notes.csv) · day 8 ran 01-07 or 01-12 by cage · Straub = 100 - Non_moving")
        outs.append(o)
    o = DAYS_DIR / f"pr_relationships_daysets_{tag}.png"
    dayset_figure(R, o, f"PR relationships in every day set (current days 1-{nd} · previous days 3-18)")
    outs.append(o)
    signature_file(nd).write_text(json.dumps(input_signature(nd), indent=1), encoding="utf-8")
    for o in outs:
        print(o)
    pd.set_option("display.width", 250)
    A = R[(R.day_set == "all") & shown(R)]
    print(A[["cohort", "session_group", "measure", "n_sessions", "n_mice", "mouse_rho", "mouse_p", "sameday_r",
             "sameday_p", "sameday_n"]].round(3).to_string(index=False))
    for coh in ("current", "previous"):
        if (R.cohort == coh).any():
            print(f"FOUND {coh}: {found_text(R, coh)}")
            print(f"DAY SETS {coh}: {robust_text(R, coh)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
