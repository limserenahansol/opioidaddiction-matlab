"""Combine Step 02 ALL_mice_longitudinal.csv files of several cohorts into one file with the same columns.

    python harmonize_longitudinal.py --in prevcohort=<run_009>\\ALL_mice_longitudinal.csv
                                     --in nov2026=<run_###>\\ALL_mice_longitudinal.csv
                                     --out <combined>\\ALL_cohorts_longitudinal.csv [--keep-duplicates]

What it does (nothing else is changed):
  * columns: the union of all inputs, in the order of the first file (Step 02 only writes the TST_/HOT_/STRAUB_
    columns when such files existed, so a cohort without assay files yet has fewer columns); missing ones are empty.
  * adds 'cohort' (the name before '=') as the first column. mouse_key stays <4-digit cage>_<mouse folder>; two cohorts
    with the same mouse_key would be ambiguous, so the script stops if that happens.
  * duplicate sessions: Step 02 counts every concat_out folder as a session, so the same recording can appear twice
    (previous cohort: 6099_white day 3 x3, 0911_orange day 7 x2). Sessions of the same mouse-day with identical
    frames, licks, rewards, camera time and pupil sum are duplicates; all but the first are dropped (listed in
    <out>_duplicates.csv). --keep-duplicates keeps them.
  * writes <out>_sessions.csv: one row per kept session (cohort, mouse_key, day, session, frames, licks, rewards,
    minutes, RequirementLast, isPassive, paradigm), the same numbers the old pipeline uses (0->1 changes of
    Lick_TTL / Injector_TTL, NaN = 0).
Reads in chunks, so multi-GB files are fine.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

KEY = ["mouse_key", "day_index", "session_idx"]


def session_stats(path: Path) -> pd.DataFrame:
    acc, last = {}, {}
    for ch in pd.read_csv(path, chunksize=1_000_000, low_memory=False,
                          usecols=lambda c: c in KEY + ["Frame", "Lick_TTL", "Injector_TTL", "CamTime_rel_s", "Diameter_px",
                                                        "RequirementLast", "isPassive", "Session_Paradigm"]):
        k = ch.mouse_key.astype(str) + "|" + ch.day_index.astype(str) + "|" + ch.session_idx.astype(str)
        for col in ("Lick_TTL", "Injector_TTL"):
            v = (pd.to_numeric(ch[col], errors="coerce").fillna(0) > 0.5).astype(np.int8)
            prev = v.groupby(k, sort=False).shift(1)
            first = prev.isna()
            prev[first] = [last.get((kk, col), 0) for kk in k[first]]
            ch[col + "_on"] = ((v == 1) & (prev == 0)).astype(int)
            ch[col + "_v"] = v
        ch["_k"] = k
        g = ch.groupby("_k", sort=False)
        a = g.agg(frames=("Frame", "size"), licks=("Lick_TTL_on", "sum"), rewards=("Injector_TTL_on", "sum"),
                  t0=("CamTime_rel_s", "min"), t1=("CamTime_rel_s", "max"), pupil_sum=("Diameter_px", "sum"),
                  req_last=("RequirementLast", "max"), passive=("isPassive", "max"), paradigm=("Session_Paradigm", "first"))
        lv = g[["Lick_TTL_v", "Injector_TTL_v"]].last()
        for kk, r in a.iterrows():
            if kk in acc:
                b = acc[kk]
                for c in ("frames", "licks", "rewards", "pupil_sum"):
                    b[c] += r[c]
                b["t0"], b["t1"] = min(b["t0"], r["t0"]), max(b["t1"], r["t1"])
            else:
                acc[kk] = dict(r)
            last[(kk, "Lick_TTL")] = int(lv.at[kk, "Lick_TTL_v"]); last[(kk, "Injector_TTL")] = int(lv.at[kk, "Injector_TTL_v"])
    S = pd.DataFrame.from_dict(acc, orient="index")
    S[KEY] = S.index.to_series().str.split("|", expand=True).values
    S["day_index"] = S.day_index.astype(float).astype(int); S["session_idx"] = S.session_idx.astype(float).astype(int)
    S["minutes"] = (S.t1 - S.t0) / 60
    return S.reset_index(drop=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inputs", action="append", required=True, help="name=path (repeat)")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--keep-duplicates", action="store_true")
    a = ap.parse_args(argv)
    ins = []
    for s in a.inputs:
        name, _, p = s.partition("=")
        if not p:
            raise SystemExit(f"--in needs name=path, got {s}")
        ins.append((name, Path(p)))
    cols = []
    for _, p in ins:
        for c in pd.read_csv(p, nrows=0).columns:
            if c not in cols:
                cols.append(c)
    stats, dups = [], []
    for name, p in ins:
        S = session_stats(p)
        S.insert(0, "cohort", name)
        sig = ["mouse_key", "day_index", "frames", "licks", "rewards", "t0", "t1", "pupil_sum"]
        S["duplicate"] = S.sort_values("session_idx").duplicated(subset=sig, keep="first")
        stats.append(S)
        dups.append(S[S.duplicate])
        print(f"{name}: {len(S)} sessions, {S.mouse_key.nunique()} mice, days {S.day_index.min()}-{S.day_index.max()}, "
              f"{int(S.duplicate.sum())} duplicate session(s)")
    allS = pd.concat(stats, ignore_index=True)
    clash = allS.groupby("mouse_key").cohort.nunique()
    if (clash > 1).any():
        raise SystemExit(f"the same mouse_key is in two cohorts: {list(clash[clash > 1].index)}")
    drop = set() if a.keep_duplicates else set(map(tuple, allS[allS.duplicate][["cohort"] + KEY].astype(str).values))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    first = True
    n_rows = 0
    for name, p in ins:
        for ch in pd.read_csv(p, chunksize=1_000_000, low_memory=False):
            if drop:
                kk = list(zip([name] * len(ch), ch.mouse_key.astype(str), ch.day_index.astype(str), ch.session_idx.astype(str)))
                ch = ch[[k not in drop for k in kk]]
            ch = ch.reindex(columns=cols)
            ch.insert(0, "cohort", name)
            ch.to_csv(a.out, mode="w" if first else "a", header=first, index=False)
            first = False
            n_rows += len(ch)
    keep = allS[~allS.duplicate] if not a.keep_duplicates else allS
    keep.drop(columns=["pupil_sum", "t0", "t1"]).to_csv(a.out.with_name(a.out.stem + "_sessions.csv"), index=False)
    pd.concat(dups).to_csv(a.out.with_name(a.out.stem + "_duplicates.csv"), index=False)
    print(f"wrote {a.out} ({n_rows:,} rows, {len(cols) + 1} columns), {len(keep)} sessions kept")
    return 0


if __name__ == "__main__":
    sys.exit(main())
