#!/usr/bin/env python3
"""
pupil_check_examples.py - one picture per session of the pupil the tables use, with the frames that are "out"
(added 2026-11-18; user: "show me png files as example of 8 mice x 5 days = 40 examples of pupil detection you did and
your confirmation this session is okay if it is noisy how many frames are out").

Report only: reads the pupil tree the tables use (concat_choice.active_pupil_root) and the videos; changes nothing.

Per session -> analysis/pupil_check/day<N>/day<N>_<cage>_<mouse>.png
  row 1   4 frames spread over the session (1/8, 3/8, 5/8, 7/8) with the table's pupil
  row 2   4 frames spread over the session's "out" frames (none -> 4 more spread frames)
  strip   the table's diameter over the session; red ticks = held frames, orange = jumps / off eye
  header  where the table's pupil comes from, the counts, and Claude's check (verdicts csv, see below)
"Out" frames - the table still has a value there (filled from the neighbours), but it is not a fresh measure:
  pupil (package, or the dark-contrast repair of a listed session)
    held   the tracker kept its last value: Diameter_pixels equal to the frame before in the CLEAN file.
           The package's qc.json splits them into "no pupil near the last place" (fracNoBlob) and "refused: a jump
           of > 18 px or a size change of > 18 %" (fracRejected).
    jump   the tracker took it, but the pupil is > 0.5 x the session's median diameter from the centre of the
           surrounding 2 s, or its size differs > 30 % from there: a short wrong blob (whisker, lid, reflection)
    odd size  the table's diameter is > 1.5 x or < 0.5 x the session's median (and not a jump): a longer stretch on
           a wrong blob, e.g. the dark area around the pupil (9669m white / 9670f white days 1 and 3). A real, large
           pupil change would also be counted, so look at the frames (morphine days!)
  whole eye (9669m black, 9670f black: eye damage, approved 11-17)
    not found  the eye model found no eye
    off eye    the mask centre is > 20 px from the eye centre of the surrounding minute (the rule of qc_summary.py):
               the spout / whiskers / fur when the damaged eye is closed (user 11-18: keep those frames)
Claude's check: analysis/pupil_check/pupil_check_verdicts.csv (day, cage, mouse, verdict, note, checked) is written
after looking at the pictures; its verdict and note are printed in the header and in the table.
Also writes analysis/pupil_check/pupil_check_days<a>-<b>.csv (one row per session) and
pupil_check_overview_days<a>-<b>.png (one tile per session: the spread frame nearest the median size; one picture per
6 days). Verdict rule (2026-11-18): ok = the circle is on the pupil (whole eye: on the eye opening) and < 10 % of the
frames are out; noisy = >= 10 % out, or wrong spikes / blob switches seen in the pictures. Both stay in every table.

usage: python pupil_check_examples.py [--days 1 2 3 4 5] [--only-overview] [--workers 6]
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

T = Path(__file__).resolve().parent
sys.path.insert(0, str(T))
import concat_choice as CC  # noqa: E402
import pupil_example_sheets as PES  # noqa: E402  (frame grab at the right time, contrast stretch, fonts)

P = T.parent
OUT = P / "analysis" / "pupil_check"
MICE = [(c, m) for c in ("9669m", "9670f") for m in ("black", "green1", "green2", "white")]
FPS = 30.26                       # real camera rate (frames are camera pulses)
TW, TH, GAP = 340, 272, 10        # tile (px)
GREEN, ORANGE, PINK, RED, CYAN = (40, 220, 60), (255, 150, 0), (255, 60, 200), (230, 40, 40), (0, 200, 255)
PKG_CSV = "pupil_diameter_FINAL_U-Net_results_CLEAN.csv"
JUMP_POS, JUMP_SIZE, WIN = 0.5, 0.30, 61          # x median diameter; fraction; frames (about 2 s)
ODD_HI, ODD_LO = 1.5, 0.5                         # x the session's median diameter (odd size, longer than a jump)
OFF_EYE_PX, OFF_WIN = 20, 1800                    # qc_summary.wholeeye_offeye
OVERVIEW_DAYS = 6                                 # days per overview picture


def verdicts():
    f = OUT / "pupil_check_verdicts.csv"
    if not f.exists():
        return {}
    V = pd.read_csv(f, dtype=str, encoding="utf-8-sig").fillna("")
    return {(int(r.day), r.cage.strip().lower(), r.mouse.strip().lower()): (r.verdict.strip(), r.note.strip())
            for r in V.itertuples()}


def grades(root):
    f = root / "pupil_session_quality.csv"
    if not f.exists():
        return {}
    Q = pd.read_csv(f, dtype={"cage": str, "mouse": str})
    return {(int(r.day), r.cage, r.mouse): r for r in Q.itertuples()}


def video_of(sd: Path, day, cage, mouse, qc_name="qc.json", sub="pupil_outputs_general"):
    """the video the tracker used (its qc json), else the largest video of the session folder"""
    q = sd / sub / qc_name if sd is not None else None
    if q is not None and q.exists():
        v = Path(json.loads(q.read_text(encoding="utf-8")).get("video", ""))
        if v.exists():
            return v
    vd = P / "video" / f"day{day}" / cage / mouse
    vids = sorted(list(vd.glob("*.mp4")) + list(vd.glob("*.avi")), key=lambda p: p.stat().st_size) if vd.is_dir() else []
    return vids[-1] if vids else None


def rolling_median(x, w, min_periods=None):
    mp = max(5, w // 6) if min_periods is None else min_periods
    return pd.Series(x).rolling(w, center=True, min_periods=mp).median().to_numpy()


# ---------------------------------------------------------------- what the table has, per session
def pupil_session(day, cage, mouse, root, sd):
    """package pupil or the dark-contrast repair of a listed session"""
    dk = CC.darkcontrast_pupil_csv(sd, root)
    if dk is not None:
        src, clean = "dark-contrast repair (approved 11-18)", dk
        F = pd.read_csv(sd / CC.DARK_DIR / "pupil_dark_frames.csv",
                        usecols=["frame", "pkg_d", "pkg_cx", "pkg_cy", "merged_d"])
        raw = pd.DataFrame(dict(d=F.pkg_d, cx=F.pkg_cx, cy=F.pkg_cy))       # what the package's detector saw
        nodet = F.merged_d.isna().to_numpy()
    else:
        hits = sorted((sd / "pupil_outputs_general").glob(PKG_CSV), key=lambda p: p.stat().st_mtime)
        if not hits:
            return None
        src, clean = "package", hits[-1]
        H = pd.read_csv(sd / "pupil_outputs_general" / "hybrid_detections.csv")
        raw = pd.DataFrame(dict(d=H.best_diam_px, cx=H.best_cx_px, cy=H.best_cy_px))
        nodet = (H["accepted_src_1unet_2dark_0none"].to_numpy() == 0)
    C = pd.read_csv(clean)
    q = json.loads((sd / "pupil_outputs_general" / "qc.json").read_text(encoding="utf-8"))
    d_raw = pd.to_numeric(C.Diameter_pixels, errors="coerce").to_numpy(float)
    d = pd.to_numeric(C.SmoothedDiameter_px, errors="coerce").to_numpy(float)
    cx = pd.to_numeric(C.CenterX_px, errors="coerce").to_numpy(float)
    cy = pd.to_numeric(C.CenterY_px, errors="coerce").to_numpy(float)
    n = len(C)
    held = np.r_[False, d_raw[1:] == d_raw[:-1]] | ~np.isfinite(d_raw)
    med = float(np.nanmedian(d))
    mx = rolling_median(np.where(held, np.nan, cx), WIN)
    my = rolling_median(np.where(held, np.nan, cy), WIN)
    md = rolling_median(np.where(held, np.nan, d_raw), WIN)
    with np.errstate(invalid="ignore", divide="ignore"):
        jump = ~held & ((np.hypot(cx - mx, cy - my) > JUMP_POS * med) | (np.abs(d_raw / md - 1) > JUMP_SIZE))
        odd = ~held & ~jump & ((d > ODD_HI * med) | (d < ODD_LO * med))
    roi = q.get("roi")
    crop = PES.roi_crop([float(v) for v in roi]) if roi else None
    out = dict(kind="pupil", src=src, clean=clean, n=n, med=med, d=d, cx=cx, cy=cy, raw=raw, crop=crop,
               held=held, jump=jump, odd=odd, out=held | jump | odd, nodet=nodet,
               pkg_noblob=float(q.get("fracNoBlob", np.nan)), pkg_rej=float(q.get("fracRejected", np.nan)),
               video=video_of(sd, day, cage, mouse))
    if src.startswith("dark"):
        S = json.loads((sd / CC.DARK_DIR / "source.json").read_text(encoding="utf-8"))
        out["dark"] = S
    return out


def wholeeye_session(day, cage, mouse, root, sd):
    we = CC.wholeeye_pupil_csv(sd, root)
    wd = P / "analysis" / "wholeeye" / f"day{day}" / cage / mouse
    if we is None or not (wd / "wholeeye_frames.csv").exists():
        return None
    W = pd.read_csv(wd / "wholeeye_frames.csv")
    J = json.loads((wd / "wholeeye_qc.json").read_text(encoding="utf-8"))
    C = pd.read_csv(we)
    fnd = W.found.to_numpy() == 1
    cx, cy = np.where(fnd, W.cx, np.nan), np.where(fnd, W.cy, np.nan)
    with np.errstate(invalid="ignore"):
        off = fnd & (np.hypot(cx - rolling_median(cx, OFF_WIN, 100), cy - rolling_median(cy, OFF_WIN, 100))
                     > OFF_EYE_PX)
    x0, y0, w, h = [int(v) for v in J["crop"]]
    v = Path(J.get("video", ""))
    qc = {int(p.stem[1:].split("_")[0]): p for p in (wd / "qc").glob("f*_mask.png")}
    return dict(kind="eye", src="whole eye (eye damage, approved 11-17)", clean=we, n=len(W),
                med=float(np.nanmedian(pd.to_numeric(C.SmoothedDiameter_px, errors="coerce"))),
                d=pd.to_numeric(C.SmoothedDiameter_px, errors="coerce").to_numpy(float), W=W, crop=(x0, y0, w, h),
                notfound=~fnd, offeye=off, held=~fnd, jump=off, odd=np.zeros(len(W), bool), out=(~fnd) | off, qc=qc,
                video=v if v.exists() else video_of(None, day, cage, mouse))


def session(day, cage, mouse, root, G):
    sd = CC.session_dir(root, day, cage, mouse)
    base = dict(day=day, cage=cage, mouse=mouse)
    g = G.get((day, cage, mouse))
    base["grade"] = str(g.grade) if g is not None else "-"
    # 2026-11-20: a mouse whose session of that day has not run yet (or is not filed yet) is not "no video"
    pl = P / "python" / f"day{day}" / cage / mouse
    if not (pl.is_dir() and any(pl.glob("*.jsonl"))):
        return dict(base, kind="none", table_run="", note="no session yet (not run or not filed)")
    vd = P / "video" / f"day{day}" / cage / mouse
    if not (vd.is_dir() and (any(vd.glob("*.mp4")) or any(vd.glob("*.avi")))):
        runs = CC.usable_runs(sd, root) if sd is not None else []
        return dict(base, kind="none", table_run=runs[-1].name if runs else "",
                    note="no video (IC Capture did not record) - behaviour-only table, pupil left empty (approved)")
    if sd is None:
        return dict(base, kind="none", note="no pupil output")
    runs = CC.usable_runs(sd, root)
    base["table_run"] = runs[-1].name if runs else ""
    if base["grade"] == "whole_eye" or CC.wholeeye_pupil_csv(sd, root) is not None:
        s = wholeeye_session(day, cage, mouse, root, sd)
    else:
        s = pupil_session(day, cage, mouse, root, sd)
    if s is None:
        return dict(base, kind="none", note="no pupil trace")
    return dict(base, **s)


# ---------------------------------------------------------------- frames to show
def spread(n, k=4, off=0.5):
    return [min(n, max(1, int(round((j + off) / k * n)) + 1)) for j in range(k)]


def _spread_of(arr, m):
    if not len(arr) or m <= 0:
        return []
    idx = np.unique(np.clip(np.round((np.arange(m) + 0.5) / m * len(arr) - 0.5).astype(int), 0, len(arr) - 1))
    return [int(arr[i]) for i in idx]


def out_frames(s, k=4):
    """k frames spread over the out frames; the slots shared round-robin by the kinds that have frames
    (held / no eye first, then odd size, then jumps / off eye)"""
    kinds = [np.flatnonzero(s[x]) + 1 for x in ("held", "odd", "jump")]
    kinds = [a for a in kinds if len(a)]
    share = [0] * len(kinds)
    for j in range(k if kinds else 0):
        share[j % len(kinds)] += 1
    pick = []
    for a, m in zip(kinds, share):
        pick += _spread_of(a, m)
    if len(pick) < k:                                  # one kind has fewer frames than its share: take more out frames
        rest = np.array([x for x in np.flatnonzero(s["out"]) + 1 if x not in set(pick)])
        pick += _spread_of(rest, k - len(pick))
    return sorted(set(pick))[:k]


def grab_eye(video, crop, k):
    x0, y0, w, h = crop
    return PES.grab(video, (x0, y0, w, h), k)


def draw_pupil_tile(img, s, k, size=(TW, TH)):
    row = dict(d=s["d"][k - 1], cx=s["cx"][k - 1], cy=s["cy"][k - 1])
    r2 = s["raw"].iloc[k - 1]
    raw = dict(d=r2.d, cx=r2.cx, cy=r2.cy)
    show_raw = bool(s["out"][k - 1]) or (np.isfinite(raw["d"]) and np.isfinite(row["d"]) and
                                         np.hypot(raw["cx"] - row["cx"], raw["cy"] - row["cy"]) > 3)
    t = PES.tile(img, row, size, raw if show_raw and np.isfinite(raw["d"]) else None, width=3)
    return t


def draw_eye_tile(img, s, k, size=(TW, TH)):
    tw, th = size
    if img is None:
        return PES.tile(None, None, size)
    sc = min(tw / img.shape[1], th / img.shape[0])
    w, h = int(img.shape[1] * sc), int(img.shape[0] * sc)
    ox, oy = (tw - w) // 2, (th - h) // 2
    t = Image.new("RGB", (tw, th), (0, 0, 0))
    t.paste(Image.fromarray(PES.stretch(img)).convert("RGB").resize((w, h), Image.LANCZOS), (ox, oy))
    dr = ImageDraw.Draw(t)
    r = s["W"].iloc[k - 1]
    x0, y0 = s["crop"][0], s["crop"][1]
    if k in s["qc"]:                                   # the model's real mask (saved for 24 frames per session)
        import cv2
        m = np.array(Image.open(s["qc"][k]).convert("L")) > 0
        mm = cv2.resize(m.astype(np.uint8), (w, h), interpolation=cv2.INTER_NEAREST)
        cs, _ = cv2.findContours(mm, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        for c in cs:
            pts = [(ox + int(x), oy + int(y)) for x, y in c[:, 0, :]]
            if len(pts) > 2:
                dr.line(pts + [pts[0]], fill=PINK, width=3)
    if r.found == 1:                                   # box of the mask + centre (every frame)
        cx, cy = ox + (r.cx - x0 + 0.5) * sc - 0.5, oy + (r.cy - y0 + 0.5) * sc - 0.5
        bw, bh = r.width_px * sc / 2, r.height_px * sc / 2
        col = ORANGE if s["offeye"][k - 1] else PINK
        if k not in s["qc"]:
            for a0 in range(0, 360, 20):
                dr.arc([cx - bw, cy - bh, cx + bw, cy + bh], a0, a0 + 10, fill=col, width=3)
        dr.line([cx - 6, cy, cx + 6, cy], fill=col, width=2)
        dr.line([cx, cy - 6, cx, cy + 6], fill=col, width=2)
    else:
        dr.rectangle([0, 0, 190, 30], fill=(0, 0, 0))
        dr.text((6, 5), "no eye found", fill=ORANGE, font=PES.font(18, True))
    return t


def tile_of(s, img, k, size=(TW, TH)):
    return draw_eye_tile(img, s, k, size) if s["kind"] == "eye" else draw_pupil_tile(img, s, k, size)


def why_out(s, k):
    i = k - 1
    if s["kind"] == "eye":
        return "no eye found" if s["notfound"][i] else ("off eye" if s["offeye"][i] else "")
    if s["held"][i]:
        return "held (no pupil)" if s["nodet"][i] else "held (refused jump/size)"
    if s["odd"][i]:
        return "odd size (> 1.5x median)" if s["d"][i] > s["med"] else "odd size (< 0.5x median)"
    return "jump" if s["jump"][i] else ""


# ---------------------------------------------------------------- numbers
def numbers(s):
    n = s["n"]
    o = dict(n_frames=n, median_px=round(s["med"], 2), out_n=int(s["out"].sum()), out_pct=round(100 * s["out"].mean(), 2))
    mins = np.array_split(s["out"], max(1, int(round(n / FPS / 60))))
    o["worst_minute_out_pct"] = round(100 * max(m.mean() for m in mins), 1)
    k910 = min(n, int(round(910 * FPS)))            # every table and metric uses the first 15 min 10 s
    o["out_first910s_pct"] = round(100 * s["out"][:k910].mean(), 2)
    # 2026-11-27: the counts and the median of the same first 15:10, so a window day's tile shows one time base
    o["out_first910s_n"], o["n_first910s"] = int(s["out"][:k910].sum()), int(k910)
    dd = np.asarray(s.get("d", []), dtype=float)[:k910]
    o["median_first910s_px"] = round(float(np.nanmedian(dd)), 2) if np.isfinite(dd).any() else np.nan
    if s["kind"] == "eye":
        o.update(held_n=0, jump_n=0, odd_n=0, notfound_n=int(s["notfound"].sum()), offeye_n=int(s["offeye"].sum()))
    else:
        o.update(held_n=int(s["held"].sum()), jump_n=int(s["jump"].sum()), odd_n=int(s["odd"].sum()), notfound_n=0,
                 offeye_n=0,
                 held_nodet_n=int((s["held"] & s["nodet"]).sum()),
                 pkg_noblob_pct=round(100 * s["pkg_noblob"], 2), pkg_rejected_pct=round(100 * s["pkg_rej"], 2))
    return o


def count_line(s, o):
    n = o["n_frames"]
    pct = lambda x: f"{x:,} ({100 * x / n:.1f} %)"
    line = _count(s, o, pct, n)
    if n > 910 * FPS + 300:                         # longer session: the tables use its first 15 min 10 s
        line += f"  ·  first 15:10: {o['out_first910s_pct']:.1f} %"
    return line


def _count(s, o, pct, n):
    if s["kind"] == "eye":
        return (f"out {pct(o['out_n'])} of {n:,} frames  =  no eye found {pct(o['notfound_n'])}  +  off eye (spout / "
                f"whiskers / fur) {pct(o['offeye_n'])}  ·  worst minute {o['worst_minute_out_pct']:g} %")
    return (f"out {pct(o['out_n'])} of {n:,} frames  =  held {pct(o['held_n'])}  +  jumps {pct(o['jump_n'])}"
            f"  +  odd size {pct(o['odd_n'])}  ·  worst minute {o['worst_minute_out_pct']:g} %")


def off_eye_decision(day, cage, mouse):
    """date of the user's 'keep the off-eye frames' answer for this session (analysis/days/qc_decisions.csv, cat
    off_eye; empty day = every day of that mouse), or None. Added 2026-11-19: the picture said 'your decision 11-18'
    for every whole-eye session, also for 9670f black, whose day-6 answer came on 11-19."""
    for r in CC.read_rows(P / "analysis" / "days" / "qc_decisions.csv", warn=None,
                          need=("cage", "mouse", "cat", "decision")):
        if (r["cat"].strip().lower() == "off_eye" and r["cage"].strip().lower() == str(cage).lower()
                and r["mouse"].strip().lower() == str(mouse).lower()
                and (not str(r.get("day", "")).strip() or int(r["day"]) == int(day))):
            w = str(r.get("decided", "")).strip()
            return w[5:] if len(w) == 10 else (w or "on file")
    return None


def detail_line(s, o):
    if s["kind"] == "eye":
        dec = off_eye_decision(s["day"], s["cage"], s["mouse"])
        return ("held = the eye model found no eye; off eye = the mask centre > 20 px from the eye centre of that "
                "minute (" + (f"kept as they are, your decision {dec}" if dec else "kept in the tables") + ")")
    if s["src"].startswith("dark"):
        S = s["dark"]
        a = S.get("actions", {})
        return (f"repair: package kept where it is on the pupil ({a.get('package_pct', 0):.0f} %), replaced where it "
                f"sat on another blob ({a.get('replaced_pct', 0):.1f} %); held = no pupil that passed the size check "
                f"({100 * o['held_nodet_n'] / o['n_frames']:.1f} %) or refused by the tracker")
    return (f"held = the tracker kept its last value: no pupil near the last place {o['pkg_noblob_pct']:.1f} % · refused "
            f"(jump > 18 px or size > 18 %) {o['pkg_rejected_pct']:.1f} % (package qc.json)")


# ---------------------------------------------------------------- pictures
def session_png(s, frames, V):
    """frames: {k: image}"""
    ks1, ks2 = s["ks1"], s["ks2"]
    W = 4 * (TW + GAP) + GAP
    HEAD, CAP, STRIP = 150, 44, 210
    v, note = V.get((s["day"], s["cage"], s["mouse"]), ("", ""))
    # 2026-11-24 (review): the note ran off the right edge; wrap it to the picture width (at most 3 lines)
    chk = (f"Claude's check: {v.upper()}" if v else "Claude's check: (not looked at yet)") + (f" - {note}" if note else "")
    chk = " ".join(chk.split())         # a line break in a note would stop Pillow's textlength (review 11-24)
    fnt = PES.font(18, True)
    probe = ImageDraw.Draw(Image.new("RGB", (8, 8)))
    wrapped, cur = [], ""
    for word in chk.split(" "):
        t = (cur + " " + word).strip()
        if cur and probe.textlength(t, font=fnt) > W - 2 * GAP:
            wrapped.append(cur); cur = word
        else:
            cur = t
    wrapped.append(cur)
    if len(wrapped) > 3:
        wrapped = wrapped[:3]; wrapped[2] = wrapped[2].rsplit(" ", 1)[0] + " …"
    extra = 24 * (len(wrapped) - 1)
    HEAD += extra
    H = HEAD + 2 * (TH + CAP + GAP) + STRIP + 40
    S = Image.new("RGB", (W, H), "white")
    dr = ImageDraw.Draw(S)
    o = s["nums"]
    col = {"ok": (30, 110, 50), "noisy": (150, 105, 0)}.get(v.split()[0].lower() if v else "", (90, 90, 90))
    dr.text((GAP, 8), f"Day {s['day']} · {s['cage']} {s['mouse']} · table = {s['src']} · median {s['med']:.1f} px",
            fill=(0, 0, 0), font=PES.font(24, True))
    dr.text((GAP, 42), count_line(s, o), fill=(0, 0, 0), font=PES.font(18, True))
    dr.text((GAP, 68), detail_line(s, o), fill=(70, 70, 70), font=PES.font(15))
    for i_, line_ in enumerate(wrapped):
        dr.text((GAP, 92 + 24 * i_), line_, fill=col, font=fnt)
    leg = ("pink = the eye model's mask (dotted: the box of the mask; orange = off eye)" if s["kind"] == "eye" else
           "green = the pupil in the table · orange dotted = what the detector found in that frame (shown when "
           "out or different)")
    dr.text((GAP, 120 + extra), "row 1: spread over the session · row 2: out frames · " + leg, fill=(90, 90, 90),
            font=PES.font(14))
    for r_, ks in enumerate((ks1, ks2)):
        y = HEAD + r_ * (TH + CAP + GAP)
        for j, k in enumerate(ks):
            x = GAP + j * (TW + GAP)
            S.paste(tile_of(s, frames.get(k), k), (x, y))
            mm, ss = divmod((k - 1) / FPS, 60)
            why = why_out(s, k)
            dv = s["d"][k - 1]
            txt = f"{int(mm):02d}:{int(ss):02d} · f{k} · " + (f"{dv:.0f} px" if np.isfinite(dv) else "no value")
            dr.text((x + 3, y + TH + 4), txt, fill=(0, 0, 0), font=PES.font(17, True))
            if why:
                dr.text((x + 3, y + TH + 24), "OUT: " + why, fill=RED if why.startswith(("held", "no eye")) else
                        (150, 60, 200) if why.startswith("odd") else (200, 110, 0), font=PES.font(15, True))
        if r_ == 1 and not s["ks2_out"]:
            dr.text((GAP, y + TH + 24), "no out frames: 4 more spread frames", fill=(30, 110, 50), font=PES.font(15, True))
    # strip: diameter over time, out ticks, example markers
    y0 = HEAD + 2 * (TH + CAP + GAP)
    x0, x1 = GAP + 46, W - GAP
    yb, yt = y0 + STRIP - 34, y0 + 40
    dr.text((GAP, y0), "table diameter (px) over the session · red ticks = held / no eye · purple = odd size · "
                       "orange = jumps / off eye · triangles = the frames above (blue row 1, red row 2)",
            fill=(60, 60, 60), font=PES.font(14))
    n = s["n"]
    d = s["d"]
    lo, hi = np.nanpercentile(d, [0.5, 99.5]) if np.isfinite(d).any() else (0, 1)
    lo, hi = lo - 2, hi + 2
    X = lambda k: x0 + (x1 - x0) * (k - 1) / max(1, n - 1)
    Y = lambda v: yb - 16 - (yb - 16 - yt) * (v - lo) / max(1e-6, hi - lo)
    dr.rectangle([x0, yt, x1, yb], outline=(200, 200, 200))
    for v_ in np.linspace(lo, hi, 3):
        dr.text((GAP, Y(v_) - 8), f"{v_:.0f}", fill=(120, 120, 120), font=PES.font(12))
    step = max(1, n // (x1 - x0))
    pts = [(X(k), Y(d[k - 1])) for k in range(1, n + 1, step) if np.isfinite(d[k - 1])]
    if len(pts) > 1:
        dr.line(pts, fill=(40, 40, 40), width=1)
    for arr, colr, yy in ((s["held"], RED, yb - 12), (s["odd"], (150, 60, 200), yb - 9), (s["jump"], ORANGE, yb - 6)):
        idx = np.flatnonzero(arr)
        xs = np.unique(np.round(x0 + (x1 - x0) * idx / max(1, n - 1)).astype(int))
        for xx in xs:
            dr.line([(xx, yy), (xx, yy + 5)], fill=colr, width=1)
    for j, k in enumerate(ks1 + ks2):
        xx = X(k)
        dr.polygon([(xx - 5, yt - 9), (xx + 5, yt - 9), (xx, yt - 1)], fill=(0, 90, 200) if j < 4 else RED)
    for mnt in range(0, int(n / FPS / 60) + 1, 2):
        xx = X(mnt * 60 * FPS + 1)
        if xx <= x1:
            dr.text((xx - 6, yb + 4), f"{mnt}", fill=(120, 120, 120), font=PES.font(12))
    dr.text((x1 - 60, yb + 4), "minutes", fill=(120, 120, 120), font=PES.font(12))
    dr.text((GAP, H - 26), f"pupil tree {CC.active_pupil_root(P).name} · table run {s.get('table_run', '')} · "
                           f"{s['clean'].name} · tools\\pupil_check_examples.py · {datetime.now():%Y-%m-%d %H:%M}",
            fill=(140, 140, 140), font=PES.font(12))
    return S


def none_png(s):
    W, H = 4 * (TW + GAP) + GAP, 300
    S = Image.new("RGB", (W, H), "white")
    dr = ImageDraw.Draw(S)
    dr.text((GAP, 10), f"Day {s['day']} · {s['cage']} {s['mouse']}", fill=(0, 0, 0), font=PES.font(24, True))
    dr.text((GAP, 60), s["note"], fill=(90, 90, 90), font=PES.font(20, True))
    return S


_WINDOW = None


def window_days() -> dict:
    """{day: seconds} of tools/step02_window_days.csv (approved rows; day 13 = 20-min sessions, 2026-11-26): days
    whose verdicts and tables use the first N s. Read once; missing file = {}."""
    global _WINDOW
    if _WINDOW is None:
        _WINDOW = {}
        f = T / "step02_window_days.csv"
        if f.exists():
            for _, r in pd.read_csv(f, dtype=str, encoding="utf-8-sig").fillna("").iterrows():
                if str(r.get("approved", "")).strip():
                    try:
                        _WINDOW[int(r["day"])] = float(r["seconds"])
                    except (KeyError, ValueError):
                        pass
    return _WINDOW


def overview(recs, frames, V, days, out):
    tw, th, cap = 300, 240, 58
    LAB, HEAD = 92, 130             # 92 (was 70): 'day 10' / 'day 11' lost their last digit (review 2026-11-24)
    W = LAB + len(MICE) * (tw + 8)
    H = HEAD + len(days) * (th + cap + 8) + 30
    S = Image.new("RGB", (W, H), (250, 250, 250))
    dr = ImageDraw.Draw(S)
    dr.text((12, 10), f"Pupil the tables use, days {days[0]}-{days[-1]}: one frame per session + Claude's check",
            fill=(20, 20, 20), font=PES.font(26, True))
    # 2026-11-18: the user asked "what is yellow and green which one you used?" - say what each colour is
    dr.text((12, 46), "CIRCLES: green = the pupil in the tables (USED; package, 9670f green2 day 5 = dark-contrast "
                      "repair) · orange dotted = the detector's blob in that frame, NOT used (refused by the tracker, "
                      "before its smoothing, or the blob the repair replaced) · pink = whole eye (eye damage)",
            fill=(60, 60, 60), font=PES.font(16))
    dr.text((12, 68), "BORDER + caption = Claude's check: green = OK, yellow = NOISY, grey = no video · out % = frames "
                      "held by the tracker (its last value) + jumps + odd sizes (these stay in the table) / no eye + off eye",
            fill=(60, 60, 60), font=PES.font(16))
    for j, (c, m) in enumerate(MICE):
        dr.text((LAB + j * (tw + 8) + 6, 100), f"{c} {m}", fill=(20, 20, 20), font=PES.font(19, True))
    for i, d in enumerate(days):
        y = HEAD + i * (th + cap + 8)
        dr.text((10, y + th // 2 - 12), f"day {d}", fill=(20, 20, 20), font=PES.font(21, True))
        for j, (c, m) in enumerate(MICE):
            x = LAB + j * (tw + 8)
            s = recs.get((d, c, m))
            if s is None or s["kind"] == "none":
                t = Image.new("RGB", (tw, th), (225, 225, 225))
                ImageDraw.Draw(t).text((10, th // 2 - 20), "no video\n(behaviour only)" if s else "no data",
                                       fill=(90, 90, 90), font=PES.font(18, True))
                S.paste(t, (x, y))
                continue
            k = s["k_best"]
            S.paste(tile_of(s, frames.get((d, c, m), {}).get(k), k, (tw, th)), (x, y))
            v, _ = V.get((d, c, m), ("", ""))
            vv = v.split()[0].lower() if v else ""
            colr = {"ok": (40, 150, 60), "noisy": (215, 160, 0)}.get(vv, (150, 150, 150))
            dr.rectangle([x - 3, y - 3, x + tw + 2, y + th + 2], outline=colr, width=4)
            o = s["nums"]
            op = 100 * o["out_n"] / o["n_frames"] if o.get("n_frames") else o["out_pct"]   # one rounding (11-24)
            wt = ""
            if window_days().get(int(d)) == 910.0 and o.get("out_first910s_pct") is not None:
                op, wt = o["out_first910s_pct"], " (15:10)"    # 2026-11-27: a window day's verdict window (day 13)
            dr.text((x + 2, y + th + 4), (v.upper() if v else "not checked") + f" · out {op:.1f} %{wt}",
                    fill=colr if vv else (90, 90, 90), font=PES.font(17, True))
            if wt and o.get("n_first910s"):     # 2026-11-27: counts + median of the first 15:10 too (one time base)
                cnt = (f"{o['out_first910s_n']:,} of {o['n_first910s']:,} frames (15:10) · median "
                       f"{o['median_first910s_px']:.0f} px")
            else:
                cnt = f"{o['out_n']:,} of {o['n_frames']:,} frames · median {s['med']:.0f} px"
            dr.text((x + 2, y + th + 26), cnt,
                    fill=(60, 60, 60), font=PES.font(14))
    dr.text((12, H - 24), "per-session pictures: analysis\\pupil_check\\dayN\\ · numbers: pupil_check_days*.csv · "
                          "tools\\pupil_check_examples.py", fill=(120, 120, 120), font=PES.font(13))
    S.save(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--days", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    ap.add_argument("--only-overview", action="store_true", help="numbers + overview only (no per-session pictures)")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args(argv)
    root = CC.active_pupil_root(P)
    G, V = grades(root), verdicts()
    recs, jobs = {}, []
    for d in a.days:
        for c, m in MICE:
            s = session(d, c, m, root, G)
            recs[(d, c, m)] = s
            if s["kind"] == "none":
                continue
            s["nums"] = numbers(s)
            s["ks1"] = spread(s["n"], 4)
            ko = out_frames(s, 4)
            s["ks2_out"] = bool(ko)
            s["ks2"] = ko if ko else spread(s["n"], 4, off=0.0)[:4]
            if len(s["ks2"]) < 4:
                s["ks2"] += [k for k in spread(s["n"], 8, 0.25) if k not in s["ks2"]][: 4 - len(s["ks2"])]
            dd = s["d"]
            cand = s["ks1"]
            s["k_best"] = min(cand, key=lambda k: abs(dd[k - 1] - s["med"]) if np.isfinite(dd[k - 1]) else 1e9)
            for k in sorted(set(s["ks1"] + s["ks2"])):
                jobs.append((s, k))
    print(f"{len(recs)} sessions, {len(jobs)} frames to read")

    def grab(job):
        s, k = job
        if s["video"] is None or s["crop"] is None:
            return None
        return PES.grab(s["video"], s["crop"], k)

    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        imgs = list(ex.map(grab, jobs))
    frames = {}
    for (s, k), img in zip(jobs, imgs):
        frames.setdefault((s["day"], s["cage"], s["mouse"]), {})[k] = img
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for (d, c, m), s in recs.items():
        v, note = V.get((d, c, m), ("", ""))
        row = dict(day=d, cage=c, mouse=m, grade=s.get("grade", "-"), source=s.get("src", s.get("note", "")),
                   table_run=s.get("table_run", ""), verdict=v, note=note)
        if s["kind"] != "none":
            row.update(s["nums"])
        rows.append(row)
        if a.only_overview:
            continue
        o = OUT / f"day{d}" / f"day{d}_{c}_{m}.png"
        o.parent.mkdir(parents=True, exist_ok=True)
        (none_png(s) if s["kind"] == "none" else session_png(s, frames.get((d, c, m), {}), V)).save(o)
    tag = f"days{a.days[0]}-{a.days[-1]}" if len(a.days) > 1 else f"day{a.days[0]}"
    T_ = pd.DataFrame(rows)
    for c in ("n_frames", "median_px", "out_n", "out_pct", "worst_minute_out_pct", "out_first910s_pct", "held_n",
              "jump_n", "odd_n", "notfound_n", "offeye_n", "out_first910s_n", "n_first910s", "median_first910s_px"):     # same columns when no session has a pupil yet
        if c not in T_:
            T_[c] = np.nan
    T_.to_csv(OUT / f"pupil_check_{tag}.csv", index=False, encoding="utf-8-sig")
    outs = []
    for i in range(0, len(a.days), OVERVIEW_DAYS):          # one overview picture per 6 days (it stays readable)
        ch = a.days[i:i + OVERVIEW_DAYS]
        o = OUT / f"pupil_check_overview_days{ch[0]}-{ch[-1]}.png"
        overview(recs, frames, V, ch, o)
        outs.append(o.name)
    print(f"-> {OUT}\\pupil_check_{tag}.csv, {', '.join(outs)}, day*/ ({len(rows)} sessions)")
    return 0


def overview_files(days):
    """the overview pictures main() writes for these days (one per 6 days)"""
    return [OUT / f"pupil_check_overview_days{days[i]}-{days[min(i + OVERVIEW_DAYS, len(days)) - 1]}.png"
            for i in range(0, len(days), OVERVIEW_DAYS)]


if __name__ == "__main__":
    raise SystemExit(main())
