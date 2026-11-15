"""Fully dilated pupils: screen sessions for them, and measure the full pupil in the sessions the user approved.

Rule (from the user, 2026-11-15): when the whole eye is dark, the pupil is fully dilated and the pupil is the WHOLE
dark area inside the eyelids. Most sessions have a small pupil; only some mice show this. The package's tracker cannot
measure it (its U-Net takes a small dark part, and its step 3 keeps only 30 px around the last centre, so no pupil can
be wider than 60 px; a full pupil on this rig is ~120 px).

This is an OPTION for rare sessions (the user, 2026-11-15: "rare event or rare session - do not make this rule for
all the others, just make options to be used in this circumstance"): nothing changes for a session that is not
listed in fullpupil_sessions.csv, and the screen below runs only when asked (--screen); it is not part of the daily run.

1) SCREEN (--screen; sessions of --screen-days, default all): 40 frames per session, the whole dark area inside the lids and the
   tracker's pupil (hybrid_detections.csv). In a normal eye the thin ring just outside the tracked pupil is lighter
   iris; when the pupil fills the eye the ring is as dark as the pupil. Writes <pupil_root>/fullpupil_screen.csv; a
   session with ring_dark >= 0.40 and ring_contrast <= 0.17 (medians) is listed as "look at it" - nothing is changed.
   (On the first cohort days it flagged every fully dilated session, but also some normal eyes: a look decides.)
1b) MOUSE OPTION (automatic; <pupil_root>/fullpupil_mice.csv: cage, mouse, approved, check_60px, reason): every
   session of a listed mouse is looked at on every run. With check_60px = yes (default) a session gets the full-pupil
   measure only when the black region around the tracker's pupil is larger than 60 px - the largest pupil the
   package's tracker can report - as a median over 40 frames (days 1-2: normal days 18-56 px, dilated days 69-93 px);
   otherwise the package's pupil is kept. The decisions go to fullpupil_auto.csv / fullpupil_check.csv. An explicit
   row in fullpupil_sessions.csv that is not approved ('no ...') overrides an automatic yes.
2) LISTED sessions (<pupil_root>/fullpupil_sessions.csv: day, cage, mouse, approved, reason; a row counts only when
   'approved' starts with 'yes' or a date YYYY-MM-DD; see concat_choice.py), plus the sessions of 1b:
   a) current_cohort_workflow/matlab/pupil_fullpupil_session.m -> <session>/pupil_outputs_fullpupil/
      (pupil_fullpupil_CLEAN.csv has the same columns as the package's CLEAN file; the daily analysis uses it for
      sessions graded full_pupil, see pupil_session_quality.py). fullpupil_source.json records the inputs (sha1 of
      qc.json, of the package's CLEAN file and of the code) and the outputs (sha1 and rows of the CLEAN file).
   b) Step 02 copy concat_out_<stamp>_fullpupil: the session's table (concat_choice.py: the package's, or its approved
      sync-aligned copy) with ONLY the pupil columns replaced (PupilTimestamp_s, Diameter_px <- SmoothedDiameter_px,
      CenterX_px, CenterY_px, UsedScoreFallback); every other cell is copied as text. Not published when the trace is
      poor (more than half of the frames without a value, or a median under 10 px). Written to a temporary folder with
      its fix_source.json, checked, then renamed in one step. prepare_step02_base.py copies the session's table.
   Redone only when its inputs changed (package re-tracked the session, new ROI, new code, new base table).
   Removing a session from the list (or 'approved' not yes / a date) takes it back to the package's table everywhere.

usage: python fullpupil_sessions.py --project <project> [--dry-run]          (listed sessions only)
       python fullpupil_sessions.py --project <project> --screen-only [--screen-days 2]   (option: screen)
"""
import argparse, csv, hashlib, json, re, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import concat_choice as CC  # noqa: E402  (which concat_out is the session's table; approvals)

OUT_DIR, CLEAN = "pupil_outputs_fullpupil", "pupil_fullpupil_CLEAN.csv"
DERIVED = CC.DERIVED
SUFFIX, TMP_PREFIX = CC.FULL, "_fullpupil_tmp_"
MATLAB_TIMEOUT_S = 3600
MAX_NAN_FRAC, MIN_MEDIAN_PX = 0.5, 10.0     # quality gate before a trace goes into Step 02
MICE_LIST, CHECK_FILE = "fullpupil_mice.csv", "fullpupil_check.csv"
CHECK_PX = 60.0     # the package's tracker cannot report a pupil wider than 60 px (step 3: 30-px search radius)
FLAG_RING_DARK, FLAG_RING_CONTRAST = 0.40, 0.17
COLMAP = {"PupilTimestamp_s": "Timestamp_seconds", "Diameter_px": "SmoothedDiameter_px", "CenterX_px": "CenterX_px",
          "CenterY_px": "CenterY_px", "UsedScoreFallback": "UsedScoreFallback"}


def matlab_exe():
    return shutil.which("matlab") or next((str(p) for p in sorted(Path("C:/Program Files/MATLAB").glob("R*/bin/matlab.exe"))), None)


def workflow_dir():
    for c in (HERE / "opioidaddiction-matlab" / "current_cohort_workflow", HERE.parent.parent):
        if (c / "matlab" / "pupil_fullpupil_session.m").exists():
            return c
    return None


def sha1(p: Path) -> str:
    h = hashlib.sha1()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def is_approved(v) -> bool:
    return CC.is_approved(v)          # 'yes...' or a date YYYY-MM-DD


def valid_derived(run: Path) -> bool:
    return run.name.endswith(DERIVED) and (run / "fix_source.json").exists() and (run / "combined_pupil_digital.csv").exists()


# ----------------------------------------------------------------------------- 1) screen
def _otsu(v):
    hist, edges = np.histogram(v, bins=128)
    c = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist); w1 = w0[-1] - w0
    s = np.cumsum(hist * c)
    m0 = s / np.maximum(w0, 1); m1 = (s[-1] - s) / np.maximum(w1, 1)
    return c[np.argmax(w0 * w1 * (m0 - m1) ** 2)]


def _disk(r):
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    return x * x + y * y <= r * r


def dark_area(C):
    """python version of pupil_finetune/fullpupil_mask.m (for the screen only)"""
    from scipy import ndimage as ndi
    S = ndi.gaussian_filter(C.astype(float), 1.5)
    G = ndi.binary_dilation(S >= np.percentile(S, 99), _disk(3))
    A = ndi.binary_opening(S < _otsu(S[~G]), _disk(4))
    lab, n = ndi.label(A)
    h, w = A.shape
    if n == 0:
        return np.zeros_like(A), S, G
    cen = np.zeros_like(A); cen[h // 4:3 * h // 4, w // 4:3 * w // 4] = True
    best = np.argmax(ndi.sum(cen, lab, index=np.arange(1, n + 1))) + 1
    m = ndi.binary_fill_holes(lab == best)
    mc = ndi.binary_fill_holes(ndi.binary_closing(m, _disk(6)))
    return ndi.binary_fill_holes(m | (mc & G)), S, G


def darkest_part(C, A, S, G, seed_yx):
    """diameter (px) of the black region around the tracker's pupil: inside the dark area A, pixels darker than halfway
    between its darkest 5 % and its 75th percentile, opened (disk 3), the piece at the tracker's centre, holes filled.
    In a normal eye this is the pupil (about 20-30 px on this rig); when the pupil is fully dilated it is the large black
    region (60-95 px on days 1-2), i.e. larger than the 60 px the package's tracker can report."""
    from scipy import ndimage as ndi
    v = S[A & ~G]
    if v.size < 50:
        return np.nan
    dark, mid = np.percentile(v, 5), np.percentile(v, 75)
    B = ((S < dark + 0.5 * (mid - dark)) | G) & A
    B = ndi.binary_opening(B, _disk(3))
    lab, n = ndi.label(B)
    if n == 0:
        return np.nan
    sy, sx = int(np.clip(round(seed_yx[0]), 0, B.shape[0] - 1)), int(np.clip(round(seed_yx[1]), 0, B.shape[1] - 1))
    l = lab[sy, sx]
    if l == 0:   # nearest piece to the tracker's centre
        idx = ndi.distance_transform_edt(lab == 0, return_indices=True)[1]
        l = lab[idx[0][sy, sx], idx[1][sy, sx]]
    return 2 * np.sqrt(ndi.binary_fill_holes(lab == l).sum() / np.pi)


def screen_session(sd: Path, nf=40):
    import imageio_ffmpeg
    from scipy import ndimage as ndi
    Q = json.loads((sd / "pupil_outputs_general" / "qc.json").read_text(encoding="utf-8"))
    x, y, w, h = [int(round(v)) for v in Q["roi"]]
    D = pd.read_csv(sd / "pupil_outputs_general" / "hybrid_detections.csv")
    ff = imageio_ffmpeg.get_ffmpeg_exe(); W, H = 1920, 1080
    rows = []
    for k in np.linspace(1, len(D), nf + 2)[1:-1].astype(int):
        r = D.iloc[k - 1]
        if not (r["accepted_src_1unet_2dark_0none"] > 0 and np.isfinite(r["best_diam_px"])):
            continue
        raw = subprocess.run([ff, "-hide_banner", "-loglevel", "error", "-ss", f"{(k - 1) / 60.0:.4f}", "-i", Q["video"],
                              "-frames:v", "1", "-f", "rawvideo", "-pix_fmt", "gray", "-"], capture_output=True).stdout
        if len(raw) < W * H:
            continue
        C = np.frombuffer(raw[: W * H], np.uint8).reshape(H, W)[y - 1:y + h, x - 1:x + w]   # imcrop: h+1 x w+1
        A, S, G = dark_area(C)
        if not A.any():
            continue
        yy, xx = np.mgrid[0:C.shape[0], 0:C.shape[1]]
        rad = r["best_diam_px"] / 2
        dist = np.hypot(xx - (r["best_cx_px"] - 1), yy - (r["best_cy_px"] - 1))
        blob = dist <= rad
        ring = (dist > rad + 2) & (dist <= rad + 6) & A & ~G
        out = ~ndi.binary_dilation(A, _disk(3)) & ~G
        if blob.sum() < 10 or ring.sum() < 20 or out.sum() < 50:
            continue
        pl, sk = np.median(S[blob & ~G]) if (blob & ~G).any() else np.median(S[blob]), np.median(S[out])
        span = max(1.0, sk - pl)
        rows.append(dict(full_px=2 * np.sqrt(A.sum() / np.pi), unet_px=r["best_diam_px"],
                         darkest_px=darkest_part(C, A, S, G, (r["best_cy_px"] - 1, r["best_cx_px"] - 1)),
                         ring_contrast=(np.median(S[ring]) - pl) / span, ring_dark=float(np.mean(S[ring] < pl + 0.15 * span))))
    if not rows:
        return None
    R = pd.DataFrame(rows)
    return dict(frames=len(R), unet_px=R.unet_px.median(), full_px=R.full_px.median(), darkest_px=R.darkest_px.median(),
                darkest_q25=R.darkest_px.quantile(0.25), darkest_q75=R.darkest_px.quantile(0.75),
                ring_contrast=R.ring_contrast.median(), ring_dark=R.ring_dark.median())


# ----------------------------------------------------------------------------- 1b) mouse option (automatic)
def mouse_option(root: Path, dry: bool) -> list:
    """Mice listed in fullpupil_mice.csv (cage, mouse, approved, check_60px, reason): EVERY session of such a mouse is
    looked at automatically. With check_60px = yes (default) the session gets the full-pupil measure only when the
    black region around the tracker's pupil is larger than 60 px (the largest pupil the package can report), as a
    median over 40 frames; otherwise the package's pupil is kept. With check_60px = no every session of the mouse gets
    it. Writes fullpupil_auto.csv (the sessions it approves; read by concat_choice.py like fullpupil_sessions.csv, an
    explicit 'no' row there wins) and fullpupil_check.csv (every check with its numbers). A check is re-done only when
    the package re-tracked the session."""
    opts = [o for o in CC.read_rows(root / MICE_LIST, need=("cage", "mouse")) if CC.is_approved(o.get("approved"))]
    cache = {}
    for r in CC.read_rows(root / CHECK_FILE, warn=None):
        if str(r["day"]).strip().isdigit():
            cache[CC.norm_key(r["day"], r["cage"], r["mouse"])] = r
    auto_rows, check_rows, msgs = [], [], []
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    for o in opts:
        need_check = str(o.get("check_60px", "yes")).strip().lower() not in ("no", "n", "false", "0")
        cage, mouse = o["cage"].strip(), o["mouse"].strip()
        for ddir in sorted(p for p in root.glob("day*") if p.is_dir() and p.name[3:].isdigit()):
            day = int(ddir.name[3:])
            sd = CC.session_dir(root, day, cage, mouse)
            pc = sd / "pupil_outputs_general" / "pupil_diameter_FINAL_U-Net_results_CLEAN.csv" if sd else None
            if pc is None or not pc.exists():
                continue
            key, pkg = CC.norm_key(day, cage, mouse), CC.sha1(pc)
            dk, nfr = np.nan, 0
            if need_check:
                c = cache.get(key)
                if c and c.get("package_clean_sha1") == pkg and str(c.get("darkest_px", "")).strip() not in ("", "nan"):
                    dk, nfr = float(c["darkest_px"]), int(float(c.get("frames") or 0))
                else:
                    try:
                        s = screen_session(sd)
                        dk, nfr = (float(s["darkest_px"]), int(s["frames"])) if s else (np.nan, 0)
                    except Exception as exc:   # noqa: BLE001
                        msgs.append(f"  !! mouse option day{day} {cage}_{mouse}: size check failed ({exc}) -> package pupil")
                passed = bool(np.isfinite(dk) and dk > CHECK_PX)
                why = f"black pupil region {dk:.0f} px {'>' if passed else '<='} {CHECK_PX:.0f} px (median of {nfr} frames)"
            else:
                passed, why = True, "mouse option without the size check (check_60px = no)"
            check_rows.append(dict(day=day, cage=cage, mouse=mouse, package_clean_sha1=pkg,
                                   darkest_px="" if not np.isfinite(dk) else round(dk, 1), frames=nfr,
                                   decision="full pupil" if passed else "package pupil", reason=why, checked=now))
            if passed:
                auto_rows.append(dict(day=day, cage=cage, mouse=mouse, approved=f"yes (auto: {MICE_LIST})", reason=why))
            msgs.append(f"  mouse option day{day} {cage}_{mouse}: {why} -> {'full pupil' if passed else 'package pupil'}")
    if not dry and (opts or (root / CC.AUTO_FULL).exists()):
        for name, rows, cols in ((CC.AUTO_FULL, auto_rows, ["day", "cage", "mouse", "approved", "reason"]),
                                 (CHECK_FILE, check_rows, ["day", "cage", "mouse", "decision", "darkest_px", "frames",
                                                           "reason", "package_clean_sha1", "checked"])):
            with open(root / name, "w", newline="", encoding="utf-8") as f:   # rewritten: a removed option has no effect
                w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    return msgs


# ----------------------------------------------------------------------------- 2) approved sessions
def code_hash(wf: Path) -> str:
    return hashlib.sha1(b"".join((wf / p).read_bytes() for p in
                                 ("matlab/pupil_fullpupil_session.m", "pupil_finetune/fullpupil_mask.m"))).hexdigest()


def trace_inputs(sd: Path, wf: Path) -> dict:
    g = sd / "pupil_outputs_general"
    return {"qc_json": sha1(g / "qc.json"), "package_clean": sha1(g / "pupil_diameter_FINAL_U-Net_results_CLEAN.csv"),
            "code": code_hash(wf)}


def make_trace(sd: Path, wf: Path, mexe, dry: bool) -> str:
    od = sd / OUT_DIR
    src = od / "fullpupil_source.json"
    want = trace_inputs(sd, wf)
    if CC.fullpupil_trace_valid(sd) and json.loads(src.read_text(encoding="utf-8")).get("inputs", {}).get("code") == want["code"]:
        return "up to date"
    if dry:
        return "would make the trace (MATLAB)"
    if not mexe:
        raise RuntimeError("MATLAB not found")
    if src.exists():
        src.unlink()                  # the old trace is no longer valid while the new one is made
    q = lambda s: str(s).replace("'", "''")  # noqa: E731
    cmd = f"cd('{q(wf / 'matlab')}'); pupil_fullpupil_session('{q(sd)}');"
    res = subprocess.run([mexe, "-batch", cmd], capture_output=True, text=True, errors="replace", timeout=MATLAB_TIMEOUT_S)
    if res.returncode != 0 or not (od / CLEAN).exists():
        raise RuntimeError("MATLAB: " + re.sub(r"[\x00-\x1f]", " ", (res.stderr or res.stdout)[-500:]))
    n = sum(1 for _ in open(od / CLEAN, encoding="utf-8", errors="replace")) - 1
    src.write_text(json.dumps({"inputs": want, "outputs": {"clean_sha1": sha1(od / CLEAN), "rows": n},
                               "made": datetime.now().isoformat(timespec="seconds")}, indent=1), encoding="utf-8")
    return "made the trace"


def base_table(sd: Path, root: Path):
    """the session's table without full-pupil copies (concat_choice.py: the package's, or its approved sync-aligned copy)"""
    ok = [p for p in CC.usable_runs(sd, root) if not p.name.endswith(SUFFIX)]
    return ok[-1] if ok else None


def raster_png(path: Path, T: pd.DataFrame, title: str):
    """the package's 3-panel figure (pupil, lick / injector tracks, onset raster) for the copy's own table"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = pd.to_numeric(T["PlotTime_s_30fps"], errors="coerce").to_numpy()
    lick = pd.to_numeric(T["Lick_TTL"], errors="coerce").fillna(0).to_numpy() > 0.5
    inj = pd.to_numeric(T["Injector_TTL"], errors="coerce").fillna(0).to_numpy() > 0.5
    on = lambda b: np.flatnonzero(b & ~np.r_[False, b[:-1]])  # noqa: E731
    fig, ax = plt.subplots(3, 1, figsize=(13, 7), dpi=110, sharex=True)
    ax[0].plot(x, pd.to_numeric(T["Diameter_px"], errors="coerce"), lw=0.8); ax[0].set_ylabel("Pupil diameter (px)")
    ax[0].set_title(title, fontsize=9)
    ax[1].step(x, lick.astype(float), lw=0.8, label="Lick"); ax[1].step(x, inj.astype(float) + 1, lw=0.8, ls="--", label="Injector")
    ax[1].set_yticks([0, 1, 2]); ax[1].set_yticklabels(["0", "Lick", "Injector"]); ax[1].legend(loc="upper right", fontsize=8)
    ax[2].vlines(x[on(lick)], 0.6, 1.4, lw=0.5); ax[2].vlines(x[on(inj)], 1.6, 2.4, lw=0.7, color="tab:orange")
    ax[2].set_yticks([1, 2]); ax[2].set_yticklabels(["Lick", "Injector"]); ax[2].set_xlabel("Time (s)")
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def read_text_csv(p: Path):
    raw = p.read_bytes()
    nl = "\r\n" if b"\r\n" in raw[:4096] else "\n"
    rows = list(csv.reader(raw.decode("utf-8-sig").splitlines()))
    return rows[0], rows[1:], nl


def make_step02_copy(sd: Path, root: Path, dry: bool, needs_sync: bool = False) -> str:
    B = base_table(sd, root)
    if B is None:
        raise RuntimeError("no concat_out to start from")
    if needs_sync and not B.name.endswith("_syncaligned"):
        # the session is also approved for the alignment fix: build on its checked sync-aligned copy only
        raise RuntimeError(f"the alignment fix is approved but the newest table is {B.name}: run fix_concat_alignment.py first")
    if not CC.fullpupil_trace_valid(sd):
        if dry:
            return f"would write a full-pupil Step 02 copy of {B.name} once the trace is (re)made"
        raise RuntimeError("the full-pupil trace is missing or no longer matches the package's output: not published")
    fclean = sd / OUT_DIR / CLEAN
    qc = json.loads((sd / OUT_DIR / "fullpupil_qc.json").read_text(encoding="utf-8"))
    nan_frac, med = float(qc.get("fracNaN", 0) or 0), float(qc.get("medianDiameter", np.nan))
    if not (nan_frac <= MAX_NAN_FRAC and np.isfinite(med) and med >= MIN_MEDIAN_PX):
        raise RuntimeError(f"poor trace (frames without a value {100 * nan_frac:.0f} %, median {med:.1f} px): not published "
                           "to Step 02 - look at pupil_outputs_fullpupil/fullpupil_sheet.png")
    want = {"base_run": B.name, "fullpupil_clean": sha1(fclean)}
    have = sorted(p for p in sd.glob(f"concat_out_*{SUFFIX}") if valid_derived(p))
    if have:
        s = json.loads((have[-1] / "fix_source.json").read_text(encoding="utf-8"))
        if s.get("base_run") == want["base_run"] and s.get("fullpupil_clean") == want["fullpupil_clean"]:
            return f"Step 02 copy up to date ({have[-1].name})"
    if dry:
        return f"would write a full-pupil Step 02 copy of {B.name}"
    hdr, rows, nl = read_text_csv(B / "combined_pupil_digital.csv")
    fh, frows, _ = read_text_csv(fclean)
    if len(frows) != len(rows):
        raise RuntimeError(f"full-pupil trace has {len(frows)} rows, {B.name} has {len(rows)}")
    ci = {c: i for i, c in enumerate(hdr)}; fi = {c: i for i, c in enumerate(fh)}
    for c, f in COLMAP.items():
        if c not in ci or f not in fi:
            raise RuntimeError(f"column {c} / {f} missing")
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    tmp, final = sd / f"{TMP_PREFIX}{stamp}", sd / f"concat_out_{stamp}{SUFFIX}"
    tmp.mkdir()
    try:
        new = []
        for r, fr in zip(rows, frows):
            r = list(r)
            for c, f in COLMAP.items():
                r[ci[c]] = fr[fi[f]]
            new.append(r)
        with open(tmp / "combined_pupil_digital.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, lineterminator=nl); w.writerow(hdr); w.writerows(new)
        for fn in ("lick_onsets.csv", "injector_onsets.csv"):   # made from the digital columns: unchanged
            if (B / fn).exists():
                shutil.copy2(B / fn, tmp / fn)
        if (B / "log.txt").exists():
            (tmp / "log.txt").write_bytes(
                f"[fullpupil_sessions.py {datetime.now():%Y-%m-%d %H:%M}] copy of {B.name} with the pupil columns from "
                f"{OUT_DIR}\\{CLEAN}; the lines below are the log of {B.name}\r\n".encode("utf-8") + (B / "log.txt").read_bytes())
        # checks: same rows and frames; every column except the pupil ones identical as text; pupil = the trace
        H2, R2, _ = read_text_csv(tmp / "combined_pupil_digital.csv")
        keep = [i for c, i in ci.items() if c not in COLMAP]
        if H2 != hdr or len(R2) != len(rows) or any(a[i] != b[i] for a, b in zip(R2, rows) for i in keep):
            raise RuntimeError("non-pupil columns changed")
        if any(a[ci["Diameter_px"]] != b[fi["SmoothedDiameter_px"]] for a, b in zip(R2, frows)):
            raise RuntimeError("Diameter_px is not the full-pupil SmoothedDiameter_px")
        pt_old = np.array([float(r[ci["PupilTimestamp_s"]]) for r in rows])
        pt_new = np.array([float(r[ci["PupilTimestamp_s"]]) for r in R2])
        if np.nanmax(np.abs(pt_old - pt_new)) > 1e-6:
            raise RuntimeError("the full-pupil trace has other frame times than the package (another video?)")
        summ = (B / "summary.txt").read_text(encoding="utf-8", errors="replace") if (B / "summary.txt").exists() else ""
        lines = []
        for ln in summ.splitlines():
            if ln.startswith("Pupil file:"):
                ln = f"Pupil file: {fclean} (full pupil; the base table {B.name} used {ln.split(':', 1)[1].strip()})"
            elif ln.startswith("Figure FIG:"):
                continue                                         # no .fig in the copy (the PNG is made again below)
            lines.append(ln.replace(B.name, final.name))
        lines.append(f"Pupil columns replaced ({datetime.now():%Y-%m-%d %H:%M}) by the full-pupil measure (whole dark "
                     f"area inside the lids, a session listed in fullpupil_sessions.csv); all other columns are those of "
                     f"{B.name}. Frames without a value: {100 * nan_frac:.1f} % (gaps > 1 s are left empty).")
        (tmp / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
        T2 = pd.read_csv(tmp / "combined_pupil_digital.csv")
        raster_png(tmp / "plot_pupil_lick_injector_STACKED_RASTER.png", T2,
                   f"{'/'.join(sd.parts[-3:])}: full pupil (whole dark area inside the lids)")
        d_old = np.nanmedian([float(r[ci["Diameter_px"]]) for r in rows])
        d_new = float(np.nanmedian(pd.to_numeric(T2["Diameter_px"], errors="coerce")))
        # fix_source.json goes in BEFORE the rename: the copy is published in one step
        (tmp / "fix_source.json").write_text(json.dumps(dict(want, made=datetime.now().isoformat(timespec="seconds"),
            made_by="fullpupil_sessions.py", median_diameter_before_px=d_old, median_diameter_after_px=d_new,
            frac_without_value=nan_frac, frac_accepted=qc.get("fracAccepted"), frac_edge_leak=qc.get("fracEdgeLeak"),
            longest_gap_frames=qc.get("longestNaNRun")), indent=1), encoding="utf-8")
        tmp.rename(final)
        return f"Step 02 copy {final.name} (median pupil {d_old:.1f} -> {d_new:.1f} px; without a value {100 * nan_frac:.1f} %)"
    except Exception:
        if tmp.exists():
            try:
                tmp.rename(sd / f"FAILED_fullpupil_{stamp}")
            except OSError:
                pass
        raise


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True, type=Path)
    ap.add_argument("--pupil-root", type=Path)
    ap.add_argument("--screen", action="store_true",
                    help="OPTION (off by default, not run daily): screen sessions for a fully dilated pupil")
    ap.add_argument("--screen-days", type=int, nargs="*", default=None, help="days to screen (default: all)")
    ap.add_argument("--screen-only", action="store_true", help="only the screen, no traces or Step 02 copies")
    ap.add_argument("--no-screen", action="store_true", help="(old flag; the screen is off unless --screen)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    root = a.pupil_root or a.project / "pupil_tst_final"
    bad = 0
    if (a.screen or a.screen_only) and not a.no_screen:
        sp = root / "fullpupil_screen.csv"
        old = pd.read_csv(sp, dtype={"cage": str, "mouse": str}) if sp.exists() else pd.DataFrame()
        rows = []
        for q in sorted(root.glob("day*/*/*/pupil_outputs_general/hybrid_detections.csv")):
            sd = q.parent.parent
            day = int(sd.parts[-3][3:]) if sd.parts[-3][3:].isdigit() else None
            if day is None or (a.screen_days and day not in a.screen_days):
                continue
            try:
                s = screen_session(sd)
            except Exception as exc:   # a screen problem must not stop the daily routine
                print(f"  screen {sd.parts[-3]} {sd.parts[-2]} {sd.parts[-1]}: {exc}")
                continue
            if s:
                flag = s["ring_dark"] >= FLAG_RING_DARK and s["ring_contrast"] <= FLAG_RING_CONTRAST
                rows.append(dict(day=day, cage=sd.parts[-2], mouse=sd.parts[-1], **{k: round(v, 3) for k, v in s.items()},
                                 flag="look: whole eye dark? (fully dilated pupil)" if flag else ""))
        T = pd.DataFrame(rows)
        if len(old) and len(T):
            done = set(zip(T.day, T.cage, T.mouse))
            old = old[[k not in done for k in zip(old.day.astype(int), old.cage, old.mouse)]]
            T = pd.concat([old, T], ignore_index=True)
        if len(T):
            T.sort_values(["day", "cage", "mouse"]).to_csv(sp, index=False)
            for _, r in T[T["flag"].astype(str).str.len() > 0].iterrows():
                if not a.screen_days or int(r["day"]) in a.screen_days:
                    print(f"  !! day{int(r['day'])} {r['cage']} {r['mouse']}: {r['flag']} (ring dark {r['ring_dark']:.2f}, "
                          f"ring contrast {r['ring_contrast']:.2f}; tracked {r['unet_px']:.0f} px, dark area "
                          f"{r['full_px']:.0f} px): look at the frames; if the pupil fills the eye, ask the user and add "
                          f"the session to fullpupil_sessions.csv")
    if a.screen_only:
        return 0
    for m in mouse_option(root, a.dry_run):                         # fullpupil_mice.csv -> fullpupil_auto.csv
        print(m)
    listed = CC.approved_sessions(root, CC.LISTS[CC.FULL])          # explicit + automatic; (day, cage, mouse) lower case
    sync_ok = CC.approved_sessions(root, CC.LISTS[CC.SYNC])
    wf, mexe = workflow_dir(), matlab_exe()
    for key in sorted(listed):
        sd = CC.session_dir(root, *key)
        name = f"day{key[0]} {key[1]}_{key[2]}"
        try:
            if sd is None or not (sd / "pupil_outputs_general" / "qc.json").exists():
                raise RuntimeError("no package output for this session")
            if not wf:
                raise RuntimeError("current_cohort_workflow (pupil_fullpupil_session.m) not found")
            m1 = make_trace(sd, wf, mexe, a.dry_run)
            m2 = make_step02_copy(sd, root, a.dry_run, key in sync_ok) if (sd / OUT_DIR / CLEAN).exists() else "no trace yet"
            print(f"  full pupil {name}: {m1}; {m2}")
        except Exception as exc:
            bad += 1
            print(f"  !! full pupil {name}: FAILED ({exc})")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
