"""Straub tail angle, v5 (user 11-21: "you can just crop only mouse part and then also increase contrast").
1. find the mouse on the 640x360 frame (darker than the background by > 40; largest blob);
2. CROP the full-resolution frame (1280x720) around it (box + 120 px), with the full-resolution background;
3. CONTRAST: the darkness d = background - frame is scaled by the mouse's own darkness (median d of its body), so the
   tail (lighter than the body) is kept: mouse mask = hysteresis, seeds d > 40, grown through pixels with
   d > max(14, 0.22 x body darkness) that connect to the seeds (isolated noise is not taken);
4. body = mask opened with a disk r 10; tail = the longest thin piece (mask minus dilated body) starting <= 32 px from
   the body, area / length <= 16 px, length from the body edge >= 50 px; at the rear end of an elongated body;
5. angle of the tail's first part (12-56 px from its base) above the body axis (PCA of the body; a round body ->
   horizontal). Straub = angle > 25 deg AND the tail rises above the horizontal by > 10 deg (a tail lying on the floor,
   e.g. while the mouse rears, is not Straub); a 'tail' pointing forward (|angle| > 135 deg) is not a tail.
All lengths in full-resolution pixels."""
import ctypes

import cv2
import numpy as np

try:
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x40)
except Exception:
    pass

SMALL = (640, 360)
SEED = 40
WEAK_MIN, WEAK_FRAC = 14, 0.15
BRIDGE = 24          # weak pieces within this distance of the mouse join it (light tail base)
UPRIGHT_DEG = 35     # body axis within this angle of vertical = rearing / upright -> horizontal reference
PAD = 120
OPEN_R = 10
START_MAX = 32
MIN_LEN = 50
MAX_WIDTH = 16.0
ELONG_MIN = 1.45
REAR_FRAC = 0.4
SEG = (12, 56)
STRAUB_DEG, RISE_DEG, FORWARD_DEG = 25.0, 10.0, 135.0
HIGH_FRAC = 0.25     # tail leaving the body this far above its centre = raised over the back (v6b)
DISK = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * OPEN_R + 1, 2 * OPEN_R + 1))
DISK2 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * OPEN_R + 9, 2 * OPEN_R + 9))


def gray_full(frame):
    return cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)


def background(cap, n, k=60):
    fr = []
    for i in np.linspace(0, n - 1, k).astype(int):
        cap.set(cv2.CAP_PROP_POS_FRAMES, int(i))
        ok, f = cap.read()
        if ok:
            fr.append(gray_full(f))
    return np.median(np.stack(fr), axis=0).astype(np.float32)


def measure(gf, bgf):
    Hf, Wf = gf.shape
    sx, sy = Wf / SMALL[0], Hf / SMALL[1]
    gs = cv2.resize(gf, SMALL, interpolation=cv2.INTER_AREA).astype(np.float32)
    bs = cv2.resize(bgf, SMALL, interpolation=cv2.INTER_AREA)
    ms = ((bs - gs) > SEED).astype(np.uint8)
    n, lab, st, _ = cv2.connectedComponentsWithStats(ms, 8)
    if n < 2:
        return dict(ok=False, why="no mouse")
    big = 1 + int(np.argmax(st[1:, cv2.CC_STAT_AREA]))
    if st[big, cv2.CC_STAT_AREA] < 400:
        return dict(ok=False, why="no mouse")
    x, y, w, h = st[big, :4]
    X0, Y0 = max(0, int(x * sx) - PAD), max(0, int(y * sy) - PAD)
    X1, Y1 = min(Wf, int((x + w) * sx) + PAD), min(Hf, int((y + h) * sy) + PAD)
    d = bgf[Y0:Y1, X0:X1] - gf[Y0:Y1, X0:X1].astype(np.float32)
    seed = (d > SEED).astype(np.uint8)
    ns, ls, ss, _ = cv2.connectedComponentsWithStats(seed, 8)
    if ns < 2:
        return dict(ok=False, why="no mouse")
    sb = 1 + int(np.argmax(ss[1:, cv2.CC_STAT_AREA]))
    core = (ls == sb)
    dark = float(np.median(d[core]))
    weak_thr = max(WEAK_MIN, WEAK_FRAC * dark)
    weak = (d > weak_thr).astype(np.uint8)
    nw, lw, _, _ = cv2.connectedComponentsWithStats(weak, 8)
    keep = np.unique(lw[core & (lw > 0)])
    m = np.isin(lw, keep).astype(np.uint8)                 # hysteresis: weak pixels connected to the mouse's core
    dm = cv2.distanceTransform((1 - m).astype(np.uint8), cv2.DIST_L2, 3)
    for j in range(1, nw):                                 # pieces cut off by a light gap (the tail base)
        if j in keep:
            continue
        pj = lw == j
        if pj.sum() >= 30 and dm[pj].min() <= BRIDGE:
            m[pj] = 1
    body = cv2.morphologyEx(m, cv2.MORPH_OPEN, DISK)
    nb, lb, sbb, _ = cv2.connectedComponentsWithStats(body, 8)
    if nb < 2:
        return dict(ok=False, why="no body")
    body = (lb == 1 + int(np.argmax(sbb[1:, cv2.CC_STAT_AREA]))).astype(np.uint8)
    ys, xs = np.nonzero(body)
    c = np.array([xs.mean(), ys.mean()])
    ev, evec = np.linalg.eigh(np.cov(np.stack([xs, ys])))
    axis, elong = evec[:, 1], float(np.sqrt(ev[1] / max(ev[0], 1e-6)))
    half = float(np.abs((np.stack([xs, ys], 1) - c) @ axis).max())
    res = dict(ok=False, c=c, axis=axis, elong=elong, why="no tail", crop=(X0, Y0, X1, Y1), weak_thr=weak_thr)
    thin = m & (1 - cv2.dilate(body, DISK2))
    nt, lt, stt, _ = cv2.connectedComponentsWithStats(thin, 8)
    db = cv2.distanceTransform((1 - body).astype(np.uint8), cv2.DIST_L2, 3)
    best, blen, base = None, 0.0, None
    for j in range(1, nt):
        a = stt[j, cv2.CC_STAT_AREA]
        if a < 40:
            continue
        py, px = np.nonzero(lt == j)
        dv = db[py, px]
        k0 = int(np.argmin(dv))
        if dv[k0] > START_MAX:
            continue
        b0 = np.array([px[k0], py[k0]], float)
        ln = float(np.hypot(px - b0[0], py - b0[1]).max())
        if ln < 1 or a / ln > MAX_WIDTH:
            continue
        if elong >= ELONG_MIN and abs(float((b0 - c) @ axis)) < REAR_FRAC * half:
            continue
        tot = ln + float(dv[k0])
        if tot > blen:
            best, blen, base = j, tot, b0
    if best is None:
        return res
    if blen < MIN_LEN:
        res["why"] = "tail too short"
        return res
    ty, tx = np.nonzero(lt == best)
    r = np.hypot(tx - base[0], ty - base[1])
    sel = (r >= SEG[0]) & (r <= SEG[1])
    if sel.sum() < 8:
        res["why"] = "tail too short"
        return res
    tip = np.array([tx[sel].mean(), ty[sel].mean()])
    tdir = (tip - base) / (np.linalg.norm(tip - base) + 1e-9)
    if np.dot(base - c, axis) < 0:
        axis = -axis
    back = 1.0 if base[0] >= c[0] else -1.0
    rise = float(np.degrees(np.arctan2(-tdir[1], abs(tdir[0]) + 1e-9)))   # elevation above the horizontal
    upright = abs(float(np.degrees(np.arctan2(abs(axis[0]), abs(axis[1]))))) < UPRIGHT_DEG
    if elong < ELONG_MIN or upright:
        ang, ref = float(np.degrees(np.arctan2(-tdir[1], tdir[0] * back))), ("upright" if upright and elong >= ELONG_MIN else "horizontal")
    else:
        a_ax = np.degrees(np.arctan2(-axis[1], axis[0]))
        a_t = np.degrees(np.arctan2(-tdir[1], tdir[0]))
        ang = (a_t - a_ax + 180) % 360 - 180
        if axis[0] < 0:
            ang = -ang
        ang, ref = float(ang), "body axis"
    if abs(ang) > FORWARD_DEG:
        res["why"] = "piece points forward (not a tail)"
        return res
    straub = bool(ang > STRAUB_DEG and rise > RISE_DEG)
    # v6b (user drawing 2026-11-21, 84.2 s): a tail that leaves the body outline HIGH on the body (above the centre
    # by > 25 % of the body height) is raised from its root behind / over the back, even if its visible part runs
    # nearly level; count it as Straub unless it points down
    ytop, ybot = float(ys.min()), float(ys.max())
    high = bool(base[1] < c[1] - HIGH_FRAC * (ybot - ytop))
    if not straub and high and rise > -10:
        straub = True
    res.update(ok=True, base=base, tdir=tdir, axis=axis, angle=ang, rise=rise, ref=ref, tail_len=blen, why="",
               straub=straub, high_start=high, tail_mask=(lt == best))
    return res


def draw(gf, r, out_w=540, out_h=320):
    """crop view around the mouse (or the middle) with the tail in magenta, body axis blue, tail line red/green"""
    im = cv2.cvtColor(gf, cv2.COLOR_GRAY2BGR)
    if "crop" in r:
        X0, Y0, X1, Y1 = r["crop"]
        sub = im[Y0:Y1, X0:X1]
        if r.get("ok"):
            sub[r["tail_mask"]] = (0.4 * sub[r["tail_mask"]] + np.array([150, 60, 150])).astype(np.uint8)
        c = r["c"]
        cv2.line(sub, tuple(int(v) for v in c - r["axis"] * 110), tuple(int(v) for v in c + r["axis"] * 110), (255, 160, 0), 2)
        if r.get("ok"):
            b, t = r["base"], r["base"] + r["tdir"] * 80
            cv2.line(sub, tuple(int(v) for v in b), tuple(int(v) for v in t), (0, 0, 255) if r["straub"] else (0, 170, 0), 3)
            cv2.circle(sub, tuple(int(v) for v in b), 5, (0, 255, 255), -1)
        cx, cy = X0 + c[0], Y0 + c[1]
    else:
        cx, cy = im.shape[1] / 2, im.shape[0] / 2
    x0 = int(np.clip(cx - out_w, 0, im.shape[1] - 2 * out_w))
    y0 = int(np.clip(cy - out_h * 1.2, 0, im.shape[0] - 2 * out_h))
    return cv2.resize(im[y0:y0 + 2 * out_h, x0:x0 + 2 * out_w], (out_w, out_h), interpolation=cv2.INTER_AREA)
