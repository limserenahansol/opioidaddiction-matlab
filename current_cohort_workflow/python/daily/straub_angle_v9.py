"""Straub tail angle, v9.1 (2026-11-23): v8h below + the v9 features (lift, floor gap, tilt) and rules A-C from the
user's answers of 2026-11-23, with the independent-review fixes (v9.1, frozen). History:
Straub tail angle, v8 (2026-11-22). User on the day-9 videos: "when it is still straub you said no as mouse watching back side..... you should make correct line of body"; "when mouse is rearing or grooming body is vertical we can not use that frame so not consider that frame for counting for straub or no straub or total frame"; "always try to use Body and tail angle not floor"; 10 frames answered by the user (2 checks). Changes to v6b (straub_angle.py, unchanged):
  a) BODY = the mouse pixels inside the convex hull of its dark fur (darkness > BODY_FRAC x its own darkness), then opened: a grey reflection / mouse behind the acrylic wall is not body; the lighter belly stays;
  b) tail THINNESS = 2 x area / perimeter (an arched tail was dropped as too thick); a long piece from the upper half of the body is preferred to a floor / foot piece; a piece leaving the bottom and pointing steeply down (< -50 deg) is a leg (no decision);
  c) a VERTICAL body (elongation >= VERT_ELONG 1.25, axis within EXCLUDE_DEG 45 of vertical = rearing / grooming) is NOT COUNTED (not Straub, not no-Straub, not in the total; straub_auto_v8.py never fills such frames);
  d) BODY LINE = from the head end (the body pixels farthest from the tail root) to the tail root; tail angle = between the tail's first part and that line continued past the root (up = positive); Straub = > 25 deg and the tail rising above the floor. The floor line only when the body line is near vertical / short (seen from the front or behind): Straub = > 25 deg above the floor;
  e) the tail root hidden behind the back: a tail leaving the TOP 10 % of the body = Straub; a tail leaving the upper quarter = Straub when its visible part rises > 10 deg (or goes forward over the back), not Straub when level; a tail from the upper quarter standing > 20 px above the body top and rising > 30 deg = Straub.
Checked on 31 frames (10 user answers, 2 user drawings of 11-21, 19 eye checks): 28 agree, 1 wrong, 2 unmeasured.
Everything else as v6b:
Straub tail angle, v5 (user 11-21: "you can just crop only mouse part and then also increase contrast").
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
VIEW_DEG = 35.0      # v8f: body line (head -> tail root) within this of vertical = seen from the front / behind
LINE_MIN = 0.6       # v8f: body line shorter than this x the body length = not measurable
VERT_ELONG = 1.25    # v8h: elongation from which a near-vertical body counts as rearing / grooming
TOP_FRAC = 0.10      # v8e: tail leaving the body in its top 10 % (user 11-22, frames 2/3 vs 5) = raised behind the back
EXCLUDE_DEG = 45.0   # v8: body axis closer to vertical than this (and elongated) = rearing / grooming, not counted
BODY_FRAC = 0.6      # v8: body = darkness > this x the mouse's darkness (reflections are 0.3-0.5)
HOOK_PX, HOOK_LIFT = 10.0, 0.15   # v9: hook / loop = highest point > 10 px above both ends and >= 0.15 body heights up
FLOOR_PX, FLOOR_FRAC = 12.0, 0.08  # v9: tail lowest point this close to the feet line = on the floor
V9 = True             # v9 rules on (False = v8h decisions with the v9 features, for comparison)
FLOOR_MAX_FRAC = 0.30  # v9 rule A: tail lowest point within this x body height of the feet line = on / near the floor
TOP_ARCH_FRAC, TOP_ARCH_LEN = 0.12, 120.0   # v9 rule B: long tail leaving the very top of the body, any direction
TILT_MAX = 25.0       # v9 rule C: body line tilted more than this = front / back view -> floor reference
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
    # v8: body = the mouse pixels inside the convex hull of its DARK fur (darkness > BODY_FRAC x its own darkness): a
    # grey reflection / a mouse behind the acrylic wall (0.3-0.5 x) lies outside the hull, the lighter belly of a
    # sitting mouse lies inside it and stays; then opened as before (removes the tail)
    fur = (m & (d > BODY_FRAC * dark)).astype(np.uint8)
    pts = cv2.findNonZero(fur)
    if pts is None:
        return dict(ok=False, why="no body")
    hull = np.zeros_like(fur)
    cv2.fillConvexPoly(hull, cv2.convexHull(pts), 1)
    body = cv2.morphologyEx((m & hull).astype(np.uint8), cv2.MORPH_OPEN, DISK)
    nb, lb, sbb, _ = cv2.connectedComponentsWithStats(body, 8)
    if nb < 2:
        return dict(ok=False, why="no body")
    body = (lb == 1 + int(np.argmax(sbb[1:, cv2.CC_STAT_AREA]))).astype(np.uint8)
    ys, xs = np.nonzero(body)
    c = np.array([xs.mean(), ys.mean()])
    ev, evec = np.linalg.eigh(np.cov(np.stack([xs, ys])))
    axis, elong = evec[:, 1], float(np.sqrt(ev[1] / max(ev[0], 1e-6)))
    half = float(np.abs((np.stack([xs, ys], 1) - c) @ axis).max())
    up0 = abs(float(np.degrees(np.arctan2(abs(axis[0]), abs(axis[1]))))) < UPRIGHT_DEG
    res = dict(ok=False, c=c, axis=axis if (elong >= ELONG_MIN and not up0) else np.array([1.0, 0.0]), body_axis=axis,
               elong=elong, why="no tail", crop=(X0, Y0, X1, Y1), weak_thr=weak_thr, excluded=False)
    # v8 (user 2026-11-22: "when mouse is rearing or grooming body is vertical we can not use that frame so not consider
    # that frame for counting for straub or no straub or total frame"): an elongated body within UPRIGHT_DEG of vertical
    # is not judged at all - not Straub, not no-Straub, not in the total
    vert = abs(float(np.degrees(np.arctan2(abs(axis[0]), abs(axis[1]))))) < EXCLUDE_DEG
    is_vert = bool(elong >= VERT_ELONG and vert)          # v8h: 1.25 (user 11-22, 2nd check frame 2: sitting up, 1.32)

    def not_counted():                                    # rearing / grooming without a clearly raised tail
        res.update(ok=False, axis=axis, why="body vertical (rearing / grooming) - not counted", excluded=True)
        return res
    thin = m & (1 - cv2.dilate(body, DISK2))
    nt, lt, stt, _ = cv2.connectedComponentsWithStats(thin, 8)
    db = cv2.distanceTransform((1 - body).astype(np.uint8), cv2.DIST_L2, 3)
    best, blen, base = None, 0.0, None
    cands = []
    ytop0, ybot0 = float(ys.min()), float(ys.max())
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
        pc, _ = cv2.findContours((lt == j).astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        per = max(float(sum(cv2.arcLength(q, True) for q in pc)), 1.0)
        if ln < 1 or 2.0 * a / per > MAX_WIDTH:                # v8: ribbon width, not area / straight length
            continue
        if elong >= ELONG_MIN and abs(float((b0 - c) @ axis)) < REAR_FRAC * half:
            continue
        tot = ln + float(dv[k0])
        cands.append((tot, j, b0))
    long_ = [q for q in cands if q[0] >= MIN_LEN]
    upper = [q for q in long_ if (q[2][1] - ytop0) / max(ybot0 - ytop0, 1.0) < 0.5]
    if cands:                                             # v8e: a long piece from the upper half first
        blen, best, base = max(upper or long_ or cands, key=lambda q: q[0])
    if best is None:
        return not_counted() if is_vert else res
    if blen < MIN_LEN:
        if is_vert:
            return not_counted()
        res["why"] = "tail too short"
        return res
    ty, tx = np.nonzero(lt == best)
    r = np.hypot(tx - base[0], ty - base[1])
    sel = (r >= SEG[0]) & (r <= SEG[1])
    if sel.sum() < 8:
        if is_vert:
            return not_counted()
        res["why"] = "tail too short"
        return res
    tip = np.array([tx[sel].mean(), ty[sel].mean()])
    tdir = (tip - base) / (np.linalg.norm(tip - base) + 1e-9)
    if np.dot(base - c, axis) < 0:
        axis = -axis
    rise = float(np.degrees(np.arctan2(-tdir[1], abs(tdir[0]) + 1e-9)))   # elevation above the floor line
    ytop, ybot = float(ys.min()), float(ys.max())
    high = bool(base[1] < c[1] - HIGH_FRAC * (ybot - ytop))
    relh = float((base[1] - ytop) / max(ybot - ytop, 1.0))  # 0 = top of the body, 1 = bottom
    # v8f (user 2026-11-22: "always try to use Body and tail angle not floor"; the body line as the user drew it on
    # 2026-11-21: from the head to the tail ROOT): head end = the body pixels farthest from where the tail leaves the
    # body; body line = head -> root; tail angle = between the tail's first part and that line continued past the root,
    # positive towards the back (up). The floor line is used ONLY when the body line cannot be measured in 2D: it lies
    # within VIEW_DEG of vertical or is short (mouse seen from the front / from behind)
    d2 = (xs - base[0]) ** 2 + (ys - base[1]) ** 2
    far = d2 >= np.quantile(d2, 0.97)
    head = np.array([xs[far].mean(), ys[far].mean()], float)
    bl = base - head
    bl_len = float(np.linalg.norm(bl))
    u = bl / (bl_len + 1e-9)                               # head -> root, continued backwards past the root
    nrm = np.array([u[1], -u[0]])
    if nrm[1] > 0:
        nrm = -nrm                                         # the normal pointing up (towards the back)
    body_ang = float(np.degrees(np.arctan2(float(tdir @ nrm), float(tdir @ u))))
    line_vert = abs(float(np.degrees(np.arctan2(abs(u[0]), abs(u[1]))))) < VIEW_DEG
    line_ok = bool((not line_vert) and bl_len >= LINE_MIN * 2 * half)
    res.update(head=head, body_line_ok=line_ok)
    # v9 features (2026-11-23, user answers on 9669m): how high the tail goes above its root (lift, in body heights),
    # a hook / loop (the tail's highest point above BOTH its ends), how close its lowest point comes to the floor (the
    # mouse's lowest pixel = the feet line), the elevation from the root to the tail's highest point, the body-line tilt
    Hb = max(ybot - ytop, 1.0)
    kf, kt = int(np.argmax(r)), int(np.argmin(ty))
    far_pt, top_pt = np.array([tx[kf], ty[kf]], float), np.array([tx[kt], ty[kt]], float)
    lift = float((base[1] - top_pt[1]) / Hb)
    hook = bool(top_pt[1] < min(base[1], far_pt[1]) - HOOK_PX and lift >= HOOK_LIFT)
    feet_y = float(np.nonzero(m)[0].max())
    floor_gap = float(feet_y - ty.max())
    floor_touch = bool(floor_gap <= max(FLOOR_PX, FLOOR_FRAC * Hb))
    elev_top = float(np.degrees(np.arctan2(base[1] - top_pt[1], abs(top_pt[0] - base[0]) + 1e-9)))
    tilt = float(np.degrees(np.arctan2(abs(u[1]), abs(u[0]) + 1e-9)))
    below = float(ty.max() - ybot)                # > 0: the tail goes lower than the body's lowest pixel
    res.update(lift=lift, hook=hook, floor_gap=floor_gap, floor_touch=floor_touch, elev_top=elev_top, tilt=tilt,
               Hb=Hb, relh=relh, tail_len=blen, below=below)
    # v9.1 (2026-11-23, independent review, confirmed): rule A must not take (1) a clearly RAISED tail whose root is
    # low (hunched mouse: the lowest tail pixel is the root itself, but the tail rises > 25 deg and stands >= 0.3 body
    # heights above the root - the user's 'no counting' tails all had lift <= 0.21), nor (2) a piece the v8 checks
    # call 'not a tail' (a leg pointing down from the bottom, a piece pointing forward), which v8 leaves unmeasured
    is_leg = bool(relh > 0.8 and rise < -50.0)
    is_fwd = bool(line_ok and abs(body_ang) > FORWARD_DEG and not (body_ang > 0 and rise > 45.0))
    raised = bool(rise > STRAUB_DEG and lift >= 0.3)
    if V9 and not is_vert and not is_leg and not is_fwd and not raised and floor_gap <= FLOOR_MAX_FRAC * Hb \
            and relh > 0.12:
        # v9 rule A (user 2026-11-23, triangle answers on 9669m green1 207.1 s, black 233.2 s / 254.0 s: tail low or on
        # the floor while hunched or walking = "no counting"): the tail's lowest point near the feet line = not counted
        res.update(ok=False, why="tail on / near the floor - not counted", excluded=True, base=base, tdir=tdir,
                   rise=rise, angle=float(body_ang if line_ok else rise))
        return res
    if is_vert:
        # v8h (user 11-22, 2nd check frame 2: vertical body with the tail raised = "do not count"; the user's rule:
        # "when mouse is rearing or grooming body is vertical ... not consider that frame"): a vertical body is never
        # counted, whatever the tail does (the v8e exception for a raised tail, from the 1st check frame 4 "rearing but
        # I can see its straub", is withdrawn - the two answers are the same case)
        return not_counted()
    if V9 and relh <= TOP_ARCH_FRAC and blen >= TOP_ARCH_LEN and rise <= -10:
        # v9 rule B (user 2026-11-23, 9669m green1 250.2 s = Straub: the tail arches over the back and hangs down in
        # front; v8h took the hanging part, -61 deg): a LONG tail leaving the very top of the body = over the back
        res.update(ok=True, base=base, tdir=tdir, axis=u if line_ok else np.array([1.0, 0.0]), body_axis=axis,
                   angle=float(body_ang if line_ok else rise), rise=rise, ref="long tail over the top of the body (arch)",
                   tail_len=blen, why="", straub=True, high_start=True, tail_mask=(lt == best), relh=relh)
        return res
    if relh <= TOP_FRAC and rise > -10:
        # v8e (user 11-22, frames 2 / 3 Straub vs 5 not): a tail leaving the TOP of the body comes from behind / over
        # the back = Straub even when its visible part runs level
        res.update(ok=True, base=base, tdir=tdir, axis=u if line_ok else np.array([1.0, 0.0]), body_axis=axis,
                   angle=float(body_ang if line_ok else rise), rise=rise, ref="tail leaves the top of the body",
                   tail_len=blen, why="", straub=True, high_start=True, tail_mask=(lt == best), relh=relh)
        return res
    if relh > 0.8 and rise < -50.0:
        # v8g: a piece leaving the bottom of the body and pointing steeply down is a leg / foot, not the tail (day-9
        # check, green1 19.4 s: the blurred tail was not found and a foot was measured) - no decision for this frame
        res["why"] = "piece points down from the bottom (a leg, not the tail)"
        return res
    if relh < HIGH_FRAC and float(ty.min()) < ytop - 20 and rise > 30.0:
        # v8g: a tail leaving the upper body and standing above the top of the body = raised above the back (day-9
        # check, green1 52.1 s: mouse seen from behind, tail straight up, the body line was not the real one)
        res.update(ok=True, base=base, tdir=tdir, axis=u if line_ok else np.array([1.0, 0.0]), body_axis=axis,
                   angle=float(body_ang if line_ok else rise), rise=rise, ref="tail raised above the back",
                   tail_len=blen, why="", straub=True, high_start=True, tail_mask=(lt == best), relh=relh)
        return res
    if relh < HIGH_FRAC and rise > -10:
        # v8h (user answers 11-22: 1st check frame 5 = level tail from the upper body NOT Straub; 2nd check frame 4 =
        # tail from the upper body rising 12 deg = Straub): a tail leaving the UPPER quarter of the body has its root
        # hidden behind the back, so the angle to the body line cannot be read from the visible part; the visible
        # part rising (> 10 deg above the floor) or going forward over the back = lifted from the root = Straub
        fwd = float((tdir @ u)) < 0 if line_ok else False
        straub = bool(rise > RISE_DEG or fwd)
        res.update(ok=True, base=base, tdir=tdir, axis=u if line_ok else np.array([1.0, 0.0]), body_axis=axis,
                   angle=float(rise), rise=rise, ref="tail leaves the upper body (root hidden): rise", tail_len=blen,
                   why="", straub=straub, high_start=True, tail_mask=(lt == best), relh=relh)
        return res
    if line_ok:
        # body line and tail angle (the user's measure): Straub = the tail raised > 25 deg above the body line continued
        # past the root, and rising above the floor (a tail lying on the floor is not Straub)
        ang, ref = body_ang, "body line (head to tail root)"
        if abs(ang) > FORWARD_DEG and not (ang > 0 and rise > 45.0):
            res["why"] = "piece points forward (not a tail)"
            return res
        straub = bool(ang > STRAUB_DEG and rise > RISE_DEG)
        draw_axis = u
        if V9 and tilt >= TILT_MAX:
            # v9 rule C (user 2026-11-23, 9669m white 259.1 s = Straub: mouse half from the front, the 'body line' ran
            # 34 deg up to the root and the tail read -3 deg): a body line tilted this much is not the side-view body
            # axis, so the tail is judged against the floor (> 25 deg), as for front / back views
            ang, ref = rise, f"body line tilted {tilt:.0f} deg (front / back view): floor"
            straub = bool(rise > STRAUB_DEG)
            draw_axis = np.array([1.0, 0.0])
    else:
        # no body line in 2D (seen from the front / behind): the floor line; a steep tail = more than 25 deg (user)
        ang, ref = rise, "floor (no body line: seen from the front / behind)"
        straub = bool(rise > STRAUB_DEG)
        draw_axis = np.array([1.0, 0.0])
    res.update(ok=True, base=base, tdir=tdir, axis=draw_axis, body_axis=axis, angle=float(ang), rise=rise, ref=ref,
               tail_len=blen, why="", straub=straub, high_start=high, tail_mask=(lt == best), relh=relh)
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
