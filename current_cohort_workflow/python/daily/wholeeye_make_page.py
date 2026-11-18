"""wholeeye_make_page.py - builds tools/wholeeye_teach_page.html (drawing the WHOLE EYE outline of 9670f black,
user 2026-11-17) from the pupil drawing page tools/pupil_teach_page.html, which is not changed.
Changes: texts; no green/orange overlays; a free-hand outline (polygon) instead of the ellipse; a click without a
drag saves nothing; the 'none' answer = the eye is fully closed; its own browser-storage keys; port 8766."""
from pathlib import Path

T = Path(__file__).resolve().parent
src = (T / "pupil_teach_page.html").read_text(encoding="utf-8")

R = [
    ("<title>Pupil Supervision</title>", "<title>Eye Outline Drawing</title>"),
    ("<h1>동공 감독</h1>", "<h1>9670f black 눈 전체 그리기</h1>"),
    ('aria-label="동공 표시와 그리기 영역"', 'aria-label="눈 윤곽 그리기 영역"'),
    ('<span><i class="sw track"></i>초록 = 지금 표의 동공</span>', ''),
    ('<span><i class="sw be"></i><span id="legend-be">주황 점선 = 다른 추적</span></span>',
     '<span hidden><i class="sw be"></i><span id="legend-be"></span></span>'),
    ('<span><i class="sw mine"></i>분홍 = 내가 그린 동공</span>', '<span><i class="sw mine"></i>분홍 = 내가 그린 눈 윤곽</span>'),
    ('<dt>표의 동공</dt><dd id="m-table">–</dd>', '<dt hidden>표</dt><dd id="m-table" hidden>–</dd>'),
    ('<button id="a-table" type="button">', '<button id="a-table" type="button" hidden>'),
    ('<span class="grow">동공 안 보임 (감김·가려짐)</span>', '<span class="grow">눈이 완전히 감김 (보이는 눈 없음)</span>'),
    ('<p class="help"><b>틀렸으면 그리세요:</b> 그림 위에서 진짜 동공(눈 안의 가장 검은 타원)을 꼭 맞게 끌어서 감싸면 저장됩니다. '
     '아직 답이 없는 프레임은 한 번 누르면 그 자리를 동공 위치로 둡니다(크기는 정하지 않음). 다시 끌어 그리면 바뀝니다.</p>',
     '<p class="help"><b>그리는 법:</b> 눈꺼풀 사이로 보이는 눈 전체(어둡든 뿌옇게 하얗든)의 윤곽을 따라, 누른 채로 '
     '한 바퀴 그리고 손을 떼면 저장됩니다. 다시 그리면 바뀝니다. 밝기·대비를 올리면 잘 보입니다. '
     '눈이 완전히 감겼으면 N을 누르세요.</p>'),
    (r"analysis\\pupil_teach\\marks.json", r"analysis\\wholeeye\\teach\\marks.json"),
    ("http://127.0.0.1:8765/", "http://127.0.0.1:8766/"),
    ('none: "동공 안 보임"', 'none: "눈 감김"'),
    ('$("#m-session").textContent = s.session ? `못 찾음 ${s.not_found_pct}% · 자리가 먼 것 ${s.wrong_pct}% (전체 ${s.frames} 프레임)` : "–";',
     '$("#m-session").textContent = s.session ? `${s.session} (전체 ${s.frames} 프레임)` : "–";'),
    ('add(s.session, `${s.day}일 ${s.cage} ${s.mouse} · 못 찾음 ${s.not_found_pct}% · ${n}장`);',
     'add(s.session, `${s.day}일 ${s.cage} ${s.mouse} · ${n}장`);'),
    ('const LS = "pupil-supervision-marks-v1";', 'const LS = "wholeeye-marks-v1";'),
    ('const LSP = "pupil-supervision-pending-v1";', 'const LSP = "wholeeye-pending-v1";'),
    ("window.__pupilSupervision", "window.__wholeEye"),
    ('$("#hide").textContent = showOv ? "원 숨기기 " : "원 보이기 ";',
     '$("#hide").textContent = showOv ? "윤곽 숨기기 " : "윤곽 보이기 ";'),
    ('<button id="hide" type="button">원 숨기기 <kbd>H</kbd></button>',
     '<button id="hide" type="button">윤곽 숨기기 <kbd>H</kbd></button>'),
    # ---- drawing: a free-hand outline
    ("""    const m = marks[f.id];
    const e = drag && drag.moved ? box(drag) : (m && m.choice === "draw" && m.ellipse ? m.ellipse : null);
    if (e) {
      svg.appendChild(el("ellipse", {cx: e.cx + 0.5, cy: e.cy + 0.5, rx: Math.max(0.5, e.w / 2), ry: Math.max(0.5, e.h / 2)},
                         "ov-mine" + (drag || (m && m.ellipse && m.ellipse.from_click) ? " temp" : "")));
      const k = Math.max(1.5, Math.min(e.w, e.h) / 6);
      svg.appendChild(el("path", {d: `M${e.cx + 0.5 - k} ${e.cy + 0.5}h${2 * k}M${e.cx + 0.5} ${e.cy + 0.5 - k}v${2 * k}`}, "ov-mine"));
    }
  }""",
     """    const m = marks[f.id];
    if (drag && drag.moved && drag.pts.length > 1) {   // the outline being drawn (SVG units)
      svg.appendChild(el("polyline", {points: drag.pts.map(([x, y]) => `${x},${y}`).join(" ")}, "ov-mine temp"));
      return;
    }
    const e = m && m.choice === "draw" && m.ellipse ? m.ellipse : null;
    if (e && e.poly) {   // a stored outline: crop pixels with centres at integers, drawn at +0.5
      svg.appendChild(el("polygon", {points: e.poly.map(([x, y]) => `${x + 0.5},${y + 0.5}`).join(" ")}, "ov-mine"));
    } else if (e) {
      svg.appendChild(el("ellipse", {cx: e.cx + 0.5, cy: e.cy + 0.5, rx: Math.max(0.5, e.w / 2), ry: Math.max(0.5, e.h / 2)}, "ov-mine"));
    }
  }"""),
    ("    drag = {x0: p.x, y0: p.y, x1: p.x, y1: p.y, id: e.pointerId, sx: e.clientX, sy: e.clientY, moved: false};",
     "    drag = {x0: p.x, y0: p.y, x1: p.x, y1: p.y, id: e.pointerId, sx: e.clientX, sy: e.clientY, moved: false, pts: [[p.x, p.y]]};"),
    ("""    const p = pt(e); drag.x1 = p.x; drag.y1 = p.y;
    if (Math.hypot(e.clientX - drag.sx, e.clientY - drag.sy) >= 6) drag.moved = true;   // click vs drag: screen px
    if (drag.moved) drawOverlays();""",
     """    const p = pt(e); drag.x1 = p.x; drag.y1 = p.y;
    const q = drag.pts[drag.pts.length - 1];
    if (Math.hypot(p.x - q[0], p.y - q[1]) >= 0.6) drag.pts.push([p.x, p.y]);
    if (Math.hypot(e.clientX - drag.sx, e.clientY - drag.sy) >= 6) drag.moved = true;   // click vs drag: screen px
    if (drag.moved) drawOverlays();"""),
    ("""    if (!d.moved) {   // a click
      if (marks[f.id]) { drawOverlays(); toast("이미 답이 있습니다. 다시 그리려면 끌어서 그리세요."); return; }
      const dia = (f.table && f.table.d) || 20;   // shown only; from_click = the position is the answer, not the size
      answer("draw", {cx: r2(d.x0 - 0.5), cy: r2(d.y0 - 0.5), w: r2(dia), h: r2(dia), from_click: true});
      return;
    }
    const ell = box(d);
    if (Math.min(ell.w, ell.h) < 3) { drawOverlays(); toast("너무 가늘어서 저장하지 않았습니다. 동공을 감싸도록 끌어 주세요."); return; }
    answer("draw", ell);""",
     """    if (!d.moved) { drawOverlays(); toast("누른 채로 끌어서 눈 윤곽을 따라 한 바퀴 그려 주세요."); return; }
    const P = d.pts.concat([[p.x, p.y]]);
    let A = 0;
    for (let i = 0; i < P.length; i++) { const [xa, ya] = P[i], [xb, yb] = P[(i + 1) % P.length]; A += xa * yb - xb * ya; }
    A = Math.abs(A) / 2;
    const xs = P.map(v => v[0]), ys = P.map(v => v[1]);
    const X0 = Math.min(...xs), X1 = Math.max(...xs), Y0 = Math.min(...ys), Y1 = Math.max(...ys);
    if (A < 6 || P.length < 4) { drawOverlays(); toast("너무 작아서 저장하지 않았습니다. 눈 윤곽을 따라 한 바퀴 그려 주세요."); return; }
    answer("draw", {poly: P.map(([x, y]) => [r2(x - 0.5), r2(y - 0.5)]), cx: r2((X0 + X1) / 2 - 0.5), cy: r2((Y0 + Y1) / 2 - 0.5),
                    w: r2(X1 - X0), h: r2(Y1 - Y0), area: r2(A)});"""),
]
out = src
for a, b in R:
    n = out.count(a)
    if n < 1:
        raise SystemExit(f"not found in the pupil page:\n{a[:160]}")
    out = out.replace(a, b)
(T / "wholeeye_teach_page.html").write_text(out, encoding="utf-8")
print(f"wrote {T / 'wholeeye_teach_page.html'} ({len(R)} changes)")
