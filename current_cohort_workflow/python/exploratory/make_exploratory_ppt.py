#!/usr/bin/env python3
r"""
make_exploratory_ppt.py - a short Korean deck that explains the exploratory day 1-5 figures (E1-E7) of
explore_days1to5.py. Made 2026-11-19 (Claude); the user: "explain me well by making short brief ppt with figure in kr
of your new plots". Reads analysis/exploratory/<date>/ (figures + results.json) and writes
analysis/ppt/exploratory_days1-5_<date>.pptx. Nothing else is changed.
usage: python make_exploratory_ppt.py [--date 20261119]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Emu, Inches, Pt

P = Path(__file__).resolve().parent.parent.parent
FONT = "Malgun Gothic"
INK, SUB, ACC, BOX, LINE = "1B2A41", "5D6D7E", "0E7C86", "E8F4F5", "C9D6DF"
W, H = 13.333, 7.5


def font(run, size, bold=False, color=INK):
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.color.rgb = RGBColor.from_string(color)
    rPr = run._r.get_or_add_rPr()
    rPr.set("lang", "ko-KR")            # Korean text: PowerPoint then wraps at spaces, not inside a word
    rPr.set("altLang", "en-US")
    for tag in ("a:latin", "a:ea", "a:cs"):
        el = rPr.find(qn(tag))
        if el is None:
            el = rPr.makeelement(qn(tag), {})
            rPr.append(el)
        el.set("typeface", FONT)


def bullet(p, on=True, level=0):
    pPr = p._p.get_or_add_pPr()
    for tag in ("a:buNone", "a:buChar", "a:buFont"):
        for el in pPr.findall(qn(tag)):
            pPr.remove(el)
    if on:
        pPr.set("marL", str(int(Emu(Inches(0.22 + 0.22 * level)))))
        pPr.set("indent", str(int(-Emu(Inches(0.18)))))
        bf = pPr.makeelement(qn("a:buFont"), {"typeface": "Arial"})
        bc = pPr.makeelement(qn("a:buChar"), {"char": "•" if level == 0 else "–"})
        pPr.append(bf)
        pPr.append(bc)
    else:
        pPr.set("marL", "0")
        pPr.set("indent", "0")
        pPr.append(pPr.makeelement(qn("a:buNone"), {}))


def textbox(slide, x, y, w, h, paras, size=14, anchor=MSO_ANCHOR.TOP, fill=None, name=None, pad=0.08):
    """paras: list of (text | [(text, bold, color)], kind) with kind 'b' bullet, 'b2' sub-bullet, 'p' plain,
    'h' small heading"""
    shp = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h)) if fill else \
        slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    if fill:
        shp.fill.solid()
        shp.fill.fore_color.rgb = RGBColor.from_string(fill)
        shp.line.fill.background()
        shp.shadow.inherit = False
    if name:
        shp.name = name
    tf = shp.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    for side in ("margin_left", "margin_right", "margin_top", "margin_bottom"):
        setattr(tf, side, Inches(pad))
    first = True
    for item, kind in paras:
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        runs = item if isinstance(item, list) else [(item, kind == "h", ACC if kind == "h" else INK)]
        for txt, b, c in runs:
            r = p.add_run()
            r.text = txt
            font(r, size - 1 if kind == "b2" else (size - 1 if kind == "h" else size), b, c)
        bullet(p, kind in ("b", "b2"), 1 if kind == "b2" else 0)
        p.alignment = PP_ALIGN.LEFT
        p._p.get_or_add_pPr().set("eaLnBrk", "0")      # never break a Korean word in the middle
        p.space_after = Pt(5 if kind != "h" else 2)
        p.line_spacing = 1.08
    return shp


def title(slide, text, sub=None):
    textbox(slide, 0.45, 0.28, 12.4, 0.7, [([(text, True, INK)], "p")], size=26, name="Title")
    if sub:
        textbox(slide, 0.47, 0.92, 12.4, 0.4, [([(sub, False, SUB)], "p")], size=13, name="Subtitle")
    ln = slide.shapes.add_connector(1, Inches(0.5), Inches(1.33), Inches(W - 0.5), Inches(1.33))
    ln.line.color.rgb = RGBColor.from_string(LINE)
    ln.line.width = Pt(1)


def footer(slide, n, date):
    textbox(slide, 0.45, 7.1, 10.5, 0.3, [([(f"탐색 분석 (exploratory) · 1–5일 · 8마리 · {date[:4]}-{date[4:6]}-{date[6:]}"
                                              f" · analysis/exploratory/{date}", False, SUB)], "p")], size=9,
            name="Footer")
    textbox(slide, 12.2, 7.1, 0.7, 0.3, [([(str(n), False, SUB)], "p")], size=9, name="Page")


def picture(slide, img: Path, x, y, max_w, max_h):
    w, h = Image.open(img).size
    s = min(max_w / (w / 100), max_h / (h / 100))
    ww, hh = w / 100 * s, h / 100 * s
    pic = slide.shapes.add_picture(str(img), Inches(x + (max_w - ww) / 2), Inches(y), Inches(ww), Inches(hh))
    pic.name = img.stem
    return ww, hh


def wide_slide(prs, img, ttl, sub, read, concl, n, date, fig_h=3.25, size=13):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    title(s, ttl, sub)
    ww, hh = picture(s, img, 0.5, 1.45, W - 1.0, fig_h)
    y = 1.45 + hh + 0.12
    hb = 7.0 - y
    textbox(s, 0.5, y, 6.05, hb, [("그림 읽는 법", "h")] + read, size=size, name="How to read")
    textbox(s, 6.75, y, 6.08, hb, [("결론", "h")] + concl, size=size, fill=BOX, name="Conclusion", pad=0.14)
    footer(s, n, date)
    return s


def side_slide(prs, img, ttl, sub, read, concl, n, date, fig_w=7.6, read_h=2.6, size=14, concl_h=None):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    title(s, ttl, sub)
    ww, hh = picture(s, img, 0.45, 1.5, fig_w, 5.5)
    x = 0.45 + fig_w + 0.25
    w = W - x - 0.45
    textbox(s, x, 1.5, w, read_h, [("그림 읽는 법", "h")] + read, size=size, name="How to read")
    yc = 1.5 + read_h + 0.12
    textbox(s, x, yc, w, concl_h or (7.0 - yc), [("결론", "h")] + concl, size=size, fill=BOX, name="Conclusion",
            pad=0.14)
    footer(s, n, date)
    return s


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="20261119")
    a = ap.parse_args(argv)
    D = P / "analysis" / "exploratory" / a.date
    R = json.loads((D / "results.json").read_text(encoding="utf-8"))["results"]
    e1, e2, e3, e4, e5, e6, e7 = (R[k] for k in ("E1", "E2", "E3", "E4", "E5", "E6", "E7"))
    icc = {k.split(" [")[0]: v for k, v in e5.items() if k.endswith("[all]")}
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(W), Inches(H)
    date = a.date

    # 1 title
    s = prs.slides.add_slide(prs.slide_layouts[6])
    textbox(s, 0.8, 2.2, 11.7, 1.1, [([("1–5일 자료로 시험한 탐색 분석", True, INK)], "p")], size=38, name="Title")
    textbox(s, 0.82, 3.35, 11.7, 0.6, [([("새 그림 7개를 쉽게 설명합니다", False, ACC)], "p")], size=22,
            name="Subtitle")
    textbox(s, 0.82, 4.35, 11.7, 1.2, [
        ("모르핀 전 (1–5일), 8마리 (9669m 4 · 9670f 4), 매일 저장한 자료만 사용", "b"),
        ("모두 탐색 분석입니다: 확정 결론이 아니라, 18일 뒤 분석을 어떻게 할지 정하는 데 씁니다", "b")], size=16,
        name="Note")
    footer(s, 1, date)

    # 2 summary cards
    s = prs.slides.add_slide(prs.slide_layouts[6])
    title(s, "한눈에 보기", "숫자 6개로 요약")
    cards = [
        (f"{e1['rule_holds']}/{e1['n_sessions']}", "세션에서 PR 이 lick 수로 정확히 정해짐", "E1"),
        (f"{e1['censored']}/{e1['n_sessions']}", "세션이 15분 10초 끝까지 핥고 있었음", "E1"),
        (f"{e3['rho_clusters_rewards']:.2f} vs {e3['rho_cluster_size_rewards']:.2f}",
         "보상과의 상관: lick 묶음 개수 vs 묶음 크기", "E3"),
        (f"{100 * (1 - e7['median_fraction_left']):.0f} %", "보상 뒤 동공 확장 중 lick 이 설명하는 부분", "E7"),
        (f"{icc['rewards']['icc31']:.2f}", f"보상 수의 날짜 간 재현성 (동공 크기는 {icc['pupil_median_px']['icc31']:.2f})",
         "E5"),
        (f"{e4['rho_strict_vs_TI']:.2f}", "hot plate 와 꼬리 담그기의 상관 (서로 반대)", "E4"),
    ]
    for i, (big, lbl, tag) in enumerate(cards):
        cx, cy = 0.55 + (i % 3) * 4.1, 1.65 + (i // 3) * 2.65
        box = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(cx), Inches(cy), Inches(3.85), Inches(2.4))
        box.fill.solid()
        box.fill.fore_color.rgb = RGBColor.from_string(BOX)
        box.line.fill.background()
        box.shadow.inherit = False
        box.adjustments[0] = 0.08
        box.name = f"Card {i + 1}"
        textbox(s, cx + 0.15, cy + 0.2, 3.55, 0.95, [([(big, True, ACC)], "p")], size=34, name=f"Card {i + 1} number")
        textbox(s, cx + 0.15, cy + 1.15, 3.55, 0.9, [([(lbl, False, INK)], "p")], size=14, name=f"Card {i + 1} label")
        textbox(s, cx + 0.15, cy + 1.95, 3.55, 0.35, [([(tag, False, SUB)], "p")], size=10, name=f"Card {i + 1} tag")
    footer(s, 2, date)

    # 3 E1
    wide_slide(prs, D / "E1_pr_licks" / "E1_pr_vs_licks.png", "PR 점수는 lick 수로 정해집니다 (E1)",
               "3–5일 PR 세션 24개", [
                   ("왼쪽: 점 하나 = 한 세션. 가로 = 센 lick 수, 세로 = 보상 수", "b"),
                   ("회색 계단 = 규칙 \"보상 n 개에는 lick 1+2+…+n 개\"", "b"),
                   ("오른쪽: 쥐가 그만둔 시각. 세모 = 15분 10초 끝까지 핥고 있었음", "b")], [
                   (f"{e1['rule_holds']}/{e1['n_sessions']} 세션이 계단 위 → PR 과 lick 수는 같은 정보", "b"),
                   ("EFA 에 둘 다 넣으면 같은 것을 두 번 셈", "b2"),
                   (f"{e1['censored']}/{e1['n_sessions']} 세션은 끝까지 핥음, 그만둔 쥐도 11–14.6분", "b"),
                   ("지금 PR = '15분 동안 얼마나 핥나' (포기 시점 아님)", "b2")], 3, date)

    # 4 E3
    wide_slide(prs, D / "E3_lick_clusters" / "E3_lick_clusters.png",
               "보상은 lick '묶음의 개수' 를 따라갑니다 (E3)", "묶음 = 0.5초 안에 이어진 lick (Dwyer 2012)", [
                   ("묶음 개수 = 얼마나 자주 핥기 시작하나 ('찾기')", "b"),
                   ("묶음 크기 = 한 번에 몇 번 핥나 ('좋아함')", "b"),
                   ("왼쪽: 날짜별 개수 · 가운데: 날짜별 크기", "b"),
                   ("빨강 = 6일부터 Active, 파랑 = 6일부터 Passive", "b2"),
                   ("오른쪽: 보상 vs 개수 (점) · 크기 (×)", "b")], [
                   (f"보상과의 상관: 개수 {e3['rho_clusters_rewards']:.2f}, 크기 {e3['rho_cluster_size_rewards']:.2f}"
                    " (관련 없음)", "b"),
                   ("쥐 차이는 '좋아함' 이 아니라 '얼마나 자주 시작하나'", "b"),
                   ("6–10일에 모르핀이 어느 쪽을 바꾸는지 보면 됨", "b")], 4, date)

    # 5 E7
    side_slide(prs, D / "E7_pupil_glm" / "E7_pupil_glm.png",
               f"보상 뒤 동공 확장의 약 {100 * (1 - e7['median_fraction_left']):.0f} % 는 lick 때문입니다 (E7)",
               "3–5일, 동공 6마리 (whole eye 2마리는 따로)", [
                   ("a: 보상 시각에 맞춘 평범한 동공 평균 (쥐마다 선 하나)", "b"),
                   ("b: 모델로 나눈 '보상만의' 반응. 회색 = lick 을 모델에 안 넣음, 빨강 = lick 도 넣음", "b"),
                   ("c: lick 한 번에 따른 동공 반응 · d: 쥐별 비교", "b")], [
                   (f"lick 을 넣으면 보상 반응이 중앙값 {100 * e7['median_fraction_left']:.0f} % 만 남음", "b"),
                   ("Passive 는 lick 이 없음 → Active vs Passive 동공 비교에는 lick 을 꼭 빼야 함", "b"),
                   ("lick 과 보상은 같이 일어나서 나눈 값은 근사치", "b2")], 5, date, fig_w=8.3, read_h=2.75, size=13,
               concl_h=2.45)

    # 6 E5
    side_slide(prs, D / "E5_reliability" / "E5_reliability.png", "날마다 얼마나 같은가: 재현성 (E5)",
               "3–5일, 8마리 (동공 6마리)", [
                   ("점 = ICC: 1 = 매일 같은 순서, 0 = 매일 바뀜", "b"),
                   ("막대 = 95 % 범위 (쥐가 적어서 넓음)", "b"),
                   ("점선: 0.5 = 중간, 0.75 = 좋음", "b")], [
                   (f"보상 {icc['rewards']['icc31']:.2f}, lick {icc['licks']['icc31']:.2f} → 하루 값은 흔들림", "b"),
                   ("기준선은 3–5일 평균으로 (3일 평균이면 약 0.6)", "b2"),
                   (f"동공 크기 {icc['pupil_median_px']['icc31']:.2f} 높음 (카메라·눈 위치 영향일 수도),"
                    f" 보상 동공 반응 {icc['pupil_reward_response_pct']['icc31']:.2f} 불안정", "b"),
                   ("이전 덱의 \"기준선이 높을수록 변화가 작다 (−0.7 ~ −0.95)\" 는 평균으로의 회귀만으로도 생길 수 있음"
                    " → 다른 방법으로 다시 확인", "b")], 6, date, fig_w=6.6, read_h=1.95, size=14, concl_h=3.45)

    # 7 E4
    strict = sorted(v for v in e4["latency_strict_s"].values() if v is not None)
    side_slide(prs, D / "E4_hp_latency" / "E4_hp_latency_vs_ti.png", "두 통증 검사가 반대 방향으로 나왔습니다 (E4)",
               "5일 (모르핀 전), 8마리", [
                   ("점 하나 = 쥐 하나", "b"),
                   ("가로 = 꼬리 담그기 지연, 세로 = hot plate 첫 반응 (떨기 · 뒷발 핥기 · 점프) 까지 시간", "b"),
                   ("빨강 = Active, 파랑 = Passive (6일부터)", "b")], [
                   (f"8마리 모두 30초 안에 반응 ({strict[0]:.1f}–{strict[-1]:.1f}초)", "b"),
                   (f"상관 {e4['rho_strict_vs_TI']:.2f} (p {e4['p_strict']:.3f}): 꼬리가 늦은 쥐가 hot plate 는 빨랐음",
                    "b"),
                   ("두 검사가 같은 것을 재지 않을 수 있음 → 금단 통증 점수에서는 따로 봄", "b"),
                   ("8마리라 우연일 수도 있음", "b2")], 7, date, fig_w=7.2, read_h=2.3, size=14, concl_h=2.75)

    # 8 E6
    wide_slide(prs, D / "E6_engagement_hmm" / "E6_engagement_hmm.png", "참여 상태 모델 (HMM): 아직 쓸모 없음 (E6)",
               "1초 단위 lick 수에 2-상태 모델 (Ashwood 2022 방식)", [
                   ("위: 한 세션의 1초당 lick (회색) 과 '참여 중' 확률 (주황)", "b"),
                   ("왼쪽 아래: 참여 시간 vs 보상", "b"),
                   ("오른쪽 아래: 참여 상태가 이어지는 시간", "b")], [
                   (f"1초 단위로는 몇 초짜리 lick 묶음 (평균 {e6['median_dwell_engaged_s']:.1f}초) 만 잡힘", "b"),
                   ("몇 분짜리 참여 기간은 못 잡음", "b2"),
                   (f"보상과의 관계: 30초 쉼 기준 {e6['rho_not_paused30_rewards']:.2f} > 모델 "
                    f"{e6['rho_engaged_min_rewards']:.2f}", "b"),
                   ("10초 단위나 3상태로 다시 시험할 것", "b")], 8, date)

    # 9 E2
    wide_slide(prs, D / "E2_time_of_day" / "E2_time_of_day.png", "시간대: 성별과 섞여 있어 나눌 수 없습니다 (E2)",
               "3–5일", [
                   ("왼쪽: 세션 시작 시각 vs 보상 (선 = 쥐 하나의 3일)", "b"),
                   ("오른쪽: 쥐마다 평소보다 늦게 / 일찍 시작한 정도 vs 보상", "b")], [
                   (f"9669m (수컷) 이 {e2['cage_mean_start_h']['9669m'] - e2['cage_mean_start_h']['9670f']:.1f}시간"
                    f" 늦게 돔 · 보상은 비슷 ({e2['cage_mean_rewards']['9669m']} vs {e2['cage_mean_rewards']['9670f']})", "b"),
                   (f"쥐 안에서 시작 시각과 보상은 관련 없음 ({e2['within_mouse_rho_time_rewards']:.2f})", "b"),
                   (f"동공은 늦은 날 컸지만 ({e2['within_mouse_rho_time_pupil']:.2f}) 날짜 효과와 섞임 (4일이 가장 늦음)",
                    "b")], 9, date)

    # 10 next
    s = prs.slides.add_slide(prs.slide_layouts[6])
    title(s, "그래서 18일 이후 분석에서 바꿀 것", "모두 계획 (§1–§5) 을 먼저 한 다음에")
    rows = [("기준선은 3–5일 평균 (하루 값 쓰지 않음)", "E5"),
            ("PR 과 lick 수는 하나로 묶고, '그만둔 시각' 을 함께 봄", "E1"),
            ("동공 보상 비교 (Active vs Passive) 에는 lick 을 꼭 모델에 넣음", "E7"),
            ("lick 묶음 개수 (찾기) 와 크기 (좋아함) 를 6–10일에 따로 봄", "E3"),
            ("통증 검사 두 개 (꼬리 담그기, hot plate) 는 따로 봄", "E4"),
            ("참여 상태 모델은 10초 단위나 3상태로 다시 시험", "E6"),
            ("성별 결과는 시간대와 겹치므로 설명용으로만", "E2")]
    for i, (txt, tag) in enumerate(rows):
        y = 1.6 + i * 0.75
        chip = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(0.6), Inches(y), Inches(0.85), Inches(0.52))
        chip.fill.solid()
        chip.fill.fore_color.rgb = RGBColor.from_string(ACC)
        chip.line.fill.background()
        chip.shadow.inherit = False
        chip.name = f"Tag {tag}"
        tf = chip.text_frame
        tf.vertical_anchor = MSO_ANCHOR.MIDDLE
        r = tf.paragraphs[0].add_run()
        r.text = tag
        tf.paragraphs[0].alignment = PP_ALIGN.CENTER
        font(r, 14, True, "FFFFFF")
        textbox(s, 1.65, y - 0.02, 11.0, 0.6, [([(txt, False, INK)], "p")], size=17, anchor=MSO_ANCHOR.MIDDLE,
                name=f"Change {tag}")
    footer(s, 10, date)

    out = P / "analysis" / "ppt" / f"exploratory_days1-5_{date}.pptx"
    out.parent.mkdir(parents=True, exist_ok=True)
    prs.save(out)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
