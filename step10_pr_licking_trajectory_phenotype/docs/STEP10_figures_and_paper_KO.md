# Step 10 그림 설명 (한국어) + 논문(레버 자가투여 궤적)과의 비교

## 실험 설계와 `day_index >= 4`

- **Day 3**: PR(Progressive Ratio) 절차가 **시작**되는 날로 가정.
- **Step 10**에서는 **day_index ≥ 4**인 행만 사용합니다.  
  → PR이 처음 도입된 day 3의 “적응/시작” 날을 빼고, **PR이 본격적으로 쌓인 뒤**의 궤적만 보려는 목적에 맞습니다.
- 스크립트 상수: `run_step10_pr_lick_trajectory_phenotype.m` 안의 `DAY_INDEX_MIN = 4` (필요 시 숫자만 조정).

---

## 출력 그림별 의미 (out_step10/)

### 1. `pr_feature_correlation_heatmap.png`

| 항목 | 설명 (한국어) |
|------|----------------|
| 무엇을 그리나 | **쥐×날짜** 단위로 집계한 세 변수 — **PR breakpoint(그날 도달한 최대 PR trial)**, **within-trial 핥기율(Hz)**, **between-trial 핥기율(Hz)** — 사이의 **피어슨 상관**을 색으로 표시합니다. |
| 어떻게 읽나 | 색이 **진할수록** (±에 가깝게) 두 변수가 같이 움직입니다. **0.8 이상**이면 (스크립트가 콘솔에 경고) 둘 중 **하나만** 군집에 써도 될 수 있습니다. |
| CSV | `pr_feature_corr_mouse_day.csv` 동일 수치. |

**논문과 비교**: 약물 자가투여 논문은 보통 **active lever, futile press, inactive** 등 **서로 다른 construct**의 상관·공선성을 봅니다. 여기서는 **lever 대신** PR breakpoint + licking rate들의 상관을 봅니다. **같은 정보가 중복**되지 않게 하려는 목적은 동일합니다.

---

### 2. `pr_cluster_mean_trajectory_shape.png`

| 항목 | 설명 (한국어) |
|------|----------------|
| 무엇을 그리나 | k-means로 나눈 **군집별**로, (z-score **이전**의) **특징 벡터 평균 ± SE**를 **에러바 선**으로 그립니다. |
| 가로축 | **특징 인덱스**: (1) **3개만** 쓰는 경우 → 대략 [max PR, mean within, mean between]. **여러 날을 이어 붙인 궤적** 모드면 `특징1×날1, …, 특징3×날D` 순으로 인덱스가 길어집니다. |
| 해석 | 군집마다 **동기(motivation)**·**섭취(consummatory)**·**지속 탐색(persistent seeking)** 쪽 패턴이 어떻게 다른지 **요약 실루엣**으로 봅니다. |

**논문과 비교**: GBTM / “escalation trajectory”류 논문은 **시간에 따른 군집 평균 궤적**을 제시하는 경우가 많습니다. 여기서는 **PR/핥기 지표**로 같은 **“서로 다른 longitudinal 패턴”**을 **요약 plot**으로 맞춘 것입니다. (정확한 GBTM은 별도 모델.)

---

### 3. `pr_breakpoint_spaghetti_by_cluster.png`  *(day-궤적 모드에서만)*

| 항목 | 설명 (한국어) |
|------|----------------|
| 무엇을 그리나 | **세로**: 그날의 **PR breakpoint** (해당 day×session에서 관측된 `Trial`의 최댓값 = 그날 PR을 얼마나 밀었는지). **가로**: `day_index` (**4 이상**만 포함). **선**: 쥐 한 마리씩. **색**: Step10에서 부여한 **궤적 표현형 군집**. |
| 해석 | “시간이 지나면서 **PR 동기(노력)**가 어떻게 변하는지” **쥐마다** 보고, **같은 군집**끼리 비슷한 **상승/하강/안정** 패턴이 있는지 눈으로 확인합니다. |

**논문과 비교**: lever-press 논문의 **“drug intake / escalation over time”** 곡선과 **역할이 비슷**합니다. 다만 **약물 주입 횟수** 대신 **PR에서 요구한 난이도( trial 수 )**로 **effort-based motivation**을 봅니다.

---

## CSV

| 파일 | 설명 |
|------|------|
| `pr_licking_mouse_by_day.csv` | (day≥4) **쥐×날** 요약: breakpoint, within, between. |
| `pr_licking_trajectory_clusters.csv` | 쥐별 **궤적 표현형 군집** ID. |
| `pr_feature_corr_mouse_day.csv` | 위 3변수 상관행렬. |

---

## 한 문단 Methods (논문용 영어 예시, 수정 가능)

Because the paradigm does not use lever-press self-administration, we adapted trajectory-based phenotyping to PR-derived licking measures. From day 4 onward (day 3 = first PR), we summarized each mouse and day by PR breakpoint (maximum PR trial), within-trial lick rate, and between-trial lick rate, then used k-means on z-scored feature vectors to identify subgroups with distinct longitudinal patterns of motivation-related and persistent seeking-related licking.

---

## 요약: 논문 vs 이 분석

| | **전형적 논문 (lever SA)** | **Step 10 (네 데이터)** |
|---|---------------------------|-------------------------|
| 핵심 질문 | 누가 **약물 / infusion** 횡적으로 올리는가 | 누가 **PR/핥기 “동기–섭취–지속탐색”** 패턴이 다른가 |
| 측정 | reinforced press, futile press, … | **PR breakpoint**, **within / between lick rate** (day≥4) |
| 군집 | trajectory / GBTM 등 | z-score 후 **k-means** (가벼운 벤치마크) |

`DAY_INDEX_MIN`이나 변수 목록을 바꾸면, **같은 “궤적 phenotyping” 틀**을 유지한 채 **다른 behavioral stream**에도 쉽게 맞출 수 있습니다.
