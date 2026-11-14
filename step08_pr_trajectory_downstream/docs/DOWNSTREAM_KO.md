# Step 08 다운스트림: 무엇을 추가했는가 (한국어)

## 배경

기존 **opioidaddiction-matlab** 흐름(Step 02 `ALL_mice_longitudinal.csv` → Step 04 lick/PR → Step 07 EFA·고급 통계)은 “요일/페이즈별 PR·행동 지표” 중심입니다.  
여기에 **리뷰 논문(펜타닐 escalation / trajectory phenotyping)**에서 쓰는 아이디어를 **PR-only·licking 기반**으로 옮긴 **추가 다운스트림**을 Step 08로 두었습니다.

## 새로 넣은 접근 (2갈래)

### 1) PR trial별 licking trajectory (필수 플롯)

- **PR1, PR2, …** 각 breakpoint 단계(= trial)마다  
  - **trials 내부** lick rate = *within-trial*  
  - **이전 trial 종료 ~ 다음 trial 시작** 사이 = *between-trial*  
- 논문의 lever reinforced / futile / inactive와 **개념적으로** 대응할 수 있음:  
  - within-trial ≈ 보상 구간의 consummatory vigor  
  - between-trial ≈ 보상이 없는 간격의 persistent / anticipatory licking  

→ 기존 파이프라인 **끝에** “effort demand가 올라갈수록 lick dynamics가 어떻게 휘는지”를 한 눈에 보는 레이어.

### 2) Trajectory phenotype (선택 클러스터)

- 쥐별로 **PR trial × within-rate(between-rate)** 곡선을 벡터로 쌓고, z-score 후 **k-means** (SAS `PROC TRAJ` 대응 **근사**; 표본이 작으면 k=2~3 권장).  
- Step 07의 **EFA addiction index 곡선**과는 달리, “먼저 시간( trial )으로 subgroup을 잡는” 쪽에 가깝습니다. EFA는 “먼저 변수를 요인으로 줄인 뒤 시간”이므로 **둘을 나란히** 두고 overlap을 보는 전략이 맞습니다.

## `run_pr_licking_full_analysis` (권장)

- **`USE_DEMO = false`**: 최신 `run_*/figs/lick_MEGA/per_trial_features.csv` + `ALL_mice_longitudinal.csv`로 **전체** 세션/일을 사용한다.  
- **All** + **Pre / During / Post / Withdrawal / Reexposure** 각각에 대해, **PR trial(= breakpoint 단계)별**  
  - **Within-trial** lick rate (`lick_rate` = time-to-reward 구간)  
  - **Between-trial** lick rate(직전 reward ~ 다음 trial 시작, Step08 enrich)  
  를 그린다.  
- **회색 = 개체별** 곡선(같은 마우스·세션을 trial 평균으로 합친 뒤), **색 = 그룹** 평균±SEM(Active/Passive 등 표에 있는 그대로).  
- `PR_lick_stats_trialwise.csv`: PR 단계마다 **두 그룹** 비교(기본 Active vs Passive가 없으면 **알파벳 순 첫 두 그룹**). `PR_lick_slopes_by_mouse.csv`: 마우스마다 lick~PR 기울기(추세).

## 기존 파이프라인과의 위치

```text
... Step 04 (lick / per_trial_*.csv) ─┐
     Step 07 (EFA, LME, ...)         ─┼─► Step 08: PR lick trajectory + optional cluster
     longitudinal table               ┘
```

- Step 04에서 `Trial`이 있는 per-trial 표가 있으면, 그걸 **정리해** `WithinTrialLickRate` / `BetweenTrialLickRate`를 넣은 CSV로 내보낸 뒤 Step 08에 넣을 수 있습니다.  
- **당장** 돌리려면 `demo_data/demo_pr_trial_table.csv`로 **시뮬**이 돌아갑니다.

## 데이터 요구(요약)

- 필수: `mouse_key`, `day_index`, `session_idx`, `Trial`, `WithinTrialLickRate`  
- 권장: `BetweenTrialLickRate` (없으면 between 플롯은 스킵)  
- 선택: `Group` (집단별 색)  

세부는 `README.md` 참고.
