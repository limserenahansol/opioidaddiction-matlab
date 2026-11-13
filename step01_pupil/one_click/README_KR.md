# Pupil one-click — 공용 U-Net 하나 + 한 번 실행 (2026-11-12)

English: [README.md](README.md)

예전 Step 01a(`step1_2dec29_nove222222222.m` → `step3_dec291123342.m` → `alignment_pupil_salae_aug20.m`)는 **영상마다** 세 가지를 반복해야 했습니다.
- 경로를 직접 고치기 (3개 스크립트)
- 프레임 2–3장에 클릭하기
- U-Net을 새로 학습하기

Dec-2025 cohort에서는 학습을 394번 했습니다. 이 폴더는 그 학습들이 남긴 라벨 약 20만 장으로 **모든 마우스·모든 날에 쓰는 U-Net 하나**를 만들고, 세션 찾기 → 추적 → TTL 정렬을 한 번에 실행합니다. 원래 스크립트 3개는 바꾸지 않았습니다.

## 실행 (MATLAB)

```matlab
% 처음 한 번: one_click_paths.example.json을 one_click_paths.json으로 복사하고 폴더를 적기 (git에는 안 올라감)
run_pupil_one_click('DryRun', true)   % 찾은 세션, TTL, 고른 영상, ROI만 확인 (아무것도 쓰지 않음)
run_pupil_one_click                   % 전부 실행 (끝난 세션은 건너뜀)
run_pupil_one_click('VideoRoot', V, 'SessionRoot', S, 'OutRoot', O)   % 경로를 직접 줘도 됨
run_pupil_one_click('Days', [12 13], 'Mice', {'0910_red'})          % 일부만
run_pupil_one_click('ROIMode', 'previous')                          % ROI가 없으면 묻지 않고 전날 ROI 사용
```

**손으로 하는 일은 ROI뿐입니다.** `roi.mat`이 없는 세션만 실행을 시작할 때 한꺼번에 묻습니다.
- 전날 ROI가 움직일 수 있는 상자로 떠 있으니, 눈 위에 맞추고 더블클릭하면 됩니다.
- 창을 닫으면 전날 ROI를 그대로 씁니다.
- 그 뒤로는 자동으로 돌아갑니다.

**필요한 것:** MATLAB R2024b(시험한 버전) + Deep Learning, Image Processing, Computer Vision Toolbox. U-Net 층이 Computer Vision Toolbox에 들어 있어서 모델을 불러올 때도 필요합니다. Signal Processing은 선택입니다. GPU가 있으면 세션당 약 3분 걸리고, 없으면 CPU로 더 느리게 돌아갑니다.

## 폴더 구조 (Dec-2025와 같음)

```
VideoRoot\dayN\<cage>\<color>\*.avi  (+ roi.mat)     카메라 2대 영상이 있어도 됨
SessionRoot\dayN\<cage>\<color>\digital.csv           Saleae export (Ch6 camera, Ch7 lick, Ch4 injector)
```
`digital.csv`가 `SessionRoot`에 없으면 `VideoRoot\dayN\` 아래 어디든 `...\<cage>\<color>\digital.csv`를 찾습니다. Dec-2025의 TTL 폴더 이름은 `ttl`, `tt`, `day4Jan2` 등으로 날마다 달랐습니다.

## 세션마다 하는 일

1. `digital.csv`에서 카메라 TTL rising edge 수를 셉니다.
2. **pupil 영상:** TTL 수와 프레임 수가 가장 가까운 `.avi`를 고릅니다. Dec-2025 세션 244개 중 243개에서 맞았고, 나머지 1개는 두 영상 길이가 같아서 경고로 표시됩니다.
3. ROI는 `roi.mat`을 씁니다.
4. **추적 (`pupil_track_session.m`):** 공용 U-Net을 적용합니다.
   - 전처리는 step 1이 학습 이미지를 만든 방식과 같습니다: ROI 자르기 → 회색조 → glint 제거 → CLAHE.
   - 그 다음은 **step 3의 추적 규칙을 그대로** 씁니다.
   - 결과 파일 이름과 열도 예전과 같고, `qc_frames.png`와 `qc.json`이 추가됩니다.
   - **예전과 다른 점:** 예전 step 3은 첫 프레임을 건너뛰었습니다. 새 코드는 모든 프레임을 처리해서 행 수 = 영상 프레임 = 카메라 TTL이 되고, 정렬이 `start_1to1`로 맞습니다.
5. **정렬 (`pupil_align_ttl.m`):** `alignment_pupil_salae_aug20.m`의 로직을 **그대로 복사**한 함수입니다. 결과는 예전과 같은 `concat_out_<시각>\`입니다.
6. **보고서:** `OutRoot\one_click_report_<시각>.csv`에 세션마다 한 줄씩, `flags` 열에 문제를 적습니다.

**주의:** Step 02는 세션 폴더의 `concat_out*`을 하나하나 별도 세션으로 셉니다. 그래서 `OutRoot`는 따로 떨어진 폴더로 두세요. `concat_out`이 이미 있는 폴더에는 두 번째를 만들지 않고 멈춥니다.

## 공용 모델

- `model\general_pupil_unet.mat`은 git에 없습니다 (미발표 영상으로 학습했기 때문). [model/README.md](model/README.md)를 보세요.
- 구성: 128×128 입력, U-Net depth 4, generalized Dice loss, 12 epoch
- 다시 만들기: `inventory_existing_pupil.py` → `make_training_list.py` → MATLAB `pupil_train_general(...)` (GPU에서 약 40분)

## 빠른 시험 (toy 데이터)

`make_toy_dataset.py` → `run_pupil_one_click(...toy 경로..., 'ROIMode','previous')` → `test_toy_compare.py`. 자세한 명령은 영어판에 있습니다.

## 확인한 것 (2026-11-12, toy만 — 전체 데이터로는 돌리지 않음)

- **학습에서 뺀 3마리** (0910_red, 0911_white, 6100_orange, run 64개):
  - Dice 중앙값 **0.952** (IQR 0.872–0.969)
  - 영상마다 학습한 모델이 자기 영상에서 잰 값은 0.962입니다 (실제보다 좋게 나오는 값).
- **toy 원클릭** (세션 3개 × 900 프레임):
  - 3/3 성공, 다시 실행하면 3/3 건너뜀
  - 영상 고르기, TTL 찾기, ROI 재사용, `start_1to1` 정렬 모두 예상대로
  - 세션당 약 3분
- **예전 결과와 비교:** 새 모델이 6–12 % 작습니다.
  - 동공이 줄어든 프레임에서, 새 모델은 홍채 안의 어두운 동공(26–30 px)을 잡고, 예전 모델은 홍채까지 포함한 눈 전체(38–40 px)를 잡습니다.
  - 동공이 클 때는 둘이 같습니다.
  - 그래서 **예전 pupil 데이터가 수축을 작게 잡았을 수 있습니다.** 다만 짧은 toy 세션 3개에서만 본 것입니다.
- **알려진 점:** glint를 지운 자리 때문에 마스크 가장자리가 조금 패일 수 있습니다. 지름 정의는 예전 step 3과 같게 두었습니다.
