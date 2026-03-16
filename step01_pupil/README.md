# Step 01: Pupil tracking and alignment (run first)

Run this step **first** in the pipeline. Per-session pupil processing: pseudo-labeling, U-Net tracking, and alignment to Saleae digital TTLs.

## Scripts (run in order for each video/session)

1. **`step1_2dec29_nove222222222.m`** — ROI + point seeds → pseudo-labels → U-Net training data. Edit at top: `BASE`, `VIDEO`, `ROI_FILE`, `RUN_ROOT`.
2. **(Step 2)** — Train U-Net (external or script not in this repo).
3. **`step3_dec291123342.m`** — Load config from Step 1/2 run folder, apply trained U-Net, write pupil CSV + video. Edit: `BASE`, `VIDEO`, `RUN_FOLDER_NAME`, `MODEL_DIR`.
4. **`alignment_pupil_salae_aug20.m`** — Align Saleae `digital.csv` to pupil CSV. Set `DIGITAL_FILE` and `PUPIL_FILE` at top.

## Next

After Step 01 (and 02–03), run **Step 04** to build `longitudinal_outputs` and `ALL_mice_longitudinal.csv`.
