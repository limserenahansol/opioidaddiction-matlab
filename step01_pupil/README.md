# Step 01a: Pupil tracking and alignment (per-session)

**Per-session** work — run for every session/video before building the longitudinal file. Pupil pseudo-labeling, U-Net tracking, and alignment to Saleae digital TTLs.

## Scripts (run in order for each video/session)

1. **`step1_2dec29_nove222222222.m`** — ROI + point seeds → pseudo-labels → U-Net training data. Edit at top: `BASE`, `VIDEO`, `ROI_FILE`, `RUN_ROOT`.
2. **(Step 2)** — Train U-Net (external or script not in this repo).
3. **`step3_dec291123342.m`** — Load config from Step 1/2 run folder, apply trained U-Net, write pupil CSV + video. Edit: `BASE`, `VIDEO`, `RUN_FOLDER_NAME`, `MODEL_DIR`.
4. **`alignment_pupil_salae_aug20.m`** — Align Saleae `digital.csv` to pupil CSV. Set `DIGITAL_FILE` and `PUPIL_FILE` at top.

## One-click alternative (2026-11-12)

**`one_click/run_pupil_one_click.m`** does steps 1–4 for every mouse x day in one call, with one general U-Net (no per-video clicks or
training): finds sessions and TTL files, picks the pupil video by its frame count, tracks, aligns, and writes a QC report. The scripts above
are unchanged. See [one_click/README.md](one_click/README.md).

## Next

After all sessions have pupil (and **Step 01b** manual assays) done, run **Step 02** to build the longitudinal folder (all mice × all days).
