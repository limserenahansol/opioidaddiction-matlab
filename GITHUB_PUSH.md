# Push this repo to GitHub

## 1. Create the repo on GitHub

1. Go to [github.com/new](https://github.com/new).
2. **Repository name:** e.g. `opioidaddiction-matlab`.
3. **Description (optional):**  
   `MATLAB pipeline for morphine PR experiment: per-session (pupil, assays) → longitudinal (all mice × days) → downstream (motivation, lick, QC, EFA, Straub).`
4. Choose **Public** (or Private).
5. **Do not** check "Add a README", "Add .gitignore", or "Choose a license".
6. Click **Create repository**.

## 2. Push from your machine

```bash
cd C:\Users\hsollim\behavior_task\opioidaddiction-matlab

git remote add origin https://github.com/YOUR_USERNAME/opioidaddiction-matlab.git
git branch -M main
git push -u origin main
```

Replace **YOUR_USERNAME** with your GitHub username (e.g. `limserenahansol`).  
If prompted for a password, use a [Personal Access Token](https://github.com/settings/tokens).

## 3. What people will see

- **README** shows the pipeline schematic image and a short table (Phase 1 → 2 → 3).
- **PIPELINE.md** has the run order and "run one of" summary for the advanced step.
- Each **step folder** has a README with a "run one of" table so users pick one script per purpose.

Clone for others:  
`git clone https://github.com/YOUR_USERNAME/opioidaddiction-matlab.git`
