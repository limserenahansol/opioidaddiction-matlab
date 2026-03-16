# Push this repo to GitHub

Follow these steps to publish **opioidaddiction-matlab** to GitHub.

## 1. Create the repo on GitHub

1. Go to [github.com/new](https://github.com/new).
2. **Repository name:** e.g. `opioidaddiction-matlab`.
3. **Description (optional):** e.g. "MATLAB pipeline for morphine PR: per-session (pupil, assays) → longitudinal → downstream (motivation, lick, QC, EFA, modules 5–12)".
4. Choose **Public** (or Private).
5. **Do not** check "Add a README", "Add .gitignore", or "Choose a license" — this repo already has them.
6. Click **Create repository**.

## 2. Push from your machine

In a terminal, from the **opioidaddiction-matlab** folder:

```bash
cd C:\Users\hsollim\behavior_task\opioidaddiction-matlab

# Add GitHub as remote (replace YOUR_USERNAME with your GitHub username)
git remote add origin https://github.com/YOUR_USERNAME/opioidaddiction-matlab.git

# Use main as default branch (optional)
git branch -M main

# Push
git push -u origin main
```

If GitHub asks for a password, use a **Personal Access Token** (Settings → Developer settings → Personal access tokens), not your account password.

## 3. Check the repo on GitHub

- The **README** will show the **pipeline schematic** (Phase 1: per-session → Phase 2: build longitudinal → Phase 3: downstream).
- **PIPELINE.md** has the flow diagram and the advanced pipeline (Plan / Module → MATLAB) table.

Done. Others can clone with:
`git clone https://github.com/YOUR_USERNAME/opioidaddiction-matlab.git`
