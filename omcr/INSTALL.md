# OMCR — opioidaddiction-matlab

Same plugin as wiki repo. **Install once** in Claude Code globally.

## Install

```
/plugin marketplace add https://github.com/youngeun1209/oh-my-claudecode-research
```

```
/plugin install oh-my-claudecode-research
```

## This repo

```bash
cd opioidaddiction-matlab
# Claude Code → /omcr-setup → /start-research
```

Root `CLAUDE.md` is pre-filled for the MATLAB pipeline.

## Preset memories

Copy from sibling repo:

```bash
git clone https://github.com/limserenahansol/opioid-model_wiki.git /tmp/wiki
for agent in supervisor analysis-implementer paper-writer figure-descriptor reviewer literature-curator; do
  cp /tmp/wiki/examples/opioid-pr-behavior/memory-templates/$agent/MEMORY.md \
     .claude/agent-memory/$agent/MEMORY.md
done
```

Or run `/start-research` with preset hint `opioid-pr-behavior` if supported.

## Two-repo workflow

| Repo | You use OMCR for |
|------|------------------|
| **This repo** | Run step01–07, QC, figures |
| **opioid-model_wiki** | M0–M1 fit, equations, Eli collaboration |

Full guide: https://github.com/limserenahansol/opioid-model_wiki/blob/main/omcr/INSTALL.md

Korean: https://github.com/limserenahansol/opioid-model_wiki/blob/main/omcr/INSTALL.ko.md
