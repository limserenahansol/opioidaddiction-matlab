# OMCR preset: opioid PR / AMM

Field overlay for [oh-my-claudecode-research](https://github.com/youngeun1209/oh-my-claudecode-research).

## Apply after `/omcr-setup`

```bash
cd opioid-model_wiki   # or opioidaddiction-matlab

for agent in supervisor analysis-implementer paper-writer figure-descriptor reviewer literature-curator; do
  cp examples/opioid-pr-behavior/memory-templates/$agent/MEMORY.md \
     .claude/agent-memory/$agent/MEMORY.md
done
```

Only overwrites if you intend to replace template memories.

## Canonical model

https://github.com/limserenahansol/opioid-model_wiki/tree/main/model
