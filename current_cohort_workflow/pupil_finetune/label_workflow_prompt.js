// Labelling workflow used on 2026-11-15 (Claude Code Workflow tool): 2 independent labels per frame + adjudication.
// Set DIR and the session list S for a new cohort. Output: votes + adjudications -> fit_pupil_rule3.py / finetune_unet_cv.m
export const meta = {
  name: 'pupil-candidate-labels-v2',
  description: 'Label which numbered contour is the pupil (16 sessions x 24 frames, 2 votes each) and adjudicate disagreements',
  phases: [
    { title: 'Label', detail: '8 labellers, 4 sessions each; every session labelled twice by different labellers' },
    { title: 'Adjudicate', detail: 'a third look at frames where the two labels disagree' },
  ],
}
const DIR = '<FOLDER WITH THE CANDIDATE SHEETS>'   // set before running
const S = ['day1_9669m_black','day1_9669m_green1','day1_9669m_green2','day1_9669m_white','day1_9670f_black','day1_9670f_green1','day1_9670f_green2','day1_9670f_white',
           'day2_9669m_black','day2_9669m_green1','day2_9669m_green2','day2_9669m_white','day2_9670f_black','day2_9670f_green1','day2_9670f_green2','day2_9670f_white']
const SHEETS = ['A','B','C','D','E','F']
const groups = [[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15],[0,4,8,12],[1,5,9,13],[2,6,10,14],[3,7,11,15]]
const FRAME = { type: 'object', properties: {
  frame: { type: 'integer' },
  choice: { type: 'integer', description: 'number of the contour that outlines the pupil; 0 = pupil visible but no contour outlines it; -1 = pupil cannot be seen' },
  fit: { type: 'string', enum: ['good', 'too_big', 'too_small', 'na'] },
  confidence: { type: 'string', enum: ['high', 'medium', 'low'] },
  note: { type: 'string' } }, required: ['frame', 'choice', 'fit', 'confidence'] }
const SCHEMA = { type: 'object', properties: { sessions: { type: 'array', items: { type: 'object', properties: {
  session: { type: 'string' }, eye: { type: 'string', description: 'one line: what this eye looks like (open/narrowed/inflamed, iris grey or dark, pupil visible?)' },
  frames: { type: 'array', items: FRAME } }, required: ['session', 'frames'] } } }, required: ['sessions'] }
const ANATOMY = `Anatomy: the eye opening is an almond- or dome-shaped area, usually bordered by a dark rim (eyelid margin). Inside it the iris looks grey (in some mice dark grey) and the PUPIL is a round or oval dark/black disc inside the iris. Its size changes between frames (small = constricted, large = dilated) and it moves with the eye. Bright white spots are reflections of the infrared lights on the cornea (glints); they often lie on or next to the pupil and can cut a notch into it. Dark areas that are NOT the pupil: the dark eyelid rim along the edge of the eye opening, the corners of the eye opening (a dark wedge at its left or right end, often next to a glint), fur, shadows, whiskers. A pupil is usually bordered by grey iris on most of its sides. Some mice have eye inflammation: the eye can be narrowed to a slit, swollen or cloudy. If you cannot see the pupil, answer -1; if you can make a reasonable guess, give it with low confidence.`
const labelPrompt = (keys) => `You are labelling infrared video frames of mouse eyes for a pupil tracker (side view, head-fixed mouse).

Sheets (PNG). Read each one ONCE with the Read tool. Do not run code, do not zoom, do not use any tool other than Read.
${keys.map(k => SHEETS.map(s => `- ${DIR}\\${k}__cand${s}.png`).join('\n')).join('\n')}

Each session has 6 sheets (A-F) with 4 frames each (frames 0-23; the sheet title gives the range). Layout: 2 rows x 4 tiles; every frame is shown twice side by side: LEFT tile = the eye crop (contrast stretched on the dark range, enlarged), RIGHT tile = the same crop with numbered coloured contours. Each contour is one candidate region an algorithm found; its number is in the small black box next to it, written in the contour's colour. Tile titles give the frame number.

${ANATOMY}

For EVERY frame (24 per session) answer:
- choice: the number of the contour that outlines the pupil (the best one if several overlap it). 0 = the pupil is visible but no contour outlines it. -1 = the pupil cannot be seen (eye closed or almost closed, blurred, not an eye, too dark to tell).
- fit: good (the contour roughly follows the pupil boundary), too_big (contains the pupil but clearly more, e.g. merged with a dark eye corner or the eyelid rim), too_small (only part of the pupil), na (choice 0 or -1).
- confidence: high / medium / low.
- note: a few words when useful.
Also give one line per session describing the eye (field "eye"). Judge each frame on its own, but you may compare frames of the same mouse to learn what its pupil looks like. Be careful and honest; do not pick a contour just because it is dark. Return all ${keys.length * 24} frames.`

phase('Label')
const votes = (await parallel(groups.map((g, i) => () =>
  agent(labelPrompt(g.map(j => S[j])), { label: `label:${i < 4 ? 'v1' : 'v2'}-${i % 4}`, phase: 'Label', schema: SCHEMA })
    .then(r => r ? ({ labeller: i, vote: i < 4 ? 1 : 2, keys: g.map(j => S[j]), result: r }) : null)))).filter(Boolean)
log(`labellers returned: ${votes.length}/8`)

// collect per session/frame choices
const byKey = {}
for (const v of votes) for (const s of (v.result.sessions || [])) {
  const key = S.find(k => (s.session || '').includes(k)); if (!key) continue
  for (const f of (s.frames || [])) {
    const id = `${key}|${f.frame}`; (byKey[id] = byKey[id] || []).push({ labeller: v.labeller, ...f })
  }
}
const disputed = Object.entries(byKey).filter(([id, a]) => a.length >= 2 && new Set(a.map(x => x.choice)).size > 1)
const single = Object.entries(byKey).filter(([id, a]) => a.length === 1).map(([id]) => id)
log(`frames with 2 labels: ${Object.values(byKey).filter(a => a.length >= 2).length}; disputed: ${disputed.length}; only one label: ${single.length}`)

phase('Adjudicate')
const bySession = {}
for (const [id, a] of disputed) { const [k, f] = id.split('|'); (bySession[k] = bySession[k] || []).push({ frame: +f, votes: a }) }
const sessKeys = Object.keys(bySession)
const chunks = [[], [], [], []]
sessKeys.forEach((k, i) => chunks[i % 4].push(k))
const ADJ = { type: 'object', properties: { decisions: { type: 'array', items: { type: 'object', properties: {
  session: { type: 'string' }, frame: { type: 'integer' }, choice: { type: 'integer' }, fit: { type: 'string', enum: ['good', 'too_big', 'too_small', 'na'] },
  confidence: { type: 'string', enum: ['high', 'medium', 'low'] }, note: { type: 'string' } }, required: ['session', 'frame', 'choice', 'confidence'] } } }, required: ['decisions'] }
const adjPrompt = (keys) => `Two independent labellers disagreed about which numbered contour outlines the PUPIL in some infrared mouse-eye frames. Decide each disputed frame yourself.

${ANATOMY}

Sheet layout: each session has sheets ${DIR}\\<session>__candA.png ... __candF.png with 4 frames each (A = frames 0-3, B = 4-7, C = 8-11, D = 12-15, E = 16-19, F = 20-23). 2 rows x 4 tiles; each frame twice: LEFT = eye crop, RIGHT = the same with numbered candidate contours (number in a small black box in the contour's colour). Read only the sheets you need, once each; do not run code and do not use any tool other than Read.

Disputed frames (with the two labels: choice, fit, confidence, note; choice 0 = pupil visible but not outlined, -1 = pupil not visible):
${keys.map(k => bySession[k].map(d => `- ${k} frame ${d.frame} (sheet ${SHEETS[Math.floor(d.frame / 4)]}): ` + d.votes.map(v => `[choice ${v.choice}, ${v.fit}, ${v.confidence}${v.note ? ', "' + v.note + '"' : ''}]`).join(' vs ')).join('\n')).join('\n')}

For every disputed frame give your own choice (a contour number, 0 or -1), fit and confidence. Do not just pick one of the two labels by default; look at the frame.`
const adj = (await parallel(chunks.filter(c => c.length).map((c, i) => () =>
  agent(adjPrompt(c), { label: `adjudicate-${i}`, phase: 'Adjudicate', schema: ADJ })))).filter(Boolean)
return { votes, adjudications: adj.flatMap(a => a.decisions || []), disputed: disputed.length, single }
