# Panel as a Workflow (token-optimal default)

Copy the script below into the `Workflow` tool's `script` argument, fill the four
slots, and launch. Each expert + the judge runs on `model:'opus', effort:'high'`.
Only the judge's verdict returns to the main conversation, so the debate never
bloats the user's context. Watch progress with `/workflows`.

**Fill these before launching:**
- `QUESTION` — the sharp one/two-line question.
- `BRIEF` — the shared evidence brief. **No conclusion of your own.** Numbers,
  units, the data table, reference values, the user's standing observations.
- `FILE_POINTERS` — absolute paths E4 must read (handover, results CSVs, campaign
  docs, run dirs).
- Adjust `EXPERTS[].identity` only if swapping a lens for a non-jellium question;
  always keep E4.

```javascript
export const meta = {
  name: 'scientific-panel',
  description: 'Four independent domain experts debate a hard result; a judge synthesises a provisional verdict',
  phases: [
    { title: 'Openings', detail: '4 experts give independent positions', model: 'opus' },
    { title: 'Rebuttal', detail: 'each expert sees the others and revises', model: 'opus' },
    { title: 'Judge', detail: 'synthesis verdict for the user', model: 'opus' },
  ],
}

const QUESTION = `<<one or two lines>>`
const BRIEF = `<<evidence brief — no conclusion>>`
const FILE_POINTERS = `<<absolute paths for E4>>`

const PREAMBLE = `You are one of four independent experts deliberating on a hard result from a first-principles TDDFT stopping-power campaign. Give your HONEST, INDEPENDENT read. Your final message IS your position (read by a judge and the other experts, not the user) — return reasoning and numbers, not a polished essay.
HOUSE STYLE: Start from the SIMPLEST model that could explain the observation (free particle, single harmonic mode, linear response / Lindhard, two-level system, dimensional analysis); state the number it predicts; add complexity only when the simple model demonstrably fails, showing each step. Ground claims in known results/limits and cite them. Separate evidence from inference (label inferences "Inference:"). Carry units; round reported numbers to 2 significant figures. Saying "the evidence does not decide this, and here is what would" is encouraged.`

const OPENING_TASK = `YOUR TASK (opening position), concise and structured: (1) the simplest reduced model that captures the phenomenon and the number it predicts; (2) your leading explanation of what is really happening, through your lens; (3) the ONE test that would most cleanly settle it, ideally runnable on existing data or a cheap new run; (4) your confidence and what would change your mind.`

const EXPERTS = [
  { key: 'E1', label: 'TDDFT / stopping-power methodologist', identity: `You compute electronic stopping from real-time TDDFT. Your expertise is WHAT THE CALCULATION ACTUALLY MEASURES and where method artifacts hide: the energy-ledger (deposited = E_total(t_f) - E_GS) vs force- or momentum-based extraction; what an absorbing boundary (CAP) removes vs leaves; self-interaction error of one added electron; the first moment <p> (drift = stopping) vs the second moment Var(p) (spread / zero-point) of a wavepacket; convergence in grid, dt, run length. Your instinct on a surprising number is to ask whether it is the method, not the physics.` },
  { key: 'E2', label: 'condensed-matter generalist', identity: `You are a broad condensed-matter theorist. You see this through many-body response: epsilon(q,omega), the loss function -Im[1/epsilon], the f-sum rule, the plasmon pole, electron-hole pair creation, Pauli/exchange, Fermi-liquid intuition, EELS analogies. You check claimed results against what linear response and the sum rules permit, and notice when a number breaks a limit it cannot break.` },
  { key: 'E3', label: 'jellium-stopping specialist', identity: `You have spent years on electronic stopping in the homogeneous electron gas: Lindhard/RPA S(v), the Bragg peak near v_F, the high-v Bethe ln(v)/v^2 tail, the low-v linear friction Q, nonlinear screening (Echenique-Ritchie-Brandt, DFT-based stopping), Barkas / Z-sign effects, heavy-ion vs light-electron projectiles, and how a finite/localised slab departs from the bulk. You place a result on the known stopping curve and judge it against established limits and density (r_s) scaling.` },
  { key: 'E4', label: 'run-data custodian', identity: `You are the custodian of this campaign's actual simulations. You have read every run spec, the handover, the campaign docs, and the results files, and you ground every statement in the real numbers. You know each run's grid, box, slab, CAP, sigma, dt, GS energy, and convergence flags, and which points are trusted/partial/unconverged/aliased/mis-referenced. Keep the debate honest: say what the files report, flag where a proposed explanation is contradicted by a run the others forgot, and identify which explanation holds across ALL runs. READ the files before answering — do not reconstruct from memory.`, files: true },
]

phase('Openings')
const openings = (await parallel(EXPERTS.map(e => () =>
  agent(
    PREAMBLE + '\n\nYOUR IDENTITY: ' + e.identity +
    (e.files ? '\n\nFILES TO READ (read them before answering):\n' + FILE_POINTERS : '') +
    '\n\nSHARED EVIDENCE BRIEF:\n' + BRIEF + '\n\n' + OPENING_TASK,
    { label: 'open:' + e.key, phase: 'Openings', model: 'opus', effort: 'high' }
  ).then(text => ({ key: e.key, label: e.label, text }))
))).filter(Boolean)

const openText = openings.map(o => '### ' + o.key + ' -- ' + o.label + '\n' + o.text).join('\n\n')

phase('Rebuttal')
const rebuttals = (await parallel(EXPERTS.map(e => () => {
  const mine = openings.find(o => o.key === e.key)
  const others = openings.filter(o => o.key !== e.key)
    .map(o => '### ' + o.key + ' -- ' + o.label + '\n' + o.text).join('\n\n')
  return agent(
    PREAMBLE + '\n\nYOUR IDENTITY: ' + e.identity +
    '\n\nSHARED EVIDENCE BRIEF:\n' + BRIEF +
    '\n\nYOUR OWN OPENING POSITION:\n' + (mine ? mine.text : '(not recorded)') +
    '\n\nTHE OTHER EXPERTS OPENING POSITIONS:\n' + others +
    '\n\nNow, in your role: (1) name the single strongest point another expert made that challenges your view; (2) concede anything you over-claimed; (3) sharpen or revise. Do not manufacture disagreement, and do not collapse into agreement to be agreeable. End with (a) your revised one-line verdict and (b) the single most decisive test.',
    { label: 'rebut:' + e.key, phase: 'Rebuttal', model: 'opus', effort: 'high' }
  ).then(text => ({ key: e.key, label: e.label, text }))
}))).filter(Boolean)

const transcript = '# OPENINGS\n\n' + openText + '\n\n# REBUTTALS\n\n' +
  rebuttals.map(r => '### ' + r.key + ' -- ' + r.label + '\n' + r.text).join('\n\n')

phase('Judge')
const verdict = await agent(
  'You are the judge for a four-expert scientific panel. You did NOT participate; weigh ALL four views and give the user a single, grounded, PROVISIONAL verdict they can push back on -- not a closed ruling.\n\nQUESTION:\n' + QUESTION +
  '\n\nEVIDENCE BRIEF:\n' + BRIEF +
  '\n\nFULL TRANSCRIPT (openings + rebuttals):\n' + transcript +
  '\n\nHOW TO JUDGE -- weigh on three axes: (a) consistency with known physics and limits; (b) parsimony (does the simplest model that actually fits win?); (c) decision-usefulness (does it tell the user what to do next?). Actively check for groupthink; if all four converged without stress-testing, say so and treat the consensus with suspicion.\n\nOUTPUT -- exactly this structure (this is what the user reads):\n1. Question -- restated in one or two lines.\n2. Reduced-model ladder -- the simplest model the panel settled on, the number it predicts, and where it breaks.\n3. Consensus -- what the experts agree is happening.\n4. Live disagreements -- the 1-3 genuine forks, what distinguishes them, the evidence each leans on.\n5. Best current answer -- grounded synthesis with a confidence level; label inferences "Inference:"; flag anything unverified; state plainly this is provisional and the user owns the verdict.\n6. Decisive next test -- the single most informative next step, concrete and where possible runnable on existing data or a cheap run.\n7. Open questions for the user.\nCarry units; round reported numbers to 2 significant figures.',
  { label: 'judge', phase: 'Judge', model: 'opus', effort: 'high' }
)

return { verdict, openings: openings.length, rebuttals: rebuttals.length }
```
