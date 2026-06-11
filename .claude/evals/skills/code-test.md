# Eval: code-test skill (LOCKED 2026-06-11)

Target skill (to be built in subtask 3): **code-test** = merge of
`physics-correctness` ∪ `development-feedback-loop`. Fires when finishing a code
function/feature (inqkit / inqview / utility) — before "done" or commit.

Evaluator types: **LLM-as-judge** (trigger + functional) + **human gate**.

## Trigger test

**Positive (should activate code-test):**
- "I added a new observable kernel to inqkit, is it done?"
- "finished the vec3 operator, let's commit"
- "wrote a new inqview analysis function — ready?"

**Negative (should activate a DIFFERENT skill — routing check):**
- "run a jellium sim at E=50" → `tddft-simulations`
- "write a handover" → `handover-update`
- "plot the stopping curve" → `report-figures`

Evaluator: LLM-judge decides which skill fired. Pass = all positives fire
code-test, no negative does.

## Functional test

Scenario: a new inqkit function has just been written and the user asks "is it
done?" / moves to commit.

Rubric (output must — all four are HARD pass/fail):
1. Refuses "compiles ⇒ works"; demands a known-case test before "done"/commit.
2. Proposes a known-case test **matched to the function type** (e.g. constant
   write → GPU reduce = n_pts; norm → ∈ [0.97, 1.03]; COD → analytic centre).
3. States the expected value **before** running it (anti-circularity — the
   expected value is independent of the code's own output).
4. Records what was tested + the observed result in the handover.

Evaluator: LLM-judge against the four criteria + human gate. PASS = 4/4.

## `<thinking>/<answer>` exemplar

Situation: just finished `center_of_density(n)`.

<thinking>
COD is defined as ∫ r·n(r) dV / ∫ n(r) dV. A known case: a uniform density in a
cubic box centred on the origin must have its centre of density exactly at the
box centre (0,0,0) by symmetry. So inject a constant n over the grid and the
result must satisfy |COD| < 1e-9 (round-off only). This value is fixed by the
definition, not read off from the function — so it is a genuine known answer.
</thinking>
<answer>
Test: fill the grid with n = 1.0 everywhere on an L-cubic box centred at origin;
assert |cod_x|, |cod_y|, |cod_z| < 1e-9. Run BEFORE declaring center_of_density
done; record "uniform-box COD = (0,0,0), |COD|<1e-9 ✓" in the handover.
</answer>

## Programmatic companion (in .claude/evals/programmatic/)

None — code-test is a reasoning skill. Its deterministic sliver (a new feature
has an accompanying test file; output has no NaN) is covered by the commit hook
eval, not here.
