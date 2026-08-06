# Handover — Cylindrical proximity ladder (weak → strong coupling)

Plan: `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/cylindrical-proximity-ladder.md`
Predecessor: `docs/handovers/cylindrical-channeling-ks-stopping.md` (rung r10, complete)

---

## 2026-08-03 — results read-out: S(T1) and S(T2) vs coupling (SUPERSEDES the earlier S table)

### What changed and why it matters

The earlier ladder table quoted S over the channeling twin's inherited windows
(`T1 9-25`, `T2 21-30`, `T2 5-20`) and T2 looked ERRATIC. That verdict was WRONG
and is retracted. Scanning the LOCAL slope -dE/ds(t) at every rung shows:

- S is ~0 at t=2 and peaks at t~8-11 in every estimator/rung — **wake build-up**
  (r_s=3 -> omega_p = sqrt(4 pi n) = 0.333 a.u., quarter period 4.7 a.u.). No
  steady-state S exists before t ~ 10.
- **S(T2) is NEGATIVE early** and crosses zero at t = 5.1 / 14.1 / 11.5 / 5.2
  (r10/r08/r06/r04); at r00 it never goes negative. `T2 5-20` straddles that sign
  change at every rung; `T2 21-30` sits in the late decay. Both averaged across a
  sign change — hence the apparent erraticism.

**One window is now used for all three estimators and all five rungs:
t = [11, 20] a.u.** After wake build-up; before the light-projectile velocity
criterion fails anywhere (classical drops below 0.85 v0 at t=20.6 at r00).
VERIFIED over the window: min v/v0 = 0.855 (classical), 0.917 (WP).

### The result (S in eV/Bohr, fit t=11-20)

| rung | R_in | f_wall | S(T1) | S(T2) | S_cl | T1/cl | T2/cl | r2(T2) |
|---|---|---|---|---|---|---|---|---|
| r10 | 2.5 sigma | 0.12 | 0.091 | 0.014 | 0.108 | 0.84 | 0.13 | 0.92 |
| r08 | 2.0 sigma | 0.37 | 0.123 | 0.026 | 0.179 | 0.69 | 0.15 | 0.66 |
| r06 | 1.5 sigma | 0.71 | 0.143 | 0.093 | 0.268 | 0.53 | 0.35 | 0.89 |
| r04 | 1.0 sigma | 0.93 | 0.162 | 0.175 | 0.366 | 0.44 | 0.48 | 0.99 |
| r00 | filled    | 1.00 | 0.222 | 0.218 | 0.483 | 0.46 | 0.45 | 0.99 |

**The headline: the two definitions BRACKET the classical answer at weak coupling
(0.84 vs 0.13, a factor 6.7) and CONVERGE at strong coupling (0.45 vs 0.46, 3%).**
Both stay below classical: the converged deficit is ~0.46, i.e. the fully
immersed wavepacket loses energy at less than half the classical rate, and the
agreement of two independent definitions rules out an estimator artefact.

**Mechanism, and it is linear:** the gap is exactly var(p)/2m. Fraction of the
drift loss diverted into momentum spread = 54 / 43 / 21 / 7.0 / 1.8 % across the
ladder. At weak coupling the projectile shears its own tail in momentum space
more than it decelerates; fully immersed the force acts on the whole packet.

Caveat kept in view: r08's T2 fit has r2 = 0.66 because its zero crossing (14.1)
is inside the window — read its 0.15 as "small", not as a value.

### Files (all committed-ready, none committed yet)

- `ResearchProject/systems/cylindrical_jellium/hypotheses/proximity_ladder/build_results_figures.py`
  — 5 figures + `results_summary.csv` + `window.json` into `figures/results/`.
- `.../hypotheses/proximity_ladder/build_results_notebook.py` -> `results.ipynb`
  (21 cells, 0 errors, 5 PNG + 2 embedded density GIFs, executed 2026-08-03).
- Figures: `F1_energy_loss_vs_path`, `F2_local_stopping_vs_time`,
  `F3_stopping_vs_coupling`, `F4_ratio_convergence`, `F5_variance_mechanism`.

`build_ladder_figures.py` and `phase_analysis.ipynb` are UNCHANGED and still
carry the old windows — they are the per-rung/campaign-internal record. The
results read-out is `results.ipynb`. If the old S table is ever quoted, quote
this one instead.

### Still not done

- `r04n160` same-N control (~2 GPU-h) — separates "wall is closer" from "more
  wall" (electron count runs 160 -> 326 to hold density fixed).
- SIC run at r00 to bracket the self-interaction bound against r10's 20.9%.
- `docs/validation/test-catalogue.md` still uncommitted (mixed with another
  session's 169 insertions) — user decision.

---

## 2026-08-02 — campaign designed, validated and launched unattended

### State: the whole campaign is queued as one dependency-chained graph

Submitted 2026-08-02 23:2x. Nothing further is required for it to complete.

| stage | job(s) | depends on |
|---|---|---|
| ground states r08/r06/r04/r00 | **32667942** (running) | — |
| RT build + guard check | **32668267** | — |
| smoke, 8 × (half, rung), 20 steps | 32668270/72/74/76/78/80/82/84 | `afterok` GS + buildcheck |
| production, 8 × (half, rung), 1500 steps | 32668271/73/75/77/79/81/83/85 | `afterok` its own smoke |
| figures + notebooks + cross-rung comparison | **32668286** | `afterany` all 8 production |

Watch: `squeue -u skcb2 -o '%.10i %.22j %.9T %.8M %.16R'`
Re-submit the whole graph: `shared/bin/submit-ladder.sh` (idempotent — GS skips populated checkpoints).

**Design choices that make it survive unattended:**
- Smoke is an `afterok` parent of production, so a rung with a bad ground state or
  mis-set geometry never burns 1.5 GPU-h to find out. Cost ~2 min per pair.
- Notebooks use `afterany`, so one failed rung does not suppress the write-up of
  the other three. `build_ladder_figures.py` reports missing rungs explicitly
  rather than silently shortening the ladder.
- Disk WARNS, never blocks (`.claude/rules/checkpoint-dont-block.md`; user
  confirmed space would be available). Only a filesystem too full to write a
  checkpoint without truncating is a hard stop.
- Every run is resumable (`.claude/rules/final-timestep-checkpoint.md`):
  `CH_RESUME=1` + larger `CH_N_STEPS` extends rather than recomputes.

### ⚠ THE RT BINARIES ARE FROZEN WHILE THE CHAIN IS QUEUED

`run-ladder-rt.slurm` deliberately **refuses to build** — it runs the prebuilt
`scripts/proximity_ladder/{wp,classical}/run`, produced once by the build check.
That is what lets 16 jobs fan out without racing on a shared build tree.

**Consequence: do NOT edit `proximity_ladder/{wp,classical}/run.cpp` or
`shared/configs/proximity_ladder_rs3.hpp` until the chain finishes.** The queued
jobs would keep using the old binary while any later rebuild used the new source,
so rungs run before and after the edit would silently differ — the one failure
mode this campaign cannot detect from its own outputs, because every rung would
still pass every gate.

If a change is genuinely needed mid-flight: cancel the pending RT jobs, edit,
re-run `run-ladder-buildcheck.slurm`, then re-run `submit-ladder.sh`. Completed
rungs are NOT reusable across such an edit and must be re-run.

Same reasoning applies to `inqkit/jellium/localised_background.hpp` — it is
header-only and compiled INTO both binaries.

### Scheduling note (2026-08-02 23:45)

The 8 smoke jobs sat at `(Priority)`: the partition was loaded (129 running / 526
pending) and the account was at its GPU-minutes cap
(`32667766_3 → AssocGrpGRESMinutes`), with 7 unrelated `s56-*` jobs holding
**24 h** requested wall time each.

SLURM charges `AssocGrpGRESMinutes` against **requested**, not used, wall time,
and backfill favours short jobs. The original 12 h request on every RT job — a
mistake, since smoke needs ~4 min — was therefore both burning the shared budget
and losing backfill slots. Corrected live via `scontrol` and at source:
**smoke 30 min, production 6 h** (≈7× and ≈3.3× margins on measured/estimated
runtimes). Keep those limits on any re-submission.

### The ladder

Fixed at r_s = 3.000000, 40×40×60 Bohr, dx 0.5, 50 eV projectile (v/v_F = 3.00),
σ_WP = 4. **R_out is DERIVED, never transcribed** — it is solved so n₀ hits the
r_s = 3 target exactly. Verified: max density error over the whole table = **2.2e-16**.

| rung | R_in | R_in/σ | N_e | R_out | states | WP charge in wall @t=0 | shape |
|---|---|---|---|---|---|---|---|
| r10 | 10 | 2.50 | 160 | 14.000 | 104 | 0.19 % | annulus — **DONE, not re-run** |
| r08 | 8 | 2.00 | 220 | 14.000 | 143 | 1.83 % | annulus |
| r06 | 6 | 1.50 | 266 | 13.986 | 173 | 10.54 % | annulus |
| r04 | 4 | 1.00 | 300 | 14.000 | 195 | 36.79 % | annulus |
| r00 | 0 | — | 326 | 13.986 | 212 | 100 % | **cylinder** (filled) |

Control (`r04n160`, not yet submitted): R_in = 4 at N = 160, R_out = √112 = 10.583 —
same jellium volume as r10, so it separates "the wall is closer" from "there is
more wall". One pair, ~2 GPU-h.

### VERIFIED

- **inqkit test suite: pure 13/13, engine 29/29** (job 32667729). New cases
  T0.7/T0.8/T0.9 confirmed compiled into the binary (built 23:11, after the 23:03
  source edit) via `strings`.
- **Rung table**: r_s = 3.000000000 at every rung, density error 2.2e-16.
  Control's jellium volume equals r10's to 6 d.p.
- **ALL FOUR GROUND STATES COMPLETE**, every gate PASS, `rungs with a non-zero
  exit : 0` (job 32667942, 1568 s total). Saved under
  `shared_gs/ladder_rs3_{r08,r06,r04,r00}_L60_dx0p5`.

  | rung | N_e | E_GS (Ha) | bore depletion | n₊/n₀ on axis | s |
  |---|---|---|---|---|---|
  | r08 | 220 | −120.24 | 0.150 | 0 | 454 |
  | r06 | 266 | −189.90 | 0.194 | 0 | 316 |
  | r04 | 300 | −261.89 | 0.259 | 4.6e-11 | 377 |
  | r00 | 326 | −342.79 | — (filled) | **1** | 421 |

  Two things worth keeping:
  * **The filled-cylinder gate is the fix, confirmed in a live SCF.** r00 reports
    n₊ = 0.0088419412828831 on the axis against a target of 0.0088419412828831 —
    ratio exactly 1. The old `annulus`-at-R_in=0 path would have put 0.00442
    there, a 50 % density hole centred on the projectile's flight path.
  * **The hollow axis gate measures a real erfc tail, not a hardcoded zero.**
    r04 reads 4.6e-11 because the r < 2 Bohr probe rim sits 4σ from the R_in = 4
    edge (½erfc(−4) ≈ 7.7e-9); at R_in ≥ 6 the argument reaches erfc(−8) ≈ 1e-29
    and underflows. A clean 0 at every rung would have been the suspicious result.
  * **Bore depletion climbs monotonically** 0.129 → 0.150 → 0.194 → 0.259 as the
    wall closes. The channel fills in smoothly: the coupling increase the campaign
    exists to measure is already visible before a projectile has been fired.
- **Both RT binaries compile** (job 32668084, wp 244 s + classical). Guard 1
  (unset `CJ_RUNG` refuses) passes on both.
- **The entire analysis stage, end-to-end against real r10 data**: 47 per-rung
  figures + 3 comparison panels + `ladder_summary.csv` + 2 notebooks. It
  **reproduces the known result exactly** — S_wp/S_cl = 0.0877/0.1096 = **0.80**,
  matching this study's previously recorded 0.801 over the 9–25 window. Measured
  f_wall(0) = 0.204 % against the analytic 0.19 %.

### The bug this campaign uncovered (fixed)

`background_shape::annulus` evaluated at R_in = 0 returns **n₀/2 exactly on the
tube axis**, relaxing to n₀ only by d ≈ 2w. The erfc step is centred ON its
nominal edge, so `background_mask(0, 0, w) = ½` for every w > 0.

Silent, and maximal precisely where a channeling projectile flies. Existing test
T0.5 could **not** have caught it: it probes R_in = 0 at w = 0, where
`background_mask(d,0,0) = 0` for all physical d ≥ 0 and the composition is
accidentally right. The defect existed only in the softened branch — i.e. only in
production, where w = 0.5.

Fix (user's call, 2026-08-02): a filled tube is **its own shape**,
`background_shape::cylinder` / `cylinder_mask`, not an annulus with a degenerate
inner edge. Purely additive — nothing else switches on the enum and no existing
run passes `inner_radius = 0`. Ladder ground states also gate on it at runtime
(gate 3b: background n₊/n₀ on the axis must read ~1 filled, ~0 hollow).

### 2026-08-03 — the f_bore gate was wrong, not the run (r04 false abort)

**Symptom.** `wp/r04/smoke` (32668278) aborted at t=0:
`[FAIL] f_bore(0) (Rayleigh): 0.62075 (expect 0.63212, dev -1.80 %, tol +/-1 %)`.
r08 and r06 passed the same gate; r00 uses the filled branch and was unaffected.

**The run was fine.** Every other t=0 gate passed to 1e-9 or better, and decisively
`<r_perp>(0) = 3.544338` vs analytic `3.544908` — **0.016 %**. `<r_perp>` is a
SMOOTH moment of the same density; `f_bore` is that density integrated over a
SHARP cylinder. A genuinely mis-injected packet moves both. Only the sharp cut
moved, so the cut was grid-limited, not the packet.
(`max_overlap = 3.0e-07` also confirms Pauli blocking is negligible as predicted.)

**Root cause — a fixed relative tolerance is not a fixed standard along this ladder.**
`f_bore = 1 − exp(−R²/2σ_d²)` has relative sensitivity `(df/dR)/f` per Bohr of:

| rung | R_in | f_bore | (df/dR)/f | ±1 % implies dR < |
|---|---|---|---|---|
| r10 | 10 | 0.9981 | 0.0024 | 4.1 Bohr (vacuous) |
| r08 | 8 | 0.9817 | 0.0187 | 0.54 Bohr |
| r06 | 6 | 0.8946 | 0.0884 | 0.11 Bohr |
| r04 | 4 | 0.6321 | 0.2910 | **0.034 Bohr = 7 % of a grid cell** |

A factor of **120** across the ladder. One constant effective-radius error of
**dR = 0.062 Bohr (12 % of a 0.5-Bohr cell**, i.e. ordinary staircase
discretisation of a cylindrical boundary on a Cartesian grid) reproduces the
observed deviation at EVERY rung: 0.01 / 0.12 / 0.55 / **1.80** %. That is exactly
the pass/pass/pass/FAIL pattern seen. The tolerance was silently tightening as the
wall closed in.

**Fix** (`proximity_ladder/wp/run.cpp`, `rayleigh_tol_pc`): 1 % floor in quadrature
with the grid-geometry term `100·(df/dR)/f·(dx/2)` → **1.0 / 1.1 / 2.4 / 7.3 %** at
r10/r08/r06/r04, and 1.0 % for the filled branch (R_out is deep in the flat tail).
Still a real gate everywhere — it asks for grid-achievable accuracy instead of
sub-cell accuracy.

**r06 passed by 8 % — the gate was marginal, not just wrong at one rung.**
Measured across the smoke stage:

| rung | measured dev | old tol | margin | implied dR (Bohr) | new tol | new margin |
|---|---|---|---|---|---|---|
| r08 | −0.087 % | 1 % | 11× | 0.047 | 1.10 % | 12× |
| r06 | **−0.919 %** | 1 % | **1.09×** | 0.104 | 2.42 % | 2.6× |
| r04 | −1.799 % | 1 % | **FAILED** | 0.062 | 7.34 % | 4.1× |
| r00 | −0.00002 % | 1 % | — | — | 1.00 % | — |

Implied radius errors span 0.047–0.104 Bohr (9–21 % of a 0.5-Bohr cell) — the
staircase error is not constant, because it depends on how each radius commensurates
with the grid. The `dx/2 = 0.25` Bohr term covers the worst observed case with 2.4×
margin. Note r06 would have aborted on any minor grid-alignment difference, so this
was a latent failure across the ladder rather than an r04 peculiarity.

**Repair chain (autonomous).** `wp/run` could not be rebuilt while wp jobs were
executing it (Linux returns ETXTBSY), so: **32669422** rebuild (`afterany` on all 9
live RT jobs) → **32669423** r04 wp smoke → **32669424** r04 wp prod. Notebook job
32668286 had its dependency extended to include 32669424, so the write-up still
covers the full ladder.

**Binary-consistency note.** r08/r06/r00 ran the pre-fix `wp/run`, r04 the post-fix
one. The change is confined to a t=0 assertion's tolerance — it touches no
wavefunction, Hamiltonian or observable, so propagated physics is bit-identical.
This is the one justified exception to the freeze rule above; a change with any
physics content would instead require re-running every rung.

### 2026-08-03 — first results (r10 vs r08, both twins complete)

| rung | R/σ | mean f_wall over fit (range) | energy lost | S_wp/S_cl (T₁ 9–25) | var(p_z) growth |
|---|---|---|---|---|---|
| r10 | 2.5 | 0.123 (0.033–0.164, 4.9×) | 10.3 % | **0.801** | 44.5 % |
| r08 | 2.0 | 0.345 (0.157–0.423, 2.7×) | 16.6 % | **0.644** | 59.7 % |

1. **The weak→strong axis is real and measured.** Classical loss 5.13 → 8.29 eV
   (10.3 % → 16.6 % of 50 eV); measured coupling roughly triples. The rungs' coupling
   ranges barely overlap (0.033–0.164 vs 0.157–0.423), i.e. they TILE the coupling
   axis rather than resampling it — the ladder working as designed.
2. **T₁ degrades with coupling**: 0.801 → 0.644.
3. **T₂ collapses.** Over t = 5–20 the ratio is **0.011** at r08 (vs 0.132 at r10):
   var(p) growth of 59.7 % now almost exactly cancels the drift loss, so total
   kinetic energy reports a projectile that has stopped losing energy while its
   classical twin loses 8.3 eV. This is the bulk contamination the channeling
   geometry was built to escape, returning under load. It is the clearest single
   result so far.
4. Cost model validated: 4381 s measured vs 4248 s predicted for wp/r08 (3 %).
5. `max_overlap` after a full 1500-step propagation = 3.7e-08 — the WP stayed
   orthogonal to the occupied manifold throughout; nothing leaned on Pauli blocking.

### 2026-08-03 04:58 — CAMPAIGN COMPLETE, all 5 rungs, `failures: 0`

Artefacts: `hypotheses/proximity_ladder/` — `rung_{r10,r08,r06,r04,r00}.ipynb`
(31 figures each), `ladder_comparison.ipynb` (4 panels), `figures/ladder_summary.csv`.

| rung | R/σ | mean f_wall (drift) | loss | S_wp (T₁) | S_cl | **ratio** | var(p) |
|---|---|---|---|---|---|---|---|
| r10 | 2.5 | 0.123 (4.9×) | 10.3 % | 0.0877 | 0.1096 | **0.801** | 44.5 % |
| r08 | 2.0 | 0.345 (2.7×) | 16.6 % | 0.1173 | 0.1822 | **0.644** | 59.7 % |
| r06 | 1.5 | 0.640 (1.7×) | 23.7 % | 0.1353 | 0.2708 | **0.499** | 68.4 % |
| r04 | 1.0 | 0.873 (1.4×) | 31.8 % | 0.1524 | 0.3682 | **0.414** | 81.4 % |
| r00 | filled | 0.989 (1.04×) | 41.6 % | 0.2017 | 0.4859 | **0.415** | 115.0 % |

**THE HEADLINE: the T₁ ratio SATURATES at ~0.41.** r04 and r00 agree to 0.3 %
despite r00 carrying 1.3× the classical stopping power and higher coupling. The
drift estimator's under-reporting bottoms out at ~59 % rather than continuing to
degrade. (Supersedes the 4-rung note below, which read the r06→r00 step as a
"turnover" — with r04 filled in it is a genuine plateau from r04 onward.)

Two reasons the plateau is trustworthy:
* **Coupling drift within the fit window SHRINKS toward the plateau** — 4.9× (r10)
  → 1.04× (r00). The strong rungs are the ladder's BEST-conditioned points, so the
  saturation is measured where the coupling is most sharply defined. The
  coupling-averaging caveat applies mainly to the WEAK end, not the plateau.
* **var(p) growth keeps climbing through it** (81.4 → 115.0 % from r04 to r00)
  while the ratio holds flat — the mechanism does not switch off; the drift channel
  simply stops losing further ground.

**Usable conclusion (final):** a KS-orbital DRIFT definition of stopping power is
good to ~20 % in the channeling limit, degrades through intermediate coupling, and
saturates at ~59 % under-reporting once the projectile is substantially immersed
(f_wall ≳ 0.87). It is not a substitute for the classical ΔE/Δs definition outside
weak coupling, but its error is BOUNDED rather than unbounded.

Also: r04 shows the most extreme spreading of any rung, f_bore 0.621 → 0.0046 —
at R_in = 1σ the channel barely exists as a confining structure by the run's end.

### 2026-08-03 — LADDER RESULT (4 of 5 rungs; r04 wp in repair) [superseded above]

| rung | R/σ | mean f_wall | loss | S_wp (T₁) | S_cl | **ratio** | var(p_z) growth |
|---|---|---|---|---|---|---|---|
| r10 | 2.5 | 0.123 | 10.3 % | 0.0877 | 0.1096 | **0.801** | 44.5 % |
| r08 | 2.0 | 0.345 | 16.6 % | 0.1173 | 0.1822 | **0.644** | 59.7 % |
| r06 | 1.5 | 0.640 | 23.7 % | 0.1353 | 0.2708 | **0.499** | 68.4 % |
| r00 | filled | 0.989 | 41.6 % | 0.2017 | 0.4859 | **0.415** | 115.0 % |

(S in eV/Bohr, T₁ = drift channel, window t = 9–25.)

**1. The campaign's aim is met.** Fractional energy loss 10.3 % → 41.6 % —
perturbative to strongly non-linear — at fixed r_s = 3, fixed 50 eV, fixed
σ_WP = 4, fixed cell/grid. Coupling f_wall spans 0.12 → 0.99.

**2. The T₁ ratio falls monotonically and SATURATES**: 0.801 → 0.644 → 0.499 →
0.415. Step sizes 0.157, 0.145, 0.084 — the last step is roughly half the
previous, so it is turning over, not continuing linearly (a 3-point extrapolation
would have predicted ~0.35 at the endpoint; measured 0.415).

**3. The mechanism, stated quantitatively.**

    S_classical:  0.110 -> 0.182 -> 0.271 -> 0.486     x4.4
    S_wp (T1):    0.088 -> 0.117 -> 0.135 -> 0.202     x2.3

Classical stopping tracks the sampled density (f_wall x8, S x4.4). The drift
channel captures only ~half that response, and the shortfall (0.199 → 0.356 →
0.501 → 0.585) rises in lockstep with var(p) growth (44.5 → 59.7 → 68.4 →
115.0 %). The missing stopping power is momentum transferred into the packet's
WIDTH rather than into decelerating its centroid — the bulk contamination the
channeling geometry was built to avoid, returning progressively as the wall closes.

**Usable conclusion:** a KS-orbital DRIFT definition of stopping power is good to
~20 % in the channeling limit and under-reports by ~60 % when the projectile is
immersed. It is not a drop-in substitute for the classical ΔE/Δs definition
outside weak coupling.

**4. T₂ is erratic, not systematically wrong.** Ratios over t = 5–20 read 0.132,
0.011, 0.134 across r10/r08/r06. The near-zero at r08 was a coincidental
cancellation in that window, NOT a trend — an earlier note in this handover called
it "the clearest single result"; that was wrong on three points of data. The honest
statement is that T₂ is unstable across rungs and unusable at any coupling.

**5. Pauli blocking never mattered**, even fully immersed: max_overlap = 3.7e-08
(r08) and 4.7e-07 (r00) after full 1500-step propagations, ~1e4 under threshold.
The k₀ = 1.92 vs k_F = 0.64 momentum separation holds throughout, as predicted.

**6. Cost model held**: predicted/measured wall time 4248/4381 (r08), 6296/5996
(r00) — within 5 %.

#### CORRECTION to the plan's fit-window reasoning (was wrong; verified against data)

The plan (§3.1, §5) argued that a fixed-TIME window cannot compare rungs because
"the rungs merge in time". **That is wrong for this ladder** and acting on it would
have discarded a valid comparison:

* every rung shares σ_WP = 4, dt, n_steps and v₀, so the packet spreads
  IDENTICALLY. At fixed t all rungs have the same packet size and differ only in
  where the wall sits — the coupling difference IS the independent variable.
* velocity is not the issue: measured v/v₀ = 0.99 → 0.93 over t = 9–25 at r08,
  well inside the `light-projectile-stopping` criterion of v ≥ 0.85 v₀.

The REAL limitation is coupling drift **within** the window — 4.9× (r10), 2.7× (r08)
— so each S is an average over a coupling RANGE that differs by rung. That cannot be
fixed by a better window: constraining the drift to 1.5× needs t < 2.4 a.u. (~120
steps), too short for a stable fit. It is intrinsic to a spreading wavepacket.

Resolution implemented: `rung_summary_row` reports `fw_fit_mean/lo/hi/drift`, and
`L04_ratio_vs_measured_coupling` draws S(coupling) with horizontal error bars. The
deliverable is a curve over coupling ranges, not points at single values.

### NOT verified — what could still stop a rung

1. **The RT binaries have never run ladder physics.** They compile and their
   guards fire, but the first real execution is the smoke stage. Most exposed:
   the `FILLED` branch of the t=0 gate and `radial_occupancy` at r_inner = 0
   (f_bore is identically 0 there — hence the explicit filled branch, since
   `gate_rel` treats `want == 0` as an automatic pass and would have silently
   stopped testing on the newest geometry).
2. **SCF convergence at r04 / r00.** A filled 14-Bohr nanowire has dense
   near-degenerate subbands and the smearing is 0.00862 eV, which is cold for
   that. If SCF struggles, raise `extra_states` before raising T.
3. **⟨n_bath⟩_WP is not implemented.** The plan names it the ideal coupling
   coordinate. The campaign currently uses **measured f_wall(t)** instead, which
   is written every step, exists at every rung including the filled one, and
   validated against the analytic Rayleigh value. Adding ⟨n⟩_WP later is an
   improvement, not a prerequisite.
4. **Buildcheck guards 2/3** were re-submitted as 32668267 after a harness bug
   (GNU `env` stops option parsing at the first NAME=VALUE, so
   `env CJ_RUNG=x -u CH_GS_DIR ./run` executes `-u` and exits 127). The binaries
   were never at fault. Confirm 32668267 reports `failures : 0`.

### The scientific caveat that must survive into the write-up

**This ladder cannot reach the textbook "full stopping power" regime, and no
choice of R_in changes that.** The projectile's Gaussian form factor
exp(−q²σ_pot²/2) is 0.37 at q = 0.5, 0.018 at q = 1, and 3e-26 at q = 2v₀ = 3.83.
The plasmon pole sits at q_min = ω_p/v = 0.174; the electron–hole continuum runs
to q = 2v. **The projectile couples to the collective response and essentially
nothing else, at every rung.** Shrinking R_in scales how much medium responds; it
does not harden the projectile.

The ladder is therefore **weak-collective → strong-collective**. Reaching the pair
channel is a **σ_WP axis**, not an R_in axis — and the vacuum σ-sweep
(`systems/vacuum/hypotheses/wp_selfinteraction/sigma_sweep.py`) already maps its
self-interaction cost.

Three further interpretation points, all recorded in the run.cpp headers so the
notebooks repeat them:
- **The rungs merge in time.** Wall-overlap spans a factor of 190 across the
  ladder at t = 0 but only 3 by t = 30. S must be fitted in a window keyed to a
  common *measured* coupling, never a common time.
- **Three variables move together** — proximity, N_e (160→326), and the target's
  mode spectrum (thin annulus → solid nanowire). Inseparable in this geometry
  because the electrons added *are* the close ones. The `r04n160` control is what
  disentangles them.
- **Pauli blocking is negligible; xc coupling is not.** k₀ = 1.917 sits 7.2 σ_p
  above k_F = 0.640, so the Gram–Schmidt overlap is ~1e-11 even at 100 % spatial
  overlap (monitor `report.max_overlap`; it should stay ≲1e-10). But the WP adds
  ~32 % to the on-axis density and LDA xc is non-linear in n, while the classical
  twin contributes nothing to xc at all. Part of the WP−classical gap at the
  strong rungs is that, not dispersion.

### Files

| path | role |
|---|---|
| `inq-stack/include/inqkit/jellium/localised_background.hpp` | `background_shape::cylinder` + `cylinder_mask` (**changed**) |
| `inq-stack/tests/include/inqkit/jellium/test_localised_background_engine.cpp` | T0.7/T0.8/T0.9 (**changed**) |
| `.../cylindrical_jellium/shared/configs/proximity_ladder_rs3.hpp` | rung table, derived R_out, self-check |
| `.../scripts/proximity_ladder/{gs,wp,classical}/run.cpp` | build-once, `CJ_RUNG` at runtime, no default |
| `.../hypotheses/proximity_ladder/build_ladder_figures.py` | reuses channeling_twin's drawing engine; shared axis limits across rungs |
| `.../hypotheses/proximity_ladder/build_ladder_notebooks.py` | per-rung + comparison notebooks |
| `shared/bin/run-ladder-{gs,buildcheck,rt,notebooks}.slurm` | the four stages |
| `shared/bin/submit-ladder.sh` | the dependency graph |

Nothing committed to git.

### Next

1. Confirm 32668267 reports `failures : 0`.
2. Watch the first smoke pair (r08) — it is the first ladder physics ever run.
3. When the chain finishes: read `hypotheses/proximity_ladder/figures/ladder_summary.csv`
   and `ladder_comparison.ipynb`.
4. Optional: submit the `r04n160` same-N control, and a SIC run at r00 to bracket
   the self-interaction bound (r10 measured 20.9 %).
