# Handover: cylindrical (annular) channeling twin — KS-orbital stopping power

**Rolling file. Latest milestone at top.**
**Repo:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research` (branch `quantum-stopping-power`)
**Machine:** CSD3, login node `login-q-2`, GPU partition `ampere` (A100, sm_80)
**Plan:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/cylindrical-channeling-ks-stopping.md`

---

## 2026-08-02 (7th) — SIC-PZ channeling run CONFIRMS the 6th entry, a width-definition error in it is FIXED, sigma sweep launched, report figures built

Three things: a defect fix, a confirmation, and two deliverables.

### 1. DEFECT (mine, 6th entry) — two different width definitions were divided into each other

`selfinteraction.py` compared the vacuum **3-D geometric-mean** excess against
`CHANNELING_EXCESS = 1.467`, which is the **transverse** `<r_perp>/free` ratio
from the 5th entry. The channeling packet is strongly anisotropic (transverse
1.52, longitudinal 1.13), so the two definitions genuinely differ. Consequence:
the reported fraction was **16.9 %** where the matched-definition answer is
**20.9 %**. Nothing else in the 6th entry is affected — the vacuum numbers
(1.436 / 1.079, 82 % xc cancellation) are all internally 3-D and stand.

Fixed: `CHANNELING_EXCESS = 1.3778` (3-D iso, uncorrected),
`CHANNELING_EXCESS_SIC = 1.2990`, `CHANNELING_EXCESS_TRANSVERSE = 1.467` kept
for cross-reference, all with provenance. New regression test
`test_channeling_constants_share_one_width_definition`. **14 tests pass.**

### 2. CONFIRMATION — the direct SIC test agrees with the vacuum prediction to 0.4 %

Job **32615191** (`channeling_sic/wp/results/wp_sic`, `sic_mode = sic-pzrun`,
1500 steps, exit 0) is identical to `channeling_twin/wp/results/wp` in every
physical parameter — same GS, box, sigma, v0 — differing only in `sic_mode`.

| | excess vs free (3-D iso) | var(p_z) growth |
|---|---|---|
| uncorrected LDA | 1.3778 | +44.5 % |
| SIC-PZ | 1.2990 | +29.0 % |
| removed by SIC | **20.9 % of the excess** | |
| PREDICTED from vacuum | **20.9 %** | |

The vacuum self-interaction transfers into the bore essentially unchanged. This
**closes candidate 2 of the 6th entry** (environment-dependent xc cancellation):
it is not significant here. So ~79 % of the excess spreading is real physics —
independently consistent with the SIC handover's stopping-power finding that
~91 % of the T1 deficit is genuine.

Vacuum validation that licenses this: `sic_pzrun` returns the width ratio to
**1.00000** (var(p) drift 0.002 %) — the correction is exact where the answer is
known. `sic_h` gives **0.5549**, i.e. 44 % UNDER-spread, confirming the 6th
entry's prediction that Hartree-only correction over-corrects badly.

### 3. sigma_WP sweep — LAUNCHED (array 32625669, 6 tasks x 4 configs = 24 runs)

`shared/bin/run-wp-si-sweep.slurm`. sigma in {1, 2, 3, 4, 6, 8}; levels
noninteracting / hartree / lda / lda+sic_pzrun.

**Scale-exact protocol** — free evolution is invariant under `r -> lambda r`,
`t -> lambda^2 t`, so everything is scaled off sigma: `L = 18 sigma`,
`h = 0.125 sigma`, `dt = 0.00125 sigma^2`, 1500 steps. Every sigma is then the
SAME discrete problem — always 144^3, always `sigma_dens(0)/h = 5.657`, always
box half-width `= 5.99 sigma_dens(t_end)` — so cost is sigma-independent
(~895 s/run, measured) and **sigma = 4 reduces exactly to the completed run**
(L=72, h=0.5, dt=0.02, t_end=30), making it a free reproduction check.

What it measures: under that rescaling only the **dimensionless coupling**
changes. Kinetic ~ 1/sigma^2 while both E_PP ~ 1/sigma and LDA `v_x ~ n^(1/3) ~
1/sigma`, so the self-interaction-to-kinetic ratio ~ sigma. Since Hartree and LDA
exchange carry the SAME 1/sigma scaling, an exchange-only cancellation would be
sigma-INDEPENDENT; LDA correlation does not follow the similarity. The sweep
therefore tests whether the 82 % cancellation is a constant of the functional or
an accident of that density. (Inference from the functional forms; the run decides.)

**Fit E_PP = A/sigma + C, not a pure power law.** Measured `E_PP(0) = 2.1814 eV`
at sigma = 4 vs free-space `1/(sigma sqrt(2 pi)) = 2.714 eV`; the 0.0196 Ha gap
matches the Madelung term `xi/2L = 2.837/144 = 0.0197 Ha` from INQ dropping G=0
in a charged periodic cell. The offset is sigma-independent, so a naive power-law
fit would return a fake exponent.

### 4. Report-ready figures — BUILT

`hypotheses/channeling_twin/build_report_figures.py` ->
`hypotheses/channeling_twin/report_figures/{setup,twin,sic}/`, **86 PNGs at 600
dpi**, one panel per file (never pre-composed), identical filenames across
`twin/` and `sic/` so they diff panel-for-panel.

- `setup/` — GS density xz + yz slices, linear and log. Hard geometry assertion:
  density peaks at |x| = 11.75 Bohr, inside the 10-14 annulus (a centre<->edge
  VTI swap would put it on the axis). `n_max = 8.27e-3` vs the r_s=3 target
  8.84e-3.
- Standard applied: canonical `inqview.visualisation.style`, no on-canvas titles,
  time in **fs** with a secondary a.u. top axis (the windows were chosen in a.u.),
  **axis limits shared across twin and sic via a two-pass driver**, one shared
  momentum clim, linear+log variants for every map, mathtext x10^n colorbars,
  fixed canvas (`bbox_inches=None`).
- `style.figure_one_col()` could NOT be used for line panels: its rect leaves
  ~3.5 % headroom and clips the secondary axis. Local `figure_line()` (3.5 x 3.3
  in, fixed inch margins) replaces it; `figure_map()` gives 2-D maps equal aspect
  and a non-clipping colorbar rect.

### 5. ~~OPEN DISCREPANCY~~ — **RETRACTED. The alarm below was MY analysis error; all numbers stand.**

> **RESOLVED, same day.** The sigma = 4 sweep point reproduces the original
> fixed-box runs **exactly**: hartree 1.43643 vs 1.43643, lda 1.07880 vs 1.07880,
> difference **0.00000** in both. The binary change altered nothing, and
> **"LDA excess 1.079 / xc cancels 82 % / SIE explains 20.9 %" all STAND.**
>
> The apparent discrepancy came from analysing the sweep WHILE THE `lda` RUNS WERE
> STILL WRITING. `_excess_at` used `argmin` over tau, which silently returned each
> partial run's last row — at a SMALLER tau — and divided it by a COMPLETED
> reference at full tau. The resulting LDA excess was ~40 % low, and (the reason it
> fooled me) perfectly smooth and monotonic in sigma, so it looked like physics.
> The self-consistent hartree pair was complete and therefore correct, which is
> exactly why only the lda curve looked "off-trend".
>
> Fixed: `_excess_at` now returns NaN when either run stops short of the requested
> tau — a partial run is a refusal, never a nearest-neighbour match. Regression
> test `test_a_run_that_stops_short_of_tau_gives_NaN_not_a_nearest_match`.
> **23 tests pass** in `hypotheses/wp_selfinteraction/tests/`.
>
> Lesson for the next session: this file's own rule — never analyse a run set
> before every member reports `run_completed = true`. The gate existed
> (`protocol_gate` carries a `complete` column) and I did not consult it.

The original (now void) alarm is kept below for the record only.

At sigma = 4 the three original runs were NOT all produced by the same binary:

| run | sic.csv | `sic_mode` line | binary |
|---|---|---|---|
| `noninteracting` | no | no | OLD |
| `hartree` | no | no | OLD |
| `lda` | **yes** | **yes** | **NEW** (SIC added) |

So `hartree - noninteracting` is a self-consistent pair, but
`lda - noninteracting` differences a NEW-binary run against an OLD-binary
reference. The sweep (all four configs on one current binary) says the mixed pair
is the outlier:

| excess - 1 | sigma=1 | sigma=2 | sigma=3 | extrapolated sigma=4 | measured sigma=4 |
|---|---|---|---|---|---|
| hartree (consistent pair) | 0.1176 | 0.2289 | 0.3349 | **0.436** | **0.4364** — agrees to 0.1 % |
| lda (MIXED pair) | 0.0210 | 0.0375 | 0.0487 | **~0.057** | **0.0788** — off by ~40 % |

Also: the sweep's xc cancellation RISES with sigma (82.2 / 83.6 / 85.5 % at
sigma = 1/2/3), whereas the sigma = 4 mixed pair gives 81.9 % — breaking
monotonicity in the same direction.

If the sweep's sigma = 4 point reproduces ~1.057 rather than 1.079, then at
sigma = 4: xc cancellation ~86 % (not 82 %) and the vacuum-predicted share of the
channeling excess ~15 % (not 20.9 %). **That would break the headline agreement
with the directly-measured SIC fraction (20.9 %)** and would mean the vacuum
number does NOT transfer as cleanly as entry-7 §2 claims — reopening candidate 2.
The direct SIC channeling measurement (1.3778 -> 1.2990) is unaffected either
way: both of those runs are same-binary and same-everything.

Resolution: sweep task 3 (sigma = 4) re-runs all four configs on one binary. It
was included as a reproduction check for exactly this class of error. Compare
`sweep_s4p0_*` against `results/{noninteracting,hartree,lda}`; a clean
reproduction of 1.4364 with a DIFFERENT lda number localises the fault to the
binary change.

Analysis machinery is ready: `hypotheses/wp_selfinteraction/sigma_sweep.py`
(+ `tests/test_sigma_sweep.py`, **8 tests pass**). `protocol_gate()` confirms
`E_PP(0) * sigma = 0.320667 Ha.Bohr` identical at sigma = 1/2/3 and 1.0017x the
analytic `[1/sqrt(2pi) - xi/36]/sigma = 0.32013` — the scaled protocol holds.
NOTE this supersedes entry-7 §3's "fit E_PP = A/sigma + C": with `L = 18 sigma`
the Madelung term scales too, so E_PP is a PURE 1/sigma with no additive constant.

### 6. SWEEP COMPLETE — all 24 runs, and the xc cancellation is NOT a constant

Array 32625669 finished, all six sigma x four levels. Protocol gate passes at
every sigma: `E_PP(0)*sigma = 0.320667 Ha.Bohr` identical to six digits,
reference width vs closed form <= 5e-12, var(p) drift ~0.

| sigma | hartree excess | lda excess | xc cancellation | sic_pzrun residual |
|---|---|---|---|---|
| 1 | 1.1176 | 1.0314 | 73.3 % | **0.0** |
| 2 | 1.2289 | 1.0527 | 77.0 % | **0.0** |
| 3 | 1.3349 | 1.0679 | 79.7 % | **0.0** |
| 4 | 1.4364 | 1.0788 | 81.9 % | **0.0** |
| 6 | 1.6285 | 1.0919 | 85.4 % | **0.0** |
| 8 | 1.8091 | 1.0970 | 88.0 % | **0.0** |

(at fixed DIMENSIONLESS time tau = 1.875; `sigma_table()` in
`hypotheses/wp_selfinteraction/sigma_sweep.py`)

Three results:

1. **The bare Hartree self-interaction grows ~linearly in sigma** (excess-1 =
   0.118 -> 0.809), confirming the dimensionless-coupling argument: coupling ~
   sigma, so a wider packet is more self-interacting per unit kinetic energy.
2. **The xc cancellation RISES monotonically with sigma, 73 % -> 88 %.** It is
   therefore NOT a constant of the functional. Since LDA *exchange* carries the
   same 1/sigma scaling as Hartree (which would give a sigma-independent
   cancellation), the drift is attributable to LDA **correlation**, whose density
   dependence does not follow the similarity. The sweep's stated purpose was to
   decide this, and it decided it.
3. **The NET LDA error SATURATES near ~10 %** (excess-1 = 0.031, 0.053, 0.068,
   0.079, 0.092, 0.097): xc cancellation improves with sigma faster than the bare
   coupling grows. Practical reading — going to a wider packet does not make the
   self-interaction artefact much worse; going narrower makes it markedly better
   (3 % at sigma = 1).
4. **SIC-PZ returns the width ratio to EXACTLY 1.0 at every sigma** (residual 0.0
   across a 73-88 % cancellation range). The correction is exact over the whole
   coupling range, not just at the one sigma it was validated at.

Caveat on reading the table: this is fixed *dimensionless* time. At fixed
PHYSICAL duration use `at_fixed_physical_time()`, which covers only
`sigma >= sqrt(T/1.875)` (sigma >= 4 for T = 30 a.u.) and returns the dropped
sigma list explicitly rather than truncating silently.

### 7. Report figures — EXTENDED (98 channeling + 10 vacuum)

`hypotheses/channeling_twin/build_report_figures.py` -> `report_figures/`
(47 twin + 47 sic + 4 setup). New since entry 6:

- `05c/05d_energy_loss_*` — classical, T1 and T2 energy LOSS on one positive
  axis (raw deltas have opposite natural signs, which made slopes read backwards
  against S). End of run: classical 5.1 eV, T1 4.0 eV, T2 1.8 eV; T2 goes
  NEGATIVE for the first ~15 Bohr (var(p) growth outruns the drift loss).
- `09d_interactions_both` / `09e_..._projectile` — BOTH halves' pairwise terms on
  one axes. Hue = term, line style = representation (the REVERSE of the
  stopping panels: five distinct terms cannot collapse to two colours). The
  zoomed variant exists because E_SS/E_SB span +-15 eV and compress the
  projectile terms (+-3 eV) that actually differ. Classical dE_PP and dE_PB are
  absent BY CONSTRUCTION (rigid cloud, z-uniform background) — the absence is the
  physics, not a gap.
- `13d_classical_fit_vs_path` / `13e_classical_energy_vs_time` — the classical
  stopping power fitted over the same three windows. **13e deliberately draws NO
  regression line**: the fit is dE/ds and the projectile decelerates, so on a
  time axis the same fit is a curve; a straight line there would misrepresent it.
- `13c_stopping_bar` rebuilt. Three real defects fixed: raw LaTeX leaked into the
  tick labels (splitting the formula on "=" stripped the `$` delimiters); the fit
  sigmas are ~1e-4 so 3 decimals printed every uncertainty as "0.000"; the legend
  collided with a bar annotation.
- **Two-colour scheme**: classical = blue, wavepacket = red, everywhere the two
  projectiles appear; windows/estimators by LINE STYLE. Two palettes are
  deliberately exempt because hue encodes something else there — the momentum
  time-slices and the five pairwise interaction terms.
- **Legend geometry assertion** (`_assert_legend_fits`) now fails the build if any
  legend would be clipped. It caught `09b` (three of five entries were being cut
  off the canvas — invisible in code, matplotlib clips silently) and `12c`.
- Momentum maps: transverse binning refined to dk_perp = 0.08 (half a grid
  spacing) — MEASURED as the finest fully-populated binning (45 distinct lattice
  radii below 1.6 Bohr^-1; at 0.05 four bins are empty and the map combs).
  `kperp_max = 2.0` deliberately exceeds the 1.6 display limit because
  `kz_kperp_map` CLIPS overflow into the top bin, which at 1.6 would have made a
  fake spike inside the panel. k_z is untouched: it is the exact FFT grid,
  dk_z = 2pi/L_z = 0.105, a hard limit of the 60 Bohr box. Zero-padding was
  REJECTED — 2.8 % of the norm sits on a box face at t=30, so padding would
  impose a discontinuity. The finer binning resolves a NEW feature: a gain stripe
  at k_perp < 0.08 previously merged into the loss region.

Vacuum: `hypotheses/wp_selfinteraction/build_report_figures.py` ->
`report_figures/{sigma4,sweep}/` (6 + 4). Sweep panels use tau = t/sigma^2, the
only axis on which six runs of different duration are comparable.

### Open

- Sweep is analysed but has NO notebook yet (module + tests only).
- `at_fixed_physical_time` for sigma < 4 at T = 30 a.u. would need longer runs.
- Nothing committed to git.
- `wp_selfinteraction` notebook and `docs/handovers/wp-self-interaction-correction.md`
  are being edited by ANOTHER session concurrently — a notebook rebuild of mine was
  clobbered once. Coordinate before rebuilding it again.
- Nothing committed to git.

---

## 2026-08-02 (6th) — VACUUM MEASUREMENT DONE. It **REFUTES** the 5th entry's attribution.

> **Amended by the 7th entry:** the "16.9 %" below mixed a 3-D vacuum ratio with a
> transverse channeling ratio; the matched-definition value is **20.9 %**, since
> confirmed by a direct SIC run. The vacuum numbers themselves are unaffected.

**Read this before using the 2026-08-02 (5th) conclusions.** That entry
attributed the channeling packet's 1.467x excess expansion, and ~two thirds of
the stopping deficit, to LDA self-interaction. A direct vacuum measurement now
says **self-interaction accounts for only ~17 % of the excess spreading.**

Runs: `ResearchProject/systems/vacuum/scripts/wp_selfinteraction/results/{noninteracting,hartree,lda}`
(job 32615079, 854/898/975 s). Notebook:
`ResearchProject/systems/vacuum/hypotheses/wp_selfinteraction/selfinteraction.ipynb`
(13 cells, 0 errors). Analysis `selfinteraction.py`, 10 tests passing.

### Method — no SIC implementation needed

One electron alone in vacuum has EXACTLY no self-interaction, so the exact answer
is free dispersion. Propagating the SAME injected packet at three theory levels
(`non_interacting` / `hartree` / `lda`, all already in INQ) measures the error by
difference AND splits it into Hartree and xc parts.
sigma_WP = 4.0 (matched to channeling), k0 = 0, 72^3 Bohr box at h = 0.5,
1500 x 0.02 = 30 a.u.

### Gates — the reference is exact, so the difference is trustworthy

| gate | value |
|---|---|
| reference sigma vs analytic | max rel error **1.7e-12** |
| reference var(p) drift | **exactly 0** |
| reference E_total drift | **exactly 0** |
| closure `E_PP` vs INQ `energy_hartree` | **0.00e+00 Ha** (both interacting runs) |
| wrap indicator | 2.6e-9 (packet never reached the box face) |

### RESULT — the xc term cancels most of the Hartree self-interaction

| | sigma at t=30 | excess vs free | var(p_z) growth |
|---|---|---|---|
| non-interacting (exact) | 6.010 | 1.000 | **0.00 %** |
| Hartree self-interaction only | 8.634 | **1.436** (+43.6 %) | +141.7 % |
| LDA (Hartree + xc) | 6.484 | **1.079** (+7.9 %) | +17.2 % |
| **xc part alone** | | **-35.8 percentage points** | |

**LDA exchange-correlation cancels ~82 % of the Hartree self-repulsion.**

Exact energy accounting (a closed identity, not a fit): for the Hartree run with
k0 = 0, `E_total = var(p)/2m + E_PP` is conserved, so `d[var(p)/2m] = -dE_PP`.
Measured: E_PP released **1.80732 eV**, excess var(p)/2m gained **1.80732 eV** —
equal to all five figures. E_total drift -3.6e-8 eV.

### CONSEQUENCE 1 — the SIC variant question is settled, and the earlier lean was WRONG

Plan §2 offered SIC-H (Hartree only) and SIC-PZ (Hartree + xc). **SIC-H would be
badly wrong**: removing the Hartree self-term while leaving the xc self-term
would leave a large uncancelled ATTRACTIVE self-binding and make the packet
severely UNDER-spread. Full Perdew-Zunger is required. (An early smoke reading
that "hartree spreads more than lda" was the first sign of this; the production
run confirms it at 43.6 % vs 7.9 %.)

### CONSEQUENCE 2 — the channeling attribution in the 5th entry does NOT hold

| | excess spreading vs free |
|---|---|
| channeling production run (LDA, in a bore) | **1.467** |
| vacuum LDA self-interaction, same sigma and duration | **1.079** |
| fraction of the channeling excess explained | **16.9 %** |

So ~83 % of the channeling packet's excess expansion is NOT self-interaction.
Two candidates, NOT yet separated:
1. the environment — the bath and the tube wall acting on the packet;
2. **environment-dependent xc cancellation** — LDA's xc self-cancellation is
   evaluated at the TOTAL local density, which in the bore is packet + bath, not
   packet alone. The 82 % cancellation measured in vacuum need not hold there.

Candidate 2 means the vacuum number is **not** directly transferable and the
17 % must be read as an indication, not a decomposition. It is reported as a
ratio in `channeling_comparison()` for exactly this reason.

**Do not repeat the 5th entry's claim that SIE drives the channeling result until
this is settled.** The honest current statement: at sigma = 4 in vacuum, LDA
self-interaction costs ~8 % excess width and ~0.22 eV; the channeling run's 47 %
excess is mostly something else.

**Note (2026-08-02):** a separate in-flight job (32615189, `wp-si`, plus pending
`chan-sic` 32615190/1) is extending this run set with `sic_h` and `sic_pzrun`
variants — not launched from this thread; results not yet assessed here.

### Highest-value next step

The `sigma_WP` sweep, now cheap and unambiguous: `E_PP ~ 1/sigma`, so repeating
these three levels across sigma maps how the error scales and tests candidate 2
by varying the packet's own density relative to a fixed environment. Same binary,
one env var.

---

## 2026-08-02 (5th) — MECHANISM IDENTIFIED: the T1 deficit is dominated by LDA self-interaction error

Analysis only, no new runs. Answers the user's questions on the var(p) shape and
the classical/WP difference. **This changes how the headline number should be
reported** and is the most consequential finding of the study so far.

### var(p) has TWO components with DIFFERENT functional forms, and they mean different things

| component | shape | fit | growth |
|---|---|---|---|
| `var(p_z)/2m` | **LINEAR** | r2 = 0.9941 vs time, 0.9934 vs path; 2.75e-4 per Bohr | +0.189 eV (8.8 % of total) |
| `var(p_perp)/2m` | **SIGMOID**, peaks t = 24.28 | — | +1.950 eV (91.2 % of total) |

- **Linear** = diffusive, no finite reservoir. Wavepacket analogue of energy-loss
  STRAGGLING: the packet is 2.83 Bohr wide, samples a spatially varying wake, and
  the spread in impulse across it accumulates additively with path.
- **Sigmoid** = a FINITE reservoir being drained to exhaustion. The reservoir is
  the packet's own self-Hartree energy.

### The reservoir is E_PP — LDA self-interaction error. Three independent lines of evidence

1. **Energy balance.** E_PP released 1.6439 eV (1.9361 -> 0.2921 eV, 85 %
   consumed); `var(p_perp)/2m` gained 1.9496 eV. Ratio 1.19.
2. **Excess expansion.** Measured against the EXACT free-Gaussian Rayleigh
   baseline (`<r> = sd sqrt(pi/2)`, `sigma_r = sd sqrt(2-pi/2)`), which matches
   the t=0 packet to 0.1 %:

   | t | `<r_perp>` meas | free | ratio | `sigma_r` meas | free | ratio |
   |---|---|---|---|---|---|---|
   | 0 | 3.544 | 3.545 | 1.000 | 1.854 | 1.853 | 1.001 |
   | 20 | 7.511 | 5.675 | 1.324 | 4.533 | 2.966 | 1.528 |
   | 30 | 11.053 | 7.533 | **1.467** | 5.835 | 3.938 | **1.482** |

   The packet expands ~1.47x faster than free dispersion. Something pushes it apart.
3. **No classical counterpart.** The classical twin's E_PP is EXACTLY constant
   (spread 0.0 Ha) — a rigid cloud has a fixed self-energy.

**REJECTED alternative: the wall attracting the packet.** `f_outside` reaches
**0.414** by t=30 — 41 % of the packet passes THROUGH the wall into the outer
vacuum. It is not being trapped in a well; it is being blown through.
(`f_wall` peaks ~0.16 at t=20 then falls: transit, not capture.)

### Why var(p_perp) turns over near the end — two causes, one real, one not

- **Plateau is REAL**: E_PP is 85 % consumed, so the driving force has died.
- **The 3 % DECLINE after t ~ 24 is BOUNDARY CONTAMINATION**: `<r>+2 sigma_r`
  exceeds the 20 Bohr transverse half-box at **t = 25.06**; the peak is at
  **t = 24.28**. User independently confirmed transverse wrapping in the density
  GIF. Do not interpret anything transverse beyond t ~ 24.

### var(p_z) growth does NOT explain the T1 deficit — category error plus numbers

`T1 = <p>^2/2m` contains NO var(p), so var growth cannot be a component of the
T1 deficit by construction. Numerically: `d[var(p_z)/2m]` = 0.189 eV against a
1.183 eV gap = **16 %**; and the TOTAL var growth (2.139 eV) is **181 %** of the
gap, i.e. larger than the quantity it would supposedly decompose. The T1 deficit
is an IMPULSE deficit (76.4 % of classical impulse), established earlier.

### Energy actually delivered to the bath — the physically decisive comparison

| | bath electron KE gain |
|---|---|
| classical | +3.9422 eV |
| wavepacket | +2.5099 eV |
| ratio | **0.637** |

WP kinetic ledger: drift T1 -3.9425 eV, internal var(p)/2m +2.1387 eV, net
orbital KE (T2) -1.8038 eV. **54.2 % of the WP's drift energy loss was
re-absorbed by its own spreading instead of reaching the bath.**

### CONSEQUENCE FOR THE STUDY (important)

The ~20 % T1 deficit (S_wp/S_cl = 0.801 over 9-25) is **NOT a clean quantum
correction to stopping power**. It is substantially an artefact of LDA
self-interaction driving the packet apart, which then weakens its coupling
(impulse ratio 0.920 compact -> 0.764 spread, r = +0.98 vs f_bore).

**A KS-orbital stopping power measured this way carries an SIE contamination at
the tens-of-per-cent level.** This should be stated in any writeup.

**Discriminating test (NOT yet run):** a sigma_WP sweep. `E_PP ~ 1/sigma`, so if
the deficit scales as 1/sigma_WP it is SIE; if sigma-independent, it is genuine
quantum kinematics. This is the single highest-value follow-up and is cheap
(same three binaries, one env var).

---

## 2026-08-02 (4th) — WINDOWS CHOSEN BY USER. Stopping powers measured.

`refined_analysis.ipynb`: **29 cells, 2.58 MB, 0 errors**. 45 tests pass.
Results also in `refined_stopping_summary.csv`.

### The three user-chosen windows and their S (eV/Bohr)

Classical fitted over the SAME window as each WP row (both decelerate, so
cross-window comparison compares different velocities).

| window | est. | S_wp | +- | r2 | S_classical | +- | ratio | n |
|---|---|---|---|---|---|---|---|---|
| **9-25** | T1 | **0.08773** | 0.00009 | 0.9991 | 0.10956 | 0.00011 | **0.801** | 801 |
| **21-30** | T2 | **0.08386** | 0.00030 | 0.9943 | 0.13217 | 0.00021 | 0.634 | 451 |
| **5-20** | T2 | **0.01289** | 0.00007 | 0.9754 | 0.09791 | 0.00030 | 0.132 | 751 |

Path definition makes no difference (`int<p>dt` vs `centroid` agree to 0.2 %),
which is the Ehrenfest consistency check holding over all three windows.

Classical cross-check: fitting the projectile KE and the bath `E_total` gives the
same S to **5e-7 eV/Bohr** in every window — the energy-budget closure carried
through to the fitted slope.

### Reading

- **T1 over 9-25 is the trustworthy number.** r2 = 0.9991, and the ratio 0.801
  sits between the early-plateau impulse ratio (0.920) and the end-of-run value
  (0.764), as it must for a window spanning both regimes.
- **T2 is window-dependent by a factor of 6.5** (0.0129 vs 0.0839). This is the
  quantitative statement of why T2 is not a usable stopping estimator here: over
  5-20 the var(p) growth is at its steepest and cancels most of the drift loss;
  over 21-30 var(p_perp) has begun to saturate so the cancellation weakens. The
  measured spread of T2 across windows is itself the result.
- **T1 varies by only ~9 %** across comparable windows, vs 650 % for T2.

### Caveats printed in the notebook alongside the numbers

- `9-25` extends 2.9 a.u. past own-wake re-entry (t = 22.1); f_bore at t1 = 0.580.
- `21-30` extends 7.9 a.u. past it; f_bore at t1 = 0.457. The classical reference
  in this window (0.132) is inflated by the same own-wake effect, so the 0.634
  ratio is the least trustworthy of the three.
- `5-20` starts inside the launch transient (local S still rising until ~10).

### Also added — var(p)/2m on its own axes (section 2b, `03b_var_p_term.png`)

Total, plus split into longitudinal and transverse, plus growth from t=0.

| term | t=0 | t=30 | growth |
|---|---|---|---|
| var(p)/2m total | 1.2755 eV | 3.4142 eV | +2.1387 eV (+167.7 %) |
| var(p_z)/2m | 0.4252 eV | 0.6143 eV | +0.1891 eV (+44.5 %) |
| var(p_perp)/2m | 0.8504 eV | 2.8000 eV | +1.9496 eV (+229.3 %) |

**Only 8.8 % of the growth is longitudinal; 91.2 % is transverse.** So the
T2 - T1 growth is overwhelmingly the packet spreading SIDEWAYS (the channeling
premise failing) rather than longitudinal straggling along the stopping
direction (which is what would directly spoil -dT2/ds).

---

## 2026-08-02 (3rd) — Fit-target plots added (section 6a). Notebook 27 cells, 2.25 MB, 0 errors.

User asked for the two quantities they will actually regress, plotted alone, to
choose the window from: `delta E_total` (classical) and `T_1` (WP).

`refined_figs/08_fit_targets.png`, section **6a**, immediately above the window
parameter cell:
- (a) both raw, with signs as measured (`dE_total` positive = bath gains,
  `dT_1` negative = projectile loses);
- (b) both as ENERGY LOST BY THE PROJECTILE vs PATH — the actual fit;
- (c) **LOCAL slope `-dE/ds`**, centred +/-1 a.u. stencil. This is the
  window-selection tool: a good window is where the local slope is FLAT.
  Edges are left NaN rather than filled with a one-sided stencil.

### What the local slope shows (eV/Bohr)

| t | S_classical | S_wp(T1) | ratio | f_bore |
|---|---|---|---|---|
| 2 | 0.0038 | 0.0035 | 0.919 | 0.997 |
| 6 | 0.0477 | 0.0441 | 0.923 | 0.988 |
| 10 | 0.0935 | 0.0852 | 0.911 | 0.953 |
| 12 | 0.1038 | 0.0929 | 0.896 | 0.921 |
| 15 | 0.1084 | 0.0928 | 0.856 | 0.852 |
| 18 | 0.1091 | 0.0863 | 0.791 | 0.770 |
| 20 | 0.1114 | 0.0822 | 0.738 | 0.714 |
| 24 | 0.1264 | 0.0824 | 0.652 | 0.606 |
| 28 | 0.1410 | 0.0888 | 0.629 | 0.506 |

- **t < ~10 is a launch transient** — S rises from ~0, it is not a plateau.
- **classical plateau: t ~ 12-21, S = 0.108-0.111 eV/Bohr.**
- **WP plateau: t ~ 11-16, S = 0.092-0.094 eV/Bohr**, then it DECLINES.
- Neither is flat over the whole run.

### NEW HARD UPPER BOUND on any window — own-wake re-entry

The projectile re-approaches its OWN launch-time disturbance through the
periodic image. Measured distances (launch z = -28, i.e. image at +32):

| t | z | distance to own launch point via image |
|---|---|---|
| 20 | 9.97 | 22.03 |
| 24 | 17.39 | 14.61 |
| 28 | 24.74 | 7.26 |
| 30 | 28.38 | **3.62** |

The wake is `lambda_p = 2 pi v / omega_p = 36.1 Bohr` long in a **60 Bohr** box,
so it fills more than half the cell. Crossing `lambda_p/2` happens at
**t = 22.1 a.u.**; within `3 sigma_pot = 8.5 Bohr` at **t = 27.4 a.u.**

**This is the most likely cause of the late-time RISE in the classical local S
(0.111 -> 0.141 after t ~ 21), which is NOT physical stopping.** Shaded red on
panel (c) and printed in the window-reference block. Any window should end
before ~22 a.u.

Inference, not established: the attribution of the rise to own-wake re-entry is
consistent with the timing but has not been proven — a longer box (L_z >= 100)
would settle it and is not in scope here.

---

## 2026-08-02 (later) — Refined notebook EXTENDED: 2-D momentum map, combined coupling, and the T1-vs-classical explanation

**Status: DELIVERED.** No new simulation. `refined_analysis.ipynb` is now
**26 cells, 2.03 MB, 0 errors**. Tests: **44 passed** (channeling_twin) +
**7 passed** (new inqview kz_kperp_map suite).

### New library capability

`inqview.visualisation.field_io.kz_kperp_map(field)` -> `(k_z, k_perp, P)`, the
2-D momentum density of one orbital, reconstructed from the COMPLEX orbital VTI
dumps (`raw/vti/wavefunction_wp/`, 11 frames every 150 steps).
Tests: `inq-stack/tests/python/inqview/visualisation/test_kz_kperp_map.py`.

- `k_z` is the NATIVE FFT axis and is **not binned** -> its moments are exact.
  Only `k_perp` is binned, at one transverse grid spacing.
- **Do not take transverse moments from the map.** Bin-centre assignment biases
  `<k_perp^2>` HIGH by +6.3 % at `sigma_p/dk = 3.0`, +9.4 % at 2.4. Verified the
  RAW unbinned moment is exact to 0.00 % on four grids, so this is bin-centre
  assignment alone. RMS bin centres were tried and made it WORSE (+7.6 %).
  Use `wp_momentum_stats.csv` for numbers; use the map for SHAPE and DIFFERENCES.
- Validated: the t=0 dump round-trips to `<k_z>=1.917011`, `var(k_z)=0.031250`,
  `var(k_perp)=0.062500`, matching `wp_momentum_stats.csv` exactly.

### Q1 ANSWERED — the momentum-histogram spikes are an ARTEFACT, not physics

Built the analytic launched Gaussian on the identical grid with ZERO
interaction, put it through the identical binning:

| | roughness |
|---|---|
| analytic non-interacting Gaussian | 2.685 |
| measured t=0 | 2.648 |

correlation of the two histograms **r = 0.99982**. Mechanism is NOT empty bins
(every bin holds >= 368 grid points; that hypothesis was tested and rejected).
It is that the packet occupies only ~**95 effective k-points** of 768 000
(participation ratio; top 50 carry 63 % of the norm), binned into 128 shells.
**Fix for future runs: `n_bins ~ 40`, not 128.**

### Q2 ANSWERED — the impulse is NOT the same for every momentum channel

A rigid shift preserves the SHAPE. It does not:

| t (a.u.) | `<k_z>` | `sigma_kz` | skewness | norm above mean |
|---|---|---|---|---|
| 0 | 1.91701 | 0.17678 | -0.00000 | 0.4987 |
| 15 | 1.88706 | 0.19138 | +0.01709 | 0.4990 |
| 30 | 1.83988 | 0.21248 | +0.15916 | 0.4881 |

`sigma_kz` +20.2 %, skewness 0 -> +0.159. Norm above the mean DOES fall,
0.4987 -> 0.4881. Norm transferred across `k0`: **0.1418** moved below.

**Metric caveat that had to be fixed first:** a hard `k_z > mean` count returns
**0.454** for the t=0 packet, which is EXACTLY symmetric, because only ~8
resolved `k_z` points carry it. `refined.kz_asymmetry` interpolates the CDF
instead and returns 0.4987. Pinned by a test.

### Q3 ANSWERED — the T1-vs-classical gap is an IMPULSE gap, in TWO regimes

`T1 = <p>^2/2m` exactly and both halves start at the same p, so the energy gap
is algebraically a `delta-p` gap. WP received **76.4 %** of the classical
impulse (`dp`: -0.0771 vs -0.1009).

| regime | t (a.u.) | impulse ratio | cause |
|---|---|---|---|
| EARLY plateau | 2-10 | **0.9204 +- 0.0014** | NOT spreading (sigma_z 2.86->3.45, f_bore >= 0.95). ~8 % deficit present from the start. |
| LATE | >10 | falls to **0.7643** | tracks delocalisation |

Correlation of the impulse ratio with spreading proxies over t in [2,28]:
`f_bore` **r = +0.9805**, `sigma_z(0)/sigma_z(t)` +0.9015, `E_PP/E_PP(0)` +0.8804.

The ~8 % early deficit is NOT explained by spreading and is **unattributed**.
Candidates (inference, not established): XC self-interaction (the WP contributes
to `n` so it feels its own `v_xc`; the classical external potential does not);
Pauli/orthogonality against the 80 occupied bath orbitals. **Discriminating test
not yet run** — a `sigma_WP` sweep would separate a form-factor effect from an
orbital-nature effect.

**The gap growth rate is NOT constant** (contrary to a first reading of the
plot) — it accelerates 9x:

| t | d(gap)/dt (eV per a.u.) |
|---|---|
| 5-10 | 0.0103 |
| 10-15 | 0.0219 |
| 15-20 | 0.0403 |
| 20-25 | 0.0704 |
| 25-30 | 0.0923 |

### Combined coupling delta(E_PS + E_PB) — user's suspicion CONFIRMED

`E_PS + E_PB = integral n_P (phi_S - phi_+) = integral n_P phi_(S+B)` — the
projectile's coupling to the NET (neutral) charge density, hence gauge-clean.
Split apart, each term is the potential of a CHARGED subsystem.

| quantity | worst WP-vs-classical disagreement |
|---|---|
| `dE_PS` alone | 2.6753 eV |
| `dE_PB` alone | 2.9124 eV |
| **`dE_PS + dE_PB`** | **0.2416 eV** (11.1x smaller) |

Classical `dE_PB` is IDENTICALLY ZERO: the tube is uniform in z and the cloud is
rigid, so a rigid cloud translating in z sees a constant background potential.
**WP `dE_PB` is therefore a pure TRANSVERSE-spreading signal.**

### Corrections to earlier readings (recorded so they are not repeated)

- "E_SB / E_SS almost identical in both halves" — same SHAPE, but max
  |wp - classical| is **4.06 eV** (`dE_SS`) and **4.45 eV** (`dE_SB`), i.e.
  **~29 % of amplitude**, and larger than the entire 5.1 eV stopping signal.
  Not negligible. `E_BB` IS exactly constant (52.081389 Ha) in both.
- "self-Hartree becomes more negative" — `dE_PP` is negative (-1.64 eV) but
  `E_PP` itself stays **positive throughout** (1.9361 -> 0.2921 eV). The
  self-Hartree of a positive-definite density cannot be negative.

### Files changed

| Path | Change |
|---|---|
| `inq-stack/python/inqview/visualisation/field_io.py` | + `kz_kperp_map` |
| `inq-stack/tests/python/inqview/visualisation/test_kz_kperp_map.py` | NEW, 7 tests |
| `.../channeling_twin/refined.py` | + `momentum_map`, `available_wf_steps`, `kz_asymmetry`, `impulse_comparison`, `combined_projectile_coupling` |
| `.../channeling_twin/tests/test_refined.py` | 15 -> 23 tests |
| `.../channeling_twin/build_refined_notebook.py` | + sections 3b, 4b, 5; window renumbered to 6 |
| `.../channeling_twin/refined_figs/` | 6 -> 10 PNGs |

**Still awaiting the user's window choice** (section 6, `T_WIN_CL` / `T_WIN_WP`).

---

## 2026-08-02 — Refined-analysis notebook (lightweight). User chooses the window.

**Status: DELIVERED, awaiting the user's window choice.** No new simulation was
run; this is analysis of the completed 2026-08-01 twin.

### Why

The user's judgement is that the fit window should be read off the diagnostics BY
EYE and chosen SEPARATELY for each half, rather than derived by the analysis
(`f_bore >= 0.95`, which gave 1.5-10.2 a.u.). This notebook shows the diagnostics
and renders NO verdict.

### Files (all absolute)

| Path | What |
|---|---|
| `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/ResearchProject/systems/cylindrical_jellium/hypotheses/channeling_twin/refined.py` | data layer, user T1/T2 convention |
| `.../channeling_twin/build_refined_notebook.py` | builder (18 cells) |
| `.../channeling_twin/refined_analysis.ipynb` | **the deliverable, 1.21 MB, 0 errors** |
| `.../channeling_twin/refined_figs/` | 6 PNGs |
| `.../channeling_twin/tests/test_refined.py` | 15 tests |
| `.../channeling_twin/tests/test_refined_notebook_cells.py` | 8 static guards |

`venv/bin/python -m pytest tests/ -q` from the `channeling_twin` dir: **36 passed**.

Rebuild:
```bash
cd .../hypotheses/channeling_twin
PYTHONPATH=<repo>/inq-stack/python <repo>/venv/bin/python build_refined_notebook.py
```

### THE LABEL SWAP — read before using this notebook

The user's convention, used throughout `refined.py` and the notebook, is the
REVERSE of `ks_stopping.py`:

| here | definition | `ks_stopping.py` |
|---|---|---|
| `T1_drift_ev` | `<p>^2/2m` (drift only) | `T2` |
| `T2_total_ev` | `<p>^2/2m + var(p)/2m = <p^2>/2m` | `T1` |

Same quantities, swapped names. Reading a curve under the wrong convention
INVERTS the study's conclusion. `T2 - T1 = var(p)/2m` is exact by construction
and is asserted in the tests.

### Three defects found and fixed (all pre-existing, none invalidates a result)

1. **`momentum_distribution.csv` was being destroyed on load.**
   `ks_stopping._concat_segments` ends with
   `drop_duplicates(subset="step", keep="last")` — right for a scalar-per-step
   observable, catastrophic for this LONG-FORMAT file which has 128 rows per step
   (one per k-bin). It kept ONE bin per step and silently discarded 127. Symptom
   was not an exception: a plausible table where every distribution collapsed to
   the Nyquist bin and integrated to 0. Fixed by a dedicated loader keying on
   `(step, k)` plus a raggedness guard that raises. Regression tests on both
   synthetic and production data.

2. **`projectile.csv` column `proj_z_unwrapped` is written one step ahead.**
   Measured: `proj_z_unwrapped[i] == proj_z[i+1]` EXACTLY for all 1501 rows, and
   `proj_z_unwrapped[0] = -27.9617` while the launch point is `-28.0`.
   Velocity-Verlet consistency against the recorded `proj_vz`: `proj_z` closes to
   1.2e-9, `proj_z_unwrapped` to only 2.1e-6. It is the POST-update position
   written against the PRE-update step index.
   **Impact on published results: NONE.** The offset is a CONSTANT 0.038 Bohr
   (= v dt), so it moves a linear fit's intercept, not its slope; and this run has
   `n_wraps == 0` throughout so both columns carry the same trajectory.
   `refined.cl_frame` now derives the path from `proj_z` via `unwrap_periodic`.
   **NOT fixed in `run.cpp`** — flagged here for the next touch of
   `scripts/channeling_twin/classical/run.cpp`.

3. **My own test suite had an isolation leak** (worth recording because it nearly
   hid #2): `CHAN_WP_RESULTS`/`CHAN_CL_RESULTS` are read ONCE into module
   constants at `channeling_stopping` import time, so a cached module keeps
   whichever tree it first saw and `monkeypatch` restoring the env at teardown
   does NOT undo it. The production tests were silently running against the
   synthetic tmp_path tree. Fixed with a `_fresh_refined()` helper that evicts
   both modules from `sys.modules`.

### New physics surfaced by this notebook

**1. The twin contract is verified electrostatically, not just by config parity.**
At t = 0 EVERY pairwise term is identical between the halves:

| term | wp (Ha) | classical (Ha) | \|diff\| |
|---|---|---|---|
| E_SS | +44.1893954886 | +44.1893954886 | 0.0 |
| E_PS | +0.2883777942 | +0.2883777942 | 6.0e-12 |
| E_PP | +0.0711489975 | +0.0711489975 | 2.1e-12 |
| E_SB | -95.4759736723 | -95.4759736723 | 0.0 |
| E_PB | -0.2651274260 | -0.2651274260 | 8.0e-12 |

Worst disagreement 8.0e-12 Ha = 2.2e-10 eV. This is the direct evidence that
`sigma_pot = sigma_WP/sqrt(2) = 2.828427` Bohr makes the classical Gaussian cloud
reproduce the wavepacket's own charge density — i.e. the answer to the user's
original "ensure they have same potentials created". Config parity
(`twin_manifest.json`) checks the INPUTS; this checks the resulting FIELDS.

**2. The classical energy budget closes to 2.2e-5 eV over 1501 steps.**
`energy_ion` and `energy_ion_kinetic` are identically 0 (the projectile is an
external perturbation, not an INQ ion), so INQ's `energy_total` is the BATH
alone. dE_bath = +5.1256 eV, d(1/2 m v^2) = -5.1256 eV.

**3. The var(p) growth is mostly TRANSVERSE, which changes its meaning.**

| quantity | t=0 | t=30 | growth |
|---|---|---|---|
| var(p_z) | 0.03125 | 0.04515 | +44.5 % |
| var(p_perp) | 0.06250 | 0.20579 | +229.3 % |
| var(p) 3D | 0.09375 | 0.25094 | +167.7 % |
| T2 - T1 | 1.2755 eV | 3.4142 eV | +167.7 % |

Longitudinal straggling contaminates `-dT2/ds` directly; transverse spreading is
instead a failure of the CHANNELING premise (the packet reaching for the wall).
Most of the growth is the latter. `var(p_perp)` also PEAKS near t ~ 24 a.u. and
then declines slightly.

**4. The two WP position definitions agree far longer than expected.**
`|s_centroid - s_pintegral|` is 4.4e-4 Bohr at t=10, 7.4e-3 at t=20, 4.3e-2 at
t=30. The Ehrenfest picture holds well past the `f_bore` breach at 10.24 a.u.
WP path 56.566 Bohr vs classical 56.380 Bohr — the WP ends 0.186 Bohr AHEAD.

**5. Whole-run energy loss**: classical 5.1256 eV, WP `T1` 3.9425 eV, WP `T2`
1.8038 eV. These are whole-run averages over a decelerating projectile, NOT
S(v0) (`.claude/rules/light-projectile-stopping.md`).

**6. The momentum histogram is resolution-limited and only qualitative.**
sigma_p = 0.177 a.u. vs k-grid spacing 2*pi/60 = 0.105 (z) and 2*pi/40 = 0.157
(xy) — the distribution is 1-2 grid points wide, so the radial |k| histogram is
spiky by construction. Confirmed NOT shell-binning geometry (correlation of the
jitter with shell-occupancy counts r = -0.04). The MOMENTS in
`wp_momentum_stats.csv` are exact grid expectation values and remain
quantitative; the notebook says so and directs the reader to the log panel.

### Window-choice reference points (printed by the notebook)

| marker | t (a.u.) |
|---|---|
| run length | 30.00 |
| plasma period 2*pi/omega_p | 18.85 |
| f_bore < 0.99 | 5.42 |
| var(p_z) > 1.05x free | 6.98 |
| f_bore < 0.95 | 10.24 |
| \|Ehrenfest residual\| > 1e-3 Bohr | 12.24 |

### NEXT STEP — the user's

Edit the single parameter cell in section 5 of `refined_analysis.ipynb`:

```python
T_WIN_CL = (t0, t1)   # classical, a.u.
T_WIN_WP = (t0, t1)   # wavepacket, a.u.
```

Both default to `None`, and while unset the notebook fits NOTHING. Once set it
fits classical `dKE/ds` and the four WP combinations (T1/T2 x
s_centroid/s_pintegral) through one code path (`refined.fit_in_window`) and
writes `refined_stopping_summary.csv`.

### Deliberate departure from a rule

`.claude/rules/notebook-density-gif.md` requires a density-matrix GIF on every
run/analysis notebook. This notebook LINKS to the existing GIFs
(`comparison_figs/density_matrix/`, already embedded in
`channeling_twin_comparison.ipynb`) instead of re-embedding them. Re-embedding
would add ~222 MB and no information, and "lightweight" was the explicit request.
Stated in the notebook header, not silent.

### Still outstanding from 2026-08-01 (unchanged)

- **/rds is nearly full.** ~5.0 GB of throwaway smoke checkpoints under
  `scripts/channeling_twin/{wp,classical}/results/smoke/` are safe to delete.
  **Nothing has been deleted; awaiting the user.**
- `channeling_check`'s `t_breach > 0.5 * t_end` criterion is run-length relative
  and flips as a run is extended. Should be redefined before any follow-up run.
- Possible follow-up at `R_in ~ 15`, `R_out ~ 19`, `L_xy >= 50` (not requested,
  needs approval).

---

## 2026-08-01 — COMPLETE. Verdict: **AIM PARTLY MET**. The bore was too narrow.

### Everything ran; all deliverables exist

`twin_manifest.json`: **valid = true**, dynamic, parity checked on
`periodicity, lz, spacing, n, sigma_wp, launch_z, gs_dir`; projectile differs
(`wavepacket_orbital` vs `gaussian_charge_perturbation`). Comparison notebook:
23 cells, **0 errors**, 9 inline figures, 9 embedded density GIFs. Plus the two
per-run notebooks and `stopping_summary.csv`.

### Result

| estimator | half | S (eV/Bohr) | uncert. | r² | n | window (a.u.) |
|---|---|---|---|---|---|---|
| S_13 ⟨p²⟩/2m vs centroid | wp | 0.0045 | 0.0086 | 0.430 | 436 | 1.5–10.2 |
| S_14 ⟨p²⟩/2m vs ∫⟨p⟩dt | wp | 0.0045 | 0.0086 | 0.430 | 436 | 1.5–10.2 |
| S_23 ⟨p⟩²/2m vs centroid | wp | 0.0427 | 0.0276 | 0.912 | 436 | 1.5–10.2 |
| **S_24 ⟨p⟩²/2m vs ∫⟨p⟩dt** | wp | **0.0427** | 0.0275 | 0.912 | 436 | 1.5–10.2 |
| **S_cl same window** | classical | **0.0464** | 0.0303 | 0.912 | 436 | 1.5–10.2 |
| S_cl initial drag | classical | 0.0971 | 0.0024 | 0.970 | 1501 | 0–30 |

**Correctness gates (excellent):** WP energy drift **−7.8e-07 eV** over 1500 steps
(gate 1e-3); WP norm drift 2.9e-09; classical conserved quantity flat to 0.0000 eV;
classical off-axis excursion max **4.3e-04 Bohr** over 56 Bohr travelled, |F_x| ≤
8.7e-07 Ha/Bohr — the tube's symmetry holds on a square grid, MEASURED not assumed.

### Why "partly"

- **RESULT — consistent but weak.** S_24 vs S_cl differ by 8 %, but the
  uncertainties are ~65 % of the values. This is *consistent with* the aim, not
  evidence for it.
- **PREMISE — FAILED.** f_bore 0.998 → **0.457**; ⟨r⊥⟩ 3.54 → **11.05 Bohr**, past
  R_in = 10. Breach at t = 10.2 a.u. vs **16.7 a.u.** predicted by free dispersion
  for the same f_bore = 0.95 criterion (the 23.3 in §3.2 is the looser 2σ_d = R_in
  test) — so the wall is ACTIVELY BROADENING the packet, not merely letting it
  disperse.
- **MECHANISM — partial.** var(p_z) +9.0 % inside the window (inside the 10 %
  tolerance, but only just) and **+44.5 % over the full run**.

**The one clean physics result:** S_13/S_14 = 0.0045 ± 0.0086 — consistent with
ZERO, r² = 0.43 — while S_23/S_24 fit at r² = 0.912. The T₁-based definitions
return no stopping at all because the growing var(p) cancels the drift loss. That
reproduces the bulk contamination in a new geometry and confirms ⟨p⟩²/2m as the
defensible definition.

### CAVEAT THAT LIMITS THE WHOLE RUN (under-weighted when the design was locked)

The channeling window is 8.7 a.u. = **0.46 plasma periods** (T_p = 18.85 a.u.).
The classical S is 0.046 in that window but 0.097 over the full run — a factor of
2. *Inference:* the wake has not formed within one plasma period, so NEITHER number
is a converged stopping power. The comparison is like-for-like (same window, same
transient) but compares two TRANSIENTS.

### THE DESIGN FIX, quantitatively (my derivation — check before relying on it)

σ_d(t) = √(σ²/2 + t²/2σ²) is minimised over σ at σ² = t, giving σ_d(T) = √T.
f_bore ≥ 0.95 requires σ_d ≤ R_in/2.448. Therefore the **maximum channeling time
for a bore is T_max = (R_in/2.448)², at the best possible σ_WP**:

| R_in (Bohr) | T_max (a.u.) | plasma periods |
|---|---|---|
| 10 (this run) | 16.7 | **0.89 — impossible at ANY σ_WP** |
| 10.6 | 18.8 | 1.0 |
| **15** | **37.6** | **2.0** |

**σ_WP = 4 was not the problem; the 10-Bohr bore was.** It could never have held a
packet for one plasma period. A follow-up wants **R_in ≈ 15, R_out ≈ 19**, keeping
L_xy ≥ 50 to preserve transverse-image clearance. Everything else (binaries,
analysis, notebooks) is parameterised and needs no code change.

### Storage — and a defect worth fixing before any re-run

Production output came in at **11.0 GB, exactly as budgeted** (1.2 GS + 4.9 WP +
4.0 classical + 0.9 notebooks). But the total was 17 GB: the **20-step SMOKE stage
wrote full-size checkpoints** (`CH_CKPT_EVERY=10`, 1.2 GB each), leaving **5.1 GB
of throwaway checkpoints** across the two smokes — 30 % of the study's footprint.
FIXED in `run-chan-twin.slurm`: the smoke now writes no interior checkpoint and
deletes its checkpoints once the gates pass.

**/rds is now at 8 GB free of a 1.1 TB quota.** The 5.0 GB in the existing
`results/smoke/` checkpoint dirs is recoverable and safe to delete (the smokes were
gates; the production runs are complete and independent) — **NOT deleted, awaiting
the user**. The comparison notebook is 222 MB because the 9 density GIFs are
base64-embedded per `.claude/rules/notebook-density-gif.md`.

### Not done

- No journal entry yet (needs the user's own observations).
- Run catalogue (`tddft-run-catalogue` skill) not updated.
- The `channeling_check` criterion (`t_breach > 0.5 × t_end`) is run-length
  relative, so it flips as a run is extended even though the measurement does not
  change. Poor definition; **deliberately NOT changed after seeing the data**.
  Fix it before the next run, not after this one.

---

## 2026-08-01 (night, later) — ALL GATES GREEN; BOTH production halves running

### Full library suite: PASS (job 32579688)

```
PURE   test_projectile ............................. 10144 assertions, 7 cases
ENGINE test_radial_occupancy_engine ................ Passed  2.96 s
ENGINE test_inq_field_export_engine ................ Passed  2.55 s
ENGINE test_wp_minimum_image_engine ................ Passed  3.37 s
ENGINE test_projectile_force_minimum_image_engine .. Passed  2.34 s
ENGINE test_gaussian_minimum_image_engine .......... Passed  2.09 s
ENGINE test_slab_occupancy_engine .................. Passed  2.56 s
```

### The WP injection fix, verified two independent ways

**(1) In isolation** — `test_wp_minimum_image_engine` passes all four cases,
including the deliberate contrast (flag OFF must FAIL the analytic moments) and
the no-op-at-box-centre backwards-compatibility check.

**(2) In production** — all NINE t=0 gates pass, most to 8–9 significant figures:

| gate | before (clipped) | after | expected |
|---|---|---|---|
| ⟨p_z⟩ | 1.882086 (−1.8 %) | **1.91701127** (−1.8e-9 %) | 1.91701127 |
| var(p_z) | 0.472693 (**+1413 %**) | **0.03125000007** (+2.2e-7 %) | 0.03125 |
| T₁ | 2.038748 (+8.2 %) | **1.884341105** (−1.6e-9 %) | 1.884341105 |
| T₁−T₂ | 0.267624 (+471 %) | **0.04687500004** (+7.5e-8 %) | 0.046875 |
| centroid z (circ) | −26.966476 | **−28** exactly | −28 |
| σ_z (circ) | 2.132770 (−24.6 %) | **2.828427125** (+1.9e-8 %) | 2.828427125 |

### Defect 6 (test-only) — my contrast assertion was wrong, and the reason matters

`test_wp_minimum_image_engine` first failed on **my own** claim
`CHECK(off.norm == Approx(1.0))` → measured **0.856**. Reading `wavepacket.hpp`:
renormalisation happens ONLY inside the `do_ortho_` branch ("No orthogonalisation
requested: … the packet is left exactly as constructed (not even renormalised)").

So: **without** ortho a truncated packet keeps a visibly depressed norm; **with**
ortho — which every production run uses — `scale = 1/√(res[0])` divides the missing
17 % away and the norm reads 1.000 while every momentum observable stays corrupted.
**The norm is structurally incapable of catching this failure; var(p_z) is not.**
The test now pins BOTH behaviours, because the second is what makes the first
dangerous.

### 20-step smoke readings (NOT the result — the run is 1500 steps)

```
[PASS] energy conservation: E_total drift = -1.76e-08 eV over 20 steps (gate 1e-3)
[info] var(p_z): 0.03125000007 -> 0.03125706433  (+0.023 %)
[info] f_bore:   0.99795       -> 0.99793        (min 0.99793)
```
Both mechanism indicators point the right way, but 20 steps is 0.4 a.u. of 30 and
dispersion is cumulative. Do not quote these as the finding.

### Running now

```
32575811 chan-cl  RUNNING (classical production, 19 min in)
32577647 chan-wp  RUNNING (wavepacket production)
32577648 chan-nb  PENDING afterany(both)
```

Both halves concurrent, on the SAME ground state, with every library kernel they
depend on tested. The next results are the science.

---

## 2026-08-01 (night) — GS PASSED, classical half RUNNING, and the t=0 gates caught a real library bug in the WP injection

### Ground state: PASSED all three gates in 6 min (job 32575807)

```
E_GS            = -61.9108203054 Ha
integral n dV   = 160.000     [PASS] exact neutrality  (159.99999999997)
num_states      = 104         [PASS] matches the RT binaries
bore depletion  = 0.1294      [PASS] a channel exists  (gate < 0.5)
```
| region | electrons | volume (Bohr³) | n̄ |
|---|---|---|---|
| bore r⊥ < 10 | 16.02 | 18675 | 8.58e-4 |
| wall 10 ≤ r⊥ < 14 | 119.73 | 18060 | 6.63e-3 |

The bore sits at **13 % of the wall density** — a real vacuum channel. r_s = 3.000000.
Checkpoint 1.2 GB (estimate was 1.28 GB, so the disk budget holds). **6 minutes, not
the 8 h allotted** — the earlier cost estimates were badly pessimistic.

### Classical half: smoke PASSED, production RUNNING (32575811)

```
[PASS] Ehrenfest conservation: drift of (E_electronic + KE_proj + U_proj_bg)
       = 1.24e-07 eV over 20 steps   (gate < 0.05 eV)
```
~400 000× inside tolerance: the Hellmann–Feynman force and the moving-Gaussian
perturbation agree, i.e. the minimum-image FORCE fix works in production.

### Defect 5 — THE WAVEPACKET WAS BEING TRUNCATED AT INJECTION (library bug)

The WP smoke aborted on its own t=0 gates. The pattern is the diagnosis:

```
[PASS] norm (real space)            1.000
[PASS] transverse std = s/sqrt2     2.828594  vs 2.828427   (+0.006 %)   <- x,y PERFECT
[PASS] f_bore(0) (Rayleigh)         0.997908  vs 0.998070
[PASS] <r_perp>(0)                  3.544236  vs 3.544908
[FAIL] <p_z> = k0                   1.882086  vs 1.917011   (-1.8 %)
[FAIL] sigma_pz^2 = 1/(2 s^2)       0.472693  vs 0.031250   (+1413 %)
[FAIL] T1                           2.038748  vs 1.884341   (+8.2 %)
[FAIL] T1 - T2 = 3/(4 s^2)          0.267624  vs 0.046875   (+471 %)
[FAIL] centroid z (CIRCULAR)      -26.966476  vs -28        (+1.03 Bohr)
[FAIL] spread sigma_z (CIRCULAR)    2.132770  vs 2.828427   (-24.6 %)
```

**Every failure is in z; x and y are perfect to 0.006 %.** That is the signature of
a packet CLIPPED at the −z face, and every number matches truncation quantitatively:
cutting the tail nearest the face pulls the centroid *into* the box (+1.03 Bohr;
the truncated-Gaussian prediction is +1.16), leaves a narrower packet (−24.6 %),
and — because a truncated Gaussian has a STEP DISCONTINUITY in real space, which is
broadband in momentum — inflates var(p_z) **fifteen-fold**.

**Root cause, and it is a claim I asserted repeatedly and wrongly.** Both run
headers and the plan said "the wavepacket has no such problem — a KS orbital lives
on a plain 3-D FFT basis and wraps exactly". That is true of the **propagation**
and false of the **injection**: `wavepacket.hpp` built its Gaussian from a plain
Cartesian displacement (`double dz_ = rz - bz;`), and `inject_into_last_extra_state`
NORMALISES afterwards, so the norm stayed 1.0 and hid the truncation completely.

**Fix — new library option `WavePacket::minimum_image(bool)`** (defaulted false, so
no published run changes). It folds the separation into [−L/2, L/2) per lattice
direction. The PHASE uses the same wrapped displacement: wrapping the amplitude
while leaving `exp(i k·r)` puts a jump of `exp(i k·L)` across the seam whenever
k·L is not a multiple of 2π (here k₀L_z = 115.0 rad = 18.3 × 2π → a 1.9 rad
discontinuity). `exp(i k·d)` differs from the old form only by the global constant
`exp(i k·b)`.

**This also matters for the TWIN specifically:** the classical half already wraps
its charge (`gaussian_density_minimum_image`), so a clipped packet is not its twin
at exactly the boundary the study introduces on purpose.

**Test** (`inq-stack/tests/include/inqkit/wavepacket/test_wp_minimum_image_engine.cpp`):
analytic moments with the flag ON; the contrast case with it OFF asserted to FAIL
them (a test that only checks the fixed path would pass if the flag did nothing);
and a no-op check at the box centre, which is the backwards-compatibility guarantee.

### Chain v5

```
32575811 chan-cl      RUNNING (production, classical)
32575976 chan-tests2  RUNNING (from_inq_field)
32577644 chan-tests3  (adds test_wp_minimum_image_engine)
32577645 chan-wp-build -> 32577646 chan-wp-smoke -> 32577647 chan-wp
32577648 chan-nb      afterany(chan-wp, chan-cl)
```

### For anyone reusing `inqkit::WavePacket`

Any run launching a packet within ~2 density σ of a cell face and NOT passing
`.minimum_image(true)` has been measuring a truncated packet with an inflated
var(p). Runs that launch near the box centre are unaffected (the flag is a no-op
there, asserted by the test).

---

## 2026-08-01 (late) — all library tests GREEN; first real compile error found and fixed; chain v4 running

### Library gate: ALL PASS (job 32575316, 2 min)

```
PURE   test_projectile ............................. 10144 assertions, 7 cases
ENGINE test_radial_occupancy_engine ................ Passed  2.78 s
ENGINE test_projectile_force_minimum_image_engine .. Passed  2.20 s   <- the L=24 fix
ENGINE test_gaussian_minimum_image_engine .......... Passed  1.98 s
ENGINE test_slab_occupancy_engine .................. Passed  2.48 s
exit status : 0  (ALL GATES PASSED)
```

### Defect 3 — the smoke stage could never have compiled anything (my design bug)

The smokes were scheduled `afterok(tests)` only, on the stated reasoning that they
"build the binaries so a compile error surfaces immediately". But
`run-chan-twin.slurm` has its OWN ground-state guard **above** the `inq-run` call,
so both smokes exited in **5 seconds without compiling** and their `afterok`
cancelled both production halves. The comment beside the guard asserted the
opposite of what the code did.

**Fix:** a third stage, `build`, which skips the GS guard, runs `inq-run` (whose
binary then exits 2 on the missing checkpoint — expected and ignored), and
succeeds iff `./run` exists and is newer than `run.cpp`. It runs **concurrently
with the ground state**, so a compile error costs minutes, not an 8-hour SCF.
Stages are now `build` → `smoke` (20 steps + t=0 gates, needs the GS) → `prod`.

### Defect 4 — the FIRST REAL COMPILE ERROR, and a genuine library gap behind it

```
gs/run.cpp(222): error: no instance of overloaded function
  "inqkit::io::RealField3DWriter::write" matches the argument list
  argument types are: (inq::basis::field<inq::basis::real_space, double>, const char [19])
```

`RealField3DWriter::write` takes an `inqkit::fields::RealField3D`, and the only
routes into that type were `fields::density::total(electrons)` and
`fields::orbital(electrons, i)` — both requiring an `electrons` object. **Every
other real-space field the wrapper builds** (the jellium background n₊, φ₊, a
projectile charge blob, a drag field) had **no way to be visualised at all**. The
GS hit it trying to dump its own annular background beside the electron density.

**Fix — NEW library function** `inqkit::fields::from_inq_field`
(`inq-stack/include/inqkit/fields/inq_field.hpp`): exports an arbitrary INQ
real-space field to a RealField3D, doing the FFT-natural → physical shift in ONE
place so a run script can never re-derive it wrongly (the bug class
`.claude/rules/vti-coordinate-mapping.md` exists for). `density::total` is
deliberately left untouched for now; the two should be unified later.

**Test** (`inq-stack/tests/include/inqkit/fields/test_inq_field_export_engine.cpp`,
registered + submitted as job 32575976): metadata vs the basis; a bump at the box
CENTRE must land at the MIDDLE index and not index 0 (the "slab at the edges"
picture); an off-centre bump at its physical index; and the decisive one —
`from_inq_field(electrons.density())` must be **bit-identical** to
`density::total(electrons)`, pinning the new code against the implementation every
production run already exercises.

### Chain v4 (running)

```
32575390 chan-wp-build  32575391 chan-cl-build   (compiling, concurrent with the GS)
32575807 chan-gs        (resubmitted with the from_inq_field fix)
32575808 chan-wp-smoke  32575809 chan-cl-smoke   afterok(build, gs)
32575810 chan-wp        32575811 chan-cl         afterok(smoke)
32575812 chan-nb                                 afterany(prod)
32575976 chan-tests2    (validates from_inq_field; concurrent)
```

`submit-channeling-twin.sh` now encodes the 4-stage structure for future use.

---

## 2026-08-01 (evening) — first submission FAILED at the test gate; two real defects found and fixed; resubmitted as 32575316–22

### What happened

Chain 1 (32573554–60) was submitted and **`chan-tests` FAILED after 10 min**,
which cancelled every downstream job. That is the gate working: **zero GPU hours
were spent on production.** Two independent defects, neither in the physics code.

### Defect 1 — MY TEST was wrong, the library was right

`test_projectile_force_minimum_image_engine` failed its FIRST case ("deep inside
the cell, clipped and minimum-image agree"):

```
CHECK( e_plain == Approx(e_mini).epsilon(1e-6) )
with expansion:  0.3292947976 == Approx( 0.3290072435 )     -> 5.3e-4 off
```

**The assertion that mattered PASSED**: `e_mini == exact_energy(3.0)` to 1e-6 —
the minimum-image kernel hits the closed-form Fourier coefficient exactly. Only
the "the two kernels agree away from the face" premise was false, because at
L = 16 the test position Z = 3 is **only 3.57 σ from the +z face**. The plain
Gaussian loses 1.8e-4 of its charge there, and the minimum-image version *wraps*
that tail to z ≈ −8 where cos(2πz/L) = −1 instead of dropping it — so the two
differ by 5.3e-4 (energy) / 1.5e-3 (force), a thousand times the tolerance.

**Fix:** box 16 → **24 Bohr**, putting Z = 3 at 6.43 σ. Parameters DERIVED by
quadrature before editing, not tuned until green:

| | L = 16 (failed) | L = 24 (fixed) |
|---|---|---|
| Z = 3 distance to face | 3.57 σ | **6.43 σ** |
| plain vs mini, energy | 5.3e-4 | **9.7e-11** |
| plain vs mini, force | 1.5e-3 | **1.8e-9** |
| straddle: plain force error | — | **328 %** (assertion needs > 5 %) |

The residual 2.9e-5 between the min-image force and the closed form is **not
slack**: it is the finite-difference sinc error (kδ)²/6 = 2.86e-5 for k = 2π/24,
δ = 0.05. Landing exactly there is itself a check.

`test_radial_occupancy_engine` **PASSED** (4.58 s) — the new observable compiles
under CUDA and satisfies its analytic Rayleigh assertions. `test_projectile`
(pure) passed 10144 assertions.

### Defect 2 — /rds IS 98 % FULL AND THE RUN AS CONFIGURED DID NOT FIT

`chan-nb` reported `disk headroom: 26 GB` and refused. **The guard was correct**:

```
/rds-d6/user/skcb2   1073.5 GB used of 1099.5 GB quota   (df: 25G avail, 98% full)
ResearchProject alone = 988 GB
```

Measured cost on THIS grid (80×80×120 = 768k points): a VTI frame is
**8.7 MB** (11.31 bytes/point, measured from an existing slab_ks_wrap frame) and
an RT checkpoint of 104 states is **1.28 GB**.

| configuration | fields WP | fields cl | checkpoints | TOTAL | vs 25 GB free |
|---|---|---|---|---|---|
| **as originally submitted** (SAVE 5 / WF 75 / MAX_CKPT 3) | 8.2 | 5.2 | 11.5 | **24.9 GB** | **would have filled the disk mid-run** |
| **lean — now the default** (SAVE 15 / WF 150 / MAX_CKPT 1) | 2.8 | 1.8 | 6.4 | **11.0 GB** | fits |
| minimal (SAVE 25 / no wavefns / MAX_CKPT 1) | 1.6 | 1.1 | 6.4 | 9.0 GB | fits |

**Nothing scientific is lost at "lean":** every scalar observable is still written
EVERY step, and 101 frames is well above the 30 the GIF builders actually sample.

Three fixes, all recorded in the scripts with their reasoning:
1. **Lean cadence is the default** in `run-chan-twin.slurm`.
2. **Pre-flight disk guards added to the RUN jobs** (`run-chan-twin.slurm`,
   `run-chan-gs.slurm`). Only the notebook job had one, so a full filesystem was
   previously discovered *after* the GPU hours were spent. This is the failure
   that killed 11 of 16 slab runs on 2026-07-31.
3. **Notebook guard 30 GB → 3 GB.** The 30 GB was inherited from a SIXTEEN-run
   campaign on a 1.65M-point grid; two runs here need well under 1 GB. It was
   simply the wrong number carried across.

### STANDING PROBLEM FOR THE USER (not actioned — it is their data)

**/rds is at 98 % of a 1 TB quota, with 988 GB in `ResearchProject`.** The lean
twin fits in the remaining 25 GB, but with ~14 GB to spare there is no room for
another campaign. Old run output (VTI frames and retained checkpoints from
completed studies) is the obvious candidate. **Nothing has been deleted.**

### Chain 2 — resubmitted 2026-08-01

```
32575316 chan-tests     32575317 chan-gs
32575318 chan-wp-smoke  32575319 chan-cl-smoke
32575320 chan-wp        32575321 chan-cl        32575322 chan-nb
```
Also fixed: `submit-channeling-twin.sh` now passes `--job-name` per half, so the
classical and wavepacket jobs are distinguishable in `squeue` (previously all four
appeared as `chan-twin`, and neither `squeue` nor `scontrol` shows a job's
arguments — the user could not tell the classical half had been submitted).

---

## 2026-08-01 — implementation complete; nothing has been RUN yet

### State in one line

Every artefact of the channeling twin is written, statically verified and
smoke-tested against synthetic data. **No GPU job has been submitted.** The whole
chain is one command: `./shared/bin/submit-channeling-twin.sh` from the repo root.

### What the study is

A matched **classical + wavepacket twin** shot on-axis down the hollow bore of a
periodic r_s = 3 annular jellium tube, to validate a **KS-orbital definition of
stopping power** against the classical ΔE/ds one. The claim has three parts and
all three are measured: the stopping powers agree (**result**), the packet stayed
in the bore (**premise**), and that froze var(p) (**mechanism**). See the plan §1.

### Done — code

| Artefact | Path | State |
|---|---|---|
| Locked geometry/physics config | `ResearchProject/systems/cylindrical_jellium/shared/configs/channeling_tube_rs3.hpp` | written; every derived number verified numerically |
| Ground state | `.../scripts/channeling_twin/gs/run.cpp` | written, **not compiled** |
| Classical twin | `.../scripts/channeling_twin/classical/run.cpp` | written, **not compiled** |
| Wavepacket twin | `.../scripts/channeling_twin/wp/run.cpp` | written, **not compiled** |
| NEW observable `radial_occupancy` | `inq-stack/include/inqkit/observables/radial_occupancy.hpp` | written, **not compiled** |
| minimum-image option on the HF force | `inq-stack/include/inqkit/dynamics/projectile_force.hpp` | edited (defaults preserve every existing run) |
| Engine test: radial occupancy | `inq-stack/tests/include/inqkit/observables/test_radial_occupancy_engine.cpp` | written, **not compiled**; registered in CMake |
| Engine test: min-image force | `inq-stack/tests/include/inqkit/dynamics/test_projectile_force_minimum_image_engine.cpp` | written, **not compiled**; registered in CMake |
| Dispatch chain | `shared/bin/run-chan-{tests,gs,twin,notebooks}.slurm`, `shared/bin/submit-channeling-twin.sh` | written, `bash -n` clean, executable |
| Analysis engine | `.../hypotheses/channeling_twin/channeling_stopping.py` | written, **13 tests PASS** |
| Per-run notebook driver | `.../hypotheses/channeling_twin/build_run_notebooks.py` | written |
| **Phase (comparison) notebook builder** | `.../hypotheses/channeling_twin/build_comparison_notebook.py` | written, **executes end-to-end on synthetic data, 0 errors, 9 inline figures** |

### Done — verification (what was actually checked, and how)

| Check | Result |
|---|---|
| `cutoff_guard.py` (mandatory pre-run gate) | **PASS** both halves. WP aliased tail **0.00 %** (σ_p = 0.177 vs k_Nyq = 6.283); classical E_cut = 537 eV ≥ 1.10 × 50 eV |
| Derived physics numbers | verified numerically: V_annulus = 18095.5737 Bohr³, n0 = 8.84194e-3, **r_s = 3.000000**, ħω_p = 9.0705 eV, v_F = 0.63972, **v/v_F = 2.997**, λ_p = 36.14 Bohr |
| Run sizing | 1500 × 0.02 = 30 a.u. → **57.51 Bohr**, z = −28 → **+29.51**: one traversal, no wrap |
| Dispersion budget | 2σ_d reaches R_in at **t = 23.32 a.u.**; 6σ_d = L_xy at **t = 34.15 a.u. (after the run ends)** — no transverse image overlap at any time |
| CMake reconfigure (login node, required — compute nodes have no network) | **done**, both new engine test targets registered |
| SLURM scripts | `bash -n` clean |
| `channeling_stopping.py` | **13/13 tests pass** (`hypotheses/channeling_twin/tests/`) |
| Comparison notebook | built + **executed** against a synthetic twin: 23 cells, **0 errors**, 9 inline figures, verdict machinery returned `AIM MET` and recovered S = 0.200 eV/Bohr on both halves to 1e-8 |

**NOT verified (cannot be, without a GPU node):** nothing C++ has been compiled.
Compile errors in the three `run.cpp` and the two new headers will surface at
stages 1 and 3 of the chain. That is where to look first if a job fails.

### The tests, and why each exists

| Test | Tier | Pins |
|---|---|---|
| `test_radial_occupancy_engine` | engine | f_bore against the **Rayleigh law** for an on-axis Gaussian (analytic, not captured); the exact 3-shell partition; the minimum-image case where a packet's tail wraps a transverse face (a non-periodic implementation returns 0.6 instead of 0.98); an off-axis packet where ⟨r⊥²⟩ = μ² + 2σ_d² is exact |
| `test_projectile_force_minimum_image_engine` | engine | the HF force against the **closed-form Fourier coefficient** of a cosine drag field, at the cell centre (both kernels agree) and straddling the face (only the min-image one is right — asserted, because that failure is the reason the flag exists) |
| `tests/test_channeling_stopping.py` | analysis | constant-S has a closed-form trajectory (`dp/dt = −S`, so p(t) is LINEAR and the trapezoid rebuild of s₄ is exact) → all four S_ij and the classical S must return the input to **1e-8**. Plus: the window follows the measured f_bore; freeze detection; **all three verdict branches**, including "clean channeling but S still differs"; resume-segment concatenation |
| `tests/test_comparison_notebook_cells.py` | static | the builder writes Python source as strings, so its failure modes are invisible to import. Guards non-raw backslash literals (`"$\approx$"` silently becomes a BEL character — this bit the build once), mathtext commands matplotlib lacks (`\le` vs `\leq`), unbalanced `$`, and that the three verdict figure blocks + the mandatory density GIF are still emitted |

### Decisions that must survive compaction

1. **No `s5` in-medium path correction, deliberately.** The tube is uniform along
   z, so the medium fills every z the projectile visits and the path IS the
   in-medium path. The slab study needed s5 only because 25 of its 85 Bohr were
   vacuum. Do not "port s5 across" — it would be wrong here.
2. **The fit window comes from the MEASURED `f_bore(t)`**, not from the
   free-dispersion formula. Window ends at the first breach of f_bore < 0.95.
3. **Minimum image everywhere, in both halves.** Launch is 0.71 σ_pot from the
   −z face; a clipped Gaussian loses 24 % of its charge ASYMMETRICALLY across the
   ±δ finite difference, so the error does NOT cancel out of the force. This is
   also why `projectile_force.hpp` gained a `minimum_image` flag.
4. **The circular centroid is mandatory**, for the same reason: ~24 % of the WP
   is across the periodic face at t = 0. ⟨p_z⟩/T1/T2 are unaffected (momentum-space
   expectation values), so the primary measurement is clean regardless — but s3 is
   not, and the naive ⟨z⟩ fails silently rather than loudly.
5. **`S_24` is the headline estimator** (drift energy vs drift path). It is built
   from ⟨p_z⟩ on both sides, so it is a stopping power whether or not var(p) is
   frozen; the other three are cross-checks.
6. **AIM MET requires all three of result/premise/mechanism.** The interesting
   failure mode — clean channeling, frozen var(p), and the stopping powers STILL
   differ — has its own verdict branch pointing at **E_PP, the WP self-Hartree**
   (SIE is a property of the orbital, not its environment, so channeling should
   NOT remove it). That branch is unit-tested.
7. **N = 160 is derived, not chosen.** It is the electron count for which
   n0 = N/V_annulus lands on r_s = 3.000000 for this exact geometry.
8. **Results layout** follows `slab_ks_wrap`
   (`scripts/channeling_twin/{half}/results/{name}/`), not the
   `<sweep>/<run>/` form of the file-placement rule — that is what the run-notebook
   builder, `check_twin.py` and `ks_stopping.py` all expect on this device.

### Not done / next actions

1. **Submit the chain** (the only remaining step):
   ```bash
   cd /rds/user/skcb2/hpc-work/tddft/inq-tddft-research
   ./shared/bin/submit-channeling-twin.sh
   ```
   Stages: `chan-tests` → `chan-gs` → 2 × `chan-twin smoke` → 2 × `chan-twin`
   (concurrent) → `chan-nb`.
   **Cost estimate (UNMEASURED — no timing run has been done):** GS up to ~8 h
   wall allotted; each production half ~2–4 GPU-h on a 768k-point grid with 104
   states for 1500 steps, both concurrent. Treat these as guesses until the smoke
   stage gives a per-step time.
2. **Kill/resume recipe** (everything checkpoints every 500 steps):
   ```bash
   scancel <jobid>
   CH_RESUME=1 sbatch shared/bin/run-chan-twin.slurm wp          # resume
   CH_RESUME=1 CH_N_STEPS=3000 sbatch shared/bin/run-chan-twin.slurm wp   # extend
   ```
3. **Catalogue the runs** (`tddft-run-catalogue` skill) once they complete.
4. **Journal entry** with the user's own observations once the notebook exists.

### Traps already paid for (do not re-discover)

- A **non-raw** Python string containing `\approx` in the notebook builder becomes
  a BEL character and kills the notebook 200 lines into execution with a mathtext
  error pointing at "pprox". Guarded now by
  `tests/test_comparison_notebook_cells.py`.
- `\le` is not a matplotlib mathtext command; `\leq` is. Same guard.
- The comparison notebook writes `comparison_figs/` and `stopping_summary.csv`
  into its own directory. A synthetic smoke build therefore drops fake artefacts
  into the repo — **they were deleted after the 2026-08-01 smoke test**; if you
  re-smoke, delete them again or they will masquerade as results.
- `check_twin.py` parses `run_summary.txt` with a `(\w+)\s*[=:]\s*(\S+)` regex and
  reads ONE token, so both halves write a single-token `projectile = ...` value
  (`gaussian_charge_perturbation` / `wavepacket_orbital`) for its
  "the twins actually differ" check. Do not turn that back into prose.
