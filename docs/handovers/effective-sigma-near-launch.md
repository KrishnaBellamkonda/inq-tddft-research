# Handover — effective-σ hypothesis: near-launch σ_WP = 0.5 wavepacket sweep

Plan: `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/effective-sigma-near-launch.md`
Parent campaign: `docs/campaigns/localised_jellium/classical-highdensity-sv-benchmark.md`
(id `classical-highdensity-sv`); parent handover
`docs/handovers/wavepacket-highdensity-sv-twin.md`.
Branch `quantum-stopping-power`. Machine **CSD3**, `ampere`,
account `mphil-nikiforakis-skcb2-sl2-gpu`. Started 2026-08-01.

---

## STATUS 2026-08-01 — verification PASSED at the floor; production sweep LAUNCHED

### Goal (one line)

Hold σ_WP = 0.5 and move the launch point from z = −24 to **z = −14** (1.5 Bohr
outside the slab face at −12.5), so the packet arrives essentially undispersed
instead of 4.7–8.1 Bohr wide. If S(v) changes, the *launch* σ was never the
parameter controlling the wavepacket–slab interaction — the **arrival width** is.

### The verification the user asked for — PASSED, at the 1.5 Bohr floor

Scan job **32528019**. Accepted on the FIRST trial, so no retreat was needed.

| | z = −24 (regression) | **z = −14 (accepted)** |
|---|---|---|
| standoff from face | 11.5 Bohr | **1.5 Bohr** |
| weight removed by orthogonalisation | 2.4e-5 % | **0.109 %** (criterion < 3 %) |
| max single overlap | 3.691564855e-4 | 0.01240 |
| closure residual (rel) | 2.6e-10 | 1.1e-13 |
| centroid z | −23.99999 | −13.9946 |
| density std z, FULL profile | 0.3537 | 0.4684 (+32.5 %) |
| **density std z, CORE** | 0.3536 | **0.3559 (+0.65 %)** |
| ⟨p_z⟩ | 1.99979 | 2.00191 |
| σ_pz² | 2.00073 | 1.99898 |
| σ_kz (k_z marginal) | 1.414473 | 1.4139 (vs 1.414214) |
| skew / excess kurtosis | −0.005 / +0.022 | −0.007 / +0.028 |
| **R² vs ANALYTIC N(k₀,σ_p²)** | 1.000000 | **0.99995** |

**The +32.5 % full-profile width is NOT a deformed packet** — resolve this before
re-deriving it. The packet CORE (|z−z_c| < 3 Bohr, holding 99.93 % of the weight)
is 0.3559 Bohr against the Gaussian's σ/√2 = 0.3536, i.e. **+0.65 %**. The
inflation comes from **0.071 % of the weight sitting 10–25 Bohr away, inside the
slab** (0.112 % inside |z|<12.5): variance weights by (z−z_c)², so a 0.1 % tail at
~15 Bohr adds ≈0.22 against a core variance of 0.125 and contributes **42 % of the
total variance**. A second moment is the wrong statistic for a distribution with a
small far tail.

That tail is **physics, not an artefact**: an electron 1.5 Bohr from a metal
surface must be Pauli-orthogonal to the occupied states, which extend through the
slab. The far-launch packet has none of it because at 11.5 Bohr there is nothing
to be orthogonal to. It does mean the near-launch packet starts with ~0.1 % of
itself already inside the slab — worth a caveat line in the notebooks, not a
blocker.

### A second, independent reason the near launch is better posed (found 2026-08-01)

Not the argument that motivated the campaign, and stronger than it. Two clocks
constrain when S can be measured at all:

- the packet's CENTROID must be inside the slab — `slab_window(v, z0)`;
- the packet must not yet overlap its own transverse periodic images —
  `transverse_overlap_time(0.5) = 4.12 a.u.`

At σ = 0.5 these two windows **barely intersect for the far launch**:

| v | FAR slab window | usable | NEAR slab window | usable |
|---|---|---|---|---|
| 2.0 | [5.75, 18.25] | **0.00** | [0.75, 13.25] | **3.37** |
| 2.5 | [4.60, 14.60] | **0.00** | [0.60, 10.60] | **3.52** |
| 3.0 | [3.83, 12.17] | 0.28 | [0.50, 8.83] | **3.62** |
| 3.5 | [3.29, 10.43] | 0.83 | [0.43, 7.57] | **3.62** |

(a.u.; "usable" = overlap of the slab window with the transversely-clean window
[0.50, 4.12].)

At v = 2.0 and 2.5 the far-launch σ = 0.5 campaign has **literally zero time**
during which the packet is both inside the slab and transversely clean. This is
NOT a new discovery — `wp_hd_stopping.slab_window`'s own docstring records it
("starts AFTER the localised window closes … which is why it returns a negative
S"). What is new is that the near launch **repairs** it: 3.4–3.6 a.u. of genuinely
usable data at every velocity.

**Scope this correctly.** It does NOT invalidate the published σ-sweep figure,
which uses the DEPOSIT estimator S = (E_total(t_f) − E_GS)/L_slab and never
touches the fit window. It does mean the localised slope estimators
(`fit_all`, T1/T2/s3/s4) were structurally unusable at far launch for σ = 0.5 and
become usable near launch — so the two campaigns should be compared on the
deposit estimator first, with the slope estimators as a near-launch-only bonus.

### Locked user decisions (2026-08-01)

| Decision | Value |
|---|---|
| Launch z | −14.0 (1.5 Bohr beyond the face); **floor — scan retreats only** |
| Scan step | 0.5 Bohr outward |
| Criterion | orthogonalisation-removed weight < 3 % |
| Scan velocity | v = 2.0 only (worst case), result applied to all four |
| Gaussianity | reported, **never vetoes** |
| σ | 0.5 only, v ∈ {2.0, 2.5, 3.0, 3.5} |

Decided by me, stated not asked: N_STEPS **identical** to far-launch
(3623/2898/2415/2070 — shorter path ⇒ more post-exit plateau, same time budget);
four **new vacuum CAP controls** at the same launch z (a baseline is only
subtractable step-for-step if the launch matches).

### The first chain ABORTED on a t=0 gate — read this before touching the gates

First submission (32528175–8): smoke **FAILED exit 4**, sweep 32528176
**CANCELLED** by its `afterok`. Cause:

    [FAIL] density std = s/sqrt2: 0.4684 (expect 0.3536, dev 32.47 %, tol +/-5 %)

That is the raw second moment discussed above — inflated by the PHYSICAL
orthogonality tail, with the packet core at +0.65 %. The gate was correct
behaviour from a safety net that had never met a near-slab launch; the packet was
not wrong. Every σ-sensitive gate passed (σ_pz² −0.05 %, T1−T2 +0.05 %).

**Fix (wp/run.cpp) — do not "simplify" this back to a wider tolerance.**
Widening would weaken a genuine blunder-catcher for every campaign. Instead:

1. NEW always-strict gate: `ortho removed weight < 3 %` (the user's criterion,
   `LJ_ORTHO_TOL_PC`). It previously had **no enforcement inside the production
   run at all**.
2. The raw real-space width gate stays strict only when
   `removed_weight < 1e-5` (far launch). Otherwise it becomes INFO plus a
   CONSISTENCY gate: the excess variance must be explainable by the measured
   loss sitting at a physically sensible distance,
   `d = sqrt((var_meas - sigma^2/2) / removed_weight)`, required to satisfy
   `3*sigma_d < d < Lz/2`.
3. Nothing is lost: σ is still probed strictly in both regimes by `sigma_pz^2`
   and `T1-T2 = 3/(4 sigma^2)`, which are momentum moments and therefore NOT
   distance-weighted — a factor-2 σ blunder is still caught.

**Re-run 32528286 — ALL t=0 GATES PASSED**, `d_implied = 9.323 Bohr` (launch −14,
slab face −12.5), i.e. the excess variance is exactly the distance to the slab
orbitals. Predicted 9.3 Bohr by hand beforehand; measured 9.323.

### RESULTS — all 8 runs complete (2026-08-01). NO VERDICT RENDERED YET.

All four near-launch runs `run_completed = true`, `steps_done == steps_target`,
`complete = True` (the `deposit_stopping` completeness flag that once silently
returned a plausible S from a 86/3623-step run — it was checked).

| v | S_raw FAR | S_raw NEAR | S_corr FAR | S_corr NEAR | norm_f FAR | norm_f NEAR |
|---|---|---|---|---|---|---|
| 2.0 | 2.4713 | 2.4626 | 0.2392 | **0.3194** | 0.0638 | **0.1061** |
| 2.5 | 2.4475 | 2.4748 | 0.2002 | **0.2385** | 0.0462 | **0.0693** |
| 3.0 | 2.4463 | 2.4359 | 0.1673 | 0.1717 | 0.0321 | 0.0398 |
| 3.5 | 2.4335 | 2.4150 | 0.1309 | 0.1267 | 0.0182 | 0.0210 |

(S in eV/Bohr, deposit definition S = (E_total(t_f) − E_GS)/L_slab, L_slab = 25.)

**Read the two columns differently — they are not raw-vs-CAP-corrected.**
`S_deposit_corrected` removes a known INQ ledger artefact: INQ reports the
kinetic term as occ·⟨ψ|T|ψ⟩/⟨ψ|ψ⟩ (`energy.hpp:50-55`), so under a CAP the
decaying WP orbital keeps contributing its per-particle MEAN and inflates
E_total. The correction `E − T1·(1 − norm)` is exact and IS the better estimate.
It is NOT a vacuum-control subtraction.

**What the numbers say.** RAW deposit is far ≈ near to within ~1 % at every
velocity, and essentially flat in v (~2.44). CORRECTED deposit separates at low
v: **+34 % at v = 2.0, +19 % at v = 2.5**, converging by v = 3.0–3.5.

**Why no verdict.** The corrected split tracks `norm_final`, and the near-launch
packet retains **66 % more norm** at v = 2.0. That must be explained before the
split is attributed to the slab interaction. The direction is, however, the
OPPOSITE of the naive artefact: the near-launch packet starts 10 Bohr CLOSER to
the +z CAP, so it reaches the absorber sooner and should be absorbed MORE, not
less. That it survives more is a genuine signal — plausibly the slab stopping or
reflecting a narrower arrival packet more strongly, which is what the hypothesis
predicts — but "plausibly" is not a result. Resolve via:
1. the per-run notebooks (job 32528288) and their density GIFs — is the packet
   being reflected, captured, or transmitted?
2. the LOCALISED slope estimators, now usable near-launch (3.4–3.6 a.u. of clean
   in-slab data vs 0.00 far-launch at v = 2.0/2.5) — an independent estimator
   that does not depend on norm_final at all. **This is the decisive check.**
3. the `interactions.csv` pairwise ledger (E_PS is the term that stops it).

### Production chain — RUNNING

| # | stage | job | state |
|---|---|---|---|
| 1 | smoke (BUILDS binary + t=0 gates) | **32528286** | **DONE — ALL GATES PASSED** |
| 2 | sweep array 0–3 (v = 2.0/2.5/3.0/3.5) | **32528287** | RUNNING (0,1 on GPU; 2,3 queued) |
| 3 | vac CAP baselines at z = −14 | **32528177** | RUNNING (`vac_nl_v2p0` … ) |
| 4 | notebooks + synthesis | **32528288** | COMPLETED exit 0 — but built **NO near-launch notebooks** (see below) |
| 5 | near-launch notebooks (after builder fix) | **32563217** | submitted |

**Trap the chained notebook stage fell into.** Job 32528288 exited 0 and looked
successful, but `build_run_notebooks.py` enumerated `W.SIGMAS` only, so it built
notebooks for the three FAR-launch campaigns and silently produced nothing for
`nl_*`. A clean exit code is not evidence the deliverable exists — check for the
files. Fixed by making a campaign a **(launch_z, sigma) pair**:

- `build_run(v, sigma, launch_z)` / `build_synthesis(sigma, launch_z)` /
  `_synth_name(sigma, launch_z)` / new `_campaign_tag(sigma, launch_z)`;
- `main()` enumerates
  `[(lz, s) for lz in (FAR, NEAR) for s in SIGMAS if has_campaign(s, lz)]`;
- new CLI target prefix `nl:` (`nl:0.5:3.0`, `nl:3.0`).

Verified enumeration (existing names unchanged):

    z=-24 sigma=0.5 -> v2p0 …            synthesis.ipynb
    z=-24 sigma=2   -> s2p0_v2p0 …       synthesis_s2p0.ipynb
    z=-24 sigma=3   -> s3p0_v2p0 …       synthesis_s3p0.ipynb
    z=-14 sigma=0.5 -> nl_v2p0 …         synthesis_nl.ipynb

**The subtle half of that fix:** the GENERATED notebooks embed `W.set_campaign(σ)`
in their own code cells. Without an accompanying `W.set_launch(...)` a re-executed
near-launch notebook would resolve the FAR-launch run of the same σ — the
near-launch tag defaults to `""` — and produce wrong numbers under a correct
title. Both emitted cells now set the launch explicitly.

Kill with `scancel 32528286 32528287 32528288 32528177`.
Run names are prefixed **`nl_`** (`nl_v2p0`, …, `vac_nl_v2p0`), so every existing
far-launch run, notebook and summary CSV keeps resolving unchanged.
Expect ~3 h (v = 2.0 is the long pole: 3623 steps at ~2.75 s/step).
`submit-wp-hd-nearlaunch.sh` submits the whole chain from scratch; the stages
above were resubmitted individually after the gate fix.

---

## What changed in code

### inqkit — the one library change (additive, back-compatible)

`inq-stack/include/inqkit/wavepacket/injection_report.hpp` +
`.../wavepacket/wavepacket.hpp`: new `norm_pre_ortho`, `norm_pre_renorm`,
`removed_weight`, `sum_overlap_sq`, `overlap_by_state`,
`ortho_closure_residual()`.

**Why it was needed:** `norm_after` is measured *after* the post-Gram-Schmidt
renormalisation (`wavepacket.hpp:394-405`) so it is ≈1 by construction and cannot
express the loss; `max_overlap` is only the single largest overlap. The total
Σᵢ|⟨ψᵢ|ψ_WP⟩|² was computed and thrown away.

**Subtlety worth not re-deriving:** `removed_weight` is a RATIO against the
measured raw-Gaussian norm, not against a hard 1.0. The injector's
`norm_fac = (πσ²)^{−3/4}` normalises the *continuum* Gaussian; at dx = 0.40 with
density std 0.354 Bohr the discrete norm departs from 1 by ~1 %, and comparing
against 1.0 would report that discretisation error as orthogonalisation loss —
enough to distort a 3 % budget.

### Test — `test_wp_ortho_loss_engine.cpp` (18 assertions, 4 cases, ALL PASS, job 32528040)

Registered in `inq-stack/tests/include/engine/CMakeLists.txt`. Oracle is the
closed form for a **constant** occupied state φ₀ = 1/√V:

    removed_weight = 8 π^{3/2} σ³ exp(−σ² k₀²) / V

| case | measured | analytic | rel |
|---|---|---|---|
| k₀ = 0 | 0.0445464505 | 0.044546624 | 4e-6 |
| k₀ = 1 | 8.159026e-4 | 8.158999e-4 | 3e-6 |
| closure residual | **0.0** bit-exact | 0 | — |

State 0 is overwritten by hand, so the test depends on no SCF result.
**nvcc gotcha:** the overwrite must live in a FREE function — an extended
`__device__` lambda inside a constructor fails to compile ("enclosing parent
function must allow its address to be taken"). Cost me one build cycle.

### New run machinery

- `ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/inject_scan/run.cpp`
  — GS load + injection, **propagates nothing**. Exit code IS the decision
  (0 accept / 3 reject / 5 closure failure), so the driver needs no parsing.
- `shared/bin/run-wp-hd-scan.slurm` — the outward 0.5 Bohr scan; its build step
  doubles as the z = −24 regression trial.
- `shared/bin/submit-wp-hd-nearlaunch.sh` — the 4-stage chain above.
- `run-wp-hd-wp.slurm` / `run-wp-hd-vac.slurm` — accept `LJ_LAUNCH_Z`
  (and `LJ_LAUNCH_Z=auto`, which reads the scan's `accepted_launch_z.txt` so no
  number is relayed by hand between jobs); `nl_` prefix. Default −24.0 leaves
  every existing name unchanged.
- `wp/run.cpp` — now records `removed_weight` / `norm_pre_ortho` /
  `ortho_closure_residual` in `wp_config.txt` for every production run.

### inqview

`visualisation/field_io.py`: `load_vti(array=...)`, new `load_complex_vti`,
`kz_marginal`, `gaussian_fit_quality`.

**Why the k_z marginal and not inqkit's radial `MomentumDistribution`:** for a
Gaussian WP the k_z MARGINAL is exactly N(k₀, σ_p²); the RADIAL n(|k|) of a
*drifting* packet is not Gaussian, so it cannot answer "is it still Gaussian".

**Ordering trap:** VTIs are in physical order, `np.fft.fftn` wants FFT-natural
order, so `kz_marginal` `ifftshift`s BEFORE transforming. This does not violate
the never-fftshift-a-VTI rule (that governs display of real-space data).

Analysis: `hypotheses/wp_highdensity_sv/scan_gaussianity.py` →
`scan_gaussianity.csv` + `scan_gaussianity.png`.

---

## Verified vs NOT verified

**VERIFIED**
- Slab geometry / 1.0 Bohr erfc softening (`slab_n100_L35x35x85.hpp:47,49`).
- Orthonormalisation is injection-only: ETRS (the default, what these runs use)
  has no `orthogonalize` call; only `crank_nicolson.hpp:139,162` does. So the
  whole risk is a t=0 property — no propagation needed to measure it.
- `removed_weight` against two independent analytic values + a bit-exact closure.
- The scan program reproduces the completed campaign's `max_overlap` at z = −24
  to 12 significant figures.
- The k_z extraction reproduces the known undeformed far-launch packet
  (R² = 1.000000 vs analytic).
- Real-space core width at z = −14 is +0.65 % of the Gaussian.

**NOT VERIFIED**
- Any near-launch **physics**. The sweep was submitted, not completed.
- Whether S(v) actually differs from far-launch — that IS the hypothesis.
- The ~0.1 % in-slab t=0 tail's effect on the deposit estimator S =
  (E_total(t_f) − E_GS)/L_slab. Expected small; not quantified.

---

## Known issue found in passing — NOT fixed, needs a user decision

`inqview.load_vti` builds its coordinate axes as `origin + (i + 0.5)·spacing`
(`field_io.py:109-111`, "cell-centred"), but inqkit writes samples at **node**
positions `origin + i·spacing`. Measured: the far-launch packet, unambiguously at
z = −24.0 by the C++ observable, reads back at **−23.803** from its VTI — exactly
+dz/2 = +0.197 Bohr.

Consequence: every real-space coordinate read off a VTI in this project carries a
**+½-cell offset**. It does NOT affect this campaign (the k_z marginal uses only
`spacing`; `removed_weight` and all C++ observables use INQ's node convention),
and it cancels in any *difference* of two VTI coordinates.

Left alone deliberately: correcting it would shift coordinates in existing
figures across the repo. **Surface to the user before changing.**

---

## Exact next steps

1. Read the smoke stage (32528175) t=0 gates and its `wp_config.txt`
   `removed_percent` — should reproduce 0.109 %.
2. When the sweep finishes, build the decisive figure: **S(v) far-launch vs
   near-launch at identical σ_WP = 0.5**, with σ = 2/3 far-launch traces for
   context. Extend `hypotheses/wp_highdensity_sv/wp_hd_stopping.py` with a
   near-launch campaign selector (it already has `set_campaign(σ)` / `sigma_tag`
   / `vac_name_for` — the `nl_` prefix needs the same treatment).
3. `deposit_stopping()` returns `complete` — **check it**. It once returned a
   perfectly plausible S from a still-propagating run (86 of 3623 steps).
4. Notebooks need the mandatory density-matrix GIF
   (`.claude/rules/notebook-density-gif.md`).
5. Caveat line in the notebooks: the near-launch packet starts with ~0.1 % of its
   weight already inside the slab (Pauli-orthogonality tail).

## Open / deliberately excluded

- No classical twin at the near-launch point. The classical projectile does not
  disperse, so its σ_pot = 0.354 already IS its arrival width — the existing
  classical curve is the correct narrow-arrival reference. Stated, not re-run.
- σ = 2 and 3 near-launch twins not run (user: σ = 0.5 only). Would upgrade the
  falsification test to a full 2-D (σ × launch distance) collapse test.
