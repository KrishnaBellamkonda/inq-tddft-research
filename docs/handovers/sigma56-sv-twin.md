# Handover — σ_WP = 5 and 6 Bohr classical+wavepacket S(v) twins

**Rolling file. Latest milestone at top.**
Plan: `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/sigma56-sv-twin.md`
Sweep name: **`sigma56_sv`**. Branch `quantum-stopping-power`.
Machine **CSD3**, `ampere`, account `mphil-nikiforakis-skcb2-sl2-gpu`. Started 2026-08-02.

Parents: `docs/handovers/wavepacket-highdensity-sv-twin.md` (σ = 0.5/2/3 WP sweep),
`docs/handovers/classical-highdensity-sv-benchmark.md` (σ = 0.5 classical benchmark),
`docs/handovers/slab-ks-orbital-stopping-wrap.md` (the E_absorbed machinery).

---

## STATUS 2026-08-05 (late) — 16/16 COMPLETE. The v = 2.5 peak is SYSTEMATIC.

Job 32880125 (`s5p0_v2p0`, fresh, gpu-q-41) **COMPLETED**: 4360/4360 steps,
4 h 11 m, `run_completed = true`. The fresh-restart decision was right — no hang,
one unsegmented dataset. **All 16 σ = 5/6 production runs are now complete.**

### The result it was launched to settle

| σ_WP | v=2.0 | 2.5 | 3.0 | 3.5 |
|---|---|---|---|---|
| 5 | **0.390** (new) | 0.459 | 0.396 | 0.305 |
| 6 | 0.374 | 0.442 | 0.380 | 0.293 |

(WP, eV/Bohr.) **σ = 5 dips at v = 2.0 exactly as σ = 6 does.** Both curves peak
at v = 2.5 with near-identical shape, so the non-monotonicity is a SYSTEMATIC
feature of the wavepacket half at these widths — not a bad run. The earlier
hypothesis that σ = 6 v = 2.0 was an outlier is FALSIFIED; do not re-open it.

Two supports: the new run's `norm_final = 3.1e-6` is the LARGEST in the σ = 5 set,
so it is the least norm-corrected point of that series and still lands low; and
its ⟨σ_r⟩ = 7.76 is unchanged from the value computed at 45 % completion,
confirming the width is fully determined long before the run ends.

**The v = 2.5 peak is now the live open question for this campaign.** The
norm-correction-dominance concern below still applies at v = 3.5 (94 % of the raw
value removed) but cannot explain a peak at 2.5 flanked by two well-corrected
points.

### Figures rebuilt (all four, standalone + report + panel slots)

`build_sv_effective_width_s6.py` coverage audit now prints **16/16**. Rebuilt:
`S_of_v_effective_width_s56.png`, `norm_loss_vs_sigma.png`, and their
`slab_*` report copies + `slab_panel/` slot copies.

### Also produced this session (all mirrored to the report draft)

Per the standing rule (`.claude` memory `feedback-slab-figures-to-report2`), every
slab figure is written to
`/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/reports/report2/drafts/draft1/figures/jellium_slab/`
by its builder, plus a `slab_panel/` copy at the LaTeX slot width:

| builder (in `hypotheses/sigma56_sv/`) | figure |
|---|---|
| `build_sv_effective_width_s6.py` | S(v), WP relabelled by ⟨σ_r⟩ |
| `make_norm_loss.py` | bath norm lost to CAP vs σ (+ `norm_loss_table.csv`) |
| `make_s6_v3_diagnostics.py` | pairwise interaction ledger; Δk line plot |
| `make_s6_v3_momentum_map.py` | Δ|ψ(k)|² heatmap in (k_z, k_⊥) |
| `make_s6_v3_energy.py` | ΔE_dep(t), classical vs WP, plateau evidence |

### ⚠ ESTIMATOR: only the CORRECTED pair plateaus — this decides the open question

Measured on the σ = 6, v = 3.0 case study (`make_s6_v3_energy.py`):

| curve | final | drift over last 10 % of run |
|---|---|---|
| classical, E_PS-subtracted | 2.77 eV | **+0.005 eV** ✅ plateaued |
| classical, raw | 11.32 eV | −1.04 eV ❌ still falling |
| WP, norm-corrected | 9.49 eV | **+0.009 eV** ✅ plateaued |
| WP, raw | 14.84 eV | −38.8 eV ❌ |

The classical raw curve never settles because the projectile's E_PS monopole tail
stays in the ledger out to z = 321 Bohr. So `S_deposit_eV_per_Bohr` is the only
estimator on which BOTH halves plateau — confirming draft1/CLAUDE.md landmine 1
by direct measurement rather than assertion.

**CONSEQUENCE, NOT YET APPLIED (needs the user):** `S_of_v_effective_width_s56.png`
still plots `S_eV_per_Bohr` for both halves. Switching it to `S_deposit` lowers
every classical curve substantially and FLIPS the headline — corrected, the WP
deposits MORE than its classical twin (9.49 vs 2.77 eV at σ = 6, v = 3.0), which
is the direction the 2026-08-03 entry below reports. Panels (a) and (c) of the
report figure therefore currently use DIFFERENT estimators. Do not ship the panel
until this is resolved.

### Mechanism finding — why classical and WP converge at σ = 6

Measured across all five widths (WP half; legacy σ = 0.5/2/3 have `interactions.csv`
but no classical twins):

| σ_WP | E_PP(0) | T_int(0) = var(p)/2m | ΔE_PP over transit | Δσ_k⊥ at exit |
|---|---|---|---|---|
| 0.5 | 20.2 eV | 81.6 eV | −21 eV | 0.068 |
| 2 | 3.9 | 5.10 | −4.4 | 0.067–0.085 |
| 3 | 2.1 | 2.27 | −2.2 | 0.061–0.068 |
| 5 | 0.65 | 0.82 | −0.66 | 0.022–0.030 |
| 6 | 0.29 | 0.57 | −0.36 | 0.013–0.018 |

`T_int(0)` reproduces the analytic `3/(4σ²)` to the printed precision at every
width — the channel is correctly identified.

1. **The two quantum-only channels shrink with width.** At σ = 0.5 the spurious
   self-Hartree release (21 eV) is ~3× the deposit; at σ = 6 it is 4 % of it.
   That ratio, not any separate mechanism, is why the twins converge.
2. **Transverse transfer has a THRESHOLD, not a trend.** Δσ_k⊥ ≈ 0.07 for every
   width up to σ = 3, then collapses 4–5×. Explanation (inference, consistent
   with the data): transverse coupling is capped by whichever is smaller — the
   medium's plasmon cutoff `q_c ≈ ω_p/v_F = 0.44` a.u. (ω_p = 0.203, v_F = 0.459
   at r_s = 4.18) or the projectile form factor's `1/σ_pot`. They cross at
   σ_WP ≈ 3.2 Bohr, exactly where the plateau ends. Longitudinal coupling needs
   only `q_z = ω_p/v = 0.068` a.u. and is untouched → quasi-1D collision.
   **TESTABLE:** q_c ∝ n^(1/6), so a denser bath moves the threshold to smaller
   σ — the r_s = 2.5 slab of `nazarov-gross-slab-mass-ladder.md` is the test.
3. **UNVERIFIED:** measured E_PP(0) matches the free-space Gaussian self-energy
   `1/(2 σ_pot √π)` at σ = 0.5 (20.2 vs 21.7 eV) but not at σ = 6 (0.29 vs 1.81).
   Suspected cause: at σ_pot = 4.24 the packet is a large fraction of L_xy = 35,
   so periodic images and the G = 0 removal suppress it. Not checked.
4. **Do NOT use σ = 0.5 Δ⟨k_z⟩ for momentum balance.** It loses 5× more momentum
   than any other width yet deposits the least; by the slab exit that packet has
   ⟨σ_r⟩ ≈ 10 Bohr and is already losing norm, so ⟨k_z⟩ is a moment of a
   truncated distribution. σ = 5/6 stay intact through transit and are safe.

---

## STATUS 2026-08-05 — EFFECTIVE-WIDTH FIGURE REMADE (legend σ_WP convention; ⟨σ_r⟩ exclusion)

**Scope: figure only. No run was launched, re-run, or re-analysed; no S value changed.**

Rebuilt `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv/S_of_v_effective_width_s56.png`
via its builder `.../sigma56_sv/build_sv_effective_width_s6.py`. Three changes,
all user-directed (2026-08-05):

1. **Classical legend entries moved to σ_WP** (`classical σ_WP = 5`), replacing
   the upstream `cl. σ_pot = 3.54`. Enforces `.claude/rules/sigma-wp-convention.md`
   — a classical twin is reported at the σ_WP it is matched to, never at the
   derived σ_pot = σ_WP/√2. WP entries now read `WP σ_WP = 5 (⟨σ_r⟩ = 7.8)`, so a
   twin pair shares one label AND keeps the effective width this figure exists to
   show. σ_pot survives only as a methods footnote in the builder docstring
   (σ = 5 → 3.54, σ = 6 → 4.24 Bohr).
2. **⟨σ_r⟩ for σ_WP = 6 now EXCLUDES the v = 2.0 run** (user: its value is not to
   be used). Implemented as the table `EFF_EXCLUDE_V = {6.0: {2.0}}`, not a
   hard-coded number, because the label is a MAX over velocities and one run
   therefore defines it silently. Effect: **σ = 6 label 8.45 → 8.15**, the max
   passing to v = 2.5. All other labels unchanged.
3. **Legend moved below the axes** (2 columns). The σ_WP labels are ~2× wider than
   the σ_pot ones and the upper-right box covered the σ = 0.5 classical point at
   v = 2.0 (S_B = 1.30), the highest point on the figure.

### Run coverage — 15/16 σ = 5/6 runs are on the figure

Verified against every run directory, not against the cached CSV. The builder now
PRINTS this audit each build (`USED`/`SKIP` + reason) instead of silently
dropping runs:

| σ_WP | half | v = 2.0 | 2.5 | 3.0 | 3.5 |
|---|---|---|---|---|---|
| 5 | classical | ✅ | ✅ | ✅ | ✅ |
| 5 | WP | ❌ **incomplete** | ✅ | ✅ | ✅ |
| 6 | classical | ✅ | ✅ | ✅ | ✅ |
| 6 | WP | ✅ | ✅ | ✅ | ✅ |

**The single gap:** `scripts/sigma56_sv/wp/results/s5p0_v2p0` stopped at
**1973/4360 steps (45%)**. Its `E_absorbed` is a mid-transit excitation, NOT a
stopping value (S_raw = 1.19 against ~0.39 for its finished siblings), so it is
excluded by the `complete` gate — user-confirmed 2026-08-05.
**σ = 6 v = 2.0 is NOT affected** — it is complete (4360/4360, `run_completed = true`,
4.1 h wall, single unsegmented dataset) and is already on the figure. Only the
σ = 5 low-velocity run is short; the two are not symmetric.

### 🚀 RELAUNCHED 2026-08-05 — job 32880125, FRESH (not resumed)

`sbatch --exclude=gpu-q-2,gpu-q-25 --time=08:00:00 --export=ALL,LJ_SIGMA=5.0,LJ_RESUME=0 shared/bin/run-s56-wp.slurm 0`

**Why fresh and not resumed**, despite a valid checkpoint at step 1744 (~2600
steps, ~2.3 h, vs 4360 steps / ~4.1 h fresh): this run has failed 2 of 3 resume
attempts, both by HANGING on the checkpoint read, and the second hang was not on
a node the bad-node diagnosis covers.

| job | node | outcome |
|---|---|---|
| 32691775 | gpu-q-4 | crashed at step 1750 — orphan VTI frame `density_t001750.vti`, `overwrite=false`. Cause since FIXED (`prune_orphan_frames.sh`). |
| 32692519 | **gpu-q-2** | hung 1 h 38 m at `loading electrons from …/checkpoint/`; cancelled. Covered by the bad-node finding below. |
| 32707545 | gpu-q-9 | hung **24 h** at the same line, killed by the time limit. **gpu-q-9 is listed below as a node that "progressed normally" — so the bad-node explanation does not fully cover this one.** |

A fresh start also avoids the resume path structurally: `run.cpp:379` sets
`overwrite = (START == 0)`, so a fresh run may overwrite the existing frames while
a resume may not, and it yields ONE unsegmented dataset — the same shape as the
σ = 6 v = 2.0 run it will be compared against, with no segment concatenation as a
confound. Wall limit set to 8 h (2× the σ = 6 reference) so a hang costs 8 h, not 24.

**Pre-launch hygiene:** the four stale `*.from1744.csv` files were renamed
`*.aborted` per the stale-segment rule below, so `_concat` cannot mix the aborted
resume into the new dataset. The base CSVs are truncated by the fresh run
(`ObservablesWriter` / `ix.open()` default to truncate) and the VTI frames are
overwritten.

**When it lands:** re-run `build_sv_effective_width_s6.py` (the coverage audit
should print 16/16) and `finalize.py`; the σ = 5 WP curve then spans v = 2.0–3.5
like σ = 6, and the σ = 5 ⟨σ_r⟩ label question below can be settled on equal
footing.

### ⚠ OPEN CONCERN — the norm correction dominates the WP points at high v

Raised while auditing the above; recorded because it bears on how much any WP
point is worth, not just the missing one. From `s56_S_summary.csv`:

| σ_WP | v | norm_final | S_raw | S (corrected) | correction removes |
|---|---|---|---|---|---|
| 5 | 2.5 | 1.5e-08 | 0.82 | 0.46 | 44 % |
| 5 | 3.0 | 5.3e-10 | 0.77 | 0.40 | 48 % |
| 5 | 3.5 | 3.7e-10 | **6.10** | 0.31 | **95 %** |
| 6 | 2.0 | 8.1e-07 | 0.62 | 0.37 | 40 % |
| 6 | 3.5 | 3.9e-10 | **4.87** | 0.29 | **94 %** |

By t_final the packet norm is 1e-7…1e-10 — the CAP has absorbed essentially all
of it. The norm correction exists because INQ reports the orbital kinetic term as
`occ*<psi|T|psi>/<psi|psi>` (`inq/src/hamiltonian/energy.hpp:50-55`), so it divides
by that vanishing norm. At v = 3.5 the reported S is what survives a ~95 %
cancellation between two large numbers. **Inference (not yet tested):** the v = 2.5
"bump" in both WP curves may be where this correction stops being small rather than
physical structure. The run just launched (v = 2.0) will have the LARGEST final
norm of the set, so it is the least corrected — and therefore the most useful point
for testing whether the curve shape is real.

### Per-velocity ⟨σ_r⟩ (1 %-norm-loss window), for auditing the labels

| σ_WP | v=2.0 | v=2.5 | v=3.0 | v=3.5 | label used |
|---|---|---|---|---|---|
| 0.5 | 7.96 | 10.20 | **10.44** | 9.90 | 10.4 (v=3.0) |
| 2 | **7.21** | 6.40 | 5.77 | 5.28 | 7.2 (v=2.0) |
| 3 | **6.52** | 5.90 | 5.45 | 5.13 | 6.5 (v=2.0) |
| 5 | **7.76** | 7.34 | 7.05 | 6.85 | 7.8 (v=2.0) |
| 6 | *(excluded 8.45)* | **8.15** | 7.95 | 7.81 | 8.2 (v=2.5) |

### OPEN — two label questions this exposed, neither resolved

- **σ = 5 parity.** Its label (7.8) still comes from its **v = 2.0** run — the very
  run whose S point is excluded as incomplete. So the σ = 5 and σ = 6 labels are
  now built from different velocity sets ({2.0…3.5} vs {2.5…3.5}). If the intent
  behind excluding σ = 6 v = 2.0 was parity with the plotted points, σ = 5 should
  also drop v = 2.0 → **7.34**. Awaiting the user; one line in `EFF_EXCLUDE_V`.
- **The builder docstring's stated rationale is falsified by the σ = 0.5 row.** It
  claims the MAX "is invariably the SLOWEST run — it spends longest in flight". At
  σ = 0.5 the max is at v = 3.0, not v = 2.0, because the packet disperses so fast
  that the 1 % window closes in 6–9 a.u. and the slowest run hits the threshold
  EARLIEST (6.40 a.u.). The MAX rule still runs, but its justification does not
  hold across all five widths. Not corrected in the docstring pending a decision.

### Validation status

- Builder runs clean; figure regenerated and visually checked (no legend/data
  collision, 10 entries, twin pairs share colour + σ_WP label).
- Plotted σ = 5/6 values cross-checked against `s56_S_summary.csv` (e.g. σ = 5
  classical v = 2.0 → 0.5984; σ = 6 WP v = 3.5 → 0.2926) — agree.
- **Not** a code-test-bearing change: no estimator, kernel, or numeric was
  touched — labelling, run-coverage reporting, and legend placement only. The
  underlying `s56_stopping.measure` / `e_absorbed.measure_dir` path is unchanged
  and retains its existing validation.

---

## STATUS 2026-08-03 11:10 — σ = 6 TWIN SET COMPLETE + FIRST RESULT; σ = 5 BLOCKED ON DISK

### The result — matched width does NOT close the classical/quantum gap

All 8 σ = 6 production runs COMPLETE (4 WP + 4 classical twins). All 8 vacuum
CAP baselines COMPLETE. 17/26 runs, 9/16 production.

| v (a.u.) | classical S | WP S | WP/cl |
|---|---|---|---|
| 2.0 | 0.195 | 0.374 | 1.9 |
| 2.5 | 0.147 | 0.442 | 3.0 |
| 3.0 | 0.111 | 0.380 | 3.4 |
| 3.5 | 0.085 | 0.293 | 3.4 |

(eV/Bohr, `S = [E_total(t_f) − E_GS − E_PS(t_f)]/25`, WP norm-corrected.)

**Finding:** at σ_WP = 6 the packet spreads only ×1.12 in transit (σ_d 4.33 →
4.85 Bohr), so dispersion cannot explain the difference — yet the WP deposits
**1.9–3.4× MORE** than its width-matched classical twin. Leading suspect: the
projectile self-Hartree E_PP (see `.claude/rules/decomposed-interaction-energies.md`
and `docs/handovers/bulk-jellium-ks-stopping.md`, where a ratio ~2.2 is open —
note this campaign's 1.9–3.4 BRACKETS it, and in the same direction).

### ✅ E_PS ARTEFACT INDEPENDENTLY CONFIRMED by the merged upstream σ-sweep (2026-08-03)

Merge `92d1a5f` brought in `1afe2da` + `16382f0` from the other device, which
carry a classical CAP σ-sweep out to **σ_WP = 20** scored by BOTH estimators —
`S_A_keloss` (projectile KE loss) and `S_B_Eabs` (E_absorbed/L). They disagree,
and the disagreement is exactly the monopole tail.

Subtracting `100/z_final` (the bare monopole, in Ha → eV) from `E_absorbed_eV`:

| σ_WP | v | z_final | E_abs (eV) | monopole (eV) | S corrected | S_A_keloss |
|---|---|---|---|---|---|---|
| 0.5 | 2.0 | 201.9 | 32.48 | 13.48 | **0.760** | **0.763** |
| 0.5 | 2.5 | 240.8 | 24.45 | 11.30 | 0.526 | 0.515 |
| 3.0 | 2.0 | 242.7 | 20.54 | 11.21 | 0.373 | 0.270 |
| 17  | 3.5 | 265.6 | 10.52 | 10.24 | 0.011 | 0.004 |
| 20  | 3.5 | 265.7 | **10.41** | **10.24** | **0.007** | **0.002** |

**The σ = 20 row is the proof.** A projectile that diffuse barely interacts, so
KE loss → 0.002 — yet raw `S_B_Eabs` still reports **0.42 eV/Bohr**, essentially
σ-INDEPENDENT. That floor is `100/z_final` to **1.6 %**: at σ = 20 the raw
estimator is **98.5 % artefact**. At σ = 0.5, v = 2.0 the corrected value agrees
with KE loss to **0.3 %**.

**Consequence: `S_B_Eabs` in `S_of_v_cap.csv`, `S_of_v_cap_sigma.csv`,
`S_of_v_cap_sigma_wide.csv` all carry this**, as does
`wp_highdensity_sv/S_of_v_effective_width.png` (its y-axis is literally S_B).
Use `S_A_keloss`, or subtract `100/z_final` from `E_absorbed_eV`. The WP half is
unaffected — the CAP annihilates the packet, so its E_PS(t_f) is already ~0.

This is the same artefact found independently in `sigma56_sv` earlier today; the
σ = 20 sweep just makes it impossible to miss.

### ⚠ THE E_PS TAIL ARTEFACT — this reversed the result, read before trusting any older classical S

The first pass of this analysis reported classical 0.43–0.55, i.e. WP stopping
LESS. That was wrong, by a factor ~4 on the classical half only. Cause:

`S = [E_total(t_f) − E_GS]/L` assumes the projectile–bath interaction has decayed
to zero by t_f.
- **WP half — true.** The CAP annihilates the packet (norm_wp → 4e-10), so
  E_PS(t_f) = 1e-5 eV. The WP numbers NEVER changed. This is also why the WP half
  was `settled` and the classical half was not.
- **Classical half — false, permanently.** The projectile is a real moving charge
  that keeps going; its monopole tail falls off only as N_e/z. At t_f it sits at
  **z = 321 Bohr and still contributes E_PS = 8.5 eV** — 62–80 % of the raw
  classical "deposit" of 10.6–13.7 eV.

**Verified two independent ways:**
1. E_PS(t_f) matches the bare monopole `100/z` to **0.6–4.4 %** at all four v.
2. Against the projectile's OWN kinetic-energy loss (free Ehrenfest back-reaction,
   `projectile.csv` — these runs DO have back-reaction, v 3.000 → 2.967):

   | v | S from KE loss | S from field − E_PS | ratio | S uncorrected |
   |---|---|---|---|---|
   | 2.0 | 0.191 | 0.195 | 0.98 | 0.548 |
   | 2.5 | 0.143 | 0.147 | 0.97 | 0.492 |
   | 3.0 | 0.107 | 0.111 | 0.96 | 0.453 |
   | 3.5 | 0.081 | 0.085 | 0.95 | 0.426 |

   Two estimators sharing no machinery agree to 2–5 %. The uncorrected column is
   4.2–4.9× too high.

**The v-independent −1.04 eV "drift" was its fingerprint.** `N_STEPS` was sized so
`v·t_f = 4.36·(|z0| + L_z/2) = 349 Bohr` is CONSTANT across the sweep — so every
run ends at the same z and carries the same tail. A drift identical at every
velocity was never physics; it was geometry.

Implemented as `Point.S_deposit_eV_per_Bohr` (+ `e_ps_final_eV`, `z_proj_final`)
in `s56_stopping.py`; `build_sv_figure.py` plots it.

**ACTION FOR OLDER CAMPAIGNS — NOT YET DONE.** Any classical run scored with
`E_absorbed/L` and a projectile still inside the box carries this. Check
`slab_ks_wrap`, `wp_highdensity_sv`, `classical_highdensity_sv`, and the bulk
`bulk_ks_stopping_*` set. Cheap test: `e_ps` at the last row of
`interactions.csv`, or `N_e/z_proj(t_f)` where that file predates 2026-08-01.
NOTE the reference figure `S_of_v_v2_timeavg_sigmar.png` (`dyn_direct`) used a
DIFFERENT estimator — projectile KE loss (`S_of_v_direct.csv` has `v_final`,
`v_mean_slab`, `deposit_eV`) — so it is NOT contaminated, but it is also not the
same measurement as the WP curves it was plotted beside.

Artefacts (all in `hypotheses/sigma56_sv/`): `s56_S_summary.csv` (15 rows),
`s56_cap_cost.csv`, `S_of_v_sigma56.png`, `S_of_sigma_eq.png`.

### Validation status — QUALITATIVE result solid, QUANTITATIVE provisional

Three open issues, recorded so they are not silently dropped:

1. **Classical "drift" — DIAGNOSED AND FIXED.** Was the main open issue; it is
   the E_PS monopole tail documented above, not a leak and not ringing. The
   `settled = False` flag on the classical half is a TRUE POSITIVE against the
   raw column and is expected to persist — the tail decays as 1/z forever, so no
   affordable run length makes it settle. Judge the classical half on
   `S_deposit_eV_per_Bohr` (and its agreement with KE loss) instead.
2. **CAP-cost — FIXED 2026-08-03, and the answer reversed.** The original
   −6.6 % differenced complete `cl_s6p0_v3p0` (2907 steps) against INCOMPLETE
   `cl_nocap_s6p0_v3p0` (2680), folding 227 missing steps into what was reported
   as the absorber's cost. `s56_stopping.cap_cost()` now matches both traces at
   their largest COMMON step (2680), averaged over a 20-sample window:

   | | S (eV/Bohr) |
   |---|---|
   | CAP-on @ 2680 | 0.4859 |
   | CAP-free @ 2680 | 0.4862 |
   | **Δ** | **−0.00027 (−0.055 %)** |

   **The CAP is essentially free on the classical half** — it validates the
   "CAP on both halves" decision. The endpoint values (0.4529 / 0.4848) are kept
   in the CSV as `S_cap_*_endpoint` so the two conventions stay visible.
3. **Norm correction dominates at high v.** v = 3.5: norm_final = 3.9e-10,
   S_raw 4.87 → 0.29 (94 % cancellation of two large numbers). The v = 2.5 WP
   bump (0.44 vs 0.37/0.38 either side) is UNRESOLVED — real structure or
   correction noise. Bound it with the vacuum controls before quoting it.

### BLOCKER — disk quota exhausted

`1099.4 GB used / 1099.5 GB quota`, flagged `!*`. All failures are
`VTIImageDataWriter: failed while writing file`. `sigma56_sv` alone holds 364 GB
(188 WP + 176 classical); per WP run 40 GB = 25 GB `raw/vti` + 6 × 2.5 GB ckpt.

**Repo-wide there are 117 `ckpt_step*` dirs totalling 248 GB** — sigma56_sv 65,
slab_ks_wrap 26, cylindrical_jellium/proximity_ladder 16, wp_highdensity_sv 5,
channeling_twin 2. Interior checkpoints of a FINISHED run are kill-insurance
only; the final checkpoint is what preserves extendability.

**AWAITING USER DECISION** on pruning scope (narrow ~96 GB / medium ~150 GB /
wide ~200 GB). NOTHING DELETED. Do not delete without an explicit answer.

σ = 5 remaining: `cl_s5p0_v2p0` 2702/4360, `cl_s5p0_v2p5` 2700/3488,
`cl_s5p0_v3p0` 2680/2907, `s5p0_v2p0` 1973/4360 (all resumable);
`s5p0_v2p5` 131/3488 (no ckpt), `s5p0_v3p0` + `s5p0_v3p5` never started,
`cl_nocap_s5p0_v3p0` missing, `cl_nocap_s6p0_v3p0` 2680/2907.

### ⚠ RESUME BUG — orphan VTI frames abort every checkpoint resume (FIXED)

Every one of the 5 resume jobs submitted 2026-08-03 ~12:40 died ~45 s in with:

```
what():  VTIImageDataWriter: file already exists and overwrite=false:
         results/s5p0_v2p0/raw/vti/density_total/density_t001750.vti
```

**Cause.** Checkpoints are written every CKPT steps (581–872) but VTI frames every
SAVE steps (8–14). A run killed BETWEEN checkpoints therefore leaves frames for
steps the checkpoint does not cover. On resume the writer — deliberately
`overwrite=false` so segment output never clobbers earlier data — hits the first
such frame and throws. Orphan counts here: `s5p0_v2p0` 70, `cl_s5p0_v2p5` 51,
`cl_s5p0_v3p0` 36, `cl_nocap_s6p0_v3p0` 36, `cl_s5p0_v2p0` 7.

**Fix.** `shared/bin/prune_orphan_frames.sh` — reads `last_step` from
`rt_state.txt` and deletes `raw/vti/**/*_t<step>.vti` with step > last_step.
Wired into `run-s56-wp.slurm` and `run-s56-cl.slurm`, guarded on `LJ_RESUME=1`.
Safe because propagation from a checkpoint is deterministic: the resumed segment
recomputes exactly those steps. **The CSV side needed nothing** — `e_absorbed._concat`
already resolves overlap via `drop_duplicates(subset="step", keep="last")`.

**Generalise this**: any run definition following
`.claude/rules/final-timestep-checkpoint.md` with frame cadence finer than
checkpoint cadence has the same latent bug. Check before relying on any resume.

### ✅ THE "RESUME HANG" WAS A BAD NODE — gpu-q-2. Earlier diagnosis below was WRONG.

**Corrected 2026-08-03 16:45.** Every `loading electrons from …` hang traces to the
NODE, not to checkpoints and not to the classical half:

| job | node | hung at | duration |
|---|---|---|---|
| finalizer 32668217 | **gpu-q-2** | `cl_s5p0_v3p0/checkpoint` | 5 h 17 m |
| 32691777_2 | **gpu-q-2** | `cl_s5p0_v3p0/checkpoint` | 2 h 21 m |
| 32691777_0 | gpu-q-25 | `cl_s5p0_v2p0/checkpoint` | 2 h 23 m |
| 32699977_2 | **gpu-q-2** | the GROUND STATE | 1 h 09 m |

**The decisive one is 32699977_2**: a FRESH run (`resume: 0`) that hung loading the
*ground state* — no checkpoint involved at all. So the common factor is a large
read on gpu-q-2, not the resume path. Concurrently, every job on gpu-q-1/3/9/10/24
progressed normally.

**Fix: `sbatch --exclude=gpu-q-2,gpu-q-25`.** Use it on every submission in this
campaign until the nodes are known good. Classical resumes are NOT broken; the
`LJ_RESUME=0` workaround recorded below was treating a symptom.

### ⚠ STALE SEGMENT FILES CONTAMINATE A FRESH RESTART

A resume that dies still leaves `observables.from<N>.csv` / `interactions.*` /
`projectile.*`. `e_absorbed._concat` merges base + segments and dedups
`keep="last"`, so those stale rows SURVIVE a later fresh restart and sit at step
numbers the new run has not reached yet. Two consequences, both silent:

1. Progress is over-reported — `cl_s5p0_v2p0` read 2618/4360 (60 %) while the
   fresh run was actually at step 1729.
2. The finished dataset would be a MIX of the aborted run and the new one.

Quarantined 2026-08-03 by renaming to `*.aborted` (reversible) in
`cl_s5p0_v2p0`, `cl_s5p0_v2p5`, `cl_nocap_s6p0_v3p0`.
**Rule: before restarting a run FRESH, move its `*.from*.csv` aside.**

### ⚠ SUPERSEDED — "classical checkpoint resume is broken" (kept: the ruling-out still stands)

Three independent hangs, all on `loading electrons from '<classical run>/checkpoint/'`:
the finalizer (5 h 17 min, cl_s5p0_v3p0), then `32691777_0` (cl_s5p0_v2p0) and
`32691777_2` (cl_s5p0_v3p0) at **2 h 21 min with zero output**. The WP half loads
the same way in **14 s**. Not the disk (400+ GB free), not corruption (a suspect
and a known-good checkpoint are byte-identical in structure: 76 files,
592 / 17496000 / 74×34992000, total 2,606,925,072 B).

**WORKAROUND ADOPTED 2026-08-03:** classical runs are RESTARTED FRESH
(`LJ_RESUME=0`, jobs 32699977/32699978) rather than resumed. Fresh starts
demonstrably work, including into a directory that already holds frames (the
σ = 5 WP fresh restarts overwrote 131 steps of prior output without complaint).
Costs 2–3.6 h per run; a hang costs the whole allocation.

**Do NOT resume a classical run until this is diagnosed.** WP resumes are fine.
Prime suspect: the rolling `checkpoint/` is STALE vs `rt_state.txt` — for
cl_s5p0_v3p0 `checkpoint/` is 03:39 while `ckpt_step002324` is 05:10 and rt_state
says last_step=2324, so the loaded state and the START offset disagree. Next step
would be to point the resume at the NUMBERED ckpt_step dir instead.

### ⚠ SUPERSEDED — earlier note on the same stall (kept for the ruling-out work)

`32691777_0` (cl_s5p0_v2p0) and `32691777_2` (cl_s5p0_v3p0) sat >20 min on
`loading electrons from 'results/<run>/checkpoint/'` with no further output.
**Same signature as the finalizer that hung 5 h 17 min this morning — and on the
same run, cl_s5p0_v3p0.** So it is NOT the full filesystem (423 GB free now).

Ruled out: checkpoint corruption. `cl_s5p0_v3p0/checkpoint` and the known-good
`cl_s6p0_v3p0/checkpoint` are structurally identical — 76 files, sizes
592 / 17496000 / 74×34992000, total 2,606,925,072 bytes each, no zero-byte files.

Unexplained, and note the WP resume loaded the same way in **14 s**. Candidates:
(a) RDS contention (5 jobs hammering the FS), (b) the rolling `checkpoint/` being
STALE relative to `rt_state.txt` — for cl_s5p0_v3p0 `checkpoint/` is 03:39 while
`ckpt_step002324` is 05:10 and rt_state says last_step=2324, so the loaded state
and the START offset may disagree. **If it recurs, try resuming from the NUMBERED
ckpt_step dir instead of the rolling one**, or fall back to a fresh run
(LJ_RESUME=0) — fresh starts demonstrably work (3 σ=5 WP runs at 300+ steps).

### Notebooks — energy ledger added; 5/8 σ=6 built, rebuild re-chained

**All run notebooks now carry the energy curves** (user request 2026-08-03):
* §5 (a) ΔE_total raw + norm-corrected, (b) INQ components kinetic/Hartree/xc/
  external as deltas from t=0, (c) pairwise E_SS/E_PP/E_PS/E_SB/E_PB.
* §6 (d) the E_PS tail with the analytic `N_e/z` monopole overlaid — the plot that
  shows WHY the classical half needs correcting; (e) closure residuals on a log
  axis with `assert < 1e-6 Ha`.
* Classical notebooks also print the independent projectile-KE-loss cross-check.

Validated against real data before committing rebuilds to it: closure residuals
**5.0e-10 Ha (classical) / 5.4e-10 Ha (WP)**; E_PS(t_f) 8.55 eV classical vs
−0.00002 eV WP.

`skip_reason` now also treats **this builder and `s56_stopping.py` as inputs**, so
adding a section invalidates existing notebooks — otherwise the freshness guard
would report every notebook "up to date" and ship the old layout forever.

### Notebooks — 5/8 σ=6 built (OLD layout, will be rebuilt), earlier job 32691161

Built and verified (parsed: 6 code cells each, 0 unexecuted, 0 errors), with the
mandatory density GIF base64-embedded so it survives VTI deletion:
`run_{wp,classical}_s6_v2.0`, `run_{wp,classical}_s6_v2.5`, `run_wp_s6_v3.0`.
WP notebooks ~136 MB, classical ~60 MB (mostly GIF bytes); ~13 min each.

Remaining 4: `run_classical_s6_v3.0`, `run_{wp,classical}_s6_v3.5`,
`run_classical_s5_v3.5` — plus `twin_s6` and `synthesis`.

**`build_run_notebooks.py` made IDEMPOTENT + INTERRUPT-SAFE** (`skip_reason`):
skips INCOMPLETE runs (a notebook off a half-finished run renders a GIF that
stops mid-flight and an S off an unplateaued trace — it LOOKS finished, which is
worse than absent), and skips notebooks already fresh vs the run's observables.
Freshness alone is not enough, so `_is_executed()` also requires the file to
PARSE with every code cell carrying non-error outputs — otherwise a notebook
truncated by a walltime kill would keep its fresh mtime and be skipped forever.
`S56_NB_FORCE=1` overrides. Twin notebooks are skipped when a σ has no complete
runs.

**LESSON: run notebook builds under SLURM, not as a session background shell.**
Two builds were lost to session teardown before this was submitted properly.

### Finalizer cancelled — and a design flaw to fix before resubmitting

Job 32668217 hung **5 h 17 min** at `loading electrons from
'results/cl_s5p0_v3p0/checkpoint/'` — zero log output, zero bytes written. It
could not have succeeded: the resumed run would have died on its first VTI write
like all the others. Cancelled with its queued twin 32668218; queue now empty.

**FLAW:** `finalize.py::repair_one` uses `timeout = max(60, remaining)` where
`remaining` is the ENTIRE repair budget, so ONE wedged repair swallows the whole
36 h window instead of failing over. FIX BEFORE RESUBMIT: per-run cap at ~2× the
run's expected wall time.

---

## STATUS 2026-08-02 (later) — GS PASSED, both binaries COMPILE, chain made SELF-HEALING

### Ground state — DONE, and it validates the whole box change

| gate | value | verdict |
|---|---|---|
| E_GS (L_z = 105, dx = 0.40) | **207.183221873 Ha** | — |
| E_GS(L_z = 85, dx = 0.50) | 207.18322156141 Ha | **Δ = 3.1e-7 Ha** |
| ∫n dV | 100.0000000001 | PASS |
| r_s | 4.18147 | PASS |
| num_states | 74 | matches |

The 3.1e-7 Ha agreement was deliberately made INFO-only (a different box is a
different calculation), and it came back essentially identical — direct evidence
that the added 20 Bohr really is pure vacuum and the bath is the same physical
system as the 85-Bohr campaigns. Checkpoint:
`shared_gs/slab_n100_L35x35x105_dx0p4_per2`.

### Both new binaries COMPILE and propagate

- **WP smoke (32667764): all t=0 gates PASSED.** density std 4.242640578 vs the
  analytic 4.242640687 (dev 2.6e-6 %); max overlap with the occupied manifold
  3.894e-6 (want < 1e-3); ALIASING 0 % of z-momentum weight beyond Nyquist
  (k_Nyq−k0)/σ_p = 49.7.
- **Classical smoke (32667765): compiles and runs.** The 4-term
  `perturbations::sum(sum(bg, proj), sum(cap_lo, cap_hi))` — the flagged risk —
  is fine. Energy climbs 210.5757 → 210.6535 Ha over 10 steps as the projectile
  does work, as it should.

**MEASURED cost: 3.15 s/step (WP), 3.00 s/step (classical)** on an A100, against
2.75 s/step at L_z = 85. Wall clock per sweep is the longest single array task,
4360 × 3.15 s = **3.8 h**; whole campaign ~10–14 h, ~60 GPU-h.

### Autonomy hardening (user instruction: "ensure all the tasks happen autonomously")

The original stage 13 was a one-shot notebook job. That is not autonomous: a
SLURM chain gets runs LAUNCHED, not FINISHED, and a one-shot post-processor would
build a figure out of whatever happened to be on disk. Replaced (job 32667775
cancelled) with a bounded self-healing finalizer.

**New: `hypotheses/sigma56_sv/finalize.py` + `shared/bin/run-s56-finalize.slurm`**
(jobs **32668217** attempt 1, **32668218** attempt 2, chained `afterany`).

Each attempt: status → repair → build → report.
- **status** — all 26 expected runs (16 production, 2 CAP-free controls, 8 vacuum)
  checked against `STEPS_TARGET`, across resume segments.
- **repair** — missing or short runs are (re)run IN PLACE by invoking their own
  dispatcher as plain bash (`bash shared/bin/run-s56-wp.slurm 2`). Deliberately
  **no sbatch from inside a job** (no precedent in this repo) and **no duplicated
  launch logic** (the dispatcher stays the single source of truth). Resumes where
  a checkpoint exists, starts fresh otherwise. `ensure_binary()` rebuilds a
  missing binary via the smoke stage — but REFUSES to proceed if the binary
  builds and its t=0 gates then fail, because that is a defect, not a hiccup.
  Priority: production → controls → vacuum (the deposit does not need vacuum).
- **bounded, always** — `S56_REPAIR_BUDGET_S` (30 h) caps repair, leaving the
  rest of the 36 h job to build and report. This is the explicit fix for the
  2026-06-28 failure where a 9.5-hour finalizer polled for production that was
  never going to run.
- **report on EVERY path** — `hypotheses/sigma56_sv/CAMPAIGN_REPORT.md` (run
  status table, S table, CAP-cost table, artefact list, and resume instructions
  if it gave up short). Written even on total failure.
- **email is best-effort** — `send_run_email` is wrapped; a failure is logged,
  never raised.

Inspect at any time, changing nothing:

    cd ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv
    python finalize.py --status-only

**Email is NOT configured on this machine** — `~/.config/inqview/gmail_credentials.json`
is absent, so the campaign will report to disk only. One interactive command fixes
it for all future runs: `python -m inqview.email setup`. The disk report is the
authoritative record either way.

**Verified for the post-processing stage:** nbformat, nbclient, ipykernel,
matplotlib, numpy 2.4.6 (`np.trapezoid` present), pandas, imageio, PIL, scipy,
inqview and its density_gifs module all import in the venv, and the `python3`
Jupyter kernelspec is registered. `finalize.py --status-only` runs clean and
enumerates all 26 runs.

### Analysis path EXECUTED against real data — three bugs found and fixed

Syntax-checking the notebook builder was not enough. The whole analysis path was
run against the COMPLETED `s3p0_v3p0` run (same file layout as the new campaign),
which found three defects that would each have reached a notebook:

1. **`z_std` does not exist.** `wp_real_space_stats.csv` writes VARIANCES
   (`sigma_z2`) plus a CIRCULAR std (`sigma_z_circ`) — no column contains the
   substring "std", so the fallback `[c for c in pos.columns if "std" in c][0]`
   raised `IndexError`. Now uses `sigma_z_circ` (correct across the periodic z
   face) with `sqrt(sigma_z2)` as the fallback.
2. **`norm_slab` is classical-only.** The WP ledger writes `norm_wp`/`norm_total`
   instead, so the bath-count line crashed on every WP notebook. Now picks the
   column present and additionally reports `norm_wp` when it exists.
3. **`e_hartree_check` is a VALUE, not a residual** (it is
   `0.5*∫n_total·φ_total` ≈ 233 Ha, `interaction_energies.hpp:139`). The cell
   printed it as "closure residual", which would have read as a catastrophic
   ledger failure. The residual is the DIFFERENCE against the observables scalar;
   measured on real data it is **5.4e-10 / 5.2e-10 Ha**. The cell now computes the
   difference per half and asserts < 1e-6 Ha.

**Estimator validated end to end:** `s56_stopping.measure()` reproduces the
published `sigma_sweep_S_deposit.csv` value for σ=3, v=3.0
(`S_deposit_corrected = 0.3052865739`) to **2.8e-9**, and T₁−T₂ at t=0 comes out
2.2676 eV against the analytic 3/(4σ²) = 2.2676 eV.

**Vacuum-run layout confirmed:** vac runs write NO `observables.csv` (only
`wp_momentum_stats`, `wp_real_space_stats`, `cap_profile`, `energies`), which is
why `finalize._steps_done` falls back to `wp_momentum_stats` — verified against a
real `vac_s3p0_v3p0` directory rather than assumed.

---

## STATUS 2026-08-02 — machinery COMPLETE, 13-job chain SUBMITTED, GS running

### Goal (one line)

Matched classical + wavepacket twin pairs at σ_WP = 5 and 6 Bohr over
v = 2.0/2.5/3.0/3.5, combined with the existing σ = 0.5/2/3 results into one S(v)
figure, to find the width at which the classical and quantum projectiles stop
being distinguishable.

### Locked user decisions (2026-08-02)

| Decision | Value |
|---|---|
| σ_WP | **5 and 6**, both halves |
| Classical binary | **`dyn_direct` lineage** — direct erf/r potential, not the Poisson perturbation |
| Velocity grid | **4 points**, v = 2.0/2.5/3.0/3.5 |
| CAP | **ON in BOTH halves** + one CAP-free classical control per σ at v = 3.0 |
| Final figure | new σ = 5/6 twins + existing σ = 0.5/2/3 WP as-is, L_z difference in the caption |
| Box | L_z 85 → **105** (+20 Bohr vacuum), launch z −24 → **−27.5** |
| Storage | **not a constraint** (user, 2026-08-02) — vacuum controls therefore also save density frames |

### Why the box had to grow — do not re-derive this

σ_d(0) = σ_WP/√2 = 4.243 Bohr at σ = 6. The 85-Bohr box has **17.5 Bohr** between
the slab face (−12.5) and the CAP inner edge (−30); 3σ_d of clearance on each side
needs 25.5. There is **no launch point in the old box** where a σ = 6 packet is
clear of both — the ceiling there is σ_WP ≈ 3.0, which is exactly the documented
"σ = 3 is the practical ceiling". At L_z = 105 the CAP edge moves to −40, the gap
becomes 27.5 Bohr, and launch z = −27.5 leaves:

| σ_WP | σ_d(0) | to CAP | in σ_d | in CAP at t=0 | to slab | in σ_d | in slab at t=0 |
|---|---|---|---|---|---|---|---|
| 5 | 3.536 | 12.5 | 3.54 | **0.020 %** | 15.0 | 4.24 | **0.001 %** |
| 6 | 4.243 | 12.5 | 2.95 | **0.16 %** | 15.0 | 3.54 | **0.020 %** |

Both at or below the 0.23 % t=0 CAP loss already accepted at σ = 3.
n0, r_s = 4.183, slab thickness, L_xy, dx, dt and the CAP are all unchanged — the
extra 20 Bohr is pure vacuum.

### Why σ = 5/6 and not more of σ = 2/3

σ_d(t) = √(σ²/2 + t²/2σ²), so σ_WP sets both the width and the spreading rate.
Over the in-slab transit at the new launch:

| σ_WP | growth ×, v 2.0→3.5 | σ_eq = √2⟨σ_d⟩ |
|---|---|---|
| 0.5 | ×3.2 | 24 → 14 |
| 2 | ×2.7 → 2.2 | 6.4 → 4.0 |
| 3 | ×1.9 → 1.4 | 5.1 → 3.8 |
| **5** | **×1.23 → 1.08** | **5.74 → 5.26** |
| **6** | **×1.12 → 1.04** | **6.45 → 6.15** |

σ = 5/6 are the first widths where the packet is effectively constant-width AND
its label agrees with its time-average — the condition under which "classical and
quantum agree at width σ" is even a well-posed statement.

**Free collapse test built into the campaign:** σ = 6 at v = 2.0 has σ_eq = 6.45;
the existing σ = 2 at v = 2.0 has σ_eq = 6.35. Same time-averaged width reached
completely differently. If their S agree, time-averaged σ is a valid collapse
variable. Wired into `synthesis.ipynb`.

### Step counts — the formula, so they are reproducible

    N_STEPS = round( 4.36 * (|launch_z| + L_z/2) / (v * dt) ),  dt = 0.04

Calibrated on the recorded 3623 at v=2.0/z0=−24/L_z=85 (this gives 3624). New box:

| idx | v | N_steps | t (a.u.) | save/ | wf/ | ckpt/ |
|---|---|---|---|---|---|---|
| 0 | 2.0 | 4360 | 174.4 | 14 | 43 | 872 |
| 1 | 2.5 | 3488 | 139.5 | 12 | 35 | 698 |
| 2 | 3.0 | 2907 | 116.3 | 10 | 29 | 581 |
| 3 | 3.5 | 2491 | 99.6 | 8 | 25 | 498 |

Packet clears the slab by t = 20 a.u., reaches the CAP by t ≈ 34; the rest is
plateau time. Transverse images do not overlap until t_ov = 32.8 (σ=5) / 34.0
(σ=6) a.u., so the **entire transit is transversely clean** — unlike σ = 0.5,
where the slab and clean windows did not intersect at all at v = 2.0/2.5.

### Files written (all new)

**Config** — `ResearchProject/systems/localised_jellium/shared/configs/`
- `slab_n100_L35x35x105.hpp` — `SlabN100_L35x35x105`. Carries the CAP fractions
  (`CAP_WIDTH_FRAC = 0.119047619048`, `CAP_MID_FRAC = 0.440476190476`,
  `CAP_Z_INNER = 40.0`) and `LAUNCH_Z_BOHR = -27.5` as constants.

**Run machinery** — `ResearchProject/systems/localised_jellium/scripts/sigma56_sv/`
- `gs/run.cpp` — GS. E_GS gate REMOVED (a different box is a different
  calculation); gates instead on ∫n dV = 100 and r_s = 4.183.
- `wp/run.cpp` — clone of the validated `wp_highdensity_sv/wp/run.cpp` with only
  Cfg + launch + step defaults changed. All t=0 gates, the orthogonalisation
  budget, the pairwise ledger and checkpoint/resume are byte-identical logic, so
  nothing validated for the 85-Bohr campaigns needs re-validating.
- `classical/run.cpp` — from `classical_highdensity_sv/dyn_direct/run.cpp` plus
  **four** changes: (1) two absorbing bands, (2) Cfg-backed L_z=105 defaults,
  (3) frames to the CANONICAL `raw/vti/density_total` with `frames/total` as a
  symlink, (4) interior + retained checkpoints.
- `vac/run.cpp` — from `wp_highdensity_sv/cap_check/run.cpp`, with **L_z made an
  env knob** (it was hard-coded 85).

**Dispatchers** — `shared/bin/`
`run-s56-gs.slurm`, `run-s56-wp.slurm`, `run-s56-cl.slurm`, `run-s56-vac.slurm`,
`run-s56-notebooks.slurm`, `submit-sigma56-sv.sh`.

**Analysis** — `ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv/`
- `s56_stopping.py` — adapter over `slab_ks_wrap/e_absorbed.py::measure_dir` (the
  validated engine; NOT reimplemented). Run naming, E_GS read from the GS
  run_summary, `complete` check, and the dispersion geometry (`sigma_eq`).
- `build_sv_figure.py` — the two figures + `s56_S_summary.csv` + `s56_cap_cost.csv`.
- `build_run_notebooks.py` — 16 per-run notebooks, 2 twin notebooks, 1 synthesis.

### Two design choices worth not re-litigating

1. **The CAP-free classical control uses the SAME binary**, via `LJ_CAP_ETA=0`
   (a zero-η absorbing band is an exact identity). A separate CAP-free binary
   would drift away from production; this cannot.
2. **The classical frame-path fix is at the source.** The ancestor wrote only
   `frames/total/` while the notebook builders read `raw/vti/density_total/` —
   that mismatch silently produced 8 classical notebooks with NO density GIFs in
   the slab_ks_wrap campaign, and nothing warned. The new binary writes the
   canonical path and symlinks the legacy one.

### The chain — SUBMITTED 2026-08-02, jobs 32667763–32667775

| # | stage | job | state at write |
|---|---|---|---|
| 1 | gs (L_z=105, dx=0.40) | **32667763** | **RUNNING** (SCF converging, 74 states) |
| 2 | wp smoke σ6 (BUILDS WP binary + t=0 gates) | 32667764 | pending afterok 1 |
| 3 | cl smoke σ6 (BUILDS classical binary) | 32667765 | pending afterok 1 |
| 4 | wp sweep σ6, array 0–3 | 32667766 | pending afterok 2 |
| 5 | cl sweep σ6, array 0–3 | 32667767 | pending afterok 3 |
| 6 | vac σ6 (4 CAP-only baselines) | 32667768 | pending afterany 4 |
| 7 | cl nocap σ6, v = 3.0 | 32667769 | pending afterany 5 |
| 8 | wp smoke σ5 (t=0 gates at the other width) | 32667770 | pending afterany 4 |
| 9 | wp sweep σ5, array 0–3 | 32667771 | pending afterok 8 |
| 10 | cl sweep σ5, array 0–3 | 32667772 | pending afterany 5 |
| 11 | vac σ5 | 32667773 | pending afterany 6 |
| 12 | cl nocap σ5, v = 3.0 | 32667774 | pending afterany 7 |
| 13 | notebooks + figures | 32667775 | pending afterany 9,10,11,12 |

`scancel 32667763 32667764 32667765 32667766 32667767 32667768 32667769 32667770 32667771 32667772 32667773 32667774 32667775`

**Cost projection (WARN, not a gate — `.claude/rules/checkpoint-dont-block.md`):**
~60 GPU-h, ~15–20 h wall clock. Grid grew ×1.24, so ~3.4 s/step against the
measured 2.75 s/step at L_z = 85. Every run checkpoints every N/5 steps; a kill
costs at most one interval. Extend or resume with:

    sbatch --export=ALL,LJ_SIGMA=6.0,LJ_RESUME=1 shared/bin/run-s56-wp.slurm 0

### Verified vs NOT verified

**VERIFIED (this session)**
- Dispersion / clearance arithmetic in the plan, reproduced independently by
  `s56_stopping.py`'s `__main__` self-test (σ_eq, t_ov, transit windows).
- The step-count formula reproduces the recorded 3623 to within one step.
- All six dispatchers pass `bash -n`.
- All 14 inqkit headers the classical binary includes exist;
  `compute_coulomb_direct` really returns `norm_slab` / `norm_p`.
- All four notebook templates validate under `nbformat` and **every generated code
  cell compiles** (0 syntax errors).
- `build_sv_figure.legacy_wp()` reproduces the σ_eq table by hand-check
  (σ = 2, v = 2.0 → 6.354).

**NOT VERIFIED**
- **Nothing has compiled or run yet.** The three new binaries have never been
  built; the smoke stages (2, 3) are the first compile. `perturbations::sum`
  nesting to four terms in the classical binary is the most likely failure point.
  A compile error costs one short job — `afterok` stops the sweeps.
- Any physics at σ = 5 or 6.
- E_GS at L_z = 105 (job 1 in flight).
- Whether the classical CAP costs a little or a lot (stages 7, 12 measure it).

### BLOCKER — needs the user

`S_of_v_v2_timeavg_sigmar.png` **and its plotting script are not on this machine.**
`hypotheses/classical_highdensity_sv/dyn_direct/` holds only `S_of_v_direct.csv`
and two notebook builders; `scripts/classical_highdensity_sv/{dyn,dyn_direct}/results/`
are empty. Searched by filename across `/rds/user/skcb2/hpc-work` and
`/home/skcb2`, and through git history — absent. The user offered to transfer it.

`build_sv_figure.py` therefore follows the PROJECT standard plus the one thing the
reference's filename states unambiguously (points on a time-averaged-σ axis), and
says so in its module docstring. **Reconcile the design when the file arrives.**
This blocks only the final figure's styling — every run and every other
deliverable proceeds.

### Exact next steps

1. Watch job 32667763 → confirm `run_completed = true`, ∫n dV = 100.000,
   r_s = 4.183, and note E_GS (it is NOT the 85-Bohr 207.18323 Ha).
2. Read the σ = 6 WP smoke (32667764) t=0 gate block — expect
   `ortho removed weight` well under 3 %, `T1-T2 = 3/(4σ²) = 0.0208 Ha`,
   `sigma_pz^2 = 1/(2σ²) = 0.0139`.
3. Read the classical smoke (32667765) — this is the first compile of the
   4-term `perturbations::sum`.
4. When the sweeps land, check `complete` on every point before quoting any S.
5. Get the reference PNG + script from the user; reconcile `build_sv_figure.py`.

### Out of scope (user-decided)

- Re-running σ = 0.5/2/3 at L_z = 105.
- Classical twins at σ = 2 and 3.
- v = 4.0/4.5 (aliasing-free at these widths, so recoverable later at no risk).
- Reproducing the σ = 0.5 classical benchmark — raw data lost with the
  `/local/data/public` machine.
