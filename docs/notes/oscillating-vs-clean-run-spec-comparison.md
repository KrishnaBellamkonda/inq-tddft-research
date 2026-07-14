# Spec comparison: oscillating family (slab_n52) vs clean run (qsp_phase3 p3_wp)

**Purpose.** Complete specification differencing between the energy-oscillation run
family (diagnosis + cap_fix + pbc_open_z campaigns; slab_n52 / effmass lineage) and
the clean-decay existence proof `p3_wp` (qsp_phase3; flat to t = 100 with a CAP on),
to isolate what gates the drain-then-rise artifact. Produced 2026-07-14 by an
independent advisor agent (TDDFT/stopping methodologist charter) reading the run
definitions, config headers, run summaries, and engine source; reviewed and filed
by the main session. Every spec claim carries file refs; inferences are labelled.

**Standing facts this analysis builds on (do not re-derive):**
- Diagnosis campaign (conf 0.90): the artifact is CAP-gated bookkeeping; dominant
  channel = the norm-divided kinetic ledger term (`inq/src/hamiltonian/energy.hpp:55`)
  rising as the CAP filters below-average-KE density.
- cap_fix: t_min (first turn) drifts with CAP config (21.6 → 27.8 → ~33 → 36.4 →
  >48); every sufficiently long two-sided run turns; wrap η=−2 w40 shows no turn
  to t = 48 (period-lengthening reading: "no turn observed", not "no oscillation").
- pbc_open_z (Arm B): periodicity 2 (open-z Poisson, matched p2 GS) reproduces the
  oscillation with t_min shifted 0.2 a.u.; amplitudes slightly larger. Electrostatic
  z-periodicity is neither cause nor clock. p2 GS tail under the CAP = 2e-11 (6
  orders tighter) with oscillation unchanged → static GS tail exonerated; the
  feeder is dynamically spilled slow density.
- `p3_wp` ran 2500 steps (t = 100), fully periodic, CAP on (η = −0.7), and decays
  to a stable plateau — a genuine no-oscillation existence proof.

Sources (all under `ResearchProject/systems/localised_jellium/` unless noted):
oscillating = `scripts/cap_fix/run.cpp`, `scripts/muon_mass_fork/effmass_sigma1/wp/run.cpp`,
`shared/configs/slab_n52_L40x40x80.hpp`, `scripts/cap_fix/results/run01_baseline_two_eta0p2/run_summary.txt`;
clean = `scripts/qsp_phase3/wp/run.cpp`, `shared/configs/slab_n82_L50x50x90.hpp`,
`scripts/qsp_phase3/wp/results/p3_wp/run_summary.txt`.

---

## 1. Complete spec-difference table

| # | Parameter | Oscillating | Clean | Same? |
|---|---|---|---|---|
| 1 | Cell | 40×40×80 Bohr (`slab_n52...hpp:34-36`) | 50×50×90 Bohr (`slab_n82...hpp:46-48`) | DIFF |
| 2 | Periodicity | 3 (default; env `EM_PERIODICITY`, `cap_fix/run.cpp:86,102`); Arm B also ran p2 | 3 (`.periodic()`, `qsp.../run.cpp:80`) | same (established non-discriminator) |
| 3 | Spacing / E_cut | 0.33333 Bohr → E_cut 44.4 Ha, k_Nyq 9.42 (`slab_n52:37`) | 0.50 Bohr → E_cut 19.7 Ha, k_Nyq 6.28 (`slab_n82:49`) | DIFF |
| 4 | Grid points | 120×120×240 ≈ 3.46 M | 100×100×180 = 1.8 M | DIFF (follows 1,3) |
| 5 | Slab | half-width 12.5, axis z, centre 0 (`slab_n52:41-43`) | identical (`slab_n82:52-54`) | same |
| 6 | Slab edge profile | erfc-softened, EDGE_WIDTH 1.0 (`slab_n52:44`) | sharp Θ, EDGE_WIDTH 0.0 (`slab_n82:55`) | DIFF |
| 7 | N electrons / n0 / r_s | 52 / 1.300e-3 / r_s 5.68 (`slab_n52:47-49`) | 82 / 1.312e-3 / r_s 5.67 (`slab_n82:58-60`) | DIFF in N; density ~same (ω_p 0.128 both, T_p ≈ 49 a.u.) |
| 8 | Extra states / total | 10 → 36 states, wp_idx 35 | 20 → 61 states, wp_idx 60 | DIFF |
| 9 | Temperature | 0.00862 eV | 0.00862 eV | same |
| 10 | GS checkpoint | `shared_gs/slab_n52_L40x40x80_dx0p333` | `shared_gs/slab_n82_L50x50x90` | DIFF (necessarily; SCF-discipline parity UNVERIFIED — see §5) |
| 11 | WP σ | 1.0 (`cap_fix:64`) | 0.5 (`slab_n82:69`) | DIFF |
| 12 | WP k0 | 5.693 (`cap_fix:65`) | 2.711 (`slab_n82:71`) | DIFF |
| 13 | WP mass | m_eff 2.10 via `electrons.inverse_mass()[0][wp_idx]=0.476190` (`cap_fix:66,126`) | m = 1 (inverse_mass never touched) | DIFF — **the mass fork** |
| 14 | WP velocity | k0/m = 2.711 a.u. | k0 = 2.711 a.u. | **same (matched by design)** |
| 15 | WP kinetic energy | k0²/(2m) = 7.72 Ha ≈ **210 eV** | **100 eV** | DIFF (2.1×) |
| 16 | Velocity spread | σ_p = 1/(2σ) = 0.5; Δv/v ≈ 0.09 (narrowband) | σ_p = 1.0; Δv/v ≈ 0.37 (dispersing) | DIFF (inference from σ, m) |
| 17 | Chirped focus | `focus_z(4.0, 2.10)` (`cap_fix:117,122`) | none (`qsp:110-112`) | DIFF |
| 18 | Launch z / standoff | −16.5 = 4 Bohr (4σ) from face; 3σ front reaches +22 vs CAP inner 25 at exit | −23.75 = 11.25 Bohr (22.5σ) from face, equidistant from CAP | DIFF |
| 19 | Orthogonalisation / injection / occ 1.0 | `orthogonalise_against_occupied` → last extra state (`cap_fix:118-124`; `wavepacket.hpp:408`) | identical mechanics (`qsp:110-113`) | same |
| 20 | Net cell charge after injection | −1 (53 e on 52-e bg) | −1 (83 on 82) | same |
| 21 | CAP profile family | two-sided sin² (stock) or wrap cos² (inqkit), env-selectable | two-sided sin² only (`qsp:124-125`) | same family |
| 22 | CAP η | default −1.0; swept −0.2…−4 | −0.7 hardcoded (`qsp:69`) | DIFF |
| 23 | CAP region / width per side | centre ±32.5, width 15 → region **[25,40]** | region **[35,45]**, width 10 | DIFF |
| 24 | Slab-face→CAP-inner gap | **12.5 Bohr** | **22.5 Bohr** | DIFF |
| 25 | CAP outer edge vs boundary | W→0 exactly at z = ±40 = L/2 | W→0 exactly at z = ±45 = L/2 | same (seam-transparent in both; wrap mode removes it in the family only) |
| 26 | CAP fraction of half-z | 15/40 = 37.5% | 10/45 = 22% | DIFF |
| 27 | dt | 0.04 | 0.04 (env override of 0.02 code default) | same |
| 28 | n_steps / horizon | 700–1200 (t = 28–48) | 2500 (t = 100) | DIFF (clean = true t=100 proof) |
| 29 | Propagator / theory / kpts | ETRS, LDA, Γ | ETRS, LDA, Γ | same |
| 30 | Engine code path | inq-study, per-state kinetic_factor ACTIVE | inq-study, kinetic_factor all-ones no-op guard (`inq-study/src/hamiltonian/ks_hamiltonian.hpp:156-167`) | DIFF in active path, same binary lineage |
| 31 | RT engine flags | plain real_time | + `.observables_current().observables_dipole()` (`qsp:197-198`) | DIFF (observation only, no back-action) |
| 32 | Write cadence / VTI | every 5, no VTIs | every 8, heavy VTI/wavefunction I/O | DIFF (I/O only) |
| 33 | SCF params | tol 1e-4 / mix_ndim 8 / mix 0.1 | identical | same |

**Record correction:** the `p3wp_run_notebook.ipynb` cell 9 caption "CAP inner faces
(|z|=25)" is a stale builder-template caption; the authoritative values are region
±35..±45 (run_summary + `qsp:69`).

## 2. Mass-fork ledger question — RESOLVED (exonerated as corruption, retained as amplifier)

The reported kinetic energy uses the per-state mass **consistently**:
`kinetic_expectation_value` (inq-study `ks_hamiltonian.hpp:313-323`) uses the same
per-state `kinetic_factor_ = −0.5·inverse_mass[ist]` as the apply kernels (set in
`real_time/propagate.hpp:73-75`; no-op guard when all masses = 1), and
`hamiltonian/energy.hpp` is byte-identical between `inq/` and `inq-study/`
(verified by diff). The norm division (`energy.hpp:50-57`, summed at `:83`)
applies uniformly to every state. So the fork does not independently corrupt the
ledger — but the forked WP state contributes a **7.7 Ha (210 eV) norm-divided
kinetic term**: a ~15% filtering shift of its surviving density is ~31 eV, exactly
the observed overshoot scale. Since classical/truncated runs reportedly also
oscillate, the WP term is an *amplifier*, not the sole channel.

## 3. Per-difference mechanistic assessment

(vs the norm-divided-kinetic mechanism with dynamically-spilled slow density as feeder)

| Rows | Difference | Assessment | Rank |
|---|---|---|---|
| 23/24/26 | CAP gap 12.5 vs 22.5 Bohr; region [25,40] vs [35,45] | sets when/how strongly spilled slow density meets the filter; v_F ≈ 0.34 a.u. → gap transit ≈ 37 vs 66 a.u. | **PRIME SUSPECT** |
| 22 | η −1 vs −0.7 | the established t_min clock knob, but −0.7 sits inside the family's swept range (−0.2 and −2 both oscillate) → cannot gate alone | plausible (co-factor) |
| 11–16 | projectile: 210 eV narrowband mass-forked vs 100 eV broadband | sets spill-pool size (4.0 vs 1.2 eV/electron, 3.3×) and CAP survival of the pulse | **PRIME SUSPECT** |
| 13 | mass fork as ledger corruption | refuted (§2); survives as amplitude amplifier | unlikely standalone |
| 17 | chirped focus | clean packet's Δv is 4× broader yet clean; classical twins oscillate without focus | unlikely |
| 18 | launch standoff 4σ vs 22.5σ | family's 3σ front (+22 vs CAP 25) puts ~0.1–1% tail into the CAP mid-run; WP-only channel, can't gate classical | plausible (secondary) |
| 3 | dx 0.333 vs 0.5 | adversarial: coarser (clean) grid would alias high-k ejecta into slow-looking k and make it WORSE — clean is clean | unlikely |
| 6 | edge width 1.0 vs 0.0 | anti-correlated (sharper edge = clean run) | exonerated |
| 8 | extra states 10 vs 20 | unoccupied states carry occ ≈ 0 in the ledger | irrelevant |
| 7 | N 52 vs 82 | only via per-electron deposition (folded into projectile suspect) | folded |
| 1 | box 80 vs 90 | only via gap/transits | folded |
| 25 | seam hole | present in BOTH (not the family/clean discriminator); a real amplitude/clock knob within the family (wrap results) | co-factor |
| others | periodicity, slab, T, v, injection, charge, dt, ETRS/LDA, SCF, I/O, flags | identical or no back-action | irrelevant |

## 4. Ranked verdict (advisor, adversarially self-checked)

**Suspect 1 (PRIME): CAP standoff geometry.** The absorber sits on the near-slab
dynamic-spill zone (gap 12.5 vs 22.5 Bohr). Spill arrival ≈ gap/v_F: ≈ 37 vs
≈ 66 a.u. — and the clean config's header *explicitly designed* the later arrival
(`slab_n82...hpp:10-11`). Consistent with: t_min tracking CAP config; classical
runs oscillating; p2≈p3; GS-tail exoneration. Strains (so geometry is necessary
but not sufficient): the clean run is flat to 1.5× its own arrival time, and wrap
w40 has a smaller gap yet turns latest → flux amplitude (S2) and filter character
(S3) co-determine visibility.

**Suspect 2 (PRIME): perturbation strength.** 3.3× more deposited energy per
electron feeds a larger slow-spill pool, and the family's narrowband packet
(Δv/v ≈ 0.09) partially SURVIVES a weak CAP (dwell ≈ 5.5 a.u.; survival at
η = −0.2 ≈ 33%) and reaches the boundary seam at t ≈ 20.8 — matching the earliest
t_min = 21.6. Explains early t_min rungs; the v_F spill arrival (S1) the late
ones. All transit arithmetic is labelled inference, not measured.

**Suspect 3 (plausible): CAP transfer function** (η·width reflectivity ×
seam-hole interplay): repeated slow-density filtering/reflection passes between
slab face and CAP = the oscillation; |η| sets the reflect/absorb split = clock.
A modulator given sufficient incident flux, not the root gate (clean run has an
in-range η and stays flat).

**Demoted explicitly:** mass-fork ledger corruption (§2), spacing (argument runs
the wrong way), edge softness (anti-correlated), extra states (occ ≈ 0),
periodicity + GS quality (established facts), chirp (classical twins oscillate).

## 5. Cheapest decisive experiments (existing infrastructure, env-only)

| Target | Run | Prediction if suspect gates | Prediction if not |
|---|---|---|---|
| S1 geometry | `EM_CAP_ETA=-0.7 EM_CAP_CENTER_BOHR=35 EM_CAP_WIDTH_BOHR=10` (region [30,40] = clean profile at the 80-box's max standoff), 1200 steps | t_min pushed far past ≈25 (η-interpolation expectation), amplitude collapses | t_min ≈ 25, family amplitude |
| S2 strength | `EM_INV_MASS=1.0 EM_K0=2.711 EM_SIGMA_WP=0.5` (clean projectile in the oscillating box; family CAP η=−1 [25,40]), 1200 steps | oscillation collapses/leaves window | oscillation persists at family amplitude (v unchanged) |
| S3 filter | zero-GPU: absorbed fraction at t_min from existing family `charge.csv` + clean `electron_number.csv`; then one run at `EM_CAP_ETA=-0.35` if needed | common absorbed-slow-fraction threshold at every turn (pure clock) | t_min at different absorbed fractions (recirculation doing real work) |

**Sharpest single discriminator: S2** — one run splits the two prime suspects
(S1 predicts persistence, S2 predicts collapse) and constrains S3 either way.
(~1–2 h per 1200-step run on one GPU; 3035 s / 700 steps measured.)

## 6. Unverified items (explicit)

- GS SCF-discipline parity between the two checkpoints (GS-generation scripts not
  read this pass); materiality low given the Arm-B matched-GS result.
- Whether the clean run's plateau would resolve a few-eV artifact at t ≈ 70–100
  (the expected clean-geometry amplitude is not independently bounded).
- All transit/survival arithmetic in §4 (v_F = 0.34 a.u., dwell/attenuation,
  arrival times) is inference from specs, not measured from run data.
- "Classical/truncated runs also oscillate" is taken from the original user note;
  those runs' summaries were not re-located in this pass.

---

**Related documents:** diagnosis `docs/campaigns/localised_jellium/energy-oscillation-diagnosis.md`
(+ notebook Parts II–IV); setup search `docs/campaigns/localised_jellium/cap-fix-experimentation.md`;
Arm B `docs/campaigns/localised_jellium/pbc-open-z-oscillation.md` +
`hypotheses/pbc_open_z/pbc_open_z_study.ipynb`; original phenomenon note
`docs/notes/localised-jellium-energy-oscillation-investigation.md`; glossary
entries in `CONTEXT.md` ("Period-lengthening reading", "PBC-vs-open-z channels").
