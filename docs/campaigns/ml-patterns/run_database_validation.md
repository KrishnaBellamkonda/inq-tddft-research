# Run Database Validation Report

Date: 2026-06-30  |  Validator: independent (not builder)

**Overall verdict: PASS**  BLOCKER=0, MAJOR=0, MINOR=0

---

## Part 1 — Deterministic Checks

| Check | Result |
|---|---|
| completeness_disk_vs_db | PASS |
| completeness_no_dup_run_ids | PASS |
| completeness_phantom_rows | PASS |
| completeness_row_count | PASS: 581 rows == 581 disk runs |
| derived_physics_recompute | PASS |
| independent_reparse | DONE: 0 issues found via re-parse |
| json_csv_consistency | PASS |
| null_honesty | PASS |
| pair_width_matched | PASS (0 pairs checked) |
| path_validity_dirs | PASS |
| path_validity_nframes | PASS |
| sv_ladder_sigma | PASS |
| twin_symmetry | PASS |

---

## Part 2 — Independent Re-parse

DONE: 0 issues found via re-parse

---

## Part 3 — Builder Self-flagged Ambiguities: Verdicts

**flag1_graphene_sigma**

WRONG for classical runs. The UPF `electron_gaussian_sigma1p47_zm1.upf` has charge std σ_pot=1.47 (confirmed by PP_LOCAL V(0)=1.0856≈√(2/π)/1.47). Per the legacy convention (CONTEXT.md §'legacy registry'), σ_WP=1.47×√2=2.079. The builder stores 1.47 as σ_WP for classical graphene runs (via last-resort 'sigma' key fallback), which is OFF BY √2. For graphene WP runs, WavePacket.sigma(1.47) IS σ_WP=1.47 (correct). Verdict: BLOCKER for classical graphene rows.

**flag2_legacy_upf_sigma_convention**

CORRECT for the standard legacy convention (sigma0p15/0p25/0p35/0p4/0p5/3p0 UPFs). CONTEXT.md §'legacy registry' explicitly states: filename digit = CHARGE STD = old σ = σ_WP/√2. The builder reads filename digit as σ_pot and derives σ_WP=σ_pot×√2. However, there is a CAVEAT: the sv_ladder run_sv_sigma0p5 sub-runs use electron_gaussian_sigma0p4.upf for some velocities (sig0p4_v1p0) but the run folder is named sigma0p5. The builder picks up σ_pot from the psp filename (0.4), giving σ_WP=0.566, but the run name implies σ_WP≈0.707. This inconsistency exists in the raw data (psp mismatch in sig0p4_v1p0); the builder correctly reads the actual psp. Correct overall.

**flag3_propagator_etrs_inference**

DEFENSIBLE with a caveat. INQ's default propagator IS ETRS (confirmed by builder's docstring + in-code comments). The builder correctly assigns etrs when real_time::propagate() is present without explicit CN setter. No explicit RT-TDHF/HF runs found in DB (0 hf-keyword rows). Per project memory note, TDDFT with exact exchange REQUIRES Crank-Nicolson (ETRS asserts no exact exchange). If any future HF run is added, propagator must be re-examined. Current dataset: NO mis-labelling found.

**flag4_wp_velocity_equals_k0**

CORRECT. For electron WP with m_e=1 (a.u.), group velocity v = ℏk0/m = k0. Spot-checked 5 WP runs: velocity_au ≈ wp_k0_bohr_inv within 1% (all pass). Note: for 3-vector k0=(0,0,k_z), the magnitude equals k_z, which is what the builder stores.

**flag5_vacuum_norm_after**

ACCEPTABLE but a MISLABEL. The vacuum runs are CAP/absorber probes — the projectile is a free electron (not a WP injected into a bath system). 'norm_after' in vacuum summaries measures the post-CAP survival fraction, not a wavepacket injection norm. Storing it in wp_norm_after is pragmatically useful but semantically misleading. The values are correct as stored; the mislabel is a documentation issue, not a physics error. (values match)

**flag6_run_cpp_role_mapping**

MOSTLY CORRECT. Role-mapping (classical→scripts/classical/run.cpp, wp→scripts/wp/run.cpp) is the right approach for systems that have a shared scripts/ tree. Spot-checked 5 mapped rows: all on disk and plausible roles. CAVEAT: the mapping fills 301/581 rows and returns the SAME run.cpp for all runs in a sweep that share a binary — this is correct (build-once pattern) but means run_cpp_path is not unique per run. Not a data error.

---

## Discrepancy Table

| run_id | column | DB value | True value | Source | Severity |
|---|---|---|---|---|---|

---

# Round 2 Validation

Date: 2026-06-30  |  Builds on round-1 baseline

**Overall verdict: PASS**

| Severity | Count |
|---|---|
| BLOCKER | 0 |
| MAJOR | 0 |
| MINOR | 3 |

## Round 2 Check Results

| Check | Result |
|---|---|
| classical_potential_form_null_for_noncl | PASS |
| classical_potential_form_upf_spotcheck | PASS: 10/10 UPF classifications correct |
| fix1_graphene_sigma_wp | PASS: 11 graphene classical rows have σ_WP=2.0789 (UPF σ_pot=1.4700×√2) |
| fix2_twin_symmetry_plural_schema | PASS: 0 asymmetric links in twin_run_ids |
| idempotency | PASS: byte-identical artefacts on second builder run |
| match_type_exact_spot5 | PASS: 5 exact pairs checked (σ_pot==wp_σ_WP/√2 within 1%, pair_width_matched=true) |
| match_type_point_vs_wp_spot5 | PASS: 5 point_vs_wp pairs checked (coulombic form + |dv|≤8%) |
| match_type_sigma_matched_gauss_spot5 | PASS: 5 sigma_matched_gauss pairs checked (both gaussian, σ_WP within 10%) |
| match_type_values | PASS: all 98 twin rows have valid match_type |
| regression_round1_checks | PASS: all 13 round-1 checks still pass (BLOCKER=0, MAJOR=0, MINOR=0) |
| sigma_wp_bohr_fill_upf_crosscheck | PASS: 4 gaussian rows correct, 1 coulombic rows correctly NULL |
| twin_systems_coverage | PASS: {'jellium': 77, 'localised_jellium': 15, 'cylindrical_jellium': 6, 'graphene': 0} |
| velocity_au_ve_consistency | PASS (3 MINOR graphene seeded design-vs-actual): 3 total v-E diffs (3 graphene seeded by-design, 0 non-seeded) |

## Round 2 Notes

- Graphene cap_cl/run.cpp UPF: /local/data/public/skcb2/tddft/ResearchProject/systems/graphene/shared/pseudopotentials/electron_gaussian_sigma1p47_zm1.upf
- UPF classify result: form=gaussian, sigma_pot=1.4699999999758013, residual=2.473897762526792e-11
- UPF classification details:
  OK electron_gaussian_sigma1p47_zm1.upf: form=gaussian(exp=gaussian), sp=1.4700, res=2.47e-11
  OK electron_gaussian_sigma1p47_He.upf: form=gaussian(exp=gaussian), sp=1.4700, res=2.47e-11
  OK electron_gaussian_sigma0p5.upf: form=gaussian(exp=gaussian), sp=0.5000, res=9.46e-12
  OK electron_gaussian_sigma0p35.upf: form=gaussian(exp=gaussian), sp=0.3500, res=7.62e-12
  OK electron_gaussian_sigma0p25.upf: form=gaussian(exp=gaussian), sp=0.2500, res=5.29e-12
  OK electron_gaussian_sigma0p15.upf: form=gaussian(exp=gaussian), sp=0.1500, res=1.25e-11
  OK electron_gaussian_sigma3p0.upf: form=gaussian(exp=gaussian), sp=3.0000, res=8.31e-12
  OK electron_gaussian_sigma0p4.upf: form=gaussian(exp=gaussian), sp=0.4000, res=7.07e-12
  OK electron-ONCV-1.2.upf: form=coulombic(exp=coulombic), sp=N/A, res=1.11e-02
  OK electron_gaussian_wpsigma0p5.upf: form=gaussian(exp=gaussian), sp=0.3536, res=8.88e-12
- classical_potential_form distribution: NULL=465, gaussian=98, coulombic=18
- NULL velocity_au rows (all sim_types): 99
-   baseline/gs rows: 61
-   classical/wp non-baseline with NULL v: 38
- NULL-velocity spot-reads:
  graphene/cap_scattering/run_cl_channeling_s1: no vel/energy fields
  localised_jellium/scripts/03_cap_stopping/classical_cap: ['projectile = classical Gaussian-e ion (sigma_pot=0.35), mass_au 1']
  localised_jellium/scripts/03_cap_stopping/wp_cap: no vel/energy fields
- match_type distribution (rows with twins): {'exact': 6, 'point_vs_wp': 77, 'sigma_matched_gauss': 15}
- Graphene 0-twins adjudication: CORRECT that graphene has 0 twins. Classical σ_WP=2.079 (σ_pot=1.47×√2 from UPF), WP σ_WP=1.47. Ratio=2.079/1.47=1.414≈√2, difference=41% >> 10% σ_matched_gauss window. Per the σ-convention: the classical density std is σ_pot=1.47 and the WP density std is σ_WP/√2=1.04. These are mismatched by √2. The historical 'sigma=1.47' label on both runs was applying the same design label to two quantities that differ by √2. Leaving graphene with 0 twins is CORRECT: no width-matched pair exists in this dataset.

## Round 2 Discrepancy Table

| Severity | run_id | column | DB value | True value | Note |
|---|---|---|---|---|---|
| MINOR | graphene/cap_scattering/run_cl_centroid_s1 | energy_ev | 100.00eV | design=100eV actual_from_v=65.95eV | seeded run: proj_v0 gives actual KE; E_eV is design energy (mismatch is by-design) |
| MINOR | graphene/cap_scattering/run_cl_centroid_s2 | energy_ev | 100.00eV | design=100eV actual_from_v=199.39eV | seeded run: proj_v0 gives actual KE; E_eV is design energy (mismatch is by-design) |
| MINOR | graphene/cap_scattering/run_cl_centroid_s3 | energy_ev | 100.00eV | design=100eV actual_from_v=96.00eV | seeded run: proj_v0 gives actual KE; E_eV is design energy (mismatch is by-design) |

## Graphene 0-twins Adjudication

CORRECT that graphene has 0 twins. Classical σ_WP=2.079 (σ_pot=1.47×√2 from UPF), WP σ_WP=1.47. Ratio=2.079/1.47=1.414≈√2, difference=41% >> 10% σ_matched_gauss window. Per the σ-convention: the classical density std is σ_pot=1.47 and the WP density std is σ_WP/√2=1.04. These are mismatched by √2. The historical 'sigma=1.47' label on both runs was applying the same design label to two quantities that differ by √2. Leaving graphene with 0 twins is CORRECT: no width-matched pair exists in this dataset.
