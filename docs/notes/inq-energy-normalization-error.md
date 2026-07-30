# INQ reports per-particle (norm-normalized) energy — the CAP "energy shoots up" artifact

Date: 2026-07-28. Status: root cause located in INQ source and verified in vacuum +
localised jellium. This note is the standing record; running hypothesis + validation
suite below.

## Symptom

In a real-time run with a Complex Absorbing Potential (CAP), as the CAP absorbs a
wavepacket (WP) — norm(t) → ~0 — the reported `energy_total` does **not** fall to ~0.
It stays near the initial value and often **rises above it** ("the energy shoots up").
Seen in vacuum (single WP) and in the localised-jellium CAP campaign (the plateau).

## Root cause (in the INQ source)

`inq/src/hamiltonian/energy.hpp:50–55`, the reduction used for the kinetic energy:

```cpp
static double occ_sum(Occupations const&, Array const& arr, Norms const& nor) {
    return gpu::run(gpu::reduce(arr.size()), 0.0,
        [occ, arr, nor] GPU_LAMBDA (auto ip) {
            return occ[ip] * real(arr[ip]) / real(nor[ip]);   // <-- divides by the orbital norm
        });
}
```
called as `kinetic_ += occ_sum(occupations, kinetic_expectation_value(phi), norms)`.

So INQ reports `sum_i occ_i * <psi_i|T|psi_i> / <psi_i|psi_i>` — the **intensive
(per-particle MEAN) kinetic energy**. For a normalized orbital `<psi|psi>=1` this is a
no-op. Under a CAP the WP orbital's norm decays, and dividing the (also-decaying)
numerator by the shrinking denominator **pins the reported energy at the mean** — which
even rises, because the CAP removes lower-energy/edge amplitude first, leaving a
higher-mean-energy remnant. The Hartree/external/xc terms are density-based (extensive,
built from `sum occ|psi|^2`) and DO track the absorbed WP correctly; only the kinetic
term carries the /norm. (`total()` = kinetic + hartree + external + non_local + xc +
exact_exchange + ion + ion_kinetic — only `kinetic` is normalized.)

## The fix

The extensive (physically-captured) energy un-normalizes the WP kinetic term:

```
E_extensive(t) = E_reported(t) * norm(t)                     (single WP, occ=1)
             or  E_reported(t) - occ_WP * e_kin_ha(t) * (1/norm - 1)   (multi-orbital)
             == E_reported(t) - occ_WP * e_kin_ha(t) * (1 - norm)      (equivalent)
```
where `e_kin_ha` = the WP per-particle mean KE (already logged per step in
`wp_momentum_stats.csv`) and `norm` from `wp_real_space_stats.csv` (or the momentum-
stats norm, normalized to its t0 value). `captured(t) = E0*norm0 - E_extensive(t)`.

## Evidence (all verified this campaign)

- **Not reflection:** vacuum cap runs at eta=-0.7 (W=25), -1.0 (W=30, full +z half),
  -3.5 (W=15) ALL give the same ~400 eV "residual" (E_reported - norm*E0). A reflection
  effect would depend strongly on eta and width. It does not.
- **Wrap is innocent:** a longer no-CAP run (t=16, ~2 periodic wraps) conserves
  `energy_total` to 0.15 meV — periodic re-entry is a clean unitary op (rules out
  "re-entry treated as a new projectile").
- **e_kin_ha == wavefunction:** the logged WP mean KE matches ½∫|∇ψ|²/∫|ψ|² computed
  from the saved complex `wavefunction_wp` frames to the digit (t0=120.4 eV, exact
  k0²/2 + 3/4σ²; tF=22.5 eV in jellium).
- **Fix decays cleanly:** vacuum cap_fulllen E_ext 402 → 1.65 eV (captured 400 eV =
  99.6% = norm absorbed); cap_better captured 386 eV (96%, 4% leaked). Physical.
- **Jellium plateau:** reported nocap-cap gap 93.5 eV → corrected 115.9 eV (the
  normalization ADDS ~22 eV; peak 44 eV mid-absorption). Only ~22 eV, not the full
  120 eV WP KE, because the WP SLOWS 120→22 eV in the bath (stopping; that kinetic is
  deposited into the bath and booked correctly by the density terms). So the jellium
  plateau is LARGELY PHYSICAL; the normalization is a real but sub-dominant systematic.

## Running hypothesis (to validate with the experiment suite)

> The CAP "energy shoots up" entirely because INQ reports the per-particle (norm-
> normalized) kinetic energy; the extensive energy `E_reported * norm` decays smoothly
> and physically as the CAP removes the packet. The effect is independent of CAP
> geometry (one/two-sided), strength, and width, and of periodic wrap; it depends only
> on the orbital norm decaying under a non-Hermitian potential.

Falsifiers: (a) if `E_ext = E_reported*norm` did NOT decay smoothly with the absorbed
norm; (b) if the "residual" scaled with eta/W (→ reflection, not normalization); (c) if
a norm-preserving absorber (e.g. mask that renormalizes) still showed the rise; (d) if
the rise appeared without any norm loss.

## Engine: this is in inq-study (the engine the CAP runs actually use)

Verified 2026-07-28: `inq-study/src/hamiltonian/energy.hpp` `occ_sum` is
BYTE-IDENTICAL to stock `inq` (the `occ[ip]*real(arr[ip])/real(nor[ip])` /norm line
is present in both; `diff` of the energy block is empty). Every CAP run in this
campaign builds against **inq-study** (`INQ_SOURCE=$ROOT/inq-study`; the vacuum
binary's CMakeCache links inq-study), because stock inq cannot compile a real-time
absorbing run (`double += complex`; see [[reference_stock_inq_cannot_compile_cap]]).
So the finding is not a stock-inq curiosity — it is in the exact engine that produces
our results, and the entire investigation (CAP runs, mask runs, jellium replica) is
run in inq-study. The normalization is NOT a modification introduced by inq-study;
inq-study only complexifies the scalar potential so the CAP's imaginary part reaches
the orbitals — the /norm in `occ_sum` is upstream INQ and is present in both trees.

## Validation status

The running hypothesis is being tested by the suite in
`docs/plans/cap-energy-normalization-validation.md`. Confirmed so far:
- Phase 6 (post-processing, 2026-07-28): `energies.csv:kinetic == e_kin_ha` to
  0.000 eV and differs from the extensive `e_kin_ha*norm` by 416 eV -> INQ prints the
  per-particle (norm-divided) value. Direct proof of the mechanism.
- The decisive Phase 3 (mask+ETRS vs mask+CN, matched geometry) is being run in
  inq-study; `WP_ABS`/`WP_PROP` switches added to
  `ResearchProject/systems/vacuum/scripts/wp_traversal_energy/run.cpp`.

### Cross-run results (2026-07-28, all 13 vacuum runs complete)

Aggregation `.../wp_traversal_energy/aggregate_investigation.py` ->
`results/investigation_summary.csv`. Key column is `E_ext = E_reported*norm`.

- **Phase 1 (eta-sweep, 5 runs).** ΔE_reported grows with |η| (+1.3 -> +14.0 eV) but
  ONLY because stronger CAPs absorb more (norm_T 0.44 -> 0.00). `E_ext/E0 == norm_T`
  to 3 d.p. in every run. The rise tracks the absorbed fraction, not η per se ->
  falsifier (b) (reflection) FAILS.
- **Phase 2 (partial-absorption, 3 runs).** `E_ext/E0` vs `norm` is the identity line
  y=x: 0.708/0.549/0.317 at norm 0.707/0.548/0.316. Cleanest confirmation.
- **Phase 3 (DECISIVE — the surprise, RESOLVED).** Both mask absorbers show ~+19 eV
  in `E_reported` (ETRS +18.86, CN +19.15). On the reported energy alone, falsifier
  (c) (a norm-preserving absorber still rises) appears to fire. BUT `E_ext` separates
  them: **ETRS -> E_ext/E0 = 0.000** (pure normalization artifact, energy correctly
  absorbed); **CN -> E_ext/E0 = 1.048** (extensive energy genuinely ROSE 4.8%).
  Crank-Nicolson does not merely preserve norm — it RENORMALIZES the sin²-clipped
  orbital every step (`crank_nicolson.hpp:139-165`), and clip-then-renormalize injects
  real high-k content. So **mask+CN is not a clean norm-preserving control**; it pumps
  real energy, a DIFFERENT mechanism from the CAP artifact. The physics CAP runs use
  ETRS + `perturbations::absorbing` (norm-losing), where `E_ext = E_reported*norm` is
  the correct extensive energy. **The normalization hypothesis STANDS; the red flag
  was a Crank-Nicolson renormalization confound, exposed by the `E_ext` diagnostic.**

Falsifier scorecard: (a) E_ext tracks norm — HELD (Phases 1,2); (b) residual scales
with η→reflection — FAILS (collapses on absorbed-fraction); (c) norm-preserving
absorber still rises — APPARENT-yes on E_reported but E_ext reveals it as energy
pumping, not the artifact; (d) rise without norm loss — only under CN's active
renormalization, itself a norm-changing operation. Phases 4 (numerics), 5 (spectral
width), 1b (W-sweep), 1c (two-sided) not yet run.

Notebooks: per-run deep dives at `results/<run>/report/run_report.ipynb`; phase
study notebooks + index at
`ResearchProject/systems/vacuum/hypotheses/cap_norm_investigation/`.

### IN-RUN FIX VALIDATED (2026-07-29, double-sided-CAP vacuum test)

`inqkit::observables::OrbitalKineticStats` (per-orbital BARE kinetic + norm,
one set-FFT per recorded step, physical units — INQ's to_fourier is an
unnormalized DFT, scale dV/N_grid) validated on `dcap_extkin` vs
`dcap_baseline` (LZ=60, CAP_L=15 BOTH ends, launch z=0, η=−3.5, 700 steps):

- **Identity EXACT:** max |Σocc·T_i/norm_i − energies.csv:kinetic| = 0.0 Ha at
  ALL 701 steps — we compute the same per-orbital T_i INQ reduces, minus /norm.
- **Artifact + fix:** E_reported(final) = 383 eV (pinned at the remnant mean,
  norm 3.5e-6) vs E_corrected = total − kinetic + kin_bare → 0.00 eV,
  == E0·norm throughout (captured 100.0% of the 402 eV).
- **Post-hoc route equivalent:** e_kin_ha·norm matches kin_bare to 2.5e-9 eV
  (single-orbital case).
- **Cost negligible:** observable self-timing 0.42 ms/step vs ~300 ms/step
  run cost (0.14%, 1 orbital, 844k grid); run-level ON−OFF wall Δ = −15 ms/step
  (i.e. within run-to-run noise, 30× the observable's own cost). Jellium-162
  estimate: ~one extra forward set-FFT per recorded step vs ETRS's several
  FFT passes per step → a few % at every-step cadence; measure in pilot,
  reducible via WP_EXTKIN_EVERY.

Artifacts: `.../hypotheses/cap_norm_investigation/extensive_kinetic/`
(fig_extkin_energies.png, fig_extkin_identity_timing.png, extkin_summary.txt).

## Related

Memory [[reference_inq_reports_normalized_energy]]; diagnostics
`ResearchProject/systems/vacuum/scripts/wp_traversal_energy/energy_diagnostics.py`;
jellium correction `.../wp_cap_energy_plateau/wp_kinetic_normalization_fix.py`;
stopping schemes `docs/notes/refined-stopping-two-schemes.md`.
