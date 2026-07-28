# Refined stopping power — the two measurement schemes (working notes)

Date: 2026-07-28. Status: both schemes formalised; Scheme 1 implemented and
validated on qsp_phase5 + wp_cap_energy_plateau; Scheme 2 partially implemented
(data already recorded; standalone S extraction not yet run).

## Problem these schemes solve

The plateau/deposit estimate `S = [E(0) − E_plateau-ish]/L` over-counts stopping
because a wavepacket, unlike a classical projectile, carries **localisation
(zero-point) energy** that is banked into the target ledger on absorption or
capture:

```
E_kin,wp = <p>²/2m  +  Var(p)/2m          ("split A", momentum space)
         = drift       localisation/spread
T_loc(0) = 3/(8·sigma_r²)     [sigma_r = sigma_WP/sqrt(2), per-axis density std]
         = 3.0 Ha = 82 eV  (sigma_WP = 0.5)   |   0.75 Ha = 20 eV  (sigma_WP = 1)
```

Measured over-count: 1.9–7.9x across qsp_phase5; plateau dissection of the σ=1
pair shows D ≈ T_loc + WP self-Hartree + ~5 eV genuine deposit.

**Refined stopping = the momentum-dependent (drift) KE lost by the projectile:**

```
S_drift = − d(<p>²/2m N) / ds        headline; excludes localisation by construction
```

Both schemes below measure this; they differ in what they call "the projectile".

---

## Scheme 1 — density-based, orbital-free (T_W and T_v)

**Projectile definition:** whatever excess density `n(r,t) − n_gs(r)` occupies
the vacuum corridors (regions where ownership is unambiguous — the bath is
exponentially localised to the slab). No Kohn–Sham label used anywhere.

**Formalism** (Madelung / field split, "split B"; exact where the density in the
region is one coherent lump):

```
<p²>/2m = T_W + T_v
T_W = ∫ |∇n|²/(8n) d³r          shape / localisation   — DENSITY ONLY, exact
T_v = ∫ |j|²/(2n)  d³r          flow (drift + spread)  — needs current j
P   = ∫ j d³r                   →  T_drift = |P|²/(2N),   N = ∫ n d³r
```

**The j problem — direct measurement is PREFERRED, reconstruction is the
fallback.** The right way to get j is to save it from the engine
(`observables::current_density`, see the dedicated section below): exact 3D
vector field, no assumptions. The continuity route below exists ONLY because
the existing runs saved no current frames; any new run must save j directly,
and any analysis on a run that has j frames must use them instead of this
reconstruction. For legacy runs: in 3D, continuity (∂n/∂t = −∇·j) fixes only
div j, not j. Fix: integrate transversally → 1D continuity fully determines
the longitudinal flux (boundary J=0 at the CAP edges):

```
∂ρ/∂t + ∂J/∂z = −2W(z)ρ    →    J(z,t) = −∫_{−L/2}^{z} (∂ρ/∂t + 2Wρ) dz′
ρ(z,t) = ∫∫ (n − n_gs) dx dy ;   W(z) = |η| sin²(π(|z|−z_cap)/w)  [CAP profile]
```

Consequences: P_z and the longitudinal flow surrogate `T_v,z = ∫ J²/(2ρ) dz`
(lower bound) are recoverable; transverse flow is NOT (free-dispersion analytic
estimate, labelled). Side-adaptive evaluation (entrance plane from left edge,
exit from right) is mandatory — long integration paths accumulate error.

**Two operational variants:**

1. **Flux / TOF picture** (robust, headline): fixed detector planes; passing
   velocity u(t) = J/ρ; exceedance curves N(>u) in/out; equal-rank
   (monotone-transport) matching, robust to capture removing the slow end:
   ```
   S(u_in(q)) = ½[u_in²(q) − u_out²(q)] / L
   ```
2. **Snapshot picture**: at fixed post-interaction t*, region integrals of the
   corridor density; quotes coverage N(t*) < 1; brackets the flux answer
   (ensemble-selection drift across t* is explicit, not hidden).

**Validation battery (all passed):** CAP-sink closure −dN/dt = A₊+A₋ (few %);
free-packet NULL test through the identical pipeline (S ≡ 0 within a trusted
rank window; systematic ±0.4 eV/Bohr for qsp5 geometry, ±1.6 for the
short-entrance σ=1 geometry); analytic-Gaussian-on-grid control for discrete
T_W (launch width 0.7·dx inflates T_W(0) 1.81x — machinery, not physics; from
sigma_z ≳ 1.5·dx cached T_W tracks the free-dispersion law to ≲3%).

**Assumption stack, ranked:** (1) transverse-flow estimate; (2) T_v,z lower
bound; (3) rank matching assumes velocity-ordering preserved; (4) bath
polarisation negligible beyond detector buffers; (5) coverage < 1 (snapshot).

**Implemented in:**
- `hypotheses/qsp_phase5/qsp_phase5_momentum_stopping.ipynb` (+ builder, cache extractor)
- `hypotheses/wp_cap_energy_plateau/wp_cap_energy_plateau_momentum_stopping.ipynb`
- `hypotheses/qsp_phase5/p5_wp_v1p3_snapshot_kinematics.ipynb`

---

## Scheme 2 — KS-orbital momentum moments (remove the localisation part explicitly)

**Projectile definition:** the wavepacket's KS orbital ψ_wp (the injected extra
state). Compute momentum moments directly from the wavefunction:

```
<p>  = ∫ ψ_wp* (−i∇) ψ_wp d³r          (per unit orbital norm)
<p²> = ∫ |∇ψ_wp|² d³r                   (= 2·E_kin,orb)
Var(p) = <p²> − <p>²
```

Split and remove the localisation part explicitly:

```
E_kin,orb = <p>²/2  +  Var(p)/2
            └ keep ┘    └ remove ┘  (localisation + spread — not stopping)
S_orb = [ <p>²/2 |_pre  −  <p>²/2 |_post ] / L
```

evaluated in the pre-contact and asymptotic post-interaction windows (the
orbital is trusted as "the projectile" outside the interaction zone; during the
collision it is a label — mid-collision moments are diagnostics only).

**Data — already recorded per step, nothing to re-run:**
- `wp_momentum_stats.csv`: px/py/pz means, p² per component, sigma_p² per
  component, e_kin_ha — i.e. the full split A per timestep.
- `momentum_distribution.csv`: |ψ̃_wp|² binned in |k| (velocity-resolved checks).
- `wavefunction_wp` VTIs: recompute any moment, incl. spatially windowed ones.

**Known failure modes (why Scheme 1 exists):**
- Fragmentation: after 50/50 reflection+transmission, <p> averages over lobes
  → drift underestimates transport. Fix: window the moments by region
  (compute from wavefunction_wp restricted to the transmitted lobe).
- CAP norm drain: moments are over the surviving |ψ_wp|²; quote the norm.
- Hybridisation/indistinguishability during contact: orbital ≠ particle.

**Status (updated 2026-07-28): standalone Scheme-2 extractions now exist as two
method notebooks on p5_wp_v1p3:**

- **Method #3 — Ehrenfest drag (local, time-resolved):**
  `hypotheses/qsp_phase5/p5_wp_v1p3_ehrenfest_drag.ipynb`.
  S = −d(⟨p_z⟩²/2)/dz_c from the per-step moments. Discovery: the electron
  packet is *accelerated into* the jellium (image attraction, pz 1.30→1.45), so
  the light-projectile early-window rule fails for attractive projectiles — the
  correct "S at v0" window is POST-PEAK (v ≥ 0.85·v_peak): S(v0) = 1.2 eV/Bohr
  at v̄ = 1.37. Full-deceleration slope 3.0 ± 0.6 at v̄ = 0.47 (sweep average,
  not a v0 number). CAP selection biases drag numbers low (lower bounds).
- **Method #4 — momentum-space distribution (asymptotic, spectral):**
  `hypotheses/qsp_phase5/p5_wp_v1p3_momentum_space.ipynb`.
  Signed P(k_z) marginal from the complex wavefunction frames; rank-matched
  against the corridor-windowed transmitted envelope: S = 0.38 ± 0.27 at
  u ≈ 1.7 — agrees with Scheme 1's TOF (0.30 ± 0.42), cross-validating the
  continuity reconstruction AND the asymptotic orbital-identity assumption.

Open question: local drag (≈1.2 at v̄ 1.37) sits above the asymptotic estimates
(≈0.3–0.4 at u ≈ 1.7) — candidates: elastic momentum return on exit, CAP
selection, genuine u-dependence. Unresolved; see handover.

---

## Relation between the schemes

```
t = 0 (min-uncertainty launch):  T_W = Var(p)/2m,   T_v = <p>²/2m   (splits coincide)
vacuum, any t:  split A constant per term;  split B converts T_W → T_v (dispersion)
asymptotically:  Scheme1 T_drift/N  ==  Scheme2 <p>²/2 (lobe-windowed)  — the test
```

Complementary failures: Scheme 1 is immune to orbital-identity questions but
needs reconstruction assumptions for j; Scheme 2 has exact j (from ψ) but needs
the orbital-identity assumption. Where both are valid they must agree.

---

## Direct j(r,t) from the simulation — the PREFERRED route (available, just never saved)

**Policy (user decision, 2026-07-28): measuring j directly from the simulation
is preferred over the continuity-equation reconstruction wherever possible.**
The reconstruction in Scheme 1 is a legacy-run fallback, not the method of
record. INQ computes the full vector current-density field natively:

```
inq/src/observables/current.hpp:
  observables::current_density(ions, electrons, ham)
      → basis::field<real_space, vector3<double, covariant>>   (j(r), all orbitals)
  observables::current(...) = integral of the above            (what runs already log)
inq/src/observables/kinetic_energy_density.hpp:
  observables::kinetic_energy_density(electrons)               (t(r) field)
```

**Future-run requirement:** in the `real_time::propagate` callback, wrap
`current_density` through an inqkit `RealField3DWriter` (one VTI per Cartesian
component, physical-order convention as for density) at the density cadence
(`LJ_SAVE_EVERY`); optionally `kinetic_energy_density` too. Per-orbital current
for the WP state is equally accessible from the orbital (j_wp = Im ψ*∇ψ).

What saved j frames buy, immediately:
- Scheme 1's T_v becomes exact 3D (no continuity reconstruction, no side-adaptive
  machinery, no transverse-flow estimate — assumption stack items 1–2 deleted).
- P from ∫j directly; detector planes read real flux; null-test systematic
  shrinks to time-sampling only.
- With kinetic_energy_density: T_full exact for the total system → the
  region-resolved energy ledger closes without any one-lump assumption.

Cost: 3 extra scalar fields at density cadence ≈ 3x density I/O — acceptable at
the 250–330-frame cadences used to date.

---

## Artefact index

| artefact | path |
|---|---|
| TOF/rank-matched sweep notebook | `ResearchProject/systems/localised_jellium/hypotheses/qsp_phase5/qsp_phase5_momentum_stopping.ipynb` |
| plateau dissection notebook (σ=1) | `.../hypotheses/wp_cap_energy_plateau/wp_cap_energy_plateau_momentum_stopping.ipynb` |
| snapshot kinematics notebook (v1p3) | `.../hypotheses/qsp_phase5/p5_wp_v1p3_snapshot_kinematics.ipynb` |
| headline numbers | `momentum_stopping_summary.json`, `plateau_dissection_summary.json`, `snapshot_kinematics_summary.json` (beside their notebooks) |
| rolling handover | `docs/handovers/qsp5-momentum-stopping.md` |
