---
id: cap-jellium-baselines
area: cap_in_jellium
title: "Jellium CAP baselines (B0-B3)"
status: done
hypothesis: "A two-sided sin^2 CAP in interacting jellium yields well-characterised drainage (B1) and projectile (B2/B3) baselines for later S(v) work."
handover: docs/handovers/cap-in-jellium-baselines.md
tasks:
  - { name: "new observables code (J-field VTI, flux reducer, region N, per-orbital E, n(k))", done: true }
  - { name: "B0 GS reference", done: true }
  - { name: "B1 drainage", done: true }
  - { name: "B2 classical E100", done: true }
  - { name: "B3 WP E100", done: true }
blocked_reason: ""
---

# Baseline runs of Jellium with CAPs

<identity>
You are a scientific computing researcher working on first principles simulations. You have a good understanding of first principles domain and are excellent at writing scientific standard code. You adhere to the rules, principles, workflows established in this repository. 
</identity>

<description>
In this task, we are going to use complex absorbing potentials in jellium and establish credible baselines. These baselines are essential as in the next tasks, we are going sweep the different energies of projectiles moving in jellium (classical and wavepacket) . However, up until now, we did not have complex absorbing potentials. So, I need to carefully understand the baslines to be able to judge the later results effectively. The idea of the simulation is as given below - 
1. Complex absorbing boundaries (sin^2) implemented in a twosided manner. I have chosen the width of 20 bohr at eta = -0.5 as the optimal configuration for the CAP
2. The box is cubic jellium with 50 bohr side, and N = 162 electrons (the same box we've have been using for a bit)

In all the follwing runs, we are going to use an extensive set of observables. Firstly, we are going to max out, meaning use all the observables that can be observed in this simulation (ALL that are possible for each baseline run). Also, we are going to have a bunch of screens at evenly spaced intervals from the start of the box to the end. These are going to be LEED screens for now. Especially, I want these screens at the edge of the CAPs, meaning, on the boundary of the free region and the CAP. This is helpful because it tells me the density being absorbed by the CAP. 

Baseline 0 (added in grill 2026-06-17):
- The converged Hermitian ground state itself, reported as the t=0 equilibrium
  reference (no propagation, no CAP). Reuse the validated checkpoint
  `gs_L50_cubic_N162_dx0p40`. This is cheap and gives the reference density /
  energy components / eigenvalues / occupations that Baselines 1-3 are measured
  against.

Baseline 1: 
- We run the simulation without any classical or quantum particle for a long duration (consider 5 * time taken for a classical particle with 100 eV energy, with its mean momentum going from the initial position to the end of the box). 

Baseline 2: 
- Now, we use a classical particle with radial potential modelled as a gaussian with width of 0.5 bohr and run it for the same time as the previous simulation. 

Baseline 3: 
- Now, we use a wavepacket with a gaussian width of 0.5 instead of the classical particle. Again same timeperiod. 

Email me the results at the end of each steps.

I want all the results to be presented with context as described in the skill in a ipynb notebook. 
</description>


<observables_set>
Have the maximal set of observables. Then, have a bunch of LEED screens for all the runs at equal intervals away from the target system either way, along the xy direction, perpendicular to the direction of projectile propagation.  
</observables_set>

<observables_resolved>
<!-- Grill 2026-06-17. Maximal set = ADR-0006 universal core + per-run-type
     (classical: projectile track; WP: momentum/real-space stats + |psi(k)|^2)
     + the two new observables below + plane screens. -->

<new_observable id="efield" kind="derived-postprocessing">
Electric-field (E-field) kernel, NEW in inqview.analysis (pure numpy, deps-clean).
FFT Poisson solver, periodic BC: rho=-n, FFT, phi(G)=rho(G)/|G|^2 with G=0->0,
E(G)=-iG phi(G), iFFT. DECISIONS (locked):
  - Native ATOMIC units (4pi eps0 = e = 1); optional units="SI" conversion flag.
  - Produce BOTH fields: electron-only (induced/screening/wake field) AND
    electron+projectile-Gaussian (total field incl. bare driver). Baseline 1 is
    electron-only (no projectile).
  - Note: setting G=0->0 imposes the neutralizing jellium background, so the
    kernel returns the field of the fluctuation dn = n - n_mean. In a.u. this
    equals -grad(v_Hartree), which INQ already produces -> free cross-check.
  - Formula-bearing -> code-test loop + formula-validation agent. Locked
    known-cases: (i) uniform n -> E=0; (ii) sinusoid dn=cos(G0 z) -> E ~ sin(G0 z);
    (iii) Gaussian charge sigma<<L -> analytic erf field.
  - Visualisation (quiver / streamlines / |E| heatmap) -> inqview.visualisation.
</new_observable>

<new_observable id="current_density_field" kind="primary">
Current-density field J(r,t), NEW primary observable. DO NOT reconstruct from
orbitals in post: INQ already has observables::current_density
(inq/src/observables/current.hpp:26), GPU-computed from all KS orbitals WITH the
nonlocal-pseudopotential position-commutator term [r,V_nl] that the bare
J=(hbar/m)Im(psi* grad psi) omits (the projectile is a nonlocal pseudo-ion, so
the bare form violates continuity). Plan: new inqkit VTI writer for the covariant
vector field + one wiring line in run.cpp; dump at the SAME reduced cadence as
density. Perturbation current = J_full - J_baseline1 (post-processing).
</new_observable>

<continuity_diagnostic kind="derived">
"How much density flows away each step" done rigorously: with a CAP the
continuity eq carries a sink, d_t n + div J = -2 eta sin^2(.) n. For the free
region: dN_free/dt (from density) vs -boundary flux ∮J·dA (from J at the
CAP-edge screens). The residual = the CAP absorption inside the slab. This
isolates "CAP-removed" from "flowed across the boundary" and cross-checks
density, current, and screens against each other.
</continuity_diagnostic>
</observables_resolved>

<guard_rails>
</guard_rails>

<resolved_decisions>
<!-- Crystallised during the grill-with-docs session 2026-06-17. -->

<gs_cap_independence>
The CAP does NOT affect the ground state. Verified in inq-study source:
- `ground_state::calculator` defaults to `Perturbation = perturbations::none`
  (calculator.hpp:53); the CAP is only an argument to `real_time::propagate`
  (propagate.hpp:32, pert.potential -> vscalar at self_consistency.hpp:195).
- The complexify types `vscalar` as `HamiltonianType::potential_type`: `double`
  in GS (bit-identical to existing N=162 GS), complex only in RT
  (self_consistency.hpp:182).
Therefore Baseline 0 reuses the existing validated GS checkpoint
`gs_L50_cubic_N162_dx0p40`; no GS is recomputed.
Consequence to remember: the CAP switches on ABRUPTLY at t=0 (step in time, no
adiabatic ramp). The sudden turn-on launches the bath-drainage transient that
Baseline 1 characterises -- this is expected, not a defect.
</gs_cap_independence>

<drainage_goal>
The aim of Baseline 1 is NOT to show negligible bath perturbation. It is to
CHARACTERISE how much the CAP drains the equilibrium bath. Required diagnostics:
- bath depletion: electron number N(t) in the free (non-CAP) region + density
  profile n(z,t);
- probability current j(z,t), especially the inward current at the CAP edge;
- energy change: total energy E(t) AND per-orbital energies/eigenvalues over time;
- everything else the extensive observable suite yields.
Baseline 1 then serves as the subtraction reference for Baselines 2-3
(wake/projectile signal = full run - Baseline 1).
</drainage_goal>

<geometry kind="locked">
Box: cubic L=50 Bohr, z in [-25,+25], N=162, r_s=5.69 (unchanged; CAP carves
slabs out of the SAME box so GS density/r_s stay valid -> reuse
gs_L50_cubic_N162_dx0p40).
CAP: 20 Bohr TOTAL, two-sided, 10 Bohr per side. Slabs [-25,-15] u [+15,+25];
free region [-15,+15] = 30 Bohr. eta = -0.5 Ha, sin^2 (inq-study built-in
perturbations::absorbing, composed two-sided via perturbations::sum).
Launch (Baselines 2-3): z0 ~ -13 (4 sigma = 2 Bohr inside the -z CAP edge),
moving +z; projectile EXITS through the far (+z) CAP. No "stop before wall" path
cap any more -> run for the fixed duration clock. The 4 sigma/1 sigma boundary
rule now references the CAP EDGES (+/-15), not the box faces (+/-25).
NOTE for the future sweep (not these baselines): 30-Bohr free traversal is
shorter than the old ~50-Bohr runs -> shorter clean-force window for S(v).
</geometry>

<propagator kind="locked">
Electronic propagator = ETRS (INQ default, real_time.hpp:183). NEVER
crank_nicolson for CAP runs (CN renormalises -> undoes absorption). .ehrenfest()
is the ionic propagator for the classical projectile, kept.
</propagator>

<duration_and_energy kind="locked">
Baselines 2 & 3 run at a single representative energy E = 100 eV (matches the
E100 configs + reused GS; the energy sweep is a later task).
Clock: 100 eV electron, v0 = 2.711 a.u., from launch z0=-13 to far face +25
(38 Bohr) at constant v0 => t_cross = 14.0 a.u. Duration = 10 x t_cross = 140 a.u.
At dt=0.020 a.u. => ~7000 steps. ALL THREE baselines (1,2,3) use this SAME 140 a.u.
window (Option 1: common window so Baseline 1 is the exact subtraction reference).
10x (vs the original 5x) chosen so the bath drainage (~44 a.u. Fermi transit of
the 15-Bohr free half-width) gets ~3 transits and can approach saturation.
dt kept at 0.020 (CAP makes results dt-sensitive; matches existing E100 runs;
dt*E_max=0.62, dt*v=0.054 both comfortable).
</duration_and_energy>

<pilot_and_io kind="locked">
PILOT FIRST: ~100-step Baseline-1 pilot before the long runs, to gate (i) real
s/step on the current GPU, (ii) CAP behaves in INTERACTING jellium RT (energy
stays real, total energy + norm decrease smoothly = absorption signature, no
NaN), (iii) GS loads + slabs absorb. This is the first-ever CAP run in
interacting jellium (all prior CAP validation was non-interacting vacuum free-WP)
and Task #7 (inq-study engine regression) is still open -> ALL absorption/eps
numbers PROVISIONAL until Task #7.
VTI cadence: 300 frames target -> WRITE_EVERY ~ 23 at 7000 steps, for BOTH field
series (density system/wp/total AND the new current-density field). Scalar
observables every step. GPU scheduling via cudaMemGetInfo probe (NVML/nvidia-smi
broken); warn if a GPU is occupied by another user.
</pilot_and_io>

<screens kind="locked">
Term "LEED screen" RETIRED for jellium (LEED = coronene diffraction; wrong
primitive). Use PLANE / FLUX SCREENS: z-normal (xy) planes, each emitting
(a) planar density rho(x,y;z,t) via inqkit::screens::PlaneScreen and
(b) integrated z-flux ∮J_z dA via a NEW flux reducer over
observables::current_density. (b) feeds the continuity / CAP-sink check.
Layout (9 planes, 5-Bohr spacing): z in {-20,-15,-10,-5,0,+5,+10,+15,+20}.
 -15/+15 = CAP edges (key absorption monitors); -20/+20 = inside the slabs
(watch decay through the absorber); interior = free-region profile.
TimeAveragedScreen <rho> also available.
</screens>

<file_placement kind="locked">
Grouped-by-sweep (ADR-0007 amendment), even though jellium is otherwise flat:
- scripts/cap_baselines/      : ONE env-driven run.cpp (modes b1/b2/b3), built
                                 ONCE against inq-study (CAP engine); dispatch.py
                                 (2-GPU, emails per baseline); analyse.py.
- cap_baselines/run_b1_drainage|run_b2_classical_E100|run_b3_wp_E100/ : outputs.
                                 (Baseline 0 = reused GS checkpoint, no run dir.)
- hypotheses/cap_baselines/   : study notebook (auto-built by dispatcher tail),
                                 combined CSVs, tests/ (task-specific checks).
New code: current-density VTI writer -> inqkit/io/; flux reducer (∮J_z) ->
inqkit/screens/; E-field kernel -> inqview/analysis/ + viz in
inqview/visualisation/. Library-generic tests -> inq-stack/tests/.
</file_placement>

<observable_enumeration kind="locked">
"Max out" per baseline (B0=GS ref, B1=drain, B2=classical, B3=WP):
- Energy components (tot/kin/Hartree/xc/ext): all, per step (B0 at t=0).
- Total current ∫J, dipole, density_l2: B1/B2/B3 per step.
- Density field VTI @300: system (all); +wp,+total for B3.
- Current-density field VTI @300 (NEW writer): B1/B2/B3.
- 9 plane/flux screens (rho + ∮J_z) @frame: B1/B2/B3 (B0 rho only).
- Region N(t) free + per-slab (NEW reducer, per step): all.
- Per-orbital energy <psi_i|H|psi_i>(t) (NEW, @300 cadence): all (B0 = eps_i).
- Eigenvalues + occupations: B0 full; B1/B2/B3 final.
- TOTAL-SYSTEM momentum distribution n(k)=sum_i f_i |FFT(psi_i)|^2 (NEW;
  generalises the single-orbital momentum_distribution to ALL occupied orbitals):
  for B1, B2 AND B3. Probes background e-e scattering (peak broadening/shift of
  the filled jellium shells) + CAP absorption (total ∫n(k) drop). Store REDUCED
  per frame (radial n(|k|) + axial n(k_z)); FULL 3D n(k) only at t=0 and t=end
  (the "before/after" comparison the user wants).
- Projectile track: B2 = ions (z,v,F) per step (F = stopping force). B3 = COD
  track via inqview center_of_density on the wp density (derived/post, position +
  d/dt = mean velocity); no in-run code.
- E-field (NEW post kernel): electron-only for all; ALSO electron+projectile for
  B2/B3.
NOTE: each NEW primary observable (current-density VTI, flux reducer, region N(t),
per-orbital energy, total-system n(k)) is substantive C++ -> code-test loop +
catalogue row before the long runs.
</observable_enumeration>
</resolved_decisions>



<tasks>
</tasks>



<rules>

</rules>
