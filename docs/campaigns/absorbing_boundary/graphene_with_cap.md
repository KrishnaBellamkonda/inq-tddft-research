---
id: ab-graphene-cap
area: absorbing_boundary
title: "Graphene/coronene with CAP - perpendicular + grazing scattering"
status: paused
hypothesis: "A CAP-bounded box captures projectile scattering off a finite flake for both perpendicular and grazing impact, giving stopping / charge-transfer observables."
handover: docs/handovers/graphene-cap.md
tasks:
  - { name: "Coronene GSes built+validated (perp x-y + grazing y-z, same 20x22x60 box/grid)", done: true }
  - { name: "Perpendicular classical head-on cl_b0", done: true }
  - { name: "Perpendicular WP head-on wp_b0 (build perp/wp binary + run)", done: false }
  - { name: "Grazing classical cl_b{1,3,6} (cl_b1 done; cl_b3, cl_b6 remain)", done: false }
  - { name: "Grazing WP wp_b{1,3,6}", done: false }
  - { name: "Comparison study notebook (with first-pass stopping fix)", done: false }
blocked_reason: "PAUSED 2026-06-21 21:55 by user; GPU 0 freed; 2/8 runs done; resume per <resume_state> below + handover docs/handovers/graphene-cap.md"
---

# Graphene with CAP

<resume_state>
⚠️ THIS CAMPAIGN IS MID-FLIGHT — PAUSED 2026-06-21 21:55, 2 of 8 runs done.
READ THIS BLOCK FIRST; it is AUTHORITATIVE and supersedes any conflicting numbers
in the original `<description>`/`<paper>` below (the design evolved). Full detail
+ rationale: handover `docs/handovers/graphene-cap.md` (top "⏸ PAUSED" section).

WHAT THIS CAMPAIGN IS NOW (locked): a perpendicular-vs-grazing **impact-parameter
comparison** on ONE common target — a finite **coronene C24H12** flake — to see
how impact geometry changes the WP-vs-classical interaction. Both arms share the
SAME box (20×22×60 Bohr), grid (50 Ha → ~65×74×192), sim-time (N=1319,
τ≈26.4 a.u.), CAP (two-sided sin², L=20, η=−0.5 Ha, |z|∈[20,30]), projectile
(E=100 eV, k0=2.711, σ=1.47, ETRS, dt=0.02). ONLY the flake orientation differs:
- PERP = coronene in x-y plane (⊥ beam), beam +z head-on. GS
  `shared_gs/gs_perp_coronene_50ha`, geom `shared/geometry/coronene_flake_perp.xyz`.
- GRAZING = coronene in y-z plane (∥ beam), beam +z at impact parameter b=x-offset.
  GS `shared_gs/gs_grazing_coronene_50ha`, geom `coronene_flake_grazing.xyz`.
Both GSes validated (closed-shell, gap ≈2.7 eV, E≈−150.77 Ha, 108 e). The CLASSICAL
projectile uses the He-symbol z_valence=−1 UPF (`electron_gaussian_sigma1p47_He.upf`)
+ `.extra_electrons(+1)` — the species-collision + Ewald fixes; DO NOT revert.
All under `ResearchProject/systems/graphene/`.

CHECKPOINT (run-level):
- DONE (results preserved): grazing `cl_b1`, perp `cl_b0`.
- REMAINING: grazing `cl_b3`, `cl_b6`, `wp_b1`, `wp_b3`, `wp_b6`; perp `wp_b0`.
- PRESERVED (no rebuild): both GSes; binaries `scripts/grazing/{cl,wp}/run`,
  `scripts/perp/cl/run`. MISSING: `scripts/perp/wp/run` (build on resume).
- INQ RT has NO mid-run checkpoint → a killed run restarts from step 0.

RESUME (pick GPU via `cudaMemGetInfo` probe `systems/vacuum/gpu_probe`, NOT
nvidia-smi which is NVML-dead; GPU 1 was ~4.6× faster than GPU 0 for this load):
1. Grazing (dispatcher is RESUMABLE — auto-skips runs with run_completed=true):
   `cd ResearchProject/systems/graphene/scripts/grazing && GPU=<g> setsid nohup bash dispatch.sh > dispatch_resume.out 2>&1 < /dev/null &`
   → reruns cl_b3, cl_b6, wp_b1, wp_b3, wp_b6, then auto-builds the notebook.
2. Perp WP (needs a fresh ~15 min build first): build+run `scripts/perp/wp/run.cpp`
   vs inq-study with GR_CX=0 GR_CY=0 GR_CAP=1 GR_E_EV=100, GR_OUTDIR=
   `…/graphene/perp/run_wp_b0/results` (exact command in the handover).

KNOWN ISSUES to carry (both post-processing, no GPU):
- Classical KE_loss is z-WRAP-contaminated: the CAP absorbs electrons not the ion,
  so over τ the ion wraps the 60-Bohr z-box ~1.2× → init-vs-final KE is spurious
  (perp cl_b0 shows −2.5 eV "gain"). Real stopping = KE(t) drop over the FIRST
  flake pass only (cross ~step 233, clear into vacuum ~step 400–500). Fix the
  classical extraction in `hypotheses/grazing/build_grazing_report.py`. WP
  unaffected (CAP absorbs the WP).
- GPU 0 ran ~4.6× slower than GPU 1 for the identical workload (likely throttled /
  slower card; NVML dead so unconfirmed). Prefer GPU 1.
- All CAP results PROVISIONAL until the inq-study engine regression (Task #7).
</resume_state>

<identity>
You are a scientific computing researcher working on first principles simulations. You have a good understanding of first principles domain and are excellent at writing scientific standard code. You adhere to the rules, principles, workflows established in this repository. 
</identity>

<description>
In this task, we are going to place a graphene target at the center of a box. The box is going to be a cubic box with the length of 50 bohr. The wave packet is initialised at a distance of 10 bohr from the target (along the z axis on one end). The wave packet propagates towards the other end. On the other end, there is a Complex Absorbing Potential (CAP) placed. Here, we are going to follow roughly the setup mentioned in the Yao & Schliefe paper. I want to confirm that the CAP is placed between 40 Bohr and -40 Bohr, looping aroung the wrapped boundary. Does this give a total of 20 Bohr of length? If so, we do the same thing in our simulation. Have a CAP that loops around the wrapped boundary conditions. Then, before launching the simulation, we need to again find the value of W. To do this, we do not have the target system. Instead, we initialise a wave packet at the location, and look for the value of W that ensures that the least amount of electrons are reflected. This is done by considering a range of W values, running the simulation, and finding the best W to use. Then, we use this value of W for the graphene case too. Then, we run the simulation. 

A main objective of this run to make these runs. 

Standard trajectory: 
1. The run where the wave packet is initialsed at 10 bohr from target along the z direction (at x, y midpoint) with its momentum moving towards the target graphene system. Let's consider the energy of 100 eV. The gaussian width sigma would be 1.1 bohr (the same as the yao & schliefe paper).  
2. We then consider an ensemble of classical particles with a psedopotential with the radial potential with a gaussian width of 0.5 bohr of gaussian width. Now, the key part is produce the ensemble the way it was produced in the paper. From the initial position, we simulate these points and a taken an average of them as the answer. Let's consider 5 randomly chosen configs and runs. 

Different trajectory: 
In this trajectory, we are going to change the start position of the projectiles. Now, in this case, we are going to move the wave packet in the x direction, close to the graphene target. There is a set impact parameter. We will consider an impact parameter, or the distance away from the xy plane to be 0.5 bohr. Initially, the projectile starts far enough away from the graphene system. Then, it move towards it in a grazing manner.
1. Make a wave packet run with the same parameters as before. However, we need to ensure that the CAp boundaries are now placed along the opposite end that the projectile is heading to. 
2.  Same CAP as the previous case. Again, have an ensemble of classical particles to run this. 

We are going going to do this in steps. 
1. Make the graphene configuration. Confirm it converges. Validate by comparing with literature. 
2. Make the setup, send me plots showing me the setup for both trajectories with visualisation of the start position, graphene. Now run without the aborbing boundary conditions, both the wave packet traectories to have a baseline to compare against. 
3. Add the absorbing boundary conditions. Remake the plot and send it to me to have a look. I need to ensure that the setups are valid before you run them.  
4. Run standard trajectory. Start with the classical particles. Ensure that everything is as expected. Use the relevant hooks and observables from the run to do this. After each run, run post  processing. 
5. Run using the wavepacket.  Run post processing. 
6. Send the user the relevant plots. You can use the plots made in report1/draft5/ where the LEED types of plots can be useful for the user.  
7. Redo the previous steps with the different trajectory. 

You are going to make a ipynb file documenting the resutls of each step, with the required plots, equations etc that are required. Ensure that in the file, the setup and key quantities and their values are clearly dicussed. 

Email me the results at the end of each steps.  
</description>


<observables_set>
Have all the minimal observables stipulated by the claude ecosystem and the deterministic code. Then, have a bunch of LEED screens for all the runs at equal intervals away from the target system either way, along the xy direction, perpendicular to the direction of projectile propagation. Ensure that the 
</observables_set>

<guard_rails>
1. Ensure that the wave packet is orthonormalised, and hence must be initialised far enough from the target system. 
2. The wave packet initialisatio must be far away from the absorbing boundary conditions.   
</guard_rails>


<ensemble>
We use an ensemble of classical trajectories to approximate the quantum mechanical behavior. As the Gaussian wave
packet remains the Gaussian distribution in both position and momentum, we use std::mt19937 in C++ to sample
Gaussian random numbers in both position and momentum space, generated from the same mean and standard
deviation as the corresponding Gaussian wave packet. The input files are included in the Yifan: mdf . More than
100 classical trajectories in the ensemble at each kinetic energy.S4Figure S3. Remaining electrons in the vacuum after the wave packet reaches the complex absorbing potenti
</ensemble>

<paper>
Yao & Schliefe paper Complex Absorbing Potential setup: 

Figure 1. A schematic of the simulation cell and the planar
integrated charge density difference after the electron projectile impacts the graphene (dark green circles represent C
atoms) at 0.28 fs after the simulation starts. The initial wave
packet is located at ⟨z⟩ = −20 a0 and approaches the graphene
perpendicularly. The solid green line denotes the region of the
complex absorbing potential, Eq. S2. Two traces of the wave
packet are visible, before and after encountering the graphene
layer, since this plot shows the electron density difference.

As the wave packet traverses the host materials, electrons are emitted from both sides of the host materials and
propagate outwards. We employ a complex absorbing potential to avoid the non-physical re-approach of the wave
packet and the emitted electrons to the graphene due to the periodic boundary condition. It has a form of [67, 68]
VCAP(z) = −i · W · sin2 (z − zs) · π2 · dz, (S2)
for zs < z < zs + 2dz, where zs is the front boundary of CAP, and dz is the half width of the CAP. In our simulation,
we set zs = 0.5, fix the width, dz = 0.1 of the absorbing potential, and tune the height, W.
To find out the most suitable W for each incident wave packet, we compare the remaining electrons of the wave
packet in a vacuum, without graphene, after impacting the absorbing potential. We select the W that gives the
smallest remaining electrons in the vacuum. In the actual simulation of irradiated graphene, the wave packet will be
distorted after interacting with the graphene, which may alter the optimal W. However, the test here still provides a
reasonable estimation of the value of W.
</paper>


<rules>

</rules>