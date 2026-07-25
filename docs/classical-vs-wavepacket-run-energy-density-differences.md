The purpose of these runs made was to compare the wavepacket and the classical simulations runs. In these comparisons, hitherto, we had the ability to analyse difference in densities and induced densities. However, we add a power new tool to this using pairwise energy decomposition and total energy comparison. Meaning, we are able to track physical energy metrics and make analysis of the physics effects that arise. Cheifly, we'd be looking at the differences in the wavepacket and the classical cases. This helps us understand what is that wavepacket run has that classical doesn't and helps us put a physical lens on it. 

For reference, a lot of the work for the energy decomposition was done - docs/notes/energy-decomposition-skill.md
and the work for comparing the perutbation and the wavepacket runs is done - docs/notes/gaussian-pertubation-for-classical-simulation.md


### Description of all the runs made
For each configuration, a classical and an analogous wavepacket was run. The runs made are - 

No.     sigma       k0        Notes
1.      0.5         1.0       This is the baseline run
1.      1.0         1.1       In this case, the width is higher, with same (similar) to baseline velocity
1.      2.1         1.1       Even higher sigma than the second run
1.      2.0         4.2       wider and faster run. Meaning, we expect least spreading. 
1.      2.0         0.4       Slow and wide.
1.      2.0         0.5       A small change from the previous run. 



### Ground State Energy check

Ground-state energy accounting — all 5 pairs

┌───────────────┬───────────┬───────────┬───────────┬───────────┬───────────┐
│     Check     │  p5_null  │ p1_reflec │ p4_captur │ p2_tunnel │ p6_ladder │
│               │           │     t     │     e     │           │           │
├───────────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
│ Pairwise      │           │           │           │           │           │
│ slab/bg Δ≈0   │ 0.0 ✅    │ 0.0 ✅    │ 0.0 ✅    │ 0.0 ✅    │ 0.0 ✅    │
│ (no gauge)    │           │           │           │           │           │
├───────────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
│ Pairwise proj │           │           │           │           │           │
│ ectile-term Δ │ 0.000 ✅  │ 0.12 ✅   │ 0.002 ✅  │ 0.70 ⚠    │ 0.009 ✅  │
│  (t=0 density │           │           │           │           │           │
│  match)       │           │           │           │           │           │
├───────────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
│ ΔE_total =    │           │           │           │           │           │
│ ΣΔcomponents  │ 1e-13 ✅  │ 1e-13 ✅  │ 2e-13 ✅  │ 3e-13 ✅  │ 2e-13 ✅  │
│ (closure)     │           │           │           │           │           │
├───────────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
│ ΔE_kin =      │ 9.007/9.0 │ 0.269/0.2 │ 0.793/0.7 │ 0.342/0.3 │ 1.355/1.3 │
│ ½k₀²+3/4σ²    │ 08 ✅     │ 68 ✅     │ 93 ✅     │ 12 ⚠      │ 55 ✅     │
├───────────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
│ E_H,E_ext     │           │           │           │           │           │
│ from pairwise │ 1e-10 ✅  │ 1e-10 ✅  │ 1e-10 ✅  │ 1e-10 ✅  │ 1e-10 ✅  │
│  (Poisson)    │           │           │           │           │           │
├───────────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
│ Residual R    │           │           │           │           │           │
│ (WP           │ 4.37 eV   │ 4.34 eV   │ 4.37 eV   │ 3.89 eV   │ 9.79 eV   │
│ self-Hartree) │           │           │           │           │           │
├───────────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
│ SIE = R +     │ −0.25 eV  │ −0.32 eV  │ −0.30 eV  │ −1.39 eV  │ +1.13 eV  │
│ ΔE_xc         │           │           │           │           │           │


So, this confirms that, more or less, every run has almost all of the energy accounted for in the ground state. The only problem occurs in the run where the wavepacket was initialised in the center of the slab. 


### Notes
1. Note: It is important to check that for each pair, the energies agree with each other. Meaning, there is no unexplained energy at the first timestep. This ensures that we are able to track each of the run comfortably. 
