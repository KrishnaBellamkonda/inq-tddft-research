# Why the electrostatic residual reads 7.4 eV, not the clean 21.5 eV self-Hartree

Decisive verification (2026-07-13), campaign localised-jellium-dynamics-analysis.
Eval binary: scripts/localised_jellium_dynamics/eval_projpot/run.cpp (single-point,
INQ p2 open-z convention, no dynamics). Densities: p5 at-rest (r=12, Lz=120).

## Exact decomposition (reproduced to +-0.05 eV vs INQ)
residual = self_Hartree  -  INT (n_slab - n_bg) . (v_ion - V_proj_ideal)
                            |__ slab spillout __|  |__ pseudopot error dv __|

| r_cut | self_Hartree | INT(spillout).dv | residual (recon) | INQ measured |
|-------|-------------:|-----------------:|-----------------:|-------------:|
| 50    | 20.83        | 6.82             | 14.0             | 13.99        |
| 120   | 20.83        | 13.45            | 7.4              | 7.36         |

## Three independent proofs it is a pseudopotential representation artifact (NOT a gauge)
1. r_cut dependence: error term doubles 6.8 -> 13.5 eV as the erf/r tail wraps further.
2. Grid pathology (dx sweep, r_cut=120 impl): -524.5 -> -268.7 -> +223.7 -> +34.3 eV
   for dx=0.5,0.4,0.3,0.25 -- swings SIGN, does not converge. r_cut=50 impl is grid-stable
   (~ -140 eV, +-3%). ideal term grid-stable (~135 eV).
3. Consistent-ideal recomputation (consistent_ideal_residual.py): residual = self-Hartree
   = 21.49 eV exactly (distortion term 0.000).

## Mechanism
UPF ghost potential v(r)=erf(r/0.5)/r ~ 1/r (long-range Coulomb tail; Z_valence=0 so INQ
places the whole tail as "local potential", no reciprocal long/short split). Tabulated to
r_cut and gridded, this conditionally-convergent tail ALIASES; r_cut=120 (>> Lx=Ly=50)
wraps 2.4x and is worse than r_cut=50. Bigger r_cut is WORSE, not better.

## Conclusion for the balance sheet
Physical r-independent residual = WP Hartree self-energy ~= 20.8 eV (INQ p2) / 21.5 free-space.
The 7.4 eV is self-Hartree minus the aliasing artifact -- not a physical quantity.
Charged-cell gauge is <1 eV (self-Hartree channel only): 21.71 free -> 21.49 periodic ->
20.83 p2-open-z.
