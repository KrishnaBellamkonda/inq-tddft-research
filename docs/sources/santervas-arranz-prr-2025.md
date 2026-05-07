# Santervás-Arranz, Stengel, Artacho (PRR 7, 033292, 2025)

**Citation.** Santervás-Arranz, A.; Stengel, M.; Artacho, E.
*Excess energy in atomic kicks of solids: countercurrent and the Mv² rule.*
Phys. Rev. Research 7, 033292 (2025).

## Methodology used by QBall analyse.py

This paper is the methodological reference for the QBall Li/Al/C analysis in
`QuantumKickExtension/qball-codebase/`. The key claim:
when a uniform velocity is impulsively applied to all ions of a periodic
crystal, the long-time excess electronic energy per unit cell satisfies

    ΔE_plateau / N_uc  =  α (v) × (M_uc v²),
    M_uc = N_e_per_uc × m_e

where `α(v) ≤ 1` quantifies the **countercurrent** effect (electrons moving
opposite to ions reduces the energy below the rigid Mv² estimate). For Li
ONCV (3 valence e/atom × 2 atoms/uc → N_e_per_uc = 6) this reference line is
hard-coded in
`qball-codebase/Li/td_kicks/analyse.py:106-108` and used as the "Mv² reference"
in every plot.

## Diagnostic protocol replicated here

The post-processing in `analyse.py:245-271` is what the new INQ runs must
reproduce verbatim:

1. ΔE(t) = (E_total − E_GS) × Ha→eV / N_uc.
2. Detrend by subtracting the second-half mean (the "plateau").
3. Hann window.
4. 8× zero-padded `np.fft.rfft`.
5. Power spectrum `|FFT|²`, normalised within the 0–20 eV range.
6. Plot frequency on a hbar omega axis (eV).

For multiple velocities, the paper additionally tests `α(v) ≈ const` in the
linear regime — the "linear-response" check we replicate after the Phase 6
velocity sweep.

## What we are testing

The QBall configuration that produced the published Li plots used:
- Γ-only k-point sampling on a 54-atom supercell;
- 1000 K Fermi smearing;
- ~3.87 fs propagation (4000 steps × 0.04 a.u.).

The INQ replication asks whether the FFT peaks in the published plots survive:
- 2×2×2 shifted MP k-grid (≈ 6×6×6 primitive);
- 200 K Fermi smearing;
- ~15 fs propagation (better frequency resolution).

If the peaks shift, broaden, or vanish, that is evidence the published shapes
are partly artefacts of the QBall configuration rather than physics. Either
outcome is reportable.

## Attribution rule

When citing this paper in plots, comments, or report copy, prefer the short
form `(Santervás-Arranz et al., PRR 7, 033292 (2025))` to keep margins clean.
The paper's Figure 1 (Mv² ratio vs v) is the figure the QBall sweep targets;
do not reproduce it without acknowledgement.
