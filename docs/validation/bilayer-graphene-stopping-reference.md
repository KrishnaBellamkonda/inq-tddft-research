# Bilayer graphene — frozen experimental reference pack (BLINDED)

**Created:** 2026-08-05. **Status:** verified extractable, FROZEN.

## Blinding protocol (user decision, 2026-08-05)

The reference data below was acquired and frozen BEFORE any run was
designed or executed, so it cannot influence the experiment. **Do not
consult, plot, or quote the pack contents during Phases 0–3.** It is
unboxed ONLY in Phase 4 (comparison plots) of
`docs/plans/real-material-stopping-comparison.md`.

Pack location:
`/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/validation/reference_data_bilayer_graphene/`

## What was verified (the extractability check)

1. **ESTAR (high-energy anchor):** NIST ESTAR queried directly
   (`e_table-t.pl`, matno 906 = CARBON GRAPHITE, ρ = 1.7 g/cm³,
   I = 78 eV). Full table 1 keV–1 GeV parsed to
   `estar_graphite_electrons.csv` (91 rows, collision/radiative/total
   stopping + CSDA range, plus derived eV/Å column). Sanity:
   S_col(1 keV) = 1.8 eV/Å. **Caveat:** ESTAR floor is 1 keV — our
   15–300 eV window is BELOW it; ESTAR is the ≥1 keV anchor only.
2. **Geelen 2019 (low-energy anchor, overlaps our grid):** quantitative
   values quoted verbatim from PRL 123, 086802 →
   `geelen2019_anchor_values.csv`: λ_inel ≈ 3 layers near 0 eV → ≈ 1
   layer at 25 eV; λ_el up to ≈ 80 layers (2LG, 0–5 eV); 2LG+ shows a
   5–15 eV reflectivity max ABSENT in 1LG (interlayer bandgap);
   n−1 reflection minima. Extraction formulae (theirs, and computable
   from our runs): T = e^(−d/λ_tot), T + R = e^(−d/λ_inel),
   1/λ_tot = 1/λ_el + 1/λ_inel.
3. **ELF bridge route (15 eV–1 keV):** CXRO atomic scattering factors
   for carbon fetched (`cxro_carbon_f1f2.nff`, 10 eV–30 keV). **Caveat:**
   independent-atom approximation; f1 undefined below ~30 eV — NOT a
   substitute for the measured graphite optical ELF (Palik / EELS) at
   our energies. Proper low-E ELF-derived S remains a Phase 4 TODO
   (options: Palik graphite optical data; measured EELS loss function;
   or our own linear-response run). Labelled uncertain.
4. **Proof-of-extraction figure:** `reference_stopping_figure.png`
   (ESTAR curve + our-window shading; Geelen anchors + overlap band),
   built by `make_reference_pack.py` (reproducible; re-run to rebuild).

## Comparison observables this pack supports (Phase 4)

| Our run observable | Reference counterpart |
|---|---|
| transmission/attenuation per layer (WP norm) | Geelen λ_inel, T(E₀), R(E₀) |
| mono vs bilayer difference (5–15 eV feature) | Geelen interlayer-bandgap signature |
| S(E) trend at grid top → keV extrapolation | ESTAR S_col graphite ≥ 1 keV |
| S(E) absolute at 15–300 eV | ELF-derived S (route 3, TODO) |

## Provenance

- NIST ESTAR: https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html
  (ICRU 37 methodology), fetched 2026-08-05.
- Geelen et al., PRL 123, 086802 (2019); source note:
  `docs/sources/geelen-2019-evtem-graphene.md`.
- CXRO: https://henke.lbl.gov/optical_constants/ (Henke tables).
