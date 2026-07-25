# Phase 4 — muon exchange–correlation: grounded candidates for the user's pick

**Status:** research complete, awaiting user decision (checkpoint).
**Decision owner:** user (2026-07-06: "Agent researches, then asks me to pick").
**To resume Phase 5:** write the chosen option to
`ResearchProject/systems/localised_jellium/scripts/muon_mass_fork/muon_xc_pick.json`
and re-run `orchestrate.py`.

Sources: [[heg-mass-scaling-xc]], [[kreibich-gross-multicomponent-dft]],
[[car-parrinello-fictitious-mass]].

---

## The physics in one paragraph

For a **one-component all-muon jellium**, the Hamiltonian `−(ℏ²/2m)∇² + Coulomb`
is *exactly* mass-scalable: measuring lengths in the muon Bohr radius
`a_μ = a₀/206.77` and energies in the muon Hartree `E_μ = 206.77·Ha`, the muon
problem is identical to the electron problem at the **same** Wigner–Seitz radius.
Consequences (grounded, [[heg-mass-scaling-xc]]):
- **Exchange is mass-independent** at fixed physical density — naive electron-LDA
  exchange is already exact for the muon.
- **Correlation carries the mass**: `ε_xc^μ(n) = m·ε_xc^HEG(r_s^μ)`, a function of
  `M·r_s`, with `r_s^μ = m·r_s^(a₀)`.
- At the LOCKED physical `r_s^(a₀)=5.69`, the muon effective coupling is
  `r_s^μ ≈ 1177` — deep in the strongly-correlated / near-Wigner regime where the
  electron-fit LDA (PZ81/PW92, fit to QMC for `r_s ≲ 100`) is a long extrapolation.
  This is exactly why the comparison is informative.

---

## Candidate prescriptions

### A — Mass-rescaled LDA  *(recommended default for the all-muon run)*
Evaluate stock LDA at the muon's own units: `ε_xc^μ(n) = m·ε_xc^LDA(r_s^μ)`,
`r_s^μ = m·r_s^(a₀)`. For a pure muon HEG this is **exact within LDA** — no new
functional, just a density-argument rescale + an energy prefactor `m` in a wrapper
around libxc. *Cost:* small (a wrapper in `inqkit`, or a rescaled effective density
fed to the existing LDA). *What it tests vs naive LDA:* the pure effect of
evaluating XC at the correct effective coupling `r_s^μ=1177` instead of `5.69`.

### B — Naive electron-LDA at physical r_s  *(the control / baseline)*
INQ's stock `options::theory{}.lda()` applied to the muon density as-is: XC read
at `r_s^(a₀)=5.69`, no mass awareness. Physically this **mis-locates** the muon's
effective coupling but is what a naive user would run. *Cost:* zero (already the
default). *Role:* the baseline the campaign compares against.

### C — Multicomponent / NEO electron–muon correlation
Coupled KS fields + an explicit electron–muon correlation functional
[[kreibich-gross-multicomponent-dft]]. *Only* needed if the muon is a **distinct
species** from the bath (muon-in-electron-jellium, handover Q3) — NOT the locked
all-muon target. *Cost:* high (second KS field + non-libxc cross functional).
*Role:* the rigorous route if the model is later made two-species.

### (Excluded) Car–Parrinello fictitious mass
Explicitly **not** a physical-XC option — a numerical MD parameter, unrelated to
the functional [[car-parrinello-fictitious-mass]]. Listed only to prevent a
category error.

---

## Recommendation for Phase 5

Run the LOCKED all-muon `r_s=5.69` bath **twice**: **B (naive electron-LDA)** as the
baseline vs **A (mass-rescaled LDA)** as the physically-correct one-component XC.
The difference is a clean, grounded measure of XC sensitivity dominated by
correlation (exchange being mass-invariant), with the SIE floor bounded by a
vacuum-WP control. **C is deferred** unless you want the muon treated as a species
distinct from the bath.

## Pick format (`muon_xc_pick.json`)
```json
{ "functional": "A_mass_rescaled_lda",
  "baseline":   "B_naive_electron_lda",
  "rationale":  "one-component all-muon HEG is exactly mass-scalable; compare correct r_s^mu vs naive r_s",
  "implement":  "inqkit LDA wrapper: feed r_s^mu = m*r_s, multiply e_xc by m" }
```
