# Skill: Simulation Validation

Use this skill when proposing benchmark systems, defining convergence tests, or writing a validation plan for a new INQ calculation.

---

## Protocol

### 1. Propose the validation menu (before running anything)

Present the user with a tiered menu of tests. The user approves which ones to run.

**Tier A — Fast, always run (~seconds to minutes)**
- SCF convergence reached (energy change < tolerance, printed in log)
- No NaN or Inf in energies
- CPU vs GPU consistency (run same input, compare total energy to 6+ significant figures)
- Restart check (save GS, reload, run a few RT steps, compare to original)

**Tier B — Medium (~minutes)**
- Known small-molecule benchmark (H2, N2, CO — compare total energy to literature or previous run)
- Energy conservation during RT propagation (plot total energy vs time; drift < 0.1% over full run)
- Forces sum to zero for an isolated finite system

**Tier C — Expensive (~hours; user must explicitly approve)**
- Full convergence study: cutoff energy (20, 30, 40, 60 Ry — plot energy vs cutoff)
- k-point convergence (Gamma, 2×2×2, 4×4×4 — plot energy vs k density)
- Time-step convergence (dt = 0.08, 0.04, 0.02 atomictime — compare spectra)
- Comparison of ETRS vs Crank-Nicolson (same system, compare dipole signal)
- Comparison of propagator with published spectrum (e.g. N2 optical response)

### 2. Record validation status

In the handover or a `docs/validation/<system>.md` file:

```md
## Validation status: <system> <calculation type>

### Tier A
- [x] SCF converged to 1e-6 Ha after 23 iterations
- [x] CPU/GPU: energies agree to 8 significant figures
- [ ] Restart check: not yet run

### Tier B
- [x] H2 total energy: -1.1745 Ha (literature: -1.1744 Ha, PBE/40 Ry)
- [ ] Energy conservation: not yet checked

### Tier C
- [ ] Cutoff convergence: not yet run (approved: no)
- [ ] k-point convergence: not yet run (approved: no)

### Remaining gaps
<anything still unverified>
```

---

## INQ-specific validation patterns

### Ground state checklist
- Total energy: compare to literature or previous run to 5+ significant figures
- Band gap / HOMO-LUMO: qualitatively correct (metal vs insulator)
- Magnetic moment (for spin-polarised): matches expected value
- Forces: computed and compared with finite-difference forces (optional, expensive)

### TDDFT / real-time checklist
- Initial dipole at t=0: should be near zero (unless system is polar)
- Dipole oscillation: correct direction, reasonable amplitude
- Long-time behaviour: signal should decay / ring naturally (no exponential growth)
- Energy conservation: total energy drift < 0.1% over run duration
- Spectrum peak positions: compare to literature or TDDFT reference values
- Current conservation: ∇·J + ∂ρ/∂t ≈ 0 (optional, numerical check)

### Ionic dynamics checklist
- Total energy (electronic + ionic) conserved in Ehrenfest run
- Ion positions do not drift unphysically
- Impulsive run: initial velocity applied correctly, forces ignored as expected

---

## Reference systems

| System | Purpose | Reference energy |
|---|---|---|
| H2 molecule (finite box) | Minimal GS benchmark | ~-1.1744 Ha (PBE) |
| N2 molecule (finite box) | Dipole response benchmark | See INQ tutorial docs |
| Li BCC (4-atom cell, k-grid) | Metal GS + TDDFT kick | See Tutorial/li-bcc/ |
| Jellium slab | Free-electron benchmark | Analytic limits known |
