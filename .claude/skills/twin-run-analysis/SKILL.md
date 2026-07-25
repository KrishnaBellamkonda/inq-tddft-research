---
name: twin-run-analysis
description: Analyse a matched classical(perturbation)+wavepacket TWIN PAIR of localised-jellium runs — compute the per-timestep energy decomposition, the residual d(E_H+E_ext)-U_proj_bg (WP self-Hartree) and the one-electron SIE, then narrate the physical (quantum) differences. Deterministic engine (twin_decompose.py) does the arithmetic; the agent writes the narrative from the interpretation rules here. Consumes the twin-run-generation contract.
---

# Twin-Run Analysis

Given a **twin pair** — one classical (Gaussian-charge perturbation) run and one
wavepacket (WP) run, identical in every physical parameter except the projectile
representation — this skill quantifies and *physically interprets* the energy
differences between them, timestep by timestep. The differences are the **quantum
effects** of treating the projectile as a wavepacket rather than a classical point
charge.

The division of labour is fixed:

- **`twin_decompose.py` (skill-local, deterministic)** — does everything with a
  known formula: parses both runs, asserts config parity, computes per-term
  `d(·)=WP-classical` at every step, the residual, the SIE, and the known
  attributions; emits a structured **findings table**. Never guesses physics.
- **The agent (this SKILL.md)** — reads the findings table and writes the
  **narrative**: which terms moved, why, and what it means. The interpretation
  rules below are the agent's knowledge base.

## When to use

- A twin pair has been produced (by `twin-run-generation`) and validated
  (`check_twin.py` → `twin_manifest.json: valid=true`).
- The user wants the energy book-keeping / quantum-effect comparison of a
  classical-vs-WP pair, at rest or dynamic.
- NOT for a single run (use `run-notebook`), and NOT for a run-SET narrative
  (use `notebook-making`).

## The two-rung ladder

- **Rung 1 — static / known-answer (validated).** At-rest pairs (k0=0). The
  decomposition is a constant; the engine must reproduce the documented golden
  numbers (below). This is the regression test that proves the engine is correct.
- **Rung 2 — dynamic.** Moving projectile: `U_proj_bg` becomes per-step (needs the
  `energy_proj_bg_ideal` column, a generation-side addition), the residual evolves,
  and centroid/velocity tracking enters the narrative. The engine is already
  step-resolved; Rung 2 only adds the per-step `U_proj_bg` source and the
  centroid overlay.

## The energy stores (both runs share INQ's ledger)

`E_total = E_kinetic + E_hartree + E_xc + E_external (+ E_nonlocal + E_ion)`

The classical projectile is an **external potential** (like a pseudopotential), so
it lands in exactly the same stores as any external field — the decomposition is
structurally identical between the twins. What differs is *what physics sits in
each store*:

| Store | Wavepacket run | Classical (perturbation) run |
|---|---|---|
| `E_kinetic` | orbitals **+ the WP** (carries its localisation + `<p>²/2m`) | orbitals only |
| `E_hartree` | electron–electron **+ WP–electron + WP–WP self** | electron–electron only |
| `E_external` | electrons ↔ background | electrons ↔ background **+ electrons ↔ projectile potential** |
| `E_xc` | whole system **+ the WP's xc** | whole system |
| `U_proj_bg` | — (WP is an electron; already in E_H) | projectile ↔ background — **not in INQ's total**, tracked separately |

## Interpretation rules (the agent's knowledge base)

**Gauge caveat — read first.** WP insertion makes the cell net −1 charged; the
classical cell is neutral. Individual `E_hartree` and `E_external` differences are
therefore **Poisson-gauge-dependent** (the G=0 convention) and are NOT individually
physical. Only the **combination** `d(E_H+E_ext)` is gauge-clean. Never interpret
`dHartree` or `dExt` in isolation — always the sum. (See
`reference_charged_cell_hartree_convention` and the handover.)

1. **`dKin` = WP localisation kinetic.** At rest, `dKin ≈ 3/(4·σ_WP²)` (Ha) — the
   zero-point energy of confining the electron to width σ_WP, which a classical
   point charge does not pay. With k0≠0 add `k0²/2`. This is the cleanest quantum
   signature. (σ=0.5 → 81.6 eV.)

2. **`residual R = d(E_H+E_ext) − U_proj_bg` = WP self-Hartree `E_H[WP–WP]`.**
   `U_proj_bg` is what the difference *would* be if the WP were an ideal classical
   charge (no self-energy). The leftover is the WP's spurious self-repulsion.
   - Expected ≈ free-space `1/(σ_WP·√2π)` (Ha) as a **reference**; the
     boundary-matched (open-z) value is ~0.9 eV lower — the engine reports this
     `unexplained` gap, which is the open-z-vs-free-space gauge, not missing physics.
   - R is **r-independent** at rest → confirms it is a self-energy, not an
     interaction. If R drifts in a dynamic run, that is the WP self-Hartree
     changing as the WP **spreads** (σ grows) — a genuine quantum effect.

3. **`dXC` = the WP's own xc.** xc is local, so with the WP far from the slab this
   is the WP-alone xc and is **r-independent**. (σ=0.5 → −16.5 eV.)

4. **`SIE = R + dXC` = LDA one-electron self-interaction error (Perdew–Zunger).**
   For one electron, exact xc cancels the self-Hartree exactly; LDA under-cancels.
   `SIE` is the irreducible, physically-meaningful residue of the whole
   decomposition (σ=0.5 → 4.34 eV). It is density-dependent (shrinks with wider σ).

5. **Representation-awareness (which classical twin) — including the U_proj_bg SIGN.**
   The residual isolates the WP self-Hartree, but the **sign of `U_proj_bg` is
   representation-dependent** (`reference_ghost_u_proj_bg_sign`):
   - **`perturbation`** (Gaussian charge) → `R = d(E_H+E_ext) − U_proj_bg` ≈
     **20.81 eV**; the ~0.9 eV shortfall vs the free-space ref is the open-z gauge. Clean.
   - **`pseudopotential`** (ghost UPF) → `R = d(E_H+E_ext) + U_proj_bg` ≈ **8.85 eV**
     (ADD — INQ omits the z_valence=0 projectile↔background compensation term; using
     −U_proj_bg gives a spurious ~−260 eV). The residual sits ~12 eV **below** the
     clean perturbation value; that ~12–14 eV gap is the **ghost-UPF tail aliasing**
     (`reference_ghost_upf_tail_aliasing`), NOT missing physics — use the *ideal*
     `U_proj_bg`, never the impl term. Consequently **SIE = R + dXC is NOT clean for
     the pseudopotential** (R is aliasing-corrupted); take the physical 4.34 eV SIE
     from the perturbation twin. The engine applies the sign automatically from the
     manifest `representation`. Comparing the two representations (why perturbation
     closes the gap the ghost leaves) is a first-class output.

6. **Stopping power — classical projectile ONLY.** `−d(energy_proj_ke)/ds` over the
   initial near-constant-velocity window gives the drag on the **classical**
   projectile → hand `proj_ke_classical` to the `stopping-power-extraction` skill.
   **NEVER apply this to the WP.** The WP orbital is not identifiable as "the
   projectile", so its kinetic energy is not projectile KE
   (`feedback_quantum_stopping_not_from_projectile_ke`). The **total quantum
   stopping power** — the reason the localised-jellium system exists — is read from
   the **total electronic energy deposited** (`E_deposited_wp`), not any projectile
   track.

7. **Narrating a dynamic run.** Interpret the **first few timesteps explicitly**
   (which term moves first, by how much), then state the **general trend** — early
   induced changes cascade (a butterfly effect) and become hard to attribute after
   a while. Combine a rising/falling `residual` with the **WP centroid vs
   projectile position** (`separation_z`) and σ(t) (`wp_sigma_z`): a falling R with
   a spreading WP = self-Hartree dropping as the WP broadens; a residual that tracks
   `separation_z` = interaction leakage. **Check energy conservation first**
   (`res.conservation`): a drifting `E_conserved_*` means an integrator bug, not
   physics. Label the attributable early differences the "quantum effects"; flag the
   rest as downstream cascade.

## Pairwise Coulomb decomposition — every difference is physically attributable

When both runs emit `interactions.csv` (the generation runs now do), the engine
resolves the lumped `E_hartree`/`E_external` into the **six pairwise Coulomb terms**
of the three charge groups **P** (projectile), **S** (slab electrons), **B**
(background), plus the kinetic split and `E_xc`. This is the primary lens: **look at
each individual energy component, compare it between classical and WP, and attribute
every difference to a named physical cause** (`reference_twin_pairwise_decomposition`).

- **The terms** (Hartree in the CSV; engine reports eV): `E_SS` slab–slab,
  `E_PP` projectile self-Hartree, `E_PS` projectile–slab, `E_SB` slab–bg,
  `E_PB` projectile–bg, `E_BB` background self (constant). Kinetic: `KE_total`
  (INQ), `KE_proj` (classical ½mV² / WP orbital KE), `KE_slab` = diff. Plus `E_xc`.
- **Closure gates** (must hold to ~1e-9 — the master sanity check):
  classical `E_hartree`=E_SS, `E_external`=E_SB+E_PS ; WP `E_hartree`=E_SS+E_PS+E_PP,
  `E_external`=E_SB+E_PB. The engine's `pairwise` table + `e_hartree_check`/
  `e_external_check` columns verify this.
- **The gauge test** (the engine runs it automatically): the physically-identical
  terms E_SS, E_SB, E_BB (same slab & background in both runs) must have **Δ≈0**.
  Verified: ΔE_SS=ΔE_SB=ΔE_BB=0.0000 eV → **no inter-run gauge, zero-points agree**.
  So `res.gauge["no_gauge"]` should be true; if not, a real gauge is present and the
  per-term comparison must be gauge-corrected before interpreting.
- **Attribution**: with no gauge, EVERY non-zero Δ is physical. The only differing
  terms are the projectile ones — e.g. ΔE_PP = −0.48 eV at step 1 = WP self-Hartree
  falling as the packet disperses (E_PP∝1/σ). Read `res.pairwise_table(step)` and
  name each Δ (dispersion, polarisation, trajectory).

The absolute E_SB/E_PB/E_BB still carry the p2 open-z convention, but it is the SAME
in both runs, so it cancels in every Δ (that is what the gauge test confirms).

## Workflow

1. **Validate the pair.** Confirm `twin_manifest.json: valid=true` (or run
   `twin-run-generation/check_twin.py`). Do not analyse an unvalidated pair.
2. **Run the engine:**
   ```bash
   /local/data/public/skcb2/tddft/venv/bin/python3 \
     .claude/skills/twin-run-analysis/twin_decompose.py <pair_dir> \
     --json out.json --csv steps.csv
   ```
   or `--wp DIR --classical DIR`. Read `res.report()` / the findings table.
3. **Narrate** using the interpretation rules. State each term's value, its
   attribution, and the unexplained remainder. Call out the gauge caveat.
4. **Build the analysis notebook** with `twin_notebook_builder.py` (skill-local):
   it embeds the findings table, the per-step decomposition plot, the residual/SIE,
   the pairwise ledger + gauge test, the n(z,t) density carpets + Δn, the WP−classical
   bar plot, the pairwise-energy GIF, and — MANDATORY (rule `notebook-density-gif`) —
   the **`density_evolution.gif`** (n(z,t) classical vs WP + Δn, animated). A twin
   notebook without the density GIF is incomplete. Output lives in the run-set's
   `hypotheses/<sweep>/` folder (ADR 0007), never in the skill. Runs feeding a
   notebook must save density frames (`LJ_SAVE_EVERY>0`).

## Golden numbers (Rung-1 regression — σ_WP=0.5, r=12, Lz120, p2)

`dKin 81.74 (loc 81.63) | dXC −16.47 | residual 20.81 | SIE 4.34 eV`

The test `tests/test_twin_decompose.py` asserts these against both a synthetic
fixture and the on-disk golden pair
(`…/proj_perturbation/{results/proj_pert_dx0p5, stress_scratch/s0p5_r12_lz120_p2/results/wp}`).
Run before trusting the engine on new data:
```bash
/local/data/public/skcb2/tddft/venv/bin/python3 -m pytest \
  .claude/skills/twin-run-analysis/tests/test_twin_decompose.py -q
```

## Files (all skill-local, shippable)

| File | Role |
|---|---|
| `twin_decompose.py` | deterministic decomposition engine + CLI |
| `twin_notebook_builder.py` | assembles the executed analysis `.ipynb` |
| `tests/test_twin_decompose.py` | synthetic + golden-pair known-answer tests |

Reference: `docs/handovers/localised-jellium-energy-book-keeping.md`,
`docs/notes/gaussian-pertubation-for-classical-simul`.
