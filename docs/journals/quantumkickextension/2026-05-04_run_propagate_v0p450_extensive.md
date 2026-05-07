# 2026-05-04_run_propagate_v0p450_extensive

**Title:** run_propagate_v0p450_extensive (li_54_atom_bcc, v=0.450)
**Run path:** `/local/data/public/skcb2/tddft/QuantumKickExtension/inq-codebase/Li/run_propagate_v0p450_extensive`
**Linked results:** `/local/data/public/skcb2/tddft/QuantumKickExtension/inq-codebase/Li/run_propagate_v0p450_extensive/results`
**Status:** running

## Config snapshot

(From `run_propagate_v0p0123_extensive`'s sibling cpp; will be replaced by
the live `run_summary.txt` once the propagation finishes. Differences from
that run are only in `kick_velocity_au` and the run name.)

| Field | Value |
|---|---|
| run | run_propagate_v0p450_extensive |
| system | li_54_atom_bcc_supercell |
| cell_angstrom | 10.53 10.53 10.53 |
| boundary | periodic |
| n_atoms | 54 |
| k_grid | 2 2 2 shifted |
| smearing | fermi_dirac |
| smearing_temperature_kelvin | 400 |
| xc | pbe |
| cutoff_ry | 74 |
| extra_states | 20 |
| checkpoint_dir | ../checkpoints/li_54_2x2x2_T200 |
| num_states | 101 |
| num_electrons | 162 |
| kick_velocity_au | 0.450 |
| kick_direction | +x |
| atoms_dynamics | impulsive |
| dt_au | 0.04 |
| n_steps | 15500 |
| total_time_fs | 14.997 |
| write_every | 100 |
| GPU | CUDA_VISIBLE_DEVICES=1 |

## Observations

### Where does v = 0.450 a.u. sit relative to the Fermi velocity of Li?

For Li in the BCC primitive cell (a = 3.51 Å, 2 atoms / primitive cell,
**1 valence-2s electron per atom contributing to the conduction
manifold**), the free-electron Fermi level computed from first principles
is

> n_val = 1 / V_prim = 6.85 × 10⁻³ Bohr⁻³ → k_F = (3π² n_val)^{1/3} = **0.588 a.u.**
> E_F = ½ k_F² = **4.70 eV** → **v_F = k_F = 0.588 a.u. = 1.286 × 10⁶ m/s**

This matches the textbook value (Ashcroft & Mermin Table 2.1 gives
v_F = 1.29 × 10⁶ m/s for Li). The 54-atom GS run produces eigenvalues
consistent with this picture — the 1s² core sits at ≈ −47 eV
(states 0–53), and the partially-filled 2s manifold runs from ≈ −5 eV up
to the Fermi level around state 79 with a small Fermi-Dirac tail
extending into bands 80–88 at 400 K.

| | v_kick | v_kick (SI) | v_kick / v_F | Mv²/N_uc (eV) |
|---|---:|---:|---:|---:|
| previous run (low) | **0.0123 a.u.** | **2.69 × 10⁴ m/s** | **0.021** | 0.025 |
| this run (high) | **0.450 a.u.** | **9.85 × 10⁵ m/s** | **0.766** | 33.06 |
| Li Fermi velocity | 0.588 a.u. | 1.286 × 10⁶ m/s | 1.000 | — |

### Categorisation

Both runs sit **below v_F**, so the perturbation is in the **electronic-
stopping regime** in both cases — the energy lost by the moving ion
lattice goes into electronic excitations across the Fermi surface, not
into phonons. The ions follow `atoms_dyn = impulsive` so they cannot
absorb energy themselves; everything observed in ΔE(t) is an electronic
response.

The two runs explore different sub-regimes of electronic stopping:

- **v = 0.0123 a.u. (deep linear-response, v/v_F ≈ 0.02).** Lindhard
  regime. The instantaneous stopping power scales as `dE/dt ∝ v²`
  (Mv²-rule integrated form), and the QBall paper's countercurrent
  prediction gives an asymptotic plateau at α(v) Mv² with α ≤ 1. The
  kick energy budget per unit cell is small (≈ 25 meV), comparable to
  the 400 K Fermi smearing window (kT ≈ 35 meV), so the relative
  fluctuation around the plateau is large — explaining the wonky-looking
  ΔE(t) trace in the previous entry.
- **v = 0.450 a.u. (near-peak linear-response, v/v_F ≈ 0.77).** Still
  Lindhard, but at the velocity where the electronic stopping power is
  near its peak (the Bragg peak for slow ions sits around v ≈ v_F in
  free-electron-like metals). The kick energy budget per unit cell is
  ≈ 33 eV, **a factor of ~1340× larger than the v = 0.0123 run** — well
  above any thermal-broadening scale, so the plateau should be sharply
  defined and the FFT signal will be much larger.

### Inference

In SI: the v = 0.450 kick gives the lattice a velocity of **985 km/s**
(0.0033 c). For comparison, the Bohr velocity v_Bohr = αc ≈ 2188 km/s
sets the natural scale of stopping-power physics; v_F for Li is at
0.59 v_Bohr, and our kick is at 0.45 v_Bohr — squarely in the slow-ion /
Lindhard regime, not in the Bethe-Bloch fast-ion regime.

The headline expectation for this run is therefore:

1. ΔE(t) builds up to a plateau **near 33 eV per unit cell**
   (≈ 0.6 eV per electron), with the countercurrent-reduction factor
   α(v) measured by `analyse_inq.py`.
2. The FFT peak position is set by the metal's plasmon/Lindhard response
   at this density. From the v = 0.0123 run we already saw a sharp peak
   at 5.72 eV; the v = 0.450 run should reproduce this position
   (linear-response peak) but with a much larger amplitude.
3. The Mermin-baseline negative-excursion artefact from the low-v run
   should be invisible here — the kick energy budget is 1340× the
   thermal smearing scale, so any negative excursions of ΔE will be
   tiny relative to the plateau.

(Pending — to be filled in once propagation completes.)

## Open questions / next steps

- Pending: live ΔE(t), FFT spectrum, and INQ-vs-QBall comparison plot
  once the run finishes (~20 h on GPU 1, ETA mid-day 2026-05-05).
- Pending: confirm whether the FFT peak position at v = 0.450 matches
  the v = 0.0123 result (linear-response → same plasmon-like
  resonance) or shifts (signature of approach to Bragg peak).
- Pending: extract α(v) = ΔE_plateau / (Mv²) and update the linear-
  response check from the QBall analyse.py.
