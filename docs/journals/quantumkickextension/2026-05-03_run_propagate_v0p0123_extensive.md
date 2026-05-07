# 2026-05-03_run_propagate_v0p0123_extensive

**Title:** run_propagate_v0p0123_extensive (li_54_atom_bcc, v=0.0123)
**Run path:** `/local/data/public/skcb2/tddft/QuantumKickExtension/inq-codebase/Li/run_propagate_v0p0123_extensive`
**Linked results:** `/local/data/public/skcb2/tddft/QuantumKickExtension/inq-codebase/Li/run_propagate_v0p0123_extensive/results`
**Status:** complete

## Config snapshot

| Field | Value |
|---|---|
| run | run_propagate_v0p0123_extensive |
| system | li_54_atom_bcc_supercell |
| cell_angstrom | 10.53 10.53 10.53 |
| boundary | periodic |
| n_atoms | 54 |
| k_grid | 2 2 2 shifted |
| smearing | fermi_dirac |
| smearing_temperature_kelvin | 400 |
| smearing_temperature_target_kelvin | 200 (fell back to 400 in GS — see GS run_summary.txt) |
| xc | pbe |
| cutoff_ry | 74 |
| extra_states | 20 |
| checkpoint_dir | ../checkpoints/li_54_2x2x2_T200 |
| num_states | 101 |
| num_electrons | 162 |
| kick_velocity_au | 0.0123 |
| kick_direction | +x |
| atoms_dynamics | impulsive |
| dt_au | 0.04 |
| n_steps | 15500 |
| total_time_fs | 14.9970808 |
| write_every | 100 |
| vti_format | binary |
| qball_kick_dir | ../../qball-codebase/Li/td_kicks/results/kick_v0.0123 |

## Observations

The excess-energy-per-unit-cell time-series for v = 0.0123 a.u. shows a
wonky low-frequency wobble: the second-half plateau is poorly defined,
and the negative excursions of ΔE(t) (signature of the Mermin baseline
at finite T) interact with the small kick energy budget so that the
oscillation around the plateau is large in relative terms.

![Excess energy per unit cell vs time](attachments/2026-05-03_run_propagate_v0p0123_extensive/excess_energy_per_uc_vs_time.png)

This same pattern was observed in the previous attempt at v = 0.0123 a.u.
with a different INQ configuration, so it is not a one-off numerical
artefact of the present GS/checkpoint pair — it appears to be a
low-velocity feature of the simulation that needs deeper investigation.

The FFT of ΔE(t) (Hann + 8× zero-pad + plateau-detrend) gives a sharp
peak at 5.72 eV with a second harmonic at 11.4 eV; the QBall reference
(Γ-only, 1000 K, 3.87 fs) gave 6.5 eV with no resolved second harmonic.
The shift of −0.78 eV under the configuration upgrade is the headline
result.

![FFT of excess energy vs ω](attachments/2026-05-03_run_propagate_v0p0123_extensive/fft_excess_energy_vs_omega.png)

![INQ vs QBall iter-5 comparison](attachments/2026-05-03_run_propagate_v0p0123_extensive/comparison_inq_vs_qball_iter5.png)

The primary aim of this section of the project is to lock down a
reliable data-production pipeline. To avoid spending more wall-time
chasing a wonky low-v case before the pipeline is validated, we move
on to the highest QBall kick velocity (v = 0.450 a.u.) next. The
low-velocity anomaly returns as a follow-up question once the
high-velocity end of the curve is well-characterised.

## Open questions / next steps

- Why does ΔE(t) at very low v have such a noisy plateau — is it
  Mermin-baseline noise, k-grid coarseness, or finite simulation length?
- Re-run with smaller dt (0.02 a.u.) and double-length propagation once
  the pipeline at v = 0.450 a.u. is validated.
- Quantify the negative-excursion contribution: is it the Mermin
  free-energy baseline (F = E − TS) or numerical drift?
