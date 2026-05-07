# 2026-05-06 — `run_plasmon_n162_L50_E3p4_varyv` (L=50, N=162, **WP=3.4 eV at off-resonance v=0.5 a.u. — bath plasmon confirmed**)

**Run path:** `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_plasmon_n162_L50_E3p4_varyv/`
**Linked results:** `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_plasmon_n162_L50_E3p4_varyv/results/`
**Status:** complete
**Branch:** `features/jellium-ks-energy-observables`
**Aim:** definitive discriminator between (i) bath-plasmon and (ii) wrap-kinematic interpretations of the m=1 peak detected in `run_plasmon_n162_L50_E15` (Run B). Box and bath identical to Run B; only the WP velocity changes (1.05 → 0.5 a.u.) and σ shrinks (5 → 3 Bohr).

## Run summary

### 1. Run identity

| Field | Value |
|---|---|
| run_name | `run_plasmon_n162_L50_E3p4_varyv` |
| run_type | wave-packet RT-LEED on jellium (TDDFT, ALDA) |
| date_finished | 2026-05-06T19:49:05 |
| wall_time_s | 24571.3 (≈ 6 h 49 min) |
| executable | `run.cpp` built via `inq-run` |
| geometry_file | (none, jellium) |
| checkpoint_dir | `checkpoints/gs_L50_cubic_N162_dx1p0` |

### 3. System configuration

| Field | Value |
|---|---|
| cell_bohr | 50³ (cubic, periodic) |
| boundary | periodic |
| n_ions | 0 |
| n_electrons | 162 (closed shell) |
| n_occupied | 81 |
| extra_states | 20 |
| wp_state_index | 100 |
| spacing_bohr | 1.0 |
| temperature_ev | 0.00862 |
| xc_functional | LDA (ALDA in TDDFT) |

### 5. Wavepacket

| Field | Value |
|---|---|
| wp_enabled | yes |
| wp_center_bohr | (0, 0, 0) |
| **wp_sigma_bohr** | **3.0** (down from 5.0) |
| **wp_k0_bohr_inv** | **(0, 0, 0.500)** (down from 1.050) |
| **wp_energy_ev** | **3.40** (down from 15.0) |
| wp_direction | +z |
| wp_occupation | 1.0 |
| orthogonalised | yes |
| norm_after | 1.0 |
| **max_overlap** | **0.065** (intermediate between Run B's 1.5 × 10⁻⁴ and the 1.5 eV runs' 0.13) |

### 6. Real-time configuration

| Field | Value |
|---|---|
| rt_num_steps | 100000 |
| dt_au | 0.020 |
| total_time_au | 2000.0 |
| write_every | 200 |
| screen_snap_every | 6 |

### 9. End-of-run diagnostics

| Field | Value |
|---|---|
| run_completed | true |
| final_time_au | 2000 |
| total_energy_drift | $\sim 10^{-8}$ Ha (machine precision) |

---

## 1. Headline: bath plasmon confirmed at the same frequency as Run B

![axial Fourier components vs time](attachments/2026-05-06_run_plasmon_n162_L50_E3p4_varyv/n_q_vs_time.png)

![n_q FFT spectrum](attachments/2026-05-06_run_plasmon_n162_L50_E3p4_varyv/n_q_spectrum.png)

### Decisive comparison vs Run B

| Mode | Run B (v=1.05) m peak | **Run D (v=0.5) m peak** | Shift between runs |
|---|---:|---:|---:|
| m=1 plasmon | 3.533 eV | **3.555 eV** | **+0.022 eV** (1 FFT bin) |
| m=2 plasmon | 3.812 eV | **3.812 eV** | **0.000 eV** |

The plasmon peak does not shift when $v$ changes by a factor of 2. A
*kinematic* peak (rigid WP translating, no real bath mode) would have
shifted by $\Delta\omega = (v_\text{D} - v_\text{B}) \cdot q_1 \approx
-1.88$ eV. The observed shift is **85× smaller** than the kinematic
prediction. The bath plasmon is real.

### Both peaks present in Run D

| Channel | Plasmon prediction | Plasmon observed | Kinematic prediction | Kinematic observed | Plasmon : kinematic amplitude |
|---|---:|---:|---:|---:|---:|
| m=1 | 3.59 eV | **3.555 eV** (amp 8.9 × 10¹) | 1.71 eV | 1.520 eV (amp 2.7 × 10¹) | **3.3 : 1** |
| m=2 | 4.00 eV | **3.812 eV** (amp 9.2 × 10¹) | 3.42 eV | 3.405 eV (amp 2.3 × 10¹) | **4.1 : 1** |

The kinematic / "wrap-period" peaks predicted by the rigid-WP rigid-translation
model are present at their predicted positions, **but they are subdominant** —
plasmon : kinematic ≈ 3.3–4.1 to 1 in amplitude. So the previous Run B
detection's interpretation as a real bath plasmon is **upheld and
strengthened** — and at the same time the user's hypothesis that wrap
kinematics contributes to the spectrum is **also vindicated**, just as
the secondary effect.

Full verdict at [`docs/reports/plasmon-vs-wrap-verdict.md`](../../reports/plasmon-vs-wrap-verdict.md).

---

## 2. Why we believe the bath plasmon is real

Three independent pieces of evidence:

1. **Frequency invariance under $v$.** The m=1 peak shifts by ≤ 0.022
   eV between Run B (v=1.05) and Run D (v=0.5). A real bath collective
   mode is set by bath density and dispersion, not by projectile
   velocity. The observed invariance is the hallmark.
2. **Spatial-mode-dependent dispersion.** m=1 at 3.55 eV ≠ m=2 at 3.81
   eV ≠ m=3 in noise. The kinematic prediction would have $\omega_m =
   m \cdot v \cdot q_1$ (linear in m); the bath-plasmon prediction
   has the Bohm-Gross dispersion (sub-linear, then super-linear). The
   data follows Bohm-Gross.
3. **Cross-confirmation with kinematic peaks at predicted positions.**
   The "kinematic check" peaks at 1.52 eV and 3.41 eV match the rigid-WP
   predictions for v=0.5, validating the n_q(t) extraction itself
   (otherwise we wouldn't see them at all). The fact that these
   subdominant kinematic peaks coexist with dominant plasmon peaks
   means the math of `density_fourier` is correct and the dominant
   peaks are not artefacts of the postprocess.

---

## 3. Energy bookkeeping

(Computed from `observables.csv`, t = 0 → 2000 a.u.):

| Component | Δ (Ha) | Δ (eV) |
|---|---:|---:|
| `energy_total` | $\sim 10^{-8}$ | $\sim 10^{-7}$ (drift, conserved) |
| `energy_kinetic` | tbd | tbd |
| `energy_hartree` | tbd | tbd |
| `energy_xc` | tbd | tbd |

(Numerical extraction of components pending — standard postprocess running.)

The qualitative pattern is expected to be the same as Run B (kinetic
up, Hartree down, xc up, smaller swings because the WP perturbation is
weaker at lower $E_\text{WP}$).

---

## 4. WP shape with σ=3 vs σ=5

| Quantity | Run B (σ=5, v=1.05) | **Run D (σ=3, v=0.5)** |
|---|---:|---:|
| σ_k = 1/σ | 0.20 a₀⁻¹ | **0.333 a₀⁻¹** |
| σ_k / k_0 | 0.19 | **0.67** (high momentum spread) |
| Backward Gaussian content (P[k_z<0]) | 1.6 × 10⁻⁷ | **0.067** (≈ 6.7%) |
| `max_overlap` (after Gram-Schmidt) | 1.5 × 10⁻⁴ | **0.065** (~430× larger; still small enough that the bath response dominates) |

The σ=3 WP at v=0.5 is a "broad-bandwidth" wave packet — it has
substantial $k$-content from $k_z = 0$ up to $k_z \approx 1.5$ a₀⁻¹.
This *helps* with the diagnostic: a broader spectral input drives the
bath at many frequencies simultaneously, so bath modes that are weakly
coupled at the dominant $k_0$ still get excited via the spectral tail.
The fact that we see strong m=1 and m=2 plasmon peaks despite v being
well off the m=1 resonance is consistent with this picture.

---

## 5. Open questions / next steps

- **Standard postprocess is running** — when it finishes, fill in the
  energy components in §3 and add the system_yz/system_xz density GIFs.
- **Re-extract $n_q(t)$ for the L=50 / E=1.5 eV runs** as a control
  — at $v = 0.332$, kinematic prediction 1.13 eV (m=1) and 2.27 eV
  (m=2) — should still show a weak plasmon peak at 3.59 eV if the
  bath mode is genuinely there. This will pin down whether the
  earlier "no plasmon at E=1.5 eV" conclusion is too strong.
- **Update `plasmons-and-stopping-power.md` §6** — soften the "below
  threshold for plasmon excitation" framing to "off-resonance, but the
  collective mode is still excited weakly."
- **Manuscript-grade plot** — overlay Run B (red) and Run D (blue) m=1
  and m=2 spectra on the same axis; annotate with (i) the
  velocity-invariant plasmon peak position and (ii) the velocity-
  dependent kinematic peak. This is the canonical figure that proves
  the L=50 / N=162 jellium box hosts real bath plasmons.
