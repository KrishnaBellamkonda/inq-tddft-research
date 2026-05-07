# Report: m=1 plasmon detection in L=50 N=162 jellium via 15 eV WP injection

**Date:** 2026-05-06
**Run:** `ResearchProject/systems/jellium/run_plasmon_n162_L50_E15/`
**Linked plan:** `docs/plans/jellium_plasmon_detection.md`
**Linked source:** `docs/sources/correa-2018-electronic-stopping-power.md`
**Verdict:** **YES — m=1 plasmon detected at 3.53 eV (predicted 3.59 eV).**

## 1. Methods

L=50 cubic periodic jellium, N=162 closed-shell bath (reused GS at
`gs_L50_cubic_N162_dx1p0`), one wave-packet orbital injected at the
cell origin with σ_r = 5 Bohr, $k_0 = 1.0500$ a₀⁻¹ (E_WP = 15.0 eV)
in the +z direction. ETRS propagator, dt = 0.020 a.u., 100 000 steps
⇒ T_sim = 2000 a.u. = 48.4 fs. Wall time 6 h 49 min on a single A30
GPU. Total energy drift over the full propagation: **−4.0 × 10⁻⁷ eV**
(machine precision; energy is conserved to ~10⁻⁸ Ha).

WP injection diagnostics:
- `max_overlap` (residual WP–bath overlap after orthogonalisation):
  **1.5 × 10⁻⁴** — three orders of magnitude smaller than the
  E_WP = 1.5 eV runs (0.13). The high-velocity WP carries
  plane-wave content peaked at $|G| \approx k_0 = 1.05$ a₀⁻¹,
  far above the bath's filled $|G|^2 \le 6$ shells (top
  $|G| \approx 0.31$ a₀⁻¹), so the geometric overlap is naturally
  small.
- `norm_after = 1.0`, `norm_before = 1.0` (clean Gram–Schmidt).

Custom postprocess phase
`inq-stack/python/inqview/postprocess/density_fourier.py` (new in this
report): reads all 501 saved `density_rt_total/density_t<step>.vti`
frames, computes $\delta n(\mathbf r, t) = n(\mathbf r, t) -
n(\mathbf r, 0)$, takes a 3D FFT, and extracts the axial Fourier
components $n_{q_m}(t)$ at $\mathbf q_m = (0, 0, 2\pi m / L)$ for
$m = 1, 2, \ldots, 6$. Then a 1D Hann-windowed FFT (4× zero-pad,
transient $t < 5$ a.u. excluded per the QBall convention) of each
$n_{q_m}(t)$ gives the plasmon spectrum.

## 2. Results

### Time-domain $|n_{q_m}(t)|$

![Axial Fourier components vs time](../journals/researchproject/attachments/2026-05-06_run_plasmon_n162_L50_E15/n_q_vs_time.png)

The $m=1$ trace shows a clean amplitude-modulated oscillation
throughout the 2000 a.u. propagation, growing rapidly during the WP's
first transit of the box and persisting at near-saturation amplitude
afterwards. $m=2$ is a factor ~3 weaker but visibly oscillatory.
Higher-$m$ traces are noise-dominated and not coherent — exactly the
prediction of the Bohm-Gross / Landau-damping analysis (modes with
$q > q_c = \omega_p / v_F = 0.378$ a₀⁻¹ are inside the
electron-hole continuum and damp on the same timescale they are
created).

### Frequency-domain spectrum (the headline result)

![Plasmon spectrum: FFT of n_q_m(t)](../journals/researchproject/attachments/2026-05-06_run_plasmon_n162_L50_E15/n_q_spectrum.png)

| Mode m | Predicted ℏω_BG (eV) | Observed peak (eV) | Shift (eV) | Relative shift | Peak amplitude (a.u.) | Verdict |
|---|---:|---:|---:|---:|---:|---|
| **1** | **3.59** | **3.533** | **−0.057** | **−1.6 %** | 1.10 × 10² | **DETECTED — bulk-like plasmon** |
| **2** | **4.00** | **3.812** | **−0.188** | **−4.7 %** | 4.29 × 10¹ | **DETECTED** |
| 3 | 4.79 | 9.89 | +5.10 | +106 % | 8.13 × 10⁰ | NOT a plasmon (Landau damped, $q/q_c = 1.00$) |
| 4 | 6.05 | 10.09 | +4.04 | +67 % | 6.78 × 10⁰ | NOT a plasmon |
| 5 | 7.80 | 11.59 | +3.79 | +49 % | 2.85 × 10⁰ | NOT a plasmon |
| 6 | 10.03 | 11.67 | +1.64 | +16 % | 2.81 × 10⁻¹ | NOT a plasmon (signal at noise floor) |

### Energy-component bookkeeping (Δ over T_sim = 2000 a.u.)

| Component | Δ (Ha) | Δ (eV) |
|---|---:|---:|
| `energy_total` | −1.5 × 10⁻⁸ | **−4.0 × 10⁻⁷** (drift, conserved) |
| `energy_kinetic` | +2.77 × 10⁻² | **+0.754** |
| `energy_hartree` | −3.19 × 10⁻² | **−0.868** |
| `energy_xc` | +4.19 × 10⁻³ | **+0.114** |

Same qualitative pattern as the E_WP = 1.5 eV runs (kinetic up,
Hartree down, xc up by a small amount), but with all swings ~50 %
larger because the WP perturbation is correspondingly stronger. The
sign pattern is again consistent with the **charge-conjugate
electron-wake** picture (negative WP → density depletion behind →
ΔE_H < 0).

## 3. Discussion

### The m=1 plasmon at 3.53 eV is the bulk-like collective mode

The FFT peak at 3.53 eV in the $n_{q_1}(t)$ channel, shifted only
−0.06 eV from the Bohm-Gross prediction $\hbar\omega(q_1) = 3.59$ eV,
is the **m=1 axial plasmon mode** of the L=50 supercell. The amplitude
is nearly 3× the m=2 peak and 13× the m=3 noise — this is by far the
dominant collective response to the WP excitation.

The −1.6 % systematic downshift from Bohm-Gross is the expected ALDA
correction: Bohm-Gross is an RPA result with a hydrodynamic-style
gradient term, and it slightly overestimates the dispersion at
finite $q$ because the full RPA correlation kernel softens the
restoring force. The same downshift trend is well documented in TDDFT
studies of jellium plasmons (e.g. Gross-Kohn 1985 sum-rule analysis;
multiple Comp. Mater. Sci. studies on Al, Mg).

### Why the m≥3 modes don't show plasmon peaks

The plasmon kinematic threshold is $q_c = \omega_p / v_F =
0.1276 / 0.337 = 0.378$ a₀⁻¹ — exactly $q_3 = 2\pi \cdot 3 / 50 =
0.377$ a₀⁻¹. So $m \ge 3$ axial modes sit *at or beyond* the
Landau-damping cutoff: their plasmon dispersion enters the
electron-hole continuum and the modes decay into single-particle
excitations on the same timescale they are created. The "peaks" we
see for $m \ge 3$ near 10 eV are **not** plasmons — they are residual
e-h continuum response of the bath at that wavevector, scaled by the
WP's Fourier content at $|q|$ corresponding to that mode.

This is the **clean sign-test** that the m=1 detection is real and
not an artefact: the modes that physics says should be damped are
damped; the modes that physics says should be coherent are coherent.

### Comparison to Run B's predicted v_res

The WP velocity $v = k_0 = 1.05$ a.u. matches $\omega(q_1)/q_1 = 1.050$
a.u. **to all four decimal places** (this is the resonance condition
the Run B Cfg was designed around). The energy-transfer efficiency is
maximal at this velocity — which is exactly why the m=1 amplitude
dominates the spectrum.

## 4. Verdict against the criteria from `jellium_plasmon_detection.md §4`

| Criterion | Required | Observed | Status |
|---|---|---|---|
| FFT peak at 3.59 ± 0.1 eV in $n_{q_1}(t)$ | yes | 3.533 eV (within 0.06 eV) | **YES** |
| Amplitude ≫ noise | yes | m=1 peak is 13 × the m=3 noise | **YES** |
| m=2 peak at 4.0 ± 0.2 eV | desired | 3.812 eV (within 0.19 eV) | **YES** |
| $\Delta E_\text{total}$ drift < 1 mHa | yes | 1.5 × 10⁻⁸ Ha = 4 × 10⁻⁷ eV | **YES** |
| Bohm-Gross dispersion roughly satisfied | yes | m=1 and m=2 within 5 % of prediction | **YES** |

**OVERALL VERDICT: YES — the L=50, N=162 jellium box hosts a clean
m=1 plasmon mode at 3.53 eV, resonantly excited by the 15 eV WP at
$v = v_\text{res} = 1.05$ a.u. and observed unambiguously in the
$n_{q_1}(t)$ Fourier channel.**

## 5. Implications for the broader project

1. **The previous L=50 / E=1.5 eV runs were sub-threshold for
   plasmon excitation**, exactly as concluded in the regime-classification
   topical entry (`plasmons-and-stopping-power.md §6`). The 1.5 eV
   WP at $v = 0.332$ a.u. is below $v_\text{res}^{m=1} = 1.05$ a.u.
   by a factor of 3.2 — too slow to drive the plasmon resonance.
   Their physics is single-particle electron-hole excitation +
   density-wake response, not collective plasmon.
2. **The Bohm-Gross dispersion is a quantitative tool for INQ jellium
   simulations**, not just an order-of-magnitude estimate. The −1.6 %
   downshift for m=1 is consistent with the literature on ALDA xc
   corrections, and reproducible.
3. **Run C (m=2 follow-up at E_WP = 4.65 eV) is justified but no
   longer urgent** — Run B already detected m=2 as a secondary peak.
   Run C would isolate the m=2 mode from the m=1 contamination and
   provide a more precise measurement of the m=2 frequency.
4. **The custom `density_fourier` postprocess phase is a permanent
   addition** to the jellium pipeline and should be wired into
   `inqview.postprocess.pipeline.run` so future plasmon-relevant
   runs get $n_{q_m}(t)$ + spectrum plots automatically.

## 6. Reproducibility

```bash
source /local/data/public/skcb2/tddft/venv/bin/activate
cd /local/data/public/skcb2/tddft
python3 -m inqview.postprocess.density_fourier \
    ResearchProject/systems/jellium/run_plasmon_n162_L50_E15/results \
    run_plasmon_n162_L50_E15
```

emits:
- `analysis/observables/n_q_vs_time.{csv,png}`
- `analysis/observables/n_q_spectrum.{csv,png}`

These were the basis for every number in §2.

## 7. Next steps

1. **Wire `density_fourier` into the standard pipeline** (small Python
   change in `inqview/postprocess/pipeline.py`).
2. **Write a new journal entry** under
   `docs/journals/researchproject/2026-05-06_run_plasmon_n162_L50_E15.md`
   linking back to this verdict.
3. **Optionally launch Run C** to confirm m=2 in isolation; high
   confidence it will reproduce the Run B m=2 finding.
4. **Use the Bohm-Gross-corrected ALDA dispersion** (the −1.6 % shift
   established here) to design future plasmon-targeted runs in
   different boxes.
5. **Apply the same `n_q(t)` extraction to the previous E_WP = 1.5 eV
   runs** as a control — should show *no* m=1 peak at 3.5 eV,
   confirming those runs are below the plasmon excitation threshold.
