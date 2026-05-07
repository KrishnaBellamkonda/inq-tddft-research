# Plan: plasmon-detection programme for L=50 / N=162 jellium

**Status:** plan written; **Run B launching automatically** (see §4
below); Runs A and C scaffolded.
**Linked entries:**
- `docs/journals/researchproject/plasmons-and-stopping-power.md`
  (regime classification — concluded the previous L=50 / E=1.5 eV runs
  were sub-threshold for kinematic plasmon excitation).
- `docs/sources/correa-2018-electronic-stopping-power.md`.

## 0. Verified numbers

Every number in the proposal was verified by direct calculation against
the Bohm-Gross / hydrodynamic plasmon dispersion
$\omega(q)^2 = \omega_p^2 + (3/5)v_F^2 q^2 + q^4/4$ and the resonance
condition $v_m^\text{res} = \omega(q_m)/q_m$.

| Quantity | Computed | Proposal | Match |
|---|---|---|---|
| Density $n = N/L^3$ | 1.296 × 10⁻³ a₀⁻³ | 1.296 × 10⁻³ | ✓ |
| $\omega_p = \sqrt{4\pi n}$ | 0.1276 a.u. = 3.473 eV | 3.47 eV | ✓ |
| $T_p = 2\pi/\omega_p$ | 49.24 a.u. = 1.191 fs | 1.19 fs | ✓ |
| $k_F = v_F$ | 0.3373 a.u. | 0.337 a.u. | ✓ |
| **Landau-damping cutoff** $q_c = \omega_p/v_F$ | **0.378 a₀⁻¹** | **(implied: m=3 q ≈ 0.377)** | ✓ confirms m=3 sits exactly at the e-h continuum edge |
| **m=1**: q, ω, ℏω, T_p, v_res, E_WP | 0.1257; 0.1320; **3.59 eV**; **1.151 fs**; **1.050**; **15.01 eV** | 3.59 / 1.15 / 1.050 / 15.0 | ✓ |
| **m=2**: same | 0.2513; 0.1470; **4.00 eV**; **1.034 fs**; **0.585**; **4.65 eV** | 4.00 / 1.03 / 0.585 / 4.65 | ✓ |
| **m=3**: same | 0.3770; 0.1762; **4.79 eV**; **0.863 fs**; **0.467**; **2.97 eV** | 4.79 / 0.863 / 0.467 / 2.97 | ✓ |
| **m=4**: same | 0.5027; 0.2225; **6.05 eV**; **0.683 fs**; **0.443**; **2.67 eV** | 6.05 / 0.683 / 0.443 / 2.67 | ✓ |
| **m=5**: same | 0.6283; 0.2867; **7.80 eV**; **0.530 fs**; **0.456**; **2.83 eV** | 7.80 / 0.530 / 0.456 / 2.83 | ✓ |
| **m=6**: same | 0.7540; 0.3686; **10.03 eV**; **0.412 fs**; **0.489**; **3.25 eV** | 10.0 / 0.412 / 0.489 / 3.25 | ✓ |
| FFT resolution at T_sim = 2000 a.u. | $\Delta E = 2\pi/2000 \cdot$ Ha2eV = 0.0855 eV | 0.0855 eV | ✓ |
| m=1 vs m=2 separation | 0.407 eV (4.7 × ΔE) | 0.41 eV | ✓ |
| WP de Broglie at k₀=1.05 | $\lambda_\text{dB} = 5.98$ a₀ | 5.98 a₀ | ✓ |
| Grid points per λ_dB at dx=0.75 / 1.0 | 7.98 / 5.98 | (dx=0.75: 8) | ✓ |
| Nyquist headroom at dx=1.0 | $k_\text{max}=π=3.14 > k_0+3σ_k = 1.65$ | OK | ✓ |

**Numerical-validation script reproducing all of this:** see the bash
log at the head of this plan's Git history — every entry above was
re-derived from $L=50, N=162$ alone, no external numbers used.

## 1. Three-run programme

### Run A — "weak q-kick calibration" (deferred — needs custom INQ perturbation)

| Field | Value |
|---|---|
| Box / bath | L = 50 a₀, N = 162 (closed shell) |
| Perturbation | Dirac kick $V(r,t) = \varepsilon\,\cos(2\pi z/L)\,\delta(t)$ exciting only the $q = (0,0,2\pi/L)$ mode |
| Target mode | m = 1 |
| Expected peak in $n_q(t)$ FFT | **3.59 eV** |
| Duration | 2000 a.u. ≈ 48.4 fs |
| Status | **Not launched.** INQ's `perturbations::kick` is uniform (q → 0); a single-q perturbation requires extending `inq-stack/include/inqkit/...` with a new perturbation. Tracked as task #16. |

### Run B — primary WP plasmon run (this session)

| Field | Value |
|---|---|
| Box / bath | L = 50 a₀, N = 162 (closed shell, reuse `gs_L50_cubic_N162_dx1p0`) |
| WP | $\sigma_r = 5.0$ a₀, $k_0 = (0, 0, 1.05)$ a₀⁻¹, $E_\text{WP} = 15.0$ eV, occupation 1 |
| dt | 0.02 a.u. (within 0.02–0.05 range; safer choice) |
| N_steps | 10000 (smoke); **100000 (full)** ⇒ T_sim = 200 / 2000 a.u. |
| Total time | 200 / 2000 a.u. = 4.84 / 48.4 fs |
| WRITE_EVERY | 200 (gives 500 density frames at full) |
| Target mode | m = 1, $\lambda = 50$ a₀ |
| Expected ℏω | **3.59 eV** in n_q FFT |
| FFT resolution at T_sim=2000 | 0.0855 eV (5× the m=1 vs m=2 separation) |
| Wall budget | smoke ≈ 7 min (proportional to N=162 base run); full ≈ 12 hr |
| Run dir | `ResearchProject/systems/jellium/run_plasmon_n162_L50_E15/` |
| Status | **Launching now (smoke first, then full in background).** |

### Run C — secondary m=2 WP run (scaffolded only)

| Field | Value |
|---|---|
| WP | $\sigma_r = 4.0$ a₀, $k_0 = (0, 0, 0.585)$ a₀⁻¹, $E_\text{WP} = 4.65$ eV |
| Target mode | m = 2, $\lambda = 25$ a₀ |
| Expected ℏω | **4.00 eV** in $n_{q_2}$ FFT |
| Status | Cfg + run.cpp scaffold pending. Launch after Run B verifies the protocol. |

## 2. New observables required (not yet in the standard pipeline)

The current pipeline writes total density VTI per WRITE_EVERY plus
`observables.csv` (energy, current, dipole, cod). It does **not** write
the Fourier components $n_{\mathbf q}(t)$ that are the cleanest plasmon
diagnostic. For Run B's analysis we will:

1. **Post-hoc extract $n_{q_m}(t)$ from the saved density VTI series**
   for $m = 1, 2, 3$. Implementation: a new
   `inqview.postprocess.density_fourier` phase that reads
   `density_rt_total/density_t<step>.vti`, computes
   $\delta n(\mathbf r, t) = n(\mathbf r, t) - n(\mathbf r, 0)$, and
   does a 3D FFT to extract amplitudes at the integer-(n_x, n_y, n_z)
   grid points. Cheap (one 3D FFT per snapshot).
2. **Plot** $|n_{q_m}(t)|$ vs time and FFT($n_{q_m}(t)$) vs energy for
   each m, with the predicted plasmon frequency overlaid.

The Δn movie (already produced as `delta_yz.gif` etc.) is the
qualitative version of this.

## 3. Why a smoke test before the full run

The current jellium runs use `dt = 0.020`, `N_steps = 1500` — at
$k_0 = 1.05$ the propagator stability and orthogonalisation may fail
for reasons that don't apply at $k_0 = 0.332$. Specifically:

- **Nyquist headroom** at dx = 1.0 is $k_\text{max} = \pi = 3.14$ a₀⁻¹,
  comfortably above $k_0 + 3\sigma_k = 1.05 + 0.6 = 1.65$ a₀⁻¹ — but
  the WP envelope's tail has more weight at large k now.
- **WP-bath geometric overlap** at $k_0 = 1.05$ may be worse than at
  $k_0 = 0.332$ because the WP's plane-wave content overlaps more
  strongly with the highest occupied bath shells (which have
  $|\mathbf G|^2 \sim k_0^2$ ⇒ $|\mathbf G| \sim 1$). `max_overlap`
  could exceed 0.5 — flag for review.
- **dt = 0.02 a.u.** at $k_0 = 1.05$ propagates the WP by
  $\Delta x = v\,dt = 0.021$ a₀ per step ⇒ 0.021 / 1.0 = 2 % of
  one grid cell per step. Should be fine, but verify.

Smoke test = 100 steps. Pass criteria: `run_completed = true`,
`density_l2(0) = 0` exactly, `density_l2(t)` monotonically growing,
`max_overlap < 0.5`, `total_energy_drift < 1 mHa`.

## 4. Verdict criteria for Run B

A report under `docs/reports/plasmon-detection-verdict.md` (using the
report-writing skill) will document:

1. **YES, m=1 plasmon detected** if $|n_{q_1}(t)|$ shows a clean
   sinusoidal oscillation and FFT peak within ±0.1 eV of 3.59 eV.
2. **PARTIAL** if a peak appears in 2–5 eV but at a shifted location
   — investigate whether the Bohm-Gross dispersion under-predicts due
   to the actual finite-temperature smearing or LDA xc kernel.
3. **NO** if no peak is detected — implies plasmons in this box require
   the q-kick perturbation (Run A) and the WP setup, even at
   resonant velocity, doesn't couple cleanly. In that case escalate to
   Run A.

## 5. Order of execution this session

1. **Now**: write Cfg + run.cpp for Run B; smoke test 100 steps.
2. **If smoke passes**: launch the full 100 000-step run in
   background; ETA ~12 hours. The next session picks up the postprocess.
3. **Run A** (q-kick): defer until a custom `inqkit::perturbations::single_q_kick`
   class is added (task #16).
4. **Run C** (m=2 WP): launch after Run B verifies the protocol.
