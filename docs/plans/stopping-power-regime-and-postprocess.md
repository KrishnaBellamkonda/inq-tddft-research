# Plan: stopping-power regime diagram, S(v) trend, and dynamical post-process

**Audience:** the Claude Code agent driving the next few sessions, with the user as scientific reviewer. Output landing zone: `/local/data/public/skcb2/tddft/docs/reports/14-05-2026-meeting-emilio/` (figures, slide markdown, future PowerPoint).

**Scope:** this is the **first and most important thread** of the 2026-05-14 meeting with Emilio. The other threads (coronene LEED, Li-54 paper reproduction, plasmon-hunt) are documented in `docs/reports/presentation-2026-05-14.md` and are not part of this plan.

---

## Context

Why this exists. The 19-day campaign that ended on 2026-05-14 produced six classical+WP jellium runs at matched energies (E = 50, 100, 300, 600, 1500 eV) plus seven WP-only scoping runs. The data are there but the **story is not yet visualised** in a form that lets a viewer locate each run inside a textbook regime classification of electronic stopping. The meeting requires a coherent first-pass narrative that

1. shows *which textbook regime* each run lives in, with the boundary in (v/v_F, κ) parameter space drawn explicitly;
2. shows the *measured stopping power* across all runs, with classical and WP points overlaid on the regime map;
3. shows whether *the WP and classical results converge* in the classical limit — i.e. whether the correspondence principle holds in this set of TDDFT simulations.

A second post-processing phase (dynamical structure factor + Lindhard function) will let us read out the *physical mechanisms* underlying the qualitative features seen in the density evolution. That phase is rigorously test-driven before deployment.

A third phase (case study + microscopic decomposition) is sketched but deferred until the user reviews the Phase 1 outputs and selects a target run.

The intended outcome of *this plan* is: four publication-ready PNG figures and three production-ready post-processing modules, all under `docs/reports/14-05-2026-meeting-emilio/` and `inq-stack/python/inqview/postprocess/`.

---

## Theoretical grounding

The regime classification rests on **two dimensionless ratios** that exhaust the configuration space of stopping experiments for a fixed projectile species and a fixed jellium target:

- **v / v_F** — the host-response axis. Below 1: only e-h pairs at the Fermi surface excitable; Pauli blocking restricts phase space; stopping is *frictional* (S ∝ v). Above ~5: full loss-function exhausted; Bethe asymptote with S ∝ ln(v)/v². Between: Bragg peak, full RPA / Lindhard linear-response dielectric.

- **κ = 2|Z₁|/v** (a.u.) — the projectile-scattering axis. The Sommerfeld parameter, equal to the ratio of the classical Coulomb collision diameter b = 2|Z₁|e²/m_e v² to the de Broglie wavelength λ_dB = ℏ/m_e v. When κ ≪ 1: Born approximation valid (Bethe). When κ ≫ 1: classical impact-parameter scattering (Bohr). When κ ~ 1: Bloch interpolation needed.

These two parameters are sufficient because of dimensional analysis: for fixed Z₁, m_p, n the stopping power is a function only of v, and any two independent dimensionless ratios formed from (v, v_F, κ-scale = 2|Z₁|/v) span the regime space. Mass enters only through E_kin = m_p v² / 2; the projectile's relevant property is its charge and its speed. The framework is laid out in `docs/plans/jellium-regime-constrained-simulations.md` §3 and is consistent with Sigmund 2006 (Vol. 1, chapter 6 onwards) and the literature survey in `mental-models-for-jellium-and-tddft-analysis-claude.md` topic 1.

For an electron projectile (Z₁ = 1, m_p = m_e), the configuration space collapses to a single line κ · (v/v_F) = 2 / v_F = 5.93. **Every electron-projectile experiment in our jellium therefore lives on this 45° diagonal in the (v/v_F, κ) plane on log-log axes.**

### Regime equations (numerical curves drawn in the plots)

| Regime | Equation for S(v) | Boundary | Source |
|---|---|---|---|
| Friction | S = Q(r_s) · v with Q ≈ 0.5–1.5 a.u. for r_s = 5.69, Z₁ = 1 | v/v_F ≲ 0.5 | Fermi-Teller 1947; ENRA 1981/1986 |
| Bragg peak / Lindhard | S = (2 Z₁²/π v²) ∫ (dq/q) ∫₀^{qv} dω ω Im[−1/ε_RPA(q,ω)] | 1 ≲ v/v_F ≲ 5 | Lindhard 1954; Mahan ch. 5 |
| Bethe (pure) | S = (4π Z₁² n / v²) · ln(2 v² / ω_p) | v/v_F ≳ 5, κ ≪ 1 | Bethe 1930 |
| Bethe-Lindhard | S_Bethe with stopping number L = ln(2v²/ω_p) − 1/2 | v/v_F ≳ 5, κ ≪ 1 | Lindhard 1954 (equipartition) |
| Bohr classical | S = (4π Z₁² n / v²) · ln(C m_e v³ / |Z₁| ω_p), C = 2 e^{−γ} | κ ≫ 1 | Bohr 1948; Jackson ch. 13 |
| Bloch interp. | L_Bloch = ln(2v²/ω_p) − Re[ψ(1+iκ/2) − ψ(1)] | κ ~ 1 | Bloch 1933 |
| Nonlinear DFT / rt-TDDFT | full self-consistent response (no closed form) | **anywhere** the perturbative theories fail | Echenique-Nieminen-Ritchie; this work |

For r_s = 5.69 jellium: n = 1.296 × 10⁻³, k_F = v_F = 0.337 a.u., ω_p = 0.1276 Ha = 3.47 eV, E_F = 1.55 eV. These are the constants used in every curve in the plots.

---

## Phase 0 — four scientific figures (immediate, today)

All four PNGs land under `docs/reports/14-05-2026-meeting-emilio/figures/`. They are generated by the new module `inq-stack/python/inqview/postprocess/regime_diagram.py` and a small driver script `docs/reports/14-05-2026-meeting-emilio/build_figures.py` that the user can re-run on demand.

### Plot 1 — General regime diagram

**File:** `figures/01_regime_diagram_general.png`.

**Axes:** log(v/v_F) horizontal from 0.1 to 100; log(κ) vertical from 0.01 to 100.

**Coloured regions** (filled translucent rectangles or polygons with one-line labels):

| Region | (v/v_F, κ) range | Colour |
|---|---|---|
| Friction (Fermi-Teller / ENRA / nonlinear DFT) | v/v_F ≤ 0.5, any κ | red (S ∝ v) |
| Bragg peak (Lindhard / RPA, full ε(q,ω)) | 0.5 ≤ v/v_F ≤ 5, κ ≤ 5 | orange |
| Bohr classical | κ ≥ 5, v/v_F ≥ 0.5 | brown |
| Bloch interpolation | 0.5 ≤ κ ≤ 5, v/v_F ≥ 3 | purple |
| Bethe (perturbative + linear-response) | v/v_F ≥ 5, κ ≤ 0.5 | blue |
| Nonlinear rt-TDDFT (this work) | overlay annotation: "captures all of the above non-perturbatively" | thin black hatched contour |

A horizontal black line at κ = 1 and a vertical black line at v/v_F = 1 mark the canonical boundaries. The κ = 1 line carries the inline annotation "Born approximation onset". The v/v_F = 1 line carries "Pauli boundary".

**Reference:** Sigmund 2006 Vol. 1 Fig. 6.1 and `jellium-regime-constrained-simulations.md` §3.4 ASCII layout. Slight extension: include explicit "nonlinear DFT/TDDFT" annotation because that is the bridge to our work.

### Plot 2 — Same regime diagram for r_s = 5.69 jellium with run points

**File:** `figures/02_regime_diagram_jellium_rs5p69.png`.

Same axes, same coloured regions as Plot 1. **Additions:**

1. The electron-projectile constraint line κ · (v/v_F) = 5.93 drawn as a dashed black diagonal across the whole plot. Every electron-in-this-jellium experiment must live on this line.

2. Marked points along the line at the seven existing energies:

| E (eV) | v (a.u.) | v/v_F | κ | Existing run dir |
|---|---|---|---|---|
| 50 | 1.92 | 5.69 | 1.04 | `run_classical_n162_L50_E50_attempt2/`, `run_wp_n162_L50_E50_attempt2/` (running) |
| 100 | 2.71 | 8.04 | 0.74 | `run_classical_n162_L50_E100/`, `run_wp_n162_L50_E100/` |
| 300 | 4.69 | 13.9 | 0.43 | `run_classical_n162_L50_E300/`, `run_wp_n162_L50_E300/` |
| 600 | 6.64 | 19.7 | 0.30 | `run_classical_n162_L50_E600/`, `run_wp_n162_L50_E600/` |
| 1500 | 10.5 | 31.2 | 0.19 | `run_classical_e1500_L50_cubic/` (WP abandoned) |

Each point gets an inline label `"<E> eV"`. Classical-complete points are filled markers; WP-complete points are open markers stacked alongside (a tiny vertical offset for legibility). The 1500 eV WP point is drawn as a red "×" to mark it as not feasible at dx = 0.40 Bohr (see §1.2 of the prior plan — Nyquist over-aliased by 41 %).

3. Inline annotation: "All runs sit in the lower-Bethe / upper-Bragg-peak corner — none probe the friction regime or the κ > 1 quantum-scattering regime."

### Plot 3 — Measured S(v) across all jellium runs with regime overlay

**File:** `figures/03_stopping_power_measured.png`.

**Axes:** log E (eV) horizontal from 30 to 2000; log S (eV/Bohr) vertical from 0.005 to 1.

**Foreground regions** (lightly tinted vertical bands using the same colour scheme as Plot 1):

- v/v_F ≤ 0.5 → E ≤ 0.78 eV (off-plot, only a left arrow marker labelled "friction")
- 0.5 ≤ v/v_F ≤ 5 → 0.78 ≤ E ≤ 38.6 eV ("Bragg peak / Lindhard")
- v/v_F ≥ 5 → E ≥ 38.6 eV ("Bethe asymptote", split at κ = 1 into "Bloch interp." and "Bethe")

**Data points:**

Two distinct marker species:

- **Classical** (filled circle, blue) — computed by the bath-energy method with proper windowing per `jellium-regime-constrained-simulations.md` §6.1. Read from the existing `results/analysis/REPORT.md` of each completed classical run; for runs without a REPORT, compute S = ΔE_bath / Δz inside the windowed range using `bath_energy_vs_time.csv`. For the 50 eV pair which is still running, use the latest in-flight `bath_energy_vs_time.csv` if Δz ≥ 10 Bohr has been reached, otherwise omit and note as "running" in the legend.

- **Wave packet** (open square, red) — computed by S_WP,bath ≡ ΔE_bath / Δ⟨z_proj⟩ on the WP run's observables.csv. Where `⟨z_proj⟩(t)` is not yet extracted into a CSV, fall back to the **crude metric** the user specified: total bath energy gain over total projectile centroid displacement during the run (read centroid from `observables.csv` columns `dipole_z` if `N_IONS=0` — for a single-orbital WP the dipole is the centroid times the WP charge).

Each point gets an inline label "E eV" near the marker; error bars come from the regression covariance for classical points and from a simple finite-difference estimate for WP points.

**Note on the existing Phase-0 Plot 3 data:** the values to plot for the three completed classical runs are already in their respective `results/analysis/REPORT.md`:

| Run | v (Bohr/a.u.) | KE loss (eV) | S(v) (eV/Bohr) | Source |
|---|---|---|---|---|
| `run_classical_n162_L50_E100` | 2.711 | 12.3703 | 0.3629 | `results/analysis/REPORT.md` line "Stopping power S(v)" |
| `run_classical_n162_L50_E600` | 6.640 | 2.1243 | 0.0606 | `results/analysis/REPORT.md` |
| `run_classical_e1500_L50_cubic` | 10.500 | 0.9371 | 0.02076 | `results/analysis/REPORT.md` |

The WP S values require computing the centroid trajectory from `observables.csv`; this is done in the same plot script.

### Plot 4 — Plot 3 with analytical theory curves overlaid

**File:** `figures/04_stopping_power_measured_vs_theory.png`.

Same as Plot 3 but with three analytical curves drawn through the regime regions:

1. **Bethe-pure** (solid line) — S = (4π n / v²) · ln(2 v² / ω_p) over the range where v/v_F ≥ 5.
2. **Bethe-Lindhard** (dashed line) — Bethe with stopping number L = ln(2v²/ω_p) − 1/2. The "equipartition" correction from Lindhard 1954.
3. **Bloch-corrected** (dotted line) — Bethe-Lindhard with the additional Re[ψ(1+iκ/2) − ψ(1)] subtraction. Shows the small turn-over near κ ~ 1.

Optional (only if compute-cheap): a **Lindhard-RPA** dotted line for v/v_F ≲ 5 obtained from numerical integration of the full Lindhard ε_RPA(q,ω). This would let the figure span both the Bragg-peak region and the Bethe asymptote with a single analytical curve. Implementation lives in the same Lindhard module that Phase 1 will produce.

The legend specifies each curve with its applicable κ / (v/v_F) condition.

### Acceptance criteria for Phase 0

- Four PNGs at 300 dpi, transparent or white background, ~12 cm wide (slide-friendly).
- Axes labels in TeX form ($v/v_F$, $\kappa$, $E$ [eV], $S(v)$ [eV/Bohr]).
- Regime labels readable at slide scale.
- All numeric values traceable to a source file path written into the script's docstring.
- Plot 3 and 4 use the same axes and legend (Plot 4 is the same figure with curves enabled — a flag in the script).

---

## Phase 1 — dynamical structure factor and Lindhard function (next session)

Two new postprocess modules. Both follow the project's `development-feedback-loop` rule: every function gets a known-case test before being applied to a production run.

### 1.1 Module `inq-stack/python/inqview/postprocess/lindhard.py`

**Purpose.** Closed-form analytical χ⁰(q, ω) for a non-interacting electron gas at r_s = 5.69, used as the reference against which the rt-TDDFT density response is benchmarked.

**Key function:** `chi0_lindhard(q, omega, kF, eta=1e-3) -> complex`. Implements the Lindhard susceptibility from chapter 8 (textbook excerpt in `ResearchProject/literature/chapter_7.pdf` Eq. 8.18) and Sigmund 2006 Vol. 1 chapter 9. Returns a complex array shaped like `q × omega`.

**Derived quantities:**

- `loss_function(q, omega)` = Im[−1/ε_RPA(q, ω)] where ε_RPA = 1 − v(q) χ⁰. The integrand of the analytical Lindhard stopping-power formula.
- `plasmon_dispersion(q)` — Bohm-Gross to leading order, full RPA root-finding for higher accuracy. Used as the analytical reference for the m=1 axial plasmon observed in `run_plasmon_n162_L50_E15` (3.533 eV measured vs Bohm-Gross 3.59 eV expected).
- `stopping_power_lindhard(v, qmin, qmax)` — numerical integration of S(v) = (2 Z²/π v²) ∫ (dq/q) ∫₀^{qv} dω ω Im[−1/ε(q,ω)]. With `qmin = 2π/L = 0.126 a.u.` we directly compute the box-truncated Lindhard prediction that should match our measurements better than the textbook qmin → 0 form.

**Known-case tests** (in `inq-stack/python/inqview/postprocess/test_lindhard.py`):

1. Static limit: χ⁰(q, 0) at q → 0 must equal −D(E_F) = −e D(E_F) with D(E_F) = 3 n / 2 E_F = the textbook Thomas-Fermi value (chapter 8 Eq. 8.2). Numerical tolerance 1 %.
2. f-sum rule: ∫₀^∞ ω Im[−1/ε(q,ω)] dω = π ω_p² / 2 for any fixed q (chapter 8 §8.3.2; ω_p² = 4π n in a.u.). Numerical tolerance 1 %.
3. High-ω limit: χ⁰(q, ω → ∞) → −n q² / m ω² (chapter 8 Eq. 8.23). Numerical tolerance 0.1 %.
4. Plasmon dispersion: ω_pl(q → 0) → ω_p = 0.1276 Ha. Numerical match to 3 decimals.
5. Pure-Bethe limit: stopping_power_lindhard at v = 10.5 a.u. with qmin → 0 must match the Bethe-Lindhard prediction (0.0279 eV/Bohr per `jellium-regime-constrained-simulations.md` §4.3 table) to 1 %.

**Acceptance:** all five tests pass. Then the module is wired into Plot 4.

### 1.2 Module `inq-stack/python/inqview/postprocess/dynamical_structure_factor.py`

**Purpose.** Compute S(q, ω) directly from the simulation's time-dependent density n(r, t). S(q, ω) is the Fourier transform of the density-density correlation function and is — by the fluctuation-dissipation theorem — proportional to Im χ(q, ω). It is the *measured* analogue of the Lindhard loss function.

**Key function:** `compute_dsf(density_series, dt, box, omega_grid, q_grid) -> array`. Two-step pipeline:

1. **Spatial FFT** of δn(r, t) = n(r, t) − n(r, 0) at each saved frame, giving δñ(q, t). Use `numpy.fft.fftn` with the centred-Cartesian convention used by inqview's VTI series.
2. **Temporal FFT** of δñ(q, t) over the clean window (same windowing rule as the stopping-power extraction: drop the first ~5 Bohr of projectile travel as a transient).
3. Return S(q, ω) = |δñ(q, ω)|² / N (with N a normalisation that maps to the textbook S(q,ω) via the fluctuation-dissipation relation).

**Known-case tests** (in `test_dynamical_structure_factor.py`):

1. **f-sum rule.** ∫₀^∞ ω S(q, ω) dω = (n q² / 2 m) · N for every q. Numerical tolerance 5 % (worse than for Lindhard because of windowing and finite-T effects).
2. **Plasmon at the resonance run.** Apply to `run_plasmon_n162_L50_E15/results/raw/vti/density_delta/`. The output S(q, ω) at q = 2π/L · ẑ must have a peak at ω = 3.533 eV — the m=1 plasmon already detected and quoted in the journal entry. Pass = peak found within ±0.05 eV.
3. **No-perturbation control.** Apply to a 100-frame fake series with constant density (no perturbation). S(q, ω) must be uniformly zero to machine precision.
4. **Single-mode injection.** Inject a synthetic time-dependent density n(r, t) = n_0 + A cos(q₀ · r − ω₀ t). S(q, ω) must show a single peak at (q₀, ω₀) of expected amplitude.

**Acceptance:** all four tests pass. Test 2 is the most diagnostic: it tells us the production-data pipeline reproduces the *known* plasmon at the *known* energy from existing simulation output.

### 1.3 Where to apply S(q, ω) once validated

The dynamical structure factor is **most meaningful** in runs where:

1. The simulation is long enough to give frequency resolution Δω = 2π / t_total ≤ 0.5 eV. → runs with N_STEPS ≥ 1500 at dt = 0.020 a.u. **The plasmon-hunt runs (`run_plasmon_n162_L50_E{15, 3p4_varyv}`) and the 1.5 eV scoping runs satisfy this** (N_STEPS = 100 000 and 1500 respectively).

2. The bath is in a well-defined regime (Bragg peak or plasmon-resonance, not deep Bethe where the projectile is too fast for the bath to respond coherently).

**Primary targets** for the production S(q, ω) deployment:

- `run_plasmon_n162_L50_E15` — already known to host the m=1 plasmon; serves as the **calibration** of the postprocess.
- `run_plasmon_n162_L50_E3p4_varyv` — the velocity-discriminator run; S(q, ω) at the two velocities should give the same plasmon ω(q) and identify the kinematic peak as a separate, velocity-dependent feature.
- `run_classical_n162_L50_E100` and `run_wp_n162_L50_E100` — matched pair at the lower-Bethe / Bragg-peak boundary. S(q, ω) here lets us check whether the classical and WP projectiles excite the *same* (q, ω) spectrum of the bath.
- `run_base_n162_L50_E1p5` — the canonical sub-threshold reference; the density hole behind the WP is expected to show up at small q, low ω (quasi-static screening).

**Secondary targets** (lower priority, run only if S(q, ω) on the primary targets gives clean signals):

- `run_classical_n162_L50_E600` and `run_wp_n162_L50_E600` — deep into the Bethe regime. The DSF should be dominated by single-particle e-h pair production with no plasmon peak.
- `run_classical_e1500_L50_cubic` — high-velocity classical reference; DSF should show only the Bethe-ridge-like band of e-h excitations.

The output: for each primary target, one figure `figures/05_dsf_<run>.png` showing |S(q, ω)| as a heatmap in (q, ω) with the Bohm-Gross plasmon dispersion ω_pl(q) and the Bethe ridge ω = q²/2 + q · k₀ overlaid.

### 1.4 Connection to the correspondence-principle question

Comparing classical S(q, ω) and WP S(q, ω) at matched energy is the **direct experimental test of the correspondence principle in this campaign**. The WP at high k₀σ (deep in the classical-packet limit per `jellium-regime-constrained-simulations.md` §1 Q2: k₀σ ranges 9.6 → 33.2 across our energies) should produce the same bath response as the classical projectile, modulo the host-side asymmetries:

- exchange / indistinguishability of the WP electron with the bath electrons;
- finite spatial extent of the source (σ = 5 Bohr Gaussian vs point charge).

If the two DSFs match within a few percent across the (q, ω) plane → correspondence principle confirmed for our setup. If they differ systematically (especially at small q where the Gaussian cutoff matters) → we can quantify the **non-classical contribution** from the WP description, which is itself a publishable result.

---

## Phase 2 — Bethe-regime metric study (added 2026-05-14)

**Question** (the user's central ask): *Given that closed-system WP TDDFT does not have a single obvious analog of the classical "projectile KE loss", which observable best captures the stopping power for the WP in the Bethe regime?*

### Case-study target

`run_classical_n162_L50_E100` ↔ `run_wp_n162_L50_E100` at v/v_F = 8.04, κ = 0.74 — lower-Bethe regime at the edge of the Bloch interpolation band. Both runs complete.

**References for the "right" S(v) at this point:**

- Classical measured (REPORT.md): **S = 0.363 eV/Bohr**
- Bethe (pure log) prediction: **S = 0.286 eV/Bohr** at v = 2.711 a.u.
- Classical via "bath state-sum" method applied to the classical run: **S = 0.347 eV/Bohr** (matches REPORT.md to within 5%, confirming the metric works for classical)

### Six candidate metrics — implementation in `docs/reports/14-05-2026-meeting-emilio/metric_comparison.py`

| # | Metric | Formula | Observable source |
|---|---|---|---|
| 1 | Bath state-sum (current) | ΔE_bath(t) / Δz | state_energies.csv minus WP slot |
| 2 | WP slot loss | \|ΔE_WP_slot(t)\| / Δz | state_energies.csv WP slot |
| 3 | System kinetic gain | \|ΔE_kinetic(t)\| / Δz | observables.csv energy_kinetic |
| 4 | Hartree anti-wake | \|ΔE_Hartree(t)\| / Δz | observables.csv energy_hartree |
| 5 | WP KE loss | (½ k₀² − ½ v_f²) / Δz, v_f from density_delta positive-centroid fit | VTI snapshots |
| 6 | Closed-system sanity | \|ΔE_total\| / Δz | observables.csv (should be ≈ 0) |

All metrics use Δz = v₀ · t (analytic). Window: Δz ∈ [3, 18] Bohr (clean window per §1.1 / `audit_wp_methods.py`).

### Measured values at E = 100 eV (run on 2026-05-14)

| Metric | S (eV/Bohr) | Ratio to classical 0.363 |
|---|---|---|
| 1 — bath state-sum | **0.0209** | ÷ 17 |
| 2 — \|WP slot\| | **0.0436** | ÷ 8 |
| 3 — system kinetic | **0.0208** | ÷ 17 |
| 4 — \|Hartree anti-wake\| | **0.0236** | ÷ 15 |
| 5 — WP KE loss (VTI) | **~2.5** | × 7 (overshoots; wake-bias) |
| 6 — closed-system sanity | **~0** | ÷ ∞ (passes as null) |

**Findings**:

1. Metrics 1–4 cluster at 0.02–0.04 eV/Bohr — about **10× below** the classical reference. None of them captures the missing energy that lives in Hartree+xc rearrangement.
2. Metric 5 (WP velocity-decay from VTI) gives **2.5 eV/Bohr** — overshoots by 7× because the density_delta positive-centroid is biased by the *co-moving polarization wake* that drifts forward with the WP. The wake adds an effective forward shift that mimics extra WP velocity decay when fit naively.
3. **No simple metric on the existing observables reproduces the classical Bethe stopping power for the WP.** The correct extraction of WP velocity requires either:
   - A bath-subtracted density (density_total − density_GS_bath), giving the WP density alone for centroid extraction.
   - A separate per-orbital ⟨ψ_WP | T̂_kinetic | ψ_WP⟩ observable, which is not currently saved.
   - Equivalently, the WP's centroid trajectory tracked outside the system-centroid contamination.

### Recommended next step (Phase 2 follow-up)

Add a new observable to the simulation pipeline (`inq-stack/include/inqkit/`): per-orbital **kinetic-energy expectation** ⟨ψ_i | T̂ | ψ_i⟩ written to `state_kinetics.csv` at every WRITE_EVERY step. Then Metric 5 becomes a direct lookup: −Δ⟨T̂⟩_WP / Δz.

This requires a small inqkit observable extension (≈ 50 lines C++) and rerunning the matched pair. With that observable, the WP stopping power could be measured to ~1% precision and the correspondence-principle question resolved cleanly.

### Outputs

- `figures/08_metric_comparison_E100eV.png` — bar chart with all six metrics + reference lines
- `figures/08_metric_comparison_E100eV.csv` — numerical table

### Original deferred case-study (still pending)

A deeper density-evolution case study (snapshots at t = 0.5, 1, 2, 3 a.u.; wake comparison; physical origins) remains pending. Likely target: the same E=100 eV matched pair, after the per-orbital ⟨T̂⟩ observable lands. Deliverable: `docs/reports/14-05-2026-meeting-emilio/case-study-E100eV.md`.

---

## Phase 3 — deferred microscopic diagnostics (remember, do not implement now)

Per the user's explicit "not to be implemented just yet" note in the 2026-05-14 conversation. These are recorded here so the next agent doing implementation work knows the queue exists.

| Diagnostic | What it computes | Why it matters | Module home |
|---|---|---|---|
| z-profile δn(z, t) GIFs | density change vs z along the projectile axis, animated | Visualises the wake formation and the orthogonalisation hole | `inqview/postprocess/density_z_profile.py` |
| KS eigenenergy time evolution | ε_i(t) = ⟨φ_i(t) \| H_KS(t) \| φ_i(t)⟩ | Excitation/de-excitation of individual KS states | `inqview/postprocess/state_energy_spectra.py` (already exists, may need extension) |
| GS-basis decomposition of evolved KS orbitals | c_ij(t) = ⟨φ_i^GS \| φ_j(t)⟩; effective occupations n_i^eff(t) = Σ_j f_j \|c_ij(t)\|² ; if clean, transition rates from time-derivative | Effective excitation spectrum of the bath; turns rt-TDDFT into linear-response language | new `inqview/postprocess/gs_basis_decomposition.py` |
| WP centroid and centre-of-density | ⟨z⟩(t) and z_cod = ∫ z n(r,t) d³r / N | Quantifies the projectile slowdown and the bath polarisation co-moving with it | extension to `wavepacket_observables.py` |
| Local density fluctuations | σ_n(r, t) = √⟨(n − ⟨n⟩)²⟩ over a small spatial window | Spatially resolved measure of bath stirring | new `inqview/postprocess/local_density_fluctuations.py` |

These get prioritised against ongoing simulation work by the user when the time comes.

---

## Files and paths

| Artefact | Path |
|---|---|
| This plan | `/local/data/public/skcb2/tddft/docs/plans/stopping-power-regime-and-postprocess.md` |
| Figures (Phase 0) | `/local/data/public/skcb2/tddft/docs/reports/14-05-2026-meeting-emilio/figures/{01..04}_*.png` |
| Future PowerPoint | `/local/data/public/skcb2/tddft/docs/reports/14-05-2026-meeting-emilio/<deck-name>.pptx` (Phase 2 / post-review) |
| Figure-builder | `/local/data/public/skcb2/tddft/docs/reports/14-05-2026-meeting-emilio/build_figures.py` |
| Regime-diagram module | `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/regime_diagram.py` |
| Lindhard module (Phase 1) | `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/lindhard.py` |
| Lindhard tests | `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/test_lindhard.py` |
| DSF module (Phase 1) | `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/dynamical_structure_factor.py` |
| DSF tests | `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/test_dynamical_structure_factor.py` |
| Case-study report (Phase 2) | `/local/data/public/skcb2/tddft/docs/reports/14-05-2026-meeting-emilio/case-study-<run>.md` |

## Existing code to reuse (do not re-implement)

- `inq-stack/python/inqview/postprocess/observables.py` — observables.csv readers.
- `inq-stack/python/inqview/postprocess/occupations.py` — KS occupation extraction.
- `inq-stack/python/inqview/postprocess/state_energy_spectra.py` — ε_i(t) extraction.
- `inq-stack/python/inqview/postprocess/density_fourier.py` — the spatial-FFT-of-density pipeline already used in the plasmon-hunt journal entries; DSF module extends this.
- `inq-stack/python/inqview/fields.py`, `data.py`, `fourier.py` — generic field IO and FFT helpers.

## Trusted sources

- Sigmund, *Particle Penetration and Radiation Effects*, Vol. 1, Springer Series in Solid-State Sciences 151 (2006) — fetched from Google Drive `1rPUoaJEaLzSiKdDdAk3XI8fTVM57VlCn`. Cite for every regime boundary.
- `mental-models-for-jellium-and-tddft-analysis-claude.md` (Google Drive `1B6PJtveryGmPXc9lcUTS1fzM2FO74K7O`) — literature-survey companion; section 1 grounds the regime classification; section 7 grounds the rt-TDDFT methodology.
- `ResearchProject/literature/chapter_7.pdf` — many-body textbook chapters 7 (interacting fermions, Hartree-Fock) and 8 (Thomas-Fermi, Kubo, Lindhard, RPA, plasmons); ground the closed-form χ⁰ used in `lindhard.py`.
- `docs/plans/jellium-regime-constrained-simulations.md` — the prior plan; contains the canonical regime classification table, the per-regime stopping-power formulas, and the operational definitions of S_classical and S_WP.
- `docs/reports/classical-vs-wp-case-study.md` — the existing case-study report; contains end-of-trajectory comparisons for the E = 100 and E = 1500 eV pairs.
- Per-run `results/analysis/REPORT.md` files for the numerical S(v) values.

## Verification

1. **Phase 0 figure-quality check.** Open each PNG; confirm regime boundaries visible at slide scale; confirm legend includes each curve; confirm every data point has an inline E-label.
2. **Phase 0 number cross-check.** For Plot 3, the three completed classical S(v) values must match the per-run REPORT.md verbatim (0.3629, 0.0606, 0.02076 eV/Bohr at E = 100, 600, 1500 eV).
3. **Phase 1 test suite.** `pytest inq-stack/python/inqview/postprocess/test_lindhard.py inq-stack/python/inqview/postprocess/test_dynamical_structure_factor.py` — all 9 tests pass.
4. **Phase 1 calibration on a known plasmon.** S(q, ω) module applied to `run_plasmon_n162_L50_E15` must reproduce the m=1 peak at 3.53 eV ± 0.05 eV.
5. **Phase 2 readiness gate.** User reviews Phase 0 + Phase 1 outputs, selects the case-study target, then case-study work begins.

## Out of scope for this plan

- The remaining three threads (coronene, plasmon, Li-54) — covered in `docs/reports/presentation-2026-05-14.md`.
- Slide-deck assembly — happens in a follow-up session after the user picks the storyline.
- New simulations — none launched in this plan; data drawn entirely from completed runs.
- Phase 3 microscopic diagnostics — recorded but deferred.
