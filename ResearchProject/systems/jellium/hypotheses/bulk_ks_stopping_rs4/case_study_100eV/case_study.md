# Case study — 100 eV projectile in high-density bulk jellium

**Run pair:** `bulk_ks_stopping_rs4` (r_s = 3.99 Bohr, σ_WP = 2 Bohr, 100 eV, ALDA,
dt = 0.04 a.u., 646 steps, 40 × 40 × 80 Bohr, 482 electrons, ω_p = 5.9 eV)

Two runs identical in every physical parameter, differing only in how the
projectile is represented — a classical point charge with a Gaussian UPF, versus
an occupied Kohn–Sham orbital. Everything below asks one question: **where does
the factor of ~5.6 in the stopping power come from?**

Regenerate with:

```bash
venv/bin/python ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping_rs4/case_study_100eV/make_case_study.py
```

---

## 1. Choosing what to compare

The two ledgers are not the same object, so "the projectile's kinetic energy"
has to be picked deliberately:

| symbol | definition | meaning |
|---|---|---|
| T_cl | ½ m v_z² | classical kinetic energy |
| **T₂** | **⟨p⟩²/2m** | **drift — the classical analogue** |
| T₁ | ⟨p²⟩/2m | total KS-orbital kinetic energy |
| T_var | var(p)/2m = T₁ − T₂ | momentum spread — quantum only |

**T₂ is the comparable quantity.** T₁ additionally contains the packet's internal
momentum spread, which at t = 0 is the pure zero-point value 3/(4σ_ψ²) = **5.10 eV**
at σ_ψ = 2 Bohr — energy the packet carries before it has moved. Charging that to
stopping would be wrong, and the data confirms the identity to 6 digits.

T_var is a channel the classical projectile does not possess: the packet can shed
drift energy into its own spreading without any of it reaching the bath.

**Position.** Classical z(t) is read directly. The wavepacket centroid is obtained
by integrating the mean momentum, z_c(t) = z₀ + (1/m)∫⟨p_z⟩dt′ — immune to the
periodic wrap at the cell face, and licensed by Ehrenfest's theorem. The density
centroid ⟨z⟩ is plotted alongside as an independent check: the two agree to
**0.079 Bohr over 68.5 Bohr** of travel (figure 06, lower panel).

---

## 2. Results

S = −dT/ds by OLS over the run's own transient-excluded window,
**t = 4.0 – 18.4 a.u.**, read from `wp_config.txt` rather than retyped.

| quantity | S (eV/Bohr) | r² |
|---|---|---|
| classical, ½mv² | **0.88 ± 0.13** | 0.991 |
| wavepacket, T₂ (drift) | **0.16 ± 0.01** | 0.997 |
| wavepacket, T₁ (total) | 0.12 ± 0.02 | 0.993 |
| wavepacket, T_var (spread) | −0.04 ± 0.02 | 0.902 |

**Ratio classical / wavepacket(T₂) = 5.6 ± 0.9.**

The three wavepacket rows are additive by construction (S(T₁) = S(T₂) + S(T_var)),
which the test suite asserts — a mismatch would mean the fits ran over different
windows or path coordinates. The negative T_var row means the spread is *growing*.

Uncertainty is the OLS slope standard error and a window-sensitivity systematic
(both edges moved ±3 a.u.) in quadrature. **The systematic dominates entirely** —
the statistical error rounds to 0.00 in every row. That is the honest ordering:
where the transient ends is a judgement call, not a noise problem.

### Energy budget over the full run (t = 0 → 25.84 a.u.)

| | start | end | change |
|---|---|---|---|
| classical ½mv² | 100.0 | 45.6 | **−54.4 eV** |
| wavepacket T₂ | 100.0 | 91.3 | **−8.7 eV** |
| wavepacket T₁ | 105.1 | 99.4 | −5.7 eV |
| wavepacket T_var | 5.10 | 8.09 | **+3.0 eV** |

Of the packet's 8.7 eV of drift loss, **3.0 eV goes into its own momentum spread**
and never reaches the bath. Only 5.7 eV actually leaves the projectile.

### Interaction energies — all six pairwise terms (change over the run)

| term | classical | wavepacket | WP − cl |
|---|---|---|---|
| E_SS (bath–bath, the wake) | +26.1 | +5.3 | **−20.8 eV** |
| E_PP (projectile self-Hartree) | 0.000 | −4.6 | −4.6 eV |
| E_PS (projectile–bath) | −5.5 | −0.4 | +5.1 eV |
| E_SB (bath–background) | 0 | 0 | 0 |
| E_PB (projectile–background) | 0 | 0 | 0 |
| E_BB (background self) | 0 | 0 | 0 |

The last three are **bitwise exactly zero at every step in both halves** — not
small, zero. Bulk has a uniform background, so poisson(n₊) is pure G = 0, which
INQ drops, making φ₊ ≡ 0. The columns are written anyway so the schema matches
the slab systems. Figures 14 and 15 plot them on the same axis as the terms that
do move, with a lower panel magnified 50× that shows them still flat — which is
what distinguishes "structurally zero" from "we forgot to compute them". Asserted
bitwise in the test suite.

**E_PP(0) is identical for both halves** (0.176996 Ha). The classical Gaussian UPF
is generated at σ_pot = σ_WP/√2 precisely so its cloud carries the packet's t = 0
density — so this equality is a *check on* the σ-matching convention, not an input
to it. It is now asserted in the test suite. Thereafter the classical E_PP is
exactly constant (a rigid cloud cannot spread) while the packet's collapses by
4.6 eV as it disperses.

---

## 3. What this says about the factor of ~5.6

Reported as inference, not as a settled result.

**E_PP alone does not account for the gap.** The classical-minus-wavepacket drift
loss is 54.4 − 8.7 = 45.7 eV over the run; the packet's E_PP releases 4.6 eV. That
is ~10% of the gap, so the projectile self-Hartree is not on its own the missing
factor that the campaign hypothesis proposed. The much larger single term is the
**wake**: the classical projectile drives 26.1 eV into E_SS against the packet's
5.3 eV, a 4.9× ratio that tracks the 5.6× stopping ratio closely.

**Inference:** the packet stops less because it is spatially extended and disperses
— it couples to the bath more weakly and drives a correspondingly weaker wake —
rather than because its self-interaction is doing hidden work. E_PP is real and
measurable, but on these numbers it is a minority channel.

**Caveats.**
- No full ledger closure was attempted here. E_SS, E_PP, E_PS do not exhaust where
  the energy goes (bath kinetic energy and XC are not in this decomposition), so
  the 20.8 eV wake gap and the 45.7 eV drift gap should not be differenced naively.
- Absolute E_PP carries the charged-cell G = 0 gauge. Only the WP-minus-classical
  differences quoted above are gauge-clean.
- Bulk jellium has a uniform background, so φ₊ ≡ 0 and E_SB = E_PB = E_BB = 0 by
  construction — the columns are written as zeros so the schema matches the slab
  systems.
- Single pair, single density, single σ. The σ = 3 and r_s = 5.70 pairs are needed
  before any of this generalises; `--family` takes them today.

---

## 4. Figure index

All 600 dpi, canonical theme (ADR 0004), fixed canvas, time axes in fs with the
fit window annotated in a.u.

| file | shows |
|---|---|
| `01_classical_KE_vs_time.png` | classical ½mv² vs t |
| `02_wp_T1_T2_vs_time.png` | T₁ and T₂ vs t |
| `03_wp_Tvar_vs_time.png` | var(p)/2m vs t, with the 5.10 eV zero-point line |
| `04_wp_kinetic_decomposition.png` | all three + classical, stacked |
| `05_classical_z_vs_time.png` | classical z(t) |
| `06_wp_centroid_vs_time.png` | Ehrenfest centroid + density centroid, residual panel |
| `07_position_vs_time_both.png` | both trajectories |
| `08_wp_T1_T2_vs_position.png` | T₁, T₂ vs centroid position |
| `09_T_vs_position_both.png` | all three vs position |
| `10_fit_classical_T_vs_position.png` | + fit window, OLS line, ±band, S |
| `11_fit_wp_T1_T2_vs_position.png` | + fit window, OLS lines, ±band, S |
| `12_fit_classical_KE_vs_time.png` | same fit shown on the time axis |
| `13_fit_wp_kinetic_vs_time.png` | same fit shown on the time axis |
| `14_interactions_classical.png` | **all six** terms vs t + magnified background panel (two-column width) |
| `15_interactions_wp.png` | **all six** terms vs t + magnified background panel (two-column width) |
| `16_interactions_difference.png` | **all six** differences, WP − classical (two-column width) |
| `stopping_power.txt` | all values, uncertainties, provenance, validation |

---

## 5. Validation

| check | result |
|---|---|
| Hartree closure vs INQ, classical | 4.99e-13 Ha |
| Hartree closure vs INQ, wavepacket | 4.98e-13 Ha |
| Ehrenfest residual, ⟨z⟩ vs ∫⟨p⟩dt | 0.079 Bohr over 68.5 Bohr |
| WP orbital norm | 0.999989595 – 1.000000000 |
| classical cloud clipping | never (window ends 18.4 a.u.) |
| S additivity, S(T₁) = S(T₂) + S(T_var) | holds to 1e-10 |
| figure margin check | all 16 clear (≥ 8 px) |
| test suite | 21 new, 123 total, all pass |

Tests: `ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/tests/test_case_study.py`

---

## 6. Production notes

Kept because two of these cost real time.

**Fixed-canvas figures crop silently.** The standard forbids `bbox_inches="tight"`
so on-page font size is exact — but that means a title wider or taller than its
reserved strip is simply *absent* from the PNG, with the build reporting success.
First cut had the two-line title flush against row 0 on 15 of 16 figures. Fixes:
reserve the full 0.42 in the title actually occupies; auto-shrink any title that
overruns horizontally; and `verify_margins()` now measures the ink bounding box of
every saved PNG and fails the build. A stacked panel's signed decimal ticks
(`−0.10`) later pushed a y-label off the left edge — caught by that same check.

**Don't cover data to label it.** The theme's frameless legend let curves run
through the label text; backing it with white then *hid* the curve. The fix is to
expand the axis limit on the legend's side so the legend sits in genuinely empty
space. Same reasoning retired the inset in figure 06 in favour of a second panel.

**Mathtext is not TeX.** `$\frac12 mv^2$` is valid Python and reads as valid TeX
but raises in matplotlib. No blocklist finds this — only the real parser does. All
34 math spans in the builder are parser-checked in CI, with a negative self-test
proving the check is not vacuous.
