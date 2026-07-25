---
name: stopping-power-extraction
description: Use when extracting electronic stopping power S from a classical-projectile rt-TDDFT/Ehrenfest run. Picks the run-geometry-appropriate method — a fixed 20%-of-time transient cut + slope fit for continuous-traversal (bulk) runs, or ΔE_total/L_z (slab thickness) for localised slabs — and treats the kinetic / ∫F·v / equal-potential-face channels as sanity checks only.
---

# Stopping-power extraction (classical projectile)

Electronic stopping power **S = energy the projectile deposits in the electrons
per unit path length**. Grounded in **Correa 2018** (*Comput. Mater. Sci.* 150,
291): Eq.(10) `S = ⟨dE/dt⟩/v` ≡ `dE/dx`; transient-vs-steady-state and the rule
to discard the transient is Fig.8 (`docs/sources/correa-2018-electronic-stopping-power.md`).

All kernels live skill-locally in **`stopping_power.py`** (this folder). Load run
data yourself, pass numpy arrays. Keep `E` and `x` in one consistent unit system
(Ha & Bohr → Ha/Bohr, or eV & Bohr → eV/Bohr).

## 0 — Before anything: the guards

1. **N(t) ≈ const** (`conservation_guard`). If a CAP drains the bath (N falls),
   raw `energy_total` is dominated by CAP energy injection, **not** the deposited
   energy — *both* methods are invalid. (Hard-won: bulk `b2_classical_E100` drained
   97% of N; its `ΔE_total` "peaks" at +20 Ha ≫ the projectile's 3.7 Ha KE — garbage.)
2. **Energy conservation:** `ΔE_total ≈ ΔKE_ion` (the projectile's KE loss). They
   come from independent channels (electronic `observables.csv` vs the classical
   track); agreement to a few % means the signal is trustworthy.

## 1 — Choose the method by run geometry (locked default, user 2026-06-25)

One deterministic branch on geometry. **Pick exactly one PRIMARY; the other
channels are sanity checks only — never average them into the answer.**

| run geometry | PRIMARY method | formula |
|---|---|---|
| **Localised slab** — finite slab, projectile *enters → deposits → exits* | **B: deposit / thickness** | `S = [E_total(t_final) − E_total(t0)] / L_z` |
| **Continuous traversal** — bulk jellium, projectile ploughs through a homogeneous medium | **A: fixed 20%-of-time cut + slope fit** | discard first **20% of the simulation time**, then `S = slope of ΔE_total(x)` (free-intercept) |

The stricter agent slope-plateau detector (`detect_x0_and_stopping_power`) is now an
**optional diagnostic**, not the default — use it only to *probe* whether a bulk run
reached steady state, not to produce the headline number.

## 2 — Method A (bulk default): fixed 20%-of-time cut + slope fit

**The locked default for every continuous-traversal run.** One deterministic rule,
no tuning:

- **Signal:** `ΔE_total(x) = E_total(t) − E_total(t0)`; `x` = projectile
  displacement (`s = z − z0`) from the track, interpolated to the energy-sample times.
- **Transient cut — fixed 20% of the SIMULATION TIME** (`fixed_time_fraction(t, x, E,
  frac=0.20)`): discard the first 20% of `t`, set `x0` = the displacement reached at
  that time, fit the remainder. (For a ~constant-velocity bulk run this ≈ 20% of the
  path, but the rule is on time.) `fixed_fraction_window` is the same idea cut in `x`
  if you have no time axis.
- **Fit:** **free-intercept** `ΔE = S·x + E0` (both free). *Never* force through the
  origin — the transient deposits a fixed `E0`; forcing it biases `S` high.
  **S = slope; error bar = the regression standard error.**
- **xT (upper bound)** — default `x_max`. Lower it only to exclude a known end-of-run
  artefact (periodic-image re-entry: truncate the series *before* the projectile wraps
  the cell — see §3a).
- **Optional probe only:** `detect_x0_and_stopping_power` (agent slope-plateau +
  `endpoint_status` + 40% gate) tells you *whether* the run reached steady state
  (`no_plateau` = it didn't). Use it to judge run length, **not** as the headline S.

## 3 — Method B (slab PRIMARY): deposit / thickness

For a **localised jellium slab** this is **the answer** — `dE_total / L_slab`. The
other channels (§4) are only sanity checks; do not let them displace this number.

```
S = [E_total(t_final) − E_total(t0)] / L_z
```

`L_z` = slab thickness = traversal length (e.g. 25 Bohr). `slab_stopping_power(...)`.

**Convergence gate.** `E_total(t_final)` should have converged — the deposit complete
(projectile cleared the slab *and* the electronic excitation settled). The gate
requires the energy change over the final `converge_frac` (default 15%) of the
(possibly truncated) window to be ≤ `converge_tol` (default 5%) of the total deposit.
Otherwise `status='not_converged'`: report S as a **lower bound** and say so.

**Reading the flag — caveat (hard-won, p2_classical 2026-06-25).** `not_converged`
does **not** always mean "extend the run." If a *charged* projectile's KE shows a
large **reversible** excursion across the slab (it slows climbing the mean-field
feature, speeds back up leaving it), the deposit signal keeps moving even though the
transit finished — and a longer run only brings the periodic image back sooner. When
that is the cause, the gate is flagging a *window* problem, not a *length* problem.
A useful sanity cross-check then is the **equal-potential slab-face window** (net
`ΔKE` between the two slab faces `z=±L_z/2`, where the reversible term cancels).

## 3a — Periodic-wrap truncation (before either method)

If the cell is periodic and the run is long enough for the projectile to **cross the
far boundary and re-enter** (check: is the *unwrapped* final `z` greater than the cell
length?), the wrapped image re-ploughs the medium and ruins `E_total` — full-run S can
be many× too high. **Truncate the time series before the wrap.** Quick estimate:
`t_wrap ≈ (z_exit − z0)/v0` where `z_exit` is the unwrapped position at which the
projectile clears the far boundary (and the first CAP beyond it); but read the
**actual** track for the step (the projectile decelerates, so the real wrap is later
than the `v0` estimate). Keep only steps ≤ that cut, then apply Method A or B.
(Worked example: p2_classical, cut at unwrapped `z=+45` → step 1572.)

## 4 — Sanity channels (run every time; cross-check ONLY, never the headline)

These confirm the primary number; they do **not** replace it. If one deviates from the
primary by more than ~10%, **report it to the user** for investigation — do not
silently average.

- **(a) kinetic** (`kinetic_channel`): `−dKE_ion/dx` over the same window — independent
  of the electronic `E_total` channel; their agreement *is* energy conservation.
- **(b) ∫F·v** (`force_power_channel`): cumulative `∫(−F·v)dt` with `F = m·dv/dt`.
  **CAVEAT (state it):** with `F` from the track this equals `ΔKE` analytically, so it
  is a deposition-*profile* / discretisation check, **not** a third independent
  channel (the track stores no force column).
- **(c) equal-potential slab-face** (slab runs): net `ΔKE`/`ΔE_total` between `z=±L_z/2`
  — the reversible mean-field term cancels, so it isolates the dissipative loss.

## 5 — Reporting

**The energy-deposit method is THE stopping power; the kinetic-energy channel is
ONLY a sanity check — state this every time (user, 2026-06-30).** The defined
stopping power is the energy the projectile transfers to the electrons per unit
path: Method B `ΔE_total/L_z` (slab) or Method A slope of `ΔE_total(x)` (continuous
traversal). `−dKE_ion/dx` is **never** the headline — it is the independent
conservation cross-check. Do not let the KE number stand in for the answer.

**Mandatory deliverables for every classical-projectile run:**
1. **The primary-method plot, with the result stated on it** — the `ΔE_total`
   deposit vs path (Method A) or `E_total(t)` with the deposit (Method B), the fit/
   window drawn, and the headline `S = … (±…)` annotated ON the figure.
2. **The KE sanity metric, provided alongside** — the `−dKE_ion/dx` value and its
   agreement ratio with the primary (energy conservation). Show it as a clearly
   labelled *sanity check*, not a co-headline.

The PRIMARY number is fixed by geometry (§1) — Method B for slabs, fixed-20%-time
Method A otherwise — so report *that* as the headline, not a menu. Alongside it give
its error/uncertainty, the window `[x0, xT]` (Method A) or `L_z` and the convergence
verdict (Method B), the guard results (N-conservation, ΔE_total vs ΔKE), the
sanity-channel cross-checks, and any flags (`not_converged`, `range_too_short`,
`endpoint_contaminated`, `no_plateau`). When a flag fires or a sanity channel
diverges, surface it (report **both** the primary and the KE numbers so the
discrepancy is visible); the verdict on accepting a flagged number is the **user's**.

## Validation status (worked examples)

Reproduced known cases via the deterministic builders in
`docs/validation/stopping-power-extraction/`:
- **bulk** — `transient_method_comparison.ipynb`: fixed-20% ≈ point-charge
  Lindhard on fast σ-sweep runs; agent probe returns `no_plateau` (signals not-steady).
- **slab** — `p5_classical_transient_comparison.ipynb`: Method B reproduces the
  campaign `S = ΔE_bath/25 = 0.93 eV/Bohr`; convergence gate flags marginal.
- **slab + periodic-wrap truncation** — `p2_classical_truncated_stopping.ipynb`
  (qsp_phase2): truncate at unwrapped `z=+45` (step 1572) → full-run Method B 0.104
  Ha/Bohr (wrap garbage) drops to 0.0125 Ha/Bohr `not_converged`; the reversible-well
  caveat (§3) is shown via the equal-potential slab-face cross-check (0.95–0.99 eV/Å,
  matches the independent 0.0186 reference).
- `stopping_power.py` ships a portable `_selftest()` (`python3 stopping_power.py`):
  known-slope recovery (incl. `fixed_time_fraction`), slab converged vs not-converged,
  the N-drainage guard.

## References
- Correa 2018 — `docs/sources/correa-2018-electronic-stopping-power.md` (Eq.10, Fig.8).
- Lindhard linear-response reference: `inqview.analysis.lindhard_elf`
  (`stopping_power_point`, `stopping_power_sigma`).
