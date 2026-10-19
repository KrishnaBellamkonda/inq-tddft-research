---
name: slab-stopping-power
description: Use when extracting electronic stopping power S from a LOCALISED-JELLIUM SLAB run (projectile enters → deposits → exits, CAP on both halves) — for the classical(direct-potential) or wavepacket half of a twin pair. Defines S = E_absorbed / L_z and, decisively, HOW E_absorbed is measured from the data: anchored to the projectile-free ground state E_GS (never E_total(t=0)), with +E_PB restored on the classical half and the norm correction on the WP half. Invoke whenever turning a completed slab run's energy ledger into S, or auditing a slab S value.
---

# Slab stopping power — S = E_absorbed / L_z

Electronic stopping power for a **localised-jellium slab**: the projectile enters
a finite slab, deposits energy, and exits (or is absorbed). The stopping power is
the energy the slab retains per unit path through it:

```
S = E_absorbed / L_z            L_z = slab thickness (Bohr), default 25
```

`E_absorbed` = the energy left in the jellium slab **once the projectile is gone**.
This skill is about measuring `E_absorbed` **correctly** from the run's energy
ledger — the whole difficulty is that INQ's `energy_total` means a *different*
thing on the classical and the wavepacket half, and neither raw endpoint is the
deposit until it is corrected.

Kernels live skill-locally in **`slab_stopping.py`** (validated below). Scope: the
localised-jellium slab twin runs (e.g. `hypotheses/sigma56_sv`,
`wp_highdensity_sv`, `classical_highdensity_sv`) with a **CAP on both halves** and
the **pairwise decomposition** (`interactions.csv`) switched on. For bulk /
continuous-traversal runs use `stopping-power-extraction` instead (slope fit).

---

## 0 — The one anchor: E_GS, never E_total(t=0)

**`E_absorbed` is measured against the projectile-free ground state `E_GS`** (the
bath-only slab GS, read from the GS run's `run_summary.txt` — never hard-coded, it
changes with the box):

```
E_absorbed = E_total_complete(t_final) − E_GS
```

**Do NOT use `E_total(t_final) − E_total(0)`.** `E_total(0)` carries the
projectile's launch energy, which `E_GS` does not:

| | `E_total(0) − E_GS` at launch | why |
|---|---|---|
| classical | **+103 eV** (σ=3, v=3) | the static projectile↔bath term `E_PS(0)` |
| wavepacket | **+124 eV** | the wavepacket's **formation energy** (its kinetic + interaction energy) |

At `t_final` that launch energy is gone — the CAP has carried the projectile (and
its energy) out of the box — so `E_total(t_f) − E_GS` isolates the slab's retained
excitation. Referenced to `E_total(0)` instead you would subtract off the
projectile's launch energy and get a large **negative** number (WP σ=3 v=3:
`E_total(t_f) − E_total(0) = −116 eV`, versus the true deposit `+7.6 eV`). The two
references agree *only* for the classical half after +E_PB (see §1), and never for
the WP.

---

## 1 — Classical half: restore E_PB (add it)

The classical projectile is a **perturbation hook** — a moving direct erf/r
external potential (`inqkit::dynamics::moving_gaussian_projectile_potential`),
summed with the background well + CAP and handed to `real_time::propagate`. It is
**not** an electron density, so INQ integrates it against the bath electrons
(→ `E_PS`, inside `energy_external`) but **never against the +background**
(→ `E_PB`). `E_PB` is therefore **absent** from `energy_total`:

```
classical:   energy_external == E_SB + E_PS      (E_PB ABSENT — verified <5e-10 Ha)
```

Restore it — **add** the ledger term (`interactions.csv` column `e_pb`, Hartree):

```
E_absorbed(classical) = [ E_total_INQ(t_f) + E_PB(t_f) ] − E_GS
```

**Why `+E_PB` (not `−`).** `E_PB` is a *missing* energy term, not a contaminant to
remove: the complete total is `E_total_INQ + E_PB`, so you add it back with its own
stored sign. `E_PB` is itself **negative** (electron-like projectile attracted to
the + background), so adding it *reduces* the inflated raw deposit. The decisive
test: at launch nothing is deposited, so the deposit must be ≈ 0 —

```
  E_total_INQ(0) − E_GS            = +103.5 eV   (spurious: this is just E_PS(0))
  E_total_INQ(0) + E_PB(0) − E_GS  =   +0.1 eV   ← PLUS gives ~0   (correct)
  E_total_INQ(0) − E_PB(0) − E_GS  = +206.8 eV   ← MINUS is absurd
```

Physically: the slab is **neutral**, so the net projectile↔slab interaction
`E_PS + E_PB ≈ 0` at all times. INQ kept the `+E_PS` half; adding the `+E_PB` half
restores the cancellation, leaving only the genuine deposit. (Numerically `+E_PB`
equals "subtract the `E_PS` monopole tail" because `E_PB(t_f) ≈ −E_PS(t_f)`, but
the principled statement is the additive one.)

---

## 2 — Wavepacket half: norm-correct; E_PB is already inside

The WP projectile **is** a real Kohn–Sham orbital, so INQ integrates it against the
background automatically:

```
wavepacket:  energy_external == E_SB + E_PB      (E_PB PRESENT — verified <5e-10 Ha)
```

Nothing to add for `E_PB`. Two things the WP half **does** need:

1. **Norm correction.** INQ reports the orbital energy as the per-particle mean
   `occ·⟨ψ|H|ψ⟩/⟨ψ|ψ⟩` (`inq/src/hamiltonian/energy.hpp`), so as the CAP drains the
   packet its kinetic mean keeps inflating `energy_total`. Subtract it:
   ```
   E_total_corrected(t) = E_total_INQ(t) − T1·(1 − norm)      T1 = e_kin_ha
   ```
2. **Full absorption is a precondition.** `E_absorbed = E_total_corrected(t_f) − E_GS`
   is the deposit **only if the WP is fully gone** (`norm_final → 0`, e.g. ~1e-8).
   Otherwise `E_total(t_f)` still holds the WP formation energy and you are reading a
   mid-transit excitation, not the deposit. **Gate on `norm_final`; if it is not
   ~0, extend the run** (`LJ_RESUME=1` with a larger `LJ_N_STEPS`), do not report S.

```
E_absorbed(wavepacket) = [ E_total_INQ(t_f) − T1·(1−norm) ]_{t_f} − E_GS
```

The WP deposit is the **total electronic energy deposition** — the right quantum
measure. Do *not* try to read it as `−d(E_projectile_KE)/ds`: the WP orbital is not
identifiable as "the projectile" (see `stopping-power-extraction` for the classical
KE channel, which does not transfer to the quantum case).

---

## 3 — Gates (run them; do not assume)

| gate | check | if it fails |
|---|---|---|
| classical closure | `max|E_ext − (E_SB+E_PS)| < 1e-6 Ha` | E_PB may already be in `energy_total` (WP run?) or the ledger convention differs — `+E_PB` invalid |
| WP absorption | `norm_final ≲ 1e-3` (aim ~1e-8) | WP not gone; `E_absorbed` holds formation energy — extend the run |
| ledger closure (optional, strong) | classical `E_SS==E_hartree`; WP `E_SS+E_PS+E_PP==E_hartree` | ledger terms untrustworthy — see `decomposed-interaction-energies.md` |
| segments concatenated | `observables*.csv` + `interactions*.csv` merged in step order | truncating a resumed run drops the tail where the deposit lives |

`slab_stopping.py` enforces the first two automatically.

---

## 4 — Quickstart

```python
import sys; sys.path.insert(0, ".claude/skills/slab-stopping-power")
import slab_stopping as SS

E_GS_HA = 207.18323030158          # from the GS run's run_summary.txt (per box)
L_Z     = 25.0                     # slab thickness (Bohr)

# classical half — needs observables.csv + interactions.csv
cl = SS.measure(".../classical/results/cl_s5p0_v3p0/raw/observables", "classical", E_GS_HA, L_Z)
# wavepacket half — needs observables.csv + wp_momentum_stats.csv + wp_real_space_stats.csv
wp = SS.measure(".../wp/results/s5p0_v3p0/raw/observables",           "wp",        E_GS_HA, L_Z)

print(cl["S_eV_per_Bohr"], wp["S_eV_per_Bohr"], "norm_final:", wp["norm_final"])
```

`obs_dir` is whatever directory holds the CSVs (legacy sweeps keep them flat in the
run dir; campaign runs put them under `raw/observables/`).

**Validated (2026-08-16):** classical σ=3,v=3 → S=0.191; WP σ=3,v=3 → S=0.305
(reproduces the published `S_deposit_corrected` to <1%); WP σ=5,v=3 → S=0.396.

---

## 5 — Common pitfalls

| symptom | cause | fix |
|---|---|---|
| classical S inflated ~2–5× | E_PB omitted (raw `E_total`) | add `E_PB(t_f)` (§1) |
| classical S grows wrong at high v | uncancelled `E_PS` tail; runs end at fixed downstream z | +E_PB removes it (it scales with the tail) |
| WP S huge / negative | referenced to `E_total(0)` not `E_GS`; or WP not absorbed | anchor to `E_GS`; check `norm_final` |
| WP S drifts as run lengthens | reading a mid-transit excitation | extend until `norm_final → 0` |
| classical & WP not comparable | comparing net vs gross ledgers | both anchored to `E_GS`, both corrected (§1, §2) |

---

## Related

- `decomposed-interaction-energies.md` (rule) — the P/S/B pairwise ledger and its
  closure gates that make `E_PB`/`E_PS` trustworthy; every projectile run emits it.
- `twin-run-analysis` (skill) — the full per-timestep energy decomposition and the
  physical narrative of the classical-vs-WP difference.
- `stopping-power-extraction` (skill) — the classical **bulk** case (slope fit) and
  the projectile-KE / ∫F·v sanity channels; explicitly NOT for the quantum WP.
- `.claude/rules/sigma-wp-convention.md` — σ labels are σ_WP; classical σ_pot = σ_WP/√2.
- `.claude/rules/final-timestep-checkpoint.md` — resume/extend a too-short run
  (needed when the WP-absorption gate fails).
```

Consumers: `ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv/`
(`build_sv_effective_width_s6.py::classical_S_epb`, `make_etotal_traces.py`,
`s56_stopping.py`).
