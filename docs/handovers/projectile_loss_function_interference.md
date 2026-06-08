# Handover: interference-free projectile loss function

Task: obtain the loss function L(q,ω) **from a moving projectile (NOT a kick)**
that is free of interference effects, in the L=50 jellium box, without
increasing the cell size.

## User comments / position (this session, 2026-06-01 — the user's voice)

- The density visualisations show that for the ~15 eV plasmon run the loss
  function is "so populated" because of **interference effects**.
- "I want to have a run that produced the loss function without any interference."
- Skeptical that wrapping does not corrupt the loss function: "there is no reason
  why this should not affect. After all, the density enters into a steady state,
  whose frequency will be captured by the loss function." Asked: "Does the
  plasmon and the interference effects have similar period?"
- "I don't want a kicked system. Instead, I want a projectile. I want to obtain
  the loss function of a projectile."
- After inspecting the E3.4 (3.4 eV) run: **"even the 3.4 eV run has interference
  effects."** → off-resonance frequency-separation of the wrap line from the
  plasmon is NOT sufficient; the real-space wrapping itself is disqualifying.
- **The interference patterns are very evidently visible in the E3.4 run, and
  they visually look like plasmons.** This is the critical problem: the
  interference is not a subtle background — it is plainly visible AND it
  resembles/mimics genuine plasmon oscillations, so it cannot be trusted or
  disentangled by eye. **A different approach is therefore required** (the
  off-resonance projectile route is abandoned).

## Current status

Open / unresolved — **a different approach is required** (per the user). No
accepted interference-free projectile run yet. The off-resonance-projectile
approach (E3.4, E25) is **rejected by the user**: the projectile still wraps the
periodic box, the resulting interference patterns are **very evidently visible
and visually mimic plasmons**, so they cannot be trusted or separated by eye —
even though the wrap frequency is nominally separated from the plasmon. All
periodic-box projectile routes are therefore considered unsuitable; the path
forward must eliminate wrapping itself (see Exact next steps).

## Established physics (verified numbers, this session)

System: L=50, N=162, r_s=5.69. Plasmon ω_p=√(4πn)=3.47 eV, period T_p=2π/ω_p=49.2 a.u.
Wrap (single box transit) period T_wrap=L/v; wrap/kinematic frequency ω_kin(m)=m·v·q₁,
q₁=2π/L=0.1257 Bohr⁻¹.

| run | v (a.u.) | T_wrap (a.u.) | ω_kin(m=1) | |ω_kin−ω_p| | wraps in 2000 a.u. |
|---|---|---|---|---|---|
| E3.4 | 0.500 | 100.0 | 1.71 eV | 1.76 eV | ~20 |
| E15 (resonant) | 1.050 | 47.6 | 3.59 eV | **0.12 eV** | ~42 |
| E25 | 1.356 | 36.9 | 4.64 eV | 1.16 eV | ~54 |

- FFT bin at T=2000 a.u.: Δω=0.086 eV. At E15, wrap and plasmon are ~1.4 bins
  apart → **not separable** (E15 was tuned to v_res where ω_kin≈ω_p by design).
- Off-resonance (E3.4, E25): wrap line separated from the plasmon in FREQUENCY,
  but the projectile STILL wraps ~20–54× in REAL SPACE. User considers the
  real-space wrapping itself an interference effect → off-resonance insufficient.
- A periodic-box projectile ALWAYS carries a wrap line; it cannot be removed,
  only frequency-shifted. ⇒ no periodic-box projectile run is interference-free.

## What changed

- Found + fixed a separate post-processing bug (moving-WP dipole residual in the
  v2 density-wake GIFs; see presentation_evidence_RESUME handover). Not this task.
- Built diagnostic `scripts/projectile_loss_clean.py` →
  `docs/presentations/storyline/tasks/batch2_figures/projectile_loss_clean.png`
  (E3.4/E15 projectile L(q,ω) with wrap vs plasmon lines labelled).

## Files touched
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/scripts/projectile_loss_clean.py`
- `/local/data/public/skcb2/tddft/docs/presentations/storyline/tasks/batch2_figures/projectile_loss_clean.png`
- E3.4 visualisations the user inspected:
  `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_plasmon_n162_L50_E3p4_varyv/results/analysis/`
  (density/wp_xz.gif, density/wp_z_profile.gif, observables/n_q_vs_time.png,
  observables/n_q_spectrum.png).

## Tests and validation
- Period/frequency arithmetic verified numerically (table above).
- projectile_loss_clean.py known-case: E3.4 m=1 peak 3.50 eV (=ω_p), wrap 1.71 eV.

## Trusted sources used
- Dielectric stopping / loss function: Lindhard 1954; Ritchie 1959.
- RT-TDDFT spectral resolution Δω=2π/T (Yabana–Bertsch convention).

## Attribution notes
- Bohm-Gross plasmon dispersion; r_s/ω_p from standard free-electron-gas relations.

## Known issues / blockers
- **Core blocker:** in a periodic box a projectile must wrap; the wrap is
  interference the user rejects. Cannot be fixed by energy choice or dt.
- A spatial absorbing potential (CAP) CANNOT be used to absorb the WP at the
  boundary here, because the uniform jellium bath fills the boundary too — a CAP
  would drain the bath, not just the projectile.

## Assumptions still in play
- User wants a PROJECTILE-excited loss function (kick rejected).
- Box size fixed at L=50 (no enlargement).
- σ=1, launch at boundary+4σ were earlier-stated preferences (may be revisited).

## Exact next steps (options for the user to choose)
1. **Remove-the-WP-before-it-wraps (custom INQ code):** inject WP, single transit,
   remove the WP orbital at the wrap boundary, keep propagating the 162-orbital
   bath long → bath rings at ω_p, no wrap. The only projectile route with NO wrap
   line. Needs custom mid-run orbital removal in a new run.cpp (not in the current
   jellium run_template). One-time transit kinematic remains (broad, weak).
2. Reconsider enlarging the box (user has so far declined) — a single pass long
   enough to resolve ω_p needs T>49 a.u. → L≈100 (new GS).
3. Accept the kick route (user has so far declined) for a clean medium L(q,ω).
Pending the user's choice; do NOT launch a new run until chosen.
