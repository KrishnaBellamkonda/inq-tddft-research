# Notebook index — classical vs wavepacket twin analysis

All paths relative to `ResearchProject/systems/localised_jellium/hypotheses/`.
Read top-to-bottom: overview → method validation → per-pair physics.

## 0. START HERE — campaign overview

| Notebook | What it shows |
|---|---|
| **`twin_dynamics/SYNTHESIS_cross_pair.ipynb`** | The whole campaign in one place: the **universal gauge test** (no gauge in any pair → every Δ is physical), the **σ-ladder power laws** (ZPE ∝ 1/σ², self-Hartree ∝ 1/σ, one-electron SIE → 0 as σ grows), and the **phenomenon comparison** (residual/dKin/E_deposited vs regime). The master table of all pairs lives in cell 1. |

## 1. Method validation — read to TRUST the decomposition

| Notebook | What it shows |
|---|---|
| `perturbation_method/perturbation_method_study.ipynb` | Foundational energy book-keeping: the residual = WP self-Hartree, the empirical (boundary-matched) self-Hartree, the LDA one-electron **SIE = 4.34 eV**, and robustness stress tests (σ, r, Lz, grid, p2-vs-p3). Why the Gaussian *perturbation* is used for the classical projectile. |
| `twin_dynamics/pdyn_ix_pairwise_study.ipynb` | The **pairwise decomposition validated**: classical & WP each decompose EXACTLY into E_SS/E_PP/E_PS/E_SB/E_PB (closure ~1e-10), the gauge test = 0.0000 (no inter-run gauge), ΔE_PP is physical WP dispersion. This is the proof the per-pair terms below are trustworthy. |

## 2. Per-pair deep dives — the physics (each covers BOTH twins)

Each `study.ipynb` describes the classical AND wavepacket run together, with: parity/provenance, findings + closure gates, **pairwise Coulomb decomposition + gauge test**, **n(z,t) density carpets (classical, WP, Δn=WP−classical)**, **WP−classical energy-budget bar plot**, **pairwise-energy GIF** (`pairwise_evolution.gif`), conservation gate, WP dispersion σ(t).

| Notebook | Pair | Phenomenon to look for |
|---|---|---|
| `twin_dynamics/p5_null_s2_k4/study.ipynb` | σ=2, k₀=4.2 | **null control** — twins agree (Δ→0 except SIE); the falsifier. Read first to calibrate "no effect". |
| `twin_dynamics/p1_reflect_s2_k04/study.ipynb` | σ=2, k₀=0.4 | **quantum reflection** — WP retains upstream charge at the attractive surface; classical never reflects. Watch Δn upstream + E_PP splitting. |
| `twin_dynamics/p4_capture_s2_k11/study.ipynb` | σ=2, k₀=1.1 | **capture vs escape** — classical trapped + rings the slab; WP partly escapes. Watch the transmitted lobe in Δn. |
| `twin_dynamics/p2_tunnel_s2_k05/study.ipynb` | σ=2, k₀=0.5, launched INSIDE | **tunnelling** — classical bounces forever; WP leaks through the surface barrier. Largest immediate divergence (WP adapts to the slab; ΔE_PS/ΔE_PB from step 0). |
| `twin_dynamics/p6_ladder_s1_k11/study.ipynb` | σ=1, k₀=1.1 | **σ-ladder rung** — the σ=1 point of the localisation power laws (ZPE 20.4, R 9.8, SIE 1.13 eV). |

## 3. Superseded / reference

| Notebook | Note |
|---|---|
| `twin_dynamics/pdyn_k1_study.ipynb` | The first dynamic pair (σ=0.5, k₀=1, 50 steps) — the **dispersion-dominated** regime (WP spreads 10× before the slab). Kept as the σ=0.5 anchor; the 200-step version feeds the synthesis ladder. |

## Suggested analysis path
1. `SYNTHESIS_cross_pair.ipynb` — get the campaign-wide picture + the power laws.
2. `pdyn_ix_pairwise_study.ipynb` (+ `perturbation_method_study.ipynb`) — confirm the decomposition is exact and gauge-free.
3. The five per-pair `study.ipynb` — the actual quantum effects, in space (Δn carpets) and time (pairwise GIF), for each regime.

Data provenance: 5-pair overnight campaign (`orchestrate_overnight.py`, 2026-07-15) + the
σ=0.5 200-step pair. Every figure regenerates deterministically from the run data.
