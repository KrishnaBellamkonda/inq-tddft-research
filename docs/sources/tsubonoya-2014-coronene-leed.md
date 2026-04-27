# Source: Tsubonoya, Hu & Watanabe (2014) — Coronene LEED via TDDFT

## Full citation

Tsubonoya, K.; Hu, C.; Watanabe, K. *Time-dependent density-functional theory
simulation of electron wave-packet scattering with nanoflakes.* Phys. Rev. B
**90**, 035416 (2014). DOI: 10.1103/PhysRevB.90.035416.

## Relevance to this project

Reference paper for every coronene wave-packet scattering simulation in
`ResearchProject/systems/coronene/`. Defines the model system (coronene C24H12
in a 18.4 × 18.4 × 31.7 Å³ box = 35 × 35 × 60 Bohr³), the wave-packet form
(Gaussian, σ = 0.53 Å = 1.0 Bohr, b = 6.35 Å = 12 Bohr above the molecule plane,
E = 200 eV directed along −z), the time-stepping (Δt = 0.020 a.u., LEED window
0.077–0.25 fs), and the LEED-screen detector convention.

## Key claims used

- §II — Gaussian wave-packet form (Eq. 1): real-space envelope and complex
  phase. Adopted in `inqkit::WavePacket` (`inq-stack/include/inqkit/wavepacket/wavepacket.hpp`).
- §III — Cell, time step, and LEED accumulator window (Eq. 5). Constants
  encoded verbatim in `inq-stack/include/inqkit/config/tsubonoya_2014_coronene.hpp`
  and mirrored to `ResearchProject/systems/coronene/shared/configs/tsubonoya_2014_base.hpp`.
- Fig. 2 — Reference LEED pattern at b = 6.35 Å, σ = 0.53 Å, E = 200 eV. The
  qualitative target for `run_base/` in the replication framework.
- ALDA functional (paper §III) — paper uses ALDA; INQ uses LDA + adiabatic
  TDDFT, which is the standard ALDA implementation (the "A" in ALDA is the
  adiabatic approximation made automatically by real-time LDA-TDDFT in INQ).

## Limitations / uncertainties

- Paper uses Troullier–Martins norm-conserving pseudopotentials. INQ's
  default norm-conserving PSPs are from the pseudo-dojo library — not
  identical, so absolute energies are not expected to match. The convergence
  sweep in `ResearchProject/systems/coronene/03_ecut_convergence/` settled on
  40 Ha as the energy-minimising cutoff for INQ + dojo-PBE; the paper used
  ≈ 54 Ha (108 Ry) which on dojo PSPs introduced SCF instability (see
  `04_leed_simulation/run_002` legacy notes). The replication framework
  therefore runs at 40 Ha, not 54 Ha.
- Paper uses fixed ions; we adopt the same — no Ehrenfest dynamics in any
  replicated run.
- The paper's electron-electron interaction beyond ALDA (e.g. memory effects,
  many-body correlation outside the adiabatic approximation) is not captured;
  this is intrinsic to ALDA-TDDFT and accepted.

## Cross-references

- Plan: `/local/data/public/skcb2/tddft/docs/plans/coronene-replication.md`
- Plan source-of-truth: `/home/raid/skcb2/.claude/plans/the-main-aim-of-mutable-taco.md`
- Code:
  - `inq-stack/include/inqkit/config/tsubonoya_2014_coronene.hpp`
  - `ResearchProject/systems/coronene/shared/configs/tsubonoya_2014_base.hpp`
  - `inq-stack/include/inqkit/wavepacket/wavepacket.hpp`
- Handover: `docs/handovers/coronene-replication.md`
