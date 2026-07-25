# Handover — Localised-jellium energetics: extend E_total(0)−E_GS vs r

Task: add data points to the EXISTING H0-style energetics plot (for the 2026-07-03
Emilio meeting, `docs/reports/03-07-2026-meting-emilio/presentation_plan.md` slide 2)
at larger projectile–slab distance r. Follow-on: the U_ext↔U_H comparison formula
(TODO 2), queued after these runs land.

## 2026-07-03 — designed + LAUNCHED

**Design (user-approved via /brainstorming Q&A):**
- Same system/density as the existing runs: a=12.5, N=82, n0=1.31e-3, σ_WP=0.5,
  spacing 0.5, w=0 — new points sit on the SAME curve.
- **Box L_z=160, slab CENTERED** (user chose centered+enlarge over off-center).
  Projectile at z=−(12.5+r); r=60 → z=−72.5, margin 7.5 Bohr (=15σ) to the ±80
  boundary (same margin the existing r=40 point had). Zero-leak OK.
- **New r set (no dup with existing {4,12,20,28,36,40}):** lower {8,16,24,32},
  higher {44,48,52,56,60}; **+ r=40 overlap** (continuity check L_z=160 vs the
  existing L_z=120 r=40 point). periodicity **2 AND 3**.
- Per r: **classical ghost** (E_GS + E_total(0)) and **WP** (N_STEPS=2 → valid
  E_total). GS built ONCE per periodicity at L_z=160 and reused (amortises the
  big-box GS cost — answers the GS-convergence-time worry).

**Machinery:** `scripts/campaign_autorun/extend_r_lz160.py` (NEW; reuses the
env-driven `gs/wp/classical` binaries already built). Detached (`setsid`,
start_new_session — SIGHUP-immune per the P0b lesson). GS build has a 40-min
liveness guard + gate; r-sweep runs 2-wide (GPU0/GPU1) with one-shot retry; then
builds `runs/extend_r160/extend_r160_excess_vs_r.png` and emails (4-part).

**LAUNCHED 2026-07-03 08:33** (driver pid 2742755→2742756). Both GS builds running
(p3 GPU0 pid 2742839, p2 GPU1 pid 2742840), entering SCF. GPUs verified free
(ml2218's 48 `lmp` are CPU-only). Output: `runs/extend_r160/`.

**Flagged for the analysis step (NOT a run blocker):** periodicity-2 WP energetics
carry the open-z net-charge G=0 bias (2D Coulomb kernel term absent from the neutral
GS) — the plot marks p2-WP as RAW; trust periodicity-3 absolute values. Correction to
be applied in the TODO-2 notebook.

## 2026-07-03 — COMPLETE (all 40 runs) + continuity PASSED

Both L_z=160 GS converged (p3 E_GS=−160.93 Ha, on the L_z trend; p2=+60.32 Ha, open-z
offset). p3 had a transient charge-sloshing episode (dn 0.02→0.4 around iter 82) that
Broyden mixing self-damped back to 2e-3 by iter 125 — NO intervention needed. All 40
projectile runs completed (~1 h wall under load-49 CPU contention).

**Plot bug found + fixed:** the wp/classical binary prepends `results/` to `LJ_OUT`, so
output nests at `<rd>/results/results/raw/observables/observables.csv` (double results);
`etotal0()` read the wrong path → the auto-emailed plot was EMPTY. Fixed `etotal0()` to
glob `**/observables.csv`; regenerated `runs/extend_r160/extend_r160_excess_vs_r.png`.
(The driver's completion email carries the OLD empty plot — resend pending user ask.)

**Continuity check PASSED:** L_z=160 r=40 WP excess 86.5 eV vs existing L_z=120 85.9 eV
(0.7%). New points are on the same curve.

**Data (E_total(0)−E_GS, eV, periodicity 3):** WP ≈ flat 86.5 (ZP ~81.6 + SIE ~5,
r-independent); classical decays 164 (r=8) → 12 (r=40) → ~1 (r=60) (unscreened ghost–slab
Coulomb). WP−classical gap → ~85 eV at large r (persistent quantum offset). Periodicity-2
WP ~78 eV (8–9 eV below p3 = the raw open-z net-charge bias, correction pending).

## 2026-07-03 — study notebook built + corrected email resent

**Notebook:** `hypotheses/extend_r160/extend_r160_study.ipynb` (builder
`build_extend_r160_report.py`, uses `_nbreport.py` house helpers; figs in
`extend_r160_figs/`). Executes 0 errors. House narrative: question → conventions/
formulas → setup → sources → §1 excess vs r (+ continuity) → §2 WP−cl gap →
**§3 the right formula (U_ext↔U_H)** → takeaway. Driver now ends with the auto-build
tail call (notebook-making convention).

**KEY PHYSICS FINDING (TODO-2, from the component decomposition U=E−T−U_H−E_xc):**
- **Classical energy is ENTIRELY in U_ext** — U_H is exactly constant across r (bath
  frozen; ghost is pure external potential). ΔU_ext,cl 163→0 eV.
- **WP splits into U_ext (−120→0, WP–background attraction) and U_H (+121→0, WP–bath
  repulsion) that CANCEL to ~1%** (neutral slab) → that is WHY the WP total is flat.
- **The right comparison = classical ΔU_ext vs WP ΔU_H** (both the projectile–bath
  electrostatic term). Same order, both decay, but **CROSS near r≈22** — a real
  quantum-vs-classical difference (polarizable WP vs rigid ghost).
- **Cleaner formula (next):** isolate the explicit Hartree cross-term
  ∫∫ n_bath·n_wp/|r−r'| from the DENSITY FIELDS (VTIs), removing the self-Hartree
  cancellation assumption. Needs density outputs, not just scalar energies.

**Corrected email RESENT** (replaces the blank-plot one) with excess_vs_r.png +
right_formula_compare.png, 4-part structure. Sent OK.

### Next
1. DONE: r=40 L_z=160 vs L_z=120 overlap verified (86.5 vs 85.9 eV, 0.7%).
2. Cross-term from density fields (needs the projectile runs to write density VTIs —
   the current extend-r runs wrote scalars only; a re-run or a dedicated density run
   is needed for the clean cross-term).
2. Merge the new points into the presentation plot (same theme).
3. TODO 2: the U_ext↔U_H comparison formula. Physics (already scoped): classical
   projectile couples via **U_ext** (∫ n_bath V_ghost); WP couples via the **Hartree
   cross-term** (∫∫ n_bath n_wp/|r−r'|, buried in U_H with the spurious WP
   self-Hartree/SIE + zero-point). Matched ghost (σ_pot=σ_WP/√2) ⇒ the two SHOULD
   coincide; the "better formula" is isolating the cross-term. Data source (user):
   static-r energetics — WP ≥1 timestep for E_total, classical E_GS + E_total(0).
   Deliver as ipynb notebooks.

### Verified vs unverified
- Driver syntax (`py_compile` OK), binaries present, GS builds started + stepping:
  **verified**. GS convergence at L_z=160, r-sweep completion, the r=40 overlap,
  and the p2 net-charge correction: **unverified** (in flight / pending).
