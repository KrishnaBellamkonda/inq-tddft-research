# Handover — jellium electron-gas analytical reference notebook

**Task.** Build a pedagogical, run-independent `.ipynb` deriving the analytical
properties of the homogeneous electron gas (HEG) realised by the jellium bath
(N=162, L=50 Bohr → r_s=5.69), driven by a single density knob, from scratch,
cross-checked against the production Python/C++ code.

**Branch.** `overnight-gaussian-classical`. **Date.** 2026-06-17.

## Status: COMPLETE

| Item | State |
|---|---|
| Grill session (design tree resolved) | done — 6 decisions, see below |
| `CONTEXT.md` glossary updated | done (new "Jellium electron-gas analytics" section) |
| Builder script | done — `build_jellium_reference_report.py` |
| Notebook built + executed (0 errors) | done — 3 embedded figures, all cross-checks green |
| Handover | this file |

## Files (absolute paths)

- Notebook: `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/hypotheses/00_jellium_reference/jellium_reference.ipynb`
- Builder:  `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/hypotheses/00_jellium_reference/build_jellium_reference_report.py`
- Figures (also embedded): `fig_box_levels.png`, `fig_loss_map.png`, `fig_screening.png` in that folder.
- Glossary: `/local/data/public/skcb2/tddft/CONTEXT.md` (§ "Jellium electron-gas analytics (2026-06-17)").

## Design decisions (from the grill session)

1. **Placement** — `ResearchProject/systems/jellium/hypotheses/00_jellium_reference/`
   (user chose hypotheses/ over docs/; folder name `00_jellium_reference` chosen
   by user, coexists with existing `00_base/`).
2. **The knob** — `RS` (Wigner-Seitz r_s); N=162 fixed (closed shell); L derived
   `L=(N/n)^(1/3)`. r_s=5.69 → L=50.0 exactly.
3. **Eigen-energies** — analytical box spectrum only (no run-data overlay), even
   though `gs_L50_cubic_N162_dx0p40/eigenvalues.csv` exists and visibly shows the
   degeneracies (a future enhancement if wanted).
4. **Loss function** — continuous RPA L(q,ω) map + Bohm-Gross dispersion as the
   backbone, with discrete box modes q_m=2πm/L overlaid + per-mode plasmon table.
5. **Extra quantities** — all four bundles: A Fermi/density, B screening,
   C energy/electron (PZ81 correlation), D stopping/dynamics (e-h continuum,
   Landau cutoff q_c, f-sum rule).
6. **Verification** — in-notebook assertion cells (no separate pytest).

## Verification (run, passed)

Rebuilt with `PYTHONPATH=.../inq-stack/python venv/bin/python3
build_jellium_reference_report.py`. Executed notebook: **0 errors, 3 PNGs**.
In-notebook asserts ALL PASS:
- [1] from-scratch χ0/ε/loss `allclose` `inqview.analysis.lindhard_elf` (rtol 1e-10).
- [2] kF, n, ω_p, k_TF match the package helpers.
- [3] from-scratch shell degeneracies + magic numbers == `shells.hpp` table
  `[2,14,38,54,66,114,162]`; Legendre exclusions (7,15,23,28) empty.
- [4] f-sum rule ∫ωL dω = (π/2)ω_p² to <1e-2 across q∈{0.2..0.8}.

Key numbers (r_s=5.69): kF=0.3373, E_F=0.0569 Ha=1.55 eV, ω_p=0.1276 Ha=3.473 eV,
k_TF=0.6553, q_c≈0.307 a₀⁻¹, e_tot=-0.0727 Ha/electron.

## How to regenerate / change density

Edit `RS` in §2 of the notebook *or* in the builder, then:
```bash
cd ResearchProject/systems/jellium/hypotheses/00_jellium_reference
PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
/local/data/public/skcb2/tddft/venv/bin/python3 build_jellium_reference_report.py
```

## 2026-06-17 update — formula-at-point-of-use restructure

User requested the notebook be reformatted so **each formula sits immediately
above the single cell that computes that quantity**, with quantities derived
**one at a time** in dependency order (k_F → v_F → E_F → λ_F → T_F → 2k_F → …),
rather than an up-front formula dump + batch calculation cells. Done:

- Builder rewritten to interleave `md(formula)` → `code(one quantity)` pairs.
  Now 41 markdown + 33 code cells; rebuilt to **0 errors, 3 figures**, all
  cross-checks still green (same physics, same values).
- Convention added to the **notebook-making skill**
  (`/home/raid/skcb2/skcb2/tddft/.claude/skills/notebook-making/SKILL.md`):
  new "Formula placement — at point of use, one quantity at a time" rule under
  the house narrative, plus an amended section 2 and a definition-of-done bullet.

## 2026-06-17 update — Rayleigh minimum-resolution-time section (§6.10)

Added a simulation-planning section computing the **minimum propagation time** to
resolve the loss-function frequencies via the Rayleigh criterion
`T_min = 1/min(Δf)`. Decision (user unsure which set → I recommended & used):

- **e-h excitation set = continuum edges** `ω±(q_m)` at box modes `q_m=2πm/L`
  (Option A), because the per-mode FFT (`plasmon_spectrum.py`) shows the plasmon
  line + the e-h band edges — those are the actual resolvable features. Options B
  (every discrete e-h microline) and C (inter-shell gaps) rejected: B over-resolves
  the continuum (meaningless huge T_min), C is not what the momentum-resolved loss
  spectrum shows.
- **plasmon set =** `ω_pl(q_m)` (Bohm-Gross).
- Combined, pooled, global min Δf (conservative). Frequencies are angular (Ha,
  ħ=1); ordinary `f=ω/2π`, so `T_min = 1/min(Δf) = 2π/min(Δω)`.

Result at r_s=5.69 (M_MODES=6): smallest gap is **ω₊(q₁)=0.0503 Ha vs
ω₋(q₂)=0.0532 Ha**, ΔE=0.00291 Ha=0.0791 eV → **T_min ≈ 2161 a.u. ≈ 52.3 fs**
(~43k steps at dt=0.05). The binding pair is two continuum onsets (ω₋ bunches near
q=2k_F), not the plasmon — the code reports the actual pair dynamically. Adds a 4th
figure `fig_frequency_comb.png`. `M_MODES` is a parameter in the cell.

Rebuilt: **0 errors, 4 figures**, all cross-checks still pass.

## Notes / not done

- No dispatcher auto-build tail (this notebook analyses no run-set, so the
  notebook-making auto-build convention does not apply).
- Real KS-eigenvalue overlay deliberately omitted (decision 3); the csv path is
  recorded above if a future session wants to add it (conditional on RS≈5.69).
- `shells.hpp` lists (|G|²=8, deg 6); from-scratch counting gives 12. Irrelevant
  here (162 closes at |G|²=6), not chased — but worth a look if shells ≥8 are ever
  used. `shells.hpp` not modified.
- Not committed (no user request to commit).
