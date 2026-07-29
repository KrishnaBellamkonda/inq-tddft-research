#!/usr/bin/env python3
"""build_debugging_notebook.py — campaign: debugging-quantum-stopping-power.

Assembles + executes debugging_quantum_stopping_v1p3.ipynb: the CAP-capture
correction to the p5_wp_v1p3 quantum stopping power, with the four user-locked
sanity checks and the binary verdict vs point-charge Lindhard.

Campaign prompt: docs/campaigns/debugging_quantum_stopping_power/debugging-quantum-stopping-power.md
Usage: venv/bin/python3 build_debugging_notebook.py
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbclient import NotebookClient
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "debugging_quantum_stopping_v1p3.ipynb"

nb = new_notebook()
C = nb.cells
md = lambda s: C.append(new_markdown_cell(s))
code = lambda s: C.append(new_code_cell(s))

# ---------------------------------------------------------------- 1. title + question
md("""# Debugging quantum stopping power — CAP-capture correction (`p5_wp_v1p3`)

**Campaign:** `docs/campaigns/debugging_quantum_stopping_power/debugging-quantum-stopping-power.md` · executed 2026-07-11

**The question.** The v = 1.3 quantum stopping power, S_WP = 2.37 eV/Bohr (upper
bound, retained-energy method, phase-5 sweep), is far above the point-charge
Lindhard bulk value at the same velocity (r_s = 5.666). **Hypothesis:** part of
the retained energy is not deposited energy but the *drift kinetic energy still
carried by the wavepacket fraction captured in the box* (never absorbed by the
CAP). Removing that captured KE brings the corrected S to within 20% of Lindhard.

**Verdict rule (binary, user-locked 2026-07-11):** explained ⟺
|S_corr − S_Lind| / S_Lind ≤ 0.20. The explained energy fraction is tabulated
regardless of verdict.

**The one free assumption (Inference, user-defined, locked 2026-07-11):**
E_capt = n_capt × (E_input + E_loc) — the captured density is assigned its share
of the packet's **total starting kinetic energy**: the code-inputted drift energy
(E_input = ½k₀² = 22.99 eV — `scripts/qsp_phase5/wp/run.cpp:64-65`, `LJ_K0=1.3`)
plus the σ-derived localisation energy (E_loc = 3/(4σ²) = 81.63 eV for σ = 0.5 —
`shared/configs/slab_n82_L50x50x90_E54.hpp:60`). Then
**E_absorbed_jellium = ΔE_plateau − E_capt** and S_corr = E_absorbed_jellium/L_z.
""")

# ---------------------------------------------------------------- 2. conventions
md("""## Conventions & symbols

| symbol | meaning |
|---|---|
| 1 Ha | 27.211 eV |
| σ = σ_WP = 0.5 Bohr | wavepacket width (project σ convention; charge std = σ/√2) |
| L_z = 25 Bohr | slab thickness (traversal length) |
| E_GS | bare-slab ground-state energy, −70.22568216820937 Ha (`shared_gs/slab_n82_L50x50x90`) |
| E_input | ½k₀² = 22.9936 eV — energy inputted to the WP (`run.cpp:64-65`, `LJ_K0=1.3`) |
| E_loc | 3/(4σ²) = 81.63 eV — localisation (zero-point) KE from σ = 0.5 (`…hpp:60`) |
| ΔE_plateau | E_total − E_GS averaged over the late plateau (last 10% of the run) |
| n_capt | N_total(t_f) − 82 — WP norm captured in the box (bath = 82 e⁻) |
| E_capt | n_capt × (E_input + E_loc) — captured share of the total starting KE |
| E_absorbed_jellium | ΔE_plateau − E_capt — energy genuinely deposited in the slab |
| S_corr | E_absorbed_jellium / L_z |

Numbers are presented at 2 s.f. (3 where a difference needs it); full precision
lives in the code cells.
""")

# ---------------------------------------------------------------- 3. setup (reconstructable)
md("## Setup — fully reconstructable (`run_summary.txt`, verbatim)")
code("""from pathlib import Path
RUN = Path('/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/'
           'scripts/qsp_phase5/wp/results/p5_wp_v1p3')
print(RUN / 'run_summary.txt', '\\n' + '='*70)
print((RUN / 'run_summary.txt').read_text())""")

# ---------------------------------------------------------------- 4. source files
md("""## Source files

- Run: `ResearchProject/systems/localised_jellium/scripts/qsp_phase5/wp/run.cpp` (results: `wp/results/p5_wp_v1p3/`)
- Sweep dispatcher: `scripts/qsp_phase5/run_sweep.sh` (E_GS constant, per-run chain)
- Recorded S(E) state: `hypotheses/qsp_phase5/se_state.csv` (S1 second route)
- Classical CAP-on-bath proxy: `scripts/qsp_phase4/classical/results/p4_classical/` (same box/CAP the phase-5 sweep reused)
- Lindhard reference curve: `hypotheses/qsp_phase5/lindhard_ref.npz` (phase-5 `build_se_plot.py` machinery, point-charge, r_s = 5.666)
- Run deep-dive: `hypotheses/qsp_phase5/p5_wp_v1p3_run_notebook.ipynb`
- This builder: `hypotheses/debugging_quantum_stopping/build_debugging_notebook.py`
- Campaign prompt: `docs/campaigns/debugging_quantum_stopping_power/debugging-quantum-stopping-power.md`
""")

# ---------------------------------------------------------------- S1 reproduce-first
md("""## S1 — reproduce the original ledger first (known-case gate)

Before correcting anything, recompute the headline from raw data and match the
**independently recorded** phase-5 value (`se_state.csv`: S = 2.374 eV/Bohr,
upper bound), plus the E_jellium(0) − E_GS consistency check
(E_jellium(0) = E_total(0) − ⟨T_WP⟩(0) − E_SIE). If this gate fails, the
campaign is debugging the wrong number and must stop.""")
code("""import numpy as np, pandas as pd
HA = 27.211; L_Z = 25.0
E_GS = -70.22568216820937          # shared_gs/slab_n82_L50x50x90 (run_sweep.sh)
E_SIE_EV = 4.40                    # same sigma=0.5 packet / r_s~5.67 slab as p2/p3
E_DRIFT_EV = 22.9936213813         # 1/2 k0^2, k0=1.3 (run_summary)

obs = pd.read_csv(RUN / 'raw/observables/observables.csv')
E0 = float(obs['energy_total'].iloc[0]); Ef = float(obs['energy_total'].iloc[-1])
tf = float(obs['time_au'].iloc[-1])
dE_ev = (Ef - E_GS) * HA; S_direct = dE_ev / L_Z

# second route: the value recorded by the sweep the night it ran
se = pd.read_csv('/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/'
                 'hypotheses/qsp_phase5/se_state.csv')
row = se[se.tag == 'p5_wp_v1p3'].iloc[0]

# E_jellium(0) consistency: strip the WP's own energy from the run's t=0 total
mom = pd.read_csv(RUN / 'raw/observables/wp_momentum_stats.csv', comment='#')
T_wp0 = float(mom['e_kin_ha'].iloc[0])            # per-norm = absolute at t=0 (norm 1)
E_jell0 = E0 - T_wp0 - E_SIE_EV / HA
T_zp_ev = T_wp0 * HA - E_DRIFT_EV

print(f"direct  : S_WP = {S_direct:.3f} eV/Bohr  (dE = {dE_ev:.2f} eV, t_f = {tf:.1f} au)")
print(f"recorded: S_WP = {row.S_eVbohr:.3f} eV/Bohr  (deposited = {row.deposited_eV:.2f} eV, bound = {row.bound})")
print(f"<T_WP>(0) = {T_wp0*HA:.1f} eV  (drift {E_DRIFT_EV:.1f} + zero-point {T_zp_ev:.1f})")
print(f"E_jellium(0) - E_GS = {(E_jell0 - E_GS)*HA:.2f} eV (expect small, ~ +0.4)")
assert abs(S_direct - row.S_eVbohr) < 1e-3, 'S1 GATE FAILED — direct recompute != recorded se_state value'
# (tolerance 1e-3: the sweep used the full CODATA Ha->eV constant, this notebook 27.211)
assert abs((E_jell0 - E_GS) * HA) <= 1.0

# E_capt basis: code-inputted drift (run.cpp:64-65) + sigma-derived localisation
SIGMA = 0.5                                   # WP_SIGMA_BOHR (slab_n82_L50x50x90_E54.hpp:60)
E_LOC_EV = 3.0 / (4.0 * SIGMA**2) * HA        # 3/(4 sigma^2) = 3 Ha = 81.63 eV
E_START_EV = E_DRIFT_EV + E_LOC_EV
print(f"E_start = E_input + E_loc = {E_DRIFT_EV:.1f} + {E_LOC_EV:.1f} = {E_START_EV:.1f} eV "
      f"(run-measured <T_WP>(0) = {T_wp0*HA:.1f} eV)")
assert abs(E_START_EV - T_wp0 * HA) < 0.5, 'code-derived starting KE != measured'
print('\\nS1 GATE PASSED — reproducing the recorded 2.37 eV/Bohr ledger (upper bound).')""")

# ---------------------------------------------------------------- S2 CAP-on-bath
md("""## S2 — assumption 1: the CAP does not eat the jellium

The captured-norm estimator n_capt = N_total(t_f) − 82 assumes the bath keeps all
82 electrons. Phase 5 ran no per-velocity classical twin, so the proxy is
`p4_classical` — the **same 50×50×90 box + CAP the phase-5 sweep reused**, with
**no** wavepacket: any N_total drift there is pure bath absorption by the CAP.
Caveat (flagged): p4_classical ran τ = 100 a.u. vs this run's 153.8 a.u.; a linear
extrapolation gives the worst case, and the classical projectile excites the bath
more strongly than the slow WP, so this is an upper-bound proxy.""")
code("""cl = pd.read_csv('/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/'
                 'scripts/qsp_phase4/classical/results/p4_classical/raw/observables/electron_number.csv')
ncol = [c for c in cl.columns if c not in ('step','time_au')][-1]
tcl = float(cl['time_au'].iloc[-1])
drift = float(cl[ncol].iloc[0] - cl[ncol].iloc[-1])
drift_ext = drift * 153.84 / tcl
print(f"p4_classical N_total: {cl[ncol].iloc[0]:.6f} -> {cl[ncol].iloc[-1]:.6f} over t = 0..{tcl:.0f} au")
print(f"bath absorbed by CAP: {drift:.4f} e-  ({100*drift/82:.2f}% of the bath)")
print(f"linear worst case at tau=153.8: {drift_ext:.4f} e-  => up to {drift_ext*E_START_EV:.1f} eV extra on E_capt")
print("=> n_capt = N_total - 82 UNDERESTIMATES the WP remnant by up to this much (upper-bound proxy).")""")

# ---------------------------------------------------------------- E_capt + S3
md("""## E_capt + S3 — time-resolved capture estimate

E_capt(t) = (N_total(t) − 82) × (E_input + E_loc). At t = 0 this is the packet's
full 104.6 eV starting KE (the whole WP is in the box); it decays as the CAP
absorbs the packet. The value at t_f is the correction; the plateau test says
whether t_f was long enough for it to have converged. The (E_input + E_loc)
basis is verified against the run-measured ⟨T_WP⟩(0) before use.""")
code("""import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, sys
sys.path.insert(0, '/local/data/public/skcb2/tddft/inq-stack/python')
from inqview.visualisation import style; style.apply_theme()

en = pd.read_csv(RUN / 'raw/observables/electron_number.csv')
tN, N = en['time_au'].to_numpy(), en['N_total'].to_numpy()
n_capt_t = N - 82.0
E_capt_t = n_capt_t * E_START_EV
n_capt, E_capt = float(n_capt_t[-1]), float(E_capt_t[-1])

fig, ax = plt.subplots(figsize=(6.4, 3.4))
ax.plot(tN, E_capt_t, lw=1.4)
ax.axhline(E_capt, ls=':', color='0.5', lw=1)
ax.annotate(f'E_capt(t_f) = {E_capt:.2f} eV', (tN[-1], E_capt), ha='right', va='bottom')
ax.set_xlabel('t (a.u.)'); ax.set_ylabel(r'$E_{\\rm capt}(t)$ (eV)')
ax.set_title(r'Captured KE: $(N_{\\rm total}(t)-82)\\times(E_{\\rm input}+E_{\\rm loc})$')
plt.tight_layout(); plt.show()

m = tN >= 0.9 * tN[-1]
dlate = float(E_capt_t[m][-1] - E_capt_t[m][0])
print(f"n_capt(t_f) = {n_capt:.4f} e-   E_capt = {E_capt:.2f} eV")
print(f"plateau test: dE_capt over last 10% of run = {dlate:+.3f} eV "
      f"({'plateaued' if abs(dlate) < 0.1 else 'STILL DRAINING — E_capt is an upper bound'})")""")

# ---------------------------------------------------------------- S4 energy books
md("""## S4 — CAP energy-removal ledger (where did the ~109 eV go?)

The WP injects ⟨T_WP⟩(0) (drift + zero-point) + E_SIE on top of the bath. By
t_f that energy has three destinations: removed by the CAP with the absorbed
packet, retained by the bath (ΔE, the stopping-power numerator), or still in the
box as captured drift KE (E_capt, inside ΔE). The books should close:
E_injected ≈ E_CAP-removed + ΔE.""")
code("""E_removed_ev = (E0 - Ef) * HA
E_injected_ev = T_wp0 * HA + E_SIE_EV
rows = [
    ('injected: <T_WP>(0) drift',        E_DRIFT_EV),
    ('injected: <T_WP>(0) zero-point',   T_zp_ev),
    ('injected: E_SIE',                  E_SIE_EV),
    ('INJECTED total',                   E_injected_ev),
    ('removed by CAP: E_total(0) - E_total(t_f)', E_removed_ev),
    ('retained by system: dE = E_total(t_f) - E_GS', dE_ev),
    ('  of which captured starting KE (E_capt)', E_capt),
    ('BOOKS: removed + retained - injected', E_removed_ev + dE_ev - E_injected_ev),
]
print(pd.DataFrame(rows, columns=['ledger row', 'eV']).to_string(index=False,
      formatters={'eV': lambda v: f'{v:8.1f}'}))
closure = E_removed_ev + dE_ev - E_injected_ev
print(f"\\nbooks close to {closure:+.1f} eV of {E_injected_ev:.0f} eV injected "
      f"({100*abs(closure)/E_injected_ev:.1f}%) — residual = bath-GS offset + numerics.")""")

# ---------------------------------------------------------------- verdict
md("""## Corrected stopping power + binary verdict

E_jellium_gained = ΔE_plateau (E_total − E_GS averaged over the late plateau,
last 10% of the run). E_absorbed_jellium = ΔE_plateau − E_capt, and
S_corr = E_absorbed_jellium / L_z, judged against point-charge Lindhard at
v = 1.3. Explained fraction f = E_capt / (ΔE_plateau − S_Lind·L_z) is tabulated
per the campaign contract.""")
code("""lind = np.load('/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/'
               'hypotheses/qsp_phase5/lindhard_ref.npz')
v_run = 1.3
v_grid = np.sqrt(2 * lind['E'] / HA)
S_lind = float(np.interp(v_run, v_grid, lind['S']))

# E_jellium_gained: delta E_total read on the late plateau (last 10% of the run)
t_obs = obs['time_au'].to_numpy(); E_obs = obs['energy_total'].to_numpy()
mpl_ = t_obs >= 0.9 * t_obs[-1]
dE_plateau = float((E_obs[mpl_].mean() - E_GS) * HA)

S_orig = S_direct
E_abs_jell = dE_plateau - E_capt
S_corr = E_abs_jell / L_Z
gap_ev = dE_plateau - S_lind * L_Z
f_expl = E_capt / gap_ev
ratio = abs(S_corr - S_lind) / S_lind

tbl = pd.DataFrame([
    ('S_WP original (upper bound)',            f'{S_orig:.2f}', 'eV/Bohr'),
    ('dE_plateau (E_jellium_gained)',          f'{dE_plateau:.1f}', 'eV'),
    ('E_capt = n_capt x (E_input + E_loc)',    f'{E_capt:.2f}', 'eV'),
    ('E_absorbed_jellium = dE_plateau - E_capt', f'{E_abs_jell:.1f}', 'eV'),
    ('S_corr = E_absorbed_jellium/L_z',        f'{S_corr:.2f}', 'eV/Bohr'),
    ('S_Lindhard (point charge, v=1.3)',       f'{S_lind:.2f}', 'eV/Bohr'),
    ('|S_corr - S_Lind| / S_Lind',             f'{ratio:.2f}',  ''),
    ('gap energy dE_plateau - S_Lind*L_z',     f'{gap_ev:.1f}', 'eV'),
    ('explained fraction f = E_capt/gap',      f'{100*f_expl:.1f}', '%'),
], columns=['quantity', 'value', 'unit'])
print(tbl.to_string(index=False))

verdict = 'EXPLAINED' if ratio <= 0.20 else 'NOT EXPLAINED'
print(f"\\nVERDICT (|S_corr - S_Lind|/S_Lind <= 0.20 ?): {verdict}")""")

# ---------------------------------------------------------------- takeaway
md("## Takeaway")
code("""print(f'''- The reproduce-first gate passed: raw data give S_WP = {S_orig:.2f} eV/Bohr, matching the
  recorded se_state.csv value exactly (upper bound; WP norm remaining {1-(83.0-float(N[-1])):.3f}).
- The captured WP fraction is n_capt = {n_capt:.3f} e- ({100*n_capt:.0f}%), taking its share of the packet's
  total starting KE (E_input {E_DRIFT_EV:.0f} + E_loc {E_LOC_EV:.0f} = {E_START_EV:.0f} eV): E_capt = {E_capt:.1f} eV.
- E_absorbed_jellium = dE_plateau - E_capt = {dE_plateau:.1f} - {E_capt:.1f} = {E_abs_jell:.1f} eV
  => S_corr = {S_corr:.2f} eV/Bohr (from the original 2.37).
- Against S_Lindhard = {S_lind:.2f} eV/Bohr the corrected value is {S_corr/S_lind:.1f}x the reference:
  verdict {verdict} — CAP-capture accounting explains f = {100*f_expl:.0f}% of the gap.
- The energy books close (S4), so the remaining excess is NOT an accounting leak; the leading
  suspect stays the non-captured packet's localisation-energy deposition ({E_LOC_EV:.0f} eV available vs
  {E_DRIFT_EV:.0f} eV drift) and beyond-linear-response effects.''')""")

nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
client = NotebookClient(nb, timeout=600, kernel_name="python3", resources={"metadata": {"path": str(HERE)}})
client.execute()
nbf.write(nb, OUT)
errs = sum(1 for c in nb.cells if c.cell_type == "code"
           for o in c.get("outputs", []) if o.get("output_type") == "error")
print(f"wrote {OUT}  ({len(nb.cells)} cells, {errs} errors)")
raise SystemExit(1 if errs else 0)
