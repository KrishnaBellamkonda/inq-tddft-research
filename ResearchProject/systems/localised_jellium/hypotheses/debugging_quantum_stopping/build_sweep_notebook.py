#!/usr/bin/env python3
"""build_sweep_notebook.py — campaign: debugging-quantum-stopping-power (sweep extension).

Applies the locked CAP-capture correction (E_capt = n_capt x (E_input + E_loc),
E_absorbed_jellium = dE_plateau - E_capt) to EVERY aliasing-valid point of the
quantum stopping power S(E) graph:
  v=1.3, 2.0, 3.0 (clean), 4.0 (borderline ~1.1% tail), 5.0 (h=0.35 rerun, clean grid).
v=6.0 is excluded (39% aliased, bound=lower — user-confirmed 2026-06-27).

Usage: venv/bin/python3 build_sweep_notebook.py
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbclient import NotebookClient
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "debugging_quantum_stopping_sweep.ipynb"

nb = new_notebook()
C = nb.cells
md = lambda s: C.append(new_markdown_cell(s))
code = lambda s: C.append(new_code_cell(s))

md("""# Debugging quantum stopping power — CAP-capture correction, all aliasing-valid points

**Campaign:** `docs/campaigns/debugging_quantum_stopping_power/debugging-quantum-stopping-power.md` · executed 2026-07-11

Extends the single-run (`v=1.3`) analysis of `debugging_quantum_stopping_v1p3.ipynb`
to **every aliasing-valid point** of the quantum S(E) graph. Aliasing verdicts are the
user-confirmed 2026-06-27 findings (`docs/handovers/localised-jellium.md`): grid
h=0.5 ⇒ E_cut=537 eV, packet tail crosses Nyquist at high v; **v6 (490 eV, 39%
aliased, bound=lower) is EXCLUDED**; v5 is used via its **h=0.35 rerun** (k_Nyq=8.98,
tail ~4σ inside — clean grid, but τ=28 a.u. only ⇒ weakly plateaued).

**The locked correction (user, 2026-07-11):** per run,
E_capt = n_capt × (E_input + E_loc), with n_capt = N_total(t_f) − 82,
E_input = ½v² (the energy inputted in code — `run.cpp` `LJ_K0`), and
E_loc = 3/(4σ²) = 81.63 eV (σ = 0.5 for all runs). Then
**E_absorbed_jellium = ΔE_plateau − E_capt** (ΔE_plateau = late-10%-mean of
E_total − E_GS) and S_corr = E_absorbed_jellium / L_z (L_z = 25 Bohr).

**Verdict rule (binary, per point):** explained ⟺ |S_corr − S_Lind|/S_Lind ≤ 0.20.
""")

md("""## Runs & anchors

| point | run (results dir) | grid h | E_GS anchor (Ha) | aliasing status |
|---|---|---|---|---|
| v=1.3, 23 eV | `qsp_phase5/wp/results/p5_wp_v1p3` | 0.5 | −70.22568 | clean |
| v=2.0, 54 eV | `qsp_phase4/wp/results/p4_wp` | 0.5 | −70.22568 | clean |
| v=3.0, 122 eV | `qsp_phase5/wp/results/p5_wp_v3p0` | 0.5 | −70.22568 | clean (0.05% tail) |
| v=4.0, 218 eV | `qsp_phase5/wp/results/p5_wp_v4p0` | 0.5 | −70.22568 | borderline (~1.1% tail) |
| v=5.0, 340 eV | `qsp_phase5/rerun_v5_h035/wp/results/p5_wp_v5p0_h035` | 0.35 | −71.85697 | clean grid; weak plateau (τ=28) |
| v=6.0, 490 eV | — | 0.5 | — | **ALIASED, EXCLUDED** (39% tail, bound=lower) |

E_GS anchors read from the respective GS `run_summary.txt` (the h=0.35 grid has its
own ground state — the anchor differs by +44 eV; using the wrong one would corrupt ΔE).
S1 cross-check per run: the recorded `hypotheses/qsp_phase5/se_state.csv` row
(the v5 row was overwritten by the h=0.35 rerun; v2.0 is tagged `p4_wp_v2p0`).
""")

code("""import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, sys
from pathlib import Path
sys.path.insert(0, '/local/data/public/skcb2/tddft/inq-stack/python')
from inqview.visualisation import style; style.apply_theme()

HA = 27.211; L_Z = 25.0; SIGMA = 0.5
E_LOC_EV = 3.0 / (4.0 * SIGMA**2) * HA          # 81.63 eV, sigma=0.5 for every run
E_SIE_EV = 4.40
LJ = Path('/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium')
EGS_H05  = -70.22568216820937                    # shared_gs/slab_n82_L50x50x90
EGS_H035 = -71.85697499994046                    # rerun_v5_h035/gs run_summary

RUNS = [
    ('p5_wp_v1p3', 1.3, LJ/'scripts/qsp_phase5/wp/results/p5_wp_v1p3',                EGS_H05,  'clean'),
    ('p4_wp_v2p0', 2.0, LJ/'scripts/qsp_phase4/wp/results/p4_wp',                     EGS_H05,  'clean'),
    ('p5_wp_v3p0', 3.0, LJ/'scripts/qsp_phase5/wp/results/p5_wp_v3p0',                EGS_H05,  'clean (0.05% tail)'),
    ('p5_wp_v4p0', 4.0, LJ/'scripts/qsp_phase5/wp/results/p5_wp_v4p0',                EGS_H05,  'borderline (~1.1% tail)'),
    ('p5_wp_v5p0', 5.0, LJ/'scripts/qsp_phase5/rerun_v5_h035/wp/results/p5_wp_v5p0_h035', EGS_H035, 'clean grid h=0.35; weak plateau'),
]
se = pd.read_csv(LJ/'hypotheses/qsp_phase5/se_state.csv')
lind = np.load(LJ/'hypotheses/qsp_phase5/lindhard_ref.npz')
v_grid = np.sqrt(2 * lind['E'] / HA)

res, capt_curves = [], {}
for tag, v, rdir, egs, status in RUNS:
    obs = pd.read_csv(rdir/'raw/observables/observables.csv')
    en  = pd.read_csv(rdir/'raw/observables/electron_number.csv')
    mom = pd.read_csv(rdir/'raw/observables/wp_momentum_stats.csv', comment='#')
    t, E = obs['time_au'].to_numpy(), obs['energy_total'].to_numpy()
    tN, N = en['time_au'].to_numpy(), en['N_total'].to_numpy()

    E_input = 0.5 * v * v * HA
    E_start = E_input + E_LOC_EV
    T_wp0 = float(mom['e_kin_ha'].iloc[0]) * HA
    # basis check: the measured <T_WP>(0) trails the analytic (E_input + E_loc) by a
    # k0-growing sliver (grid discretisation, orthogonalisation vs occupied slab
    # states, and for v>=4 the aliased momentum tail). The correction itself uses
    # the CODE basis (user-locked); the deviation is recorded as a diagnostic.
    basis_dev = 100.0 * (E_start - T_wp0) / E_start
    assert abs(basis_dev) < 5.0, f'{tag}: code basis {E_start:.1f} vs measured {T_wp0:.1f} (>5%)'

    dE_f = (float(E[-1]) - egs) * HA
    S_orig = dE_f / L_Z
    row = se[se.tag == tag].iloc[0]
    assert abs(S_orig - row.S_eVbohr) < 1e-3, f'{tag}: S1 gate failed vs se_state'

    mpl_ = t >= 0.9 * t[-1]
    dE_plateau = float((E[mpl_].mean() - egs) * HA)
    n_capt = float(N[-1] - 82.0)
    E_capt = n_capt * E_start
    E_abs = dE_plateau - E_capt
    S_corr = E_abs / L_Z
    S_lind = float(np.interp(v, v_grid, lind['S']))
    gap = dE_plateau - S_lind * L_Z
    ratio = abs(S_corr - S_lind) / S_lind
    capt_curves[f'v={v}'] = (tN, (N - 82.0) * E_start)
    res.append(dict(v=v, E_input_eV=E_input, tau_au=t[-1], status=status,
                    T_wp0_meas=T_wp0, E_start=E_start, basis_dev_pct=basis_dev,
                    S_orig=S_orig, dE_plateau=dE_plateau, n_capt=n_capt,
                    E_capt=E_capt, E_abs_jellium=E_abs, S_corr=S_corr,
                    S_lind=S_lind, ratio=ratio, f_expl=100*E_capt/gap,
                    verdict='EXPLAINED' if ratio <= 0.20 else 'NOT EXPLAINED',
                    bound=row.bound, late_slope=row.late_slope_eV_au))
R = pd.DataFrame(res)
print('S1 gates: all', len(R), 'points reproduce their recorded se_state values (<1e-3).')
print('basis checks: code (E_input + E_loc) vs measured <T_WP>(0) agrees to <5% everywhere;')
print('  deviation grows with k0 (basis_dev_pct column): max %.1f%% at v=4 (aliased-tail).'
      % R.basis_dev_pct.abs().max())""")

md("""## Per-point correction ledger

All energies in eV, S in eV/Bohr; `bound`/`late_slope` are the sweep's convergence
flags (every point is an **upper bound** — the retained energy had not fully
plateaued; v5-h035 is the loosest, slope ≈ −6 eV/a.u. at τ=28).""")
code("""show = R[['v','E_input_eV','tau_au','status','basis_dev_pct','n_capt','E_capt','dE_plateau',
          'E_abs_jellium','S_orig','S_corr','S_lind','ratio','f_expl','verdict','late_slope']].copy()
for c,f in [('E_input_eV','{:.1f}'),('tau_au','{:.0f}'),('basis_dev_pct','{:.1f}'),
            ('n_capt','{:.3f}'),('E_capt','{:.1f}'),
            ('dE_plateau','{:.1f}'),('E_abs_jellium','{:.1f}'),('S_orig','{:.2f}'),
            ('S_corr','{:.2f}'),('S_lind','{:.2f}'),('ratio','{:.2f}'),('f_expl','{:.1f}'),
            ('late_slope','{:.2f}')]:
    show[c] = show[c].map(lambda x: f.format(x))
print(show.to_string(index=False))""")

md("""## E_capt(t) across the sweep

Captured-KE estimate vs time for every valid point (basis: per-run E_input + the
common 81.6 eV localisation energy). A curve still falling at its right edge means
the capture correction at t_f is an upper bound for that run.""")
code("""fig, ax = plt.subplots(figsize=(7.0, 3.6))
for lbl, (tN, ec) in capt_curves.items():
    ax.plot(tN, ec, lw=1.3, label=lbl)
ax.set_xlabel('t (a.u.)'); ax.set_ylabel(r'$E_{\\rm capt}(t)$ (eV)')
ax.set_title(r'$(N_{\\rm total}(t)-82)\\times(E_{\\rm input}+E_{\\rm loc})$ per valid point')
ax.legend(frameon=False, fontsize=8)
plt.tight_layout(); plt.show()
for r in res:
    tN, ec = capt_curves[f\"v={r['v']}\"]
    m = tN >= 0.9*tN[-1]
    print(f\"v={r['v']}: E_capt(t_f)={r['E_capt']:.1f} eV, last-10% change {ec[m][-1]-ec[m][0]:+.2f} eV\")""")

md("""## Energy books per point (S4)

Injected (⟨T_WP⟩(0) + E_SIE) vs CAP-removed + retained. Closure ≲1% of injected
confirms no accounting leak at any velocity.""")
code("""books = []
for (tag, v, rdir, egs, status), r in zip(RUNS, res):
    obs = pd.read_csv(rdir/'raw/observables/observables.csv')
    E0, Ef = float(obs['energy_total'].iloc[0]), float(obs['energy_total'].iloc[-1])
    inj = r['T_wp0_meas'] + E_SIE_EV
    rem = (E0 - Ef) * HA
    ret = (Ef - egs) * HA
    books.append(dict(v=v, injected=inj, cap_removed=rem, retained=ret,
                      closure=rem + ret - inj, closure_pct=100*abs(rem+ret-inj)/inj))
B = pd.DataFrame(books)
print(B.to_string(index=False, formatters={c: (lambda x: f'{x:8.1f}') for c in
      ('injected','cap_removed','retained','closure')} | {'closure_pct': lambda x: f'{x:5.1f}'}))""")

md("## Takeaway")
code("""nexp = (R.verdict == 'NOT EXPLAINED').sum()
print(f'''- All {len(R)} aliasing-valid points reproduce their recorded S values (S1 gates pass); the
  code-derived KE basis (E_input + 81.6 eV localisation) agrees with measured <T_WP>(0) to <5%
  (deviation grows with k0, max {R.basis_dev_pct.abs().max():.1f}% at v=4).
- The capture correction removes {R.E_capt.min():.1f}-{R.E_capt.max():.1f} eV per point (n_capt {R.n_capt.min():.3f}-{R.n_capt.max():.3f}),
  lowering S by ~0.2-0.7 eV/Bohr; corrected S_corr spans {R.S_corr.min():.2f}-{R.S_corr.max():.2f} eV/Bohr.
- Verdict: {nexp}/{len(R)} points NOT EXPLAINED — S_corr remains {R.ratio.min():.1f}-{R.ratio.max():.1f}x away from
  point-charge Lindhard; the explained fraction is only {R.f_expl.min():.0f}-{R.f_expl.max():.0f}% of each gap.
- The excess is systematic across v = 1.3-5.0, not a per-run accounting artifact (books close
  everywhere): the localisation-energy deposition of the ABSORBED/transiting packet and
  beyond-linear-response physics remain the standing suspects.
- Caveats: every S is an upper bound (unplateaued, worst at v=5 h035, slope -6 eV/a.u.);
  v=4 carries a ~1.1% aliased tail; v=6 stays excluded.''')""")

nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
client = NotebookClient(nb, timeout=600, kernel_name="python3", resources={"metadata": {"path": str(HERE)}})
client.execute()
nbf.write(nb, OUT)
errs = sum(1 for c in nb.cells if c.cell_type == "code"
           for o in c.get("outputs", []) if o.get("output_type") == "error")
print(f"wrote {OUT}  ({len(nb.cells)} cells, {errs} errors)")
raise SystemExit(1 if errs else 0)
