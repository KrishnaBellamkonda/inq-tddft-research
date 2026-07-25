#!/usr/bin/env python3
"""Builder for stopping_comparison_study.ipynb — quantum (wavepacket) vs classical
electronic stopping power in the localised jellium slab (qsp_phase2, p2), with the
linear-response and bulk-jellium references.

Executes to an .ipynb (nbformat + ExecutePreprocessor). Figures -> _figs/.

Run:
  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
  /local/data/public/skcb2/tddft/venv/bin/python3 build_stopping_comparison_report.py
"""
from __future__ import annotations

import os
from pathlib import Path

import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = Path(__file__).resolve().parent
OUT = HERE / "stopping_comparison_study.ipynb"
FIGS = HERE / "stopping_comparison_study_figs"
FIGS.mkdir(exist_ok=True)

cells = []
def md(s): cells.append(new_markdown_cell(s))
def co(s): cells.append(new_code_cell(s))

# ── 1. Title + question ──────────────────────────────────────────────────────
md(r"""# Localised jellium slab — quantum vs classical electronic stopping power

**p2 · N=82 · r_s≈5.67 · σ_WP = 0.5 Bohr · E = 100 eV · τ = 40 a.u.**

**Question.** A 100 eV projectile crosses a finite jellium slab. How does the
**quantum wavepacket** stopping power compare with the **classical (Ehrenfest)
projectile** in the *same* slab, and with the **bulk-jellium** classical S(v) and
the **point-charge linear-response** reference at this velocity?

The earlier comparison wrongly read the classical slab stopping as a low number
(0.5 eV/Bohr, the *equal-potential-faces* slab-transit loss). That is a
**lower bound**: the classical ion actually loses **73 eV** of kinetic energy
(100 → 27 eV) to the electron gas — measured *consistently* with the wavepacket
(energy deposited / slab thickness) the classical stopping is **also high**
(~2.9 eV/Bohr), comparable to the wavepacket's 2.7 eV/Bohr upper bound.
""")

# ── 2. Method / conventions ──────────────────────────────────────────────────
md(r"""## 1. Conventions & symbols

Atomic units throughout; energies reported in eV (1 Ha = 27.211 eV). Stopping
power $S$ in eV/Bohr.

| symbol | meaning |
|---|---|
| $E_{\rm total}(t)$ | total Kohn–Sham energy of the run at time $t$ |
| $E_{\rm GS}$ | ground-state energy of the **bare** slab (independent run) = −45.7588 Ha |
| $L_z$ | slab thickness = 25 Bohr (half-width 12.5) |
| $\mathrm{KE}_{\rm ion}(z)$ | classical projectile kinetic energy vs position |
| $v$ | projectile velocity (a.u.); $E=100$ eV $\Rightarrow v=k_0=2.711$ |
| $\sigma_{\rm WP}$ | wavepacket width (= 0.5 Bohr here) |

The wavepacket and bath are one inseparable Kohn–Sham system (no projectile
force/trajectory), so the WP stopping is read from the **energy balance**; the
classical Ehrenfest ion *does* have a trajectory, so its stopping is read from the
**ion kinetic-energy loss**.
""")

co(r"""import sys
sys.path.insert(0, '/local/data/public/skcb2/tddft/inq-stack/python')
import numpy as np, pandas as pd
from inqview.visualisation import style
from inqview.analysis import lindhard_elf as LR
style.apply_theme()
HA = 27.211386245988
LZ = 25.0
E_GS_HA = -45.75884855005157
V_100 = float(np.sqrt(2*100.0/HA))          # 2.711 a.u.
SYS = '/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium'
WP = f'{SYS}/scripts/qsp_phase2/wp/results/p2_wp/raw/observables'
CL = f'{SYS}/scripts/qsp_phase2/classical/results/p2_classical/raw/observables'
print('v(100 eV) =', round(V_100,3), 'a.u.;  E_GS =', round(E_GS_HA*HA,1), 'eV')""")

# ── 3. Setup ─────────────────────────────────────────────────────────────────
md(r"""## 2. Simulation setup (reconstructable)

- **Cell:** 50 × 50 × 70 Bohr, spacing 0.5; slab full 50×50 face, |z| < 12.5
  (25 Bohr thick), centred z=0. Region layout: slab · free · CAP |z|∈[25,35].
- **Electronic structure:** LDA, N=82 electrons, n₀=1.312e-3 a₀⁻³ (r_s≈5.67).
- **Dynamics:** real-time TDDFT, dt=0.02 a.u., 2000 steps → τ=40 a.u.;
  two-sided sin² CAP (η=−0.7 Ha) at |z|∈[25,35].
- **Projectile:** E=100 eV → k₀=2.711; WP σ_WP=0.5 (UPF electron_gaussian_sigma0p35,
  σ_pot=σ_WP/√2). Classical twin = matched Gaussian-electron Ehrenfest ion.
""")

md(r"""## 3. Source files

| role | path |
|---|---|
| WP run | `scripts/qsp_phase2/wp/results/p2_wp/` |
| classical run | `scripts/qsp_phase2/classical/results/p2_classical/` |
| ledger notebook (formula source) | `hypotheses/qsp_phase2/quantum_stopping_ledger_26-6-26.ipynb` |
| linear response | `inqview.analysis.lindhard_elf.stopping_power_point` |
| bulk S(v) extractor | `docs/reports/26-06-2026-meeting-emilio/build/build_section1.py` |
| this builder | `hypotheses/qsp_phase2/build_stopping_comparison_report.py` |
""")

# ── 4. WP energy-balance ─────────────────────────────────────────────────────
md(r"""## 4. Quantum (wavepacket) stopping — energy balance

Once the CAP has absorbed the wavepacket remnants, the retained energy is the
slab's; the stopping power is the energy *gained by the electron gas* per unit
path:

$$S_{\rm WP} = \frac{E_{\rm total}(t_f) - E_{\rm GS}}{L_z}.$$

Reported as an **UPPER BOUND**: the σ=0.5 packet carries ~82 eV of zero-point
energy whose fate (removed with the absorbed packet) is the one free assumption
(`quantum_stopping_ledger_26-6-26.ipynb`).""")

co(r"""o = pd.read_csv(f'{WP}/observables.csv')
E0, Ef = o['energy_total'].iloc[0]*HA, o['energy_total'].iloc[-1]*HA
S_WP = (Ef - E_GS_HA*HA) / LZ
print(f'E_total(0)={E0:.1f}  E_total(tf)={Ef:.1f}  E_GS={E_GS_HA*HA:.1f} eV')
print(f'dE = E_total(tf)-E_GS = {Ef - E_GS_HA*HA:+.1f} eV')
print(f'S_WP = {S_WP:.2f} eV/Bohr  (UPPER BOUND)')""")

# ── 5. Classical direct KE ───────────────────────────────────────────────────
md(r"""## 5. Classical projectile — ion kinetic-energy loss

The Ehrenfest ion is **not** absorbed (the CAP removes wavefunctions, not the
point ion), so its stopping is read directly from the ion KE. Two measures:

$$S_{\rm faces} = \frac{\mathrm{KE}(z{=}{-}12.5) - \mathrm{KE}(z{=}{+}12.5)}{L_z}
  \quad(\text{slab-transit, equal-potential faces — lower bound})$$
$$S_{\rm dep} = \frac{\mathrm{KE}_{\rm launch} - \mathrm{KE}_{\rm final}}{L_z}
  \quad(\text{total energy deposited / }L_z\text{ — matches the WP convention}).$$""")

co(r"""t = pd.read_csv(f'{CL}/electron_track.csv').drop_duplicates('step').reset_index(drop=True)
z = t['z'].to_numpy(); ke = t['ke_ion_ha'].to_numpy()*HA
def ke_cross(zt):
    for i in range(1,len(z)):
        if z[i-1] < zt <= z[i]: return ke[i]
    return np.nan
ke_in, ke_out = ke_cross(-12.5), ke_cross(12.5)
S_faces = (ke_in - ke_out)/LZ
S_dep   = (ke[0] - ke[-1])/LZ
print(f'launch KE={ke[0]:.1f}  final KE={ke[-1]:.1f} eV  (z {z[0]:.0f}->{z[-1]:.0f})')
print(f'faces: KE {ke_in:.1f}->{ke_out:.1f}  S_faces={S_faces:.2f} eV/Bohr  (lower bound)')
print(f'total KE loss = {ke[0]-ke[-1]:.1f} eV  ->  S_dep={S_dep:.2f} eV/Bohr  (deposited)')""")

# ── 6. References: linear response + bulk ────────────────────────────────────
md(r"""## 6. References — linear response & bulk-jellium classical

Point-charge **linear response** (Lindhard, $k_F=0.337$, r_s=5.69) and the
**bulk-jellium** classical S(v) at $\sigma_{\rm WP}=0.5$ (continuous-traversal run
`run_sv_sigma0p5`, Method-A KE-loss slope).""")

co(r"""KF = 0.33729
vg = np.logspace(np.log10(0.18), np.log10(3.2), 50)
S_lr = np.array([LR.stopping_power_point(float(v), KF) for v in vg]) * HA   # eV/Bohr
S_lr_100 = float(np.interp(V_100, vg, S_lr))
print(f'linear response @ v={V_100:.2f}: {S_lr_100:.2f} eV/Bohr')
# bulk classical sigma_WP=0.5 S(v) via the section-1 extractor
sys.path.insert(0, '/local/data/public/skcb2/tddft/docs/reports/26-06-2026-meeting-emilio/build')
import build_section1 as B1
bulk_all = B1._collect_sigma_rows()
bulk = bulk_all[min(bulk_all, key=lambda s: abs(s-0.5))].sort_values('v_au')
print('bulk sigma_WP=0.5 S(v):', list(np.round(bulk['S'].to_numpy(),2)))""")

# ── 7. Results: scientific comparison plot ──────────────────────────────────
md(r"""## 7. Result — stopping power comparison

Both localised-slab projectiles deposit ~70 eV at 100 eV → **both are high**
(WP 2.7, classical 2.9 eV/Bohr), far above the bulk/linear-response (~0.3) at this
velocity. The classical *slab-transit* lower bound (0.5) is shown for reference.""")

co(r"""import matplotlib.pyplot as plt
fig, ax = style.figure_one_col()
ax.plot(vg, S_lr, '-', color='0.15', lw=1.8, zorder=2, label='linear response')
ax.errorbar(bulk['v_au'], bulk['S'], yerr=bulk['se'], fmt='o', ms=4.5,
            color='#1f77b4', mec='k', mew=0.4, capsize=2, elinewidth=0.8,
            zorder=3, label=r'classical, bulk ($\sigma_{\rm WP}$=0.5)')
ax.plot(V_100, S_dep, 's', ms=8, color='#2e8b57', mec='k', mew=0.6, zorder=6,
        label='classical, slab (deposited)')
ax.plot(V_100, S_faces, 'v', ms=7, color='#2e8b57', mec='k', mew=0.6, mfc='white',
        zorder=6, label='classical, slab (faces, lower bound)')
ax.plot(V_100, S_WP, 'D', ms=8, color='#7d3c98', mec='k', mew=0.6, zorder=6,
        label='wavepacket, slab (upper bound)')
ax.axvline(V_100, color='0.7', lw=0.8, ls=':', zorder=1)
ax.set_xscale('log'); ax.set_xlim(0.18, 3.2); ax.set_ylim(0, max(S_WP, S_dep)*1.18)
ax.set_xlabel(r'velocity  $v$  (a.u.)')
ax.set_ylabel(r'stopping power  $S$  (eV/Bohr)')
ax.set_title('Stopping power — quantum vs classical (localised slab, E=100 eV)')
ax.legend(fontsize=6.6, frameon=False, loc='upper left', ncol=1)
fig.savefig('stopping_comparison_study_figs/stopping_compare.png', dpi=600,
            bbox_inches='tight', pad_inches=0.02)
plt.show()""")

md(r"""### Classical ion KE(z) — where the energy goes""")
co(r"""fig, ax = style.figure_one_col()
ax.plot(z, ke, '-', color='#2e8b57', lw=1.3)
for zz in (-12.5, 12.5):
    ax.axvline(zz, color='0.6', ls='--', lw=0.8)
ax.set_xlabel(r'ion position  $z$  (Bohr)')
ax.set_ylabel(r'ion kinetic energy  (eV)')
ax.set_title('Classical projectile KE vs position (slab edges dashed)')
fig.savefig('stopping_comparison_study_figs/classical_ke_z.png', dpi=600,
            bbox_inches='tight', pad_inches=0.02)
plt.show()""")

# ── 8. Summary table + takeaway ──────────────────────────────────────────────
md(r"""## 8. Summary""")
co(r"""rows = [
  ('linear response (point charge)', f'{S_lr_100:.2f}', 'Lindhard @ v=2.71'),
  ('classical, bulk jellium sigma_WP=0.5', f'{float(np.interp(V_100, bulk["v_au"], bulk["S"])):.2f}', 'Method-A S(v)'),
  ('classical, slab — faces (lower bound)', f'{S_faces:.2f}', 'equal-potential faces'),
  ('classical, slab — deposited', f'{S_dep:.2f}', 'total KE loss / Lz'),
  ('wavepacket, slab (upper bound)', f'{S_WP:.2f}', '(E_total(tf)-E_GS)/Lz'),
]
print(f'{"method":42s} {"S (eV/Bohr)":>12s}  note')
for a,b,c in rows: print(f'{a:42s} {b:>12s}  {c}')""")

md(r"""## 9. Takeaway

- **At 100 eV the classical projectile is NOT weakly stopped:** it loses ~73 eV
  of kinetic energy to the gas → $S_{\rm dep}\approx2.9$ eV/Bohr, comparable to
  the wavepacket's $S_{\rm WP}\approx2.7$ (upper bound). The earlier "0.5"
  (equal-potential faces) is only the slab-transit **lower bound**.
- **Localised ≫ bulk at this velocity:** both slab projectiles (~2.7–2.9) sit far
  above the bulk-jellium classical S(v) and the point-charge linear response
  (~0.3) at v=2.71 — a finite-slab / wake effect, not a bulk stopping.
- **Open:** the WP number is a zero-point-inflated **upper bound**; a converged
  value needs the packet fully absorbed (residual WP KE ≈ ΔE at τ=40).
""")

nb = new_notebook(cells=cells)
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3',
                             'language': 'python'}
ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
ep.preprocess(nb, {'metadata': {'path': str(HERE)}})
nbf.write(nb, str(OUT))
print(f'[wrote] {OUT}')
