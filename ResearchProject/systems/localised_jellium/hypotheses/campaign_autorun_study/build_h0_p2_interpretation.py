#!/usr/bin/env python3
"""build_h0_p2_interpretation.py — interpretation-aid notebook for the H0
periodicity-2 full-decomposition re-run (2026-07-07).

Two aids the user asked for while interpreting the results:
  1. Individual-run energy decomposition as a WATERFALL bar chart — every logged
     component in a logical order, the running sum landing exactly on E_total.
  2. Ground-state charge distribution: the electron density n_-, the positive
     jellium background n_+, and their difference n_- - n_+ (the surface dipole),
     in xz / yz / xy slices — to see whether the two distributions coincide.

Neutral: shows the data, no verdict. Run (venv + stack on path):
    PYTHONPATH=.../inq-stack/python .../venv/bin/python3 build_h0_p2_interpretation.py
    python3 -m nbconvert --to notebook --execute --inplace H0_p2_interpretation.ipynb
"""
from __future__ import annotations
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

REPO = "/local/data/public/skcb2/tddft"
LJ = f"{REPO}/ResearchProject/systems/localised_jellium"
CA = f"{LJ}/scripts/campaign_autorun"
OUT = Path(f"{LJ}/hypotheses/campaign_autorun_study")

PRE = f"""import sys, glob, csv
import numpy as np
import matplotlib.pyplot as plt  # kernel inline backend captures plt.show()
sys.path.insert(0, {REPO+'/inq-stack/python'!r})
try:
    from inqview.visualisation import style as _st; _st.apply()
except Exception: pass
from inqview import load_vti
from pathlib import Path
RUNS = Path({CA!r}) / 'runs'
HA_EV = 27.211386
# logical component order; the eight that sum to total (see energy.hpp::total())
ORDER = [('kinetic','T\\nkinetic'), ('external','E_ext\\nexternal/local'),
         ('nonlocal','E_nl\\nnon-local PP'), ('hartree','U_H\\nHartree'),
         ('xc','E_xc\\nexch-corr'), ('ion','E_ion\\nion-ion'),
         ('ion_kinetic','E_ion,kin'), ('exact_exchange','E_xx\\nexact exch')]
KEYS = [k for k,_ in ORDER] + ['total']
def comp(run, root='h0_p2'):
    f = next(iter(glob.glob(str(RUNS/root/run/'**/observables.csv'), recursive=True)))
    rr = list(csv.reader(open(f))); h,d = rr[0], rr[1]; g = lambda c: float(d[h.index(c)])
    return {{k: g('energy_'+k) for k in KEYS}}
print('helpers ready — RUNS =', RUNS)"""

WATERFALL = """# Individual-run energy decomposition — waterfall: each component adds onto the
# running sum, the last bar lands on E_total (segments add to the total, in Ha).
def waterfall(ax, run, title):
    c = comp(run); running = 0.0; xs = []
    for i,(key,lab) in enumerate(ORDER):
        v = c[key]
        ax.bar(i, v, bottom=running, width=0.62,
               color=('#c0392b' if v < 0 else '#1b6ca8'), edgecolor='k', linewidth=.5)
        if abs(v) > 1e-9:
            ax.text(i, running + v + (2 if v>=0 else -2), f'{v:+.1f}',
                    ha='center', va=('bottom' if v>=0 else 'top'), fontsize=7)
        running += v; xs.append(lab)
    # closing 'total' bar from 0
    ax.bar(len(ORDER), c['total'], width=0.62, color='#555', edgecolor='k', linewidth=.5)
    ax.text(len(ORDER), c['total'], f"{c['total']:+.1f}", ha='center',
            va=('bottom' if c['total']>=0 else 'top'), fontsize=7, fontweight='bold')
    xs.append('E_total')
    ax.axhline(0, color='.5', lw=.8)
    ax.set_xticks(range(len(xs))); ax.set_xticklabels(xs, fontsize=7)
    ax.set_ylabel('energy (Ha)'); ax.set_title(title)
    print(f'{run}: running sum = {running:.6f} Ha,  E_total = {c["total"]:.6f} Ha,'
          f'  diff = {abs(running-c["total"]):.1e}')

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
waterfall(axes[0], 'wp_r4_p2', 'wavepacket  wp_r4_p2  (quantum electron)')
waterfall(axes[1], 'cl_r4_p2', 'classical  cl_r4_p2  (ghost, z_valence 0)')
fig.suptitle('H0 periodicity-2: single-run energy decomposition (components sum to E_total)')
fig.tight_layout(); plt.show()
print('NOTE: E_nl, E_ion, E_ion,kin, E_xx are identically 0 for these LDA jellium runs'
      ' (no non-local PP, background is not an Ewald ion, pure LDA).')"""

CHARGE = """# Ground-state charge distribution: electron density n_- vs positive background n_+.
# n_- from the periodicity-2 GS VTI (physical order via load_vti; NEVER fftshift).
# n_+ constructed analytically: slab n0 for |z|<half, uniform in x,y (edge_width 0).
GS_VTI = RUNS/'h2/gs_p2_lz120/results/density_gs_system/density_gs_system.vti'
N0, HALF = 1.312e-3, 12.5
v = load_vti(str(GS_VTI), expect_centered_axis='z')
ne = v.data                          # (nx,ny,nz) electron density, e/Bohr^3
X, Y, Z = v.x, v.y, v.z
nz_prof = np.where(np.abs(Z) < HALF, N0, 0.0)
nplus = np.broadcast_to(nz_prof, ne.shape)   # uniform in x,y
diff = ne - nplus
print('int n_- =', f'{ne.sum()*np.prod(v.spacing):.2f}', ' int n_+ =',
      f'{nplus.sum()*np.prod(v.spacing):.2f}', '(both ~82 e; neutral)')

ix0, iy0, iz0 = len(X)//2, len(Y)//2, len(Z)//2
ZLIM, XYLIM = 30, 25
def crop_z(a, ax):  # a is (n_row, n_col); rows=z for xz/yz
    return a
slices = [
    ('xz  (y=0)', ne[:,iy0,:].T, nplus[:,iy0,:].T, diff[:,iy0,:].T,
     [X[0],X[-1],Z[0],Z[-1]], 'x (Bohr)', 'z (Bohr)', (-XYLIM,XYLIM,-ZLIM,ZLIM)),
    ('yz  (x=0)', ne[ix0,:,:].T, nplus[ix0,:,:].T, diff[ix0,:,:].T,
     [Y[0],Y[-1],Z[0],Z[-1]], 'y (Bohr)', 'z (Bohr)', (-XYLIM,XYLIM,-ZLIM,ZLIM)),
    ('xy  (z=0, slab interior)', ne[:,:,iz0].T, nplus[:,:,iz0].T, diff[:,:,iz0].T,
     [X[0],X[-1],Y[0],Y[-1]], 'x (Bohr)', 'y (Bohr)', (-XYLIM,XYLIM,-XYLIM,XYLIM)),
]
vmax = float(max(ne.max(), N0)); dmax = float(np.abs(diff).max())
fig, axes = plt.subplots(3, 3, figsize=(12.5, 11))
for row,(name, a_ne, a_np, a_df, ext, xl, yl, lim) in enumerate(slices):
    for col,(dat, ttl, cmap, vlim) in enumerate([
        (a_ne, f'n_-  ({name})', 'viridis', (0, vmax)),
        (a_np, f'n_+  ({name})', 'viridis', (0, vmax)),
        (a_df, f'n_- - n_+  ({name})', 'RdBu_r', (-dmax, dmax))]):
        ax = axes[row, col]
        im = ax.imshow(dat, origin='lower', extent=ext, aspect='auto',
                       cmap=cmap, vmin=vlim[0], vmax=vlim[1])
        ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3])
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(ttl, fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
fig.suptitle('Ground state (periodicity 2): electron density n_-, positive background '
             'n_+, and their difference (surface dipole)', y=1.005)
fig.tight_layout(); plt.show()
# 1-D planar profile along z: the cleanest view of where they differ
fig, ax = plt.subplots(figsize=(7,4))
ax.plot(Z, ne.mean(axis=(0,1)), color='#1b6ca8', label='n_- (electrons, planar mean)')
ax.plot(Z, nz_prof, color='#c0392b', ls='--', label='n_+ (background)')
ax.axvspan(-HALF, HALF, color='.92', zorder=0)
ax.set_xlim(-30, 30); ax.set_xlabel('z (Bohr)'); ax.set_ylabel('planar-mean density (e/Bohr^3)')
ax.set_title('Planar-mean n_-(z) vs n_+(z): spill-out + Friedel vs sharp slab')
ax.legend(frameon=False); fig.tight_layout(); plt.show()"""

EXTENDED = """# Extended-r sweep in a bigger box (Lz=200, its own periodicity-2 GS): does the
# classical excess actually reach zero? r pushed from 40 out to 76 Bohr.
import sys; sys.path.insert(0, str(RUNS.parent))   # scripts/campaign_autorun -> analyse_phase
from analyse_phase import gs_energy
EGS = gs_energy(RUNS/'h0_p2_far/gs_p2_lz200/results'); ZP = 3/(4*0.5**2)*HA_EV
def comp_far(run):
    f = next(iter(glob.glob(str(RUNS/'h0_p2_far'/run/'**/observables.csv'), recursive=True)))
    rr = list(csv.reader(open(f))); h,d = rr[0], rr[1]; g = lambda c: float(d[h.index(c)])
    return {k: g('energy_'+k) for k in KEYS}
rs = [4,12,20,28,36,44,52,60,68,76]
wp = [(comp_far(f'wp_r{r}_p2')['total']-EGS)*HA_EV for r in rs]
cl = [(comp_far(f'cl_r{r}_p2')['total']-EGS)*HA_EV for r in rs]
print('E_GS(p2, Lz=200) =', f'{EGS:.4f} Ha')
print('WP  excess (eV):', [f'{x:.1f}' for x in wp])
print('cl  excess (eV):', [f'{x:.1f}' for x in cl])
print('classical minimum:', f'{min(cl):.1f} eV at r={rs[cl.index(min(cl))]}')
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.plot(rs, wp, 'o-', color='#1b6ca8', label='wavepacket (quantum)')
ax.plot(rs, cl, 's-', color='#c0392b', label='classical ghost')
ax.axhline(0, color='.5', lw=.8)
ax.axhline(ZP, ls=':', color='#1b6ca8', lw=1, label=f'WP zero-point {ZP:.0f} eV')
ax.set_xlabel('r (Bohr from near slab face)'); ax.set_ylabel('E_tot(0) - E_GS  (eV)')
ax.set_title('H0 periodicity-2, extended r (Lz=200): excess vs distance')
ax.legend(frameon=False); fig.tight_layout(); plt.show()"""


def nb(cells):
    n = new_notebook(); n.cells = cells
    n.metadata.kernelspec = {"name": "python3", "display_name": "Python 3"}
    return n

def build():
    cells = [
        new_markdown_cell(
            "# H0 periodicity-2 — interpretation aids\n\n"
            "*Companion to `H0_base_difference.ipynb`. Two views requested while "
            "interpreting the full-decomposition re-run (2026-07-07): (1) a single-run "
            "energy decomposition whose segments sum to `E_total`, and (2) the ground-"
            "state electron density vs the positive background and their difference. "
            "Numbers/fields read from the run files; no interpretation added.*"),
        new_code_cell(PRE),
        new_markdown_cell(
            "## 1. Individual-run energy decomposition (segments add to E_total)\n"
            "Waterfall bars in a logical order — kinetic, external (local pseudopotential "
            "/ background), non-local PP, Hartree, xc, ion-ion, ionic-kinetic, exact "
            "exchange — each stacked onto the running sum so the final position lands on "
            "`E_total`. Left = wavepacket, right = classical ghost, both at r=4 (closest, "
            "largest interaction). The printed `diff` confirms Σcomponents == E_total."),
        new_code_cell(WATERFALL),
        new_markdown_cell(
            "## 2. Ground-state charge distribution — electrons vs positive background\n"
            "Do the negative (electron) and positive (jellium background) distributions "
            "coincide? Rows: xz (y=0), yz (x=0), xy (z=0, slab interior). Columns: n_- "
            "(electrons), n_+ (background), and n_- - n_+ (the difference / surface "
            "dipole). The background is a sharp slab (uniform in x,y, |z|<12.5); the "
            "electron density relaxes. Both integrate to ~82 e (neutral). Loaded via the "
            "canonical `load_vti` (physical order, centre-checked — no fftshift)."),
        new_code_cell(CHARGE),
        new_markdown_cell(
            "## 3. Does the classical excess reach zero? — extended-r sweep (Lz=200)\n"
            "The Lz=120 sweep only reached r=40 (classical excess still ~12 eV there). "
            "This extends r to 76 Bohr in a larger open-z box with its own converged GS "
            "(`runs/h0_p2_far/gs_p2_lz200`, E_GS = 60.25 Ha), full components streamed and "
            "sum-checked. Excess = E_tot(0) - E_GS(Lz=200) for both projectiles; the dotted "
            "line marks the WP zero-point 3/4σ². The printed line reports where the classical "
            "excess is minimal."),
        new_code_cell(EXTENDED),
    ]
    p = OUT / "H0_p2_interpretation.ipynb"
    nbf.write(nb(cells), str(p)); print("wrote", p.name); return p

if __name__ == "__main__":
    build()
    print("execute: python3 -m nbconvert --to notebook --execute --inplace H0_p2_interpretation.ipynb (venv)")
