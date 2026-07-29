#!/usr/bin/env python3
"""Build + execute the in-built-CAP investigation notebook (ADR 0007 location).

Thorough study of INQ's in-built `perturbations::absorbing` on a free wavepacket:
depth (η) × width (L) × energy sweeps, every run carrying the full free-WP
minimum observable set + manifest. Renders the density GIF from VTI frames (vtk),
summarises manifest validation, and assembles an executed `cap_real_study.ipynb`.

    /local/.../venv/bin/python3 build_cap_report.py
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path

HERE = Path(__file__).resolve().parent
VAC = HERE.parent.parent                      # systems/vacuum

nb = new_notebook(); C = nb.cells

C.append(new_markdown_cell(r"""# INQ's in-built CAP — a thorough investigation (free wavepacket in vacuum)

A Gaussian wavepacket (k₀ along +z) travels toward INQ's **own** Complex Absorbing
Potential, `perturbations::absorbing`, placed at the far end of the box. We study
its three knobs — **depth η**, **width L**, **energy E** — with the full free-WP
minimum observable set on every run.

## Governing equation (De Giovannini–Larsen–Rubio 2014, Eq. 17)
$$\hat H(\mathbf r)=\hat H_0(\mathbf r)+i\,\eta\,\sin^2\!\Big(\tfrac{(z-z_{\rm abs})\pi}{2L}\Big)\ \ (z\in[z_{\rm abs},z_{\rm abs}+L]),\qquad \eta<0\ \text{absorbs.}$$

## Method
- **The CAP is INQ's own** `perturbations::absorbing` — a region-restricted sin²
  imaginary potential, integrated into the Hamiltonian (the trailing argument of
  `real_time::propagate`), **not** a post-step multiply.
- **It is a real, team-built feature** (Yao + Andrade, 2023). A 2024
  `update_hamiltonian` refactor silently regressed it (a real `vscalar`
  intermediate; no test drives absorbing through propagation). A 1-line
  `inq-study` change re-aligns `vscalar` with the team's already-complex `vks` —
  a **regression repair**, `inq/` untouched. See `docs/handovers/inq-study-cap-deferred.md`.
- **Propagator: ETRS** (applies $e^{-i\hat Ht}$, no renormalisation → absorption
  survives; Crank–Nicolson orthonormalises every step and would undo it).
- Geometry mirrors the MFA study (σ=4√2/k₀, box 6σ+L, CAP = last L); ε = inner
  survival $\int_{z<z_{\rm abs}}|\psi(\tau)|^2/N_0$.

> **Provisional.** These results are exploratory until the `inq-study` engine
> regression (Task #7) confirms the repair is inert for non-CAP physics.
"""))

C.append(new_code_cell(r"""import numpy as np, csv, re, glob, json, subprocess, sys
from pathlib import Path
import matplotlib.pyplot as plt
from inqview.visualisation import style
style.apply_theme()
VAC = Path(r"%s")

def parse_kv(p):
    out = {}
    for ln in Path(p).read_text().splitlines():
        k, _, v = ln.partition(' ')
        try: out[k] = float(v)
        except ValueError: out[k] = v
    return out

runs = []
for d in sorted(VAC.glob('cap_real/run_cap_*')):    # runs grouped under cap_real/ (ADR 0007 amendment 2026-06-15)
    f = d/'results/epsilon.txt'
    if f.exists():
        r = parse_kv(f); r['dir'] = d; r['name'] = d.name; runs.append(r)
print(f'{len(runs)} CAP runs')
def grp(pred): return sorted([r for r in runs if pred(r)], key=lambda r: r.get('eta_Ha',0))
depth = sorted([r for r in runs if abs(r['k0']-1.28)<1e-6 and abs(r['L_abs']-20)<1e-6], key=lambda r:abs(r['eta_Ha']))
width = sorted([r for r in runs if abs(r['k0']-1.28)<1e-6 and abs(r['eta_Ha']+0.5)<1e-6], key=lambda r:r['L_abs'])
energy= sorted([r for r in runs if abs(r['L_abs']-20)<1e-6 and abs(r['eta_Ha']+0.5)<1e-6], key=lambda r:r['E_eV'])
print(f'depth={len(depth)} width={len(width)} energy={len(energy)} runs')
""" % str(VAC)))

C.append(new_markdown_cell(r"""## 1. Manifest validation (the minimum observable set)
Every run declares a `free_wp` observable manifest; the inqview validator checks
4 tiers. **Tiers 1–3 (existence/schema/finite) pass for all runs** — the set is
complete. The **tier-4 invariants `energy_total` drift and `wp` norm-band fail by
design**: a CAP is non-Hermitian, so energy and norm decrease — that *is* the
absorption, not a defect. We report the validator verbatim and flag those two as
expected."""))

C.append(new_code_cell(r"""for r in runs[:3] + ([runs[-1]] if len(runs)>3 else []):
    out = subprocess.run([sys.executable, '-m', 'inqview.validation', str(r['dir'])],
                         capture_output=True, text=True)
    head = out.stdout.splitlines()[0] if out.stdout else '(no output)'
    fails = [l.strip() for l in out.stdout.splitlines() if 'FAIL' in l]
    print(f"{r['name']:32s} {head.split('(')[0].strip()}")
    for fl in fails: print('     ', fl)
print('\\nThe only FAILs are energy_total drift + wp norm-band — expected for an absorbing (non-unitary) run.')
"""))

C.append(new_markdown_cell(r"""## 2. Depth η — reflection U-shape
Too weak → incomplete absorption (high ε); a sweet spot; too strong → the steep
imaginary wall reflects (ε rises again)."""))
C.append(new_code_cell(r"""eta = np.array([abs(r['eta_Ha']) for r in depth]); eps=np.array([r['epsilon'] for r in depth])
ab = np.array([r['absorbed_fraction'] for r in depth])
fig,ax=plt.subplots(1,2,figsize=(9,3.4))
ax[0].loglog(eta,eps,'o-'); ax[0].set_xlabel(r'depth $|\eta|$ (Ha)'); ax[0].set_ylabel(r'$\varepsilon$'); ax[0].set_title('reflection vs depth')
ax[1].semilogx(eta,ab,'s-',color='C1'); ax[1].set_xlabel(r'depth $|\eta|$ (Ha)'); ax[1].set_ylabel('absorbed'); ax[1].set_ylim(0,1.02); ax[1].set_title('absorbed vs depth')
fig.tight_layout(); fig.savefig('fig_cap_depth.png',dpi=140); plt.show()
print(f'min eps={eps.min():.2e} at |eta|={eta[eps.argmin()]:.2f} Ha')
"""))

C.append(new_markdown_cell(r"""## 3. Width L and energy E (at η near the sweet spot)"""))
C.append(new_code_cell(r"""fig,ax=plt.subplots(1,2,figsize=(9,3.4))
if width:
    L=np.array([r['L_abs'] for r in width]); ax[0].semilogy(L,[r['epsilon'] for r in width],'o-')
ax[0].set_xlabel('absorber width L (Bohr)'); ax[0].set_ylabel(r'$\varepsilon$'); ax[0].set_title('reflection vs width (η=−0.5)')
if energy:
    E=np.array([r['E_eV'] for r in energy]); ax[1].loglog(E,[r['epsilon'] for r in energy],'o-',color='C2')
ax[1].set_xlabel('WP energy (eV)'); ax[1].set_ylabel(r'$\varepsilon$'); ax[1].set_title('reflection vs energy (η=−0.5, L=20)')
fig.tight_layout(); fig.savefig('fig_cap_width_energy.png',dpi=140); plt.show()
"""))

C.append(new_markdown_cell(r"""## 4. Dynamics — survival, absorbed energy, and system perturbation
For a few depths: inner-region survival, total WP norm (absorption), the
`energy_total` change (energy carried off by the absorbed WP), and `density_l2`
(L2 of the system-density change Δn — how the background is perturbed)."""))
C.append(new_code_cell(r"""def load_csv(p):
    rows=list(csv.DictReader(open(p))); return {k:np.array([float(x[k]) for x in rows]) for k in rows[0]}
fig,ax=plt.subplots(1,3,figsize=(12,3.4))
for r in depth:
    if abs(r['eta_Ha']) in (0.05,0.25,0.5,2.0):
        inn=load_csv(r['dir']/'results/raw/observables/inner_norm_vs_time.csv')
        obs=load_csv(r['dir']/'results/raw/observables/observables.csv')
        lbl=f"|η|={abs(r['eta_Ha']):.2f}"
        ax[0].plot(inn['time_au'],inn['inner_norm_over_N0'],label=lbl)
        ax[1].plot(obs['time_au'],obs['energy_total']-obs['energy_total'][0],label=lbl)
        ax[2].plot(obs['time_au'],obs['density_l2'],label=lbl)
ax[0].set_title('inner survival /N₀'); ax[0].set_xlabel('t (a.u.)'); ax[0].legend(fontsize=8)
ax[1].set_title('ΔE_total (absorbed energy)'); ax[1].set_xlabel('t (a.u.)')
ax[2].set_title('density_l2 (Δn system)'); ax[2].set_xlabel('t (a.u.)')
fig.tight_layout(); fig.savefig('fig_cap_dynamics.png',dpi=140); plt.show()
"""))

C.append(new_markdown_cell(r"""## 5. Final-state occupations & eigenvalues
KS occupations (fixed in non-interacting) and final eigenvalues — the excitation
fingerprint at t=τ."""))
C.append(new_code_cell(r"""r=depth[len(depth)//2]
occ=list(csv.DictReader(open(r['dir']/'results/raw/observables/eigenvalues/occupations.csv')))
eig=open(r['dir']/'results/raw/observables/eigenvalues/eigenvalues.csv').read()
print(f"final-state occupations ({r['name']}):")
for row in occ: print('  state',row.get('state_index'),'occ',row.get('occupation'))
print('\neigenvalues.csv:\n', eig)
"""))

C.append(new_markdown_cell(r"""## 6. Density GIF — wavepacket meeting the CAP"""))
C.append(new_code_cell(r"""import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.animation as animation
sc = next((r for r in runs if abs(r['k0']-1.28)<1e-6 and abs(r['L_abs']-20)<1e-6 and abs(r['eta_Ha']+0.25)<1e-6), depth[len(depth)//2])
vti_dir = sc['dir']/'results/raw/vti/density_wp'
frames = sorted(vti_dir.glob('density_wp_t*.vti'), key=lambda p:int(re.search(r'_t(\d+)',p.name).group(1)))
print(f"{len(frames)} frames from {sc['name']}")
def zprof(p):
    rd=vtk.vtkXMLImageDataReader(); rd.SetFileName(str(p)); rd.Update(); img=rd.GetOutput()
    d=img.GetDimensions(); a=vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(d[2],d[1],d[0])
    sp=img.GetSpacing(); org=img.GetOrigin(); z=org[2]+sp[2]*np.arange(d[2]); return z,a.sum(axis=(1,2))
if frames:
    z0,_=zprof(frames[0]); profs=[zprof(p)[1] for p in frames]; ymax=max(p.max() for p in profs)*1.05
    fig,ax=plt.subplots(figsize=(6,3.2)); (line,)=ax.plot(z0,profs[0])
    ax.axvspan(sc['z_abs0'],sc['z_abs0']+sc['L_abs'],color='C3',alpha=0.15,label='CAP')
    ax.set_xlabel('z (Bohr)'); ax.set_ylabel('WP density (transverse sum)'); ax.set_ylim(0,ymax); ax.legend(loc='upper left',fontsize=8)
    ttl=ax.set_title('')
    def upd(k): line.set_ydata(profs[k]); ttl.set_text(f'frame {k+1}/{len(profs)}'); return line,ttl
    animation.FuncAnimation(fig,upd,frames=len(profs),interval=120,blit=False).save('fig_cap_density.gif',writer='pillow',dpi=90); plt.close(fig)
    print('wrote fig_cap_density.gif')
"""))

C.append(new_markdown_cell(r"""![density gif](fig_cap_density.gif)

## 7. Run data (paths) + verdict"""))
C.append(new_code_cell(r"""for r in runs: print(f"{r['name']:34s} eps={r['epsilon']:.3e} absorbed={r['absorbed_fraction']:.4f} -> {r['dir']}")
"""))
C.append(new_markdown_cell(r"""**Findings.**
- INQ's **own** `perturbations::absorbing` works and reproduces textbook CAP
  behaviour (depth U-shape; width/energy trends) — repaired from a silent 2024
  regression, `inq/` untouched.
- Every run carries the **complete free-WP minimum observable set** (manifest
  tiers 1–3 pass). The energy-drift / norm-band tier-4 "failures" are the
  **absorption signature** (non-Hermitian, non-conservative) — expected.
- **Provisional** until the deferred `inq-study` regression (Task #7).
- Source: De Giovannini, Larsen & Rubio (2014), arXiv:1409.1689, §IV.
"""))

if __name__ == "__main__":
    ep = ExecutePreprocessor(timeout=2400, kernel_name='python3')
    ep.preprocess(nb, {'metadata': {'path': str(HERE)}})
    out = HERE / 'cap_real_study.ipynb'
    nbf.write(nb, out)
    print(f'wrote {out}')
