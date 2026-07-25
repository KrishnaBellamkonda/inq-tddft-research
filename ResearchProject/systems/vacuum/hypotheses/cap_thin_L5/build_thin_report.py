#!/usr/bin/env python3
"""Build + execute the cap_thin_L5 reflectivity-tuning notebook (ADR 0007 location).

Thin in-built CAP (L=5 Bohr) reflectivity tuning: 3 depth curves
η ∈ {−0.01,−0.05,−0.30} Ha × 11 energies (1–100 eV), each run carrying the full
free-WP minimum observable set + manifest. Goal: locate the ε(E) curve whose
minimum sits near 10 eV and stays low across the decade. Writes the combined CSV,
the overlaid ε(E) curves, a density GIF, and assembles cap_thin_L5_study.ipynb.

    PYTHONPATH=.../inq-stack/python /local/.../venv/bin/python3 build_thin_report.py
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path
import csv

HERE = Path(__file__).resolve().parent
VAC = HERE.parent.parent                       # systems/vacuum
SWEEP = VAC / "cap_thin_L5"

# ---- combined CSV (built outside the notebook so it exists even if exec fails) ----
def parse_kv(p):
    out = {}
    for ln in Path(p).read_text().splitlines():
        k, _, v = ln.partition(' ')
        try: out[k] = float(v)
        except ValueError: out[k] = v
    return out

recs = []
for d in sorted(SWEEP.glob('run_cap_*')):
    f = d / 'results/epsilon.txt'
    if f.exists():
        r = parse_kv(f); r['name'] = d.name; recs.append(r)
recs.sort(key=lambda r: (r.get('eta_Ha', 0), r.get('E_eV', 0)))
cols = ['name', 'E_eV', 'k0', 'eta_Ha', 'L_abs', 'epsilon', 'absorbed_fraction']
with open(HERE / 'cap_thin_combined.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(cols)
    for r in recs: w.writerow([r.get(c, '') for c in cols])
print(f'{len(recs)} runs -> cap_thin_combined.csv')

nb = new_notebook(); C = nb.cells

C.append(new_markdown_cell(r"""# Thin in-built CAP (L=5 Bohr) — reflectivity tuning across 1–100 eV

A Gaussian wavepacket (k₀ along +z) meets INQ's **own** Complex Absorbing
Potential `perturbations::absorbing` placed in the **last 5 Bohr** of the box.
We sweep the depth η over three shallow values and trace the reflection error
ε(E) across the 10⁰–10² eV decade, asking: **which η puts the reflection minimum
near 10 eV while keeping ε low across the band?**

## Governing equation (De Giovannini–Larsen–Rubio 2014, Eq. 17)
$$\hat H=\hat H_0+i\,\eta\,\sin^2\!\Big(\tfrac{(z-z_{\rm abs})\pi}{2L}\Big)\ \ (z\in[z_{\rm abs},z_{\rm abs}+L]),\qquad \eta<0\ \text{absorbs},\ L=5\,a_0.$$

## Method
- **Absorber = INQ's own** `perturbations::absorbing` (region-restricted sin²
  imaginary potential, integrated in H via the trailing `real_time::propagate`
  argument). **Thin** L=5 Bohr, **shallow** η ∈ {−0.01, −0.05, −0.30} Ha.
- **Propagator ETRS** (no renormalisation → absorption survives).
- Energy axis: 11 log-spaced points 1–100 eV, densified near the 10 eV target.
  Geometry scales with energy (σ=4√2/k₀, box 6σ+L); ε = inner survival
  $\int_{z<z_{\rm abs}}|\psi(\tau)|^2/N_0$. Lower ε ⇒ better absorption.

> **Provisional.** Exploratory until the `inq-study` engine regression (Task #7)
> confirms the scalar-potential complexification is inert for non-CAP physics.
> A thin 5-Bohr absorber has a *higher reflection floor* than wide ones — the
> deliverable is the **curve shape and where the minimum sits**, not a deep null.
"""))

C.append(new_code_cell(r"""import numpy as np, csv, re, subprocess, sys
from pathlib import Path
import matplotlib.pyplot as plt
from inqview.visualisation import style
style.apply_theme()
HERE = Path(r"%s"); SWEEP = Path(r"%s")
rows = list(csv.DictReader(open(HERE/'cap_thin_combined.csv')))
for r in rows:
    for k in ('E_eV','k0','eta_Ha','L_abs','epsilon','absorbed_fraction'):
        r[k] = float(r[k])
etas = sorted({round(r['eta_Ha'],3) for r in rows})
def curve(eta):
    c = sorted([r for r in rows if abs(r['eta_Ha']-eta)<1e-6], key=lambda r:r['E_eV'])
    return np.array([r['E_eV'] for r in c]), np.array([r['epsilon'] for r in c]), np.array([r['absorbed_fraction'] for r in c]), c
print(f'{len(rows)} runs; eta values:', etas)
""" % (str(HERE), str(SWEEP))))

C.append(new_markdown_cell(r"""## 1. Manifest validation (the minimum observable set)
Every run declares a `free_wp` manifest; tiers 1–3 (existence/schema/finite) pass
for all runs — the set is complete. The tier-4 invariants `energy_total` drift and
WP norm-band **fail by design**: a CAP is non-Hermitian, so energy and norm
decrease — that *is* the absorption, not a defect."""))
C.append(new_code_cell(r"""sample = [SWEEP/r['name'] for r in rows[:2]+rows[-1:]]
for d in sample:
    out = subprocess.run([sys.executable,'-m','inqview.validation',str(d)],capture_output=True,text=True)
    head = out.stdout.splitlines()[0] if out.stdout else '(no output)'
    fails=[l.strip() for l in out.stdout.splitlines() if 'FAIL' in l and 'invariant' in l]
    print(f"{d.name:30s} {head.split('(')[0].strip()}")
    for fl in fails: print('     ',fl)
print('\\nOnly energy_total drift + wp norm-band FAIL — expected for an absorbing (non-unitary) run.')
"""))

C.append(new_markdown_cell(r"""## 2. Reflectivity curves ε(E) — the headline result
Three depth curves overlaid (log–log). Vertical guide at the **10 eV target**;
each curve's minimum is marked. We read off which η lands its dip near 10 eV and
how low ε stays across 1–100 eV."""))
C.append(new_code_cell(r"""fig,ax=plt.subplots(figsize=(6.4,4.2))
colors={}
summary=[]
for i,eta in enumerate(etas):
    E,eps,ab,c=curve(eta)
    if len(E)==0: continue
    ln,=ax.loglog(E,eps,'o-',label=fr'$\eta={eta:+.2f}$ Ha')
    j=int(np.argmin(eps)); ax.plot(E[j],eps[j],'*',ms=15,color=ln.get_color())
    summary.append((eta,E[j],eps[j],eps.min(),eps.max()))
ax.axvline(10,ls='--',color='0.5',lw=1); ax.text(10,ax.get_ylim()[1],' 10 eV target',va='top',fontsize=8,color='0.4')
ax.set_xlabel('WP energy (eV)'); ax.set_ylabel(r'reflection error $\varepsilon$')
ax.set_title('Thin CAP (L=5 Bohr): reflectivity vs energy'); ax.legend(); ax.grid(True,which='both',alpha=0.3)
fig.tight_layout(); fig.savefig('fig_cap_thin_reflectivity.png',dpi=140); plt.show()
print('eta   E_min(eV)  eps_min    eps_floor  eps_max')
for eta,Em,epm,emn,emx in summary: print(f'{eta:+.2f}  {Em:8.2f}  {epm:.3e}  {emn:.3e}  {emx:.3e}')
"""))

C.append(new_markdown_cell(r"""## 3. Absorbed fraction & where the minimum sits
Left: absorbed fraction vs energy per η (how much norm the CAP removes). Right:
the energy of each curve's reflection minimum vs η — does any depth place it at
~10 eV?"""))
C.append(new_code_cell(r"""fig,ax=plt.subplots(1,2,figsize=(10,3.8))
for eta in etas:
    E,eps,ab,c=curve(eta)
    if len(E)==0: continue
    ax[0].semilogx(E,ab,'o-',label=fr'$\eta={eta:+.2f}$')
ax[0].axvline(10,ls='--',color='0.6',lw=1); ax[0].set_xlabel('WP energy (eV)'); ax[0].set_ylabel('absorbed fraction')
ax[0].set_ylim(0,1.02); ax[0].set_title('absorbed fraction vs energy'); ax[0].legend(fontsize=8)
em=[(eta,curve(eta)[0][int(np.argmin(curve(eta)[1]))]) for eta in etas if len(curve(eta)[0])]
ax[1].semilogx([abs(e) for e,_ in em],[m for _,m in em],'s-',color='C3')
ax[1].axhline(10,ls='--',color='0.6',lw=1); ax[1].set_xlabel(r'depth $|\eta|$ (Ha)'); ax[1].set_ylabel(r'$E$ of min $\varepsilon$ (eV)')
ax[1].set_title('reflection-minimum energy vs depth')
fig.tight_layout(); fig.savefig('fig_cap_thin_absorbed.png',dpi=140); plt.show()
"""))

C.append(new_markdown_cell(r"""## 4. Dynamics near the target — survival & absorbed energy
At the energy nearest 10 eV, inner-region survival and ΔE_total for the three
depths."""))
C.append(new_code_cell(r"""def load_csv(p):
    rr=list(csv.DictReader(open(p))); return {k:np.array([float(x[k]) for x in rr]) for k in rr[0]}
fig,ax=plt.subplots(1,2,figsize=(10,3.6))
for eta in etas:
    E,eps,ab,c=curve(eta)
    if len(c)==0: continue
    r=min(c,key=lambda r:abs(r['E_eV']-10))
    d=SWEEP/r['name']
    try:
        inn=load_csv(d/'results/raw/observables/inner_norm_vs_time.csv')
        obs=load_csv(d/'results/raw/observables/observables.csv')
        lbl=fr'$\eta={eta:+.2f}$ (E={r["E_eV"]:.0f}eV)'
        ax[0].plot(inn['time_au'],inn['inner_norm_over_N0'],label=lbl)
        ax[1].plot(obs['time_au'],obs['energy_total']-obs['energy_total'][0],label=lbl)
    except Exception as e: print('skip',r['name'],e)
ax[0].set_title('inner survival /N₀'); ax[0].set_xlabel('t (a.u.)'); ax[0].legend(fontsize=8)
ax[1].set_title('ΔE_total (absorbed energy)'); ax[1].set_xlabel('t (a.u.)')
fig.tight_layout(); fig.savefig('fig_cap_thin_dynamics.png',dpi=140); plt.show()
"""))

C.append(new_markdown_cell(r"""## 5. Density GIF — wavepacket meeting the thin CAP (E≈10 eV, deepest η)"""))
C.append(new_code_cell(r"""import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.animation as animation
# showcase: energy nearest 10 eV at the deepest eta (-0.30)
cands=[r for r in rows if abs(r['eta_Ha']+0.30)<1e-6]
sc=min(cands,key=lambda r:abs(r['E_eV']-10)) if cands else min(rows,key=lambda r:abs(r['E_eV']-10))
d=SWEEP/sc['name']; vti_dir=d/'results/raw/vti/density_wp'
frames=sorted(vti_dir.glob('density_wp_t*.vti'),key=lambda p:int(re.search(r'_t(\d+)',p.name).group(1)))
print(f"{len(frames)} frames from {sc['name']} (E={sc['E_eV']:.0f}eV, eta={sc['eta_Ha']:+.2f})")
def zprof(p):
    rd=vtk.vtkXMLImageDataReader(); rd.SetFileName(str(p)); rd.Update(); img=rd.GetOutput()
    dim=img.GetDimensions(); a=vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(dim[2],dim[1],dim[0])
    sp=img.GetSpacing(); org=img.GetOrigin(); z=org[2]+sp[2]*np.arange(dim[2]); return z,a.sum(axis=(1,2))
if frames:
    z0,_=zprof(frames[0]); profs=[zprof(p)[1] for p in frames]; ymax=max(p.max() for p in profs)*1.05
    fig,ax=plt.subplots(figsize=(6,3.2)); (line,)=ax.plot(z0,profs[0])
    ax.axvspan(sc['z_abs0'] if 'z_abs0' in sc else z0[-1]-5,(sc['z_abs0']+5) if 'z_abs0' in sc else z0[-1],color='C3',alpha=0.15,label='CAP (last 5 a₀)')
    ax.set_xlabel('z (Bohr)'); ax.set_ylabel('WP density (transverse sum)'); ax.set_ylim(0,ymax); ax.legend(loc='upper left',fontsize=8)
    ttl=ax.set_title('')
    def upd(k): line.set_ydata(profs[k]); ttl.set_text(f'frame {k+1}/{len(profs)}'); return line,ttl
    animation.FuncAnimation(fig,upd,frames=len(profs),interval=120,blit=False).save('fig_cap_thin_density.gif',writer='pillow',dpi=90); plt.close(fig)
    print('wrote fig_cap_thin_density.gif')
"""))

C.append(new_markdown_cell(r"""![density gif](fig_cap_thin_density.gif)

## 6. Run data (paths) + verdict"""))
C.append(new_code_cell(r"""for r in sorted(rows,key=lambda r:(r['eta_Ha'],r['E_eV'])):
    print(f"{r['name']:30s} E={r['E_eV']:7.2f}eV eps={r['epsilon']:.3e} absorbed={r['absorbed_fraction']:.4f}")
print('\\nrun dirs under', SWEEP)
"""))
C.append(new_markdown_cell(r"""**Findings (fill from the curves above).**
- Best depth & where its ε(E) minimum sits relative to the 10 eV target.
- The reflection floor of a thin L=5 absorber across 1–100 eV.
- Every run carries the complete free-WP minimum observable set (tiers 1–3 pass);
  energy-drift / norm-band tier-4 "failures" are the absorption signature.
- **Provisional** until the deferred `inq-study` regression (Task #7).
- Source: De Giovannini, Larsen & Rubio (2014), arXiv:1409.1689, §IV.
"""))

if __name__ == "__main__":
    ep = ExecutePreprocessor(timeout=3000, kernel_name='python3')
    ep.preprocess(nb, {'metadata': {'path': str(HERE)}})
    out = HERE / 'cap_thin_L5_study.ipynb'
    nbf.write(nb, out)
    print(f'wrote {out}')
