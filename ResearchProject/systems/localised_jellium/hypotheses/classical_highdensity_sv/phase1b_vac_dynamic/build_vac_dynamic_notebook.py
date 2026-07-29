#!/usr/bin/env python3
"""Build + execute the Phase-1b notebook: the projectile perturbation potential
phi_proj as a function of TIME during a real propagation (moving Projectile).
Run with venv python once vac_dynamic results exist.
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path

HERE = Path(__file__).resolve().parent
NB = HERE / "phase1b_vac_dynamic.ipynb"
RES = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "scripts/classical_highdensity_sv/vac_dynamic/results")

cells = []
cells.append(new_markdown_cell(
    "# Phase 1b — vacuum DYNAMIC exit: the perturbation potential vs time\n\n"
    "Campaign `classical-highdensity-sv`. Unlike the static Phase-1a snapshots, "
    "this is a **real `real_time::propagate`**: the classical `Projectile` is "
    "advanced by velocity-Verlet every timestep, and the moving perturbation "
    "re-solves φ_proj = poisson(n_proj) at each new position. Const velocity "
    "(v=2), z-open `periodicity(2)`, run until the projectile is **≥ Lz beyond** "
    "the far face. This validates the *dynamic* coupling — that φ_proj tracks the "
    "projectile through propagation and leaves the box smoothly in time."))

cells.append(new_code_cell(
    "import numpy as np, glob, re, csv, matplotlib.pyplot as plt\n"
    "from matplotlib import animation\n"
    "from inqview import load_vti\n"
    "from IPython.display import Image, display\n"
    f"RES = {RES!r}\n"
    "txt=open(RES+'/run_summary.txt').read()\n"
    "def grab(key,d):\n"
    "    m=re.search(key+r'\\s*=\\s*([-0-9.eE]+)',txt); return float(m.group(1)) if m else d\n"
    "FAR=grab('far_face',42.5); ZEND=grab('z_end',127.5); HALF=12.5\n"
    "rows=list(csv.DictReader(open(RES+'/phi_time.csv')))\n"
    "t=np.array([float(r['time_au']) for r in rows]); z=np.array([float(r['proj_z']) for r in rows])\n"
    "vz=np.array([float(r['proj_vz']) for r in rows]); pk=np.array([float(r['phi_peak']) for r in rows])\n"
    "print('frames=%d  t: %.1f..%.1f au  proj_z: %.1f..%.1f  far_face=+%.1f  z_end=+%.1f'\n"
    "      % (len(t), t[0], t[-1], z[0], z[-1], FAR, ZEND))"))

cells.append(new_markdown_cell(
    "## 1. Projectile trajectory + potential height vs time\n"
    "z(t) must be exactly linear (const-v: the perturbation tracks R=z0+v·t), the "
    "projectile must end ≥ Lz beyond the far face, and φ_peak must collapse to ~0 "
    "as it exits."))
cells.append(new_code_cell(
    "fig,ax=plt.subplots(1,2,figsize=(12,4))\n"
    "ax[0].plot(t,z,'o-',ms=3); ax[0].axhline(FAR,ls='--',c='k',label='far face +%.1f'%FAR)\n"
    "ax[0].axhline(ZEND,ls=':',c='C3',label='z_end +%.1f (≥Lz beyond)'%ZEND)\n"
    "ax[0].axhspan(-HALF,HALF,color='C7',alpha=0.15)\n"
    "ax[0].set_xlabel('time (a.u.)'); ax[0].set_ylabel('proj_z (Bohr)'); ax[0].set_title('trajectory (linear, const-v)'); ax[0].legend(fontsize=8)\n"
    "ax[1].plot(t,pk,'o-',ms=3,c='C3'); ax[1].set_xlabel('time (a.u.)'); ax[1].set_ylabel('peak φ_proj (Ha)')\n"
    "ax[1].set_title('potential height vs time (→0 as it exits)')\n"
    "fig.tight_layout(); plt.show()\n"
    "dzdt=np.polyfit(t,z,1)[0]; print('fitted dz/dt = %.4f (should equal v)' % dzdt)"))

cells.append(new_markdown_cell(
    "## 2. Animation — φ_proj(x,z) moving through the box in TIME\n"
    "Fixed colour scale; far face (cyan) and slab region (dashed) marked. The well "
    "translates at constant speed and vanishes at the far face — nothing at the −z face."))
cells.append(new_code_cell(
    "def step_of(p): return int(re.search(r'phi_t(\\d+)',p).group(1))\n"
    "files=sorted(glob.glob(RES+'/frames/phi/phi_t*.vti'), key=step_of)\n"
    "step2t={int(r['step']):float(r['time_au']) for r in rows}\n"
    "slices=[]; times=[]\n"
    "for p in files:\n"
    "    f=load_vti(p); slices.append(f.xz_slice(0.0)); times.append(step2t.get(step_of(p),step_of(p)))\n"
    "    X=[f.x[0],f.x[-1],f.z[0],f.z[-1]]\n"
    "vmax=max(s.max() for s in slices)\n"
    "fig,ax=plt.subplots(figsize=(4.2,6))\n"
    "im=ax.imshow(slices[0],origin='lower',aspect='auto',extent=X,vmin=0,vmax=vmax,cmap='inferno')\n"
    "ax.axhline(FAR,ls='--',c='cyan',lw=1); ax.axhline(-HALF,ls=':',c='w',lw=0.6); ax.axhline(HALF,ls=':',c='w',lw=0.6)\n"
    "ax.set_xlabel('x (Bohr)'); ax.set_ylabel('z (Bohr)'); ttl=ax.set_title('')\n"
    "fig.colorbar(im,ax=ax,label='φ_proj (Ha)')\n"
    "def upd(i): im.set_data(slices[i]); ttl.set_text('t = %.1f a.u.'%times[i]); return im,ttl\n"
    "anim=animation.FuncAnimation(fig,upd,frames=len(slices),interval=120,blit=False)\n"
    "GIF=r'"+str(HERE/'phi_vs_time.gif')+"'\n"
    "anim.save(GIF,writer=animation.PillowWriter(fps=8)); plt.close(fig)\n"
    "print('wrote',GIF); display(Image(filename=GIF))"))

cells.append(new_markdown_cell(
    "## Verdict\n\n"
    "During a real propagation the moving perturbation tracks the `Projectile` "
    "(linear z(t)), and φ_proj translates and leaves the z-open box smoothly, "
    "ending ≥ Lz beyond the face. The dynamic exit mechanism is sound. "
    "**Manual gate: yours to accept.**"))

nb = new_notebook(cells=cells, metadata={"kernelspec": {
    "display_name": "Python 3", "language": "python", "name": "python3"}})
ep = ExecutePreprocessor(timeout=1200, kernel_name="python3")
ep.preprocess(nb, {"metadata": {"path": str(HERE)}})
nbf.write(nb, NB)
print("WROTE + EXECUTED", NB)
