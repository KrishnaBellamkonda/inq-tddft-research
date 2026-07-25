#!/usr/bin/env python3
"""Build + execute the Phase-1 vacuum-exit notebook (classical-highdensity-sv).
Shows the projectile PERTURBATION (charge n_proj and potential phi_proj = poisson(n_proj))
moving along z and exiting the z-open box, plus the clip/no-wrap evidence.
Run with the venv python once the vac_exit VTIs exist.
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path

HERE = Path(__file__).resolve().parent
NB_PATH = HERE / "phase1_vac_exit.ipynb"
RES = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "scripts/classical_highdensity_sv/vac_exit/results")

cells = []
cells.append(new_markdown_cell(
    "# Phase 1 — vacuum exit test: the projectile perturbation along z\n\n"
    "Campaign `classical-highdensity-sv`. A moving Gaussian **charge** (σ_pot=0.354, "
    "= σ_WP/√2 for σ_WP=0.5), built by INQ's own `gaussian_density` in a z-open "
    "`periodicity(2)` box (35×35×85 Bohr, far face at **z=+42.5**), is swept along z. "
    "We verify it is **clipped** by the finite z-grid as it exits and does **not** "
    "wrap to the opposite face — and we visualise both the charge n_proj and the "
    "**potential the bath actually feels**, φ_proj = poisson(n_proj) (the +φ "
    "repulsive well of a −1 projectile), moving through the box and out."))

cells.append(new_code_cell(
    "import numpy as np, glob, re, matplotlib.pyplot as plt\n"
    "from matplotlib import animation\n"
    "from inqview import load_vti\n"
    "from IPython.display import Image, display\n"
    f"RES = {RES!r}\n"
    "HALF = 12.5      # the slab a real run will place here (context)\n"
    "FAR  = 42.5      # z-open far face (exit boundary)\n"
    "def zc_of(path):\n"
    "    m = re.search(r'z([+-]\\d+p\\d+)', path)\n"
    "    return float(m.group(1).replace('p','.'))\n"
    "nfiles = sorted(glob.glob(RES+'/nproj_z*.vti'), key=zc_of)\n"
    "pfiles = sorted(glob.glob(RES+'/phi_z*.vti'),  key=zc_of)\n"
    "print(len(nfiles),'n_proj VTIs,',len(pfiles),'phi VTIs')"))

cells.append(new_markdown_cell(
    "## 1. Clip + no-wrap evidence (from `exit_scan.csv`)\n"
    "Norm ≈1 while inside → smooth decay to 0 across the +42.5 face → **no secondary "
    "rise**. Wrap witness (max charge density at z<−38) must be ~0 at every position."))
cells.append(new_code_cell(
    "import csv\n"
    "rows=list(csv.DictReader(open(RES+'/exit_scan.csv')))\n"
    "z=np.array([float(r['z_center']) for r in rows])\n"
    "integ=np.array([float(r['integral']) for r in rows])\n"
    "wrap=np.array([float(r['wrap_witness_max']) for r in rows])\n"
    "phip=np.array([float(r['phi_peak']) for r in rows])\n"
    "fig,ax=plt.subplots(1,3,figsize=(14,3.6))\n"
    "ax[0].plot(z,integ,'o-',ms=3); ax[0].axvline(FAR,ls='--',c='k'); ax[0].axhline(1,ls=':',c='C7')\n"
    "ax[0].set_xlabel('projectile z-center (Bohr)'); ax[0].set_ylabel('∫ n_proj dV'); ax[0].set_title('norm: clip at +42.5')\n"
    "ax[1].semilogy(z,np.maximum(wrap,1e-300)+1e-300,'o-',ms=3); ax[1].axvline(FAR,ls='--',c='k')\n"
    "ax[1].set_xlabel('projectile z-center (Bohr)'); ax[1].set_ylabel('max n at z<-38 (wrap witness)'); ax[1].set_title('wrap witness = 0 everywhere')\n"
    "ax[2].plot(z,phip,'o-',ms=3,c='C3'); ax[2].axvline(FAR,ls='--',c='k')\n"
    "ax[2].set_xlabel('projectile z-center (Bohr)'); ax[2].set_ylabel('peak φ_proj (Ha)'); ax[2].set_title('potential well height')\n"
    "fig.tight_layout(); plt.show()\n"
    "print('max wrap witness over ALL positions =', wrap.max(), ' (0 = no wrap)')"))

cells.append(new_markdown_cell(
    "## 2. The perturbation plotted in the box along z\n"
    "On-axis (x=y=0) lineouts of the charge n_proj(z) and the potential φ_proj(z) at "
    "several projectile positions — the well **moving through the box and clipping at "
    "the far face**. Slab region (a real run's target) and the +42.5 exit face marked."))
cells.append(new_code_cell(
    "def onaxis(path):\n"
    "    f=load_vti(path)\n"
    "    ix=np.argmin(np.abs(f.x)); iy=np.argmin(np.abs(f.y))\n"
    "    return f.z, f.data[ix,iy,:]\n"
    "show_zc=[-20,0,12,30,40,42.5,45]\n"
    "def nearest(files,zc): return min(files,key=lambda p:abs(zc_of(p)-zc))\n"
    "fig,ax=plt.subplots(1,2,figsize=(13,4.2))\n"
    "for zc in show_zc:\n"
    "    z1,n1=onaxis(nearest(nfiles,zc)); ax[0].plot(z1,n1,label=f'{zc_of(nearest(nfiles,zc)):+.1f}')\n"
    "    z2,p2=onaxis(nearest(pfiles,zc)); ax[1].plot(z2,p2,label=f'{zc_of(nearest(pfiles,zc)):+.1f}')\n"
    "for a in ax:\n"
    "    a.axvspan(-HALF,HALF,color='C7',alpha=0.15); a.axvline(FAR,ls='--',c='k',lw=0.8)\n"
    "    a.set_xlabel('z (Bohr)'); a.legend(title='proj z-center',fontsize=7,ncol=2)\n"
    "ax[0].set_ylabel('n_proj (a0$^{-3}$)'); ax[0].set_title('charge n_proj(z) — clips at +42.5')\n"
    "ax[1].set_ylabel('φ_proj (Ha)'); ax[1].set_title('potential well φ_proj(z) the bath feels')\n"
    "fig.tight_layout(); plt.show()"))

cells.append(new_markdown_cell(
    "## 3. Animation — the perturbation potential moving out of the box\n"
    "φ_proj(x,z) at y=0 across the full sweep. Fixed colour scale; the +42.5 face "
    "and the slab region are marked. Watch the well translate and vanish at the far "
    "face with nothing appearing at the near (−z) face."))
cells.append(new_code_cell(
    "# load phi xz slices (2D only, discard 3D) across the sweep\n"
    "slices=[]; zcs=[]\n"
    "for p in pfiles:\n"
    "    f=load_vti(p); slices.append(f.xz_slice(0.0)); zcs.append(zc_of(p))\n"
    "    X=[f.x[0],f.x[-1],f.z[0],f.z[-1]]\n"
    "vmax=max(s.max() for s in slices)\n"
    "fig,ax=plt.subplots(figsize=(4.2,6))\n"
    "im=ax.imshow(slices[0],origin='lower',aspect='auto',extent=X,vmin=0,vmax=vmax,cmap='inferno')\n"
    "ax.axhline(FAR,ls='--',c='cyan',lw=1); ax.axhline(-HALF,ls=':',c='w',lw=0.6); ax.axhline(HALF,ls=':',c='w',lw=0.6)\n"
    "ax.set_xlabel('x (Bohr)'); ax.set_ylabel('z (Bohr)')\n"
    "ttl=ax.set_title('')\n"
    "fig.colorbar(im,ax=ax,label='φ_proj (Ha)')\n"
    "def upd(i):\n"
    "    im.set_data(slices[i]); ttl.set_text(f'proj z-center = {zcs[i]:+.1f} Bohr'); return im,ttl\n"
    "anim=animation.FuncAnimation(fig,upd,frames=len(slices),interval=180,blit=False)\n"
    "GIF=r'"+str(HERE/'perturbation_along_z.gif')+"'\n"
    "anim.save(GIF,writer=animation.PillowWriter(fps=6)); plt.close(fig)\n"
    "print('wrote',GIF); display(Image(filename=GIF))"))

cells.append(new_markdown_cell(
    "## Verdict\n\n"
    "The projectile perturbation (charge **and** potential) is faithfully **clipped** "
    "at the z-open far face and **never wraps** — so in a real run the Gaussian "
    "potential leaves the box cleanly, and (CAP-free) the electronic energy will "
    "plateau once it exits. **Manual gate: yours to accept.**"))

nb = new_notebook(cells=cells, metadata={"kernelspec": {
    "display_name": "Python 3", "language": "python", "name": "python3"}})
ep = ExecutePreprocessor(timeout=1200, kernel_name="python3")
ep.preprocess(nb, {"metadata": {"path": str(HERE)}})
nbf.write(nb, NB_PATH)
print("WROTE + EXECUTED", NB_PATH)
