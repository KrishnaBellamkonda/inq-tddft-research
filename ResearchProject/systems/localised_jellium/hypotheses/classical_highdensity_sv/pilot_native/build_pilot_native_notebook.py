#!/usr/bin/env python3
"""Phase-3 pilot Run B notebook: native-Ehrenfest ghost-UPF ion (contrast to Run A).
Documents that native Ehrenfest MOVES the z_valence=0 ghost, but the ghost's unscreened
1/r tail (+ no background compensation in INQ's ion-force, since the jellium background is
a perturbation not ions) makes the trajectory UNPHYSICAL in this geometry — so the
perturbation projectile is required. Plus a native-vs-perturbation z(t) overlay.
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path
import glob

HERE = Path(__file__).resolve().parent
NB = HERE / "pilot_native_run_notebook.ipynb"
NAT = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "scripts/classical_highdensity_sv/pilot_native/results/pilot_native/native.csv")
PERT = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "scripts/classical_highdensity_sv/pilot/results/pilot/raw/observables/projectile.csv")
frames = sorted(glob.glob("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
                          "scripts/classical_highdensity_sv/pilot_native/results/**/frames/**/*.vti", recursive=True))
FRAMEDIR = str(Path(frames[0]).parent) if frames else ""

cells = []
cells.append(new_markdown_cell(
    "# Phase 3 pilot — Run B: native-Ehrenfest ghost-UPF ion (contrast)\n\n"
    "Same slab as Run A, but a REAL mass-1 ghost-UPF ion moved by INQ's OWN native "
    "Ehrenfest (`ion_dynamics=EHRENFEST`), own GS with the ghost present (E_GS=217.3 Ha). "
    "Question: does native Ehrenfest move a z_valence=0 ghost, and does it reproduce Run A?"))

cells.append(new_code_cell(
    "import numpy as np, glob, re, csv, matplotlib.pyplot as plt\n"
    "from matplotlib import animation\n"
    "from inqview import load_vti\n"
    "from IPython.display import Image, display\n"
    f"NAT={NAT!r}; PERT={PERT!r}; FRAMEDIR={FRAMEDIR!r}\n"
    "HA=27.211386; HALF=12.5; FARFACE=42.5\n"
    "r=[x for x in csv.DictReader(open(NAT)) if all(x.get(k) not in (None,'') for k in ('time','z','vz','E_total'))]\n"
    "C=lambda k:np.array([float(x[k]) for x in r])\n"
    "t=C('time'); z=C('z'); vz=C('vz'); Et=C('E_total')\n"
    "print('native ghost: z %.1f -> min %.1f max %.1f -> %.1f ; vz %.2f->%.2f ; E swing %.0f eV'%(z[0],z.min(),z.max(),z[-1],vz[0],vz[-1],(Et.max()-Et.min())*HA))"))

cells.append(new_markdown_cell(
    "## Density evolution — the ghost never reaches the slab\n"
    "n(x,z,t), mid-y, physical order; slab faces dashed. The ghost stalls and oscillates "
    "in vacuum on the launch side."))
cells.append(new_code_cell(
    "if FRAMEDIR:\n"
    "    so=lambda p:int(re.search(r'_t?(\\d+)\\.vti',p).group(1)) if re.search(r'(\\d+)\\.vti',p) else 0\n"
    "    ff=sorted(glob.glob(FRAMEDIR+'/*.vti'),key=so); ff=ff[::max(1,len(ff)//50)]\n"
    "    f0=load_vti(ff[0]); base=f0.xz_slice(0.0); X=[f0.x[0],f0.x[-1],f0.z[0],f0.z[-1]]\n"
    "    sl=[load_vti(f).xz_slice(0.0) for f in ff]; vmax=np.percentile(np.abs(np.array(sl)),99.5)\n"
    "    fig,ax=plt.subplots(figsize=(4.5,6)); im=ax.imshow(sl[0],origin='lower',aspect='auto',extent=X,vmin=0,vmax=vmax,cmap='viridis')\n"
    "    ax.axhline(-HALF,ls='--',c='w'); ax.axhline(HALF,ls='--',c='w'); ax.set_xlabel('x'); ax.set_ylabel('z (Bohr)'); tt=ax.set_title('')\n"
    "    def upd(i): im.set_data(sl[i]); tt.set_text('frame %d/%d'%(i,len(ff))); return im,tt\n"
    "    an=animation.FuncAnimation(fig,upd,frames=len(ff),interval=120); GIF=r'"+str(HERE/'native_density.gif')+"'\n"
    "    an.save(GIF,writer=animation.PillowWriter(fps=8)); plt.close(fig); print('wrote',GIF); display(Image(filename=GIF))\n"
    "else: print('no native density frames')"))

cells.append(new_markdown_cell(
    "## Native-vs-perturbation trajectory — the decisive contrast\n"
    "Run A (perturbation, analytic force) transits the slab; Run B (native ghost) stalls "
    "in vacuum, reverses, and oscillates — it never reaches the slab, and E is not conserved."))
cells.append(new_code_cell(
    "pa=list(csv.DictReader(open(PERT))); paz=np.array([float(x['proj_z']) for x in pa]); pat=np.array([float(x['time_au']) for x in pa])\n"
    "fig,ax=plt.subplots(1,2,figsize=(13,4.5))\n"
    "ax[0].plot(pat,paz,label='Run A perturbation (analytic force)',c='C0')\n"
    "ax[0].plot(t,z,label='Run B native ghost-UPF ion',c='C3')\n"
    "ax[0].axhspan(-HALF,HALF,color='C7',alpha=.2,label='slab'); ax[0].axhline(FARFACE,ls=':',c='k')\n"
    "ax[0].set_xlabel('t (a.u.)'); ax[0].set_ylabel('proj_z (Bohr)'); ax[0].legend(fontsize=8); ax[0].set_title('trajectory: perturbation transits, native stalls in vacuum')\n"
    "ax[1].plot(t,(Et-Et[0])*HA,c='C3'); ax[1].set_xlabel('t'); ax[1].set_ylabel('ΔE_total (eV)'); ax[1].set_title('native E_total NOT conserved (%.0f eV swing)'%((Et.max()-Et.min())*HA))\n"
    "fig.tight_layout(); plt.show()\n"
    "print('native reached slab? %s (max z=%.1f, near face -12.5)'%(z.max()>-HALF, z.max()))"))

cells.append(new_markdown_cell(
    "## Verdict\n\n"
    "- **Native Ehrenfest DOES move the z_valence=0 ghost ion** (it feels its local HF "
    "force) — settles that open question.\n"
    "- **But the native ghost is DISQUALIFIED for this jellium geometry.** Its unscreened "
    "1/r Coulomb tail + the fact that INQ's ion-force sees the slab ELECTRONS but not the "
    "positive background (a perturbation, not ions) gives a large spurious static vacuum "
    "force → the ghost stalls ~14 Bohr short of the slab, reverses, oscillates, and E is "
    "not conserved (229 eV swing).\n"
    "- **Run A (perturbation, analytic force) is correct and REQUIRED:** its drag is "
    "(electrons − background), so the neutral slab exerts ~no static force in vacuum and "
    "the projectile transits. In a *compatible* setup (Test C: finite box, no background "
    "perturbation) native and perturbation agreed to 0.11% — so this is a setup "
    "incompatibility, not a native-Ehrenfest flaw. The perturbation projectile stands."))

nb = new_notebook(cells=cells, metadata={"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"}})
ep = ExecutePreprocessor(timeout=1800, kernel_name="python3"); ep.preprocess(nb, {"metadata":{"path":str(HERE)}})
nbf.write(nb, NB); print("WROTE + EXECUTED", NB)
