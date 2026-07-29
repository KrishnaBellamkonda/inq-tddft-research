#!/usr/bin/env python3
"""Total-density GIFs for the DIRECT-potential pilot: n(x,z,t) linear + log, and induced
Δn = n(t)−n(0). Mid-y x–z slice, physical order (load_vti), slab faces dashed. Embedded in
a small notebook + saved as standalone .gif."""
import nbformat as nbf, glob, re
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path

HERE = Path(__file__).resolve().parent
NB = HERE / "pilot_direct_density.ipynb"
FRAMES = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
          "scripts/classical_highdensity_sv/pilot_direct/results/pilot_direct/frames/total")
PRJ = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "scripts/classical_highdensity_sv/pilot_direct/results/pilot_direct/raw/observables/projectile.csv")

cells = []
cells.append(new_markdown_cell(
    "# Direct-potential pilot — total density n(x,z,t)\n\n"
    "Total electron density (mid-y x–z slice, physical order) as the mass-1 projectile "
    "(direct erf/r potential, v=2) transits the r_s=4.18 slab and exits the z-open box. "
    "Slab faces dashed, far face (+42.5) dotted."))
cells.append(new_code_cell(
    "import numpy as np, glob, re, csv, matplotlib.pyplot as plt\n"
    "from matplotlib import animation\n"
    "from matplotlib.colors import LogNorm\n"
    "from inqview import load_vti\n"
    "from IPython.display import Image, display\n"
    f"FRAMES={FRAMES!r}; PRJ={PRJ!r}\n"
    "HALF=12.5; FARFACE=42.5; dt=0.04\n"
    "so=lambda p:int(re.search(r'density_t(\\d+)',p).group(1))\n"
    "allf=sorted(glob.glob(FRAMES+'/*.vti'),key=so)\n"
    "files=allf[::max(1,len(allf)//60)]\n"
    "p=list(csv.DictReader(open(PRJ))); pt=np.array([float(x['time_au']) for x in p]); pz=np.array([float(x['proj_z']) for x in p])\n"
    "f0=load_vti(files[0]); X=[f0.x[0],f0.x[-1],f0.z[0],f0.z[-1]]\n"
    "tot=[load_vti(f).xz_slice(0.0) for f in files]; base=tot[0]\n"
    "times=[so(f)*dt for f in files]\n"
    "print('%d frames, t=%.1f..%.1f au'%(len(files),times[0],times[-1]))"))

cells.append(new_markdown_cell("## Total density — linear"))
cells.append(new_code_cell(
    "vmax=np.percentile(np.array(tot),99.8)\n"
    "fig,ax=plt.subplots(figsize=(4.6,6))\n"
    "im=ax.imshow(tot[0],origin='lower',aspect='auto',extent=X,vmin=0,vmax=vmax,cmap='viridis')\n"
    "ax.axhline(-HALF,ls='--',c='w',lw=.7); ax.axhline(HALF,ls='--',c='w',lw=.7); ax.axhline(FARFACE,ls=':',c='c',lw=.8)\n"
    "ax.set_xlabel('x (Bohr)'); ax.set_ylabel('z (Bohr)'); tt=ax.set_title(''); fig.colorbar(im,ax=ax,label='n')\n"
    "def upd(i): im.set_data(tot[i]); tt.set_text('t=%.1f au  proj_z=%.1f'%(times[i],np.interp(times[i],pt,pz))); return im,tt\n"
    "an=animation.FuncAnimation(fig,upd,frames=len(tot),interval=120); GIF=r'"+str(HERE/'total_density.gif')+"'\n"
    "an.save(GIF,writer=animation.PillowWriter(fps=8)); plt.close(fig); print('wrote',GIF); display(Image(filename=GIF))"))

cells.append(new_markdown_cell("## Total density — log scale (surface + spill-out visible)"))
cells.append(new_code_cell(
    "arr=np.array(tot); floor=max(arr[arr>0].min(),1e-6)\n"
    "fig,ax=plt.subplots(figsize=(4.6,6))\n"
    "im=ax.imshow(np.clip(tot[0],floor,None),origin='lower',aspect='auto',extent=X,norm=LogNorm(vmin=floor,vmax=vmax),cmap='magma')\n"
    "ax.axhline(-HALF,ls='--',c='w',lw=.7); ax.axhline(HALF,ls='--',c='w',lw=.7); ax.axhline(FARFACE,ls=':',c='c',lw=.8)\n"
    "ax.set_xlabel('x'); ax.set_ylabel('z (Bohr)'); tt=ax.set_title(''); fig.colorbar(im,ax=ax,label='n (log)')\n"
    "def upd(i): im.set_data(np.clip(tot[i],floor,None)); tt.set_text('t=%.1f au  proj_z=%.1f'%(times[i],np.interp(times[i],pt,pz))); return im,tt\n"
    "an=animation.FuncAnimation(fig,upd,frames=len(tot),interval=120); GIF=r'"+str(HERE/'total_density_log.gif')+"'\n"
    "an.save(GIF,writer=animation.PillowWriter(fps=8)); plt.close(fig); print('wrote',GIF); display(Image(filename=GIF))"))

cells.append(new_markdown_cell("## Induced density Δn = n(t) − n(0) (the bath response / anti-wake)"))
cells.append(new_code_cell(
    "ind=[s-base for s in tot]; vi=np.percentile(np.abs(np.array(ind)),99.5)\n"
    "fig,ax=plt.subplots(figsize=(4.6,6))\n"
    "im=ax.imshow(ind[0],origin='lower',aspect='auto',extent=X,vmin=-vi,vmax=vi,cmap='RdBu_r')\n"
    "ax.axhline(-HALF,ls='--',c='k',lw=.7); ax.axhline(HALF,ls='--',c='k',lw=.7); ax.axhline(FARFACE,ls=':',c='c',lw=.8)\n"
    "ax.set_xlabel('x'); ax.set_ylabel('z (Bohr)'); tt=ax.set_title(''); fig.colorbar(im,ax=ax,label='Δn')\n"
    "def upd(i): im.set_data(ind[i]); tt.set_text('t=%.1f au  proj_z=%.1f'%(times[i],np.interp(times[i],pt,pz))); return im,tt\n"
    "an=animation.FuncAnimation(fig,upd,frames=len(ind),interval=120); GIF=r'"+str(HERE/'induced_density.gif')+"'\n"
    "an.save(GIF,writer=animation.PillowWriter(fps=8)); plt.close(fig); print('wrote',GIF); display(Image(filename=GIF))"))

nb=new_notebook(cells=cells,metadata={"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"}})
ep=ExecutePreprocessor(timeout=1200,kernel_name="python3"); ep.preprocess(nb,{"metadata":{"path":str(HERE)}})
nbf.write(nb,NB); print("WROTE + EXECUTED", NB)
