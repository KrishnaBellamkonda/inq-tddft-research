#!/usr/bin/env python3
"""Two focused run-notebooks centred on the raw INQ total-energy curve — Poisson pilot
vs direct-potential pilot. Fast (no density GIF). Run with venv python."""
import nbformat as nbf, glob
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path

BASE = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/classical_highdensity_sv"
RUNS = [
  dict(label="POISSON perturbation (charge → Poisson → potential; old)",
       tag="pilot",
       run=BASE+"/pilot",
       out="/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/classical_highdensity_sv/pilot/pilot_energy_notebook.ipynb"),
  dict(label="DIRECT potential (erf/r added directly; no charge/Poisson/background; new)",
       tag="pilot_direct",
       run=BASE+"/pilot_direct",
       out="/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/classical_highdensity_sv/pilot_direct/pilot_direct_energy_notebook.ipynb"),
]
EGS = 207.18322156141

for R in RUNS:
    obs = sorted(glob.glob(R["run"]+"/results/**/observables.csv", recursive=True))[0]
    prj = sorted(glob.glob(R["run"]+"/results/**/projectile.csv", recursive=True))[0]
    cells = []
    cells.append(new_markdown_cell(
        f"# Total-energy curve — {R['label']}\n\n"
        f"Raw INQ `energy_total` (= `data.energy().total()`) for the v=2 pilot through the "
        f"r_s=4.18 slab. Mass-1 Gaussian projectile, σ_WP=0.5, z-open periodicity(2), no CAP.\n\n"
        f"Run: `{R['run']}`"))
    cells.append(new_code_cell(
        "import numpy as np, csv, matplotlib.pyplot as plt\n"
        f"OBS={obs!r}; PRJ={prj!r}; EGS={EGS}\n"
        "HA=27.211386; HALF=12.5; FARFACE=42.5\n"
        "o=list(csv.DictReader(open(OBS))); p=list(csv.DictReader(open(PRJ)))\n"
        "C=lambda r,k:np.array([float(x[k]) for x in r])\n"
        "t=C(o,'time_au'); Et=C(o,'energy_total')\n"
        "Ek=C(o,'energy_kinetic'); Eh=C(o,'energy_hartree'); Ex=C(o,'energy_xc'); Eext=C(o,'energy_external')\n"
        "pz=C(p,'proj_z'); pv=C(p,'proj_vz'); ke=C(p,'energy_proj_ke')*HA\n"
        "# map projectile z to the observable time grid\n"
        "pzt=np.interp(t,C(p,'time_au'),pz)\n"
        "t_enter=t[np.argmax(pzt>-HALF)]; t_slabx=t[np.argmax(pzt>HALF)]; t_exit=t[np.argmax(pzt>FARFACE)]\n"
        "print('proj_z: -30 -> %.1f ; vz 2.00 -> %.3f'%(pz[-1],pv[-1]))"))

    cells.append(new_markdown_cell("## Raw `energy_total(t)` — the curve you asked for (Ha, as INQ outputs)"))
    cells.append(new_code_cell(
        "fig,ax=plt.subplots(figsize=(9,4.5))\n"
        "ax.plot(t,Et,lw=1.5,c='C0')\n"
        "ax.axhline(EGS,ls='--',c='C7',lw=1,label='E_GS (projectile-absent) = %.3f Ha'%EGS)\n"
        "for tv,lab,c in [(t_enter,'enters slab','C2'),(t_slabx,'exits slab','C1'),(t_exit,'exits box (+42.5)','k')]:\n"
        "    ax.axvline(tv,ls=':',c=c,lw=1,label=lab)\n"
        "ax.set_xlabel('time (a.u.)'); ax.set_ylabel('energy_total (Ha)')\n"
        "plateau=Et[int(0.9*len(Et)):].mean()\n"
        "ax.set_title('energy_total: E(0)=%.3f -> plateau=%.3f Ha  (Δ vs E_GS = %.2f eV)'%(Et[0],plateau,(plateau-EGS)*HA))\n"
        "ax.legend(fontsize=8); fig.tight_layout(); plt.show()\n"
        "print('E_total(0)=%.4f Ha ; plateau=%.4f Ha ; plateau-E_GS=%.2f eV ; -dKE_proj(total)=%.2f eV'%(Et[0],plateau,(plateau-EGS)*HA,ke[0]-ke[-1]))"))

    cells.append(new_markdown_cell("## Per-step ΔE_total — is there an abrupt change at box exit?"))
    cells.append(new_code_cell(
        "dEt=np.diff(Et)*HA\n"
        "fig,ax=plt.subplots(1,2,figsize=(13,4))\n"
        "ax[0].plot(t[1:],dEt,lw=.8); ax[0].axvline(t_exit,ls=':',c='k',label='box exit')\n"
        "ax[0].set_xlabel('t'); ax[0].set_ylabel('ΔE_total per step (eV)'); ax[0].set_title('per-step jump (full run)'); ax[0].legend(fontsize=8)\n"
        "m=(pzt[1:]>39)&(pzt[1:]<46)\n"
        "ax[1].plot(pzt[1:][m],dEt[m],'.-',ms=3); ax[1].axvline(FARFACE,ls=':',c='k')\n"
        "ax[1].set_xlabel('proj_z (Bohr)'); ax[1].set_ylabel('ΔE_total per step (eV)'); ax[1].set_title('zoom on the far-face crossing')\n"
        "fig.tight_layout(); plt.show()\n"
        "print('max |ΔE_total/step| in the exit window (proj_z 39-46) = %.2f eV/step'%(np.max(np.abs(dEt[m])) if m.any() else 0))"))

    cells.append(new_markdown_cell("## Energy components (Ha) — where the total moves"))
    cells.append(new_code_cell(
        "fig,ax=plt.subplots(figsize=(9,4.5))\n"
        "for arr,lab in [(Ek,'kinetic'),(Eh,'hartree'),(Ex,'xc'),(Eext,'external'),(Et,'TOTAL')]:\n"
        "    ax.plot(t,arr-arr[0],label=lab)\n"
        "ax.axvline(t_exit,ls=':',c='k'); ax.set_xlabel('t'); ax.set_ylabel('component − its t=0 value (Ha)')\n"
        "ax.set_title('energy components (Δ from t=0)'); ax.legend(fontsize=8); fig.tight_layout(); plt.show()"))

    nb=new_notebook(cells=cells,metadata={"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"}})
    ep=ExecutePreprocessor(timeout=300,kernel_name="python3"); ep.preprocess(nb,{"metadata":{"path":str(Path(R['out']).parent)}})
    nbf.write(nb,R['out']); print("WROTE",R['out'])
