#!/usr/bin/env python3
"""Build + execute the per-velocity RUN notebooks and the PHASE notebook for the
high-density classical S(v) sweep, plus a tabulation. venv python.

Run notebook  : density GIF on top + trajectory + step-by-step S calc (deposit +
                KE-loss cross-check) + full pairwise ledger + conservation check.
Phase notebook: results table + S(v) curve (both methods) + Bethe-tail power-law
                fit + component-ledger deltas (Definition-1 staging) + caveats.
"""
import nbformat as nbf, json, os, math
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path

SYS = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
HYP = f"{SYS}/hypotheses/classical_highdensity_sv/sv_sweep"
RUNS = f"{SYS}/scripts/classical_highdensity_sv/dyn/results"
E_GS = 207.18322156141
HA = 27.2114; LSLAB = 25.0; FAR = 42.5
VELS = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5]

def vname(v): return f"v{v:.1f}".replace(".", "p")

# ---------- per-velocity RUN notebook ----------
def run_notebook(v):
    n = vname(v); rdir = f"{RUNS}/{n}"; adir = f"{HYP}/{n}"
    gif = f"{adir}/density_evolution.gif"
    cells = []
    cells.append(new_markdown_cell(
        f"# Run notebook — classical projectile v={v} (r_s=4.18 slab)\n\n"
        f"Campaign `classical-highdensity-sv`. Mass-1 Gaussian-charge electron "
        f"(σ_WP=0.5), Ehrenfest, launched at z=−24 through the 25-Bohr r_s=4.18 "
        f"slab, z-open `periodicity(2)`, **no CAP**. This notebook re-derives the "
        f"stopping power from the raw run and shows the full energy ledger."))
    # density GIF at TOP (rule: notebook-density-gif)
    cells.append(new_markdown_cell("## Density evolution n(x,z,t) — visual intuition (top of notebook)"))
    cells.append(new_code_cell(
        "from IPython.display import Image, display\n"
        f"import os\n"
        f"g = {gif!r}\n"
        "display(Image(filename=g)) if os.path.exists(g) else print('no density GIF')"))
    cells.append(new_code_cell(
        "import numpy as np, pandas as pd, matplotlib.pyplot as plt\n"
        f"R = {rdir!r}; E_GS={E_GS}; HA={HA}; LSLAB={LSLAB}; FAR={FAR}; V={v}\n"
        "obs=pd.read_csv(R+'/raw/observables/observables.csv')\n"
        "proj=pd.read_csv(R+'/raw/observables/projectile.csv')\n"
        "ix=pd.read_csv(R+'/raw/observables/interactions.csv')\n"
        "t=obs['time_au'].values; E=obs['energy_total'].values\n"
        "z=proj['proj_z'].values; vz=proj['proj_vz'].values; ke=proj['energy_proj_ke'].values\n"
        "ubg=proj['energy_proj_bg_ideal'].values\n"
        "print('frames',len(t),'  z:',z[0],'->',round(z[-1],1),'  v:',vz[0],'->',round(vz[-1],3))"))
    cells.append(new_markdown_cell(
        "## Trajectory + deposit\n"
        "The projectile transits the slab (shaded) and exits past the far face "
        "(dashed); once it is fully out, E_electronic is flat = the plateau."))
    cells.append(new_code_cell(
        "fig,ax=plt.subplots(1,2,figsize=(12,4))\n"
        "ax[0].plot(t,z); ax[0].axhspan(-12.5,12.5,color='C7',alpha=.15); ax[0].axhline(FAR,ls='--',c='k')\n"
        "ax[0].set_xlabel('t (au)'); ax[0].set_ylabel('proj_z (Bohr)'); ax[0].set_title('trajectory')\n"
        "ax[1].plot(t,(E-E_GS)*HA); ax[1].axhline(0,ls=':',c='k')\n"
        "ax[1].set_xlabel('t (au)'); ax[1].set_ylabel('E_elec − E_GS (eV)'); ax[1].set_title('slab excitation')\n"
        "fig.tight_layout(); plt.show()"))
    cells.append(new_markdown_cell(
        "## Stopping power — step by step\n"
        "**Definition 2 (headline):** S = E_absorbed / L_slab, with E_absorbed the "
        "slab excitation once the projectile has left (referenced to the "
        "projectile-free GS, E_GS). **Cross-check:** the projectile's KE loss "
        "(energy conservation, CAP-free) must give the same S."))
    cells.append(new_code_cell(
        "# 1) exit frame: projectile centre fully past the far face (+2 Bohr tail)\n"
        "i0=int(np.argmax(z>FAR+2.0))\n"
        "print('1) exit at frame',i0,' t=%.2f au, z=%.1f Bohr'%(t[i0],z[i0]))\n"
        "# 2) plateau = mean E_total after exit\n"
        "plateau=E[i0:].mean(); flat=E[i0:].std()*HA\n"
        "print('2) plateau E_total = %.5f Ha  (flatness %.2e eV)'%(plateau,flat))\n"
        "# 3) E_absorbed = plateau - E_GS (projectile-free reference)\n"
        "Eabs=(plateau-E_GS)*HA; print('3) E_absorbed = (%.5f - %.5f)*27.2114 = %.2f eV'%(plateau,E_GS,Eabs))\n"
        "# 4) S = E_absorbed / L_slab\n"
        "S=Eabs/LSLAB; print('4) S = %.2f / %.0f = %.3f eV/Bohr'%(Eabs,LSLAB,S))\n"
        "# 5) cross-check: projectile KE loss to exit\n"
        "keloss=(ke[0]-ke[i0])*HA; Sk=keloss/LSLAB\n"
        "print('5) KE loss = %.2f eV -> S_keloss = %.3f eV/Bohr'%(keloss,Sk))\n"
        "print('   agreement: S/S_keloss = %.4f'%(S/Sk))\n"
        "# effective (mean) velocity over the slab\n"
        "inslab=np.abs(z)<12.5; vmean=vz[inslab].mean() if inslab.any() else np.nan\n"
        "print('   mean v in slab = %.3f (launch %.1f -> exit %.3f)'%(vmean,vz[0],vz[i0]))"))
    cells.append(new_markdown_cell(
        "## Full pairwise Coulomb ledger (Definition-1 staging)\n"
        "Every pairwise term vs time — the raw material for the (still-TBD) "
        "energy-decomposition stopping definition. P=projectile, S=slab electrons, "
        "B=background."))
    cells.append(new_code_cell(
        "tt=ix['time_au'].values\n"
        "fig,ax=plt.subplots(figsize=(9,4.5))\n"
        "for c,lab in [('e_pp','E_PP proj self'),('e_ps','E_PS proj–slab'),('e_ss','E_SS slab self'),\n"
        "              ('e_sb','E_SB slab–bg'),('e_pb','E_PB proj–bg'),('e_bb','E_BB bg self')]:\n"
        "    ax.plot(tt,(ix[c]-ix[c].iloc[0])*HA,label=lab)\n"
        "ax.set_xlabel('t (au)'); ax.set_ylabel('Δ term (eV)'); ax.legend(fontsize=8); ax.set_title('pairwise ledger (Δ vs t=0)')\n"
        "fig.tight_layout(); plt.show()"))
    cells.append(new_markdown_cell(
        "## Conservation check (correctness)\n"
        "E_electronic + KE_proj + U_proj_bg must stay flat (CAP-free ⇒ no energy "
        "sink). Its constancy is what makes E_absorbed trustworthy."))
    cells.append(new_code_cell(
        "cons=(E+ke+ubg)  # Ha ; note U_proj_bg sign convention per run.cpp\n"
        "fig,ax=plt.subplots(figsize=(8,3.5))\n"
        "ax.plot(t,(cons-cons[0])*HA); ax.set_xlabel('t (au)'); ax.set_ylabel('Δ(E_elec+KE+U_bg) (eV)')\n"
        "ax.set_title('conservation drift'); fig.tight_layout(); plt.show()\n"
        "print('max |drift| = %.3f eV'%(np.abs((cons-cons[0])*HA).max()))"))
    nb = new_notebook(cells=cells, metadata={"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"}})
    ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
    ep.preprocess(nb, {"metadata": {"path": adir}})
    nbf.write(nb, f"{adir}/run_{n}.ipynb")
    print("wrote", f"{adir}/run_{n}.ipynb")

# ---------- PHASE notebook ----------
def phase_notebook():
    cells = []
    cells.append(new_markdown_cell(
        "# Phase notebook — high-density classical S(v) benchmark\n\n"
        "Campaign `classical-highdensity-sv`. Mass-1 Gaussian-charge electron "
        "(σ_WP=0.5) Ehrenfest through a 25-Bohr r_s=4.18 jellium slab, z-open "
        "`periodicity(2)`, CAP-free. Six launch velocities. Headline S = "
        "E_absorbed/L_slab (slab excitation vs the projectile-free GS); "
        "cross-checked by the projectile's KE loss (energy conservation)."))
    cells.append(new_code_cell(
        "import numpy as np, pandas as pd, matplotlib.pyplot as plt, json\n"
        f"HYP={HYP!r}; RUNS={RUNS!r}; E_GS={E_GS}; HA={HA}; LSLAB={LSLAB}; FAR={FAR}\n"
        f"VELS={VELS}\n"
        "def vname(v): return ('v%.1f'%v).replace('.','p')\n"
        "rows=[]\n"
        "for v in VELS:\n"
        "    r=json.load(open(f'{HYP}/{vname(v)}/result.json'))\n"
        "    proj=pd.read_csv(f'{RUNS}/{vname(v)}/raw/observables/projectile.csv')\n"
        "    z=proj['proj_z'].values; vz=proj['proj_vz'].values\n"
        "    inslab=np.abs(z)<12.5; vmean=float(vz[inslab].mean()) if inslab.any() else np.nan\n"
        "    rows.append(dict(v=v, v_final=r['v_final'], v_mean=vmean, S=r['S'],\n"
        "                     S_keloss=r.get('S_keloss'), E_abs_eV=r['E_absorbed_eV'],\n"
        "                     flat_eV=r.get('plateau_flatness_eV')))\n"
        "df=pd.DataFrame(rows); df.round(3)"))
    cells.append(new_markdown_cell("## Results table"))
    cells.append(new_code_cell("df.round({'v':1,'v_final':2,'v_mean':2,'S':3,'S_keloss':3,'E_abs_eV':2,'flat_eV':5})"))
    cells.append(new_markdown_cell(
        "## S(v) curve — two independent methods\n"
        "The two channels (slab deposit vs projectile KE loss) agree to ~0.1%, "
        "confirming energy conservation and a clean plateau. Plotted vs launch v "
        "and vs the mean in-slab v (the electron decelerates, so the effective "
        "velocity is below the launch value)."))
    cells.append(new_code_cell(
        "fig,ax=plt.subplots(1,2,figsize=(12,4.5))\n"
        "ax[0].plot(df.v,df.S,'o-',label='S (deposit/L)'); ax[0].plot(df.v,df.S_keloss,'x--',label='S (KE-loss/L)')\n"
        "ax[0].set_xlabel('launch v (au)'); ax[0].set_ylabel('S (eV/Bohr)'); ax[0].set_title('S vs launch v'); ax[0].legend()\n"
        "ax[1].plot(df.v_mean,df.S,'s-',c='C2'); ax[1].set_xlabel('mean in-slab v (au)'); ax[1].set_ylabel('S (eV/Bohr)')\n"
        "ax[1].set_title('S vs effective (mean) v'); fig.tight_layout(); plt.show()"))
    cells.append(new_markdown_cell(
        "## Bethe-tail analysis\n"
        "Above the Lindhard peak (v_F≈0.46; peak ~v≈0.5–1), stopping falls with v. "
        "Fit S ∝ v^(−n) to characterise the tail."))
    cells.append(new_code_cell(
        "p=np.polyfit(np.log(df.v.values), np.log(df.S.values), 1)\n"
        "print('power-law fit: S ∝ v^%.2f  (n=%.2f)'%(p[0],-p[0]))\n"
        "fig,ax=plt.subplots(figsize=(6,4))\n"
        "ax.loglog(df.v,df.S,'o-'); ax.loglog(df.v, np.exp(p[1])*df.v**p[0],'k:',label=f'v^{p[0]:.2f}')\n"
        "ax.set_xlabel('v (au)'); ax.set_ylabel('S (eV/Bohr)'); ax.legend(); ax.set_title('log–log tail')\n"
        "fig.tight_layout(); plt.show()"))
    cells.append(new_markdown_cell(
        "## Component-ledger deltas across the transit (Definition-1 staging)\n"
        "Net change of each pairwise Coulomb term from t=0 to the plateau, per v — "
        "the raw material for the energy-decomposition stopping definition."))
    cells.append(new_code_cell(
        "led=[]\n"
        "for v in VELS:\n"
        "    ix=pd.read_csv(f'{RUNS}/{vname(v)}/raw/observables/interactions.csv')\n"
        "    d={c:(ix[c].iloc[-1]-ix[c].iloc[0])*HA for c in ['e_pp','e_ps','e_ss','e_sb','e_pb']}\n"
        "    d['v']=v; led.append(d)\n"
        "ldf=pd.DataFrame(led)[['v','e_pp','e_ps','e_ss','e_sb','e_pb']]\n"
        "print('Δ pairwise terms (eV), t=0 -> final:'); display(ldf.round(2))"))
    cells.append(new_markdown_cell(
        "## Caveats\n"
        "- **Deceleration:** each S is an average over v∈[v_final, v_launch]; use "
        "the mean-v column for theory overlays.\n"
        "- **High-v tail only:** a transiting mass-1 electron cannot reach the "
        "Lindhard peak (it would stop inside), so this curve is the Bethe tail.\n"
        "- **dx=0.5:** the strict 2v-response cutoff bound is exceeded for v>3.1 — "
        "the fast points (v≥3.5) carry a mild resolution caveat (finer-dx GS would "
        "confirm).\n"
        "- **Definition-1 formula still TBD** — the ledger is collected, not yet "
        "reduced to a closed S."))
    nb = new_notebook(cells=cells, metadata={"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"}})
    ep = ExecutePreprocessor(timeout=900, kernel_name="python3")
    ep.preprocess(nb, {"metadata": {"path": HYP}})
    nbf.write(nb, f"{HYP}/phase_notebook.ipynb")
    print("wrote", f"{HYP}/phase_notebook.ipynb")

if __name__ == "__main__":
    for v in VELS:
        run_notebook(v)
    phase_notebook()
    print("DONE")
