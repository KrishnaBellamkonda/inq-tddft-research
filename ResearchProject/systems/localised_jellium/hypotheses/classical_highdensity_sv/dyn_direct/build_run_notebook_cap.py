#!/usr/bin/env python3
"""
Generalized run notebook for a DIRECT-potential classical projectile run (any velocity).
Usage: build_run_notebook.py <vtag>   e.g.  v3p0
Produces run_<vtag>_direct.ipynb beside this file. Compact deep-dive: density GIFs
(total + induced, linear|log), energetics (dE_total + N(t)), no-kink OLD-vs-NEW overlay,
clean pairwise ledger, projectile transport + stopping, takeaway. Auto-gates on data.
"""
import os, glob, sys, re
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

TAG = sys.argv[1] if len(sys.argv) > 1 else "v4p5"   # v3p0 (sigma_WP=0.5 family) OR s2p0_v2p0 (sigma sweep)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/local/data/public/skcb2/tddft"
SYS  = f"{ROOT}/ResearchProject/systems/localised_jellium"
NEW  = f"{SYS}/scripts/classical_highdensity_sv/dyn_direct_cap/results/{TAG}_cap"
OLD  = f"{SYS}/scripts/classical_highdensity_sv/dyn_direct/results/{TAG}_direct"  # overlay: no-CAP direct twin (sigma_WP=0.5 family only)
E_GS = 207.18322156141
NB   = f"{HERE}/run_{TAG}_cap.ipynb"
sys.path.insert(0, f"{ROOT}/inq-stack/python")
# --- reconstruct config from run_summary.txt so one builder serves both run families ---
def _rs(key, default, cast=float):
    m = re.search(rf"\b{key}\s*=\s*([-\d.eE+]+)", open(NEW + "/run_summary.txt").read())
    return cast(m.group(1)) if m else default
V         = _rs("v0", 0.0)
SIGMA_WP  = _rs("sigma_wp", 0.5)
SIGMA_POT = SIGMA_WP / 2**0.5
FACE      = _rs("slab_half", 12.5)
WALL      = _rs("Lz", 85.0) / 2.0
DT        = _rs("dt", 0.04)
VTAG = TAG   # back-compat alias for gif filenames / titles below

def build_gifs():
    import numpy as np, matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, matplotlib.animation as anim
    from matplotlib.colors import LogNorm, SymLogNorm
    from inqview import load_vti; import pandas as pd
    frames = sorted(glob.glob(f"{NEW}/frames/total/*.vti"), key=lambda f: int(f.split('_t')[-1].split('.')[0]))
    if not frames: print("no frames"); return {}
    pj = pd.read_csv(f"{NEW}/raw/observables/projectile.csv"); zmap = dict(zip(pj.step, pj.proj_z))
    sel = frames[::max(1, len(frames)//48)]
    d0 = load_vti(frames[0], expect_centered_axis="z"); iy = np.argmin(abs(d0.y))
    x, z = d0.x, d0.z; ext = [x[0], x[-1], z[0], z[-1]]; n0 = d0.data[:, iy, :].T
    vmax = float(np.percentile(load_vti(sel[len(sel)//3], expect_centered_axis="z").data[:, iy, :], 99.9))
    ind = 0.0
    for f in sel:
        ind = max(ind, float(np.percentile(np.abs(load_vti(f, expect_centered_axis="z").data[:, iy, :].T - n0), 99.5)))
    ind = ind or 1e-6; out = {}
    def faces(a):
        for zf in (-FACE, FACE): a.axhline(zf, ls="--", c="w", lw=0.8, alpha=0.7)
        a.axhline(WALL, ls=":", c="r", lw=1.0); a.set_xlabel("x"); a.set_ylabel("z")
    def mk(kind, fn, ttl):
        fig, ax = plt.subplots(1, 2, figsize=(9, 5), constrained_layout=True); ref = {"d": n0.copy()}
        def upd(f):
            n = load_vti(f, expect_centered_axis="z").data[:, iy, :].T
            step = int(f.split('_t')[-1].split('.')[0]); zc = zmap.get(step, np.nan)
            for a in ax: a.clear()
            if kind == "total":
                ax[0].imshow(n, origin="lower", extent=ext, aspect="auto", cmap="viridis", vmin=0, vmax=vmax); ax[0].set_title("n linear")
                ax[1].imshow(np.maximum(n, 1e-8), origin="lower", extent=ext, aspect="auto", cmap="viridis", norm=LogNorm(vmax*1e-4, vmax)); ax[1].set_title("n log")
            else:
                d = n - n0
                ax[0].imshow(d, origin="lower", extent=ext, aspect="auto", cmap="RdBu_r", vmin=-ind, vmax=ind); ax[0].set_title("Δn linear")
                ax[1].imshow(d, origin="lower", extent=ext, aspect="auto", cmap="RdBu_r", norm=SymLogNorm(ind/100, vmin=-ind, vmax=ind)); ax[1].set_title("Δn symlog")
            for a in ax:
                faces(a)
                if np.isfinite(zc): a.axhline(zc, ls="-", c="k", lw=0.6, alpha=0.5)
            fig.suptitle(f"{VTAG} direct {ttl} — step {step}, proj_z={zc:+.1f}")
        anim.FuncAnimation(fig, upd, frames=sel, interval=120).save(f"{HERE}/{fn}", writer=anim.PillowWriter(fps=8)); plt.close(fig)
    mk("total", f"gif_total_{VTAG}_cap.gif", "total density"); out["total"] = f"{HERE}/gif_total_{VTAG}_cap.gif"
    mk("induced", f"gif_induced_{VTAG}_cap.gif", "induced Δn"); out["induced"] = f"{HERE}/gif_induced_{VTAG}_cap.gif"
    return out
GIFS = build_gifs()

cells = []
def md(s): cells.append(new_markdown_cell(s))
def code(s): cells.append(new_code_cell(s))

md(f"# Classical projectile in jellium slab — {TAG} (v={V:g}, σ_WP={SIGMA_WP:g}), DIRECT potential + two-sided CAP\n"
   f"Direct erf/r potential (charge-free) classical projectile. System: r_s≈4.18, N=100, box 35×35×85, "
   f"slab ±{FACE:g}, `periodicity(2)`. σ_WP={SIGMA_WP:g} (σ_pot=σ_WP/√2={SIGMA_POT:.3f}), launch z=−24, "
   "mass 1, free Ehrenfest, two-sided sin² CAP (η=−1 Ha, inner ±30, width 12.5). "
   "Shows: density evolution, physical ledger, projectile transport, stopping (Def-A KE-loss & Def-B E_absorbed).")
md("## Setup (run_summary.txt)")
code(f"import pandas as pd,numpy as np,matplotlib.pyplot as plt,os\nfrom matplotlib.colors import LogNorm\n"
     f"NEW=r'{NEW}';OLD=r'{OLD}';HA=27.211386;E_GS={E_GS!r};FACE={FACE};WALL={WALL};DT={DT};V={V}\n"
     "print(open(NEW+'/run_summary.txt').read())\n"
     "import glob as _g\n"
     "def _cat(b,stem):\n"
     "    fs=sorted(_g.glob(b+stem+'*.csv'))\n"
     "    return pd.concat([pd.read_csv(f) for f in fs]).drop_duplicates('step').sort_values('step').reset_index(drop=True)\n"
     "def load(r):\n"
     "    b=r+'/raw/observables/'\n"
     "    return _cat(b,'observables').merge(_cat(b,'projectile'),on=['step','time_au']).merge(_cat(b,'interactions'),on=['step','time_au'])\n"
     "N=load(NEW); O=load(OLD) if os.path.exists(OLD+'/raw/observables/observables.csv') else None\n"
     "print('NEW steps: %d -> %d (segments concatenated)'%(N.step.min(),N.step.max()))")
md("## Density evolution (mid-y x–z, LINEAR | LOG) — total then induced Δn")
for k, cap in [("total", "Total n(x,z,t)"), ("induced", "Induced Δn=n(t)−n(0)")]:
    g = GIFS.get(k)
    if g:
        md(f"**{cap}**"); code(f"from IPython.display import Image,display\ndisplay(Image(filename=r'{g}'))")
md("## Energetics — ΔE_total(t) with N(t) beside; conservation")
code("from inqview import load_vti\nimport glob\n"
     "fr=sorted(glob.glob(NEW+'/frames/total/*.vti'),key=lambda f:int(f.split('_t')[-1].split('.')[0]))\n"
     "sel=fr[::max(1,len(fr)//120)];Nt=[];tt=[]\n"
     "for f in sel:\n"
     "    d=load_vti(f,expect_centered_axis='z');dv=(d.x[1]-d.x[0])*(d.y[1]-d.y[0])*(d.z[1]-d.z[0]);Nt.append(d.data.sum()*dv);tt.append(int(f.split('_t')[-1].split('.')[0])*DT)\n"
     "fig,ax=plt.subplots(1,2,figsize=(13,4))\n"
     "ax[0].plot(N.time_au,(N.energy_total-E_GS)*HA);ax[0].set_title('ΔE_total=E-E_GS (eV)');ax[0].set_xlabel('t')\n"
     "ax[1].plot(tt,Nt);ax[1].set_title('N(t)=∫n dV (conserved)');ax[1].set_xlabel('t');plt.tight_layout();plt.show()\n"
     "cons=N.energy_total+N.energy_proj_ke+N.energy_proj_bg_ideal\n"
     "print('conservation drift %.2e eV; N(t) %.4f->%.4f'%((cons.max()-cons.min())*HA,Nt[0],Nt[-1]))")
md("## No-kink demonstration — no-CAP direct vs NEW direct (curvature max should be IN-SLAB)")
code("if O is not None:\n"
     "    fig,ax=plt.subplots(2,3,figsize=(15,8))\n"
     "    def pan(a,c,t,ref=0,vz=False):\n"
     "        a.plot(O.proj_z,O[c] if vz else (O[c]-ref)*HA,c='tab:red',alpha=0.7,label='no-CAP direct');a.plot(N.proj_z,N[c] if vz else (N[c]-ref)*HA,c='tab:blue',label='NEW direct')\n"
     "        a.axvline(WALL,ls=':',c='k');a.axvline(FACE,ls='--',c='grey');a.set_title(t);a.set_xlabel('proj_z');a.legend(fontsize=8)\n"
     "    pan(ax[0,0],'energy_total','E-E_GS (eV)',ref=E_GS);pan(ax[0,1],'e_ps','E_PS (eV)');pan(ax[0,2],'e_pb','E_PB (eV)')\n"
     "    pan(ax[1,0],'e_pp','E_PP (eV)');pan(ax[1,1],'energy_proj_bg_ideal','U_proj_bg (eV)');pan(ax[1,2],'proj_vz','vz',vz=True)\n"
     "    plt.tight_layout();plt.show()\n"
     "    mc=lambda d,c:d.proj_z.values[np.argmax(np.abs(np.gradient(np.gradient(d[c].values))))]\n"
     "    print('%-20s %8s %8s'%('series','OLD z*','NEW z*'))\n"
     "    for c in ['energy_total','e_ps','e_pb','e_pp']:print('%-20s %8.1f %8.1f'%(c,mc(O,c),mc(N,c)))\n"
     "else: print('no OLD run to overlay')")
md("## Direct pairwise ledger — E_PS>0→0, E_PP const, E_PB<0→0")
code("fig,ax=plt.subplots(figsize=(8,4))\n"
     "for c,l in [('e_ss','E_SS'),('e_pp','E_PP'),('e_ps','E_PS'),('e_sb','E_SB'),('e_pb','E_PB')]:ax.plot(N.proj_z,N[c],label=l)\n"
     "ax.axvline(FACE,ls='--',c='grey');ax.axvline(WALL,ls=':',c='r');ax.axhline(0,c='k',lw=.5);ax.legend(fontsize=8);ax.set_xlabel('proj_z');ax.set_ylabel('Ha');plt.show()\n"
     "print('E_PP mean %.4f std %.1e'%(N.e_pp.mean(),N.e_pp.std()))")
md("## Projectile transport & stopping — TWO definitions\n"
   "**Def-A** (KE-loss): S_A = [KE(−12.5)−KE(+12.5)]/L_slab — net projectile KE lost across the "
   "equal-potential slab window, per 25 Bohr (gauge-free).  \n"
   "**Def-B** (E_absorbed): S_B = [E_total(t_final)−E_GS]/L_slab — total electronic energy gain at "
   "the final timestep, per 25 Bohr. NOTE: with the direct 1/r potential E_total does not fully "
   "plateau (residual e_ps offset), so S_B is an upper bound that drifts down as the run is extended.")
code("fig,ax=plt.subplots(1,3,figsize=(16,4))\n"
     "ax[0].plot(N.time_au,N.proj_z);ax[0].axhline(FACE,ls='--',c='grey');ax[0].axhline(-FACE,ls='--',c='grey');ax[0].set_title('z(t)')\n"
     "ax[1].plot(N.proj_z,N.proj_vz);ax[1].axvline(FACE,ls='--',c='grey');ax[1].axvline(-FACE,ls='--',c='grey');ax[1].set_title('v(z)')\n"
     "ax[2].plot(N.proj_z,N.energy_proj_ke*HA);ax[2].axvline(FACE,ls='--',c='grey');ax[2].axvline(-FACE,ls='--',c='grey');ax[2].set_title('KE(z) eV');plt.tight_layout();plt.show()\n"
     "ke=lambda d,z:d.loc[(d.proj_z-z).abs().idxmin(),'energy_proj_ke']\n"
     "S_A=(ke(N,-FACE)-ke(N,FACE))*HA/(2*FACE)\n"
     "Eabs=(N.energy_total.iloc[-1]-E_GS)*HA; S_B=Eabs/(2*FACE)\n"
     "print('v_launch=%.3f v_final=%.3f  (t_final=%.1f a.u., z_final=%.1f Bohr)'%(N.proj_vz.iloc[0],N.proj_vz.iloc[-1],N.time_au.iloc[-1],N.proj_z.iloc[-1]))\n"
     "print('Def-A  S_A = KEloss(-12.5..+12.5)/25       = %.3f eV/Bohr'%S_A)\n"
     "print('Def-B  S_B = [E_total(t_final)-E_GS]/25    = %.3f eV/Bohr   (E_absorbed=%.2f eV)'%(S_B,Eabs))\n"
     "if O is not None: print('  (no-CAP direct Def-A = %.3f eV/Bohr -> no-CAP)'%((ke(O,-FACE)-ke(O,FACE))*HA/(2*FACE)))")
md(f"## Takeaway\n- Direct erf/r potential: no exit kink, physical ledger (E_PS>0→0, E_PP const, E_PB<0→0), energy conserved.\n"
   f"- Two stopping definitions reported: **Def-A** (KE-loss across slab /25, gauge-free) and **Def-B** "
   f"([E_total(t_final)−E_GS]/25). Def-B > Def-A because E_total retains the 1/r e_ps offset (does not fully plateau).\n"
   "- Extended run (+25% steps) concatenated from segments; corrected classical baseline at r_s=4.18.")

nb = new_notebook(cells=cells)
nb.metadata.kernelspec = {"name": "python3", "display_name": "Python 3", "language": "python"}
ExecutePreprocessor(timeout=2400, kernel_name="python3").preprocess(nb, {"metadata": {"path": HERE}})
with open(NB, "w") as f: nbf.write(nb, f)
nerr = sum(1 for c in nb.cells if c.cell_type == "code" for o in c.get("outputs", []) if o.get("output_type") == "error")
print(f"wrote {NB} | errors {nerr}")
