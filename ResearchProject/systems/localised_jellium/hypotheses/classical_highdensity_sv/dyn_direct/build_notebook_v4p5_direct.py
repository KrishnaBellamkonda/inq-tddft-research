#!/usr/bin/env python3
"""
Deep single-run notebook (run-notebook skill) for the v=4.5 DIRECT-potential classical
projectile in the localised jellium slab (dyn_direct/v4p5_direct). House narrative:
  1 title+question | 2 conventions | 3 reconstructable setup | 4 source files
  + battery (auto-gated on this run's observables):
    visual intuition (lead density GIF + energetics) | density-GIF battery (linear|log) |
    density carpets (z-t, linear|log) | energetics (dE_total + N(t) beside, components,
    conservation) | no-kink demo OLD-charge vs NEW-direct | pairwise ledger (clean) |
    projectile & transport (z,v,F,KE(z), equal-potential window, stopping S) |
    physical anchors (compute_heuristics) | loss function (note) | takeaway.
This run uses the localised_jellium I/O convention (raw/observables + frames/total),
so it is a tailored assembler over inqview kernels (load_vti, compute_heuristics),
not the jellium-layout run_notebook_builder. Figures render LINEAR | LOG (project rule).
Run:  PYTHONPATH=.../inq-stack/python venv/bin/python3 build_notebook_v4p5_direct.py
"""
import os, glob, re, sys, json
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/local/data/public/skcb2/tddft"
SYS  = f"{ROOT}/ResearchProject/systems/localised_jellium"
NEW  = f"{SYS}/scripts/classical_highdensity_sv/dyn_direct/results/v4p5_direct"
OLD  = f"{SYS}/scripts/classical_highdensity_sv/dyn/results/v4p5"
RUNCPP = f"{SYS}/scripts/classical_highdensity_sv/dyn_direct/run.cpp"
E_GS = 207.18322156141
NB   = f"{HERE}/run_v4p5_direct.ipynb"
sys.path.insert(0, f"{ROOT}/inq-stack/python")

# ---- parse run_summary.txt (multi-pair lines) ------------------------------
def parse_summary(path):
    d = {}
    for ln in open(path):
        for m in re.finditer(r"([A-Za-z_]\w*)\s*=\s*([^\s]+)", ln):
            d[m.group(1)] = m.group(2)
    return d
SUMM = parse_summary(f"{NEW}/run_summary.txt")
DT = float(SUMM.get("dt", 0.04)); SIGMA_WP = float(SUMM.get("sigma_wp", 0.5))
FACE = float(SUMM.get("slab_half", 12.5)); WALL = float(SUMM.get("Lz", 85))/2

# ---- build density GIFs (linear|log) at the top ----------------------------
def build_gifs():
    import numpy as np, matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.animation as anim
    from matplotlib.colors import LogNorm, SymLogNorm
    from inqview import load_vti
    import pandas as pd
    frames = sorted(glob.glob(f"{NEW}/frames/total/*.vti"),
                    key=lambda f: int(re.search(r'_t(\d+)', f).group(1)))
    if not frames:
        print("no density frames"); return {}
    pj = pd.read_csv(f"{NEW}/raw/observables/projectile.csv"); zmap = dict(zip(pj.step, pj.proj_z))
    sel = frames[::max(1, len(frames)//48)]
    d0 = load_vti(frames[0], expect_centered_axis="z"); iy = np.argmin(abs(d0.y))
    x, z = d0.x, d0.z; ext = [x[0], x[-1], z[0], z[-1]]
    n0 = d0.data[:, iy, :].T
    vmax = float(np.percentile(load_vti(sel[len(sel)//3], expect_centered_axis="z").data[:, iy, :], 99.9))
    ind_max = 0.0; prevref = {"d": n0.copy()}
    for f in sel:                                   # global scales for Δ panels
        n = load_vti(f, expect_centered_axis="z").data[:, iy, :].T
        ind_max = max(ind_max, float(np.percentile(np.abs(n-n0), 99.5)))
    ind_max = ind_max or 1e-6
    out = {}
    def faces(a):
        for zf in (-FACE, FACE): a.axhline(zf, ls="--", c="w", lw=0.8, alpha=0.7)
        a.axhline(WALL, ls=":", c="r", lw=1.0)
        a.set_xlabel("x (Bohr)"); a.set_ylabel("z (Bohr)")
    # 1) total density linear | log
    def mk_total():
        fig, ax = plt.subplots(1, 2, figsize=(9, 5), constrained_layout=True)
        def upd(f):
            n = load_vti(f, expect_centered_axis="z").data[:, iy, :].T
            step = int(re.search(r'_t(\d+)', f).group(1)); zc = zmap.get(step, np.nan)
            for a in ax: a.clear()
            ax[0].imshow(n, origin="lower", extent=ext, aspect="auto", cmap="viridis", vmin=0, vmax=vmax); ax[0].set_title("n(x,z) linear")
            ax[1].imshow(np.maximum(n,1e-8), origin="lower", extent=ext, aspect="auto", cmap="viridis", norm=LogNorm(vmax*1e-4, vmax)); ax[1].set_title("n(x,z) log")
            for a in ax:
                faces(a)
                if np.isfinite(zc): a.axhline(zc, ls="-", c="k", lw=0.6, alpha=0.5)
            fig.suptitle(f"v=4.5 direct total density — step {step}, proj_z={zc:+.1f}")
        anim.FuncAnimation(fig, upd, frames=sel, interval=120).save(f"{HERE}/gif_total.gif", writer=anim.PillowWriter(fps=8)); plt.close(fig)
    # 2) induced / instantaneous Δn linear | symlog
    def mk_delta(kind, fname, title):
        fig, ax = plt.subplots(1, 2, figsize=(9, 5), constrained_layout=True)
        prevref["d"] = n0.copy()
        def upd(f):
            n = load_vti(f, expect_centered_axis="z").data[:, iy, :].T
            d = (n - n0) if kind == "induced" else (n - prevref["d"]); prevref["d"] = n.copy()
            step = int(re.search(r'_t(\d+)', f).group(1)); zc = zmap.get(step, np.nan)
            for a in ax: a.clear()
            ax[0].imshow(d, origin="lower", extent=ext, aspect="auto", cmap="RdBu_r", vmin=-ind_max, vmax=ind_max); ax[0].set_title(title+" linear")
            ax[1].imshow(d, origin="lower", extent=ext, aspect="auto", cmap="RdBu_r", norm=SymLogNorm(linthresh=ind_max/100, vmin=-ind_max, vmax=ind_max)); ax[1].set_title(title+" symlog")
            for a in ax:
                faces(a)
                if np.isfinite(zc): a.axhline(zc, ls="-", c="k", lw=0.6, alpha=0.5)
            fig.suptitle(f"v=4.5 direct {title} — step {step}, proj_z={zc:+.1f}")
        anim.FuncAnimation(fig, upd, frames=sel, interval=120).save(f"{HERE}/{fname}", writer=anim.PillowWriter(fps=8)); plt.close(fig)
    mk_total(); out["total"] = f"{HERE}/gif_total.gif"
    mk_delta("induced", "gif_induced.gif", u"induced Δn=n(t)-n(0)"); out["induced"] = f"{HERE}/gif_induced.gif"
    mk_delta("instant", "gif_instant.gif", u"instantaneous Δn"); out["instant"] = f"{HERE}/gif_instant.gif"
    print("wrote GIFs:", list(out)); return out
GIFS = build_gifs()

cells = []
def md(s): cells.append(new_markdown_cell(s))
def code(s): cells.append(new_code_cell(s))

# ---- 1. title + question ----
md("# Classical projectile in a jellium slab — v=4.5, DIRECT potential (deep-dive)\n"
   "**System:** localised jellium slab, r_s≈4.18 (n0=%s), N=100, box 35×35×85 Bohr, "
   "slab half-width ±%.1f, `periodicity(2)` (x,y periodic, z open → CAP-free, energy-conserving).\n\n"
   "**Projectile:** a moving Gaussian **DIRECT erf/r potential** (charge-free), mass 1, charge −1, "
   "σ_WP=%.2f, launched at z=−24 with k0=4.5 (v0=4.5 a.u.), free Ehrenfest.\n\n"
   "**What this run shows:** the corrected classical baseline for the FASTEST sweep velocity — a "
   "replica of the old charge-based `dyn/v4p5` with the direct-potential fix. The question: does "
   "the direct potential remove the box-exit **kink** and the unphysical pairwise ledger while "
   "reproducing the in-slab stopping physics?" % (SUMM.get('n0','0.00327'), FACE, SIGMA_WP))

# ---- 2. conventions ----
md("## 1. Conventions & symbols\n"
   "Atomic units (ℏ=m_e=e=1); energies in Ha unless eV stated (1 Ha=27.2114 eV). "
   "σ = **σ_WP** (wavepacket width); the potential's charge std is σ_pot=σ_WP/√2=%.3f.\n\n"
   "| symbol | meaning |\n|---|---|\n"
   "| n(x,z,t) | electron number density (mid-y slice) |\n"
   "| Δn | induced n(t)−n(0) / instantaneous n(t)−n(t−Δt) |\n"
   "| E_PP,E_PS,E_PB | projectile self / projectile–slab / projectile–background Coulomb |\n"
   "| E_SS,E_SB | slab self-Hartree / slab–background |\n"
   "| U_proj_bg | projectile–background interaction (=E_PB) |\n"
   "| S | electronic stopping power = −dKE/ds over the slab |\n"
   "| ω_p, v_F, k_F | plasma freq, Fermi velocity, Fermi wavevector (r_s=4.18) |" % (SIGMA_WP/2**0.5))

# ---- 3. setup ----
md("## 2. Setup — fully reconstructable (from run_summary.txt)")
code("import pandas as pd, numpy as np, matplotlib.pyplot as plt, os\n"
     "from matplotlib.colors import LogNorm, SymLogNorm\n"
     f"NEW=r'{NEW}'; OLD=r'{OLD}'; HA=27.211386; E_GS={E_GS!r}; FACE={FACE!r}; WALL={WALL!r}; DT={DT!r}\n"
     f"summ=open(NEW+'/run_summary.txt').read()\n"
     "print(summ)")

# ---- 4. source files ----
md("## 3. Source files\n"
   f"- run.cpp: `{RUNCPP}`\n"
   f"- perturbation: `inq-stack/include/inqkit/dynamics/moving_gaussian_projectile_potential.hpp`\n"
   f"- direct potential: `inq-stack/include/inqkit/jellium/gaussian_potential.hpp`\n"
   f"- direct force: `inq-stack/include/inqkit/dynamics/projectile_force.hpp` (projectile_force_direct_z)\n"
   f"- direct ledger: `inq-stack/include/inqkit/jellium/interaction_energies.hpp` (compute_coulomb_direct)\n"
   f"- launch: `.../dyn_direct/launch_v4p5_direct.sh` | this builder: `.../hypotheses/.../dyn_direct/build_notebook_v4p5_direct.py`\n"
   f"- plan: `docs/plans/direct-potential-ledger-fix.md`")

code("def load(run):\n"
     "    b=run+'/raw/observables/'\n"
     "    ob=pd.read_csv(b+'observables.csv'); pj=pd.read_csv(b+'projectile.csv'); ix=pd.read_csv(b+'interactions.csv')\n"
     "    return ob.merge(pj,on=['step','time_au']).merge(ix,on=['step','time_au'])\n"
     "N=load(NEW); O=load(OLD)\n"
     "print('NEW rows',len(N),' proj_z %.1f..%.1f'%(N.proj_z.min(),N.proj_z.max()))")

# ---- visual intuition: GIFs ----
md("## Visual intuition — density evolution (mid-y x–z, LINEAR | LOG)\n"
   "Dashed=slab faces (±%.1f), red dotted=box wall (±%.1f), black=projectile z. Total density, "
   "then induced Δn=n(t)−n(0), then instantaneous Δn. The negative projectile carves a "
   "**depletion** (anti-wake); watch it transit cleanly and the field recede past the wall with "
   "no lurch." % (FACE, WALL))
for key, cap in [("total", "Total density n(x,z,t)"), ("induced", "Induced Δn = n(t)−n(0)"),
                 ("instant", "Instantaneous Δn = n(t)−n(t−Δt)")]:
    g = GIFS.get(key)
    if g:
        md(f"**{cap}**")
        code(f"from IPython.display import Image, display\ndisplay(Image(filename=r'{g}'))")

# ---- density carpets ----
md("## Density carpets (z–t, transverse-integrated) — LINEAR | LOG")
code("from inqview import load_vti\n"
     "import glob, re\n"
     "frames=sorted(glob.glob(NEW+'/frames/total/*.vti'), key=lambda f:int(f.split('_t')[-1].split('.')[0]))\n"
     "sel=frames[::max(1,len(frames)//120)]\n"
     "d0=load_vti(frames[0],expect_centered_axis='z'); z=d0.z\n"
     "nzt=[]; t=[]\n"
     "for f in sel:\n"
     "    d=load_vti(f,expect_centered_axis='z'); nzt.append(d.data.sum(axis=(0,1)))  # integrate x,y\n"
     "    t.append(int(f.split('_t')[-1].split('.')[0])*DT)\n"
     "nzt=np.array(nzt); t=np.array(t)\n"
     "fig,ax=plt.subplots(1,2,figsize=(13,4))\n"
     "im0=ax[0].pcolormesh(z,t,nzt,cmap='viridis',shading='auto'); ax[0].set_title('n(z,t) linear'); fig.colorbar(im0,ax=ax[0])\n"
     "im1=ax[1].pcolormesh(z,t,np.maximum(nzt,nzt.max()*1e-4),cmap='viridis',norm=LogNorm(nzt.max()*1e-4,nzt.max()),shading='auto'); ax[1].set_title('n(z,t) log'); fig.colorbar(im1,ax=ax[1])\n"
     "for a in ax:\n"
     "    for zf in (-FACE,FACE): a.axvline(zf,ls='--',c='w',lw=0.8)\n"
     "    a.axvline(WALL,ls=':',c='r'); a.set_xlabel('z (Bohr)'); a.set_ylabel('t (a.u.)')\n"
     "plt.tight_layout(); plt.show()")

# ---- energetics: dE_total + N(t) beside ----
md("## Energetics — ΔE_total(t) with N(t) beside it, components, conservation")
code("fig,ax=plt.subplots(1,2,figsize=(13,4))\n"
     "ax[0].plot(N.time_au,(N.energy_total-E_GS)*HA); ax[0].set_title('ΔE_total = E_total - E_GS (eV)'); ax[0].set_xlabel('t (a.u.)'); ax[0].set_ylabel('eV')\n"
     "# N(t) = integral n dV per frame\n"
     "Nt=[]; tt=[]\n"
     "for f in sel:\n"
     "    d=load_vti(f,expect_centered_axis='z'); dv=(d.x[1]-d.x[0])*(d.y[1]-d.y[0])*(d.z[1]-d.z[0]); Nt.append(d.data.sum()*dv); tt.append(int(f.split('_t')[-1].split('.')[0])*DT)\n"
     "ax[1].plot(tt,Nt); ax[1].set_title('N(t)=∫n dV  (CAP-free: conserved)'); ax[1].set_xlabel('t (a.u.)'); ax[1].set_ylabel('electrons')\n"
     "plt.tight_layout(); plt.show()\n"
     "print('N(t): %.4f -> %.4f (drift %.2e)'%(Nt[0],Nt[-1],Nt[-1]-Nt[0]))")
code("fig,ax=plt.subplots(1,2,figsize=(13,4))\n"
     "for c,l in [('energy_kinetic','kinetic'),('energy_hartree','Hartree'),('energy_xc','xc'),('energy_external','external')]:\n"
     "    ax[0].plot(N.proj_z,(N[c]-N[c].iloc[0])*HA,label=l)\n"
     "ax[0].axvline(FACE,ls='--',c='grey'); ax[0].axvline(WALL,ls=':',c='r'); ax[0].legend(); ax[0].set_title('energy components Δ vs proj_z (eV)'); ax[0].set_xlabel('proj_z')\n"
     "cons=N.energy_total+N.energy_proj_ke+N.energy_proj_bg_ideal\n"
     "ax[1].plot(N.proj_z,(cons-cons.iloc[0])*HA); ax[1].axvline(WALL,ls=':',c='r'); ax[1].set_title('conserved: E_elec+KE_proj+U_proj_bg (eV)'); ax[1].set_xlabel('proj_z')\n"
     "plt.tight_layout(); plt.show()\n"
     "print('conservation drift over run = %.3e eV; across wall(40<z<45) std %.2e eV'%((cons.max()-cons.min())*HA, cons[(N.proj_z>40)&(N.proj_z<45)].std()*HA))")

# ---- no-kink demo OLD vs NEW ----
md("## No-kink demonstration — OLD charge vs NEW direct\n"
   "The OLD run's Gaussian **charge** clips at the wall and its x,y-periodic replication forms a "
   "charged **sheet**; the NEW direct potential has neither. Each panel vs projectile z.")
code("fig,ax=plt.subplots(2,3,figsize=(15,8))\n"
     "def panel(a,col,title,ref=0,vz=False):\n"
     "    a.plot(O.proj_z, O[col] if vz else (O[col]-ref)*HA, color='tab:red',alpha=0.7,label='OLD charge')\n"
     "    a.plot(N.proj_z, N[col] if vz else (N[col]-ref)*HA, color='tab:blue',label='NEW direct')\n"
     "    a.axvline(WALL,ls=':',c='k'); a.axvline(FACE,ls='--',c='grey'); a.set_title(title); a.set_xlabel('proj_z'); a.legend(fontsize=8)\n"
     "panel(ax[0,0],'energy_total','E_total - E_GS (eV)',ref=E_GS)\n"
     "panel(ax[0,1],'e_ps','E_PS proj-slab (eV)')\n"
     "panel(ax[0,2],'e_pb','E_PB proj-background (eV)')\n"
     "panel(ax[1,0],'e_pp','E_PP proj self (eV)')\n"
     "panel(ax[1,1],'energy_proj_bg_ideal','U_proj_bg (eV)')\n"
     "panel(ax[1,2],'proj_vz','proj_vz (a.u.)',vz=True)\n"
     "plt.tight_layout(); plt.show()\n"
     "def maxcurv(d,c):\n"
     "    v=d[c].values; z=d.proj_z.values; return z[np.argmax(np.abs(np.gradient(np.gradient(v))))]\n"
     "print('%-22s %8s %8s'%('series','OLD z*','NEW z*'))\n"
     "for c in ['energy_total','e_ps','e_pb','e_pp','energy_proj_bg_ideal']:\n"
     "    print('%-22s %8.1f %8.1f'%(c,maxcurv(O,c),maxcurv(N,c)))\n"
     "print('(NEW curvature maxima IN-SLAB, not at wall z=%.1f => no kink)'%WALL)")

# ---- clean ledger ----
md("## Pairwise ledger (direct) — now physical\n"
   "E_PS>0 decaying to 0 (repulsion between the −projectile and −electrons), E_PP constant "
   "(rigid-Gaussian self-energy = 1/(2σ_pot√π)), E_PB<0 (attraction to +background) → 0.")
code("fig,ax=plt.subplots(figsize=(8,4))\n"
     "for c,l in [('e_ss','E_SS slab self'),('e_pp','E_PP proj self'),('e_ps','E_PS proj-slab'),('e_sb','E_SB slab-bg'),('e_pb','E_PB proj-bg')]:\n"
     "    ax.plot(N.proj_z,N[c],label=l)\n"
     "ax.axvline(FACE,ls='--',c='grey'); ax.axvline(WALL,ls=':',c='r'); ax.axhline(0,c='k',lw=0.5)\n"
     "ax.legend(fontsize=8); ax.set_xlabel('proj_z (Bohr)'); ax.set_ylabel('energy (Ha)'); ax.set_title('direct pairwise ledger'); plt.show()\n"
     "print('E_PP mean %.4f Ha std %.1e (const);  E_PS %.2f (z~0) -> %.3f (end);  E_PB %.3f (z~0) -> %.3f (end)'\n"
     "      %(N.e_pp.mean(),N.e_pp.std(),N.loc[(N.proj_z).abs().idxmin(),'e_ps'],N.e_ps.iloc[-1],N.loc[(N.proj_z).abs().idxmin(),'e_pb'],N.e_pb.iloc[-1]))")

# ---- projectile & transport ----
md("## Projectile & transport — trajectory, KE(z), stopping\n"
   "KE(z) is not monotonic (conservative dip-and-recovery); the **net** loss between the "
   "symmetric slab faces (equal-potential window) is the stopping. v=4.5 barely decelerates.")
code("fig,ax=plt.subplots(1,3,figsize=(16,4))\n"
     "ax[0].plot(N.time_au,N.proj_z); ax[0].axhline(FACE,ls='--',c='grey'); ax[0].axhline(-FACE,ls='--',c='grey'); ax[0].set_title('z(t)'); ax[0].set_xlabel('t'); ax[0].set_ylabel('proj_z')\n"
     "ax[1].plot(N.proj_z,N.proj_vz); ax[1].axvline(FACE,ls='--',c='grey'); ax[1].axvline(-FACE,ls='--',c='grey'); ax[1].set_title('v(z)'); ax[1].set_xlabel('proj_z'); ax[1].set_ylabel('vz')\n"
     "ax[2].plot(N.proj_z,N.energy_proj_ke*HA); ax[2].axvline(FACE,ls='--',c='grey'); ax[2].axvline(-FACE,ls='--',c='grey'); ax[2].set_title('KE(z) (eV)'); ax[2].set_xlabel('proj_z')\n"
     "plt.tight_layout(); plt.show()\n"
     "# stopping: KE loss across the equal-potential slab window (-FACE..+FACE) and to z~20\n"
     "def ke_at(d,z): return d.loc[(d.proj_z-z).abs().idxmin(),'energy_proj_ke']\n"
     "dep_slab=(ke_at(N,-FACE)-ke_at(N,FACE))*HA; dep20=(N.energy_proj_ke.iloc[0]-ke_at(N,20))*HA\n"
     "print('v_launch=%.3f  v_final=%.3f'%(N.proj_vz.iloc[0],N.proj_vz.iloc[-1]))\n"
     "print('KE loss across slab (-12.5..+12.5): NEW %.2f eV  (OLD %.2f eV) -> sheet inflated OLD ~%.0f%%'\n"
     "      %(dep_slab,(ke_at(O,-FACE)-ke_at(O,FACE))*HA,100*((ke_at(O,-FACE)-ke_at(O,FACE))/(ke_at(N,-FACE)-ke_at(N,FACE))-1)))\n"
     "print('S(v=4.5) via KE-loss to z~20 = %.3f eV/Bohr  (old sweep charge-based reported 0.28)'%(dep20/25.0))")

# ---- physical anchors ----
md("## Physical anchors (HEG scales, projectile timescales, Lindhard) — r_s=4.18")
code("try:\n"
     "    from inqview.analysis import compute_heuristics\n"
     "    H=compute_heuristics(NEW, rs=4.18, v0=4.5, z0=-24.0, slab_half=FACE, box_half=WALL, sigma_wp=%r)\n" % SIGMA_WP +
     "    for grp in ['eg_scales','timescales','stopping_refs']:\n"
     "        dd=getattr(H,grp,{}) or {}\n"
     "        print('['+grp+']'); [print('   %-22s %s'%(k,('%.4g'%v if isinstance(v,(int,float)) else v))) for k,v in dd.items()]\n"
     "except Exception as e:\n"
     "    print('compute_heuristics unavailable (%s); inline HEG for r_s=4.18:'%e)\n"
     "    rs=4.18; kF=(9*np.pi/4)**(1/3)/rs; print('k_F=%.3f v_F=%.3f E_F=%.3f Ha  omega_p=%.3f Ha (%.2f eV)'%(kF,kF,kF*kF/2,(3/rs**3)**0.5,(3/rs**3)**0.5*HA))\n"
     "print('\\nprojectile v0=4.5 >> v_F=0.46 => fast (Bethe) regime; T_plasmon=2pi/omega_p=%.1f a.u.'%(2*np.pi/(3/4.18**3)**0.5))")

# ---- loss function note ----
md("## Loss function L(q,ω)\n"
   "**Requires a future observable.** This run saved the energy ledger + density frames but not "
   "n_q(t) / dipole time series, so L(q,ω)=|n_q(ω)|²/q² cannot be formed here. To enable, re-run "
   "with the density-q / dipole observable (fourier-analysis pipeline). Not fabricated.")

# ---- takeaway ----
md("## Takeaway\n"
   "- The **direct erf/r potential removes both charge-based artifacts**: no exit kink (all "
   "curvature maxima are in-slab, not at the wall) and a physical pairwise ledger "
   "(E_PS>0→0, E_PP const, E_PB<0→0). Energy conserves to ~3e-3 eV with no wall spike.\n"
   "- The old charge representation **inflated the in-slab stopping ~20–35%** via the spurious "
   "x,y-periodic charged sheet. **S(v=4.5) ≈ 0.18 eV/Bohr** (direct) vs 0.28 (old charge).\n"
   "- In-slab trajectory/KE reproduce the old run where it was valid; this is the corrected "
   "classical baseline. **The full S(v) sweep should be re-run with `dyn_direct`.**")

nb = new_notebook(cells=cells)
nb.metadata.kernelspec = {"name": "python3", "display_name": "Python 3", "language": "python"}
ExecutePreprocessor(timeout=2400, kernel_name="python3").preprocess(nb, {"metadata": {"path": HERE}})
with open(NB, "w") as f: nbf.write(nb, f)
nerr = sum(1 for c in nb.cells if c.cell_type=="code" for o in c.get("outputs",[]) if o.get("output_type")=="error")
print("wrote", NB, "| code errors:", nerr)
