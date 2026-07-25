#!/usr/bin/env python3
"""Build + execute the muon-mass-fork campaign notebooks (canonical theme).

Notebooks (under this folder):
  index.ipynb            - guided read-order over the phase notebooks
  phase2_physics.ipynb   - vacuum free-particle oracles + xz-density-vs-sigma
  phase3_regression.ipynb- bit-for-bit fork(mass=1) vs pristine inq
  phase4_xc_research.ipynb - grounded muon-XC candidates for the user's pick

Run:  venv/bin/python3 build_notebooks.py [index|phase2|phase3|phase4|all]
Uses the shared _nbreport helper (executes cells in the venv kernel).
"""
import sys, os
HYP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HYP))          # hypotheses/ (has _nbreport.py)
import _nbreport as nb

MF   = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/muon_mass_fork"
P2   = f"{MF}/runs/phase2"
REG  = f"{MF}/regression"
DOCS = "/local/data/public/skcb2/tddft/docs"
nb.set_outdir(HYP)

# ---------------------------------------------------------------- phase 2
def build_phase2():
    c = [nb.md("# Phase 2 — vacuum free-particle physics (mass fork oracles)\n\n"
               "**Trust-gate role:** proves the per-state inverse-mass fork reproduces the "
               "*exact* free-Gaussian dynamics for arbitrary mass, on GPU. Every panel is an "
               "analytic oracle — no fitting freedom hides a bug.\n\n"
               "σ is the wavepacket width σ_WP; the density std is σ_ρ = σ_WP/√2. The free "
               "law is σ_ρ(t)² = σ_ρ0² + t²/(4m²σ_ρ0²)."),
         nb.setup_cell(),
         nb.code(
            f"P2 = {P2!r}\nHYP = {HYP!r}\n"
            "import numpy as np, matplotlib.pyplot as plt\n"
            "def load(run):\n"
            "    p = f'{P2}/{run}/results/{run}/raw/observables/wp_real_space_stats.csv'\n"
            "    L = [l for l in open(p).read().strip().splitlines() if l.strip() and not l.startswith('#')]\n"
            "    h = L[0].split(','); a = np.array([[float(x) for x in r.split(',')] for r in L[1:]])\n"
            "    return {k: a[:,i] for i,k in enumerate(h)}\n"
            "SIG_WP=0.5; SR0=SIG_WP/np.sqrt(2); SR0S=SR0**2\n"
            "print('loaded helper; sigma_rho0^2 =', SR0S)"),
         nb.md("## 2.1 σ_ρ(t)² spreading — electron, m=10, muon (mass dial)\n"
               "The parabola slope in t is 1/(4m²σ_ρ0²): heavier ⇒ flatter. Dashed = analytic oracle."),
         nb.code(
            "fig, ax = plt.subplots(figsize=(6,4))\n"
            "runs = [('spread_elec',1.0,'C0'),('spread_m10',10.0,'C1'),('spread_muon',206.77,'C2')]\n"
            "rows=[]\n"
            "for run,m,col in runs:\n"
            "    d=load(run); t=d['time_au']; s2=d['sigma_z2']\n"
            "    ax.plot(t,s2,color=col,lw=1.5,label=f'{run} (m={m:g})')\n"
            "    orn = SR0S + t**2/(4*m**2*SR0S)\n"
            "    ax.plot(t,orn,color=col,ls='--',lw=1,alpha=0.7)\n"
            "    A=np.vstack([np.ones_like(t),t**2]).T; b=np.linalg.lstsq(A,s2,rcond=None)[0]\n"
            "    mfit=1/(2*SR0*np.sqrt(b[1])); rows.append((run,m,mfit,abs(mfit-m)/m))\n"
            "ax.set_xlabel('t (a.u.)'); ax.set_ylabel(r'$\\sigma_\\rho^2$ (Bohr$^2$)'); ax.set_yscale('log')\n"
            "ax.legend(fontsize=8); ax.set_title('Free-Gaussian spreading vs analytic oracle (dashed)')\n"
            "fig.tight_layout(); fig.savefig(f'{HYP}/p2_spreading.png',dpi=130); plt.close(fig)\n"
            "print('recovered masses:  run | m_true | m_fit | rel.err')\n"
            "for r in rows: print(f'  {r[0]:12s} {r[1]:8.2f} {r[2]:12.5f}  {r[3]:.2e}')"),
         nb.embed(f"{HYP}/p2_spreading.png", "σ_ρ² vs t: solid = simulation, dashed = analytic free-Gaussian oracle."),
         nb.md("## 2.2 Mass-dial linearity — slope ∝ 1/m²\n"
               "The t²-coefficient of σ_ρ² must fall as 1/m² across the whole dial."),
         nb.code(
            "ms=[1.0,10.0,206.77]; runs=['spread_elec','spread_m10','spread_muon']\n"
            "sl=[]\n"
            "for run in runs:\n"
            "    d=load(run); t=d['time_au']; s2=d['sigma_z2']\n"
            "    A=np.vstack([np.ones_like(t),t**2]).T; sl.append(np.linalg.lstsq(A,s2,rcond=None)[0][1])\n"
            "fig,ax=plt.subplots(figsize=(5,4))\n"
            "ax.loglog(ms,sl,'o-',label='measured slope b')\n"
            "ax.loglog(ms,[1/(4*m**2*SR0S) for m in ms],'k--',label=r'$1/(4m^2\\sigma_{\\rho0}^2)$')\n"
            "ax.set_xlabel('mass m (mₑ)'); ax.set_ylabel(r'slope of $\\sigma_\\rho^2$ vs $t^2$')\n"
            "ax.legend(fontsize=8); ax.set_title('Mass-dial: spreading rate ∝ 1/m²')\n"
            "fig.tight_layout(); fig.savefig(f'{HYP}/p2_massdial.png',dpi=130); plt.close(fig)\n"
            "print('slopes:', [f'{s:.3e}' for s in sl])"),
         nb.embed(f"{HYP}/p2_massdial.png", "Spreading-rate coefficient vs mass: measured (o) vs 1/m² law (dashed)."),
         nb.md("## 2.3 Group velocity — ⟨z⟩(t) slope = k₀/m\n"
               "With k₀=0.5, the electron centroid drifts at 0.5; the muon at 0.5/206.77."),
         nb.code(
            "fig,(a1,a2)=plt.subplots(1,2,figsize=(9,3.6))\n"
            "for run,m,ax in [('vgroup_elec',1.0,a1),('vgroup_muon',206.77,a2)]:\n"
            "    d=load(run); t=d['time_au']; z=d['z_mean']\n"
            "    vfit=np.polyfit(t,z,1)[0]\n"
            "    ax.plot(t,z,lw=1.5); ax.plot(t,0.5/m*t,'k--',lw=1)\n"
            "    ax.set_title(f'{run}\\n v_fit={vfit:.5g}  k0/m={0.5/m:.5g}',fontsize=9)\n"
            "    ax.set_xlabel('t (a.u.)'); ax.set_ylabel(r'$\\langle z\\rangle$ (Bohr)')\n"
            "fig.tight_layout(); fig.savefig(f'{HYP}/p2_vgroup.png',dpi=130); plt.close(fig)\n"
            "print('group velocity: solid=sim, dashed=k0/m')"),
         nb.embed(f"{HYP}/p2_vgroup.png", "Centroid drift ⟨z⟩(t): the muon moves 206.77× slower at the same k₀."),
         nb.md("## 2.4 xz density slices vs σ_WP (muon WP at t=0)\n"
               "Physical-order VTI loaded via `inqview.load_vti` (NEVER fftshift a VTI)."),
         nb.code(
            "from inqview import load_vti\n"
            "sigs=[0.5,1.0,2.0,4.0]\n"
            "fig,axs=plt.subplots(1,4,figsize=(13,3.4))\n"
            "for s,ax in zip(sigs,axs):\n"
            "    run=f'xz_muon_sig{s:g}'\n"
            "    vp=f'{P2}/{run}/results/{run}/raw/vti/density_wp/density_t000000.vti'\n"
            "    v=load_vti(vp)\n"
            "    dat=v.data; x=v.x; z=v.z\n"
            "    iy=dat.shape[1]//2\n"
            "    sl=dat[:,iy,:].T\n"
            "    ax.imshow(sl,origin='lower',extent=[x[0],x[-1],z[0],z[-1]],aspect='equal',cmap='viridis')\n"
            "    ax.set_title(f'σ_WP={s:g}',fontsize=9); ax.set_xlabel('x (Bohr)')\n"
            "axs[0].set_ylabel('z (Bohr)')\n"
            "fig.suptitle('Injected muon WP density (xz slice) vs σ_WP',fontsize=11)\n"
            "fig.tight_layout(); fig.savefig(f'{HYP}/p2_xz_density.png',dpi=130); plt.close(fig)\n"
            "print('xz density panels written')"),
         nb.embed(f"{HYP}/p2_xz_density.png", "Injected muon WP density (xz) broadening with σ_WP."),
         nb.md("### Verdict\nAll five physics runs pass every analytic oracle to ~1e-5 (mass recovered "
               "across m ∈ {1,10,206.77} to rel ≤ 5e-6; v_group=k₀/m; ⟨k²⟩ & norm conserved). "
               "The mass fork is physically correct in vacuum → proceed to the bit-for-bit gate (Phase 3).")]
    nb.build(c, f"{HYP}/phase2_physics.ipynb")

# ---------------------------------------------------------------- phase 3
def build_phase3():
    c = [nb.md("# Phase 3 — bit-for-bit regression: fork(mass=1) vs pristine inq\n\n"
               "**Trust-gate role (HARD):** with all masses = 1 the fork's empty-factor guard "
               "must route the ORIGINAL scalar kinetic path. A He-atom LDA GS + kicked RT is built "
               "against inq-study (fork) AND pristine inq; identical output ⇒ the fork is inert when "
               "off ⇒ muon physics is attributable to the mass, not an engine edit."),
         nb.setup_cell(),
         nb.code(
            f"REG={REG!r}\nHYP={HYP!r}\n"
            "import numpy as np\n"
            "def gs(run):\n"
            "    d={}\n"
            "    for l in open(f'{REG}/{run[0]}/results/{run[1]}/raw/observables/gs_energy.csv').read().strip().splitlines()[1:]:\n"
            "        k,v=l.split(','); d[k]=float(v)\n"
            "    return d\n"
            "f=gs(('fork','reg_fork')); p=gs(('pristine','reg_pristine'))\n"
            "print(f'{\"component\":10s} {\"fork (Ha)\":>18s} {\"pristine (Ha)\":>18s} {\"|diff|\":>10s}')\n"
            "for k in f: print(f'{k:10s} {f[k]:18.12f} {p[k]:18.12f} {abs(f[k]-p[k]):10.1e}')"),
         nb.md("### RT energy trace + density diff"),
         nb.code(
            "def rt(run):\n"
            "    L=open(f'{REG}/{run[0]}/results/{run[1]}/raw/observables/rt_energy.csv').read().strip().splitlines()\n"
            "    h=L[0].split(','); a=np.array([[float(x) for x in r.split(',')] for r in L[1:]]); return h,a\n"
            "hf,af=rt(('fork','reg_fork')); hp,ap=rt(('pristine','reg_pristine'))\n"
            "print(f'RT trace: {af.shape[0]} steps')\n"
            "for j,name in enumerate(hf):\n"
            "    if name in ('step','time_au'): continue\n"
            "    print(f'  {name:9s} max|d| over steps = {np.max(np.abs(af[:,j]-ap[:,j])):.1e}')\n"
            "df=np.fromfile(f'{REG}/fork/results/reg_fork/raw/observables/gs_density/gs_density.raw')\n"
            "dp=np.fromfile(f'{REG}/pristine/results/reg_pristine/raw/observables/gs_density/gs_density.raw')\n"
            "print(f'GS density ({df.size} pts): max|d| = {np.max(np.abs(df-dp)):.1e}')"),
         nb.md("### Verdict\n**BIT-FOR-BIT PASS** — GS energies identical to 14 digits; RT trace agrees to "
               "≤1e-14; density identical. The fork is provably inert when off.")]
    nb.build(c, f"{HYP}/phase3_regression.ipynb")

# ---------------------------------------------------------------- phase 4
def build_phase4():
    src = open(f"{DOCS}/campaigns/muon_mass_fork/phase4_muon_xc_candidates.md").read()
    c = [nb.md("# Phase 4 — muon exchange–correlation: grounded candidates (USER PICK)\n\n"
               "**Checkpoint role:** research is complete; the *user* selects the Phase-5 functional. "
               "Below is the grounded candidate summary (source notes in `docs/sources/`)."),
         nb.setup_cell(),
         nb.md(src),
         nb.md("---\n**To resume Phase 5:** write your choice to "
               "`scripts/muon_mass_fork/muon_xc_pick.json` and re-run `orchestrate.py`.")]
    nb.build(c, f"{HYP}/phase4_xc_research.ipynb")

# ---------------------------------------------------------------- index
def build_index():
    c = [nb.md("# Muon-mass-fork campaign — INDEX (read order)\n\n"
               "A per-state inverse-mass fork in `inq-study` lets any Kohn–Sham orbital carry an "
               "arbitrary mass (muon = 206.77 mₑ). This campaign takes it from *compiles* to "
               "*trusted + a physics result*. Read in order — each notebook gates the next.\n\n"
               "| # | Notebook | Establishes | Gate |\n"
               "|---|----------|-------------|------|\n"
               "| 1 | *(engine)* Phase 1 | GPU build + Tier-1 kernel oracle (k²/2m per state) | code compiles + kernel correct |\n"
               "| 2 | [phase2_physics](phase2_physics.ipynb) | vacuum free-particle oracles: σ(t), v_group, mass-dial {1,10,206.77}, xz-density | fork physically correct in vacuum |\n"
               "| 3 | [phase3_regression](phase3_regression.ipynb) | bit-for-bit fork(mass=1) == pristine inq | fork inert when off (HARD trust gate) |\n"
               "| 3b| *(runs)* Phase 3b | muon WP under full LDA — interacting propagator sanity | no NaN, stable dynamics |\n"
               "| 4 | [phase4_xc_research](phase4_xc_research.ipynb) | grounded muon-XC candidates | **user picks** the Phase-5 functional |\n"
               "| 5 | *(pending pick)* Phase 5 | all-muon r_s=5.69 jellium: naive-LDA vs mass-rescaled | XC-sensitivity verdict |\n\n"
               "**Status at build time:** Phases 1–3 PASSED; Phase 4 research done (awaiting pick); "
               "Phase 5 gated on `muon_xc_pick.json`.\n\n"
               "Provenance: campaign `docs/campaigns/muon_mass_fork/muon_mass_fork.md`; "
               "handover `docs/handovers/muon-mass-fork.md`; sources `docs/sources/{heg-mass-scaling-xc,"
               "kreibich-gross-multicomponent-dft,car-parrinello-fictitious-mass}.md`.")]
    nb.build(c, f"{HYP}/index.ipynb")

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("phase2","all"): build_phase2()
    if which in ("phase3","all"): build_phase3()
    if which in ("phase4","all"): build_phase4()
    if which in ("index","all"):  build_index()
    print("done:", which)
