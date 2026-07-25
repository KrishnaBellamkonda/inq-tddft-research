#!/usr/bin/env python3
"""Build hypotheses/perturbation_method/perturbation_method_study.ipynb — the
energy-ledger + stress-test analysis of the Gaussian-charge PERTURBATION method
for representing a classical projectile (campaign localised-jellium-dynamics-analysis).

Reads the CSVs produced by proj_perturbation/perturbation_stress.py (stress_*.csv,
this folder) + proj_perturbation/grid_sweep.csv. Partial-tolerant. Executes to 0
errors via ExecutePreprocessor and writes the notebook with outputs.
"""
import sys, os
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = Path(__file__).resolve().parent
GRID = HERE.parent.parent/"scripts/localised_jellium_dynamics/proj_perturbation/grid_sweep.csv"
OUT  = HERE/"perturbation_method_study.ipynb"
cells = []
def md(s, anchor=None):
    c = new_markdown_cell(s); c.metadata["gen"]="builder"
    if anchor: c.metadata["anchor"]=anchor
    cells.append(c)
def code(s):
    c = new_code_cell(s); c.metadata["gen"]="builder"; cells.append(c)

# ── 1. Title + question ──────────────────────────────────────────────────────
md("# Localised jellium — the Gaussian-charge **perturbation** projectile: energy ledger + stress tests\n\n"
   "**Question.** Can we represent a *classical* projectile not as a UPF pseudopotential ion (whose "
   "long-range erf/r ghost tail aliases with the radial cutoff) but as a **stationary Gaussian charge "
   "added to the KS potential via its Poisson potential** — and is that representation *accurate*?\n\n"
   "**Accuracy criterion (locked with the user).** (1) a fully-accounted energy ledger in which the "
   "measured `d(E_H+E_ext) − U_proj_bg` sits close to the analytic wavepacket Hartree self-energy; "
   "(2) the LDA self-interaction error falls out of it; (3) the result is stable against grid spacing and "
   "(for the pseudopotential contrast) radial cutoff. We confirm this against a WP reference and then "
   "stress-test it across σ, projectile position r, box length Lz, and periodicity (open-z vs full PBC).\n\n"
   "Sibling notebook `../localised_jellium_dynamics/ledger_rcut.ipynb` established that the *pseudopotential* "
   "residual (7.4 eV at r_cut=120) was a representation artifact; this notebook validates the clean method.")

# ── 2. Method / formulas (conventions + symbol table up front) ───────────────
md("## Method — conventions and symbols\n\n"
   "Atomic units (Hartree). 1 Ha = 27.211386 eV. INQ total energy = "
   "$E_\\text{tot}=E_\\text{kin}+E_\\text{H}+E_\\text{xc}+E_\\text{ext}$ (8-term; nonlocal/ion vanish here).\n\n"
   "| symbol | meaning | value / range |\n|---|---|---|\n"
   "| $\\sigma_{WP}$ | wavepacket (amplitude) width | swept {0.35, 0.5, 0.7, 1.0} Bohr |\n"
   "| $\\sigma_\\rho=\\sigma_{WP}/\\sqrt2$ | projectile **charge/density** std | = σ_pot |\n"
   "| $r$ | projectile distance from slab face | swept {4,12,20,28} Bohr |\n"
   "| $L_z$ | box length (open or periodic z) | swept {90,120,160,240} Bohr |\n"
   "| periodicity | 2 = periodic x,y + open z; 3 = full PBC | {2,3} |\n"
   "| $n_+$ | positive jellium background (∫=N=82) | slab \\|z\\|<12.5, n₀=1.312e-3 |\n"
   "| $n_\\text{proj}$ | projectile Gaussian charge (∫=1) | std σ_ρ, at (0,0,−(12.5+r)) |\n")

md("### The perturbation potential (at its point of use)\n"
   "The projectile is a −1 electron ⇒ *repulsive* to electrons. Its potential energy for an electron is "
   "$v_\\text{proj}=+\\,\\text{poisson}(n_\\text{proj})$, added to the KS potential alongside the background "
   "well $v_\\text{bg}=-\\text{poisson}(n_+)$ via `perturbations::sum`. INQ then captures "
   "$\\int n_e\\,v_\\text{proj}$ inside $E_\\text{ext}$ automatically. No UPF, no ion, no r_cut.")
md("### The projectile↔background diagnostic $U_\\text{proj,bg}$\n"
   "Absent from INQ's total (the projectile is not an Ewald-coupled ion). Computed cleanly as "
   "$U_\\text{proj,bg}=-\\int n_\\text{proj}\\,\\varphi_+$ with $\\varphi_+=\\text{poisson}(n_+)$ — the "
   "r_cut-invariant 'ideal' value.")
md("### The ledger residual and the self-Hartree target\n"
   "$$R \\equiv d(E_\\text{H}+E_\\text{ext}) - U_\\text{proj,bg} = "
   "(E_\\text{H}+E_\\text{ext})_\\text{WP} - (E_\\text{H}+E_\\text{ext})_\\text{pert} - U_\\text{proj,bg}$$\n"
   "Claim: $R$ equals the wavepacket's **Hartree self-energy**. Free-space analytic value of a normalised "
   "Gaussian charge of density-std $\\sigma_\\rho$:\n"
   "$$E_\\text{self}=\\frac{1}{2\\,\\sigma_\\rho\\sqrt\\pi}=\\frac{1}{\\sigma_{WP}\\sqrt{2\\pi}}\\ \\text{Ha}.$$\n"
   "The in-cell (periodic, G=0-dropped) value differs by the **charged-cell gauge** (the WP cell is net −1). "
   "Localisation kinetic term (for reference): $E_\\text{loc}=\\tfrac{3}{4\\sigma_{WP}^2}$ Ha. "
   "Self-interaction error: $\\text{SIE}=R+d E_\\text{xc}$ (self-Hartree partly cancelled by self-XC).")

# theory helpers
code("import pandas as pd, numpy as np\nimport matplotlib.pyplot as plt\n"
     "import sys; sys.path.insert(0,'/local/data/public/skcb2/tddft/inq-stack/python')\n"
     "try:\n    from inqview.visualisation import style as _st; _st.apply()\nexcept Exception: pass\n"
     "HA=27.211386\n"
     "def self_hartree_free(sig_wp):\n"
     "    return 1.0/(sig_wp*np.sqrt(2*np.pi))*HA           # eV, free-space analytic\n"
     "def self_hartree_incell(sig_wp, Lz, per, dx=0.5, Lx=50.0, Ly=50.0):\n"
     "    '''periodic FFT self-Hartree of the WP Gaussian, G=0 dropped (INQ convention).'''\n"
     "    sr=sig_wp/np.sqrt(2.0)\n"
     "    nx,ny,nz=int(round(Lx/dx)),int(round(Ly/dx)),int(round(Lz/dx))\n"
     "    x=(np.arange(nx)-nx//2)*(Lx/nx); y=(np.arange(ny)-ny//2)*(Ly/ny); z=(np.arange(nz)-nz//2)*(Lz/nz)\n"
     "    X,Y,Z=np.meshgrid(x,y,z,indexing='ij'); dV=(Lx/nx)*(Ly/ny)*(Lz/nz)\n"
     "    n=np.exp(-(X*X+Y*Y+Z*Z)/(2*sr*sr)); n/= (n.sum()*dV)\n"
     "    kx=2*np.pi*np.fft.fftfreq(nx,d=Lx/nx); ky=2*np.pi*np.fft.fftfreq(ny,d=Ly/ny); kz=2*np.pi*np.fft.fftfreq(nz,d=Lz/nz)\n"
     "    KX,KY,KZ=np.meshgrid(kx,ky,kz,indexing='ij'); G2=KX*KX+KY*KY+KZ*KZ; G2[0,0,0]=np.inf\n"
     "    nk=np.fft.fftn(n)*dV; V=Lx*Ly*Lz\n"
     "    return float(0.5/V*np.sum(4*np.pi*np.abs(nk)**2/G2))*HA   # eV\n"
     "loc=lambda s: 3/(4*s*s)*HA\n"
     "print('checks (sigma_wp=0.5): self_hartree_free=%.2f eV  incell(Lz120,p2)=%.2f eV  E_loc=%.2f eV'\n"
     "      % (self_hartree_free(0.5), self_hartree_incell(0.5,120,2), loc(0.5)))")

# ── 3. Simulation setup ──────────────────────────────────────────────────────
md("## Simulation setup (reconstructable)\n"
   "- **Cell/geometry:** orthorhombic Lx=Ly=50, Lz∈{90,120,160,240} Bohr; slab \\|z\\|<12.5; dx=0.5 (grid "
   "sweep also 0.4, 0.3). Periodicity 2 (open z) baseline; 3 (full PBC) for the gauge test.\n"
   "- **Electronic structure:** LDA; N=82 electrons; jellium background as a perturbation (no ions in the "
   "perturbation runs). GS reused (bare slab, `campaign_autorun/.../gs_p2_lz120` and "
   "`semiempirical_spillout/runs/lz*`, `p3_lz120`).\n"
   "- **Projectile:** stationary Gaussian charge, std σ_ρ=σ_WP/√2, at (0,0,−(12.5+r)); v_proj=+poisson(n_proj).\n"
   "- **WP reference:** wavepacket at rest (k₀=0), σ_WP, same position; the −1 is a real extra electron "
   "(cell net −1 — the source of the gauge).\n"
   "- **Dynamics:** 2 RT steps (dt=0.01 a.u.), energies read at t=0 (density = GS, unrelaxed — matches the "
   "ledger methodology). Diagnostic only; no production propagation.")

# ── 4. Source files ──────────────────────────────────────────────────────────
md("## Source files\n"
   "| role | path |\n|---|---|\n"
   "| perturbation C++ | `inq-stack/include/inqkit/jellium/gaussian_projectile_perturbation.hpp` |\n"
   "| U_proj_bg / gaussian_density | `inq-stack/include/inqkit/jellium/projectile_background_energy.hpp` |\n"
   "| perturbation run | `scripts/localised_jellium_dynamics/proj_perturbation/run.cpp` |\n"
   "| WP run | `scripts/localised_jellium_dynamics/phase5_wp/run.cpp` |\n"
   "| stress orchestrator | `scripts/localised_jellium_dynamics/proj_perturbation/perturbation_stress.py` |\n"
   "| grid sweep | `scripts/localised_jellium_dynamics/proj_perturbation/grid_sweep.py` |\n"
   "| this builder | `hypotheses/perturbation_method/build_perturbation_report.py` |\n")

# ── 4b. The projectile potential — shape comparison ──────────────────────────
md("## The projectile potential — Gaussian-charge vs pseudopotential vs Coulomb\n"
   "Before the energy ledger, characterise *the potential the projectile makes*. Three radial curves for a "
   "unit ($|Z|=1$) projectile:\n"
   "1. **Gaussian-charge (our perturbation):** $\\varphi_\\text{proj}(r)=\\text{poisson}(n_\\text{proj})$. The "
   "free-space Poisson solution of a normalised Gaussian charge of density-std $\\sigma_\\rho=\\sigma_{WP}/\\sqrt2$ "
   "is the **error-function-smoothed Coulomb**\n"
   "$$\\varphi_\\text{proj}(r)=\\frac{1}{r}\\,\\text{erf}\\!\\left(\\frac{r}{\\sqrt2\\,\\sigma_\\rho}\\right)"
   "=\\frac{\\text{erf}(r/\\sigma_{WP})}{r}\\ \\text{Ha}\\quad(\\sigma_{WP}=0.5\\Rightarrow\\text{erf}(r/0.5)/r).$$\n"
   "2. **Pseudopotential (UPF ghost):** the local potential $V_\\text{loc}(r)$ read verbatim from "
   "`electron_gaussian_wpsigma0p5.upf` (PP_LOCAL, Rydberg → Ha).\n"
   "3. **Exact Coulomb:** $1/r$ Ha.\n\n"
   "Limits of the smoothed form: $r\\gg\\sigma$ → $1/r$ (pure Coulomb); $r\\to0$ → finite plateau "
   "$\\varphi(0)=\\sqrt{2/\\pi}/\\sigma_\\rho=2/(\\sqrt\\pi\\,\\sigma_{WP})=2.26$ Ha (the Gaussian rounds off the "
   "$1/r$ cusp). We also overlay the **actual** periodic-Poisson potential INQ computes (FFT, $4\\pi n_k/G^2$, "
   "$G{=}0$ dropped) to confirm the solver reproduces the analytic form.")
code("import numpy as np, re\n"
     "from scipy.special import erf\n"
     "UPF='/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf'\n"
     "SIGWP=0.5; SIGR=SIGWP/np.sqrt(2)\n"
     "def _blk(tag,txt):\n"
     "    m=re.search(r'<%s[^>]*>(.*?)</%s>'%(tag,tag),txt,re.S); return np.array(m.group(1).split(),float)\n"
     "txt=open(UPF).read(); rU=_blk('PP_R',txt); vU=_blk('PP_LOCAL',txt)/2.0  # Ry->Ha\n"
     "m=rU>0\n"
     "phi_g=lambda r: erf(r/SIGWP)/r                       # our Gaussian-charge potential (Ha)\n"
     "plateau=2.0/(np.sqrt(np.pi)*SIGWP)\n"
     "# actual periodic-Poisson potential INQ computes (radial +z cut through centre)\n"
     "Lx=Ly=50.0; Lz=120.0; dx=0.5\n"
     "nx,ny,nz=[int(round(L/dx)) for L in (Lx,Ly,Lz)]\n"
     "gx=(np.arange(nx)-nx//2)*(Lx/nx); gy=(np.arange(ny)-ny//2)*(Ly/ny); gz=(np.arange(nz)-nz//2)*(Lz/nz)\n"
     "X,Y,Z=np.meshgrid(gx,gy,gz,indexing='ij'); dV=(Lx/nx)*(Ly/ny)*(Lz/nz)\n"
     "ncharge=np.exp(-(X*X+Y*Y+Z*Z)/(2*SIGR*SIGR)); ncharge/=(ncharge.sum()*dV)\n"
     "kx=2*np.pi*np.fft.fftfreq(nx,d=Lx/nx); ky=2*np.pi*np.fft.fftfreq(ny,d=Ly/ny); kz=2*np.pi*np.fft.fftfreq(nz,d=Lz/nz)\n"
     "KX,KY,KZ=np.meshgrid(kx,ky,kz,indexing='ij'); G2=KX*KX+KY*KY+KZ*KZ; G2[0,0,0]=np.inf\n"
     "phi3=np.fft.ifftn(4*np.pi*np.fft.fftn(ncharge)/G2).real\n"
     "z_fft=gz[nz//2:]; phi_fft=phi3[nx//2,ny//2,nz//2:]\n"
     "# quantitative agreement\n"
     "resid_upf=np.max(np.abs(vU[m]-phi_g(rU[m])))\n"
     "r1=rU[m]; onepct=r1[np.argmax(np.abs(phi_g(r1)-1.0/r1)/(1.0/r1)<0.01)]\n"
     "print(f'plateau phi(0)          = {plateau:.4f} Ha  (UPF r->0: {vU[m][0]:.4f})')\n"
     "print(f'max |V_upf - erf/r|     = {resid_upf:.2e} Ha   (identical shapes)')\n"
     "print(f'erf/r within 1%% of 1/r beyond r ~ {onepct:.2f} Bohr')\n"
     "for rr in [0.25,0.5,1,2,5,20]:\n"
     "    i=np.argmin(abs(rU-rr))\n"
     "    print(f'  r={rU[i]:5.2f}  our erf/r={phi_g(rU[i]):+.4f}  V_upf={vU[i]:+.4f}  1/r={1.0/rU[i]:+.4f} Ha')")
code("fig,(a1,a2)=plt.subplots(1,2,figsize=(11.5,4.6))\n"
     "rr=np.linspace(0.02,6,400)\n"
     "a1.plot(rr,1.0/rr,'--',color='0.55',label='exact Coulomb 1/r')\n"
     "a1.plot(rr,phi_g(rr),'-',color='#2e8b57',lw=2,label='Gaussian-charge erf(r/0.5)/r (ours)')\n"
     "a1.plot(rU[m][::120],vU[m][::120],'o',color='#c0392b',ms=5,label='UPF pseudopotential V(r)')\n"
     "a1.plot(z_fft[z_fft<=6],phi_fft[z_fft<=6],':',color='#1b6ca8',lw=1.6,label='INQ periodic Poisson (+G=0 const)')\n"
     "a1.axhline(plateau,ls=':',color='#2e8b57',alpha=.5); a1.text(3.2,plateau+.15,f'plateau {plateau:.2f} Ha',color='#2e8b57',fontsize=8)\n"
     "a1.set_xlabel('r (Bohr)'); a1.set_ylabel('V(r) (Ha)'); a1.set_ylim(0,6); a1.set_xlim(0,6)\n"
     "a1.legend(frameon=False,fontsize=8); a1.set_title('Core: 1/r cusp vs erf-smoothed plateau')\n"
     "rl=np.logspace(np.log10(0.05),np.log10(50),400)\n"
     "a2.loglog(rl,1.0/rl,'--',color='0.55',label='1/r')\n"
     "a2.loglog(rl,phi_g(rl),'-',color='#2e8b57',lw=2,label='ours erf/r')\n"
     "a2.loglog(rU[m][::60],vU[m][::60],'o',color='#c0392b',ms=4,label='UPF V(r)')\n"
     "a2.axvline(SIGWP,ls=':',color='0.7'); a2.text(SIGWP*1.05,3e-2,'σ_WP',fontsize=8,color='0.4')\n"
     "a2.set_xlabel('r (Bohr)'); a2.set_ylabel('V(r) (Ha)'); a2.legend(frameon=False,fontsize=8)\n"
     "a2.set_title('Tail: all three collapse onto 1/r')\n"
     "fig.tight_layout(); fig.savefig('proj_potential.png',dpi=140); plt.show()")
md("**Reading it.** The pseudopotential's radial potential is **exactly** our Gaussian-charge potential — "
   "$V_\\text{loc}(r)=\\text{erf}(r/0.5)/r$ Ha pointwise (max difference ~$10^{-16}$). Both cap the core at the "
   "finite plateau $2.26$ Ha where the bare $1/r$ diverges, and both merge with Coulomb beyond ~1.5 Bohr. INQ's "
   "periodic Poisson solver reproduces this up to the standard $G{=}0$ constant (a uniform ~0.03 Ha shift of the "
   "potential zero, physically irrelevant). **Key point:** the projectile potential's *shape is correct* in both "
   "representations — it is **not** the source of the 7.4 eV pseudopotential residual. That artifact comes from "
   "how INQ *treats* this identical potential: with $Z_\\text{valence}=0$ it places the whole erf/r tail as a "
   "*truncated short-range local* potential, conditionally-convergent against the laterally-infinite background "
   "(see `reference_ghost_upf_tail_aliasing`). The perturbation feeds the same potential through the Poisson "
   "channel with no truncation, so the ledger is clean.")

# ── 5a. Case study — WP vs perturbation decomposition ────────────────────────
md("## Case study — WP vs perturbation, full energy decomposition (σ=0.5, r=12, Lz=120, p2)", anchor="case-study")
code("b = pd.read_csv('stress_baseline.csv').iloc[0]\n"
     "tab = pd.DataFrame({\n"
     "  'component (eV)': ['E_kinetic','E_hartree','E_xc','E_external','E_H+E_ext'],\n"
     "  'WP':    [b.Ekin_WP,b.EH_WP,b.Exc_WP,b.Eext_WP,b.EH_WP+b.Eext_WP],\n"
     "  'perturbation': [b.Ekin_pert,b.EH_pert,b.Exc_pert,b.Eext_pert,b.EH_pert+b.Eext_pert],\n"
     "  'WP − pert': [b.Ekin_WP-b.Ekin_pert,b.EH_WP-b.EH_pert,b.Exc_WP-b.Exc_pert,\n"
     "                b.Eext_WP-b.Eext_pert,(b.EH_WP+b.Eext_WP)-(b.EH_pert+b.Eext_pert)]})\n"
     "print(tab.to_string(index=False, float_format=lambda v:f'{v:9.2f}'))\n"
     "print(f'\\nU_proj_bg (clean ideal)      = {b.U_proj_bg:8.2f} eV')\n"
     "print(f'residual = d(E_H+E_ext) - U_proj_bg = {b.residual:7.2f} eV')\n"
     "print(f'self-Hartree: free {self_hartree_free(0.5):.2f} | in-cell {self_hartree_incell(0.5,120,2):.2f} eV'\n"
     "      f'  -> gauge = {self_hartree_free(0.5)-b.residual:.2f} eV')\n"
     "print(f'dKin={b.dKin:.2f} (loc {loc(0.5):.2f}) | dXC={b.dXC:.2f} | SIE=R+dXC={b.SIE:.2f} eV')")
md("**Reading it.** Every WP−pert difference is a *wavepacket* contribution: dKin ≈ the localisation "
   "$3/4\\sigma^2$; dXC = the WP self-XC; and $R=d(E_H+E_ext)-U_\\text{proj,bg}$ = the WP Hartree self-energy "
   "(≈ analytic $1/\\sigma_{WP}\\sqrt{2\\pi}$ minus the ~1 eV charged-cell gauge). The net H+XC self-interaction "
   "$R+dE_\\text{xc}$ is the LDA SIE — the only energy not attributable to localisation.")

# ── 5b. σ-sweep vs analytic self-Hartree ─────────────────────────────────────
md("## Stress test 1 — σ-sweep: does the residual track the analytic self-Hartree?", anchor="sigma")
code("s = pd.read_csv('stress_sigma.csv').sort_values('sigma')\n"
     "s['SH_free']=s.sigma.map(self_hartree_free)\n"
     "s['SH_incell']=[self_hartree_incell(sg,120,2) for sg in s.sigma]\n"
     "s['loc']=s.sigma.map(loc)\n"
     "print(s[['sigma','residual','SH_incell','SH_free','dKin','loc','dXC','SIE']].to_string(index=False,float_format=lambda v:f'{v:8.2f}'))\n"
     "fig,(a1,a2)=plt.subplots(1,2,figsize=(11,4.4))\n"
     "a1.plot(s.sigma,s.SH_free,'--',color='0.5',label='analytic 1/(σ√2π)')\n"
     "a1.plot(s.sigma,s.SH_incell,'-',color='#2e8b57',label='FFT in-cell (dx=0.5)')\n"
     "a1.plot(s.sigma,s.residual,'o',color='#c0392b',ms=8,label='measured residual R')\n"
     "a1.set_xlabel('σ_WP (Bohr)'); a1.set_ylabel('self-Hartree (eV)'); a1.legend(frameon=False,fontsize=8)\n"
     "a1.set_title('R tracks the self-Hartree across σ')\n"
     "a2.plot(s.sigma,s['loc'],'--',color='0.5',label='3/(4σ²)')\n"  # s.loc is the .loc indexer!
     "a2.plot(s.sigma,s.dKin,'s',color='#1b6ca8',ms=8,label='measured dKin')\n"
     "a2.set_xlabel('σ_WP (Bohr)'); a2.set_ylabel('localisation (eV)'); a2.legend(frameon=False,fontsize=8)\n"
     "a2.set_title('dKin tracks 3/(4σ²)')\n"
     "fig.tight_layout(); fig.savefig('sigma_sweep.png',dpi=140); plt.show()\n"
     "print('\\nmax |R - SH_incell| = %.2f eV' % (s.residual-s.SH_incell).abs().max())")

# ── 5c. r-independence ───────────────────────────────────────────────────────
md("## Stress test 2 — position r: the residual is a self-energy (must be flat in r)", anchor="rindep")
code("rr = pd.read_csv('stress_r.csv').sort_values('r')\n"
     "print(rr[['r','residual','dKin','dXC','SIE']].to_string(index=False,float_format=lambda v:f'{v:8.2f}'))\n"
     "fig,ax=plt.subplots(figsize=(7,4.4))\n"
     "ax.plot(rr.r,rr.residual,'o-',color='#c0392b',label='residual R')\n"
     "ax.axhline(rr.residual.mean(),ls='--',color='0.5',label=f'mean {rr.residual.mean():.2f} eV')\n"
     "ax.set_xlabel('r (Bohr)'); ax.set_ylabel('residual (eV)'); ax.legend(frameon=False,fontsize=8)\n"
     "ax.set_title(f'R flat in r (spread {rr.residual.max()-rr.residual.min():.2f} eV)')\n"
     "fig.tight_layout(); fig.savefig('r_independence.png',dpi=140); plt.show()")

# ── 5d. Lz gauge growth ──────────────────────────────────────────────────────
md("## Stress test 3 — box length Lz: the open-z charged-cell gauge", anchor="lz")
code("lz = pd.read_csv('stress_lz.csv').sort_values('Lz')\n"
     "lz['SH_incell']=[self_hartree_incell(0.5,L,2) for L in lz.Lz]\n"
     "lz['gauge_free_minus_R']=self_hartree_free(0.5)-lz.residual\n"
     "print(lz[['Lz','residual','SH_incell','gauge_free_minus_R']].to_string(index=False,float_format=lambda v:f'{v:8.2f}'))\n"
     "fig,ax=plt.subplots(figsize=(7,4.4))\n"
     "ax.plot(lz.Lz,lz.residual,'o-',color='#c0392b',label='measured residual R')\n"
     "ax.plot(lz.Lz,lz.SH_incell,'s--',color='#2e8b57',label='FFT in-cell self-Hartree')\n"
     "ax.axhline(self_hartree_free(0.5),ls=':',color='0.5',label=f'free-space {self_hartree_free(0.5):.1f} eV')\n"
     "ax.set_xlabel('Lz (Bohr)'); ax.set_ylabel('energy (eV)'); ax.legend(frameon=False,fontsize=8)\n"
     "ax.set_title('open-z gauge: R vs Lz'); fig.tight_layout(); fig.savefig('lz_gauge.png',dpi=140); plt.show()")

# ── 5e. p3 vs p2 ─────────────────────────────────────────────────────────────
md("## Stress test 4 — full PBC (p3) vs open-z (p2): the gauge convention", anchor="p3vp2")
code("p = pd.read_csv('stress_p3vp2.csv')\n"
     "print(p[['per','Ekin_WP','EH_WP','Exc_WP','Eext_WP','dEH_Eext','U_proj_bg','residual','dXC','SIE']]\n"
     "      .to_string(index=False,float_format=lambda v:f'{v:9.2f}'))\n"
     "r2=p[p.per==2].residual.iloc[0]; r3=p[p.per==3].residual.iloc[0]\n"
     "print(f'\\nresidual p2 (open-z) = {r2:.2f} eV ; p3 (full PBC) = {r3:.2f} eV ; Δ(gauge conv.) = {r3-r2:.2f} eV')\n"
     "print(f'FFT in-cell self-Hartree: p2 {self_hartree_incell(0.5,120,2):.2f} | p3 {self_hartree_incell(0.5,120,3):.2f} eV')")
md("**Reading it.** p3 (fully periodic, Makov–Payne) and p2 (open-z) place the net-charge gauge differently; "
   "the residual difference is the gauge-convention shift, not physics. A fresh WP decomposition is used for "
   "p3 (its GS/energies differ from p2). The self-Hartree itself is convention-robust to ~1 eV.")

# ── 5f. grid sweep (from proj_perturbation/grid_sweep.csv) ────────────────────
md("## Stress test 5 — grid spacing: the clean residual is grid-stable", anchor="grid")
code("import os\n"
     "gp='"+str(GRID)+"'\n"
     "if os.path.exists(gp):\n"
     "    g=pd.read_csv(gp).sort_values('dx',ascending=False)\n"
     "    print(g[['dx','HE_WP','HE_pert','U_proj_bg','residual']].to_string(index=False,float_format=lambda v:f'{v:9.2f}'))\n"
     "    fig,ax=plt.subplots(figsize=(7,4.4))\n"
     "    ax.plot(g.dx,g.residual,'o-',color='#2e8b57',ms=8,label='perturbation (clean)')\n"
     "    ax.axhline(7.36,ls='--',color='#c0392b',label='pseudopotential r_cut=120: 7.4 eV')\n"
     "    ax.set_xlabel('dx (Bohr)'); ax.set_ylabel('residual (eV)'); ax.invert_xaxis()\n"
     "    ax.legend(frameon=False,fontsize=8); ax.set_title('clean residual grid-stable near self-Hartree')\n"
     "    fig.tight_layout(); fig.savefig('grid_sweep.png',dpi=140); plt.show()\n"
     "else:\n    print('grid_sweep.csv not found yet:', gp)")

# ── 5g. Empirical boundary-matched self-Hartree ──────────────────────────────
md("## Stress test 6 — the EMPIRICAL self-Hartree: the residual is fully accounted, no analytic gauge\n"
   "The earlier sections compared $R$ to the *analytic* free-space self-Hartree $1/(\\sigma_{WP}\\sqrt{2\\pi})=21.71$ eV "
   "and absorbed the ~0.9 eV difference into a 'charged-cell gauge'. We can do better: compute the WP self-Hartree "
   "**empirically**, $E_\\text{self}=\\tfrac12\\int n_{WP}\\,\\text{poisson}(n_{WP})$, with **INQ's own Poisson solver "
   "in the actual run cell**. INQ picks the boundary-matched kernel from the cell periodicity — "
   "periodicity 3 → fully-periodic FFT (Makov–Payne), periodicity 2 → **Rozzi et al. (2006) 2D Coulomb-cutoff = "
   "open-z** (`inq/src/solvers/poisson.hpp:190`). So for the p2 production runs $E_\\text{self}$ carries the exact "
   "open-z electrostatics, no analytic correction. Source: `proj_perturbation/self_hartree.cpp` + `self_hartree_sweep.py`.")
code("emp = pd.read_csv('self_hartree_empirical.csv')\n"
     "meas = {}\n"
     "b = pd.read_csv('stress_baseline.csv').iloc[0]; meas[(2,0.5,120,0.5)] = b.residual\n"
     "for _,x in pd.read_csv('stress_p3vp2.csv').iterrows(): meas[(int(x.per),0.5,120,0.5)] = x.residual\n"
     "gp='"+str(GRID)+"'\n"
     "if os.path.exists(gp):\n"
     "    for _,x in pd.read_csv(gp).iterrows(): meas[(2,round(x.dx,2),120,0.5)] = x.residual\n"
     "for _,x in pd.read_csv('stress_lz.csv').iterrows(): meas[(2,0.5,int(x.Lz),0.5)] = x.residual\n"
     "for _,x in pd.read_csv('stress_sigma.csv').iterrows(): meas[(2,0.5,120,round(x.sigma,2))] = x.residual\n"
     "mk=lambda r:(int(r.per),round(r.dx,2),int(r.Lz),round(r.sigma,2))\n"
     "emp['R_measured']=[meas.get(mk(r),np.nan) for _,r in emp.iterrows()]\n"
     "emp['R_minus_Eself']=emp.R_measured-emp.E_self_ev\n"
     "emp['bc']=emp.per.map({2:'open-z (p2)',3:'full PBC (p3)'})\n"
     "print(emp[['bc','dx','Lz','sigma','E_self_ev','R_measured','R_minus_Eself']].to_string(index=False,float_format=lambda v:f'{v:8.3f}'))\n"
     "free=1.0/(0.5*np.sqrt(2*np.pi))*27.211386\n"
     "m=emp.dropna(subset=['R_measured'])\n"
     "print(f'\\nfree-space analytic self-Hartree = {free:.2f} eV -> vs measured p2 leaves a {free-b.residual:+.2f} eV \"gauge\" gap')\n"
     "print(f'EMPIRICAL open-z self-Hartree    = {emp[(emp.per==2)&(emp.dx==0.5)&(emp.Lz==120)&(emp.sigma==0.5)].E_self_ev.iloc[0]:.2f} eV -> vs measured p2 leaves {b.residual-emp[(emp.per==2)&(emp.dx==0.5)&(emp.Lz==120)&(emp.sigma==0.5)].E_self_ev.iloc[0]:+.2f} eV  (gap closed)')\n"
     "print(f'max |R - E_self| over σ>=0.5, both BC, grid-matched = {m[m.sigma>=0.5].R_minus_Eself.abs().max():.3f} eV')")
code("fig,(a1,a2)=plt.subplots(1,2,figsize=(11.5,4.6))\n"
     "mm=emp.dropna(subset=['R_measured'])\n"
     "for bc,c in [('open-z (p2)','#2e8b57'),('full PBC (p3)','#8e44ad')]:\n"
     "    s=mm[mm.bc==bc]; a1.scatter(s.E_self_ev,s.R_measured,s=70,color=c,label=bc,zorder=3)\n"
     "lim=[8,40]; a1.plot(lim,lim,'--',color='0.6',label='R = E_self (parity)')\n"
     "a1.set_xlim(lim); a1.set_ylim(lim); a1.set_xlabel('empirical E_self (eV)'); a1.set_ylabel('measured residual R (eV)')\n"
     "a1.legend(frameon=False,fontsize=8); a1.set_title('R vs empirical self-Hartree (all axes)')\n"
     "a1.annotate('σ=0.35 (under-resolved\\nσ_ρ=0.25 at dx=0.5)',xy=(38.1,34.9),xytext=(20,37),fontsize=7,\n"
     "            arrowprops=dict(arrowstyle='->',color='0.5'))\n"
     "# gap-closure bar at the baseline (σ=0.5, p2, dx0.5)\n"
     "labels=['free-space\\nanalytic','empirical\\nopen-z']\n"
     "e0=emp[(emp.per==2)&(emp.dx==0.5)&(emp.Lz==120)&(emp.sigma==0.5)].E_self_ev.iloc[0]\n"
     "gaps=[free-b.residual, b.residual-e0]\n"
     "a2.bar(labels,[abs(g) for g in gaps],color=['#c0392b','#2e8b57'])\n"
     "for i,g in enumerate(gaps): a2.text(i,abs(g)+0.02,f'{g:+.2f} eV',ha='center',fontsize=9)\n"
     "a2.set_ylabel('|R - self-Hartree reference| (eV)'); a2.set_ylim(0,1.05)\n"
     "a2.set_title('Gap to explain: analytic vs empirical reference')\n"
     "fig.tight_layout(); fig.savefig('empirical_self_hartree.png',dpi=140); plt.show()")
md("**Reading it.** With INQ's boundary-matched Poisson, the empirical self-Hartree reproduces the measured "
   "residual to **~0.01 eV** for every production case (σ≥0.5, both p2 and p3, grid-for-grid): open-z 20.82 (dx0.5) "
   "/ 20.65 (dx0.3), full-PBC 21.50 (dx0.5); Lz-flat at 20.82 across Lz={90..240} exactly as the measured residual "
   "is. So $R$ **is** the WP self-Hartree, with *no* analytic gauge fudge — the 0.9 eV 'gauge' was simply the "
   "open-z-vs-free-space mismatch of the wrong reference. (σ=0.35 is the lone outlier, 38.1 vs 34.9: σ_ρ=0.25 is "
   "only 0.5 grid points/σ at dx=0.5, badly under-resolved — a discretisation artifact, not physics.) "
   "**Net:** the residual is now fully explained; the only genuinely unaccounted energy is the LDA "
   "self-interaction error $R+dE_\\text{xc}\\approx4.3$ eV.")

# ── 6. Takeaway ──────────────────────────────────────────────────────────────
md("## Takeaway\n"
   "- The **perturbation (Gaussian-charge) projectile** yields a clean, r_cut-free ledger: "
   "$R=d(E_H+E_ext)-U_\\text{proj,bg}$ equals the wavepacket Hartree self-energy (≈20.8 eV open-z, σ=0.5).\n"
   "- **σ-sweep:** $R$ tracks the analytic $1/\\sigma_{WP}\\sqrt{2\\pi}$ across σ — the method measures the right "
   "physical quantity at every scale.\n"
   "- **r-independence:** $R$ is flat in projectile position — it is a genuine self-energy.\n"
   "- **Empirical closure:** with INQ's boundary-matched Poisson (Rozzi 2D for open-z), the *empirical* self-Hartree "
   "matches $R$ to ~0.01 eV in both p2 and p3, grid-for-grid — the residual is fully accounted with **no** analytic "
   "gauge correction. The p3−p2 = 0.68 eV shift is just the periodic-image (Makov–Payne) term, reproduced exactly.\n"
   "- **Grid:** $R$ is grid-stable (contrast the pseudopotential's sign-swinging 7.4 eV).\n"
   "- The only genuinely unaccounted energy is the **LDA self-interaction error** $R+dE_\\text{xc}\\approx4.3$ eV — "
   "a known functional artifact, not missing physics.")

nb = new_notebook(); nb.cells = cells
nb.metadata.kernelspec = {"name":"python3","display_name":"Python 3"}
try:
    from nbconvert.preprocessors import ExecutePreprocessor
    ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
    ep.preprocess(nb, {"metadata": {"path": str(HERE)}})
    print("executed OK")
except Exception as e:
    print("execution error (writing unexecuted):", e)
nbf.write(nb, str(OUT)); print("wrote", OUT)
