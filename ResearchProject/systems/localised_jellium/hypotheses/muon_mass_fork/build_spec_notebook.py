#!/usr/bin/env python3
"""Builder for muon_effmass_spec.ipynb — Phase-4 (momentum-matched effective-mass
projectile) specification & feasibility.

v3 (2026-07-07): grid budget relaxed to 1.7× (dx≥0.294); scan refocused on
CONCENTRATED packets with <1% spread up to impact.
Run:  python build_spec_notebook.py   (auto-executes via nbconvert afterwards)
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT  = HERE / "muon_effmass_spec.ipynb"
nb = new_notebook(); C = []

C.append(new_markdown_cell(r"""
# Phase 4 — Momentum-matched effective-mass projectile on the $r_s{=}5.69$ localised jellium

**Specification & feasibility.** The `inq-study` mass fork (Phases 1–3 validated) adds a *third*
projectile parameter — **mass** — to (energy, width). We match the projectile's **initial
momentum** to the "good" **100 eV electron** localised-jellium run so it lives on a near-standard
grid, then use mass to bring velocity/energy to a comparable regime while suppressing spreading.

**This version:** grid budget relaxed to **≤1.7× finer** (dx≥0.294), and the scan targets the most
**concentrated** packet that still spreads **<1% up to impact**.

**Sources / provenance**
- Reference geometry: `shared/configs/slab_n82_L50x50x90.hpp` (slab 25 Bohr, $N{=}82$,
  $n_0{=}1.312\times10^{-3}$, $r_s{=}5.665$, dx=0.50, $\sigma_{WP}{=}0.5$, 100 eV).
- Spreading oracle validated Phase 2 (`sigma_z2(0)=0.125` for $\sigma_{WP}{=}0.5$).
- GPU calibration: measured `qsp_phase4` (3.75 s/step, dx=0.50, 61 states).
- Concept note: `docs/misc/thoughts/mass-as-a-knob-in-simulations.md`.
""".strip()))

C.append(new_markdown_cell(r"""
## 1. Physics & formulas (every term defined)

Density-width spreading ($\sigma_\rho=$ std of $|\psi|^2$, $\sigma_{\rho,0}=\sigma_{WP}/\sqrt2$):
$$\sigma_\rho(t)=\sigma_{\rho,0}\sqrt{1+(t/\tau_s)^2},\qquad \tau_s=m\,\sigma_{WP}^2 .$$
Kinematics for mean wavevector $k_0$:  $v=k_0/m$, $E=k_0^2/2m=\tfrac12 mv^2$, $\sigma_p=1/(\sqrt2\sigma_{WP})$.
**Grid ceiling (≤1.7× refinement of dx=0.50):** $k_0+3\sigma_p\le k_{max}=\pi/0.294=10.68\ a_0^{-1}$.

**Two spreading measures.** *Impact* spread (launch→slab face) $\propto D_{launch}/(k_0\sigma_{WP}^2)$
(mass-free) — beaten <1% by widening $\sigma_{WP}$, shrinking $D_{launch}$, or the larger $k_0$ the
1.7× grid now allows. *Through-slab* growth over the fixed 25 Bohr $\propto 25/(k_0\sigma_{WP}^2)$ — also
momentum-set; its floor is a few % (sub-1% would need $\sigma_{WP}{>}5$, wider than the slab). So **<1%
is the *impact* target**; through-slab we minimise. The extra k-space from 1.7× lets the *concentrated*
end (smaller $\sigma_{WP}$) clear <1% at impact — the whole point of relaxing the grid.
""".strip()))

C.append(new_code_cell(r"""
%matplotlib inline
import numpy as np, pandas as pd, math
import matplotlib.pyplot as plt
try:
    from inqview.visualisation.style import apply_theme; apply_theme()
except Exception as e:
    print(f"[theme fallback] {e}")
    plt.rcParams.update({"figure.dpi":120,"axes.grid":True,"grid.alpha":.3,"font.size":10})

HA_EV=27.211386; DX_NOW=0.50; REFINE_MAX=1.7
DX_MIN=DX_NOW/REFINE_MAX; KMAX=np.pi/DX_MIN
V_REF=np.sqrt(2*100/HA_EV); SLAB_HALF, LZ = 12.5, 90.0
sig0=lambda s: s/np.sqrt(2)
sigp=lambda s: 1/(np.sqrt(2)*s)
tail=lambda n: 0.5*math.erfc(n/np.sqrt(2))*100      # % WP density beyond n*sigma (one tail)

def config(sWP, sep):
    "Grid-max-k0, velocity-matched config; launch sep*sigma_rho0 before the slab face."
    k0=KMAX-3*sigp(sWP); m=k0/V_REF; E=0.5*m*V_REF**2*HA_EV
    dx=np.pi/(k0+3*sigp(sWP)); s0=sig0(sWP); Dl=sep*s0
    z0=-SLAB_HALF-Dl; Dcell=LZ/2-z0; tau_tr=Dcell/V_REF
    sr=lambda d: s0*np.sqrt(1+(d/(2*k0*s0**2))**2)
    imp=(sr(Dl)/s0-1)*100; grow=(sr(Dl+2*SLAB_HALF)/sr(Dl)-1)*100
    return dict(sigma_WP=sWP, sep=sep, k0=k0, mass=m, E_eV=E, dx=dx, sig0=s0, D_launch=Dl,
                z0=z0, tau_trav=tau_tr, T=3.5*tau_tr, impact=imp, grow=grow,
                overlap=tail(sep), w_impact=sr(Dl), w_exit=sr(Dl+2*SLAB_HALF))
print(f"grid budget {REFINE_MAX}x -> dx_min={DX_MIN:.3f}  k_max={KMAX:.3f}  v_ref={V_REF:.4f}")
""".strip()))

C.append(new_markdown_cell("## 2. Scan over $\\sigma_{WP}$ × launch separation (impact < 1%, overlap < 1%)"))

C.append(new_code_cell(r"""
rows=[]
for sWP in [1.0,1.25,1.5,1.6,1.75,2.0]:
    for sep in [2.0,2.5,2.75,3.0]:
        c=config(sWP,sep)
        rows.append(dict(sigma_WP=sWP, launch_sep=f"{sep:.2f}σ", mass=round(c['mass'],2),
                         E_eV=round(c['E_eV']), dx=round(c['dx'],3),
                         D_launch=round(c['D_launch'],2), overlap_pct=round(c['overlap'],2),
                         impact_pct=round(c['impact'],2), slab_growth_pct=round(c['grow'],0),
                         ok=("✓" if (c['impact']<1 and c['overlap']<1) else "")))
scan=pd.DataFrame(rows)
scan.style.apply(lambda r:['background-color:#d6efd6' if r['ok']=="✓" else '' for _ in r],axis=1)
""".strip()))

C.append(new_markdown_cell(r"""
**Reading the scan.** With the 1.7× budget, a *clean* <1% impact (overlap <1%, launch ≥ ~2.7σ) now
reaches down to $\sigma_{WP}\approx1.5$ — vs $\approx2.0$ at 1.5×. So the extra grid buys a genuinely
more concentrated packet ($\sigma_{\rho,0}$ from 1.4 → 1.06 Bohr). Going below $\sigma_{WP}{\approx}1.25$
still forces a launch inside ~2σ (the WP tail overlaps the slab at $t{=}0$).
""".strip()))

C.append(new_markdown_cell("## 3. Contenders (context) — **B selected**"))

C.append(new_code_cell(r"""
TOP = {
 "A · tight":         config(1.25, 2.00),
 "B · concentrated ★":config(1.50, 2.75),   # <-- SELECTED
 "C · clean":         config(1.75, 3.00),
 "D · balanced":      config(2.00, 3.00),
}
tab=pd.DataFrame({k:{
    "σ_WP (Bohr)":f"{v['sigma_WP']:.2f}", "σ_ρ,0 (Bohr)":f"{v['sig0']:.2f}",
    "mass":f"{v['mass']:.2f}", "E (eV)":f"{v['E_eV']:.0f}", "k0":f"{v['k0']:.2f}",
    "launch":f"{v['sep']:.2f}σ", "overlap%":f"{v['overlap']:.2f}",
    "impact%":f"{v['impact']:.2f}", "slab-growth%":f"{v['grow']:.0f}",
} for k,v in TOP.items()})
tab
""".strip()))

C.append(new_markdown_cell(r"""
## 4. Deep dive — the selected **B** packet ($\sigma_{WP}=1.5,\ m=3.4,\ 342$ eV)

Geometry from `slab_n82_L50x50x90`: box $z\in[-45,45]$, slab $[-12.5,12.5]$, two-sided CAP
$[\pm35,\pm45]$. B launches at $z_0=-15.4$ (2.75σ before the slab face). Below: the widths ledger,
then four diagnostic views.
""".strip()))

C.append(new_code_cell(r"""
B = config(1.50, 2.75)
s0,k0,m,Dl,z0 = B['sig0'],B['k0'],B['mass'],B['D_launch'],B['z0']
tau_s = m*B['sigma_WP']**2
CAP_IN = 35.0                                          # CAP inner face (reference config)
def z_of_t(t): return z0 + V_REF*t
def sr_t(t):   return s0*np.sqrt(1+(t/tau_s)**2)
# key events (by z crossing)
events = [("launch",   z0),   ("impact (slab face)", -SLAB_HALF),
          ("slab centre", 0.0),("exit slab", SLAB_HALF), ("reach CAP", CAP_IN)]
led=[]
for lab,z in events:
    t=(z-z0)/V_REF; led.append(dict(event=lab, z=round(z,1), t_au=round(t,2),
        sigma_rho=round(sr_t(t),2), spread_vs_init_pct=round((sr_t(t)/s0-1)*100,1)))
ledger=pd.DataFrame(led); print("through-slab growth (impact→exit): "
      f"{(sr_t((SLAB_HALF-z0)/V_REF)/sr_t((-SLAB_HALF-z0)/V_REF)-1)*100:.0f}%"); ledger
""".strip()))

C.append(new_code_cell(r"""
fig,axs=plt.subplots(2,2,figsize=(11.5,8.4)); (a,b),(c,d)=axs
tcap=(CAP_IN-z0)/V_REF                                  # packet absorbed ~here
ev_c={"launch":"#555","impact (slab face)":"#2c7a2c","slab centre":"#1f4e79",
      "exit slab":"#b5651d","reach CAP":"#8b1a1a"}

# (a) sigma_rho(t) annotated
t=np.linspace(0,tcap*1.02,500); a.plot(t,sr_t(t),color="#1f4e79",lw=2.4)
for lab,z in events:
    tt=(z-z0)/V_REF; a.axvline(tt,color=ev_c[lab],ls="--",lw=1,alpha=.7)
    a.plot([tt],[sr_t(tt)],'o',color=ev_c[lab],ms=5)
    a.annotate(f"{lab}\n{sr_t(tt):.2f} $a_0$",(tt,sr_t(tt)),(tt+0.3,sr_t(tt)-0.25),fontsize=7.5,color=ev_c[lab])
a.set_xlabel("time [a.u.]"); a.set_ylabel(r"$\sigma_\rho(t)$ [Bohr]")
a.set_title("(a) width vs time (free/vacuum estimate)")

# (b) real-space |psi|^2(z) snapshots vs slab
zz=np.linspace(-20,20,800)
for lab,z in events[:4]:
    tt=(z-z0)/V_REF; sg=sr_t(tt); zc=z_of_t(tt)
    b.plot(zz, np.exp(-(zz-zc)**2/(2*sg**2))/(np.sqrt(2*np.pi)*sg), color=ev_c[lab], lw=2, label=f"{lab} (t={tt:.1f})")
b.axvspan(-SLAB_HALF,SLAB_HALF,color="grey",alpha=.15); b.text(0,b.get_ylim()[1]*0.9,"slab (25 Bohr)",ha="center",fontsize=8,color="grey")
b.set_xlabel("z [Bohr]"); b.set_ylabel(r"$|\psi(z)|^2$ (norm.)"); b.set_title("(b) packet vs slab — real space")
b.legend(fontsize=7,frameon=False)

# (c) momentum space vs grid Nyquist
kk=np.linspace(0,11.6,600); sp=sigp(B['sigma_WP'])
c.plot(kk, np.exp(-(kk-k0)**2/(2*sp**2)), color="#1f4e79", lw=2.4)
c.axvline(k0,color="#1f4e79",ls=":",lw=1); c.text(k0,1.02,f"$k_0$={k0:.2f}",ha="center",fontsize=8,color="#1f4e79")
c.axvspan(KMAX,11.6,color="#8b1a1a",alpha=.15); c.axvline(KMAX,color="#8b1a1a",lw=1.5)
c.text(KMAX,0.5,f" $k_{{max}}$={KMAX:.2f}\n (aliasing →)",fontsize=8,color="#8b1a1a")
c.annotate("",(k0+3*sp,0.15),(k0,0.15),arrowprops=dict(arrowstyle="<->",color="grey"))
c.text(k0+1.5*sp,0.18,f"$3\\sigma_p$={3*sp:.2f}",fontsize=7.5,color="grey",ha="center")
c.set_xlabel("k [1/Bohr]"); c.set_ylabel(r"$|\psi(k)|^2$ (norm.)"); c.set_title("(c) momentum vs Nyquist (0.13% aliased)")

# (d) run geometry schematic
d.set_xlim(-46,46); d.set_ylim(0,1); d.set_yticks([])
d.axvspan(-SLAB_HALF,SLAB_HALF,0.25,0.75,color="#8899bb",alpha=.6)
d.text(0,0.85,"jellium slab\nr_s=5.665, N=82",ha="center",fontsize=8)
for lo,hi in [(-45,-35),(35,45)]:
    d.axvspan(lo,hi,0.25,0.75,facecolor="none",edgecolor="#8b1a1a",hatch="///",lw=1)
d.text(40,0.85,"CAP",ha="center",fontsize=8,color="#8b1a1a")
zc=np.linspace(z0-4,z0+4,100); d.plot(zc,0.5+0.18*np.exp(-(zc-z0)**2/(2*s0**2)),color="#1f4e79",lw=2)
d.annotate("",(z0+7,0.5),(z0+2.5,0.5),arrowprops=dict(arrowstyle="->",color="#1f4e79",lw=2))
d.text(z0,0.28,f"launch $z_0$={z0:.1f}\n$k_0$→ +z",ha="center",fontsize=8,color="#1f4e79")
d.set_xlabel("z [Bohr]"); d.set_title("(d) run geometry (50×50×90 box)")
fig.tight_layout(); plt.show()
""".strip()))

C.append(new_code_cell(r"""
# (e) the win: B vs the old 100 eV electron, over the packet's lifetime in the box
fig,ax=plt.subplots(figsize=(7.2,4.2)); t=np.linspace(0,tcap*1.02,500)
ax.plot(t, sig0(0.5)*np.sqrt(1+(t/(1.0*0.5**2))**2), lw=2.4, color="#8b1a1a", label="old 100 eV e⁻ (σ=0.5, m=1)")
ax.plot(t, sr_t(t), lw=2.4, color="#1f4e79", label="B (σ=1.5, m=3.4, 342 eV)")
ax.axvspan((-SLAB_HALF-z0)/V_REF,(SLAB_HALF-z0)/V_REF,color="grey",alpha=.13)
ax.text(((-SLAB_HALF-z0)/V_REF+(SLAB_HALF-z0)/V_REF)/2, ax.get_ylim()[1]*0.9,"inside slab",ha="center",fontsize=8,color="grey")
ax.set_xlabel("time [a.u.]"); ax.set_ylabel(r"$\sigma_\rho(t)$ [Bohr]"); ax.set_ylim(0,8)
ax.set_title("(e) B stays a defined projectile; the old electron does not")
ax.legend(fontsize=9,frameon=False); fig.tight_layout(); plt.show()
""".strip()))

C.append(new_markdown_cell(r"""
**What the panels say.** (a) B spreads 0.97% by impact, reaching the slab at 1.07 Bohr and leaving it
at 1.77 Bohr. (b) even at exit the packet ($\approx$1.8 Bohr) is far narrower than the 25 Bohr slab — a
well-defined projectile throughout. (c) the momentum packet sits with its $3\sigma_p$ edge exactly at the
grid Nyquist (0.13% aliased) — the 1.7× grid is *used*, not wasted. (d) launch/slab/CAP geometry. (e) B
vs the old electron, which balloons past 8 Bohr before it even reaches the slab. *(a,b,e are the free
vacuum estimate; the medium modifies the true widths, measured in the run.)*
""".strip()))

C.append(new_markdown_cell(r"""
## 5. Time step & stability — **MEASURED** (smoke tests, 2026-07-07)

dt stability is set by the **mass-1 bath electrons** at the plane-wave cutoff $E_{cut}=\tfrac12 k_{max}^2$
(the heavy projectile's max is ~6× smaller — the *bath* limits dt). Rather than trust the ETRS Taylor
estimate, we ran the built regression probe (He, LDA, **ETRS** — same propagator) at the production
cutoff and swept dt. **The instability cliff is $H\,dt=E_{cut}\,dt\approx2.2$, consistent across grids.**
For B's 1.7× grid ($E_{cut}=57$ Ha), **dt=0.08 and 0.04 both diverge (NaN); dt=0.02 is safe.** dt=0.03
survives but sits on the cliff (drift ~5× larger) — too risky for a ~2600-step run.
""".strip()))

C.append(new_code_cell(r"""
# Measured smoke-test results (regression He/LDA/ETRS probe at each cutoff), 2026-07-07.
smoke = pd.DataFrame([
 ("1.7× (B)", 0.294, 57.0, 0.08, "NaN @ step 3"),
 ("1.7× (B)", 0.294, 57.0, 0.04, "NaN @ step 9"),
 ("1.7× (B)", 0.294, 57.0, 0.03, "stable (edge, drift 8e-3)"),
 ("1.7× (B)", 0.294, 57.0, 0.025,"stable (drift 4e-3)"),
 ("1.7× (B)", 0.294, 57.0, 0.02, "STABLE ✓ (drift 1.5e-3)"),
 ("1.5× (D)", 0.333, 44.4, 0.05, "NaN @ step 9"),
 ("1.5× (D)", 0.333, 44.4, 0.04, "STABLE ✓ (drift 7e-3)"),
 ("coarse ctl",0.57, 19.7, 0.08, "stable (grid-dependence check)"),
], columns=["grid","dx","E_cut_Ha","dt","verdict"])
smoke["H_dt"]=(smoke.E_cut_Ha*smoke.dt).round(2)
smoke[["grid","dx","E_cut_Ha","dt","H_dt","verdict"]]
""".strip()))

C.append(new_markdown_cell("## 6. Expected GPU wall time"))

C.append(new_code_cell(r"""
def ngrid(dx,L=(50,50,90)): return int(np.prod([np.ceil(l/dx) for l in L]))
def sstep(dx): g=ngrid(dx); return 3.754*(g/ngrid(0.50))*(np.log2(g)/np.log2(ngrid(0.50)))  # meas. calib
T=TOP["B · concentrated ★"]['T']
rows=[]
for name,dx,dt in [("B · concentrated (1.7×)",0.294,0.02),("B · concentrated (1.7×)",0.294,0.025),
                   ("D · balanced (1.5×)",0.333,0.04)]:
    s=sstep(dx); n=int(round(T/dt))
    rows.append(dict(config=name, dx=dx, grid_M=round(ngrid(dx)/1e6,2), dt=dt, n_steps=n,
                     s_per_step=round(s,1), wall_h_1gpu=round(s*n/3600,1), wall_h_2gpu=round(s*n/3600/1.8,1)))
pd.DataFrame(rows)
""".strip()))

C.append(new_markdown_cell(r"""
The concentrated packet (B) is expensive on **two** fronts: the 1.7× grid has ~1.5× more points *and*
forces **dt=0.02** (measured). Net **~22 h on one GPU (~12 h on two)**. The coarser 1.5× grid of D allows
**dt=0.04**, giving **~7 h** — ~3× faster — for a slightly wider launch packet. `s/step` calibrated to the
measured `qsp_phase4` run; GS convergence on the fine grid is the remaining pre-run gate.
""".strip()))

C.append(new_markdown_cell("## 7. Suggested run specification (recommend **B · concentrated**)"))

C.append(new_code_cell(r"""
b=TOP["B · concentrated ★"]
spec=pd.DataFrame([
 ("Projectile width σ_WP", f"{b['sigma_WP']:.2f} Bohr", "concentrated & clean <1% impact"),
 ("Effective mass m",      f"{b['mass']:.1f} m_e",   "tuned via inq-study mass fork"),
 ("Velocity v",            f"{V_REF:.2f} a.u.",      "= 100 eV electron ⇒ same S(v)"),
 ("Energy E",              f"{b['E_eV']:.0f} eV",    "= ½ m v²"),
 ("Momentum k0",           f"{b['k0']:.2f} 1/Bohr",  "sets grid"),
 ("Grid dx",               f"{b['dx']:.3f} Bohr",    f"{DX_NOW/b['dx']:.2f}× finer (≤1.7×)"),
 ("Cell / slab",           "50×50×90 Bohr / 25 Bohr", "reuse slab_n82_L50x50x90"),
 ("Density",               "r_s=5.665, N=82 (even)", "matches reference exactly"),
 ("Launch z0",             f"{b['z0']:.1f} Bohr",    f"{b['sep']:.2f}σ before slab ({b['overlap']:.2f}% overlap)"),
 ("Spread @ impact",       f"{b['impact']:.2f} %",   "target <1% ✓"),
 ("Growth through slab",   f"{b['grow']:.0f} %",     f"width {b['sig0']:.2f}→{b['w_exit']:.2f} Bohr"),
 ("Time step dt",          "0.02 a.u. (MEASURED)",   "smoke test: 0.08 & 0.04 diverge on 1.7× grid"),
 ("Total time T",          f"{b['T']:.0f} a.u.",     "= 3.5× traversal (<80)"),
 ("Theory",                "LDA, ETRS, γ-only",      "free Ehrenfest; bath = mass-1 electrons"),
 ("Est. GPU wall",         "~22 h (1 GPU, dt=0.02)", "~12 h on 2 GPUs; measured s/step calib"),
 ("Stopping readout",      "n(k,t) coherent peak",   "robust to residual real-space spread"),
], columns=["parameter","value","note"])
spec
""".strip()))

C.append(new_markdown_cell(r"""
## 8. Takeaway

Relaxing the grid to **1.7×** buys the concentrated packet we wanted: **B — $\sigma_{WP}{=}1.5$,
$m{=}3.4$, $v{=}2.71$ (=100 eV electron), $E{\approx}342$ eV, dx=0.294** — launches at $\sigma_{\rho,0}{=}1.06$
Bohr (vs 1.41 at 1.5×), stays **<1% up to impact (0.97%)** with a clean 2.75σ launch (0.30% $t{=}0$
overlap), at the same $r_s{=}5.665$ density. The cost is a heavier grid (~9 M pts ⇒ ~20–24 h/GPU) and
more through-slab growth than the wider options. Tighter still (A, $\sigma_{WP}{=}1.25$) is available if
you accept ~2% launch overlap. Through-slab <1% remains out of reach for a localized projectile — read
stopping from the $n(k,t)$ coherent peak.

**Open:** confirm B (or A for max concentration); dt 0.02 vs 0.04; run the 300 eV electron twin on the
same dx for a one-grid/two-mass comparison?
""".strip()))

nb["cells"]=C; nbf.write(nb, OUT); print("wrote", OUT)
