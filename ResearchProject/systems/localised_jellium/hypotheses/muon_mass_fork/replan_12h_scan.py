#!/usr/bin/env python3
"""Re-plan the effective-mass muon WP run for a <=12 h wall time.

Root cause of the old miss: the cost model was calibrated on a *different* run
(qsp_phase4, 3.75 s/step @ dx=0.50) with an n*log n scaling, predicting ~14 s/step
at dx=0.333. The live run MEASURES ~200 s/step at dx=0.333 (contended). We
re-anchor the cost model on that MEASURED point and scan coarser grids / smaller
cells / lighter masses to land <=12 h.

Locked physics anchors (unchanged):
  * v_ref = 2.711 a.u.  (= 100 eV electron velocity  =>  same S(v) comparison)
  * r_s = 5.665 localised slab (N0 = 1.312e-3), 25 Bohr thick slab
  * effective-mass fork: mass m = k0 / v_ref  (want m meaningfully > 1 for the demo)
  * T_total = 3.5 x traversal, dt from the ETRS cliff (H*dt = E_cut*dt ~ 2.2)

Relaxed this round:
  * impact spread threshold: <1%  ->  <2%
  * grid + cell + mass are FREE levers to cut cost
"""
import numpy as np

HA_EV = 27.211386
V_REF = 2.711                      # 100 eV electron velocity (LOCKED)

# ---- MEASURED cost anchor (live run, contended) -------------------------------
DX_MEAS   = 1.0/3.0                # dx of the live run
CELL_MEAS = (50.0, 50.0, 90.0)
SPS_MEAS  = 200.0                  # s/step measured (run-avg 188, recent 217) -> 200 conservative

def ngrid(dx, cell):
    return int(np.prod([np.ceil(l/dx) for l in cell]))

N_REF = ngrid(DX_MEAS, CELL_MEAS)

def sstep(dx, cell):
    """s/step, n*log n scaling anchored on the MEASURED point (contended)."""
    g = ngrid(dx, cell)
    return SPS_MEAS * (g*np.log2(g)) / (N_REF*np.log2(N_REF))

# ---- kinematics / geometry ----------------------------------------------------
sigp   = lambda sWP: 1.0/(np.sqrt(2)*sWP)          # momentum-space std
sig0   = lambda sWP: sWP/np.sqrt(2)                 # density (|psi|^2) std at t=0
SLAB_HALF = 12.5                                    # 25 Bohr slab

def config(dx, sWP, sep=2.75, cell=(50,50,90)):
    """Max out k0 to the grid ceiling for this dx; derive mass; free-space spread."""
    kmax = np.pi/dx
    k0   = kmax - 3*sigp(sWP)                        # 3-sigma_p edge at Nyquist
    m    = k0/V_REF
    E    = 0.5*m*V_REF**2 * HA_EV
    s0   = sig0(sWP)
    D_launch = sep*s0                                # launch this far before slab face
    z0   = -(SLAB_HALF + D_launch)                   # launch z
    # free/vacuum spread law: sigma_rho(t) = s0*sqrt(1+(t/tau)^2), tau = m*sWP^2 = k0*sWP^2/v
    # spread vs DISTANCE d (mass-free): sigma_rho(d)=s0*sqrt(1+(d/(2 k0 s0^2))^2)
    sr_d = lambda d: s0*np.sqrt(1+(d/(2*k0*s0**2))**2)
    impact_pct = (sr_d(D_launch)/s0 - 1)*100         # launch -> slab face
    exit_pct   = (sr_d(D_launch + 2*SLAB_HALF)/s0 - 1)*100
    # aliased fraction beyond Nyquist
    aliased_pct = (1 - 0.5*(1+np_erf((kmax-k0)/(np.sqrt(2)*sigp(sWP)))))*100
    # dt from cliff: E_cut = 0.5 kmax^2 ; safe H*dt ~ 1.8 (margin below 2.2)
    Ecut = 0.5*kmax**2
    dt_raw = 1.8/Ecut
    dt = max([x for x in (0.10,0.08,0.06,0.05,0.04,0.03,0.025,0.02,0.015,0.01) if x <= dt_raw],
             default=0.01)
    # traversal + total time: crossing = approach + slab + a little exit past face
    cross = (D_launch + 2*SLAB_HALF + sep*s0)/V_REF   # launch -> fully past slab
    T = 3.5*cross
    n = int(round(T/dt))
    s = sstep(dx, cell)
    wall_h = s*n/3600.0
    return dict(dx=round(dx,3), sWP=sWP, cell=cell, grid_M=round(ngrid(dx,cell)/1e6,2),
                k0=round(k0,2), mass=round(m,2), E_eV=round(E), z0=round(z0,1),
                impact_pct=round(impact_pct,2), exit_pct=round(exit_pct,1),
                aliased_pct=round(aliased_pct,3), Ecut=round(Ecut,1), dt=dt,
                T=round(T), n_steps=n, s_step=round(s,1),
                wall_1gpu_h=round(wall_h,1), wall_2gpu_h=round(wall_h/1.8,1))

# math.erf without importing scipy
def np_erf(x):
    # Abramowitz-Stegun 7.1.26
    t = 1.0/(1.0+0.3275911*abs(x)); y = 1 - (((((1.061405429*t-1.453152027)*t)+1.421413741)*t-0.284496736)*t+0.254829592)*t*np.exp(-x*x)
    return np.sign(x)*y

print(f"MEASURED anchor: dx={DX_MEAS:.3f} cell={CELL_MEAS} grid={N_REF/1e6:.2f}M -> {SPS_MEAS:.0f} s/step (contended)\n")
print("Old model predicted {:.1f} s/step here -> that is the ~15x miss.\n".format(
      3.754*(ngrid(DX_MEAS,(50,50,90))/ngrid(0.5,(50,50,90)))*(np.log2(ngrid(DX_MEAS,(50,50,90)))/np.log2(ngrid(0.5,(50,50,90))))))

cols = ["dx","sWP","cell","grid_M","k0","mass","E_eV","impact_pct","exit_pct",
        "aliased_pct","dt","n_steps","s_step","wall_1gpu_h","wall_2gpu_h"]
hdr = " ".join(f"{c:>10}" for c in cols); print(hdr); print("-"*len(hdr))
rows = []
for cell in [(50,50,90),(36,36,90),(36,36,80)]:
    for sWP in [1.5, 2.0]:
        for dx in [1/3, 0.40, 0.45, 0.50]:
            r = config(dx, sWP, cell=cell); rows.append(r)
            print(" ".join(f"{str(r[c]):>10}" for c in cols))
    print()

print("\nFILTER: impact<2% AND mass>=1.8 AND wall_1gpu_h<=12  (sorted by wall):")
ok = sorted([r for r in rows if r['impact_pct']<2 and r['mass']>=1.8 and r['wall_1gpu_h']<=12],
            key=lambda r: r['wall_1gpu_h'])
for r in ok:
    print(f"  dx={r['dx']} sWP={r['sWP']} cell={r['cell']} m={r['mass']} E={r['E_eV']}eV "
          f"impact={r['impact_pct']}% dt={r['dt']} n={r['n_steps']} "
          f"-> {r['wall_1gpu_h']}h/1gpu ({r['wall_2gpu_h']}h/2gpu)")
