#!/usr/bin/env python3
"""Pre-run estimate: free-Gaussian broadening of a momentum-matched muon WP.

Design question (2026-07-07): to keep the SAME grid as the electron localised
jellium (r_s=5.69 slab), match the INITIAL MOMENTUM of the muon WP to a 300 eV
electron WP (same k0 -> same de Broglie wavelength -> same grid). This script:
  * computes the muon energy/velocity at matched momentum,
  * plots the free (vacuum) density-width sigma_rho(t) with the slab
    enter / exit / CAP timestamps annotated,
  * shows the mass-independence of spreading-per-DISTANCE at fixed k0.

Free-particle law (matches the Phase-2 oracle, sigma_z2(0)=0.125 for sigma_WP=0.5):
    sigma_rho(t)^2 = sigma_rho0^2 + t^2 / (4 m^2 sigma_rho0^2),   sigma_rho0 = sigma_WP/sqrt(2)
Recast vs distance d = v t (with m v = k0):
    sigma_rho(d)^2 = sigma_rho0^2 + d^2 / (4 k0^2 sigma_rho0^2)   <- NO mass.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from inqview.visualisation.style import apply_theme
    apply_theme()
except Exception as e:
    print(f"[warn] canonical theme unavailable ({e}); using fallback rcParams")
    plt.rcParams.update({"figure.dpi": 140, "font.size": 10, "axes.grid": True,
                         "grid.alpha": 0.3, "axes.sppath" if False else "axes.axisbelow": True})

HA_EV   = 27.211386
M_MU    = 206.77          # muon / electron mass ratio
SIG_WP  = 0.5             # wavefunction std (INJECTED width); density std = SIG_WP/sqrt2
SIG0    = SIG_WP / np.sqrt(2.0)   # density position std, 0.35355
SIG0_2  = SIG0**2                 # 0.125  (== Phase-2 sigma_z2(0))

# --- momentum matching to a 300 eV electron ---------------------------------
KE_E_EV = 300.0
KE_E_HA = KE_E_EV / HA_EV
K0      = np.sqrt(2.0 * 1.0 * KE_E_HA)     # electron k0 = muon k0 (matched)
V_E     = K0 / 1.0
V_MU    = K0 / M_MU
KE_MU_HA = K0**2 / (2.0 * M_MU)
KE_MU_EV = KE_MU_HA * HA_EV

# --- slab geometry (reuse slab_n82_L50x50x90: launch -23.75, faces +/-12.5, CAP +/-35) ---
D_ENTER, D_EXIT, D_CAP = 11.25, 36.25, 58.75     # distances from launch z0=-23.75
t_enter, t_exit, t_cap = D_ENTER/V_MU, D_EXIT/V_MU, D_CAP/V_MU

def sigma_rho_t(t, mass):
    return np.sqrt(SIG0_2 + t**2 / (4.0 * mass**2 * SIG0_2))

def sigma_rho_d(d):                     # mass-independent at fixed k0
    return np.sqrt(SIG0_2 + d**2 / (4.0 * K0**2 * SIG0_2))

print(f"matched k0        = {K0:.4f} Bohr^-1")
print(f"electron: E=300 eV v={V_E:.4f} a.u.")
print(f"muon    : E={KE_MU_EV:.4f} eV  v={V_MU:.5f} a.u.")
print(f"t_enter={t_enter:.0f}  t_exit={t_exit:.0f}  t_cap={t_cap:.0f} a.u.")
for label, t in [("enter", t_enter), ("exit", t_exit), ("CAP", t_cap)]:
    print(f"  sigma_rho(t_{label}) = {sigma_rho_t(t, M_MU):.2f} Bohr")

# --- figure -----------------------------------------------------------------
fig, (axA, axB) = plt.subplots(1, 2, figsize=(10.5, 4.2))

# Panel A: sigma_rho(t) for the muon, annotated
t = np.linspace(0.0, t_cap * 1.03, 2000)
axA.plot(t, sigma_rho_t(t, M_MU), color="#1f4e79", lw=2.2,
         label=r"muon WP, $E_\mu$=1.45 eV")
for label, tt, col in [("enter slab", t_enter, "#2c7a2c"),
                       ("exit slab",  t_exit,  "#b5651d"),
                       ("reach CAP",  t_cap,   "#8b1a1a")]:
    s = sigma_rho_t(tt, M_MU)
    axA.axvline(tt, color=col, ls="--", lw=1.1, alpha=0.8)
    axA.plot([tt], [s], "o", color=col, ms=5, zorder=5)
    axA.annotate(f"{label}\n"+r"$\sigma_\rho$="+f"{s:.1f} $a_0$\nt={tt:.0f} a.u.",
                 xy=(tt, s), xytext=(tt-160, s+2.6), fontsize=8, color=col,
                 ha="right", va="bottom")
axA.axhline(SIG0, color="grey", ls=":", lw=1)
axA.text(40, SIG0+0.2, r"$\sigma_{\rho,0}$="+f"{SIG0:.2f} $a_0$", fontsize=8, color="grey")
axA.set_xlabel("time  t  [a.u.]")
axA.set_ylabel(r"density width  $\sigma_\rho(t)$  [Bohr]")
axA.set_title("Free (vacuum) muon-WP broadening vs time", fontsize=10)
axA.set_xlim(0, t_cap*1.03); axA.set_ylim(0, sigma_rho_t(t_cap, M_MU)*1.12)
axA.legend(loc="upper left", fontsize=8, frameon=False)

# Panel B: sigma_rho vs DISTANCE -- electron and muon coincide
d = np.linspace(0.0, D_CAP*1.03, 800)
axB.plot(d, sigma_rho_d(d), color="#5a3d7a", lw=2.4, label="muon  &  electron (coincide)")
for label, dd, col in [("slab entry", D_ENTER, "#2c7a2c"),
                       ("slab exit",  D_EXIT,  "#b5651d"),
                       ("CAP",        D_CAP,   "#8b1a1a")]:
    axB.axvline(dd, color=col, ls="--", lw=1.1, alpha=0.8)
    axB.plot([dd], [sigma_rho_d(dd)], "o", color=col, ms=5, zorder=5)
    axB.text(dd-1.0, sigma_rho_d(dd)+0.35, label, fontsize=8, color=col, ha="right")
axB.set_xlabel("distance travelled  d = v t  [Bohr]")
axB.set_ylabel(r"$\sigma_\rho$  [Bohr]")
axB.set_title(r"At matched $k_0$, spread-per-distance is mass-independent", fontsize=10)
axB.set_xlim(0, D_CAP*1.03); axB.set_ylim(0, sigma_rho_d(D_CAP)*1.12)
axB.legend(loc="upper left", fontsize=8, frameon=False)

fig.suptitle(r"Momentum-matched muon WP  ($k_0$=4.696 $a_0^{-1}$; $E_\mu$=1.45 eV $\leftrightarrow$ 300 eV electron), $\sigma_{WP}$=0.5",
             fontsize=10.5)
fig.tight_layout(rect=(0, 0, 1, 0.96))
out = __file__.rsplit("/", 1)[0] + "/muon_wp_broadening.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"wrote {out}")
