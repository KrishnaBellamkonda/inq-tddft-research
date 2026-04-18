#!/usr/bin/env python3
"""
plot_shell_structure.py — Free-electron shell structure for the jellium ground state.

Produces a side-by-side level diagram:
  Left  panel: free-electron energies  ε_k = k²/2  (kinetic only)
  Right panel: KS eigenvalues          ε_k = k²/2 + V_xc(n₀)

Shells are colour-coded by occupancy and labelled with |n|², electron count,
and degeneracy. The constant V_xc shift (LDA, Perdew-Zunger 1981) is shown as
a bracketed annotation between the panels.

Usage:
    python3 plot_shell_structure.py

Output:
    shell_structure.pdf, shell_structure.png

Reference for V_xc: Perdew & Zunger, PRB 23, 5048 (1981).
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

plt.rcParams.update({
    'font.family':      'serif',
    'font.size':        10,
    'axes.linewidth':   0.8,
    'xtick.direction':  'in',
    'ytick.direction':  'in',
})

# ── System parameters ─────────────────────────────────────────────────────────

N  = 40
L  = 40.0           # cell side (bohr) — updated to match run.cpp
HA_TO_EV = 27.211386245988

n0  = N / L**3                        # mean electron density (bohr⁻³)
rs  = (3.0 / (4.0 * np.pi * n0))**(1.0/3.0)   # Wigner-Seitz radius (bohr)
kF  = (3.0 * np.pi**2 * n0)**(1.0/3.0)         # Fermi wavevector (bohr⁻¹)
EF  = 0.5 * kF**2                    # Fermi energy (Ha)
k0  = 2.0 * np.pi / L               # primitive reciprocal vector (bohr⁻¹)

# ── PZ81 LDA exchange-correlation potential ───────────────────────────────────
# Perdew & Zunger, PRB 23, 5048 (1981), metallic regime (r_s ≥ 1).

def vxc_pz81(rs):
    ex  = -0.4582 / rs
    vx  = (4.0/3.0) * ex
    gamma, beta1, beta2 = -0.1423, 1.0529, 0.3334
    sqrtrs = np.sqrt(rs)
    denom  = 1.0 + beta1*sqrtrs + beta2*rs
    ec     = gamma / denom
    vc     = ec * (1.0 + (7.0/6.0)*beta1*sqrtrs + (4.0/3.0)*beta2*rs) / denom
    return vx + vc

Vxc_Ha = vxc_pz81(rs)   # in Hartree (negative)
Vxc_eV = Vxc_Ha * HA_TO_EV

# ── Shell structure ───────────────────────────────────────────────────────────
# Enumerate all k = (2π/L)n for n ∈ ℤ³ up to |n|²_max.

N2_MAX = 6
nmax   = int(np.sqrt(N2_MAX)) + 1
degen  = {}
for nx in range(-nmax, nmax+1):
    for ny in range(-nmax, nmax+1):
        for nz in range(-nmax, nmax+1):
            n2 = nx**2 + ny**2 + nz**2
            if n2 <= N2_MAX:
                degen[n2] = degen.get(n2, 0) + 1

# Build shell list, compute occupancy by filling from bottom up.
shells = []
remaining = N
for n2 in sorted(degen.keys()):
    deg  = degen[n2]
    cap  = 2 * deg          # electrons this shell can hold (×2 spin)
    fill = min(remaining, cap)
    remaining -= fill
    shells.append({
        'n2':          n2,
        'deg':         deg,
        'Ek_Ha':       0.5 * n2 * k0**2,
        'electrons':   fill,
        'capacity':    cap,
        'frac':        fill / cap,
    })

# ── Layout ────────────────────────────────────────────────────────────────────

fig, (ax_left, ax_right) = plt.subplots(
    1, 2, figsize=(9, 6), sharey=True,
    gridspec_kw={'width_ratios': [1, 1], 'wspace': 0.04}
)

# Colour scheme
COL_FULL    = '#2166ac'   # blue  — fully occupied
COL_PARTIAL = '#e6711f'   # amber — partially occupied
COL_EMPTY   = '#aaaaaa'   # grey  — unoccupied

LWIDTH_SCALE = 12         # linewidth (pt) for the most degenerate shell

def level_lw(deg):
    """Line width proportional to shell degeneracy (visual density of states)."""
    max_deg = max(sh['deg'] for sh in shells)
    return max(2.5, LWIDTH_SCALE * deg / max_deg)

def level_color(frac):
    if frac >= 0.999:   return COL_FULL,    0.90
    elif frac > 0.001:  return COL_PARTIAL, 0.90
    else:               return COL_EMPTY,   0.55

# ── Draw both panels ──────────────────────────────────────────────────────────

for ax, offset_Ha, title in [
    (ax_left,  0.0,      r'Free-electron energies   $\varepsilon_k = k^2/2$'),
    (ax_right, Vxc_Ha,   r'KS eigenvalues   $\varepsilon_k = k^2/2 + V_{xc}$'),
]:
    for sh in shells:
        E_eV = (sh['Ek_Ha'] + offset_Ha) * HA_TO_EV
        col, alpha = level_color(sh['frac'])
        lw = level_lw(sh['deg'])

        # Horizontal energy level (x from 0.15 to 0.85)
        ax.plot([0.15, 0.85], [E_eV, E_eV],
                color=col, linewidth=lw, alpha=alpha,
                solid_capstyle='round', zorder=3)

        # Label on the right side: |n|², electron count, degeneracy
        occ_str = (f'{sh["electrons"]}/{sh["capacity"]} e⁻'
                   if sh['frac'] < 0.999
                   else f'{sh["capacity"]} e⁻')
        ax.text(0.88, E_eV,
                f'$|n|^2={sh["n2"]}$  {occ_str}  ×{sh["deg"]}',
                va='center', ha='left', fontsize=8.5, color='#333333')

    # Fermi energy
    EF_plot = (EF + offset_Ha) * HA_TO_EV
    ax.axhline(EF_plot, color='#d73027', ls='--', lw=1.4, zorder=4,
               label=f'$E_F = {EF*HA_TO_EV:.2f}$ eV')

    ax.set_xlim(0, 1.85)
    ax.set_xticks([])
    ax.set_xlabel(title, fontsize=10, labelpad=10)
    ax.set_title(title.split('  ')[0], fontsize=11, pad=8)
    ax.spines[['top', 'right', 'bottom']].set_visible(False)
    ax.grid(axis='y', ls=':', alpha=0.35, zorder=0)
    ax.tick_params(axis='y', labelsize=10)

ax_left.set_ylabel('Energy  (eV)', fontsize=11)

# ── V_xc shift annotation between panels ─────────────────────────────────────
# Draw a double-headed arrow on the figure spanning between the two panels.
# We annotate at the |n|²=0 level (E_k=0 → most visible gap).

E0_left_eV  = 0.0
E0_right_eV = Vxc_eV   # negative (downward shift)

# Use figure-level annotation via ax_right's coordinate system
ax_right.annotate(
    '',
    xy    =(0.12, E0_right_eV),
    xytext=(0.12, E0_left_eV),
    xycoords='data', textcoords='data',
    arrowprops=dict(arrowstyle='<->', color='#555555', lw=1.2),
    zorder=5
)
ax_right.text(0.08, 0.5*(E0_left_eV + E0_right_eV),
              f'$V_{{xc}} = {Vxc_eV:.2f}$ eV',
              ha='right', va='center', fontsize=9, color='#555555',
              rotation=90)

# ── Legend ────────────────────────────────────────────────────────────────────

legend_handles = [
    mpatches.Patch(color=COL_FULL,    alpha=0.90, label='Fully occupied'),
    mpatches.Patch(color=COL_PARTIAL, alpha=0.90, label='Partially occupied'),
    mpatches.Patch(color=COL_EMPTY,   alpha=0.55, label='Empty'),
    Line2D([0], [0], color='#d73027', ls='--', lw=1.4,
           label=f'$E_F$ (continuum) $= {EF*HA_TO_EV:.2f}$ eV'),
]
fig.legend(handles=legend_handles, loc='lower center', ncol=4,
           fontsize=9, frameon=True, framealpha=0.95,
           bbox_to_anchor=(0.5, -0.06))

# ── System info subtitle ──────────────────────────────────────────────────────

fig.suptitle(
    f'Jellium free-electron shell structure\n'
    f'$N={N}$,  $L={L}\\,a_0$,  $r_s = {rs:.3f}\\,a_0$  '
    f'(PZ81 LDA,  $V_{{xc}} = {Vxc_eV:.2f}$ eV)',
    fontsize=12, y=1.02
)

plt.tight_layout()
fig.savefig('shell_structure.pdf', bbox_inches='tight', dpi=150)
fig.savefig('shell_structure.png', bbox_inches='tight', dpi=150)
print(f'Saved shell_structure.pdf / .png')
print(f'V_xc = {Vxc_eV:.4f} eV  ({Vxc_Ha:.6f} Ha)')
print(f'E_F  = {EF*HA_TO_EV:.4f} eV')
