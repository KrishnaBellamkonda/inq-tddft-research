#!/usr/bin/env python3
"""
plot_results.py — Visualisation suite for the jellium ground-state.

Run from the 01_ground_state/ directory after running inq-run:

    python3 plot_results.py

Produces three figures, all saved as PNG to results/:
  1. shell_structure.png   — free-electron vs KS energy level diagram
  2. xc_offset.png         — scatter plot: ε_i vs k²/2, checks slope=1 intercept=V_xc
  3. orbitals.png          — 2D slices of Re[ψ_k], Im[ψ_k], and |ψ_k|² for each shell
                             (reads separate _real.txt and _imag.txt orbital files)

All plots use a consistent serif style suitable for publication or reports.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

plt.rcParams.update({
    'font.family':      'serif',
    'font.size':        11,
    'axes.linewidth':   0.8,
    'xtick.direction':  'in',
    'ytick.direction':  'in',
    'legend.framealpha': 0.95,
})

os.makedirs('results', exist_ok=True)

# ── System parameters ─────────────────────────────────────────────────────────

N  = 40
L  = 40.0       # bohr
HA_TO_EV = 27.211386245988

n0  = N / L**3
rs  = (3.0 / (4.0 * np.pi * n0))**(1.0/3.0)
kF  = (3.0 * np.pi**2 * n0)**(1.0/3.0)
EF  = 0.5 * kF**2
k0  = 2.0 * np.pi / L

# ── PZ81 LDA (Perdew & Zunger, PRB 23, 5048, 1981) ───────────────────────────

def vxc_pz81(rs):
    ex  = -0.4582 / rs;  vx = (4.0/3.0) * ex
    g, b1, b2 = -0.1423, 1.0529, 0.3334
    srs = np.sqrt(rs)
    d   = 1.0 + b1*srs + b2*rs
    ec  = g / d
    vc  = ec * (1.0 + (7.0/6.0)*b1*srs + (4.0/3.0)*b2*rs) / d
    return vx + vc

Vxc_Ha = vxc_pz81(rs)
Vxc_eV = Vxc_Ha * HA_TO_EV

# ── Shell structure ───────────────────────────────────────────────────────────

def compute_shells(L, N, n2_max=6):
    k0 = 2.0 * np.pi / L
    degen = {}
    nmax = int(np.sqrt(n2_max)) + 1
    for nx in range(-nmax, nmax+1):
        for ny in range(-nmax, nmax+1):
            for nz in range(-nmax, nmax+1):
                n2 = nx**2 + ny**2 + nz**2
                if n2 <= n2_max:
                    degen[n2] = degen.get(n2, 0) + 1
    shells = []
    remaining = N
    for n2 in sorted(degen):
        deg = degen[n2]
        cap = 2 * deg
        fill = min(remaining, cap)
        remaining -= fill
        shells.append({'n2': n2, 'deg': deg,
                       'Ek_Ha': 0.5 * n2 * k0**2,
                       'electrons': fill, 'cap': cap,
                       'frac': fill / cap})
    return shells

shells = compute_shells(L, N)

COL_FULL    = '#2166ac'
COL_PARTIAL = '#e6711f'
COL_EMPTY   = '#aaaaaa'

def shell_style(frac):
    if frac >= 0.999:   return COL_FULL,    0.90
    elif frac > 0.001:  return COL_PARTIAL, 0.90
    else:               return COL_EMPTY,   0.55

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Shell structure level diagram
# ═════════════════════════════════════════════════════════════════════════════

def plot_shell_structure():
    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(9, 6), sharey=True,
                                      gridspec_kw={'width_ratios': [1,1], 'wspace': 0.04})

    LW_MAX = 10
    max_deg = max(sh['deg'] for sh in shells)

    for ax, offset_Ha in [(ax_l, 0.0), (ax_r, Vxc_Ha)]:
        for sh in shells:
            E_eV = (sh['Ek_Ha'] + offset_Ha) * HA_TO_EV
            col, alpha = shell_style(sh['frac'])
            lw = max(2.5, LW_MAX * sh['deg'] / max_deg)

            ax.plot([0.15, 0.85], [E_eV, E_eV],
                    color=col, lw=lw, alpha=alpha, solid_capstyle='round', zorder=3)

            occ = (f'{sh["electrons"]}/{sh["cap"]} e⁻' if sh['frac'] < 0.999
                   else f'{sh["cap"]} e⁻')
            ax.text(0.88, E_eV,
                    f'$|n|^2={sh["n2"]}$  {occ}  ×{sh["deg"]}',
                    va='center', ha='left', fontsize=8.5, color='#333333')

        EF_plot = (EF + offset_Ha) * HA_TO_EV
        ax.axhline(EF_plot, color='#d73027', ls='--', lw=1.4, zorder=4)
        ax.set_xlim(0, 1.85)
        ax.set_xticks([])
        ax.spines[['top', 'right', 'bottom']].set_visible(False)
        ax.grid(axis='y', ls=':', alpha=0.35)
        ax.tick_params(axis='y', labelsize=10)

    ax_l.set_ylabel('Energy  (eV)', fontsize=12)
    ax_l.set_xlabel(r'Free-electron  $\varepsilon_k = k^2/2$', fontsize=10, labelpad=10)
    ax_r.set_xlabel(r'KS eigenvalue  $\varepsilon_k = k^2/2 + V_{xc}$', fontsize=10, labelpad=10)

    # V_xc bracket annotation
    ax_r.annotate('', xy=(0.12, Vxc_eV), xytext=(0.12, 0.0),
                  arrowprops=dict(arrowstyle='<->', color='#555555', lw=1.2))
    ax_r.text(0.08, 0.5*Vxc_eV,
              f'$V_{{xc}} = {Vxc_eV:.2f}$ eV',
              ha='right', va='center', fontsize=9, color='#555555', rotation=90)

    legend_handles = [
        mpatches.Patch(color=COL_FULL, alpha=0.9, label='Fully occupied'),
        mpatches.Patch(color=COL_PARTIAL, alpha=0.9, label='Partially occupied'),
        mpatches.Patch(color=COL_EMPTY, alpha=0.55, label='Empty'),
        Line2D([0],[0], color='#d73027', ls='--', lw=1.4,
               label=f'$E_F = {EF*HA_TO_EV:.3f}$ eV'),
    ]
    fig.legend(handles=legend_handles, loc='lower center', ncol=4,
               fontsize=9, frameon=True, bbox_to_anchor=(0.5, -0.05))

    fig.suptitle(
        f'Jellium free-electron shell structure\n'
        f'$N={N}$,  $L={L}\\,a_0$,  $r_s={rs:.3f}\\,a_0$  '
        f'(PZ81 LDA,  $V_{{xc}}={Vxc_eV:.2f}$ eV)',
        fontsize=12, y=1.02
    )
    plt.tight_layout()
    fig.savefig('results/shell_structure.png', bbox_inches='tight', dpi=150)
    print('Saved results/shell_structure.png')
    plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — XC offset verification (ε_i vs k²/2)
# ═════════════════════════════════════════════════════════════════════════════

def plot_xc_offset():
    evfile = 'results/eigenvalues.txt'
    if not os.path.exists(evfile):
        print(f'  [skip] {evfile} not found — run inq-run first')
        return

    data = np.loadtxt(evfile, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # Columns: state_idx  shell_n2  k2_over_2_Ha  eigenvalue_Ha  predicted_Ha  residual_Ha
    n2_vals  = data[:, 1].astype(int)
    Ek_Ha    = data[:, 2]
    ev_num   = data[:, 3]
    ev_pred  = data[:, 4]
    resid    = data[:, 5]

    unique_n2 = sorted(set(n2_vals))
    cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(unique_n2)))
    col_map = {n2: cmap[i] for i, n2 in enumerate(unique_n2)}

    fig, (ax_main, ax_res) = plt.subplots(2, 1, figsize=(7, 7), sharex=True,
                                           gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.06})

    # Main panel: ε_i vs k²/2
    for n2 in unique_n2:
        mask = n2_vals == n2
        ax_main.scatter(Ek_Ha[mask] * HA_TO_EV, ev_num[mask] * HA_TO_EV,
                        color=col_map[n2], s=60, zorder=4, label=f'$|n|^2={n2}$')

    # Reference line: slope 1, intercept V_xc
    x_range = np.array([Ek_Ha.min() - 0.005, Ek_Ha.max() + 0.005])
    ax_main.plot(x_range * HA_TO_EV, (x_range + Vxc_Ha) * HA_TO_EV,
                 'r--', lw=1.5, zorder=3,
                 label=f'$\\varepsilon = k^2/2 + V_{{xc}}$  ($V_{{xc}}={Vxc_eV:.3f}$ eV)')

    ax_main.set_ylabel(r'KS eigenvalue $\varepsilon_i$  (eV)', fontsize=12)
    ax_main.legend(fontsize=9, loc='upper left')
    ax_main.grid(ls=':', alpha=0.4)
    ax_main.set_title(
        r'XC offset verification: $\varepsilon_i = k_i^2/2 + V_{xc}(n_0)$',
        fontsize=12, pad=8
    )

    # Fit the line to check intercept
    coeffs = np.polyfit(Ek_Ha, ev_num, 1)
    slope, intercept = coeffs
    ax_main.text(0.05, 0.95,
                 f'Fit: slope = {slope:.4f}  (expected 1.0000)\n'
                 f'     intercept = {intercept*HA_TO_EV:.4f} eV  '
                 f'(V_xc predicted = {Vxc_eV:.4f} eV)',
                 transform=ax_main.transAxes, fontsize=9, va='top',
                 bbox=dict(boxstyle='round', fc='white', alpha=0.85, ec='#cccccc'))

    # Residual panel
    for n2 in unique_n2:
        mask = n2_vals == n2
        ax_res.scatter(Ek_Ha[mask] * HA_TO_EV, resid[mask] * 1000,
                       color=col_map[n2], s=40, zorder=4)
    ax_res.axhline(0, color='red', ls='--', lw=1.2)
    ax_res.set_xlabel(r'Free-electron energy $k^2/2$  (eV)', fontsize=12)
    ax_res.set_ylabel(r'Residual $\varepsilon_i - (k^2/2 + V_{xc})$  (mHa)', fontsize=10)
    ax_res.grid(ls=':', alpha=0.4)

    plt.tight_layout()
    fig.savefig('results/xc_offset.png', bbox_inches='tight', dpi=150)
    print('Saved results/xc_offset.png')
    plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — 2D orbital slices (Re[ψ_k] for each shell)
# ═════════════════════════════════════════════════════════════════════════════

def load_orbital(orbital_idx):
    """Load a 2D slice from separate _real.txt and _imag.txt files.

    Returns (Re_psi, Im_psi, density, x_vals, n2) or (None,)*5 if missing.
    density = Re² + Im²  (computed here; always 1/Ω for plane waves).
    """
    import glob
    re_files = glob.glob(f'results/orbitals/orbital_{orbital_idx}_n2_*_real.txt')
    im_files = glob.glob(f'results/orbitals/orbital_{orbital_idx}_n2_*_imag.txt')
    if not re_files or not im_files:
        return None, None, None, None, None
    fname_re = re_files[0]
    n2 = int(fname_re.split('_n2_')[1].split('_real')[0])

    # Columns in each file: ix  iy  x_bohr  y_bohr  psi_value
    raw_re = np.loadtxt(fname_re, comments='#')
    raw_im = np.loadtxt(im_files[0], comments='#')
    N_g    = int(raw_re[:, 0].max()) + 1
    Re_psi  = raw_re[:, 4].reshape(N_g, N_g)
    Im_psi  = raw_im[:, 4].reshape(N_g, N_g)
    density = Re_psi**2 + Im_psi**2   # = 1/Ω everywhere for plane waves
    x_vals  = raw_re[:N_g, 2]
    return Re_psi, Im_psi, density, x_vals, n2


def plot_orbitals():
    import glob

    # Count by unique orbital indices (look for _real.txt files)
    re_files = glob.glob('results/orbitals/orbital_*_real.txt')
    n_orbs = len(re_files)
    if n_orbs == 0:
        print('  [skip] No orbital files found — run inq-run first')
        return

    # 3 rows: Re[ψ], Im[ψ], density
    fig, axes = plt.subplots(3, n_orbs, figsize=(3*n_orbs, 9),
                              gridspec_kw={'hspace': 0.35, 'wspace': 0.08})
    if n_orbs == 1:
        axes = axes.reshape(3, 1)

    row_labels = [r'Re[$\psi_k$]', r'Im[$\psi_k$]', r'$|\psi_k|^2$']
    row_cmaps  = ['RdBu_r', 'PuOr_r', 'viridis']

    for oi in range(n_orbs):
        result = load_orbital(oi)
        if result[0] is None:
            continue
        Re_psi, Im_psi, density, x_vals, n2 = result
        fields = [Re_psi, Im_psi, density]
        extent = [0, L, 0, L]

        for row, (field, cmap) in enumerate(zip(fields, row_cmaps)):
            ax = axes[row, oi]
            if row < 2:
                # Centre the diverging colormap on zero
                lim = np.abs(field).max()
                im = ax.imshow(field.T, origin='lower', extent=extent,
                               cmap=cmap, vmin=-lim, vmax=lim, aspect='equal')
            else:
                # Density — tight colour scale to show uniformity
                vmin = field.min() * 0.999
                vmax = field.max() * 1.001
                im = ax.imshow(field.T, origin='lower', extent=extent,
                               cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal')
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            ax.set_title(f'$|n|^2={n2}$' if row == 0 else '', fontsize=10)
            if oi == 0:
                ax.set_ylabel(r'$y$  (bohr)', fontsize=9)
            if row == 2:
                ax.set_xlabel(r'$x$  (bohr)', fontsize=9)

    # Row labels on left side
    for row, label in enumerate(row_labels):
        axes[row, 0].set_ylabel(f'{label}\n$y$  (bohr)', fontsize=9)

    fig.suptitle(
        f'Jellium KS orbitals — analytical plane waves\n'
        f'$N={N}$, $L={L}\\,a_0$, $r_s={rs:.3f}\\,a_0$,  slice $z={L/2:.1f}$ bohr',
        fontsize=12, y=1.01
    )
    plt.tight_layout()
    fig.savefig('results/orbitals.png', bbox_inches='tight', dpi=120)
    print('Saved results/orbitals.png')
    plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print(f'System: N={N}, L={L} bohr, r_s={rs:.4f} bohr')
    print(f'V_xc (PZ81) = {Vxc_eV:.4f} eV,  E_F = {EF*HA_TO_EV:.4f} eV\n')

    print('Plotting shell structure...')
    plot_shell_structure()

    print('Plotting XC offset...')
    plot_xc_offset()

    print('Plotting orbital slices...')
    plot_orbitals()

    print('\nDone. Figures written to results/.')
