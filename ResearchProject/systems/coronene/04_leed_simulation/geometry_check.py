#!/usr/bin/env python3
"""
geometry_check.py
=================
Plots:
  1. Coronene atoms (xy plane) with WP Gaussian footprint superimposed —
     confirms geometric alignment of WP centre with molecule centre.
  2. WP orbital density |psi^WP(x,y,z,t)|^2 at several (t, z) slices,
     using the free-space analytical formula (valid before collision).
     Saves frames to results/wp_orbital/.

Free-space WP formula (paper Eq. 1, propagated in free space):
  psi(r,t) = N(t) * exp(-|r - b(t)|^2 / (2 sigma(t)^2)) * exp(i*(k.r - omega*t))
where
  b(t) = (0, 0, D - k0*t)          WP centre (moving in -z)
  sigma(t) = d * sqrt(1 + (t/d^2)^2)   spreading width (atomic units)
  N(t) = normalisation factor

For the heatmaps we only need |psi|^2, which on a z=z0 plane is:
  n(x,y; z0, t) = |N(t)|^2 * exp(-(x^2 + y^2 + (z0-bz(t))^2) / sigma(t)^2)

All units atomic units unless stated.
"""

import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import LogNorm

# ── Parameters (must match config.hpp) ───────────────────────────────────────
ANG_TO_BOHR   = 1.8897259886
BOHR_TO_ANG   = 0.529177210903

LX_ANG = 18.4;  LY_ANG = 18.4
WP_D_ANG      = 0.53          # current WP width (Å)
WP_D_BOHR     = WP_D_ANG * ANG_TO_BOHR
WP_D_IMPACT_BOHR = 11.9998    # D (bohr)
WP_EKIN_EV    = 200.0
WP_EKIN_HA    = WP_EKIN_EV / 27.21138625
k0            = np.sqrt(2.0 * WP_EKIN_HA)   # bohr^-1
DT_AU         = 0.0200092
T2_AU         = 10.3353
t1_au         = WP_D_IMPACT_BOHR / k0        # WP arrival time

# Coronene atom positions (from coronene.xyz, Angstrom)
CARBONS = np.array([
    [ 1.421,  0.000], [ 2.842,  0.000],
    [ 0.7105, 1.2306], [ 1.421,  2.4612],
    [-0.7105, 1.2306], [-1.421,  2.4612],
    [-1.421,  0.000], [-2.842,  0.000],
    [-0.7105,-1.2306], [-1.421, -2.4612],
    [ 0.7105,-1.2306], [ 1.421, -2.4612],
    [ 3.5525, 1.2306], [ 2.842,  2.4612],
    [ 0.7105, 3.6919], [-0.7105, 3.6919],
    [-2.842,  2.4612], [-3.5525, 1.2306],
    [-3.5525,-1.2306], [-2.842, -2.4612],
    [-0.7105,-3.6919], [ 0.7105,-3.6919],
    [ 2.842, -2.4612], [ 3.5525,-1.2306],
])
HYDROGENS = np.array([
    [ 4.5787, 1.5861], [ 3.6629, 3.1722],
    [ 0.9157, 4.7583], [-0.9157, 4.7583],
    [-3.6629, 3.1722], [-4.5787, 1.5861],
    [-4.5787,-1.5861], [-3.6629,-3.1722],
    [-0.9157,-4.7583], [ 0.9157,-4.7583],
    [ 3.6629,-3.1722], [ 4.5787,-1.5861],
])

RESULTS    = "results"
WP_RESULTS = os.path.join(RESULTS, "wp_orbital")
os.makedirs(WP_RESULTS, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# 1. Geometry alignment plot
# ─────────────────────────────────────────────────────────────────────────────
def plot_geometry_alignment(d_ang=WP_D_ANG):
    fig, ax = plt.subplots(figsize=(6, 6))

    # Draw bonds (nearest-neighbour C-C, threshold 1.5 Å)
    for i, ci in enumerate(CARBONS):
        for j, cj in enumerate(CARBONS):
            if j <= i: continue
            if np.linalg.norm(ci - cj) < 1.50:
                ax.plot([ci[0], cj[0]], [ci[1], cj[1]], 'k-', lw=1.5, zorder=1)
    for i, ci in enumerate(CARBONS):
        for j, hj in enumerate(HYDROGENS):
            if np.linalg.norm(ci - hj) < 1.20:
                ax.plot([ci[0], hj[0]], [ci[1], hj[1]], 'k-', lw=1.0, zorder=1)

    # Atoms
    ax.scatter(CARBONS[:,0],   CARBONS[:,1],   s=120, c='dimgray',  zorder=3, label='C')
    ax.scatter(HYDROGENS[:,0], HYDROGENS[:,1], s=50,  c='lightblue', zorder=3, label='H')
    ax.scatter([0], [0], s=60, c='red', marker='+', zorder=4, lw=2, label='WP centre (0,0)')

    # Gaussian footprint circles: 1σ, 2σ, 3σ
    for nsig, alpha, ls in [(1, 0.8, '-'), (2, 0.4, '--'), (3, 0.2, ':')]:
        c = Circle((0, 0), nsig * d_ang, fill=False, edgecolor='red',
                   linestyle=ls, linewidth=1.8, alpha=alpha, zorder=5,
                   label=f'{nsig}σ  ({nsig*d_ang:.2f} Å)')
    ax.add_patch(Circle((0, 0), 1*d_ang, fill=False, edgecolor='red', lw=1.8, ls='-',  alpha=0.8, zorder=5))
    ax.add_patch(Circle((0, 0), 2*d_ang, fill=False, edgecolor='red', lw=1.8, ls='--', alpha=0.4, zorder=5))
    ax.add_patch(Circle((0, 0), 3*d_ang, fill=False, edgecolor='red', lw=1.8, ls=':',  alpha=0.2, zorder=5))

    # Ring radii annotation
    for r, label in [(1.421, 'inner ring\n(1.421 Å)'), (2.842, 'middle ring\n(2.842 Å)'),
                     (3.760, 'outer ring\n(3.760 Å)')]:
        ax.add_patch(Circle((0, 0), r, fill=False, edgecolor='steelblue',
                            linestyle=':', lw=0.8, alpha=0.5))

    ax.set_xlim(-6, 6); ax.set_ylim(-6, 6)
    ax.set_aspect('equal')
    ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
    ax.set_title(f'Coronene + WP Gaussian footprint\nd = {d_ang:.2f} Å '
                 f'(red circles: 1σ, 2σ, 3σ)\nBlue dashed: ring radii at 1.421, 2.842, 3.760 Å',
                 fontsize=11)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.2)

    out = os.path.join(RESULTS, f"geometry_alignment_d{d_ang:.2f}A.png")
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
# 2. WP orbital density heatmaps at (t, z) slices (free-space approximation)
# ─────────────────────────────────────────────────────────────────────────────
def wp_density_xy(x2d, y2d, bz_bohr, z_bohr, d_bohr, sigma_bohr):
    """
    |psi^WP(x,y,z)|^2 on a 2D (x,y) grid at fixed z.
    Uses free-space Gaussian formula; ignores time-dependent normalisation
    (we normalise to peak=1 for display).
    """
    dz = z_bohr - bz_bohr
    density = np.exp(-(x2d**2 + y2d**2 + dz**2) / sigma_bohr**2)
    return density


def plot_wp_slices():
    """
    For a set of simulation times t, plot |psi^WP|^2 in the xy plane at:
      - z = z_center(t)       (WP centre plane)
      - z = z_center(t) + d   (1 sigma above centre)
      - z = z_center(t) - d   (1 sigma below centre, towards molecule)
    """
    d = WP_D_BOHR

    # Time points: before, near arrival, during, after
    t_labels = [
        (0.0,     'initial (t=0)'),
        (t1_au/2, f'halfway (t={t1_au/2:.1f} a.u.)'),
        (t1_au,   f'arrival (t={t1_au:.1f} a.u.)'),
        (t1_au + 1.0, f'during (t={t1_au+1.0:.1f} a.u.)'),
        (T2_AU,   f'end (t={T2_AU:.1f} a.u.)'),
    ]

    # xy grid in bohr (-Lx/2 to Lx/2)
    Lx_bohr = LX_ANG * ANG_TO_BOHR
    Ly_bohr = LY_ANG * ANG_TO_BOHR
    N = 200
    xs = np.linspace(-Lx_bohr/2, Lx_bohr/2, N)
    ys = np.linspace(-Ly_bohr/2, Ly_bohr/2, N)
    x2d, y2d = np.meshgrid(xs, ys)
    extent_ang = [-Lx_bohr/2*BOHR_TO_ANG, Lx_bohr/2*BOHR_TO_ANG,
                  -Ly_bohr/2*BOHR_TO_ANG, Ly_bohr/2*BOHR_TO_ANG]

    for t_au, t_label in t_labels:
        bz = WP_D_IMPACT_BOHR - k0 * t_au          # WP centre z(t)
        sigma = d * np.sqrt(1 + (t_au / d**2)**2)  # spreading width

        z_slices = [
            (bz,     'z=centre',    'centre plane (peak density)'),
            (bz + d, 'z=ctr+1sig',  'z = centre + 1σ'),
            (bz - d, 'z=ctr-1sig',  'z = centre − 1σ (towards molecule)'),
        ]

        fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
        for ax, (z0, z_tag, z_title) in zip(axes, z_slices):
            dens = wp_density_xy(x2d, y2d, bz, z0, d, sigma)
            im = ax.imshow(dens, origin='lower', cmap='inferno',
                           vmin=0, vmax=1.0, extent=extent_ang, aspect='equal')
            # Overlay coronene atoms
            ax.scatter(CARBONS[:,0],   CARBONS[:,1],   s=15, c='white',  zorder=3, alpha=0.7)
            ax.scatter(HYDROGENS[:,0], HYDROGENS[:,1], s=8,  c='cyan',   zorder=3, alpha=0.5)
            ax.set_title(f'{z_title}\nz_WP={bz*BOHR_TO_ANG:.2f} Å, z_slice={z0*BOHR_TO_ANG:.2f} Å',
                         fontsize=9)
            ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
            plt.colorbar(im, ax=ax, label='|ψ|² (norm. to peak)')

        bz_ang = bz * BOHR_TO_ANG
        sigma_ang = sigma * BOHR_TO_ANG
        fig.suptitle(f'WP orbital density — {t_label}\n'
                     f'σ(t) = {sigma_ang:.3f} Å  |  z_centre = {bz_ang:.2f} Å  '
                     f'(free-space approx.)',
                     fontsize=11)
        plt.tight_layout()

        fname = f"wp_t{t_au:.2f}au".replace('.', 'p') + ".png"
        out = os.path.join(WP_RESULTS, fname)
        fig.savefig(out, dpi=130, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
# 3. Summary: WP overlap vs d for each ring
# ─────────────────────────────────────────────────────────────────────────────
def plot_overlap_vs_d():
    d_vals = np.linspace(0.3, 3.0, 300)
    r_central = 1.421; r_middle = 2.842; r_outer = 3.760

    fig, ax = plt.subplots(figsize=(7, 4))
    for r, label, color in [
        (r_central, f'Central ring  r={r_central} Å', 'steelblue'),
        (r_middle,  f'Middle ring   r={r_middle} Å',  'darkorange'),
        (r_outer,   f'Outer ring    r={r_outer} Å',   'green'),
    ]:
        ax.plot(d_vals, np.exp(-(r/d_vals)**2), label=label, color=color, lw=2)

    ax.axvline(WP_D_ANG, color='red', ls='--', lw=1.5, label=f'Current d = {WP_D_ANG} Å')
    ax.axhline(0.5, color='gray', ls=':', lw=1, alpha=0.7, label='50% level')
    ax.axhline(0.1, color='gray', ls=':', lw=0.8, alpha=0.5, label='10% level')

    ax.set_xlabel('WP width d (Å)', fontsize=12)
    ax.set_ylabel('Density at ring radius (norm.)', fontsize=12)
    ax.set_title('WP illumination of each coronene ring vs. width d\n'
                 'n(r) = exp(−r²/d²)', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.3, 3.0); ax.set_ylim(0, 1)

    out = os.path.join(RESULTS, "wp_ring_overlap_vs_d.png")
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print('=== Geometry alignment plot ===')
    plot_geometry_alignment(WP_D_ANG)

    print('\n=== WP overlap vs d sweep ===')
    plot_overlap_vs_d()

    print('\n=== WP orbital slices (free-space) ===')
    plot_wp_slices()

    print('\nDone.')
