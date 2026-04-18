#!/usr/bin/env python3
"""
geometry_check.py — Validate coronene geometry and produce geometry_check.png.

Checks:
  1. All atom coordinates inside [0, Lx] × [0, Ly] × [0, Lz]
  2. Molecular centroid at (Lx/2, Ly/2, Lz/2) ± 0.01 Å
  3. Molecule flat: all z within Lz/2 ± 0.01 Å
  4. WP injection point bz = Lz/2 + D is inside the cell

INQ cell convention: (0,0,0) is the CELL CORNER. Grid runs from 0 → L.
This was confirmed by run_001 failure (molecule at corner created 4 copies).

Reference: Tsubonoya, Hu, Watanabe PRB 90, 035416 (2014)
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle

# ── Cell and WP parameters (must match config.hpp) ───────────────────────────
LX_ANG = 18.4
LY_ANG = 18.4
LZ_ANG = 31.7
WP_D_ANG = 1.4        # WP Gaussian width σ
WP_D_IMPACT_ANG = 6.35  # WP-to-molecule distance D
ANG_TO_BOHR = 1.8897259886

Z_FLAKE  = LZ_ANG / 2.0                         # 15.850 Å
Z_START  = LZ_ANG / 2.0 + WP_D_IMPACT_ANG      # 22.200 Å (WP initial centre)
Z_OBS    = Z_START                               # LEED observation plane = z_start
Z_MID    = (Z_FLAKE + Z_OBS) / 2.0              # 19.025 Å

# ── Parse XYZ file ────────────────────────────────────────────────────────────
def parse_xyz(filename):
    atoms, xs, ys, zs = [], [], [], []
    with open(filename) as f:
        n = int(f.readline())
        _ = f.readline()  # comment
        for _ in range(n):
            parts = f.readline().split()
            atoms.append(parts[0])
            xs.append(float(parts[1]))
            ys.append(float(parts[2]))
            zs.append(float(parts[3]))
    return atoms, np.array(xs), np.array(ys), np.array(zs)

atoms, xs, ys, zs = parse_xyz('coronene_centered.xyz')
n_atoms = len(atoms)

# ── Print key geometry values ─────────────────────────────────────────────────
print("=" * 60)
print("Geometry check — coronene_centered.xyz")
print("=" * 60)
print(f"Cell (INQ corner-origin convention: origin at (0,0,0)):")
print(f"  Lx = {LX_ANG:.3f} Å  Ly = {LY_ANG:.3f} Å  Lz = {LZ_ANG:.3f} Å")
print()
print(f"Key z planes:")
print(f"  z_flake = {Z_FLAKE:.3f} Å = {Z_FLAKE*ANG_TO_BOHR:.3f} bohr  [molecule plane = Lz/2]")
print(f"  z_start = {Z_START:.3f} Å = {Z_START*ANG_TO_BOHR:.3f} bohr  [WP initial centre = Lz/2 + D]")
print(f"  z_obs   = {Z_OBS:.3f} Å = {Z_OBS*ANG_TO_BOHR:.3f} bohr  [LEED observation = z_start]")
print(f"  z_mid   = {Z_MID:.3f} Å = {Z_MID*ANG_TO_BOHR:.3f} bohr  [midpoint]")
print(f"  D       = {WP_D_IMPACT_ANG:.3f} Å = {WP_D_IMPACT_ANG*ANG_TO_BOHR:.3f} bohr  [WP impact distance]")
print()

cx = np.mean(xs[np.array(atoms) == 'C'])
cy = np.mean(ys[np.array(atoms) == 'C'])
cz = np.mean(zs[np.array(atoms) == 'C'])
print(f"Molecular centroid (C atoms only): ({cx:.3f}, {cy:.3f}, {cz:.3f}) Å")
print(f"Expected centroid:                 ({LX_ANG/2:.3f}, {LY_ANG/2:.3f}, {LZ_ANG/2:.3f}) Å")
print()

# ── Validation checks ─────────────────────────────────────────────────────────
print("Validation checks:")
failures = []

def check(cond, msg_ok, msg_fail):
    if cond:
        print(f"  [PASS] {msg_ok}")
    else:
        print(f"  [FAIL] {msg_fail}")
        failures.append(msg_fail)

check(np.all(xs >= 0) and np.all(xs <= LX_ANG),
      f"All x-coords in [0, {LX_ANG}] Å",
      f"x-coords out of range: min={xs.min():.3f} max={xs.max():.3f}")

check(np.all(ys >= 0) and np.all(ys <= LY_ANG),
      f"All y-coords in [0, {LY_ANG}] Å",
      f"y-coords out of range: min={ys.min():.3f} max={ys.max():.3f}")

check(np.all(zs >= 0) and np.all(zs <= LZ_ANG),
      f"All z-coords in [0, {LZ_ANG}] Å",
      f"z-coords out of range: min={zs.min():.3f} max={zs.max():.3f}")

tol = 0.01
check(abs(cx - LX_ANG/2) < tol and abs(cy - LY_ANG/2) < tol,
      f"Centroid at (Lx/2, Ly/2) = ({LX_ANG/2:.3f}, {LY_ANG/2:.3f}) Å ± {tol}",
      f"Centroid offset from cell centre: Δx={cx-LX_ANG/2:.4f} Δy={cy-LY_ANG/2:.4f} Å")

c_atoms = np.array(atoms) == 'C'
check(np.all(np.abs(zs[c_atoms] - LZ_ANG/2) < tol),
      f"All C atoms flat at z = Lz/2 = {LZ_ANG/2:.3f} Å ± {tol}",
      f"C atoms not flat: z range [{zs[c_atoms].min():.4f}, {zs[c_atoms].max():.4f}] Å")

check(Z_START < LZ_ANG,
      f"WP start z={Z_START:.3f} Å is inside cell (< Lz={LZ_ANG:.3f} Å)",
      f"WP start z={Z_START:.3f} Å is OUTSIDE cell!")

check(n_atoms == 36,
      f"36 atoms (C24H12 coronene)",
      f"Wrong atom count: {n_atoms}")

n_C = np.sum(np.array(atoms) == 'C')
n_H = np.sum(np.array(atoms) == 'H')
check(n_C == 24 and n_H == 12,
      f"24 C + 12 H = coronene C24H12",
      f"Wrong composition: {n_C} C + {n_H} H")

print()
if failures:
    print(f"GEOMETRY CHECK FAILED: {len(failures)} failure(s). Do not proceed with simulation.")
    sys.exit(1)
else:
    print("All checks PASSED. Geometry is correct for INQ run_003.")

# ── Build bond connectivity for plot ─────────────────────────────────────────
def bond_pairs(atoms, xs, ys, zs, C_C=1.55, C_H=1.20):
    """Return list of (i,j) bonded pairs based on distance threshold."""
    pairs = []
    n = len(atoms)
    for i in range(n):
        for j in range(i+1, n):
            d = np.sqrt((xs[i]-xs[j])**2 + (ys[i]-ys[j])**2 + (zs[i]-zs[j])**2)
            thr = C_C if (atoms[i]=='C' and atoms[j]=='C') else C_H
            if d < thr:
                pairs.append((i,j))
    return pairs

bonds = bond_pairs(atoms, xs, ys, zs)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("run_003: Coronene LEED geometry check\n"
             "Cell origin at (0,0,0) — INQ corner convention", fontsize=11)

# Atom colours
col = {'C': '#404040', 'H': '#aaaaaa'}
sz  = {'C': 60, 'H': 25}

# ── Left panel: side view (x vs z) ───────────────────────────────────────────
ax = axes[0]
ax.set_title("Side view (x–z plane, y at cell centre)")

# Bonds (project onto x-z)
for i, j in bonds:
    ax.plot([xs[i], xs[j]], [zs[i], zs[j]], 'k-', lw=0.8, alpha=0.4)

# Atoms
for s, x, z in zip(atoms, xs, zs):
    ax.scatter(x, z, c=col[s], s=sz[s], zorder=3)

# Cell boundaries
for z_val, lbl in [(0, 'z=0'), (LZ_ANG, f'z=Lz={LZ_ANG} Å')]:
    ax.axhline(z_val, color='black', lw=1.0, ls='-')
    ax.text(0.5, z_val + 0.3, lbl, fontsize=7, color='black')

# Key z planes
plane_style = dict(lw=1.2, alpha=0.8)
ax.axhline(Z_FLAKE, color='steelblue', ls='--', label=f'z_flake={Z_FLAKE:.2f} Å (molecule)', **plane_style)
ax.axhline(Z_OBS,   color='firebrick',  ls='--', label=f'z_obs={Z_OBS:.2f} Å (LEED screen)', **plane_style)
ax.axhline(Z_MID,   color='darkorange', ls=':',  label=f'z_mid={Z_MID:.2f} Å', **plane_style)

# WP Gaussian sketch at z_start
wp_z0 = Z_START
sigma = WP_D_ANG
wp_x = np.linspace(0, LX_ANG, 300)
wp_z = wp_z0 + 1.5 * np.exp(-0.5 * ((wp_x - LX_ANG/2) / sigma)**2)  # offset for visibility
ax.plot(wp_x, wp_z, 'r-', lw=1.5, label=f'WP Gaussian (d={WP_D_ANG} Å)', zorder=4)
ax.scatter(LX_ANG/2, wp_z0, c='red', s=50, zorder=5, marker='x')

ax.set_xlabel('x (Å)')
ax.set_ylabel('z (Å)')
ax.set_xlim(-0.5, LX_ANG + 0.5)
ax.set_ylim(-1.0, LZ_ANG + 1.5)
ax.legend(fontsize=7, loc='upper right')

# ── Right panel: top view (x vs y) ───────────────────────────────────────────
ax = axes[1]
ax.set_title("Top view (x–y plane, coronene at z=Lz/2)")

# Bonds (project onto x-y)
for i, j in bonds:
    ax.plot([xs[i], xs[j]], [ys[i], ys[j]], 'k-', lw=0.8, alpha=0.4)

# Atoms
for s, x, y in zip(atoms, xs, ys):
    ax.scatter(x, y, c=col[s], s=sz[s], zorder=3, edgecolors='none')

# Centroid cross-hair
cx_plot = LX_ANG / 2
cy_plot = LY_ANG / 2
ax.axhline(cy_plot, color='grey', lw=0.6, ls=':')
ax.axvline(cx_plot, color='grey', lw=0.6, ls=':')
ax.scatter(cx_plot, cy_plot, c='red', s=80, zorder=6, marker='+', linewidths=1.5,
           label=f'Centroid ({cx_plot:.1f}, {cy_plot:.1f})')

# WP beam footprint (1σ circle in red)
circle_1sig = Circle((cx_plot, cy_plot), WP_D_ANG, fill=False, edgecolor='red',
                      lw=1.5, ls='-', label=f'WP beam 1σ (d={WP_D_ANG} Å)')
ax.add_patch(circle_1sig)

# Cell boundary
rect = plt.Rectangle((0, 0), LX_ANG, LY_ANG, fill=False, edgecolor='black', lw=1.2)
ax.add_patch(rect)

# Inner ring radius label
r_inner = 1.421  # C-C nearest distance (inner ring)
circle_inner = Circle((cx_plot, cy_plot), r_inner, fill=False, edgecolor='steelblue',
                       lw=1.0, ls='--', label=f'Inner ring r={r_inner} Å')
ax.add_patch(circle_inner)

ax.set_xlabel('x (Å)')
ax.set_ylabel('y (Å)')
ax.set_xlim(-0.5, LX_ANG + 0.5)
ax.set_ylim(-0.5, LY_ANG + 0.5)
ax.set_aspect('equal')
ax.legend(fontsize=7, loc='upper right')

# Colour legend for atoms
c_patch = mpatches.Patch(color='#404040', label='C')
h_patch = mpatches.Patch(color='#aaaaaa', label='H')
axes[1].legend(handles=[c_patch, h_patch,
                         mpatches.Patch(color='red', label=f'WP beam 1σ (d={WP_D_ANG} Å)'),
                         mpatches.Patch(color='steelblue', label=f'Inner ring r={r_inner} Å')],
               fontsize=7, loc='upper right')

plt.tight_layout()
plt.savefig('geometry_check.png', dpi=150, bbox_inches='tight')
print("\nSaved: geometry_check.png")
