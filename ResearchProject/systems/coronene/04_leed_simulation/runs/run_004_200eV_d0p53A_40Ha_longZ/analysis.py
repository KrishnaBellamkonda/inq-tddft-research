#!/usr/bin/env python3
"""
analysis.py — Post-processing for run_004 (coronene TDDFT LEED)
Ref: Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ──────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────
RUN_DIR  = os.path.dirname(os.path.abspath(__file__))
RES_DIR  = os.path.join(RUN_DIR, 'results')
FIG_DIR  = os.path.join(RES_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

BOHR_TO_ANG      = 0.529177   # bohr → Å
BOHR_INV_TO_ANG_INV = 1.0 / BOHR_TO_ANG   # bohr⁻¹ → Å⁻¹

# ──────────────────────────────────────────────────────────────
# 1. Load simulation summary
# ──────────────────────────────────────────────────────────────
def load_summary(path):
    params = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                params[parts[0]] = parts[1]
    return params

summary = load_summary(os.path.join(RES_DIR, 'sim_summary.txt'))
print("=== Simulation summary ===")
for k, v in summary.items():
    print(f"  {k:<25} {v}")

# ──────────────────────────────────────────────────────────────
# 2. Load grid metadata
# ──────────────────────────────────────────────────────────────
with open(os.path.join(RES_DIR, 'grid', 'grid_metadata.txt')) as f:
    lines = [l for l in f if not l.startswith('#')]
meta = lines[0].split()
Nx, Ny, Nz = int(meta[0]), int(meta[1]), int(meta[2])
dx, dy, dz  = float(meta[3]), float(meta[4]), float(meta[5])   # bohr
Lx, Ly, Lz  = float(meta[6]), float(meta[7]), float(meta[8])   # bohr

print(f"\n=== Grid ===")
print(f"  Nx={Nx}  Ny={Ny}  Nz={Nz}")
print(f"  dx={dx:.5f} bohr  Lx={Lx:.5f} bohr = {Lx*BOHR_TO_ANG:.4f} Å")
print(f"  dz={dz:.5f} bohr  Lz={Lz:.5f} bohr = {Lz*BOHR_TO_ANG:.4f} Å")

# Spatial axes in Å
x_ang = np.arange(Nx) * dx * BOHR_TO_ANG
y_ang = np.arange(Ny) * dy * BOHR_TO_ANG
z_ang = np.arange(Nz) * dz * BOHR_TO_ANG

# Key parameters from summary
z_screen0_bohr = float(summary.get('z_screen0_bohr', '0'))
z_screen1_bohr = float(summary.get('z_screen1_bohr', '0'))
z_screen2_bohr = float(summary.get('z_screen2_bohr', '0'))
z_screen0_ang  = z_screen0_bohr * BOHR_TO_ANG
z_screen1_ang  = z_screen1_bohr * BOHR_TO_ANG
z_screen2_ang  = z_screen2_bohr * BOHR_TO_ANG
t_leed_start   = float(summary.get('t_leed_start_au', '0'))
t2             = float(summary.get('t2_au', '0'))
t1_au          = float(summary.get('t1_arrival_au', '0'))
k0             = float(summary.get('WP_k0_bohr_inv', '3.834'))   # bohr⁻¹
wp_d_bohr      = float(summary.get('WP_d_bohr', '1.0'))
wp_d_ang       = wp_d_bohr * BOHR_TO_ANG
wp_bz_bohr     = float(summary.get('WP_bz_bohr', '0'))
wp_D_bohr      = float(summary.get('WP_D_bohr', '0'))
z_flake_ang    = (wp_bz_bohr - wp_D_bohr) * BOHR_TO_ANG

print(f"\n  z_flake  = {z_flake_ang:.3f} Å")
print(f"  Screen 0 = {z_screen0_ang:.3f} Å  (= z_obs)")
print(f"  Screen 1 = {z_screen1_ang:.3f} Å")
print(f"  Screen 2 = {z_screen2_ang:.3f} Å")
print(f"  t_leed_start = {t_leed_start:.3f} a.u.  (= 10σ/k₀)")
print(f"  T1 = {t1_au:.3f} a.u.,  T2 = {t2:.3f} a.u.")

# ──────────────────────────────────────────────────────────────
# 3. Load LEED patterns
# ──────────────────────────────────────────────────────────────
def load_leed(filename):
    """Load Ny×Nx LEED pattern from space-separated text (5 header lines)."""
    data = []
    with open(filename) as f:
        for i, line in enumerate(f):
            if i < 5:
                continue
            vals = list(map(float, line.split()))
            if vals:
                data.append(vals)
    arr = np.array(data)   # shape (Ny, Nx)
    print(f"  {os.path.basename(filename)}: shape={arr.shape}  "
          f"min={arr.min():.3e}  max={arr.max():.3e}  sum={arr.sum():.3e}")
    return arr

print("\n=== Loading LEED patterns ===")
leed0 = load_leed(os.path.join(RES_DIR, 'leed_pattern', 'leed_screen0.txt'))
leed1 = load_leed(os.path.join(RES_DIR, 'leed_pattern', 'leed_screen1.txt'))
leed2 = load_leed(os.path.join(RES_DIR, 'leed_pattern', 'leed_screen2.txt'))

leed_list   = [leed0, leed1, leed2]
screen_labels = [
    f'Screen 0  z={z_screen0_ang:.2f} Å\n(= z_obs)',
    f'Screen 1  z={z_screen1_ang:.2f} Å',
    f'Screen 2  z={z_screen2_ang:.2f} Å',
]

# Centred spatial coordinates (molecule centroid at origin)
cx_ang = Lx * BOHR_TO_ANG / 2.0
cy_ang = Ly * BOHR_TO_ANG / 2.0
xc = x_ang - cx_ang
yc = y_ang - cy_ang

# ──────────────────────────────────────────────────────────────
# 4. Real-space LEED (Figure 1)
# ──────────────────────────────────────────────────────────────
print("\n=== Plotting real-space LEED ===")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
vmax = max(np.abs(L).max() for L in leed_list)

for ax, leed, label in zip(axes, leed_list, screen_labels):
    im = ax.pcolormesh(xc, yc, leed,
                       cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                       shading='auto')
    ax.set_xlabel('x − x_c  (Å)', fontsize=11)
    ax.set_ylabel('y − y_c  (Å)', fontsize=11)
    ax.set_title(label, fontsize=11)
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='∫Δn dt  (bohr⁻³·a.u.)')

fig.suptitle(
    f'Coronene TDDFT LEED — real-space (background subtracted)\n'
    f'200 eV | d={wp_d_ang:.3f} Å | t_leed=[{t_leed_start:.2f}, {t2:.2f}] a.u.',
    fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'leed_real_space.png'), dpi=150)
plt.close(fig)
print("  Saved: leed_real_space.png")

# ──────────────────────────────────────────────────────────────
# 5. k-space LEED (2D FFT) — Figures 2 & 3
# ──────────────────────────────────────────────────────────────
# k-space axes in Å⁻¹
kx_ang = np.fft.fftshift(np.fft.fftfreq(Nx, d=dx) * 2 * np.pi) * BOHR_INV_TO_ANG_INV
ky_ang = np.fft.fftshift(np.fft.fftfreq(Ny, d=dy) * 2 * np.pi) * BOHR_INV_TO_ANG_INV
KX, KY = np.meshgrid(kx_ang, ky_ang)

# Graphene/coronene reference ring radii (Å⁻¹)
# a = 2.46 Å, |g₁| = 4π/(√3·a), |g₂| = √3·|g₁|, ...
a_graph  = 2.46  # Å
g1 = 4 * np.pi / (np.sqrt(3) * a_graph)   # ≈ 2.95 Å⁻¹  (first BZ boundary)
g2 = np.sqrt(3) * g1                        # ≈ 5.12 Å⁻¹
g3 = 2.0 * g1                               # ≈ 5.90 Å⁻¹
g_refs = [(g1, '--', f'|g₁|={g1:.2f} Å⁻¹'),
          (g2, ':',  f'|g₂|={g2:.2f} Å⁻¹'),
          (g3, '-.', f'|g₃|={g3:.2f} Å⁻¹')]
theta = np.linspace(0, 2*np.pi, 300)

def plot_kspace(leed_list, labels, filename, klim):
    fig, axes = plt.subplots(1, 3, figsize=(5*3, 5))
    Ik_list = []
    for L in leed_list:
        F = np.fft.fftshift(np.fft.fft2(L))
        Ik_list.append(np.abs(F)**2)
    vmax = max(I.max() for I in Ik_list)

    for ax, Ik, label in zip(axes, Ik_list, labels):
        Ik_disp = Ik.copy()
        cy, cx = Ny // 2, Nx // 2
        Ik_disp[cy-2:cy+3, cx-2:cx+3] = 0   # suppress DC

        im = ax.pcolormesh(KX, KY, Ik_disp,
                           cmap='hot', vmin=0, vmax=vmax,
                           shading='auto')
        ax.set_xlim(-klim, klim)
        ax.set_ylim(-klim, klim)
        ax.set_xlabel('k_x  (Å⁻¹)', fontsize=11)
        ax.set_ylabel('k_y  (Å⁻¹)', fontsize=11)
        ax.set_title(label, fontsize=11)
        ax.set_aspect('equal')
        plt.colorbar(im, ax=ax, label='|FFT|²  (arb.)')

        for g_mag, ls, lbl in g_refs:
            ax.plot(g_mag*np.cos(theta), g_mag*np.sin(theta),
                    color='cyan', lw=0.8, ls=ls, alpha=0.8, label=lbl)
        ax.legend(fontsize=7, loc='upper right')

    fig.suptitle(
        f'Coronene TDDFT LEED — k-space (2D FFT)\n'
        f'200 eV | d={wp_d_ang:.3f} Å | Tsubonoya PRB 90, 035416 (2014)',
        fontsize=12)
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"  Saved: {os.path.basename(filename)}")

print("\n=== Plotting k-space LEED ===")
plot_kspace(leed_list, screen_labels,
            os.path.join(FIG_DIR, 'leed_kspace.png'), klim=15.0)
plot_kspace(leed_list, screen_labels,
            os.path.join(FIG_DIR, 'leed_kspace_zoomed.png'), klim=6.0)

# ──────────────────────────────────────────────────────────────
# 6. Energy vs time
# ──────────────────────────────────────────────────────────────
print("\n=== Energy conservation ===")
energy_data = np.genfromtxt(
    os.path.join(RES_DIR, 'energy', 'energy_vs_time.csv'),
    delimiter=',', comments='#')
t_e, E_e = energy_data[:, 1], energy_data[:, 2]
dE       = E_e - E_e[0]
drift    = (E_e[-1] - E_e[0]) / (t_e[-1] - t_e[0])
print(f"  E(0)   = {E_e[0]:.8f} Ha")
print(f"  E(end) = {E_e[-1]:.8f} Ha")
print(f"  Drift  = {drift:.3e} Ha/a.u.")

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_e, dE * 1000, 'b-', lw=1)
ax.axhline(0, color='k', lw=0.5, ls='--')
ax.axvline(t1_au, color='r', lw=1, ls='--', label=f'T1={t1_au:.2f} a.u.')
ax.axvline(t_leed_start, color='g', lw=1, ls=':', label=f't_leed={t_leed_start:.2f} a.u.')
ax.set_xlabel('t  (a.u.)', fontsize=12)
ax.set_ylabel('ΔE  (mHa)', fontsize=12)
ax.set_title(f'Energy conservation | drift = {drift:.2e} Ha/a.u.', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'energy_vs_time.png'), dpi=150)
plt.close(fig)
print("  Saved: energy_vs_time.png")

# ──────────────────────────────────────────────────────────────
# 7. Momentum vs time
# ──────────────────────────────────────────────────────────────
print("\n=== Momentum ===")
mom_data = np.genfromtxt(
    os.path.join(RES_DIR, 'momentum', 'momentum_vs_time.csv'),
    delimiter=',', comments='#')
t_m, Jx, Jy, Jz = mom_data[:,1], mom_data[:,2], mom_data[:,3], mom_data[:,4]
print(f"  Jz(t=0) = {Jz[0]:.4f} a.u.  k₀={k0:.4f} a.u.  "
      f"rel err={(abs(Jz[0])-k0)/k0*100:.2f}%")

fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
axes[0].plot(t_m, Jz, 'b-', lw=1, label='Jz')
axes[0].axhline( k0, color='gray', lw=0.8, ls='--', label=f'+k₀={k0:.3f}')
axes[0].axhline(-k0, color='gray', lw=0.8, ls='--', label=f'-k₀={-k0:.3f}')
axes[0].axvline(t1_au, color='r', lw=1, ls='--', label=f'T1={t1_au:.2f} a.u.')
axes[0].axvline(t_leed_start, color='g', lw=1, ls=':', label='t_leed')
axes[0].set_ylabel('Jz  (a.u.)', fontsize=11)
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.3)
axes[1].plot(t_m, Jx, 'r-', lw=0.8, label='Jx', alpha=0.7)
axes[1].plot(t_m, Jy, 'g-', lw=0.8, label='Jy', alpha=0.7)
axes[1].set_xlabel('t  (a.u.)', fontsize=12)
axes[1].set_ylabel('Jx, Jy  (a.u.)', fontsize=11)
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)
fig.suptitle('Total current (momentum) vs time', fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'momentum_vs_time.png'), dpi=150)
plt.close(fig)
print("  Saved: momentum_vs_time.png")

# ──────────────────────────────────────────────────────────────
# 8. WP z-centroid from z-profile heatmap
# ──────────────────────────────────────────────────────────────
print("\n=== z-profile heatmap and WP centroid ===")
zprof_data = np.genfromtxt(
    os.path.join(RES_DIR, 'wp_trajectory', 'density_z_profile_vs_time.csv'),
    delimiter=',', comments='#')
t_zp    = zprof_data[:, 1]       # a.u.
profile = zprof_data[:, 2:]      # (n_snap, Nz)

# Background: mean of last 5 snapshots (WP has left the grid or dispersed)
n_baseline = np.mean(profile[-5:, :], axis=0)
dprof = profile - n_baseline[np.newaxis, :]

# WP z-centroid: centroid of positive Δn
z_centroid = np.array([
    np.sum(z_ang * np.maximum(dp, 0)) / (np.sum(np.maximum(dp, 0)) + 1e-30)
    for dp in dprof
])

# Expected free-particle: starting z depends on WP direction
# Jz>0 at t=0 → WP initially moves in +z; from z_screen0 it wraps
# around and arrives at z_flake after traveling (Lz−z_screen0+z_flake)/k₀_ang_per_au
k0_ang = k0 * BOHR_TO_ANG   # Å per a.u.
z_expect_direct = z_screen0_ang + k0_ang * t_zp   # +z direction (mod Lz_ang)
Lz_ang = Lz * BOHR_TO_ANG
z_expect_mod = z_expect_direct % Lz_ang

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(t_zp, z_centroid, 'b-', lw=1.5, label='WP z-centroid (Δn>0)')
ax.plot(t_zp, z_expect_mod, 'k--', lw=1, alpha=0.6, label='Free-particle (+z, PBC)')
ax.axhline(z_screen0_ang, color='purple', lw=1, ls='--',
           label=f'Screen 0: {z_screen0_ang:.2f} Å')
ax.axhline(z_flake_ang, color='orange', lw=1.5, ls='--',
           label=f'z_flake: {z_flake_ang:.2f} Å')
ax.axvline(t1_au, color='r', lw=1, ls='--', label=f'T1={t1_au:.2f} a.u.')
ax.axvline(t_leed_start, color='g', lw=1, ls=':', label=f't_leed')
ax.set_xlabel('t  (a.u.)', fontsize=12)
ax.set_ylabel('z  (Å)', fontsize=12)
ax.set_title('WP z-centroid vs time', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'wp_trajectory_z.png'), dpi=150)
plt.close(fig)
print("  Saved: wp_trajectory_z.png")

# 2D heatmap of z-profile
vmax_prof = np.percentile(np.abs(dprof), 99)
fig, ax = plt.subplots(figsize=(10, 5))
im = ax.pcolormesh(t_zp, z_ang, dprof.T,
                   cmap='RdBu_r', vmin=-vmax_prof, vmax=vmax_prof,
                   shading='auto')
ax.axhline(z_screen0_ang, color='purple', lw=1, ls='--', label='Screen 0')
ax.axhline(z_flake_ang,   color='orange', lw=1.5, ls='--', label='z_flake')
ax.axvline(t1_au,         color='r', lw=1, ls='--', label='T1')
ax.axvline(t_leed_start,  color='g', lw=1, ls=':', label='t_leed_start')
ax.set_xlabel('t  (a.u.)', fontsize=12)
ax.set_ylabel('z  (Å)', fontsize=12)
ax.set_title('Δn(z,t) = n(z,t) − n_baseline', fontsize=12)
ax.legend(fontsize=9, loc='upper right')
plt.colorbar(im, ax=ax, label='Δn  (bohr⁻³)')
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'zprofile_heatmap.png'), dpi=150)
plt.close(fig)
print("  Saved: zprofile_heatmap.png")

# ──────────────────────────────────────────────────────────────
# 9. Density snapshots
# ──────────────────────────────────────────────────────────────
def load_snapshot(path):
    data = []
    with open(path) as f:
        for i, line in enumerate(f):
            if line.startswith('#'):
                continue
            vals = list(map(float, line.split()))
            if vals:
                data.append(vals)
    return np.array(data)

def plot_snapshot_grid(snap_dir_path, title_prefix, out_name, n_cols=5, dt_au=0.02):
    files = sorted([f for f in os.listdir(snap_dir_path) if f.endswith('.txt')])
    n_files = len(files)
    if n_files == 0:
        print(f"  No snapshots in {snap_dir_path}")
        return
    n_show = min(10, n_files)
    indices = np.round(np.linspace(0, n_files-1, n_show)).astype(int)
    selected = [files[i] for i in indices]

    n_rows = (n_show + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows), squeeze=False)

    snaps = [load_snapshot(os.path.join(snap_dir_path, f)) for f in selected]
    vmax = max(np.max(s) for s in snaps)

    for idx, (fname, snap) in enumerate(zip(selected, snaps)):
        ax = axes.flat[idx]
        step_str = fname.replace('snapshot_t', '').replace('.txt', '')
        try:
            step_val = int(step_str)
        except:
            step_val = -1
        t_val = step_val * dt_au
        im = ax.pcolormesh(xc, yc, snap,
                           cmap='inferno', vmin=0, vmax=vmax, shading='auto')
        ax.set_title(f'step={step_val}  t={t_val:.2f} a.u.', fontsize=8)
        ax.set_aspect('equal')

    for ax in list(axes.flat)[n_show:]:
        ax.set_visible(False)

    fig.suptitle(title_prefix, fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, out_name), dpi=120)
    plt.close(fig)
    print(f"  Saved: {out_name}")

print("\n=== Density snapshots ===")
dt_au_val = float(summary.get('dt_au', '0.02'))
plot_snapshot_grid(os.path.join(RES_DIR, 'density_snapshots'),
                   'Density at z_flake (molecular plane)',
                   'density_snapshots_flake.png', dt_au=dt_au_val)
plot_snapshot_grid(os.path.join(RES_DIR, 'density_obs_snapshots'),
                   'Density at z_obs (screen 0)',
                   'density_snapshots_obs.png', dt_au=dt_au_val)
plot_snapshot_grid(os.path.join(RES_DIR, 'density_mid_snapshots'),
                   'Density at z_mid',
                   'density_snapshots_mid.png', dt_au=dt_au_val)

# ──────────────────────────────────────────────────────────────
# 10. Validation summary
# ──────────────────────────────────────────────────────────────
print("\n=== Validation summary ===")
print(f"  GS energy:     {float(summary.get('GS_energy_Ha','0')):.6f} Ha")
print(f"  SCF steps:     {summary.get('SCF_steps','?')}")
print(f"  WP norm:       {float(summary.get('WP_norm','0')):.6f}  (target ≈ 1.0)")
print(f"  Jz(t=0):       {Jz[0]:.4f} a.u.  k₀={k0:.4f}  "
      f"rel error={(abs(Jz[0])-k0)/k0*100:.2f}%")
print(f"  Energy drift:  {drift:.3e} Ha/a.u.")
print(f"  Wall time:     {float(summary.get('wall_time_sec','0'))/3600:.2f} h")
print(f"\nAll figures in: {FIG_DIR}")
print("Analysis complete.")
