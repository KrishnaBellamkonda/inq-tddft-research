#!/usr/bin/env python3
"""
analysis.py — Post-processing for run_005 (coronene TDDFT LEED)
Paper parameters: D=6.35 Å, Lz=31.7 Å, t1=0.077 fs, t2=0.25 fs
LEED formula: I(r) = integral_{t1}^{t2} n_total(r,t) dt  (paper Eq. 5)
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
RUN_DIR = os.path.dirname(os.path.abspath(__file__))
RES_DIR = os.path.join(RUN_DIR, 'results')
FIG_DIR = os.path.join(RES_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

BOHR_TO_ANG         = 0.529177
BOHR_INV_TO_ANG_INV = 1.0 / BOHR_TO_ANG

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
k0          = float(summary.get('WP_k0_bohr_inv', '3.834'))
wp_d_bohr   = float(summary.get('WP_d_bohr', '1.0'))
wp_d_ang    = wp_d_bohr * BOHR_TO_ANG
wp_bz_bohr  = float(summary.get('WP_bz_bohr', '0'))
wp_D_bohr   = float(summary.get('WP_D_bohr', '0'))
z_screen_bohr = float(summary.get('z_screen_bohr', '0'))
z_screen_ang  = z_screen_bohr * BOHR_TO_ANG
t1_au         = float(summary.get('t1_au', '0'))
t2_au         = float(summary.get('t2_au', '0'))
z_flake_ang   = float(summary.get('WP_bz_bohr', '0')) * BOHR_TO_ANG - float(summary.get('WP_D_ang', '6.35'))
dt_au_val     = float(summary.get('dt_au', '0.02'))

# z_flake: WP_bz = z_flake + D  → z_flake = WP_bz - D
z_flake_bohr = float(summary.get('z_flake_bohr', str(wp_bz_bohr - wp_D_bohr)))
z_flake_ang  = z_flake_bohr * BOHR_TO_ANG

print(f"\n  z_flake  = {z_flake_ang:.3f} Å")
print(f"  z_screen = {z_screen_ang:.3f} Å  (= z_obs)")
print(f"  t1       = {t1_au:.4f} a.u.  ({t1_au*0.024189:.4f} fs)")
print(f"  t2       = {t2_au:.4f} a.u.  ({t2_au*0.024189:.4f} fs)")
print(f"  k0       = {k0:.4f} bohr⁻¹")
print(f"  d        = {wp_d_ang:.4f} Å")

# Centred spatial coordinates (molecule centroid at origin)
cx_ang = Lx * BOHR_TO_ANG / 2.0
cy_ang = Ly * BOHR_TO_ANG / 2.0
xc = x_ang - cx_ang
yc = y_ang - cy_ang

# ──────────────────────────────────────────────────────────────
# 3. Load LEED pattern and GS baseline
# ──────────────────────────────────────────────────────────────
def load_leed(filename):
    """Load Ny×Nx 2D array from space-separated text with # header lines."""
    data = []
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            vals = list(map(float, line.split()))
            if vals:
                data.append(vals)
    arr = np.array(data)
    print(f"  {os.path.basename(filename)}: shape={arr.shape}  "
          f"min={arr.min():.3e}  max={arr.max():.3e}  sum={arr.sum():.3e}")
    return arr

print("\n=== Loading LEED data ===")
leed_total = load_leed(os.path.join(RES_DIR, 'leed_pattern', 'leed_screen.txt'))
gs_baseline = load_leed(os.path.join(RES_DIR, 'leed_pattern', 'gs_baseline_z_obs.txt'))

# Background-subtracted LEED: I_total - n_GS * (t2 - t1)
dt_leed = t2_au - t1_au
leed_sub = leed_total - gs_baseline * dt_leed
print(f"\n  dt_leed = t2-t1 = {dt_leed:.4f} a.u.")
print(f"  leed_total  sum = {leed_total.sum():.3e}")
print(f"  gs_contrib  sum = {(gs_baseline * dt_leed).sum():.3e}")
print(f"  leed_sub    sum = {leed_sub.sum():.3e}")

# ──────────────────────────────────────────────────────────────
# 4. Real-space LEED — total and background-subtracted (Figure 1)
# ──────────────────────────────────────────────────────────────
print("\n=== Plotting real-space LEED ===")
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# Panel 1: total (paper formula, no subtraction)
vmax0 = leed_total.max()
im0 = axes[0].pcolormesh(xc, yc, leed_total,
                          cmap='hot', vmin=0, vmax=vmax0, shading='auto')
axes[0].set_xlabel('x − x_c  (Å)', fontsize=11)
axes[0].set_ylabel('y − y_c  (Å)', fontsize=11)
axes[0].set_title(f'∫n_total dt  (paper Eq.5)\nz_obs={z_screen_ang:.2f} Å  '
                  f't=[{t1_au:.2f},{t2_au:.2f}] a.u.', fontsize=10)
axes[0].set_aspect('equal')
plt.colorbar(im0, ax=axes[0], label='∫n dt  (bohr⁻³·a.u.)')

# Panel 2: background-subtracted
vmax1 = np.abs(leed_sub).max()
im1 = axes[1].pcolormesh(xc, yc, leed_sub,
                          cmap='RdBu_r', vmin=-vmax1, vmax=vmax1, shading='auto')
axes[1].set_xlabel('x − x_c  (Å)', fontsize=11)
axes[1].set_ylabel('y − y_c  (Å)', fontsize=11)
axes[1].set_title(f'∫(n−n_GS) dt  (BG-subtracted)\nz_obs={z_screen_ang:.2f} Å  '
                  f't=[{t1_au:.2f},{t2_au:.2f}] a.u.', fontsize=10)
axes[1].set_aspect('equal')
plt.colorbar(im1, ax=axes[1], label='∫Δn dt  (bohr⁻³·a.u.)')

fig.suptitle(
    f'Coronene TDDFT LEED — real-space\n'
    f'200 eV | d={wp_d_ang:.3f} Å | D=6.35 Å | Tsubonoya PRB 90, 035416 (2014)',
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
# a = 2.46 Å  |g₁| = 4π/(√3·a) ≈ 2.95 Å⁻¹
a_graph = 2.46
g1 = 4 * np.pi / (np.sqrt(3) * a_graph)
g2 = np.sqrt(3) * g1
g3 = 2.0 * g1
g_refs = [(g1, '--', f'|g₁|={g1:.2f} Å⁻¹'),
          (g2, ':',  f'|g₂|={g2:.2f} Å⁻¹'),
          (g3, '-.', f'|g₃|={g3:.2f} Å⁻¹')]
theta = np.linspace(0, 2*np.pi, 300)

def plot_kspace_single(leed, title, filename, klim):
    F  = np.fft.fftshift(np.fft.fft2(leed))
    Ik = np.abs(F)**2
    Ik_disp = Ik.copy()
    cy, cx = Ny // 2, Nx // 2
    Ik_disp[cy-2:cy+3, cx-2:cx+3] = 0   # suppress DC spike

    fig, ax = plt.subplots(figsize=(6, 5.5))
    im = ax.pcolormesh(KX, KY, Ik_disp,
                       cmap='hot', vmin=0, vmax=Ik_disp.max(),
                       shading='auto')
    ax.set_xlim(-klim, klim)
    ax.set_ylim(-klim, klim)
    ax.set_xlabel('k_x  (Å⁻¹)', fontsize=12)
    ax.set_ylabel('k_y  (Å⁻¹)', fontsize=12)
    ax.set_title(title, fontsize=11)
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='|FFT|²  (arb.)')
    for g_mag, ls, lbl in g_refs:
        ax.plot(g_mag*np.cos(theta), g_mag*np.sin(theta),
                color='cyan', lw=0.8, ls=ls, alpha=0.8, label=lbl)
    ax.legend(fontsize=8, loc='upper right')
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"  Saved: {os.path.basename(filename)}")

print("\n=== Plotting k-space LEED ===")
# FFT of background-subtracted LEED (isolates WP scattering signal)
plot_kspace_single(
    leed_sub,
    f'k-space LEED (BG-subtracted)  klim={15} Å⁻¹\n200 eV | d={wp_d_ang:.3f} Å | D=6.35 Å',
    os.path.join(FIG_DIR, 'leed_kspace.png'), klim=15.0)
plot_kspace_single(
    leed_sub,
    f'k-space LEED (BG-subtracted)  klim={6} Å⁻¹\n200 eV | d={wp_d_ang:.3f} Å | D=6.35 Å',
    os.path.join(FIG_DIR, 'leed_kspace_zoomed.png'), klim=6.0)

# Also FFT of raw total LEED (paper formula, no subtraction)
plot_kspace_single(
    leed_total,
    f'k-space LEED (total, paper Eq.5)  klim={6} Å⁻¹\n200 eV | d={wp_d_ang:.3f} Å | D=6.35 Å',
    os.path.join(FIG_DIR, 'leed_kspace_total_zoomed.png'), klim=6.0)

# ──────────────────────────────────────────────────────────────
# 6. Energy vs time
# ──────────────────────────────────────────────────────────────
print("\n=== Energy conservation ===")
energy_data = np.genfromtxt(
    os.path.join(RES_DIR, 'energy', 'energy_vs_time.csv'),
    delimiter=',', comments='#')
t_e, E_e = energy_data[:, 1], energy_data[:, 2]
dE    = E_e - E_e[0]
drift = (E_e[-1] - E_e[0]) / (t_e[-1] - t_e[0])
print(f"  E(0)   = {E_e[0]:.8f} Ha")
print(f"  E(end) = {E_e[-1]:.8f} Ha")
print(f"  Drift  = {drift:.3e} Ha/a.u.")

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_e, dE * 1000, 'b-', lw=1)
ax.axhline(0, color='k', lw=0.5, ls='--')
ax.axvline(t1_au, color='r', lw=1.2, ls='--', label=f't1={t1_au:.2f} a.u. (WP at flake)')
ax.axvline(t2_au, color='g', lw=1.2, ls='--', label=f't2={t2_au:.2f} a.u. (WP at boundary)')
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
print(f"  Jz(t=0) = {Jz[0]:.4f} a.u.  k₀={k0:.4f}  "
      f"rel err={(abs(Jz[0])-k0)/k0*100:.2f}%")

fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
axes[0].plot(t_m, Jz, 'b-', lw=1, label='Jz')
axes[0].axhline( k0, color='gray', lw=0.8, ls='--', label=f'+k₀={k0:.3f}')
axes[0].axhline(-k0, color='gray', lw=0.8, ls='--', label=f'-k₀')
axes[0].axvline(t1_au, color='r', lw=1.2, ls='--', label=f't1={t1_au:.2f} a.u.')
axes[0].axvline(t2_au, color='g', lw=1.2, ls='--', label=f't2={t2_au:.2f} a.u.')
axes[0].set_ylabel('Jz  (a.u.)', fontsize=11)
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.3)
axes[1].plot(t_m, Jx, 'r-', lw=0.8, label='Jx', alpha=0.7)
axes[1].plot(t_m, Jy, 'g-', lw=0.8, label='Jy', alpha=0.7)
axes[1].axvline(t1_au, color='r', lw=1.2, ls='--')
axes[1].axvline(t2_au, color='g', lw=1.2, ls='--')
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
# 8. WP z-profile heatmap and centroid
# ──────────────────────────────────────────────────────────────
print("\n=== z-profile heatmap and WP centroid ===")
zprof_data = np.genfromtxt(
    os.path.join(RES_DIR, 'wp_trajectory', 'density_z_profile_vs_time.csv'),
    delimiter=',', comments='#')
t_zp    = zprof_data[:, 1]
profile = zprof_data[:, 2:]   # (n_snap, Nz)

# Background: mean of first 3 snapshots (before WP arrives at flake)
n_baseline = np.mean(profile[:3, :], axis=0)
dprof = profile - n_baseline[np.newaxis, :]

# WP z-centroid: centroid of positive Δn
z_centroid = np.array([
    np.sum(z_ang * np.maximum(dp, 0)) / (np.sum(np.maximum(dp, 0)) + 1e-30)
    for dp in dprof
])

# Expected: WP starts at z_screen_ang, moves in -z at k0 (Å/a.u.)
k0_ang = k0 * BOHR_TO_ANG
z_expect = z_screen_ang - k0_ang * t_zp

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(t_zp, z_centroid, 'b-', lw=1.5, label='WP z-centroid (Δn>0)')
ax.plot(t_zp, z_expect,   'k--', lw=1, alpha=0.6, label='Free-particle (−z direction)')
ax.axhline(z_screen_ang, color='purple', lw=1, ls='--',
           label=f'z_obs: {z_screen_ang:.2f} Å')
ax.axhline(z_flake_ang,  color='orange', lw=1.5, ls='--',
           label=f'z_flake: {z_flake_ang:.2f} Å')
ax.axvline(t1_au, color='r', lw=1.2, ls='--', label=f't1={t1_au:.2f} a.u.')
ax.axvline(t2_au, color='g', lw=1.2, ls='--', label=f't2={t2_au:.2f} a.u.')
ax.set_xlabel('t  (a.u.)', fontsize=12)
ax.set_ylabel('z  (Å)', fontsize=12)
ax.set_title('WP z-centroid vs time', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'wp_trajectory_z.png'), dpi=150)
plt.close(fig)
print("  Saved: wp_trajectory_z.png")

# 2D heatmap Δn(z,t)
vmax_prof = np.percentile(np.abs(dprof), 99)
fig, ax = plt.subplots(figsize=(10, 5))
im = ax.pcolormesh(t_zp, z_ang, dprof.T,
                   cmap='RdBu_r', vmin=-vmax_prof, vmax=vmax_prof,
                   shading='auto')
ax.axhline(z_screen_ang, color='purple', lw=1, ls='--', label='z_obs')
ax.axhline(z_flake_ang,  color='orange', lw=1.5, ls='--', label='z_flake')
ax.axvline(t1_au, color='r', lw=1.2, ls='--', label='t1')
ax.axvline(t2_au, color='g', lw=1.2, ls='--', label='t2')
ax.set_xlabel('t  (a.u.)', fontsize=12)
ax.set_ylabel('z  (Å)', fontsize=12)
ax.set_title('Δn(z,t) — z-profile heatmap', fontsize=12)
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
        for line in f:
            if line.startswith('#'):
                continue
            vals = list(map(float, line.split()))
            if vals:
                data.append(vals)
    return np.array(data)

def plot_snapshot_grid(snap_dir_path, title_prefix, out_name, n_cols=5):
    if not os.path.isdir(snap_dir_path):
        print(f"  Skipping {out_name}: directory not found")
        return
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
        t_val = step_val * dt_au_val
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
plot_snapshot_grid(os.path.join(RES_DIR, 'density_snapshots'),
                   'Density at z_flake = 15.85 Å (molecular plane)',
                   'density_snapshots_flake.png')
plot_snapshot_grid(os.path.join(RES_DIR, 'density_obs_snapshots'),
                   'Density at z_obs = 22.20 Å (LEED screen)',
                   'density_snapshots_obs.png')

# ──────────────────────────────────────────────────────────────
# 10. Validation summary
# ──────────────────────────────────────────────────────────────
print("\n=== Validation summary ===")
gs_e  = float(summary.get('GS_energy_Ha', '0'))
norm  = float(summary.get('WP_norm', '0'))
scf_s = summary.get('SCF_steps', '?')
wtime = float(summary.get('wall_time_sec', '0'))

print(f"  GS energy:     {gs_e:.6f} Ha  (run_004: 371.551 Ha — same expected)")
print(f"  SCF steps:     {scf_s}")
print(f"  WP norm:       {norm:.6f}  (target ≈ 1.0)  {'PASS' if abs(norm-1)<0.03 else 'FAIL'}")
print(f"  Jz(t=0):       {Jz[0]:.4f} a.u.  k₀={k0:.4f}  "
      f"rel err={(abs(Jz[0])-k0)/k0*100:.2f}%  "
      f"{'PASS' if abs(abs(Jz[0])-k0)/k0 < 0.05 else 'FAIL'}")
print(f"  Energy drift:  {drift:.3e} Ha/a.u.  {'PASS' if abs(drift) < 0.01 else 'WARN'}")
print(f"  Wall time:     {wtime/3600:.2f} h")
print(f"\nAll figures in: {FIG_DIR}")
print("Analysis complete.")
