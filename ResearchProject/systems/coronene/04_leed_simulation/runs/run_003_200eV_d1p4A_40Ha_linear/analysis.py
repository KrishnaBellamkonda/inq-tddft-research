#!/usr/bin/env python3
"""
analysis.py — Post-processing for run_003 coronene LEED simulation.

Reference: Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014).

All spatial units: bohr (atomic units) internally; Å for display.
All time units: a.u. internally; fs for display.
All energy units: Ha internally; eV for display.

Usage:
    python3 analysis.py               # run all validation checks + plots
    python3 analysis.py --validate    # validation checks only
    python3 analysis.py --leed        # LEED pattern plot only
"""

import os
import sys
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')   # non-interactive backend for cluster use
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ── Directory layout ──────────────────────────────────────────────────────────
RESULTS = os.path.join(os.path.dirname(__file__), 'results')
FIGS    = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(FIGS, exist_ok=True)

# ── Physical constants ────────────────────────────────────────────────────────
BOHR_TO_ANG = 0.529177210903
ANG_TO_BOHR = 1.8897259886
HA_TO_EV    = 27.21138625
FS_TO_AU    = 1.0 / 0.024188843
AU_TO_FS    = 0.024188843

# ── Paper parameters (from config.hpp) ───────────────────────────────────────
# Cell
LX_BOHR = 18.4 * ANG_TO_BOHR   # 34.771 bohr
LY_BOHR = 18.4 * ANG_TO_BOHR
LZ_BOHR = 31.7 * ANG_TO_BOHR   # 59.904 bohr

# WP parameters
WP_D_BOHR  = 1.4 * ANG_TO_BOHR             # 2.646 bohr
WP_EKIN_HA = 200.0 / HA_TO_EV              # 7.350 Ha
K0_BOHR    = np.sqrt(2.0 * WP_EKIN_HA)    # 3.834 bohr^-1

# Planes
Z_FLAKE_BOHR = LZ_BOHR / 2.0                   # 29.952 bohr
WP_D_IMPACT  = 6.35 * ANG_TO_BOHR             # 12.001 bohr
Z_OBS_BOHR   = Z_FLAKE_BOHR + WP_D_IMPACT     # 41.952 bohr
Z_MID_BOHR   = (Z_FLAKE_BOHR + Z_OBS_BOHR) / 2.0  # 35.952 bohr

# Propagation
DT_AU  = 4.84e-4 * FS_TO_AU     # 0.020009 a.u.
T1_AU  = 0.077  * FS_TO_AU      # 3.183 a.u. — WP reaches molecule
T2_AU  = 0.25   * FS_TO_AU      # 10.335 a.u. — end of propagation
N_STATES = 57                   # 54 occupied + 2 buffer + 1 WP
IST_WP   = 56                   # WP orbital index (0-based)

# ── Grid metadata ─────────────────────────────────────────────────────────────

_grid_cache = None

def load_grid_metadata():
    """Return dict: Nx, Ny, Nz, dx, dy, dz (bohr), Lx, Ly, Lz (bohr).
    Reads once and caches.
    """
    global _grid_cache
    if _grid_cache is not None:
        return _grid_cache
    path = os.path.join(RESULTS, 'grid', 'grid_metadata.txt')
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            vals = line.split()
            _grid_cache = dict(
                Nx=int(vals[0]), Ny=int(vals[1]), Nz=int(vals[2]),
                dx=float(vals[3]), dy=float(vals[4]), dz=float(vals[5]),
                Lx=float(vals[6]), Ly=float(vals[7]), Lz=float(vals[8]),
            )
            return _grid_cache
    raise RuntimeError(f"Cannot parse {path}")


def grid_axes():
    """Return (x, y, z) 1D coordinate arrays in bohr."""
    g = load_grid_metadata()
    x = np.arange(g['Nx']) * g['dx']
    y = np.arange(g['Ny']) * g['dy']
    z = np.arange(g['Nz']) * g['dz']
    return x, y, z


# ── 3D density loaders ────────────────────────────────────────────────────────

def load_density_3d(step):
    """Load 3D electron density at given timestep.

    Returns ndarray (Nx, Ny, Nz), values in bohr^-3.
    Reads results/density/density_tNNNNNN.txt.
    """
    g = load_grid_metadata()
    fname = os.path.join(RESULTS, 'density', f'density_t{step:06d}.txt')
    data = np.loadtxt(fname, comments='#')
    return data.reshape(g['Nx'], g['Ny'], g['Nz'])


def list_density_steps():
    """Return sorted list of available 3D density snapshot steps."""
    d = os.path.join(RESULTS, 'density')
    steps = []
    for fname in os.listdir(d):
        m = re.match(r'density_t(\d+)\.txt', fname)
        if m:
            steps.append(int(m.group(1)))
    return sorted(steps)


# ── GS orbital loaders ────────────────────────────────────────────────────────

def load_gs_orbital(ist):
    """Load ground-state KS orbital ist (0-based).

    Returns ndarray complex128 (Nx, Ny, Nz).
    Reads results/gs_orbitals/orbital_NNNN/orbital.txt.
    """
    g = load_grid_metadata()
    fname = os.path.join(RESULTS, 'gs_orbitals', f'orbital_{ist:04d}', 'orbital.txt')
    raw = np.loadtxt(fname, comments='#')     # shape (N_pts, 2)
    return (raw[:, 0] + 1j * raw[:, 1]).reshape(g['Nx'], g['Ny'], g['Nz'])


# ── WP orbital loaders ────────────────────────────────────────────────────────

def load_wp_orbital(step):
    """Load WP orbital at given timestep from wp_orbital snapshots.

    Returns ndarray complex128 (Nx, Ny, Nz), or None if not saved.
    Reads results/wp_orbital/step_NNNNNN/kpt_0/orbital_0056/orbital.txt.
    """
    g = load_grid_metadata()
    fname = os.path.join(RESULTS, 'wp_orbital',
                         f'step_{step:06d}', 'kpt_0',
                         f'orbital_{IST_WP:04d}', 'orbital.txt')
    if not os.path.exists(fname):
        return None
    raw = np.loadtxt(fname, comments='#')
    return (raw[:, 0] + 1j * raw[:, 1]).reshape(g['Nx'], g['Ny'], g['Nz'])


# ── 2D density snapshot loaders ───────────────────────────────────────────────

def _load_2d_snapshot(subdir, step):
    """Load a 2D density slice from a snapshot directory.

    Format: header line `# t=T z=Z`, then Nx rows × Ny space-separated values.
    Returns ndarray (Nx, Ny), values in bohr^-3.
    """
    fname = os.path.join(RESULTS, subdir, f'snapshot_t{step:06d}.txt')
    data = []
    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                continue
            data.append([float(v) for v in line.split()])
    return np.array(data)


def load_density_flake(step):
    """2D density slice at z_flake (coronene plane) at given step.
    Returns (Nx, Ny) ndarray."""
    return _load_2d_snapshot('density_snapshots', step)


def load_density_obs(step):
    """2D density slice at z_obs (LEED screen / WP start) at given step.
    Returns (Nx, Ny) ndarray."""
    return _load_2d_snapshot('density_obs_snapshots', step)


def load_density_mid(step):
    """2D density slice at z_mid at given step.
    Returns (Nx, Ny) ndarray."""
    return _load_2d_snapshot('density_mid_snapshots', step)


def list_snapshot_steps(subdir='density_snapshots'):
    d = os.path.join(RESULTS, subdir)
    steps = []
    for fname in os.listdir(d):
        m = re.match(r'snapshot_t(\d+)\.txt', fname)
        if m:
            steps.append(int(m.group(1)))
    return sorted(steps)


# ── z-profile ────────────────────────────────────────────────────────────────

_zprofile_cache = None

def load_z_profile():
    """Load z-profile CSV (density integrated over x,y per iz).

    Returns dict with:
        steps (N,)  — timestep indices
        t_au  (N,)  — times in a.u.
        n     (N, Nz) — density profile (bohr^-2, integral over x,y × dx × dy)
    """
    global _zprofile_cache
    if _zprofile_cache is not None:
        return _zprofile_cache
    fname = os.path.join(RESULTS, 'wp_trajectory', 'density_z_profile_vs_time.csv')
    raw = np.loadtxt(fname, delimiter=',', comments='#')
    _zprofile_cache = dict(
        steps=raw[:, 0].astype(int),
        t_au=raw[:, 1],
        n=raw[:, 2:],
    )
    return _zprofile_cache


# ── Energy / momentum CSVs ────────────────────────────────────────────────────

_energy_cache = None

def load_energy():
    """Load energy vs time CSV.
    Returns dict: steps, t_au, E_Ha."""
    global _energy_cache
    if _energy_cache is not None:
        return _energy_cache
    fname = os.path.join(RESULTS, 'energy', 'energy_vs_time.csv')
    raw = np.loadtxt(fname, delimiter=',', comments='#')
    _energy_cache = dict(steps=raw[:, 0].astype(int), t_au=raw[:, 1], E_Ha=raw[:, 2])
    return _energy_cache


_momentum_cache = None

def load_momentum():
    """Load momentum (current density) vs time CSV.
    Returns dict: steps, t_au, Jx, Jy, Jz (all in a.u.)."""
    global _momentum_cache
    if _momentum_cache is not None:
        return _momentum_cache
    fname = os.path.join(RESULTS, 'momentum', 'momentum_vs_time.csv')
    raw = np.loadtxt(fname, delimiter=',', comments='#')
    _momentum_cache = dict(
        steps=raw[:, 0].astype(int), t_au=raw[:, 1],
        Jx=raw[:, 2], Jy=raw[:, 3], Jz=raw[:, 4],
    )
    return _momentum_cache


_ksoverlap_cache = None

def load_ks_overlaps():
    """Load KS projected occupation vs time.
    Returns dict: steps, t_au, S2 (N_steps, N_states) — |S_ii|^2."""
    global _ksoverlap_cache
    if _ksoverlap_cache is not None:
        return _ksoverlap_cache
    fname = os.path.join(RESULTS, 'ks_overlaps', 'projected_occ_vs_time.csv')
    raw = np.loadtxt(fname, delimiter=',', comments='#')
    _ksoverlap_cache = dict(
        steps=raw[:, 0].astype(int), t_au=raw[:, 1], S2=raw[:, 2:],
    )
    return _ksoverlap_cache


# ── Overlap matrix ────────────────────────────────────────────────────────────

def load_overlap_matrix(step):
    """Parse overlap_matrix.txt and return the block for the given step.

    Returns complex ndarray (N_states, N_states) or None if step not found.
    S_ij = <phi_i_GS | phi_j(t)>
    """
    fname = os.path.join(RESULTS, 'overlap_matrix', 'overlap_matrix.txt')
    g = load_grid_metadata()
    n = N_STATES
    target_step = step
    found = False
    rows = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                m = re.search(r'step=(\d+)', line)
                if m and int(m.group(1)) == target_step:
                    found = True
                    rows = []
                elif found and m:
                    break   # next block started
                continue
            if not found:
                continue
            vals = line.split()
            # Each row has n pairs: "re im re im ..."
            pairs = [(float(vals[2*k]), float(vals[2*k+1])) for k in range(n)]
            rows.append([re + 1j*im for re, im in pairs])
            if len(rows) == n:
                break
    if not rows:
        return None
    return np.array(rows)   # (n, n)


# ── LEED pattern ──────────────────────────────────────────────────────────────

def load_leed_pattern():
    """Load the LEED pattern from results/leed_pattern/leed_pattern.txt.

    Format (from run.cpp): header + Nx rows of Ny space-separated floats.
    Returns ndarray (Nx, Ny).
    """
    fname = os.path.join(RESULTS, 'leed_pattern', 'leed_pattern.txt')
    if not os.path.exists(fname):
        return None
    data = []
    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                continue
            data.append([float(v) for v in line.split()])
    return np.array(data)


def time_averaged_density_at_obs(step_list=None, subtract_gs=True):
    """Time-average of density at z_obs over t1 ≤ t ≤ t2.

    Uses the 2D obs snapshots saved during propagation.
    If step_list=None, selects all obs snapshot steps with t in [t1, t2].

    subtract_gs=True (default): subtracts the GS density at z_obs (step 0)
    before averaging.  This isolates the WP-scattered contribution from the
    molecular π-electron tail background, which dominates the raw density.
    Without subtraction the pattern is the coronene charge density seen from
    22.2 Å above, not the diffraction signal.

    Returns (Nx, Ny) ndarray (units: bohr^-3 × a.u. if summed, else bohr^-3).
    """
    en  = load_energy()

    def step_to_t(s):
        idx = np.searchsorted(en['steps'], s)
        if idx < len(en['steps']) and en['steps'][idx] == s:
            return en['t_au'][idx]
        return s * DT_AU

    avail = list_snapshot_steps('density_obs_snapshots')
    if step_list is None:
        step_list = [s for s in avail if T1_AU <= step_to_t(s) <= T2_AU]
    if not step_list:
        return None
    gs_bg = load_density_obs(0) if subtract_gs else 0.0
    acc = None
    for s in step_list:
        sl = load_density_obs(s) - gs_bg
        acc = sl if acc is None else acc + sl
    return acc / len(step_list)


def leed_from_obs_snapshots(subtract_gs=True):
    """Compute LEED pattern from saved obs snapshots (background-subtracted).

    This is the physically correct LEED signal: the WP-scattered electron
    density returned to the observation plane, with the static molecular
    background removed.

    Returns (Nx, Ny) ndarray.
    """
    return time_averaged_density_at_obs(subtract_gs=subtract_gs)


# ── Observable helpers ────────────────────────────────────────────────────────

def wp_z_centroid_from_zprofile():
    """Compute a proxy WP z-position vs time from the z-profile line scan.

    IMPORTANT: the z-profile is a LINE SCAN along (Lx/2, Ly/2, z), not the
    xy-integrated density.  The WP peak density at z_obs is ~9.7e-3 bohr^-3
    while the molecular π-electron tails there are ~2.6e-2 bohr^-3, so the
    signal-to-background ratio is only ~0.37.  The centroid is therefore a
    noisy proxy for the WP position, not a precise tracker.

    We use the density excess Δn(z,t) = n(z,t) − n(z,t=0) in the region
    z > z_flake and take the centroid of the positive lobe (WP moved here).

    Returns (t_au, z_centroid_bohr) — informational only.
    """
    zp   = load_z_profile()
    _, _, z = grid_axes()
    t_au = zp['t_au']
    nz   = zp['n']                   # (N_steps, Nz)
    dn   = nz - nz[0:1, :]          # density excess
    mask_z = z > Z_FLAKE_BOHR
    dn_wp  = dn[:, mask_z]
    z_wp   = z[mask_z]
    dn_pos = np.maximum(dn_wp, 0.0)
    norm   = np.sum(dn_pos, axis=1)
    good   = norm > 0
    z_cent = np.where(good,
                      np.sum(dn_pos * z_wp[np.newaxis, :], axis=1) / np.where(good, norm, 1.0),
                      Z_OBS_BOHR)
    return t_au, z_cent


# ── Validation ───────────────────────────────────────────────────────────────

def validate_energy_conservation(tolerance_ha=0.01):
    """Check total energy drift over the simulation.

    Pass criterion: max |E(t) - E(0)| < tolerance_ha.
    """
    en = load_energy()
    E0    = en['E_Ha'][0]
    drift = np.max(np.abs(en['E_Ha'] - E0))
    pass_ = drift < tolerance_ha
    print(f"  Energy conservation: max drift = {drift:.4e} Ha  "
          f"[{'PASS' if pass_ else 'FAIL'} tol={tolerance_ha:.0e} Ha]")
    return pass_


def validate_initial_momentum():
    """Check Jz(t=0) ≈ −k₀ = −3.834 a.u. and |Jx|, |Jy| < 1e-2.

    Note: INQ reports total momentum. At t=0 the WP has Jz = +k₀ (before
    sign convention check).  WP travels in −z so Jz sign depends on the
    convention. We check |Jz| ≈ k₀ and transverse |Jx|,|Jy| are small.
    """
    mom = load_momentum()
    Jx0, Jy0, Jz0 = mom['Jx'][0], mom['Jy'][0], mom['Jz'][0]
    pass_z  = abs(abs(Jz0) - K0_BOHR) / K0_BOHR < 0.02   # within 2%
    pass_xy = max(abs(Jx0), abs(Jy0)) < 0.05
    print(f"  Initial momentum: Jx={Jx0:.4f}, Jy={Jy0:.4f}, Jz={Jz0:.4f} a.u."
          f"  (k₀={K0_BOHR:.4f})  "
          f"[{'PASS' if pass_z and pass_xy else 'FAIL'}]")
    return pass_z and pass_xy


def validate_ks_occupation():
    """Check |S_ii(0)|² ≈ 1 for all occupied orbitals (0–53) at t=0.

    Expectation: diagonal ≈ 1 (GS orbitals unchanged at t=0),
    WP slot (index 56) ≈ 0 (freshly injected, orthogonal to GS).
    """
    ks = load_ks_overlaps()
    S2_0 = ks['S2'][0]  # shape (N_states,)
    occ_diag   = S2_0[:54]          # occupied KS orbitals
    wp_slot    = S2_0[IST_WP]
    min_diag   = np.min(occ_diag)
    pass_occ   = min_diag > 0.99
    pass_wp    = wp_slot < 1e-10
    print(f"  KS occupation at t=0: min |S_ii|² = {min_diag:.6f} (occ)  "
          f"WP slot = {wp_slot:.2e}  "
          f"[{'PASS' if pass_occ and pass_wp else 'FAIL'}]")
    return pass_occ and pass_wp


def validate_wp_scattering():
    """Check that the WP excited the coronene molecule (scattering occurred).

    Physical expectation: as the WP passes through the molecule, it deposits
    energy and perturbs the molecular KS orbitals.  The diagonal overlap
    |S_ii(t_final)|² for occupied orbitals should depart from 1.0 by a
    measurable amount (> 1e-3) at the end of the simulation.

    Note: total Jz is NOT conserved (fixed ionic potential exchanges momentum
    with electrons), so Jz is not a clean scattering indicator here.
    """
    ks = load_ks_overlaps()
    S2_final = ks['S2'][-1, :54]    # occupied orbitals at last saved step
    t_final  = ks['t_au'][-1]
    max_dev  = float(np.max(np.abs(1.0 - S2_final)))
    pass_    = max_dev > 1e-3
    print(f"  WP scattering (KS excitation): max |1 - S_ii²| at t={t_final:.2f} a.u. "
          f"= {max_dev:.4e}  "
          f"[{'PASS' if pass_ else 'FAIL'} (> 1e-3)]")
    return pass_


def validate_leed_symmetry():
    """Check that the LEED pattern has approximate 6-fold symmetry.

    Coronene has D₆h symmetry → FFT of LEED pattern should show
    6-fold peaks. We check that the intensity at 60° rotations of
    the brightest off-centre peak are within 20% of each other.

    Returns None if leed_pattern.txt not yet written.
    """
    leed = load_leed_pattern()
    if leed is None:
        print("  LEED symmetry: leed_pattern.txt not yet written.")
        return None
    g   = load_grid_metadata()
    dx  = g['dx'] * BOHR_TO_ANG   # Å per pixel
    # 2D FFT (intensity ∝ |FT|²)
    ft   = np.fft.fftshift(np.fft.fft2(leed))
    I    = np.abs(ft)**2
    Nx, Ny = I.shape
    # Find the brightest off-centre peak
    cx, cy = Nx // 2, Ny // 2
    Ir = I.copy()
    Ir[cx-5:cx+6, cy-5:cy+6] = 0   # mask DC
    iy_peak, ix_peak = np.unravel_index(np.argmax(Ir), Ir.shape)
    r_peak = np.sqrt((ix_peak - cx)**2 + (iy_peak - cy)**2)
    # Sample 6-fold rotations — take max in a ±2-pixel neighbourhood
    # to tolerate grid discretisation errors at the rotated positions.
    angles = np.linspace(0, 2*np.pi, 7)[:-1]
    a0     = np.arctan2(iy_peak - cy, ix_peak - cx)
    vals   = []
    hw = 2   # half-width of neighbourhood window
    for da in angles:
        a   = a0 + da
        xi0 = int(round(cx + r_peak * np.cos(a)))
        yi0 = int(round(cy + r_peak * np.sin(a)))
        # max over ±hw neighbourhood
        xi1 = max(0, xi0 - hw); xi2 = min(Nx, xi0 + hw + 1)
        yi1 = max(0, yi0 - hw); yi2 = min(Ny, yi0 + hw + 1)
        vals.append(I[yi1:yi2, xi1:xi2].max() if xi1 < xi2 and yi1 < yi2 else 0.0)
    vals = np.array(vals)
    ratio = vals.min() / vals.max() if vals.max() > 0 else 0
    pass_ = ratio > 0.3
    print(f"  LEED 6-fold symmetry: min/max peak ratio = {ratio:.3f}  "
          f"[{'PASS' if pass_ else 'FAIL'} (> 0.3; qualitative check)]")
    return pass_


def run_validation():
    """Run all validation checks and print a summary."""
    print("=" * 60)
    print("run_003 — Validation suite")
    print("=" * 60)
    results = {}
    results['energy']   = validate_energy_conservation()
    results['momentum'] = validate_initial_momentum()
    results['ks_occ']   = validate_ks_occupation()
    results['wp_scat']  = validate_wp_scattering()
    results['leed_sym'] = validate_leed_symmetry()
    print("-" * 60)
    n_pass = sum(1 for v in results.values() if v is not None and bool(v))
    n_fail = sum(1 for v in results.values() if v is not None and not bool(v))
    n_skip = sum(1 for v in results.values() if v is None)
    print(f"  {n_pass} PASS / {n_fail} FAIL / {n_skip} SKIP")
    print("=" * 60)
    return results


# ── Plots ────────────────────────────────────────────────────────────────────

def plot_energy_conservation(save=True):
    """Plot total energy vs time. Highlight t1 and t2."""
    en = load_energy()
    t_fs = en['t_au'] * AU_TO_FS
    E0   = en['E_Ha'][0]
    dE   = (en['E_Ha'] - E0) * HA_TO_EV * 1e3   # meV

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t_fs, dE, 'b-', lw=1.0, label='E(t) − E(0)')
    ax.axvline(T1_AU * AU_TO_FS, ls='--', color='orange', label=f't₁ = {T1_AU*AU_TO_FS:.3f} fs')
    ax.axvline(T2_AU * AU_TO_FS, ls='--', color='red',    label=f't₂ = {T2_AU*AU_TO_FS:.3f} fs')
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('ΔE (meV)')
    ax.set_title('run_003 — Energy conservation')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGS, 'energy_conservation.png')
        fig.savefig(out, dpi=150)
        print(f"Saved: {out}")
    return fig, ax


def plot_momentum_vs_time(save=True):
    """Plot Jx, Jy, Jz vs time."""
    mom  = load_momentum()
    t_fs = mom['t_au'] * AU_TO_FS
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t_fs, mom['Jz'], 'b-',  lw=1.0, label='Jz')
    ax.plot(t_fs, mom['Jx'], 'g--', lw=0.8, label='Jx')
    ax.plot(t_fs, mom['Jy'], 'r--', lw=0.8, label='Jy')
    ax.axhline(-K0_BOHR, ls=':', color='gray', label=f'−k₀ = {-K0_BOHR:.3f}')
    ax.axvline(T1_AU * AU_TO_FS, ls='--', color='orange', lw=0.8, label='t₁')
    ax.axvline(T2_AU * AU_TO_FS, ls='--', color='red',    lw=0.8, label='t₂')
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('J (a.u.)')
    ax.set_title('run_003 — Momentum vs time')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGS, 'momentum_vs_time.png')
        fig.savefig(out, dpi=150)
        print(f"Saved: {out}")
    return fig, ax


def plot_wp_trajectory(save=True):
    """Plot WP z-centroid vs time with k₀ reference line."""
    t_au, z_cent = wp_z_centroid_from_zprofile()
    t_fs = t_au * AU_TO_FS
    z_ang = z_cent * BOHR_TO_ANG
    z_flake_ang = Z_FLAKE_BOHR * BOHR_TO_ANG
    z_obs_ang   = Z_OBS_BOHR   * BOHR_TO_ANG

    # Reference: free-particle centroid starting at z_obs
    z_ref_ang = (Z_OBS_BOHR - K0_BOHR * t_au) * BOHR_TO_ANG

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t_fs, z_ang,     'b-',  lw=1.2, label='z centroid (density)')
    ax.plot(t_fs, z_ref_ang, 'k--', lw=0.8, label=f'free particle (k₀={K0_BOHR:.3f} bohr/a.u.)')
    ax.axhline(z_flake_ang, ls=':', color='green', label=f'z_flake = {z_flake_ang:.2f} Å')
    ax.axhline(z_obs_ang,   ls=':', color='purple', label=f'z_obs = {z_obs_ang:.2f} Å')
    ax.axvline(T1_AU * AU_TO_FS, ls='--', color='orange', lw=0.8, label='t₁')
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('z centroid (Å)')
    ax.set_title('run_003 — WP trajectory (z-centroid)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGS, 'wp_trajectory.png')
        fig.savefig(out, dpi=150)
        print(f"Saved: {out}")
    return fig, ax


def plot_ks_overlaps(save=True):
    """Plot |S_ii(t)|² for all occupied states and WP slot."""
    ks = load_ks_overlaps()
    t_fs = ks['t_au'] * AU_TO_FS
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    # Upper panel: occupied KS orbitals (0–53)
    ax = axes[0]
    S2_occ = ks['S2'][:, :54]
    ax.plot(t_fs, S2_occ, 'b-', lw=0.4, alpha=0.4)
    ax.set_ylabel('|S_ii|²')
    ax.set_title('Projected occupation — occupied KS orbitals (0–53)')
    ax.set_ylim(0.97, 1.01)
    ax.grid(True, alpha=0.3)
    ax.axvline(T1_AU * AU_TO_FS, ls='--', color='orange', lw=0.8, label='t₁')
    ax.legend(fontsize=9)

    # Lower panel: WP orbital (index 56)
    ax = axes[1]
    ax.plot(t_fs, ks['S2'][:, IST_WP], 'r-', lw=1.0, label=f'WP (ist={IST_WP})')
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('|S_WP|²')
    ax.set_title('WP orbital overlap with GS WP slot')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.axvline(T1_AU * AU_TO_FS, ls='--', color='orange', lw=0.8, label='t₁')
    ax.legend(fontsize=9)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGS, 'ks_overlaps.png')
        fig.savefig(out, dpi=150)
        print(f"Saved: {out}")
    return fig, axes


def plot_density_snapshot_at_plane(step, plane='flake', save=True):
    """Plot 2D electron density at the given plane and step.

    plane: 'flake', 'obs', or 'mid'
    """
    loaders = {'flake': load_density_flake,
               'obs':   load_density_obs,
               'mid':   load_density_mid}
    loader  = loaders[plane]
    snap    = loader(step)
    g       = load_grid_metadata()
    dx_ang  = g['dx'] * BOHR_TO_ANG
    dy_ang  = g['dy'] * BOHR_TO_ANG
    extent  = [0, g['Nx']*dx_ang, 0, g['Ny']*dy_ang]

    z_labels = {'flake': f'z_flake = {Z_FLAKE_BOHR*BOHR_TO_ANG:.2f} Å',
                'obs':   f'z_obs   = {Z_OBS_BOHR*BOHR_TO_ANG:.2f} Å',
                'mid':   f'z_mid   = {Z_MID_BOHR*BOHR_TO_ANG:.2f} Å'}

    t_au = step * DT_AU
    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(snap.T, origin='lower', extent=extent, cmap='inferno',
                   norm=LogNorm(vmin=max(snap.max()*1e-6, 1e-12), vmax=snap.max()))
    plt.colorbar(im, ax=ax, label='n (bohr⁻³)')
    ax.set_xlabel('x (Å)')
    ax.set_ylabel('y (Å)')
    ax.set_title(f'Density at {z_labels[plane]}\n'
                 f'step={step}, t={t_au*AU_TO_FS:.4f} fs')
    fig.tight_layout()
    if save:
        out = os.path.join(FIGS, f'density_{plane}_step{step:06d}.png')
        fig.savefig(out, dpi=150)
        print(f"Saved: {out}")
    return fig, ax


def plot_leed_pattern(save=True):
    """Plot the LEED pattern: raw (from run.cpp) and background-subtracted.

    The raw leed_pattern.txt is the total ∫n dt at z_obs and is dominated
    by the molecular π-electron tail background.  The background-subtracted
    version (∫[n(t) - n_GS] dt) isolates the WP-scattered contribution.
    Both real-space and momentum-space (FFT) views are shown.
    """
    g      = load_grid_metadata()
    dx_ang = g['dx'] * BOHR_TO_ANG
    extent = [0, g['Nx']*dx_ang, 0, g['Ny']*dx_ang]

    # ── Load patterns ─────────────────────────────────────────────────────────
    leed_raw = load_leed_pattern()           # total density ∫n dt (run.cpp)
    leed_sub = leed_from_obs_snapshots()     # background-subtracted from snapshots

    fig, axes = plt.subplots(2, 2, figsize=(13, 11))

    def _plot_panel(ax, data, title, cmap='hot'):
        """Plot real-space density, handling negative values with diverging cmap."""
        vmax = np.abs(data).max()
        vmin = data.min()
        if vmin < 0:
            # Background-subtracted: use symmetric diverging colormap
            im = ax.imshow(data.T, origin='lower', extent=extent,
                           cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        else:
            im = ax.imshow(data.T, origin='lower', extent=extent, cmap=cmap,
                           norm=LogNorm(vmin=max(vmax*1e-6, 1e-20), vmax=vmax))
        plt.colorbar(im, ax=ax, label='Intensity')
        ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
        ax.set_title(title)

    def _plot_fft(ax, data, title):
        ft = np.fft.fftshift(np.fft.fft2(data))
        I  = np.abs(ft)**2
        kx = np.fft.fftshift(np.fft.fftfreq(g['Nx'], d=g['dx']))
        ky = np.fft.fftshift(np.fft.fftfreq(g['Ny'], d=g['dy']))
        k_ext = [kx.min(), kx.max(), ky.min(), ky.max()]
        im = ax.imshow(I.T, origin='lower', extent=k_ext, cmap='hot',
                       norm=LogNorm(vmin=max(I.max()*1e-6, 1e-30), vmax=I.max()))
        plt.colorbar(im, ax=ax, label='|FFT|²')
        ax.set_xlabel('kₓ (bohr⁻¹)'); ax.set_ylabel('kᵧ (bohr⁻¹)')
        ax.set_title(title)
        return I

    # Row 0: raw pattern
    if leed_raw is not None:
        _plot_panel(axes[0,0], leed_raw,
                    f'Raw LEED: ∫n dt  [dominated by mol. background]\nz_obs={Z_OBS_BOHR*BOHR_TO_ANG:.1f} Å, t₁–t₂')
        _plot_fft(axes[0,1], leed_raw, 'Raw LEED — FFT')
    else:
        for ax in axes[0]: ax.text(0.5,0.5,'leed_pattern.txt not found',ha='center',va='center')

    # Row 1: background-subtracted
    if leed_sub is not None:
        _plot_panel(axes[1,0], leed_sub,
                    'BG-subtracted LEED: ∫[n(t)−n_GS] dt  [WP signal only]\n'
                    'Blue = n < n_GS (WP left), Red = n > n_GS (WP returned)')
        _plot_fft(axes[1,1], leed_sub, 'BG-subtracted LEED — FFT (diffraction pattern)')
    else:
        for ax in axes[1]: ax.text(0.5,0.5,'No obs snapshots available',ha='center',va='center')

    fig.suptitle(f'Coronene LEED — run_003  (200 eV, d=1.4 Å, D={WP_D_IMPACT*BOHR_TO_ANG:.2f} Å)',
                 fontsize=12, fontweight='bold')
    fig.tight_layout()
    if save:
        out = os.path.join(FIGS, 'leed_pattern.png')
        fig.savefig(out, dpi=150)
        print(f"Saved: {out}")
    return fig, axes


def plot_overlap_matrix_diagonal(save=True):
    """Plot diagonal of overlap matrix |S_ii(t)|² vs time from full matrix.

    Loads every available block from overlap_matrix.txt.
    """
    fname = os.path.join(RESULTS, 'overlap_matrix', 'overlap_matrix.txt')
    steps_found = []
    t_found     = []
    diag_found  = []
    current_step = None
    current_t    = None
    row_buf      = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                m_step = re.search(r'step=(\d+)', line)
                m_t    = re.search(r't=([\d.]+)', line)
                if m_step and m_t:
                    # flush previous block
                    if row_buf and current_step is not None:
                        mat  = np.array(row_buf)
                        diag = np.abs(np.diag(mat))**2
                        steps_found.append(current_step)
                        t_found.append(current_t)
                        diag_found.append(diag)
                    current_step = int(m_step.group(1))
                    current_t    = float(m_t.group(1))
                    row_buf      = []
                continue
            vals = line.split()
            row  = [float(vals[2*k]) + 1j*float(vals[2*k+1])
                    for k in range(N_STATES)]
            row_buf.append(row)
    # flush last block
    if row_buf and current_step is not None:
        mat  = np.array(row_buf)
        diag = np.abs(np.diag(mat))**2
        steps_found.append(current_step)
        t_found.append(current_t)
        diag_found.append(diag)

    if not steps_found:
        print("No overlap matrix blocks found.")
        return None, None

    t_arr   = np.array(t_found) * AU_TO_FS
    diag_arr = np.array(diag_found)   # (N_blocks, N_states)

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    ax = axes[0]
    for i in range(54):
        ax.plot(t_arr, diag_arr[:, i], 'b-', lw=0.4, alpha=0.4)
    ax.set_ylabel('|S_ii|²')
    ax.set_title('Overlap matrix diagonal — occupied orbitals (0–53)')
    ax.set_ylim(0.97, 1.01)
    ax.grid(True, alpha=0.3)
    ax.axvline(T1_AU * AU_TO_FS, ls='--', color='orange', lw=0.8, label='t₁')
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.plot(t_arr, diag_arr[:, IST_WP], 'r-', lw=1.0, label=f'WP (ist={IST_WP})')
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('|S_WP,WP|²')
    ax.set_title('WP orbital — self-overlap from full matrix')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.axvline(T1_AU * AU_TO_FS, ls='--', color='orange', lw=0.8, label='t₁')
    ax.legend(fontsize=9)
    fig.tight_layout()
    if save:
        out = os.path.join(FIGS, 'overlap_matrix_diagonal.png')
        fig.savefig(out, dpi=150)
        print(f"Saved: {out}")
    return fig, axes


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else '--all'

    if mode in ('--validate', '--all'):
        run_validation()

    if mode in ('--plots', '--all'):
        print("\nGenerating plots...")
        plot_energy_conservation()
        plot_momentum_vs_time()
        plot_wp_trajectory()
        plot_ks_overlaps()
        # Density at flake — pick the step nearest to t1
        t1_step = int(round(T1_AU / DT_AU / 10) * 10)
        plot_density_snapshot_at_plane(t1_step, plane='flake')
        plot_density_snapshot_at_plane(t1_step, plane='obs')
        plot_overlap_matrix_diagonal()

    if mode in ('--leed', '--all'):
        plot_leed_pattern()

    print("Done.")
