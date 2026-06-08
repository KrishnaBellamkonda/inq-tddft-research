"""spectral_weight.py — dynamic structure factor from density perturbation.

Full 5-stage pipeline:
  Stage 1: δn(r,t) → 3D FFT → δn(q,t)
  Stage 2: temporal FFT with Hann window → δn(q,ω)
  Stage 3: W(q_z,ω) = |δn(q_z,ω)|²  (raw spectral weight)
  Stage 4a: Compute free-WP density n_WP(r,t) analytically on the grid
  Stage 4b: δn_resp = δn - δn_WP  (subtract WP's own density)
  Stage 4c: χ(q,ω) = δn_resp(q,ω) / V_ext(q,ω)
  Stage 4d: L(q,ω) = -(4π/q²) Im[χ(q,ω)]

On-axis (q_⊥ = 0) longitudinal response.

Reference curves: plasmon dispersion, P-H boundaries, kinematic line.
All quantities in atomic units internally; axes in Bohr⁻¹ and eV.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

HA_TO_EV = 27.211386245988


def _parse_run_summary(results_dir: Path) -> dict:
    rs_path = results_dir / "run_summary.txt"
    if not rs_path.exists():
        return {}
    d = {}
    for line in rs_path.read_text().splitlines():
        m = re.match(r"\s*(\w+)\s*=\s*(.+)", line)
        if m:
            key, val = m.group(1).strip(), m.group(2).strip()
            try:
                d[key] = float(val)
            except ValueError:
                d[key] = val
    return d


def _load_vti_as_array(vti_path: Path) -> np.ndarray:
    """Load a VTI file and return the scalar field as a 3D numpy array (nx,ny,nz)."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(vti_path))
    reader.Update()
    img = reader.GetOutput()
    nx, ny, nz = img.GetDimensions()
    flat = vtk_to_numpy(img.GetPointData().GetArray(0)).astype(np.float64, copy=False)
    return flat.reshape((nz, ny, nx)).transpose(2, 1, 0)


def _free_wp_density_on_grid(N: int, L: float, dx: float,
                              sigma: float, z_c: float,
                              x_c: float, y_c: float,
                              t: float) -> np.ndarray:
    """Compute the free-WP density n_WP(r, t) on a cubic N³ grid.

    Uses the analytical Gaussian spreading formula:
      n_WP(r,t) = (2π s_t²)^(-3/2) exp(-|r - r_c(t)|² / (2 s_t²))
    where s_t² = σ² + t²/(4σ²) and r_c(t) is the centroid.
    Periodic minimum-image convention applied.
    """
    s_t2 = sigma**2 + t**2 / (4 * sigma**2)
    norm = (2 * np.pi * s_t2) ** (-1.5)

    coords = np.arange(N) * dx - L / 2
    dz = coords - z_c
    dx_arr = coords - x_c
    dy_arr = coords - y_c

    # Minimum image convention for periodic boundaries
    dz = dz - L * np.round(dz / L)
    dx_arr = dx_arr - L * np.round(dx_arr / L)
    dy_arr = dy_arr - L * np.round(dy_arr / L)

    r2 = (dx_arr[:, None, None]**2 +
          dy_arr[None, :, None]**2 +
          dz[None, None, :]**2)

    return norm * np.exp(-r2 / (2 * s_t2))


def _reference_curves(q_ref, omega_p, k_F, v0):
    """Compute plasmon, P-H, and kinematic reference curves."""
    omega_pl = omega_p * np.sqrt(1 + 3 * q_ref**2 / (5 * k_F**2))
    omega_plus = q_ref**2 / 2 + q_ref * k_F
    omega_minus = np.abs(q_ref**2 / 2 - q_ref * k_F)
    omega_kin = q_ref * v0
    return omega_pl, omega_plus, omega_minus, omega_kin


def _plot_map(q_plot, omega_eV, Z, q_ref, curves_eV, run_name,
              title_extra, cbar_label, out_path,
              omega_max_eV, q_max, use_symlog=False):
    """Plot a (q, ω) colourmap with reference curves."""
    omega_pl_eV, omega_plus_eV, omega_minus_eV, omega_kin_eV = curves_eV

    fig, ax = plt.subplots(figsize=(10, 7))

    Z_safe = np.where(np.isfinite(Z) & (Z != 0), Z, np.nan)
    pos_vals = Z_safe[np.isfinite(Z_safe) & (Z_safe > 0)]

    if use_symlog or len(pos_vals) == 0:
        vabs = np.nanmax(np.abs(Z_safe)) if np.any(np.isfinite(Z_safe)) else 1
        norm = mcolors.SymLogNorm(linthresh=vabs * 1e-3, vmin=-vabs, vmax=vabs)
        cmap = "RdBu_r"
    else:
        vmin = np.percentile(pos_vals, 5)
        vmax = np.percentile(pos_vals, 99)
        norm = mcolors.LogNorm(vmin=max(vmin, 1e-30), vmax=vmax)
        cmap = "inferno"

    im = ax.pcolormesh(q_plot, omega_eV, Z.T, shading="auto",
                       norm=norm, cmap=cmap, rasterized=True)
    fig.colorbar(im, ax=ax, label=cbar_label)

    ax.plot(q_ref, omega_pl_eV, "c-", lw=1.5, label=r"Plasmon $\omega_{pl}(q)$")
    ax.plot(q_ref, omega_plus_eV, "w--", lw=1.0, alpha=0.7, label=r"P-H $\omega_+$")
    ax.plot(q_ref, omega_minus_eV, "w--", lw=1.0, alpha=0.7, label=r"P-H $\omega_-$")
    ax.plot(q_ref, omega_kin_eV, "lime", lw=1.5, ls=":",
            label=r"$\omega = q v_0$")

    ax.set_xlabel(r"$q_z$ / Bohr$^{-1}$", fontsize=13)
    ax.set_ylabel(r"$\omega$ / eV", fontsize=13)
    ax.set_title(f"{run_name}: {title_extra}", fontsize=11)
    ax.set_xlim(0, q_max)
    ax.set_ylim(0, omega_max_eV)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def run(results_dir: Path | str, run_name: str = "",
        zero_pad_factor: int = 16,
        eta_factor: float = 3.0,
        omega_max_eV: float = 20.0,
        q_max_inv_bohr: float = 1.0,
        deconv_threshold: float = 1e-3) -> Path | None:
    """Compute spectral weight W, response W_resp, and loss function L.

    Parameters
    ----------
    results_dir : path to the run's results/ directory
    run_name : label for plot titles
    zero_pad_factor : temporal FFT zero-padding (default 16)
    eta_factor : damping η = eta_factor / T (default 3.0)
    omega_max_eV : upper ω for plots (default 20)
    q_max_inv_bohr : upper q for plots (default 1.0)
    deconv_threshold : |V_ext| threshold for safe division (default 1e-3)
    """
    results_dir = Path(results_dir)
    rs = _parse_run_summary(results_dir)

    dd_dir = results_dir / "raw" / "vti" / "density_delta"
    if not dd_dir.exists():
        dd_dir = results_dir / "raw" / "vti" / "density_rt_delta"
    if not dd_dir.exists():
        print(f"spectral_weight: no density_delta VTI dir in {results_dir}")
        return None

    vti_files = sorted(dd_dir.glob("*.vti"))
    if len(vti_files) < 10:
        print(f"spectral_weight: only {len(vti_files)} frames — need ≥10")
        return None

    # --- Timing ---
    dt_au = float(rs.get("dt_au", 0.02))
    write_every = int(rs.get("write_every", 2))
    N_frames = len(vti_files)

    steps = []
    for f in vti_files:
        m = re.search(r"t(\d+)", f.stem)
        if m:
            steps.append(int(m.group(1)))
    steps = np.array(sorted(steps))
    times = steps * dt_au
    t_rel = times - times[0]
    T_total = t_rel[-1]
    eta = eta_factor / T_total

    # --- Physical parameters ---
    L_str = str(rs.get("cell_bohr", "50"))
    L = float(re.search(r"([\d.]+)", L_str).group(1))
    dx = float(rs.get("spacing_bohr", 0.40))
    N_grid = int(round(L / dx))
    N_electrons = int(rs.get("n_electrons", 162))

    n0 = N_electrons / L**3
    omega_p = np.sqrt(4 * np.pi * n0)
    k_F = (3 * np.pi**2 * n0) ** (1.0 / 3)
    v_F = k_F

    sigma = float(rs.get("wp_sigma_bohr", 5.0))
    k0_str = str(rs.get("wp_k0_bohr_inv", "0 0 2.711"))
    k0_parts = k0_str.split()
    k0_z = float(k0_parts[-1]) if len(k0_parts) >= 3 else float(k0_parts[0])
    v0 = abs(k0_z)

    wp_center_str = str(rs.get("wp_center_bohr", "0 0 -21"))
    wp_parts = wp_center_str.split()
    x_c0 = float(wp_parts[0]) if len(wp_parts) >= 3 else 0.0
    y_c0 = float(wp_parts[1]) if len(wp_parts) >= 3 else 0.0
    z_c0 = float(wp_parts[2]) if len(wp_parts) >= 3 else float(wp_parts[-1])

    # WP direction: positive k0_z means +z travel
    vz_sign = 1.0 if k0_z >= 0 else -1.0

    print(f"spectral_weight: {run_name}")
    print(f"  {N_frames} frames, T = {T_total:.2f} a.u., η = {eta:.4f} a.u.⁻¹")
    print(f"  L = {L}, dx = {dx}, N = {N_grid}³")
    print(f"  n₀ = {n0:.4e}, ω_p = {omega_p*HA_TO_EV:.2f} eV, k_F = {k_F:.4f}")
    print(f"  σ = {sigma} Bohr, v₀ = {v0:.4f}, WP start = ({x_c0},{y_c0},{z_c0})")

    # =====================================================================
    # Stage 1: Load δn VTI frames, 3D FFT, extract on-axis q_z profiles
    # Stage 4a: Simultaneously compute free-WP analytical δn on same grid
    # =====================================================================
    N_qz = N_grid
    dn_total_qz_t = np.zeros((N_qz, N_frames), dtype=np.complex128)
    dn_wp_qz_t = np.zeros((N_qz, N_frames), dtype=np.complex128)

    for j, vti_path in enumerate(vti_files):
        if j % 50 == 0:
            print(f"  Frame {j}/{N_frames} ...")
        t = t_rel[j]

        # Measured δn from simulation
        dn_real = _load_vti_as_array(vti_path)
        if j == 0 and dn_real.shape[0] != N_grid:
            N_grid = dn_real.shape[0]
            N_qz = N_grid
            dn_total_qz_t = np.zeros((N_qz, N_frames), dtype=np.complex128)
            dn_wp_qz_t = np.zeros((N_qz, N_frames), dtype=np.complex128)

        dn_q = np.fft.fftn(dn_real)
        dn_total_qz_t[:, j] = dn_q[0, 0, :]

        # Analytical free-WP contribution: n_WP(t) - n_WP(0)
        z_c_t = z_c0 + vz_sign * v0 * t
        wp_t = _free_wp_density_on_grid(N_grid, L, dx, sigma, z_c_t, x_c0, y_c0, t)
        wp_0 = _free_wp_density_on_grid(N_grid, L, dx, sigma, z_c0, x_c0, y_c0, 0.0)
        dn_wp = wp_t - wp_0

        dn_wp_q = np.fft.fftn(dn_wp)
        dn_wp_qz_t[:, j] = dn_wp_q[0, 0, :]

    # =====================================================================
    # Stage 4b: Response density = total - free WP
    # =====================================================================
    dn_resp_qz_t = dn_total_qz_t - dn_wp_qz_t

    print(f"  WP subtraction: max|δn_WP|/max|δn_total| = "
          f"{np.max(np.abs(dn_wp_qz_t)) / max(np.max(np.abs(dn_total_qz_t)), 1e-30):.2f}")

    # =====================================================================
    # Stage 2: Temporal FFT with Hann window
    # =====================================================================
    N_omega = N_frames * zero_pad_factor
    hann = np.hanning(N_frames)

    def temporal_fft(data_qz_t):
        out = np.zeros((N_qz, N_omega), dtype=np.complex128)
        for iq in range(N_qz):
            sig = data_qz_t[iq, :].copy()
            sig -= np.mean(sig)
            out[iq, :] = np.fft.fft(sig * hann, n=N_omega)
        return out

    dn_total_qz_w = temporal_fft(dn_total_qz_t)
    dn_resp_qz_w = temporal_fft(dn_resp_qz_t)

    # Also FFT the WP itself for V_ext computation
    dn_wp_qz_w = temporal_fft(dn_wp_qz_t)

    # =====================================================================
    # Stage 3: Spectral weight maps
    # =====================================================================
    W_total = np.abs(dn_total_qz_w) ** 2
    W_resp = np.abs(dn_resp_qz_w) ** 2

    # =====================================================================
    # Stage 4c-d: Deconvolution → χ(q,ω) → L(q,ω)
    # =====================================================================
    # V_ext(q,ω) = -(4π/q²) n_WP(q,ω)
    # But we need n_WP(q,ω) not δn_WP(q,ω). Since δn_WP = n_WP(t) - n_WP(0),
    # and the DC-subtracted temporal FT of n_WP(0) is zero, the DC-subtracted
    # FT of δn_WP equals the DC-subtracted FT of n_WP. So we can use dn_wp_qz_w.
    q_z_all = np.fft.fftfreq(N_qz, d=dx) * 2 * np.pi
    q_z_col = q_z_all[:, None]
    q2 = q_z_col**2
    q2_safe = np.where(q2 > 1e-10, q2, 1e-10)

    V_ext_qz_w = -(4 * np.pi / q2_safe) * dn_wp_qz_w

    # The q=0 (uniform) mode carries no momentum transfer and L(q,ω)=-(4π/q²)Im[χ]
    # is undefined there; its 4π/q² is only finite because q² was clamped to 1e-10,
    # which inflates |V_ext| at q=0 by ~10³× and poisons the relative threshold below
    # (every physical mode then fails |V_ext|>thresh → χ≡0 → L≡0). Zero the q=0 row so
    # the threshold is set by physical modes and χ/L stay 0 at q=0 (excluded from plot).
    V_ext_qz_w[np.abs(q_z_all) < 1e-5, :] = 0.0
    V_ext_max = np.max(np.abs(V_ext_qz_w))
    thresh = deconv_threshold * V_ext_max

    chi_qz_w = np.zeros_like(dn_resp_qz_w)
    safe_mask = np.abs(V_ext_qz_w) > thresh
    chi_qz_w[safe_mask] = dn_resp_qz_w[safe_mask] / V_ext_qz_w[safe_mask]

    # L(q,ω) = -(4π/q²) Im[χ]
    L_qz_w = -(4 * np.pi / q2_safe) * np.imag(chi_qz_w)

    safe_frac = safe_mask.sum() / safe_mask.size
    print(f"  Deconvolution: {safe_frac*100:.1f}% of (q,ω) above threshold")

    # =====================================================================
    # Build axes and clip to plot region
    # =====================================================================
    omega_all = np.fft.fftfreq(N_omega, d=(times[1] - times[0])) * 2 * np.pi
    q_z_pos = q_z_all[1:N_qz // 2]  # skip q=0 (divergent)
    omega_pos = omega_all[:N_omega // 2]
    omega_eV = omega_pos * HA_TO_EV

    q_mask = q_z_pos <= q_max_inv_bohr
    omega_mask = omega_eV <= omega_max_eV
    q_plot = q_z_pos[q_mask]
    omega_plot = omega_eV[omega_mask]

    def clip(arr):
        return arr[1:N_qz // 2, :N_omega // 2][np.ix_(q_mask, omega_mask)]

    W_total_p = clip(W_total)
    W_resp_p = clip(W_resp)
    L_p = clip(L_qz_w)

    # =====================================================================
    # Reference curves
    # =====================================================================
    q_ref = np.linspace(0.01, q_max_inv_bohr, 200)
    omega_pl, omega_plus, omega_minus, omega_kin = _reference_curves(
        q_ref, omega_p, k_F, v0)
    curves_eV = (omega_pl * HA_TO_EV, omega_plus * HA_TO_EV,
                 omega_minus * HA_TO_EV, omega_kin * HA_TO_EV)

    out_dir = results_dir / "analysis" / "observables"
    out_dir.mkdir(parents=True, exist_ok=True)

    # =====================================================================
    # Plot 1: Raw W(q,ω) — includes WP's own signal
    # =====================================================================
    _plot_map(q_plot, omega_plot, W_total_p, q_ref, curves_eV, run_name,
              r"raw $W(q_z, \omega) = |\delta n_{total}|^2$ (includes WP)",
              r"$|\delta n|^2$",
              out_dir / "spectral_weight_raw.png",
              omega_max_eV, q_max_inv_bohr)
    print(f"  Saved: spectral_weight_raw.png")

    # =====================================================================
    # Plot 2: Response W_resp — WP subtracted
    # =====================================================================
    _plot_map(q_plot, omega_plot, W_resp_p, q_ref, curves_eV, run_name,
              r"response $W_{resp}(q_z, \omega) = |\delta n_{resp}|^2$ (WP subtracted)",
              r"$|\delta n_{resp}|^2$",
              out_dir / "spectral_weight_response.png",
              omega_max_eV, q_max_inv_bohr)
    print(f"  Saved: spectral_weight_response.png")

    # =====================================================================
    # Plot 3: Loss function L(q,ω)
    # =====================================================================
    _plot_map(q_plot, omega_plot, L_p, q_ref, curves_eV, run_name,
              r"loss function $L(q_z, \omega) = -\frac{4\pi}{q^2}\mathrm{Im}[\chi]$",
              r"$L(q, \omega)$",
              out_dir / "loss_function_qz_omega.png",
              omega_max_eV, q_max_inv_bohr, use_symlog=True)
    print(f"  Saved: loss_function_qz_omega.png")

    # =====================================================================
    # Diagnostics
    # =====================================================================
    dw_eV = (omega_all[1] - omega_all[0]) * HA_TO_EV
    print(f"\n  === Diagnostics ===")
    print(f"  ω_p = {omega_p*HA_TO_EV:.2f} eV, dω = {dw_eV:.2f} eV")
    print(f"  Kinematic × plasmon crossing: q = {omega_p/v0:.3f} Bohr⁻¹")
    print(f"  Landau onset: q_c = {omega_p/v_F:.3f} Bohr⁻¹")

    for label, W_arr in [("W_total", W_total_p), ("W_resp", W_resp_p)]:
        if W_arr.size > 0 and np.any(W_arr > 0):
            pk = np.unravel_index(np.argmax(W_arr), W_arr.shape)
            print(f"  Peak {label}: q = {q_plot[pk[0]]:.3f}, ω = {omega_plot[pk[1]]:.2f} eV")

    return out_dir / "spectral_weight_response.png"
