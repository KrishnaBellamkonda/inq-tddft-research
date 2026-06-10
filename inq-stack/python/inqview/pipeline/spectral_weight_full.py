"""spectral_weight_full.py — full q_⊥-integrated dynamic structure factor.

Unlike the on-axis (q_⊥=0) version in spectral_weight.py, this module
computes the FULL q_⊥-integrated spectral weight and loss function:

  W_∥(q_z, ω) = Σ_{q_x, q_y} |δn_resp(q_x, q_y, q_z, ω)|²

Memory strategy: store δn_resp_q(q_z, q_x, q_y, t) for all frames in a
single buffer (~2.7 GB for 125³ × 344 frames in complex64). Single pass
through the VTI files. Temporal FFT vectorised over (q_x, q_y) per q_z slice.

Full pipeline:
  Stage 1: δn(r,t) → 3D FFT each frame
  Stage 4a: subtract analytical free-WP density (per frame, on grid)
  Stage 4b: δn_resp = δn - δn_WP_free (in q-space)
  Stage 2: temporal FFT with Hann window per (q_x, q_y, q_z)
  Stage 3: W_∥(q_z, ω) = Σ_{q_⊥} |δn_resp(q_⊥, q_z, ω)|²
  Stage 4c: χ = δn_resp / V_ext (where V_ext is appreciable)
  Stage 4d: L(q_z, ω) = -(4π/q²) Im[χ], q_⊥-integrated
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
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(vti_path))
    reader.Update()
    img = reader.GetOutput()
    nx, ny, nz = img.GetDimensions()
    flat = vtk_to_numpy(img.GetPointData().GetArray(0)).astype(np.float64, copy=False)
    return flat.reshape((nz, ny, nx)).transpose(2, 1, 0)


def _free_wp_density_on_grid(N, L, dx, sigma, z_c, x_c, y_c, t):
    s_t2 = sigma**2 + t**2 / (4 * sigma**2)
    norm = (2 * np.pi * s_t2) ** (-1.5)
    coords = np.arange(N) * dx - L / 2
    dz = coords - z_c;  dz -= L * np.round(dz / L)
    dxx = coords - x_c;  dxx -= L * np.round(dxx / L)
    dyy = coords - y_c;  dyy -= L * np.round(dyy / L)
    r2 = dxx[:, None, None]**2 + dyy[None, :, None]**2 + dz[None, None, :]**2
    return norm * np.exp(-r2 / (2 * s_t2))


def _reference_curves(q_ref, omega_p, k_F, v0):
    omega_pl = omega_p * np.sqrt(1 + 3 * q_ref**2 / (5 * k_F**2))
    omega_plus = q_ref**2 / 2 + q_ref * k_F
    omega_minus = np.abs(q_ref**2 / 2 - q_ref * k_F)
    omega_kin = q_ref * v0
    return omega_pl, omega_plus, omega_minus, omega_kin


def _plot_map(q_plot, omega_eV, Z, q_ref, curves_eV, run_name,
              title_extra, cbar_label, out_path,
              omega_max_eV, q_max, use_symlog=False):
    omega_pl_eV, omega_plus_eV, omega_minus_eV, omega_kin_eV = curves_eV
    fig, ax = plt.subplots(figsize=(10, 7))

    Z_safe = np.where(np.isfinite(Z) & (Z != 0), Z, np.nan)
    pos_vals = Z_safe[np.isfinite(Z_safe) & (Z_safe > 0)]

    if use_symlog or len(pos_vals) == 0:
        vabs = np.nanmax(np.abs(Z_safe)) if np.any(np.isfinite(Z_safe)) else 1
        if vabs == 0:
            vabs = 1
        norm = mcolors.SymLogNorm(linthresh=vabs * 1e-3, vmin=-vabs, vmax=vabs)
        cmap = "RdBu_r"
    else:
        vmin = np.percentile(pos_vals, 2)
        vmax = np.percentile(pos_vals, 99.5)
        norm = mcolors.LogNorm(vmin=max(vmin, 1e-30), vmax=vmax)
        cmap = "inferno"

    im = ax.pcolormesh(q_plot, omega_eV, Z.T, shading="auto",
                       norm=norm, cmap=cmap, rasterized=True)
    fig.colorbar(im, ax=ax, label=cbar_label)

    ax.plot(q_ref, omega_pl_eV, "c-", lw=2, label=r"Plasmon $\omega_{pl}(q)$")
    ax.plot(q_ref, omega_plus_eV, "w--", lw=1.2, alpha=0.8, label=r"P-H $\omega_+$")
    ax.plot(q_ref, omega_minus_eV, "w--", lw=1.2, alpha=0.8, label=r"P-H $\omega_-$")
    ax.plot(q_ref, omega_kin_eV, "lime", lw=2, ls=":", label=r"$\omega = q v_0$")

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
        omega_max_eV: float = 15.0,
        q_max_inv_bohr: float = 0.8,
        deconv_threshold: float = 1e-2) -> Path | None:
    """Full q_⊥-integrated spectral weight and loss function.

    Memory: ~3 GB for the main buffer (complex64). Single pass through VTIs.
    """
    results_dir = Path(results_dir)
    rs = _parse_run_summary(results_dir)

    dd_dir = results_dir / "raw" / "vti" / "density_delta"
    if not dd_dir.exists():
        dd_dir = results_dir / "raw" / "vti" / "density_rt_delta"
    if not dd_dir.exists():
        print(f"spectral_weight_full: no density_delta dir in {results_dir}")
        return None

    vti_files = sorted(dd_dir.glob("*.vti"))
    if len(vti_files) < 10:
        print(f"spectral_weight_full: only {len(vti_files)} frames")
        return None

    # --- Timing ---
    dt_au = float(rs.get("dt_au", 0.02))
    N_frames = len(vti_files)
    steps = np.array(sorted(
        int(re.search(r"t(\d+)", f.stem).group(1)) for f in vti_files))
    times = steps * dt_au
    t_rel = times - times[0]
    T_total = t_rel[-1]

    # --- Physics ---
    L = float(re.search(r"([\d.]+)", str(rs.get("cell_bohr", "50"))).group(1))
    dx = float(rs.get("spacing_bohr", 0.40))
    N_grid = int(round(L / dx))
    N_electrons = int(rs.get("n_electrons", 162))

    n0 = N_electrons / L**3
    omega_p = np.sqrt(4 * np.pi * n0)
    k_F = (3 * np.pi**2 * n0) ** (1.0 / 3)
    v_F = k_F

    sigma = float(rs.get("wp_sigma_bohr", 5.0))
    k0_parts = str(rs.get("wp_k0_bohr_inv", "0 0 2.711")).split()
    k0_z = float(k0_parts[-1]) if len(k0_parts) >= 3 else float(k0_parts[0])
    v0 = abs(k0_z)
    vz_sign = 1.0 if k0_z >= 0 else -1.0

    wp_parts = str(rs.get("wp_center_bohr", "0 0 -5")).split()
    x_c0 = float(wp_parts[0]) if len(wp_parts) >= 3 else 0.0
    y_c0 = float(wp_parts[1]) if len(wp_parts) >= 3 else 0.0
    z_c0 = float(wp_parts[2]) if len(wp_parts) >= 3 else float(wp_parts[-1])

    N_qz_pos = N_grid // 2  # positive q_z values (skip DC)
    dq = 2 * np.pi / L

    print(f"spectral_weight_full: {run_name}")
    print(f"  {N_frames} frames, T = {T_total:.2f} a.u.")
    print(f"  L = {L}, dx = {dx}, N = {N_grid}³, N_qz_pos = {N_qz_pos}")
    print(f"  n₀ = {n0:.4e}, ω_p = {omega_p*HA_TO_EV:.2f} eV, k_F = {k_F:.4f}")
    print(f"  σ = {sigma}, v₀ = {v0:.4f}, WP start = ({x_c0},{y_c0},{z_c0})")

    buf_gb = N_qz_pos * N_grid * N_grid * N_frames * 8 / 1e9
    print(f"  Buffer: {N_qz_pos}×{N_grid}×{N_grid}×{N_frames} complex64 = {buf_gb:.1f} GB")

    # =====================================================================
    # Stage 1 + 4a+4b: Load all frames, 3D FFT, subtract free WP,
    # store δn_resp(q_z, q_x, q_y, t) for positive q_z
    # =====================================================================
    # Buffer: [N_qz_pos, N_grid, N_grid, N_frames] complex64
    dn_resp_buf = np.zeros((N_qz_pos, N_grid, N_grid, N_frames),
                           dtype=np.complex64)
    # Also store WP contribution for V_ext computation
    dn_wp_buf = np.zeros((N_qz_pos, N_grid, N_grid, N_frames),
                         dtype=np.complex64)

    for j in range(N_frames):
        if j % 50 == 0:
            print(f"  Frame {j}/{N_frames} ...")
        t = t_rel[j]

        # Load measured δn
        dn_real = _load_vti_as_array(vti_files[j])

        # Analytical free-WP δn
        z_c_t = z_c0 + vz_sign * v0 * t
        wp_t = _free_wp_density_on_grid(N_grid, L, dx, sigma, z_c_t, x_c0, y_c0, t)
        wp_0 = _free_wp_density_on_grid(N_grid, L, dx, sigma, z_c0, x_c0, y_c0, 0.0)
        dn_wp_real = wp_t - wp_0

        # Response = measured - free WP
        dn_resp_real = dn_real - dn_wp_real

        # 3D FFT
        dn_resp_q = np.fft.fftn(dn_resp_real)
        dn_wp_q = np.fft.fftn(dn_wp_real)

        # Store positive q_z half (indices 1 to N//2)
        dn_resp_buf[:, :, :, j] = dn_resp_q[:, :, 1:N_qz_pos+1].transpose(2, 0, 1).astype(np.complex64)
        dn_wp_buf[:, :, :, j] = dn_wp_q[:, :, 1:N_qz_pos+1].transpose(2, 0, 1).astype(np.complex64)

    print(f"  All frames loaded. Computing temporal FFT ...")

    # =====================================================================
    # Stage 2+3: Temporal FFT per q_z slice, then sum |.|² over (q_x, q_y)
    # =====================================================================
    N_omega = N_frames * zero_pad_factor
    hann = np.hanning(N_frames).astype(np.float32)

    q_z_all = np.fft.fftfreq(N_grid, d=dx) * 2 * np.pi
    q_z_pos = q_z_all[1:N_qz_pos+1]
    omega_all = np.fft.fftfreq(N_omega, d=(times[1]-times[0])) * 2 * np.pi
    omega_pos = omega_all[:N_omega // 2]
    omega_eV = omega_pos * HA_TO_EV
    N_omega_pos = N_omega // 2

    # q_x, q_y grids for |q|² computation
    q_x_all = np.fft.fftfreq(N_grid, d=dx) * 2 * np.pi
    q_y_all = np.fft.fftfreq(N_grid, d=dx) * 2 * np.pi
    qx2 = q_x_all[:, None]**2  # (Nx, 1)
    qy2 = q_y_all[None, :]**2  # (1, Ny)

    W_resp = np.zeros((N_qz_pos, N_omega_pos), dtype=np.float64)
    W_raw_total = np.zeros((N_qz_pos, N_omega_pos), dtype=np.float64)
    L_integrated = np.zeros((N_qz_pos, N_omega_pos), dtype=np.float64)
    deconv_count = np.zeros((N_qz_pos, N_omega_pos), dtype=np.int32)

    for k in range(N_qz_pos):
        if k % 10 == 0:
            print(f"  q_z slice {k}/{N_qz_pos} (q_z = {q_z_pos[k]:.3f} Bohr⁻¹) ...")

        qz2 = q_z_pos[k]**2

        # Response: temporal FFT over (q_x, q_y) simultaneously
        # data shape: (Nx, Ny, N_frames) → apply hann → FFT → (Nx, Ny, N_omega)
        data_resp = dn_resp_buf[k, :, :, :] * hann[None, None, :]
        spectra_resp = np.fft.fft(data_resp, n=N_omega, axis=2)[:, :, :N_omega_pos]
        W_resp[k, :] = np.sum(np.abs(spectra_resp)**2, axis=(0, 1))

        # WP contribution for V_ext
        data_wp = dn_wp_buf[k, :, :, :] * hann[None, None, :]
        spectra_wp = np.fft.fft(data_wp, n=N_omega, axis=2)[:, :, :N_omega_pos]

        # Raw total (resp + WP, for comparison)
        data_total = (dn_resp_buf[k] + dn_wp_buf[k]) * hann[None, None, :]
        spectra_total = np.fft.fft(data_total, n=N_omega, axis=2)[:, :, :N_omega_pos]
        W_raw_total[k, :] = np.sum(np.abs(spectra_total)**2, axis=(0, 1))

        # Deconvolution: for each (q_x, q_y), compute χ and L, then sum
        # V_ext(q, ω) = -(4π/|q|²) × n_WP(q, ω)
        # L(q, ω) = -(4π/|q|²) Im[χ(q, ω)]
        # χ(q, ω) = δn_resp(q, ω) / V_ext(q, ω)
        q_sq = qx2 + qy2 + qz2  # (Nx, Ny) — |q|² for this q_z slice
        q_sq_safe = np.where(q_sq > 1e-10, q_sq, 1e-10)
        coulomb = 4 * np.pi / q_sq_safe  # (Nx, Ny)

        # V_ext(q_x, q_y, q_z, ω) = -coulomb(q_x,q_y) × n_WP(q_x,q_y,q_z,ω)
        V_ext = -coulomb[:, :, None] * spectra_wp  # (Nx, Ny, N_omega_pos)
        V_ext_abs = np.abs(V_ext)
        V_ext_max = V_ext_abs.max()
        thresh = deconv_threshold * V_ext_max if V_ext_max > 0 else 1e30

        safe = V_ext_abs > thresh
        chi = np.zeros_like(spectra_resp)
        chi[safe] = spectra_resp[safe] / V_ext[safe]

        # L contribution: -(4π/q²) Im[χ], summed over (q_x, q_y)
        L_per_qperp = -coulomb[:, :, None] * np.imag(chi)
        L_integrated[k, :] = np.sum(L_per_qperp, axis=(0, 1))
        deconv_count[k, :] = np.sum(safe, axis=(0, 1))

    # Free the big buffers
    del dn_resp_buf, dn_wp_buf

    total_deconv = deconv_count.sum()
    total_points = N_qz_pos * N_grid * N_grid * N_omega_pos
    print(f"  Deconvolution: {total_deconv}/{total_points} points "
          f"({100*total_deconv/total_points:.1f}%) above threshold")

    # =====================================================================
    # Plot region
    # =====================================================================
    q_mask = q_z_pos <= q_max_inv_bohr
    omega_mask = omega_eV <= omega_max_eV
    q_plot = q_z_pos[q_mask]
    omega_plot = omega_eV[omega_mask]

    W_resp_p = W_resp[np.ix_(q_mask, omega_mask)]
    W_raw_p = W_raw_total[np.ix_(q_mask, omega_mask)]
    L_p = L_integrated[np.ix_(q_mask, omega_mask)]

    q_ref = np.linspace(0.01, q_max_inv_bohr, 200)
    omega_pl, omega_plus, omega_minus, omega_kin = _reference_curves(
        q_ref, omega_p, k_F, v0)
    curves_eV = (omega_pl * HA_TO_EV, omega_plus * HA_TO_EV,
                 omega_minus * HA_TO_EV, omega_kin * HA_TO_EV)

    out_dir = results_dir / "analysis" / "observables"
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Plot 1: Raw W (includes WP) ---
    _plot_map(q_plot, omega_plot, W_raw_p, q_ref, curves_eV, run_name,
              r"raw $W_\parallel(q_z, \omega)$ — full $q_\perp$ integration (includes WP)",
              r"$\sum_{q_\perp} |\delta n|^2$",
              out_dir / "spectral_weight_full_raw.png",
              omega_max_eV, q_max_inv_bohr)
    print(f"  Saved: spectral_weight_full_raw.png")

    # --- Plot 2: Response W (WP subtracted) ---
    _plot_map(q_plot, omega_plot, W_resp_p, q_ref, curves_eV, run_name,
              r"response $W_\parallel(q_z, \omega)$ — full $q_\perp$ (WP subtracted)",
              r"$\sum_{q_\perp} |\delta n_{resp}|^2$",
              out_dir / "spectral_weight_full_response.png",
              omega_max_eV, q_max_inv_bohr)
    print(f"  Saved: spectral_weight_full_response.png")

    # --- Plot 3: Loss function L ---
    _plot_map(q_plot, omega_plot, L_p, q_ref, curves_eV, run_name,
              r"loss function $L_\parallel(q_z, \omega)$ — full $q_\perp$ integration",
              r"$\sum_{q_\perp} -\frac{4\pi}{q^2} \mathrm{Im}[\chi]$",
              out_dir / "loss_function_full.png",
              omega_max_eV, q_max_inv_bohr, use_symlog=True)
    print(f"  Saved: loss_function_full.png")

    # --- Diagnostics ---
    dw_eV = (omega_all[1] - omega_all[0]) * HA_TO_EV
    print(f"\n  === Diagnostics ===")
    print(f"  ω_p = {omega_p*HA_TO_EV:.2f} eV, dω = {dw_eV:.2f} eV")
    print(f"  Kinematic × plasmon: q = {omega_p/v0:.3f} Bohr⁻¹")
    print(f"  Landau onset: q_c = {omega_p/v_F:.3f} Bohr⁻¹")

    for label, W_arr in [("W_raw", W_raw_p), ("W_resp", W_resp_p)]:
        if W_arr.size > 0 and np.any(W_arr > 0):
            pk = np.unravel_index(np.argmax(W_arr), W_arr.shape)
            print(f"  Peak {label}: q = {q_plot[pk[0]]:.3f}, ω = {omega_plot[pk[1]]:.2f} eV")

    return out_dir / "spectral_weight_full_response.png"
