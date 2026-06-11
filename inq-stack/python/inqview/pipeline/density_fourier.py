"""density_fourier — Fourier components n_q(t) of the time-dependent density.

n_q - density related to the momentum transfer q 

For plasmon detection: read every VTI in `raw/vti/density_rt_total/`,
compute delta_n(r, t) = n(r, t) - n(r, 0), 3D-FFT, and pick out the
axial Fourier components n_q_m(t) with q_m = (0, 0, 2*pi*m/L_z). Then
1D-FFT each n_q_m(t) over time to localise the plasmon resonance.


Usage:
    from inqview.pipeline.density_fourier import run
    run(results_dir, run_name="...", m_max=6, dt_au=0.02)

Outputs:
    analysis/observables/n_q_vs_time.csv          long format
    analysis/observables/n_q_vs_time.png          time-domain |n_q_m(t)|
    analysis/observables/n_q_spectrum.csv         long format (omega, m, |FFT|)
    analysis/observables/n_q_spectrum.png         frequency-domain
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

# TODO: Shouldn't the loss function calculation be in 3D and not be restricted to the 
# z axis where the propagation happens? Need to research this. I understand that we
# can split the momentum transfer q into transverse and longitudinal oscillations. 
# I need to check this thoroughly. 

# TODO: Are these the axial modes that ar eploted in the loss function, or the 1D loss function
# for a given q? what is the equation for this? 

# Reuse the VTI loader from the density phase (avoids duplicating the
# vtkXMLImageDataReader plumbing).
from .density import _load_vti_array

Ha2eV = 27.211386245988

_TIME_RE = re.compile(r"_t(\d{6})\.vti$")


def _step_from_filename(p: Path) -> int:
    m = _TIME_RE.search(p.name)
    return int(m.group(1)) if m else -1


def _bohm_gross_omega(q: float, omega_p: float, vF: float) -> float:
    """omega(q)^2 = omega_p^2 + (3/5) vF^2 q^2 + q^4 / 4   (atomic units)."""
    return float(np.sqrt(omega_p**2 + (3.0 / 5.0) * vF**2 * q**2 + q**4 / 4.0))


def run(
    results_dir: Path,
    run_name: str,
    m_max: int = 6,
    dt_au: float = 0.02,
    omega_p_au: float | None = None,
    vF_au: float | None = None,
) -> dict:
    """Extract n_q_m(t) for m = 1..m_max axial modes and FFT over time.

    Parameters
    ----------
    results_dir : `<run>/results/` directory.
    run_name    : run name (for plot titles).
    m_max       : highest axial mode to extract (default 6).
    dt_au       : time step a.u. (default 0.02 — matches Plasmon_N162_L50_E15).
    omega_p_au, vF_au : optional overlays on the FFT plot for the Bohm-Gross
                  predictions; default to the L=50, N=162 values
                  (omega_p = sqrt(4*pi*162/50^3) = 0.1276; vF = kF = 0.337).
    """
    results_dir = Path(results_dir)
    vti_dir = results_dir / "raw" / "vti" / "density_rt_total"
    files = sorted(vti_dir.glob("density_t*.vti"), key=_step_from_filename)
    if not files:
        raise FileNotFoundError(f"No density VTI files in {vti_dir}")

    # First pass: load t=0 to get reference n(r, 0) and grid metadata.
    n0, meta = _load_vti_array(files[0])
    Lx = meta["nx"] * meta["spacing"][0]
    Ly = meta["ny"] * meta["spacing"][1]
    Lz = meta["nz"] * meta["spacing"][2]
    print(f"  Grid: {meta['nx']}x{meta['ny']}x{meta['nz']}  spacing={meta['spacing']}")
    print(f"  Cell: {Lx:.2f} x {Ly:.2f} x {Lz:.2f} bohr")

    # Defaults for L=50, N=162
    if omega_p_au is None:
        n_e_per_vol = 162.0 / (50.0**3)
        omega_p_au = float(np.sqrt(4 * np.pi * n_e_per_vol))
    if vF_au is None:
        vF_au = float((3 * np.pi**2 * 162.0 / (50.0**3))**(1/3))
    print(f"  omega_p = {omega_p_au:.4f} a.u. = {omega_p_au*Ha2eV:.3f} eV")
    print(f"  v_F     = {vF_au:.4f} a.u.")

    # Pre-allocate n_q_m(t) for axial m=1..m_max (complex).
    nframes = len(files)
    n_q = np.zeros((nframes, m_max), dtype=np.complex128)
    times = np.zeros(nframes)

    print(f"  Reading {nframes} VTI frames + FFT-3D each...")
    for i, fp in enumerate(files):
        if i == 0:
            cube = n0
        else:
            cube, _ = _load_vti_array(fp)
        step = _step_from_filename(fp)
        times[i] = step * dt_au
        # delta_n(r, t) (the m=(0,0,0) DC term drops automatically when we look
        # at any nonzero index, but we keep this subtraction for clarity).
        dn = cube - n0
        # 3D FFT (numpy convention: F[k1,k2,k3] = sum dn(r) exp(-i 2pi k.n / N))
        F = np.fft.fftn(dn) * (meta["spacing"][0] * meta["spacing"][1] *
                                meta["spacing"][2])
        # Axial mode index (0, 0, m) → q_m = 2*pi*m/L_z
        for m in range(1, m_max + 1):
            n_q[i, m - 1] = F[0, 0, m]
        if (i + 1) % 50 == 0 or i == nframes - 1:
            print(f"    frame {i+1}/{nframes}  t={times[i]:.1f} a.u.")

    # Wavevectors and predicted Bohm-Gross omega(q_m)
    q_vals = np.array([2 * np.pi * m / Lz for m in range(1, m_max + 1)])
    omega_pred = np.array([_bohm_gross_omega(q, omega_p_au, vF_au) for q in q_vals])

    # Write n_q(t) CSV
    out_dir = results_dir / "analysis" / "observables"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "n_q_vs_time.csv"
    with open(csv_path, "w") as fh:
        fh.write("time_au,m,q_au,re_n_q,im_n_q,abs_n_q\n")
        for i in range(nframes):
            for m in range(1, m_max + 1):
                z = n_q[i, m - 1]
                fh.write(f"{times[i]:.6f},{m},{q_vals[m-1]:.6f},"
                         f"{z.real:.6e},{z.imag:.6e},{abs(z):.6e}\n")
    print(f"  Wrote {csv_path}")

    # Plots: time-domain
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))
    cmap = plt.get_cmap("tab10")
    for m in range(1, m_max + 1):
        ax.plot(times, np.abs(n_q[:, m - 1]), color=cmap(m - 1), lw=0.8,
                label=f"m={m} (q={q_vals[m-1]:.3f}, "
                      f"$\\hbar\\omega_p^{{BG}}$={omega_pred[m-1]*Ha2eV:.2f} eV)")
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"$|n_{q_m}(t)|$  (e $\cdot$ bohr$^{-3}$ $\cdot$ bohr$^{3}$)")
    ax.set_title(f"{run_name} — axial Fourier components of $\\delta n(r,t)$")
    ax.legend(fontsize=8, ncol=2, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "n_q_vs_time.png", dpi=140)
    plt.close(fig)
    print(f"  Wrote {out_dir / 'n_q_vs_time.png'}")

    # FFT each n_q_m(t) over the time series. Detrend by subtracting initial
    # value (= 0 since dn(r, 0) = 0 by construction, so n_q(0) = 0). Window
    # with Hann to suppress endpoint leakage. Skip the first 5 a.u. as
    # injection-shake-up transient (per docs/observables_reference.md §13.6).
    t_start_au = 5.0
    mask = times >= t_start_au
    t_used = times[mask]
    nused = mask.sum()
    print(f"  FFT input: {nused} samples after t_start={t_start_au} a.u. cut")

    # Frequency axis: dt_sample = WRITE_EVERY * dt = (200) * 0.02 = 4 a.u.
    dt_sample = times[1] - times[0]
    print(f"  dt_sample = {dt_sample:.3f} a.u.")
    win = np.hanning(nused)
    n_pad = nused * 4  # zero-pad x4 for smoother spectrum
    spectra = np.zeros((n_pad // 2 + 1, m_max))
    for m in range(1, m_max + 1):
        sig = n_q[mask, m - 1] * win
        # We can FFT real and imaginary parts separately or use complex FFT.
        # Plasmon n_q(t) oscillates as a complex phasor; complex FFT gives
        # the right spectrum.
        full = np.fft.fft(sig.real, n=n_pad)
        spectra[:, m - 1] = np.abs(full[:n_pad // 2 + 1])

    freq_au = np.fft.rfftfreq(n_pad, d=dt_sample)
    omega_au = 2 * np.pi * freq_au
    omega_eV = omega_au * Ha2eV
    print(f"  Frequency grid: 0 to {omega_eV[-1]:.2f} eV in {len(omega_eV)} bins, dE={omega_eV[1]:.4f} eV")

    # Write spectrum CSV
    spec_csv = out_dir / "n_q_spectrum.csv"
    with open(spec_csv, "w") as fh:
        fh.write("omega_au,omega_eV,m,abs_FFT_n_q\n")
        for k, om in enumerate(omega_au):
            for m in range(1, m_max + 1):
                fh.write(f"{om:.6f},{omega_eV[k]:.6f},{m},{spectra[k, m-1]:.6e}\n")
    print(f"  Wrote {spec_csv}")

    # Spectrum plot
    fig, ax = plt.subplots(figsize=(10, 5))
    energy_window = (omega_eV >= 0.1) & (omega_eV <= 15.0)
    for m in range(1, m_max + 1):
        ax.plot(omega_eV[energy_window], spectra[energy_window, m - 1],
                color=cmap(m - 1), lw=1.0,
                label=f"m={m} (predict $\\hbar\\omega$={omega_pred[m-1]*Ha2eV:.2f} eV)")
        ax.axvline(omega_pred[m - 1] * Ha2eV, color=cmap(m - 1), ls=":",
                   lw=0.8, alpha=0.7)
    ax.set_xlabel(r"$\hbar\omega$ (eV)")
    ax.set_ylabel(r"$|\mathrm{FFT}[n_{q_m}(t)]|$")
    ax.set_title(f"{run_name} — plasmon spectrum (Hann + 4x zero-pad, "
                 f"transient t<{t_start_au} a.u. cut)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "n_q_spectrum.png", dpi=140)
    plt.close(fig)
    print(f"  Wrote {out_dir / 'n_q_spectrum.png'}")

    return {
        "times_au": times,
        "n_q": n_q,
        "q_vals": q_vals,
        "omega_pred": omega_pred,
        "omega_eV": omega_eV,
        "spectra": spectra,
        "out_dir": out_dir,
    }


if __name__ == "__main__":
    import sys
    results_dir = Path(sys.argv[1])
    run_name = sys.argv[2] if len(sys.argv) > 2 else results_dir.parent.name
    run(results_dir, run_name)
