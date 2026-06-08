"""Run-2b: Python Schrödinger split-step FFT toy.

Independent free-particle propagation of a Gaussian WP at the same
parameters as Run-2 (the INQ non-interacting run in this directory).
The three traces of sigma_r²(t) we then have:

  1. INQ non-interacting       — results/raw/observables/wp_real_space_stats.csv
  2. Python split-step (this)  — results/raw/observables/python_toy_sigma.csv
  3. Analytical Gaussian       — written here as a third column

All three must agree to ≤ 1% over the IFW window for the WP injector
+ inq propagator to be considered validated against free-particle
physics (plan §"Per-family motivations" Run-2+2b+3).

Physics: free-particle Schrödinger equation in 1D along z (the 3D
problem separates for a free Gaussian, so only z propagates the
launch momentum; transverse axes spread by the same analytic rule).
Split-step method:

    psi(t+dt) = F^-1[ exp(-i k² dt / 2) * F[psi(t)] ]   (kinetic-only)

where F is the FFT, exp(-i k² dt / 2) is the kinetic propagator in
k-space, and m_e = 1 in atomic units.

Convention reminder: the inqkit::WavePacket injector writes
psi(z) ~ exp(-z² / (2 sigma²)), so the *density* sigma is
sigma / sqrt(2). Both the INQ run and this toy use the same psi
convention, so the wavefunction-sigma == sigma comparison is direct.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np


# ----- Match Run-2 physics parameters ------------------------------------
# The split-step FFT uses periodic boundaries; if the Gaussian density tail
# reaches the box edge it wraps and pollutes the second moment. So we use a
# much wider box than Run-2 (which uses .finite() walls). The physics is
# box-size invariant for the free particle, so this is just headroom.
L_BOHR        = 200.0
DX_BOHR       = 0.40
SIGMA_W       = 5.0
WP_K0         = math.sqrt(2.0 * 100.0 / 27.21138625)   # 2.7110633...
LAUNCH_Z      = -5.0   # boundary_rule launch_z(σ=5, L=50)
DT_AU         = 0.020
N_STEPS       = 462
WRITE_EVERY   = 2

OUT_CSV = Path(__file__).parent / "results" / "raw" / "observables" \
                                  / "python_toy_sigma.csv"


def _free_particle_sigma_density_au(t_au: np.ndarray,
                                     sigma_w: float,
                                     m_au: float = 1.0) -> np.ndarray:
    """Analytic density σ_r(t) for a free Gaussian.

    For psi(r,0) ~ exp(-r² / (2 σ²)) (injector convention),
        |psi(r,t)|² is Gaussian with std-dev
        σ_density(t) = (σ/√2) * sqrt(1 + (t / (m σ²))²)
    """
    s0 = sigma_w / math.sqrt(2.0)              # density σ at t=0
    return s0 * np.sqrt(1.0 + (t_au / (m_au * sigma_w * sigma_w)) ** 2)


def split_step_1d(z: np.ndarray, dt_au: float, n_steps: int,
                  sigma_w: float, k0: float, launch_z: float,
                  write_every: int):
    """Free-particle split-step FFT propagation of a 1D Gaussian wp.

    Returns (steps, times_au, sigmas_density) where sigmas_density is
    the time series of std-dev of |psi(z, t)|².
    """
    dz = z[1] - z[0]
    nz = z.size
    # k-grid for FFT, consistent with numpy's fftfreq convention.
    k = 2.0 * math.pi * np.fft.fftfreq(nz, d=dz)
    kinetic_phase = np.exp(-0.5j * (k * k) * dt_au)        # m = 1

    # Initial WP: psi(z, 0) = norm * exp(-(z - z0)² / (2σ²)) * exp(i k0 z)
    psi = np.exp(-((z - launch_z) ** 2) / (2.0 * sigma_w * sigma_w)) \
          * np.exp(1.0j * k0 * z)
    psi /= math.sqrt(float((np.abs(psi) ** 2 * dz).sum()))  # normalise to 1

    steps = []
    times = []
    sigmas = []
    for step in range(n_steps + 1):
        if step % write_every == 0:
            rho = np.abs(psi) ** 2
            zmean = float((z * rho * dz).sum())
            z2mean = float((z * z * rho * dz).sum())
            steps.append(step)
            times.append(step * dt_au)
            sigmas.append(math.sqrt(max(0.0, z2mean - zmean * zmean)))

        if step == n_steps:
            break
        # split-step: full kinetic step (no potential, free particle)
        psi = np.fft.ifft(kinetic_phase * np.fft.fft(psi))
    return np.array(steps), np.array(times), np.array(sigmas)


def main() -> None:
    nz = int(round(L_BOHR / DX_BOHR))
    z  = np.linspace(-L_BOHR / 2.0, L_BOHR / 2.0, nz, endpoint=False) \
         + DX_BOHR * 0.5
    print(f"Run-2b toy: nz={nz}, dz={DX_BOHR}, k0={WP_K0:.4f}, "
          f"sigma_w={SIGMA_W}")

    steps, t_au, sigma_split = split_step_1d(
        z, DT_AU, N_STEPS, SIGMA_W, WP_K0, LAUNCH_Z, WRITE_EVERY)
    sigma_analytic = _free_particle_sigma_density_au(t_au, SIGMA_W)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w") as f:
        f.write("step,time_au,sigma_density_split_step_au,"
                "sigma_density_analytic_au\n")
        for s, t, a, b in zip(steps, t_au, sigma_split, sigma_analytic):
            f.write(f"{int(s)},{t:.6f},{a:.10f},{b:.10f}\n")

    diff = np.abs(sigma_split - sigma_analytic)
    rel  = diff / np.maximum(sigma_analytic, 1e-12)
    print(f"  wrote {len(steps)} rows to {OUT_CSV}")
    print(f"  σ_density (t=0)        = {sigma_split[0]:.5f}  "
          f"(analytic: {sigma_analytic[0]:.5f})")
    print(f"  σ_density (t={t_au[-1]:.2f}) = {sigma_split[-1]:.5f}  "
          f"(analytic: {sigma_analytic[-1]:.5f})")
    print(f"  max |split − analytic| = {diff.max():.2e}  "
          f"(rel: {rel.max():.2e})")
    assert rel.max() < 1e-3, \
        f"split-step disagrees with analytic by {rel.max():.2e}"
    print("PASS")


if __name__ == "__main__":
    main()
