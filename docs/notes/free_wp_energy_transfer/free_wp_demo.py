#!/usr/bin/env python3
"""Free-wavepacket energy-store transfer — quick exact experiment.

Question (user, 2026-07-16): during FREE propagation of a wavepacket the packet
spreads and one expects the "localisation energy" to fall. Which component of the
kinetic energy changes so that TOTAL energy is conserved?

Answer being demonstrated: total KE (= <p^2>/2m) is a frozen constant of a free
particle (momentum distribution is conserved). What changes is the REAL-SPACE
partition of that KE (Madelung / hydrodynamic decomposition):

    KE_total = T_W(localisation) + T_flow(current)
      T_W    = 1/2 . integral |grad sqrt(n)|^2   (von Weizsaecker; density SHAPE)
      T_flow = 1/2 . integral n |v|^2 , v = j/n  (probability-current KE)

T_W falls as the packet delocalises; T_flow rises as it develops outward flow;
their sum is pinned to the conserved momentum-space KE.

Method: exact split-operator (Strang) propagation of a free 1D Gaussian, hbar=m=1.
1D is sufficient — the free Gaussian is separable, so a 3D isotropic packet is a
product of three of these and the energy stores simply add per axis. We also print
the analytic 3D Gaussian self-Hartree (the pairwise e_pp a KS run would carry) to
show how it would decay as ~1/sigma(t) — but note a genuinely FREE particle has NO
Coulomb interaction; that term only appears once the WP is embedded in a mean field
(and is exactly the self-interaction error we track in the twin runs).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).resolve().parent

# ------------------------------------------------------------------ parameters
hbar = 1.0
m = 1.0
sigma0 = 0.5          # initial DENSITY std (Bohr) — matches the sigma_WP=0.5 anchor
k0 = 0.0              # mean momentum: 0 => packet at rest, purely spreading
L = 200.0             # box (Bohr) — large so the packet never wraps
N = 8192              # grid points
T_END = 6.0           # a.u. of time
N_OUT = 13            # number of report rows

x = (np.arange(N) - N // 2) * (L / N)
dx = L / N
k = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)   # angular wavenumber grid

# ------------------------------------------------------------- initial packet
# amplitude std sigma_a = sqrt(2)*sigma0 so that |psi|^2 has std = sigma0
sigma_a = np.sqrt(2.0) * sigma0
psi = np.exp(-x**2 / (4.0 * sigma_a**2 / 2.0)) * np.exp(1j * k0 * x)
# ^ |psi|^2 = exp(-x^2/(2 sigma0^2)) => density std sigma0. (4*sigma_a^2/2 = 2*sigma0^2*2? verify below)
# Build density-std directly and renormalise, to avoid convention slips:
psi = np.exp(-x**2 / (4.0 * sigma0**2)) * np.exp(1j * k0 * x)   # |psi|^2 ~ exp(-x^2/(2 sigma0^2))
psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx)

# --------------------------------------------------------------- propagators
dt = 1e-3
Kfac = np.exp(-1j * (hbar * k**2 / (2.0 * m)) * dt)         # full kinetic step
# (no potential — free particle — so Strang reduces to pure kinetic phase in k-space)

def energies(psi):
    """Return dict of energy components for the current state."""
    n = np.abs(psi)**2
    norm = np.sum(n) * dx
    # --- momentum-space (total KE and its frozen split) ---
    psik = np.fft.fft(psi) * dx
    Pk = np.abs(psik)**2
    Pk /= np.sum(Pk) * (k[1]-k[0] if False else 1.0)  # not used for norm; use direct expectation below
    # expectation values via k-grid
    dk = k[1] - k[0]
    Wk = np.abs(np.fft.fft(psi))**2
    Wk /= np.sum(Wk)                       # normalised momentum pdf on the k grid
    p_mean = np.sum(k * Wk)
    p2_mean = np.sum(k**2 * Wk)
    KE_tot = 0.5 * p2_mean / m
    KE_trans = 0.5 * p_mean**2 / m
    KE_spread = 0.5 * (p2_mean - p_mean**2) / m
    # --- real-space Madelung split ---
    sqrtn = np.sqrt(np.maximum(n, 0.0))
    dsqrtn = np.gradient(sqrtn, dx)
    T_W = 0.5 * hbar**2 / m * np.sum(dsqrtn**2) * dx           # localisation
    # current j = hbar/m Im(psi* dpsi/dx)
    dpsi = np.gradient(psi, dx)
    j = (hbar / m) * np.imag(np.conj(psi) * dpsi)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(n > 1e-14, j / n, 0.0)
    T_flow = 0.5 * m * np.sum(n * v**2) * dx                   # current KE
    # density std
    xm = np.sum(x * n) * dx / norm
    s = np.sqrt(np.sum((x - xm)**2 * n) * dx / norm)
    return dict(norm=norm, KE_tot=KE_tot, KE_trans=KE_trans, KE_spread=KE_spread,
                T_W=T_W, T_flow=T_flow, T_sum=T_W + T_flow, s=s)

# analytic 3D isotropic Gaussian self-Hartree (charge -1, density std s): 1/(sigma_rho*sqrt(pi))
# with sigma_rho the 1D std of the 3D isotropic density. U_H = 1/(2 sigma_rho sqrt(pi)) * ...
# memory: 1/(2 sigma_rho sqrt(pi)) = 0.80 Ha for sigma_rho=... use closed form U_H = 1/(sqrt(pi) s) * (1/2)
def self_hartree_3d(s):
    return 1.0 / (2.0 * s * np.sqrt(np.pi))   # Ha, for a unit-charge 3D isotropic Gaussian, density std s

# ----------------------------------------------------------------- integrate
n_steps = int(round(T_END / dt))
out_every = n_steps // (N_OUT - 1)
rows = []
E0 = None
for step in range(n_steps + 1):
    if step % out_every == 0 or step == n_steps:
        e = energies(psi)
        t = step * dt
        e["t"] = t
        e["U_H"] = self_hartree_3d(e["s"])
        if E0 is None:
            E0 = e["KE_tot"]
        e["dE_rel"] = (e["KE_tot"] - E0) / E0
        rows.append(e)
    if step < n_steps:
        psi = np.fft.ifft(np.fft.fft(psi) * Kfac)

# ------------------------------------------------------------------- report
hdr = f"{'t':>6} {'std':>7} {'KE_tot':>9} {'KE_spread':>10} {'KE_trans':>9} " \
      f"{'T_W(loc)':>9} {'T_flow':>9} {'T_W+T_flow':>10} {'dE/E':>10} {'U_H(3D)':>8}"
print(hdr)
print("-" * len(hdr))
for e in rows:
    print(f"{e['t']:6.2f} {e['s']:7.3f} {e['KE_tot']:9.5f} {e['KE_spread']:10.5f} "
          f"{e['KE_trans']:9.5f} {e['T_W']:9.5f} {e['T_flow']:9.5f} {e['T_sum']:10.5f} "
          f"{e['dE_rel']:10.2e} {e['U_H']:8.4f}")

# analytic cross-checks
print("\n--- analytic cross-checks (hbar=m=1) ---")
sp = 1.0 / (2.0 * sigma0)                      # momentum std of amplitude
print(f"analytic total KE (k0=0) = 1/2 * sigma_p^2 = {0.5*sp**2:.5f} Ha  "
      f"(sigma_p = 1/(2 sigma0) = {sp:.4f})")
print(f"numeric  total KE (t=0)  = {rows[0]['KE_tot']:.5f} Ha")
print(f"T_W(t=0) should equal total KE (no flow at t=0): "
      f"T_W={rows[0]['T_W']:.5f}, T_flow={rows[0]['T_flow']:.2e}")

# ------------------------------------------------------------------- figures
ts = np.array([e["t"] for e in rows])
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(ts, [e["KE_tot"] for e in rows], "k-", lw=2.2, label="KE total (frozen)")
ax[0].plot(ts, [e["T_W"] for e in rows], "C3-o", ms=4, label="T_W  localisation")
ax[0].plot(ts, [e["T_flow"] for e in rows], "C0-s", ms=4, label="T_flow  current")
ax[0].plot(ts, [e["T_sum"] for e in rows], "C2--", lw=1.2, label="T_W + T_flow")
ax[0].set_xlabel("time  (a.u.)")
ax[0].set_ylabel("energy  (Ha)")
ax[0].set_title("Kinetic-energy stores of a free spreading WP")
ax[0].legend(fontsize=8, frameon=False)

ax2 = ax[1]
ax2.plot(ts, [e["s"] for e in rows], "C4-o", ms=4, label="density std  σ(t)")
ax2.set_xlabel("time  (a.u.)")
ax2.set_ylabel("density std  (Bohr)", color="C4")
ax2b = ax2.twinx()
ax2b.plot(ts, [e["U_H"] for e in rows], "C1-s", ms=4, label="self-Hartree e_pp  (~1/σ)")
ax2b.set_ylabel("self-Hartree  (Ha)", color="C1")
ax2.set_title("Spreading and the self-Hartree (pairwise e_pp) it would carry")
fig.tight_layout()
fig.savefig(OUT / "free_wp_energy_stores.png", dpi=140)
print(f"\nwrote {OUT/'free_wp_energy_stores.png'}")
