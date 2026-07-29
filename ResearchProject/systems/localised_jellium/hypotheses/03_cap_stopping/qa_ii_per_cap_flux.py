#!/usr/bin/env python3
"""Q1(ii): per-CAP cumulative flux = transmission vs reflection.

WP run: exact, from the saved complex wavefunction psi_WP. Probability current
  J_z = Im(conj(psi) d psi/dz); flux(z0,t) = integral of J_z over the xy-plane at z0.
  Forward flux across z=+17.5 (the +z CAP inner edge) -> TRANSMISSION; backward flux
  across z=-17.5 -> REFLECTION. Cumulate in time.
Classical: the projectile is a rigid 1-unit Gaussian moving at v_z(t); its flux across
  a plane is v_z * (Gaussian marginal at the plane). It transits forward (T~1) and
  never moves -z (R=0) — the ballistic reference.

Validation: integral|psi|^2 must match N_wp(density_wp); and T_end+R_end must equal the
WP that has LEFT the inner region [-17.5,17.5] (= 1 - inner WP).
Outputs: qa_ii_per_cap_flux.png + .csv in this directory.
"""
import glob
import re
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.stats import norm as gauss
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from inqview.visualisation import style

style.apply_theme()
DT = 0.02
SIG_CHARGE = 0.350
Z_PLUS, Z_MINUS = 17.5, -17.5
ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
WP = f"{ROOT}/scripts/fullsuite_wp/results/p5_wp/raw"
CL = f"{ROOT}/scripts/fullsuite_classical/results/p5_classical/raw"
OUT = f"{ROOT}/hypotheses/03_cap_stopping"


def ftime(p):
    return int(re.search(r"_t(\d+)\.vti", p).group(1)) * DT


def load_psi(path):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(path); r.Update()
    img = r.GetOutput()
    nx, ny, nz = img.GetDimensions()
    pdt = img.GetPointData()
    re_ = vtk_to_numpy(pdt.GetArray("wavefunction_real")).reshape((nz, ny, nx)).transpose(2, 1, 0)
    im_ = vtk_to_numpy(pdt.GetArray("wavefunction_imag")).reshape((nz, ny, nx)).transpose(2, 1, 0)
    oz = img.GetOrigin()[2]; sx, sy, sz = img.GetSpacing()
    z = oz + (np.arange(nz) + 0.5) * sz
    return (re_ + 1j * im_), z, (sx, sy, sz)


# --- WP probability current flux at the two CAP inner edges ---
wff = sorted(glob.glob(f"{WP}/vti/wavefunction_wp/wavefunction_t*.vti"), key=ftime)
t_f = np.array([ftime(p) for p in wff])
flux_p = np.zeros(len(wff))   # at +17.5
flux_m = np.zeros(len(wff))   # at -17.5
norm_psi = np.zeros(len(wff))
def dpsi_dz_spectral(psi, sz):
    """Exact z-derivative on the periodic box (FFT) — avoids the finite-difference
    underestimate of high-k current (k0=2.71 is only ~4.6 pts/wavelength at dx=0.5)."""
    nz = psi.shape[2]
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=sz)
    return np.fft.ifft(1j * kz[None, None, :] * np.fft.fft(psi, axis=2), axis=2)

for i, p in enumerate(wff):
    psi, z, (sx, sy, sz) = load_psi(p)
    Jz = np.imag(np.conj(psi) * dpsi_dz_spectral(psi, sz))
    ip = int(np.argmin(np.abs(z - Z_PLUS)))
    im = int(np.argmin(np.abs(z - Z_MINUS)))
    flux_p[i] = Jz[:, :, ip].sum() * sx * sy
    flux_m[i] = Jz[:, :, im].sum() * sx * sy
    norm_psi[i] = (np.abs(psi) ** 2).sum() * sx * sy * sz

# cumulative: forward across +17.5 = transmission; backward across -17.5 = reflection
T = cumulative_trapezoid(flux_p, t_f, initial=0.0)
R = cumulative_trapezoid(-flux_m, t_f, initial=0.0)

# --- classical ballistic flux (rigid Gaussian * v_z) ---
trk = pd.read_csv(f"{CL}/observables/electron_track.csv").drop_duplicates("step")
tt = trk.time_au.values
zion = ((trk.z.values + 25.0) % 50.0) - 25.0
vz = trk.vz.values
fcl_p = vz * gauss.pdf(Z_PLUS, loc=zion, scale=SIG_CHARGE)
fcl_m = vz * gauss.pdf(Z_MINUS, loc=zion, scale=SIG_CHARGE)
Tcl = cumulative_trapezoid(fcl_p, tt, initial=0.0)
Rcl = cumulative_trapezoid(-fcl_m, tt, initial=0.0)

# --- validation cross-checks ---
# WP inner-region norm at the last frame (slab+left+right bands), from density_wp:
inner_wp = pd.read_csv(f"{OUT}/qa_i_region_densities.csv")
inner_last = float(inner_wp.iloc[-1][["wp_leftfree", "wp_slab", "wp_rightfree"]].sum())
chk_TR = 1.0 - inner_last

# --- plot ---
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(t_f, T, "C0-", lw=1.8, label=f"WP transmission (+z CAP): {T[-1]:.3f}")
ax.plot(t_f, R, "C3-", lw=1.8, label=f"WP reflection (-z CAP): {R[-1]:.3f}")
ax.plot(tt, Tcl, "C0--", lw=1.3, label=f"classical transmission: {Tcl[-1]:.3f}")
ax.plot(tt, Rcl, "C3--", lw=1.3, label=f"classical reflection: {Rcl[-1]:.3f}")
ax.set_xlabel("time (a.u.)")
ax.set_ylabel("cumulative norm across CAP edge  (electrons)")
ax.set_title("Per-CAP flux: transmission vs reflection — WP vs classical (slab, 100 eV)",
             fontsize=9)
ax.legend(fontsize=7, frameon=False, loc="upper left")
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(f"{OUT}/qa_ii_per_cap_flux.png", dpi=200)
plt.close(fig)

pd.DataFrame({"time_au": t_f, "flux_plus": flux_p, "flux_minus": flux_m,
              "T_cum": T, "R_cum": R, "norm_psi": norm_psi}).to_csv(
    f"{OUT}/qa_ii_per_cap_flux.csv", index=False)

print(f"|psi|^2 norm check: t0={norm_psi[0]:.4f} tEnd={norm_psi[-1]:.4f} "
      f"(density_wp gave 1.000 / 0.378)")
print(f"WP : transmission(+z)={T[-1]:.3f}  reflection(-z)={R[-1]:.3f}  "
      f"T+R={T[-1]+R[-1]:.3f}")
print(f"     cross-check 1-inner_WP = {chk_TR:.3f}  (should ~= T+R)")
print(f"     reflection fraction R/(T+R) = {R[-1]/(T[-1]+R[-1]):.1%}")
print(f"classical: transmission={Tcl[-1]:.3f}  reflection={Rcl[-1]:.3f} (ballistic)")
print("wrote qa_ii_per_cap_flux.png + .csv")
