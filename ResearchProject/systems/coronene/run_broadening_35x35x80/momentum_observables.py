#!/usr/bin/env python3
"""Higher-order momentum-space observables for run_broadening_35x35x80.

Builds everything that the standard analyse.py could NOT (no momentum_distribution
was saved), directly from the 170 complex WP wavefunction frames in
results/raw/vti/wavefunction_wp_rt/ via FFT:

  1. momentum_distribution.gif   - P(k_z, t) marginal animation
  2. momentum_kz_heatmap.png     - P(k_z) vs time heatmap (collision marked)
  3. momentum_radial.png         - radial P(|k|) before/after
  4. momentum_before_after.png   - P(k_z) at t=0 vs final (deceleration/broadening)
  5. momentum_scatter_map_2d.png - Delta P(k_z, k_perp) after-before, elastic ring overlay
  6. momentum_trajectory.png     - <k_z>(t), sigma_kz(t)  (also the k-space free-evolution check)
  7. energy_loss_spectrum.png    - P(E) before/after, E=|k|^2/2 (EELS-type projectile loss)

psi(r) = real + i*imag (RealField writer applied fft_shift -> physical order, so the
loaded array is physically centred; we FFT and fftshift back to centred k-space).

Known-case test (t=0): P(k_z) must peak at |k_z| = k0 = 3.834 rad/Bohr. ASSERTED.
"""
import glob, os, re, numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
from matplotlib import animation

RUN = "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/run_broadening_35x35x80"
WDIR = f"{RUN}/results/raw/vti/wavefunction_wp_rt"
OUT = f"{RUN}/results/analysis/momentum"
os.makedirs(OUT, exist_ok=True)
DT = 0.02; K0 = 3.83402254402536; T_COLL = 30.0 / K0
HA_EV = 27.211386245988

files = sorted(glob.glob(os.path.join(WDIR, "*.vti")),
               key=lambda p: int(re.search(r"_t(\d+)\.vti$", p).group(1)))
print(f"{len(files)} wavefunction frames")

# --- grid / k-axes from first frame ---
r = vtk.vtkXMLImageDataReader(); r.SetFileName(files[0]); r.Update()
img = r.GetOutput(); nx, ny, nz = img.GetDimensions(); sx, sy, sz = img.GetSpacing()
kx = np.fft.fftshift(2 * np.pi * np.fft.fftfreq(nx, d=sx))
ky = np.fft.fftshift(2 * np.pi * np.fft.fftfreq(ny, d=sy))
kz = np.fft.fftshift(2 * np.pi * np.fft.fftfreq(nz, d=sz))
KX, KY = np.meshgrid(kx, ky, indexing="xy")           # (ny,nx)
KPERP = np.sqrt(KX**2 + KY**2)                          # (ny,nx)

def psik2(path):
    rr = vtk.vtkXMLImageDataReader(); rr.SetFileName(path); rr.Update()
    im = rr.GetOutput().GetPointData()
    re_ = vtk_to_numpy(im.GetArray(0)).reshape(nz, ny, nx)
    im_ = vtk_to_numpy(im.GetArray(1)).reshape(nz, ny, nx)
    psi = re_ + 1j * im_
    pk = np.fft.fftshift(np.fft.fftn(psi))             # centred k-space
    d = np.abs(pk)**2
    return d / d.sum()                                  # normalised |psi(k)|^2, (nz,ny,nx)

# k_perp bins for the 2D map
nperp = 60
pb = np.linspace(0, min(KPERP.max(), K0 * 2.2), nperp + 1)
pcen = 0.5 * (pb[:-1] + pb[1:])

def kz_kperp_map(d):
    """Collapse |psi(k)|^2 (nz,ny,nx) to P(k_z, k_perp) by azimuthal binning."""
    M = np.zeros((nz, nperp))
    idx = np.clip(np.digitize(KPERP.ravel(), pb) - 1, 0, nperp - 1)
    for iz in range(nz):
        M[iz] = np.bincount(idx, weights=d[iz].ravel(), minlength=nperp)
    return M  # (nz, nperp)

ts, mean_kz, sig_kz, Pkz = [], [], [], []
rad_bins = np.linspace(0, K0 * 2.2, 120); rcen = 0.5 * (rad_bins[:-1] + rad_bins[1:])
KMAG = None
first_map = last_map = None; rad_first = rad_last = None
for i, f in enumerate(files):
    step = int(re.search(r"_t(\d+)\.vti$", f).group(1)); t = step * DT
    d = psik2(f)
    pkz = d.sum(axis=(1, 2))                            # P(k_z) marginal, len nz
    ts.append(t); Pkz.append(pkz)
    mkz = (kz * pkz).sum(); vkz = (kz * kz * pkz).sum() - mkz**2
    mean_kz.append(mkz); sig_kz.append(np.sqrt(max(vkz, 0)))
    if i == 0 or i == len(files) - 1:
        if KMAG is None:
            KZ3 = kz[:, None, None]
            KMAG = np.sqrt(KZ3**2 + KPERP[None]**2)     # (nz,ny,nx)
        rad = np.histogram(KMAG.ravel(), bins=rad_bins, weights=d.ravel())[0]
        m = kz_kperp_map(d)
        if i == 0: first_map, rad_first, P0 = m, rad, pkz
        else:      last_map, rad_last, Pend = m, rad, pkz
ts = np.array(ts); Pkz = np.array(Pkz); mean_kz = np.array(mean_kz); sig_kz = np.array(sig_kz)

# ---- KNOWN-CASE TEST: t=0 peak at |k_z| = K0 ----
kz_peak0 = kz[np.argmax(Pkz[0])]
assert abs(abs(kz_peak0) - K0) < 0.25, f"t=0 k_z peak {kz_peak0:.3f} != +/-{K0:.3f}"
print(f"[known-case PASS] t=0 P(k_z) peak = {kz_peak0:.3f} rad/Bohr (|k0|={K0:.3f}); "
      f"<k_z>(0)={mean_kz[0]:.3f}, sigma_kz(0)={sig_kz[0]:.3f}")
SIGN = np.sign(kz_peak0)                                # +1 or -1 travel direction

# =================== PLOTS ===================
# (2) heatmap P(k_z) vs t
fig, ax = plt.subplots(figsize=(6.6, 4.2))
im0 = ax.pcolormesh(ts, kz, Pkz.T, shading="auto", cmap="viridis")
ax.axhline(SIGN * K0, color="w", ls=":", lw=1, label=r"$k_0$")
ax.axvline(T_COLL, color="r", ls="--", lw=1.2, label="collision")
ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$k_z$ (rad/Bohr)")
ax.set_ylim(SIGN * K0 - 6, SIGN * K0 + 6)
ax.set_title("WP longitudinal momentum $P(k_z,t)$"); ax.legend(fontsize=8, loc="upper right")
fig.colorbar(im0, ax=ax, label="probability"); fig.savefig(f"{OUT}/momentum_kz_heatmap.png", dpi=150)
plt.close(fig)

# (4) before/after marginal
fig, ax = plt.subplots(figsize=(6.4, 4.0))
ax.plot(kz, P0, color="#0072B2", lw=2, label=f"t=0 (before)")
ax.plot(kz, Pend, color="#D55E00", lw=2, label=f"t={ts[-1]:.1f} (after)")
ax.axvline(SIGN * K0, color="k", ls=":", lw=1)
ax.set_xlim(SIGN * K0 - 6, SIGN * K0 + 6)
ax.set_xlabel(r"$k_z$ (rad/Bohr)"); ax.set_ylabel("probability")
ax.set_title("WP $P(k_z)$ before vs after collision"); ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.savefig(f"{OUT}/momentum_before_after.png", dpi=150); plt.close(fig)

# (3) radial |k| before/after
fig, ax = plt.subplots(figsize=(6.4, 4.0))
ax.plot(rcen, rad_first, color="#0072B2", lw=2, label="t=0")
ax.plot(rcen, rad_last, color="#D55E00", lw=2, label=f"t={ts[-1]:.1f}")
ax.axvline(K0, color="k", ls=":", lw=1, label=r"$|k_0|$")
ax.set_xlabel(r"$|k|$ (rad/Bohr)"); ax.set_ylabel("probability")
ax.set_title("WP radial momentum $P(|k|)$ before vs after"); ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.savefig(f"{OUT}/momentum_radial.png", dpi=150); plt.close(fig)

# (5) 2D scattering difference map  (k_z, k_perp), after-before, elastic ring
fig, ax = plt.subplots(figsize=(6.8, 4.6))
diff = (last_map - first_map).T                        # (nperp, nz)
vmax = np.abs(diff).max()
im2 = ax.pcolormesh(kz, pcen, diff, shading="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
th = np.linspace(0, np.pi, 200)
ax.plot(SIGN * K0 * np.cos(th), K0 * np.sin(th), "k--", lw=1.2, label=r"elastic ring $|k|=k_0$")
ax.set_xlabel(r"$k_z$ (rad/Bohr)"); ax.set_ylabel(r"$k_\perp$ (rad/Bohr)")
ax.set_xlim(SIGN * K0 - 6, SIGN * K0 + 6); ax.set_ylim(0, K0 * 2.0)
ax.set_title(r"WP scattering map $\Delta P(k_z,k_\perp)$ (after $-$ before)")
ax.legend(fontsize=8, loc="upper right"); fig.colorbar(im2, ax=ax, label=r"$\Delta$ probability")
fig.savefig(f"{OUT}/momentum_scatter_map_2d.png", dpi=150); plt.close(fig)

# (6) momentum trajectory + k-space free-evolution check
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.6, 5.6), sharex=True)
a1.plot(ts, mean_kz, color="#0072B2", lw=1.8); a1.axhline(SIGN * K0, color="k", ls=":", lw=1, label=r"$k_0$")
a1.axvline(T_COLL, color="r", ls="--", lw=1); a1.set_ylabel(r"$\langle k_z\rangle$ (rad/Bohr)")
a1.legend(fontsize=8); a1.grid(alpha=0.3); a1.set_title("WP momentum evolution (pre-collision = free-evolution check)")
a2.plot(ts, sig_kz, color="#009E73", lw=1.8); a2.axvline(T_COLL, color="r", ls="--", lw=1)
a2.set_ylabel(r"$\sigma_{k_z}$ (rad/Bohr)"); a2.set_xlabel("time (a.u.)"); a2.grid(alpha=0.3)
fig.savefig(f"{OUT}/momentum_trajectory.png", dpi=150); plt.close(fig)

# (7) energy-loss (EELS-type) spectrum from |k|^2/2
Ebins = np.linspace(0, (K0 * 1.8)**2 / 2 * HA_EV, 160); Ecen = 0.5 * (Ebins[:-1] + Ebins[1:])
Emag = (KMAG.ravel()**2 / 2) * HA_EV
# recompute first/last full |psik|^2 weights for energy histogram
d0 = psik2(files[0]); dN = psik2(files[-1])
PE0 = np.histogram(Emag, bins=Ebins, weights=d0.ravel())[0]
PEN = np.histogram(Emag, bins=Ebins, weights=dN.ravel())[0]
Emean0 = (Emag * d0.ravel()).sum(); EmeanN = (Emag * dN.ravel()).sum()
fig, ax = plt.subplots(figsize=(6.6, 4.0))
ax.plot(Ecen, PE0, color="#0072B2", lw=2, label=f"before  $\\langle E\\rangle$={Emean0:.1f} eV")
ax.plot(Ecen, PEN, color="#D55E00", lw=2, label=f"after   $\\langle E\\rangle$={EmeanN:.1f} eV")
ax.axvline(K0**2 / 2 * HA_EV, color="k", ls=":", lw=1, label=f"$E_0$={K0**2/2*HA_EV:.0f} eV")
ax.set_xlabel("kinetic energy (eV)"); ax.set_ylabel("probability")
ax.set_title(f"WP energy-loss spectrum  ($\\Delta\\langle E\\rangle$ = {EmeanN-Emean0:+.2f} eV)")
ax.legend(fontsize=8); ax.grid(alpha=0.3); fig.savefig(f"{OUT}/energy_loss_spectrum.png", dpi=150); plt.close(fig)
print(f"energy: <E>_before={Emean0:.2f} eV, <E>_after={EmeanN:.2f} eV, loss={Emean0-EmeanN:+.3f} eV")

# (1) momentum distribution GIF
kzlo, kzhi = SIGN * K0 - 6, SIGN * K0 + 6
sel = (kz >= min(kzlo, kzhi)) & (kz <= max(kzlo, kzhi))
ymax = Pkz[:, sel].max() * 1.1
fig, ax = plt.subplots(figsize=(6.4, 4.0))
(line,) = ax.plot([], [], color="#0072B2", lw=2)
ax.axvline(SIGN * K0, color="k", ls=":", lw=1, label=r"$k_0$")
ax.set_xlim(min(kzlo, kzhi), max(kzlo, kzhi)); ax.set_ylim(0, ymax)
ax.set_xlabel(r"$k_z$ (rad/Bohr)"); ax.set_ylabel("probability"); ax.legend(fontsize=8, loc="upper right")
ax.grid(alpha=0.3); ttl = ax.set_title("")
def upd(i):
    line.set_data(kz, Pkz[i])
    phase = "pre-collision" if ts[i] < T_COLL else "post-collision"
    ttl.set_text(f"WP $P(k_z)$  t={ts[i]:.2f} a.u.  ({phase})")
    line.set_color("#0072B2" if ts[i] < T_COLL else "#D55E00")
    return line, ttl
anim = animation.FuncAnimation(fig, upd, frames=len(ts), blit=False)
anim.save(f"{OUT}/momentum_distribution.gif", writer=animation.PillowWriter(fps=20), dpi=110)
plt.close(fig)
print(f"wrote 7 momentum products to {OUT}")
