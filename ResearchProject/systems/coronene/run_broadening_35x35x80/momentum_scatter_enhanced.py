#!/usr/bin/env python3
"""Enhanced momentum scatter map + 3D directional scattering analysis.

Recomputes the WP momentum-difference DP(k) = |psi_after(k)|^2 - |psi_before(k)|^2
from the first/last complex wavefunction frames (one FFT each), saves the arrays,
and produces:

  momentum_scatter_map_enhanced.{pdf,png}  - annotated 2D DP(k_z,k_perp) map
  momentum_3d_directions.{pdf,png}         - 3D scattering directions on coronene
  momentum_scatter_arrays.npz              - k grids + 2D map + 3D DP

Convention: WP drifts in -z (<k_z> = -k0), so transmission = k_z<0 (forward),
backscatter = k_z>0.  k0 = sqrt(2*200/27.211) = 3.834 Bohr^-1.

Known-case test: both distributions normalised to sum=1, so DP must integrate to 0.
"""
import glob, os, re, numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from mpl_toolkits.mplot3d import Axes3D  # noqa

RUN = "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/run_broadening_35x35x80"
WDIR = f"{RUN}/results/raw/vti/wavefunction_wp_rt"
OUT = f"{RUN}/results/analysis/momentum"
XYZ = "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/shared/geometry/coronene.xyz"
os.makedirs(OUT, exist_ok=True)
K0 = np.sqrt(2 * 200 / 27.211)        # 3.834 Bohr^-1
SHELL = 0.15                           # elastic-band half-width (Bohr^-1)

# ---------- load first & last wavefunction frames, FFT ----------
files = sorted(glob.glob(os.path.join(WDIR, "*.vti")),
               key=lambda p: int(re.search(r"_t(\d+)\.vti$", p).group(1)))

def psik2(path):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(path); r.Update()
    img = r.GetOutput(); nx, ny, nz = img.GetDimensions(); sp = img.GetSpacing()
    pd = img.GetPointData()
    psi = (vtk_to_numpy(pd.GetArray(0)) + 1j * vtk_to_numpy(pd.GetArray(1))).reshape(nz, ny, nx)
    d = np.abs(np.fft.fftshift(np.fft.fftn(psi)))**2
    return d / d.sum(), (nx, ny, nz), sp

Pb, (nx, ny, nz), sp = psik2(files[0])
Pa, _, _ = psik2(files[-1])
dP3 = Pa - Pb                                            # (nz,ny,nx), sums to 0

# known-case test
assert abs(dP3.sum()) < 1e-9, f"DP does not integrate to 0: {dP3.sum():.2e}"
print(f"[known-case PASS] sum(DP) = {dP3.sum():.2e} (~0); frames {len(files)}")

kx = np.fft.fftshift(2 * np.pi * np.fft.fftfreq(nx, d=sp[0]))
ky = np.fft.fftshift(2 * np.pi * np.fft.fftfreq(ny, d=sp[1]))
kz = np.fft.fftshift(2 * np.pi * np.fft.fftfreq(nz, d=sp[2]))
KZ, KY, KX = np.meshgrid(kz, ky, kx, indexing="ij")       # (nz,ny,nx)
KMAG = np.sqrt(KX**2 + KY**2 + KZ**2)
KPERP = np.sqrt(KX**2 + KY**2)

# ---------- 3D probabilities (exact, per-voxel probability) ----------
def msum(m): return float(dP3[m].sum())
inel = KMAG < (K0 - SHELL)
elas = np.abs(KMAG - K0) < SHELL
trans = KZ < 0
back = KZ > 0
I = {
    "total":        dP3.sum(),
    "inelastic":    msum(inel),
    "elastic_band": msum(elas),
    "transmission": msum(trans),
    "backscatter":  msum(back),
}
# absolute forward/back probability of the AFTER distribution
P_trans_after = float(Pa[KZ < 0].sum()); P_back_after = float(Pa[KZ > 0].sum())
print("=== 3D DP integrals (change in probability) ===")
for k, v in I.items(): print(f"  {k:14s} = {v:+.4e}")
print(f"  P(transmit) after = {P_trans_after:.4f}   P(backscatter) after = {P_back_after:.4f}")

# ---------- 2D azimuthal-mean map M(k_z, k_perp) ----------
nperp = 70
pb = np.linspace(0, K0 * 2.0, nperp + 1); pcen = 0.5 * (pb[:-1] + pb[1:])
idx = np.clip(np.digitize(KPERP.reshape(nz, -1), pb) - 1, 0, nperp - 1)  # (nz, ny*nx)
flat = dP3.reshape(nz, -1)
M = np.zeros((nz, nperp))
for iz in range(nz):
    s = np.bincount(idx[iz], weights=flat[iz], minlength=nperp)
    c = np.bincount(idx[iz], minlength=nperp).astype(float); c[c == 0] = 1
    M[iz] = s / c                                          # azimuthal MEAN (density)

np.savez_compressed(f"{OUT}/momentum_scatter_arrays.npz",
                    kz=kz, kperp=pcen, M=M, kx=kx, ky=ky, dP3=dP3.astype(np.float32),
                    k0=K0, integrals=np.array([I[k] for k in I]))

# =================== PLOT 1: enhanced 2D map ===================
vmax = 0.8 * np.percentile(np.abs(M), 99)
fig, ax = plt.subplots(figsize=(10, 7), dpi=150)
im = ax.pcolormesh(kz, pcen, M.T, shading="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
cb = fig.colorbar(im, ax=ax, label=r"$\Delta P$")

# (4) region shading
ax.add_patch(Circle((0, 0), K0, facecolor="blue", alpha=0.05, edgecolor="none",
                    clip_on=True, zorder=0))
ax.add_patch(Rectangle((0, 0), kz.max() + 1, pcen.max() + 1, facecolor="orange",
                       alpha=0.05, edgecolor="none", clip_on=True, zorder=0))
ax.text(-K0 * 0.45, K0 * 0.30, "inelastic\n$(|k|<k_0)$", fontsize=8, color="navy",
        ha="center", va="center")
ax.text(K0 * 0.55, K0 * 1.55, "backscatter", fontsize=8, color="#b35900", ha="center")

# (1) transmission / backscatter divide
ax.axvline(0.0, color="#444", lw=1.5, ls="--")
ax.text(0.02, 0.97, "Transmission ($k_z<0$)", transform=ax.transAxes, ha="left",
        va="top", fontsize=9, color="#222")
ax.text(0.98, 0.97, "Backscatter ($k_z>0$)", transform=ax.transAxes, ha="right",
        va="top", fontsize=9, color="#222")
# (6) beam axis
ax.axhline(0.0, color="#444", lw=0.8, ls="--", alpha=0.7)

# (2) elastic ring (upper semicircle)
kzr = np.linspace(-K0, K0, 500); kpr = np.sqrt(np.clip(K0**2 - kzr**2, 0, None))
ax.plot(kzr, kpr, "k--", lw=1.5, label=r"$|k|=k_0=3.83\ \mathrm{Bohr}^{-1}$")

# (3) scattering-angle annotations on the ring
for th in [0, 30, 60, 90, 120, 150, 180]:
    a = np.deg2rad(th); zc = -K0 * np.cos(a); pc = K0 * np.sin(a)
    off = 0.45
    ax.annotate(rf"$\theta={th}^\circ$", (zc, pc),
                xytext=(zc * (1 + off * 0.15) + 0.1, pc + off), fontsize=8,
                ha="center", color="k")
    ax.plot([zc], [pc], "k.", ms=3)

# (5) integral textbox (3D-exact values)
box = (f"$\\Delta P$ integrals (3D)\n"
       f"total = {I['total']:+.1e}\n"
       f"inelastic = {I['inelastic']:+.3f}\n"
       f"elastic band = {I['elastic_band']:+.3f}\n"
       f"transmit = {I['transmission']:+.3f}\n"
       f"backscatter = {I['backscatter']:+.3f}")
ax.text(0.015, 0.03, box, transform=ax.transAxes, fontsize=8, va="bottom", ha="left",
        bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9))

ax.set_xlim(kz.min() * 0.55, kz.max() * 0.55)
ax.set_ylim(0, K0 * 1.9)
ax.set_xlabel(r"$k_z$  (Bohr$^{-1}$)"); ax.set_ylabel(r"$k_\perp$  (Bohr$^{-1}$)")
ax.set_title(r"$\Delta P(k_z, k_\perp)$ | E = 200 eV | e$^-$ + coronene")
ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
fig.tight_layout()
fig.savefig(f"{OUT}/momentum_scatter_map_enhanced.pdf")
fig.savefig(f"{OUT}/momentum_scatter_map_enhanced.png", dpi=150)
plt.close(fig)
print("wrote momentum_scatter_map_enhanced.{pdf,png}")

# =================== PLOT 2: 3D directional scattering ===================
# dominant positive-DP directions on the elastic shell -> real-space directions (k_hat)
cand = elas & (dP3 > 0)
ii = np.where(cand.ravel())[0]
vals = dP3.ravel()[ii]
order = ii[np.argsort(vals)[::-1]]
kxf, kyf, kzf, kmf = KX.ravel(), KY.ravel(), KZ.ravel(), KMAG.ravel()
dirs, weights = [], []
for j in order:
    u = np.array([kxf[j], kyf[j], kzf[j]]) / kmf[j]
    if all(np.dot(u, d) < np.cos(np.deg2rad(18)) for d in dirs):  # >18 deg apart
        dirs.append(u); weights.append(dP3.ravel()[j])
    if len(dirs) >= 12:
        break
dirs = np.array(dirs); weights = np.array(weights)
print(f"{len(dirs)} dominant scattering channels found")

# coronene atoms
atoms = []
if os.path.exists(XYZ):
    with open(XYZ) as f:
        lines = f.read().splitlines()
    for ln in lines[2:]:
        p = ln.split()
        if len(p) >= 4:
            atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
A = np.array([[a[1], a[2], a[3]] for a in atoms]) if atoms else np.zeros((0, 3))
sym = [a[0] for a in atoms]

fig = plt.figure(figsize=(10, 7), dpi=150)
ax = fig.add_subplot(111, projection="3d")
# molecule
if len(A):
    cc = ["#222" if s == "C" else "#bbb" for s in sym]
    ax.scatter(A[:, 0], A[:, 1], A[:, 2], c=cc, s=40, depthshade=True, edgecolor="k", lw=0.3)
# incident beam (from +z toward -z)
ax.quiver(0, 0, 12, 0, 0, -6, color="green", lw=2.5, arrow_length_ratio=0.25)
ax.text(0, 0, 13, "incident", color="green", fontsize=9, ha="center")
# scattered arrows: direction = k_hat, colour transmit(blue)/back(red), intensity by DP
L = 8.0
wn = weights / weights.max()
for u, w in zip(dirs, wn):
    col = plt.cm.Blues(0.4 + 0.6 * w) if u[2] < 0 else plt.cm.Reds(0.4 + 0.6 * w)
    ax.quiver(0, 0, 0, u[0] * L, u[1] * L, u[2] * L, color=col, lw=1.8,
              arrow_length_ratio=0.18)
import matplotlib.lines as mlines
ax.legend(handles=[mlines.Line2D([], [], color=plt.cm.Blues(0.8), lw=2, label="transmitted ($k_z<0$)"),
                   mlines.Line2D([], [], color=plt.cm.Reds(0.8), lw=2, label="backscattered ($k_z>0$)"),
                   mlines.Line2D([], [], color="green", lw=2, label="incident")],
          fontsize=8, loc="upper left")
ax.set_xlabel("x (Bohr)"); ax.set_ylabel("y (Bohr)"); ax.set_zlabel("z (Bohr)")
ax.set_title(f"Scattering directions from $\\Delta P$ on $|k|=k_0$  "
             f"(transmit {P_trans_after:.2f} / back {P_back_after:.3f})")
lim = 10; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
ax.view_init(elev=12, azim=-60)
fig.tight_layout()
fig.savefig(f"{OUT}/momentum_3d_directions.pdf")
fig.savefig(f"{OUT}/momentum_3d_directions.png", dpi=150)
plt.close(fig)
print("wrote momentum_3d_directions.{pdf,png}")

# 2D cylindrical cross-check of the integrals (per user's formula, on the mean map)
dkz = np.diff(kz).mean(); dkp = np.diff(pcen).mean()
integ2d = (M * (2 * np.pi * pcen[None, :]) * dkz * dkp)
print(f"[2D cylindrical cross-check] total = {integ2d.sum():+.2e} (density form; "
      f"3D per-voxel is the authoritative one above)")
