#!/usr/bin/env python3
"""Q&A ix: signed WP probability-current flux across ALL internal planes vs time.

Generalises qa_ii (which used only the two CAP inner edges) to every internal
boundary of the region partition, to expose the DIRECTION of charge transport at
each interface and test the "fast components reflect/transmit first" intuition.

Planes (z, Bohr):  -17.5 (-z CAP inner edge) | -12.5 (near face) |
                   +12.5 (far face)          | +17.5 (+z CAP inner edge).
At each plane: flux(z0,t) = int_xy J_z(z0) dx dy, with J_z = Im(conj(psi) d psi/dz)
(spectral z-derivative). Sign convention: + = forward (+z), - = backward (-z).
  Panel A: cumulative int_0^t flux dt across each plane (the transport budget).
  Panel B: instantaneous flux at the near face (-12.5) and -z CAP edge (-17.5),
           zoomed early — a backward (negative) excursion here = genuine reflection,
           and its TIMING tests whether fast components turn back first.

Wrap caveat: beyond t ~ 14.9 the transmitted WP periodic-wraps (+25 -> -25) and
re-enters moving +z, contaminating the -z-side planes. The wrap-free window is
t < ~14.9 (annotated). Frame cadence is dt=0.2 a.u. (91 frames) — a brief early
reflection pulse may be undersampled; flagged, not hidden.

Outputs: qa_ix_cumulative_current_regions.png + .csv in this directory.
"""
import glob
import re
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from inqview.visualisation import style

style.apply_theme()
DT = 0.02
ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
WP = f"{ROOT}/scripts/fullsuite_wp/results/p5_wp/raw"
OUT = f"{ROOT}/hypotheses/03_cap_stopping"

PLANES = [(-17.5, "-z CAP edge (-17.5)", "C3"),
          (-12.5, "near face (-12.5)", "C1"),
          (12.5,  "far face (+12.5)",  "C0"),
          (17.5,  "+z CAP edge (+17.5)", "C4")]
# kinematic anchors (from qa_viii / wp_real_space_stats), annotated as vlines
EVENTS = [(0.4, "lead enters near face"), (1.2, "centroid enters near face"),
          (4.0, "lead reaches far face"), (9.0, "centroid stalls"),
          (14.9, "periodic wrap begins")]
T_WRAP = 14.9


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


def dpsi_dz_spectral(psi, sz):
    nz = psi.shape[2]
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=sz)
    return np.fft.ifft(1j * kz[None, None, :] * np.fft.fft(psi, axis=2), axis=2)


wff = sorted(glob.glob(f"{WP}/vti/wavefunction_wp/wavefunction_t*.vti"), key=ftime)
t = np.array([ftime(p) for p in wff])
flux = {zp: np.zeros(len(wff)) for zp, _, _ in PLANES}
for i, p in enumerate(wff):
    psi, z, (sx, sy, sz) = load_psi(p)
    Jz = np.imag(np.conj(psi) * dpsi_dz_spectral(psi, sz))
    for zp, _, _ in PLANES:
        iz = int(np.argmin(np.abs(z - zp)))
        flux[zp][i] = Jz[:, :, iz].sum() * sx * sy

cum = {zp: cumulative_trapezoid(flux[zp], t, initial=0.0) for zp, _, _ in PLANES}

# --- report: does the -z CAP edge ever run net-backward in the wrap-free window? ---
m = t <= T_WRAP
fminus = flux[-17.5]
print("Wrap-free window (t <= 14.9):")
print(f"  -z CAP edge (-17.5) flux: min={fminus[m].min():+.4f}  max={fminus[m].max():+.4f}  "
      f"(negative => backward/reflection)")
print(f"  near face  (-12.5) flux: min={flux[-12.5][m].min():+.4f}  max={flux[-12.5][m].max():+.4f}")
for zp, lab, _ in PLANES:
    print(f"  cum flux across {lab:>22}: t_wrap={cum[zp][m][-1]:+.3f}  final={cum[zp][-1]:+.3f}")

# --- plot ---
fig, (axA, axB) = plt.subplots(2, 1, figsize=(7.0, 7.6))


def mark(ax):
    for te, _ in EVENTS:
        ax.axvline(te, ls=":", lw=0.9, color="0.55", alpha=0.8)
    ax.axvspan(T_WRAP, t[-1], color="0.9", alpha=0.5, zorder=0)


# (A) cumulative flux across each plane
mark(axA)
for zp, lab, c in PLANES:
    axA.plot(t, cum[zp], c, lw=1.8, label=f"{lab}: {cum[zp][-1]:+.3f}")
axA.axhline(0, ls="-", lw=0.8, color="0.5")
ymin = min(cum[zp].min() for zp, _, _ in PLANES)
for te, lbl in EVENTS:
    axA.text(te, ymin, " " + lbl, rotation=90, va="bottom", ha="right",
             fontsize=6.5, color="0.4")
axA.text(T_WRAP + 0.1, axA.get_ylim()[1], " wrap-contaminated", fontsize=6.5,
         color="0.5", va="top")
axA.set_ylabel("cumulative ∫J_z dt across plane  (electrons; + forward)")
axA.set_title("Cumulative WP current across every internal plane — direction of transport",
              fontsize=9)
axA.legend(fontsize=7, frameon=False, loc="upper left")
axA.grid(alpha=0.25)

# (B) instantaneous backward-side flux, early window
mark(axB)
axB.axhline(0, ls="-", lw=0.8, color="0.5")
axB.plot(t, flux[-12.5], "C1-", lw=1.8, label="near face (-12.5)")
axB.plot(t, flux[-17.5], "C3-", lw=1.8, label="-z CAP edge (-17.5)")
axB.set_xlim(0, T_WRAP)
axB.set_xlabel("time (a.u.)")
axB.set_ylabel("instantaneous flux J_z  (+ forward / − backward)")
axB.set_title("Backward-side flux (wrap-free): a NEGATIVE excursion = genuine reflection",
              fontsize=9)
axB.legend(fontsize=7, frameon=False, loc="upper right")
axB.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(f"{OUT}/qa_ix_cumulative_current_regions.png", dpi=200)
plt.close(fig)

df = pd.DataFrame({"time_au": t})
for zp, lab, _ in PLANES:
    df[f"flux_{zp}"] = flux[zp]
    df[f"cum_{zp}"] = cum[zp]
df.to_csv(f"{OUT}/qa_ix_cumulative_current_regions.csv", index=False)
print("wrote qa_ix_cumulative_current_regions.png + .csv")
