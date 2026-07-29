#!/usr/bin/env python3
"""T-a — coronene MOLECULAR response map S(q_z,ω) (labelled; NOT a dielectric loss fn).

From the molecular (WP-excluded) electron density density_rt_system of the WP-scattering
run. For each frame: z-line density n(z,t)=∫∫ n_system dx dy; response δn(z,t)=n(z,t)-n(z,0);
FFT in z → δn(q_z,t); FFT in time (Hann) → δn(q_z,ω); map S=|δn(q_z,ω)|².

RESOLUTION CAVEAT: run T=16.96 a.u. → native Δω=2π/T≈10 eV. Sharp molecular lines are NOT
resolved; only BROAD electronic features survive (PAH π-plasmon ~5-6 eV, σ+π plasmon
~15-18 eV are several eV wide). A properly resolved map needs a long kick-response run
(cf. the jellium E15 2000-a.u. run). Labelled with the q=0 dipole-spectrum_z peaks.

Known-case (printed): system electron count ~const in t (WP excluded); δn(t0)=0.
Output: results/analysis/observables/ta_molecular_response_map.png
"""
import glob, re, numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path
import pandas as pd

HA = 27.211386245988
RUN = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/run_broadening_35x35x80")
SYS = RUN / "results/raw/vti/density_rt_system"
TOT = RUN / "results/raw/vti/density_rt_total"
OUT = RUN / "results/analysis/observables"
DT = 0.02
STRIDE = 2                       # subsample (preserves Δω; lowers Nyquist, still ample)


def frames(d):
    out = []
    for f in glob.glob(f"{d}/*.vti"):
        m = re.search(r"_t(\d+)\.vti$", f)
        if m:
            out.append((int(m.group(1)), f))
    return sorted(out)


def zline(path):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(path); r.Update()
    img = r.GetOutput(); nx, ny, nz = img.GetDimensions()
    oz, sz = img.GetOrigin()[2], img.GetSpacing()[2]
    sx, sy = img.GetSpacing()[0], img.GetSpacing()[1]
    a = vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx)
    return oz + sz * np.arange(nz), a.sum(axis=(1, 2)) * sx * sy, sz


fr = frames(SYS)[::STRIDE]
frt = {s: f for s, f in frames(TOT)}
steps = np.array([s for s, _ in fr]); times = steps * DT; T = times[-1] - times[0]
z, n0, dz = zline(fr[0][1])
Nz = len(z); Nf = len(fr)
print(f"T-a: {Nf} frames (stride {STRIDE}), Nz={Nz}, T={T:.2f} a.u., "
      f"native Δω=2π/T={2*np.pi/T*HA:.1f} eV")

# known-case: system count vs total count at t0 (system must exclude WP → count < total)
sys_e = n0.sum() * dz
ztt, ntt, _ = zline(frt[fr[0][0]]); tot_e = ntt.sum() * dz
print(f"[known-case] t0 system e={sys_e:.2f}, total e={tot_e:.2f} "
      f"(system<total ⇒ WP excluded: {sys_e < tot_e - 0.3})")

nz_t = np.zeros((Nz, Nf))
e_t = np.zeros(Nf)
for j, (s, f) in enumerate(fr):
    _, n, _ = zline(f); nz_t[:, j] = n; e_t[j] = n.sum() * dz
print(f"[known-case] system electron count drift over t: "
      f"min={e_t.min():.2f} max={e_t.max():.2f} (≈const ⇒ molecule only)")

dn = nz_t - nz_t[:, [0]]                       # response δn(z,t), δn(t0)=0 by construction
print(f"[known-case] max|δn(t0)| = {np.abs(dn[:, 0]).max():.2e} (==0)")

# FFT z then t
dn_qz_t = np.fft.fft(dn, axis=0)               # (q_z, t)
hann = np.hanning(Nf)
ZP = 8
No = Nf * ZP
dn_qz_w = np.fft.fft((dn_qz_t - dn_qz_t.mean(axis=1, keepdims=True)) * hann, n=No, axis=1)
S = np.abs(dn_qz_w) ** 2

qz = np.fft.fftfreq(Nz, d=dz) * 2 * np.pi
omega = np.fft.fftfreq(No, d=(times[1] - times[0])) * 2 * np.pi * HA   # eV
# positive halves, trim to physical window
qm = (qz > 0) & (qz <= 2.0)
om = (omega >= 0) & (omega <= 25.0)
Qp, Op = qz[qm], omega[om]
Sp = S[np.ix_(qm, om)]

# dipole_z peaks for labels (q=0 response)
dip = pd.read_csv(RUN / "results/raw/observables/dipole_spectrum_z.csv")
wcol, acol = dip.columns[0], dip.columns[1]
wev = dip[wcol].values * HA                      # frequency_au → eV
amp = dip[acol].values
# search for the strongest broad features in the molecular window 2–22 eV (skip DC)
win = (wev > 2.0) & (wev < 22.0)
amp_w = amp.copy(); amp_w[~win] = 0.0
pk_idx = [i for i in range(2, len(amp) - 2)
          if amp_w[i] > amp_w[i-1] and amp_w[i] > amp_w[i+1] and amp_w[i] > 0.25 * amp_w.max()]
pk_w = sorted(set(round(wev[i], 1) for i in pk_idx))[:6]

fig, ax = plt.subplots(figsize=(7, 5))
vmax = np.percentile(Sp[Sp > 0], 99.5)
im = ax.pcolormesh(Op, Qp, Sp, shading="auto", cmap="inferno",
                   norm=LogNorm(vmin=max(vmax * 1e-4, Sp[Sp > 0].min()), vmax=vmax))
fig.colorbar(im, ax=ax, label=r"$S(q_z,\omega)=|\delta n(q_z,\omega)|^2$ (arb.)")
for w in pk_w:
    ax.axvline(w, color="#66CCFF", ls=":", lw=0.9, alpha=0.8)
ax.set_xlabel(r"$\omega$ (eV)"); ax.set_ylabel(r"$q_z$ (Bohr$^{-1}$)")
ax.set_title("T-a — coronene molecular response map S(q$_z$,ω)  "
             f"(native Δω≈{2*np.pi/T*HA:.0f} eV: broad features only)")
ax.text(0.98, 0.02, "dotted = dipole-z (q=0) peaks", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=7, color="#66CCFF")
fig.tight_layout()
OUT.mkdir(parents=True, exist_ok=True)
fp = OUT / "ta_molecular_response_map.png"
fig.savefig(fp, dpi=150); plt.close(fig)
print(f"dipole-z peak labels (eV): {pk_w}")
print(f"wrote {fp}")
