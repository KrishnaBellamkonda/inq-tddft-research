#!/usr/bin/env python3
"""Loss function on LINEAR colour scale + 1D cuts — channel-excitation check.

The pipeline's loss_function_qz_omega.png uses SymLogNorm, which compresses the
dynamic range and can make every (q,ω) channel look equally bright. This re-plots
the SAME L(q_z,ω) on a LINEAR scale (and a percentile-clipped linear scale), plus
1D cuts L(ω) at fixed q, so the user can judge whether channels are genuinely
uniformly excited or whether structure was hidden by the log scale.

Reuses inqview.postprocess.spectral_weight internals (no recomputation drift).
Output (in this run's results/analysis/observables/):
  loss_function_linear.png        L(q_z,ω) linear + clipped-linear (2 panels)
  loss_function_1d_cuts.png       L(ω) at several fixed q_z
  spectral_weight_response_linear.png  W_resp linear (WP-subtracted weight)
"""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from pathlib import Path
from inqview.postprocess import spectral_weight as sw

HA_TO_EV = 27.211386245988
RUN = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_wp_n162_L50_E20_sigma1_v2")
RES = RUN / "results"
OUT = RES / "analysis" / "observables"
OUT.mkdir(parents=True, exist_ok=True)

# ---- replicate run() up to L_qz_w (omega_max/q_max generous for inspection) ----
omega_max_eV, q_max_inv_bohr, zero_pad_factor, eta_factor = 20.0, 1.0, 16, 3.0
rs = sw._parse_run_summary(RES)
import glob, re
dd = RES / "raw" / "vti" / "density_rt_delta"
if not dd.exists(): dd = RES / "raw" / "vti" / "density_delta"
vti = sorted(dd.glob("*.vti"))
dt_au = float(rs.get("dt_au", 0.01))
steps = np.array(sorted(int(re.search(r"t(\d+)", f.stem).group(1)) for f in vti))
times = steps * dt_au; t_rel = times - times[0]; T = t_rel[-1]
L = float(re.search(r"([\d.]+)", str(rs.get("cell_bohr", "50"))).group(1))
dx = float(rs.get("spacing_bohr", 0.40)); Ng = int(round(L / dx))
Ne = int(rs.get("n_electrons", 162)); n0 = Ne / L**3
omega_p = np.sqrt(4 * np.pi * n0); k_F = (3 * np.pi**2 * n0)**(1/3)
sigma = float(rs.get("wp_sigma_bohr", 1.0))
k0p = str(rs.get("wp_k0_bohr_inv", "0 0 1.21")).split()
k0_z = float(k0p[-1]); v0 = abs(k0_z); vz = 1.0 if k0_z >= 0 else -1.0
wc = str(rs.get("wp_center_bohr", "0 0 -21")).split()
x_c0, y_c0, z_c0 = float(wc[0]), float(wc[1]), float(wc[-1])
print(f"L={L} dx={dx} Ng={Ng} omega_p={omega_p*HA_TO_EV:.2f} eV k_F={k_F:.3f} v0={v0:.3f} sigma={sigma}")

# determine true grid from the first frame (density_rt_delta may be Ng-1)
Ng = sw._load_vti_as_array(vti[0]).shape[0]
print(f"true grid from VTI = {Ng}")
Nf = len(vti); Nq = Ng
dn_tot = np.zeros((Nq, Nf), complex); dn_wp = np.zeros((Nq, Nf), complex)
for j, p in enumerate(vti):
    dn = sw._load_vti_as_array(p)
    dn_tot[:, j] = np.fft.fftn(dn)[0, 0, :]
    t = t_rel[j]; zc = z_c0 + vz * v0 * t
    w_t = sw._free_wp_density_on_grid(Ng, L, dx, sigma, zc, x_c0, y_c0, t)
    w_0 = sw._free_wp_density_on_grid(Ng, L, dx, sigma, z_c0, x_c0, y_c0, 0.0)
    dn_wp[:, j] = np.fft.fftn(w_t - w_0)[0, 0, :]
dn_resp = dn_tot - dn_wp

No = Nf * zero_pad_factor; hann = np.hanning(Nf)
def tfft(d):
    o = np.zeros((Nq, No), complex)
    for i in range(Nq):
        s = d[i].copy() - d[i].mean(); o[i] = np.fft.fft(s * hann, n=No)
    return o
dn_resp_w = tfft(dn_resp); dn_wp_w = tfft(dn_wp)
W_resp = np.abs(dn_resp_w)**2
qz = np.fft.fftfreq(Nq, d=dx) * 2 * np.pi
q2 = np.where(qz[:, None]**2 > 1e-10, qz[:, None]**2, 1e-10)
Vext = -(4 * np.pi / q2) * dn_wp_w
# Zero the unphysical q=0 row: its 4π/q² is only finite via the 1e-10 clamp, which
# inflates |Vext| ~10^3x and poisons thr below → all physical modes rejected → L≡0.
Vext[np.abs(qz) < 1e-5, :] = 0.0
thr = 1e-3 * np.max(np.abs(Vext))
chi = np.zeros_like(dn_resp_w); m = np.abs(Vext) > thr
chi[m] = dn_resp_w[m] / Vext[m]
Lqw = -(4 * np.pi / q2) * np.imag(chi)

omega = np.fft.fftfreq(No, d=(times[1] - times[0])) * 2 * np.pi
qpos = qz[1:Nq // 2]; opos = omega[:No // 2]; oeV = opos * HA_TO_EV
qm = qpos <= q_max_inv_bohr; om = oeV <= omega_max_eV
qpl = qpos[qm]; opl = oeV[om]
def clip(a): return a[1:Nq // 2, :No // 2][np.ix_(qm, om)]
Lp = clip(Lqw); Wp = clip(W_resp)
qref = np.linspace(0.01, q_max_inv_bohr, 200)
opl_pl, opl_plus, opl_minus, opl_kin = sw._reference_curves(qref, omega_p, k_F, v0)
cur = [c * HA_TO_EV for c in (opl_pl, opl_plus, opl_minus, opl_kin)]

def overlay(ax):
    ax.plot(qref, cur[0], "c-", lw=1.5, label=r"plasmon $\omega_{pl}$")
    ax.plot(qref, cur[1], "w--", lw=0.9, alpha=0.7)
    ax.plot(qref, cur[2], "w--", lw=0.9, alpha=0.7, label=r"P-H edges")
    ax.plot(qref, cur[3], "lime", lw=1.3, ls=":", label=r"$\omega=qv_0$")
    ax.set_xlim(0, q_max_inv_bohr); ax.set_ylim(0, omega_max_eV)
    ax.set_xlabel(r"$q_z$ (Bohr$^{-1}$)"); ax.set_ylabel(r"$\omega$ (eV)")

# ---- linear maps ----
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
im1 = a1.pcolormesh(qpl, opl, Lp.T, shading="auto", cmap="inferno",
                    vmin=0, vmax=np.nanmax(Lp))
fig.colorbar(im1, ax=a1, label=r"$L$"); overlay(a1)
a1.set_title("L(q,ω) LINEAR (full range)"); a1.legend(fontsize=7, loc="upper left")
vmax99 = np.nanpercentile(Lp[Lp > 0], 99) if np.any(Lp > 0) else 1
im2 = a2.pcolormesh(qpl, opl, Lp.T, shading="auto", cmap="inferno", vmin=0, vmax=vmax99)
fig.colorbar(im2, ax=a2, label=r"$L$"); overlay(a2)
a2.set_title("L(q,ω) LINEAR (clipped at 99th pct)"); a2.legend(fontsize=7, loc="upper left")
fig.tight_layout(); fig.savefig(OUT / "loss_function_linear.png", dpi=150); plt.close(fig)
print("wrote loss_function_linear.png")

# ---- W_resp linear ----
fig, ax = plt.subplots(figsize=(7, 5))
im = ax.pcolormesh(qpl, opl, Wp.T, shading="auto", cmap="inferno",
                   vmin=0, vmax=np.nanpercentile(Wp[Wp > 0], 99) if np.any(Wp > 0) else 1)
fig.colorbar(im, ax=ax, label=r"$|\delta n_{resp}|^2$"); overlay(ax)
ax.set_title("Response weight (linear, clipped)"); ax.legend(fontsize=7, loc="upper left")
fig.tight_layout(); fig.savefig(OUT / "spectral_weight_response_linear.png", dpi=150); plt.close(fig)
print("wrote spectral_weight_response_linear.png")

# ---- 1D cuts at fixed q ----
fig, ax = plt.subplots(figsize=(8, 5))
qcuts = [0.2, 0.4, 0.6, 0.8]
for qc in qcuts:
    iq = int(np.argmin(np.abs(qpl - qc)))
    ax.plot(opl, Lp[iq], lw=1.4, label=f"$q_z$={qpl[iq]:.2f}")
ax.axvline(omega_p * HA_TO_EV, color="c", ls="--", lw=1, label=r"$\omega_p$=%.1f eV" % (omega_p*HA_TO_EV))
ax.set_xlabel(r"$\omega$ (eV)"); ax.set_ylabel(r"$L(q_z,\omega)$")
ax.set_title("L(ω) cuts at fixed $q_z$ (linear) — channel structure")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(OUT / "loss_function_1d_cuts.png", dpi=150); plt.close(fig)
print("wrote loss_function_1d_cuts.png")

# numeric summary: peak location per q-cut
print("=== L(ω) peak per q-cut ===")
for qc in qcuts:
    iq = int(np.argmin(np.abs(qpl - qc)))
    ip = int(np.argmax(Lp[iq]));
    print(f"  q={qpl[iq]:.2f}: peak L={Lp[iq][ip]:.3e} at ω={opl[ip]:.2f} eV  "
          f"(ω_p={omega_p*HA_TO_EV:.2f}, ω=qv0={qpl[iq]*v0*HA_TO_EV:.2f})")
print(f"L range: min={np.nanmin(Lp):.2e} max={np.nanmax(Lp):.2e} "
      f"dynamic-range={np.nanmax(np.abs(Lp))/max(np.nanpercentile(np.abs(Lp),50),1e-30):.1f}x median")
