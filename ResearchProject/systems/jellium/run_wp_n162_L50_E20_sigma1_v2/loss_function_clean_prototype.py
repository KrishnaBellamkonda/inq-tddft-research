#!/usr/bin/env python3
"""Prototype: regularised loss-function extraction (cleaner χ deconvolution).

The naive χ(q,ω) = δn_resp/V_ext divides by V_ext, which decays at high ω, so it
amplifies noise there (raw E20 L: 48% negative, 65% of weight piled at the ω-edge).
This prototype adds two grounded levers and sweeps them, reporting known-case
metrics so we can pick values that recover a physical (≥0, plasmon-peaked) L BEFORE
porting into inqview/postprocess/spectral_weight.py.

Levers (both standard, currently absent from the pipeline):
  (1) Lorentzian damping exp(-η t) before the temporal FFT — RT-TDDFT broadening
      convention (Yabana & Bertsch PRB 54, 4484 (1996); standard absorption-spectrum
      apodisation). η = eta_factor / T.
  (2) Tikhonov/Wiener regularised inversion:
          χ = δn_resp · conj(V_ext) / (|V_ext|² + λ²),   λ = reg_factor·max|V_ext_phys|
      → as |V_ext|→0 the estimate →0 instead of blowing up (Numerical Recipes §13.3;
      standard ill-posed-inverse regularisation).

Known-case acceptance targets for a physical loss function:
  - fraction of L<0 should drop well below the naive 48% (ideally <~10%),
  - weight should move OUT of the ω>8 eV noise band toward the plasmon ω_p,
  - the small-q peak should sit near ω_p (Bohm-Gross), not at the window edge.
"""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from pathlib import Path
from inqview.pipeline import spectral_weight as sw

HA_TO_EV = 27.211386245988
RUN = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_wp_n162_L50_E20_sigma1_v2")
RES = RUN / "results"
OUT = RES / "analysis" / "observables"
OUT.mkdir(parents=True, exist_ok=True)

omega_max_eV, q_max_inv_bohr = 20.0, 1.0
zero_pad_factor = 16
rs = sw._parse_run_summary(RES)
import re
dd = RES / "raw" / "vti" / "density_rt_delta"
if not dd.exists():
    dd = RES / "raw" / "vti" / "density_delta"
vti = sorted(dd.glob("*.vti"))
dt_au = float(rs.get("dt_au", 0.01))
steps = np.array(sorted(int(re.search(r"t(\d+)", f.stem).group(1)) for f in vti))
times = steps * dt_au; t_rel = times - times[0]; T = t_rel[-1]
L = float(re.search(r"([\d.]+)", str(rs.get("cell_bohr", "50"))).group(1))
dx = float(rs.get("spacing_bohr", 0.40))
Ne = int(rs.get("n_electrons", 162)); n0 = Ne / L**3
omega_p = np.sqrt(4 * np.pi * n0); k_F = (3 * np.pi**2 * n0)**(1/3)
sigma = float(rs.get("wp_sigma_bohr", 1.0))
k0p = str(rs.get("wp_k0_bohr_inv", "0 0 1.21")).split()
k0_z = float(k0p[-1]); v0 = abs(k0_z); vz = 1.0 if k0_z >= 0 else -1.0
wc = str(rs.get("wp_center_bohr", "0 0 -21")).split()
x_c0, y_c0, z_c0 = float(wc[0]), float(wc[1]), float(wc[-1])

# ---- load δn(q_z,t) and the analytic WP density once ----
Ng = sw._load_vti_as_array(vti[0]).shape[0]
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
print(f"loaded {Nf} frames, T={T:.2f} a.u., omega_p={omega_p*HA_TO_EV:.2f} eV, "
      f"k_F={k_F:.3f}, v0={v0:.3f}")

qz = np.fft.fftfreq(Nq, d=dx) * 2 * np.pi
q2 = np.where(qz[:, None]**2 > 1e-10, qz[:, None]**2, 1e-10)
No = Nf * zero_pad_factor
hann = np.hanning(Nf)
omega = np.fft.fftfreq(No, d=(times[1] - times[0])) * 2 * np.pi
qpos = qz[1:Nq // 2]; oeV = (omega[:No // 2]) * HA_TO_EV
qm = qpos <= q_max_inv_bohr; om = oeV <= omega_max_eV
qpl = qpos[qm]; opl = oeV[om]
def clip(a): return a[1:Nq // 2, :No // 2][np.ix_(qm, om)]


def tfft(d, eta):
    """temporal FFT with mean removal + Lorentzian damping exp(-eta t) + Hann."""
    damp = np.exp(-eta * t_rel) * hann
    o = np.zeros((Nq, No), complex)
    for i in range(Nq):
        s = d[i].copy() - d[i].mean()
        o[i] = np.fft.fft(s * damp, n=No)
    return o


def extract_L(eta_factor, reg_factor):
    eta = eta_factor / T
    dn_resp_w = tfft(dn_resp, eta); dn_wp_w = tfft(dn_wp, eta)
    Vext = -(4 * np.pi / q2) * dn_wp_w
    Vext[np.abs(qz) < 1e-5, :] = 0.0                      # drop unphysical q=0 row
    Vmax = np.max(np.abs(Vext))
    lam = reg_factor * Vmax                               # Tikhonov noise floor
    # Wiener/Tikhonov inversion: χ = δn_resp·conj(Vext)/(|Vext|²+λ²)
    chi = dn_resp_w * np.conj(Vext) / (np.abs(Vext)**2 + lam**2)
    Lqw = -(4 * np.pi / q2) * np.imag(chi)
    return clip(Lqw)


def metrics(Lp):
    olow = opl <= 8.0; ohigh = opl > 8.0
    wlow = float(np.sum(np.abs(Lp[:, olow]))); whigh = float(np.sum(np.abs(Lp[:, ohigh])))
    frac_neg = 100 * np.mean(Lp < 0)
    high_frac = 100 * whigh / max(wlow + whigh, 1e-30)
    # small-q peak omega
    iq = int(np.argmin(np.abs(qpl - 0.25)))
    ip = int(np.argmax(Lp[iq])); peak_w = opl[ip]
    return frac_neg, high_frac, peak_w, float(np.max(Lp)), float(np.min(Lp))


print(f"\n{'eta_f':>6} {'reg':>6} | {'%neg':>6} {'%w>8eV':>7} {'peakw@q.25':>11} "
      f"{'maxL':>9} {'minL':>9}")
print("-" * 62)
for ef in [3.0, 6.0, 10.0]:
    for rf in [1e-2, 3e-2, 1e-1, 3e-1]:
        Lp = extract_L(ef, rf)
        fn, hf, pw, mx, mn = metrics(Lp)
        print(f"{ef:6.1f} {rf:6.3f} | {fn:6.1f} {hf:7.1f} {pw:11.2f} {mx:9.2e} {mn:9.2e}")
print(f"\n(plasmon target omega_p = {omega_p*HA_TO_EV:.2f} eV; want %neg low, "
      f"%w>8eV low, peak near omega_p)")
