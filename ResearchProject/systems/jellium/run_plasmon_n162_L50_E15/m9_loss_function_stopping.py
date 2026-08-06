#!/usr/bin/env python3
"""M9 — loss-function stopping power S(v) from the E15 long-run loss function.

Uses the SAME loss function as report1 draft5 (`make_fig_loss_function.py`):
  L(q,ω) = |FFT(n_q(t)·hann)|² / q²   (positive-definite density power spectrum / q²)
from run_plasmon_n162_L50_E15 (T=2000 a.u.=48.4 fs, Δω≈0.09 eV — properly resolved).

A loss function is a property of the MEDIUM, so one L(q,ω) yields S(v) for ANY v via
the dielectric stopping formula with the kinematic limit ω≤qv:
    S(v) = (2/πv²) ∫_qmin^qmax (dq/q) ∫_0^{qv} ω L(q,ω) dω      [Lindhard 1954; Ritchie 1959]

We compare the E15-loss-function S(v) (arb. units — power spectrum, not Im[-1/ε], so
NOT absolute without FDT calibration) against:
  • analytic Lindhard/RPA S(v)  (absolute Ha/Bohr; lindhard.stopping_power), box-matched
    q-window AND full kinematic q-window,
  • classical Ehrenfest S(v) points from the L50 classical runs (absolute eV/Bohr).
The E15 curve is normalised to the box-matched analytic curve by total area over the v
range (shape/trend comparison — the 'usable metric' test). All caveats printed.

Known-case test (printed): m=1 peak must sit near ω_p≈3.47 eV and disperse up ~Bohm-Gross.

Run: venv/bin/python3 <thisdir>/m9_loss_function_stopping.py
"""
from __future__ import annotations
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
from pathlib import Path
from inqview.pipeline import lindhard

HA = 27.211386245988
ROOT = Path("/local/data/public/skcb2/tddft")
E15 = ROOT / "ResearchProject/systems/jellium/run_plasmon_n162_L50_E15"
NQ = E15 / "results/analysis/observables/n_q_vs_time.csv"
OUT = E15 / "results/analysis/observables"
L_BOHR = 50.0
N_DENS = 162 / L_BOHR**3
kF = (3 * np.pi**2 * N_DENS) ** (1 / 3)
omega_p = np.sqrt(4 * np.pi * N_DENS)
vF = kF
v0_E15 = 1.04999  # wp_k0_bohr_inv z-component (E=15 eV)
print(f"system: r_s={(3/(4*np.pi*N_DENS))**(1/3):.2f}  k_F={kF:.4f}  "
      f"omega_p={omega_p*HA:.2f} eV  v0(E15)={v0_E15:.3f} a.u.")


def bohm_gross(q):
    return np.sqrt(omega_p**2 + 0.6 * vF**2 * q**2)


# ---- Stage 1: L(q,ω) per mode (draft5 method) ----
df = pd.read_csv(NQ)
modes = sorted(df["m"].unique())
q_of_mode, omega_au, L_modes = [], None, []
print("\n=== known-case: loss-function peak per mode vs Bohm-Gross ===")
for m in modes:
    sub = df[df["m"] == m].sort_values("time_au")
    t = sub["time_au"].values
    nq = sub["re_n_q"].values + 1j * sub["im_n_q"].values
    nq = nq - nq.mean()                            # remove static screening (DC) component
    N = len(t); q = sub["q_au"].values[0]
    fft = np.fft.fft(nq * np.hanning(N))
    fr = np.fft.fftfreq(N, d=t[1] - t[0])
    pos = fr >= 0
    om = fr[pos] * 2 * np.pi                       # a.u. (Ha)
    L = (np.abs(fft[pos]) ** 2) / q**2             # arb. units
    if omega_au is None:
        omega_au = om
    q_of_mode.append(q); L_modes.append(L)
    ipk = int(np.argmax(L[1:]) + 1)                # skip ω=0 bin
    iwl = int(np.argmax((om * L)[1:]) + 1)         # ω·L peak (draft5's plotted quantity)
    print(f"  m={m} q={q:.3f}: L-peak ω={om[ipk]*HA:6.2f}  ωL-peak ω={om[iwl]*HA:6.2f}  "
          f"Bohm-Gross={bohm_gross(q)*HA:6.2f}  ω_p={omega_p*HA:.2f} eV")
q_modes = np.array(q_of_mode)
# L(q,ω) matrix on the shared ω grid (modes share N, dt → same grid)
Lqw = np.vstack(L_modes)                            # (n_q, n_omega), arb units


# physical ω support cap (eV→a.u.): plasmon ~3.5 eV + e-h continuum top ~15 eV all lie
# below this; it excludes the dt=4 a.u. Nyquist-edge artifact (~19-21 eV) in high-q modes.
# The analytic Im[-1/ε] is already 0 above the e-h top, so this only affects the noisy
# numerical E15 L → comparison stays fair.
OMEGA_PHYS_CAP = 16.0 / HA


def stopping_from_L(v, q_arr, om_arr, Lmat):
    """S(v) = (2/πv²) ∫(dq/q) ∫_0^{min(qv, cap)} ω L dω, trapz on the mode q-grid."""
    g_over_q = np.zeros_like(q_arr)
    for i, q in enumerate(q_arr):
        wmax = min(q * v, OMEGA_PHYS_CAP)
        sel = (om_arr > 0) & (om_arr <= wmax)
        if sel.sum() < 2:
            continue
        inner = np.trapezoid(om_arr[sel] * Lmat[i, sel], om_arr[sel])  # ∫ ω L dω
        g_over_q[i] = inner / q                                    # (1/q)·g(q)
    return (2.0 / (np.pi * v**2)) * np.trapezoid(g_over_q, q_arr)      # ∫(g/q)dq


# ---- Stage 2: S(v) curves ----
v_grid = np.linspace(0.2, 4.0, 60)
q_min_box, q_max_mode = q_modes[0], q_modes[-1]
S_LF = np.array([stopping_from_L(v, q_modes, omega_au, Lqw) for v in v_grid])  # arb
S_lin_box = np.array([lindhard.stopping_power(v, kF, qmin=q_min_box, qmax=q_max_mode)
                      for v in v_grid]) * HA  # eV/Bohr, box-matched q window
S_lin_full = np.array([lindhard.stopping_power(v, kF) for v in v_grid]) * HA   # full kinematic

# normalise the arb-unit LF curve to the box-matched analytic curve by total area
mask = S_lin_box > 0
scale = np.trapezoid(S_lin_box[mask], v_grid[mask]) / max(np.trapezoid(S_LF[mask], v_grid[mask]), 1e-30)
S_LF_norm = S_LF * scale

# ---- Stage 3: classical S(v) points (absolute eV/Bohr) ----
cls_v, cls_S = [], []
try:
    from applications.report1 import stopping_power_data as spd
    for spec in spd.get_L50_classical_runs():
        try:
            p = ROOT / spec.run_dir; params = spd.parse_run_summary(p)
            win = spd.compute_time_window(params)
            S, _ = spd.compute_classical_S(p, win)            # already eV/Bohr
            v = np.sqrt(2 * spec.energy_eV / HA)
            if np.isfinite(S):
                cls_v.append(v); cls_S.append(abs(S))
                print(f"  classical E={spec.energy_eV}eV v={v:.3f}: S={abs(S):.4f} eV/Bohr")
        except Exception as e:
            print(f"  classical {spec.run_dir}: skip ({e})")
except Exception as e:
    print(f"classical data unavailable: {e}")

# ---- Stage 4: figure ----
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(v_grid, S_lin_full, "-", color="#1f77b4", lw=2,
        label="analytic Lindhard (full kinematic q)")
ax.plot(v_grid, S_lin_box, "--", color="#1f77b4", lw=1.5,
        label=f"analytic Lindhard (box q: {q_min_box:.2f}–{q_max_mode:.2f})")
ax.plot(v_grid, S_LF_norm, "-", color="#d62728", lw=2,
        label="E15 loss-function S(v) (arb, area-normalised to box)")
if cls_v:
    ax.plot(cls_v, cls_S, "ks", ms=7, label="classical Ehrenfest runs")
ax.axvline(v0_E15, color="green", ls=":", lw=1.2, label=f"E15 WP v₀={v0_E15:.2f}")
ax.set_xlabel("v (a.u.)"); ax.set_ylabel("S(v) (eV/Bohr)")
ax.set_title("M9 — loss-function stopping power vs analytic & classical")
ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.set_xlim(0, 4); ax.set_ylim(bottom=0)
fig.tight_layout()
OUT.mkdir(parents=True, exist_ok=True)
fp = OUT / "m9_loss_function_stopping.png"
fig.savefig(fp, dpi=150); plt.close(fig)
print(f"\nwrote {fp}")

# numeric summary
iv = int(np.argmin(np.abs(v_grid - v0_E15)))
print(f"\n=== S(v) at E15 v0={v0_E15:.2f} ===")
print(f"  analytic Lindhard (full)   = {S_lin_full[iv]:.4f} eV/Bohr")
print(f"  analytic Lindhard (box q)  = {S_lin_box[iv]:.4f} eV/Bohr")
print(f"  E15 loss-fn (normalised)   = {S_LF_norm[iv]:.4f} eV/Bohr  (scale={scale:.3e})")
iv_pk_lf = int(np.argmax(S_LF_norm)); iv_pk_lin = int(np.argmax(S_lin_full))
print(f"  peak-v: loss-fn={v_grid[iv_pk_lf]:.2f}  analytic-full={v_grid[iv_pk_lin]:.2f} a.u.")
