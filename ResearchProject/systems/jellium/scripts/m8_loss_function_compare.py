#!/usr/bin/env python3
"""M8 — is the loss function L(q,ω) a projectile-INDEPENDENT medium property?

Extracts L(q,ω) = |FFT(n_q·hann)|²/q² (the draft5 / m9 method) from the two
fully-resolved long runs at different projectile velocities:
  * E15  (15 eV, v=1.05 a.u.)  run_plasmon_n162_L50_E15
  * E3.4 (3.4 eV, v=0.50 a.u.) run_plasmon_n162_L50_E3p4_varyv
both T=2000 a.u. (Δω=0.086 eV).  A 25 eV companion (run_plasmon_n162_L50_E25)
is added automatically when its n_q_vs_time.csv exists.

The MEDIUM plasmon peak ω_BG(q) (Bohm-Gross) is velocity-independent; a
KINEMATIC peak ω_kin(m)=m·v·q₁ scales with v.  Overlaying the runs per mode
therefore separates medium (coincide) from projectile (shift) features.

Optionally overlays the classical E100 loss function (T≈33 a.u., Δω≈5.8 eV —
resolution-starved, plotted only as a coarse cross-check with that caveat).

Outputs (batch2_figures/):
  m8_loss_function_2d_compare.png   L(q,ω) maps side-by-side + Bohm-Gross/ω_p/ω_kin
  m8_loss_function_1d_overlay.png   per-mode ωL(ω) overlays (area-normalised)

Known-case (printed): each run's m=1 peak vs ω_p and its kinematic ω_kin.
"""
from __future__ import annotations
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HA = 27.211386245988
JB = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium")
OUT = Path("/local/data/public/skcb2/tddft/docs/presentations/storyline/tasks/batch2_figures")
L_BOHR, N = 50.0, 162
N_DENS = N / L_BOHR**3
kF = (3 * np.pi**2 * N_DENS) ** (1 / 3)
omega_p = np.sqrt(4 * np.pi * N_DENS)
q1 = 2 * np.pi / L_BOHR
vF = kF


def bohm_gross(q):
    return np.sqrt(omega_p**2 + 0.6 * vF**2 * q**2 + 0.25 * q**4)


def loss_qw(csv):
    """Return (q_modes, omega_au, Lqw[n_q,n_omega]) via |FFT(n_q·hann)|²/q²."""
    df = pd.read_csv(csv)
    modes = sorted(df["m"].unique())
    qs, om0, Ls = [], None, []
    for m in modes:
        sub = df[df["m"] == m].sort_values("time_au")
        t = sub["time_au"].values
        nq = sub["re_n_q"].values + 1j * sub["im_n_q"].values
        nq = nq - nq.mean()
        Nn = len(t); q = sub["q_au"].values[0]
        fft = np.fft.fft(nq * np.hanning(Nn))
        fr = np.fft.fftfreq(Nn, d=t[1] - t[0]); pos = fr >= 0
        om = fr[pos] * 2 * np.pi
        L = (np.abs(fft[pos]) ** 2) / q**2
        if om0 is None:
            om0 = om
        qs.append(q); Ls.append(np.interp(om0, om, L) if om0 is not om else L)
    return np.array(qs), om0, np.vstack(Ls)


RUNS = [
    ("E15 (15 eV, v=1.05)", 1.050, JB / "run_plasmon_n162_L50_E15"),
    ("E3.4 (3.4 eV, v=0.50)", 0.500, JB / "run_plasmon_n162_L50_E3p4_varyv"),
]
e25 = JB / "run_plasmon_n162_L50_E25/results/analysis/observables/n_q_vs_time.csv"
if e25.exists():
    RUNS.append(("E25 (25 eV, v=1.36)", 1.356, e25.parent.parent.parent.parent))

data = []
print(f"system r_s={(3/(4*np.pi*N_DENS))**(1/3):.2f}  omega_p={omega_p*HA:.2f} eV  q1={q1:.4f}")
print("\n=== known-case: m=1 peak vs medium ω_p and kinematic ω_kin=v·q1 ===")
for label, v, run in RUNS:
    csv = run / "results/analysis/observables/n_q_vs_time.csv"
    if not csv.exists():
        print(f"  {label}: NO n_q csv, skip"); continue
    q, om, L = loss_qw(csv)
    data.append((label, v, q, om, L))
    i1 = 0; ipk = int(np.argmax(L[i1, 1:]) + 1)
    print(f"  {label}: m=1 L-peak ω={om[ipk]*HA:5.2f} eV | ω_p={omega_p*HA:.2f} | "
          f"ω_kin=v·q1={v*q1*HA:5.2f} eV | Bohm-Gross={bohm_gross(q[i1])*HA:5.2f}")

# ---- Fig 1: 2D maps side by side ----
nrun = len(data)
fig, axs = plt.subplots(1, nrun, figsize=(5.2 * nrun, 4.6), squeeze=False)
omax_ev = 12.0
for ax, (label, v, q, om, L) in zip(axs[0], data):
    sel = (om * HA <= omax_ev) & (om > 0)
    Q, W = np.meshgrid(q, om[sel] * HA, indexing="ij")
    Lp = L[:, sel]
    vmax = np.percentile(Lp[Lp > 0], 99.5)
    pc = ax.pcolormesh(Q, W, Lp, shading="auto", cmap="inferno",
                       norm=LogNorm(vmax * 1e-4, vmax))
    qq = np.linspace(q.min(), q.max(), 100)
    ax.plot(qq, bohm_gross(qq) * HA, "c-", lw=1.5, label="Bohm-Gross ω(q)")
    ax.axhline(omega_p * HA, color="w", ls=":", lw=1.0, label=f"ω_p={omega_p*HA:.2f} eV")
    ax.plot(q, np.arange(1, len(q) + 1) * v * q1 * HA, "g--", lw=1.2, label="ω_kin=m·v·q₁")
    ax.set_title(label); ax.set_xlabel("q (Bohr⁻¹)"); ax.set_ylim(0, omax_ev)
    fig.colorbar(pc, ax=ax, shrink=0.8)
axs[0][0].set_ylabel("ω (eV)"); axs[0][-1].legend(fontsize=7, loc="upper left")
fig.suptitle("M8 — loss function L(q,ω): does it depend on the projectile? "
             "(plasmon = medium → should coincide; kinematic → shifts with v)")
fig.tight_layout()
fig.savefig(OUT / "m8_loss_function_2d_compare.png", dpi=150); plt.close(fig)
print(f"\nwrote {OUT / 'm8_loss_function_2d_compare.png'}")

# ---- Fig 2: per-mode 1D overlay (area-normalised ωL) ----
nq = min(len(d[2]) for d in data)
modes_show = list(range(min(4, nq)))
fig, axs = plt.subplots(1, len(modes_show), figsize=(4 * len(modes_show), 4),
                        sharex=True, squeeze=False)
for j, mi in enumerate(modes_show):
    ax = axs[0][j]
    for label, v, q, om, L in data:
        sel = (om * HA <= 12) & (om > 0)
        wL = om[sel] * L[mi, sel]
        area = np.trapezoid(wL, om[sel] * HA)
        ax.plot(om[sel] * HA, wL / max(area, 1e-30), lw=1.6, label=label)
    qm = data[0][2][mi]
    ax.axvline(bohm_gross(qm) * HA, color="c", ls="-", lw=1.0)
    ax.axvline(omega_p * HA, color="k", ls=":", lw=0.9)
    ax.set_title(f"m={mi+1}  q={qm:.3f}"); ax.set_xlabel("ω (eV)"); ax.grid(alpha=0.3)
    ax.set_xlim(0, 12)
axs[0][0].set_ylabel("ω·L(ω)  (area-normalised)"); axs[0][-1].legend(fontsize=7)
fig.suptitle("M8 — per-mode loss-function shape across projectile velocities "
             "(cyan = Bohm-Gross ω(q); dotted = ω_p)")
fig.tight_layout()
fig.savefig(OUT / "m8_loss_function_1d_overlay.png", dpi=150); plt.close(fig)
print(f"wrote {OUT / 'm8_loss_function_1d_overlay.png'}")
