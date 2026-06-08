#!/usr/bin/env python3
"""Projectile loss function with the WRAP line and the PLASMON separated.

A projectile in a periodic box imprints a kinematic/wrap line at
ω_kin(m) = m·v·q₁ (q₁=2π/L) in its loss function L(q,ω). This only OVERLAPS the
plasmon ω_p when v ≈ v_res (the E15 resonance). Off-resonance projectiles
(E3.4: wrap below; E25: wrap above) leave the plasmon peak clean.

This makes, per available projectile run, the L(q,ω) map + the m=1 1D cut, with
the Bohm-Gross plasmon ω(q) and the kinematic line ω_kin overlaid and the
plasmon vs wrap peaks labelled — so the clean (E3.4, E25) vs contaminated (E15)
cases are visible. NO kick; these are genuine PROJECTILE loss functions.

Output: batch2_figures/projectile_loss_clean.png
Known-case (printed): per run, m=1 peak ω, ω_kin, |ω_kin-ω_p|.
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
kF = (3 * np.pi**2 * N_DENS) ** (1 / 3); omega_p = np.sqrt(4 * np.pi * N_DENS); vF = kF
q1 = 2 * np.pi / L_BOHR


def bohm_gross(q):
    return np.sqrt(omega_p**2 + 0.6 * vF**2 * q**2 + 0.25 * q**4)


def loss_qw(csv):
    df = pd.read_csv(csv); modes = sorted(df["m"].unique())
    qs, om0, Ls = [], None, []
    for m in modes:
        sub = df[df["m"] == m].sort_values("time_au"); t = sub["time_au"].values
        nq = sub["re_n_q"].values + 1j * sub["im_n_q"].values
        nq = nq - nq.mean(); Nn = len(t); q = sub["q_au"].values[0]
        fft = np.fft.fft(nq * np.hanning(Nn)); fr = np.fft.fftfreq(Nn, d=t[1] - t[0]); pos = fr >= 0
        om = fr[pos] * 2 * np.pi; L = (np.abs(fft[pos]) ** 2) / q**2
        if om0 is None: om0 = om
        qs.append(q); Ls.append(L)
    return np.array(qs), om0, np.vstack(Ls)


RUNS = [("E3.4 (v=0.50) — wrap BELOW plasmon", 0.500, JB / "run_plasmon_n162_L50_E3p4_varyv"),
        ("E15 (v=1.05) — RESONANT, wrap≈plasmon", 1.050, JB / "run_plasmon_n162_L50_E15")]
e25 = JB / "run_plasmon_n162_L50_E25/results/analysis/observables/n_q_vs_time.csv"
if e25.exists():
    RUNS.insert(1, ("E25 (v=1.36) — wrap ABOVE plasmon", 1.356, e25.parent.parent.parent.parent))

avail = [(lab, v, r) for lab, v, r in RUNS
         if (r / "results/analysis/observables/n_q_vs_time.csv").exists()]
fig, axes = plt.subplots(2, len(avail), figsize=(5.2 * len(avail), 8.4), squeeze=False)
omax = 9.0
print(f"ω_p={omega_p*HA:.2f} eV  q1={q1:.4f}")
for j, (lab, v, run) in enumerate(avail):
    q, om, L = loss_qw(run / "results/analysis/observables/n_q_vs_time.csv")
    wkin1 = v * q1 * HA
    i1 = 0; ipk = int(np.argmax(L[i1, 1:]) + 1)
    print(f"  {lab}: m=1 peak ω={om[ipk]*HA:.2f} eV | ω_kin(m1)={wkin1:.2f} | "
          f"|ω_kin-ω_p|={abs(wkin1-omega_p*HA):.2f} eV")
    # top: 2D map
    ax = axes[0][j]; sel = (om * HA <= omax) & (om > 0)
    Q, W = np.meshgrid(q, om[sel] * HA, indexing="ij"); Lp = L[:, sel]
    vmax = np.percentile(Lp[Lp > 0], 99.5)
    pc = ax.pcolormesh(Q, W, Lp, shading="auto", cmap="inferno", norm=LogNorm(vmax * 1e-4, vmax))
    qq = np.linspace(q.min(), q.max(), 80)
    ax.plot(qq, bohm_gross(qq) * HA, "c-", lw=1.6, label="plasmon ω(q) (Bohm-Gross)")
    ax.plot(q, np.arange(1, len(q) + 1) * v * q1 * HA, "g--", lw=1.5, label="wrap ω_kin=m·v·q₁")
    ax.axhline(omega_p * HA, color="w", ls=":", lw=1.0)
    ax.set_title(lab, fontsize=9.5); ax.set_xlabel("q (Bohr⁻¹)"); ax.set_ylim(0, omax)
    ax.legend(fontsize=6.5, loc="upper left"); fig.colorbar(pc, ax=ax, shrink=0.8)
    # bottom: m=1 1D cut with peaks labelled
    ax2 = axes[1][j]; s = (om * HA <= omax) & (om > 0)
    ax2.plot(om[s] * HA, om[s] * L[i1, s], "k-", lw=1.4)
    ax2.axvline(omega_p * HA, color="c", lw=1.6, label=f"plasmon ω_p={omega_p*HA:.2f}")
    ax2.axvline(wkin1, color="g", ls="--", lw=1.6, label=f"wrap ω_kin={wkin1:.2f}")
    ax2.set_xlabel("ω (eV)"); ax2.set_ylabel("ω·L(ω), m=1"); ax2.set_xlim(0, omax)
    ax2.legend(fontsize=7.5); ax2.grid(alpha=0.3)
axes[0][0].set_ylabel("ω (eV)")
fig.suptitle("PROJECTILE loss function: wrap line (green) vs plasmon (cyan). "
             "Off-resonance (E3.4, E25) → plasmon clean; E15 (resonant) → degenerate.",
             fontsize=11)
fig.tight_layout()
fp = OUT / "projectile_loss_clean.png"
fig.savefig(fp, dpi=150); plt.close(fig)
print(f"\nwrote {fp}")
