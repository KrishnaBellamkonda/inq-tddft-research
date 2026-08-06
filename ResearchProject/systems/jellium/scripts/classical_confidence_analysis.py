#!/usr/bin/env python3
"""Classical-electron momentum / scattering / loss-function analysis.

Two purposes (user, 2026-06-01):
  (1) CONFIDENCE in the classical (cusp-pseudopotential) runs: does the
      classical point projectile excite the SAME medium response (plasmon at
      ω_p, stopping ~ Lindhard) as expected?  A clean match builds confidence;
      a distortion would flag the cusp-infested pseudopotential approach.
  (2) WP-vs-classical contrast at the same energy (E=100 eV) to isolate what is
      genuinely QUANTUM: the classical projectile is a single deterministic
      trajectory (a δ in momentum); the WP carries a momentum DISTRIBUTION that
      broadens / scatters.

Panels (batch2_figures/classical_confidence_analysis.png):
  A  classical L(q,ω) map + Bohm-Gross + ω_p  (CAVEAT: classical run T≈33 a.u.
     → Δω≈5.8 eV, plasmon only marginally resolved — coarse cross-check).
  B  classical m=1 ω·L(ω) vs the resolved E15 medium m=1 (area-normalised):
     does the classical projectile peak sit at ω_p?
  C  classical projectile trajectory v_z(t) from electron_track.csv — the
     deterministic momentum loss (Δp_z) = the classical "momentum/scattering".
  D  WP momentum distribution n_wp(|k|) before vs after (WP transit run) — the
     quantum momentum spread/redistribution the classical projectile lacks.

Known-case (printed): classical loss-fn m=1 peak ω; classical Δp_z<0 (decel).
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
CL = JB / "run_classical_n162_L50_E100_v2"
E15 = JB / "run_plasmon_n162_L50_E15"
WP = JB / "run_wp_n162_L50_E100_sigma1_v2"
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


qc, omc, Lc = loss_qw(CL / "results/analysis/observables/n_q_vs_time.csv")
qe, ome, Le = loss_qw(E15 / "results/analysis/observables/n_q_vs_time.csv")
Tc = pd.read_csv(CL / "results/analysis/observables/n_q_vs_time.csv")["time_au"].max()
ip = int(np.argmax(Lc[0, 1:]) + 1)
print(f"[KC] classical T={Tc:.1f} a.u. (Δω={2*np.pi/Tc*HA:.1f} eV); "
      f"classical m=1 L-peak ω={omc[ip]*HA:.2f} eV (ω_p={omega_p*HA:.2f})")

fig, axs = plt.subplots(2, 2, figsize=(13, 9))

# A: classical L(q,w) map
ax = axs[0, 0]; sel = (omc * HA <= 12) & (omc > 0)
Q, W = np.meshgrid(qc, omc[sel] * HA, indexing="ij"); Lp = Lc[:, sel]
vmax = np.percentile(Lp[Lp > 0], 99.5)
pc = ax.pcolormesh(Q, W, Lp, shading="auto", cmap="inferno", norm=LogNorm(vmax * 1e-4, vmax))
qq = np.linspace(qc.min(), qc.max(), 80)
ax.plot(qq, bohm_gross(qq) * HA, "c-", lw=1.5, label="Bohm-Gross")
ax.axhline(omega_p * HA, color="w", ls=":", lw=1.0, label=f"ω_p={omega_p*HA:.2f} eV")
ax.set_xlabel("q (Bohr⁻¹)"); ax.set_ylabel("ω (eV)"); ax.set_ylim(0, 12)
ax.set_title(f"A — classical L(q,ω)  [CAVEAT Δω≈{2*np.pi/Tc*HA:.1f} eV: coarse]")
ax.legend(fontsize=7); fig.colorbar(pc, ax=ax, shrink=0.8)

# B: classical vs E15 medium m=1 (area-normalised)
ax = axs[0, 1]
for q, om, L, lab, c in [(qc, omc, Lc, f"classical (T={Tc:.0f}a.u.)", "C3"),
                         (qe, ome, Le, "E15 medium (T=2000a.u.)", "C0")]:
    s = (om * HA <= 12) & (om > 0); wL = om[s] * L[0, s]
    ax.plot(om[s] * HA, wL / max(np.trapezoid(wL, om[s] * HA), 1e-30), color=c, lw=1.6, label=lab)
ax.axvline(omega_p * HA, color="k", ls=":", lw=1.0, label="ω_p")
ax.axvline(bohm_gross(q1) * HA, color="c", ls="-", lw=1.0, label="Bohm-Gross(q₁)")
ax.set_xlabel("ω (eV)"); ax.set_ylabel("ω·L(ω) m=1 (area-norm)"); ax.set_xlim(0, 12)
ax.set_title("B — does the classical projectile excite the plasmon at ω_p?")
ax.legend(fontsize=7); ax.grid(alpha=0.3)

# C: classical projectile v_z(t)
ax = axs[1, 0]
tr = pd.read_csv(CL / "results/raw/observables/electron_track.csv")
ax.plot(tr["time_au"], tr["vz"], "C3-", lw=1.8)
dpz = tr["vz"].values[-1] - tr["vz"].values[0]
ax.set_xlabel("t (a.u.)"); ax.set_ylabel("projectile v_z (a.u.)")
ax.set_title(f"C — classical trajectory: deterministic momentum loss\nΔp_z={dpz:+.3f} a.u. (single trajectory, no distribution)")
ax.grid(alpha=0.3)
print(f"[KC] classical Δp_z={dpz:+.3f} a.u. (expect <0 decel)")

# D: WP momentum distribution before/after (quantum spread)
ax = axs[1, 1]
md = pd.read_csv(WP / "results/raw/observables/momentum_distribution.csv", comment="#")
t0 = md["time_au"].min(); tN = md["time_au"].max()
b = md[np.isclose(md["time_au"], t0)]; a = md[np.isclose(md["time_au"], tN)]
ax.plot(b["k_bohr_inv"], b["n_wp"], "C0-", lw=1.6, label=f"WP n(|k|) t={t0:.1f}")
ax.plot(a["k_bohr_inv"], a["n_wp"], "C2-", lw=1.6, label=f"WP n(|k|) t={tN:.1f}")
k0 = np.sqrt(2 * 100 / HA)
ax.axvline(k0, color="0.5", ls="--", lw=1.0, label=f"k₀(100eV)={k0:.2f}")
ax.set_xlabel("|k| (Bohr⁻¹)"); ax.set_ylabel("WP momentum density n_wp(|k|)")
ax.set_title("D — WP momentum distribution broadens/redistributes (quantum)\nvs classical δ-trajectory in C")
ax.legend(fontsize=7); ax.grid(alpha=0.3); ax.set_xlim(0, min(6, a["k_bohr_inv"].max()))

fig.suptitle("Classical-electron confidence (A,B) + WP-vs-classical momentum contrast (C,D), E=100 eV",
             fontsize=12)
fig.tight_layout()
fp = OUT / "classical_confidence_analysis.png"
fig.savefig(fp, dpi=150); plt.close(fig)
print(f"wrote {fp}")
