#!/usr/bin/env python3
"""Per-run loss-function stopping power: integrate EACH run's L(q,ω) -> S(v).

For each long, frequency-resolved run (E15, E3.4, and E25 when ready) extract
L(q,ω)=|FFT(n_q·hann)|²/q² and integrate via the dielectric stopping formula
    S(v) = (2/πv²) ∫(dq/q) ∫_0^{min(qv,cap)} ω L(q,ω) dω      [Lindhard/Ritchie]
to get a stopping curve FROM THAT RUN'S OWN loss function.

Consistency test: L(q,ω) is a MEDIUM property, so if the extraction is clean &
linear, every run's S(v) SHAPE should coincide (independent of the projectile
that excited it). Curves are normalised to unit area over v so shapes compare
directly (absolute scale depends on excitation amplitude + needs FDT calib).
Analytic Lindhard (full q) shown for reference. Each run's own velocity marked.

Output: batch2_figures/loss_function_stopping_per_run.png
Known-case (printed): each run's m=1 L-peak near ω_p; S(v) peak velocity.
"""
from __future__ import annotations
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
from inqview.postprocess import lindhard

HA = 27.211386245988
JB = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium")
OUT = Path("/local/data/public/skcb2/tddft/docs/presentations/storyline/tasks/batch2_figures")
L_BOHR, N = 50.0, 162
N_DENS = N / L_BOHR**3
kF = (3 * np.pi**2 * N_DENS) ** (1 / 3); omega_p = np.sqrt(4 * np.pi * N_DENS); vF = kF
OMEGA_CAP = 16.0 / HA


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


def S_from_L(v, q, om, Lmat):
    g = np.zeros_like(q)
    for i, qi in enumerate(q):
        wmax = min(qi * v, OMEGA_CAP); sel = (om > 0) & (om <= wmax)
        if sel.sum() < 2: continue
        g[i] = np.trapezoid(om[sel] * Lmat[i, sel], om[sel]) / qi
    return (2.0 / (np.pi * v**2)) * np.trapezoid(g, q)


RUNS = [
    ("E15 (15 eV, v=1.05)", 1.050, "#d62728", JB / "run_plasmon_n162_L50_E15"),
    ("E3.4 (3.4 eV, v=0.50)", 0.500, "#1f77b4", JB / "run_plasmon_n162_L50_E3p4_varyv"),
]
e25 = JB / "run_plasmon_n162_L50_E25"
if (e25 / "results/analysis/observables/n_q_vs_time.csv").exists():
    RUNS.append(("E25 (25 eV, v=1.36)", 1.356, "#2ca02c", e25))

v_grid = np.linspace(0.2, 5.0, 80)
S_full = np.array([lindhard.stopping_power(v, kF) for v in v_grid]) * HA
print(f"system r_s={(3/(4*np.pi*N_DENS))**(1/3):.2f}  ω_p={omega_p*HA:.2f} eV")

fig, (axS, axN) = plt.subplots(1, 2, figsize=(13.5, 5))
print("\n=== per-run loss-function S(v) ===")
curves = []
for label, v0, c, run in RUNS:
    csv = run / "results/analysis/observables/n_q_vs_time.csv"
    if not csv.exists():
        print(f"  {label}: no n_q csv, skip"); continue
    q, om, L = loss_qw(csv)
    ipk = int(np.argmax(L[0, 1:]) + 1)
    S = np.array([S_from_L(v, q, om, L) for v in v_grid])
    curves.append((label, v0, c, S))
    iv = int(np.argmin(np.abs(v_grid - v0)))
    print(f"  {label}: m=1 peak ω={om[ipk]*HA:.2f} eV | S(v0) (arb)={S[iv]:.3e} | "
          f"S-peak at v={v_grid[int(np.argmax(S))]:.2f}")

# Panel A: absolute (arb units) — shows excitation-amplitude differences
for label, v0, c, S in curves:
    axS.plot(v_grid, S, color=c, lw=2, label=label)
    axS.axvline(v0, color=c, ls=":", lw=1.0)
axS.set_xlabel("v (a.u.)"); axS.set_ylabel("loss-function S(v) (arb units)")
axS.set_title("Per-run loss-function S(v) — absolute (arb)\n(dotted = each run's own velocity)")
axS.legend(fontsize=8); axS.grid(alpha=0.3); axS.set_xlim(0, 5); axS.set_ylim(bottom=0)

# Panel B: unit-area normalised — SHAPE consistency test + analytic Lindhard
for label, v0, c, S in curves:
    area = np.trapezoid(S, v_grid)
    axN.plot(v_grid, S / max(area, 1e-30), color=c, lw=2, label=f"{label} (norm)")
axN.plot(v_grid, S_full / np.trapezoid(S_full, v_grid), "k--", lw=1.5, label="analytic Lindhard (norm)")
axN.set_xlabel("v (a.u.)"); axN.set_ylabel("S(v) (unit-area normalised)")
axN.set_title("SHAPE consistency: do the runs' S(v) coincide?\n(L is a medium property → should overlap)")
axN.legend(fontsize=8); axN.grid(alpha=0.3); axN.set_xlim(0, 5); axN.set_ylim(bottom=0)
fig.suptitle("Loss-function stopping power integrated PER RUN (E15, E3.4"
             + (", E25" if len(curves) > 2 else "") + ")", fontsize=12)
fig.tight_layout()
fp = OUT / "loss_function_stopping_per_run.png"
fig.savefig(fp, dpi=150); plt.close(fig)
print(f"\nwrote {fp}")
