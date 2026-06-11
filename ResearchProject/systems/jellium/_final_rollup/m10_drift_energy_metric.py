#!/usr/bin/env python3
"""M10 — drift-energy metric vs total-KE metric for the WP σ=1 stopping.

M10 question: "classical is not an upper bound in the master plot → metric or physics?"
A WP loses z-kinetic energy two ways:
  • DRIFT deceleration   E_drift(t)  = ⟨p_z⟩²/2          (genuine slowing / stopping)
  • momentum BROADENING  σ_pz²(t)/2  (spread; a FREE WP keeps σ_p CONSTANT, so any
                                      growth is interaction/metric, not real stopping)
z-total KE = ⟨p_z²⟩/2 = (⟨p_z⟩² + σ_pz²)/2 mixes both. This script separates them:
  S_drift = −Δ(⟨p_z⟩²/2)/Δz      S_ztot = −Δ(⟨p_z²⟩/2)/Δz
over the interference-free window, for E=20..300 σ=1, alongside classical S.
If S_ztot > S_drift, part of the "extra" WP stopping is momentum spreading, not drag.

Output (this dir): m10_drift_energy_metric.png  + printed table.
"""
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
from pathlib import Path
from applications.report1 import stopping_power_data as spd

HA = 27.211386245988
ROOT = Path("/local/data/public/skcb2/tddft")


def metrics_for(run_dir: Path):
    params = spd.parse_run_summary(run_dir)
    win = spd.compute_time_window(params)
    csv = run_dir / "results/raw/observables/wp_momentum_stats.csv"
    if not csv.exists():
        return None
    df = pd.read_csv(csv, comment="#")
    df = df[df["time_au"] <= win.t_end]
    if len(df) < 2:
        return None
    pz = df["pz_mean"].values
    pz2 = df["pz2_mean"].values            # ⟨p_z²⟩
    E_drift = pz**2 / 2.0                   # ⟨p_z⟩²/2  (Ha)
    E_ztot = pz2 / 2.0                      # ⟨p_z²⟩/2   (Ha)
    dz = params.wp_k0_z * win.t_end
    if abs(dz) < 0.1:
        return None
    S_drift = -(E_drift[-1] - E_drift[0]) * HA / dz
    S_ztot = -(E_ztot[-1] - E_ztot[0]) * HA / dz
    sig0 = (pz2[0] - pz[0]**2); sigE = (pz2[-1] - pz[-1]**2)
    return dict(S_drift=S_drift, S_ztot=S_ztot, dz=dz, t_end=win.t_end,
                sig_pz2_0=sig0, sig_pz2_end=sigE, v0=abs(params.wp_k0_z))


# WP σ=1 v2 runs (E=20..300)
wp_specs = spd.get_L50_wp_sigma1_runs()
rows = []
print(f"{'E(eV)':>6} {'v0':>6} {'S_drift':>9} {'S_ztot':>9} {'σpz²0':>9} {'σpz²end':>9} {'Δσpz²%':>8}")
for sp in wp_specs:
    rd = ROOT / sp.run_dir if not str(sp.run_dir).startswith("/") else Path(sp.run_dir)
    m = metrics_for(rd)
    if m is None:
        print(f"{sp.energy_eV:>6} skip (no momentum csv / window)"); continue
    dpct = 100 * (m["sig_pz2_end"] - m["sig_pz2_0"]) / max(m["sig_pz2_0"], 1e-12)
    rows.append((sp.energy_eV, m))
    print(f"{sp.energy_eV:>6} {m['v0']:>6.2f} {m['S_drift']:>9.4f} {m['S_ztot']:>9.4f} "
          f"{m['sig_pz2_0']:>9.3e} {m['sig_pz2_end']:>9.3e} {dpct:>8.1f}")

# classical S at each energy
cls = {}
for sp in spd.get_L50_classical_runs():
    rd = ROOT / sp.run_dir
    try:
        params = spd.parse_run_summary(rd); win = spd.compute_time_window(params)
        S, _ = spd.compute_classical_S(rd, win)
        if np.isfinite(S):
            cls[sp.energy_eV] = abs(S)
    except Exception:
        pass

E = np.array([r[0] for r in rows])
Sd = np.array([r[1]["S_drift"] for r in rows])
Sz = np.array([r[1]["S_ztot"] for r in rows])
v0 = np.array([r[1]["v0"] for r in rows])

fig, ax = plt.subplots(figsize=(7.5, 5))
ax.plot(v0, Sd, "o-", color="#1f77b4", lw=1.8, ms=6, label=r"WP $S_{\rm drift}=-\Delta(\langle p_z\rangle^2/2)/\Delta z$")
ax.plot(v0, Sz, "s--", color="#d62728", lw=1.8, ms=6, label=r"WP $S_{z\rm tot}=-\Delta(\langle p_z^2\rangle/2)/\Delta z$")
if cls:
    ce = sorted(cls); cv = [np.sqrt(2 * e / HA) for e in ce]
    ax.plot(cv, [cls[e] for e in ce], "k^-", lw=1.5, ms=7, label="classical Ehrenfest")
ax.axhline(0, color="0.7", lw=0.6)
ax.set_xlabel(r"$v_0$ (a.u.)"); ax.set_ylabel("S (eV/Bohr)")
ax.set_title("M10 — WP drift vs total-KE stopping metric (σ=1) vs classical")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout()
fp = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/_final_rollup/m10_drift_energy_metric.png")
fig.savefig(fp, dpi=150); plt.close(fig)
print(f"\nwrote {fp}")
