#!/usr/bin/env python3
"""Extensive-kinetic (norm-division fix) validation — vacuum, double-sided CAP.

Compares, from the dcap_extkin run (OrbitalKineticStats ON):
  * INQ out-of-the-box:  energies.csv:kinetic  (norm-divided, energy.hpp:55)
  * our schema:          orbital_kinetic_stats.csv:kin_bare_total_ha (extensive)
  * identity check:      kin_normdiv_total_ha == energies.csv:kinetic per step
  * post-hoc route:      e_kin_ha * norm (wp_momentum_stats x real-space norm)
and the timing overhead vs the dcap_baseline run (OrbitalKineticStats OFF).

Outputs (this folder): fig_extkin_energies.png, fig_extkin_identity_timing.png,
extkin_summary.txt
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, re

ROOT = Path("/local/data/public/skcb2/tddft")
sys.path.insert(0, str(ROOT / "inq-stack" / "python"))
from inqview.visualisation.style import apply_theme  # noqa: E402
apply_theme()

HA = 27.211386
RES = ROOT / "ResearchProject/systems/vacuum/scripts/wp_traversal_energy/results"
OUT = Path(__file__).resolve().parent


def load(run):
    d = RES / run / "raw/observables"
    en = pd.read_csv(d / "energies.csv")
    ek = None
    p = d / "orbital_kinetic_stats.csv"
    if p.exists():
        ek = pd.read_csv(p, comment="#")
    mom = pd.read_csv(d / "wp_momentum_stats.csv", comment="#")
    return en, ek, mom


def summary_field(run, key):
    txt = (RES / run / "run_summary.txt").read_text()
    m = re.search(rf"{key}\s*=\s*([0-9.eE+-]+)", txt)
    return float(m.group(1)) if m else np.nan


en, ek, mom = load("dcap_extkin")
m = en.merge(ek, left_on="step", right_on="step", suffixes=("", "_ek"))

E_rep = m["total"].to_numpy(float)              # reported total (Ha)
K_rep = m["kinetic"].to_numpy(float)            # reported kinetic (norm-divided)
K_bare = m["kin_bare_total_ha"].to_numpy(float) # extensive kinetic (ours)
K_recon = m["kin_normdiv_total_ha"].to_numpy(float)
norm = m["norm_total"].to_numpy(float)
t = m["time_au"].to_numpy(float)
E_corr = E_rep - K_rep + K_bare                 # corrected (extensive) total

# post-hoc route on the same run: e_kin_ha * norm
tm = mom["time_au"].to_numpy(float)
ekin_mean = np.interp(t, tm, mom["e_kin_ha"].to_numpy(float))
K_posthoc = ekin_mean * norm                    # occ=1, single orbital

ident = K_recon - K_rep                         # must be ~solver precision
E0 = E_rep[0]

# ---- figure 1: the energy comparison -----------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(t, E_rep * HA, label="E_total reported (INQ)", lw=2)
ax[0].plot(t, E_corr * HA, label="E_total corrected (bare kinetic)", lw=2)
ax[0].plot(t, E0 * norm * HA, "--", label="E0 · norm (expected extensive)", lw=1.4)
axn = ax[0].twinx()
axn.plot(t, norm, color="0.5", lw=1, alpha=0.7)
axn.set_ylabel("WP norm", color="0.5")
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("E (eV)")
ax[0].set_title("Total energy: reported vs corrected")
ax[0].legend(loc="center left", fontsize=8)

ax[1].plot(t, K_rep * HA, label="kinetic reported (norm-divided)", lw=2)
ax[1].plot(t, K_bare * HA, label="kinetic bare (extensive, ours)", lw=2)
ax[1].plot(t, K_posthoc * HA, "--", label="e_kin_ha · norm (post-hoc route)", lw=1.4)
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("E_kin (eV)")
ax[1].set_title("Kinetic channel")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "fig_extkin_energies.png", dpi=160)

# ---- figure 2: identity residual + timing -------------------------------
wall = m["wall_ms"].to_numpy(float)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].semilogy(t, np.abs(ident) * HA, lw=1)
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("|Σocc·T/N − kinetic_INQ| (eV)")
ax[0].set_title("Per-step identity residual (ours vs INQ)")
ax[1].plot(t, wall, lw=0.8)
ax[1].axhline(np.mean(wall[1:]), color="C1", ls="--",
              label=f"mean {np.mean(wall[1:]):.2f} ms")
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("observable wall (ms/step)")
ax[1].set_title("OrbitalKineticStats cost")
ax[1].legend()
fig.tight_layout()
fig.savefig(OUT / "fig_extkin_identity_timing.png", dpi=160)

# ---- summary -------------------------------------------------------------
per_step_on = summary_field("dcap_extkin", "per_step_ms")
per_step_off = summary_field("dcap_baseline", "per_step_ms")
overhead_ms = per_step_on - per_step_off
lines = [
    "extensive-kinetic validation — vacuum double-sided CAP (dcap_extkin vs dcap_baseline)",
    f"E0 = {E0*HA:.1f} eV   final norm = {norm[-1]:.2e}",
    f"E_reported(final) = {E_rep[-1]*HA:.1f} eV   (artifact: stays at the mean)",
    f"E_corrected(final) = {E_corr[-1]*HA:.2f} eV  (E0·norm = {E0*norm[-1]*HA:.2f} eV)",
    f"captured = {(E0 - E_corr[-1])/E0*100:.1f}% of E0",
    f"identity max |Σocc·T/N − kinetic_INQ| = {np.max(np.abs(ident)):.2e} Ha",
    f"post-hoc vs in-run bare kinetic: max |Δ| = {np.max(np.abs(K_posthoc-K_bare))*HA:.2e} eV",
    f"timing: per-step {per_step_on:.1f} ms (obs ON) vs {per_step_off:.1f} ms (OFF)"
    f" -> overhead {overhead_ms:.1f} ms/step = {overhead_ms/per_step_off*100:.1f}%",
    f"observable self-timing: mean {np.mean(wall[1:]):.2f} ms/step (1 orbital)",
]
(OUT / "extkin_summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
