#!/usr/bin/env python3
"""Phase 0 — mine existing localised-jellium runs' energy series (FREE, no GPU).

Characterises the ΔE_total anomaly from data already on disk: which runs show the
unphysical >0 rise, which decay cleanly, and the component gap (only
total/kin/hartree/xc recorded). Emits result.json + a phase0 plot for the Advisor.

Reference: campaign docs/campaigns/localised_jellium/energy-oscillation-diagnosis.md
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HA_EV = 27.211386
ROOT = Path("/local/data/public/skcb2/tddft")
HYP = ROOT / "ResearchProject/systems/localised_jellium/hypotheses/energy_oscillation_diagnosis"
PROBE_DIR = HYP / "probes" / "phase0_mine"
PROBE_DIR.mkdir(parents=True, exist_ok=True)

E_GS = -36.940459047122  # effmass_sigma1 GS SCF total (gs_build_run.log), Ha

RUNS = {
    "effmass_sigma1 (eta=-1.0, default)":
        "scripts/muon_mass_fork/effmass_sigma1/wp/results/sigma1/raw/observables/observables.csv",
    "cap_eta0p4 (weak eta=-0.4)":
        "scripts/muon_mass_fork/effmass_sigma1/wp/results/cap_eta0p4/raw/observables/observables.csv",
    "cap_eta2p0 (strong eta=-2.0)":
        "scripts/muon_mass_fork/effmass_sigma1/wp/results/cap_eta2p0/raw/observables/observables.csv",
    "cap_gap19p5 (wider gap)":
        "scripts/muon_mass_fork/effmass_sigma1/wp/results/cap_gap19p5/raw/observables/observables.csv",
    "classical twin":
        "scripts/muon_mass_fork/effmass_sigma1/classical/results/classical/raw/observables/observables.csv",
    "p3_wp (clean contrast)":
        "scripts/fullsuite_wp/results/p3_wp/raw/observables/observables.csv",
}
LJ = ROOT / "ResearchProject/systems/localised_jellium"

records = []
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for label, rel in RUNS.items():
    csv = LJ / rel
    df = pd.read_csv(csv)
    t = df["time_au"].to_numpy(float)
    et = df["energy_total"].to_numpy(float)
    d_rt0 = (et - et[0]) * HA_EV
    d_gs = (et - E_GS) * HA_EV  # note: E_GS is 52-electron slab; RT is 53 (WP added) except classical
    rec = {
        "run": label, "csv": str(csv), "columns": list(df.columns),
        "n_rows": int(len(df)), "t_final_au": float(t[-1]),
        "E_total0_ha": float(et[0]), "E_total_final_ha": float(et[-1]),
        "dE_vs_rt0_final_ev": float(d_rt0[-1]), "dE_vs_rt0_max_ev": float(np.max(d_rt0)),
        "dE_vs_rt0_min_ev": float(np.min(d_rt0)),
        "crosses_zero_above": bool(np.max(d_rt0) > 1e-2),
        "dE_vs_gs_final_ev": float(d_gs[-1]),
    }
    records.append(rec)
    axes[0].plot(t, d_rt0, marker=".", ms=3, label=label)
    # oscillation view: same, zoomed
    axes[1].plot(t, d_rt0, marker=".", ms=3, label=label)

axes[0].axhline(0, color="0.5", lw=1)
axes[0].set_xlabel("t (a.u.)"); axes[0].set_ylabel(r"$E_{tot}(t)-E_{tot}(0_{RT})$ (eV)")
axes[0].set_title("All runs (vs RT t=0)"); axes[0].legend(fontsize=7)
axes[1].axhline(0, color="0.5", lw=1)
axes[1].set_ylim(-8, 45)
axes[1].set_xlabel("t (a.u.)"); axes[1].set_ylabel(r"$\Delta E_{tot}$ (eV)")
axes[1].set_title("Zoom: positive-rise phenomenon"); axes[1].legend(fontsize=7)
fig.suptitle("Phase 0 — mined localised-jellium energy series (component gap: only total/kin/H/xc)")
fig.tight_layout()
plot_png = PROBE_DIR / "phase0_mine_energy.png"
fig.savefig(plot_png, dpi=130); plt.close(fig)

component_gap = ["energy_external", "energy_nonlocal", "energy_ion",
                 "energy_ion_kinetic", "energy_eigenvalues", "energy_nvxc"]
result = {
    "name": "phase0_mine",
    "aim": "Characterise the DeltaE_total anomaly from existing runs; confirm the component gap.",
    "method": "Read observables.csv of 6 localised-jellium runs; compute DeltaE_total vs RT t=0 and vs E_GS; plot.",
    "e_gs_ha": E_GS,
    "component_gap_all_runs": component_gap,
    "plot_png": str(plot_png),
    "runs": records,
}
(PROBE_DIR / "result.json").write_text(json.dumps(result, indent=2))
print(json.dumps({r["run"]: {"dE_vs_rt0_final_ev": r["dE_vs_rt0_final_ev"],
                             "dE_vs_rt0_max_ev": r["dE_vs_rt0_max_ev"],
                             "crosses_zero_above": r["crosses_zero_above"]}
                  for r in records}, indent=2))
print(f"\n[component gap] all mined runs record only: {records[0]['columns']}")
print(f"plot: {plot_png}")
