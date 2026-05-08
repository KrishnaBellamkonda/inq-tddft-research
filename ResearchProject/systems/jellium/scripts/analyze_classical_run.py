#!/usr/bin/env python3
"""Comprehensive post-run analysis for the classical-electron jellium run.

Runs the full inqview postprocess pipeline on a results/ tree, then
computes physics-summary numbers (stopping power, energy budget,
GS-projected excitation amplitude) and writes a markdown report.

Usage:
    python3 analyze_classical_run.py <run_dir>

Example:
    python3 analyze_classical_run.py \\
      /local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_classical_e1500_L50_cubic
"""
from __future__ import annotations

import argparse
import re
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

HA_TO_EV = 27.211386245988
M_E_AU   = 1.0    # electron mass in atomic units


def read_summary_value(path: Path, key: str, default=None):
    if not path.exists():
        return default
    m = re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)",
                  path.read_text(), flags=re.MULTILINE)
    if not m:
        return default
    s = m.group(1)
    try:
        return float(s)
    except ValueError:
        return s


def physics_block(run_dir: Path) -> dict:
    """Extract physics summary numbers from CSVs."""
    results = run_dir / "results"
    obs_dir = results / "raw" / "observables"

    obs   = pd.read_csv(obs_dir / "observables.csv")
    track = pd.read_csv(obs_dir / "electron_track.csv")

    # Trajectory
    z0 = float(track["z"].iloc[0]); z1 = float(track["z"].iloc[-1])
    distance = z1 - z0
    v0 = float(track["vz"].iloc[0]); v1 = float(track["vz"].iloc[-1])
    KE0 = 0.5 * M_E_AU * v0 ** 2
    KE1 = 0.5 * M_E_AU * v1 ** 2
    KE_loss_eV = (KE0 - KE1) * HA_TO_EV

    # Energies
    E0 = float(obs["energy_total"].iloc[0])
    E1 = float(obs["energy_total"].iloc[-1])
    Ek0 = float(obs["energy_kinetic"].iloc[0])
    Ek1 = float(obs["energy_kinetic"].iloc[-1])
    Eh0 = float(obs["energy_hartree"].iloc[0])
    Eh1 = float(obs["energy_hartree"].iloc[-1])
    Exc0 = float(obs["energy_xc"].iloc[0])
    Exc1 = float(obs["energy_xc"].iloc[-1])

    # Energy conservation
    proj_loss = KE_loss_eV
    bath_gain = (E1 - E0) * HA_TO_EV
    cons_residual = bath_gain + (-proj_loss)  # signed: bath gain minus proj loss

    # Stopping power
    S = (KE0 - KE1) / distance * HA_TO_EV  # eV/Bohr

    # Eigenvalue summary
    ev_path = obs_dir / "eigenvalues" / "eigenvalues.csv"
    occ_path = obs_dir / "eigenvalues" / "occupations.csv"
    homo_idx = lumo_idx = gap_eV = None
    if ev_path.exists() and occ_path.exists():
        ev = pd.read_csv(ev_path).set_index("state_index").sort_index()
        occ = pd.read_csv(occ_path).set_index("state_index").sort_index()
        joined = ev.join(occ, lsuffix="_e", rsuffix="_o")
        homo_idx = int(joined[joined["occupation"] >= 0.5].index.max())
        if homo_idx + 1 in joined.index:
            lumo_idx = homo_idx + 1
            gap_eV = float(joined.loc[lumo_idx, "eigenvalue_ha"]
                           - joined.loc[homo_idx, "eigenvalue_ha"]) * HA_TO_EV

    # GS-projected excitation (uses gs_projected_occupations output)
    gs_proj_csv = (results / "analysis" / "observables"
                   / "gs_projected_occupations" / "excitation_total_vs_time.csv")
    excitation_loss = excitation_t = None
    if gs_proj_csv.exists():
        df = pd.read_csv(gs_proj_csv)
        last = df.iloc[-1]
        excitation_loss = float(last["loss_from_occupied"])
        excitation_t = float(last["time_au"])

    return {
        "z0": z0, "z1": z1, "distance": distance,
        "v0": v0, "v1": v1,
        "KE0_eV": KE0 * HA_TO_EV, "KE1_eV": KE1 * HA_TO_EV,
        "KE_loss_eV": KE_loss_eV,
        "S_eV_per_Bohr": S,
        "E0_Ha": E0, "E1_Ha": E1,
        "dE_total_eV": (E1 - E0) * HA_TO_EV,
        "dE_kinetic_eV": (Ek1 - Ek0) * HA_TO_EV,
        "dE_hartree_eV": (Eh1 - Eh0) * HA_TO_EV,
        "dE_xc_eV":      (Exc1 - Exc0) * HA_TO_EV,
        "energy_conservation_residual_eV": cons_residual,
        "homo_idx": homo_idx,
        "lumo_idx": lumo_idx,
        "gap_eV": gap_eV,
        "excitation_loss": excitation_loss,
        "excitation_t":    excitation_t,
        "n_steps": len(obs),
        "total_time_atu": float(obs["time_au"].iloc[-1]),
    }


def render_markdown(run_dir: Path, phys: dict, run_summary: dict) -> str:
    """Markdown report skeleton ready for journal-writing skill ingestion."""
    L_BOHR     = run_summary.get("cell_bohr", "?")
    N_e        = run_summary.get("n_electrons", "?")
    spacing    = run_summary.get("spacing_bohr", "?")
    KE_eV      = run_summary.get("projectile_KE_eV", "?")
    v_target   = run_summary.get("velocity_atu", "?")
    ion_dyn    = run_summary.get("ion_dynamics", "?")
    wall_s     = run_summary.get("wall_time_s", "?")

    md = f"""# Analysis: {run_dir.name}

## Run identity

| Field | Value |
|---|---|
| Run dir | `{run_dir}` |
| Cell | {L_BOHR} Bohr cubic, periodic |
| N_electrons | {N_e} |
| Spacing | {spacing} Bohr |
| Projectile KE | {KE_eV} eV |
| Initial velocity | {v_target} bohr/atu |
| Ion dynamics | {ion_dyn} |
| Total propagation time | {phys['total_time_atu']:.3f} atu |
| Wall time | {wall_s} s |

## Trajectory

| | Value |
|---|---|
| z(t=0) | {phys['z0']:.4f} Bohr (INQ centred) |
| z(t=end) | {phys['z1']:.4f} Bohr |
| Distance traveled | {phys['distance']:.4f} Bohr |
| v_z(t=0) | {phys['v0']:.6f} bohr/atu |
| v_z(t=end) | {phys['v1']:.6f} bohr/atu |
| Velocity drop | {phys['v0'] - phys['v1']:+.6f} bohr/atu ({100*(phys['v0']-phys['v1'])/phys['v0']:+.4f} %) |

## Stopping power

| | Value |
|---|---|
| Initial KE | {phys['KE0_eV']:.4f} eV |
| Final KE | {phys['KE1_eV']:.4f} eV |
| **Projectile KE loss** | **{phys['KE_loss_eV']:.4f} eV** |
| **Stopping power S(v)** | **{phys['S_eV_per_Bohr']:.6f} eV/Bohr** |

S(v) = -dE/dz = projectile KE loss / distance traveled.

## Energy budget (system)

| Component | dE (eV) |
|---|---|
| Total | {phys['dE_total_eV']:+.4f} |
| Kinetic | {phys['dE_kinetic_eV']:+.4f} |
| Hartree | {phys['dE_hartree_eV']:+.4f} |
| XC | {phys['dE_xc_eV']:+.4f} |

**Energy conservation residual**: ΔE_proj + ΔE_bath = {phys['energy_conservation_residual_eV']:+.6f} eV (should be 0).

## Excitation analysis (GS-basis projection)

The overlap matrix at t=end was used to compute

    n_i^GS(t)  =  sum_j  f_j(0)  *  |<psi_i^GS | psi_j(t)>|^2

The total occupation transferred *out of* the initially-occupied
subspace (= excitation amount) is recorded below.

| | Value |
|---|---|
| Snapshot time | {phys['excitation_t']:.3f} atu |
| **Loss from occupied subspace** | **{phys['excitation_loss']:.6f} electrons** |
| Implied # of e-h pairs | ~{(phys['excitation_loss'] or 0)/2:.4f} |

For comparison, ω_p = 3.47 eV (jellium plasmon at this density).
If average e-h pair gap ~ω_p, the equivalent energy in excitation is
~{(phys['excitation_loss'] or 0)*3.47/2:.3f} eV — to be compared against
the bath kinetic gain of {phys['dE_kinetic_eV']:.3f} eV.

## Spectrum

GS HOMO/LUMO: index {phys['homo_idx']} / {phys['lumo_idx']},
HOMO-LUMO gap = {phys['gap_eV']:.4f} eV.

## Files of interest

```
results/analysis/observables/
├── observables_summary.png
├── total_energy_vs_time.png
├── current_components_vs_time.png
├── dipole_components_vs_time.png
├── dE_kinetic_vs_z.png            <-- stopping curve
├── stopping_force_vs_z.png        <-- F_z(z)
├── bath_energy_vs_time.png        <-- band-structure-summed bath E(t)
├── ks_energies_absolute.gif       <-- per-orbital E_i(t)
├── ks_energies_delta.gif          <-- delta from t=0
├── occupations_absolute.gif       <-- per-orbital f_i(t)
├── occupations_delta.gif
├── fft_total_energy.png
├── fft_current_*.png
├── dipole_spectrum_*.png
└── gs_projected_occupations/
    ├── n_i_gs_vs_time.csv         <-- raw GS-projected occupations
    ├── excitation_total_vs_time.{{csv,png}}
    ├── gs_projection_t000000.png  <-- bar chart at t=0 (= identity)
    ├── gs_projection_t<NN>.png    <-- bar charts at later snapshots
    └── heatmap_overlap_t<NN>.png  <-- |O_ij(t)|^2 heatmaps (full-matrix)
```

## Source data

```
results/raw/observables/
├── observables.csv             ({phys['n_steps']} rows)
├── electron_track.csv          (every step: pos, vel, F=0 placeholder)
├── state_energies.csv          (per-orbital E_i(t))
├── occupations_vs_time.csv     (per-orbital f_i(t))
├── momentum_distribution.csv   (n(|k|, t))
├── overlap_full/               (full O matrix at t=0, mid, end)
├── overlap_proxies/            (proxy O matrix at many snapshots, future runs)
└── eigenvalues/                (GS eigenvalues + occupations)
```
"""
    return textwrap.dedent(md)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("run_dir", type=Path,
                   help="Path to run directory containing results/")
    p.add_argument("--report-out", type=Path, default=None,
                   help="Where to write the markdown report "
                        "(default: <run_dir>/results/analysis/REPORT.md)")
    args = p.parse_args()

    run_dir = args.run_dir.resolve()
    if not (run_dir / "results").exists():
        print(f"FATAL: {run_dir}/results does not exist", file=sys.stderr)
        sys.exit(2)

    # Run the postprocess pipeline.
    print(f"=== Running postprocess pipeline on {run_dir.name} ===")
    from inqview.postprocess import pipeline
    pipeline_phases = ["summary", "observables", "state_energies",
                       "bath_energy", "stopping",
                       "gs_projected_occupations", "occupations"]
    res = pipeline.run(
        run_dir / "results",
        run_name=run_dir.name,
        phases=pipeline_phases,
        rebuild=False,
        skip_paraview=True,
    )

    # Read run_summary.txt for metadata.
    rs_path = run_dir / "results" / "run_summary.txt"
    rs_keys = ["cell_bohr", "n_electrons", "spacing_bohr",
               "projectile_KE_eV", "velocity_atu", "ion_dynamics",
               "wall_time_s"]
    rs = {k: read_summary_value(rs_path, k, "?") for k in rs_keys}

    # Compute physics summary.
    phys = physics_block(run_dir)

    # Render and write the markdown report.
    md = render_markdown(run_dir, phys, rs)
    out_md = args.report_out or (run_dir / "results" / "analysis" / "REPORT.md")
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(md)
    print(f"\n=== Wrote analysis report ===\n{out_md}\n")

    # Print compact physics summary to stdout.
    print("=== Physics summary ===")
    print(f"Distance traveled:   {phys['distance']:+.4f} Bohr")
    print(f"v_z drop:            {phys['v0']-phys['v1']:+.6f} bohr/atu "
          f"({100*(phys['v0']-phys['v1'])/phys['v0']:+.4f} %)")
    print(f"KE loss:             {phys['KE_loss_eV']:.4f} eV")
    print(f"Stopping power S(v): {phys['S_eV_per_Bohr']:.6f} eV/Bohr")
    print(f"Energy conservation: {phys['energy_conservation_residual_eV']:+.6f} eV (should be 0)")
    if phys['excitation_loss'] is not None:
        print(f"GS-basis excitation: {phys['excitation_loss']:.6f} electrons (= ~{phys['excitation_loss']/2:.4f} e-h pairs)")
    print(f"Pipeline phases:")
    for ph in pipeline_phases:
        if ph in res.ok:    print(f"  [ok]   {ph}")
        if ph in res.skipped:print(f"  [skip] {ph}: {res.skipped[ph]}")
        if ph in res.failed: print(f"  [fail] {ph}: {res.failed[ph]}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
