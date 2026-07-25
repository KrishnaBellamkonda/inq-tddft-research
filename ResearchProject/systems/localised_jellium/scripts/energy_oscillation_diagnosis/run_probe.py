#!/usr/bin/env python3
"""run_probe.py — mechanical dispatch of ONE tiny ablation probe for the
localised-jellium ΔE_total energy-oscillation diagnosis campaign.

The Investigator agent calls this per iteration. It is the deterministic half of
the loop: build + run a probe (or reuse an existing run's CSV), extract the full
energy decomposition, compute BOTH ΔE conventions, plot, and emit result.json for
the Advisor to read. The *reasoning* (which probe next) is the agent's job.

Contract (result.json), consumed by the Advisor and build_master_notebook.py:
{
  "name": str, "aim": str, "method": str,
  "run_dir": str, "observables_csv": str, "plot_png": str,
  "columns": [str...],            # which energy components were recorded
  "component_gap": [str...],      # requested-but-missing components
  "summary": {                    # raw numbers, faithfully reported
     "n_steps": int, "t_final_au": float,
     "dE_total_vs_gs_final": float, "dE_total_vs_gs_max": float,
     "dE_total_vs_rt0_final": float, "dE_total_vs_rt0_max": float,
     "crosses_zero_above": bool,  # does E_total(t) - E_total(0) ever go > 0
     "N_initial": float, "N_final": float,
     "eig_tracks_total": bool|null   # does Sum eps_i track E_total (hyp e)
  }
}

Usage:
  run_probe.py --name pure_gs_floor --aim "..." --run-dir <dir> [--env EM_CAP=0 ...]
  run_probe.py --name mine_effmass --mine <path/to/observables.csv> --e-gs <Ha>
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys
from pathlib import Path

CAMPAIGN = "lj-energy-oscillation-diagnosis"
HYP_ROOT = Path(__file__).resolve().parents[2] / "hypotheses" / "energy_oscillation_diagnosis"
PROBES = HYP_ROOT / "probes"
HA_EV = 27.211386

# Full energy decomposition we want recorded (see observables_writer.hpp).
WANT_COMPONENTS = [
    "energy_total", "energy_kinetic", "energy_hartree", "energy_xc",
    "energy_external", "energy_nonlocal", "energy_ion", "energy_ion_kinetic",
    "energy_eigenvalues", "energy_nvxc",
]


def _venv_python() -> str:
    return "/local/data/public/skcb2/tddft/venv/bin/python3"


def build_and_run(run_dir: Path, out_dir: Path, env: dict[str, str]) -> Path:
    """Build via inq-run (GPU) and run the probe, directing output to out_dir.

    CONTRACT for the agent: point the run.cpp output at out_dir (LJ_OUT env or the
    run's own OUT knob), pass the ablation env (e.g. EM_CAP=0), keep it tiny
    (~100-300 steps). Return the path to observables.csv. Reuse the effmass_sigma1
    GS checkpoint where possible; run cutoff_guard.py first. This skeleton shells
    out to inq-run; the agent wires the exact env per probe.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    runenv = {**os.environ, **env, "LJ_OUT": str(out_dir)}
    # NOTE: agent confirms the binary/build path; inq-run auto-detects the .cpp.
    subprocess.run(["inq-run"], cwd=run_dir, env=runenv, check=True)
    # observables path convention from the reference run.cpp:
    csv = out_dir / "raw" / "observables" / "observables.csv"
    if not csv.exists():
        # fall back to any observables*.csv the run wrote
        cand = sorted(out_dir.rglob("observables*.csv"))
        if not cand:
            raise FileNotFoundError(f"no observables csv under {out_dir}")
        csv = cand[0]
    return csv


def summarise(csv: Path, e_gs_ha: float | None, plot_png: Path) -> dict:
    import numpy as np, pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.read_csv(csv)
    cols = list(df.columns)
    gap = [c for c in WANT_COMPONENTS if c not in cols]
    et = df["energy_total"].to_numpy(dtype=float)
    t = df["time_au"].to_numpy(dtype=float) if "time_au" in cols else np.arange(len(et), dtype=float)

    d_rt0 = (et - et[0]) * HA_EV                       # vs E_total(0) of the RT run
    d_gs = (et - e_gs_ha) * HA_EV if e_gs_ha is not None else None  # vs E_GS

    eig_tracks = None
    if "energy_eigenvalues" in cols:
        eig = df["energy_eigenvalues"].to_numpy(dtype=float)
        # crude "tracks" test: correlation of increments over the window
        de, deig = np.diff(et), np.diff(eig)
        if de.std() > 0 and deig.std() > 0:
            eig_tracks = bool(np.corrcoef(de, deig)[0, 1] > 0.9)

    N = None
    for nc in ("electron_number", "n_electrons", "num_electrons"):
        if nc in cols:
            N = df[nc].to_numpy(dtype=float)
            break

    # plot: both conventions + any per-component drift
    fig, ax = plt.subplots(2, 1, figsize=(7.5, 6.5), sharex=True)
    ax[0].axhline(0, color="0.6", lw=0.8)
    ax[0].plot(t, d_rt0, label=r"$E_{tot}(t)-E_{tot}(0_{RT})$")
    if d_gs is not None:
        ax[0].plot(t, d_gs, label=r"$E_{tot}(t)-E_{GS}$", ls="--")
    ax[0].set_ylabel(r"$\Delta E$ (eV)"); ax[0].legend(fontsize=8)
    for c in ("energy_kinetic", "energy_hartree", "energy_xc", "energy_external", "energy_ion"):
        if c in cols:
            ax[1].plot(t, (df[c].to_numpy(dtype=float) - df[c].iloc[0]) * HA_EV, label=c, lw=1)
    ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel(r"$\Delta$comp (eV)"); ax[1].legend(fontsize=7)
    fig.tight_layout(); plot_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_png, dpi=130); plt.close(fig)

    def _f(x):
        return None if x is None else float(x)

    return {
        "columns": cols, "component_gap": gap,
        "summary": {
            "n_steps": int(len(et)), "t_final_au": _f(t[-1]),
            "dE_total_vs_rt0_final": _f(d_rt0[-1]),
            "dE_total_vs_rt0_max": _f(np.max(d_rt0)),
            "dE_total_vs_gs_final": _f(None if d_gs is None else d_gs[-1]),
            "dE_total_vs_gs_max": _f(None if d_gs is None else np.max(d_gs)),
            "crosses_zero_above": bool(np.max(d_rt0) > 1e-3),
            "N_initial": _f(None if N is None else N[0]),
            "N_final": _f(None if N is None else N[-1]),
            "eig_tracks_total": eig_tracks,
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--aim", default="")
    ap.add_argument("--method", default="")
    ap.add_argument("--run-dir", help="dir with the ablation run.cpp (build+run)")
    ap.add_argument("--mine", help="path to an EXISTING observables.csv (Phase 0)")
    ap.add_argument("--e-gs", type=float, default=None, help="E_GS in Ha for the ΔE-vs-GS convention")
    ap.add_argument("--env", action="append", default=[], help="KEY=VAL ablation env (repeatable)")
    a = ap.parse_args()

    probe_dir = PROBES / a.name
    probe_dir.mkdir(parents=True, exist_ok=True)
    plot_png = probe_dir / f"{a.name}_energy.png"

    if a.mine:
        csv = Path(a.mine)
    else:
        if not a.run_dir:
            ap.error("need --run-dir (to build+run) or --mine (existing csv)")
        env = dict(kv.split("=", 1) for kv in a.env)
        csv = build_and_run(Path(a.run_dir), probe_dir, env)

    res = summarise(csv, a.e_gs, plot_png)
    result = {
        "name": a.name, "aim": a.aim, "method": a.method,
        "run_dir": a.run_dir, "observables_csv": str(csv), "plot_png": str(plot_png),
        **res,
    }
    (probe_dir / "result.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result["summary"], indent=2))
    if res["component_gap"]:
        print(f"[component gap] missing: {res['component_gap']}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
