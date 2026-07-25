#!/usr/bin/env python3
"""rerun_h0_p3.py — full-component PERIODICITY-3 mirror of rerun_h0_p2.py.

Motivation
----------
The h0_p2 runs stream the full energy decomposition but sit at periodicity 2,
whose Poisson G=0 term is 0.5*rc^2 (poisson.hpp:49). For a net-charged electron
density this inflates E_hartree and E_external by a large, box-dependent constant
of opposite sign — so the INDIVIDUAL Hartree/external components are reference
artefacts and only their SUM is physical. That makes "match the analytic Gaussian
self-Hartree against E_hartree(WP) - E_hartree(classical)" impossible in p2.

Periodicity 3 sets V(G=0) = 0 (poisson.hpp:31, zeroterm defaults to 0) — no
0.5*rc^2 inflation — so each of E_hartree and E_external is individually physical.
This script reruns the SAME radius sweep as h0_p2, at periodicity 3, with the
p3 ground state, streaming ALL energy components. Output: runs/h0_p3/.

The original h0 (p3) runs (runs/h0/) predate the full-component streaming and
carry only {total,kinetic,hartree,xc}; this is their full-component remake.

Backends (physics is backend-identical; verified per-row by sum(parts)==total):
  * wp/run is a CUDA binary  -> GPU  (CUDA_VISIBLE_DEVICES)
  * classical/run is a CPU binary -> CPU (OMP/OPENBLAS threads)

Single-point runs. Idempotent: skips completed runs.

Launch:
    cd .../campaign_autorun
    GPU=1 /local/data/public/skcb2/tddft/venv/bin/python3 rerun_h0_p3.py
"""
from __future__ import annotations
import os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
LJ = ROOT / "ResearchProject/systems/localised_jellium"
CA = LJ / "scripts/campaign_autorun"
RUNS = CA / "runs"
WPBIN, CLBIN = str(CA / "wp/run"), str(CA / "classical/run")
GS_P3_DIR = str(LJ / "shared_gs/slab_n82_L50x50x120")    # periodicity-3 GS (dir IS the checkpoint)
GPU = os.environ.get("GPU", "1")
RADII = (4, 12, 20, 28, 36, 40)                          # same set as h0_p2
CL_WORKERS = 4                                            # concurrent CPU classical runs
CL_THREADS = 8                                            # threads per classical run

BASE = {**os.environ,
        "INQ_SHARE_PATH": str(ROOT / "inq/install/share"),
        "PSEUDOPOD_SHARE_PATH": str(ROOT / "inq/install/share/pseudopod"),
        "INQ_SOURCE": str(ROOT / "inq-study")}

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)

def _done(rundir: Path) -> bool:
    for rs in rundir.glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs.read_text():
                return True
        except Exception:
            pass
    return False

def run_sim(binary: str, rundir: Path, ov: dict, extra_env: dict, label: str) -> bool:
    rundir.mkdir(parents=True, exist_ok=True)
    if _done(rundir):
        log(f"  SKIP {label} (already complete)"); return True
    env = {**BASE, **extra_env, **{k: str(v) for k, v in ov.items()}}
    log(f"  RUN  {label}")
    with open(rundir / "run.log", "w") as lf:
        rc = subprocess.run([binary], cwd=str(rundir), env=env,
                            stdout=lf, stderr=subprocess.STDOUT).returncode
    ok = rc == 0 and _done(rundir)
    log(f"  {'OK  ' if ok else 'FAIL'} {label}" + ("" if ok else f" rc={rc} (see {rundir/'run.log'})"))
    return ok

def main():
    if not Path(GS_P3_DIR).exists():
        log(f"FATAL: periodicity-3 GS missing: {GS_P3_DIR}"); sys.exit(2)
    log(f"H0 periodicity-3 full-component mirror; wp->GPU {GPU}, cl->CPU; GS={GS_P3_DIR}")

    # --- WP runs: GPU, sequential (single device) ---
    wp_env = {"CUDA_VISIBLE_DEVICES": GPU}
    for r in RADII:
        z = -(12.5 + r)
        ov = dict(LJ_OUT=f"wp_r{r}_p3", LJ_LZ=120, LJ_PERIODICITY=3, LJ_LAUNCH_Z=z,
                  LJ_GS_DIR=GS_P3_DIR, LJ_K0=0, LJ_SIGMA=0.5)
        run_sim(WPBIN, RUNS / f"h0_p3/wp_r{r}_p3", ov, wp_env, f"H0-p3 wp r={r}")

    # --- classical runs: CPU, concurrent ---
    cl_env = {"OMP_NUM_THREADS": str(CL_THREADS), "OPENBLAS_NUM_THREADS": str(CL_THREADS),
              "CUDA_VISIBLE_DEVICES": ""}
    def cl_task(r):
        z = -(12.5 + r)
        ov = dict(LJ_OUT=f"cl_r{r}_p3", LJ_LZ=120, LJ_PERIODICITY=3, LJ_LAUNCH_Z=z,
                  LJ_GS_DIR=GS_P3_DIR)
        return run_sim(CLBIN, RUNS / f"h0_p3/cl_r{r}_p3", ov, cl_env, f"H0-p3 cl r={r}")
    with ThreadPoolExecutor(max_workers=CL_WORKERS) as ex:
        list(ex.map(cl_task, RADII))

    log("H0 periodicity-3 full-component mirror done")

if __name__ == "__main__":
    main()
