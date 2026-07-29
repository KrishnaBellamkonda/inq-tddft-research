#!/usr/bin/env python3
"""screening_wp_test.py — save the t=0 WP densities for the WP-potential / screening test.

Runs the WP insertion (periodicity 2, off the Lz=120 open-z GS) at two radii with
LJ_SAVE_DENSITY=1 so the rebuilt wp/run binary writes, BEFORE propagation:
  density_wp/    |psi_WP|^2   (the WP source charge, ~norm 1)
  density_total/ slab + WP
  density_bath/  slab only (screening baseline; = GS slab at t=0)

  r=12 — clean, far from the slab face  -> intrinsic orthogonalisation distortion
  r=4  — near the surface               -> slab-proximity distortion

n_WP feeds the Coulomb-potential comparison vs the classical Gaussian ghost
(Learning #2). CPU build; a couple of runs, 1 step each (+ the slow one-shot
orbital() write ~minutes). Output: runs/screening_wp/wp_r{r}_p2/ .

Launch (venv):
    cd .../campaign_autorun
    /local/data/public/skcb2/tddft/venv/bin/python3 screening_wp_test.py
"""
from __future__ import annotations
import os, subprocess, sys
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
CA = ROOT / "ResearchProject/systems/localised_jellium/scripts/campaign_autorun"
RUNS = CA / "runs"
WPBIN = str(CA / "wp/run")
GS_P2 = str(RUNS / "h2/gs_p2_lz120/checkpoint")     # periodicity-2 GS, Lz=120
RADII = [12, 4]
THREADS = int(os.environ.get("LJ_THREADS", "12"))

ENV = {**os.environ,
       "INQ_SHARE_PATH": str(ROOT / "inq/install/share"),
       "PSEUDOPOD_SHARE_PATH": str(ROOT / "inq/install/share/pseudopod"),
       "INQ_SOURCE": str(ROOT / "inq-study"),
       "OMP_NUM_THREADS": str(THREADS),
       "OPENBLAS_NUM_THREADS": str(THREADS)}

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)

def _done(rundir: Path) -> bool:
    # done = run_summary present AND the WP density VTI written
    ok_summary = any("run_completed = true" in rs.read_text()
                     for rs in rundir.glob("**/run_summary.txt") if rs.is_file())
    ok_vti = any(rundir.glob("**/density_wp/*.vti"))
    return ok_summary and ok_vti

def run_sim(rundir: Path, ov: dict, label: str) -> bool:
    rundir.mkdir(parents=True, exist_ok=True)
    if _done(rundir):
        log(f"  SKIP {label} (already complete)"); return True
    env = {**ENV, **{k: str(v) for k, v in ov.items()}}
    log(f"  RUN  {label}")
    with open(rundir / "run.log", "w") as lf:
        rc = subprocess.run([WPBIN], cwd=str(rundir), env=env,
                            stdout=lf, stderr=subprocess.STDOUT).returncode
    ok = rc == 0 and _done(rundir)
    if not ok: log(f"  FAIL {label} rc={rc} (see {rundir/'run.log'})")
    return ok

def main():
    if not Path(GS_P2).exists():
        log(f"FATAL: p2 GS missing {GS_P2}"); sys.exit(2)
    if not Path(WPBIN).exists():
        log(f"FATAL: wp binary missing {WPBIN}"); sys.exit(2)
    log(f"screening WP-potential test (CPU, {THREADS} threads); GS={GS_P2}")
    ok = 0
    for r in RADII:
        z = -(12.5 + r)
        ov = dict(LJ_OUT=f"wp_r{r}_p2", LJ_LZ=120, LJ_PERIODICITY=2, LJ_N_STEPS=1,
                  LJ_LAUNCH_Z=z, LJ_K0=0, LJ_SIGMA=0.5, LJ_GS_DIR=GS_P2,
                  LJ_SAVE_DENSITY=1)
        if run_sim(RUNS / f"screening_wp/wp_r{r}_p2", ov, f"screening wp r={r}"):
            ok += 1
    log(f"screening WP-potential test done: {ok}/{len(RADII)} complete")

if __name__ == "__main__":
    main()
