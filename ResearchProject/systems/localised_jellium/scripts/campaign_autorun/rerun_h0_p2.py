#!/usr/bin/env python3
"""rerun_h0_p2.py — re-run the H0 insertion-energy experiment at PERIODICITY 2.

H0 originally ran at full PBC (periodicity 3). This re-runs the same WP + classical
radius sweep with:
  * LJ_PERIODICITY = 2 (open-z slab; removes the full-PBC G=0/background ambiguity)
  * the periodicity-2 GS checkpoint (runs/h2/gs_p2_lz120/checkpoint)
  * the extended observables schema (all energy components — external, non_local,
    ion, ion_kinetic, exact_exchange, nvxc, eigenvalues — streamed each step)

Single-point 3-step runs → minutes on GPU. Idempotent: skips completed runs.
Output: runs/h0_p2/{tag}_r{r}_p2/ .

Launch:
    cd .../campaign_autorun
    GPU=1 /local/data/public/skcb2/tddft/venv/bin/python3 rerun_h0_p2.py
"""
from __future__ import annotations
import os, subprocess, sys
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
LJ = ROOT / "ResearchProject/systems/localised_jellium"
CA = LJ / "scripts/campaign_autorun"
RUNS = CA / "runs"
WPBIN, CLBIN = str(CA / "wp/run"), str(CA / "classical/run")
GS_P2_CKPT = str(RUNS / "h2/gs_p2_lz120/checkpoint")     # periodicity-2 GS
GPU = os.environ.get("GPU", "1")
RADII = (4, 12, 20, 28, 36, 40)                          # same set as original H0

ENV = {**os.environ,
       "INQ_SHARE_PATH": str(ROOT / "inq/install/share"),
       "PSEUDOPOD_SHARE_PATH": str(ROOT / "inq/install/share/pseudopod"),
       "INQ_SOURCE": str(ROOT / "inq-study"),
       "CUDA_VISIBLE_DEVICES": GPU}

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)

def _done(rundir: Path) -> bool:
    for rs in rundir.glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs.read_text():
                return True
        except Exception:
            pass
    return False

def run_sim(binary: str, rundir: Path, ov: dict, label: str) -> bool:
    rundir.mkdir(parents=True, exist_ok=True)
    if _done(rundir):
        log(f"  SKIP {label} (already complete)"); return True
    env = {**ENV, **{k: str(v) for k, v in ov.items()}}
    log(f"  RUN  {label}")
    with open(rundir / "run.log", "w") as lf:
        rc = subprocess.run([binary], cwd=str(rundir), env=env,
                            stdout=lf, stderr=subprocess.STDOUT).returncode
    ok = rc == 0 and _done(rundir)
    if not ok:
        log(f"  FAIL {label} rc={rc} (see {rundir/'run.log'})")
    return ok

def main():
    if not Path(GS_P2_CKPT).exists():
        log(f"FATAL: periodicity-2 GS checkpoint missing: {GS_P2_CKPT}"); sys.exit(2)
    log(f"H0 periodicity-2 re-run (GPU={GPU}); GS={GS_P2_CKPT}")
    for tag, binary in (("wp", WPBIN), ("cl", CLBIN)):
        for r in RADII:
            z = -(12.5 + r)
            ov = dict(LJ_OUT=f"{tag}_r{r}_p2", LJ_LZ=120, LJ_PERIODICITY=2,
                      LJ_LAUNCH_Z=z, LJ_GS_DIR=GS_P2_CKPT)
            if tag == "wp":
                ov.update(LJ_K0=0, LJ_SIGMA=0.5)
            run_sim(binary, RUNS / f"h0_p2/{tag}_r{r}_p2", ov, f"H0-p2 {tag} r={r}")
    log("H0 periodicity-2 re-run done")

if __name__ == "__main__":
    main()
