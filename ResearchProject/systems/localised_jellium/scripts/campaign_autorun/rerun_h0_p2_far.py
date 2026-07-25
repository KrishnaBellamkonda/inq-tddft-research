#!/usr/bin/env python3
"""rerun_h0_p2_far.py — extended-r H0 sweep at periodicity 2 in a BIGGER box (Lz=200).

The Lz=120 p2 sweep only reaches r~40 before the projectile nears the box edge (the
classical excess there is still ~12 eV, not 0). This extends r to ~76 Bohr in an
Lz=200 open-z box, off its own GS (runs/h0_p2_far/gs_p2_lz200/checkpoint), streaming
the full energy decomposition. Idempotent; output runs/h0_p2_far/{tag}_r{r}_p2/.

Launch (after the Lz=200 GS exists):
    cd .../campaign_autorun
    GPU=1 /local/data/public/skcb2/tddft/venv/bin/python3 rerun_h0_p2_far.py
"""
from __future__ import annotations
import os, subprocess, sys
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
CA = ROOT / "ResearchProject/systems/localised_jellium/scripts/campaign_autorun"
RUNS = CA / "runs"
WPBIN, CLBIN = str(CA / "wp/run"), str(CA / "classical/run")
GS_FAR = str(RUNS / "h0_p2_far/gs_p2_lz200/checkpoint")
GPU = os.environ.get("GPU", "1")
LZ = 200
RADII = (4, 12, 20, 28, 36, 44, 52, 60, 68, 76)   # extends past the Lz=120 max (40)

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

def run_sim(binary, rundir: Path, ov: dict, label: str) -> bool:
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
    if not Path(GS_FAR).exists():
        log(f"FATAL: Lz=200 p2 GS checkpoint missing: {GS_FAR}"); sys.exit(2)
    log(f"H0 periodicity-2 FAR sweep (Lz={LZ}, GPU={GPU}); GS={GS_FAR}")
    for tag, binary in (("wp", WPBIN), ("cl", CLBIN)):
        for r in RADII:
            z = -(12.5 + r)
            ov = dict(LJ_OUT=f"{tag}_r{r}_p2", LJ_LZ=LZ, LJ_PERIODICITY=2,
                      LJ_LAUNCH_Z=z, LJ_GS_DIR=GS_FAR)
            if tag == "wp":
                ov.update(LJ_K0=0, LJ_SIGMA=0.5)
            run_sim(binary, RUNS / f"h0_p2_far/{tag}_r{r}_p2", ov, f"far {tag} r={r}")
    log("H0 periodicity-2 FAR sweep done")

if __name__ == "__main__":
    main()
