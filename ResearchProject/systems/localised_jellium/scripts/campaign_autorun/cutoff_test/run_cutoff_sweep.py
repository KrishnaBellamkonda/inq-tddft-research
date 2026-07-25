#!/usr/bin/env python3
"""run_cutoff_sweep.py — classical projectile ΔE_total(r) for several UPF radial
cutoffs, to test whether the classical decay is set by the projectile potential's
finite range.

For each truncated UPF (r_cut in {10,20,30,40} Bohr, made by make_cutoff_upfs.py),
run the classical ghost at r in {2..40} Bohr, periodicity 2, off the Lz=120 open-z GS,
streaming the full energy decomposition. Output: runs/cutoff_test/rc{rc}/cl_r{r}_p2/.

Launch (after make_cutoff_upfs.py + a rebuilt classical binary):
    cd .../campaign_autorun
    GPU=1 /local/data/public/skcb2/tddft/venv/bin/python3 cutoff_test/run_cutoff_sweep.py
"""
from __future__ import annotations
import os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
CA = ROOT / "ResearchProject/systems/localised_jellium/scripts/campaign_autorun"
RUNS = CA / "runs"
CLBIN = str(CA / "classical/run")
UPFS = CA / "cutoff_test/upfs"
GS_P2 = str(RUNS / "h2/gs_p2_lz120/checkpoint")     # periodicity-2 GS, Lz=120
CUTOFFS = [10, 20, 30, 40]
RADII = [2, 8, 16, 24, 32, 40]   # 6 points: enough to resolve each cutoff's decay (CPU is slow)
# CPU run: this binary is a --cpu build. Cap threads per run and run several
# concurrently to use the idle cores (THREADS*CONCURRENCY should stay < ~55 free).
THREADS = int(os.environ.get("LJ_THREADS", "8"))
CONCURRENCY = int(os.environ.get("LJ_CONCURRENCY", "6"))

ENV = {**os.environ,
       "INQ_SHARE_PATH": str(ROOT / "inq/install/share"),
       "PSEUDOPOD_SHARE_PATH": str(ROOT / "inq/install/share/pseudopod"),
       "INQ_SOURCE": str(ROOT / "inq-study"),
       "OMP_NUM_THREADS": str(THREADS),
       "OPENBLAS_NUM_THREADS": str(THREADS)}

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)

def _done(rundir: Path) -> bool:
    for rs in rundir.glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs.read_text(): return True
        except Exception: pass
    return False

def run_sim(rundir: Path, ov: dict, label: str) -> bool:
    rundir.mkdir(parents=True, exist_ok=True)
    if _done(rundir):
        log(f"  SKIP {label}"); return True
    env = {**ENV, **{k: str(v) for k, v in ov.items()}}
    log(f"  RUN  {label}")
    with open(rundir / "run.log", "w") as lf:
        rc = subprocess.run([CLBIN], cwd=str(rundir), env=env, stdout=lf, stderr=subprocess.STDOUT).returncode
    ok = rc == 0 and _done(rundir)
    if not ok: log(f"  FAIL {label} rc={rc}")
    return ok

def main():
    if not Path(GS_P2).exists(): log(f"FATAL: p2 GS missing {GS_P2}"); sys.exit(2)
    for rc in CUTOFFS:
        upf = UPFS / f"electron_gaussian_wpsigma0p5_rc{rc}.upf"
        if not upf.exists(): log(f"FATAL: UPF missing {upf}"); sys.exit(2)
    jobs = []
    for rc in CUTOFFS:
        upf = str(UPFS / f"electron_gaussian_wpsigma0p5_rc{rc}.upf")
        for r in RADII:
            z = -(12.5 + r)
            ov = dict(LJ_OUT=f"cl_r{r}_p2", LJ_LZ=120, LJ_PERIODICITY=2, LJ_N_STEPS=1,
                      LJ_LAUNCH_Z=z, LJ_GS_DIR=GS_P2, LJ_PROJ_UPF=upf)
            jobs.append((RUNS / f"cutoff_test/rc{rc}/cl_r{r}_p2", ov, f"rc{rc} r={r}"))
    log(f"cutoff sweep on CPU: {len(jobs)} runs, {CONCURRENCY} concurrent x {THREADS} threads")
    ok = 0
    with ThreadPoolExecutor(max_workers=CONCURRENCY) as ex:
        futs = {ex.submit(run_sim, rd, ov, lbl): lbl for rd, ov, lbl in jobs}
        for f in as_completed(futs):
            if f.result(): ok += 1
    log(f"cutoff sweep done: {ok}/{len(jobs)} complete")

if __name__ == "__main__":
    main()
