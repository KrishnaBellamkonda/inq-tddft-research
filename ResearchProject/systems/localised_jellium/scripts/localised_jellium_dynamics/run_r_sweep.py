#!/usr/bin/env python3
"""Radius sweep for the energy-decomposition r-independence proof.

For each r (projectile-slab-face distance) generate a complete at-rest twin pair
(Gaussian-charge classical + WP, sigma_WP=0.5, periodicity 2, Lz=120), then run
twin_decompose to extract R (residual = WP self-Hartree) and SIE. No density
frames (SAVE_EVERY=0) and only 2 steps -> ~1 min/run. Writes sweep_R_SIE.csv.
"""
import os, re, subprocess, csv
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

DYN = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
           "scripts/localised_jellium_dynamics")
GS_P2 = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
         "scripts/campaign_autorun/runs/h2/gs_p2_lz120/checkpoint")
PY = "/local/data/public/skcb2/tddft/venv/bin/python3"
DECOMP = "/local/data/public/skcb2/tddft/.claude/skills/twin-run-analysis/twin_decompose.py"

RADII = [4, 12, 20, 28, 36, 40]
SWEEP = DYN / "runs" / "twin_ec_rsweep"
SWEEP.mkdir(parents=True, exist_ok=True)

ENV = dict(os.environ,
           INQ_SHARE_PATH="/local/data/public/skcb2/tddft/inq/install/share",
           PSEUDOPOD_SHARE_PATH="/local/data/public/skcb2/tddft/inq/install/share/pseudopod")

def z_of(r): return -(12.5 + r)

def common(r):
    return dict(LJ_LZ="120", LJ_PERIODICITY="2", LJ_LAUNCH_Z=str(z_of(r)),
                LJ_K0="0", LJ_SIGMA="0.5", LJ_GS_DIR=GS_P2,
                LJ_N_STEPS="2", LJ_DT="0.01", LJ_SAVE_EVERY="0")

def launch(r, side, gpu):
    binary = str(DYN / ("phase5_wp/run" if side == "wp" else "proj_dyn/run"))
    out = f"r{r}_{side}"
    env = dict(ENV, CUDA_VISIBLE_DEVICES=str(gpu), LJ_OUT=out, **common(r))
    if side == "classical":
        env["LJ_MASS"] = "1"
    logf = SWEEP / f"{out}.log"
    with open(logf, "w") as lf:
        rc = subprocess.run([binary], cwd=str(SWEEP), env=env,
                            stdout=lf, stderr=subprocess.STDOUT).returncode
    return r, side, rc

# Launch all 12 runs, 2 at a time (GPU0/GPU1).
jobs = [(r, side) for r in RADII for side in ("wp", "classical")]
def worker(i_job):
    i, (r, side) = i_job
    return launch(r, side, gpu=i % 2)
with ThreadPoolExecutor(max_workers=2) as ex:
    results = list(ex.map(worker, list(enumerate(jobs))))
for r, side, rc in results:
    print(f"r={r:3d} {side:9s} rc={rc}")

# Decompose each pair, parse R and SIE.
def decompose(r):
    wp = SWEEP / "results" / f"r{r}_wp"
    cl = SWEEP / "results" / f"r{r}_classical"
    out = subprocess.run([PY, DECOMP, "--wp", str(wp), "--classical", str(cl)],
                         capture_output=True, text=True, env=ENV).stdout
    def grab(pat):
        m = re.search(pat, out)
        return float(m.group(1)) if m else float("nan")
    R = grab(r"residual R = d\(E_H\+E_ext\) - U_proj_bg\s+([-\d.]+)")
    sie = grab(r"SIE = R \+ dXC\s+([-\d.]+)")
    dkin = grab(r"dKin_localisation.*?\s([-\d.]+)\s+81")
    dxc = grab(r"dXC \(xc surplus\)\s+([-\d.]+)")
    return dict(r=r, R=R, SIE=sie, dKin=dkin, dXC=dxc)

rows = [decompose(r) for r in RADII]
csv_path = SWEEP / "sweep_R_SIE.csv"
with open(csv_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["r", "R", "SIE", "dKin", "dXC"])
    w.writeheader(); w.writerows(rows)
print("\nwrote", csv_path)
for row in rows:
    print(row)
