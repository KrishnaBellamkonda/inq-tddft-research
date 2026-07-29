#!/usr/bin/env python3
"""Depth-sweep dispatcher for the engine-integrated CAP (perturbations::absorbing).

Builds the cap_sweep binary ONCE against inq-study (the fork carrying the
complexified scalar potential), then runs a small η (depth) sweep — CAP at the
far end of the box, WP travelling toward it — writing FLAT run dirs
`systems/vacuum/run_cap_<...>/` (ADR 0007), each with the full minimum observable
suite. One showcase run writes density VTI frames for the report GIF.

Run AFTER Phase-1 validation (cap_probe compiles + absorbs against inq-study).
    python3 dispatch.py            # build once + run the menu on 2 GPUs

Env that must reach inq-run (set here): INQ_SOURCE=inq-study, share paths -> inq.
"""
import os, subprocess, shutil, sys, time
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
INQ_STUDY = ROOT / "inq-study"
SHARE = ROOT / "inq/install/share"
SYS = ROOT / "ResearchProject/systems/vacuum"
BIN_SRC = SYS / "scripts/cap_sweep/run.cpp"
BUILD_DIR = SYS / "scripts/cap_sweep/build"
GPUS = [0, 1]

# --- the depth menu (editable after Phase-1 observations) -------------------
# Fixed E≈22 eV (k0≈1.28), L=20; sweep CAP depth η (Ha, <0 absorbs).
K0_MAIN, L_MAIN, ETA_MID = 1.28, 20.0, -0.50    # E≈22 eV; η near the absorption sweet spot
# Thorough investigation of the in-built CAP's three knobs (every run = full
# free-WP minimum observable set + manifest; all runs emit density_wp VTI).
ETAS  = [-0.01, -0.02, -0.05, -0.10, -0.25, -0.50, -1.00, -2.00, -4.00]  # depth, at L_MAIN
LS    = [5.0, 10.0, 20.0, 30.0, 40.0, 50.0]                              # width, at ETA_MID
K0S   = [0.86, 1.28, 2.00, 2.71]                                         # energy (~10,22,54,100 eV), at ETA_MID,L_MAIN


def env_for(gpu):
    e = dict(os.environ)
    e["INQ_SOURCE"] = str(INQ_STUDY)
    e["INQ_SHARE_PATH"] = str(SHARE)
    e["PSEUDOPOD_SHARE_PATH"] = str(SHARE / "pseudopod")
    e["CUDA_VISIBLE_DEVICES"] = str(gpu)
    return e


def jobs():
    seen, out = set(), []
    def add(k0, L, eta):
        key = (round(k0, 2), round(L, 1), round(eta, 3))
        if key not in seen:
            seen.add(key); out.append(dict(k0=k0, L=L, eta=eta))
    for eta in ETAS: add(K0_MAIN, L_MAIN, eta)     # depth sweep
    for L in LS:     add(K0_MAIN, L, ETA_MID)      # width sweep
    for k0 in K0S:   add(k0, L_MAIN, ETA_MID)      # energy sweep
    return out


def run_dir(j):
    return SYS / f"run_cap_k{j['k0']:.2f}_L{int(j['L'])}_eta{abs(j['eta']):.2f}"


BINARY = BIN_SRC.parent / "run"          # inq-run places the binary next to the .cpp


def build_once():
    """Build the cap_sweep binary against inq-study (skip if already built)."""
    if BINARY.exists():
        print(f"==> binary present ({BINARY}); skipping build", flush=True)
        return
    print("==> building cap_sweep against inq-study (first build is full)...", flush=True)
    e = env_for(GPUS[0])
    r = subprocess.run(["inq-run", "--reconfig"], cwd=BIN_SRC.parent, env={**e,
                        "CAP_K0": "1.0", "CAP_L": "5", "CAP_ETA": "-0.5", "CAP_OUTDIR": "build_smoke"})
    if r.returncode != 0:
        print("BUILD FAILED", file=sys.stderr); sys.exit(1)
    print("==> build OK", flush=True)


def launch(j, gpu):
    d = run_dir(j)
    d.mkdir(parents=True, exist_ok=True)
    binary = BINARY
    e = env_for(gpu)
    e.update(CAP_K0=str(j["k0"]), CAP_L=str(j["L"]), CAP_ETA=str(j["eta"]),
             CAP_OUTDIR=str(d / "results"))
    log = open(d / "run.log", "w")
    return subprocess.Popen([str(binary)], cwd=d, env=e, stdout=log, stderr=subprocess.STDOUT)


def main():
    build_once()
    pending = jobs()
    running = {}   # gpu -> (proc, job)
    while pending or running:
        for gpu in GPUS:
            if gpu not in running and pending:
                j = pending.pop(0)
                print(f"  launch {run_dir(j).name} on GPU{gpu}", flush=True)
                running[gpu] = (launch(j, gpu), j)
        for gpu in list(running):
            proc, j = running[gpu]
            if proc.poll() is not None:
                print(f"  done   {run_dir(j).name} rc={proc.returncode}", flush=True)
                del running[gpu]
        time.sleep(2)
    print("==> depth sweep complete", flush=True)


if __name__ == "__main__":
    main()
