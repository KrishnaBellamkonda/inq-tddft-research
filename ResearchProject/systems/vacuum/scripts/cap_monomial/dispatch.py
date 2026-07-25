#!/usr/bin/env python3
"""cap_monomial: benchmark the inq-study MONOMIAL CAP vs the built-in sin² at L=5.

Fixed thin L=5 Bohr, E=10 eV (k0≈0.857). Sweep monomial order n × depth η to ask:
does the ramp V=i·eta·s^n beat the sin² hump's ~0.20 reflection floor at L=5?
Compare against the sin² L5/E10 points already in cap_thin_L5/.

Build the cap_monomial binary ONCE against inq-study (it #includes the new
perturbations::absorbing_monomial header), then dispatch. ε PROVISIONAL until Task #7.

    python3 dispatch.py            # build once (if needed) + 16 runs on GPUs 0,1
"""
import os, subprocess, sys, time, math
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
INQ_STUDY = ROOT / "inq-study"
SHARE = ROOT / "inq/install/share"
SYS = ROOT / "ResearchProject/systems/vacuum"
SRC = SYS / "scripts/cap_monomial/run.cpp"
BINARY = SYS / "scripts/cap_monomial/run"
SWEEP_DIR = SYS / "cap_monomial"
GPUS = [0, 1]

E_FIXED_EV = 10.0
K0 = math.sqrt(2.0 * E_FIXED_EV / 27.211386245988)   # ≈0.857
L_FIXED = 5.0
ORDERS = [1, 2, 3, 4]
ETAS   = [-0.10, -0.20, -0.30, -0.50]


def env_for(gpu):
    e = dict(os.environ)
    e["INQ_SOURCE"] = str(INQ_STUDY)
    e["INQ_SHARE_PATH"] = str(SHARE)
    e["PSEUDOPOD_SHARE_PATH"] = str(SHARE / "pseudopod")
    e["CUDA_VISIBLE_DEVICES"] = str(gpu)
    return e


def build_once():
    if BINARY.exists():
        print(f"==> binary present ({BINARY}); skip build", flush=True); return
    print("==> building cap_monomial against inq-study (first build is full)...", flush=True)
    e = env_for(GPUS[0])
    r = subprocess.run(["inq-run", "--reconfig"], cwd=SRC.parent, env={**e,
                        "CAP_K0": str(K0), "CAP_L": "5", "CAP_ETA": "-0.3",
                        "CAP_ORDER": "2", "CAP_OUTDIR": "build_smoke"})
    if r.returncode != 0:
        print("BUILD FAILED", file=sys.stderr); sys.exit(1)
    print("==> build OK", flush=True)


def jobs():
    return [dict(k0=K0, L=L_FIXED, eta=eta, order=n) for n in ORDERS for eta in ETAS]


def run_dir(j):
    return SWEEP_DIR / f"run_mono_k{j['k0']:.2f}_L{int(j['L'])}_eta{abs(j['eta']):.2f}_n{j['order']}"


def launch(j, gpu):
    d = run_dir(j); d.mkdir(parents=True, exist_ok=True)
    e = env_for(gpu)
    e.update(CAP_K0=str(j["k0"]), CAP_L=str(j["L"]), CAP_ETA=str(j["eta"]),
             CAP_ORDER=str(j["order"]), CAP_OUTDIR=str(d / "results"))
    log = open(d / "run.log", "w")
    return subprocess.Popen([str(BINARY)], cwd=d, env=e, stdout=log, stderr=subprocess.STDOUT)


def main():
    build_once()
    pending = jobs()
    print(f"==> {len(pending)} monomial runs pending (L=5, E=10 eV, n×η)", flush=True)
    running = {}
    while pending or running:
        for gpu in GPUS:
            if gpu not in running and pending:
                j = pending.pop(0)
                print(f"  launch {run_dir(j).name} on GPU{gpu}", flush=True)
                running[gpu] = (launch(j, gpu), j, time.time())
        for gpu in list(running):
            proc, j, t0 = running[gpu]
            if proc.poll() is not None:
                print(f"  done   {run_dir(j).name} rc={proc.returncode} ({time.time()-t0:.0f}s)", flush=True)
                del running[gpu]
        time.sleep(2)
    print("==> cap_monomial sweep complete", flush=True)


if __name__ == "__main__":
    main()
