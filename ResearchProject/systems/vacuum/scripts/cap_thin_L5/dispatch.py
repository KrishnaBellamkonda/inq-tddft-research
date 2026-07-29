#!/usr/bin/env python3
"""cap_thin_L5: thin-absorber (L=5 Bohr) CAP reflectivity-tuning dispatcher.

Goal (docs/plans/cap-thin-absorber-tuning.md): find the in-built INQ CAP
(perturbations::absorbing) parameters giving low reflection error ε across
1-100 eV with the minimum near 10 eV, under a THIN absorber L=5 Bohr and SHALLOW
depth η ∈ [-0.01,-0.30] Ha. Three η-curves × 11 energies = 33 CAP runs.

Reuses the already-built cap_sweep binary (env-parameterised by CAP_L/CAP_ETA/
CAP_K0) — no rebuild. Each run emits the full free-WP minimum observable set +
manifest; ε is PROVISIONAL until Task #7 (inq-study engine ctest).

    python3 dispatch.py            # run the 33-run menu on GPUs 0,1

Env that must reach the binary: INQ_SOURCE=inq-study, share paths -> inq.
"""
import os, subprocess, sys, time, math
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
INQ_STUDY = ROOT / "inq-study"
SHARE = ROOT / "inq/install/share"
SYS = ROOT / "ResearchProject/systems/vacuum"
BINARY = SYS / "scripts/cap_sweep/run"          # reuse the existing cap_sweep binary
SWEEP_DIR = SYS / "cap_thin_L5"                  # runs grouped here (ADR 0007 amendment)
GPUS = [0, 1]

# --- the menu -----------------------------------------------------------------
L_FIXED = 5.0                                    # thin absorber, Bohr
ETAS = [-0.01, -0.05, -0.30]                     # depth (Ha, <0 absorbs); log-even over [-0.01,-0.30]
ENERGIES_EV = [1.01, 1.87, 3.48, 6.46, 7.0, 10.0, 15.0, 22.24, 41.28, 76.62, 100.0]
HA_TO_EV = 27.211386245988

def k0_of(E_eV):                                 # E = 0.5 k0^2 Ha  ->  k0 = sqrt(2E/Ha)
    return math.sqrt(2.0 * E_eV / HA_TO_EV)


def env_for(gpu):
    e = dict(os.environ)
    e["INQ_SOURCE"] = str(INQ_STUDY)
    e["INQ_SHARE_PATH"] = str(SHARE)
    e["PSEUDOPOD_SHARE_PATH"] = str(SHARE / "pseudopod")
    e["CUDA_VISIBLE_DEVICES"] = str(gpu)
    return e


def jobs():
    out = []
    for eta in ETAS:
        for E in ENERGIES_EV:
            k0 = k0_of(E)
            out.append(dict(k0=k0, E=E, L=L_FIXED, eta=eta))
    return out


def run_dir(j):
    return SWEEP_DIR / f"run_cap_k{j['k0']:.2f}_L{int(j['L'])}_eta{abs(j['eta']):.2f}"


def launch(j, gpu):
    d = run_dir(j)
    d.mkdir(parents=True, exist_ok=True)
    e = env_for(gpu)
    e.update(CAP_K0=str(j["k0"]), CAP_L=str(j["L"]), CAP_ETA=str(j["eta"]),
             CAP_OUTDIR=str(d / "results"))
    log = open(d / "run.log", "w")
    return subprocess.Popen([str(BINARY)], cwd=d, env=e, stdout=log, stderr=subprocess.STDOUT)


def main():
    if not BINARY.exists():
        print(f"BINARY missing: {BINARY} (build cap_sweep first)", file=sys.stderr); sys.exit(1)
    only = sys.argv[1] if len(sys.argv) > 1 else None     # optional substring filter (smoke)
    pending = [j for j in jobs() if (only is None or only in run_dir(j).name)]
    print(f"==> {len(pending)} runs pending (L=5 thin CAP tuning)", flush=True)
    running = {}
    while pending or running:
        for gpu in GPUS:
            if gpu not in running and pending:
                j = pending.pop(0)
                print(f"  launch {run_dir(j).name} (E={j['E']}eV) on GPU{gpu}", flush=True)
                running[gpu] = (launch(j, gpu), j, time.time())
        for gpu in list(running):
            proc, j, t0 = running[gpu]
            if proc.poll() is not None:
                dt = time.time() - t0
                print(f"  done   {run_dir(j).name} rc={proc.returncode} ({dt:.0f}s)", flush=True)
                del running[gpu]
        time.sleep(2)
    print("==> cap_thin_L5 sweep complete", flush=True)


if __name__ == "__main__":
    main()
