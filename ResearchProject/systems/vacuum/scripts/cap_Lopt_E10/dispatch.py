#!/usr/bin/env python3
"""cap_Lopt_E10: built-in sin² CAP — how THIN can we go at E=10 eV?

Baseline for the reflectivity-floor question (docs/plans/cap-monomial-inq-study.md):
fix E=10 eV (k0≈0.857) and 2D-sweep width L × depth η with the in-built
perturbations::absorbing, to find the smallest L (and best η) that reaches a low
reflection error ε. Reuses the cap_sweep binary (env-parameterised) — no rebuild.
ε PROVISIONAL until Task #7.

    python3 dispatch.py            # 15 runs (5 L × 3 η) on GPUs 0,1
"""
import os, subprocess, sys, time, math
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
INQ_STUDY = ROOT / "inq-study"
SHARE = ROOT / "inq/install/share"
SYS = ROOT / "ResearchProject/systems/vacuum"
BINARY = SYS / "scripts/cap_sweep/run"          # reuse built-in sin² CAP binary
SWEEP_DIR = SYS / "cap_Lopt_E10"
GPUS = [0, 1]

E_FIXED_EV = 10.0
K0 = math.sqrt(2.0 * E_FIXED_EV / 27.211386245988)   # ≈0.857
LS   = [6.0, 8.0, 10.0, 12.0, 15.0]
ETAS = [-0.30, -0.50, -1.00]


def env_for(gpu):
    e = dict(os.environ)
    e["INQ_SOURCE"] = str(INQ_STUDY)
    e["INQ_SHARE_PATH"] = str(SHARE)
    e["PSEUDOPOD_SHARE_PATH"] = str(SHARE / "pseudopod")
    e["CUDA_VISIBLE_DEVICES"] = str(gpu)
    return e


def jobs():
    return [dict(k0=K0, L=L, eta=eta) for eta in ETAS for L in LS]


def run_dir(j):
    return SWEEP_DIR / f"run_cap_k{j['k0']:.2f}_L{int(j['L'])}_eta{abs(j['eta']):.2f}"


def launch(j, gpu):
    d = run_dir(j); d.mkdir(parents=True, exist_ok=True)
    e = env_for(gpu)
    e.update(CAP_K0=str(j["k0"]), CAP_L=str(j["L"]), CAP_ETA=str(j["eta"]),
             CAP_OUTDIR=str(d / "results"))
    log = open(d / "run.log", "w")
    return subprocess.Popen([str(BINARY)], cwd=d, env=e, stdout=log, stderr=subprocess.STDOUT)


def main():
    if not BINARY.exists():
        print(f"BINARY missing: {BINARY}", file=sys.stderr); sys.exit(1)
    pending = jobs()
    print(f"==> {len(pending)} runs (E=10 eV, L×η) pending", flush=True)
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
    print("==> cap_Lopt_E10 sweep complete", flush=True)
    autobuild_notebook()


def autobuild_notebook():
    """Run-machinery auto-build (notebook-making skill): rebuild the sweep study
    notebook ONCE at end of batch — fires regardless of who launched the run."""
    VENV = "/local/data/public/skcb2/tddft/venv/bin/python3"
    STACK = "/local/data/public/skcb2/tddft/inq-stack/python"
    builder = SYS / "hypotheses/cap_Lopt_E10/build_Lopt_report.py"
    print(f"==> auto-build notebook: {builder.name}", flush=True)
    subprocess.run([VENV, str(builder)],
                   env={**os.environ, "PYTHONPATH": STACK}, check=False)


if __name__ == "__main__":
    main()
