#!/usr/bin/env python3
"""dispatch.py — loss-function feasibility 9-run matrix (campaign cap-jellium-loss-function).

Drives the env-driven ./run binary (built ONCE against inq-study; see run.cpp header)
over 3 modes x 3 energies. The binary must already be built:

    cd scripts/cap_loss_function
    INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study \
    INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share \
    PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod \
    inq-run run.cpp          # produces build/run

Then launch:
    python3 dispatch.py --smoke              # 20 steps each, all 9 (readiness check)
    python3 dispatch.py --gpu 1              # PRODUCTION: T~2000 a.u. (~100k steps) each
    python3 dispatch.py --gpu 1 --modes kick --energies E15   # subset

GPU: pass --gpu N (default auto-pick a free device via cudaMemGetInfo; NVML/nvidia-smi
is broken on this host, so we never rely on it). Production runs are HEAVY
(~12 h each per the campaign cost note) — confirm GPU is free first.
"""
from __future__ import annotations
import argparse, ctypes, os, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
HA2EV = 27.211386

# Energy -> velocity v=sqrt(2E/Ha2eV). E15/E20/E30 are the locked low-v window.
ENERGIES = {"E15": 15.0, "E20": 20.0, "E30": 30.0}
MODES = ("classical", "wp", "kick")

# Production I/O (campaign <duration_and_io>): dt=0.020 (in Cfg), T~2000 a.u.
PROD_N_STEPS = 100_000          # ~2000 a.u.
PROD_WRITE_EVERY = 200          # n_q source VTI every 200 steps = 4 a.u. (>=500 frames)
SMOKE_N_STEPS = 20
SMOKE_WRITE_EVERY = 2

ENV_BUILD = {
    "INQ_SOURCE": "/local/data/public/skcb2/tddft/inq-study",
    "INQ_SHARE_PATH": "/local/data/public/skcb2/tddft/inq/install/share",
    "PSEUDOPOD_SHARE_PATH": "/local/data/public/skcb2/tddft/inq/install/share/pseudopod",
}


def free_gpu() -> int:
    """Pick the GPU with the most free memory via cudaMemGetInfo (NVML-independent)."""
    try:
        cuda = ctypes.CDLL("libcudart.so")
        n = ctypes.c_int()
        cuda.cudaGetDeviceCount(ctypes.byref(n))
        best, best_free = 0, -1
        for i in range(n.value):
            cuda.cudaSetDevice(i)
            free, total = ctypes.c_size_t(), ctypes.c_size_t()
            cuda.cudaMemGetInfo(ctypes.byref(free), ctypes.byref(total))
            if free.value > best_free:
                best, best_free = i, free.value
        print(f"[dispatch] auto-picked GPU {best} ({best_free/1e9:.1f} GB free)")
        return best
    except Exception as e:
        print(f"[dispatch] GPU probe failed ({e}); defaulting to GPU 0")
        return 0


def v0(energy_key: str) -> float:
    return (2.0 * ENERGIES[energy_key] / HA2EV) ** 0.5


def run_one(mode: str, ekey: str, gpu: int, smoke: bool) -> int:
    binary = HERE / "run"            # inq-run emits the binary in RUN_DIR, not build/
    if not binary.exists():
        sys.exit(f"FATAL: {binary} not built. See module docstring for the build command.")
    sub = f"{'smoke_' if smoke else 'run_'}{mode}_{ekey}"
    env = {**os.environ, **ENV_BUILD,
           "CUDA_VISIBLE_DEVICES": str(gpu),
           "CAP_MODE": mode,
           "CAP_V0": f"{v0(ekey):.6f}",
           "CAP_OUT_SUBDIR": sub,
           "CAP_N_STEPS": str(SMOKE_N_STEPS if smoke else PROD_N_STEPS),
           "CAP_WRITE_EVERY": str(SMOKE_WRITE_EVERY if smoke else PROD_WRITE_EVERY)}
    log = HERE / f"{sub}.log"
    print(f"[dispatch] {sub}: v0={v0(ekey):.4f} steps={env['CAP_N_STEPS']} gpu={gpu} -> {log.name}")
    t0 = time.time()
    with open(log, "w") as fh:
        st = subprocess.run([str(binary)], cwd=HERE, env=env, stdout=fh, stderr=subprocess.STDOUT).returncode
    print(f"[dispatch]   exit={st} wall={time.time()-t0:.0f}s")
    # per-run email if the helper exists (campaign: per-phase Gmail)
    emailer = HERE / "email_run.py"
    if emailer.exists():
        subprocess.run([sys.executable, str(emailer), sub, f"loss-function {mode} {ekey}",
                        str(st), "cap_loss_function"], cwd=HERE, check=False)
    return st


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="20-step readiness check, all 9")
    ap.add_argument("--gpu", type=int, default=None, help="GPU id (default: auto-pick free)")
    ap.add_argument("--modes", nargs="+", default=list(MODES), choices=MODES)
    ap.add_argument("--energies", nargs="+", default=list(ENERGIES), choices=list(ENERGIES))
    a = ap.parse_args()
    gpu = a.gpu if a.gpu is not None else free_gpu()
    results = {}
    for ekey in a.energies:
        for mode in a.modes:
            results[f"{mode}_{ekey}"] = run_one(mode, ekey, gpu, a.smoke)
    print("\n[dispatch] summary:")
    for k, st in results.items():
        print(f"  {k:20s} {'OK' if st == 0 else f'FAIL(exit {st})'}")
    # PRODUCTION auto-build of the study notebook would go here (notebook-making
    # tail), once the gated analysis (check-stopping-power fourier skill) is ready.
    sys.exit(0 if all(s == 0 for s in results.values()) else 1)


if __name__ == "__main__":
    main()
