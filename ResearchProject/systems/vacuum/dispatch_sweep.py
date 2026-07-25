#!/usr/bin/env python3
"""Dispatch the MFA reflection ε(E,L) sweep across 2 GPUs.

One prebuilt binary (mfa_sweep/run) is invoked per (k0, L) job with env vars.
GPU availability is checked with the cudaMemGetInfo probe (NVML/nvidia-smi is
broken on this box). Gated externally: only run after gate-1 + gate-2 pass.

Grid (plan §6/D6): L ∈ {5,10,20,30,40,50} × 12 log-spaced k0 (E≈0.5–490 eV)
= 72 masked runs + 4 no-absorber anchors (ε≈1). ~6 showcase runs write density
frames + final wavefunction.

Usage:
  python3 dispatch_sweep.py --dry-run      # print job list, no runs
  python3 dispatch_sweep.py                 # run the full sweep on free GPUs
  python3 dispatch_sweep.py --jobs k0=1.5,L=20         # run a single ad-hoc job
"""
import argparse
import math
import os
import queue
import subprocess
import sys
import threading
import time
from pathlib import Path

VAC = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/vacuum")
BIN = VAC / "mfa_sweep" / "run"
RUNS = VAC / "runs"
PROBE = VAC / "gpu_probe"
HA_TO_EV = 27.211386245988

ENV_BASE = {
    **os.environ,
    "PATH": "/local/data/public/skcb2/tddft/shared/bin:" + os.environ.get("PATH", ""),
    "INQ_SHARE_PATH": "/local/data/public/skcb2/tddft/inq/install/share",
    "PSEUDOPOD_SHARE_PATH": "/local/data/public/skcb2/tddft/inq/install/share/pseudopod",
}

L_VALUES = [5, 10, 20, 30, 40, 50]
N_K0 = 12
K0_MIN, K0_MAX = 0.2, 6.0
K0_VALUES = [10 ** (math.log10(K0_MIN) + (math.log10(K0_MAX) - math.log10(K0_MIN)) * i / (N_K0 - 1))
             for i in range(N_K0)]
# Anchors dropped: INQ finite-cell hard-wall reflection is messy (ε≈0.43, not 1),
# and the paper's ε itself plateaus below 1 at low E — so a clean ε≈1 anchor is
# neither achievable nor the right reference. The low-E masked runs supply the
# high-ε end, and the independent cap_toy 1D reference validates the pipeline.
ANCHOR_K0 = []
ANCHOR_L = 20
# showcase (k0, L): low-E reflecting + high-E absorbed + mid
SHOWCASE = {(round(K0_VALUES[1], 4), 10), (round(K0_VALUES[1], 4), 50),
            (round(K0_VALUES[9], 4), 10), (round(K0_VALUES[9], 4), 50),
            (round(K0_VALUES[5], 4), 20), (round(K0_VALUES[7], 4), 30)}


def e_ev(k0):
    return 0.5 * k0 * k0 * HA_TO_EV


def free_gpus(min_free_mb=4000):
    if not PROBE.exists():
        print(f"WARN: probe {PROBE} missing; assuming GPUs 0,1", file=sys.stderr)
        return [0, 1]
    out = subprocess.run([str(PROBE)], capture_output=True, text=True).stdout
    gpus = []
    for line in out.strip().splitlines():
        parts = line.split()
        if len(parts) >= 4 and parts[0] == "GPU":
            idx, free_mb = int(parts[1]), int(parts[3])
            if free_mb >= min_free_mb:
                gpus.append(idx)
    return gpus


def build_jobs():
    jobs = []
    for k0 in K0_VALUES:
        for L in L_VALUES:
            show = (round(k0, 4), L) in SHOWCASE
            jobs.append(dict(k0=k0, L=L, anchor=0, showcase=1 if show else 0))
    for k0 in ANCHOR_K0:
        jobs.append(dict(k0=k0, L=ANCHOR_L, anchor=1, showcase=0))
    return jobs


def job_name(j):
    tag = "anchor" if j["anchor"] else "mfa"
    s = "_show" if j["showcase"] else ""
    return f"run_{tag}_E{e_ev(j['k0']):07.2f}_L{int(j['L'])}{s}"


def run_job(j, gpu, log):
    name = job_name(j)
    outdir = RUNS / name
    outdir.mkdir(parents=True, exist_ok=True)
    env = dict(ENV_BASE)
    env.update(CUDA_VISIBLE_DEVICES=str(gpu), MFA_K0=f"{j['k0']:.10g}",
               MFA_LABS=f"{j['L']:.10g}", MFA_ANCHOR=str(j["anchor"]),
               MFA_SHOWCASE=str(j["showcase"]), MFA_OUTDIR=str(outdir),
               MFA_NPERP="8")  # 8 (not 12): 2.25× faster, ε transverse-insensitive (proven)
    t0 = time.time()
    with open(outdir / "run.log", "w") as lf:
        rc = subprocess.run([str(BIN)], cwd=str(outdir), env=env,
                            stdout=lf, stderr=subprocess.STDOUT).returncode
    dt = time.time() - t0
    eps = read_eps(outdir / "epsilon.txt")
    status = "OK" if (rc == 0 and eps is not None) else "FAIL"
    line = f"[GPU{gpu}] {status:4s} {name}  E={e_ev(j['k0']):.2f}eV L={j['L']} eps={eps} ({dt:.0f}s rc={rc})"
    print(line, flush=True)
    log.append(dict(name=name, k0=j["k0"], E_eV=e_ev(j["k0"]), L=j["L"],
                    anchor=j["anchor"], showcase=j["showcase"], eps=eps,
                    rc=rc, seconds=dt, status=status))


def read_eps(path):
    try:
        for line in Path(path).read_text().splitlines():
            if line.startswith("epsilon "):
                return float(line.split()[1])
    except Exception:
        return None
    return None


_lock = threading.Lock()


def write_grid(log):
    import csv
    rows = sorted(log, key=lambda r: (r["anchor"], r["L"], r["k0"]))
    with open(RUNS / "epsilon_grid.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "k0", "E_eV", "L", "anchor",
                                          "showcase", "eps", "rc", "seconds", "status"])
        w.writeheader()
        w.writerows(rows)


def worker(gpu, q, log):
    while True:
        try:
            j = q.get_nowait()
        except queue.Empty:
            return
        try:
            run_job(j, gpu, log)
            with _lock:
                write_grid(log)   # incremental: partial results usable if interrupted
        except Exception as e:
            print(f"[GPU{gpu}] EXC {job_name(j)}: {e}", flush=True)
        finally:
            q.task_done()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--jobs", default=None, help="ad-hoc single job 'k0=1.5,L=20[,anchor=1][,showcase=1]'")
    args = ap.parse_args()

    if args.jobs:
        kv = dict(p.split("=") for p in args.jobs.split(","))
        jobs = [dict(k0=float(kv["k0"]), L=float(kv["L"]),
                     anchor=int(kv.get("anchor", 0)), showcase=int(kv.get("showcase", 0)))]
    else:
        jobs = build_jobs()

    print(f"jobs={len(jobs)}  k0={[round(k,3) for k in K0_VALUES]}")
    print(f"E_eV range = {e_ev(K0_VALUES[0]):.2f} .. {e_ev(K0_VALUES[-1]):.2f}")
    if args.dry_run:
        for j in jobs:
            print(f"  {job_name(j)}  k0={j['k0']:.4f} E={e_ev(j['k0']):.2f} L={j['L']} "
                  f"anchor={j['anchor']} showcase={j['showcase']}")
        return

    gpus = free_gpus()
    if not gpus:
        print("ERROR: no free GPU found by probe. Aborting.", file=sys.stderr)
        sys.exit(2)
    print(f"using GPUs {gpus}")

    RUNS.mkdir(parents=True, exist_ok=True)
    q = queue.Queue()
    # CHEAPEST (high k0) first: the scientifically valuable rolloff (E≈3–490 eV)
    # completes early; the expensive low-E plateau runs last. If the sweep
    # overruns, the curve is mostly complete (CSV is written incrementally).
    for j in sorted(jobs, key=lambda x: -x["k0"]):
        q.put(j)
    log = []
    threads = [threading.Thread(target=worker, args=(g, q, log)) for g in gpus]
    t0 = time.time()
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    # aggregate
    import csv
    log.sort(key=lambda r: (r["anchor"], r["L"], r["k0"]))
    with open(RUNS / "epsilon_grid.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "k0", "E_eV", "L", "anchor",
                                          "showcase", "eps", "rc", "seconds", "status"])
        w.writeheader()
        w.writerows(log)
    ok = sum(1 for r in log if r["status"] == "OK")
    print(f"\nDONE: {ok}/{len(log)} OK in {(time.time()-t0)/60:.1f} min. "
          f"Grid → {RUNS/'epsilon_grid.csv'}")


if __name__ == "__main__":
    main()
