#!/usr/bin/env python3
"""dispatch_v2_runs.py — launch 15 jellium v2 re-runs with queue-based GPU scheduling.

Queue-based: polls nvidia-smi for free GPUs, launches next run when a GPU
becomes available. Max 2 concurrent runs (one per GPU).

After each run: analyse.py with venv Python + Gmail notification.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

PROJECT_ROOT = Path("/local/data/public/skcb2/tddft")
JELLIUM_DIR  = PROJECT_ROOT / "ResearchProject" / "systems" / "jellium"
VENV_PYTHON  = str(PROJECT_ROOT / "venv" / "bin" / "python3")

ENV_BASE = {
    **os.environ,
    "PATH": str(PROJECT_ROOT / "shared" / "bin") + ":" + os.environ.get("PATH", ""),
    "INQ_SHARE_PATH": str(PROJECT_ROOT / "inq" / "install" / "share"),
    "PSEUDOPOD_SHARE_PATH": str(PROJECT_ROOT / "inq" / "install" / "share" / "pseudopod"),
}

sys.path.insert(0, str(PROJECT_ROOT / "inq-stack" / "python"))
RECIPIENT = "chiddukanna@gmail.com"

# 15 planned v2 runs (pilot E=100 σ=1 already completed)
RUNS = [
    # WP σ=1 standard density (L=50) — 5 runs
    "run_wp_n162_L50_E20_sigma1_v2",
    "run_wp_n162_L50_E25_sigma1_v2",
    "run_wp_n162_L50_E50_sigma1_v2",
    "run_wp_n162_L50_E200_sigma1_v2",
    "run_wp_n162_L50_E300_sigma1_v2",
    # WP σ=1 high-density (L=30) — 4 runs
    "run_wp_n162_L30_E50_highdens_sigma1_v2",
    "run_wp_n162_L30_E100_highdens_sigma1_v2",
    "run_wp_n162_L30_E200_highdens_sigma1_v2",
    "run_wp_n162_L30_E300_highdens_sigma1_v2",
    # Classical standard (L=50) — 3 runs
    "run_classical_n162_L50_E50_v2",
    "run_classical_n162_L50_E100_v2",
    "run_classical_n162_L50_E300_v2",
    # Classical high-density (L=30) — 3 runs
    "run_classical_n162_L30_E50_highdens_v2",
    "run_classical_n162_L30_E100_highdens_v2",
    "run_classical_n162_L30_E300_highdens_v2",
]


def find_free_gpu() -> int | None:
    result = subprocess.run(
        ["nvidia-smi", "--query-gpu=index,memory.used",
         "--format=csv,noheader,nounits"],
        capture_output=True, text=True)
    for line in result.stdout.strip().split("\n"):
        idx, mem = line.split(",")
        if int(mem.strip()) < 500:
            return int(idx.strip())
    return None


def read_wall_time(run_dir: str) -> str:
    summary = Path(run_dir) / "results" / "run_summary.txt"
    if not summary.exists():
        return "unknown"
    m = re.search(r"^\s*wall_time_s\s*=\s*(\S+)",
                  summary.read_text(), flags=re.MULTILINE)
    return m.group(1) if m else "unknown"


def send_notification(name: str, run_dir: str, status: str, wall_time: str):
    try:
        from inqview.email import send_run_email
        body = f"Run: {name}\nDir: {run_dir}\nStatus: {status}\nWall time: {wall_time} s"
        attachments = []
        obs_png = Path(run_dir) / "results" / "analysis" / "observables" / "observables_summary.png"
        if obs_png.exists():
            attachments.append(str(obs_png))
        send_run_email(subject=f"v2 complete: {name}", body=body,
                       attachments=attachments or None, to=RECIPIENT)
        print(f"    [email] sent for {name}")
    except Exception as exc:
        print(f"    [email] FAILED for {name}: {exc}")


def execute_run(name: str, gpu_id: int) -> dict:
    run_dir = str(JELLIUM_DIR / name)
    print(f"  [{name}] Starting on GPU {gpu_id} ...")
    t0 = time.time()

    env = {**ENV_BASE, "CUDA_VISIBLE_DEVICES": str(gpu_id)}

    inq_result = subprocess.run(["inq-run"], cwd=run_dir, env=env,
                                capture_output=True, text=True)
    elapsed = time.time() - t0
    inq_ok = inq_result.returncode == 0

    if not inq_ok:
        print(f"  [{name}] inq-run FAILED (rc={inq_result.returncode}) after {elapsed:.0f}s")
        send_notification(name, run_dir, "FAILED", f"{elapsed:.0f}")
        return {"name": name, "inq_ok": False, "analyse_ok": False}

    print(f"  [{name}] inq-run OK ({elapsed:.0f}s). Running analyse.py ...")

    analyse_result = subprocess.run(
        [VENV_PYTHON, "analyse.py"], cwd=run_dir, env=env,
        capture_output=True, text=True)
    analyse_ok = analyse_result.returncode == 0
    total = time.time() - t0

    wall_time = read_wall_time(run_dir)
    status = "SUCCESS" if analyse_ok else "SUCCESS (analyse partial)"
    send_notification(name, run_dir, status, wall_time)
    print(f"  [{name}] Done ({total:.0f}s total, analyse={'OK' if analyse_ok else 'FAIL'})")

    return {"name": name, "inq_ok": True, "analyse_ok": analyse_ok}


def main() -> int:
    print("=" * 72)
    print(f"dispatch_v2_runs.py — {len(RUNS)} runs, queue-based GPU scheduling")
    print("=" * 72)

    results = []
    queue = list(RUNS)
    active = {}  # gpu_id -> Future

    with ThreadPoolExecutor(max_workers=2) as executor:
        while queue or active:
            # Launch on free GPUs
            while queue:
                gpu = find_free_gpu()
                if gpu is None or gpu in active:
                    break
                name = queue.pop(0)
                future = executor.submit(execute_run, name, gpu)
                active[gpu] = (future, name)
                print(f"\n  Launched {name} on GPU {gpu} ({len(queue)} remaining in queue)")

            # Wait for any completion
            if active:
                done_gpus = []
                for gpu_id, (future, name) in active.items():
                    if future.done():
                        try:
                            result = future.result()
                            results.append(result)
                        except Exception as exc:
                            print(f"  [{name}] EXCEPTION: {exc}")
                            results.append({"name": name, "inq_ok": False, "analyse_ok": False})
                        done_gpus.append(gpu_id)

                for gpu_id in done_gpus:
                    del active[gpu_id]

                if not done_gpus:
                    time.sleep(30)

    # Summary
    print(f"\n{'='*72}")
    print("DISPATCH SUMMARY")
    print(f"{'='*72}")
    n_ok = sum(1 for r in results if r["inq_ok"])
    n_ana = sum(1 for r in results if r["analyse_ok"])
    print(f"  inq-run:    {n_ok}/{len(results)}")
    print(f"  analyse.py: {n_ana}/{len(results)}")
    for r in results:
        print(f"    {r['name']:<50} inq={'OK' if r['inq_ok'] else 'FAIL'}  "
              f"analyse={'OK' if r['analyse_ok'] else 'FAIL'}")

    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
