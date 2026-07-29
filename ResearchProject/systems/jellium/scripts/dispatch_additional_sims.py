#!/usr/bin/env python3
"""dispatch_additional_sims.py — launch 8 simulations in 4 GPU pairs.

Runs pairs sequentially, with 2 runs in parallel per pair (one per GPU).
After each run completes, runs analyse.py and sends a Gmail notification
with the results.

Pairs:
  1: GPU 0 = coronene run_cc_bond (~50 min)
     GPU 1 = WP E=50 standard (~35 min)
  2: GPU 0 = WP E=50 HD (~20 min)
     GPU 1 = Classical E=50 HD (~30 min)
  3: GPU 0 = WP E=200 HD (~15 min)
     GPU 1 = Classical E=200 HD (~15 min)
  4: GPU 0 = WP E=300 HD (~10 min)
     GPU 1 = Classical E=300 HD (~12 min)

Usage:
    python3 dispatch_additional_sims.py
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path("/local/data/public/skcb2/tddft")
JELLIUM_DIR  = PROJECT_ROOT / "ResearchProject" / "systems" / "jellium"
CORONENE_DIR = PROJECT_ROOT / "ResearchProject" / "systems" / "coronene"

ENV_BASE = {
    **os.environ,
    "PATH": str(PROJECT_ROOT / "shared" / "bin") + ":" + os.environ.get("PATH", ""),
    "INQ_SHARE_PATH": str(PROJECT_ROOT / "inq" / "install" / "share"),
    "PSEUDOPOD_SHARE_PATH": str(PROJECT_ROOT / "inq" / "install" / "share" / "pseudopod"),
}

# Ensure inqview is importable for email sending.
sys.path.insert(0, str(PROJECT_ROOT / "inq-stack" / "python"))

RECIPIENT = "chiddukanna@gmail.com"

# ---------------------------------------------------------------------------
# Run definitions
# ---------------------------------------------------------------------------
RUNS: list[dict] = [
    # Pair 1
    {
        "name": "coronene_run_cc_bond",
        "dir": str(CORONENE_DIR / "run_cc_bond"),
        "gpu": 0,
        "pair": 1,
        "type": "coronene",
    },
    {
        "name": "wp_n162_L50_E50_sigma1",
        "dir": str(JELLIUM_DIR / "run_wp_n162_L50_E50_sigma1"),
        "gpu": 1,
        "pair": 1,
        "type": "wp",
    },
    # Pair 2
    {
        "name": "wp_n162_L30_E50_highdens_sigma1",
        "dir": str(JELLIUM_DIR / "run_wp_n162_L30_E50_highdens_sigma1"),
        "gpu": 0,
        "pair": 2,
        "type": "wp",
    },
    {
        "name": "classical_n162_L30_E50_highdens",
        "dir": str(JELLIUM_DIR / "run_classical_n162_L30_E50_highdens"),
        "gpu": 1,
        "pair": 2,
        "type": "classical",
    },
    # Pair 3
    {
        "name": "wp_n162_L30_E200_highdens_sigma1",
        "dir": str(JELLIUM_DIR / "run_wp_n162_L30_E200_highdens_sigma1"),
        "gpu": 0,
        "pair": 3,
        "type": "wp",
    },
    {
        "name": "classical_n162_L30_E200_highdens",
        "dir": str(JELLIUM_DIR / "run_classical_n162_L30_E200_highdens"),
        "gpu": 1,
        "pair": 3,
        "type": "classical",
    },
    # Pair 4
    {
        "name": "wp_n162_L30_E300_highdens_sigma1",
        "dir": str(JELLIUM_DIR / "run_wp_n162_L30_E300_highdens_sigma1"),
        "gpu": 0,
        "pair": 4,
        "type": "wp",
    },
    {
        "name": "classical_n162_L30_E300_highdens",
        "dir": str(JELLIUM_DIR / "run_classical_n162_L30_E300_highdens"),
        "gpu": 1,
        "pair": 4,
        "type": "classical",
    },
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def read_wall_time(run_dir: str) -> str:
    """Extract wall_time_s from run_summary.txt, return as string."""
    summary = Path(run_dir) / "results" / "run_summary.txt"
    if not summary.exists():
        return "unknown"
    m = re.search(r"^\s*wall_time_s\s*=\s*(\S+)",
                  summary.read_text(), flags=re.MULTILINE)
    return m.group(1) if m else "unknown"


def send_notification(run: dict, status: str, wall_time: str,
                      extra_info: str = "") -> None:
    """Send a Gmail notification for a completed run."""
    try:
        from inqview.email import send_run_email

        subject = f"Simulation complete: {run['name']}"
        body_lines = [
            f"Run:       {run['name']}",
            f"Type:      {run['type']}",
            f"Directory: {run['dir']}",
            f"GPU:       {run['gpu']}",
            f"Pair:      {run['pair']}",
            f"Status:    {status}",
            f"Wall time: {wall_time} s",
        ]
        if extra_info:
            body_lines.append(f"\n{extra_info}")
        body = "\n".join(body_lines)

        attachments: list[str] = []
        obs_png = Path(run["dir"]) / "results" / "analysis" / "observables" / "observables_summary.png"
        if obs_png.exists():
            attachments.append(str(obs_png))

        send_run_email(
            subject=subject,
            body=body,
            attachments=attachments if attachments else None,
            to=RECIPIENT,
        )
        print(f"    [email] Notification sent for {run['name']}")
    except Exception as exc:
        print(f"    [email] WARNING: failed to send email for {run['name']}: {exc}")


def execute_run(run: dict) -> dict:
    """Execute inq-run + analyse.py for a single run. Returns result dict."""
    name = run["name"]
    run_dir = run["dir"]
    gpu_id = run["gpu"]

    print(f"  [{name}] Starting on GPU {gpu_id} ...")
    t0 = time.time()

    env = {**ENV_BASE, "CUDA_VISIBLE_DEVICES": str(gpu_id)}

    # Step 1: inq-run
    inq_result = subprocess.run(
        ["inq-run"],
        cwd=run_dir,
        env=env,
        capture_output=True,
        text=True,
    )

    elapsed_inq = time.time() - t0
    inq_ok = inq_result.returncode == 0

    if inq_ok:
        print(f"  [{name}] inq-run completed in {elapsed_inq:.0f}s (rc=0)")
    else:
        print(f"  [{name}] inq-run FAILED (rc={inq_result.returncode}) "
              f"after {elapsed_inq:.0f}s")
        stderr_tail = inq_result.stderr[-500:] if inq_result.stderr else ""
        send_notification(run, "FAILED (inq-run)", f"{elapsed_inq:.0f}",
                          extra_info=f"stderr tail:\n{stderr_tail}")
        return {"name": name, "inq_ok": False, "analyse_ok": False}

    # Step 2: analyse.py
    analyse_script = Path(run_dir) / "analyse.py"
    if not analyse_script.exists():
        print(f"  [{name}] WARNING: analyse.py not found, skipping analysis")
        wall_time = read_wall_time(run_dir)
        send_notification(run, "SUCCESS (no analyse.py)", wall_time)
        return {"name": name, "inq_ok": True, "analyse_ok": False}

    venv_python = str(PROJECT_ROOT / "venv" / "bin" / "python3")
    analyse_result = subprocess.run(
        [venv_python, "analyse.py"],
        cwd=run_dir,
        env=env,
        capture_output=True,
        text=True,
    )

    elapsed_total = time.time() - t0
    analyse_ok = analyse_result.returncode == 0

    if analyse_ok:
        print(f"  [{name}] analyse.py completed (total {elapsed_total:.0f}s)")
    else:
        print(f"  [{name}] analyse.py FAILED (rc={analyse_result.returncode})")
        stderr_tail = analyse_result.stderr[-500:] if analyse_result.stderr else ""
        print(f"    stderr: {stderr_tail}")

    wall_time = read_wall_time(run_dir)
    status = "SUCCESS" if analyse_ok else "SUCCESS (analyse failed)"
    send_notification(run, status, wall_time)

    return {"name": name, "inq_ok": True, "analyse_ok": analyse_ok}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("=" * 72)
    print("dispatch_additional_sims.py")
    print(f"  {len(RUNS)} runs in 4 pairs")
    print("=" * 72)

    all_results: list[dict] = []

    # Group runs by pair.
    pairs: dict[int, list[dict]] = {}
    for run in RUNS:
        pairs.setdefault(run["pair"], []).append(run)

    for pair_id in sorted(pairs.keys()):
        pair_runs = pairs[pair_id]
        names = [r["name"] for r in pair_runs]
        print(f"\n{'='*72}")
        print(f"Pair {pair_id}: {', '.join(names)}")
        print(f"{'='*72}")

        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = {executor.submit(execute_run, r): r for r in pair_runs}
            for future in as_completed(futures):
                run = futures[future]
                try:
                    result = future.result()
                    all_results.append(result)
                except Exception as exc:
                    print(f"  [{run['name']}] EXCEPTION: {exc}")
                    all_results.append({
                        "name": run["name"],
                        "inq_ok": False,
                        "analyse_ok": False,
                    })

    # Final summary.
    print(f"\n{'='*72}")
    print("DISPATCH SUMMARY")
    print(f"{'='*72}")
    n_inq_ok = sum(1 for r in all_results if r["inq_ok"])
    n_analyse_ok = sum(1 for r in all_results if r["analyse_ok"])
    print(f"  inq-run:    {n_inq_ok}/{len(all_results)} succeeded")
    print(f"  analyse.py: {n_analyse_ok}/{len(all_results)} succeeded")
    for r in all_results:
        inq_flag = "OK" if r["inq_ok"] else "FAIL"
        ana_flag = "OK" if r["analyse_ok"] else "FAIL"
        print(f"    {r['name']:<50} inq={inq_flag}  analyse={ana_flag}")

    failed = [r for r in all_results if not r["inq_ok"]]
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
