#!/usr/bin/env python3
"""dispatch_runs.py — multi-GPU dispatcher for coronene runs.

Reads a queue file (one run-directory absolute path per line, ``#`` for
comments). For each GPU in ``--gpus`` it polls ``nvidia-smi`` for activity
from any process; when a GPU is free it pops the next queue entry and
launches::

    cd <run_dir> && CUDA_VISIBLE_DEVICES=<gpu_id> inq-run run.cpp \\
        > run.log 2>&1

It then waits, re-polling every ``--poll-seconds`` seconds, refilling any
freed GPU until the queue is empty *and* every dispatched child has
finished. **It will not preempt other users' jobs**; if both GPUs are busy
on entry, it sleep-polls until at least one frees up.

Tier-A GPU-execution check (per the plan §5):

* ``inq-run`` (no ``--cpu`` flag) is used.
* The dispatcher records the child PID and the assigned GPU. After the
  child exits, the run's log is grepped for the INQ GPU initialisation
  banner; failure is reported but does not abort the queue (so an
  unrelated bug doesn't stop the rest).

CLI
---

::

    dispatch_runs.py run_queue.txt
                     [--poll-seconds 30]
                     [--gpus 0,1]
                     [--free-mem-mb 200]
                     [--dry-run]
                     [--log-file dispatch.log]

The ``--dry-run`` mode prints the launch command for each entry without
running anything; used for the user's review pass.
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

POLL_DEFAULT_S = 30.0
FREE_MEM_DEFAULT_MB = 200
GPU_BUSY_CONFIRM_POLLS = 1   # confirm idleness only once for our own children
GPU_FREE_CONFIRM_POLLS = 2   # require 2 consecutive idle reads for new launches


@dataclass
class GPUSlot:
    gpu_id: int
    child: subprocess.Popen | None = None
    run_dir: Path | None = None
    launched_at: float = 0.0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="dispatch_runs",
                                description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("queue_file", type=Path,
                   help="newline-delimited run-dir paths (# = comment)")
    p.add_argument("--poll-seconds", type=float, default=POLL_DEFAULT_S,
                   help=f"poll interval (default {POLL_DEFAULT_S})")
    p.add_argument("--gpus", type=str, default="0,1",
                   help="comma-separated GPU ids to use (default 0,1)")
    p.add_argument("--free-mem-mb", type=int, default=FREE_MEM_DEFAULT_MB,
                   help="GPU mem usage below this (MiB) counts as free "
                        f"(default {FREE_MEM_DEFAULT_MB})")
    p.add_argument("--dry-run", action="store_true",
                   help="print launch commands without executing")
    p.add_argument("--log-file", type=Path, default=Path("dispatch.log"),
                   help="path for the rolling dispatcher log")
    p.add_argument("--inq-run", type=str, default="inq-run",
                   help="inq-run executable (default: 'inq-run' on PATH)")
    p.add_argument("--cpp-file", type=str, default="run.cpp",
                   help="run.cpp basename (default: run.cpp)")
    p.add_argument("--clear-results", action="store_true",
                   help="wipe each run's results/{raw,analysis,run_summary.txt}"
                        " before launching, so a re-run starts clean.")
    return p.parse_args(argv)


def read_queue(path: Path) -> list[Path]:
    if not path.exists():
        raise FileNotFoundError(path)
    out: list[Path] = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        out.append(Path(s).resolve())
    return out


def gpu_is_free(gpu_id: int, *, free_mem_mb: int, owned_pid: int | None) -> bool:
    """A GPU is free if no compute apps run on it (other than our own
    ``owned_pid``) and used memory < free_mem_mb."""
    try:
        # Compute apps
        out = subprocess.check_output(
            ["nvidia-smi",
             f"--query-compute-apps=pid,used_memory",
             "--format=csv,noheader,nounits",
             f"-i", str(gpu_id)],
            text=True, stderr=subprocess.DEVNULL)
        for line in out.splitlines():
            parts = [s.strip() for s in line.split(",")]
            if not parts or not parts[0]:
                continue
            try:
                pid = int(parts[0])
            except ValueError:
                continue
            if owned_pid is not None and pid == owned_pid:
                continue
            return False  # someone else's job
        # Memory used (total fb)
        out2 = subprocess.check_output(
            ["nvidia-smi",
             f"--query-gpu=memory.used",
             "--format=csv,noheader,nounits",
             f"-i", str(gpu_id)],
            text=True, stderr=subprocess.DEVNULL)
        used_mib = int(out2.strip().splitlines()[0])
        if owned_pid is None and used_mib > free_mem_mb:
            return False
        return True
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        print(f"[warn] nvidia-smi failed for gpu {gpu_id}: {e}",
              file=sys.stderr)
        return False


def gpu_free_with_confirmation(gpu_id: int, *, free_mem_mb: int,
                                poll_seconds: float, n: int) -> bool:
    """Read free-state ``n`` times, ``poll_seconds / n`` apart. Conservative
    against transient one-shot probes flickering on then off."""
    interval = max(0.5, poll_seconds / max(n, 1))
    for _ in range(n):
        if not gpu_is_free(gpu_id, free_mem_mb=free_mem_mb, owned_pid=None):
            return False
        time.sleep(interval)
    return True


class DispatchLog:
    def __init__(self, path: Path):
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = path.open("a")
    def __call__(self, msg: str) -> None:
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        self.fh.write(line + "\n")
        self.fh.flush()
    def close(self):
        self.fh.close()


def _clear_results(run_dir: Path, log: DispatchLog) -> None:
    """Wipe the run's output trees so a re-run starts clean.

    Removes ``results/raw`` and ``results/analysis`` (the populated
    subtrees) but leaves ``results/`` itself intact. Build artefacts
    (build/, run binary, build_run.log, run.log) are left so
    ``inq-run`` can incrementally rebuild — only old simulation output
    is purged.
    """
    for sub in ("raw", "analysis"):
        target = run_dir / "results" / sub
        if target.exists():
            shutil.rmtree(target)
    rs = run_dir / "results" / "run_summary.txt"
    if rs.exists():
        rs.unlink()
    log(f"[gpu --] cleared {run_dir.name}/results/{{raw,analysis,run_summary.txt}}")


def launch(slot: GPUSlot, run_dir: Path, *, inq_run: str, cpp_file: str,
           dry_run: bool, log: DispatchLog,
           clear_results: bool = False) -> None:
    cmd = ["env", f"CUDA_VISIBLE_DEVICES={slot.gpu_id}", inq_run, cpp_file]
    log_str = f"cd {run_dir} && CUDA_VISIBLE_DEVICES={slot.gpu_id} {inq_run} {cpp_file}"
    if dry_run:
        if clear_results:
            log(f"[dry-run] [gpu {slot.gpu_id}] would clear "
                f"{run_dir.name}/results/{{raw,analysis,run_summary.txt}}")
        log(f"[dry-run] [gpu {slot.gpu_id}] {log_str}")
        return
    if clear_results:
        _clear_results(run_dir, log)
    run_log = run_dir / "run.log"
    fh = run_log.open("ab")
    log(f"[gpu {slot.gpu_id}] launching {run_dir.name} (log: {run_log})")
    p = subprocess.Popen(
        cmd, cwd=str(run_dir),
        stdout=fh, stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    slot.child = p
    slot.run_dir = run_dir
    slot.launched_at = time.monotonic()


def reap(slot: GPUSlot, log: DispatchLog) -> bool:
    if slot.child is None:
        return False
    rc = slot.child.poll()
    if rc is None:
        return False
    walltime = time.monotonic() - slot.launched_at
    log(f"[gpu {slot.gpu_id}] finished {slot.run_dir.name if slot.run_dir else '?'} "
        f"exit={rc} walltime={walltime:.1f}s")
    if slot.run_dir is not None and rc == 0:
        run_log = slot.run_dir / "run.log"
        if run_log.exists():
            txt = run_log.read_text(errors="replace")
            if "GPU" not in txt and "CUDA" not in txt:
                log(f"[warn] [gpu {slot.gpu_id}] no GPU/CUDA banner in "
                    f"{run_log}; tier-A GPU check FAILED")
    slot.child = None
    slot.run_dir = None
    return True


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    queue = read_queue(args.queue_file)
    gpus = [int(x) for x in args.gpus.split(",") if x.strip()]
    log = DispatchLog(args.log_file)

    if args.dry_run:
        log(f"DRY RUN: queue has {len(queue)} entries, gpus={gpus}")
    else:
        if shutil.which(args.inq_run) is None and not Path(args.inq_run).exists():
            log(f"[fatal] inq-run not found on PATH ({args.inq_run!r}). "
                "source ~/.bashrc first.")
            return 2
        log(f"queue={len(queue)} entries, gpus={gpus}, "
            f"poll={args.poll_seconds}s, free_mem_threshold={args.free_mem_mb}MiB")

    slots = [GPUSlot(g) for g in gpus]
    pending = list(queue)
    stop = False

    def _sigint(_signum, _frame):
        nonlocal stop
        log("Ctrl-C received; will not start new jobs. "
            "Existing children will finish naturally.")
        stop = True
    signal.signal(signal.SIGINT, _sigint)

    try:
        while pending or any(s.child is not None for s in slots):
            # Reap any finished children first
            for s in slots:
                reap(s, log)

            # Try to launch on idle slots
            if not stop:
                for s in slots:
                    if s.child is not None or not pending:
                        continue
                    if args.dry_run:
                        run_dir = pending.pop(0)
                        launch(s, run_dir, inq_run=args.inq_run,
                               cpp_file=args.cpp_file, dry_run=True, log=log,
                               clear_results=args.clear_results)
                        continue
                    free = gpu_free_with_confirmation(
                        s.gpu_id,
                        free_mem_mb=args.free_mem_mb,
                        poll_seconds=args.poll_seconds,
                        n=GPU_FREE_CONFIRM_POLLS,
                    )
                    if not free:
                        log(f"[gpu {s.gpu_id}] busy (someone else); waiting")
                        continue
                    run_dir = pending.pop(0)
                    launch(s, run_dir, inq_run=args.inq_run,
                           cpp_file=args.cpp_file, dry_run=False, log=log,
                           clear_results=args.clear_results)

            if args.dry_run:
                # nothing actually running; exit when queue is empty
                if not pending:
                    break
                continue

            time.sleep(args.poll_seconds)

    finally:
        log.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
