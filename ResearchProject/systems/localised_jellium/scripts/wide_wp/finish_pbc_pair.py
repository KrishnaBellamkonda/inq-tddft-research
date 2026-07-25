#!/usr/bin/env python3
"""
Waiter for the PBC WP+classical pair (2026-07-06). The original orchestrator
(run_pbc_pair.py) was retired after the classical run was relaunched with a corrected
length (1500 steps, tau~30 a.u. — the 3540-step version was ~59 h, over the guard).
This WATCHES the two already-running detached processes (never relaunches) and calls
run_pbc_pair.finish() once BOTH write run_completed=true, so the pair is never orphaned.

START (detached):
  cd .../scripts/wide_wp
  setsid nohup /local/data/public/skcb2/tddft/venv/bin/python3 finish_pbc_pair.py \
      > finish_pbc_pair.log 2>&1 &
"""
from __future__ import annotations
import time
from datetime import datetime
from pathlib import Path
import run_pbc_pair as O

WP_RES = O.WPDIR / "results" / O.WP_OUT
CL_RES = O.CLDIR / "results" / O.CL_OUT
WP_LOG = O.WPDIR / "wp_pbc.log"
CL_LOG = O.CLDIR / "classical_openz.log"
DEADLINE_H = 30
POLL_S = 120
STALL_MIN = 30          # log silent this long (no sentinel) => presumed dead


def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)


def stalled(logpath: Path) -> bool:
    """True if the run log has not been written for > STALL_MIN (proc likely dead/hung).
    Immune to MPI/pid quirks — the compute process names/pids vary, but the log mtime
    advances every step (~30-60 s)."""
    if not logpath.exists():
        return False       # not yet created — give it time, not 'dead'
    return (time.time() - logpath.stat().st_mtime) / 60 > STALL_MIN


def main():
    t_end = time.time() + DEADLINE_H * 3600
    log(f"WAIT for PBC pair. WP={WP_RES}  CL={CL_RES}  (liveness = log-staleness > {STALL_MIN}min)")
    while True:
        wp_done, cl_done = O.sentinel(WP_RES), O.sentinel(CL_RES)
        if wp_done and cl_done:
            log("BOTH complete — building notebooks + comparison."); O.finish(WP_RES, CL_RES); return
        if time.time() > t_end:
            log("DEADLINE exceeded.")
            O.email("[wide-wp PBC] waiter TIMEOUT",
                     f"Deadline {DEADLINE_H} h hit. WP done={wp_done}, CL done={cl_done}. "
                     "Building whatever finished.")
            O.finish(WP_RES if wp_done else None, CL_RES if cl_done else None); return
        wp_dead = (not wp_done) and stalled(WP_LOG)
        cl_dead = (not cl_done) and stalled(CL_LOG)
        if wp_dead or cl_dead:
            who = "WP" if wp_dead else "CLASSICAL"
            log(f"{who}: log stalled > {STALL_MIN} min, no sentinel — presumed dead.")
            O.email(f"[wide-wp PBC] {who} STALLED/DIED",
                     f"{who} log silent > {STALL_MIN} min with no completion sentinel. "
                     f"WP done={wp_done}, CL done={cl_done}. Building whatever finished.")
            O.finish(WP_RES if wp_done else None, CL_RES if cl_done else None); return
        time.sleep(POLL_S)


if __name__ == "__main__":
    main()
