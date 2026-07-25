#!/usr/bin/env python3
"""
Standalone WAITER for the already-running long open-z WP run (2026-07-04).

Why this exists: the parent orchestrator (run_long_wp_per2.py) had MAX_HOURS=14,
but the realised rate (~38 s/step under open-z 2D-Poisson z-doubling + VTI writes)
makes the 1773-step run ~18.7 h. To avoid a premature kill, the orchestrator is
retired and THIS waiter takes over: it watches the EXISTING WP process (it does NOT
relaunch it), and on genuine completion builds the run-notebook + emails the 4-part
result via the reused run_long_wp_per2.finish().

START (detached):
  cd .../scripts/wide_wp
  setsid nohup /local/data/public/skcb2/tddft/venv/bin/python3 finish_wp.py \
      > finish_wp.log 2>&1 &
"""
from __future__ import annotations
import sys, time, subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import run_long_wp_per2 as O   # reuse log/email/sentinel/finish/constants

WP_RES = O.WPDIR / "results" / O.WP_OUT   # run.cpp prepends "results/" to LJ_OUT
WP_LOG = O.WPDIR / "wp_per2_long.log"
DEADLINE_H = 24
POLL_S = 120
STALL_MIN = 30


def wp_alive() -> bool:
    try:
        subprocess.check_output(["pgrep", "-f", "wide_wp/wp/run"])
        return True
    except subprocess.CalledProcessError:
        return False


def main():
    t_deadline = time.time() + DEADLINE_H * 3600
    O.log(f"WAITER: watching existing WP (deadline {DEADLINE_H} h). NO relaunch.")
    while True:
        if O.sentinel(WP_RES):
            O.log("WAITER: WP sentinel found — running finish().")
            O.finish(WP_RES)
            O.log("WAITER: done.")
            return
        if time.time() > t_deadline:
            O.email("[wide-wavepacket] open-z WP — WAITER TIMEOUT",
                    f"WP did not complete within {DEADLINE_H} h. See {WP_LOG}.\n"
                    f"Partial data in {WP_RES}. Nothing auto-built.")
            O.log("WAITER: deadline exceeded — abort."); return
        silent_min = (time.time() - WP_LOG.stat().st_mtime) / 60 if WP_LOG.exists() else 999
        if not wp_alive() and not O.sentinel(WP_RES):
            O.email("[wide-wavepacket] open-z WP — PROCESS DIED before completion",
                    f"The WP process is gone and no completion sentinel was written.\n"
                    f"Last log activity {silent_min:.0f} min ago. See {WP_LOG}.\n"
                    f"Partial data in {WP_RES}.")
            O.log("WAITER: WP process died without sentinel — alerted, abort."); return
        if silent_min > STALL_MIN:
            O.log(f"WAITER: WP log silent {silent_min:.0f} min but process still listed — keep waiting.")
        time.sleep(POLL_S)


if __name__ == "__main__":
    main()
