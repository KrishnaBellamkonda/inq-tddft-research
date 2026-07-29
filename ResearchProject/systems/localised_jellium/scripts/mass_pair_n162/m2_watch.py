#!/usr/bin/env python3
"""Resume-safety watcher for the manually-relaunched m2 run.

The orchestrator (orchestrate.py) monitors ONLY m1; m2 crashed at step 0 (no
checkpoint) and was dropped from its monitor loop, then relaunched by hand on
GPU1. This watcher gives m2 the SAME auto-resume the orchestrator gives m1:
poll the run; if it dies before 2500 steps, relaunch from its last checkpoint
(EM_RESUME=1), or fresh if it died before the first (step-500) checkpoint.

Adopts the already-running m2 via results/m2/m2.pid. Fully detached children
(start_new_session) so they outlive this watcher. Emails on resume / done / fail.
"""
import os, subprocess, sys, time
from pathlib import Path
from datetime import datetime

ROOT   = Path("/local/data/public/skcb2/tddft")
LJ     = ROOT/"ResearchProject/systems/localised_jellium"
WP     = LJ/"scripts/mass_pair_n162/wp"
GSDIR  = LJ/"shared_gs/slab_n120_L60x60x62_dx0p40"
PIDF   = WP/"results/m2/m2.pid"
RTSTAY = WP/"results/m2/rt_ckpt/rt_state.txt"
NSTEPS = 2500
MAXTRIES = 4
TO     = "chiddukanna@gmail.com"

ENV = {**os.environ,
       "PATH":                 f"{ROOT/'shared/bin'}:" + os.environ.get("PATH",""),
       "INQ_SOURCE":           str(ROOT/"inq-study"),
       "INQ_SHARE_PATH":       str(ROOT/"inq/install/share"),
       "PSEUDOPOD_SHARE_PATH": str(ROOT/"inq/install/share/pseudopod")}
RUNENV = {"CUDA_VISIBLE_DEVICES":"1", "EM_OUT":"m2",
          "EM_DT":"0.04", "EM_N_STEPS":str(NSTEPS), "EM_WRITE_EVERY":"8",
          "EM_WF_EVERY":"40", "EM_CKPT_EVERY":"500", "EM_CAP":"1", "EM_CAP_ETA":"-1.0",
          "EM_CAP_CENTER_BOHR":"26.2", "EM_CAP_WIDTH_BOHR":"10.0",
          "EM_SIGMA_WP":"1.0", "EM_LAUNCH_Z":"-16.5", "EM_GS_DIR":str(GSDIR),
          "EM_OBS_DIPCUR":"0", "EM_OBS_MOM":"0", "EM_OBS_OVL":"0",
          "EM_OBS_WF":"0", "EM_OBS_STATE":"0",
          "EM_INV_MASS":"0.5", "EM_K0":"3.834"}

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)

def email(subj, body):
    try:
        sys.path.insert(0, str(ROOT/"inq-stack/python"))
        from inqview.email import send_run_email
        send_run_email(subject=f"[mass_pair_n162] {subj}", body=body, to=TO)
        log(f"  emailed: {subj}")
    except Exception as e:
        log(f"  email FAILED ({subj}): {e}")

def last_step():
    try:
        for ln in RTSTAY.read_text().splitlines():
            if ln.startswith("last_step="):
                return int(ln.split("=",1)[1])
    except Exception:
        pass
    return 0

def alive(pid):
    if pid <= 0: return False
    try:
        os.kill(pid, 0); return True
    except Exception:
        return False

def cur_pid():
    try:
        return int(PIDF.read_text().strip())
    except Exception:
        return 0

def launch(resume):
    env = {**ENV, **RUNENV, "EM_RESUME": "1" if resume else "0"}
    f = open(WP/"rt_m2.log", "a")
    p = subprocess.Popen([str(WP/"run")], cwd=str(WP), env=env,
                         stdout=f, stderr=subprocess.STDOUT,
                         stdin=subprocess.DEVNULL, start_new_session=True)
    PIDF.write_text(str(p.pid))
    log(f"launched m2 pid={p.pid} resume={resume}")
    return p.pid

def main():
    tries = 0
    log(f"m2 watcher start — adopting pid={cur_pid()}, target {NSTEPS} steps, max {MAXTRIES} resumes")
    # grace period so an in-flight launch can allocate before first liveness check
    time.sleep(120)
    while True:
        step = last_step()
        if step >= NSTEPS:
            email("m2 complete", f"m2 (projectile mass 2) reached {step}/{NSTEPS} steps. "
                                 f"Results: {WP}/results/m2/")
            log(f"m2 DONE at step {step}"); break
        pid = cur_pid()
        if alive(pid):
            time.sleep(120); continue
        # not alive and not done
        if tries >= MAXTRIES:
            email("m2 RUN FAILED (watcher)", f"m2 died and exhausted {tries} resume attempts; "
                  f"last checkpoint step={step}. See {WP}/rt_m2.log. Resume manually with EM_RESUME=1.")
            log(f"m2 FAILED permanently, tries={tries}, last_step={step}"); break
        tries += 1
        log(f"m2 not alive; resuming from step {step} (attempt {tries}/{MAXTRIES})")
        email("m2 auto-resume", f"m2 exited early; resuming from checkpoint step {step} "
              f"(attempt {tries}/{MAXTRIES}).")
        launch(resume=step > 0)
        time.sleep(120)

if __name__ == "__main__":
    main()
