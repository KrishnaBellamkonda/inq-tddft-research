#!/usr/bin/env python3
"""B2 finisher (energy book-keeping campaign): wait for the in-flight ghost SCFs,
launch the 83-electron illustration run, then rebuild+execute the campaign notebook.
Idempotent: skips any run whose run_summary shows run_completed = true."""
import os, subprocess, time, sys
from datetime import datetime
from pathlib import Path

GH = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/campaign_autorun/gs_ghost")
HYP = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/campaign_autorun_study")
PY = "/local/data/public/skcb2/tddft/venv/bin/python3"
ENV = {**os.environ,
       "INQ_SHARE_PATH": "/local/data/public/skcb2/tddft/inq/install/share",
       "PSEUDOPOD_SHARE_PATH": "/local/data/public/skcb2/tddft/inq/install/share/pseudopod"}

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)

def done(tag):
    rs = GH / "runs" / tag / "results" / "run_summary.txt"
    return rs.exists() and "run_completed = true" in rs.read_text()

def wait_for(tags, timeout_h=8):
    t0 = time.time()
    while time.time() - t0 < timeout_h * 3600:
        if all(done(t) for t in tags):
            return True
        time.sleep(60)
    return False

def run(tag, extra, gpu):
    if done(tag):
        log(f"SKIP {tag} (complete)"); return True
    d = GH / "runs" / tag; d.mkdir(parents=True, exist_ok=True)
    env = {**ENV, "CUDA_VISIBLE_DEVICES": gpu, "LJ_LZ": "120", "LJ_PERIODICITY": "2",
           "LJ_TAG": tag, "LJ_GS_DIR": str(d / "checkpoint"), **extra}
    log(f"RUN {tag} on GPU {gpu}")
    with open(d / "run.log", "a") as lf:
        rc = subprocess.run([str(GH / "run")], cwd=str(d), env=env,
                            stdout=lf, stderr=subprocess.STDOUT).returncode
    log(f"{tag} rc={rc} done={done(tag)}")
    return done(tag)

log("waiting for ghost_r12_p2 + ghost_r4_p2 ...")
ok = wait_for(["ghost_r12_p2", "ghost_r4_p2"])
log(f"in-flight runs complete: {ok}")
run("n83_r12", {"LJ_GHOST": "0", "LJ_N": "83", "LJ_N_BG": "82", "LJ_LAUNCH_Z": "-24.5"}, gpu="1")
log("rebuilding + executing campaign notebook ...")
rc = subprocess.run([PY, str(HYP / "build_energy_book_keeping_notebook.py"), "--execute"],
                    cwd=str(HYP)).returncode
log(f"notebook execute rc={rc}")
(GH / "B2_ALL_DONE.marker").write_text(f"{datetime.now()} rc={rc}\n")
sys.exit(rc)
