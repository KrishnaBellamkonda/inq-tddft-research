#!/usr/bin/env python3
"""Autonomous orchestrator for the effective-mass quantum-vs-classical pair.

Pipeline (each step gated; a failure emails + stops, never wastes hours):
  1. build + run GS (dx=0.333 slab) if the shared checkpoint is absent  [GPU 0]
  2. build quantum WP binary; build classical binary                    [compile gate]
  3. launch quantum [GPU 0] and classical [GPU 1] concurrently
  4. wait for both (run_summary run_completed=true)
  5. build per-run notebooks + the model-comparison notebook
  6. email at GS-done, launch, and completion

Idempotent: skips a step whose output already exists. Headless:
  cd .../effmass_pair && nohup ../../../../../venv/bin/python3 orchestrate.py > orch.log 2>&1 &
"""
from __future__ import annotations
import os, subprocess, sys, time
from datetime import datetime
from pathlib import Path

ROOT   = Path("/local/data/public/skcb2/tddft")
LJ     = ROOT/"ResearchProject/systems/localised_jellium"
HERE   = LJ/"scripts/muon_mass_fork/effmass_pair"
GS_DIR = LJ/"shared_gs/slab_n82_L50x50x90_dx0p333"
NB_DIR = LJ/"hypotheses/muon_mass_fork"
PY     = str(ROOT/"venv/bin/python3")
TO     = "chiddukanna@gmail.com"
ENV = {**os.environ,
       "INQ_SOURCE":           str(ROOT/"inq-study"),
       "INQ_SHARE_PATH":       str(ROOT/"inq/install/share"),
       "PSEUDOPOD_SHARE_PATH": str(ROOT/"inq/install/share/pseudopod")}

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)

def email(subj, body, attach=None):
    try:
        sys.path.insert(0, str(ROOT/"inq-stack/python"))
        from inqview.email import send_run_email
        send_run_email(subject=f"[effmass-pair] {subj}", body=body,
                       attachments=[a for a in (attach or []) if Path(a).exists()], to=TO)
    except Exception as e:
        log(f"  email failed ({subj}): {e}")

def done(rundir: Path) -> bool:
    for rs in rundir.glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs.read_text(): return True
        except Exception: pass
    return False

def inq_run(rundir: Path, gpu: str, extra_env: dict, logname: str) -> int:
    """Build+run via inq-run in rundir on the given GPU. Returns rc."""
    env = {**ENV, **extra_env, "CUDA_VISIBLE_DEVICES": gpu}
    with open(rundir/logname, "w") as f:
        return subprocess.run(["inq-run"], cwd=rundir, env=env, stdout=f, stderr=subprocess.STDOUT).returncode

def build_only(rundir: Path, gpu: str) -> bool:
    """Compile the binary by running it once with a bogus GS so it exits at the
    GS-missing guard (rc 2). Returns True if the binary now exists."""
    if (rundir/"run").exists(): return True
    env = {**ENV, "CUDA_VISIBLE_DEVICES": gpu, "EM_GS_DIR": "/nonexistent_gs_compile_probe"}
    with open(rundir/"build.log", "w") as f:
        subprocess.run(["inq-run"], cwd=rundir, env=env, stdout=f, stderr=subprocess.STDOUT)
    return (rundir/"run").exists()

# ---------------------------------------------------------------- pipeline
def main():
    # 1. GS -- wait for the already-launched GS; launch only if none running --
    def gs_ready(): return GS_DIR.exists() and any(GS_DIR.iterdir())
    if not gs_ready():
        # A GS was launched separately iff gs/gs_run.log exists -> just POLL for it
        # (never double-launch on the same GPU). Only start GS if nobody has.
        if not (HERE/"gs/gs_run.log").exists():
            log("no GS started — launching GS on GPU0"); inq_run(HERE/"gs", "0", {}, "gs_run.log")
        log("waiting for GS checkpoint (poll 60s, up to 4h)")
        t_end = time.time() + 4*3600
        while not gs_ready() and time.time() < t_end:
            time.sleep(60)
        if not gs_ready():
            email("GS TIMEOUT/FAILED", "GS did not produce a checkpoint in 4h; see gs/gs_run.log. STOPPED.")
            log("GS timeout — stop"); return
    log("GS present ✓")
    email("GS ready", f"Shared GS at {GS_DIR} converged. Launching quantum (GPU0) + classical (GPU1).")

    # 2. compile both -------------------------------------------------------
    for name, gpu in [("quantum","0"), ("classical","1")]:
        if not build_only(HERE/name, gpu):
            email(f"{name} BUILD FAILED", f"see {name}/build.log. STOPPED.")
            log(f"{name} build failed — stop"); return
        log(f"{name} binary built ✓")

    # 3. launch — quantum on GPU0; classical on GPU1 IF it has room, else -----
    #    sequentially on GPU0 (a dx=0.333 run needs ~15-20 GB; never OOM overnight)
    def gpu_free_mb(idx):
        exe = HERE.parent/".gpuprobe"
        if not exe.exists(): return 999999
        try:
            out = subprocess.run([str(exe)], env={**os.environ,"CUDA_VISIBLE_DEVICES":str(idx)},
                                 capture_output=True, text=True)
            return int(out.stdout.strip() or "-1")
        except Exception: return -1

    def launch(name, gpu):
        f = open(HERE/name/"rt_run.log", "w")
        p = subprocess.Popen(["inq-run"], cwd=HERE/name,
                             env={**ENV, "CUDA_VISIBLE_DEVICES": gpu}, stdout=f, stderr=subprocess.STDOUT)
        log(f"{name} launched on GPU {gpu} (pid {p.pid})"); return p

    need_q = not done(HERE/"quantum"/"results")
    need_c = not done(HERE/"classical"/"results")
    g1 = gpu_free_mb(1)
    log(f"GPU1 free = {g1} MB")
    parallel = need_c and g1 >= 18000
    if not parallel and need_c:
        log(f"GPU1 too small ({g1} MB < 18 GB) — classical will run SEQUENTIALLY on GPU0 after quantum")

    pq = launch("quantum","0") if need_q else None
    pc = launch("classical","1") if (need_c and parallel) else None
    if pq: pq.wait(); log(f"quantum rc={pq.returncode} completed={done(HERE/'quantum'/'results')}")
    if pc: pc.wait(); log(f"classical rc={pc.returncode} completed={done(HERE/'classical'/'results')}")
    # sequential fallback: classical on GPU0 now that quantum is done
    if need_c and not parallel and not done(HERE/"classical"/"results"):
        pc = launch("classical","0"); pc.wait()
        log(f"classical (seq, GPU0) rc={pc.returncode} completed={done(HERE/'classical'/'results')}")
    for name in ("quantum","classical"):
        if not done(HERE/name/"results"):
            email(f"{name} RUN FAILED", f"see {name}/rt_run.log. Other run may still be ok.")

    q_ok = done(HERE/"quantum"/"results"); c_ok = done(HERE/"classical"/"results")
    if not (q_ok and c_ok):
        email("PAIR INCOMPLETE", f"quantum={q_ok} classical={c_ok}. Notebooks skipped for missing run.")

    # 5. notebooks ----------------------------------------------------------
    log("building comparison + run notebooks")
    rc = subprocess.run([PY, str(HERE/"build_notebooks.py")], env=ENV,
                        cwd=HERE, capture_output=True, text=True)
    log(rc.stdout[-2000:] if rc.stdout else "");
    if rc.returncode != 0: log("notebook build stderr:\n"+rc.stderr[-2000:])
    nb = NB_DIR/"effmass_pair_comparison.ipynb"

    # 6. done ---------------------------------------------------------------
    email("PAIR COMPLETE — quantum vs classical",
          f"quantum(σ_WP=2, m=3.09 WP) and classical(Gaussian-charge, m=3.09) through the "
          f"r_s=5.665 slab, dt=0.04.\nquantum_ok={q_ok} classical_ok={c_ok}\n"
          f"Comparison notebook: {nb}\nRun notebooks under {NB_DIR}.",
          attach=[str(NB_DIR/'effmass_pair_stopping.png')])
    log("DONE")

if __name__ == "__main__":
    main()
