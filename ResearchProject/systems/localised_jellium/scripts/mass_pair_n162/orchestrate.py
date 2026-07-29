#!/usr/bin/env python3
"""Autonomous orchestrator — GENUINE 162-electron localised-jellium mass pair.

Two matched quantum-wavepacket runs, identical except the projectile mass:
  * m1 : EM_INV_MASS=1.0, k0=2.711  (projectile mass 1)  -> GPU 0
  * m2 : EM_INV_MASS=0.5, k0=3.834  (projectile mass 2)  -> GPU 1
Bath electrons stay mass 1 in both. sigma_WP=1, E=100 eV, CAP 10 Bohr eta=-1.0,
100 a.u. (2500 steps, dt=0.04), checkpoint every 500 steps.
Plan: docs/plans/mass-pair-n162-sigma1-cap.md.

Pipeline (each gate emails; correctness failures stop, cost does NOT — the
checkpoint-dont-block rule):
  1. wait for the N=162 GS checkpoint (build launched separately; poll up to 4h)
  2. build the WP binary once (inq-study engine)
  3. PILOT gate: 24-step m1 smoke -> memory fit + WP norm + finite energy +
     s/step projection. Blocks ONLY on crash/OOM/NaN/norm; warns on drift.
  4. launch m1 (GPU0) + m2 (GPU1) concurrently; auto-resume (EM_RESUME=1) on crash
  5. wait for both (run_completed=true); best-effort analyse.py per run
  6. email at each milestone + on any failure

Idempotent + resume-aware: a completed run is skipped; a partial run resumes from
its last 500-step checkpoint. Headless launch:
  cd .../mass_pair_n162 && nohup ../../../../venv/bin/python3 orchestrate.py > orch.log 2>&1 &
"""
from __future__ import annotations
import os, re, subprocess, sys, time
from datetime import datetime
from pathlib import Path

ROOT   = Path("/local/data/public/skcb2/tddft")
LJ     = ROOT/"ResearchProject/systems/localised_jellium"
HERE   = LJ/"scripts/mass_pair_n162"
WP     = HERE/"wp"
GSDIR  = LJ/"shared_gs/slab_n120_L60x60x62_dx0p40"
PY     = str(ROOT/"venv/bin/python3")
PROBE  = LJ/"scripts/../../vacuum/gpu_probe"   # NVML-independent free-MB probe
TO     = "chiddukanna@gmail.com"

ENV = {**os.environ,
       "PATH":                 f"{ROOT/'shared/bin'}:" + os.environ.get("PATH",""),
       "INQ_SOURCE":           str(ROOT/"inq-study"),
       "INQ_SHARE_PATH":       str(ROOT/"inq/install/share"),
       "PSEUDOPOD_SHARE_PATH": str(ROOT/"inq/install/share/pseudopod")}

# shared run parameters (both masses)
# Observable set: current/dipole OFF and momentum_distribution OFF — both allocate
# a full extra field on the 24 GB A30 and caused a step-0 CUDA illegal-access on
# this 176x176x156 / 99-state grid (diagnosed 2026-07-19). The per-step full energy
# decomposition (energy_decomp.csv) + density VTIs + wp stats + overlaps + state
# energies remain. See docs/handovers/mass-pair-n162-sigma1-cap.md.
BASE = {"EM_DT":"0.04", "EM_N_STEPS":"2500", "EM_WRITE_EVERY":"8", "EM_WF_EVERY":"40",
        "EM_CKPT_EVERY":"500", "EM_CAP":"1", "EM_CAP_ETA":"-1.0",
        "EM_CAP_CENTER_BOHR":"26.2", "EM_CAP_WIDTH_BOHR":"10.0",
        "EM_SIGMA_WP":"1.0", "EM_LAUNCH_Z":"-16.5", "EM_GS_DIR":str(GSDIR),
        "EM_OBS_DIPCUR":"0", "EM_OBS_MOM":"0", "EM_OBS_OVL":"0",
        "EM_OBS_WF":"0", "EM_OBS_STATE":"0"}
# NOTE: the base 99-state ETRS propagation sits near the 24 GB A30 ceiling, so
# EVERY all-states observable (current/dipole, momentum FFT, overlaps, state
# energies) independently OOMs at step 0 (diagnosed 2026-07-19). Only the
# always-on core fits: per-step energy_decomp, density VTIs (total/system/wp/
# delta), wp_momentum_stats, wp_real_space_stats. Dropped extras are recoverable
# in post from the saved density fields.
RUNS = {  # out-subdir -> (gpu, per-run env)
    "m1": ("0", {"EM_INV_MASS":"1.0", "EM_K0":"2.711"}),
    "m2": ("1", {"EM_INV_MASS":"0.5", "EM_K0":"3.834"}),
}

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)

def email(subj, body, attach=None):
    try:
        sys.path.insert(0, str(ROOT/"inq-stack/python"))
        from inqview.email import send_run_email
        send_run_email(subject=f"[mass_pair_n162] {subj}", body=body,
                       attachments=[a for a in (attach or []) if a and Path(a).exists()], to=TO)
        log(f"  emailed: {subj}")
    except Exception as e:
        log(f"  email FAILED ({subj}): {e}")

def gs_ready() -> bool:
    # require BOTH the checkpoint AND the summary's run_completed (written last)
    if not (GSDIR.exists() and any(GSDIR.iterdir())):
        return False
    try:
        return "run_completed = true" in (HERE/"gs/results/run_summary.txt").read_text()
    except Exception:
        return False

def gpu_free_mb(gpu: str) -> int:
    """Free MB on the given GPU via the vacuum gpu_probe (NVML-independent)."""
    try:
        out = subprocess.run([str(PROBE)], env={**os.environ, "CUDA_VISIBLE_DEVICES": str(gpu)},
                             capture_output=True, text=True, timeout=30).stdout
        import re as _re
        m = _re.search(r"free_MB\s+(\d+)", out)
        return int(m.group(1)) if m else -1
    except Exception:
        return -1

def wait_gpu_free(gpu: str, min_mb: int = 20000, timeout: int = 900) -> bool:
    """Block until the GPU has >= min_mb free (avoids the CUDA teardown race)."""
    t_end = time.time() + timeout
    while time.time() < t_end:
        f = gpu_free_mb(gpu)
        if f >= min_mb:
            return True
        log(f"  waiting for GPU{gpu} to free (now {f} MB, need {min_mb})")
        time.sleep(15)
    return False

def run_done(out: str) -> bool:
    rs = WP/"results"/out/"run_summary.txt"
    try:
        return rs.exists() and "run_completed = true" in rs.read_text()
    except Exception:
        return False

def has_ckpt(out: str) -> int:
    st = WP/"results"/out/"rt_ckpt"/"rt_state.txt"
    try:
        m = re.search(r"last_step=(\d+)", st.read_text())
        return int(m.group(1)) if m else 0
    except Exception:
        return 0

def build_wp() -> bool:
    """Compile ./run once (bogus GS -> exits at the GS-missing guard, rc 2)."""
    if (WP/"run").exists():
        return True
    env = {**ENV, "CUDA_VISIBLE_DEVICES":"0", "EM_GS_DIR":"/nonexistent_compile_probe"}
    with open(WP/"build.log","w") as f:
        subprocess.run(["inq-run"], cwd=WP, env=env, stdout=f, stderr=subprocess.STDOUT)
    return (WP/"run").exists()

def launch(out: str, gpu: str, extra: dict, resume: bool):
    if not wait_gpu_free(gpu):
        log(f"  WARN: GPU{gpu} never freed; launching {out} anyway")
    env = {**ENV, **BASE, **extra, "EM_OUT":out, "CUDA_VISIBLE_DEVICES":gpu,
           "EM_RESUME":"1" if resume else "0"}
    logf = open(WP/f"rt_{out}.log","a" if resume else "w")
    p = subprocess.Popen([str(WP/"run")], cwd=WP, env=env, stdout=logf, stderr=subprocess.STDOUT)
    log(f"  {out} launched on GPU{gpu} pid={p.pid} resume={resume}")
    return p

def parse_kv(text: str, key: str):
    m = re.search(rf"{re.escape(key)}\s*=\s*([-\d.eE+]+)", text)
    return float(m.group(1)) if m else None

# --------------------------------------------------------------------- pipeline
def main():
    log("=== mass_pair_n162 orchestrator start ===")

    # 1. GS ------------------------------------------------------------------
    if not gs_ready():
        log("waiting for GS checkpoint (poll 60s, up to 4h)")
        t_end = time.time() + 4*3600
        while not gs_ready() and time.time() < t_end:
            time.sleep(60)
    if not gs_ready():
        email("GS TIMEOUT", f"No GS checkpoint at {GSDIR} within 4h. See gs/gs_build_run.log. STOPPED.")
        log("GS timeout — stop"); return
    gs_summary = ""
    try: gs_summary = (HERE/"gs/results/run_summary.txt").read_text()
    except Exception: pass
    egs = parse_kv(gs_summary, "ground_state_energy_ha")
    if egs is None or egs != egs:  # NaN check
        email("GS INVALID", f"GS energy not finite ({egs}). See gs/results/run_summary.txt. STOPPED.")
        log("GS invalid — stop"); return
    log(f"GS ready ✓  E_gs={egs} Ha")
    email("GS ready", f"N=162 localised-jellium GS converged: E_gs={egs} Ha, r_s≈5.68.\n"
                      f"Checkpoint: {GSDIR}\nBuilding WP binary + running pilot next.")

    # 2. build WP ------------------------------------------------------------
    if not build_wp():
        email("WP BUILD FAILED", f"WP binary did not compile. See {WP}/build.log. STOPPED — no runs executed.")
        log("WP build failed — stop"); return
    log("WP binary built ✓")

    # 3. PILOT gate (correctness only) --------------------------------------
    if not (WP/"results/pilot/run_summary.txt").exists():
        log("running 24-step pilot (m1) on GPU0")
        if not wait_gpu_free("0"):
            log("  WARN: GPU0 never freed before pilot; proceeding anyway")
        t0 = time.time()
        env = {**ENV, **BASE, **RUNS["m1"][1], "EM_OUT":"pilot", "EM_N_STEPS":"24",
               "EM_WF_EVERY":"1000", "EM_CKPT_EVERY":"100000", "CUDA_VISIBLE_DEVICES":"0"}
        with open(WP/"rt_pilot.log","w") as f:
            rc = subprocess.run([str(WP/"run")], cwd=WP, env=env, stdout=f, stderr=subprocess.STDOUT).returncode
        dt_wall = time.time() - t0
    else:
        rc = 0; dt_wall = None
    psum = ""
    try: psum = (WP/"results/pilot/run_summary.txt").read_text()
    except Exception: pass
    completed = "run_completed = true" in psum
    norm = parse_kv(psum, "wp_norm_after")
    # energy sanity from per-step decomposition
    drift = None; nan_energy = False
    try:
        import csv
        rows = list(csv.DictReader(open(WP/"results/pilot/raw/observables/energy_decomp.csv")))
        if rows:
            e0 = float(rows[0]["energy_total"]); e1 = float(rows[-1]["energy_total"])
            drift = e1 - e0
            nan_energy = any(v != v for v in (e0, e1))
    except Exception as e:
        log(f"  pilot energy read failed: {e}")
    per_step = (dt_wall/24.0) if dt_wall else None
    proj_h = (per_step*2500/3600.0) if per_step else None

    hard_fail = (not completed) or nan_energy or (norm is not None and not (0.9 <= norm <= 1.1))
    pilot_report = (f"completed={completed}  wp_norm_after={norm}  dE_total(24 steps)={drift}\n"
                    f"pilot wall={dt_wall}s  ~{per_step:.1f}s/step  =>  ~{proj_h:.1f} h per 2500-step run"
                    if per_step else
                    f"completed={completed}  wp_norm_after={norm}  dE_total(24 steps)={drift}")
    if hard_fail:
        email("PILOT FAILED — production NOT launched",
              "Hard correctness failure (crash/OOM/NaN/norm). Production runs NOT started to avoid "
              f"wasting GPU hours.\n\n{pilot_report}\n\nSee {WP}/rt_pilot.log.")
        log(f"PILOT hard-fail — stop. {pilot_report}"); return
    log(f"PILOT ok ✓  {pilot_report}")
    email("PILOT passed — launching both production runs",
          f"Pilot green (memory fits, WP norm≈1, energy finite).\n\n{pilot_report}\n\n"
          f"Launching m1 (mass 1, GPU0) + m2 (mass 2, GPU1), 2500 steps, checkpoint every 500. "
          f"Projected ~{proj_h:.0f} h/run if shown." if proj_h else
          f"Pilot green.\n\n{pilot_report}\n\nLaunching m1 (GPU0) + m2 (GPU1).")

    # 4. production (parallel, auto-resume on crash) ------------------------
    MAX_TRIES = 3
    tries = {o:0 for o in RUNS}
    procs = {}
    for out,(gpu,extra) in RUNS.items():
        if run_done(out): log(f"{out} already complete — skip"); continue
        procs[out] = launch(out, gpu, extra, resume=has_ckpt(out) > 0); tries[out] += 1

    # monitor loop
    while procs:
        time.sleep(120)
        for out in list(procs):
            p = procs[out]
            if p.poll() is None:
                continue  # still running
            del procs[out]
            if run_done(out):
                log(f"{out} COMPLETE ✓ (rc={p.returncode})")
                email(f"{out} run complete", f"{out} (projectile mass {'1' if out=='m1' else '2'}) finished "
                      f"2500 steps.\nResults: {WP}/results/{out}/\nRunning analysis next.")
                continue
            # unexpected exit before completion -> resume if we can
            ls = has_ckpt(out)
            if tries[out] < MAX_TRIES and ls > 0:
                log(f"{out} died at rc={p.returncode}, last_ckpt={ls} — resuming (try {tries[out]+1})")
                email(f"{out} auto-resume", f"{out} exited early (rc={p.returncode}) at step {ls}; "
                      f"resuming from checkpoint (attempt {tries[out]+1}/{MAX_TRIES}).")
                gpu = RUNS[out][0]
                procs[out] = launch(out, gpu, RUNS[out][1], resume=True); tries[out] += 1
            else:
                log(f"{out} FAILED permanently (rc={p.returncode}, tries={tries[out]})")
                email(f"{out} RUN FAILED", f"{out} did not complete after {tries[out]} attempts "
                      f"(rc={p.returncode}). Last checkpoint step={ls}. See {WP}/rt_{out}.log. "
                      f"Data up to the last checkpoint is preserved (EM_RESUME=1 to continue).")

    # 5. analysis (best-effort; never blocks) -------------------------------
    for out in RUNS:
        if not run_done(out): continue
        try:
            if (HERE/"analyse.py").exists():
                log(f"analyse.py on {out}")
                r = subprocess.run([PY, str(HERE/"analyse.py"), "--run", str(WP/"results"/out)],
                                   cwd=HERE, env=ENV, capture_output=True, text=True, timeout=3*3600)
                log((r.stdout or "")[-1500:]);
                if r.returncode: log("analyse stderr:\n"+(r.stderr or "")[-1500:])
        except Exception as e:
            log(f"  analyse {out} failed: {e}")

    # 6. final -------------------------------------------------------------
    ok = {o: run_done(o) for o in RUNS}
    email("MASS PAIR — final status",
          f"Genuine N=162 localised-jellium mass pair (sigma_WP=1, E=100 eV, CAP eta=-1.0, 100 a.u.).\n"
          f"m1 (projectile mass 1) complete = {ok['m1']}\n"
          f"m2 (projectile mass 2) complete = {ok['m2']}\n"
          f"Results: {WP}/results/{{m1,m2}}/\n"
          f"Per-step full energy decomposition: results/<run>/raw/observables/energy_decomp.csv")
    log(f"DONE  m1={ok['m1']} m2={ok['m2']}")

if __name__ == "__main__":
    main()
