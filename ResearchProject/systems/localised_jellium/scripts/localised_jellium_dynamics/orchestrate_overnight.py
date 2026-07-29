#!/usr/bin/env python3
"""Overnight classical-vs-wavepacket TWIN campaign (localised jellium).

Runs a fable-informed, budget-sized matrix of twin pairs (classical proj_dyn +
WP phase5_wp), each identical except the projectile representation, to extract
maximally-diverse quantum (classical-vs-WP) differences. Every run emits the full
energy decomposition INCLUDING the pairwise P/S/B Coulomb terms (interactions.csv),
density frames, and a final checkpoint. Idempotent resume; ~9h soft budget
(checkpoint-and-continue, never blocks); per-pair analysis + Gmail.

Design: docs/plans/wide-wp-and-notebook-enhancements.md + the fable agent advice.
Pairs (m=1, dispatch order = validate-first then highest quantum value):
  p5 null control (falsifies framework)  | p1 quantum reflection | p4 capture-vs-escape
  p2 tunnelling from inside slab          | p6 sigma-ladder ZPE/SIE
"""
import os, subprocess, sys, time, math, traceback
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
LJD  = ROOT/"ResearchProject/systems/localised_jellium/scripts/localised_jellium_dynamics"
GS   = ROOT/"ResearchProject/systems/localised_jellium/scripts/campaign_autorun/runs/h2/gs_p2_lz120/checkpoint"
PY   = str(ROOT/"venv/bin/python3")
ENGINE = str(ROOT/".claude/skills/twin-run-analysis/twin_decompose.py")
CHECK  = str(ROOT/".claude/skills/twin-run-generation/check_twin.py")
GPU = "0"
BUDGET_S = 9*3600
START = time.time()

ENV = {**os.environ,
       "INQ_SHARE_PATH": str(ROOT/"inq/install/share"),
       "PSEUDOPOD_SHARE_PATH": str(ROOT/"inq/install/share/pseudopod"),
       "INQ_SOURCE": str(ROOT/"inq-study"),
       "TMPDIR": str(ROOT/".buildtmp"),
       "CUDA_VISIBLE_DEVICES": GPU}

PROJ_DYN, PHASE5_WP = LJD/"proj_dyn", LJD/"phase5_wp"

COMMON = dict(LJ_GS_DIR=str(GS), LJ_PERIODICITY=2, LJ_LZ=120, LJ_HALF=12.5,
              LJ_N=82, LJ_SPACING=0.5, LJ_SAVE_EVERY=25, LJ_MASS=1)

PAIRS = [
  dict(name="p5_null_s2_k4",     sigma=2.0, k0=4.2, launch=-24.5, dt=0.025, nsteps=300, phenom="null control (Δ→0 except SIE; framework falsifier)"),
  dict(name="p1_reflect_s2_k04", sigma=2.0, k0=0.4, launch=-17.0, dt=0.05,  nsteps=400, phenom="quantum reflection at the attractive surface"),
  dict(name="p4_capture_s2_k11", sigma=2.0, k0=1.1, launch=-17.0, dt=0.05,  nsteps=350, phenom="capture-vs-escape branching + plasmon ringing"),
  dict(name="p2_tunnel_s2_k05",  sigma=2.0, k0=0.5, launch=0.0,   dt=0.05,  nsteps=350, phenom="tunnelling: launched INSIDE the slab"),
  dict(name="p6_ladder_s1_k11",  sigma=1.0, k0=1.1, launch=-17.0, dt=0.05,  nsteps=350, phenom="sigma-ladder ZPE/SIE scaling (vs sigma=2, 0.5)"),
]

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)
def budget_left_h(): return (BUDGET_S-(time.time()-START))/3600.0

def email(subject, body, attach=None):
    try:
        from inqview.email import send_run_email
        send_run_email(subject=subject, body=body, attachments=attach or [])
        log(f"  emailed: {subject}")
    except Exception as e:
        log(f"  email FAILED ({e})")

def cutoff_ok(sigma, k0):
    sp = 1.0/(math.sqrt(2)*sigma); kmax = k0 + 3*sp; nyq = math.pi/0.5
    return kmax < 0.9*nyq, kmax, nyq

def _done(rundir):
    for rs in Path(rundir).glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs.read_text(): return True
        except Exception: pass
    return False

_GPU_PROBE = (
    "import ctypes,sys\n"
    "c=None\n"
    "for lib in ('libcudart.so','libcudart.so.12',"
    "'/lsc/opt/cuda-12.6.2/lib64/libcudart.so',"
    "'/lsc/opt/cuda-12.6.2/targets/x86_64-linux/lib/libcudart.so'):\n"
    " try: c=ctypes.CDLL(lib); break\n"
    " except OSError: pass\n"
    "if c is None: sys.exit(2)\n"
    "f=ctypes.c_size_t(); t=ctypes.c_size_t(); c.cudaSetDevice(0)\n"
    "sys.exit(0 if c.cudaMemGetInfo(ctypes.byref(f),ctypes.byref(t))==0 and f.value>0.85*t.value else 1)\n"
)

def wait_for_gpu(max_wait_s=3600):
    """Wait until PHYSICAL GPU0 is free — device-specific cudaMemGetInfo probe run in
    a throwaway subprocess (holds no CUDA context). Ignores jobs on other GPUs
    (a bare `/run$` grep wrongly matched a GPU1 job)."""
    t0 = time.time()
    while time.time()-t0 < max_wait_s:
        rc = subprocess.run([PY, "-c", _GPU_PROBE]).returncode
        if rc in (0, 2):   # 0 = free, 2 = cannot-probe (assume free, avoid deadlock)
            log("  GPU0 free"); return True
        log("  GPU0 in use (waiting 120s)…"); time.sleep(120)
    log("  wait_for_gpu timed out — proceeding anyway"); return False

INQ_RUN = str(ROOT/"shared/bin/inq-run")   # builds (incremental) + runs; idempotent build

def run_one(bindir, out, cfg, is_wp):
    """One twin run (classical or WP) via inq-run: incremental build + run,
    checkpointed + idempotent + resume."""
    rundir = bindir/"results"/out
    if _done(rundir):
        log(f"  SKIP {out} (done)"); return True
    ov = {**COMMON, "LJ_OUT": out, "LJ_SIGMA": cfg["sigma"], "LJ_K0": cfg["k0"],
          "LJ_LAUNCH_Z": cfg["launch"], "LJ_DT": cfg["dt"], "LJ_N_STEPS": cfg["nsteps"]}
    if (rundir/"checkpoint").exists(): ov["LJ_RESUME"]=1; log(f"  RESUME {out}")
    env = {**ENV, **{k:str(v) for k,v in ov.items()}}
    log(f"  RUN {out} (budget {budget_left_h():.1f}h left)")
    with open(bindir/f"orch_{out}.log","w") as lf:
        rc = subprocess.run([INQ_RUN], cwd=str(bindir), env=env,
                            stdout=lf, stderr=subprocess.STDOUT).returncode
    ok = rc==0 and _done(rundir)
    log(f"  {out} {'OK' if ok else 'FAIL rc='+str(rc)}")
    return ok

def analyse_pair(p):
    cl = PROJ_DYN/"results"/(p["name"]+"_cl")
    wp = PHASE5_WP/"results"/(p["name"]+"_wp")
    rep = subprocess.run([PY, ENGINE, "--wp", str(wp), "--classical", str(cl)],
                         capture_output=True, text=True, env=ENV)
    txt = rep.stdout or rep.stderr
    (LJD/f"analysis_{p['name']}.txt").write_text(txt)
    return txt

def do_pair(p):
    ok_g, kmax, nyq = cutoff_ok(p["sigma"], p["k0"])
    if not ok_g:
        log(f"PAIR {p['name']} ALIASED (kmax {kmax:.2f} vs 0.9·nyq {0.9*nyq:.2f}) — SKIP")
        email(f"[twin-campaign] {p['name']} SKIPPED (aliasing)",
              f"kmax={kmax:.2f} exceeds 0.9*Nyquist={0.9*nyq:.2f}. {p['phenom']}")
        return
    log(f"=== PAIR {p['name']}: {p['phenom']} ===")
    cl_ok = run_one(PROJ_DYN,  p["name"]+"_cl", p, is_wp=False)
    wp_ok = run_one(PHASE5_WP, p["name"]+"_wp", p, is_wp=True)
    if cl_ok and wp_ok:
        txt = analyse_pair(p)
        email(f"[twin-campaign] {p['name']} DONE — {p['phenom']}",
              f"σ={p['sigma']} k0={p['k0']} launch={p['launch']} dt={p['dt']} N={p['nsteps']}\n\n{txt}")
    else:
        email(f"[twin-campaign] {p['name']} PARTIAL (cl={cl_ok} wp={wp_ok})",
              f"One twin did not complete; checkpointed, resumable. {p['phenom']}")

def main():
    log(f"OVERNIGHT TWIN CAMPAIGN start (GPU={GPU}, budget {BUDGET_S/3600:.0f}h)")
    email("[twin-campaign] STARTED",
          "Overnight classical-vs-WP twin campaign started.\nPairs: " +
          ", ".join(p["name"] for p in PAIRS) + f"\nBudget ~{BUDGET_S/3600:.0f}h on GPU{GPU}.")
    wait_for_gpu()
    for p in PAIRS:
        if budget_left_h() < 0.3:
            log(f"BUDGET reached ({budget_left_h():.2f}h left) — stopping before {p['name']}")
            email("[twin-campaign] BUDGET REACHED — stopping (resumable)",
                  f"Stopped before {p['name']}; completed pairs are done, remainder checkpointed.")
            break
        try:
            do_pair(p)
        except Exception:
            tb = traceback.format_exc(); log(f"PAIR {p['name']} EXCEPTION:\n{tb}")
            email(f"[twin-campaign] {p['name']} EXCEPTION (chain continues)", tb)
    log(f"CAMPAIGN done ({(time.time()-START)/3600:.1f}h elapsed)")
    email("[twin-campaign] COMPLETE",
          f"Elapsed {(time.time()-START)/3600:.1f}h. See analysis_<pair>.txt + per-pair emails.\n"
          "Follow-ups: m=10 Bragg/image pairs (need WP fictitious mass), notebooks + GIFs.")

if __name__ == "__main__":
    main()
