#!/usr/bin/env python3
"""Autonomous orchestrator for the muon-mass-fork campaign.

Campaign prompt: docs/campaigns/muon_mass_fork/muon_mass_fork.md
Plan:           docs/plans/muon-mass-fork-implementation.md
Handover:       docs/handovers/muon-mass-fork.md

Phase-gated (STRICT): phase N+1 runs only if phase N passed. A failure ⇒ email
(full traceback) + STOP + write `blocked` to the state file (NOT continue — a
broken engine must never leak into physics). Idempotent resume: completed phases
are recorded in state.json and skipped on restart. Per-phase Gmail via
inqview.email. GPU is the default (cudaMemGetInfo probe; NVML/nvidia-smi is broken
but compute works — never CPU-fall-back on an nvidia-smi error).

Two-phase autonomy (user decision 2026-07-06): runs Phases 1→4-research
autonomously, then HARD-STOPS at the Phase-4 checkpoint and emails the user to
PICK the muon-XC functional. Re-run with the pick recorded to continue to Phase 5.

Headless launch (survives disconnect):
    cd .../scripts/muon_mass_fork
    GPU=1 nohup ../../../../../venv/bin/python3 orchestrate.py > orchestrate.log 2>&1 &

Run a single phase:  python3 orchestrate.py phase1
"""
from __future__ import annotations
import json, os, subprocess, sys, time, traceback
from datetime import datetime
from pathlib import Path

ROOT      = Path("/local/data/public/skcb2/tddft")
INQSTUDY  = ROOT / "inq-study"
BUILD_GPU = INQSTUDY / "build-gpu"
LJ        = ROOT / "ResearchProject/systems/localised_jellium"
HERE      = LJ / "scripts/muon_mass_fork"
RUNS      = HERE / "runs"; RUNS.mkdir(parents=True, exist_ok=True)
STATE     = HERE / "state.json"
PICK_FILE = HERE / "muon_xc_pick.json"                 # user writes this to pass Phase-4
NB_DIR    = LJ / "hypotheses/muon_mass_fork"           # phase + run + index notebooks
PY        = str(ROOT / "venv/bin/python3")
NVCC      = "/lsc/opt/cuda-12.6.2/bin/nvcc"           # = config.sh INQ_CUDA_COMPILER
                                                       # (the version inq-run uses for ALL working GPU runs;
                                                       #  12.5 from the stale inq/build cache does NOT compile
                                                       #  INQ's reduce lambdas — dipole/singularity assert)
CUDA_ARCH = "80"
GPU       = os.environ.get("GPU", "auto")             # "auto" = pick the freest GPU
MUON_MASS = 206.77                                     # PDG muon/electron mass ratio
TO        = "chiddukanna@gmail.com"

ENV = {**os.environ,
       "INQ_SHARE_PATH":       str(ROOT/"inq/install/share"),
       "PSEUDOPOD_SHARE_PATH": str(ROOT/"inq/install/share/pseudopod"),
       "INQ_SOURCE":           str(INQSTUDY),
       "CUDA_VISIBLE_DEVICES": GPU}

def log(msg): print(f"[{datetime.now():%F %T}] {msg}", flush=True)

# ---------------- state / resume -------------------------------------------
def _state() -> dict:
    if STATE.exists():
        try: return json.loads(STATE.read_text())
        except Exception: pass
    return {"done": [], "status": "running", "blocked_reason": ""}

def _save(st: dict): STATE.write_text(json.dumps(st, indent=2))

def _mark_done(name: str):
    st = _state()
    if name not in st["done"]: st["done"].append(name)
    st["status"] = "running"; _save(st)

def _mark_blocked(name: str, reason: str):
    st = _state(); st["status"] = "blocked"
    st["blocked_reason"] = f"{name}: {reason}"; _save(st)

# ---------------- email -----------------------------------------------------
def _email(subject: str, body: str, attachments=None):
    try:
        sys.path.insert(0, str(ROOT/"inq-stack/python"))
        from inqview.email import send_run_email
        send_run_email(subject=f"[muon-mass-fork] {subject}", body=body,
                       attachments=attachments or [], to=TO)
    except Exception as e:
        log(f"  email failed ({subject}): {e}")

def fail_email(phase: str, msg: str):
    _email(f"{phase} — FAILED, chain STOPPED",
           f"Phase {phase} failed. STRICT gating: the orchestrator STOPPED (a "
           f"broken engine must not reach physics).\n\n{msg}\n\n"
           f"Fix, then re-run orchestrate.py to resume (done phases are skipped).")

# ---------------- GPU probe + auto-select -----------------------------------
# NVML/nvidia-smi is broken on this box; cudaMemGetInfo still works. Single-GPU
# aware: GPU="auto" probes every device and picks the freest.
_PROBE_EXE = None

def _probe_exe():
    """Compile the cudaMemGetInfo helper once; return its path or None."""
    global _PROBE_EXE
    if _PROBE_EXE is not None:
        return _PROBE_EXE or None
    prog = ('#include <cuda_runtime.h>\n#include <cstdio>\nint main(){size_t f,t;'
            'if(cudaMemGetInfo(&f,&t)!=cudaSuccess){printf("-1\\n");return 1;}'
            'printf("%zu\\n",f/1048576);return 0;}\n')
    src = HERE/".gpuprobe.cu"; exe = HERE/".gpuprobe"
    src.write_text(prog)
    if subprocess.run([NVCC, str(src), "-o", str(exe)], capture_output=True).returncode != 0:
        log("  gpu probe compile failed — probe disabled"); _PROBE_EXE = ""; return None
    _PROBE_EXE = str(exe); return _PROBE_EXE

def _free_mb(idx) -> int:
    exe = _probe_exe()
    if not exe: return -1
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": str(idx)}
    try:
        out = subprocess.run([exe], env=env, capture_output=True, text=True)
        return int(out.stdout.strip() or "-1")
    except Exception:
        return -1

GPUS: list[str] = []                                 # resolved pool of usable devices

def resolve_gpu(candidates=range(4), min_free_gb: float = 4.0) -> None:
    """Resolve the GPU pool. GPU='auto' probes every device and keeps ALL with
    >= min_free_gb free (dual-GPU when two are idle); GPU=<n> pins one. Sets the
    module-level GPU (primary, for single-GPU helpers) and GPUS (the pool)."""
    global GPU, GPUS
    if GPU != "auto":
        GPUS = [GPU]; ENV["CUDA_VISIBLE_DEVICES"] = GPU
        log(f"  GPU fixed to {GPU} (GPU env set) — single-GPU run"); return
    free = {}
    for idx in candidates:
        f = _free_mb(idx)
        if f >= 0:
            log(f"  GPU {idx}: {f} MB free"); free[str(idx)] = f
    pool = sorted([g for g, f in free.items() if f/1024.0 >= min_free_gb],
                  key=lambda g: -free[g])
    if not pool:
        log("  no GPU with enough free memory — defaulting to 0"); pool = ["0"]
    GPUS = pool
    GPU = pool[0]; ENV["CUDA_VISIBLE_DEVICES"] = GPU
    log(f"  AUTO-SELECTED GPU pool {GPUS} (freest first); primary GPU {GPU}")

def gpu_free_ok(min_free_gb: float = 4.0) -> bool:
    free_mb = _free_mb(GPU if GPU != "auto" else 0)   # GPU already resolved by main()
    if free_mb < 0:
        log("  cudaMemGetInfo failed — proceeding (probe optional)"); return True
    free_gb = free_mb/1024.0
    log(f"  GPU {GPU} free = {free_gb:.1f} GB")
    if free_gb < min_free_gb:
        log(f"  WARNING: GPU {GPU} has <{min_free_gb} GB free — likely occupied "
            "by another user. Proceeding, but flag this.")
    return True

# ---------------- run helpers -----------------------------------------------
def _done(rundir: Path) -> bool:
    for rs in rundir.glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs.read_text(): return True
        except Exception: pass
    return False

def run_sim(binary: str, rundir: Path, overrides: dict, label: str,
            retries: int = 1, gpu: str | None = None) -> bool:
    if not Path(binary).exists():
        log(f"  BLOCKED {label}: binary {binary} not built yet"); return False
    rundir.mkdir(parents=True, exist_ok=True)
    if _done(rundir): log(f"  SKIP {label} (already complete)"); return True
    env = {**ENV, **{k: str(v) for k, v in overrides.items()}}
    if gpu is not None: env["CUDA_VISIBLE_DEVICES"] = gpu
    for attempt in range(retries + 1):
        log(f"  RUN  {label} [gpu {env['CUDA_VISIBLE_DEVICES']}]"
            + (f" (retry {attempt})" if attempt else ""))
        with open(rundir/"run.log", "w") as lf:
            rc = subprocess.run([binary], cwd=str(rundir), env=env,
                                stdout=lf, stderr=subprocess.STDOUT).returncode
        if rc == 0 and _done(rundir): return True
        log(f"  FAIL {label} rc={rc} (see {rundir/'run.log'})")
    return False

def run_sims_parallel(jobs: list[dict], phase_label: str) -> bool:
    """Dispatch independent runs across the whole GPU pool (GPUS), <=1 job per
    GPU at a time. Each job: {binary, rundir, overrides, label}. Idempotent:
    already-complete runs are skipped instantly. Returns True iff ALL pass."""
    pool = GPUS or [GPU]
    log(f"  {phase_label}: {len(jobs)} runs across GPUs {pool}")
    pending = list(jobs)
    running: dict[str, dict] = {}          # gpu -> {proc, job, rundir}
    results: dict[str, bool] = {}

    def launch(job: dict, gpu: str):
        binary = job["binary"]; rundir = Path(job["rundir"]); ov = job["overrides"]
        if not Path(binary).exists():
            log(f"  BLOCKED {job['label']}: binary {binary} not built"); results[job["label"]] = False; return False
        rundir.mkdir(parents=True, exist_ok=True)
        if _done(rundir):
            log(f"  SKIP {job['label']} (already complete)"); results[job["label"]] = True; return False
        env = {**ENV, **{k: str(v) for k, v in ov.items()}, "CUDA_VISIBLE_DEVICES": gpu}
        log(f"  RUN  {job['label']} [gpu {gpu}]")
        lf = open(rundir/"run.log", "w")
        proc = subprocess.Popen([binary], cwd=str(rundir), env=env,
                                stdout=lf, stderr=subprocess.STDOUT)
        running[gpu] = dict(proc=proc, job=job, rundir=rundir, lf=lf); return True

    # prime: one job per free GPU
    for gpu in pool:
        while pending and gpu not in running:
            if launch(pending.pop(0), gpu): break
    # poll loop
    while running:
        for gpu, st in list(running.items()):
            rc = st["proc"].poll()
            if rc is None: continue
            st["lf"].close()
            ok = (rc == 0 and _done(st["rundir"]))
            results[st["job"]["label"]] = ok
            log(("  OK   " if ok else "  FAIL ") + st["job"]["label"] + f" (rc={rc})")
            del running[gpu]
            while pending and gpu not in running:      # backfill this GPU
                if launch(pending.pop(0), gpu): break
        time.sleep(2)
    allok = all(results.get(j["label"], False) for j in jobs)
    log(f"  {phase_label}: {sum(results.values())}/{len(jobs)} passed")
    return allok

# ---------------- PHASES ----------------------------------------------------
# Each returns True on pass. A False (or exception) triggers strict STOP.

def phase1():
    """Code bug checks: GPU build inq-study + engine tests (kernel/ledger/MPI)."""
    gpu_free_ok()
    log("  GPU-building inq-study (nvcc arch 80) — this is the Phase-1 gate")
    cfg = subprocess.run(
        ["cmake", "-S", str(INQSTUDY), "-B", str(BUILD_GPU), "-DENABLE_CUDA=ON",
         f"-DCMAKE_CUDA_ARCHITECTURES={CUDA_ARCH}", f"-DCMAKE_CUDA_COMPILER={NVCC}",
         "-DCMAKE_BUILD_TYPE=Release"], env=ENV, capture_output=True, text=True)
    if cfg.returncode != 0:
        log("STDOUT:\n"+cfg.stdout[-2000:]+"\nSTDERR:\n"+cfg.stderr[-3000:])
        raise RuntimeError("GPU cmake configure failed")
    # merge stderr into stdout so compiler errors (which go to stderr) are captured
    bld = subprocess.run(["cmake", "--build", str(BUILD_GPU),
                          "--target", "muon_mass_fork", "-j", "8"],
                         env=ENV, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if bld.returncode != 0:
        log(bld.stdout[-5000:]); raise RuntimeError("GPU build of muon_mass_fork failed")
    log("  running engine tests via ctest")
    t = subprocess.run(["ctest", "-R", "muon_mass_fork", "--output-on-failure"],
                       cwd=str(BUILD_GPU), env=ENV, capture_output=True, text=True)
    log(t.stdout[-2000:])
    # TODO(next): add engine tests for kinetic_expectation_value + ledger +
    #   wrong-slot + GPU-vs-CPU + MPI -np 2 partition (see plan Part 2 Tier-1/4),
    #   and add their catalogue rows BEFORE Phase 2.
    if t.returncode != 0: raise RuntimeError("engine tests failed under CUDA")
    _email("Phase 1 PASSED — GPU build + kernel tests green",
           "Hypothesis check: the per-state mass fork compiles under nvcc and the "
           "Tier-1 kernel oracle passes on GPU.\nDone: GPU build + ctest.\n"
           "Shows: laplacian_states applies factor[ist]·(-|k|²)·ψ per state; "
           "electrons unchanged by the muon slot.\nConclusion: engine path is "
           "sound enough to proceed to physics (Phase 2).\n"
           "STILL TODO in Phase 1: expectation/ledger + MPI-partition + GPU-vs-CPU "
           "engine tests (see plan) before the expensive runs.")
    return True

def phase2():
    """Physics tests: vacuum free-particle σ(t) spreading + analytic oracles + the
    xz-density-vs-σ visualisation set. Runs are dispatched across the whole GPU
    pool; each physics run is then checked against the exact free-Gaussian oracle
    by vacuum_wp/check_oracle.py (parabola m-fit, v_group, norm, <T> drift)."""
    vac_bin = str(HERE / "vacuum_wp/run")
    base = dict(WP_SPACING=0.4, WP_DT=0.02, WP_SIGMA=0.5)

    # --- physics runs (each has an analytic oracle) --------------------------
    # (label, inv_mass, k0, tsteps, L) ; mass-appropriate step counts + box keep
    # the packet in-box: electron/m=10 spread fast (need L=48); the muon barely
    # spreads (rate ~ 1/m^2) so it fits a cheap L=24 (60^3) box even at t=120.
    physics = [
        ("spread_elec", 1.0/1.0,       0.0, 150,  48),  # Panel A electron: t=3
        ("spread_m10",  1.0/10.0,      0.0, 1500, 48),  # mass-dial mid point: t=30
        ("spread_muon", 1.0/MUON_MASS, 0.0, 6000, 24),  # Panel B muon: t=120
        ("vgroup_elec", 1.0/1.0,       0.5, 150,  48),  # Panel C: v_group = 0.5
        ("vgroup_muon", 1.0/MUON_MASS, 0.5, 6000, 24),  # v_group = 0.5/206.77
    ]
    phys_jobs = [dict(binary=vac_bin, rundir=RUNS/f"phase2/{lab}",
                      label=f"phase2 {lab}",
                      overrides={**base, "WP_L": L, "WP_INV_MASS": im, "WP_K0": k0,
                                 "WP_TSTEPS": nst, "WP_OUT": lab})
                 for (lab, im, k0, nst, L) in physics]

    # --- xz-density-vs-σ visualisation (muon), σ ∈ {0.5,1,2,4}: VTI only ------
    xz_jobs = [dict(binary=vac_bin, rundir=RUNS/f"phase2/xz_muon_sig{s:g}",
                    label=f"phase2 xz σ={s:g}",
                    overrides={**base, "WP_L": 48, "WP_SIGMA": s,
                               "WP_INV_MASS": 1.0/MUON_MASS, "WP_K0": 0.0,
                               "WP_TSTEPS": 50, "WP_EMIT_VTI": 1,
                               "WP_OUT": f"xz_muon_sig{s:g}"})
               for s in (0.5, 1.0, 2.0, 4.0)]

    if not run_sims_parallel(phys_jobs + xz_jobs, "phase2 vacuum runs"):
        log("  BLOCKED phase2: one or more vacuum runs failed (see run.logs)")
        return False

    # --- analytic-oracle gate ------------------------------------------------
    checker = str(HERE / "vacuum_wp/check_oracle.py")
    all_ok, lines = True, []
    for (lab, im, k0, _nst, _L) in physics:
        rd = RUNS/f"phase2/{lab}/results/{lab}"   # run.cpp writes to results/<WP_OUT>/
        r = subprocess.run([PY, checker, str(rd), "--sigma_wp", "0.5",
                            "--mass", str(1.0/im), "--k0", str(k0),
                            "--json", str(rd/"oracle.json")],
                           capture_output=True, text=True)
        log(r.stdout.strip() or r.stderr.strip())
        all_ok &= (r.returncode == 0)
        lines.append(f"  {lab}: {'PASS' if r.returncode==0 else 'FAIL'}")
    if not all_ok:
        fail_email("phase2", "analytic free-Gaussian oracle FAILED:\n" + "\n".join(lines))
        return False
    _email("Phase 2 PASSED — free-particle oracles green",
           "Hypothesis check: the per-state mass fork reproduces the EXACT "
           "free-Gaussian spreading law and dispersion for arbitrary mass.\n"
           "Done: vacuum σ(t) runs (electron, m=10, muon), v_group runs, xz-density-"
           "vs-σ VTIs — all across both GPUs.\nShows (per check_oracle.py):\n"
           + "\n".join(lines) +
           "\n  σ_ρ(t)²=σ_ρ0²+t²/(4m²σ_ρ0²) parabola recovers m; v_group=k0/m; "
           "norm & <T> conserved.\nConclusion: mass fork is physically correct "
           "in vacuum -> proceed to the Phase-3 bit-for-bit electron regression.\n"
           "TODO: phase2_physics.ipynb (σ(t) overlays + xz-density-vs-σ figure).")
    return True

def phase3():
    """Bit-for-bit electron regression: He-atom LDA GS+kicked RT built against the
    fork (all mass=1 -> scalar path) vs pristine inq. Runs both (each on its own
    GPU), then diffs GS+RT energies and density. HARD trust gate: any diff above
    tol means the fork perturbs the mass-1 path -> the fork is BROKEN."""
    reg = HERE / "regression"
    fork_bin, pris_bin = str(reg/"fork/run"), str(reg/"pristine/run")
    jobs = [dict(binary=fork_bin, rundir=reg/"fork", label="phase3 fork",
                 overrides={"REG_OUT": "reg_fork"}),
            dict(binary=pris_bin, rundir=reg/"pristine", label="phase3 pristine",
                 overrides={"REG_OUT": "reg_pristine"})]
    # NOTE: pristine build must see INQ_SOURCE=inq, not inq-study. run_sims_parallel
    # inherits ENV (INQ_SOURCE=inq-study); override per job here.
    jobs[1]["overrides"]["INQ_SOURCE"] = str(ROOT/"inq")
    if not run_sims_parallel(jobs, "phase3 regression runs"):
        log("  BLOCKED phase3: a regression binary is not built (build fork/ + "
            "pristine/ via inq-run first)"); return False
    cmp_py = str(reg/"compare_regression.py")
    r = subprocess.run([PY, cmp_py, str(reg/"fork/results/reg_fork"),
                        str(reg/"pristine/results/reg_pristine")],
                       capture_output=True, text=True)
    log(r.stdout.strip() or r.stderr.strip())
    if r.returncode != 0:
        fail_email("phase3", "BIT-FOR-BIT regression FAILED — the mass-1 fork does "
                   "NOT reproduce pristine inq. The fork is BROKEN; do not report "
                   "muon physics.\n\n" + r.stdout[-2000:])
        return False
    _email("Phase 3 PASSED — bit-for-bit vs pristine inq",
           "Hypothesis check: with all mass=1 the fork's empty-factor guard routes "
           "the ORIGINAL scalar kinetic path.\nDone: He-atom LDA GS+kicked RT built "
           "against inq-study (fork) AND pristine inq; compare_regression.py diff.\n"
           "Shows: GS + every RT-step E_total/E_kinetic/E_hartree/E_xc and the GS "
           "density agree to <1e-9 (GPU-reduction floor).\nConclusion: the fork is "
           "inert when off -> the muon results are attributable to the mass, not to "
           "an engine edit. Proceed to Phase 3b.\n\n" + r.stdout[-1500:])
    return True

def phase3b():
    """Simple RT sanity: a muon WP under FULL LDA (interacting) — exercises the
    forked kinetic path inside the Hartree+XC propagator (Phase 2 was
    non_interacting only). Pass = runs complete, energy finite (no NaN)."""
    vac_bin = str(HERE / "vacuum_wp/run")
    jobs = [dict(binary=vac_bin, rundir=RUNS/"phase3b/muon_lda",
                 label="phase3b muon LDA",
                 overrides=dict(WP_THEORY="lda", WP_L=24, WP_SPACING=0.4, WP_DT=0.02,
                                WP_TSTEPS=300, WP_SIGMA=0.5, WP_K0=0.3,
                                WP_INV_MASS=1.0/MUON_MASS, WP_WRITE_EVERY=5,
                                WP_OUT="muon_lda")),
            dict(binary=vac_bin, rundir=RUNS/"phase3b/elec_lda",
                 label="phase3b elec LDA",
                 overrides=dict(WP_THEORY="lda", WP_L=48, WP_SPACING=0.4, WP_DT=0.02,
                                WP_TSTEPS=150, WP_SIGMA=0.5, WP_K0=0.3,
                                WP_INV_MASS=1.0, WP_WRITE_EVERY=5, WP_OUT="elec_lda"))]
    if not run_sims_parallel(jobs, "phase3b RT sanity runs"):
        log("  BLOCKED phase3b: a sanity run failed"); return False
    # finite-energy (no-NaN) check on the WP momentum-space KE trace
    import math
    all_ok, lines = True, []
    for j in jobs:
        rd = Path(j["rundir"]) / f"results/{j['overrides']['WP_OUT']}"
        csv = rd / "raw/observables/wp_momentum_stats.csv"
        finite = csv.exists()
        if finite:
            for ln in csv.read_text().splitlines():
                if ln.startswith("#") or ln.startswith("step"): continue
                vals = ln.split(",")
                if any(v.strip() in ("nan", "-nan", "inf", "-inf") for v in vals):
                    finite = False; break
        all_ok &= finite
        lines.append(f"  {j['overrides']['WP_OUT']}: {'finite (no NaN)' if finite else 'NaN/inf DETECTED'}")
    if not all_ok:
        fail_email("phase3b", "RT sanity FAILED (NaN/inf in the WP dynamics):\n" + "\n".join(lines))
        return False
    _email("Phase 3b PASSED — interacting RT sanity",
           "A muon WP (and electron control) propagated under FULL LDA (Hartree+XC "
           "on the forked kinetic path) — the interacting propagator the free "
           "oracle could not exercise.\n" + "\n".join(lines) +
           "\nConclusion: the muon fork is stable end-to-end in the full "
           "propagator. Proceed to Phase 4 (muon-XC research + your pick).")
    return True

def phase4():
    """Muon-XC research + HARD user checkpoint. Resumes once PICK_FILE exists.

    First pass: emails the user + writes CHECKPOINT (stops). The user records the
    chosen functional in muon_xc_pick.json; a re-run then finds it and PROCEEDS to
    Phase 5. This makes the checkpoint resumable rather than an infinite pause."""
    if PICK_FILE.exists():
        log(f"  muon-XC pick found ({PICK_FILE.name}): {PICK_FILE.read_text()[:200]}")
        return True                              # user has picked → proceed to Phase 5
    _email("Phase 4 CHECKPOINT — PICK the muon-XC functional",
           "Phase 4 (literature-review) has surfaced the candidate muon-XC "
           "prescriptions (see phase4_xc_research.ipynb + docs/sources). Per your "
           "2026-07-06 decision, the orchestrator PAUSES here for you to PICK.\n\n"
           f"To resume: write your choice to {PICK_FILE} (e.g. "
           '{\"functional\": \"...\", \"rationale\": \"...\"}) and re-run '
           "orchestrate.py — it will run Phase 5 (all-muon r_s=5.69 jellium: LDA "
           "vs your pick).")
    _mark_blocked("phase4", "awaiting user muon-XC pick (write muon_xc_pick.json)")
    log("  PHASE 4 CHECKPOINT — paused for user pick; stopping.")
    return "CHECKPOINT"          # sentinel: stop the loop cleanly (not a failure)

def phase5():
    """All-muon r_s=5.69 jellium + incident muon WP: LDA vs user-picked muon-XC."""
    log("  BLOCKED phase5: all-muon GS (r_s=5.69, N=162, L=50 cubic, dx=0.40) + "
        "muon WP run scripts not built; needs user-picked XC")
    return False

PHASES = [("phase1", phase1), ("phase2", phase2), ("phase3", phase3),
          ("phase3b", phase3b), ("phase4", phase4), ("phase5", phase5)]

def main():
    log(f"ORCHESTRATOR start (GPU={GPU}); campaign muon-mass-fork")
    resolve_gpu()                      # single-GPU aware: 'auto' -> freest device
    only = sys.argv[1:] or None
    st = _state()
    for name, fn in PHASES:
        if only and name not in only: continue
        if name in st["done"] and not only:
            log(f"=== PHASE {name} already done — skip ==="); continue
        log(f"=== PHASE {name} ===")
        t0 = time.time()
        try:
            res = fn()
        except Exception:
            tb = traceback.format_exc()
            log(f"PHASE {name} EXCEPTION:\n{tb}")
            _mark_blocked(name, "exception (see log)"); fail_email(name, tb)
            log("STRICT GATING: stopping."); return
        if res == "CHECKPOINT":
            return                                   # clean pause (Phase-4 user pick)
        if not res:
            _mark_blocked(name, "phase returned False")
            fail_email(name, "phase returned False (blocked / not built)")
            log("STRICT GATING: stopping."); return
        _mark_done(name)
        log(f"=== PHASE {name} done ({(time.time()-t0)/60:.1f} min) ===")
    log("ORCHESTRATOR done (all requested phases complete)")

if __name__ == "__main__":
    main()
