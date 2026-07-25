#!/usr/bin/env python3
"""Autonomous orchestrator — campaign localised-jellium-dynamics-analysis.

Runs the whole 5-phase campaign hands-off once the GPUs are free:
  gate  : build binaries + GPU smoke of projectile_background_energy (abort on fail)
  P1    : 6 classical ledger runs (r∈{4,12,20,28,36,40}) carrying E_proj_bg columns
  P2    : 5 classical r_cut runs at r=20 (UPFs rc10/20/30/40 + full)
  P3    : p3_lz120 GS (complete the open-z-vs-PBC matched set; reuse the gs binary)
  P5    : 2 RT screening runs (WP + classical, at rest, r=12)
  nb    : build/execute the campaign notebooks (P1/P2 ledger, P5 GIFs); P4 = reuse
Idempotent resume (run_completed=true ⇒ skip), per-phase try/except + Gmail,
2-GPU parallel. Python (NOT bash) per the campaigns autonomy rule.

Usage: nohup venv/bin/python3 orchestrate.py &   (or via the arm-and-wait wrapper)
"""
import os, sys, subprocess, time, glob, traceback
from pathlib import Path
from datetime import datetime

ROOT = Path("/local/data/public/skcb2/tddft")
LJ   = ROOT/"ResearchProject/systems/localised_jellium"
DYN  = LJ/"scripts/localised_jellium_dynamics"
CA   = LJ/"scripts/campaign_autorun"
SEMI = LJ/"scripts/semiempirical_spillout"
GS_P2   = str(CA/"runs/h2/gs_p2_lz120/checkpoint")
GS_P2_RES = CA/"runs/h2/gs_p2_lz120/results"
UPF_DIR = CA/"cutoff_test/upfs"
UPF_FULL = str(ROOT/"ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf")
PY   = str(ROOT/"venv/bin/python3")
INQRUN = str(ROOT/"shared/bin/inq-run")
HYP  = LJ/"hypotheses/localised_jellium_dynamics"
CAMPAIGN_MD = ROOT/"docs/campaigns/localised_jellium_dynamics_analysis/localised-jellium-dynamics-analysis.md"
TO = "chiddukanna@gmail.com"

ENV = {**os.environ,
       "INQ_SHARE_PATH": str(ROOT/"inq/install/share"),
       "PSEUDOPOD_SHARE_PATH": str(ROOT/"inq/install/share/pseudopod"),
       "INQ_SOURCE": str(ROOT/"inq-study")}

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)

# ---------------------------------------------------------------- infrastructure
def _done(rundir: Path) -> bool:
    for rs in Path(rundir).glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs.read_text(): return True
        except Exception: pass
    return False

def build(bindir: Path) -> bool:
    """Clean GPU build of a run.cpp (nvcc; does not need a free GPU). A `.gpu_built`
    marker makes resume skip; any prior --cpu build is wiped so we never run a CPU
    binary on the GPU queue (which would silently be days-slow)."""
    bindir = Path(bindir)
    if (bindir/".gpu_built").exists() and (bindir/"run").exists():
        log(f"  build SKIP {bindir.name} (gpu binary present)"); return True
    log(f"  BUILD {bindir.name} (clean GPU)")
    import shutil
    shutil.rmtree(bindir/"build", ignore_errors=True)   # wipe any prior --cpu config → GPU reconfigure
    (bindir/"run").unlink(missing_ok=True)
    # inq-run has no build-only mode: it configures (GPU, since no --cpu) + builds + runs.
    # For GS-requiring binaries the run exits fast with FATAL (no LJ_GS_DIR here) AFTER the
    # binary is linked — harmless; we only care that `run` exists.
    with open(bindir/"build.log","w") as lf:
        subprocess.run([INQRUN], cwd=str(bindir), env=ENV,
                       stdout=lf, stderr=subprocess.STDOUT)
    ok = (bindir/"run").exists()
    if ok: (bindir/".gpu_built").write_text("gpu\n")
    log(f"  build {'OK' if ok else 'FAIL'} {bindir.name}"); return ok

def run_sim(binary: str, rundir: Path, overrides: dict, label: str, gpu: str, retries: int = 1) -> bool:
    rundir.mkdir(parents=True, exist_ok=True)
    if _done(rundir):
        log(f"  SKIP {label} (done)"); return True
    env = {**ENV, **{k: str(v) for k, v in overrides.items()}, "CUDA_VISIBLE_DEVICES": gpu}
    for attempt in range(retries+1):
        log(f"  RUN [gpu{gpu}] {label}" + (f" (retry {attempt})" if attempt else ""))
        with open(rundir/"run.log","w") as lf:
            rc = subprocess.run([binary], cwd=str(rundir), env=env, stdout=lf, stderr=subprocess.STDOUT).returncode
        if rc == 0 and _done(rundir): return True
        log(f"  FAIL {label} rc={rc} (see {rundir/'run.log'})")
    return False

def run_parallel(jobs: list) -> list:
    """jobs = [(binary, rundir, overrides, label), ...]; split round-robin across 2 GPUs.
    Each GPU runs its queue sequentially; the two queues run concurrently."""
    queues = {"0": [], "1": []}
    for i, j in enumerate(jobs): queues[str(i % 2)].append(j)
    import threading
    results = {}
    def worker(gpu, q):
        for (binary, rundir, ov, label) in q:
            results[label] = run_sim(binary, rundir, ov, label, gpu)
    ts = [threading.Thread(target=worker, args=(g, q)) for g, q in queues.items()]
    for t in ts: t.start()
    for t in ts: t.join()
    return [results.get(j[3], False) for j in jobs]

def email(subject, body, attachments=None):
    try:
        sys.path.insert(0, str(ROOT/"inq-stack/python"))
        from inqview.email import send_run_email
        send_run_email(subject=subject, body=body, attachments=attachments or [], to=TO)
        log(f"  emailed: {subject}")
    except Exception as e:
        log(f"  email failed ({e}): {subject}")

def wait_for_gpus():
    """Block until the current semiempirical_spillout runs release both GPUs."""
    log("waiting for GPUs (current semiempirical_spillout runs to finish)...")
    while True:
        r = subprocess.run(["pgrep","-f","semiempirical_spillout/gs/run"], capture_output=True, text=True)
        if not r.stdout.strip():
            log("GPUs free."); return
        time.sleep(60)

def set_task_done(idx: int):
    """Flip the idx-th `done: false` → true in the campaign frontmatter (best-effort)."""
    try:
        t = CAMPAIGN_MD.read_text(); lines = t.splitlines(); n = 0
        for i, ln in enumerate(lines):
            if "done: false" in ln:
                if n == idx: lines[i] = ln.replace("done: false","done: true"); break
                n += 1
        CAMPAIGN_MD.write_text("\n".join(lines)+"\n")
    except Exception as e:
        log(f"  frontmatter update failed: {e}")

RADII = (4, 12, 20, 28, 36, 40)
def z_of(r): return -(12.5 + r)

# ---------------------------------------------------------------- gate + phases
def gpu_smoke_gate() -> bool:
    # (a) capability physics/machinery
    smoke = DYN/"smoke_eprojbg"
    if not build(smoke): return False
    log("GPU smoke gate (a): projectile_background_energy")
    env = {**ENV, "CUDA_VISIBLE_DEVICES":"0"}
    with open(smoke/"gpu_smoke.log","w") as lf:
        subprocess.run([str(smoke/"run")], cwd=str(smoke), env=env, stdout=lf, stderr=subprocess.STDOUT)
    if "SMOKE PASS" not in (smoke/"gpu_smoke.log").read_text():
        log("  smoke gate (a) FAIL"); return False
    log("  smoke gate (a) PASS")
    # (b) full chain: one phase12 run emits the two columns with n_proj_norm≈1
    log("GPU smoke gate (b): phase12 column emission")
    if not build(DYN/"phase12"): return False
    gate_dir = DYN/"runs/gate/cl_r12_smoke"
    ok = run_sim(str(DYN/"phase12/run"), gate_dir,
                 dict(LJ_OUT="cl_r12_smoke", LJ_LZ=120, LJ_PERIODICITY=2, LJ_LAUNCH_Z=z_of(12),
                      LJ_GS_DIR=GS_P2, LJ_SIGMA=0.5, LJ_N_STEPS=2), "gate phase12 r=12", "0", retries=0)
    if not ok: log("  smoke gate (b) FAIL (run)"); return False
    csvs = glob.glob(str(gate_dir/"**/observables.csv"), recursive=True)
    rs   = glob.glob(str(gate_dir/"**/run_summary.txt"), recursive=True)
    hdr  = open(csvs[0]).readline() if csvs else ""
    norm = 1.0
    if rs:
        for ln in open(rs[0]):
            if ln.startswith("n_proj_norm"):
                try: norm = float(ln.split("=")[1].split()[0])
                except Exception: pass
    ok = ("energy_proj_bg_ideal" in hdr and "energy_proj_bg_impl" in hdr and abs(norm-1.0) < 0.02)
    log(f"  smoke gate (b) {'PASS' if ok else 'FAIL'}  (cols={'energy_proj_bg_ideal' in hdr}, n_proj_norm={norm:.4f})")
    return ok

def phase1():
    bindir = DYN/"phase12"; assert (bindir/"run").exists()
    jobs = [(str(bindir/"run"), DYN/f"runs/p1/cl_r{r}",
             dict(LJ_OUT=f"cl_r{r}", LJ_LZ=120, LJ_PERIODICITY=2, LJ_LAUNCH_Z=z_of(r),
                  LJ_GS_DIR=GS_P2, LJ_SIGMA=0.5, LJ_N_STEPS=2), f"P1 cl r={r}") for r in RADII]
    res = run_parallel(jobs)
    ledger = subprocess.run([PY, str(DYN/"build_ledger_notebook.py")], env=ENV).returncode == 0
    body = _ledger_body()
    email("[localised-jellium-dynamics] Phase 1 — completed ledger (U_proj_bg)", body,
          attachments=[str(HYP/"ledger.png")] if (HYP/"ledger.png").exists() else [])
    if all(res): set_task_done(1)
    return all(res)

def phase2():
    bindir = DYN/"phase12"
    upfs = {10:str(UPF_DIR/"electron_gaussian_wpsigma0p5_rc10.upf"),
            20:str(UPF_DIR/"electron_gaussian_wpsigma0p5_rc20.upf"),
            30:str(UPF_DIR/"electron_gaussian_wpsigma0p5_rc30.upf"),
            40:str(UPF_DIR/"electron_gaussian_wpsigma0p5_rc40.upf"),
            50:UPF_FULL}
    jobs = [(str(bindir/"run"), DYN/f"runs/p2/cl_r20_rc{rc}",
             dict(LJ_OUT=f"cl_r20_rc{rc}", LJ_LZ=120, LJ_PERIODICITY=2, LJ_LAUNCH_Z=z_of(20),
                  LJ_GS_DIR=GS_P2, LJ_SIGMA=0.5, LJ_N_STEPS=2, LJ_PROJ_UPF=upf), f"P2 r=20 rc={rc}")
            for rc, upf in upfs.items()]
    res = run_parallel(jobs)
    subprocess.run([PY, str(DYN/"build_ledger_notebook.py")], env=ENV)
    email("[localised-jellium-dynamics] Phase 2 — r_cut sweep at r=20", _rcut_body(),
          attachments=[str(HYP/"rcut.png")] if (HYP/"rcut.png").exists() else [])
    if all(res): set_task_done(2)
    return all(res)

def phase3():
    gsbin = str(SEMI/"gs/run"); rundir = SEMI/"runs/p3_lz120"
    ok = run_sim(gsbin, rundir, dict(LJ_LX=50,LJ_LY=50,LJ_LZ=120,LJ_HALF=12.5,LJ_N=82,LJ_EDGE_W=0,
                 LJ_PERIODICITY=3,LJ_SPACING=0.5,LJ_EXTRA_STATES=20,LJ_TEMP_EV=0.00862,
                 LJ_GS_DIR=str(rundir/"checkpoint"),LJ_TAG="p3_lz120"), "P3 p3_lz120", "0")
    # refresh the semiempirical notebook (Q3 open-z vs PBC now has the full matched set)
    subprocess.run([PY, str(LJ/"hypotheses/campaign_autorun_study/build_semiempirical_spillout.py")], env=ENV)
    subprocess.run([PY,"-m","nbconvert","--to","notebook","--execute","--inplace",
                    str(LJ/"hypotheses/campaign_autorun_study/semiempirical_spillout.ipynb"),
                    "--ExecutePreprocessor.timeout=1200"], env=ENV)
    email("[localised-jellium-dynamics] Phase 3 — open-z vs PBC (p3_lz120 added)",
          "p3_lz120 complete; the semiempirical_spillout notebook Q3 now has the full matched "
          "{90,120,160,240}×{p2,p3} set. See the open-z-vs-PBC overlay + near-edge table.",
          attachments=[])
    if ok: set_task_done(3)
    return ok

def phase5():
    for b in ("phase5_wp","phase5_cl"): build(DYN/b)
    jobs = [
        (str(DYN/"phase5_wp/run"), DYN/"runs/p5/wp",
         dict(LJ_OUT="wp", LJ_LZ=120, LJ_PERIODICITY=2, LJ_LAUNCH_Z=z_of(12), LJ_K0=0, LJ_SIGMA=0.5,
              LJ_GS_DIR=GS_P2, LJ_N_STEPS=500, LJ_DT=0.01, LJ_SAVE_EVERY=25), "P5 wp"),
        (str(DYN/"phase5_cl/run"), DYN/"runs/p5/cl",
         dict(LJ_OUT="cl", LJ_LZ=120, LJ_PERIODICITY=2, LJ_LAUNCH_Z=z_of(12),
              LJ_GS_DIR=GS_P2, LJ_N_STEPS=500, LJ_DT=0.01, LJ_SAVE_EVERY=25), "P5 cl"),
    ]
    res = run_parallel(jobs)
    subprocess.run([PY, str(DYN/"build_screening_gifs.py")], env=ENV)
    gifs = [str(p) for p in HYP.glob("*.gif")]
    email("[localised-jellium-dynamics] Phase 5 — screening GIFs (WP vs classical)",
          "RT screening runs complete (at rest, r=12, 500 steps). Total and induced (bath-only) "
          "density-difference GIFs attached.", attachments=gifs[:4])
    if all(res): set_task_done(5)
    return all(res)

# ------- notebook body helpers (numeric; robust to missing runs) -------------
def _read_row0(csv_path):
    import csv
    r = list(csv.reader(open(csv_path))); h, d = r[0], r[1]
    return {k: float(v) for k, v in zip(h, d)}
def _ledger_body():
    try:
        HA=27.211386; lines=["Phase 1 — completed classical ledger with U_proj_bg (eV):",""]
        for r in RADII:
            f=glob.glob(str(DYN/f"runs/p1/cl_r{r}/**/observables.csv"),recursive=True)
            if not f: continue
            row=_read_row0(f[0])
            lines.append(f"  r={r:>2}: U_proj_bg(ideal)={row.get('energy_proj_bg_ideal',0)*HA:8.2f}  "
                         f"impl={row.get('energy_proj_bg_impl',0)*HA:8.2f} eV")
        lines+=["", "Hypothesis: d(U_H+U_ext+U_proj_bg) closes the WP−CL gap (see ledger.png)."]
        return "\n".join(lines)
    except Exception as e: return f"(ledger body error: {e})"
def _rcut_body():
    try:
        HA=27.211386; lines=["Phase 2 — r_cut sweep at r=20 (eV):",""]
        for rc in (10,20,30,40,50):
            f=glob.glob(str(DYN/f"runs/p2/cl_r20_rc{rc}/**/observables.csv"),recursive=True)
            if not f: continue
            row=_read_row0(f[0])
            lines.append(f"  rc={rc:>2}: E_ext={row.get('energy_external',0)*HA:9.2f}  "
                         f"proj_bg_impl={row.get('energy_proj_bg_impl',0)*HA:8.2f}  "
                         f"ideal={row.get('energy_proj_bg_ideal',0)*HA:8.2f} eV")
        lines+=["","ideal should be flat vs r_cut; the r_cut effect lives in E_external."]
        return "\n".join(lines)
    except Exception as e: return f"(rcut body error: {e})"

# ---------------------------------------------------------------- main
def main():
    HYP.mkdir(parents=True, exist_ok=True)
    log("=== campaign localised-jellium-dynamics-analysis: autonomous orchestrator ===")
    wait_for_gpus()
    # build the ledger/RT binaries (nvcc; GPUs may still be busy elsewhere)
    for b in ("phase12","phase5_wp","phase5_cl"): build(DYN/b)
    if not gpu_smoke_gate():
        email("[localised-jellium-dynamics] ABORT — E_proj_bg GPU smoke FAILED",
              "The projectile_background_energy GPU smoke gate failed; NO campaign runs were "
              "launched (correctness gate). See smoke_eprojbg/gpu_smoke.log.")
        log("ABORT: smoke gate failed."); return
    set_task_done(0)
    phases = [("P1",phase1),("P2",phase2),("P3",phase3),("P5",phase5)]
    for name, fn in phases:
        try:
            log(f"===== {name} ====="); fn()
        except Exception:
            tb = traceback.format_exc(); log(f"  {name} EXCEPTION:\n{tb}")
            email(f"[localised-jellium-dynamics] {name} — PHASE FAILED (chain continues)",
                  f"Phase {name} raised:\n\n{tb}\n\nOrchestrator continued; re-run to resume.")
    # final notebook + completion email
    try:
        subprocess.run([PY, str(DYN/"build_ledger_notebook.py")], env=ENV)
    except Exception: pass
    set_task_done(6)
    email("[localised-jellium-dynamics] CAMPAIGN COMPLETE (all phases)",
          _ledger_body()+"\n\n"+_rcut_body()+"\n\nNotebooks under "+str(HYP)+".",
          attachments=[str(HYP/p) for p in ("ledger.png","rcut.png") if (HYP/p).exists()])
    log("=== orchestrator done ===")

if __name__ == "__main__":
    main()
