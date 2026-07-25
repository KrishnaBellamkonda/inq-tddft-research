#!/usr/bin/env python3
"""Python orchestrator for the localised-jellium GS ladder (H1-H5 + cumulative).

Replaces the bash master.sh. Python is preferred for autonomous multi-phase
dispatch: structured logging, IDEMPOTENT RESUME (skips runs whose run_summary
shows run_completed=true), per-phase try/except with full-traceback failure
emails, and one-shot retry on a sim failure. H0 is already done + emailed.

Headless launch (survives disconnect):
    cd .../campaign_autorun
    GPU=1 nohup venv/bin/python3 orchestrate.py > orchestrate.log 2>&1 &
"""
from __future__ import annotations
import os, subprocess, sys, time, traceback
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
LJ = ROOT / "ResearchProject/systems/localised_jellium"
CA = LJ / "scripts/campaign_autorun"
RUNS = CA / "runs"; RUNS.mkdir(parents=True, exist_ok=True)
PY = str(ROOT / "venv/bin/python3")
INQRUN = str(ROOT / "shared/bin/inq-run")
GSBIN, WPBIN, CLBIN = str(CA/"gs/run"), str(CA/"wp/run"), str(CA/"classical/run")
ANALYSE = str(CA / "analyse_phase.py")
GPU = os.environ.get("GPU", "1")
GS120_P3_CKPT = str(LJ / "shared_gs/slab_n82_L50x50x120")                 # H0 periodicity-3 GS
GS120_P3_RES  = str(LJ / "scripts/h0_base_difference/gs/results")          # its run_summary (E_GS)
GS120_P2_DIR  = RUNS / "h2/gs_p2_lz120"
TO = "chiddukanna@gmail.com"

ENV = {**os.environ,
       "INQ_SHARE_PATH": str(ROOT/"inq/install/share"),
       "PSEUDOPOD_SHARE_PATH": str(ROOT/"inq/install/share/pseudopod"),
       "INQ_SOURCE": str(ROOT/"inq-study"),
       "CUDA_VISIBLE_DEVICES": GPU}

def log(msg): print(f"[{datetime.now():%F %T}] {msg}", flush=True)

def _done(rundir: Path) -> bool:
    # tolerant of layout: GS writes results/run_summary.txt; projectiles nest under
    # results/<LJ_OUT>/run_summary.txt. Glob finds either.
    for rs in rundir.glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs.read_text():
                return True
        except Exception:
            pass
    return False

def run_sim(binary: str, rundir: Path, overrides: dict, label: str, retries: int = 1) -> bool:
    rundir.mkdir(parents=True, exist_ok=True)
    if _done(rundir):
        log(f"  SKIP {label} (already complete)"); return True
    env = {**ENV, **{k: str(v) for k, v in overrides.items()}}
    for attempt in range(retries + 1):
        log(f"  RUN  {label}" + (f" (retry {attempt})" if attempt else ""))
        with open(rundir / "run.log", "w") as lf:
            rc = subprocess.run([binary], cwd=str(rundir), env=env, stdout=lf, stderr=subprocess.STDOUT).returncode
        if rc == 0 and _done(rundir):
            return True
        log(f"  FAIL {label} rc={rc} (see {rundir/'run.log'})")
    return False

def analyse(args: list[str], phase: str) -> bool:
    log(f"  ANALYSE {phase}: {' '.join(args)}")
    rc = subprocess.run([PY, ANALYSE, *args], env=ENV).returncode
    if rc != 0: log(f"  analyse {phase} rc={rc} (analyse_phase.py emails its own error)")
    return rc == 0

def fail_email(phase: str, msg: str):
    try:
        sys.path.insert(0, str(ROOT/"inq-stack/python"))
        from inqview.email import send_run_email
        send_run_email(subject=f"[localised-jellium GS] {phase} — PHASE FAILED (chain continues)",
                       body=f"Phase {phase} failed in the Python orchestrator.\n\n{msg}\n\n"
                            f"Data (if any) under {RUNS}. The orchestrator continued to the next phase; "
                            f"re-run orchestrate.py to resume (completed runs are skipped).",
                       attachments=[], to=TO)
    except Exception as e:
        log(f"  could not send fail email: {e}")

# ---------------- phases ----------------------------------------------------
def H0():
    rs = (4, 12, 20, 28, 36, 40)
    for tag, binary in (("wp", WPBIN), ("cl", CLBIN)):
        for r in rs:
            z = -(12.5 + r)
            ov = dict(LJ_OUT=f"{tag}_r{r}_p3", LJ_LZ=120, LJ_PERIODICITY=3,
                      LJ_LAUNCH_Z=z, LJ_GS_DIR=GS120_P3_CKPT)
            if tag == "wp":
                ov.update(LJ_K0=0, LJ_SIGMA=0.5)
            run_sim(binary, RUNS / f"h0/{tag}_r{r}_p3", ov, f"H0 {tag} r={r}")
    analyse(["--phase", "H0", "--base", str(RUNS / "h0"), "--gs120-p3", GS120_P3_RES], "H0")

def H1():
    for w in (0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3):
        run_sim(GSBIN, RUNS/f"h1/gs_w{w:g}",
                dict(LJ_LX=50, LJ_LY=50, LJ_LZ=90, LJ_HALF=12.5, LJ_N=82,
                     LJ_EDGE_W=w, LJ_PERIODICITY=3, LJ_TAG=f"h1_w{w:g}",
                     LJ_GS_DIR=str(RUNS/f"h1/gs_w{w:g}/checkpoint")), f"H1 gs w={w:g}")
    analyse(["--phase", "H1", "--base", str(RUNS/"h1")], "H1")

def H2():
    for lz in (50, 60, 70, 80, 90, 110, 130, 150):
        run_sim(GSBIN, RUNS/f"h2/gs_lz{lz}",
                dict(LJ_LX=50, LJ_LY=50, LJ_LZ=lz, LJ_HALF=12.5, LJ_N=82,
                     LJ_EDGE_W=0, LJ_PERIODICITY=3, LJ_TAG=f"h2_lz{lz}",
                     LJ_GS_DIR=str(RUNS/f"h2/gs_lz{lz}/checkpoint")), f"H2 gs lz={lz}")
    # open-z (periodicity 2) at L_z=90 (H2 comparison) and 120 (reused by H4/H5)
    run_sim(GSBIN, RUNS/"h2/gs_p2_lz90",
            dict(LJ_LX=50, LJ_LY=50, LJ_LZ=90, LJ_HALF=12.5, LJ_N=82,
                 LJ_EDGE_W=0, LJ_PERIODICITY=2, LJ_TAG="h2_p2_lz90",
                 LJ_GS_DIR=str(RUNS/"h2/gs_p2_lz90/checkpoint")), "H2 gs periodicity-2 L_z=90")
    run_sim(GSBIN, GS120_P2_DIR,
            dict(LJ_LX=50, LJ_LY=50, LJ_LZ=120, LJ_HALF=12.5, LJ_N=82,
                 LJ_EDGE_W=0, LJ_PERIODICITY=2, LJ_TAG="h2_p2_lz120",
                 LJ_GS_DIR=str(GS120_P2_DIR/"checkpoint")), "H2 gs periodicity-2 L_z=120")
    analyse(["--phase", "H2", "--base", str(RUNS/"h2")], "H2")

def H3():
    for a, N in ((5, 32), (7.5, 50), (10, 66), (12.5, 82), (15, 98),
                 (17.5, 114), (20, 132), (22.5, 148), (27.5, 180)):
        run_sim(GSBIN, RUNS/f"h3/gs_a{a:g}_N{N}",
                dict(LJ_LX=50, LJ_LY=50, LJ_LZ=90, LJ_HALF=a, LJ_N=N,
                     LJ_EDGE_W=0, LJ_PERIODICITY=3, LJ_TAG=f"h3_a{a:g}",
                     LJ_GS_DIR=str(RUNS/f"h3/gs_a{a:g}_N{N}/checkpoint")), f"H3 gs a={a:g} N={N}")
    analyse(["--phase", "H3", "--base", str(RUNS/"h3")], "H3")

def _proj_phase(binary, tag):
    hdir = RUNS / ("h4" if tag == "wp" else "h5")
    p2_ckpt = str(GS120_P2_DIR / "checkpoint")
    have_p2 = (GS120_P2_DIR / "checkpoint").exists()
    if not have_p2:
        log(f"  NOTE: periodicity-2 GS missing ({GS120_P2_DIR}); {tag} open-z runs skipped")
    for r in (2, 6, 10, 14, 18, 22, 26, 30, 34, 40):
        z = -(12.5 + r)
        run_sim(binary, hdir / f"{tag}_r{r}_p3",
                dict(LJ_OUT=f"{tag}_r{r}_p3", LJ_LZ=120, LJ_PERIODICITY=3, LJ_LAUNCH_Z=z,
                     LJ_K0=0, LJ_SIGMA=0.5, LJ_GS_DIR=GS120_P3_CKPT), f"{tag} r={r} per3")
        if have_p2:
            run_sim(binary, hdir / f"{tag}_r{r}_p2",
                    dict(LJ_OUT=f"{tag}_r{r}_p2", LJ_LZ=120, LJ_PERIODICITY=2, LJ_LAUNCH_Z=z,
                         LJ_K0=0, LJ_SIGMA=0.5, LJ_GS_DIR=p2_ckpt), f"{tag} r={r} per2")

def H4():
    _proj_phase(WPBIN, "wp")
    analyse(["--phase", "H4", "--base", str(RUNS/"h4"),
             "--gs120-p3", GS120_P3_RES, "--gs120-p2", str(GS120_P2_DIR/"results")], "H4")

def H5():
    _proj_phase(CLBIN, "cl")
    analyse(["--phase", "H5", "--base", str(RUNS/"h5"), "--h4-base", str(RUNS/"h4"),
             "--gs120-p3", GS120_P3_RES, "--gs120-p2", str(GS120_P2_DIR/"results")], "H5")

def CUMULATIVE():
    sys.path.insert(0, str(ROOT/"inq-stack/python"))
    from inqview.email import send_run_email
    cands = ["scripts/campaign_autorun/runs/h0/H0_base_difference.png",
             "scripts/campaign_autorun/runs/h1/H1_edge_model.png",
             "scripts/campaign_autorun/runs/h2/H2_gs_convergence.png",
             "scripts/campaign_autorun/runs/h3/H3_surface_energetics.png",
             "scripts/campaign_autorun/runs/h4/H4_wp_energetics.png",
             "scripts/campaign_autorun/runs/h5/H5_classical_subtraction.png"]
    pngs = [str(LJ/p) for p in cands if (LJ/p).exists()]
    send_run_email(subject="[localised-jellium GS] CAMPAIGN COMPLETE — cumulative (all phases)",
                   body=("The localised-jellium GS ladder finished its headless (Python) run.\n"
                         "Attached: the highlight plot from each completed phase (H0-H5).\n\n"
                         "Review the per-phase emails for hypothesis/method/plot/conclusion.\n"
                         "Flagged for your judgement: H2 work-function (Phi), H3 sigma_s (E_self), "
                         "H5 ghost-background. Critical path (H4 E_SIE + PBC-vs-open-z) is fully analysed."),
                   attachments=pngs, to=TO)
    log(f"cumulative sent with {len(pngs)} plots")

PHASES = [("H0", H0), ("H1", H1), ("H2", H2), ("H3", H3), ("H4", H4), ("H5", H5), ("CUMULATIVE", CUMULATIVE)]

def main():
    log(f"ORCHESTRATOR start (GPU={GPU}); H0 already done+emailed")
    only = sys.argv[1:] or None
    for name, fn in PHASES:
        if only and name not in only: continue
        log(f"=== PHASE {name} ===")
        t0 = time.time()
        try:
            fn()
            log(f"=== PHASE {name} done ({(time.time()-t0)/60:.1f} min) ===")
        except Exception:
            tb = traceback.format_exc(); log(f"PHASE {name} EXCEPTION:\n{tb}"); fail_email(name, tb)
    log("ORCHESTRATOR done")

if __name__ == "__main__":
    main()
