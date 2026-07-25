#!/usr/bin/env python3
"""Autonomous orchestrator — classical stopping baseline for the localised jellium slab.

Campaign: docs/campaigns/localised_jellium/classical-stopping-baseline.md
Handover: docs/handovers/classical-stopping-baseline-localised-jellium.md

Runs the two phases CONCURRENTLY, one per GPU:
  Phase 1 (GPU 0) — Ehrenfest, light electron (LJ_CONST_V=0), 2000 steps.
  Phase 2 (GPU 1) — prescribed constant velocity (LJ_CONST_V=1), 1034 steps.

Both twin the WP run p5_wp_v1p3 (σ_WP=0.5, v=1.3, cell 50x50x90, per-3, N=82),
CAP-free. Python (not bash) per the campaigns skill: structured logging,
IDEMPOTENT RESUME (skip runs whose run_summary shows run_completed=true; else
LJ_RESUME=1 from the last checkpoint), per-phase try/except + full-traceback
failure email, per-phase success quick-look email with an energy/projectile plot.

Headless launch (survives disconnect):
    cd .../scripts/classical_slab_stopping
    nohup ../../../../../venv/bin/python3 orchestrate.py > orchestrate.log 2>&1 &
"""
from __future__ import annotations
import os, subprocess, sys, threading, traceback
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
SWEEP = ROOT / "ResearchProject/systems/localised_jellium/scripts/classical_slab_stopping"
BIN = SWEEP / "run"    # inq-run builds the binary as ./run in the sweep dir
RESULTS = SWEEP / "results"
GS = ROOT / "ResearchProject/systems/localised_jellium/shared_gs/slab_n82_L50x50x90"
PY = str(ROOT / "venv/bin/python3")
TO = "chiddukanna@gmail.com"
FAMILY = "[classical-slab-stopping]"

# shared physical parameters — EXACT twin of p5_wp_v1p3
COMMON = dict(
    LJ_LX="50", LJ_LY="50", LJ_LZ="90", LJ_HALF="12.5", LJ_N="82",
    LJ_PERIODICITY="3", LJ_SPACING="0.5", LJ_SIGMA="0.5",
    LJ_LAUNCH_Z="-23.75", LJ_K0="1.3", LJ_MASS="1.0", LJ_DT="0.04",
    LJ_GS_DIR=str(GS),
)
PHASES = [
    dict(name="p1_ehrenfest_v1p3", gpu="0", const_v="0", n_steps="2000", save_every="7",
         desc="Ehrenfest light electron (decelerating); S(v0)=initial-drag slope"),
    dict(name="p2_constv_v1p3",   gpu="1", const_v="1", n_steps="1034", save_every="4",
         desc="prescribed constant velocity; S=dE_deposited/L_slab, stop center at +30"),
]

def log(msg): print(f"[{datetime.now():%F %T}] {msg}", flush=True)

def _summary_done(rundir: Path) -> bool:
    s = rundir / "run_summary.txt"
    return s.exists() and "run_completed = true" in s.read_text()

def _has_checkpoint(rundir: Path) -> bool:
    return (rundir / "rt_state.txt").exists() and (rundir / "checkpoint").exists()

def run_phase(ph: dict):
    name, out = ph["name"], RESULTS / ph["name"]
    logf = SWEEP / f"run_{name}.log"
    if _summary_done(out):
        log(f"{name}: already complete (run_completed=true) — skipping.")
        return
    resume = "1" if _has_checkpoint(out) else "0"
    env = {**os.environ, **COMMON,
           "INQ_SHARE_PATH": str(ROOT / "inq/install/share"),
           "PSEUDOPOD_SHARE_PATH": str(ROOT / "inq/install/share/pseudopod"),
           "TMPDIR": str(ROOT / ".build_tmp"),
           "CUDA_VISIBLE_DEVICES": ph["gpu"],
           "LJ_CONST_V": ph["const_v"], "LJ_N_STEPS": ph["n_steps"],
           "LJ_SAVE_EVERY": ph["save_every"], "LJ_OUT": name, "LJ_RESUME": resume}
    log(f"{name}: launching on GPU {ph['gpu']} (steps={ph['n_steps']}, resume={resume}) -> {logf.name}")
    with open(logf, "a") as lf:
        rc = subprocess.call([str(BIN)], cwd=str(SWEEP), env=env, stdout=lf, stderr=subprocess.STDOUT)
    if rc != 0 or not _summary_done(out):
        raise RuntimeError(f"{name}: binary exited rc={rc} or no run_completed=true (see {logf})")
    log(f"{name}: DONE.")
    _quicklook_email(ph, out)

def _quicklook_email(ph: dict, out: Path):
    """Per-phase success email: hypothesis reminder + a quick energy/projectile plot."""
    try:
        png = _make_quicklook_plot(ph, out)
        sys.path.insert(0, str(ROOT / "inq-stack/python"))
        from inqview.email import send_run_email
        body = (
            f"Campaign: classical stopping baseline for the localised jellium slab.\n"
            f"Hypothesis: the matched localised-slab classical baseline sets S_classical(v=1.3)\n"
            f"as the benchmark the WP quantum stopping (p5_wp_v1p3, S~2.4 eV/Bohr upper bound)\n"
            f"must be judged against.\n\n"
            f"Phase {ph['name']} DONE ({ph['desc']}).\n"
            f"Twin of p5_wp_v1p3: cell 50x50x90 per-3, N=82, sigma_WP=0.5, v=1.3, CAP-free.\n"
            f"The plot shows the total-energy ledger and the projectile track (proj_z, proj_vz,\n"
            f"KE) over the run. Full S extraction + twin decomposition follow in post-processing.\n"
        )
        send_run_email(f"{FAMILY} {ph['name']} complete", body,
                       attachments=[str(png)] if png else None, to=TO)
        log(f"{ph['name']}: quick-look email sent.")
    except Exception:
        log(f"{ph['name']}: quick-look email FAILED:\n{traceback.format_exc()}")

def _make_quicklook_plot(ph: dict, out: Path):
    import numpy as np, matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    obs = out / "raw/observables/observables.csv"
    prj = out / "raw/observables/projectile.csv"
    if not obs.exists(): return None
    o = np.genfromtxt(obs, delimiter=",", names=True)
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
    HA = 27.211386
    ax[0].plot(o["time_au"], (o["energy_total"] - o["energy_total"][0]) * HA)
    ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("dE_total (eV)"); ax[0].set_title("energy ledger")
    if prj.exists():
        p = np.genfromtxt(prj, delimiter=",", names=True)
        ax[1].plot(p["time_au"], p["proj_z"]); ax[1].axhline(-12.5, ls=":", c="k"); ax[1].axhline(12.5, ls=":", c="k")
        ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("proj_z (Bohr)"); ax[1].set_title("trajectory (slab dashed)")
        ax[2].plot(p["time_au"], p["proj_vz"]); ax[2].set_xlabel("t (a.u.)"); ax[2].set_ylabel("proj_vz (a.u.)")
        ax[2].set_title("velocity")
    fig.tight_layout(); png = out / "quicklook.png"; fig.savefig(png, dpi=110); plt.close(fig)
    return png

def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    if not BIN.exists():
        log(f"FATAL: binary not built at {BIN} — build via inq-run first."); sys.exit(2)
    log(f"orchestrator start — launching {len(PHASES)} phases concurrently on GPUs "
        f"{[p['gpu'] for p in PHASES]}")
    threads, errors = [], {}
    def _wrap(ph):
        try: run_phase(ph)
        except Exception:
            errors[ph["name"]] = traceback.format_exc()
            log(f"{ph['name']}: FAILED:\n{errors[ph['name']]}")
            _failure_email(ph, errors[ph["name"]])
    for ph in PHASES:
        t = threading.Thread(target=_wrap, args=(ph,)); t.start(); threads.append(t)
    for t in threads: t.join()
    if errors:
        log(f"orchestrator finished WITH ERRORS in: {list(errors)}"); sys.exit(1)
    log("orchestrator: both phases complete. Post-processing is the next stage "
        "(analyse.py + run-notebooks + twin-decompose vs p5_wp_v1p3 + comparison figure).")

def _failure_email(ph, tb):
    try:
        sys.path.insert(0, str(ROOT / "inq-stack/python"))
        from inqview.email import send_run_email
        send_run_email(f"{FAMILY} {ph['name']} FAILED",
                       f"Phase {ph['name']} raised:\n\n{tb}\n\nOther phases continue.", to=TO)
    except Exception:
        log(f"failure email itself failed:\n{traceback.format_exc()}")

if __name__ == "__main__":
    main()
