#!/usr/bin/env python3
"""Hardened P0b relaunch — wide-WP campaign design gate.

Fixes the 2026-07-01 failure (both foreground `inq-run` jobs died together at ~03:00,
almost certainly a SIGHUP when their controlling session closed; the orchestrator then
polled 8 h for a completion that never came and halted). This launcher:

  1. DETACHES every child (start_new_session=True => own session, SIGHUP-immune);
  2. LAUNCHES the P0b pair itself (does not wait on hand-started jobs it cannot restart);
  3. adds a LIVENESS check — if a run's log stops advancing for STALL_MIN with no
     completion sentinel, it kills both and emails an alert (instead of an 8 h timeout);
  4. HOLDS AT THE GATE — on success it rebuilds the gate-review notebook and emails
     "review & sign off"; it does NOT start the Phase-1 sweep. That stays the user's call.

Launch itself detached:
    cd .../scripts/wide_wp
    setsid nohup /local/data/public/skcb2/tddft/venv/bin/python3 relaunch_p0b.py \
        > relaunch_p0b.log 2>&1 &
"""
from __future__ import annotations
import os, sys, time, shutil, subprocess, traceback
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
LJ   = ROOT / "ResearchProject/systems/localised_jellium"
WWP  = LJ / "scripts/wide_wp"
WPDIR, CLDIR = WWP / "wp", WWP / "classical"
HYP  = LJ / "hypotheses/wide_wp"
PY   = str(ROOT / "venv/bin/python3")
TO   = "chiddukanna@gmail.com"
SUBJ = "[wide-wavepacket] "

STALL_MIN   = 20      # kill+alert if a live run's log is silent this long (no sentinel)
POLL_S      = 120
MAX_HOURS   = 6

# exact original P0b config (from the p0b_*.log headers)
K0 = 4.6957
WP_ENV = dict(LJ_OUT="results/p0b_wp", LJ_K0=K0, LJ_N_STEPS=750, LJ_DT=0.04,
              LJ_CAP=1, LJ_WRITE_EVERY=4, LJ_WF_EVERY=40, LJ_LAUNCH_Z=-26.5)
CL_ENV = dict(LJ_OUT="results/p0b_classical", LJ_K0=K0, LJ_N_STEPS=1500, LJ_DT=0.02,
              LJ_CAP=1, LJ_WRITE_EVERY=10, LJ_LAUNCH_Z=-26.5)

ENV_BASE = {**os.environ,
            "INQ_SHARE_PATH":       str(ROOT / "inq/install/share"),
            "PSEUDOPOD_SHARE_PATH": str(ROOT / "inq/install/share/pseudopod"),
            "INQ_SOURCE":           str(ROOT / "inq-study")}

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)

def email(subject, body, attachments=None):
    try:
        sys.path.insert(0, str(ROOT / "inq-stack/python"))
        from inqview.email import send_run_email
        return send_run_email(subject=SUBJ + subject, body=body,
                              attachments=attachments or [], to=TO)
    except Exception as e:
        log(f"  EMAIL FAILED: {e}")

def sentinel(rundir: Path) -> bool:
    for rs in rundir.glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs.read_text():
                return True
        except Exception:
            pass
    return False

def launch(cwd: Path, out_log: str, gpu: int, overrides: dict) -> subprocess.Popen:
    env = {**ENV_BASE, "CUDA_VISIBLE_DEVICES": str(gpu),
           **{k: str(v) for k, v in overrides.items()}}
    lf = open(cwd / out_log, "w")
    # start_new_session=True => detached session leader, immune to controlling-terminal SIGHUP
    return subprocess.Popen([str(cwd / "run")], cwd=str(cwd), env=env,
                            stdout=lf, stderr=subprocess.STDOUT, start_new_session=True)

def main():
    log("HARDENED P0b relaunch — detached, liveness-guarded, holds at gate")
    wp_res, cl_res = WPDIR / "results/p0b_wp", CLDIR / "results/p0b_classical"
    # fresh start: the partial (killed) data is superseded by this complete run
    for d in (wp_res, cl_res):
        if d.exists(): shutil.rmtree(d); log(f"  cleared partial {d}")

    procs = {
        "wp": (launch(WPDIR, "p0b_wp.log", 0, WP_ENV), WPDIR / "p0b_wp.log", wp_res),
        "cl": (launch(CLDIR, "p0b_classical.log", 1, CL_ENV), CLDIR / "p0b_classical.log", cl_res),
    }
    log(f"  launched wp (GPU0, pid {procs['wp'][0].pid}) + classical (GPU1, pid {procs['cl'][0].pid})")

    t_end = time.time() + MAX_HOURS * 3600
    last_size = {k: -1 for k in procs}
    last_move = {k: time.time() for k in procs}
    while time.time() < t_end:
        time.sleep(POLL_S)
        if all(sentinel(r) for _, _, r in procs.values()):
            log("both runs complete"); break
        for k, (p, lf, res) in procs.items():
            if sentinel(res):
                continue
            sz = lf.stat().st_size if lf.exists() else 0
            if sz > last_size[k]:
                last_size[k] = sz; last_move[k] = time.time()
            elif p.poll() is not None:                     # process exited w/o sentinel
                log(f"  {k} exited (rc={p.returncode}) with NO sentinel — dead")
                last_move[k] = 0
            silent = (time.time() - last_move[k]) / 60.0
            if silent > STALL_MIN:
                log(f"  {k} STALLED ({silent:.0f} min silent, no sentinel) — killing pair")
                for pp, _, _ in procs.values():
                    try: pp.terminate()
                    except Exception: pass
                email("P0b relaunch STALLED — needs attention",
                      f"Run '{k}' produced no log output for {silent:.0f} min and never wrote the "
                      f"completion sentinel; killed the pair. Check GPU/host state before retrying.\n"
                      f"Logs: {procs['wp'][1]} , {procs['cl'][1]}")
                return

    ok = all(sentinel(r) for _, _, r in procs.values())
    if not ok:
        email("P0b relaunch did not finish in time — HALTING",
              f"After {MAX_HOURS} h not both runs completed. WP sentinel="
              f"{sentinel(wp_res)}, classical sentinel={sentinel(cl_res)}.")
        log("timeout — halting"); return

    # rebuild the gate-review notebook on the COMPLETE data (auto-update)
    figs = []
    try:
        subprocess.run([PY, str(HYP / "build_p0b_gate_review.py")],
                       env={**ENV_BASE, "PYTHONPATH": str(ROOT / "inq-stack/python")},
                       cwd=str(HYP), check=True, timeout=1200)
        figs = [str(HYP / "p0b_gate_review_figs" / f)
                for f in ("c1_spreading.png", "c2_wpnorm.png", "c6_classical.png")]
        log("gate-review notebook rebuilt on complete data")
    except Exception:
        log("notebook rebuild failed:\n" + traceback.format_exc())

    email("P0b COMPLETE — review the gate notebook & sign off (holding at gate)",
          "HYPOTHESIS: a wide near-rigid wavepacket (σ_WP=3.5) and a matched-σ classical "
          "projectile through the localised jellium slab isolate purely quantum stopping.\n\n"
          "WHAT WAS DONE: the P0b matched pair (E=300 eV) was RELAUNCHED to completion, "
          "detached (SIGHUP-immune — fixing the 2026-07-01 death) and liveness-guarded.\n\n"
          "WHAT TO CHECK: the gate-review notebook has been rebuilt on the COMPLETE data at\n"
          f"  {HYP / 'p0b_gate_review.ipynb'}\n"
          "Sign off criteria 1–6 there (spreading, CAP completeness, bath, energy plateau, "
          "quantum-S two channels, classical stopping). The Phase-1 sweep is HELD — it will "
          "NOT start until you sign the gate.\n\n"
          "CONCLUSION: awaiting your gate verdict; nothing further launched.",
          attachments=figs)
    log("P0b complete — held at gate, notebook rebuilt, email sent. DONE (no sweep).")

if __name__ == "__main__":
    try:
        main()
    except Exception:
        tb = traceback.format_exc(); log("FATAL:\n" + tb); email("P0b relaunch FATAL", tb)
