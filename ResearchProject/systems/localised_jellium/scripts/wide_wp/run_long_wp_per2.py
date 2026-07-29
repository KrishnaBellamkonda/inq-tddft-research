#!/usr/bin/env python3
"""
Hardened overnight orchestrator for the wide-WP open-z LONG run (2026-07-03).

Chain (single GPU, GPU 0):
  1. GS   — build-once binary gs/run; skipped if the per2 checkpoint already exists.
  2. WP   — long open-z run (periodicity 2, box 50x50x111, CAP eta -1.0 / 14 Bohr-side,
            tau = 3 x rigid end-to-end traversal = ~71 a.u. => N_STEPS=1775, dt=0.04).
            Aim: watch E_total(t) reach a long-time PLATEAU (deposited energy).
  3. NOTEBOOK — run-notebook skill builder (full single-run battery incl. the locked
            Fourier pipeline via fft_pipeline_panel) -> hypotheses/wide_wp/.
  4. EMAIL — four-part result email with the E_total(t) plateau + N(t) figure.

Crash-hardening (the 2026-07-01 SIGHUP lesson):
  * every child detached (start_new_session=True) -> immune to controlling-terminal SIGHUP;
  * liveness guard: if a live run's log is silent > STALL_MIN with no completion sentinel,
    kill + relaunch ONCE; a second death -> email an alert and abort (no 8 h dead poll);
  * MAX_HOURS wall cap; failure at any stage emails an alert.

START (detached, survives logout):
  cd .../scripts/wide_wp
  setsid nohup /local/data/public/skcb2/tddft/venv/bin/python3 run_long_wp_per2.py \
      > run_long_wp_per2.log 2>&1 &
"""
from __future__ import annotations
import os, sys, time, signal, subprocess, traceback
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
LJ   = ROOT / "ResearchProject/systems/localised_jellium"
WWP  = LJ / "scripts/wide_wp"
GSDIR_SRC, WPDIR = WWP / "gs", WWP / "wp"
HYP  = LJ / "hypotheses/wide_wp"
GS_CKPT = LJ / "shared_gs/slab_n82_L50x50x111_h0p40_per2"
PY   = str(ROOT / "venv/bin/python3")
BUILDER = ROOT / ".claude/skills/run-notebook/run_notebook_builder.py"
STACK = str(ROOT / "inq-stack/python")
TO   = "chiddukanna@gmail.com"

GPU        = 0
POLL_S     = 120
STALL_MIN  = 25          # log silent this long (no sentinel) => presumed dead
MAX_HOURS  = 14

# ---- run config (grilling decisions 2026-07-03) --------------------------------
K0        = 4.6957                      # E = 300 eV, v = k0 (m_e=1)
LZ        = 111.0
DT        = 0.04
TAU_AU    = 3.0 * LZ / K0               # 3 x rigid end-to-end traversal ~ 70.9 a.u.
N_STEPS   = round(TAU_AU / DT)          # ~1775
# NOTE: wp/run.cpp prepends "results/" to LJ_OUT, so LJ_OUT must be the BARE name
# (not "results/..."); the actual output tree is WPDIR/results/<WP_OUT>.
WP_OUT    = "wp_per2_E300_long"
GS_ENV    = dict(LJ_GS_DIR=str(GS_CKPT), LJ_SPACING=0.40)
WP_ENV    = dict(LJ_OUT=WP_OUT, LJ_K0=K0, LJ_N_STEPS=N_STEPS, LJ_DT=DT, LJ_CAP=1,
                 LJ_WRITE_EVERY=6, LJ_WF_EVERY=40, LJ_LAUNCH_Z=-26.5,
                 LJ_GS_DIR=str(GS_CKPT))

ENV_BASE = {**os.environ,
            "INQ_SHARE_PATH":       str(ROOT / "inq/install/share"),
            "PSEUDOPOD_SHARE_PATH": str(ROOT / "inq/install/share/pseudopod"),
            "INQ_SOURCE":           str(ROOT / "inq-study")}


def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)


def email(subject, body, attachments=None):
    try:
        sys.path.insert(0, STACK)
        from inqview.email import send_run_email
        send_run_email(subject=subject, body=body,
                       attachments=[str(a) for a in (attachments or [])], to=TO)
        log(f"  emailed: {subject}")
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


def launch(cwd: Path, out_log: str, overrides: dict) -> subprocess.Popen:
    env = {**ENV_BASE, "CUDA_VISIBLE_DEVICES": str(GPU),
           **{k: str(v) for k, v in overrides.items()}}
    lf = open(cwd / out_log, "w")
    return subprocess.Popen([str(cwd / "run")], cwd=str(cwd), env=env,
                            stdout=lf, stderr=subprocess.STDOUT, start_new_session=True)


def run_to_completion(name, cwd, out_log, overrides, sentinel_dir, t_deadline) -> bool:
    """Launch + guard a single binary; relaunch once on a silent-death. Returns True on sentinel."""
    for attempt in (1, 2):
        if sentinel(sentinel_dir):
            log(f"{name}: already complete."); return True
        log(f"{name}: launch attempt {attempt} (N_STEPS={overrides.get('LJ_N_STEPS','-')})")
        proc = launch(cwd, out_log, overrides)
        logpath = cwd / out_log
        while True:
            time.sleep(POLL_S)
            if sentinel(sentinel_dir):
                log(f"{name}: DONE (sentinel)."); return True
            if time.time() > t_deadline:
                log(f"{name}: MAX_HOURS exceeded — abort."); _kill(proc); return False
            alive = proc.poll() is None
            silent_min = (time.time() - logpath.stat().st_mtime) / 60 if logpath.exists() else 999
            if not alive and not sentinel(sentinel_dir):
                log(f"{name}: process EXITED without sentinel (rc={proc.returncode})."); break
            if silent_min > STALL_MIN:
                log(f"{name}: log silent {silent_min:.0f} min > {STALL_MIN} — presumed dead, kill+retry.")
                _kill(proc); break
        # fell through -> dead; loop retries once
    log(f"{name}: FAILED after 2 attempts."); return False


def _kill(proc):
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM); time.sleep(5)
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except Exception:
        pass


def plateau_fig(results_dir: Path, out_png: Path):
    """E_total(t) plateau (top) + N_total(t) CAP-drain (bottom). Guaranteed email plot."""
    import numpy as np, matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, STACK)
    try:
        from inqview.visualisation.style import apply_theme; apply_theme()
    except Exception:
        pass
    obs = np.genfromtxt(results_dir / "raw/observables/observables.csv", delimiter=",", names=True)
    HA = 27.21138625
    t, E = obs["time_au"], obs["energy_total"] * HA
    fig, ax = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    ax[0].plot(t, E - E[0], lw=1.3)
    ax[0].set_ylabel(r"$E_{\rm total}(t)-E_{\rm total}(0)$  (eV)")
    ax[0].set_title("Open-z long run: total-energy evolution toward plateau")
    ax[0].axvline(LZ/2 / K0, ls=":", c="0.5")          # rigid slab-centre time (context)
    try:
        ncsv = np.genfromtxt(results_dir / "raw/observables/electron_number.csv",
                             delimiter=",", names=True)
        ax[1].plot(ncsv["time_au"], ncsv["N_total"], lw=1.3, c="C3")
    except Exception:
        pass
    ax[1].set_ylabel(r"$N_{\rm total}(t)$"); ax[1].set_xlabel("time (a.u.)")
    ax[1].set_title("Total electron number (CAP absorption of the WP)")
    fig.tight_layout(); fig.savefig(out_png, dpi=150); plt.close(fig)
    return float(E[-1] - E[0]), float(t[-1])


def main():
    t_deadline = time.time() + MAX_HOURS * 3600
    log(f"START long open-z WP run. tau={TAU_AU:.1f} a.u., N_STEPS={N_STEPS}, GPU {GPU}")
    try:
        # 1. GS (skip if checkpoint present)
        if GS_CKPT.is_dir() and any(GS_CKPT.iterdir()):
            log("GS: per2 checkpoint present — skip.")
        else:
            # Drop any STALE gs/results/run_summary.txt (e.g. the old LZ=101 PBC run)
            # so the generic sentinel does not falsely report GS complete.
            (GSDIR_SRC / "results" / "run_summary.txt").unlink(missing_ok=True)
            if not run_to_completion("GS", GSDIR_SRC, "gs_per2.log", GS_ENV,
                                     GSDIR_SRC / "results", t_deadline):
                email("[wide-wavepacket] OPEN-Z RUN FAILED — GS did not complete",
                      "The open-z GS build/run failed or stalled. See gs_per2.log. Nothing downstream ran.")
                return
            log("GS: complete.")

        # 2. WP long run  (run.cpp writes under WPDIR/results/<WP_OUT>)
        wp_res = WPDIR / "results" / WP_OUT
        if not run_to_completion("WP", WPDIR, "wp_per2_long.log", WP_ENV, wp_res, t_deadline):
            email("[wide-wavepacket] OPEN-Z LONG WP RUN FAILED / STALLED",
                  f"The long open-z WP run did not complete (2 attempts). See {WPDIR/'wp_per2_long.log'}.\n"
                  f"Partial data (if any) in {wp_res}.")
            return
        log("WP: complete.")
        finish(wp_res)
    except Exception:
        log("FATAL:\n" + traceback.format_exc())
        email("[wide-wavepacket] OPEN-Z ORCHESTRATOR CRASHED",
              "The orchestrator hit an unhandled exception:\n\n" + traceback.format_exc())


def finish(wp_res: Path):
    """Post-WP: plateau figure (guaranteed) + run-notebook (battery + Fourier) + 4-part email.
    Reused by the standalone waiter (finish_wp.py) so a long run is never orphaned."""
    try:
        # 3. plateau figure (guaranteed) + run-notebook (full battery + Fourier)
        HYP.mkdir(parents=True, exist_ok=True)
        fig_png = HYP / "wp_per2_E300_long_plateau.png"
        dE, tf = plateau_fig(wp_res, fig_png)
        nb = HYP / "wp_per2_E300_long_run_notebook.ipynb"
        log("NOTEBOOK: building run-notebook ...")
        nb_ok = True
        try:
            subprocess.run([PY, str(BUILDER), str(wp_res), str(nb),
                            "--run-cpp", str(WPDIR / "run.cpp")],
                           env={**ENV_BASE, "PYTHONPATH": STACK}, check=True,
                           timeout=3600)
        except Exception as e:
            nb_ok = False; log(f"NOTEBOOK build failed: {e}")

        # 4. result email
        body = f"""HYPOTHESIS
  A wide near-rigid WP (sigma=3.5) fired through the localised jellium slab under
  open-z (periodicity 2) boundaries deposits energy that, after the CAP fully
  absorbs the packet, leaves E_total(t) at a constant long-time PLATEAU = the
  deposited energy (the basis of the energy-method stopping power).

WHAT WAS DONE
  - Localised jellium slab, box 50x50x111 Bohr, dx=0.40, periodicity 2 (open-z);
    CAP two-sided sin^2, eta=-1.0 Ha, 14 Bohr/side, region +/-41.5..+/-55.5.
  - Single WP run: sigma_WP=3.5, E=300 eV (k0={K0}, v={K0} a.u.), launch z0=-26.5,
    dt={DT}, N_STEPS={N_STEPS} (tau={TAU_AU:.1f} a.u. = 3x rigid end-to-end traversal).

PLOT (attached: {fig_png.name})
  Top: E_total(t)-E_total(0) in eV vs time — look for a flat long-time plateau.
  Bottom: N_total(t) — the CAP draining the WP (Q: 83 -> 82); the plateau forms
  once this settles.

CONCLUSION
  Run COMPLETE. dE_total(end) = {dE:.2f} eV over tau={tf:.1f} a.u.
  CAVEAT (open-z, to debug later): under periodicity 2 the plateau carries the
  G=0 monopole step (~L_z^2*dQ^2) once the CAP drains the WP, so this level is
  NOT yet a physical deposited energy — Q(t) is logged for the monopole
  correction. Judge the plateau SHAPE (did E_total go flat?), not its level.
  Run-notebook (full battery + locked Fourier panel): {'built at ' + str(nb) if nb_ok else 'BUILD FAILED — see log'}.
"""
        email("[wide-wavepacket] open-z long WP run COMPLETE — E_total plateau",
              body, attachments=[fig_png])
        log("ALL DONE.")
    except Exception:
        log("FATAL:\n" + traceback.format_exc())
        email("[wide-wavepacket] OPEN-Z ORCHESTRATOR CRASHED",
              "The orchestrator hit an unhandled exception:\n\n" + traceback.format_exc())


if __name__ == "__main__":
    main()
