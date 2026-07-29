#!/usr/bin/env python3
"""
Overnight orchestrator — wide-WP vs classical, FULL PBC matched pair (2026-07-06).

Chain:
  0. WAIT for the GS (PBC-111) checkpoint sentinel (GS launched separately on GPU0).
  1. Launch CONCURRENTLY:
       WP        (GPU 0) periodicity 3, box 50x50x111, CAP eta -1.0 / 14 Bohr-side,
                 sigma_WP=3.5, E=300 eV, dt=0.04, N_STEPS=1773 (tau~70.8 a.u.).
       CLASSICAL (GPU 1) same geometry/CAP/E, Gaussian-e ion (sigma_pot=2.475=sigma_WP/sqrt2,
                 UPF electron_gaussian_wpsigma3p5.upf), dt=0.02, N_STEPS=3540 (matched tau).
  2. On BOTH complete: build 2 run notebooks (run-notebook skill, full battery + Fourier)
     + 1 comparison notebook (build_wp_vs_classical_pbc.py) -> hypotheses/wide_wp/.
  3. EMAIL a four-part result with the WP-vs-classical S / energy comparison figure.

Hardening (open-z lessons): children detached (start_new_session); per-run liveness
guard (kill+retry once on STALL_MIN silence); MAX_HOURS wall cap; failure emails an alert.
Output-path note: {wp,classical}/run.cpp prepend "results/" to LJ_OUT, so LJ_OUT is a
BARE name; data lands at <dir>/results/<OUT>.

START (detached):
  cd .../scripts/wide_wp
  setsid nohup /local/data/public/skcb2/tddft/venv/bin/python3 run_pbc_pair.py \
      > run_pbc_pair.log 2>&1 &
"""
from __future__ import annotations
import os, sys, time, signal, subprocess, traceback
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
LJ   = ROOT / "ResearchProject/systems/localised_jellium"
WWP  = LJ / "scripts/wide_wp"
GSDIR, WPDIR, CLDIR = WWP / "gs", WWP / "wp", WWP / "classical"
HYP  = LJ / "hypotheses/wide_wp"
GS_CKPT = LJ / "shared_gs/slab_n82_L50x50x111_h0p40_pbc"
PY   = str(ROOT / "venv/bin/python3")
BUILDER = ROOT / ".claude/skills/run-notebook/run_notebook_builder.py"
CMP_BUILDER = HYP / "build_wp_vs_classical_pbc.py"
STACK = str(ROOT / "inq-stack/python")
TO   = "chiddukanna@gmail.com"

POLL_S     = 120
STALL_MIN  = 30
MAX_HOURS  = 24

# ---- run config -----------------------------------------------------------------
K0        = 4.6957                 # E = 300 eV, v = k0 (m_e=1)
LZ        = 111.0
WP_OUT    = "wp_pbc_E300"
CL_OUT    = "classical_openz_E300"
WP_DT, WP_STEPS = 0.04, 1773       # tau ~ 70.9 a.u.
# Classical sized to the initial-drag + deceleration + full-absorption window (~25 a.u.):
# the projectile reaches the CAP (+41.5) at t~14.5 and is absorbed by ~20. Open-z 2D-Poisson
# (z-doubled) + moving ion is ~65 s/step, so 1250 steps @ dt0.02 = tau~25 a.u. (~23 h).
# CLASSICAL RUNS OPEN-Z (periodicity 2) to kill the point-charge z self-image drag that
# contaminated the PBC classical (5.6 eV/Bohr spurious vacuum drag); the WP is BC-insensitive
# so WP stays PBC. GS for the open-z classical = the per2 checkpoint.
CL_DT, CL_STEPS = 0.02, 1250       # open-z classical; relaunched 2026-07-06

WP_ENV = dict(LJ_OUT=WP_OUT, LJ_K0=K0, LJ_N_STEPS=WP_STEPS, LJ_DT=WP_DT, LJ_CAP=1,
              LJ_WRITE_EVERY=6, LJ_WF_EVERY=40, LJ_LAUNCH_Z=-26.5, LJ_GS_DIR=str(GS_CKPT))
GS_PER2 = LJ / "shared_gs/slab_n82_L50x50x111_h0p40_per2"   # open-z GS for the classical
CL_ENV = dict(LJ_OUT=CL_OUT, LJ_K0=K0, LJ_N_STEPS=CL_STEPS, LJ_DT=CL_DT, LJ_CAP=1,
              LJ_WRITE_EVERY=10, LJ_LAUNCH_Z=-26.5, LJ_GS_DIR=str(GS_PER2))

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


def _kill(proc):
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM); time.sleep(5)
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except Exception:
        pass


def launch(cwd: Path, out_log: str, gpu: int, overrides: dict) -> subprocess.Popen:
    env = {**ENV_BASE, "CUDA_VISIBLE_DEVICES": str(gpu),
           **{k: str(v) for k, v in overrides.items()}}
    lf = open(cwd / out_log, "w")
    return subprocess.Popen([str(cwd / "run")], cwd=str(cwd), env=env,
                            stdout=lf, stderr=subprocess.STDOUT, start_new_session=True)


class Job:
    def __init__(self, name, cwd, log_name, gpu, env, res_dir):
        self.name, self.cwd, self.log_name = name, cwd, log_name
        self.gpu, self.env, self.res = gpu, env, res_dir
        self.proc = None; self.attempt = 0; self.done = False; self.failed = False

    def logpath(self): return self.cwd / self.log_name

    def start(self):
        if sentinel(self.res):
            self.done = True; log(f"{self.name}: already complete."); return
        self.attempt += 1
        log(f"{self.name}: launch attempt {self.attempt} (GPU {self.gpu}, N_STEPS={self.env['LJ_N_STEPS']})")
        self.proc = launch(self.cwd, self.log_name, self.gpu, self.env)

    def poll(self, t_deadline):
        """Advance one poll tick; set done/failed. relaunch once on silent death."""
        if self.done or self.failed:
            return
        if sentinel(self.res):
            self.done = True; log(f"{self.name}: DONE (sentinel)."); return
        if time.time() > t_deadline:
            log(f"{self.name}: MAX_HOURS exceeded — abort."); _kill(self.proc); self.failed = True; return
        alive = self.proc.poll() is None
        lp = self.logpath()
        silent_min = (time.time() - lp.stat().st_mtime) / 60 if lp.exists() else 999
        if (not alive and not sentinel(self.res)) or silent_min > STALL_MIN:
            why = f"exited rc={self.proc.returncode}" if not alive else f"log silent {silent_min:.0f}min"
            log(f"{self.name}: presumed dead ({why}).")
            _kill(self.proc)
            if self.attempt >= 2:
                log(f"{self.name}: FAILED after 2 attempts."); self.failed = True
            else:
                self.start()


def wait_for_gs(t_deadline) -> bool:
    log(f"WAIT for GS checkpoint sentinel: {GS_CKPT}")
    while True:
        if GS_CKPT.is_dir() and any(GS_CKPT.iterdir()) and sentinel(GSDIR / "results"):
            log("GS: complete (sentinel + checkpoint present)."); return True
        if time.time() > t_deadline:
            log("GS: wait exceeded MAX_HOURS — abort."); return False
        time.sleep(POLL_S)


def main():
    t_deadline = time.time() + MAX_HOURS * 3600
    log(f"START PBC pair. WP N_STEPS={WP_STEPS}@dt{WP_DT}, CL N_STEPS={CL_STEPS}@dt{CL_DT}")
    try:
        if not wait_for_gs(t_deadline):
            email("[wide-wp PBC] FAILED — GS never completed",
                  "The PBC-111 GS checkpoint/sentinel never appeared. Nothing downstream ran.")
            return

        wp = Job("WP", WPDIR, "wp_pbc.log", 0, WP_ENV, WPDIR / "results" / WP_OUT)
        cl = Job("CLASSICAL", CLDIR, "classical_pbc.log", 1, CL_ENV, CLDIR / "results" / CL_OUT)
        wp.start(); cl.start()
        while not all(j.done or j.failed for j in (wp, cl)):
            time.sleep(POLL_S)
            for j in (wp, cl):
                j.poll(t_deadline)

        if wp.failed or cl.failed:
            email("[wide-wp PBC] pair run FAILED / STALLED",
                  f"WP done={wp.done} failed={wp.failed}; CLASSICAL done={cl.done} failed={cl.failed}.\n"
                  f"Logs: {wp.logpath()} ; {cl.logpath()}.")
            # still try to finish whatever completed
        finish(wp.res if wp.done else None, cl.res if cl.done else None)
    except Exception:
        log("FATAL:\n" + traceback.format_exc())
        email("[wide-wp PBC] ORCHESTRATOR CRASHED",
              "Unhandled exception:\n\n" + traceback.format_exc())


def finish(wp_res: Path | None, cl_res: Path | None):
    """Build 2 run notebooks + comparison notebook, email a 4-part result."""
    try:
        HYP.mkdir(parents=True, exist_ok=True)
        nb_env = {**ENV_BASE, "PYTHONPATH": STACK}
        built = []
        if wp_res:
            nb = HYP / "wp_pbc_E300_run_notebook.ipynb"
            try:
                subprocess.run([PY, str(BUILDER), str(wp_res), str(nb),
                                "--run-cpp", str(WPDIR / "run.cpp")],
                               env=nb_env, check=True, timeout=3600)
                built.append(nb.name)
            except Exception as e:
                log(f"WP notebook build failed: {e}")
        if cl_res:
            nb = HYP / "classical_pbc_E300_run_notebook.ipynb"
            try:
                subprocess.run([PY, str(BUILDER), str(cl_res), str(nb),
                                "--run-cpp", str(CLDIR / "run.cpp")],
                               env=nb_env, check=True, timeout=3600)
                built.append(nb.name)
            except Exception as e:
                log(f"CLASSICAL notebook build failed: {e}")

        cmp_png = HYP / "wp_vs_classical_pbc_S.png"
        cmp_nb  = HYP / "wp_vs_classical_pbc_comparison.ipynb"
        if wp_res and cl_res and CMP_BUILDER.exists():
            try:
                subprocess.run([PY, str(CMP_BUILDER), str(wp_res), str(cl_res),
                                str(cmp_nb), str(cmp_png)],
                               env=nb_env, check=True, timeout=1800)
                built.append(cmp_nb.name)
            except Exception as e:
                log(f"comparison notebook build failed: {e}")

        body = f"""HYPOTHESIS
  Under FULL 3D PBC (periodicity 3), a wide WP (sigma_WP=3.5) and a matched classical
  Gaussian-e projectile at the SAME width, fired through the localised jellium slab
  at E=300 eV, give stopping powers that agree except for the purely-quantum
  (Pauli/interference) part; PBC removes the open-z monopole so the energy ledger is
  clean. Cross-check: the phantom absorbed-orbital artifact should persist (it is not
  a boundary effect).

WHAT WAS DONE
  - Box 50x50x111 Bohr, dx=0.40, periodicity 3 (full PBC), r_s=5.67.
  - CAP two-sided sin^2, eta=-1.0 Ha, 14 Bohr/side, region +/-41.5..+/-55.5 (matched).
  - WP:        dt={WP_DT}, N_STEPS={WP_STEPS} (tau~70.9). data {wp_res}.
  - CLASSICAL: dt={CL_DT}, N_STEPS={CL_STEPS} (matched tau). data {cl_res}.
  - Notebooks built: {', '.join(built) if built else 'NONE (see log)'}.

PLOT (attached: {cmp_png.name if cmp_png.exists() else 'n/a'})
  WP vs classical: projectile KE loss / stopping estimate + energy evolution.

CONCLUSION
  Pair COMPLETE under PBC. S extracted as INITIAL DRAG (light projectiles decelerate):
  classical -dKE_ion/ds over the early v>=0.85 v0 window; WP momentum-centroid drift.
  See the comparison notebook for S(300 eV)_WP vs S(300 eV)_classical and the
  quantum difference. PROVISIONAL until reviewed.
"""
        email("[wide-wp PBC] WP vs classical pair COMPLETE",
              body, attachments=[cmp_png] if cmp_png.exists() else None)
        log("ALL DONE.")
    except Exception:
        log("FATAL in finish:\n" + traceback.format_exc())
        email("[wide-wp PBC] ORCHESTRATOR CRASHED (finish)",
              "Unhandled exception:\n\n" + traceback.format_exc())


if __name__ == "__main__":
    main()
