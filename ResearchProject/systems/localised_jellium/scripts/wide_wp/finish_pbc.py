#!/usr/bin/env python3
"""
WATCH-ONLY finisher for the user-launched PBC long WP run `wp_pbc_E300`.
Does NOT launch or relaunch anything — it only watches the existing run
(user's pid, GPU 0), then on genuine completion builds the run-notebook +
sends a 4-part email. Liveness guard emails an alert if the WP dies first.

START (detached):
  cd .../scripts/wide_wp
  setsid nohup /local/data/public/skcb2/tddft/venv/bin/python3 finish_pbc.py \
      > finish_pbc.log 2>&1 &
"""
from __future__ import annotations
import sys, time, subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import run_long_wp_per2 as O   # reuse plateau_fig / email / sentinel / log / constants

RUN_NAME = "wp_pbc_E300"                          # user's LJ_OUT (bare -> results/wp_pbc_E300)
WP_RES   = O.WPDIR / "results" / RUN_NAME
WP_LOG   = O.WPDIR / "wp_pbc.log"
NB       = O.HYP / f"{RUN_NAME}_run_notebook.ipynb"
FIG      = O.HYP / f"{RUN_NAME}_plateau.png"
DEADLINE_H, POLL_S, STALL_MIN = 24, 120, 40


def wp_alive() -> bool:
    try:
        subprocess.check_output(["pgrep", "-f", "wide_wp/wp/run"]); return True
    except subprocess.CalledProcessError:
        return False


def main():
    t_deadline = time.time() + DEADLINE_H * 3600
    O.log(f"PBC WAITER: watching {WP_RES} (24 h). NO relaunch.")
    while True:
        if O.sentinel(WP_RES):
            O.log("PBC WAITER: sentinel found."); break
        if time.time() > t_deadline:
            O.email("[wide-wavepacket] PBC WP — WAITER TIMEOUT (24 h)",
                    f"{RUN_NAME} did not complete within 24 h. See {WP_LOG}. Nothing auto-built.")
            return
        if not wp_alive() and not O.sentinel(WP_RES):
            O.email("[wide-wavepacket] PBC WP — PROCESS GONE before completion",
                    f"The WP process is gone with no completion sentinel for {RUN_NAME}.\n"
                    f"See {WP_LOG}; partial data in {WP_RES}.")
            O.log("PBC WAITER: WP gone without sentinel — alerted."); return
        time.sleep(POLL_S)

    # completed -> plateau fig + notebook + email
    O.HYP.mkdir(parents=True, exist_ok=True)
    dE, tf = O.plateau_fig(WP_RES, FIG)
    nb_ok = True
    try:
        subprocess.run([O.PY, str(O.BUILDER), str(WP_RES), str(NB),
                        "--run-cpp", str(O.WPDIR / "run.cpp")],
                       env={**O.ENV_BASE, "PYTHONPATH": O.STACK}, check=True, timeout=5400)
    except Exception as e:
        nb_ok = False; O.log(f"notebook build failed: {e}")

    body = f"""HYPOTHESIS
  Repeat of the wide-WP (sigma=3.5, E=300 eV) slab run under FULL 3D PBC
  (periodicity 3) instead of open-z, to test whether the E_total oscillation
  seen in the open-z run was a G=0-monopole artifact or genuine dynamics, and
  whether E_total reaches a clean long-time plateau (basis of the energy S method).

WHAT WAS DONE
  - Localised jellium slab, box 50x50x111 Bohr, dx=0.40, FULL PBC (periodicity 3);
    CAP sin^2 eta=-1.0 Ha, 14 Bohr/side. GS E=-99.3 Ha (physical, no monopole).
  - ONE WP run, sigma_WP=3.5, E=300 eV, launch z0=-26.5, dt=0.04, N_STEPS=1773
    (tau=70.9 a.u.). COMPLETED.

PLOT (attached: {FIG.name})
  Top: E_total(t)-E_total(0) (eV). Bottom: N_total(t) (CAP absorption 83->82).

CONCLUSION
  dE_total(end) = {dE:.1f} eV at tau={tf:.1f} a.u.  Compare the E_total(t) SHAPE to
  the open-z run: partial data already showed the +/-25 eV oscillation PERSISTS under
  PBC, i.e. it is NOT the open-z monopole. Likely the CAP draining energy as density
  sloshes through the absorber (a CAP makes E_total non-conserved). Judge whether a
  clean plateau forms; if not, the energy S method needs time-averaging / rethink.
  Run-notebook (battery + Fourier panel): {'built at ' + str(NB) if nb_ok else 'BUILD FAILED'}.
"""
    O.email(f"[wide-wavepacket] PBC WP run DONE — {RUN_NAME} (plateau vs open-z compare)",
            body, attachments=[FIG])
    O.log("PBC WAITER: done.")


if __name__ == "__main__":
    main()
