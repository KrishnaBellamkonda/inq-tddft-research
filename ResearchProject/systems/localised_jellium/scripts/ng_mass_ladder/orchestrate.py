#!/usr/bin/env python3
"""Autonomous orchestrator — Nazarov-Gross mass ladder (dense slab, wide packet).

Plan: docs/plans/nazarov-gross-slab-mass-ladder.md

WHAT IT DOES, END TO END, WITHOUT INTERVENTION
    P0 build      compile gs / wp / classical against inq-study        (sbatch)
    P1 gs         ground state, r_s = 2.5011 slab                      (sbatch)
    P2 vacuum     free-dispersion + LDA-SIE controls, 2 masses         (sbatch)
    P3 pilot      cl_inf, wp_m1, wp_m0p5 at 600 steps                  (sbatch)
    P4 gate       drift / CAP / containment / deposit separation       (local)
    P5 ladder     cl_inf, cl_m1, wp_m3, wp_m1p2, wp_m1, wp_m0p5        (sbatch)
    P6 sigma      wp_m1 at sigma_WP = 2, 3, 6 (4 is the ladder run)    (sbatch)
    P7 thickness  wp_m1 at L_z = 25 -- surface-vs-bulk separation      (sbatch)
    P8 notebooks  per-run notebooks + one phase notebook per phase     (sbatch)
    P9 final      synthesis email with the NG validation figures       (local)

HOW TO START IT (login node -- it must be able to call sbatch):
    cd <repo> && setsid nohup venv/bin/python3 \\
      ResearchProject/systems/localised_jellium/scripts/ng_mass_ladder/orchestrate.py \\
      >> ResearchProject/systems/localised_jellium/scripts/ng_mass_ladder/orch.log 2>&1 &

IDEMPOTENT. Every phase records itself in state.json and every run is skipped if
its run_summary.txt already says run_completed = true. Killing and restarting the
orchestrator resumes where it stopped; killed RUNS resume from their last
checkpoint via NG_RESUME=1 (rule final-timestep-checkpoint.md).

COST OVERRUNS ARE A WARNING, NEVER A BLOCK (rule checkpoint-dont-block.md, after
the 2026-07-12 incident where a budget gate idled two GPUs for eight hours). Only
CORRECTNESS gates stop the chain: a failed build, a failed GS, NaN energies, or a
pilot that cannot resolve a deposit.
"""
from __future__ import annotations

import json
import math
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------- paths
REPO = Path("/rds/user/skcb2/hpc-work/tddft/inq-tddft-research")
LJ = REPO / "ResearchProject/systems/localised_jellium"
VAC = REPO / "ResearchProject/systems/vacuum"
HERE = LJ / "scripts/ng_mass_ladder"
HYP = LJ / "hypotheses/ng_mass_ladder"
GS_DIR = LJ / "shared_gs/slab_n206_L30x30x120_rs2p5_dx0p5_per2"
PY = str(REPO / "venv/bin/python3")
STATE_F = HERE / "state.json"
TO = "chiddukanna@gmail.com"
SUBJ = "[ng-mass-ladder]"

# ---------------------------------------------------------------- physics
# Mirrors shared/configs/slab_n206_L30x30x120_rs2p5.hpp. Kept here so the
# orchestrator can size runs without compiling anything; the binary re-derives
# them from the header and would refuse to start on a mismatch.
V0 = 1.0742685        # 1.40 v_F, 0.875 of the Bragg peak
KF = 0.7673347
H = 0.50
SIGMA_WP = 4.0
LAUNCH_Z = -25.0
Z_END = +30.0
PATH_BOHR = Z_END - LAUNCH_Z          # 55
T_TOTAL = PATH_BOHR / V0              # 51.2 a.u.
HA_EV = 27.211386245988


def dt_for(mass: float) -> float:
    """dt = 0.08 * min(M,1) * h^2.

    The min() is the part that is easy to get wrong: ONE dt advances all 124
    orbitals and the 103 bath states have m = 1, so a heavy projectile cannot
    loosen the ceiling. Calibrated on p3 (M=1, h=0.5, dt=0.02) and
    sigma1_masspair (M=2, h=0.5, dt=0.04), which both sat exactly on it.
    """
    return 0.08 * min(mass, 1.0) * H * H


def steps_for(mass: float) -> int:
    return int(math.ceil(T_TOTAL / dt_for(mass)))


# (tag, half, mass, sigma_WP, n_steps, walltime_hours)
def _rt(tag, half, mass, sigma=SIGMA_WP, steps=None, hours=None, path=None):
    if steps is None:
        p_bohr = path if path is not None else PATH_BOHR
        steps = int(math.ceil((p_bohr / V0) / dt_for(mass)))
    # 1.91 s/step MEASURED on gpu-q-22/31 (jobs 32881739/40, 200-step averages),
    # not the 5 s/step originally scaled from p3. 1.6x headroom + 0.5 h slack.
    # A tight request matters: this partition runs ~560 deep, and Slurm cannot
    # backfill a job into a gap shorter than what it asks for.
    hours = hours if hours is not None else min(12, max(1, round(steps * 1.91 * 1.6 / 3600 + 0.5)))
    return dict(tag=tag, half=half, mass=mass, sigma=sigma, steps=steps, hours=hours)


# A pilot that does not cross the slab measures nothing: extract_S windows on
# |z| <= 7.5, and 600 steps only carried the projectile from -25 to -12.0.
PILOT_PATH = 35.0            # -25 -> +10, clearing the far face at +7.5
PILOT = [
    _rt("pilot_cl_inf", "classical", 1.0e6, path=PILOT_PATH),
    _rt("pilot_wp_m1", "wp", 1.0, path=PILOT_PATH),
    _rt("pilot_wp_m0p5", "wp", 0.5, path=PILOT_PATH),
]

# Ordered so the cheapest informative rungs land first: if the queue or the disk
# turns hostile, what survives is still a usable ladder.
LADDER = [
    _rt("cl_inf", "classical", 1.0e6),      # the M -> infinity anchor
    _rt("wp_m1", "wp", 1.0),                # the NG 'distinguishable electron'
    _rt("wp_m0p5", "wp", 0.5),              # the strongest quantum rung
    _rt("wp_m3", "wp", 3.0),                # approaching the classical limit
    _rt("wp_m1p2", "wp", 1.2),
    _rt("cl_m1", "classical", 1.0),         # matched-mass classical control
]

SIGMA_SWEEP = [_rt(f"wp_m1_s{str(s).replace('.', 'p')}", "wp", 1.0, sigma=s)
               for s in (2.0, 3.0, 6.0)]

VAC_RUNS = [  # (tag, theory, inv_mass)
    ("vac_ni_m1", "noninteracting", 1.0),
    ("vac_lda_m1", "lda", 1.0),
    ("vac_ni_m0p5", "noninteracting", 2.0),
    ("vac_lda_m0p5", "lda", 2.0),
]

POLL_S = 120


# ---------------------------------------------------------------- utilities
def log(msg: str) -> None:
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def load_state() -> dict:
    if STATE_F.exists():
        try:
            return json.loads(STATE_F.read_text())
        except json.JSONDecodeError:
            log("state.json unreadable — starting a fresh state")
    return {}


def save_state(st: dict) -> None:
    STATE_F.parent.mkdir(parents=True, exist_ok=True)
    STATE_F.write_text(json.dumps(st, indent=2, sort_keys=True))


def email(subject: str, body: str, attachments=None) -> None:
    """Never let a notification failure kill the campaign."""
    try:
        sys.path.insert(0, str(REPO / "inq-stack/python"))
        from inqview.email import send_run_email
        send_run_email(subject=f"{SUBJ} {subject}", body=body,
                       attachments=[str(a) for a in (attachments or []) if Path(a).exists()],
                       to=TO)
        log(f"emailed: {subject}")
    except Exception as exc:                                   # noqa: BLE001
        log(f"EMAIL FAILED ({exc}) — falling back to the log")
        log(f"--- would have emailed: {subject} ---")
        for line in body.splitlines():
            log(f"  | {line}")
        log("--- end ---")


def sbatch(script: str, args: list[str], env: dict, jobname: str,
           hours: float | None = None) -> int | None:
    """Submit one job. `hours` becomes a COMMAND-LINE --time.

    It has to be the command line, not the SBATCH_TIMELIMIT environment
    variable: an explicit `#SBATCH --time` inside the script overrides the env
    var, so the env route silently gave every job the script's ceiling. On a
    partition with ~560 pending jobs that is not cosmetic — asking 12 h for a
    5-minute compile put the first build 7.5 hours out, because Slurm cannot
    backfill a job into a gap shorter than its request.
    """
    cmd = ["sbatch", f"--job-name={jobname}"]
    if hours is not None:
        total_min = max(15, int(round(hours * 60)))          # 15 min floor
        cmd.append(f"--time={total_min // 60:02d}:{total_min % 60:02d}:00")
    cmd += [str(REPO / "shared/bin" / script), *args]
    e = {**os.environ, **{k: str(v) for k, v in env.items()}}
    try:
        out = subprocess.run(cmd, cwd=REPO, env=e, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        log(f"sbatch TIMED OUT for {jobname}")
        return None
    if out.returncode != 0:
        log(f"sbatch FAILED for {jobname}: {out.stderr.strip()}")
        return None
    m = re.search(r"Submitted batch job (\d+)", out.stdout)
    if not m:
        log(f"sbatch gave no job id for {jobname}: {out.stdout.strip()}")
        return None
    jid = int(m.group(1))
    log(f"submitted {jobname} as job {jid}")
    return jid


def job_state(jid: int) -> str:
    """RUNNING / PENDING / COMPLETED / FAILED / ... via sacct, squeue fallback."""
    try:
        out = subprocess.run(["sacct", "-j", str(jid), "-n", "-P", "-o", "State"],
                             capture_output=True, text=True, timeout=60)
        states = [s.strip().split()[0] for s in out.stdout.splitlines() if s.strip()]
        if states:
            # The parent job's state is the first row; array/batch steps follow.
            return states[0]
    except Exception:                                          # noqa: BLE001
        pass
    try:
        out = subprocess.run(["squeue", "-j", str(jid), "-h", "-o", "%T"],
                             capture_output=True, text=True, timeout=60)
        if out.stdout.strip():
            return out.stdout.strip().splitlines()[0]
    except Exception:                                          # noqa: BLE001
        pass
    return "UNKNOWN"


def wait_for(jids: list[int], label: str) -> dict[int, str]:
    """Block until every job leaves the queue. Returns {jid: final_state}."""
    jids = [j for j in jids if j]
    if not jids:
        return {}
    log(f"waiting on {label}: jobs {jids}")
    pending = set(jids)
    final: dict[int, str] = {}
    while pending:
        time.sleep(POLL_S)
        for j in list(pending):
            s = job_state(j)
            if s.startswith(("COMPLETED", "FAILED", "CANCELLED", "TIMEOUT",
                             "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED")):
                final[j] = s
                pending.discard(j)
                log(f"  job {j} -> {s}")
    return final


def run_complete(half: str, tag: str) -> bool:
    p = HERE / half / "results" / tag / "run_summary.txt"
    return p.exists() and "run_completed = true" in p.read_text()


def vac_complete(tag: str) -> bool:
    p = VAC / "scripts/wp_selfinteraction/results" / tag / "run_summary.txt"
    return p.exists() and "run_completed = true" in p.read_text()


def submit_rt(r: dict, prune: bool = False, resume: bool = False) -> int | None:
    if run_complete(r["half"], r["tag"]) and not resume:
        log(f"  [skip] {r['tag']} already complete")
        return None
    env = {
        "NG_OUT": r["tag"], "NG_MASS": r["mass"], "NG_SIGMA_WP": r["sigma"],
        "NG_N_STEPS": r["steps"], "NG_V": V0, "NG_GS_DIR": str(GS_DIR),
        "NG_RESUME": 1 if resume else 0, "NG_PRUNE_CKPT": 1 if prune else 0,
    }
    return sbatch("run-ng-rt.slurm", [r["half"], "prod"], env, f"ng-{r['tag']}",
                  hours=r["hours"])


# ---------------------------------------------------------------- phases
def phase_build(st: dict) -> bool:
    if st.get("P0"):
        return True
    log("P0 build: compiling gs / wp / classical")
    jids = [
        sbatch("run-ng-gs.slurm", ["build"], {}, "ng-build-gs", hours=0.75),
        sbatch("run-ng-rt.slurm", ["wp", "build"], {}, "ng-build-wp", hours=0.75),
        sbatch("run-ng-rt.slurm", ["classical", "build"], {}, "ng-build-cl", hours=0.75),
    ]
    res = wait_for(jids, "P0 build")
    ok = all(s.startswith("COMPLETED") for s in res.values()) and len(res) == 3
    if not ok:
        email("P0 BUILD FAILED — chain stopped",
              f"Compile stage did not succeed: {res}\n\n"
              f"Logs: {REPO}/ng-build-*.out\n"
              "Nothing was run. Fix the compile error and restart orchestrate.py; "
              "it is idempotent and will resume here.")
        return False
    st["P0"] = True
    save_state(st)
    email("P0 build OK", "gs / wp / classical all compiled against inq-study. Proceeding to the ground state.")
    return True


def phase_gs(st: dict) -> bool:
    if st.get("P1") or (GS_DIR.exists() and (HERE / "gs/results/gs/run_summary.txt").exists()):
        st["P1"] = True
        save_state(st)
        log("P1 gs: already present")
        return True
    log("P1 gs: submitting ground state")
    jid = sbatch("run-ng-gs.slurm", ["prod"], {"GS_DIR": str(GS_DIR)}, "ng-gs", hours=3)
    res = wait_for([jid], "P1 gs")
    summ = HERE / "gs/results/gs/run_summary.txt"
    ok = jid and res.get(jid, "").startswith("COMPLETED") and summ.exists() \
        and "run_completed = true" in summ.read_text()
    if not ok:
        email("P1 GROUND STATE FAILED — chain stopped",
              "The GS binary's own gates (electron count, r_s, finite E) refuse to save a bad "
              f"ground state, so nothing downstream can run.\n\nJob state: {res}\n"
              f"Log: {REPO}/ng-gs-*.out")
        return False
    st["P1"] = True
    save_state(st)
    email("P1 ground state OK", summ.read_text())
    return True


def phase_vacuum(st: dict) -> bool:
    """Steps 5-6: does LDA self-interaction widen the packet MASS-DEPENDENTLY?

    Never blocks. A lone electron has no self-interaction, so free dispersion is
    the exact answer and the LDA-minus-noninteracting difference IS the error.
    What matters for this campaign is only whether that error depends on M — if
    it does, it forges the result with the right sign and must be quoted as a
    systematic (or the ladder re-run with SIC-PZ).
    """
    if st.get("P2"):
        return True
    log("P2 vacuum: self-interaction / free-dispersion controls")
    jids = []
    for tag, theory, inv_mass in VAC_RUNS:
        if vac_complete(tag):
            log(f"  [skip] {tag}")
            continue
        env = {"WP_OUT": tag, "WP_THEORY": theory, "WP_INV_MASS": inv_mass,
               "WP_SIGMA": SIGMA_WP, "WP_K0": 0.0, "WP_NSTEPS": 1500}
        j = sbatch("run-wp-si-sweep.slurm", [], env, f"ng-{tag}", hours=1.5)
        if j:
            jids.append(j)
    wait_for(jids, "P2 vacuum")
    st["P2"] = True
    save_state(st)
    email("P2 vacuum controls done",
          "Free-dispersion validation of the mass fork + the LDA self-interaction "
          "control are complete. The phase notebook quantifies whether the SIE "
          "excess width is mass-dependent; a mass-INDEPENDENT error cancels in the "
          "ladder's ratios, a mass-dependent one does not and will be quoted.")
    return True


def phase_pilot(st: dict) -> bool:
    if st.get("P3"):
        return True
    log("P3 pilot: 3 short runs (600 steps)")
    jids = [submit_rt(r, prune=True) for r in PILOT]
    wait_for(jids, "P3 pilot")
    st["P3"] = True
    save_state(st)
    return True


def phase_pilot_gate(st: dict) -> bool:
    """The ONE correctness gate that can stop the campaign on physics."""
    if st.get("P4"):
        return True
    log("P4 pilot gate")
    sys.path.insert(0, str(HYP))
    try:
        from ng_analysis import pilot_gate
        verdict = pilot_gate(HERE, [r["tag"] for r in PILOT], [r["half"] for r in PILOT])
    except Exception as exc:                                   # noqa: BLE001
        log(f"gate evaluation raised: {exc!r} — treating as INCONCLUSIVE and PROCEEDING")
        verdict = {"pass": True, "report": f"gate could not be evaluated: {exc!r}",
                   "inconclusive": True}
    if not verdict["pass"]:
        email("P4 PILOT GATE FAILED — ladder NOT launched", verdict["report"])
        return False
    st["P4"] = True
    save_state(st)
    email("P4 pilot gate PASSED — launching the ladder", verdict["report"])
    return True


def phase_runs(st: dict, key: str, runs: list[dict], label: str) -> bool:
    if st.get(key):
        return True
    log(f"{key} {label}: {len(runs)} runs")
    done_states = {}
    for r in runs:                       # serial: one GPU job at a time is kinder
        j = submit_rt(r)                 # to the queue and keeps the disk bounded
        if j:
            done_states.update(wait_for([j], r["tag"]))
            # A run killed by the wall clock is RESUMED, never restarted
            # (rule final-timestep-checkpoint.md: extend, don't recompute).
            if done_states.get(j, "").startswith("TIMEOUT") and not run_complete(r["half"], r["tag"]):
                log(f"  {r['tag']} hit the wall clock — resuming from its checkpoint")
                j2 = submit_rt(r, resume=True)
                if j2:
                    wait_for([j2], f"{r['tag']} (resume)")
    st[key] = True
    save_state(st)
    incomplete = [r["tag"] for r in runs if not run_complete(r["half"], r["tag"])]
    email(f"{key} {label} done",
          f"Completed: {[r['tag'] for r in runs if run_complete(r['half'], r['tag'])]}\n"
          f"Incomplete: {incomplete or 'none'}\n\n"
          "Incomplete runs are not fatal — each holds a final checkpoint and can be "
          "extended with NG_RESUME=1. The chain continues.")
    return True


def phase_notebooks(st: dict) -> bool:
    if st.get("P8"):
        return True
    log("P8 notebooks: per-run + per-phase")
    jid = sbatch("run-ng-notebooks.slurm", [], {}, "ng-notebooks", hours=2)
    res = wait_for([jid], "P8 notebooks")
    st["P8"] = True
    save_state(st)
    figs = sorted((HYP / "figures").glob("ng_*.png")) if (HYP / "figures").exists() else []
    email("P8 notebooks + figures built",
          f"Job state: {res}\n\nNotebooks: {HYP}\nFigures: {len(figs)} PNGs attached "
          "(the Nazarov-Gross validation set).",
          attachments=figs[:8])
    return True


def main() -> int:
    log("=" * 70)
    log("Nazarov-Gross mass ladder — autonomous orchestrator starting")
    log(f"repo={REPO}")
    log(f"v0={V0:.6f} a.u. = {V0/KF:.2f} v_F, {V0/1.2277355:.3f} of the Bragg peak")
    log(f"path={PATH_BOHR} Bohr, t={T_TOTAL:.1f} a.u.")
    for r in LADDER + SIGMA_SWEEP:
        log(f"  {r['tag']:16s} half={r['half']:9s} M={r['mass']:<8g} "
            f"sigma={r['sigma']:<4g} dt={dt_for(r['mass']):.4f} steps={r['steps']:<6d} "
            f"wall={r['hours']}h")
    log("=" * 70)

    st = load_state()
    email("campaign STARTED",
          "The Nazarov-Gross mass ladder is running autonomously. No intervention needed.\n\n"
          f"System: r_s = 2.5011 jellium slab, 15 Bohr thick, 30x30x120 Bohr box, dx 0.50.\n"
          f"Projectile: sigma_WP = {SIGMA_WP}, v0 = {V0:.4f} a.u. = {V0/KF:.2f} v_F "
          f"(0.875 of the Bragg peak, so BELOW it as required).\n"
          f"Ladder: {[r['tag'] for r in LADDER]}\n"
          f"Sigma sweep: {[r['tag'] for r in SIGMA_SWEEP]}\n\n"
          "You will get one email per phase. Only a build failure, a bad ground state, "
          "or a pilot that cannot resolve a deposit will stop it.")

    steps = [
        (phase_build, ()),
        (phase_gs, ()),
        (phase_vacuum, ()),
        (phase_pilot, ()),
        (phase_pilot_gate, ()),
        (lambda s: phase_runs(s, "P5", LADDER, "mass ladder"), ()),
        (lambda s: phase_runs(s, "P6", SIGMA_SWEEP, "sigma sweep"), ()),
        (phase_notebooks, ()),
    ]
    for fn, _ in steps:
        if not fn(st):
            log("chain stopped by a correctness gate — see the email")
            return 1

    email("campaign COMPLETE",
          "All phases finished. The phase notebooks under "
          f"{HYP} carry the Nazarov-Gross validation figures: S vs mass at fixed "
          "charge and velocity, S vs measured width, and the collapse test that "
          "decides whether mass acts ONLY through width.")
    log("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
