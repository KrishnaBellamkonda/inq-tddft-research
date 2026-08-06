"""
Autonomous finalizer for the sigma56_sv twin campaign: status -> repair -> report.

    python finalize.py --attempt 1 --max 2
    python finalize.py --status-only          # inspect, change nothing

WHAT IT IS FOR. A SLURM chain gets the runs launched; it does not get them
FINISHED. A run can come back short (walltime, preemption, a node fault) or never
start at all (a dependency that was never satisfied), and the one-shot notebook
stage would then quietly build a figure out of whatever happened to be there. This
driver closes that loop: it checks every expected run against its step target,
RESUMES the ones that fell short IN PLACE, and only then builds the deliverables —
and it reports regardless of outcome.

--------------------------------------------------------------------------------
DESIGN CONSTRAINTS, AND WHY THEY ARE WHAT THEY ARE
--------------------------------------------------------------------------------
* NO sbatch FROM INSIDE THE JOB. Submitting work from a compute node is not
  something this repo has ever relied on, so it is not relied on here either.
  Repairs run in-process by invoking the SAME dispatcher script the chain uses
  (`bash shared/bin/run-s56-wp.slurm 2`), which does its own environment setup and
  is therefore the single source of truth for how a run is launched. No duplicated
  launch logic means no drift.

* BOUNDED, ALWAYS. `--max` caps the repair attempts and REPAIR_BUDGET_S caps the
  wall time spent repairing, leaving room to still build and report. An unbounded
  finalizer is not autonomy -- a 9.5-hour finalizer once polled for production
  that was never going to run (docs/handovers, 2026-06-28). This one always
  reaches the reporting stage.

* REPORT EVEN ON FAILURE. Silence is the one outcome that is never acceptable
  (.claude/rules/checkpoint-dont-block.md). CAMPAIGN_REPORT.md is written on every
  path, including the one where everything went wrong; email is attempted and its
  failure is logged, never raised.

* THE DEPOSIT DOES NOT NEED THE VACUUM CONTROLS. vac runs are diagnostics (they
  quantify the t=0 CAP loss), so a missing one degrades the analysis rather than
  invalidating it: they are repaired at LOWEST priority and never block the report.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import s56_stopping as S   # noqa: E402

REPO = S.REPO
BIN = REPO / "shared/bin"
V_INDEX = {v: i for i, v in enumerate(S.VELOCITIES)}

# Leave this much of the job's wall time for building figures + notebooks, which
# is the part the user actually wants; repairs get whatever is left.
REPORT_RESERVE_S = int(os.environ.get("S56_REPORT_RESERVE_S", 3 * 3600))
REPAIR_BUDGET_S = int(os.environ.get("S56_REPAIR_BUDGET_S", 30 * 3600))


# ---------------------------------------------------------------------------
# status
# ---------------------------------------------------------------------------
def expected() -> list[dict]:
    """Every run the campaign is supposed to produce, in REPAIR PRIORITY order.

    Production first (they are the measurement), CAP-free controls next (they
    calibrate it), vacuum baselines last (they only annotate it).
    """
    out = []
    for sigma in S.SIGMAS:
        for v in S.VELOCITIES:
            for half in ("wp", "classical"):
                out.append({"sigma": sigma, "v": v, "half": half, "cap": True,
                            "kind": "production"})
    for sigma in S.SIGMAS:
        out.append({"sigma": sigma, "v": 3.0, "half": "classical", "cap": False,
                    "kind": "control"})
    for sigma in S.SIGMAS:
        for v in S.VELOCITIES:
            out.append({"sigma": sigma, "v": v, "half": "vac", "cap": True,
                        "kind": "vacuum"})
    return out


def _steps_done(d: Path) -> int:
    """Last step actually on disk, across resume segments. -1 if unreadable."""
    obs = d / "raw" / "observables"
    if not obs.exists():
        return -1
    try:
        return int(S._concat(obs, "observables")["step"].iloc[-1])
    except Exception:                                          # noqa: BLE001
        # A vacuum run writes no observables.csv; fall back to its own stats file.
        try:
            return int(S._concat(obs, "wp_momentum_stats")["step"].iloc[-1])
        except Exception:                                      # noqa: BLE001
            return -1


def status() -> list[dict]:
    rows = []
    for e in expected():
        d = S.run_dir(e["sigma"], e["v"], e["half"], e["cap"])
        done = _steps_done(d)
        target = S.STEPS_TARGET[e["v"]]
        rows.append({**e, "run": d.name, "dir": d, "exists": d.exists(),
                     "steps_done": done, "steps_target": target,
                     "complete": done >= target,
                     "resumable": (d / "checkpoint").exists()})
    return rows


def print_status(rows: list[dict]) -> None:
    print(f"{'run':28s} {'kind':11s} {'steps':>13s}  state")
    for r in rows:
        state = ("COMPLETE" if r["complete"] else
                 "MISSING" if not r["exists"] else
                 f"SHORT ({'resumable' if r['resumable'] else 'no ckpt'})")
        print(f"  {r['run']:26s} {r['kind']:11s} "
              f"{r['steps_done']:>6d}/{r['steps_target']:<6d} {state}")
    n_ok = sum(r["complete"] for r in rows)
    print(f"\n  {n_ok}/{len(rows)} complete "
          f"({sum(r['complete'] for r in rows if r['kind']=='production')}"
          f"/{sum(1 for r in rows if r['kind']=='production')} production)")


# ---------------------------------------------------------------------------
# repair
# ---------------------------------------------------------------------------
def dispatcher_for(half: str) -> Path:
    return BIN / {"wp": "run-s56-wp.slurm",
                  "classical": "run-s56-cl.slurm",
                  "vac": "run-s56-vac.slurm"}[half]


def binary_dir(half: str) -> Path:
    return S.SCRIPTS / {"wp": "wp", "classical": "classical", "vac": "vac"}[half]


def ensure_binary(half: str, sigma: float, deadline: float) -> bool:
    """Build the run binary if it is missing, by running the dispatcher's own
    `smoke` stage (which is the stage that builds).

    WHY THIS EXISTS. The chain gates each sweep on `afterok` of its smoke, so a
    failed smoke means the binary was never built and every array task exits 2
    with "not built". Without this, the finalizer would faithfully re-invoke a
    dispatcher that can only fail, 16 times. With it, a smoke failure caused by a
    transient (a bad node, a full TMPDIR, a filesystem hiccup) self-heals; a
    smoke failure caused by a real compile error or a genuinely bad packet still
    fails fast, which is correct — that is a defect, not a hiccup.
    """
    if half == "vac":
        return True          # the vac dispatcher builds itself on first use
    run_bin = binary_dir(half) / "run"
    if run_bin.exists() and os.access(run_bin, os.X_OK):
        return True
    remaining = deadline - time.time()
    if remaining <= 0:
        return False
    print(f"  BUILD {half} binary missing — running its smoke stage")
    env = dict(os.environ)
    env["LJ_SIGMA"] = f"{sigma}"
    env["SLURM_SUBMIT_DIR"] = str(REPO)
    env.pop("SLURM_ARRAY_TASK_ID", None)
    try:
        p = subprocess.run(["bash", str(dispatcher_for(half)), "smoke"],
                           cwd=REPO, env=env, timeout=max(60, remaining))
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT building {half}")
        return False
    ok = run_bin.exists() and os.access(run_bin, os.X_OK)
    print(f"  smoke exit {p.returncode}; binary present: {ok}")
    # A non-zero smoke exit with a BUILT binary means the t=0 gates failed, not
    # the compile. Report it and stop: gates exist to prevent exactly the
    # production runs this function is about to enable.
    if p.returncode != 0 and ok:
        print(f"  {half} binary built but its t=0 gates FAILED — refusing to "
              f"repair production runs from a packet the gates rejected")
        return False
    return ok


# MEASURED on the 2026-08-02 smoke stages, A100, 88x88x264 grid. Used only to
# bound a repair, so being approximate is fine — being ABSENT is not.
S_PER_STEP = {"wp": 3.15, "classical": 3.00, "vac": 3.15}
STARTUP_S = 900.0     # toolchain + optional rebuild + GS/checkpoint load
SAFETY = 2.0          # a repair gets 2x its expected cost before we give up on it


def run_budget_s(r: dict) -> float:
    """Wall-clock a single repair is allowed before the finalizer abandons it.

    Sized from the steps it actually still has to do (a resume from 92 % is not
    given a from-scratch budget), so a wedged run is detected in a couple of hours
    rather than at the end of the allocation.
    """
    todo = r["steps_target"] - (r["steps_done"] if r["resumable"] else 0)
    todo = max(todo, 0)
    return SAFETY * (STARTUP_S + todo * S_PER_STEP.get(r["half"], 3.15))


def repair_one(r: dict, deadline: float) -> bool:
    """Run (or resume) one point by invoking its own dispatcher as plain bash.

    The dispatcher is the single source of truth for how a run is launched — it
    sets up the toolchain, resolves the GS, and picks the step count from the same
    table this module uses. Calling it (rather than re-deriving the command here)
    is what keeps repair and production from drifting apart.
    """
    script = dispatcher_for(r["half"])
    env = dict(os.environ)
    env["LJ_SIGMA"] = f"{r['sigma']}"
    env["SLURM_SUBMIT_DIR"] = str(REPO)
    env.pop("SLURM_ARRAY_TASK_ID", None)

    if r["half"] == "vac":
        # cap_check-derived, no resume support: a short vacuum control is redone
        # from scratch. It is a diagnostic, so this is cheap insurance rather than
        # a requirement, and it is attempted last.
        args = [str(script)]
    else:
        env["LJ_RESUME"] = "1" if r["resumable"] else "0"
        if not r["cap"]:
            env["LJ_CAP_ETA"] = "0"
        args = [str(script), str(V_INDEX[r["v"]])]

    remaining = deadline - time.time()
    if remaining <= 0:
        print(f"  SKIP {r['run']}: repair budget exhausted")
        return False
    if not ensure_binary(r["half"], r["sigma"], deadline):
        print(f"  SKIP {r['run']}: no usable binary for the {r['half']} half")
        return False
    remaining = deadline - time.time()
    # PER-RUN CAP, not the whole remaining budget. Giving one repair the entire
    # window means a single WEDGED run swallows it and every later run is skipped
    # — which is exactly what happened 2026-08-03: a resume hung for 5 h 17 min on
    # `electrons.load` against a full filesystem (zero output, zero bytes written)
    # and burned the finalizer's whole 36 h allocation without repairing anything.
    # Budget each run 2x its measured cost and move on when it overruns; its
    # checkpoint survives, so the next attempt picks it up.
    budget = run_budget_s(r)
    cap = min(remaining, budget)
    print(f"  REPAIR {r['run']}  (resume={env.get('LJ_RESUME','n/a')}, "
          f"cap {cap/3600:.1f} h, {remaining/3600:.1f} h of budget left)")
    try:
        p = subprocess.run(["bash", *args], cwd=REPO, env=env,
                           timeout=max(60, cap))
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT {r['run']} after {cap/3600:.1f} h — its checkpoint "
              f"survives; moving on, the next finalize attempt resumes it")
        return False
    if p.returncode != 0:
        print(f"  FAILED {r['run']}: exit {p.returncode}")
    return p.returncode == 0


def repair(rows: list[dict], deadline: float) -> int:
    todo = [r for r in rows if not r["complete"]]
    if not todo:
        print("  nothing to repair — every expected run is complete")
        return 0
    fixed = 0
    for r in todo:
        if time.time() >= deadline:
            print(f"  STOP: repair budget exhausted with "
                  f"{len(todo)-fixed} run(s) still short")
            break
        if repair_one(r, deadline):
            fixed += 1
    return fixed


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def build_deliverables() -> tuple[int, int]:
    """Figures first, then notebooks: the figures are the deliverable, and the
    notebooks are long. If the job dies between them the figures already exist."""
    fig_rc, nb_rc = 1, 1
    try:
        import build_sv_figure
        fig_rc = build_sv_figure.main()
    except Exception as exc:                                   # noqa: BLE001
        print(f"  figure stage raised {type(exc).__name__}: {exc}")
    try:
        import build_run_notebooks
        nb_rc = build_run_notebooks.main()
    except Exception as exc:                                   # noqa: BLE001
        print(f"  notebook stage raised {type(exc).__name__}: {exc}")
    return fig_rc, nb_rc


def write_report(rows: list[dict], attempt: int, max_attempt: int,
                 fixed: int, fig_rc: int, nb_rc: int) -> Path:
    """The on-disk record. Written on EVERY path, including total failure —
    it is what makes the campaign readable without a live session or email."""
    prod = [r for r in rows if r["kind"] == "production"]
    done = sum(r["complete"] for r in prod)
    lines = [
        "# sigma56_sv — campaign report",
        "",
        f"Generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC} "
        f"(finalize attempt {attempt}/{max_attempt}, SLURM job "
        f"{os.environ.get('SLURM_JOB_ID','n/a')})",
        "",
        f"**Production runs complete: {done}/{len(prod)}.** "
        f"Repaired this attempt: {fixed}.",
        "",
        "Plan: `docs/plans/sigma56-sv-twin.md` · "
        "Handover: `docs/handovers/sigma56-sv-twin.md`",
        "",
        "## Run status", "",
        "| run | kind | steps | state |", "|---|---|---|---|",
    ]
    for r in rows:
        state = ("complete" if r["complete"] else
                 "MISSING" if not r["exists"] else "SHORT")
        lines.append(f"| `{r['run']}` | {r['kind']} | "
                     f"{r['steps_done']}/{r['steps_target']} | {state} |")

    lines += ["", "## Stopping power", ""]
    try:
        t = S.table()
        if t.empty:
            lines.append("_No run has produced observables yet._")
        else:
            ok = t[t.complete & t.cap]
            if ok.empty:
                lines.append("_No COMPLETE production point yet._")
            else:
                piv = ok.pivot_table(index=["sigma_wp", "half"], columns="v",
                                     values="S_eV_per_Bohr")
                lines += ["S = [E_total(t_f) − E_GS] / 25 Bohr, in eV/Bohr:", "",
                          "```", piv.round(3).to_string(), "```"]
            excluded = t[~t.complete]
            if not excluded.empty:
                lines += ["", "Excluded as incomplete: "
                          + ", ".join(f"`{r.run}`" for r in excluded.itertuples())]
        c = S.cap_cost()
        if not c.empty:
            lines += ["", "CAP cost on the classical half (v = 3.0):", "",
                      "```", c.round(4).to_string(index=False), "```"]
    except Exception as exc:                                   # noqa: BLE001
        lines.append(f"_Analysis unavailable: {type(exc).__name__}: {exc}_")

    lines += ["", "## Deliverables", "",
              f"- figure stage exit code: {fig_rc}",
              f"- notebook stage exit code: {nb_rc}", ""]
    for f in sorted(HERE.glob("*.png")) + sorted(HERE.glob("*.csv")):
        lines.append(f"- `{f.name}`")

    if done < len(prod) and attempt >= max_attempt:
        lines += ["", "## ⚠ Incomplete at the final attempt", "",
                  "The repair budget ran out with runs still short. They are "
                  "checkpointed, so nothing is lost — resume with:", "",
                  "```",
                  "sbatch --export=ALL,LJ_SIGMA=<s>,LJ_RESUME=1 "
                  "shared/bin/run-s56-wp.slurm <idx>",
                  "sbatch shared/bin/run-s56-finalize.slurm 1 2",
                  "```"]

    out = HERE / "CAMPAIGN_REPORT.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {out}")
    return out


def notify(report: Path, rows: list[dict]) -> None:
    """Best-effort email. NEVER raises: a missing credential must not turn a
    finished campaign into a failed job. The report on disk is the guaranteed
    channel; email is the convenience one."""
    prod = [r for r in rows if r["kind"] == "production"]
    done = sum(r["complete"] for r in prod)
    figs = [p for p in (HERE / "S_of_v_sigma56.png",
                        HERE / "S_of_sigma_eq.png") if p.exists()]
    subject = (f"[sigma56_sv] {done}/{len(prod)} production runs complete — "
               f"S(v) twins at sigma_WP = 5 and 6")
    body = (
        "HYPOTHESIS\n"
        "  A wavepacket disperses in flight, so a classical projectile of fixed\n"
        "  width is only a fair comparison when the packet's sigma label agrees\n"
        "  with its time-average. sigma_WP = 5 and 6 are the first widths where\n"
        "  it does (growth x1.23 / x1.12 over the in-slab transit). At what width\n"
        "  do the classical and quantum stopping powers coincide?\n\n"
        "WHAT WAS DONE\n"
        f"  {len(prod)} production runs (sigma 5 and 6 x v 2.0-3.5 x "
        "{classical, wavepacket}), matched pair by pair: identical ground state,\n"
        "  35x35x105 Bohr box, dx 0.40, dt 0.04, launch z = -27.5, CAP eta = -1 Ha\n"
        "  on BOTH halves. The classical twin uses a direct erf/r potential of\n"
        "  sigma_pot = sigma_WP/sqrt(2). Plus CAP-free controls and vacuum baselines.\n"
        f"  {done}/{len(prod)} are complete at the time of writing.\n\n"
        "WHAT THE PLOTS SHOW\n"
        "  S_of_v_sigma56.png  — S = [E_total(t_f) - E_GS]/25 Bohr against v, the\n"
        "    new twin pairs plus the existing sigma = 0.5/2/3 wavepacket traces.\n"
        "  S_of_sigma_eq.png   — the same points against the TIME-AVERAGED width\n"
        "    sqrt(2)<sigma_d>. The collapse test: a sigma = 6 run at v = 2.0 and\n"
        "    the old sigma = 2 run at v = 2.0 have the same time-averaged width\n"
        "    (6.45 vs 6.35) reached completely differently.\n\n"
        "CONCLUSION\n"
        "  See the attached report for the S table and the CAP-cost calibration.\n"
        "  Caveat that travels with these numbers: sigma = 5/6 ran at L_z = 105\n"
        "  and the legacy sigma = 0.5/2/3 traces at L_z = 85; the sigma = 0.5\n"
        "  classical benchmark is CAP-FREE and therefore a DIFFERENT estimator.\n\n"
        f"  Report: {report}\n"
        "  Plan: docs/plans/sigma56-sv-twin.md\n"
        "  Handover: docs/handovers/sigma56-sv-twin.md\n"
    )
    try:
        from inqview.email import send_run_email
        send_run_email(subject=subject, body=body,
                       attachments=[str(report), *[str(f) for f in figs]])
        print("  emailed the campaign report")
    except Exception as exc:                                   # noqa: BLE001
        print(f"  EMAIL SKIPPED ({type(exc).__name__}: {exc})")
        print("  The report on disk is complete and is the authoritative record.")
        print("  To enable email for future runs: python -m inqview.email setup")


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--attempt", type=int, default=1)
    ap.add_argument("--max", type=int, default=2)
    ap.add_argument("--status-only", action="store_true")
    ap.add_argument("--no-repair", action="store_true")
    a = ap.parse_args()

    print(f"=== sigma56_sv finalize — attempt {a.attempt}/{a.max} ===\n")
    rows = status()
    print_status(rows)

    if a.status_only:
        return 0

    fixed = 0
    if not a.no_repair:
        # Budget: whichever is smaller, the configured repair budget or what is
        # left of this job after reserving time to build and report.
        deadline = time.time() + REPAIR_BUDGET_S
        print(f"\n=== repair (budget {REPAIR_BUDGET_S/3600:.0f} h, "
              f"reserving {REPORT_RESERVE_S/3600:.0f} h to report) ===")
        fixed = repair(rows, deadline)
        rows = status()          # re-read: repairs change the picture
        print()
        print_status(rows)

    print("\n=== deliverables ===")
    fig_rc, nb_rc = build_deliverables()

    report = write_report(rows, a.attempt, a.max, fixed, fig_rc, nb_rc)
    notify(report, rows)

    prod_done = all(r["complete"] for r in rows if r["kind"] == "production")
    print(f"\nproduction complete: {prod_done}")
    # Exit 0 whenever the report was written. A non-zero code here would only
    # make the NEXT chained attempt look like it should not run — and the next
    # attempt is exactly what an incomplete campaign needs.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
