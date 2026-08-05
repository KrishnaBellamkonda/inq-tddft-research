"""
Autonomous finalizer for lz_bulk_sweep: status -> repair -> report.

    python finalize.py --attempt 1 --max 2
    python finalize.py --status-only          # inspect, change nothing

Clone of hypotheses/sigma56_sv/finalize.py (see its docstring for the design
constraints: no sbatch from inside a job, bounded always, report on every path)
adapted to this campaign's run matrix:

    production : 4 boxes x 4 velocities x {wp, classical}  (32 runs)
    vacuum     : 4 boxes x 4 velocities                    (16 CAP-only baselines)

The per-run repair timeout cap (2x measured cost) is inherited from the sigma56
fix — one wedged `electrons.load` can no longer swallow the whole allocation.
Vacuum repairs are per-BOX (the vac dispatcher runs its four velocities in one
job), deduplicated so four short vac runs of one box trigger one repair.
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

import lzb_stopping as L   # noqa: E402

REPO = L.REPO
BIN = REPO / "shared/bin"
V_INDEX = {v: i for i, v in enumerate(L.VELOCITIES)}

REPORT_RESERVE_S = int(os.environ.get("LZB_REPORT_RESERVE_S", 3 * 3600))
REPAIR_BUDGET_S = int(os.environ.get("LZB_REPAIR_BUDGET_S", 30 * 3600))


# ---------------------------------------------------------------------------
# status
# ---------------------------------------------------------------------------
def expected() -> list[dict]:
    """Every run the campaign should produce, in repair-priority order:
    production first, vacuum baselines last (they only annotate the result)."""
    out = []
    for cfg in L.CFGS:
        for v in L.VELOCITIES:
            for half in ("wp", "classical"):
                out.append({"cfg": cfg, "v": v, "half": half, "kind": "production"})
    for cfg in L.CFGS:
        for v in L.VELOCITIES:
            out.append({"cfg": cfg, "v": v, "half": "vac", "kind": "vacuum"})
    return out


def _steps_done(d: Path) -> int:
    obs = d / "raw" / "observables"
    if not obs.exists():
        return -1
    try:
        return int(L._concat(obs, "observables")["step"].iloc[-1])
    except Exception:                                          # noqa: BLE001
        # vacuum runs write no observables.csv
        try:
            return int(L._concat(obs, "wp_momentum_stats")["step"].iloc[-1])
        except Exception:                                      # noqa: BLE001
            return -1


def status() -> list[dict]:
    rows = []
    for e in expected():
        d = L.run_dir(e["cfg"], e["v"], e["half"])
        done = _steps_done(d)
        target = L.STEPS_TARGET[e["cfg"]][e["v"]]
        rows.append({**e, "run": d.name, "dir": d, "exists": d.exists(),
                     "steps_done": done, "steps_target": target,
                     "complete": done >= target,
                     "resumable": (d / "checkpoint").exists()})
    return rows


def print_status(rows: list[dict]) -> None:
    print(f"{'run':30s} {'kind':11s} {'steps':>13s}  state")
    for r in rows:
        state = ("COMPLETE" if r["complete"] else
                 "MISSING" if not r["exists"] else
                 f"SHORT ({'resumable' if r['resumable'] else 'no ckpt'})")
        print(f"  {r['run']:28s} {r['kind']:11s} "
              f"{r['steps_done']:>6d}/{r['steps_target']:<6d} {state}")
    n_ok = sum(r["complete"] for r in rows)
    print(f"\n  {n_ok}/{len(rows)} complete "
          f"({sum(r['complete'] for r in rows if r['kind']=='production')}"
          f"/{sum(1 for r in rows if r['kind']=='production')} production)")


# ---------------------------------------------------------------------------
# repair
# ---------------------------------------------------------------------------
def dispatcher_for(half: str) -> Path:
    return BIN / {"wp": "run-lzb-wp.slurm",
                  "classical": "run-lzb-cl.slurm",
                  "vac": "run-lzb-vac.slurm"}[half]


def binary_dir(half: str) -> Path:
    return L.SCRIPTS / {"wp": "wp", "classical": "classical", "vac": "vac"}[half]


def ensure_binary(half: str, cfg: str, deadline: float) -> bool:
    """Build the run binary if missing, via the dispatcher's own smoke stage.
    A smoke that builds but fails its t=0 gates is a DEFECT, not a hiccup —
    refuse to repair production from a packet the gates rejected."""
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
    env["LZB_CFG"] = cfg
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
    if p.returncode != 0 and ok:
        print(f"  {half} binary built but its t=0 gates FAILED — refusing to "
              f"repair production runs from a packet the gates rejected")
        return False
    return ok


# Per-box cost estimates (s/step, either half; vac is far cheaper so these are
# safe upper bounds for its cap too). Used only to BOUND a repair.
S_PER_STEP = {"s0p5_L15": 1.3, "s5p0_L15": 1.7, "s0p5_L35": 4.0, "s5p0_L35": 4.8}
STARTUP_S = 900.0
SAFETY = 2.0


def run_budget_s(r: dict) -> float:
    todo = r["steps_target"] - (r["steps_done"] if r["resumable"] else 0)
    todo = max(todo, 0)
    if r["half"] == "vac":
        # the vac dispatcher redoes all four velocities of the box
        todo = sum(L.STEPS_TARGET[r["cfg"]].values())
    return SAFETY * (STARTUP_S + todo * S_PER_STEP.get(r["cfg"], 4.8))


def repair_one(r: dict, deadline: float) -> bool:
    script = dispatcher_for(r["half"])
    env = dict(os.environ)
    env["LZB_CFG"] = r["cfg"]
    env["SLURM_SUBMIT_DIR"] = str(REPO)
    env.pop("SLURM_ARRAY_TASK_ID", None)

    if r["half"] == "vac":
        args = [str(script)]            # redoes the box's four baselines
    else:
        env["LJ_RESUME"] = "1" if r["resumable"] else "0"
        args = [str(script), str(V_INDEX[r["v"]])]

    remaining = deadline - time.time()
    if remaining <= 0:
        print(f"  SKIP {r['run']}: repair budget exhausted")
        return False
    if not ensure_binary(r["half"], r["cfg"], deadline):
        print(f"  SKIP {r['run']}: no usable binary for the {r['half']} half")
        return False
    remaining = deadline - time.time()
    # PER-RUN CAP (the sigma56 lesson): 2x measured cost, then move on — the
    # checkpoint survives and the next attempt resumes it.
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
    vac_done: set[str] = set()
    for r in todo:
        if time.time() >= deadline:
            print(f"  STOP: repair budget exhausted with "
                  f"{len(todo)-fixed} run(s) still short")
            break
        if r["half"] == "vac":
            if r["cfg"] in vac_done:
                continue                # one vac repair covers the whole box
            vac_done.add(r["cfg"])
        if repair_one(r, deadline):
            fixed += 1
    return fixed


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def build_deliverables() -> tuple[int, int]:
    fig_rc, nb_rc = 1, 1
    try:
        import build_lzb_figures
        fig_rc = build_lzb_figures.main()
    except Exception as exc:                                   # noqa: BLE001
        print(f"  figure stage raised {type(exc).__name__}: {exc}")
    # Run notebooks (density-GIF battery) are built by a follow-up session per
    # the notebook rule; the finalizer calls the builder IF it exists so adding
    # it later needs no dispatcher change.
    try:
        import build_run_notebooks                            # noqa: F401
        nb_rc = build_run_notebooks.main()
    except ImportError:
        print("  notebook builder not present yet (planned follow-up) — skipped")
        nb_rc = 0
    except Exception as exc:                                   # noqa: BLE001
        print(f"  notebook stage raised {type(exc).__name__}: {exc}")
    return fig_rc, nb_rc


def write_report(rows: list[dict], attempt: int, max_attempt: int,
                 fixed: int, fig_rc: int, nb_rc: int) -> Path:
    prod = [r for r in rows if r["kind"] == "production"]
    done = sum(r["complete"] for r in prod)
    lines = [
        "# lz_bulk_sweep — campaign report",
        "",
        f"Generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC} "
        f"(finalize attempt {attempt}/{max_attempt}, SLURM job "
        f"{os.environ.get('SLURM_JOB_ID','n/a')})",
        "",
        f"**Production runs complete: {done}/{len(prod)}.** "
        f"Repaired this attempt: {fixed}.",
        "",
        "Plan: `docs/plans/jellium-slab-extend-Lz.md` · "
        "Handover: `docs/handovers/jellium-slab-extend-Lz.md`",
        "",
        "## Run status", "",
        "| run | kind | steps | state |", "|---|---|---|---|",
    ]
    for r in rows:
        state = ("complete" if r["complete"] else
                 "MISSING" if not r["exists"] else "SHORT")
        lines.append(f"| `{r['run']}` | {r['kind']} | "
                     f"{r['steps_done']}/{r['steps_target']} | {state} |")

    lines += ["", "## Stopping power (corrected deposit, eV/Bohr)", ""]
    try:
        t = L.table()
        if t.empty:
            lines.append("_No run has produced observables yet._")
        else:
            ok = t[t.complete]
            if ok.empty:
                lines.append("_No COMPLETE production point yet._")
            else:
                piv = ok.pivot_table(index=["sigma_wp", "half", "L_slab"],
                                     columns="v",
                                     values="S_deposit_eV_per_Bohr")
                lines += ["S = [E_total(t_f) − E_GS − E_PS(t_f)] / L_slab:", "",
                          "```", piv.round(3).to_string(), "```"]
            excluded = t[~t.complete]
            if not excluded.empty:
                lines += ["", "Excluded as incomplete: "
                          + ", ".join(f"`{r.run}`" for r in excluded.itertuples())]
        fits = HERE / "lzb_fits.csv"
        if fits.exists():
            import pandas as pd
            lines += ["", "1/L fits (S_bulk = intercept at 1/L -> 0):", "",
                      "```", pd.read_csv(fits).round(3).to_string(index=False), "```"]
    except Exception as exc:                                   # noqa: BLE001
        lines.append(f"_Analysis unavailable: {type(exc).__name__}: {exc}_")

    lines += ["", "## Deliverables", "",
              f"- figure stage exit code: {fig_rc}",
              f"- notebook stage exit code: {nb_rc}", ""]
    for f in sorted(HERE.glob("*.png")) + sorted(HERE.glob("*.csv")):
        lines.append(f"- `{f.name}`")

    if done < len(prod) and attempt >= max_attempt:
        lines += ["", "## ⚠ Incomplete at the final attempt", "",
                  "Short runs are checkpointed — resume with:", "",
                  "```",
                  "sbatch --export=ALL,LZB_CFG=<preset>,LJ_RESUME=1 "
                  "shared/bin/run-lzb-wp.slurm <idx>",
                  "sbatch shared/bin/run-lzb-finalize.slurm 1 2",
                  "```"]

    out = HERE / "CAMPAIGN_REPORT.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {out}")
    return out


def notify(report: Path, rows: list[dict]) -> None:
    """Best-effort email; NEVER raises."""
    prod = [r for r in rows if r["kind"] == "production"]
    done = sum(r["complete"] for r in prod)
    figs = [p for p in (HERE / "S_of_invL.png",) if p.exists()]
    subject = (f"[lz_bulk_sweep] {done}/{len(prod)} production runs complete — "
               f"S(L) slab->bulk extrapolation")
    body = (
        "HYPOTHESIS\n"
        "  The slab deposit separates as S(L) = S_bulk + c/L: surface and\n"
        "  per-traversal terms (incl. the WP self-interaction overhead) scale\n"
        "  as 1/L, so extrapolating L_slab = 15/25/35 to 1/L -> 0 yields the\n"
        "  bulk stopping power — and tests whether the classical/WP gap\n"
        "  (1.9-3.4x at sigma = 6, L = 25) closes in the bulk limit.\n\n"
        "WHAT WAS DONE\n"
        f"  {len(prod)} production runs: L_slab 15/35 x sigma_WP {{0.5, 5}} x\n"
        "  v 2.0-3.5 x {classical, wavepacket}, per-family geometry matched to\n"
        "  the existing L = 25 anchors (standoff 11.5 / 15 Bohr), corrected\n"
        "  deposit estimator on both halves, pilot-gated chain.\n"
        f"  {done}/{len(prod)} complete at the time of writing.\n\n"
        "WHAT THE PLOT SHOWS\n"
        "  S_of_invL.png — S_deposit against 1/L_slab per sigma family with the\n"
        "  linear 1/L fits; the starred intercepts at 1/L = 0 are S_bulk(v).\n"
        "  The middle (L = 25) point's residual is the linearity check.\n\n"
        "CONCLUSION\n"
        "  See CAMPAIGN_REPORT.md for the S table, the fits and the run\n"
        "  status. Caveats that travel with the numbers: the sigma = 0.5 WP\n"
        "  trace is qualitative (dispersing instrument), and L = 15 holds only\n"
        "  ~2 Friedel periods (see PILOT_REPORT.md bulk-likeness).\n\n"
        f"  Report: {report}\n"
        "  Plan: docs/plans/jellium-slab-extend-Lz.md\n"
    )
    try:
        from inqview.email import send_run_email
        send_run_email(subject=subject, body=body,
                       attachments=[str(report), *[str(f) for f in figs]])
        print("  emailed the campaign report")
    except Exception as exc:                                   # noqa: BLE001
        print(f"  EMAIL SKIPPED ({type(exc).__name__}: {exc})")
        print("  The report on disk is complete and is the authoritative record.")


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--attempt", type=int, default=1)
    ap.add_argument("--max", type=int, default=2)
    ap.add_argument("--status-only", action="store_true")
    ap.add_argument("--no-repair", action="store_true")
    a = ap.parse_args()

    print(f"=== lz_bulk_sweep finalize — attempt {a.attempt}/{a.max} ===\n")
    rows = status()
    print_status(rows)

    if a.status_only:
        return 0

    fixed = 0
    if not a.no_repair:
        deadline = time.time() + REPAIR_BUDGET_S
        print(f"\n=== repair (budget {REPAIR_BUDGET_S/3600:.0f} h, "
              f"reserving {REPORT_RESERVE_S/3600:.0f} h to report) ===")
        fixed = repair(rows, deadline)
        rows = status()
        print()
        print_status(rows)

    print("\n=== deliverables ===")
    fig_rc, nb_rc = build_deliverables()

    report = write_report(rows, a.attempt, a.max, fixed, fig_rc, nb_rc)
    notify(report, rows)

    prod_done = all(r["complete"] for r in rows if r["kind"] == "production")
    print(f"\nproduction complete: {prod_done}")
    return 0    # a report was written; the next chained attempt must still run


if __name__ == "__main__":
    raise SystemExit(main())
