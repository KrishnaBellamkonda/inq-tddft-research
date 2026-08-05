"""
Pilot gate for lz_bulk_sweep: the automated check between the v = 3.0 pilot
(all four boxes, both halves, ~20 GPU-h) and the ~70 GPU-h of production.

    python pilot_gate.py            # exit 0 -> SLURM releases production
    python pilot_gate.py --report   # rebuild the report, never gate (exit 0)

User decision (2026-08-05): "do one or two runs of velocity for all the Lz.
Check if everything is alright before committing the massive number of GPU
hours." The check is AUTOMATED so the campaign stays autonomous: production
arrays hang on afterok of this job, so a non-zero exit here leaves them
DependencyNeverSatisfied — blocked, reported, never silent.

HARD gates — correctness only (checkpoint-dont-block):
  1. every pilot run exists, run_completed = true, steps_done >= target;
  2. pairwise-ledger closure against INQ's own scalars at the last common step
     (wp: e_hartree_check/e_external_check; classical: E_SS vs energy_hartree
     and E_SB + E_PS vs energy_external) — tolerance 1e-5 Ha against the
     5e-10 Ha measured on sigma56 data;
  3. S finite, and S_deposit not significantly negative.

WARN — reported (and emailed), never gating:
  * GS interior bulk-likeness: n(z~0) from the GS density VTI vs n0 (the
    Friedel-overlap concern for the 15-Bohr slab; plan section "caveats");
  * S(L) ordering vs the L = 25 anchors;
  * production cost re-projection from the MEASURED pilot s/step.
"""
from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import lzb_stopping as L   # noqa: E402

PILOT_V = L.PILOT_V
PROD_V = tuple(v for v in L.VELOCITIES if v != PILOT_V)
CLOSURE_TOL_HA = 1.0e-5
BULK_TOL_PC = 2.0          # WARN threshold on |n(0)/n0 - 1|

# rough per-step estimates, replaced by measured values when a pilot summary
# carries wall_time_s (used only for the cost WARN, never a gate)
S_PER_STEP_EST = {"s0p5_L15": 1.3, "s0p5_L35": 4.0, "s5p0_L15": 1.7, "s5p0_L35": 4.8}


def closure_residuals(cfg: str, half: str) -> dict:
    """Ledger-vs-INQ closure at the last common step, in Ha."""
    obs_dir = L.run_dir(cfg, PILOT_V, half) / "raw" / "observables"
    obs = L._concat(obs_dir, "observables")
    ix = L._concat(obs_dir, "interactions")
    m = pd.merge(obs, ix, on="step", suffixes=("", "_ix"))
    if m.empty:
        raise ValueError("no common steps between observables and interactions")
    r = m.iloc[-1]
    if half == "wp":
        return {"step": int(r.step),
                "hartree": abs(float(r.e_hartree_check) - float(r.energy_hartree)),
                "external": abs(float(r.e_external_check) - float(r.energy_external))}
    return {"step": int(r.step),
            "hartree": abs(float(r.e_ss) - float(r.energy_hartree)),
            "external": abs(float(r.e_sb) + float(r.e_ps) - float(r.energy_external))}


def check_pilot_run(cfg: str, half: str) -> tuple[list[str], list[str], dict]:
    """-> (failures, notes, info) for one pilot run. Empty failures = PASS."""
    fails, notes, info = [], [], {}
    name = L.run_name(cfg, PILOT_V, half)
    d = L.run_dir(cfg, PILOT_V, half)
    if not d.exists():
        return [f"{name}: run directory missing"], notes, info

    kv = L.summary_kv(d / "run_summary.txt")
    if kv.get("run_completed") != "true":
        fails.append(f"{name}: run_summary.txt missing or run_completed != true")

    try:
        p = L.measure(cfg, PILOT_V, half)
        info["point"] = p
        if not p.complete:
            fails.append(f"{name}: SHORT ({p.steps_done}/{p.steps_target} steps)")
        if not np.isfinite(p.S_deposit_eV_per_Bohr):
            fails.append(f"{name}: S_deposit is not finite")
        elif p.S_deposit_eV_per_Bohr < -0.05:
            fails.append(f"{name}: S_deposit = {p.S_deposit_eV_per_Bohr:.3f} eV/Bohr "
                         "(significantly negative — the deposit reference or the "
                         "E_PS correction is broken)")
        if not p.settled:
            notes.append(f"{name}: plateau not settled "
                         f"(drift {p.plateau_drift_eV:.3f} eV over the last 10 %)")
    except Exception as e:                                        # noqa: BLE001
        fails.append(f"{name}: measurement failed ({type(e).__name__}: {e})")
        return fails, notes, info

    try:
        c = closure_residuals(cfg, half)
        info["closure"] = c
        for k in ("hartree", "external"):
            if c[k] > CLOSURE_TOL_HA:
                fails.append(f"{name}: {k} closure residual {c[k]:.2e} Ha "
                             f"> {CLOSURE_TOL_HA:.0e} (step {c['step']})")
    except Exception as e:                                        # noqa: BLE001
        fails.append(f"{name}: closure unreadable ({type(e).__name__}: {e})")

    return fails, notes, info


def gs_bulk_likeness() -> list[str]:
    """WARN lines: interior density vs n0 per box, from the GS density VTI."""
    out = []
    try:
        from inqview import load_vti
    except Exception as e:                                        # noqa: BLE001
        return [f"bulk-likeness check unavailable (inqview import: {e})"]
    for cfg, b in L.CFGS.items():
        vti_dir = L.SCRIPTS / "gs" / "results" / cfg / "density_gs"
        vtis = sorted(vti_dir.glob("*.vti")) if vti_dir.exists() else []
        if not vtis:
            out.append(f"{cfg}: no GS density VTI at {vti_dir} — check skipped")
            continue
        try:
            f = load_vti(vtis[0])
            prof = f.data.mean(axis=(0, 1))                  # n(z), x/y-averaged
            n0 = b.n_e / (L.LXY * L.LXY * b.l_slab)
            iz0 = int(np.argmin(np.abs(f.z)))
            dev = 100.0 * (prof[iz0] / n0 - 1.0)
            core = np.abs(f.z) < max(b.half - 2.0, 1.0)
            ripple = 100.0 * (prof[core].max() - prof[core].min()) / n0
            flag = "WARN" if abs(dev) > BULK_TOL_PC else "ok"
            out.append(f"{cfg}: n(z=0)/n0 - 1 = {dev:+.2f} % [{flag}], "
                       f"interior peak-to-peak {ripple:.1f} % of n0")
        except Exception as e:                                    # noqa: BLE001
            out.append(f"{cfg}: bulk-likeness read failed ({type(e).__name__}: {e})")
    return out


def cost_projection(points: dict) -> tuple[list[str], float]:
    """Projected production cost from MEASURED pilot s/step where available."""
    lines, total_h = [], 0.0
    for cfg in L.CFGS:
        for half in ("wp", "classical"):
            kv = L.summary_kv(L.run_dir(cfg, PILOT_V, half) / "run_summary.txt")
            steps = L.STEPS_TARGET[cfg][PILOT_V]
            try:
                sps = float(kv["wall_time_s"]) / steps
                src = "measured"
            except (KeyError, ValueError, ZeroDivisionError):
                sps, src = S_PER_STEP_EST[cfg], "estimate"
            prod_steps = sum(L.STEPS_TARGET[cfg][v] for v in PROD_V)
            h = prod_steps * sps / 3600.0
            total_h += h
            lines.append(f"{cfg}/{half}: {sps:.2f} s/step ({src}) -> "
                         f"{prod_steps} production steps = {h:.1f} h")
    lines.append(f"TOTAL projected production: {total_h:.0f} GPU-h "
                 "(+ vacuum baselines; proceeding is the default — "
                 "kill with scancel if unwanted)")
    return lines, total_h


def ordering_info(points: dict) -> list[str]:
    """S(L15) vs anchor(L25) vs S(L35), per (sigma, half) at the pilot v."""
    out = []
    for sigma in (0.5, 5.0):
        stag = "s0p5" if sigma == 0.5 else "s5p0"
        for half in ("wp", "classical"):
            s15 = points.get((f"{stag}_L15", half))
            s35 = points.get((f"{stag}_L35", half))
            a25 = L.anchor_S(sigma, PILOT_V, half)
            v15 = s15.S_deposit_eV_per_Bohr if s15 else float("nan")
            v35 = s35.S_deposit_eV_per_Bohr if s35 else float("nan")
            line = (f"sigma={sigma} {half}: S(15)={v15:.3f}  "
                    f"S(25,anchor)={a25:.3f}  S(35)={v35:.3f} eV/Bohr")
            if np.isfinite(v15) and np.isfinite(v35) and np.isfinite(a25):
                mono = (v15 >= a25 >= v35) or (v15 <= a25 <= v35)
                line += "  [monotone in L]" if mono else \
                        "  [NOT monotone — thin-slab effects or noise; inspect]"
            out.append(line)
    return out


def write_report(all_fails, all_notes, warn_bulk, warn_cost, warn_order,
                 points, gated: bool) -> Path:
    ok = not all_fails
    lines = [
        "# lz_bulk_sweep — PILOT REPORT (v = 3.0, all four boxes, both halves)",
        "",
        f"Generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC} "
        f"(SLURM job {os.environ.get('SLURM_JOB_ID', 'n/a')})",
        "",
        f"## VERDICT: {'PASS — production released' if ok else 'FAIL — production BLOCKED'}"
        + ("" if gated else " (report-only mode, nothing gated)"),
        "",
    ]
    if all_fails:
        lines += ["### Hard-gate failures", ""] + [f"- {f}" for f in all_fails] + [""]
    if all_notes:
        lines += ["### Notes (non-gating)", ""] + [f"- {n}" for n in all_notes] + [""]

    lines += ["## Pilot S values (corrected deposit, eV/Bohr)", ""]
    if points:
        lines += ["| run | S_deposit | S (no E_PS cut) | norm_final | settled | steps |",
                  "|---|---|---|---|---|---|"]
        for (cfg, half), p in sorted(points.items()):
            lines.append(f"| `{p.run}` | {p.S_deposit_eV_per_Bohr:.3f} | "
                         f"{p.S_eV_per_Bohr:.3f} | {p.norm_final:.2e} | "
                         f"{p.settled} | {p.steps_done}/{p.steps_target} |")
    else:
        lines.append("_no pilot point measurable_")

    lines += ["", "## S(L) ordering vs the L = 25 anchors (INFO)", ""]
    lines += [f"- {x}" for x in warn_order]
    lines += ["", "## GS interior bulk-likeness (WARN-only)", ""]
    lines += [f"- {x}" for x in warn_bulk]
    lines += ["", "## Production cost projection (WARN, never a gate)", ""]
    lines += [f"- {x}" for x in warn_cost]
    lines += ["", "Plan: `docs/plans/jellium-slab-extend-Lz.md` · "
                  "Handover: `docs/handovers/jellium-slab-extend-Lz.md`", ""]

    out = HERE / "PILOT_REPORT.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {out}")
    return out


def notify(report: Path, ok: bool, total_h: float) -> None:
    """Best-effort email — never raises (no Gmail credentials on this machine
    turns into a log line, not a failed gate)."""
    subject = (f"[lz_bulk_sweep] pilot {'PASSED — production released' if ok else 'FAILED — production BLOCKED'}"
               f" (~{total_h:.0f} GPU-h at stake)")
    body = (
        "HYPOTHESIS\n"
        "  Slab deposit stopping separates as S(L) = S_bulk + c/L. The pilot\n"
        "  (v = 3.0 at L_slab = 15/25/35, sigma_WP = 0.5 and 5, classical and\n"
        "  wavepacket) tests the machinery before the full velocity grid runs.\n\n"
        "WHAT WAS DONE\n"
        "  8 pilot runs, per-family matched geometry (standoff 11.5 / 15 Bohr),\n"
        "  corrected deposit estimator, automated gate on completeness, ledger\n"
        "  closure (<= 1e-5 Ha) and S sanity. Details in PILOT_REPORT.md.\n\n"
        "WHAT THE REPORT SHOWS\n"
        "  Pilot S table, S(L) ordering vs the L = 25 anchors, GS interior\n"
        "  bulk-likeness per box, and the measured-cost projection for the\n"
        "  remaining velocities.\n\n"
        "CONCLUSION\n"
        f"  {'All hard gates passed; SLURM has released production (v = 2.0/2.5/3.5).' if ok else 'A hard gate failed; production is DependencyNeverSatisfied until resubmitted.'}\n\n"
        f"  Report: {report}\n"
        "  Plan: docs/plans/jellium-slab-extend-Lz.md\n"
    )
    try:
        from inqview.email import send_run_email
        send_run_email(subject=subject, body=body, attachments=[str(report)])
        print("  emailed the pilot report")
    except Exception as exc:                                     # noqa: BLE001
        print(f"  EMAIL SKIPPED ({type(exc).__name__}: {exc}) — the report on "
              f"disk is the authoritative record.")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true",
                    help="rebuild the report only; always exit 0")
    a = ap.parse_args()

    print("=== lz_bulk_sweep pilot gate ===\n")
    all_fails, all_notes, points = [], [], {}
    for cfg in L.CFGS:
        for half in ("wp", "classical"):
            fails, notes, info = check_pilot_run(cfg, half)
            all_fails += fails
            all_notes += notes
            if "point" in info:
                points[(cfg, half)] = info["point"]
            tag = "PASS" if not fails else "FAIL"
            print(f"  [{tag}] {L.run_name(cfg, PILOT_V, half)}"
                  + (f" — {'; '.join(fails)}" if fails else ""))
            if "closure" in info:
                c = info["closure"]
                print(f"         closure: hartree {c['hartree']:.2e} Ha, "
                      f"external {c['external']:.2e} Ha (step {c['step']})")

    print("\n--- WARN-only checks ---")
    warn_bulk = gs_bulk_likeness()
    warn_cost, total_h = cost_projection(points)
    warn_order = ordering_info(points)
    for x in warn_bulk + warn_order + warn_cost:
        print(f"  {x}")

    ok = not all_fails
    report = write_report(all_fails, all_notes, warn_bulk, warn_cost, warn_order,
                          points, gated=not a.report)
    notify(report, ok, total_h)

    if a.report:
        return 0
    print(f"\nverdict: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
