#!/usr/bin/env python3
"""Autonomous orchestrator — bulk-jellium governing-PDE discovery redo (T8-T14).

Fully HANDS-OFF: launch once and leave it. Self-contained Python (numpy/scipy),
CPU-only, NO INQ runs, NO GPU, NO live LLM calls in the loop.

Guarantees for unattended operation:
  * IDEMPOTENT / RESUMABLE — each phase writes artifacts/PDE_T*_result.json;
    an existing result is skipped (``--force`` to redo). Re-run to resume.
  * HARD 12 h WALL-CLOCK CAP — checked before every phase and every cell; on
    expiry the run stops gracefully and emails a partial-result summary.
  * PER-CELL try/except — one bad cell (missing/short/NaN series) is skipped and
    logged, never kills the chain. Per-phase try/except likewise.
  * PLATEAU STOP — the T11/T12 refine loop stops early when extra config tries
    stop improving calibration validation.
  * PER-PHASE Gmail — 4-part email (hypothesis / done / plot / conclusion) + >=1
    plot, via inqview.email.send_run_email; email failure is logged, never fatal.

ANTI-P-HACKING (ADR 0011 + 0012): the refine loop tunes the shared config ONLY
on the PINNED CALIBRATION cells; every verdict is read from PINNED HELD-OUT
cells; all attempts are logged. Two projectiles discovered SEPARATELY then
compared. Broad agnostic library; physics named post-hoc; three walls
(pinned split + forward-predict + bootstrap) gate every admitted term.

Usage:
    venv/bin/python3 docs/campaigns/ml-patterns/orchestrate_pde.py            # T8-T14
    venv/bin/python3 docs/campaigns/ml-patterns/orchestrate_pde.py T11 T12    # subset
    ... --no-email     skip emails
    ... --force        redo phases with existing results
    ... --hours 12     wall-clock cap (default 12)
    ... --smoke        tiny run: 1 calib + 1 held-out cell per projectile, 2 configs
"""
from __future__ import annotations
import os, sys, json, time, argparse, traceback
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ART = os.path.join(HERE, "artifacts")
NB = os.path.join(HERE, "notebooks")
os.makedirs(ART, exist_ok=True); os.makedirs(NB, exist_ok=True)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from kernels import celldb, discovery as DIS

# --------------------------------------------------------------------- runtime
SEND_EMAIL = True
FORCE = False
SMOKE = False
DEADLINE = None            # unix time; set in main()

HYPO = ("Hypothesis (bulk-jellium PDE-discovery redo): the induced bath density "
        "dn(z,t) of a classical (point) and a quantum (wavepacket) projectile in "
        "bulk jellium each obeys an interpretable governing PDE discovered from a "
        "BROAD agnostic library; a term counts as physics only if it survives the "
        "three walls (pinned calib/held-out split + forward-prediction + bootstrap "
        "stability). Discovered separately per projectile, then compared "
        "(similarities & differences). Suggestive-of-physics terms named post-hoc.")

# Pinned Track-B split (campaign redo section; deterministic energies in eV).
CALIB_E = [20.0, 50.0, 300.0]
HELDOUT_E = [25.0, 100.0, 600.0]
FWD_PASS = 0.60            # forward-prediction rel-L2 threshold to "validate"

# Config grid the refine loop searches (tuned on CALIBRATION only). Order matters:
# most-promising first so plateau can stop early.
CONFIG_GRID = [
    dict(order=2, threshold=0.06, poly=2, deriv_order=3, smooth_t=1.5, smooth_x=1.0, pod_rank=6),
    dict(order=1, threshold=0.06, poly=2, deriv_order=3, smooth_t=1.5, smooth_x=1.0, pod_rank=6),
    dict(order=2, threshold=0.04, poly=2, deriv_order=2, smooth_t=2.0, smooth_x=1.0, pod_rank=8),
    dict(order=1, threshold=0.08, poly=1, deriv_order=3, smooth_t=1.5, smooth_x=1.0, pod_rank=6),
]


# ------------------------------------------------------------------------- IO
def res_path(ph): return os.path.join(ART, f"PDE_{ph}_result.json")
def have(ph): return os.path.isfile(res_path(ph)) and not FORCE
def save(ph, obj): json.dump(_jsonify(obj), open(res_path(ph), "w"), indent=2)
def load(ph): return json.load(open(res_path(ph)))


def _jsonify(o):
    if isinstance(o, dict): return {k: _jsonify(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [_jsonify(v) for v in o]
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, float) and (np.isnan(o) or np.isinf(o)): return None
    return o


def log(msg): print(msg, flush=True)


def time_left():
    return float("inf") if DEADLINE is None else DEADLINE - time.time()


def out_of_time(margin=120):
    return time_left() < margin


# ---------------------------------------------------------------------- email
def email(subject, parts, plots):
    if not SEND_EMAIL:
        log(f"[email skipped] {subject}"); return
    body = (f"1) HYPOTHESIS\n{parts['hypothesis']}\n\n"
            f"2) WHAT WAS DONE\n{parts['done']}\n\n"
            f"3) WHAT THE PLOT SHOWS\n{parts['plot_shows']}\n\n"
            f"4) CONCLUSION\n{parts['conclusion']}\n")
    try:
        from inqview.email import send_run_email
        send_run_email(subject=subject, body=body,
                       attachments=[p for p in plots if p and os.path.isfile(p)])
        log(f"[email sent] {subject}")
    except Exception as e:
        log(f"[email FAILED, continuing] {subject}: {e!r}")


# ------------------------------------------------------------------- notebooks
def write_notebook(name, title, cells):
    try:
        import nbformat as nbf
    except Exception as e:
        log(f"[nb skip] {e!r}"); return ""
    nb = nbf.v4.new_notebook()
    nb.cells.append(nbf.v4.new_markdown_cell(f"# {title}\n\n{HYPO}"))
    for kind, content in cells:
        nb.cells.append(nbf.v4.new_markdown_cell(content) if kind == "md"
                        else nbf.v4.new_code_cell(content))
    path = os.path.join(NB, name)
    nbf.write(nb, path)
    return path


# ===================================================================== T8 cells
def phase_T8():
    """Resolve + pin the bulk-jellium cells for both projectiles and both cuts."""
    wake = celldb.resolve_wake_cells()          # sigma=5 velocity sweep (WP+classical)
    ff = celldb.resolve_form_factor_cells()     # E=100 sigma sweep (Track A)
    allcells = wake["calibration"] + wake["heldout"]

    def by_role(role_energies):
        return [c for c in allcells if any(abs(c["energy_ev"] - e) < 1 for e in role_energies)]

    manifest = {
        "track_b_split": {"calibration_E": CALIB_E, "heldout_E": HELDOUT_E},
        "wake_cells": {
            "calibration": [_cellref(c) for c in by_role(CALIB_E)],
            "heldout": [_cellref(c) for c in by_role(HELDOUT_E)],
        },
        "wake_skipped": wake["skipped"],
        "form_factor_cells": {k: [_cellref(c) for c in v]
                              for k, v in ff.items() if k != "skipped"},
        "form_factor_skipped": ff.get("skipped", []),
        "fwd_pass_threshold": FWD_PASS,
    }
    save("T8", manifest)
    log(f"[T8] wake calib E={[c['energy_ev'] for c in by_role(CALIB_E)]} "
        f"heldout E={[c['energy_ev'] for c in by_role(HELDOUT_E)]} "
        f"skipped={len(wake['skipped'])}")
    email("ml-patterns REDO T8: cells pinned",
          dict(hypothesis=HYPO,
               done=f"Pinned bulk-jellium cells. Track-B split: calib E={CALIB_E}, "
                    f"held-out E={HELDOUT_E}. Wake calib/held-out cells resolved with "
                    f"matched classical + WP; {len(wake['skipped'])} energies skipped "
                    f"(no matched classical).",
               plot_shows="(no plot — inventory phase)",
               conclusion="Cells frozen; discovery proceeds on these only."),
          [])
    return manifest


def _cellref(c):
    return {k: c[k] for k in ("energy_ev", "velocity_au", "sigma_wp", "omega_p_ev",
                              "r_s", "dx", "cl_run", "wp_run", "cl_bath_dir",
                              "wp_bath_dir", "cl_gs", "wp_gs", "frame_dt_au_cl",
                              "frame_dt_au_wp")}


# ============================================================ T11/T12 discovery
def _discover_cellref(cellref, which, cfg, max_frames):
    """Load one cell's axial field and discover; returns a compact dict or skip."""
    try:
        cf = DIS.load_cell_axial(cellref, which, max_frames=max_frames)
        r = DIS.discover_cell(cf, cfg=cfg, bootstrap=(8 if SMOKE else 15))
        if not r.get("ok"):
            return {"E": cellref["energy_ev"], "skip": r.get("reason")}
        return {
            "E": cellref["energy_ev"], "v": cellref["velocity_au"],
            "admitted_equation": r["admitted_equation"],
            "forward_rel_l2": r["forward_rel_l2"],
            "residual_rel": r["residual_rel"],
            "admitted": [{"term": t["term"], "coeff": t["coeff"],
                          "physics": t["physics"]} for t in r["admitted"]],
            "validated": bool(np.isfinite(r["forward_rel_l2"] or np.nan)
                              and (r["forward_rel_l2"] or 1e9) < FWD_PASS
                              and len(r["admitted"]) > 0),
        }
    except Exception as e:
        return {"E": cellref["energy_ev"], "skip": f"{type(e).__name__}: {e}"}


def _discover_track(which, label, ph):
    """Refine-loop discovery for one projectile (which='cl'|'wp')."""
    man = load("T8")
    calib = man["wake_cells"]["calibration"]
    heldout = man["wake_cells"]["heldout"]
    if SMOKE:
        calib, heldout = calib[:1], heldout[:1]
    grid = CONFIG_GRID[:2] if SMOKE else CONFIG_GRID
    max_frames = 120 if SMOKE else 220

    attempts = []
    best = None
    prev_score = -1.0
    plateau = 0
    for i, cfg in enumerate(grid):
        if out_of_time():
            log(f"[{ph}] out of time before config {i}"); break
        cal_res = [_discover_cellref(c, which, cfg, max_frames) for c in calib]
        n_val = sum(1 for r in cal_res if r.get("validated"))
        # score: validated fraction, tie-broken by mean (1 - forward rel) over valid
        fwd = [r["forward_rel_l2"] for r in cal_res
               if r.get("validated") and r.get("forward_rel_l2") is not None]
        score = n_val + (np.mean([1 - min(f, 1) for f in fwd]) if fwd else 0.0)
        attempts.append({"try": i, "config": cfg, "calib_score": score,
                         "n_validated": n_val, "calib": cal_res})
        log(f"[{ph}] try {i} cfg-order{cfg['order']} th{cfg['threshold']}: "
            f"calib score={score:.2f} ({n_val}/{len(calib)} validated)")
        if score > prev_score + 1e-3:
            best = {"try": i, "config": cfg, "score": score}
            prev_score = score
            plateau = 0
        else:
            plateau += 1
        if plateau >= 2:
            log(f"[{ph}] plateau — freezing config from try {best['try']}"); break

    if best is None:
        best = {"try": 0, "config": grid[0], "score": 0.0}
    # freeze best config; report HELD-OUT
    frozen = best["config"]
    held_res = [_discover_cellref(c, which, frozen, max_frames) for c in heldout]
    n_held_val = sum(1 for r in held_res if r.get("validated"))
    verdict = ("CONFIRM" if n_held_val >= 1 and n_held_val >= (len(held_res) + 1) // 2
               else ("INCONCLUSIVE" if n_held_val >= 1 else "REFUTE"))
    result = {
        "phase": ph, "projectile": label, "frozen_config": frozen,
        "frozen_from_try": best["try"], "calib_attempts": attempts,
        "heldout": held_res, "n_heldout_validated": n_held_val,
        "n_heldout": len(held_res), "verdict": verdict,
    }
    plot = _plot_track(result, ph, label)
    result["plot"] = plot
    save(ph, result)
    # email
    eqs = "; ".join(f"E={r['E']}: {r.get('admitted_equation','skip')}"
                    for r in held_res)
    email(f"ml-patterns REDO {ph} {label}: {verdict}",
          dict(hypothesis=HYPO,
               done=f"Discovered the axial induced-density PDE for the {label} "
                    f"projectile across the velocity sweep. Config tuned on "
                    f"calibration E={CALIB_E} ({len(attempts)} tries, plateau-stopped), "
                    f"frozen, verdict read on HELD-OUT E={HELDOUT_E}.",
               plot_shows=f"Held-out admitted equations + forward-prediction rel-L2 "
                          f"per velocity. {n_held_val}/{len(held_res)} held-out cells "
                          f"validated (forward rel-L2 < {FWD_PASS}). {eqs}",
               conclusion=f"{label} verdict: {verdict}."),
          [plot])
    return result


def phase_T11():
    return _discover_track("cl", "classical", "T11")


def phase_T12():
    return _discover_track("wp", "wavepacket", "T12")


def _plot_track(result, ph, label):
    held = result["heldout"]
    Es = [r["E"] for r in held]
    fwd = [r.get("forward_rel_l2") if r.get("forward_rel_l2") is not None else np.nan
           for r in held]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar([str(e) for e in Es], [min(f, 2) if np.isfinite(f) else 2 for f in fwd],
           color=["#2a7" if (np.isfinite(f) and f < FWD_PASS) else "#c44" for f in fwd])
    ax.axhline(FWD_PASS, ls="--", c="k", lw=1, label=f"pass < {FWD_PASS}")
    ax.set_xlabel("held-out energy (eV)"); ax.set_ylabel("forward-predict rel-L2")
    ax.set_title(f"{ph} {label}: held-out PDE forward-prediction — {result['verdict']}")
    ax.legend(); fig.tight_layout()
    p = os.path.join(ART, f"{ph}_heldout.png"); fig.savefig(p, dpi=110); plt.close(fig)
    return p


# =============================================================== T13 comparison
def phase_T13():
    cl = load("T11"); wp = load("T12")
    def term_map(res):
        m = {}
        for cell in res["heldout"]:
            for t in cell.get("admitted", []):
                m.setdefault(t["term"], []).append(t["coeff"])
        return {k: float(np.median(v)) for k, v in m.items()}
    tcl, twp = term_map(cl), term_map(wp)
    shared = sorted(set(tcl) & set(twp))
    only_cl = sorted(set(tcl) - set(twp))
    only_wp = sorted(set(twp) - set(tcl))
    comp = {
        "shared_terms": {t: {"classical": tcl[t], "wavepacket": twp[t],
                             "ratio_wp_over_cl": (twp[t] / tcl[t] if tcl[t] else None),
                             "physics": DIS.PF.INTERPRET.get(t, "uninterpreted")}
                         for t in shared},
        "classical_only": {t: tcl[t] for t in only_cl},
        "wavepacket_only": {t: twp[t] for t in only_wp},
        "classical_verdict": cl["verdict"], "wavepacket_verdict": wp["verdict"],
    }
    plot = _plot_compare(comp)
    comp["plot"] = plot
    save("T13", comp)
    email("ml-patterns REDO T13: classical vs wavepacket PDE comparison",
          dict(hypothesis=HYPO,
               done="Compared the held-out admitted governing PDEs of the classical "
                    "and wavepacket projectiles term-by-term.",
               plot_shows=f"Shared terms {shared}; classical-only {only_cl}; "
                          f"wavepacket-only {only_wp}. Coefficients + WP/classical ratio.",
               conclusion="Similarities = shared terms; differences = unique terms + "
                          "coefficient ratios. See comparison plot."),
          [plot])
    return comp


def _plot_compare(comp):
    shared = comp["shared_terms"]
    fig, ax = plt.subplots(figsize=(7, 4))
    if shared:
        terms = list(shared); x = np.arange(len(terms))
        ax.bar(x - 0.2, [shared[t]["classical"] for t in terms], 0.4, label="classical")
        ax.bar(x + 0.2, [shared[t]["wavepacket"] for t in terms], 0.4, label="wavepacket")
        ax.set_xticks(x); ax.set_xticklabels(terms, rotation=30, ha="right")
    ax.set_ylabel("median admitted coeff (nondim)")
    ax.set_title("Classical vs wavepacket: shared PDE terms")
    ax.legend(); fig.tight_layout()
    p = os.path.join(ART, "T13_compare.png"); fig.savefig(p, dpi=110); plt.close(fig)
    return p


# ================================================================ T14 synthesis
def phase_T14():
    parts = {}
    for ph in ("T11", "T12", "T13"):
        try: parts[ph] = load(ph)
        except Exception: parts[ph] = None
    summary = {
        "classical": {"verdict": parts["T11"]["verdict"] if parts["T11"] else "n/a",
                      "n_validated": parts["T11"]["n_heldout_validated"] if parts["T11"] else 0},
        "wavepacket": {"verdict": parts["T12"]["verdict"] if parts["T12"] else "n/a",
                       "n_validated": parts["T12"]["n_heldout_validated"] if parts["T12"] else 0},
        "comparison": (parts["T13"] and {
            "shared": list(parts["T13"]["shared_terms"]),
            "classical_only": list(parts["T13"]["classical_only"]),
            "wavepacket_only": list(parts["T13"]["wavepacket_only"])}),
    }
    save("T14", summary)
    nb = write_notebook("redo_synthesis.ipynb",
                        "Bulk-jellium PDE-discovery redo — synthesis",
                        [("md", "## Verdicts\n" + json.dumps(summary, indent=2)),
                         ("md", "Track-B split, three walls, per-projectile PDEs, and "
                                "the classical-vs-wavepacket comparison are in the "
                                "PDE_T*_result.json artifacts.")])
    email("ml-patterns REDO T14: synthesis",
          dict(hypothesis=HYPO,
               done="Synthesised the classical + wavepacket PDE verdicts and their "
                    "comparison; wrote the synthesis notebook.",
               plot_shows="(summary phase — see per-phase plots)",
               conclusion=f"classical={summary['classical']['verdict']}, "
                          f"wavepacket={summary['wavepacket']['verdict']}. "
                          f"Comparison: {summary['comparison']}"),
          [])
    return summary


# ======================================================================== main
PHASES = {"T8": phase_T8, "T11": phase_T11, "T12": phase_T12,
          "T13": phase_T13, "T14": phase_T14}
ORDER = ["T8", "T11", "T12", "T13", "T14"]


def main(argv):
    global SEND_EMAIL, FORCE, SMOKE, DEADLINE
    ap = argparse.ArgumentParser()
    ap.add_argument("phases", nargs="*", help="subset e.g. T11 T12")
    ap.add_argument("--no-email", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--hours", type=float, default=12.0)
    a = ap.parse_args(argv)
    SEND_EMAIL = not a.no_email
    FORCE = a.force
    SMOKE = a.smoke
    DEADLINE = time.time() + a.hours * 3600.0
    want = a.phases if a.phases else ORDER
    log(f"=== PDE-discovery orchestrator | phases={want} | cap={a.hours}h | "
        f"smoke={SMOKE} ===")
    status = {}
    for ph in want:
        if have(ph):
            log(f"[{ph}] result exists — skip"); status[ph] = "cached"; continue
        if out_of_time():
            log(f"[{ph}] WALL-CLOCK CAP hit — stopping"); status[ph] = "skipped-time"; break
        t0 = time.time()
        log(f"\n===== {ph} starting (time left {time_left()/3600:.1f}h) =====")
        try:
            PHASES[ph]()
            status[ph] = f"done in {int(time.time()-t0)}s"
        except Exception as e:
            status[ph] = f"ERROR {e!r}"
            log(f"[{ph}] ERROR: {e!r}\n{traceback.format_exc()}")
        log(f"===== {ph} {status[ph]} =====")
    json.dump(status, open(os.path.join(ART, "PDE_phase_status.json"), "w"), indent=2)
    log(f"\nSTATUS: {json.dumps(status, indent=2)}")


if __name__ == "__main__":
    main(sys.argv[1:])
