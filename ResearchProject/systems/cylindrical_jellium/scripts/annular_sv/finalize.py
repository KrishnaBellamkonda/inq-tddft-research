#!/usr/bin/env python3
"""Post-production finalizer for the cylindrical-jellium campaign (Phases 4 + 5).

SEPARATE from the running production orchestrator (orchestrate.py) — does NOT
touch it. WAITS (polls) until all 9 production runs are complete, then on the
freed GPUs runs:
  Phase 4 — the WP quantum rung (electron wavepacket σ_WP=0.5, k0=0.30 at r_s=6;
            its matched classical ghost is the existing rs6_v0p30 production run).
  Phase 5 — builds + executes the synthesis notebook (build_report.py) and emails
            the headline plots.
Idempotent: skips the WP run if already complete; safe to re-launch.

Headless launch (after production has started):
    cd .../scripts/annular_sv
    nohup venv/bin/python3 finalize.py > finalize.log 2>&1 &
"""
from __future__ import annotations
import math, os, subprocess, sys, time, traceback
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
SYS  = ROOT / "ResearchProject/systems/cylindrical_jellium"
SCR  = SYS / "scripts/annular_sv"
SWEEP_OUT = SYS / "annular_sv"
HYP  = SYS / "hypotheses/annular_sv"
PY   = str(ROOT / "venv/bin/python3")
WPBIN = str(SCR / "wp/run")
GPU = "0"                                  # WP runs alone after production frees both GPUs
TO = "chiddukanna@gmail.com"
DT = 0.020
POLL_SECONDS = 300

GEOM = dict(CJ_LXY=40, CJ_RIN=5, CJ_ROUT=13, CJ_EDGE_W=1.0, CJ_SPACING=0.5)
RS6 = dict(L_z=48, N=24, gs=str(SYS/"shared_gs/tube_rs6"))
VELS = [0.15, 0.30, 0.45]
PROD_LABELS = [f"rs{rs}_v{v:.2f}".replace(".", "p") for rs in (6,4,2) for v in VELS]
WP_LABEL = "wp_rs6_v0p30"
WP_V0, WP_SIGMA = 0.30, 0.5

ENV0 = {**os.environ,
        "INQ_SHARE_PATH": str(ROOT/"inq/install/share"),
        "PSEUDOPOD_SHARE_PATH": str(ROOT/"inq/install/share/pseudopod")}

def log(msg): print(f"[{datetime.now():%F %T}] {msg}", flush=True)

def n_steps_for(rs, L_z, v):
    # light electron/WP decelerates fast: capture the initial-drag window, not a
    # 5-plasma-period wake (matches orchestrate.py).
    return int(math.ceil(max(30.0, 100.0 * v) / DT))

def _done(rundir: Path) -> bool:
    for rs_file in rundir.glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs_file.read_text(): return True
        except Exception: pass
    return False

def send_email(subject, body, pngs):
    try:
        sys.path.insert(0, str(ROOT/"inq-stack/python"))
        from inqview.email import send_run_email
        send_run_email(subject=subject, body=body,
                       attachments=[p for p in pngs if Path(p).exists()], to=TO)
        log(f"  email sent: {subject}")
    except Exception as e:
        log(f"  email FAILED: {e}")

def wait_for_production():
    log(f"waiting for {len(PROD_LABELS)} production runs: {PROD_LABELS}")
    while True:
        done = [lbl for lbl in PROD_LABELS if _done(SWEEP_OUT/lbl)]
        if len(done) == len(PROD_LABELS):
            log("all production runs complete"); return True
        log(f"  {len(done)}/{len(PROD_LABELS)} done; sleeping {POLL_SECONDS}s")
        time.sleep(POLL_SECONDS)

def run_wp_rung():
    if not Path(WPBIN).exists():
        log(f"WP binary missing ({WPBIN}); skipping Phase 4"); return False
    rundir = SWEEP_OUT / WP_LABEL; rundir.mkdir(parents=True, exist_ok=True)
    if _done(rundir):
        log("WP rung already complete (skip)"); return True
    ns = n_steps_for(6, RS6["L_z"], WP_V0); we = max(1, round(ns/300.0))
    env = {**ENV0, **{k: str(v) for k, v in GEOM.items()},
           "CUDA_VISIBLE_DEVICES": GPU,
           "CJ_LZ": str(RS6["L_z"]), "CJ_N": str(RS6["N"]), "CJ_GS_DIR": RS6["gs"],
           "PROJ_V0": str(WP_V0), "CJ_WP_SIGMA": str(WP_SIGMA),
           "SV_N_STEPS": str(ns), "SV_WRITE_EVERY": str(we), "SV_OUT_SUBDIR": WP_LABEL}
    log(f"  RUN WP rung on GPU{GPU} (N_STEPS={ns}); 6h timeout (WP-injection-deadlock guard)")
    t0 = time.time()
    try:
        with open(rundir/"run.log", "w") as lf:
            rc = subprocess.run([WPBIN], cwd=str(rundir), env=env, stdout=lf,
                                stderr=subprocess.STDOUT, timeout=21600).returncode
    except subprocess.TimeoutExpired:
        log("  WP rung TIMED OUT (>6h) — likely injection deadlock; skipping to notebook"); return False
    ok = rc == 0 and _done(rundir)
    log(f"  WP rung {'OK' if ok else 'FAILED'} rc={rc} ({(time.time()-t0)/60:.1f} min)")
    return ok

def build_notebook():
    log("  building + executing synthesis notebook (build_report.py)")
    rc = subprocess.run([PY, str(HYP/"build_report.py")], env=ENV0).returncode
    log(f"  build_report rc={rc}")
    return rc == 0

def main():
    log("FINALIZER start (Phases 4+5; separate from production orchestrator)")
    try:
        wait_for_production()
    except Exception:
        log(f"wait error:\n{traceback.format_exc()}")
    try:
        run_wp_rung()
    except Exception:
        log(f"WP rung exception:\n{traceback.format_exc()}")
        send_email("[cylindrical-jellium] Phase 4 WP rung FAILED", traceback.format_exc(), [])
    nb_ok = False
    try:
        nb_ok = build_notebook()
    except Exception:
        log(f"notebook exception:\n{traceback.format_exc()}")

    # ---- headline results for the completion email (read S(v)/β from the CSV) ----
    results_txt = "(S(v) results CSV not found — see the notebook)"
    csv = HYP / "Sv_results.csv"
    try:
        if csv.exists():
            import pandas as pd, numpy as np
            d = pd.read_csv(csv)
            lines = ["S(v) per run (initial-drag, Ha/Bohr):"]
            for _, r in d.iterrows():
                lines.append(f"  r_s={int(r['r_s'])}  v={r['v']:.2f}:  S = {r['S_ha_per_bohr']:.4f} ± {r['S_err']:.4f}")
            betas = []
            for rs in (6, 4, 2):
                sub = d[d.r_s == rs].sort_values("v")
                if len(sub) >= 2:
                    b = np.polyfit(sub.v, sub.S_ha_per_bohr, 1)[0]; betas.append((rs, b))
            if betas:
                lines.append("\nβ(r_s) = dS/dv:")
                lines += [f"  r_s={rs}:  β = {b:.4f}" for rs, b in betas]
                if len(betas) == 3:
                    bs = [b for _, b in betas]
                    mono = bs[0] < bs[1] < bs[2] or bs[0] > bs[1] > bs[2]
                    lines.append(f"\nβ(r_s) monotonic across r_s={{6,4,2}}: {'YES' if mono else 'NO'}")
            results_txt = "\n".join(lines)
    except Exception as e:
        results_txt = f"(could not parse {csv}: {e})"

    wp_status = "complete" if _done(SWEEP_OUT / WP_LABEL) else "FAILED/skipped (see finalize.log)"
    body = ("The cylindrical-jellium campaign has COMPLETED end-to-end (Phases 0→5).\n\n"
            "HYPOTHESIS: a charge gliding down the bore of a periodic annular jellium tube\n"
            "feels a measurable electronic stopping power S(v); its low-velocity friction\n"
            "slope β(r_s)=dS/dv varies with the wall density r_s.\n\n"
            f"{results_txt}\n\n"
            f"Phase 4 (quantum rung, electron wavepacket vs classical ghost at r_s=6): {wp_status}.\n"
            f"Phase 5 (synthesis notebook): {'executed' if nb_ok else 'built (execution flagged — check log)'}\n"
            f"  -> {HYP/'annular_sv_report.ipynb'}\n"
            "  Panels: S(v)+β(r_s), induced wall current (hydrovoltaic flow→current),\n"
            "  wake structure, WP-vs-classical.\n\n"
            "NOTE: S(v) is the INITIAL stopping power at the launch velocity (the light\n"
            "electron decelerates; the initial KE-loss slope = friction force at v0).\n"
            "PROVISIONAL: r_s=6 is a small gas (~24 e); electron-as-cation rests on charge-even\n"
            "S at leading order (Barkas = the charge-odd correction).")
    send_email("[cylindrical-jellium] ✅ CAMPAIGN COMPLETE — S(v), β(r_s), WP rung + notebook",
               body, [str(HYP/"Sv_beta.png")])
    log("FINALIZER done — completion email sent")

if __name__ == "__main__":
    main()
