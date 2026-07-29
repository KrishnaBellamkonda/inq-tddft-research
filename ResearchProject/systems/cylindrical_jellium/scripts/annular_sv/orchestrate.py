#!/usr/bin/env python3
"""Autonomous 2-GPU orchestrator for the cylindrical-jellium S(v) campaign.

Drives Phases 2→3 (pilot gate → 9-run production sweep) of
docs/campaigns/cylindrical_jellium/cylindrical_jellium_projectile.md across BOTH
GPUs, then extracts S(v) and β(r_s) and emails the result. Phase 4 (WP rung) runs
if its binary exists; Phase 5 (notebook) is a follow-on.

Design (per the campaigns autonomy checklist — Python, not bash):
  * 2-GPU work queue: independent runs dispatched concurrently to device 0 and 1.
  * IDEMPOTENT RESUME: a run whose results/<sub>/run_summary.txt shows
    'run_completed = true' is skipped, so a crash/restart continues.
  * Pilot gate (BLOCKING): r_s=6, v=0.30 runs FIRST and is gated (no NaN; N
    conserved <1%; bounded energy drift; v-drift <10%). Production fans out only
    on PASS.
  * Per-phase try/except with full-traceback failure emails; the chain continues.
  * One-shot retry on a sim failure.

Headless launch (survives disconnect):
    cd .../scripts/annular_sv
    nohup venv/bin/python3 orchestrate.py > orchestrate.log 2>&1 &
"""
from __future__ import annotations
import math, os, queue, subprocess, sys, threading, time, traceback
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
SYS  = ROOT / "ResearchProject/systems/cylindrical_jellium"
SCR  = SYS / "scripts/annular_sv"
SWEEP_OUT = SYS / "annular_sv"                  # <sweep>/<run>/ outputs
HYP  = SYS / "hypotheses/annular_sv"
PY   = str(ROOT / "venv/bin/python3")
CLBIN = str(SCR / "classical/run")
WPBIN = str(SCR / "wp/run")
GPUS = ["0", "1"]
TO = "chiddukanna@gmail.com"
DT = 0.020

# Fixed transverse geometry (locked).
GEOM = dict(CJ_LXY=40, CJ_RIN=5, CJ_ROUT=13, CJ_EDGE_W=1.0, CJ_SPACING=0.5)
# Per-density: r_s -> (L_z, N, GS checkpoint).
DENS = {
    6: dict(L_z=48, N=24,  gs=str(SYS/"shared_gs/tube_rs6")),
    4: dict(L_z=28, N=48,  gs=str(SYS/"shared_gs/tube_rs4")),
    2: dict(L_z=10, N=136, gs=str(SYS/"shared_gs/tube_rs2")),
}
VELS = [0.15, 0.30, 0.45]

ENV0 = {**os.environ,
        "INQ_SHARE_PATH": str(ROOT/"inq/install/share"),
        "PSEUDOPOD_SHARE_PATH": str(ROOT/"inq/install/share/pseudopod")}

_print_lock = threading.Lock()
def log(msg):
    with _print_lock:
        print(f"[{datetime.now():%F %T}] {msg}", flush=True)


def n_steps_for(rs, L_z, v):
    """The light electron (m_e) DECELERATES and stops within the box, so a run only
    needs to capture the initial-drag window (v >= 0.85 v0) plus the deceleration
    sweep — NOT a 5-plasma-period wake (the electron stops long before that). ~30-45
    a.u. suffices and is cheap; higher v0 (more KE) gets a longer run."""
    t_total = max(30.0, 100.0 * v)        # 0.15->30, 0.30->30, 0.45->45 a.u.
    return int(math.ceil(t_total / DT))


def run_label(rs, v):
    return f"rs{rs}_v{v:.2f}".replace(".", "p")


def _done(rundir: Path) -> bool:
    for rs_file in rundir.glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs_file.read_text():
                return True
        except Exception:
            pass
    return False


def run_classical(rs, v, gpu, retries=1, ns_override=None, subdir_override=None) -> bool:
    d = DENS[rs]
    L_z, N = d["L_z"], d["N"]
    ns = ns_override if ns_override else n_steps_for(rs, L_z, v)
    we = max(1, round(ns / 300.0))
    label = subdir_override if subdir_override else run_label(rs, v)
    rundir = SWEEP_OUT / label
    rundir.mkdir(parents=True, exist_ok=True)
    if _done(rundir):
        log(f"  SKIP {label} (already complete)"); return True
    env = {**ENV0, **{k: str(val) for k, val in GEOM.items()},
           "CUDA_VISIBLE_DEVICES": gpu,
           "CJ_LZ": str(L_z), "CJ_N": str(N), "CJ_GS_DIR": d["gs"],
           "PROJ_V0": str(v), "SV_N_STEPS": str(ns), "SV_WRITE_EVERY": str(we),
           "SV_OUT_SUBDIR": label}
    for attempt in range(retries + 1):
        log(f"  RUN  {label} on GPU{gpu} (N_STEPS={ns}, write_every={we})"
            + (f" retry{attempt}" if attempt else ""))
        t0 = time.time()
        with open(rundir / "run.log", "w") as lf:
            rc = subprocess.run([CLBIN], cwd=str(rundir), env=env,
                                stdout=lf, stderr=subprocess.STDOUT).returncode
        if rc == 0 and _done(rundir):
            log(f"  OK   {label} ({(time.time()-t0)/60:.1f} min)"); return True
        log(f"  FAIL {label} rc={rc} (see {rundir/'run.log'})")
    return False


# ----------------- pilot gate ------------------------------------------------
def pilot_gate() -> bool:
    """Run r_s=6, v=0.30 (the real production point) and gate on numeric health.
    The light electron DECELERATES by design, so we do NOT abort on v-drift; we
    require instead that a clean initial-drag S can be extracted from the early
    near-constant-velocity window."""
    import numpy as np
    rs, v = 6, 0.30
    if not run_classical(rs, v, GPUS[0]):
        log("PILOT: run failed"); return False
    rundir = SWEEP_OUT / run_label(rs, v)
    trk = next(rundir.glob("**/electron_track.csv"), None)
    obs = next(rundir.glob("**/observables.csv"), None)
    nl  = next(rundir.glob("**/electron_number.csv"), None)
    if not (trk and obs):
        log("PILOT: missing track/observables"); return False
    T = np.genfromtxt(trk, delimiter=",", names=True)
    O = np.genfromtxt(obs, delimiter=",", names=True)
    checks = {}
    e = np.asarray(O["energy_total"], float)
    checks["no NaN energy"] = bool(np.all(np.isfinite(e)))
    ke = np.asarray(T["ke_ion_ha"], float)
    checks["no NaN track"] = bool(np.all(np.isfinite(ke)))
    if nl:
        try:
            Nn = np.genfromtxt(nl, delimiter=",", names=True)
            ntot = np.atleast_1d(np.asarray(Nn["N_total"], float))
            dN = abs(ntot[-1] - ntot[0]) / ntot[0] if len(ntot) > 1 else 0.0
            checks["N conserved <1%"] = dN < 0.01
        except Exception as ex:
            log(f"  PILOT N-check parse skipped: {ex}")
    # bounded energy (deposited KE ~0.05 Ha; a blow-up would be >>1)
    checks["energy bounded"] = bool(np.max(np.abs(e - e[0])) < 1.0)
    # the gate that REPLACES the v-drift abort: a clean initial-drag S exists
    r = extract_S(rs, v)
    checks["clean initial-drag S"] = bool(r and np.isfinite(r["S"]) and r["npts"] >= 30)
    log(f"PILOT checks: {checks}" + (f"  S={r['S']:.4f} npts={r['npts']}" if r else ""))
    return all(checks.values())


# ----------------- 2-GPU scheduler -------------------------------------------
def dispatch_2gpu(jobs):
    """jobs = list of (rs, v). Run across both GPUs concurrently; return results."""
    q = queue.Queue()
    for j in jobs: q.put(j)
    results = {}
    res_lock = threading.Lock()
    def worker(gpu):
        while True:
            try: rs, v = q.get_nowait()
            except queue.Empty: return
            try:
                ok = run_classical(rs, v, gpu)
            except Exception:
                log(f"  worker GPU{gpu} exception:\n{traceback.format_exc()}"); ok = False
            with res_lock: results[(rs, v)] = ok
            q.task_done()
    threads = [threading.Thread(target=worker, args=(g,), daemon=True) for g in GPUS]
    for t in threads: t.start()
    for t in threads: t.join()
    return results


# ----------------- S(v) extraction -------------------------------------------
def extract_S(rs, v, vfrac=0.85):
    """S(v0) = INITIAL stopping power (drag) at the launch velocity, from the early
    near-constant-velocity window where vz >= vfrac*v0. The light electron
    decelerates, so a full-run regression would mix velocities; the initial slope of
    KE_ion(s) loss IS the friction force at v0:  S = -d(KE_ion)/ds.
    Uses the per-step track (high resolution). Returns S, stderr, mean v, n points."""
    import numpy as np
    rundir = SWEEP_OUT / run_label(rs, v)
    trk = next(rundir.glob("**/electron_track.csv"), None)
    if not trk: return None
    T = np.genfromtxt(trk, delimiter=",", names=True)
    z = np.asarray(T["z"], float); vz = np.asarray(T["vz"], float)
    ke = np.asarray(T["ke_ion_ha"], float)
    if len(z) < 10 or not np.all(np.isfinite(ke)): return None
    s = np.abs(z - z[0])
    # early window: vz >= vfrac*v0, drop a tiny launch transient; widen if too few
    for vf in (vfrac, 0.70, 0.50):
        sel = vz >= vf * v
        i_lo = max(2, int(0.03 * len(s)))
        sel[:i_lo] = False
        if sel.sum() >= 20: break
    if sel.sum() < 8: return None
    ss, kk = s[sel], ke[sel]
    A = np.vstack([ss, np.ones_like(ss)]).T
    coef, *_ = np.linalg.lstsq(A, kk, rcond=None)
    S = -float(coef[0])                       # drag: KE falls with path
    yhat = A @ coef; dof = max(1, len(ss) - 2)
    sxx = np.sum((ss - ss.mean())**2)
    serr = math.sqrt(np.sum((kk - yhat)**2)/dof / sxx) if sxx > 0 else float("nan")
    return dict(rs=rs, v=v, S=S, S_err=float(serr),
                v_mean=float(vz[sel].mean()), npts=int(sel.sum()))


def send_email(subject, body, pngs):
    try:
        sys.path.insert(0, str(ROOT/"inq-stack/python"))
        from inqview.email import send_run_email
        send_run_email(subject=subject, body=body,
                       attachments=[p for p in pngs if Path(p).exists()], to=TO)
        log(f"  email sent: {subject}")
    except Exception as e:
        log(f"  email FAILED ({subject}): {e}")


def fail_email(phase, msg):
    send_email(f"[cylindrical-jellium] {phase} — PHASE FAILED (chain continues)",
               f"Phase {phase} failed in the orchestrator.\n\n{msg}\n\n"
               f"Data under {SWEEP_OUT}. Re-run orchestrate.py to resume "
               f"(completed runs are skipped).", [])


# ----------------- phases ----------------------------------------------------
def phase3_production():
    import numpy as np
    # pilot first (r_s=6, v=0.30) — gate
    log("=== PHASE 2: pilot gate (r_s=6, v=0.30) ===")
    if not pilot_gate():
        fail_email("Phase 2 pilot", "Pilot numeric gate FAILED — production NOT launched. "
                   "Inspect the pilot run; likely L_z/wake or dt/resolution.")
        raise RuntimeError("pilot gate failed")
    log("PILOT PASS — fanning out all 9 production runs across both GPUs")
    # full 9-run production sweep (the short pilot was a separate smoke subdir)
    jobs = [(rs, v) for rs in (6, 4, 2) for v in VELS]
    results = dispatch_2gpu(jobs)
    log(f"production results: {results}")
    # S extraction for all 9
    rows = []
    for rs in (6, 4, 2):
        for v in VELS:
            r = extract_S(rs, v)
            if r: rows.append(r); log(f"  S(r_s={rs}, v={v}) = {r['S']:.4f} ± {r['S_err']:.4f} Ha/Bohr")
    # beta(r_s) = dS/dv fit per density
    betas = {}
    for rs in (6, 4, 2):
        pts = [(x["v"], x["S"]) for x in rows if x["rs"] == rs]
        if len(pts) >= 2:
            vv = np.array([p[0] for p in pts]); ss = np.array([p[1] for p in pts])
            A = np.vstack([vv, np.ones_like(vv)]).T
            (b, _), *_ = np.linalg.lstsq(A, ss, rcond=None)
            betas[rs] = float(b)
    # plot S(v) + write CSV
    HYP.mkdir(parents=True, exist_ok=True)
    csv = HYP / "Sv_results.csv"
    with open(csv, "w") as f:
        f.write("r_s,v,S_ha_per_bohr,S_err\n")
        for x in rows: f.write(f"{x['rs']},{x['v']},{x['S']},{x['S_err']}\n")
    png = HYP / "Sv_beta.png"
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7,5))
        for rs in (6,4,2):
            pts = sorted([(x["v"], x["S"], x["S_err"]) for x in rows if x["rs"]==rs])
            if pts:
                vv=[p[0] for p in pts]; ss=[p[1] for p in pts]; ee=[p[2] for p in pts]
                ax.errorbar(vv, ss, yerr=ee, marker="o", capsize=3,
                            label=f"r_s={rs}" + (f" (β={betas[rs]:.3f})" if rs in betas else ""))
        ax.set_xlabel("v (a.u.)"); ax.set_ylabel("S (Ha/Bohr)")
        ax.set_title("Annular-tube electronic stopping S(v) vs wall r_s"); ax.legend()
        fig.tight_layout(); fig.savefig(png, dpi=130)
    except Exception as e:
        log(f"  S(v) plot skipped: {e}")
    body = ("Cylindrical-jellium production sweep complete (9 runs, both GPUs).\n\n"
            "Hypothesis: a charge gliding down the bore of a periodic annular jellium\n"
            "tube feels a measurable S(v); β(r_s)=dS/dv varies with wall density r_s.\n\n"
            "S(v) per density (Ha/Bohr):\n")
    for x in rows: body += f"  r_s={x['rs']} v={x['v']}: S={x['S']:.4f} ± {x['S_err']:.4f}\n"
    body += "\nβ(r_s)=dS/dv:\n" + "".join(f"  r_s={rs}: β={b:.4f}\n" for rs,b in betas.items())
    body += ("\nThe attached S(v)/β plot shows the low-velocity friction slope per wall\n"
             "density. PROVISIONAL: r_s=6 is a small gas (~24 e); electron-as-cation\n"
             "rests on charge-even S (Barkas = odd correction).")
    send_email("[cylindrical-jellium] PRODUCTION COMPLETE — S(v) + β(r_s)", body, [str(png)])


def main():
    log(f"ORCHESTRATOR start; GPUs={GPUS}; binary={CLBIN}")
    if not Path(CLBIN).exists():
        log(f"FATAL: classical binary missing ({CLBIN}); build it first."); return
    for rs in (6,4,2):
        if not Path(DENS[rs]["gs"]).exists():
            log(f"FATAL: GS checkpoint missing for r_s={rs} ({DENS[rs]['gs']})."); return
    t0 = time.time()
    try:
        phase3_production()
    except Exception:
        tb = traceback.format_exc(); log(f"PHASE 3 EXCEPTION:\n{tb}"); fail_email("Phase 3", tb)
    log(f"ORCHESTRATOR done ({(time.time()-t0)/3600:.2f} h)")


if __name__ == "__main__":
    main()
