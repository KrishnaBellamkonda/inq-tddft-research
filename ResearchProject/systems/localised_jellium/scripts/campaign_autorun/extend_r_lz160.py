#!/usr/bin/env python3
"""Extend the localised-jellium E_total(0)-E_GS vs r energetics plot — larger box.

Adds data points to the EXISTING H0-style plot (same density/system: a=12.5, N=82,
n0=1.31e-3, sigma_WP=0.5, spacing 0.5, w=0) at a bigger box L_z=160 so the projectile
can sit at larger r without hitting the boundary. Slab kept CENTERED (user decision).

New r set (no duplicates with the existing {4,12,20,28,36,40}):
  lower fill : 8, 16, 24, 32
  higher     : 44, 48, 52, 56, 60      (extends 40 -> 60; z=-(12.5+r), margin>=7.5 Bohr=15σ)
  overlap    : 40  (continuity check: L_z=160 vs the existing L_z=120 r=40 point)
periodicity 2 AND 3 for every point. Per r: classical ghost (E_GS + E_total(0)) and
WP (>=1 timestep -> E_total). GS built ONCE per periodicity at L_z=160 and reused.

Root-cause hardening (cf. the 2026-07-01 P0b SIGHUP death): every child is launched
with start_new_session=True (detached); the GS build (the one slow step) has a
liveness guard; the launcher itself is meant to run under `setsid nohup`.

Launch:
  cd .../campaign_autorun
  setsid nohup /local/data/public/skcb2/tddft/venv/bin/python3 extend_r_lz160.py \
      > extend_r_lz160.log 2>&1 &
"""
from __future__ import annotations
import os, sys, csv, time, subprocess, traceback
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
LJ   = ROOT / "ResearchProject/systems/localised_jellium"
CA   = LJ / "scripts/campaign_autorun"
GSBIN, WPBIN, CLBIN = str(CA/"gs/run"), str(CA/"wp/run"), str(CA/"classical/run")
RUNS = CA / "runs/extend_r160"; RUNS.mkdir(parents=True, exist_ok=True)
TO   = "chiddukanna@gmail.com"
HA   = 27.211386

LZ, HALF, N = 160, 12.5, 82
R_LOW, R_HIGH, R_OVERLAP = (8, 16, 24, 32), (44, 48, 52, 56, 60), (40,)
R_ALL = R_LOW + R_HIGH + R_OVERLAP
PERS = (3, 2)
GS_STALL_MIN = 40                     # kill+fail the GS build if its log is silent this long

ENV = {**os.environ,
       "INQ_SHARE_PATH":       str(ROOT/"inq/install/share"),
       "PSEUDOPOD_SHARE_PATH": str(ROOT/"inq/install/share/pseudopod"),
       "INQ_SOURCE":           str(ROOT/"inq-study")}

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)

def gs_ckpt(per): return RUNS / f"gs_lz160_p{per}" / "checkpoint"

def done(rundir: Path) -> bool:
    for rs in Path(rundir).glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs.read_text(): return True
        except Exception: pass
    return False

def launch(binary, rundir: Path, overrides: dict, gpu: int) -> subprocess.Popen:
    rundir.mkdir(parents=True, exist_ok=True)
    env = {**ENV, "CUDA_VISIBLE_DEVICES": str(gpu), **{k: str(v) for k, v in overrides.items()}}
    lf = open(rundir / "run.log", "w")
    return subprocess.Popen([binary], cwd=str(rundir), env=env, stdout=lf,
                            stderr=subprocess.STDOUT, start_new_session=True)

def email(subject, body, attachments=None):
    try:
        sys.path.insert(0, str(ROOT / "inq-stack/python"))
        from inqview.email import send_run_email
        return send_run_email(subject="[extend-r-lz160] " + subject, body=body,
                              attachments=attachments or [], to=TO)
    except Exception as e:
        log(f"  EMAIL FAILED: {e}")

# ----------------------------------------------------------------- GS build (gated)
def build_gs():
    ckpts = {per: gs_ckpt(per) for per in PERS}
    if all(done(RUNS / f"gs_lz160_p{per}") for per in PERS):
        log("both L_z=160 GS already complete — skipping build"); return True
    procs = {}
    for i, per in enumerate(PERS):
        rd = RUNS / f"gs_lz160_p{per}"
        if done(rd): log(f"  GS p{per} already done"); continue
        ov = dict(LJ_LX=50, LJ_LY=50, LJ_LZ=LZ, LJ_HALF=HALF, LJ_N=N, LJ_EDGE_W=0,
                  LJ_PERIODICITY=per, LJ_SPACING=0.5, LJ_GS_DIR=str(ckpts[per]),
                  LJ_TAG=f"gs_lz160_p{per}")
        procs[per] = (launch(GSBIN, rd, ov, gpu=i), rd / "run.log", rd)
        log(f"  building GS p{per} (GPU{i}, pid {procs[per][0].pid})")
    last_sz = {p: -1 for p in procs}; last_mv = {p: time.time() for p in procs}
    while procs and not all(done(RUNS / f"gs_lz160_p{per}") for per in PERS):
        time.sleep(60)
        for per, (proc, lf, rd) in list(procs.items()):
            if done(rd): continue
            sz = lf.stat().st_size if lf.exists() else 0
            if sz > last_sz[per]: last_sz[per] = sz; last_mv[per] = time.time()
            elif proc.poll() is not None: last_mv[per] = 0
            if (time.time() - last_mv[per]) / 60.0 > GS_STALL_MIN:
                log(f"  GS p{per} STALLED/died — aborting GS build")
                for pr, _, _ in procs.values():
                    try: pr.terminate()
                    except Exception: pass
                return False
    ok = all(done(RUNS / f"gs_lz160_p{per}") for per in PERS)
    log(f"GS build {'complete' if ok else 'FAILED'} for both periodicities")
    return ok

# ------------------------------------------------------------- projectile r-sweep
def sweep():
    jobs = [(t, r, per) for t in ("wp", "cl") for per in PERS for r in R_ALL]
    # run 2 at a time (one per GPU)
    i = 0
    while i < len(jobs):
        batch = jobs[i:i+2]; running = []
        for g, (t, r, per) in enumerate(batch):
            out = f"{t}_r{r}_p{per}"; rd = RUNS / out
            if done(rd): log(f"  SKIP {out} (done)"); continue
            z = -(HALF + r)
            ov = dict(LJ_OUT="results", LJ_LZ=LZ, LJ_HALF=HALF, LJ_N=N, LJ_EDGE_W=0,
                      LJ_PERIODICITY=per, LJ_SPACING=0.5, LJ_LAUNCH_Z=z,
                      LJ_GS_DIR=str(gs_ckpt(per)))
            if t == "wp": ov.update(LJ_K0=0, LJ_SIGMA=0.5)
            binary = WPBIN if t == "wp" else CLBIN
            running.append((launch(binary, rd, ov, gpu=g), out, rd))
            log(f"  RUN {out} (GPU{g})")
        for proc, out, rd in running:
            proc.wait()
            if not done(rd):                       # one-shot retry
                log(f"  {out} incomplete -> retry")
                g = 0
                z = -(HALF + int(out.split("_r")[1].split("_")[0]))
                # rebuild env from name
                t = out.split("_")[0]; per = int(out.split("_p")[1])
                r = int(out.split("_r")[1].split("_")[0]); z = -(HALF + r)
                ov = dict(LJ_OUT="results", LJ_LZ=LZ, LJ_HALF=HALF, LJ_N=N, LJ_EDGE_W=0,
                          LJ_PERIODICITY=per, LJ_SPACING=0.5, LJ_LAUNCH_Z=z,
                          LJ_GS_DIR=str(gs_ckpt(per)))
                if t == "wp": ov.update(LJ_K0=0, LJ_SIGMA=0.5)
                launch(WPBIN if t == "wp" else CLBIN, rd, ov, gpu=0).wait()
        i += 2

# ------------------------------------------------------------------ plot + report
def gs_energy(per):
    for line in (RUNS / f"gs_lz160_p{per}").glob("**/run_summary.txt"):
        for ln in line.read_text().splitlines():
            if ln.startswith("ground_state_energy_ha"):
                return float(ln.split("=")[1])
    return None

def etotal0(rd: Path):
    # the binary prepends "results/" to LJ_OUT, so the CSV nests at
    # <rd>/results/results/raw/observables/observables.csv — glob to be layout-robust.
    matches = sorted(rd.glob("**/observables.csv"))
    if not matches: return None
    with open(matches[0]) as f: rows = list(csv.reader(f))
    if len(rows) < 2: return None
    return float(rows[1][rows[0].index("energy_total")])

def build_plot():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    try:
        from inqview.visualisation.style import apply_theme; apply_theme()
    except Exception: pass
    egs = {per: gs_energy(per) for per in PERS}
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    styles = {("cl", 3): "s-", ("cl", 2): "s--", ("wp", 3): "o-", ("wp", 2): "o--"}
    data = {}
    for t in ("wp", "cl"):
        for per in PERS:
            pts = []
            for r in sorted(set(R_ALL)):
                e0 = etotal0(RUNS / f"{t}_r{r}_p{per}")
                if e0 is not None and egs[per] is not None:
                    pts.append((r, (e0 - egs[per]) * HA))
            if pts:
                xs, ys = zip(*pts); data[(t, per)] = pts
                lab = f"{'WP' if t=='wp' else 'classical'} p{per}" + (" (raw; net-Q corr pending)" if t=="wp" and per==2 else "")
                ax.plot(xs, ys, styles[(t, per)], label=lab, ms=5)
    ax.set_xlabel("projectile-slab distance r (Bohr)")
    ax.set_ylabel(r"$E_\mathrm{total}(0)-E_\mathrm{GS}$ (eV)")
    ax.set_title("Localised jellium — energetics vs r, extended to r=60 (L_z=160)")
    ax.legend(frameon=False, fontsize=8)
    png = RUNS / "extend_r160_excess_vs_r.png"
    fig.savefig(png, dpi=140, bbox_inches="tight"); plt.close(fig)
    return png, data

# ------------------------------------------------------------------------- main
def main():
    log("EXTEND-r L_z=160 — build GS (gated) -> r-sweep -> plot")
    if not build_gs():
        email("GS build FAILED — HALTING",
              "The L_z=160 neutral GS (periodicity 3 and/or 2) did not converge/complete; "
              "the r-sweep was NOT started. Check runs/extend_r160/gs_lz160_p*/run.log.")
        log("GS gate failed — halting"); return
    sweep()
    n_done = sum(done(RUNS / f"{t}_r{r}_p{per}")
                 for t in ("wp", "cl") for per in PERS for r in R_ALL)
    try:
        png, data = build_plot()
    except Exception:
        png = None; log("plot failed:\n" + traceback.format_exc())
    body = (
        "HYPOTHESIS/PURPOSE: extend the localised-jellium E_total(0)-E_GS vs r energetics "
        "plot to larger projectile-slab distance r, same system/density as before, using a "
        "bigger centered box (L_z=160) so higher r fits without charge leaking to the boundary.\n\n"
        f"WHAT WAS DONE: built the L_z=160 neutral GS for periodicity 3 and 2 (once each), then "
        f"swept r={list(R_ALL)} Bohr x periodicity{{3,2}} x {{classical ghost, WP}}. "
        f"{n_done}/{len(R_ALL)*4} projectile runs completed. Classical gives E_GS + E_total(0); "
        f"WP is propagated >=1 step for a valid E_total.\n\n"
        "WHAT THE PLOT SHOWS: E_total(0)-E_GS (eV) vs r for WP and classical at periodicity 3 and "
        "2, extending the existing curve to r=60. The r=40 point at L_z=160 is a continuity check "
        "against the existing L_z=120 r=40 point.\n\n"
        "CONCLUSION (PROVISIONAL): new r points staged for the existing plot; verify the r=40 "
        "overlap agrees across boxes before merging. NOTE: periodicity-2 WP excess is RAW — the "
        "open-z net-charge G=0 correction is still pending (trust periodicity-3 for absolute values)."
    )
    email("extend-r L_z=160 COMPLETE — energetics vs r", body,
          attachments=[str(png)] if png else None)
    # auto-build the study notebook from all completed runs (notebook-making convention)
    try:
        builder = LJ / "hypotheses/extend_r160/build_extend_r160_report.py"
        subprocess.run([sys.executable, str(builder)],
                       env={**ENV, "PYTHONPATH": str(ROOT / "inq-stack/python")},
                       cwd=str(builder.parent), check=False, timeout=1200)
        log("study notebook rebuilt")
    except Exception:
        log("notebook rebuild failed:\n" + traceback.format_exc())
    log(f"DONE — {n_done} projectile runs; plot={'ok' if png else 'failed'}")

if __name__ == "__main__":
    try: main()
    except Exception:
        tb = traceback.format_exc(); log("FATAL:\n" + tb); email("extend-r FATAL", tb)
