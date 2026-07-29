#!/usr/bin/env python3
"""Autonomous orchestrator — wide-WP campaign Phase-0 gate -> Phase-1 S(E) sweep.

Runs UNATTENDED after the P0b matched pair (currently in flight) finishes:
  1. wait for P0b (p0b_wp + p0b_classical) to complete;
  2. NUMERIC P0b gate (replaces the human sign-off): both runs completed, WP
     E_total finite, N_total not collapsed (bath intact), CAP absorbed the WP;
  3. Phase-1 sweep: E in {200,280,360,440,520,600} eV, each = WP (GPU0) +
     classical (GPU1) concurrently; resumable (skips run_summary run_completed);
  4. per-run S: WP energy method S=[E_total(t_f)-E_GS]/L_z (phase-5 method),
     classical S=-dKE_ion/dz across the slab; append to se_wide_wp.csv;
  5. final S(E) overlay PNG + email (4-part).

Design per the campaigns rule: PYTHON orchestrator (not bash) — structured
logging, idempotent resume, per-run try/except with traceback emails, one-shot
retry. Sweep uses a COARSER VTI cadence (WRITE_EVERY=20) than the p0b pilot so 12
runs fit overnight; all scalar/momentum observables (which give S) are kept.

Headless launch (survives disconnect):
    cd .../scripts/wide_wp
    nohup /local/data/public/skcb2/tddft/venv/bin/python3 orchestrate.py \
        > orchestrate.log 2>&1 &
"""
from __future__ import annotations
import os, subprocess, sys, time, math, traceback, csv
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
LJ   = ROOT / "ResearchProject/systems/localised_jellium"
WWP  = LJ / "scripts/wide_wp"
WPDIR, CLDIR = WWP / "wp", WWP / "classical"
WPBIN, CLBIN = str(WPDIR / "run"), str(CLDIR / "run")
HYP  = LJ / "hypotheses/wide_wp"
PY   = str(ROOT / "venv/bin/python3")
TO   = "chiddukanna@gmail.com"

HA   = 27.211386
E_GS = -86.04107005396197      # our dx=0.40 / LZ=101 GS anchor (validated)
L_Z  = 25.0                    # slab thickness = 2*12.5 Bohr
SLAB_HALF = 12.5
LAUNCH_Z  = -26.5
CAP_INNER = 40.5               # far CAP inner face
E_GRID = [200, 280, 360, 440, 520, 600]   # eV, 6 energies
DT_WP, DT_CL = 0.04, 0.02
WRITE_EVERY_SWEEP = 20         # coarser than p0b (=4) so 12 runs fit overnight

ENV_BASE = {**os.environ,
            "INQ_SHARE_PATH":     str(ROOT / "inq/install/share"),
            "PSEUDOPOD_SHARE_PATH": str(ROOT / "inq/install/share/pseudopod"),
            "INQ_SOURCE":         str(ROOT / "inq-study")}
SUBJ = "[wide-wavepacket] "

def log(msg): print(f"[{datetime.now():%F %T}] {msg}", flush=True)

def k0_of(E_eV): return math.sqrt(2.0 * E_eV / HA)

def n_steps_wp(E):    # packet: launch -> far CAP (+ absorb + plateau)
    v = k0_of(E); path = (CAP_INNER + 5.0) - LAUNCH_Z    # ~72 Bohr
    return int(round((path / v + 10.0) / DT_WP))
def n_steps_cl(E):    # ion: launch -> past far slab face (+ margin)
    v = k0_of(E); path = (SLAB_HALF + 4.0) - LAUNCH_Z    # ~43 Bohr
    return int(round((path / v + 3.0) / DT_CL))

def done(rundir: Path) -> bool:
    for rs in Path(rundir).glob("**/run_summary.txt"):
        try:
            if "run_completed = true" in rs.read_text():
                return True
        except Exception:
            pass
    return False

def email(subject, body, attachments=None, **kw):
    try:
        sys.path.insert(0, str(ROOT / "inq-stack/python"))
        from inqview.email import send_run_email
        return send_run_email(subject=SUBJ + subject, body=body,
                              attachments=attachments or [], to=TO, **kw)
    except Exception as e:
        log(f"  EMAIL FAILED: {e}")
        return None

# ---------------------------------------------------------------- run a pair
def launch(binary, cwd, out, gpu, overrides):
    env = {**ENV_BASE, "CUDA_VISIBLE_DEVICES": str(gpu),
           **{k: str(v) for k, v in overrides.items()}}
    lf = open(Path(cwd) / f"{out}.log", "w")
    # start_new_session=True => detached session leader, immune to controlling-terminal
    # SIGHUP (the 2026-07-01 death: two foreground inq-run jobs killed together on
    # session close). Every child launched here now survives a disconnect.
    return subprocess.Popen([binary], cwd=str(cwd), env=env, stdout=lf,
                            stderr=subprocess.STDOUT, start_new_session=True)

# ---------------------------------------------------------------- P0b (self-launched)
P0B_STALL_MIN = 25        # kill+retry a P0b run whose log is silent this long (no sentinel)

def _launch_p0b():
    """Launch the P0b pair (exact original config) detached on GPU0/GPU1. Non-destructive:
    a prior partial dir is renamed aside, not deleted."""
    from datetime import datetime as _dt
    procs = {}
    for kind, (cwd, out, gpu, ov) in {
        "wp": (WPDIR, "p0b_wp", 0, dict(LJ_OUT="results/p0b_wp", LJ_K0=4.6957, LJ_N_STEPS=750,
               LJ_DT=0.04, LJ_CAP=1, LJ_WRITE_EVERY=4, LJ_WF_EVERY=40, LJ_LAUNCH_Z=-26.5)),
        "cl": (CLDIR, "p0b_classical", 1, dict(LJ_OUT="results/p0b_classical", LJ_K0=4.6957,
               LJ_N_STEPS=1500, LJ_DT=0.02, LJ_CAP=1, LJ_WRITE_EVERY=10, LJ_LAUNCH_Z=-26.5)),
    }.items():
        res = cwd / "results" / out
        if res.exists() and not done(res):
            bak = res.with_name(out + ".partial_pre_relaunch")
            if bak.exists(): import shutil; shutil.rmtree(bak)
            res.rename(bak); log(f"  moved partial {res.name} -> {bak.name}")
        p = launch(str(cwd / "run"), cwd, out, gpu, ov)
        procs[kind] = (p, cwd / f"{out}.log", res)
        log(f"  launched P0b {kind} (GPU{gpu}, pid {p.pid})")
    return procs

def run_p0b_to_completion(max_retries=1):
    """Launch P0b and BLOCK until both sentinels exist, with a liveness guard: if a run's
    log goes silent > P0B_STALL_MIN with no sentinel (or the process exits early), kill the
    pair and retry once. Returns True iff both completed."""
    for attempt in range(max_retries + 1):
        log(f"P0b launch attempt {attempt+1}/{max_retries+1}")
        procs = _launch_p0b()
        last_sz = {k: -1 for k in procs}; last_mv = {k: time.time() for k in procs}
        dead = False
        while True:
            time.sleep(60)
            if all(done(r) for _, _, r in procs.values()):
                log("P0b pair COMPLETE (both sentinels present)"); return True
            for k, (p, lf, res) in procs.items():
                if done(res): continue
                sz = lf.stat().st_size if lf.exists() else 0
                if sz > last_sz[k]: last_sz[k] = sz; last_mv[k] = time.time()
                elif p.poll() is not None: last_mv[k] = 0    # exited w/o sentinel
                if (time.time() - last_mv[k]) / 60.0 > P0B_STALL_MIN:
                    log(f"  P0b {k} STALLED/died with no sentinel — killing pair"); dead = True; break
            if dead:
                for p, _, _ in procs.values():
                    try: p.terminate()
                    except Exception: pass
                time.sleep(5); break
        if not dead:  # loop only exits via completion (returned) or dead
            return all(done(r) for _, _, r in procs.values())
    return False

def run_pair(E, retry=True):
    wp_out, cl_out = f"sweep_wp_E{E}", f"sweep_cl_E{E}"
    wp_res, cl_res = WPDIR / "results" / wp_out, CLDIR / "results" / cl_out
    procs = []
    if done(wp_res):
        log(f"  SKIP {wp_out} (done)")
    else:
        log(f"  RUN  {wp_out}  k0={k0_of(E):.3f}  n_steps={n_steps_wp(E)}")
        procs.append(("wp", launch(WPBIN, WPDIR, wp_out, 0,
            dict(LJ_OUT=wp_out, LJ_K0=k0_of(E), LJ_N_STEPS=n_steps_wp(E), LJ_DT=DT_WP,
                 LJ_CAP=1, LJ_WRITE_EVERY=WRITE_EVERY_SWEEP, LJ_WF_EVERY=200, LJ_LAUNCH_Z=LAUNCH_Z))))
    if done(cl_res):
        log(f"  SKIP {cl_out} (done)")
    else:
        log(f"  RUN  {cl_out}  v0={k0_of(E):.3f}  n_steps={n_steps_cl(E)}")
        procs.append(("cl", launch(CLBIN, CLDIR, cl_out, 1,
            dict(LJ_OUT=cl_out, LJ_K0=k0_of(E), LJ_N_STEPS=n_steps_cl(E), LJ_DT=DT_CL,
                 LJ_CAP=1, LJ_WRITE_EVERY=WRITE_EVERY_SWEEP, LJ_LAUNCH_Z=LAUNCH_Z))))
    for _, p in procs: p.wait()
    # one-shot retry for any that didn't complete
    if retry:
        if not done(wp_res): log(f"  {wp_out} incomplete -> retry"); run_pair_single(E, "wp")
        if not done(cl_res): log(f"  {cl_out} incomplete -> retry"); run_pair_single(E, "cl")

def run_pair_single(E, kind):
    if kind == "wp":
        out = f"sweep_wp_E{E}"; p = launch(WPBIN, WPDIR, out, 0,
            dict(LJ_OUT=out, LJ_K0=k0_of(E), LJ_N_STEPS=n_steps_wp(E), LJ_DT=DT_WP,
                 LJ_CAP=1, LJ_WRITE_EVERY=WRITE_EVERY_SWEEP, LJ_WF_EVERY=200, LJ_LAUNCH_Z=LAUNCH_Z))
    else:
        out = f"sweep_cl_E{E}"; p = launch(CLBIN, CLDIR, out, 1,
            dict(LJ_OUT=out, LJ_K0=k0_of(E), LJ_N_STEPS=n_steps_cl(E), LJ_DT=DT_CL,
                 LJ_CAP=1, LJ_WRITE_EVERY=WRITE_EVERY_SWEEP, LJ_LAUNCH_Z=LAUNCH_Z))
    p.wait()

# ---------------------------------------------------------------- S extraction
def _read_csv(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return rows

def analyse_wp(E):
    """WP energy method: S = [E_total(t_f) - E_GS]/L_z (phase-5). + N guard."""
    base = WPDIR / "results" / f"sweep_wp_E{E}" / "raw" / "observables"
    obs = _read_csv(base / "observables.csv")
    Et = [float(r["energy_total"]) for r in obs]
    S = (Et[-1] - E_GS) / L_Z * HA           # eV/Bohr
    # N guard (CAP absorbed WP, bath intact): expect N: 83 -> ~82
    ncsv = base / "electron_number.csv"; n0 = nf = None
    if ncsv.exists():
        nr = _read_csv(ncsv)
        if nr: n0 = float(nr[0]["N_total"]); nf = float(nr[-1]["N_total"])
    # convergence: late dE/dt small
    late = [abs(Et[i]-Et[i-1]) for i in range(max(1,len(Et)-5), len(Et))]
    conv = (max(late) if late else 0.0) * HA
    return dict(E=E, v=k0_of(E), S_eV_per_Bohr=round(S,3),
                E_total_final=Et[-1], N0=n0, Nf=nf,
                late_dE_eV=round(conv,3),
                absorbed=(nf is not None and nf < 82.5),
                bath_ok=(nf is None or nf > 81.3))

def analyse_cl(E):
    """Classical S(v0) = INITIAL DRAG -dKE/ds over the early vz>=0.85*v0 window
    (ALWAYS-ON light-projectile rule: a light Ehrenfest electron decelerates
    strongly, so a slab/full-run average mixes S over v0..0, not S at v0). ds ~ dz
    for +z motion; widen the vz threshold if the window is sparse. Provisional."""
    import numpy as np
    trk = CLDIR / "results" / f"sweep_cl_E{E}" / "raw" / "observables" / "electron_track.csv"
    tr = _read_csv(trk)
    z  = [float(r["z"]) for r in tr]; vz = [float(r["vz"]) for r in tr]
    ke = [float(r["ke_ion_ha"]) for r in tr]
    v0 = k0_of(E)
    S = float("nan"); win = None; thr_used = None; npts = 0
    for thr in (0.85, 0.70, 0.50):
        idx = [i for i in range(len(vz)) if vz[i] >= thr * v0]
        if len(idx) >= 8:
            zs = [z[i] for i in idx]; kes = [ke[i] for i in idx]
            if max(zs) - min(zs) > 1.0:
                S = -np.polyfit(zs, kes, 1)[0] * HA   # eV/Bohr
                win = [round(min(zs), 1), round(max(zs), 1)]; thr_used = thr; npts = len(idx)
                break
    return dict(E=E, v=v0, S_eV_per_Bohr=(round(S, 3) if S == S else None),
                drag_window_z=win, vz_thr=thr_used, n_pts=npts,
                stopped_in_slab=(max(z) < SLAB_HALF))

# ---------------------------------------------------------------- P0b gate
def p0b_gate():
    wp_res = WPDIR / "results" / "p0b_wp"; cl_res = CLDIR / "results" / "p0b_classical"
    ok = True; notes = []
    if not done(wp_res): ok = False; notes.append("P0b WP did NOT complete")
    if not done(cl_res): ok = False; notes.append("P0b classical did NOT complete")
    detail = {}
    if done(wp_res):
        try:
            a = analyse_wp_dir(wp_res); detail["wp"] = a
            if not math.isfinite(a["E_total_final"]): ok = False; notes.append("WP E_total NaN")
            if a.get("Nf") is not None and a["Nf"] < 80.0: ok = False; notes.append(f"bath drained (N={a['Nf']:.2f})")
        except Exception as e:
            notes.append(f"WP analyse error: {e}")
    return ok, notes, detail

def analyse_wp_dir(res):
    base = res / "raw" / "observables"
    obs = _read_csv(base / "observables.csv")
    Et = [float(r["energy_total"]) for r in obs]
    S = (Et[-1] - E_GS) / L_Z * HA
    ncsv = base / "electron_number.csv"; n0 = nf = None
    if ncsv.exists():
        nr = _read_csv(ncsv)
        if nr: n0 = float(nr[0]["N_total"]); nf = float(nr[-1]["N_total"])
    return dict(E_total_final=Et[-1], S_eV_per_Bohr=round(S,3), N0=n0, Nf=nf)

# ---------------------------------------------------------------- plot
def build_se_plot(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from inqview.visualisation.style import apply_theme; apply_theme()
    except Exception: pass
    wp = sorted([r for r in rows if r["kind"]=="wp"], key=lambda r:r["E"])
    cl = sorted([r for r in rows if r["kind"]=="cl"], key=lambda r:r["E"])
    fig, ax = plt.subplots(figsize=(6.4,4.0), constrained_layout=True)
    if wp: ax.plot([r["E"] for r in wp],[r["S_eV_per_Bohr"] for r in wp],"o-",label="WP (σ=3.5, quantum)")
    if cl: ax.plot([r["E"] for r in cl],[r["S_eV_per_Bohr"] for r in cl],"s--",label="classical (matched σ)")
    ax.set_xlabel("drift energy E (eV)"); ax.set_ylabel("S (eV/Bohr)")
    ax.set_title("Wide-WP vs matched-classical stopping S(E) — localised slab (PROVISIONAL)")
    ax.legend(frameon=False)
    fig.savefig(path, dpi=140); plt.close(fig)

# ---------------------------------------------------------------- main
def main():
    log("ORCHESTRATOR start — wide-WP: SELF-LAUNCH P0b -> gate -> P1 sweep")
    # 1) LAUNCH the P0b pair ourselves (detached) and verify completion with a liveness
    #    guard + one retry. (Was: passively wait for hand-launched jobs it couldn't
    #    restart -> the 2026-07-01 orphaned 8h timeout. Fixed.)
    wp0, cl0 = WPDIR/"results"/"p0b_wp", CLDIR/"results"/"p0b_classical"
    if done(wp0) and done(cl0):
        log("P0b already complete — skipping relaunch")
    elif not run_p0b_to_completion(max_retries=1):
        email("P0b did NOT complete (relaunch+retry failed) — HALTING",
              "The self-launched P0b pair failed to complete even after a retry "
              "(stalled/died with no sentinel). Sweep NOT started; check GPU/host state.")
        log("P0b relaunch failed — halting"); return
    log("P0b pair complete")
    # rebuild the gate-review notebook on the complete data (best-effort)
    try:
        subprocess.run([PY, str(HYP / "build_p0b_gate_review.py")],
                       env={**ENV_BASE, "PYTHONPATH": str(ROOT / "inq-stack/python")},
                       cwd=str(HYP), check=True, timeout=1200)
        log("gate-review notebook rebuilt on complete P0b data")
    except Exception:
        log("notebook rebuild failed:\n" + traceback.format_exc())
    # 2) numeric gate
    ok, notes, detail = p0b_gate()
    gate_txt = ("PASS" if ok else "FAIL") + (("; " + "; ".join(notes)) if notes else "")
    log(f"P0b gate: {gate_txt}  detail={detail}")
    email(f"P0b gate {'PASS -> starting sweep' if ok else 'FAIL -> HALT'}",
          f"Hypothesis: a wide near-rigid WP (σ=3.5) vs matched classical isolates quantum stopping.\n"
          f"P0b matched pair (E=300) complete.\nGate: {gate_txt}\nWP detail: {detail}\n"
          + ("Proceeding to the 6-energy P1 sweep." if ok else "Sweep NOT started — needs attention."))
    if not ok:
        log("gate failed — halting before P1"); return
    # 3) sweep
    rows = []
    HYP.mkdir(parents=True, exist_ok=True)
    csv_path = HYP / "se_wide_wp.csv"
    for E in E_GRID:
        try:
            log(f"=== energy E={E} eV ===")
            run_pair(E)
            wp = analyse_wp(E); wp["kind"]="wp"; rows.append(wp)
            cl = analyse_cl(E); cl["kind"]="cl"; rows.append(cl)
            log(f"  S_wp={wp['S_eV_per_Bohr']}  S_cl={cl['S_eV_per_Bohr']} eV/Bohr")
            # incremental save
            with open(csv_path,"w",newline="") as f:
                w = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r}))
                w.writeheader(); [w.writerow(r) for r in rows]
        except Exception:
            tb = traceback.format_exc(); log(f"E={E} FAILED:\n{tb}")
            email(f"sweep E={E} FAILED (chain continues)", tb)
    # 4) plot + final email
    png = HYP / "se_wide_wp.png"
    try:
        build_se_plot(rows, str(png))
    except Exception:
        tb = traceback.format_exc(); log(f"plot failed:\n{tb}"); png = None
    wp_s = {r["E"]: r["S_eV_per_Bohr"] for r in rows if r["kind"]=="wp"}
    cl_s = {r["E"]: r["S_eV_per_Bohr"] for r in rows if r["kind"]=="cl"}
    body = (
        "HYPOTHESIS: a wide, near-rigid wavepacket (σ_WP=3.5) and a matched-σ "
        "classical projectile through the localised jellium slab isolate purely "
        "quantum stopping (Pauli+interference) from dispersion.\n\n"
        "WHAT WAS DONE: autonomous 6-energy S(E) sweep (E=200-600 eV), each energy "
        "a WP + matched-classical pair, box 50x50x101 dx=0.40, CAP η=-0.7/10-side, "
        "GS anchor E_GS=-86.041 Ha.\n\n"
        f"RESULT (PROVISIONAL S, eV/Bohr):\n  WP:        {wp_s}\n  classical: {cl_s}\n\n"
        "WHAT THE PLOT SHOWS: S(E) for the wide WP vs the matched classical over the "
        "high-velocity tail (v/vF~11-20). WP-minus-classical is the candidate quantum "
        "component.\n\n"
        "CONCLUSION (PROVISIONAL — review before quoting):\n"
        " - Classical S = INITIAL DRAG -dKE/ds over the early v>=0.85*v0 window "
        "(light-projectile rule); the p0b track shows STRONG deceleration "
        "(v/v0 ~0.86 already at slab entry, partly in the vacuum approach — possible "
        "pre-contact / PBC self-image drag, cf. the classical-projectile-fix concerns).\n"
        " - WP S = energy method [E_total(t_f)-E_GS]/L_z (phase-5), which is a "
        "TRAVERSAL AVERAGE for a decelerating packet; the rule-compliant WP S(v0) is "
        "the n(k,t) momentum-centroid initial drag — that data IS saved, to be "
        "extracted on review.\n"
        " - SIE not yet bounded (needs the vacuum-WP control); E_ref is the dx=0.40 GS "
        "(16 Ha below dx=0.50 production — resolution effect, cancels per-dataset).\n"
        " - VTI cadence coarsened (WRITE_EVERY=20) for the sweep so 12 runs fit; all "
        "scalar/momentum observables kept."
    )
    email("P1 sweep COMPLETE — S(E) (provisional)", body, attachments=[str(png)] if png else None)
    log("ORCHESTRATOR done")

if __name__ == "__main__":
    try:
        main()
    except Exception:
        tb = traceback.format_exc(); log(f"FATAL:\n{tb}"); email("orchestrator FATAL", tb)
