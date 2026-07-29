#!/usr/bin/env python3
"""Autonomous orchestrator — Nazarov-Gross fixed-velocity mass sweep, Phase 1.

Campaign: docs/campaigns/nazarov_gross_comparison/nazarov-gross-comparison.md
Budget:   14 h wall x 2 GPUs, enforced from T0 (orchestrator start).

Pipeline (each stage gated; a failure emails full context + stops):
  0. compile gates: gs + wp binaries build against inq-study
  1. GS35 [GPU0] || GS40 [GPU1]  (skip any that already exists)
  2. cutoff guard on every rung (BLOCK -> stop)
  3. smoke: m=2.2 (worst-case k0), 60 steps at h=0.35 -> completed, finite,
     |E_total drift| < 1e-3 Ha; measures the per-step wall cost
  4. null branch: (m0.5 GPU0 || m0.71 GPU1) then (m1.41 GPU0 || m2.2 GPU1),
     880 steps, v=2.711, h=0.35; one retry per run
  5. slow pilots (budget-gated sacrificial tail): m=1 GPU0 || m=10 GPU1 at
     v=0.25, h=0.40, launch_z=-13.5; N_STEPS sized to remaining wall
     (cap 1500, floor 800; below floor -> deferred + emailed)
  6. final email + flip campaign frontmatter done-flags

Idempotent: reruns skip any run whose run_summary.txt has run_completed=true.
Headless:
  cd .../scripts/nazarov_gross && setsid nohup \
    ../../../../../venv/bin/python3 orchestrate.py > orch.log 2>&1 &
"""
from __future__ import annotations
import math, os, re, subprocess, sys, time
from datetime import datetime
from pathlib import Path

ROOT  = Path("/local/data/public/skcb2/tddft")
LJ    = ROOT/"ResearchProject/systems/localised_jellium"
HERE  = LJ/"scripts/nazarov_gross"
GS35  = LJ/"shared_gs/slab_n234_L50_h0p35"
GS40  = LJ/"shared_gs/slab_n234_L50_h0p40"
CAMP  = ROOT/"docs/campaigns/nazarov_gross_comparison/nazarov-gross-comparison.md"
GUARD = Path("/home/raid/skcb2/skcb2/tddft/.claude/skills/tddft-simulations/cutoff_guard.py")
PY    = str(ROOT/"venv/bin/python3")
TO    = "chiddukanna@gmail.com"

ENV = {**os.environ,
       "PATH":                 f"{ROOT/'shared/bin'}:{os.environ.get('PATH','')}",
       "INQ_SOURCE":           str(ROOT/"inq-study"),
       "INQ_SHARE_PATH":       str(ROOT/"inq/install/share"),
       "PSEUDOPOD_SHARE_PATH": str(ROOT/"inq/install/share/pseudopod")}

V_FAST   = 2.7110633401       # p3 electron velocity (k0 for m=1)
V_SLOW   = 0.25               # 0.52 * v_F(r_s=3.996)
SIGMA_WP = 0.5
HA_EV    = 27.211386245988
BUDGET_S = 14*3600
T0       = time.time()

# RE-PLAN 2026-07-12 (user-locked): h=0.35 unusable on 24 GB GPUs (memory-thrash,
# ~260 s/step effective vs 16 predicted; see wp/smoke.log). Null branch moved to
# h=0.40 on GS40; m=2.2 dropped (aliasing BLOCK at 0.40) -> 3 rungs. Pilot m=10
# rides GPU1 in round 2; pilot m=1 is the budget-sized tail.
NULL_RUNS = [  # (name, mass, gpu, round)
    ("null_m0p5",  0.50, "0", 1), ("null_m0p71", 0.71, "1", 1),
    ("null_m1p41", 1.41, "0", 2)]
PILOT_RUNS = [("pilot_slow_m1", 1.0, "0"), ("pilot_slow_m10", 10.0, "1")]
PILOT_STEPS_TARGET = 1400

def log(m): print(f"[{datetime.now():%F %T}] [{(time.time()-T0)/3600:5.2f}h] {m}", flush=True)
def remaining(): return BUDGET_S - (time.time() - T0)

def email(subj, body, attach=None):
    try:
        sys.path.insert(0, str(ROOT/"inq-stack/python"))
        from inqview.email import send_run_email
        send_run_email(subject=f"[nazarov-gross] {subj}", body=body,
                       attachments=[a for a in (attach or []) if Path(a).exists()], to=TO)
        log(f"emailed: {subj}")
    except Exception as e:
        log(f"  email failed ({subj}): {e}")

def done(rundir: Path) -> bool:
    rs = rundir/"run_summary.txt"
    try:    return "run_completed = true" in rs.read_text()
    except Exception: return False

def flip_task(prefix: str):
    """Flip `done: false` -> `done: true` on the frontmatter task starting with prefix."""
    try:
        t = CAMP.read_text()
        pat = re.compile(r'(\{ name: "' + re.escape(prefix) + r'[^"]*", done: )false')
        CAMP.write_text(pat.sub(r"\1true", t, count=1))
    except Exception as e:
        log(f"  frontmatter flip failed ({prefix}): {e}")

# ---------------------------------------------------------------- helpers
def build(sub: str, probe_env: dict) -> bool:
    d = HERE/sub
    if (d/"run").exists(): return True
    with open(d/"build.log", "w") as f:
        subprocess.run(["inq-run"], cwd=d, env={**ENV, "CUDA_VISIBLE_DEVICES": "0", **probe_env},
                       stdout=f, stderr=subprocess.STDOUT)
    return (d/"run").exists()

def launch(sub: str, gpu: str, extra: dict, logname: str):
    d = HERE/sub
    f = open(d/logname, "w")
    p = subprocess.Popen([str(d/"run")], cwd=d,
                         env={**ENV, "CUDA_VISIBLE_DEVICES": gpu, **extra},
                         stdout=f, stderr=subprocess.STDOUT)
    log(f"{sub}:{logname} launched on GPU {gpu} (pid {p.pid})")
    return p

def guard_rung(mass: float, v: float, spacing: float) -> dict:
    """cutoff_guard with the p0-equivalent drift energy: its WP formula takes
    p0 = sqrt(2 E) (m=1 convention), so feed E_eff = (m v)^2 / 2."""
    sys.path.insert(0, str(GUARD.parent))
    from cutoff_guard import check_run
    e_eff_ev = 0.5 * (mass*v)**2 * HA_EV
    return check_run(spacing_bohr=spacing, kind="wp", energy_ev=e_eff_ev, sigma_wp_bohr=SIGMA_WP)

def obs_energy_ok(rundir: Path, tol_ha: float = 1e-3):
    """observables.csv: energy_total finite and |drift| < tol."""
    f = rundir/"raw/observables/observables.csv"
    if not f.exists(): return False, "observables.csv missing"
    et = []
    with open(f) as fh:
        hdr = fh.readline().strip().split(",")
        try: i = hdr.index("energy_total")
        except ValueError: return False, "no energy_total column"
        for line in fh:
            try: et.append(float(line.split(",")[i]))
            except Exception: pass
    if not et: return False, "no rows"
    if any(math.isnan(x) or math.isinf(x) for x in et): return False, "NaN/inf energy"
    drift = max(et) - min(et)
    return drift < tol_ha, f"E_total drift {drift:.3e} Ha over {len(et)} rows (tol {tol_ha})"

def run_summary_field(rundir: Path, key: str):
    try:
        for line in (rundir/"run_summary.txt").read_text().splitlines():
            if line.startswith(key):
                return float(line.split("=")[1].split()[0])
    except Exception: pass
    return None

# ---------------------------------------------------------------- pipeline
def main():
    log(f"START budget {BUDGET_S/3600:.0f} h x 2 GPUs")

    # 0. compile gates -------------------------------------------------------
    if not build("gs", {"NG_COMPILE_PROBE": "1"}):
        email("gs BUILD FAILED", f"see {HERE/'gs/build.log'}. STOPPED."); return
    if not build("wp", {"NG_GS_DIR": "/nonexistent_gs_compile_probe"}):
        email("wp BUILD FAILED", f"see {HERE/'wp/build.log'}. STOPPED."); return
    log("binaries built (gs, wp) [engine inq-study]")

    # 1. ground states (GS35 on GPU0 || GS40 on GPU1) ------------------------
    need = [(GS35, "0.35", "0"), (GS40, "0.40", "1")]
    procs = []
    for ckpt, h, gpu in need:
        if done(ckpt):
            log(f"GS h={h} present, skip"); continue
        procs.append((ckpt, h, launch("gs", gpu, {"NG_SPACING": h, "NG_CKPT": str(ckpt)},
                                      f"gs_h{h.replace('.','p')}.log")))
    for ckpt, h, p in procs:
        p.wait()
        if not done(ckpt):
            email(f"GS h={h} FAILED", f"rc={p.returncode}; see gs logs in {HERE/'gs'}. STOPPED."); return
    e35 = run_summary_field(GS35, "ground_state_energy_ha")
    e40 = run_summary_field(GS40, "ground_state_energy_ha")
    log(f"GS ready: E(h=0.35)={e35} Ha, E(h=0.40)={e40} Ha")
    flip_task("Phase 1a")
    email("GS ready", f"slab_n234 ground states converged.\nE(h=0.35) = {e35} Ha\n"
          f"E(h=0.40) = {e40} Ha\nNext: cutoff guards + m=2.2 smoke, then the 4 null runs.")

    # 2. cutoff guards --------------------------------------------------------
    lines = []
    for name, m, _, _ in NULL_RUNS:
        r = guard_rung(m, V_FAST, 0.40)
        lines.append(f"{name}: {r['status']} ({r.get('reason','')})")
        if r["block"]:
            email("CUTOFF GUARD BLOCK", "\n".join(lines) + "\nSTOPPED."); return
    for name, m, _ in PILOT_RUNS:
        r = guard_rung(m, V_SLOW, 0.40)
        lines.append(f"{name}: {r['status']} ({r.get('reason','')})")
        if r["block"]:
            email("CUTOFF GUARD BLOCK (pilot)", "\n".join(lines) + "\nSTOPPED."); return
    log("cutoff guards passed:\n  " + "\n  ".join(lines))

    # 3. smoke: worst-case h=0.40 rung (m=1.41), 60 steps ----------------------
    smoke = HERE/"wp/results/smoke_m1p41_h0p40"
    if not done(smoke):
        p = launch("wp", "0", {"NG_OUT": "smoke_m1p41_h0p40", "NG_MASS": "1.41",
                               "NG_N_STEPS": "60", "NG_SPACING": "0.40",
                               "NG_GS_DIR": str(GS40)}, "smoke_h0p40.log")
        p.wait()
    if not done(smoke):
        email("SMOKE FAILED", "m=1.41 60-step h=0.40 smoke did not complete; "
              "see wp/smoke_h0p40.log. STOPPED."); return
    ok, msg = obs_energy_ok(smoke)
    if not ok:
        email("SMOKE GATE FAILED", f"m=1.41 h=0.40 smoke energy check: {msg}. STOPPED."); return
    wall = run_summary_field(smoke, "wall_time_s") or 0.0
    step40 = wall/60.0
    null_est_s = 880 * step40
    log(f"smoke passed ({msg}); step cost h=0.40 ~ {step40:.1f} s/step "
        f"-> null run est {null_est_s/3600:.1f} h each")
    # Budget policy (user decision 2026-07-12): NEVER self-block on a projected
    # overrun. Runs checkpoint every NG_CKPT_EVERY steps, so the user can kill a
    # run at any time and resume with NG_RESUME=1 — losing at most one
    # checkpoint interval. We WARN with the projection and proceed at full scope.
    if 2*null_est_s + 1800 > remaining():
        email("BUDGET OVERRUN PROJECTED — proceeding (checkpointed)",
              f"Measured {step40:.0f} s/step at h=0.40 -> {null_est_s/3600:.1f} h per "
              f"880-step run; the full schedule exceeds the remaining "
              f"{remaining()/3600:.1f} h. Proceeding anyway per the checkpoint-"
              f"don't-block policy: every run saves an RT checkpoint every "
              f"200 steps (results/<run>/rt_ckpt). To stop: kill the run "
              f"process; to continue later: rerun with NG_RESUME=1.")
        log("budget overrun projected — proceeding (checkpointed)")
    flip_task("Phase 1b")

    # 4. null branch (h=0.40, GS40). Round 1: m0.5||m0.71. Round 2: m1.41 on
    #    GPU0 with pilot m=10 riding GPU1 if the budget supports its full target.
    def null_env(m):
        return {"NG_MASS": str(m), "NG_SPACING": "0.40", "NG_GS_DIR": str(GS40)}
    def pilot_env(m, n_steps):
        return {"NG_MASS": str(m), "NG_V": str(V_SLOW), "NG_SPACING": "0.40",
                "NG_GS_DIR": str(GS40), "NG_LAUNCH_Z": "-13.5", "NG_N_STEPS": str(n_steps)}

    for rnd in (1, 2):
        procs = []
        for name, m, gpu, r in NULL_RUNS:
            if r != rnd or done(HERE/"wp/results"/name): continue
            procs.append((name, null_env(m), gpu,
                          launch("wp", gpu, {"NG_OUT": name, **null_env(m)}, f"{name}.log")))
        if rnd == 2 and not done(HERE/"wp/results/pilot_slow_m10"):
            procs.append(("pilot_slow_m10", pilot_env(10.0, PILOT_STEPS_TARGET), "1",
                          launch("wp", "1", {"NG_OUT": "pilot_slow_m10",
                                             **pilot_env(10.0, PILOT_STEPS_TARGET)},
                                 "pilot_slow_m10.log")))
        for name, env_, gpu, p in procs:
            p.wait()
            if not done(HERE/"wp/results"/name):           # one-shot retry
                log(f"{name} failed (rc={p.returncode}) — retrying once")
                p2 = launch("wp", gpu, {"NG_OUT": name, **env_}, f"{name}.retry.log")
                p2.wait()
            ok = done(HERE/"wp/results"/name)
            log(f"{name} completed={ok}")
            if not ok:
                email(f"{name} FAILED twice", f"see wp/{name}.log + retry log. Continuing with the rest.")
        email(f"round {rnd}/2 done",
              "\n".join(f"{n}: completed={done(HERE/'wp/results'/n)}, "
                        f"wall={run_summary_field(HERE/'wp/results'/n,'wall_time_s')}"
                        for n, _, _, _ in procs) or "(nothing to run)")
    if all(done(HERE/"wp/results"/n) for n, _, _, _ in NULL_RUNS):
        flip_task("Phase 1c")

    # 5. remaining pilot (m=1) — full target, checkpointed (never self-blocked) --
    if not done(HERE/"wp/results/pilot_slow_m1"):
        est = PILOT_STEPS_TARGET*step40/3600
        log(f"pilot_slow_m1: {PILOT_STEPS_TARGET} steps (~{est:.1f} h; "
            f"remaining budget {remaining()/3600:.1f} h; checkpointed, user may kill)")
        p = launch("wp", "0", {"NG_OUT": "pilot_slow_m1",
                               **pilot_env(1.0, PILOT_STEPS_TARGET)}, "pilot_slow_m1.log")
        p.wait(); log(f"pilot_slow_m1 completed={done(HERE/'wp/results/pilot_slow_m1')}")
    if all(done(HERE/"wp/results"/n) for n, _, _ in PILOT_RUNS):
        flip_task("Phase 1d")
    final_email(deferred=False)

def final_email(deferred: bool):
    rows = []
    for n, m, _, _ in NULL_RUNS:
        rows.append(f"{n} (m={m}, v={V_FAST:.3f}, h=0.40): completed={done(HERE/'wp/results'/n)} "
                    f"wall={run_summary_field(HERE/'wp/results'/n, 'wall_time_s')}")
    for n, m, _ in PILOT_RUNS:
        rows.append(f"{n} (m={m}, v={V_SLOW}, h=0.40): completed={done(HERE/'wp/results'/n)} "
                    f"wall={run_summary_field(HERE/'wp/results'/n, 'wall_time_s')}")
    email("PHASE 1 RUNS FINISHED" + (" (pilots deferred)" if deferred else ""),
          "Nazarov-Gross fixed-velocity mass sweep — run stage complete.\n\n"
          "Hypothesis reminder: at v=2.711 (5.6 v_F, HIGH-v regime) NG predicts S is\n"
          "FLAT across masses (null branch); mass-dependent friction lives at v<v_F\n"
          "(slow pilots de-risk that next phase).\n\n"
          "What was done:\n" + "\n".join(rows) + "\n\n"
          f"Results under {HERE/'wp/results'}. Analysis (retained-energy ledgers, S(m)\n"
          "figure, spreading-systematic check, pilot initial-drag verdicts) is the next\n"
          "session's task 1e — notebooks not yet built, so no plot attached here.\n\n"
          "Conclusion: pending analysis; no NG verdict is claimable from runs alone.")
    log("DONE")

if __name__ == "__main__":
    try:
        main()
    except Exception:
        import traceback
        tb = traceback.format_exc()
        log(tb)
        email("ORCHESTRATOR CRASHED", tb)
        raise
