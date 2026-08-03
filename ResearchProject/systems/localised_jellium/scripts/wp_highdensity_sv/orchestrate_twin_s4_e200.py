#!/usr/bin/env python3
"""
Autonomous WP-vs-classical MATCHED TWIN at sigma_WP=4, 200 eV (density A, r_s=4.18).

One clean pairwise comparison at a width where the WP barely spreads during transit
(free-dispersion floor: +15% across the slab, +17% by far face; see
docs/plans/twin-sigma4-e200.md). Both halves share the SAME dx=0.5 GS/grid/box/CAP
=> a perfectly matched pair. Motivation: test "wide classical ~ WP" directly, after
the sigma=17/20 wide runs proved contaminated (S_A~0, S_B a static-overlap artefact).

Config (both halves):
  sigma_WP=4  E=200 eV -> v=k0=3.8340 (m=1)  launch_z=-24.5 (12 Bohr from near face)
  dt=0.04  dx=0.5  N_STEPS=1890 (t_final=75.6 a.u., 290-Bohr path)
  CAP eta=-1 Ha, 12.5 Bohr/face (+/-30..+/-42.5)  GS=slab_n100_L35x35x85_dx0p5_per2

Binaries:
  WP : scripts/wp_highdensity_sv/wp/run          (built against inq-study; LJ_CAP_L)
  CL : scripts/classical_highdensity_sv/dyn_direct_cap/run  (sweep binary; direct erf/r)

Dual-GPU, concurrent, detached children (start_new_session), checkpointed
(LJ_RESUME=1), idempotent (skips run_completed), liveness guard (kill+retry once on
stall), MAX_HOURS cap, per-run + final emails (checkpoint-dont-block).

START (detached, survives session):
  cd .../scripts/wp_highdensity_sv
  setsid nohup /local/data/public/skcb2/tddft/venv/bin/python3 orchestrate_twin_s4_e200.py \
      > orchestrate_twin_s4_e200.log 2>&1 < /dev/null &
"""
from __future__ import annotations
import os, sys, time, glob, signal, subprocess
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
LJ   = ROOT / "ResearchProject/systems/localised_jellium"
WPDIR = LJ / "scripts/wp_highdensity_sv/wp"
CLDIR = LJ / "scripts/classical_highdensity_sv/dyn_direct_cap"
GS    = LJ / "shared_gs/slab_n100_L35x35x85_dx0p5_per2"
HYP   = LJ / "hypotheses/wp_highdensity_sv/twin_s4_e200"
STACK = str(ROOT / "inq-stack/python")
TO    = "chiddukanna@gmail.com"

HA = 27.211386
E_GS = 207.18322156141      # Ha, bare slab GS (dx0p5) — same GS both halves
LSLAB = 25.0; FACE = 12.5

# ---- run config ----
SIGMA = 4.0
V     = 3.8340              # 200 eV, m=1
LAUNCH_Z = -24.5
DT    = 0.04
DX    = 0.5
N_STEPS = 1890             # ceil(7245/V)
SAVE_EVERY = 6            # ~300 density frames
WF_EVERY   = 18          # ~105 wavefunction/momentum frames
WP_OUT = "twin_s4_e200"
CL_OUT = "twin_s4_e200_cl"

POLL_S = 60; STALL_MIN = 25; MAX_HOURS = 12

SHARE = {"INQ_SHARE_PATH": str(ROOT / "inq/install/share"),
         "PSEUDOPOD_SHARE_PATH": str(ROOT / "inq/install/share/pseudopod"),
         "INQ_SOURCE": str(ROOT / "inq-study")}


def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)


def email(subject, body, attach=None):
    try:
        sys.path.insert(0, STACK)
        from inqview.email import send_run_email
        send_run_email(subject=subject, body=body,
                       attachments=[str(a) for a in (attach or []) if a and os.path.exists(a)], to=TO)
        log(f"  emailed: {subject}")
    except Exception as e:
        log(f"  EMAIL FAILED: {e}")


def wp_env(gpu):
    e = {**os.environ, **SHARE, "PATH": f"{ROOT}/shared/bin:" + os.environ.get("PATH", ""),
         "CUDA_VISIBLE_DEVICES": str(gpu),
         "LJ_SIGMA": f"{SIGMA}", "LJ_K0": f"{V}", "LJ_LAUNCH_Z": f"{LAUNCH_Z}",
         "LJ_SPACING": f"{DX}", "LJ_DT": f"{DT}", "LJ_N_STEPS": str(N_STEPS),
         "LJ_CAP_ETA": "-1.0", "LJ_CAP_L": f"{FACE}",
         "LJ_SAVE_EVERY": str(SAVE_EVERY), "LJ_WF_EVERY": str(WF_EVERY),
         "LJ_CKPT_EVERY": "0", "LJ_RESUME": "1",
         "LJ_OUT": WP_OUT, "LJ_GS_DIR": str(GS)}
    return e


def cl_env(gpu):
    e = {**os.environ, **SHARE, "PATH": f"{ROOT}/shared/bin:" + os.environ.get("PATH", ""),
         "CUDA_VISIBLE_DEVICES": str(gpu),
         "LJ_LX": "35", "LJ_LY": "35", "LJ_LZ": "85", "LJ_HALF": "12.5", "LJ_N": "100",
         "LJ_EDGE_W": "1.0", "LJ_PERIODICITY": "2", "LJ_SPACING": f"{DX}",
         "LJ_SIGMA": f"{SIGMA}", "LJ_MASS": "1.0", "LJ_DELTA": "0.1", "LJ_DT": f"{DT}",
         "LJ_CONST_V": "0", "LJ_LAUNCH_Z": f"{LAUNCH_Z}", "LJ_K0": f"{V}",
         "LJ_CAP_ETA": "-1.0", "LJ_CAP_INNER_BOHR": "30.0", "LJ_CAP_WIDTH_BOHR": f"{FACE}",
         "LJ_N_STEPS": str(N_STEPS), "LJ_SAVE_EVERY": str(SAVE_EVERY),
         "LJ_GS_DIR": str(GS), "LJ_OUT": CL_OUT, "LJ_RESUME": "1"}
    return e


def completed(rundir: Path, out: str) -> bool:
    s = rundir / "results" / out / "run_summary.txt"
    return s.exists() and "run_completed = true" in s.read_text()


def gpu_free(g):
    return not subprocess.run(["fuser", f"/dev/nvidia{g}"], capture_output=True, text=True).stdout.strip()


def launch(rundir: Path, env, label):
    lf = open(rundir / f"run_{label}.log", "a")
    log(f"launch {label} on GPU {env['CUDA_VISIBLE_DEVICES']} ({N_STEPS} steps)")
    return subprocess.Popen([str(rundir / "run")], cwd=str(rundir), env=env,
                            stdout=lf, stderr=subprocess.STDOUT, start_new_session=True)


def _cat(rundir: Path, out: str, stem: str):
    import pandas as pd
    fs = sorted(glob.glob(str(rundir / "results" / out / "raw/observables" / f"{stem}*.csv")))
    if not fs:
        return None
    return (pd.concat([pd.read_csv(f) for f in fs])
              .drop_duplicates("step").sort_values("step").reset_index(drop=True))


def extract_classical():
    pj = _cat(CLDIR, CL_OUT, "projectile"); ob = _cat(CLDIR, CL_OUT, "observables")
    ix = _cat(CLDIR, CL_OUT, "interactions")
    ke = lambda z: pj.loc[(pj.proj_z - z).abs().idxmin(), "energy_proj_ke"]
    S_A = (ke(-FACE) - ke(FACE)) * HA / LSLAB
    S_B = (ob.energy_total.iloc[-1] - E_GS) * HA / LSLAB
    return dict(S_A=S_A, S_B=S_B, v_final=float(pj.proj_vz.iloc[-1]),
                z_final=float(pj.proj_z.iloc[-1]),
                norm_final=float(ix.norm_slab.iloc[-1]) if ix is not None and "norm_slab" in ix else float("nan"))


def extract_wp():
    ob = _cat(WPDIR, WP_OUT, "observables")
    rs = _cat(WPDIR, WP_OUT, "wp_real_space_stats")
    mo = _cat(WPDIR, WP_OUT, "wp_momentum_stats")
    E_raw = float(ob.energy_total.iloc[-1])
    normWP = float(rs.norm_check.iloc[-1])      # physical WP norm fraction (=1 at t=0)
    T1 = float(mo.e_kin_ha.iloc[-1])
    E_corr = E_raw - T1 * (1.0 - normWP)
    return dict(S_deposit_raw=(E_raw - E_GS) * HA / LSLAB,
                S_deposit_corrected_naive=(E_corr - E_GS) * HA / LSLAB,
                norm_WP_final=normWP, T1_final_Ha=T1,
                sigma_z_final=float(rs.sigma_z2.iloc[-1] ** 0.5))


def wait_for_wp_binary(timeout_min=75):
    """Block until the WP binary is built AND its 5-step smoke run completed."""
    smoke = WPDIR / "results/smoke/run_summary.txt"
    t0 = time.time()
    while time.time() - t0 < timeout_min * 60:
        building = subprocess.run(["pgrep", "-f", "cmake --build"], capture_output=True).returncode == 0
        if (WPDIR / "run").exists() and smoke.exists() and "run_completed = true" in smoke.read_text():
            log("WP binary built + smoke passed"); return True
        if (WPDIR / "run").exists() and not building and smoke.exists():
            log("WP binary present, smoke sentinel found"); return True
        time.sleep(30)
    return False


def run_pair():
    HYP.mkdir(parents=True, exist_ok=True)
    if not wait_for_wp_binary():
        email("[twin-s4] ABORT: WP binary not ready", f"{WPDIR/'run'} not built / smoke did not pass in time"); return
    jobs = []  # (rundir, out, env_fn, label)
    if not completed(WPDIR, WP_OUT): jobs.append((WPDIR, WP_OUT, wp_env, "wp"))
    else: log("WP already complete — skip")
    if not completed(CLDIR, CL_OUT): jobs.append((CLDIR, CL_OUT, cl_env, "cl"))
    else: log("classical already complete — skip")

    email("[twin-s4] matched twin started (sigma_WP=4, 200 eV)",
          "One clean WP-vs-classical pairwise run at a barely-spreading width.\n\n"
          f"sigma_WP=4 (sigma_pot=2.83), E=200 eV (v=3.834), launch 12 Bohr from near face,\n"
          f"dt=0.04, dx=0.5, N_STEPS={N_STEPS} (t_final=75.6 a.u.), CAP eta=-1/12.5 per face,\n"
          f"SAME dx0.5 GS for both halves. Free-dispersion spread: +15% across the slab.\n\n"
          f"Running: {[l for *_ , l in jobs] or 'nothing (both complete)'}. Dual-GPU, checkpointed.")

    # assign a free GPU to each job (wait until enough are free)
    procs = {}     # gpu -> (proc, rundir, out, env_fn, label, last_step, last_change)
    retried = set()
    t0 = time.time()
    pending = list(jobs)
    while pending or procs:
        # fill free GPUs
        free = [g for g in (0, 1) if gpu_free(g) and g not in procs]
        while pending and free:
            rundir, out, env_fn, label = pending.pop(0)
            g = free.pop(0)
            p = launch(rundir, env_fn(g), label)
            procs[g] = [p, rundir, out, env_fn, label, -1, time.time()]
            time.sleep(5)
        if not procs:
            time.sleep(POLL_S); continue
        time.sleep(POLL_S)
        if time.time() - t0 > MAX_HOURS * 3600:
            log("MAX_HOURS exceeded — leaving runs; they are checkpointed")
            email("[twin-s4] MAX_HOURS cap hit", "Runs are checkpointed; resume with LJ_RESUME=1.")
            return
        for g, st in list(procs.items()):
            p, rundir, out, env_fn, label, last, changed = st
            if p.poll() is None:
                # liveness: track last step in the run log
                try:
                    lf = (rundir / f"run_{label}.log").read_text()
                    cur = 0
                    for line in reversed(lf.splitlines()):
                        if "step" in line:
                            import re
                            m = re.search(r"step\s+(\d+)", line)
                            if m: cur = int(m.group(1)); break
                    if cur != last:
                        st[5] = cur; st[6] = time.time()
                    elif time.time() - changed > STALL_MIN * 60:
                        log(f"STALL {label} @step {cur} — killing")
                        os.killpg(os.getpgid(p.pid), signal.SIGKILL)
                except Exception as e:
                    log(f"liveness {label}: {e}")
                continue
            # finished
            del procs[g]
            if completed(rundir, out):
                log(f"DONE {label}")
            else:
                if label not in retried:
                    retried.add(label)
                    log(f"{label} did not complete — retry once (RESUME)")
                    pending.append((rundir, out, env_fn, label))
                else:
                    log(f"{label} FAILED twice")
                    email(f"[twin-s4] {label} FAILED twice",
                          f"See {rundir}/run_{label}.log and results/{out}. Other half continues.")

    # ---- final extraction ----
    import pandas as pd
    lines = []
    cl = wp = None
    try:
        cl = extract_classical()
        lines.append(f"CLASSICAL (sigma_WP=4): S_A(KE-loss)={cl['S_A']:.3f}  S_B(E-abs)={cl['S_B']:.3f} eV/Bohr  "
                     f"v_final={cl['v_final']:.3f} (v0={V})  z_final={cl['z_final']:.0f}  norm_slab={cl['norm_final']:.3f}")
    except Exception as e:
        lines.append(f"classical extract failed: {e}")
    try:
        wp = extract_wp()
        lines.append(f"WP (sigma_WP=4): S_deposit_raw={wp['S_deposit_raw']:.3f}  "
                     f"S_deposit_corrected(naive)={wp['S_deposit_corrected_naive']:.3f} eV/Bohr  "
                     f"norm_WP_final={wp['norm_WP_final']:.3f}  T1_final={wp['T1_final_Ha']:.2f} Ha  "
                     f"sigma_z_final={wp['sigma_z_final']:.1f} Bohr")
    except Exception as e:
        lines.append(f"WP extract failed: {e}")
    if cl or wp:
        pd.DataFrame([{**({} if not cl else {f"cl_{k}": v for k, v in cl.items()}),
                       **({} if not wp else {f"wp_{k}": v for k, v in wp.items()})}]).to_csv(
            HYP / "twin_s4_e200_summary.csv", index=False)
    body = ("Matched WP-vs-classical twin complete (sigma_WP=4, 200 eV).\n\n" + "\n".join(lines) +
            "\n\nNOTE: WP corrected-S here is the NAIVE application of the sweep formula "
            "E_corr = E_raw - T1*(1-norm_WP); the definitive corrected S_deposit + the "
            "pairwise comparison figure (with n(z,t) spreading overlay) are a reviewed "
            "follow-up using the tested wp_hd_stopping machinery.\n"
            "Outputs: wp/results/twin_s4_e200/, dyn_direct_cap/results/twin_s4_e200_cl/.")
    log("FINAL:\n" + "\n".join(lines))
    email("[twin-s4] TWIN COMPLETE — sigma_WP=4, 200 eV", body,
          attach=[HYP / "twin_s4_e200_summary.csv"])


if __name__ == "__main__":
    log("=== matched twin sigma_WP=4, 200 eV (dual-GPU) ===")
    run_pair()
    log("=== done ===")
