#!/usr/bin/env python3
"""
Autonomous WP-vs-classical MATCHED TWINS (density A, r_s=4.18), dual-GPU queue.

Runs a QUEUE of twin configs; for each, WP + classical launch concurrently on the
two GPUs, both share the SAME dx=0.5 GS/grid/box/CAP (perfectly matched), both save
density frames so the total + induced density GIFs can be built.

Queue (priority order):
  1. cs_s2p0_v2p0 — CASE STUDY (user-chosen): sigma_WP=2, v=2.0, launch -24.0
     (matches the sweep run s2p0_v2p0). Re-run locally at dx=0.5 to regenerate frames
     for the case-study GIFs. n_steps=3623 (t_final=145 a.u.).
  2. s4_e200      — sigma_WP=4, 200 eV (v=3.834), launch -24.5 (12 Bohr standoff).
     The barely-spreading pairwise from the previous task. n_steps=1890.

Both: dt=0.04, dx=0.5, CAP eta=-1 Ha / 12.5 Bohr per face, GS slab_n100_L35x35x85_dx0p5_per2.

RESUME is AUTO-DETECTED (1 only if results/<out>/rt_state.txt exists) — fixes the
earlier bug where a hardcoded LJ_RESUME=1 aborted fresh runs. Checkpointed,
idempotent (skips run_completed), retry-once, liveness guard, MAX_HOURS cap,
per-run + final emails (checkpoint-dont-block).

START (detached):
  cd .../scripts/wp_highdensity_sv
  setsid nohup /local/data/public/skcb2/tddft/venv/bin/python3 orchestrate_case_twins.py \
      > orchestrate_case_twins.log 2>&1 < /dev/null &
"""
from __future__ import annotations
import os, sys, re, time, glob, signal, subprocess
from datetime import datetime
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
LJ   = ROOT / "ResearchProject/systems/localised_jellium"
WPDIR = LJ / "scripts/wp_highdensity_sv/wp"
CLDIR = LJ / "scripts/classical_highdensity_sv/dyn_direct_cap"
GS    = LJ / "shared_gs/slab_n100_L35x35x85_dx0p5_per2"
HYP   = LJ / "hypotheses/wp_highdensity_sv"
STACK = str(ROOT / "inq-stack/python")
TO    = "chiddukanna@gmail.com"

HA = 27.211386; E_GS = 207.18322156141; LSLAB = 25.0; FACE = 12.5
DT = 0.04; DX = 0.5
POLL_S = 60; STALL_MIN = 30; MAX_HOURS = 14
SHARE = {"INQ_SHARE_PATH": str(ROOT / "inq/install/share"),
         "PSEUDOPOD_SHARE_PATH": str(ROOT / "inq/install/share/pseudopod"),
         "INQ_SOURCE": str(ROOT / "inq-study")}

CONFIGS = [
    dict(tag="cs_s2p0_v2p0", sigma=2.0, v=2.0, launch_z=-24.0, n_steps=3623,
         save_every=13, wf_every=13, wp_out="cs_s2p0_v2p0", cl_out="cs_s2p0_v2p0_cl"),
    dict(tag="s4_e200", sigma=4.0, v=3.8340, launch_z=-24.5, n_steps=1890,
         save_every=6, wf_every=18, wp_out="twin_s4_e200", cl_out="twin_s4_e200_cl"),
]


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


def resume_flag(rundir: Path, out: str) -> str:
    return "1" if (rundir / "results" / out / "rt_state.txt").exists() else "0"


def wp_env(cfg, gpu):
    return {**os.environ, **SHARE, "PATH": f"{ROOT}/shared/bin:" + os.environ.get("PATH", ""),
            "CUDA_VISIBLE_DEVICES": str(gpu),
            "LJ_SIGMA": f"{cfg['sigma']}", "LJ_K0": f"{cfg['v']}", "LJ_LAUNCH_Z": f"{cfg['launch_z']}",
            "LJ_SPACING": f"{DX}", "LJ_DT": f"{DT}", "LJ_N_STEPS": str(cfg["n_steps"]),
            "LJ_CAP_ETA": "-1.0", "LJ_CAP_L": f"{FACE}",
            "LJ_SAVE_EVERY": str(cfg["save_every"]), "LJ_WF_EVERY": str(cfg["wf_every"]),
            "LJ_CKPT_EVERY": "0", "LJ_RESUME": resume_flag(WPDIR, cfg["wp_out"]),
            "LJ_OUT": cfg["wp_out"], "LJ_GS_DIR": str(GS)}


def cl_env(cfg, gpu):
    return {**os.environ, **SHARE, "PATH": f"{ROOT}/shared/bin:" + os.environ.get("PATH", ""),
            "CUDA_VISIBLE_DEVICES": str(gpu),
            "LJ_LX": "35", "LJ_LY": "35", "LJ_LZ": "85", "LJ_HALF": "12.5", "LJ_N": "100",
            "LJ_EDGE_W": "1.0", "LJ_PERIODICITY": "2", "LJ_SPACING": f"{DX}",
            "LJ_SIGMA": f"{cfg['sigma']}", "LJ_MASS": "1.0", "LJ_DELTA": "0.1", "LJ_DT": f"{DT}",
            "LJ_CONST_V": "0", "LJ_LAUNCH_Z": f"{cfg['launch_z']}", "LJ_K0": f"{cfg['v']}",
            "LJ_CAP_ETA": "-1.0", "LJ_CAP_INNER_BOHR": "30.0", "LJ_CAP_WIDTH_BOHR": f"{FACE}",
            "LJ_N_STEPS": str(cfg["n_steps"]), "LJ_SAVE_EVERY": str(cfg["save_every"]),
            "LJ_GS_DIR": str(GS), "LJ_OUT": cfg["cl_out"], "LJ_RESUME": resume_flag(CLDIR, cfg["cl_out"])}


def completed(rundir: Path, out: str) -> bool:
    s = rundir / "results" / out / "run_summary.txt"
    return s.exists() and "run_completed = true" in s.read_text()


def gpu_free(g):
    return not subprocess.run(["fuser", f"/dev/nvidia{g}"], capture_output=True, text=True).stdout.strip()


def launch(rundir: Path, env, label):
    lf = open(rundir / f"run_{label}.log", "a")
    log(f"launch {label} on GPU {env['CUDA_VISIBLE_DEVICES']} (resume={env['LJ_RESUME']})")
    return subprocess.Popen([str(rundir / "run")], cwd=str(rundir), env=env,
                            stdout=lf, stderr=subprocess.STDOUT, start_new_session=True)


def _cat(rundir: Path, out: str, stem: str):
    import pandas as pd
    fs = sorted(glob.glob(str(rundir / "results" / out / "raw/observables" / f"{stem}*.csv")))
    if not fs:
        return None
    return (pd.concat([pd.read_csv(f, comment="#") for f in fs])
              .drop_duplicates("step").sort_values("step").reset_index(drop=True))


def extract(cfg):
    lines = []
    try:
        pj = _cat(CLDIR, cfg["cl_out"], "projectile"); ob = _cat(CLDIR, cfg["cl_out"], "observables")
        ke = lambda z: pj.loc[(pj.proj_z - z).abs().idxmin(), "energy_proj_ke"]
        S_A = (ke(-FACE) - ke(FACE)) * HA / LSLAB
        S_B = (ob.energy_total.iloc[-1] - E_GS) * HA / LSLAB
        lines.append(f"CLASSICAL: S_A(KE-loss)={S_A:.3f}  S_B(E-abs)={S_B:.3f} eV/Bohr  "
                     f"v_final={float(pj.proj_vz.iloc[-1]):.3f} (v0={cfg['v']})")
    except Exception as e:
        lines.append(f"classical extract failed: {e}")
    try:
        ob = _cat(WPDIR, cfg["wp_out"], "observables")
        rs = _cat(WPDIR, cfg["wp_out"], "wp_real_space_stats")
        mo = _cat(WPDIR, cfg["wp_out"], "wp_momentum_stats")
        E_raw = float(ob.energy_total.iloc[-1]); normWP = float(rs.norm_check.iloc[-1])
        T1 = float(mo.e_kin_ha.iloc[-1]); E_corr = E_raw - T1 * (1 - normWP)
        lines.append(f"WP: S_deposit_raw={(E_raw-E_GS)*HA/LSLAB:.3f}  "
                     f"S_deposit_corr(naive)={(E_corr-E_GS)*HA/LSLAB:.3f} eV/Bohr  "
                     f"norm_WP_final={normWP:.3f}  sigma_z_final={float(rs.sigma_z2.iloc[-1])**0.5:.1f} Bohr")
    except Exception as e:
        lines.append(f"WP extract failed: {e}")
    return lines


def run_cfg(cfg):
    jobs = []
    if not completed(WPDIR, cfg["wp_out"]): jobs.append((WPDIR, cfg["wp_out"], wp_env, "wp"))
    if not completed(CLDIR, cfg["cl_out"]): jobs.append((CLDIR, cfg["cl_out"], cl_env, "cl"))
    if not jobs:
        log(f"{cfg['tag']}: both halves already complete"); return
    log(f"=== config {cfg['tag']}: running {[l for *_, l in jobs]} ===")
    email(f"[twins] {cfg['tag']} started",
          f"WP+classical matched twin, sigma_WP={cfg['sigma']}, v={cfg['v']}, launch={cfg['launch_z']}, "
          f"N_STEPS={cfg['n_steps']}, dx=0.5, CAP eta=-1/12.5, frames every {cfg['save_every']} steps. "
          f"Running: {[l for *_, l in jobs]}. Dual-GPU, checkpointed.")
    procs = {}; retried = set(); t0 = time.time(); pending = list(jobs)
    while pending or procs:
        free = [g for g in (0, 1) if gpu_free(g) and g not in procs]
        while pending and free:
            rundir, out, env_fn, label = pending.pop(0); g = free.pop(0)
            procs[g] = [launch(rundir, env_fn(cfg, g), f"{cfg['tag']}_{label}"), rundir, out, env_fn, label, -1, time.time()]
            time.sleep(5)
        if not procs:
            time.sleep(POLL_S); continue
        time.sleep(POLL_S)
        if time.time() - t0 > MAX_HOURS * 3600:
            log("MAX_HOURS — leaving (checkpointed)"); email("[twins] MAX_HOURS cap", "Checkpointed; resume auto."); return
        for g, stt in list(procs.items()):
            p, rundir, out, env_fn, label, last, changed = stt
            if p.poll() is None:
                try:
                    txt = (rundir / f"run_{cfg['tag']}_{label}.log").read_text()
                    cur = 0
                    for line in reversed(txt.splitlines()):
                        m = re.search(r"step\s+(\d+)", line)
                        if m: cur = int(m.group(1)); break
                    if cur != last: stt[5] = cur; stt[6] = time.time()
                    elif time.time() - changed > STALL_MIN * 60:
                        log(f"STALL {label} @step {cur} — kill"); os.killpg(os.getpgid(p.pid), signal.SIGKILL)
                except Exception as e:
                    log(f"liveness {label}: {e}")
                continue
            del procs[g]
            if completed(rundir, out):
                log(f"DONE {cfg['tag']}_{label}")
            elif label not in retried:
                retried.add(label); log(f"{label} incomplete — retry once (resume auto)")
                pending.append((rundir, out, env_fn, label))
            else:
                log(f"{cfg['tag']}_{label} FAILED twice")
                email(f"[twins] {cfg['tag']}_{label} FAILED twice", f"See {rundir}/run_{cfg['tag']}_{label}.log")
    lines = extract(cfg)
    log(f"FINAL {cfg['tag']}:\n" + "\n".join(lines))
    email(f"[twins] {cfg['tag']} COMPLETE", f"Matched twin done (sigma_WP={cfg['sigma']}, v={cfg['v']}).\n\n" +
          "\n".join(lines) + f"\n\nFrames saved for GIFs. WP: wp/results/{cfg['wp_out']}/  "
          f"CL: dyn_direct_cap/results/{cfg['cl_out']}/")


if __name__ == "__main__":
    log("=== case-study + s4 matched twins (dual-GPU queue) ===")
    if not (WPDIR / "run").exists():
        email("[twins] ABORT: WP binary missing", f"{WPDIR/'run'} not built")
    else:
        for cfg in CONFIGS:
            run_cfg(cfg)
    log("=== all configs done ===")
