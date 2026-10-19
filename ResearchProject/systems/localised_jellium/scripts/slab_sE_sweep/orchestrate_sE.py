#!/usr/bin/env python3
"""slab_sE_sweep — autonomous S(E) campaign orchestrator (2-GPU work-pool).

Generalises slab_bulk_limit (S vs 1/L_z at ONE energy) to a RANGE of projectile
energies, producing a final S(E) plot with a classical and a wavepacket trace,
each point being the bulk (L_z→∞) limit from a quadratic 1/L_z extrapolation of a
5-slab-length sweep.

Design (user-locked 2026-09-11):
  * density r_s=5.70 (n0=0.001288), identical to slab_bulk_limit
  * energies v0 = 1.5, 2.0, 2.71, 3.5, 4.5 a.u.  (KE ≈ 31, 54, 100, 167, 276 eV)
  * slab lengths L_z = 15, 30, 45, 60, 75 Bohr
  * σ_WP=2.0, launch run-up d=10 Bohr from the near face at every L_z
  * σ_r(t) time-averaged (in-slab + full-run) on the WP half
  * NO VTI files (LJ_SAVE_EVERY=0); RT checkpoints deleted on completion

WP-ABSORPTION HANDLING (fix 2026-09-11): a slow, dispersing WP needs far longer to
be fully CAP-absorbed (norm_wp→0) than the classical projectile needs to traverse.
So the WP half is (a) sized longer (traversal + a 60 a.u. drain margin) and (b)
AUTO-EXTENDED (LJ_RESUME=1, larger LJ_N_STEPS) if norm_wp_final > 1e-3, up to a cap.
The WP checkpoint is preserved until absorption is confirmed (then pruned); the
classical checkpoint is pruned as soon as it reaches its target. This is the
final-timestep-checkpoint "extend, don't restart" rule.

Phases:
  1. GS: build the 5 geometry-only ground states, pool of 2. A length whose GS
     fails to converge is DROPPED (≥3 points still extrapolate).
  2. Production: 5 energies × 5 lengths × {classical,wp} chained onto both GPUs,
     ENERGY BY ENERGY (ascending). A cell is MEASURED when its classical half has
     reached target AND its WP half is finalised (absorbed or extension-capped).
     When all lengths of an energy are measured, that energy's S vs 1/L_z is built
     and (S_cl_bulk, S_wp_bulk) appended to sE_results.csv.
  3. Final: S(E) plot.

RESUMABLE + ADOPTING (scans /proc), disk-guarded, per-phase Gmail.
"""
import csv, importlib, json, math, os, shutil, subprocess, sys, time
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
LJ   = ROOT / "ResearchProject/systems/localised_jellium"
BASE = LJ / "scripts/slab_sE_sweep"
GS_BIN = LJ / "scripts/campaign_autorun/gs/run"
CL_BIN = BASE / "classical/run"
WP_BIN = BASE / "wp/run"
SHARED_GS = LJ / "shared_gs/slab_sE_sweep"
RUNS = BASE / "runs"
HYP  = LJ / "hypotheses/slab_sE_sweep"
SKILL = "/home/raid/skcb2/skcb2/tddft/.claude/skills/slab-stopping-power"

# ---------------- locked design ----------------
ENERGIES = [(4.5, 276), (3.5, 167), (2.71, 100), (2.0, 54), (1.5, 31)]  # (v0 a.u., KE eV) — DESCENDING (user 2026-09-12: clean/fast high-energy first)
SWEEP    = [15, 30, 45]          # user 2026-09-14: DROP L_z=60,75 — they oversubscribe the 24GB A30 (UVM thrash, ~100× slowdown). Small boxes fit; extrapolate bulk from 3 points.
N_OF_LZ  = {lz: round(0.001288 * 50 * 50 * lz) for lz in SWEEP}         # {15:48,...,75:241}
DT       = 0.04
SIZE_FAC = 1.6                    # traversal margin (both halves)
WP_DRAIN_AU = 60.0               # extra a.u. the WP half runs to CAP-drain after traversal
WP_EXT_STEPS = 1250              # extension chunk (50 a.u.) if still not absorbed
WP_EXT_ROUNDS = 4                # max auto-extensions (cap = base + 4·1250 steps); E31 showed 2 too few
NORM_ABS = 1e-3                  # WP absorbed when norm_wp_final <= this
CAP_ETA, CAP_L = -1.0, 10.0
SPACING, SIGMA, EDGE_W, PER = 0.5, 2.0, 1.0, 2
LX = LY = 50.0
MIN_FREE_GB = 50.0
GPUS = [0, 1]
WATCH_WINDOW_S = 1200            # warmup grace before judging a run's throughput
WATCH_MIN_RATE = 0.5             # steps/min; below this over the window => thrashing/hung -> kill+fail
                                 # (healthy small-box runs do 5–25 steps/min; UVM thrash ~0.24)

RESULTS_CSV = BASE / "sweep_results.csv"
SE_CSV      = BASE / "sE_results.csv"
STATUS_JSON = BASE / "sweep_status.json"
LOG         = BASE / "orchestrate_sE.log"


def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def free_gb(path=ROOT):
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize / 1e9


def base_env():
    e = dict(os.environ)
    e.update(INQ_SOURCE=str(ROOT / "inq-study"),
             INQ_SHARE_PATH=str(ROOT / "inq/install/share"),
             PSEUDOPOD_SHARE_PATH=str(ROOT / "inq/install/share/pseudopod"))
    return e


def read_state_d(rt_state, key):
    try:
        for ln in Path(rt_state).read_text().splitlines():
            if ln.startswith(key + "="):
                return float(ln.split("=")[1])
    except FileNotFoundError:
        return None
    return None


def read_last_step(rt_state):
    v = read_state_d(rt_state, "last_step")
    return int(v) if v is not None else -1


def n_steps_classical(v0, lz):
    """Traversal + CAP arrival for the point-like projectile."""
    return int(math.ceil(SIZE_FAC * (lz + 45.0) / v0 / DT))


def n_steps_wp(v0, lz):
    """Traversal + a drain margin so the dispersing WP fully CAP-absorbs."""
    return int(math.ceil((SIZE_FAC * (lz + 45.0) / v0 + WP_DRAIN_AU) / DT))


# ================= ground state (phase 1) =================
def gs_dir(lz):
    return SHARED_GS / f"Lz{lz}" / "checkpoint"


def gs_energy_ha(lz):
    for cand in (SHARED_GS / f"Lz{lz}" / "run_summary.txt",
                 SHARED_GS / f"Lz{lz}" / "gs.log"):
        if cand.exists():
            for ln in cand.read_text().splitlines():
                if "ground_state_energy_ha" in ln:
                    return float(ln.split("=")[1])
                if "GS energy =" in ln:
                    return float(ln.split("GS energy =")[1].split("Ha")[0])
    return None


def gs_ok(lz):
    return gs_dir(lz).exists() and gs_energy_ha(lz) is not None


def gs_env(lz):
    return dict(LJ_LX=str(LX), LJ_LY=str(LY), LJ_LZ=str(lz + 60), LJ_HALF=str(lz / 2.0),
                LJ_N=str(N_OF_LZ[lz]), LJ_EDGE_W=str(EDGE_W), LJ_PERIODICITY=str(PER),
                LJ_SPACING=str(SPACING), LJ_GS_DIR=str(gs_dir(lz)), LJ_TAG=f"Lz{lz}")


def build_gs_jobs():
    jobs = []
    for lz in SWEEP:
        if gs_ok(lz):
            continue
        rundir = SHARED_GS / f"Lz{lz}"
        jobs.append(dict(id=f"gs_Lz{lz}", kind="gs", half="", lz=lz, binary=GS_BIN,
                         rundir=rundir, out="", env=gs_env(lz), target=None,
                         rt_state=rundir / "run_summary.txt"))
    return jobs


# ================= production (phase 2) =================
def slab_env(v0, lz, half):
    ns = n_steps_classical(v0, lz) if half == "classical" else n_steps_wp(v0, lz)
    e = dict(LJ_LX=str(LX), LJ_LY=str(LY), LJ_LZ=str(lz + 60), LJ_HALF=str(lz / 2.0),
             LJ_N=str(N_OF_LZ[lz]), LJ_EDGE_W=str(EDGE_W), LJ_PERIODICITY=str(PER),
             LJ_SPACING=str(SPACING), LJ_SIGMA=str(SIGMA),
             LJ_LAUNCH_Z=str(-lz / 2.0 - 10.0), LJ_K0=str(v0),
             LJ_N_STEPS=str(ns), LJ_DT=str(DT),
             LJ_CAP_ETA=str(CAP_ETA), LJ_CAP_L=str(CAP_L),
             LJ_SAVE_EVERY="0", LJ_CKPT_EVERY="200", LJ_GS_DIR=str(gs_dir(lz)),
             LJ_OUT=half)
    if half == "classical":
        e["LJ_MASS"] = "1.0"; e["LJ_DELTA"] = "0.1"
    return e, ns


def build_prod_jobs(active_lengths):
    jobs = []
    for v0, ev in ENERGIES:
        for lz in active_lengths:
            rundir = RUNS / f"E{ev}_Lz{lz}"
            for half, binary in (("classical", CL_BIN), ("wp", WP_BIN)):
                env, ns = slab_env(v0, lz, half)
                maxt = ns if half == "classical" else ns + WP_EXT_ROUNDS * WP_EXT_STEPS
                jobs.append(dict(id=f"E{ev}_Lz{lz}_{half}", kind="prod", ev=ev, v0=v0,
                                 lz=lz, half=half, binary=binary, rundir=rundir, out=half,
                                 target=ns, maxt=maxt, env=env,
                                 rt_state=rundir / "results" / half / "rt_state.txt"))
    return jobs


# ================= shared pool mechanics =================
def obs_dir_of(job):
    return job["rundir"] / "results" / job["out"] / "raw/observables"


def wp_norm_final(job):
    """Last norm_wp across all interaction segments; None if unreadable."""
    import pandas as pd
    d = obs_dir_of(job)
    parts = sorted(list(d.glob("interactions.csv")) + list(d.glob("interactions.from*.csv")))
    if not parts:
        return None
    try:
        df = pd.concat([pd.read_csv(p) for p in parts]).drop_duplicates("step").sort_values("step")
        return float(df["norm_wp"].iloc[-1]) if "norm_wp" in df else None
    except Exception:
        return None


def reached_target(job):
    if job["kind"] == "gs":
        return gs_ok(job["lz"])
    return read_last_step(job["rt_state"]) >= job["target"]


def checkpoint_dir(job):
    if job["kind"] == "gs":
        return None
    return job["rundir"] / "results" / job["out"] / "checkpoint"


def prune_checkpoint(job):
    ck = checkpoint_dir(job)
    if ck and ck.exists():
        try:
            shutil.rmtree(ck); log(f"    pruned checkpoint {ck}")
        except OSError as ex:
            log(f"    prune failed {ck}: {ex}")


def discover_running(jobs):
    key2id = {(str(j["rundir"]), j["env"].get("LJ_OUT", "")): j["id"] for j in jobs}
    live = {}
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            if not os.readlink(f"/proc/{pid}/exe").endswith("/run"):
                continue
            cwd = os.readlink(f"/proc/{pid}/cwd")
            data = open(f"/proc/{pid}/environ", "rb").read()
        except (OSError, IOError):
            continue
        env = {}
        for kv in data.split(b"\0"):
            if b"=" in kv:
                k, v = kv.split(b"=", 1)
                env[k.decode("utf-8", "ignore")] = v.decode("utf-8", "ignore")
        jid = key2id.get((cwd, env.get("LJ_OUT", "")))
        cvd = env.get("CUDA_VISIBLE_DEVICES", "")
        if jid and cvd in ("0", "1"):
            live[int(cvd)] = jid
    return live


def current_step(job):
    """Latest step from the newest observables*.csv (cheap tail); -1 if none."""
    d = obs_dir_of(job)
    files = list(d.glob("observables.csv")) + list(d.glob("observables.from*.csv"))
    if not files:
        return -1
    newest = max(files, key=lambda p: p.stat().st_mtime)
    try:
        with open(newest, "rb") as f:
            f.seek(0, 2); size = f.tell()
            f.seek(max(0, size - 4096)); tail = f.read().decode("utf-8", "ignore")
        for ln in reversed(tail.splitlines()):
            tok = ln.split(",")[0]
            if tok.isdigit():
                return int(tok)
    except OSError:
        pass
    return -1


def find_run_pid(job):
    """PID of the live /run process for this job (cwd==rundir & LJ_OUT match)."""
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            if not os.readlink(f"/proc/{pid}/exe").endswith("/run"):
                continue
            if os.readlink(f"/proc/{pid}/cwd") != str(job["rundir"]):
                continue
            data = open(f"/proc/{pid}/environ", "rb").read()
        except (OSError, IOError):
            continue
        if f"LJ_OUT={job['out']}".encode() in data.split(b"\0") or \
           any(kv == f"LJ_OUT={job['out']}".encode() for kv in data.split(b"\0")):
            return int(pid)
    return None


def kill_run(job):
    pid = find_run_pid(job)
    if pid:
        for sig in (9,):
            try:
                os.kill(pid, sig)
            except OSError:
                pass
    return pid


def launch(job, gpu):
    e = base_env(); e.update(job["env"]); e["CUDA_VISIBLE_DEVICES"] = str(gpu)
    e["LJ_N_STEPS"] = str(job["target"])          # current (possibly extended) target
    ls = read_last_step(job["rt_state"]) if job["kind"] == "prod" else -1
    if job["kind"] == "prod" and 0 <= ls < job["target"] and checkpoint_dir(job).exists():
        e["LJ_RESUME"] = "1"; how = f"RESUME {ls}->{job['target']}"
    else:
        how = f"fresh->{job['target']}"
    Path(job["rundir"]).mkdir(parents=True, exist_ok=True)
    lf = open(job["rundir"] / ((job["out"] or "gs") + ".log"), "a")
    lf.write(f"\n===== {time.strftime('%F %T')} GPU{gpu} {job['id']} ({how}) =====\n"); lf.flush()
    log(f"  launch {job['id']} -> GPU{gpu} ({how})")
    return subprocess.Popen(["setsid", str(job["binary"])], cwd=str(job["rundir"]),
                            env=e, stdout=lf, stderr=subprocess.STDOUT)


def notify(subject, body, attachments=None):
    try:
        sys.path.insert(0, str(ROOT / "inq-stack/python"))
        from inqview.email import send_run_email
        send_run_email(subject=f"[slab_sE_sweep] {subject}", body=body, attachments=attachments or [])
        log(f"  emailed: {subject}")
    except Exception as ex:
        log(f"  email failed ({subject}): {ex}")


def run_pool(jobs, phase, on_done=None):
    """2-GPU work-pool. GS/classical: done when target reached (prune immediately).
    WP: on reaching target, check norm_wp — if absorbed(≤NORM_ABS) or capped, finalise
    (prune); else EXTEND (grow target, keep checkpoint, relaunch RESUME)."""
    by_id = {j["id"]: j for j in jobs}
    done, failed = set(), set()
    for j in jobs:                       # seed: already-finalised jobs
        if j["kind"] in ("gs", "prod") and reached_target(j):
            if j["kind"] != "prod" or j["half"] == "classical":
                done.add(j["id"]); prune_checkpoint(j)
            else:                        # wp already at target — finalise only if absorbed/capped
                nf = wp_norm_final(j)
                if (nf is not None and nf <= NORM_ABS) or j["target"] >= j["maxt"]:
                    j["absorbed"] = (nf is not None and nf <= NORM_ABS)
                    done.add(j["id"]); prune_checkpoint(j)
                # else: leave pending -> will be extended below
    slots = {g: None for g in GPUS}
    for gpu, jid in discover_running(jobs).items():
        if jid not in done:
            slots[gpu] = (jid, None); log(f"ADOPTED {jid} on GPU{gpu}")
    started = {s[0] for s in slots.values() if s} | done
    prog = {g: None for g in GPUS}       # watchdog: per-GPU (jid, anchor_step, anchor_time)

    while True:
        live_ids = set(discover_running(jobs).values())
        for gpu in GPUS:                 # reap
            if not slots[gpu]:
                continue
            jid, p = slots[gpu]; j = by_id[jid]
            if reached_target(j):
                if j["kind"] != "prod" or j["half"] == "classical":
                    log(f"  {jid} DONE on GPU{gpu}"); done.add(jid); slots[gpu] = None
                    prune_checkpoint(j)
                    if on_done: on_done(jid, done)
                else:                    # WP finished a segment — absorbed?
                    nf = wp_norm_final(j)
                    if (nf is not None and nf <= NORM_ABS):
                        log(f"  {jid} ABSORBED (norm_wp={nf:.2e}) on GPU{gpu}")
                        j["absorbed"] = True; done.add(jid); slots[gpu] = None
                        prune_checkpoint(j)
                        if on_done: on_done(jid, done)
                    elif j["target"] >= j["maxt"]:
                        log(f"  {jid} EXT-CAPPED (norm_wp={nf}) at {j['target']} on GPU{gpu}")
                        j["absorbed"] = False; done.add(jid); slots[gpu] = None
                        prune_checkpoint(j)
                        if on_done: on_done(jid, done)
                    else:                # extend: keep checkpoint, relaunch RESUME
                        j["target"] = min(j["target"] + WP_EXT_STEPS, j["maxt"])
                        log(f"  {jid} not absorbed (norm_wp={nf}) -> EXTEND to {j['target']}")
                        started.discard(jid); slots[gpu] = None
            elif (p is not None and p.poll() is not None) or (p is None and jid not in live_ids):
                rc = p.returncode if p is not None else "adopted-gone"
                log(f"  {jid} EXITED-NOT-DONE GPU{gpu} (rc={rc}, {read_last_step(j['rt_state'])}/{j['target']}) -> failed")
                notify(f"{jid} FAILED", f"[{phase}] exited early. Queue continues.")
                failed.add(jid); slots[gpu] = None
        for gpu in GPUS:                 # WATCHDOG: kill runs whose throughput collapsed
            if not slots[gpu]:
                prog[gpu] = None; continue
            jid, _ = slots[gpu]; j = by_id[jid]
            if j["kind"] != "prod":
                continue
            cs = current_step(j); now = time.time()
            if prog[gpu] is None or prog[gpu][0] != jid:
                prog[gpu] = (jid, cs, now); continue
            _, a_step, a_time = prog[gpu]
            if now - a_time >= WATCH_WINDOW_S:
                rate = (cs - a_step) / ((now - a_time) / 60.0)
                if rate < WATCH_MIN_RATE:
                    log(f"  WATCHDOG: {jid} {rate:.2f} steps/min < {WATCH_MIN_RATE} over "
                        f"{now-a_time:.0f}s (step {cs}) -> KILL+fail (likely GPU oversubscription)")
                    kill_run(j); notify(f"{jid} THRASHING killed",
                        f"[{phase}] {rate:.2f} steps/min at step {cs}. Killed+skipped; queue continues.")
                    failed.add(jid); slots[gpu] = None; prog[gpu] = None
                else:
                    prog[gpu] = (jid, cs, now)
        for gpu in GPUS:                 # fill
            if slots[gpu]:
                continue
            nxt = next((j for j in jobs if j["id"] not in done and j["id"] not in failed
                        and j["id"] not in started), None)
            if nxt is None:
                continue
            if free_gb() < MIN_FREE_GB:
                log(f"DISK GUARD: free {free_gb():.0f} GB < {MIN_FREE_GB}; pausing before {nxt['id']}")
                notify("DISK GUARD hit", f"[{phase}] free={free_gb():.0f} GB before {nxt['id']}. "
                       f"Free space and re-run orchestrate_sE.py to resume.")
                sys.exit(3)
            slots[gpu] = (nxt["id"], launch(nxt, gpu)); started.add(nxt["id"])
        STATUS_JSON.write_text(json.dumps(
            {"updated": time.strftime("%F %T"), "phase": phase, "free_gb": round(free_gb(), 1),
             "done": sorted(done), "failed": sorted(failed),
             "running": {str(g): (s[0] if s else None) for g, s in slots.items()}}, indent=2))
        if all(j["id"] in done or j["id"] in failed for j in jobs) and not any(slots.values()):
            break
        time.sleep(20)
    return done, failed


# ================= measurement + σ_r =================
def sigma_r_averages(wp_obs_dir, half_width):
    import pandas as pd
    parts = sorted(wp_obs_dir.glob("wp_real_space_stats*.csv"),
                   key=lambda p: (".from" in p.name, p.name))
    frames = []
    for p in parts:
        try:
            frames.append(pd.read_csv(p, comment="#"))
        except Exception:
            pass
    if not frames:
        return {}
    df = pd.concat(frames).drop_duplicates("step").sort_values("step")
    sr = (df["sigma_x2"].clip(lower=0) + df["sigma_y2"].clip(lower=0)) ** 0.5
    sz = df["sigma_z2"].clip(lower=0) ** 0.5
    inslab = df["z_mean"].abs() < half_width
    out = {"sigma_r_full": float(sr.mean()), "sigma_z_full": float(sz.mean()),
           "sigma_r_t0": float(sr.iloc[0]), "n_wp_steps": int(len(df))}
    if inslab.any():
        out["sigma_r_inslab"] = float(sr[inslab].mean())
        out["sigma_z_inslab"] = float(sz[inslab].mean())
        out["n_inslab_steps"] = int(inslab.sum())
    else:
        out["sigma_r_inslab"] = None; out["sigma_z_inslab"] = None; out["n_inslab_steps"] = 0
    return out


def measure_ez(ev, v0, lz):
    sys.path.insert(0, SKILL)
    SS = importlib.import_module("slab_stopping"); importlib.reload(SS)
    e_gs = gs_energy_ha(lz)
    row = {"E_eV": ev, "v0": v0, "L_z": lz, "E_GS_ha": e_gs}
    for half in ("classical", "wp"):
        obs = RUNS / f"E{ev}_Lz{lz}" / "results" / half / "raw/observables"
        try:
            m = SS.measure(obs, half, e_gs, float(lz))
            row[f"S_{half}"] = m["S_eV_per_Bohr"]; row[f"E_abs_{half}"] = m["E_absorbed_eV"]
            if "norm_final" in m:
                row["norm_wp_final"] = m["norm_final"]
        except Exception as ex:
            row[f"S_{half}"] = None; row[f"err_{half}"] = str(ex)[:160]
            log(f"  E{ev} Lz{lz} {half}: measure FAILED: {ex}")
    try:
        row.update(sigma_r_averages(RUNS / f"E{ev}_Lz{lz}" / "results/wp/raw/observables", lz / 2.0))
    except Exception as ex:
        log(f"  E{ev} Lz{lz}: sigma_r FAILED: {ex}")
    try:
        rs = RUNS / f"E{ev}_Lz{lz}" / "results/classical/rt_state.txt"
        row["proj_z_final"] = read_state_d(rs, "proj_z"); row["proj_vz_final"] = read_state_d(rs, "proj_vz")
        row["traversed"] = (row["proj_z_final"] is not None and row["proj_z_final"] > lz / 2.0)
    except Exception:
        pass
    return row


RESULT_COLS = ["E_eV", "v0", "L_z", "E_GS_ha", "S_classical", "S_wp",
               "E_abs_classical", "E_abs_wp", "norm_wp_final",
               "sigma_r_inslab", "sigma_r_full", "sigma_z_inslab", "sigma_z_full",
               "sigma_r_t0", "n_inslab_steps", "n_wp_steps",
               "proj_z_final", "proj_vz_final", "traversed", "err_classical", "err_wp"]


def upsert_result(row):
    """Idempotent: replace any existing (E_eV,L_z) row, keep sorted."""
    rows = []
    if RESULTS_CSV.exists():
        rows = [r for r in csv.DictReader(open(RESULTS_CSV))
                if not (int(float(r["E_eV"])) == row["E_eV"] and int(float(r["L_z"])) == row["L_z"])]
    rows.append({k: row.get(k) for k in RESULT_COLS})
    rows.sort(key=lambda r: (float(r["E_eV"]), float(r["L_z"])))
    with open(RESULTS_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=RESULT_COLS); w.writeheader(); w.writerows(rows)


def build_energy_plot(ev):
    r = subprocess.run([str(ROOT / "venv/bin/python3"), str(HYP / "build_S_vs_invLz.py"), str(ev)],
                       capture_output=True, text=True)
    log(f"  build_S_vs_invLz E{ev}: rc={r.returncode}\n{r.stdout[-800:]}{r.stderr[-400:]}")
    return HYP / f"S_vs_invLz_E{ev}.png"


def build_final_sE():
    r = subprocess.run([str(ROOT / "venv/bin/python3"), str(HYP / "build_S_vs_E.py")],
                       capture_output=True, text=True)
    log(f"  build_S_vs_E: rc={r.returncode}\n{r.stdout[-800:]}{r.stderr[-400:]}")
    return HYP / "S_vs_E.png"


# ================= main =================
def main():
    log("=" * 72)
    log(f"slab_sE_sweep START  energies={[e[1] for e in ENERGIES]}eV  L_z={SWEEP}  free={free_gb():.0f} GB")
    notify("campaign START/RESUME",
           f"S(E) slab→bulk. Energies(eV): {[e[1] for e in ENERGIES]}  L_z: {SWEEP}\n"
           f"WP half auto-extends until absorbed (norm_wp≤{NORM_ABS}). Disk free {free_gb():.0f} GB.")

    gs_jobs = build_gs_jobs()
    if gs_jobs:
        log(f"PHASE 1 (GS): {[j['id'] for j in gs_jobs]}")
        run_pool(gs_jobs, phase="GS")
    active = [lz for lz in SWEEP if gs_ok(lz)]
    dropped = [lz for lz in SWEEP if lz not in active]
    if dropped:
        notify("GS drop", f"GS unconverged for L_z={dropped}; proceeding with {active}.")
    if len(active) < 3:
        notify("ABORT", f"Only {len(active)} GS converged ({active})."); log("ABORT <3 GS"); return
    log(f"GS ready: L_z={active}")

    prod_jobs = build_prod_jobs(active)
    energies_done = set()
    if SE_CSV.exists():
        import pandas as pd
        try:
            energies_done = set(int(x) for x in pd.read_csv(SE_CSV)["E_eV"].tolist())
        except Exception:
            pass
    measured = set()
    if RESULTS_CSV.exists():
        import pandas as pd
        try:
            for _, rr in pd.read_csv(RESULTS_CSV).iterrows():
                if pd.notna(rr.get("S_wp")):                 # only fully-measured cells count as done
                    measured.add((int(rr["E_eV"]), int(rr["L_z"])))
        except Exception:
            pass

    def on_prod_done(_jid, done):
        for v0, ev in ENERGIES:
            for lz in active:
                if (ev, lz) in measured:
                    continue
                if f"E{ev}_Lz{lz}_classical" in done and f"E{ev}_Lz{lz}_wp" in done:
                    row = measure_ez(ev, v0, lz); upsert_result(row); measured.add((ev, lz))
                    log(f"  MEASURED E{ev} Lz{lz}: S_cl={row.get('S_classical')} S_wp={row.get('S_wp')} "
                        f"σr_inslab={row.get('sigma_r_inslab')}")
                    notify(f"E{ev} Lz{lz} done",
                           f"S_classical={row.get('S_classical')} eV/Bohr\nS_wp={row.get('S_wp')} eV/Bohr\n"
                           f"norm_wp_final={row.get('norm_wp_final')}\n"
                           f"σ_r(in-slab)={row.get('sigma_r_inslab')}  σ_r(full)={row.get('sigma_r_full')} Bohr\n"
                           f"traversed={row.get('traversed')}")
            if ev not in energies_done and all((ev, lz) in measured for lz in active):
                png = build_energy_plot(ev); energies_done.add(ev)
                notify(f"E{ev} S vs 1/L_z ready", f"All {len(active)} lengths measured for {ev} eV.",
                       attachments=[str(png)] if png.exists() else [])

    run_pool(prod_jobs, phase="production", on_done=on_prod_done)
    on_prod_done(None, {j["id"] for j in prod_jobs if reached_target(j)})

    png = build_final_sE()
    log(f"CAMPAIGN COMPLETE. -> {RESULTS_CSV}, {SE_CSV}")
    notify("campaign COMPLETE", f"S(E) built. Per-energy plots + final S(E) in {HYP}.",
           attachments=[str(png)] if png.exists() else [])


if __name__ == "__main__":
    main()
