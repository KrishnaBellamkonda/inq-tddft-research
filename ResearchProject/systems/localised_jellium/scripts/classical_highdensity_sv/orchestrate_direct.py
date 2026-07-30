#!/usr/bin/env python3
"""
Autonomous orchestrator for the DIRECT-potential classical S(v) sweep (localised jellium
slab, r_s=4.18, N=100, 35x35x85, CAP-free). Re-runs the old charge-based sweep velocities
with the direct erf/r potential fix (dyn_direct/run binary: direct perturbation + direct
force + direct ledger). v=4.5 already done; this completes v=2.0,2.5,3.0,3.5,4.0.

Uses BOTH GPUs when free (one run per free GPU, refills as runs finish). Idempotent
(skips run_completed). Full energy decomposition + projectile z(t)/v(t) stored by the
binary (observables.csv, projectile.csv, interactions.csv) + density frames + checkpoint.
Per-run: compute S (gauge-free KE-loss across the slab), best-effort run notebook, email.
Final: S(v) synthesis plot + summary CSV + email. checkpoint-dont-block: never self-block.

Launch detached:
  cd .../classical_highdensity_sv
  setsid nohup .../venv/bin/python3 orchestrate_direct.py > dyn_direct/orchestrate.log 2>&1 < /dev/null &
"""
import os, sys, math, time, subprocess, glob

ROOT = "/local/data/public/skcb2/tddft"
SYS  = f"{ROOT}/ResearchProject/systems/localised_jellium"
DIR  = f"{SYS}/scripts/classical_highdensity_sv/dyn_direct"
BIN  = f"{DIR}/run"
GS   = f"{SYS}/shared_gs/slab_n100_L35x35x85_dx0p5_per2"
HYP  = f"{SYS}/hypotheses/classical_highdensity_sv/dyn_direct"
VENV = f"{ROOT}/venv/bin/python3"
EMAIL_TO = "chiddukanna@gmail.com"
HA = 27.211386; E_GS = 207.18322156141; LSLAB = 25.0; FACE = 12.5
VELOCITIES = [2.0, 2.5, 3.0, 3.5, 4.0]          # v4.5 already complete
sys.path.insert(0, f"{ROOT}/inq-stack/python")

def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def email(subject, body, attach=None):
    try:
        from inqview.email import send_run_email
        send_run_email(subject, body, attachments=[a for a in (attach or []) if a and os.path.exists(a)], to=EMAIL_TO)
        log(f"emailed: {subject}")
    except Exception as e:
        log(f"email FAILED ({subject}): {e}")

def n_steps_for(v):
    return int(math.ceil(1.4 * 69.0 / (0.5 * v * 0.04)))

def vtag(v): return f"v{v:.1f}".replace(".", "p")

def completed(out):
    s = f"{DIR}/results/{out}/run_summary.txt"
    return os.path.exists(s) and "run_completed = true" in open(s).read()

def gpu_free(g):
    r = subprocess.run(["fuser", f"/dev/nvidia{g}"], capture_output=True, text=True)
    return not r.stdout.strip()

def env_for(v, gpu):
    out = f"{vtag(v)}_direct"; ns = n_steps_for(v)
    e = dict(os.environ)
    e.update({
        "PATH": f"{ROOT}/shared/bin:" + e.get("PATH", ""),
        "INQ_SHARE_PATH": f"{ROOT}/inq/install/share",
        "PSEUDOPOD_SHARE_PATH": f"{ROOT}/inq/install/share/pseudopod",
        "CUDA_VISIBLE_DEVICES": str(gpu),
        "LJ_LX": "35", "LJ_LY": "35", "LJ_LZ": "85", "LJ_HALF": "12.5", "LJ_N": "100",
        "LJ_EDGE_W": "1.0", "LJ_PERIODICITY": "2", "LJ_SPACING": "0.5", "LJ_SIGMA": "0.5",
        "LJ_MASS": "1.0", "LJ_DELTA": "0.1", "LJ_DT": "0.04", "LJ_CONST_V": "0",
        "LJ_LAUNCH_Z": "-24.0", "LJ_K0": f"{v:.4f}",
        "LJ_N_STEPS": str(ns), "LJ_SAVE_EVERY": str(max(1, round(ns / 300))),
        "LJ_GS_DIR": GS, "LJ_OUT": out, "LJ_RESUME": "0",
    })
    return e, out, ns

def launch(v, gpu):
    e, out, ns = env_for(v, gpu)
    lf = open(f"{DIR}/run_{out}.log", "a")
    log(f"launch {out} on GPU {gpu}  n_steps={ns}")
    p = subprocess.Popen([BIN], cwd=DIR, env=e, stdout=lf, stderr=subprocess.STDOUT)
    return p

def extract_S(v):
    """gauge-free KE-loss across the equal-potential slab window (-FACE..+FACE)."""
    import pandas as pd
    out = f"{vtag(v)}_direct"
    pj = pd.read_csv(f"{DIR}/results/{out}/raw/observables/projectile.csv")
    ke = lambda z: pj.loc[(pj.proj_z - z).abs().idxmin(), "energy_proj_ke"]
    dep = (ke(-FACE) - ke(FACE)) * HA
    vmean = pj[(pj.proj_z > -FACE) & (pj.proj_z < FACE)].proj_vz.mean()
    return dict(v=v, v_final=float(pj.proj_vz.iloc[-1]), v_mean_slab=float(vmean),
                deposit_eV=float(dep), S_eV_per_Bohr=float(dep / LSLAB))

def build_notebook(v):
    try:
        r = subprocess.run([VENV, f"{HYP}/build_run_notebook.py", vtag(v)],
                           cwd=HYP, capture_output=True, text=True, timeout=1800,
                           env={**os.environ, "PYTHONPATH": f"{ROOT}/inq-stack/python"})
        log(f"notebook {vtag(v)}: {r.stdout.strip().splitlines()[-1] if r.stdout.strip() else 'built'}")
    except Exception as ex:
        log(f"notebook {vtag(v)} FAILED: {ex}")

def synthesise():
    import pandas as pd, matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = []
    for v in [2.0, 2.5, 3.0, 3.5, 4.0, 4.5]:
        try: rows.append(extract_S(v))
        except Exception as e: log(f"S({v}) skipped: {e}")
    if not rows: return None
    df = pd.DataFrame(rows).sort_values("v")
    os.makedirs(HYP, exist_ok=True)
    csv = f"{HYP}/S_of_v_direct.csv"; df.to_csv(csv, index=False)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(df.v, df.S_eV_per_Bohr, "o-", label="direct (corrected)")
    # old charge sweep for comparison (if the summary csv exists)
    old = f"{SYS}/hypotheses/classical_highdensity_sv/sv_sweep/S_summary.csv"
    if os.path.exists(old):
        try:
            o = pd.read_csv(old); vc = "v" if "v" in o else o.columns[0]
            sc = [c for c in o.columns if "S" in c][0]
            ax.plot(o[vc], o[sc], "s--", color="tab:red", alpha=0.7, label="old charge (sheet-inflated)")
        except Exception: pass
    ax.set_xlabel("launch velocity v (a.u.)"); ax.set_ylabel("S (eV/Bohr)")
    ax.set_title("Direct-potential S(v), localised jellium slab r_s=4.18"); ax.legend()
    png = f"{HYP}/S_of_v_direct.png"; fig.tight_layout(); fig.savefig(png, dpi=140); plt.close(fig)
    log(f"synthesis: {csv} + {png}")
    return df, csv, png

def main():
    log("=== DIRECT S(v) sweep orchestrator ===")
    if not os.path.exists(BIN): email("[direct-sv] ABORT: binary missing", f"{BIN} not found"); return
    queue = [v for v in VELOCITIES if not completed(f"{vtag(v)}_direct")]
    done_already = [v for v in VELOCITIES if v not in queue]
    log(f"queue={queue}  already-complete={done_already}")
    free = [g for g in (0, 1) if gpu_free(g)]
    log(f"free GPUs at start: {free}")
    if not free: free = [0]                                  # never self-block; use GPU 0
    email("[direct-sv] sweep started",
          f"Direct-potential classical S(v) sweep launched.\nQueue: {queue}\nGPUs: {free}\n"
          f"Each run stores full energy decomposition (observables.csv), projectile z(t)/v(t) "
          f"(projectile.csv), pairwise ledger (interactions.csv), density frames + checkpoint.\n"
          f"Idempotent; per-run email on completion; S(v) synthesis at the end.")
    running = {}                                             # gpu -> (Popen, v)
    while queue or running:
        while free and queue:
            v = queue.pop(0); g = free.pop(0); running[g] = (launch(v, g), v)
            time.sleep(3)
        time.sleep(15)
        for g, (p, v) in list(running.items()):
            if p.poll() is None: continue
            rc = p.returncode; out = f"{vtag(v)}_direct"; ok = completed(out)
            del running[g]; free.append(g)
            if ok:
                try:
                    S = extract_S(v)
                    log(f"DONE {out} rc={rc}: S={S['S_eV_per_Bohr']:.3f} eV/Bohr (v_mean={S['v_mean_slab']:.2f})")
                    build_notebook(v)
                    email(f"[direct-sv] {out} done: S={S['S_eV_per_Bohr']:.2f} eV/Bohr",
                          f"v={v}: deposit {S['deposit_eV']:.2f} eV over 25 Bohr -> S={S['S_eV_per_Bohr']:.3f} eV/Bohr\n"
                          f"v_final={S['v_final']:.3f}, v_mean_in_slab={S['v_mean_slab']:.3f}\n"
                          f"Full ledger + z(t)/v(t) in results/{out}/raw/observables/.\n"
                          f"Notebook: hypotheses/.../dyn_direct/run_{vtag(v)}_direct.ipynb",
                          attach=[f"{HYP}/run_{vtag(v)}_direct.ipynb"])
                except Exception as e:
                    log(f"post-processing {out} failed: {e}")
                    email(f"[direct-sv] {out} finished (post-proc issue)", f"rc={rc}; error {e}")
            else:
                log(f"FAILED {out} rc={rc} (no run_completed) -> requeue once")
                if v not in getattr(main, "_retried", set()):
                    main._retried = getattr(main, "_retried", set()) | {v}; queue.append(v)
                    email(f"[direct-sv] {out} FAILED (rc={rc}), retrying once", f"see run_{out}.log")
                else:
                    email(f"[direct-sv] {out} FAILED twice (rc={rc}) - skipping", f"see run_{out}.log")
    res = synthesise()
    if res:
        df, csv, png = res
        body = "Direct-potential S(v) sweep COMPLETE (r_s=4.18 slab).\n\n" + df.to_string(index=False) + \
               "\n\nS = gauge-free KE-loss across the equal-potential slab (-12.5..+12.5)/25 Bohr.\n" \
               "The old charge sweep was sheet-inflated ~20-35%; this is the corrected curve."
        email("[direct-sv] SWEEP COMPLETE — corrected S(v)", body, attach=[png, csv])
    log("=== orchestrator done ===")

if __name__ == "__main__":
    main()
