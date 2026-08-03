#!/usr/bin/env python3
"""
Extend every direct-potential run by +25% steps (LJ_RESUME=1), then rebuild all run
notebooks and the S(v) synthesis with BOTH stopping definitions:
  Def-A  S_A = [KE(-12.5) - KE(+12.5)] / L_slab          (KE-loss across the slab; gauge-free)
  Def-B  S_B = [E_total(t_final) - E_GS] / L_slab         (E_absorbed at the final step; user spec)
Dual-GPU (both free), idempotent (skips runs already at the extended target). Segment CSVs
(observables.from<N>.csv ...) are concatenated by the notebook builders.

Launch detached:
  cd .../classical_highdensity_sv
  setsid nohup .../venv/bin/python3 extend_and_update.py > dyn_direct/extend.log 2>&1 < /dev/null &
"""
import os, sys, math, time, subprocess, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orchestrate_direct import (n_steps_for, vtag, email, log, DIR, HYP, GS, BIN, ROOT, HA, E_GS, LSLAB, FACE)

VENV = f"{ROOT}/venv/bin/python3"
ALL = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
FACTOR = 1.25                                             # +25% steps
V4P5_BUILDER = f"{HYP}/build_notebook_v4p5_direct.py"     # full deep-dive for v4.5
GEN_BUILDER  = f"{HYP}/build_run_notebook.py"             # compact deep-dive for the rest
import pandas as pd

def target(v): return int(math.ceil(FACTOR * n_steps_for(v)))

def read_last_step(out):
    p = f"{DIR}/results/{out}/rt_state.txt"
    if not os.path.exists(p): return -1
    for ln in open(p):
        if ln.startswith("last_step="):
            try: return int(float(ln.split("=")[1]))
            except: return -1
    return -1

def gpu_free(g):
    return not subprocess.run(["fuser", f"/dev/nvidia{g}"], capture_output=True, text=True).stdout.strip()

def env_ext(v, gpu):
    out = f"{vtag(v)}_direct"; ns = target(v)
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
        "LJ_N_STEPS": str(ns), "LJ_SAVE_EVERY": str(max(1, round(n_steps_for(v) / 300))),
        "LJ_GS_DIR": GS, "LJ_OUT": out, "LJ_RESUME": "1",
    })
    return e, out, ns

def launch(v, gpu):
    e, out, ns = env_ext(v, gpu)
    lf = open(f"{DIR}/run_{out}.log", "a")
    log(f"EXTEND {out} on GPU {gpu}: resume -> {ns} steps (+25%)")
    return subprocess.Popen([BIN], cwd=DIR, env=e, stdout=lf, stderr=subprocess.STDOUT)

# ---- both-definition S extraction from concatenated segments ----
def _cat(out, stem):
    fs = sorted(glob.glob(f"{DIR}/results/{out}/raw/observables/{stem}*.csv"))
    return pd.concat([pd.read_csv(f) for f in fs]).drop_duplicates("step").sort_values("step").reset_index(drop=True)

def both_S(v):
    out = f"{vtag(v)}_direct"
    pj = _cat(out, "projectile"); ob = _cat(out, "observables")
    ke = lambda z: pj.loc[(pj.proj_z - z).abs().idxmin(), "energy_proj_ke"]
    S_A = (ke(-FACE) - ke(FACE)) * HA / LSLAB
    Eabs = (ob.energy_total.iloc[-1] - E_GS) * HA
    S_B = Eabs / LSLAB
    return dict(v=v, v_final=float(pj.proj_vz.iloc[-1]), z_final=float(pj.proj_z.iloc[-1]),
                n_steps=int(pj.step.iloc[-1]),
                S_A_keloss=S_A, S_B_Eabs=S_B, E_absorbed_eV=Eabs)

def rebuild(v):
    try:
        script = V4P5_BUILDER if abs(v - 4.5) < 1e-6 else GEN_BUILDER
        args = [VENV, script] + ([] if abs(v - 4.5) < 1e-6 else [vtag(v)])
        r = subprocess.run(args, cwd=HYP, capture_output=True, text=True, timeout=2400,
                           env={**os.environ, "PYTHONPATH": f"{ROOT}/inq-stack/python"})
        tail = (r.stdout.strip().splitlines() or ["(no output)"])[-1]
        log(f"rebuild {vtag(v)}: {tail}")
    except Exception as ex:
        log(f"rebuild {vtag(v)} FAILED: {ex}")

def synthesise():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rows = [both_S(v) for v in ALL]
    df = pd.DataFrame(rows).sort_values("v")
    csv = f"{HYP}/S_of_v_direct.csv"; df.to_csv(csv, index=False)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(df.v, df.S_A_keloss, "o-", color="tab:blue", label="Def-A: KE-loss across slab /25")
    ax.plot(df.v, df.S_B_Eabs, "^-", color="tab:green", label="Def-B: [E_total(t_final)-E_GS]/25")
    old = f"{DIR}/../../hypotheses/classical_highdensity_sv/sv_sweep/S_summary.csv"
    if os.path.exists(old):
        try:
            o = pd.read_csv(old); vc = "v" if "v" in o else o.columns[0]; sc = [c for c in o.columns if "S" in c][0]
            ax.plot(o[vc], o[sc], "s--", color="tab:red", alpha=0.6, label="old charge (sheet-inflated)")
        except Exception: pass
    ax.set_xlabel("launch velocity v (a.u.)"); ax.set_ylabel("S (eV/Bohr)")
    ax.set_title("Direct-potential S(v), r_s=4.18 slab (extended +25%) — two definitions"); ax.legend()
    png = f"{HYP}/S_of_v_direct.png"; fig.tight_layout(); fig.savefig(png, dpi=140); plt.close(fig)
    log(f"synthesis -> {csv} + {png}")
    return df, csv, png

def main():
    log("=== extend (+25%) + rebuild (both defs) ===")
    email("[direct-sv] extension started",
          f"Extending all 6 direct runs by +25% steps (LJ_RESUME=1), then rebuilding notebooks + S(v) "
          f"with BOTH definitions (Def-A KE-loss/25 ; Def-B [E_total(t_final)-E_GS]/25). Targets: " +
          ", ".join(f"{vtag(v)}->{target(v)}" for v in ALL))
    queue = [v for v in ALL if read_last_step(f"{vtag(v)}_direct") < target(v)]
    log(f"to-extend: {queue}   already-at-target: {[v for v in ALL if v not in queue]}")
    free = [g for g in (0, 1) if gpu_free(g)] or [0]
    log(f"free GPUs: {free}")
    running = {}
    while queue or running:
        while free and queue:
            v = queue.pop(0); g = free.pop(0); running[g] = (launch(v, g), v); time.sleep(3)
        time.sleep(15)
        for g, (p, v) in list(running.items()):
            if p.poll() is None: continue
            out = f"{vtag(v)}_direct"; ok = read_last_step(out) >= target(v)
            del running[g]; free.append(g)
            log(f"{out} extension {'OK' if ok else 'INCOMPLETE'} (last_step={read_last_step(out)}/{target(v)})")
            if not ok and v not in getattr(main, "_retried", set()):
                main._retried = getattr(main, "_retried", set()) | {v}; queue.append(v)
    log("all extensions done; rebuilding notebooks")
    for v in ALL: rebuild(v)
    df, csv, png = synthesise()
    email("[direct-sv] EXTENDED sweep + both definitions COMPLETE",
          "All direct runs extended +25%; notebooks + S(v) rebuilt with two definitions.\n\n" +
          df.to_string(index=False) +
          "\n\nDef-A = KE-loss across slab /25 (gauge-free). Def-B = [E_total(t_final)-E_GS]/25 "
          "(user spec; note E_total retains a 1/r e_ps offset so Def-B > Def-A and is an upper bound).",
          attach=[png, csv])
    log("=== done ===")

if __name__ == "__main__":
    main()
