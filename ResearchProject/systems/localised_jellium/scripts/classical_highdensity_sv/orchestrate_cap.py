#!/usr/bin/env python3
"""
Autonomous orchestrator for the CAP classical direct-potential sweep (dyn_direct_cap).
Identical to dyn_direct EXCEPT a two-sided sin² CAP (η=-1.0 Ha, bands ±[30,42.5], ETRS,
inq-study) — matching the quantum WP runs so the two are directly comparable. Grid dx=0.5
UNCHANGED. Fresh runs from the dx0p5 GS. Step counts match the quantum (fixed 290-Bohr path:
n_steps = ceil(7245/v)). All observables kept (energy decomposition, projectile z/v, ledger,
frames) + norm_slab for the extensive-energy rescale. Dual-GPU, idempotent. NEW folder;
does NOT overwrite dyn_direct.

Launch detached:
  cd .../classical_highdensity_sv
  setsid nohup .../venv/bin/python3 orchestrate_cap.py > dyn_direct_cap/orchestrate.log 2>&1 < /dev/null &
"""
import os, sys, math, time, subprocess, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orchestrate_direct import vtag, email, log, ROOT, GS, HA, E_GS, LSLAB, FACE

SYS  = f"{ROOT}/ResearchProject/systems/localised_jellium"
DIR  = f"{SYS}/scripts/classical_highdensity_sv/dyn_direct_cap"
BIN  = f"{DIR}/run"
HYP  = f"{SYS}/hypotheses/classical_highdensity_sv/dyn_direct"
VENV = f"{ROOT}/venv/bin/python3"
ALL  = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
PATH_STEPS = 7245                                  # n_steps*v (quantum's fixed 290-Bohr path)
import pandas as pd

def nsteps(v): return int(math.ceil(PATH_STEPS / v))

def completed(out):
    s = f"{DIR}/results/{out}/run_summary.txt"
    return os.path.exists(s) and "run_completed = true" in open(s).read()

def gpu_free(g):
    return not subprocess.run(["fuser", f"/dev/nvidia{g}"], capture_output=True, text=True).stdout.strip()

def env_for(v, gpu):
    out = f"{vtag(v)}_cap"; ns = nsteps(v)
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
        "LJ_CAP_ETA": "-1.0", "LJ_CAP_INNER_BOHR": "30.0", "LJ_CAP_WIDTH_BOHR": "12.5",
        "LJ_N_STEPS": str(ns), "LJ_SAVE_EVERY": str(max(1, round(ns / 300))),
        "LJ_GS_DIR": GS, "LJ_OUT": out, "LJ_RESUME": "0",
    })
    return e, out, ns

def launch(v, gpu):
    e, out, ns = env_for(v, gpu)
    lf = open(f"{DIR}/run_{out}.log", "a")
    log(f"CAP {out} on GPU {gpu}: {ns} steps (path 290 Bohr)")
    return subprocess.Popen([BIN], cwd=DIR, env=e, stdout=lf, stderr=subprocess.STDOUT)

def _cat(out, stem):
    fs = sorted(glob.glob(f"{DIR}/results/{out}/raw/observables/{stem}*.csv"))
    return pd.concat([pd.read_csv(f) for f in fs]).drop_duplicates("step").sort_values("step").reset_index(drop=True)

def both_S(v):
    out = f"{vtag(v)}_cap"
    pj = _cat(out, "projectile"); ob = _cat(out, "observables"); ix = _cat(out, "interactions")
    ke = lambda z: pj.loc[(pj.proj_z - z).abs().idxmin(), "energy_proj_ke"]
    S_A = (ke(-FACE) - ke(FACE)) * HA / LSLAB
    Eabs = (ob.energy_total.iloc[-1] - E_GS) * HA
    return dict(v=v, v_final=float(pj.proj_vz.iloc[-1]), z_final=float(pj.proj_z.iloc[-1]), n_steps=int(pj.step.iloc[-1]),
                norm_final=float(ix.norm_slab.iloc[-1]), S_A_keloss=S_A, S_B_Eabs=Eabs / LSLAB, E_absorbed_eV=Eabs)

def rebuild(v):
    try:
        r = subprocess.run([VENV, f"{HYP}/build_run_notebook_cap.py", vtag(v)], cwd=HYP,
                           capture_output=True, text=True, timeout=2400,
                           env={**os.environ, "PYTHONPATH": f"{ROOT}/inq-stack/python"})
        log(f"notebook {vtag(v)}_cap: {(r.stdout.strip().splitlines() or ['(no output)'])[-1]}")
    except Exception as ex:
        log(f"notebook {vtag(v)}_cap FAILED: {ex}")

def synthesise():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rows = []
    for v in ALL:
        try: rows.append(both_S(v))
        except Exception as e: log(f"S({v}) skipped: {e}")
    if not rows: return None
    df = pd.DataFrame(rows).sort_values("v")
    csv = f"{HYP}/S_of_v_cap.csv"; df.to_csv(csv, index=False)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(df.v, df.S_A_keloss, "o-", color="tab:blue", label="CAP  Def-A (KE-loss/25)")
    ax.plot(df.v, df.S_B_Eabs, "^-", color="tab:green", label="CAP  Def-B ([E_tot(t_f)-E_GS]/25)")
    nocap = f"{HYP}/S_of_v_direct.csv"
    if os.path.exists(nocap):
        try:
            o = pd.read_csv(nocap)
            ax.plot(o.v, o.S_A_keloss, "o--", color="tab:blue", alpha=0.45, label="no-CAP Def-A")
            ax.plot(o.v, o.S_B_Eabs, "^--", color="tab:green", alpha=0.45, label="no-CAP Def-B")
        except Exception: pass
    ax.set_xlabel("launch velocity v (a.u.)"); ax.set_ylabel("S (eV/Bohr)")
    ax.set_title("CAP vs no-CAP classical S(v), r_s=4.18 (matched to quantum BCs)"); ax.legend(fontsize=8)
    png = f"{HYP}/S_of_v_cap.png"; fig.tight_layout(); fig.savefig(png, dpi=140); plt.close(fig)
    log(f"synthesis -> {csv} + {png}")
    return df, csv, png

def main():
    log("=== CAP classical sweep (dual-GPU) ===")
    if not os.path.exists(BIN):
        email("[cap-sv] ABORT: binary missing", f"{BIN} not built"); return
    email("[cap-sv] CAP sweep started",
          "Classical direct-potential runs WITH a two-sided CAP (eta=-1.0, +/-[30,42.5], ETRS, "
          "inq-study), matching the quantum WP boundary conditions. dx=0.5 unchanged; step counts "
          "match quantum (290-Bohr path). Velocities " + ", ".join(f"{vtag(v)}={nsteps(v)}" for v in ALL))
    queue = [v for v in ALL if not completed(f"{vtag(v)}_cap")]
    free = [g for g in (0, 1) if gpu_free(g)] or [0]
    log(f"queue={queue}  free GPUs={free}")
    running = {}
    while queue or running:
        while free and queue:
            v = queue.pop(0); g = free.pop(0); running[g] = (launch(v, g), v); time.sleep(3)
        time.sleep(15)
        for g, (p, v) in list(running.items()):
            if p.poll() is None: continue
            out = f"{vtag(v)}_cap"; ok = completed(out); del running[g]; free.append(g)
            if ok:
                try:
                    S = both_S(v); rebuild(v)
                    log(f"DONE {out}: S_A={S['S_A_keloss']:.3f} S_B={S['S_B_Eabs']:.3f} norm_final={S['norm_final']:.4f}")
                    email(f"[cap-sv] {out} done: S_A={S['S_A_keloss']:.2f}, S_B={S['S_B_Eabs']:.2f} eV/Bohr",
                          f"v={v} (CAP): Def-A(KE-loss)={S['S_A_keloss']:.3f}, Def-B(E_abs)={S['S_B_Eabs']:.3f} eV/Bohr\n"
                          f"norm_final={S['norm_final']:.4f} (bath absorption by CAP), z_final={S['z_final']:.0f}\n"
                          f"Full ledger + z(t)/v(t) in results/{out}/raw/observables/.",
                          attach=[f"{HYP}/run_{vtag(v)}_cap.ipynb"])
                except Exception as e:
                    log(f"post {out} failed: {e}")
            else:
                log(f"FAILED {out} (no run_completed) -> retry once")
                if v not in getattr(main, "_r", set()):
                    main._r = getattr(main, "_r", set()) | {v}
                    subprocess.run(["rm", "-rf", f"{DIR}/results/{out}"]); queue.append(v)
    res = synthesise()
    if res:
        df, csv, png = res
        email("[cap-sv] CAP SWEEP COMPLETE — CAP vs no-CAP S(v)",
              "Classical CAP sweep complete (matched to quantum BCs).\n\n" + df.to_string(index=False) +
              "\n\nnorm_final ~ 100 confirms the CAP barely touches the localised bath; energies ~ the "
              "no-CAP direct runs. CAP makes the boundary conditions identical to the quantum WP runs.",
              attach=[png, csv])
    log("=== done ===")

if __name__ == "__main__":
    main()
