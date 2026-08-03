#!/usr/bin/env python3
"""
Autonomous orchestrator for the classical CAP twins of the WP high-density S(v)
sweep at sigma_WP = 2.0 and 3.0 (density A, r_s=4.18, N=100, 35x35x85 slab).

Purpose: give every wavepacket run in hypotheses/wp_highdensity_sv/sweep_data/
(sigma_WP in {2.0, 3.0} x v in {2.0, 2.5, 3.0, 3.5}) an identical-setup classical
twin. The sigma_WP=0.5 twins already exist (dyn_direct_cap/results/v*_cap); this
covers the two missing widths => 8 runs.

Twin parity (see docs/plans/classical-cap-twins.md): identical cell, periodicity,
N, r_s, slab geometry, launch_z, dt, velocity/k0, n_steps, and CAP parameters as
the WP runs. The ONE tolerated departure is grid spacing (dx=0.5 classical vs 0.4
WP; user decision 2026-08-01 - the classical charge has no momentum aliasing). The
classical Gaussian potential width is sigma_pot = sigma_WP/sqrt2, derived INSIDE
run.cpp from LJ_SIGMA, so we pass LJ_SIGMA = sigma_WP.

Analysis contract: these are CAP slab runs => the reported estimator is
  S_B = (E_total(t_final) - E_GS) / L_slab      (L_slab = 25 Bohr)
the E-absorbed definition the CAP makes well-defined (docs/plans/classical-cap-twins.md).
S_A (gauge-free KE-loss across the slab) is kept as a sanity channel.

Dual-GPU (0 and 1, one run per free GPU, refills as runs finish), idempotent
(skips run_completed), retry-once on failure, never self-blocks (checkpoint-dont-block).
Runs are checkpointed (rt_state.txt + checkpoint dir) => resumable via LJ_RESUME=1.

Launch detached (survives the session; NOT the harness-tracked background):
  cd .../scripts/classical_highdensity_sv
  setsid nohup .../venv/bin/python3 orchestrate_cap_sigma.py \
      > dyn_direct_cap/orchestrate_sigma.log 2>&1 < /dev/null &
"""
import os, sys, math, time, subprocess, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# reuse the proven helpers from the sigma=0.5 dispatchers
from orchestrate_direct import vtag, email, log, ROOT, GS, HA, E_GS, LSLAB, FACE
import pandas as pd

SYS  = f"{ROOT}/ResearchProject/systems/localised_jellium"
DIR  = f"{SYS}/scripts/classical_highdensity_sv/dyn_direct_cap"
BIN  = f"{DIR}/run"
HYP  = f"{SYS}/hypotheses/classical_highdensity_sv/dyn_direct"

SIGMAS = [2.0, 3.0]                 # the two missing WP widths (sigma_WP=0.5 twins already exist)
VELS   = [2.0, 2.5, 3.0, 3.5]      # the WP twin velocity grid (v=4.0/4.5 have no WP twin)
PATH_STEPS = 7245                   # n_steps*v: fixed 290-Bohr path, matches the WP n_steps exactly

def stag(s): return f"s{s:.1f}".replace(".", "p")          # 2.0 -> s2p0
def outname(s, v): return f"{stag(s)}_{vtag(v)}_cap"       # -> s2p0_v2p0_cap
def nsteps(v): return int(math.ceil(PATH_STEPS / v))       # v2.0->3623 2.5->2898 3.0->2415 3.5->2070

def completed(out):
    s = f"{DIR}/results/{out}/run_summary.txt"
    return os.path.exists(s) and "run_completed = true" in open(s).read()

def gpu_free(g):
    return not subprocess.run(["fuser", f"/dev/nvidia{g}"], capture_output=True, text=True).stdout.strip()

def env_for(s, v, gpu):
    out = outname(s, v); ns = nsteps(v)
    e = dict(os.environ)
    e.update({
        "PATH": f"{ROOT}/shared/bin:" + e.get("PATH", ""),
        "INQ_SHARE_PATH": f"{ROOT}/inq/install/share",
        "PSEUDOPOD_SHARE_PATH": f"{ROOT}/inq/install/share/pseudopod",
        "CUDA_VISIBLE_DEVICES": str(gpu),
        "LJ_LX": "35", "LJ_LY": "35", "LJ_LZ": "85", "LJ_HALF": "12.5", "LJ_N": "100",
        "LJ_EDGE_W": "1.0", "LJ_PERIODICITY": "2", "LJ_SPACING": "0.5",
        "LJ_SIGMA": f"{s:.4f}",                              # sigma_WP; run.cpp uses sigma_pot=this/sqrt2
        "LJ_MASS": "1.0", "LJ_DELTA": "0.1", "LJ_DT": "0.04", "LJ_CONST_V": "0",
        "LJ_LAUNCH_Z": "-24.0", "LJ_K0": f"{v:.4f}",
        "LJ_CAP_ETA": "-1.0", "LJ_CAP_INNER_BOHR": "30.0", "LJ_CAP_WIDTH_BOHR": "12.5",
        "LJ_N_STEPS": str(ns), "LJ_SAVE_EVERY": str(max(1, round(ns / 300))),
        "LJ_GS_DIR": GS, "LJ_OUT": out, "LJ_RESUME": "0",
    })
    return e, out, ns

def launch(s, v, gpu):
    e, out, ns = env_for(s, v, gpu)
    lf = open(f"{DIR}/run_{out}.log", "a")
    log(f"launch {out} (sigma_WP={s}, v={v}) on GPU {gpu}: {ns} steps (290-Bohr path)")
    return subprocess.Popen([BIN], cwd=DIR, env=e, stdout=lf, stderr=subprocess.STDOUT)

def _cat(out, stem):
    fs = sorted(glob.glob(f"{DIR}/results/{out}/raw/observables/{stem}*.csv"))
    return (pd.concat([pd.read_csv(f) for f in fs])
              .drop_duplicates("step").sort_values("step").reset_index(drop=True))

def extract_S(s, v):
    out = outname(s, v)
    pj = _cat(out, "projectile"); ob = _cat(out, "observables"); ix = _cat(out, "interactions")
    ke = lambda z: pj.loc[(pj.proj_z - z).abs().idxmin(), "energy_proj_ke"]
    S_A = (ke(-FACE) - ke(FACE)) * HA / LSLAB                # gauge-free KE-loss (sanity channel)
    Eabs = (ob.energy_total.iloc[-1] - E_GS) * HA            # E-absorbed
    return dict(sigma_WP=s, v=v,
                v_final=float(pj.proj_vz.iloc[-1]), z_final=float(pj.proj_z.iloc[-1]),
                n_steps=int(pj.step.iloc[-1]), norm_final=float(ix.norm_slab.iloc[-1]),
                S_A_keloss=S_A, S_B_Eabs=Eabs / LSLAB, E_absorbed_eV=Eabs)

def synthesise():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rows = []
    for s in SIGMAS:
        for v in VELS:
            try: rows.append(extract_S(s, v))
            except Exception as e: log(f"S(sigma={s}, v={v}) skipped: {e}")
    if not rows: return None
    df = pd.DataFrame(rows).sort_values(["sigma_WP", "v"])
    csv = f"{HYP}/S_of_v_cap_sigma.csv"; df.to_csv(csv, index=False)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    colours = {2.0: "tab:orange", 3.0: "tab:red"}
    for s in SIGMAS:
        d = df[df.sigma_WP == s]
        if d.empty: continue
        ax.plot(d.v, d.S_B_Eabs, "^-", color=colours[s], label=f"sigma_WP={s}  Def-B (E_abs/25)")
        ax.plot(d.v, d.S_A_keloss, "o--", color=colours[s], alpha=0.45, label=f"sigma_WP={s}  Def-A (KE-loss/25)")
    # sigma_WP=0.5 twins, if present, for context
    s05 = f"{HYP}/S_of_v_cap.csv"
    if os.path.exists(s05):
        try:
            o = pd.read_csv(s05)
            ax.plot(o.v, o.S_B_Eabs, "^-", color="tab:green", alpha=0.6, label="sigma_WP=0.5  Def-B")
        except Exception: pass
    ax.set_xlabel("launch velocity v (a.u.)"); ax.set_ylabel("S (eV/Bohr)")
    ax.set_title("Classical CAP twins S(v), r_s=4.18 (sigma_WP sweep)"); ax.legend(fontsize=8)
    png = f"{HYP}/S_of_v_cap_sigma.png"; fig.tight_layout(); fig.savefig(png, dpi=140); plt.close(fig)
    log(f"synthesis -> {csv} + {png}")
    return df, csv, png

def main():
    log("=== classical CAP twins, sigma_WP sweep (dual-GPU) ===")
    if not os.path.exists(BIN):
        email("[cap-twins] ABORT: binary missing", f"{BIN} not built"); return
    queue = [(s, v) for s in SIGMAS for v in VELS if not completed(outname(s, v))]
    email("[cap-twins] classical CAP twins started",
          "Hypothesis: every WP high-density run needs an identical-setup classical twin so that "
          "S_B (E-absorbed) can be compared quantum-vs-classical at r_s=4.18.\n\n"
          "What is running: 8 classical direct-potential + two-sided CAP (eta=-1, +/-[30,42.5], ETRS, "
          "inq-study) runs, sigma_WP in {2.0, 3.0} x v in {2.0, 2.5, 3.0, 3.5}, dx=0.5, dx0p5 GS. "
          "n_steps = ceil(7245/v) (matches the WP twins). Dual-GPU, checkpointed, idempotent.\n\n"
          "Queue: " + ", ".join(outname(s, v) + f"={nsteps(v)}steps" for s, v in queue))
    free = [g for g in (0, 1) if gpu_free(g)] or [0]
    log(f"queue={[outname(s,v) for s,v in queue]}  free GPUs={free}")
    running = {}; retried = set()
    while queue or running:
        while free and queue:
            s, v = queue.pop(0); g = free.pop(0)
            running[g] = (launch(s, v, g), s, v); time.sleep(3)
        time.sleep(15)
        for g, (p, s, v) in list(running.items()):
            if p.poll() is None: continue
            out = outname(s, v); ok = completed(out); del running[g]; free.append(g)
            if ok:
                try:
                    S = extract_S(s, v)
                    log(f"DONE {out}: S_B={S['S_B_Eabs']:.3f} S_A={S['S_A_keloss']:.3f} "
                        f"norm_final={S['norm_final']:.4f} E_abs={S['E_absorbed_eV']:.2f} eV")
                    email(f"[cap-twins] {out} done: S_B={S['S_B_Eabs']:.2f} eV/Bohr",
                          f"Classical CAP twin complete (sigma_WP={s}, v={v}).\n"
                          f"S_B (E-absorbed / 25 Bohr) = {S['S_B_Eabs']:.3f} eV/Bohr   "
                          f"[E_absorbed = {S['E_absorbed_eV']:.2f} eV]\n"
                          f"S_A (KE-loss / 25 Bohr, sanity) = {S['S_A_keloss']:.3f} eV/Bohr\n"
                          f"norm_slab_final = {S['norm_final']:.4f}, z_final = {S['z_final']:.0f} Bohr, "
                          f"v_final = {S['v_final']:.3f}.\nFull ledger in results/{out}/raw/observables/.")
                except Exception as e:
                    log(f"post {out} failed: {e}")
            else:
                log(f"FAILED {out} (no run_completed)")
                if (s, v) not in retried:
                    retried.add((s, v))
                    log(f"retry-once {out}: wiping and re-queueing")
                    subprocess.run(["rm", "-rf", f"{DIR}/results/{out}"]); queue.append((s, v))
                else:
                    email(f"[cap-twins] {out} FAILED twice", f"{out} did not complete after one retry. "
                          f"See results/{out} and run_{out}.log. Other runs continue.")
    res = synthesise()
    if res:
        df, csv, png = res
        email("[cap-twins] SWEEP COMPLETE - classical CAP twins S(v)",
              "Hypothesis: quantum-vs-classical S_B comparison needs matched classical twins at "
              "r_s=4.18 for sigma_WP=2.0 and 3.0.\n\n"
              "What was done: 8 classical direct-potential+CAP runs completed (see table).\n\n"
              "What the plot shows: S(v) for sigma_WP=2.0 and 3.0, the E-absorbed definition (Def-B, "
              "solid) that the CAP makes well-defined, with the KE-loss channel (Def-A, dashed) and "
              "the existing sigma_WP=0.5 twin curve for context.\n\n"
              "Conclusion: these are the classical halves of the twin pairs; ready to overlay against "
              "the WP S_B (which is a LOWER BOUND - the WP is in the ledger and the CAP removes it).\n\n"
              + df.to_string(index=False), attach=[png, csv])
    log("=== done ===")

if __name__ == "__main__":
    main()
