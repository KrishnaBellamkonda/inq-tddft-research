#!/usr/bin/env python3
"""
Autonomous orchestrator for WIDE classical CAP runs at sigma_WP = 17 and 20
(density A, r_s=4.18, N=100, 35x35x85 slab) - the classical-point-of-view check.

Motivation (user, 2026-08-02): the WP effective (time-averaged) real-space width
<sigma_r> lands at ~17-22 Bohr (hypotheses/wp_highdensity_sv/momentum_and_effective_width.ipynb).
This runs CLASSICAL projectiles at comparable widths (sigma_WP = 17, 20) to test
whether the classical S(v) at those widths tracks / brackets the WP S(v) - i.e.
whether the quantum trend "makes sense from a classical point of view" once the
classical projectile is as wide as the dispersed packet.

Identical setup to orchestrate_cap_sigma.py (the sigma_WP=2,3 twins): same cell,
periodicity, N, r_s, slab geometry, launch_z, dt, velocity/k0, n_steps, CAP, dx=0.5
GS. The classical Gaussian potential width sigma_pot = sigma_WP/sqrt2 is derived
INSIDE run.cpp from LJ_SIGMA, so we pass LJ_SIGMA = sigma_WP.

  sigma_WP=17 -> sigma_pot=12.02   ;   sigma_WP=20 -> sigma_pot=14.14

CAVEAT (surfaced in the email + report): at these widths sigma_pot exceeds the
transverse half-box (17.5 Bohr), so the projectile potential is broader than the
cell - exactly the "wider than the box" regime the dispersed WP is ALSO in, so the
comparison is fair but qualitative, not a converged infinite-medium number.

Estimator (CAP slab): S_B = (E_total(t_final) - E_GS) / L_slab  (L_slab = 25 Bohr).
S_A (KE-loss across the slab) kept as a gauge-free sanity channel.

Dual-GPU (0 and 1), idempotent, retry-once, checkpointed (checkpoint-dont-block).

Launch detached (survives the session):
  cd .../scripts/classical_highdensity_sv
  setsid nohup .../venv/bin/python3 orchestrate_cap_sigma_wide.py \
      > dyn_direct_cap/orchestrate_sigma_wide.log 2>&1 < /dev/null &
"""
import os, sys, math, time, subprocess, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orchestrate_direct import vtag, email, log, ROOT, GS, HA, E_GS, LSLAB, FACE
# reuse the identical launch/extract machinery from the sigma=2,3 orchestrator
from orchestrate_cap_sigma import (DIR, BIN, HYP, VELS, PATH_STEPS,
                                   stag, outname, nsteps, completed, gpu_free,
                                   env_for, launch, extract_S)
import pandas as pd

SIGMAS = [17.0, 20.0]                # wide widths ~ WP effective <sigma_r>

def synthesise():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rows = []
    for s in SIGMAS:
        for v in VELS:
            try: rows.append(extract_S(s, v))
            except Exception as e: log(f"S(sigma={s}, v={v}) skipped: {e}")
    if not rows: return None
    df = pd.DataFrame(rows).sort_values(["sigma_WP", "v"])
    csv = f"{HYP}/S_of_v_cap_sigma_wide.csv"; df.to_csv(csv, index=False)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    colours = {17.0: "tab:purple", 20.0: "tab:brown"}
    for s in SIGMAS:
        d = df[df.sigma_WP == s]
        if d.empty: continue
        ax.plot(d.v, d.S_B_Eabs, "^-", color=colours[s], label=f"sigma_WP={s:.0f}  Def-B (E_abs/25)")
        ax.plot(d.v, d.S_A_keloss, "o--", color=colours[s], alpha=0.45, label=f"sigma_WP={s:.0f}  Def-A")
    ax.set_xlabel("launch velocity v (a.u.)"); ax.set_ylabel("S (eV/Bohr)")
    ax.set_title("Wide classical CAP S(v), r_s=4.18 (sigma_WP=17,20)"); ax.legend(fontsize=8)
    png = f"{HYP}/S_of_v_cap_sigma_wide.png"; fig.tight_layout(); fig.savefig(png, dpi=140); plt.close(fig)
    log(f"synthesis -> {csv} + {png}")
    return df, csv, png

def main():
    log("=== WIDE classical CAP runs, sigma_WP = 17, 20 (dual-GPU) ===")
    if not os.path.exists(BIN):
        email("[cap-wide] ABORT: binary missing", f"{BIN} not built"); return
    queue = [(s, v) for s in SIGMAS for v in VELS if not completed(outname(s, v))]
    email("[cap-wide] wide classical CAP runs started",
          "Motivation: the WP effective <sigma_r> is ~17-22 Bohr; run classical "
          "projectiles at comparable widths (sigma_WP=17,20; sigma_pot=12.0,14.1) to "
          "check whether classical S(v) there tracks the WP curves.\n\n"
          "What is running: 8 classical direct-potential + two-sided CAP (eta=-1, ETRS, "
          "inq-study) runs, sigma_WP in {17,20} x v in {2.0,2.5,3.0,3.5}, dx=0.5. "
          "n_steps = ceil(7245/v). Dual-GPU, checkpointed, idempotent.\n\n"
          "CAVEAT: sigma_pot exceeds the transverse half-box (17.5 Bohr), so the "
          "potential is broader than the cell - the same 'wider than box' regime the "
          "dispersed WP is in. Comparison is fair but qualitative.\n\n"
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
                    email(f"[cap-wide] {out} done: S_B={S['S_B_Eabs']:.2f} eV/Bohr",
                          f"Wide classical CAP run complete (sigma_WP={s:.0f}, v={v}).\n"
                          f"S_B (E-absorbed / 25 Bohr) = {S['S_B_Eabs']:.3f} eV/Bohr   "
                          f"[E_absorbed = {S['E_absorbed_eV']:.2f} eV]\n"
                          f"S_A (KE-loss / 25 Bohr, sanity) = {S['S_A_keloss']:.3f} eV/Bohr\n"
                          f"norm_slab_final = {S['norm_final']:.4f}, z_final = {S['z_final']:.0f} Bohr, "
                          f"v_final = {S['v_final']:.3f}.")
                except Exception as e:
                    log(f"post {out} failed: {e}")
            else:
                log(f"FAILED {out} (no run_completed)")
                if (s, v) not in retried:
                    retried.add((s, v))
                    log(f"retry-once {out}: wiping and re-queueing")
                    subprocess.run(["rm", "-rf", f"{DIR}/results/{out}"]); queue.append((s, v))
                else:
                    email(f"[cap-wide] {out} FAILED twice", f"{out} did not complete after one retry.")
    res = synthesise()
    if res:
        df, csv, png = res
        email("[cap-wide] SWEEP COMPLETE - wide classical CAP S(v)",
              "Motivation: classical-point-of-view check at WP-comparable widths.\n\n"
              "What was done: 8 wide classical direct-potential+CAP runs (sigma_WP=17,20).\n\n"
              "What the plot shows: S(v) for sigma_WP=17 and 20, E-absorbed (Def-B, solid) "
              "and KE-loss (Def-A, dashed). Overlay these against the WP S_B and the narrow "
              "classical twins to see whether classical stopping at the dispersed-packet width "
              "reproduces the WP values.\n\n"
              "CAVEAT: sigma_pot > transverse half-box; qualitative comparison.\n\n"
              + df.to_string(index=False), attach=[png, csv])
    log("=== done ===")

if __name__ == "__main__":
    main()
