#!/usr/bin/env python3
"""
GPU-0-only sequential finisher for the direct-potential S(v) sweep. Fallback after
GPU 1 turned out to be intermittently contended (external containerized job time-slices
it: 2 s/step then 50 s/step). v=2.0 done; v=3.0 running unmanaged on GPU 0. This waits
for GPU 0 to free, then runs the remaining incomplete velocities sequentially on GPU 0,
post-processes each (S + notebook + email), builds any missing notebook, and synthesises.

Launch detached:
  cd .../classical_highdensity_sv
  setsid nohup .../venv/bin/python3 finish_gpu0.py > dyn_direct/finish.log 2>&1 < /dev/null &
"""
import os, sys, time, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orchestrate_direct import (n_steps_for, vtag, completed, launch, extract_S,
                                build_notebook, synthesise, email, log, DIR, HYP)

GPU = 0
ALL = [2.0, 2.5, 3.0, 3.5, 4.0]

def gpu0_busy():
    return bool(subprocess.run(["fuser", "/dev/nvidia0"], capture_output=True, text=True).stdout.strip())

def post(v, rc=0):
    out = f"{vtag(v)}_direct"
    try:
        S = extract_S(v)
        log(f"{out}: S={S['S_eV_per_Bohr']:.3f} eV/Bohr (v_mean={S['v_mean_slab']:.2f})")
        build_notebook(v)
        email(f"[direct-sv] {out} done: S={S['S_eV_per_Bohr']:.2f} eV/Bohr",
              f"v={v}: deposit {S['deposit_eV']:.2f} eV / 25 Bohr -> S={S['S_eV_per_Bohr']:.3f} eV/Bohr\n"
              f"v_final={S['v_final']:.3f}, v_mean_in_slab={S['v_mean_slab']:.3f}\n"
              f"Full ledger + z(t)/v(t) in results/{out}/raw/observables/.",
              attach=[f"{HYP}/run_{vtag(v)}_direct.ipynb"])
    except Exception as e:
        log(f"post-proc {out} failed: {e}")

def run_seq(v):
    out = f"{vtag(v)}_direct"
    if completed(out):
        log(f"{out} already complete; ensuring notebook");
        if not os.path.exists(f"{HYP}/run_{vtag(v)}_direct.ipynb"): build_notebook(v)
        return
    # clear any partial (crashed/killed) output so we start clean (no resume-prune issue)
    subprocess.run(["rm", "-rf", f"{DIR}/results/{out}"])
    p = launch(v, GPU)
    while p.poll() is None:
        time.sleep(15)
    if completed(out):
        post(v, p.returncode)
    else:
        log(f"{out} FAILED rc={p.returncode}; one retry")
        subprocess.run(["rm", "-rf", f"{DIR}/results/{out}"])
        p = launch(v, GPU)
        while p.poll() is None: time.sleep(15)
        if completed(out): post(v, p.returncode)
        else: email(f"[direct-sv] {out} FAILED twice", f"rc={p.returncode}; see run_{out}.log")

def main():
    log("=== GPU0-only finisher ===")
    email("[direct-sv] GPU0-only fallback",
          "GPU 1 is intermittently contended (external job: 2 s/step -> 50 s/step), so the "
          "sweep is finishing on GPU 0 only, sequentially. v=2.0 done; v=3.0 running; "
          "then v=2.5, 3.5, 4.0. ETA ~4 h. Per-run emails + final S(v) as before.")
    # 1. wait for the in-flight v3.0 (running unmanaged on GPU 0) to finish and free GPU 0
    log("waiting for GPU 0 to free (in-flight v3.0)...")
    waited = 0
    while gpu0_busy() and not completed("v3p0_direct"):
        time.sleep(30); waited += 30
        if waited % 300 == 0: log(f"  still waiting ({waited//60} min); v3p0 completed={completed('v3p0_direct')}")
    if completed("v3p0_direct"):
        log("v3.0 complete"); post(3.0)
    # 2. sequential remaining on GPU 0 (skip completed)
    for v in [2.5, 3.5, 4.0]:
        # guard: don't collide if GPU0 still busy for some reason
        while gpu0_busy(): time.sleep(20)
        run_seq(v)
    # 3. make sure v2.0 / v3.0 notebooks exist, then synthesise across all six
    for v in [2.0, 3.0]:
        if completed(f"{vtag(v)}_direct") and not os.path.exists(f"{HYP}/run_{vtag(v)}_direct.ipynb"):
            build_notebook(v)
    res = synthesise()
    if res:
        df, csv, png = res
        email("[direct-sv] SWEEP COMPLETE — corrected S(v)",
              "Direct-potential S(v) sweep COMPLETE (r_s=4.18 slab, GPU0-only).\n\n" +
              df.to_string(index=False) +
              "\n\nS = gauge-free KE-loss across the equal-potential slab (-12.5..+12.5)/25 Bohr.\n"
              "Old charge sweep was sheet-inflated ~20-35%; this is the corrected curve.",
              attach=[png, csv])
    log("=== finisher done ===")

if __name__ == "__main__":
    main()
