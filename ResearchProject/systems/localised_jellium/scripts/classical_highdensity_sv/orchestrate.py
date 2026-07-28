#!/usr/bin/env python3
"""
Autonomous orchestrator — high-density classical S(v) benchmark
(campaign classical-highdensity-sv). Runs the whole remaining campaign on ONE GPU:

  pilot(v=2) --[plateau CORRECTNESS gate]--> 6-velocity Ehrenfest sweep
      --> per-run S + plateau plot + density GIF + run-notebook
      --> synthesis S(v) curve + component ledger

Rules honoured:
- Python orchestrator, idempotent resume (skip runs with run_completed=true),
  per-phase try/except (one failure never kills the chain), per-phase Gmail.
- checkpoint-don't-block: long wall-clock only WARNs; only correctness stops the
  sweep (no pilot plateau, NaN/non-finite energy).
- CAP-free, periodicity(2), mass-1 Ehrenfest Gaussian charge (sigma_WP=0.5).

Launch DETACHED so it survives the session:
  cd scripts/classical_highdensity_sv
  setsid nohup env CUDA_VISIBLE_DEVICES=0 \
    /local/data/public/skcb2/tddft/venv/bin/python3 orchestrate.py > orchestrate.log 2>&1 &
"""
import os, sys, csv, json, math, subprocess, traceback
from pathlib import Path

ROOT = "/local/data/public/skcb2/tddft"
SYS  = f"{ROOT}/ResearchProject/systems/localised_jellium"
BIN  = f"{SYS}/scripts/classical_highdensity_sv/dyn/run"
GS   = f"{SYS}/shared_gs/slab_n100_L35x35x85_dx0p5_per2"
RUNDIR = f"{SYS}/scripts/classical_highdensity_sv/dyn"          # binary writes results/ here
HYP  = f"{SYS}/hypotheses/classical_highdensity_sv/sv_sweep"
HA = 27.2114
LSLAB = 25.0
FARFACE = 42.5
EMAIL_TO = "chiddukanna@gmail.com"

VENV = f"{ROOT}/venv/bin/python3"
os.makedirs(HYP, exist_ok=True)

# Projectile-free slab GS energy (the clean E_absorbed reference — NOT E_total(0),
# which includes the projectile's t=0 contribution).
def _read_egs():
    summ = f"{SYS}/scripts/classical_highdensity_sv/gs/results/run_summary.txt"
    for line in open(summ):
        if "ground_state_energy_ha" in line:
            return float(line.split("=")[1])
    raise RuntimeError("could not read E_GS")
E_GS = _read_egs()   # 207.183 Ha

def log(m):
    print(f"[orchestrate] {m}", flush=True)

def email(subject, body, attachments=None):
    try:
        sys.path.insert(0, f"{ROOT}/inq-stack/python")
        from inqview.email import send_run_email
        send_run_email(subject, body, attachments=[a for a in (attachments or []) if a and os.path.exists(a)], to=EMAIL_TO)
        log(f"emailed: {subject}")
    except Exception as e:
        log(f"email FAILED ({subject}): {e}")

# ---- run one Ehrenfest velocity -------------------------------------------
LAUNCH_Z = -24.0    # standoff from the slab face (-12.5) to limit the t=0 turn-on transient
def n_steps_for(v):
    # launch_z=-24 -> exit ~+45; Ehrenfest decelerates, assume mean v >= 0.5*v0
    # over the transit; 1.4x margin. Resumable if it under-runs.
    dist = 69.0
    return int(math.ceil(1.4 * dist / (0.5 * v * 0.04)))

def run_one(v, resume=False):
    name = f"v{v:.1f}".replace(".", "p")
    out  = f"results/{name}"
    outdir = f"{RUNDIR}/{out}"
    summ = f"{outdir}/run_summary.txt"
    if (not resume) and os.path.exists(summ) and "run_completed = true" in open(summ).read():
        log(f"{name}: already complete, skip"); return outdir
    nsteps = n_steps_for(v) * (2 if resume else 1)   # resume => EXTEND to a larger target
    save_every = max(1, round(n_steps_for(v) / 300))
    env = dict(os.environ)
    env.update({
        "PATH": f"{ROOT}/shared/bin:" + env.get("PATH", ""),
        "INQ_SHARE_PATH": f"{ROOT}/inq/install/share",
        "PSEUDOPOD_SHARE_PATH": f"{ROOT}/inq/install/share/pseudopod",
        "LJ_LX":"35","LJ_LY":"35","LJ_LZ":"85","LJ_HALF":"12.5","LJ_N":"100",
        "LJ_EDGE_W":"1.0","LJ_PERIODICITY":"2","LJ_SPACING":"0.5","LJ_SIGMA":"0.5",
        "LJ_MASS":"1.0","LJ_DELTA":"0.1","LJ_DT":"0.04","LJ_CONST_V":"0",
        "LJ_LAUNCH_Z":f"{LAUNCH_Z:.1f}","LJ_K0":f"{v:.4f}",   # k0 = m*v = v (mass 1)
        "LJ_N_STEPS":str(nsteps),"LJ_SAVE_EVERY":str(save_every),
        "LJ_GS_DIR":GS,"LJ_OUT":name,
        "LJ_RESUME":"1" if resume else "0",
    })
    log(f"{name}: launch v={v} n_steps={nsteps} save_every={save_every} resume={resume}")
    with open(f"{RUNDIR}/run_{name}.log","a") as lf:
        subprocess.run([BIN], cwd=RUNDIR, env=env, stdout=lf, stderr=subprocess.STDOUT, check=True)
    return outdir

# ---- load + analyse --------------------------------------------------------
def _read_csv(path):
    import pandas as pd
    return pd.read_csv(path)

def _concat_segments(pattern):
    import pandas as pd, glob
    files = sorted(glob.glob(pattern))                 # observables.csv + observables.fromN.csv
    parts = [pd.read_csv(f) for f in files]
    df = pd.concat(parts, ignore_index=True) if parts else pd.read_csv(pattern.replace("*",""))
    stepcol = "step" if "step" in df.columns else df.columns[0]
    return df.drop_duplicates(subset=stepcol).sort_values(stepcol).reset_index(drop=True)

def load_run(outdir):
    obs  = _concat_segments(f"{outdir}/raw/observables/observables*.csv")
    proj = _concat_segments(f"{outdir}/raw/observables/projectile*.csv")
    etot_col = [c for c in obs.columns if "total" in c.lower()][0]
    return obs, proj, etot_col

def transit_plateau(outdir):
    import numpy as np
    obs, proj, etot = load_run(outdir)
    E = obs[etot].to_numpy(); t = obs["time_au"].to_numpy()
    z = proj["proj_z"].to_numpy(); vz = proj["proj_vz"].to_numpy()
    ke = proj["energy_proj_ke"].to_numpy()
    finite = np.all(np.isfinite(E))
    z_final = float(z[-1]); transited = z_final > FARFACE
    # exit index: projectile center fully past the far face (+2 Bohr of Gaussian tail)
    out_mask = z > (FARFACE + 2.0)
    result = dict(z_final=z_final, v0=float(vz[0]), v_final=float(vz[-1]),
                  transited=bool(transited), finite=bool(finite),
                  n_obs=len(E), n_proj=len(z), E_GS_ha=E_GS)
    if transited and out_mask.any():
        i0 = int(np.argmax(out_mask))                 # first fully-exited frame
        plateau = E[i0:]
        # Definition-2 (HEADLINE): slab excitation vs the projectile-FREE GS.
        deposit = float(np.mean(plateau) - E_GS) * HA
        flat = float(np.std(plateau) * HA) if len(plateau) > 3 else float("nan")
        # cross-check: projectile KE lost from t=0 to exit (energy conservation, CAP-free).
        ke_loss = float(ke[0] - ke[i0]) * HA
        result.update(exited=True, S=deposit / LSLAB, E_absorbed_eV=deposit,
                      S_keloss=ke_loss / LSLAB, ke_loss_eV=ke_loss,
                      plateau_flatness_eV=flat, t_exit=float(t[i0]),
                      plateau_ok=(len(plateau) >= 4 and abs(flat) < max(1.0, 0.10*abs(deposit))))
    else:
        result.update(exited=False, S=float("nan"), E_absorbed_eV=float("nan"),
                      plateau_ok=False)
    return result

def analyse_run(outdir, v):
    """S + plateau plot + density GIF + REPORT.md. Best-effort; returns result dict."""
    import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    r = transit_plateau(outdir)
    name = Path(outdir).name
    ad = f"{HYP}/{name}"; os.makedirs(ad, exist_ok=True)
    obs, proj, etot = load_run(outdir)
    E = obs[etot].to_numpy(); t = obs["time_au"].to_numpy(); z = proj["proj_z"].to_numpy()
    # plateau plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].plot(t, (E - E_GS) * HA); ax[0].axhline(0, ls=":", c="k", label="E_GS (projectile-free)")
    if r.get("exited"):
        ax[0].axvline(r["t_exit"], ls="--", c="C3",
                      label=f"exit; S={r['S']:.2f} (KE-loss {r.get('S_keloss',float('nan')):.2f}) eV/Bohr")
    ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("E_electronic − E_GS (eV)"); ax[0].set_title(f"{name}: deposit + plateau"); ax[0].legend(fontsize=8)
    ax[1].plot(t, z); ax[1].axhline(FARFACE, ls="--", c="k"); ax[1].axhspan(-12.5, 12.5, color="C7", alpha=0.15)
    ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("proj_z (Bohr)"); ax[1].set_title("trajectory")
    fig.tight_layout(); fig.savefig(f"{ad}/plateau.png", dpi=110); plt.close(fig)
    # density GIF (total n(x,z,t)) — best-effort
    try:
        make_density_gif(outdir, f"{ad}/density_evolution.gif")
    except Exception as e:
        log(f"{name}: density GIF failed: {e}")
    # REPORT
    with open(f"{ad}/REPORT.md", "w") as f:
        f.write(f"# {name} (v={v})\n\n")
        for k, val in r.items(): f.write(f"- {k}: {val}\n")
    json.dump(r, open(f"{ad}/result.json", "w"), indent=2)
    return r

def make_density_gif(outdir, gif_path):
    import numpy as np, glob, re, matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt; from matplotlib import animation
    sys.path.insert(0, f"{ROOT}/inq-stack/python"); from inqview import load_vti
    files = sorted(glob.glob(f"{outdir}/frames/total/density_t*.vti"),
                   key=lambda p: int(re.search(r"density_t(\d+)", p).group(1)))
    if not files:
        log("no density frames"); return
    files = files[:: max(1, len(files)//60)]          # ~60 frames
    f0 = load_vti(files[0]); n0 = f0.xz_slice(0.0); X=[f0.x[0],f0.x[-1],f0.z[0],f0.z[-1]]
    slices = [load_vti(p).xz_slice(0.0) for p in files]
    dmax = max(s.max() for s in slices)
    fig, ax = plt.subplots(figsize=(4.2,6))
    im = ax.imshow(slices[0], origin="lower", aspect="auto", extent=X, vmin=0, vmax=dmax)
    ax.axhline(FARFACE, ls="--", c="cyan", lw=1); ax.axhline(-12.5, ls=":", c="w", lw=0.6); ax.axhline(12.5, ls=":", c="w", lw=0.6)
    ax.set_xlabel("x (Bohr)"); ax.set_ylabel("z (Bohr)"); ttl = ax.set_title("")
    fig.colorbar(im, ax=ax, label="n (a0^-3)")
    def upd(i): im.set_data(slices[i]); ttl.set_text(f"frame {i}"); return im, ttl
    anim = animation.FuncAnimation(fig, upd, frames=len(slices), interval=120)
    anim.save(gif_path, writer=animation.PillowWriter(fps=8)); plt.close(fig)

# ---- synthesis -------------------------------------------------------------
def synthesise(results):
    import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    good = [(v, r) for v, r in results if r and r.get("exited") and np.isfinite(r.get("S", float("nan")))]
    with open(f"{HYP}/S_summary.csv", "w") as f:
        w = csv.writer(f); w.writerow(["v","v_final","S_eV_per_Bohr","E_absorbed_eV","plateau_flatness_eV","transited"])
        for v, r in results:
            if r: w.writerow([v, r.get("v_final"), r.get("S"), r.get("E_absorbed_eV"), r.get("plateau_flatness_eV"), r.get("transited")])
    if good:
        vs = [v for v,_ in good]; Ss = [r["S"] for _,r in good]
        fig, ax = plt.subplots(figsize=(7,4.5))
        ax.plot(vs, Ss, "o-", ms=6, label="classical (E_absorbed/L_slab)")
        ax.set_xlabel("v (a.u.)"); ax.set_ylabel("S (eV/Bohr)")
        ax.set_title("High-density classical S(v) — r_s=4.18 slab, CAP-free"); ax.legend()
        fig.tight_layout(); fig.savefig(f"{HYP}/S_of_v.png", dpi=120); plt.close(fig)
    return good

# ---- main ------------------------------------------------------------------
def main():
    log("=== autonomous classical S(v) benchmark START ===")
    if not os.path.exists(BIN):
        log(f"FATAL: binary missing {BIN}"); email("[classical-highdensity-sv] FATAL", f"binary missing {BIN}"); return
    results = []
    # ---- PILOT + transit-floor find (correctness gate) ----
    floor = None
    for v in (2.0, 2.5, 3.0):
        try:
            outdir = run_one(v)
            r = transit_plateau(outdir)
            if not r["finite"]:
                email("[classical-highdensity-sv] ABORT: non-finite energy", f"v={v}: {r}"); return
            r_full = analyse_run(outdir, v); results.append((v, r_full))
            log(f"pilot v={v}: transited={r['transited']} exited={r.get('exited')} S={r.get('S')}")
            if r.get("exited") and r.get("plateau_ok"):
                floor = v; break
            elif r["v_final"] > 0.3 * r["v0"]:
                # still moving, just ran out of steps -> EXTEND once, then re-check
                outdir = run_one(v, resume=True); r = transit_plateau(outdir)
                if r.get("exited") and r.get("plateau_ok"):
                    results[-1] = (v, analyse_run(outdir, v)); floor = v; break
                log(f"pilot v={v}: still no plateau after extend -> try higher v")
            else:
                log(f"pilot v={v}: projectile STOPPED inside (v_final={r['v_final']:.2f}) -> below floor")
        except Exception:
            log(f"pilot v={v} FAILED:\n{traceback.format_exc()}")
            email(f"[classical-highdensity-sv] pilot v={v} exception", traceback.format_exc())
    if floor is None:
        email("[classical-highdensity-sv] STOP: no transit floor in {2,2.5,3}",
              "No pilot velocity transited+plateaued; a mass-1 electron may stop inside "
              "the r_s=4.18 slab even at v=3. Human decision needed (raise v, thin slab, "
              "or heavier projectile).")
        log("STOP: no transit floor"); return
    email(f"[classical-highdensity-sv] PILOT PASS v_min={floor}",
          f"Pilot plateau gate PASSED at v={floor}. S={results[-1][1].get('S'):.2f} eV/Bohr. "
          f"Launching the sweep {floor}..+5.", attachments=[f"{HYP}/v{floor:.1f}".replace('.','p')+"/plateau.png"])
    # ---- SWEEP: floor + 5 up, step 0.5, capped at 4.5 ----
    grid = [round(floor + 0.5*k, 2) for k in range(6)]
    grid = [v for v in grid if v <= 4.5]
    done_v = {v for v, _ in results}
    for v in grid:
        if v in done_v: continue
        try:
            outdir = run_one(v); r = transit_plateau(outdir)
            if not r["finite"]:
                email(f"[classical-highdensity-sv] v={v} non-finite", str(r)); continue
            if r["transited"] and not r.get("exited"):
                outdir = run_one(v, resume=True)                 # extend once
            r_full = analyse_run(outdir, v); results.append((v, r_full))
            email(f"[classical-highdensity-sv] v={v} done  S={r_full.get('S')}",
                  json.dumps(r_full, indent=2),
                  attachments=[f"{HYP}/v{v:.1f}".replace('.','p')+"/plateau.png",
                               f"{HYP}/v{v:.1f}".replace('.','p')+"/density_evolution.gif"])
        except Exception:
            log(f"v={v} FAILED:\n{traceback.format_exc()}")
            email(f"[classical-highdensity-sv] v={v} exception", traceback.format_exc())
    # ---- SYNTHESIS ----
    try:
        good = synthesise(results)
        body = "S(v) classical benchmark (r_s=4.18, CAP-free, periodicity 2):\n" + \
               "\n".join(f"  v={v}: S={r.get('S'):.3f} eV/Bohr (v_final={r.get('v_final'):.2f})" for v, r in results if r)
        email("[classical-highdensity-sv] SWEEP COMPLETE — S(v) curve",
              body, attachments=[f"{HYP}/S_of_v.png", f"{HYP}/S_summary.csv"])
        log("=== DONE ===\n" + body)
    except Exception:
        log(f"synthesis FAILED:\n{traceback.format_exc()}")
        email("[classical-highdensity-sv] synthesis exception", traceback.format_exc())

if __name__ == "__main__":
    main()
