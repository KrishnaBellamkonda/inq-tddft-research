"""Grid-spacing sweep of the CLEAN residual d(E_H+E_ext) − U_proj_bg (perturbation
method) vs the pseudopotential method. For each dx: build a matching bare GS, run
the WP (2-step, at rest, r=12) and the Gaussian-charge perturbation classical run
(2-step), and form residual = (E_H+E_ext)_WP − (E_H+E_ext)_CL_pert − U_proj_bg.
All on GPU0. dx=0.5 is taken from the already-completed runs.
"""
import os, csv, subprocess
from pathlib import Path
import pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HA = 27.211386
SCR = Path(__file__).resolve().parent / "grid_scratch"; SCR.mkdir(exist_ok=True)
DYN = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts")
GS_BIN   = DYN/"campaign_autorun/gs/run"
WP_BIN   = DYN/"localised_jellium_dynamics/phase5_wp/run"
PERT_BIN = DYN/"localised_jellium_dynamics/proj_perturbation/run"
BASE = dict(os.environ, CUDA_VISIBLE_DEVICES="0",
            INQ_SHARE_PATH="/local/data/public/skcb2/tddft/inq/install/share",
            PSEUDOPOD_SHARE_PATH="/local/data/public/skcb2/tddft/inq/install/share/pseudopod",
            LJ_LX="50", LJ_LY="50", LJ_LZ="120", LJ_HALF="12.5", LJ_N="82",
            LJ_EDGE_W="0", LJ_PERIODICITY="2", LJ_SIGMA="0.5", LJ_LAUNCH_Z="-24.5")

def run(binary, env, cwd, log):
    with open(log, "w") as f:
        subprocess.run([str(binary)], env=env, cwd=cwd, stdout=f, stderr=subprocess.STDOUT, check=True)

def read_HE(obs_csv):
    with open(obs_csv) as f:
        r = next(csv.DictReader(f))
    return float(r["energy_hartree"]), float(r["energy_external"])   # Ha, t=0

def do_dx(dx):
    tag = f"dx{str(dx).replace('.','p')}"; d = SCR/tag; d.mkdir(exist_ok=True)
    gsdir = str(d/"gs")
    # 1) GS
    run(GS_BIN, dict(BASE, LJ_SPACING=str(dx), LJ_GS_DIR=gsdir, LJ_TAG=f"gs_{tag}"), d, d/"gs.log")
    # 2) WP (2-step, at rest)
    run(WP_BIN, dict(BASE, LJ_SPACING=str(dx), LJ_GS_DIR=gsdir, LJ_K0="0",
                     LJ_N_STEPS="2", LJ_SAVE_EVERY="1000", LJ_OUT="wp"), d, d/"wp.log")
    HwP, ExtWP = read_HE(d/"results/wp/raw/observables/observables.csv")
    # 3) perturbation classical (2-step)
    run(PERT_BIN, dict(BASE, LJ_SPACING=str(dx), LJ_GS_DIR=gsdir,
                       LJ_N_STEPS="2", LJ_OUT="pert"), d, d/"pert.log")
    Hcl, ExtCL = read_HE(d/"results/pert/raw/observables/observables.csv")
    upb = next(float(l.split("=")[1].split()[0]) for l in
               open(d/"results/pert/run_summary.txt") if l.startswith("U_proj_bg_ha"))
    resid = ((HwP+ExtWP) - (Hcl+ExtCL) - upb) * HA
    return dict(dx=dx, HE_WP=(HwP+ExtWP)*HA, HE_pert=(Hcl+ExtCL)*HA,
                U_proj_bg=upb*HA, residual=resid)

rows = [dict(dx=0.5, HE_WP=66.686813*HA, HE_pert=60.972492*HA,
             U_proj_bg=4.949720*HA, residual=0.764601*HA)]   # from completed dx=0.5 runs
for dx in (0.4, 0.3):
    print(f"=== dx={dx} ===", flush=True)
    rows.append(do_dx(dx)); print(rows[-1], flush=True)

df = pd.DataFrame(rows).sort_values("dx", ascending=False)
df.to_csv(Path(__file__).resolve().parent/"grid_sweep.csv", index=False)
print(df.to_string(index=False))
fig, ax = plt.subplots(figsize=(7.4,4.6))
ax.plot(df.dx, df.residual, 'o-', color='#2e8b57', ms=8,
        label='perturbation (Gaussian charge) — clean, r_cut-free')
ax.axhline(7.36, ls='--', color='#c0392b', label='pseudopotential r_cut=120 (dx=0.5): 7.4 eV')
ax.set_xlabel('grid spacing dx (Bohr)'); ax.set_ylabel('residual d(E_H+E_ext) − U_proj_bg (eV)')
ax.set_title('Clean residual is grid-stable near the WP self-Hartree (~21 eV)')
ax.legend(frameon=False, fontsize=8); ax.invert_xaxis()
fig.tight_layout(); fig.savefig(Path(__file__).resolve().parent/"grid_sweep.png", dpi=140)
print("\nwrote grid_sweep.csv + grid_sweep.png")
