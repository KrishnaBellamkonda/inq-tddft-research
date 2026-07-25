"""Stress-test the Gaussian-charge PERTURBATION projectile method (campaign
localised-jellium-dynamics-analysis). For each config (sigma, r, Lz, periodicity)
run the WP (2-step, at rest) and the perturbation classical run (2-step), read the
full energy decomposition + clean U_proj_bg, and form the ledger:
  residual = (E_H+E_ext)_WP − (E_H+E_ext)_pert − U_proj_bg
  dKin = E_kin_WP − E_kin_pert ;  dXC = E_xc_WP − E_xc_pert ;  SIE = residual + dXC
Writes per-axis CSVs to hypotheses/perturbation_method/. All GS reused; GPU0 only.
Waits for grid_sweep.py (also GPU0) to finish first.
"""
import os, csv, time, subprocess, traceback
from pathlib import Path

HA = 27.211386
HERE = Path(__file__).resolve().parent
OUT  = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
            "hypotheses/perturbation_method")
SCR  = HERE/"stress_scratch"; SCR.mkdir(exist_ok=True)
DYN  = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts")
WP_BIN   = DYN/"localised_jellium_dynamics/phase5_wp/run"
PERT_BIN = DYN/"localised_jellium_dynamics/proj_perturbation/run"
SEMI = DYN/"semiempirical_spillout/runs"
GS = {(120,2): str(DYN/"campaign_autorun/runs/h2/gs_p2_lz120/checkpoint"),
      (90,2):  str(SEMI/"lz90/checkpoint"),  (160,2): str(SEMI/"lz160/checkpoint"),
      (240,2): str(SEMI/"lz240/checkpoint"), (120,3): str(SEMI/"p3_lz120/checkpoint")}
BASE = dict(os.environ, CUDA_VISIBLE_DEVICES="0",
            INQ_SHARE_PATH="/local/data/public/skcb2/tddft/inq/install/share",
            PSEUDOPOD_SHARE_PATH="/local/data/public/skcb2/tddft/inq/install/share/pseudopod",
            LJ_LX="50", LJ_LY="50", LJ_HALF="12.5", LJ_N="82", LJ_EDGE_W="0", LJ_SPACING="0.5")

def sh(binary, env, cwd, log):
    with open(log,"w") as f:
        subprocess.run([str(binary)], env=env, cwd=cwd, stdout=f, stderr=subprocess.STDOUT, check=True)

def energies(obs_csv):
    with open(obs_csv) as f: r = next(csv.DictReader(f))
    return {k: float(r["energy_"+k])*HA for k in ("kinetic","hartree","xc","external")}

def do(sigma, r, Lz, per):
    tag = f"s{sigma}_r{r}_lz{Lz}_p{per}".replace(".","p")
    d = SCR/tag; d.mkdir(exist_ok=True)
    gs = GS[(Lz, per)]; launch_z = -(12.5 + r)
    env = dict(BASE, LJ_LZ=str(Lz), LJ_PERIODICITY=str(per), LJ_SIGMA=str(sigma),
               LJ_LAUNCH_Z=str(launch_z), LJ_GS_DIR=gs)
    sh(WP_BIN,   dict(env, LJ_K0="0", LJ_N_STEPS="2", LJ_SAVE_EVERY="1000", LJ_OUT="wp"), d, d/"wp.log")
    sh(PERT_BIN, dict(env, LJ_N_STEPS="2", LJ_OUT="pert"), d, d/"pert.log")
    wp = energies(d/"results/wp/raw/observables/observables.csv")
    pt = energies(d/"results/pert/raw/observables/observables.csv")
    upb = next(float(l.split("=")[1].split()[0])*HA for l in
               open(d/"results/pert/run_summary.txt") if l.startswith("U_proj_bg_ha"))
    HE_wp, HE_pt = wp["hartree"]+wp["external"], pt["hartree"]+pt["external"]
    resid = HE_wp - HE_pt - upb
    dKin, dXC = wp["kinetic"]-pt["kinetic"], wp["xc"]-pt["xc"]
    return dict(sigma=sigma, r=r, Lz=Lz, per=per,
                Ekin_WP=wp["kinetic"], EH_WP=wp["hartree"], Exc_WP=wp["xc"], Eext_WP=wp["external"],
                Ekin_pert=pt["kinetic"], EH_pert=pt["hartree"], Exc_pert=pt["xc"], Eext_pert=pt["external"],
                dEH_Eext=HE_wp-HE_pt, U_proj_bg=upb, residual=resid, dKin=dKin, dXC=dXC, SIE=resid+dXC)

# ---- wait for the grid sweep (same GPU) ----
while subprocess.run(["pgrep","-f","grid_sweep.py"], capture_output=True).returncode == 0:
    print("waiting for grid_sweep.py to free GPU0 ...", flush=True); time.sleep(30)

CONFIGS = {
    "baseline": [(0.5,12,120,2)],
    "sigma":    [(0.35,12,120,2),(0.5,12,120,2),(0.7,12,120,2),(1.0,12,120,2)],
    "r":        [(0.5,4,120,2),(0.5,12,120,2),(0.5,20,120,2),(0.5,28,120,2)],
    "lz":       [(0.5,12,90,2),(0.5,12,120,2),(0.5,12,160,2),(0.5,12,240,2)],
    "p3vp2":    [(0.5,12,120,2),(0.5,12,120,3)],
}
cache, cols = {}, ["sigma","r","Lz","per","Ekin_WP","EH_WP","Exc_WP","Eext_WP",
    "Ekin_pert","EH_pert","Exc_pert","Eext_pert","dEH_Eext","U_proj_bg","residual","dKin","dXC","SIE"]
def get(cfg):
    if cfg not in cache:
        print(f"RUN {cfg}", flush=True)
        try: cache[cfg] = do(*cfg)
        except Exception: print("FAIL", cfg, "\n", traceback.format_exc(), flush=True); cache[cfg] = None
    return cache[cfg]

for axis, cfgs in CONFIGS.items():
    rows = [get(c) for c in cfgs]; rows = [r for r in rows if r]
    with open(OUT/f"stress_{axis}.csv","w",newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows: w.writerow({k:r[k] for k in cols})
    print(f"wrote stress_{axis}.csv ({len(rows)} rows)", flush=True)
print("STRESS_DONE", flush=True)
