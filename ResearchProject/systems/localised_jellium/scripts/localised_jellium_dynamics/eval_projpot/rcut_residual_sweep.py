"""Residual d(E_H+E_ext) − U_proj_bg vs r_cut for the PSEUDOPOTENTIAL method.

At t=0 the classical electron density = bare GS (2-step run), so
  (E_H+E_ext)_CL = (E_H+E_ext)_GS + e_proj(r_cut),   e_proj = ∫ n_GS·v_ion
  residual(r_cut) = (E_H+E_ext)_WP − (E_H+E_ext)_GS − e_proj − impl
                  = 15.617 eV − e_proj(r_cut) − impl(r_cut)
with 15.617 eV = (E_H+E_ext)_WP[p5] − (E_H+E_ext)_bareGS[eval_gs_xc]
   = (−85.038+151.725) − (−80.664+146.777) Ha = 0.574 Ha.
Runs eval_projpot WITH the GS loaded (so it prints e_proj) at fine r_cut on GPU0.
"""
import re, os, subprocess, tempfile
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BIN  = os.path.join(HERE, "run")
SRC  = ("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/"
        "pseudopotentials/electron_gaussian_wpsigma0p5_rc120.upf")
GS   = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "scripts/campaign_autorun/runs/h2/gs_p2_lz120/checkpoint")
CONST = 15.617   # eV
DR, T = 0.01, open(SRC).read()

def truncate_upf(r_cut, out):
    n = int(round(r_cut / DR)) + 1; t = T
    for tag in ("PP_R", "PP_RAB", "PP_LOCAL", "PP_RHOATOM"):
        m = re.search(rf'(<{tag}\b[^>]*>)(.*?)(</{tag}>)', t, re.S)
        if not m: continue
        head = re.sub(r'size="\s*\d+"', f'size="{n}"', m.group(1))
        vals = m.group(2).split()[:n]; ncol = 4 if tag in ("PP_LOCAL","PP_RHOATOM") else 8
        body = "\n".join(" ".join(vals[i:i+ncol]) for i in range(0,len(vals),ncol))
        t = t[:m.start()] + head + "\n" + body + "\n" + m.group(3) + t[m.end():]
    open(out,"w").write(re.sub(r'mesh_size="\s*\d+"', f'mesh_size="{n}"', t))

env = dict(os.environ, CUDA_VISIBLE_DEVICES="0",
           INQ_SHARE_PATH="/local/data/public/skcb2/tddft/inq/install/share",
           PSEUDOPOD_SHARE_PATH="/local/data/public/skcb2/tddft/inq/install/share/pseudopod",
           LJ_GS_DIR=GS, LJ_LAUNCH_Z="-24.5", LJ_SPACING="0.5")

def run_one(r_cut):
    with tempfile.NamedTemporaryFile("w", suffix=".upf", delete=False) as f: upf=f.name
    truncate_upf(r_cut, upf)
    out = subprocess.run([BIN], env=dict(env, LJ_PROJ_UPF=upf), capture_output=True, text=True).stdout
    os.unlink(upf)
    g = lambda k: float(re.search(rf'^  {k}\b.*?=\s*(\S+)\s*eV', out, re.M).group(1))
    return g("ideal"), g("impl"), g("e_proj")

RCUTS = [10,12,14,16,20,25,30,37,40,45,50,60,70,80,90,100,110,120]
rows = []
for rc in RCUTS:
    idl, imp, ep = run_one(rc)
    res = CONST - ep - imp
    rows.append(dict(r_cut=rc, ideal=idl, impl=imp, e_proj=ep, residual=res))
    print(f"r_cut={rc:4}  impl={imp:9.2f}  e_proj={ep:9.2f}  residual={res:7.2f} eV")
df = pd.DataFrame(rows); df.to_csv(os.path.join(HERE,"rcut_residual_sweep.csv"), index=False)

fig, ax = plt.subplots(figsize=(7.8,4.8))
ax.plot(df.r_cut, df.residual, 'o-', color='#1b6ca8',
        label='residual d(E_H+E_ext) − U_proj_bg  (pseudopotential)')
ax.axhline(20.81, ls='--', color='#c0392b', label='WP self-Hartree (clean target, 20.8 eV)')
ax.axhline(0, color='0.6', lw=0.6)
ax.set_xlabel('pseudopotential radial cutoff r_cut (Bohr)'); ax.set_ylabel('residual (eV)')
ax.set_title('The ACTUAL observable drifts with r_cut (52→35% of the clean value) — pseudopotential contamination')
ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(HERE,"rcut_residual_sweep.png"), dpi=140)
print("\nwrote rcut_residual_sweep.csv + .png")
