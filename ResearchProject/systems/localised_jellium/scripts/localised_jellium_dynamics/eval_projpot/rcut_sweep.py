"""Fine r_cut sweep of E_proj_bg (ideal & impl), at fixed geometry (r=12, z=-24.5).

Truncates the rc120 ghost UPF to each r_cut (keep mesh r<=r_cut; the erf(r/0.5)/r
Coulomb tail is cut to 0 beyond — exactly how the rc10..rc50 UPFs were made) and
runs the eval_projpot binary (no GS: ideal/impl depend only on background +
pseudopotential). Writes rcut_impl_sweep.csv + rcut_impl_sweep.png.
"""
import re, os, subprocess, tempfile
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BIN  = os.path.join(HERE, "run")
SRC  = ("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/"
        "pseudopotentials/electron_gaussian_wpsigma0p5_rc120.upf")
DR   = 0.01
T    = open(SRC).read()

def truncate_upf(r_cut, out):
    n = int(round(r_cut / DR)) + 1                      # keep r = 0 .. r_cut
    t = T
    for tag in ("PP_R", "PP_RAB", "PP_LOCAL", "PP_RHOATOM"):
        m = re.search(rf'(<{tag}\b[^>]*>)(.*?)(</{tag}>)', t, re.S)
        if not m: continue
        head, data, tail = m.group(1), m.group(2), m.group(3)
        vals = data.split()[:n]
        head = re.sub(r'size="\s*\d+"', f'size="{n}"', head)
        ncol = 4 if tag in ("PP_LOCAL", "PP_RHOATOM") else 8
        body = "\n".join(" ".join(vals[i:i+ncol]) for i in range(0, len(vals), ncol))
        t = t[:m.start()] + head + "\n" + body + "\n" + tail + t[m.end():]
    t = re.sub(r'mesh_size="\s*\d+"', f'mesh_size="{n}"', t)
    open(out, "w").write(t)

env = dict(os.environ,
           INQ_SHARE_PATH="/local/data/public/skcb2/tddft/inq/install/share",
           PSEUDOPOD_SHARE_PATH="/local/data/public/skcb2/tddft/inq/install/share/pseudopod",
           LJ_GS_DIR="", LJ_LAUNCH_Z="-24.5")

def run_one(r_cut):
    with tempfile.NamedTemporaryFile("w", suffix=".upf", delete=False) as f:
        upf = f.name
    truncate_upf(r_cut, upf)
    e = dict(env, LJ_PROJ_UPF=upf, LJ_SPACING="0.5")
    out = subprocess.run([BIN], env=e, capture_output=True, text=True).stdout
    os.unlink(upf)
    def grab(key):
        m = re.search(rf'^  {key}\b.*?=\s*(\S+)\s*eV', out, re.M)
        return float(m.group(1)) if m else np.nan
    return grab("ideal"), grab("impl"), grab("gap")

RCUTS = [5,8,10,12,14,16,20,24,25,26,30,35,40,45,50,55,60,70,80,90,100,110,120]
rows = []
for rc in RCUTS:
    idl, imp, gap = run_one(rc)
    rows.append(dict(r_cut=rc, ideal_eV=idl, impl_eV=imp, gap_eV=gap))
    print(f"r_cut={rc:4}  ideal={idl:8.2f}  impl={imp:10.2f}  gap={gap:9.2f}")
df = pd.DataFrame(rows)
df.to_csv(os.path.join(HERE, "rcut_impl_sweep.csv"), index=False)

fig, ax = plt.subplots(figsize=(7.8,4.8))
ax.plot(df.r_cut, df.impl_eV, 'o-', color='#c0392b', label='impl  = −∫ n₊·v_ion  (as-implemented pseudopotential)')
ax.axhline(df.ideal_eV.mean(), ls='--', color='#2e8b57',
           label=f'ideal = ∫ n_proj·v_bg  (true Gaussian, r_cut-invariant, {df.ideal_eV.mean():.0f} eV)')
ax.axhline(0, color='0.6', lw=0.6)
for x,lab in [(12,'proj→near slab face (12)'),(25,'Lx/2 = 25 (lateral wrap onset)'),(37,'proj→far slab face (37)')]:
    ax.axvline(x, ls=':', color='0.5', lw=0.8); ax.text(x, ax.get_ylim()[0], ' '+lab, rotation=90, va='bottom', ha='right', fontsize=7, color='0.4')
ax.set_xlabel('pseudopotential radial cutoff r_cut (Bohr)'); ax.set_ylabel('E_proj_bg (eV)')
ax.set_title('E_proj_bg impl vs r_cut (projectile at r=12, dx=0.5) — grows without bound as the erf/r tail wraps')
ax.legend(frameon=False, fontsize=8, loc='lower left')
fig.tight_layout(); fig.savefig(os.path.join(HERE, "rcut_impl_sweep.png"), dpi=140)
print("\nwrote rcut_impl_sweep.csv + rcut_impl_sweep.png")
