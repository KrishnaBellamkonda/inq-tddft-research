#!/usr/bin/env python3
"""Build the effective-mass quantum-vs-classical comparison ("phase") notebook +
per-run notebooks, from the two completed runs. Called by orchestrate.py; also
runnable standalone. Defensive: missing files/columns degrade to notes, never crash.
"""
import sys, traceback
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

ROOT = Path("/local/data/public/skcb2/tddft")
LJ   = ROOT/"ResearchProject/systems/localised_jellium"
HERE = LJ/"scripts/muon_mass_fork/effmass_pair"
NB   = LJ/"hypotheses/muon_mass_fork"; NB.mkdir(parents=True, exist_ok=True)
QDIR = HERE/"quantum/results/quantum"
CDIR = HERE/"classical/results/classical"
try:
    from inqview.visualisation.style import apply_theme; apply_theme()
except Exception:
    plt.rcParams.update({"figure.dpi":130,"axes.grid":True,"grid.alpha":.3,"font.size":10})

HA_EV = 27.211386
def rd(p):
    p = Path(p)
    return pd.read_csv(p) if p.exists() else None

# ---- load ------------------------------------------------------------------
q_obs = rd(QDIR/"raw/observables/observables.csv")
q_rs  = rd(QDIR/"raw/observables/wp_real_space_stats.csv")
q_ms  = rd(QDIR/"raw/observables/wp_momentum_stats.csv")
c_trk = rd(CDIR/"raw/observables/electron_track.csv")
c_obs = rd(CDIR/"raw/observables/observables.csv")

def col(df, *names):
    if df is None: return None
    for n in names:
        if n in df.columns: return df[n].values
    return None

# ---- comparison figure: projectile KE vs path (the stopping) ---------------
fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
a, b = axs
# classical: ion KE vs z from the track
cz, cke = col(c_trk,"z"), col(c_trk,"ke_ion_ha")
if cz is not None and cke is not None:
    a.plot(cz, cke*HA_EV, color="#b5651d", lw=2, label="classical (ion KE)")
# quantum: WP mean KE from momentum stats (k^2/2m) vs centroid z
qz = col(q_rs,"z_mean"); qk2 = col(q_ms,"k2_mean","kz2_mean")
INV_M = 0.324127; M = 1/INV_M
if qz is not None and qk2 is not None:
    n = min(len(qz), len(qk2))
    a.plot(qz[:n], (0.5*qk2[:n]*INV_M)*HA_EV, color="#1f4e79", lw=2, label="quantum (WP ⟨T⟩)")
a.axvspan(-12.5,12.5,color="grey",alpha=.13); a.text(0,a.get_ylim()[1]*0.05 if a.get_ylim()[1] else 0,"slab",ha="center",fontsize=8,color="grey")
a.set_xlabel("projectile z [Bohr]"); a.set_ylabel("projectile KE [eV]")
a.set_title("Stopping: projectile KE vs path"); a.legend(fontsize=8, frameon=False)

# system (bath) energy vs time — both
qt, qE = col(q_obs,"time_au"), col(q_obs,"energy_total")
ct, cE = col(c_obs,"time_au"), col(c_obs,"energy_total")
if qt is not None: b.plot(qt, (qE-qE[0])*HA_EV, color="#1f4e79", lw=2, label="quantum ΔE_total")
if ct is not None: b.plot(ct, (cE-cE[0])*HA_EV, color="#b5651d", lw=2, label="classical ΔE_total")
b.set_xlabel("time [a.u.]"); b.set_ylabel("ΔE_total [eV]")
b.set_title("System energy vs time"); b.legend(fontsize=8, frameon=False)
fig.suptitle("Effective-mass projectile (m=3.09, v=2.71, σ_WP=2) — quantum WP vs classical, r_s=5.665 slab", fontsize=10)
fig.tight_layout(rect=(0,0,1,0.95))
png = NB/"effmass_pair_stopping.png"; fig.savefig(png, dpi=150, bbox_inches="tight")
print("wrote", png)

# ---- comparison notebook ---------------------------------------------------
C=[]
C.append(new_markdown_cell(r"""
# Effective-mass projectile — quantum WP vs classical (model comparison)

Contender **D**: a projectile of effective mass **m=3.09 mₑ** at velocity **v=2.71 a.u.**
(= 100 eV electron ⇒ same S(v)), momentum k₀=8.36, through the **r_s=5.665** jellium slab
(N=82, 50×50×90, dx=0.333), LDA/ETRS, **dt=0.04**, CAP.

Two realisations of the *same* projectile:
- **Quantum** — a Gaussian wavepacket (σ_WP=2) in the last extra state with `inverse_mass=1/3.09`
  (the inq-study mass fork). Stopping read from ⟨T⟩=½⟨k²⟩/m.
- **Classical** — a Gaussian-charge ion (`electron_gaussian_wpsigma2p0.upf`, σ_pot=1.414=σ_WP/√2)
  with `.mass(3.09/1822.8885)` amu, launched at v=2.71 (matched momentum). Ehrenfest; stopping
  from d(KE_ion)/dz.

The pair isolates the **quantum vs classical** difference at identical mass, momentum, density.
Runs: `scripts/muon_mass_fork/effmass_pair/{quantum,classical}`.
""".strip()))
C.append(new_markdown_cell("## Stopping & energy comparison"))
C.append(new_code_cell("from IPython.display import Image\nImage('effmass_pair_stopping.png')".replace(
    "effmass_pair_stopping.png", str(png))))
C.append(new_markdown_cell(r"""
**Left:** projectile KE along its path — the slope through the slab is the stopping power. Quantum
(⟨T⟩ of the WP) vs classical (ion KE). **Right:** total-system energy vs time — the energy deposited
into the bath. Divergence between the two curves is the genuine quantum-vs-classical difference
(the WP also carries a known ~few-eV self-interaction the classical lacks — bound it with the
vacuum-WP control before over-reading small gaps).
"""))
# raw provenance
prov = []
for name, d in [("quantum", QDIR), ("classical", CDIR)]:
    rs = d/"run_summary.txt"
    prov.append(f"### {name}\n```\n{rs.read_text() if rs.exists() else 'MISSING'}\n```")
C.append(new_markdown_cell("## Run provenance\n"+"\n".join(prov)))
nb = new_notebook(); nb["cells"]=C
out = NB/"effmass_pair_comparison.ipynb"; nbf.write(nb, out); print("wrote", out)

# ---- per-run notebooks (minimal: summary + observables plot) ---------------
for name, d in [("quantum", QDIR), ("classical", CDIR)]:
    try:
        rs = (d/"run_summary.txt")
        cells=[new_markdown_cell(f"# effmass_pair — {name} run\n\n```\n{rs.read_text() if rs.exists() else 'MISSING'}\n```"),
               new_markdown_cell("## Observables\nSee the comparison notebook `effmass_pair_comparison.ipynb` for the paired analysis.")]
        n2=new_notebook(); n2["cells"]=cells
        nbf.write(n2, NB/f"effmass_pair_{name}_run.ipynb")
        print("wrote", NB/f"effmass_pair_{name}_run.ipynb")
    except Exception:
        traceback.print_exc()
print("notebooks done")
