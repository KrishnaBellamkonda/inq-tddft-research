#!/usr/bin/env python3
"""Build + execute the A/B pilot notebook: DIRECT erf/r potential (pilot_direct)
vs POISSON charge perturbation (pilot).

Headline question: does the abrupt exit transient (the ~+40 eV/step spike in per-step
ΔE_total as the clipped periodic charge crosses the far box face at proj_z≈+42.5) DISAPPEAR
when the projectile is added as a DIRECT free-space erf/r potential (no charge/Poisson/
neutralizing background)?

Panels (per prompt):
  1. Density-evolution GIF for the DIRECT run at the TOP (base64-embedded).
  2. Exit transient: per-step ΔE_total vs proj_z near the far-face crossing (proj_z 41..44)
     for BOTH runs — the headline check.
  3. Raw energy_total(t) for both, overlaid.
  4. S comparison: S_def2 (ΔE_total across slab / 25) and S_KEloss (−ΔKE_proj / 25) for both,
     with each run's OWN E_total(0) baseline (the DIRECT run has a different external-energy
     offset than POISSON, so a shared E_GS is NOT used blindly — both interpretations shown).
  5. Trajectory z(t)/vz(t) overlay (should be near-identical if the force is unchanged).

Self-contained: inqview only (load_vti / density_gifs / lindhard). VTIs via load_vti
(physical order, never fftshift'd — vti-coordinate-mapping rule). Figures PNG.

Usage: build_ab_notebook.py
"""
import sys, os, json

ROOT = "/local/data/public/skcb2/tddft"
SYS = f"{ROOT}/ResearchProject/systems/localised_jellium"
SCRIPTS = f"{SYS}/scripts/classical_highdensity_sv"
HYP = f"{SYS}/hypotheses/classical_highdensity_sv/pilot_direct"
os.makedirs(HYP, exist_ok=True)

RES_DIR = f"{SCRIPTS}/pilot_direct/results/pilot_direct"       # DIRECT run outputs
RES_POI = f"{SCRIPTS}/pilot/results/pilot"                     # POISSON run outputs
FRAMES_DIR = f"{RES_DIR}/frames/total"

import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

nb = new_notebook()
cells = []
def md(t): cells.append(new_markdown_cell(t))
def code(t): cells.append(new_code_cell(t))

md(r"""# Phase-3 pilot A/B — DIRECT erf/r potential vs POISSON charge perturbation

**Campaign:** classical-highdensity-sv (high-density localised jellium slab, r_s≈4.18).
**System:** 35×35×85 Bohr box, 25-Bohr slab (faces ±12.5), N=100, dx=0.5, periodicity 2.
**Projectile:** mass-1 electron, charge −1, σ_WP=0.5 (σ_pot=0.35355), launched z=−30, v=2
(K0=2, E_kin=2 Ha=54 eV). Ehrenfest — the light electron decelerates by design.

**The two runs (identical in EVERY physical parameter except the projectile representation):**

| | POISSON (`pilot`) | DIRECT (`pilot_direct`) |
|---|---|---|
| perturbation | `moving_gaussian_projectile_perturbation` = poisson(gaussian_density) | `moving_gaussian_projectile_potential` = erf(\|r−R\|/(√2σ))/\|r−R\| added directly |
| built via | charge density + Poisson solve + periodic G=0 neutralizing background | direct free-space potential, NO charge, NO Poisson, NO background |
| Ehrenfest force | `projectile_force_analytic_z` (−∫poisson(n_proj)·∇n) | `projectile_force_direct_z` (−∫erf/r·∇n) |

**Hypothesis.** The Poisson representation carries a periodic neutralizing background whose
G=0 offset lurches abruptly as the clipped periodic charge crosses the far *box* face
(proj_z≈+42.5=Lz/2) — the ~+40 eV/step spike in per-step ΔE_total. The DIRECT erf/r potential
has no charge in the cell to clip, so **the transient should be GONE**. The force gradient is
insensitive to that offset, so **the trajectory should be near-identical**, and **S should be
similar** (near-field dominates).

**Headline check (Section 2):** per-step ΔE_total near the far-face crossing.
""")

code(f"""%matplotlib inline
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt
ROOT = "{ROOT}"
RES_DIR = "{RES_DIR}"        # DIRECT run
RES_POI = "{RES_POI}"        # POISSON run
FRAMES_DIR = "{FRAMES_DIR}"
HYP = "{HYP}"
HA = 27.211386
SLAB_HALF = 12.5
L_SLAB = 25.0
BOX_HALF = 42.5              # Lz/2 — the far box face where the clipped charge lurches
DT = 0.04
V0 = 2.0; MASS = 1.0; RS = 4.183
E_GS = 207.18322156141       # clean-slab GS E_total (Ha), the Poisson S_def2 reference
os.makedirs(HYP, exist_ok=True)
print("DIRECT  results:", RES_DIR, "exists:", os.path.isdir(RES_DIR))
print("POISSON results:", RES_POI, "exists:", os.path.isdir(RES_POI))
print("DIRECT frames:", len(glob.glob(FRAMES_DIR + "/*.vti")))""")

# ---------------- density GIF for DIRECT at TOP ----------------
md("""## 1. Density-evolution GIF — DIRECT run (visual intuition, read first)

n(x,z,t) on the mid-y x–z plane, as **total density** and **induced Δn = n(t)−n(0)** (the
projectile-driven wake) for the DIRECT run. VTIs via `inqview.load_vti` (physical order —
never fftshift'd). Slab faces ±12.5 dashed. Base64-embedded so it animates on reopen.""")
code("""from inqview.visualisation.density_gifs import _slice_stack, _save_gif
from IPython.display import Image, display
gif_paths = []
files = sorted(glob.glob(FRAMES_DIR + "/*.vti"))
if not files:
    print("no density frames for DIRECT run — GIF skipped")
else:
    nfiles = len(files); frames_max = 30
    idx = list(range(0, nfiles, max(1, nfiles // frames_max)))
    times, tot, axes = _slice_stack(FRAMES_DIR, idx, DT)
    cap = (SLAB_HALF, -SLAB_HALF)
    p_tot = os.path.join(HYP, "pilot_direct_total_density.gif")
    _save_gif(tot, times, axes, p_tot, title="pilot_direct · total density n(x,z,t)",
              cap_lines=cap, kind="density", fps=10); gif_paths.append(p_tot)
    p_ind = os.path.join(HYP, "pilot_direct_induced_delta.gif")
    _save_gif(tot - tot[0][None], times, axes, p_ind,
              title="pilot_direct · induced Δn = n(t) − n(0)", cap_lines=cap, kind="diff", fps=10)
    gif_paths.append(p_ind)
    print("wrote:", gif_paths)""")
md("### Total density  n(x,z,t) — DIRECT")
code("""display(Image(filename=gif_paths[0])) if gif_paths else print("no frames")""")
md("### Induced wake  Δn = n(t) − n(0) — DIRECT")
code("""display(Image(filename=gif_paths[1])) if len(gif_paths) > 1 else print("no frames")""")

# ---------------- load both runs ----------------
md("""## 2. Load both runs

For each run merge observables (E_total etc.) with the projectile track (proj_z, proj_vz,
KE_proj) on `step`. E_total is INQ's native `energy_total`.""")
code("""def load_run(res):
    obs = pd.read_csv(res + "/raw/observables/observables.csv")
    proj = pd.read_csv(res + "/raw/observables/projectile.csv")
    df = obs.merge(proj, on="step", suffixes=("", "_p"))
    df["t"] = df["time_au"]
    return df
D = load_run(RES_DIR)   # DIRECT
P = load_run(RES_POI)   # POISSON
for name, df in [("DIRECT ", D), ("POISSON", P)]:
    print(f"{name}: {len(df):4d} steps | E_total(0)={df.energy_total.iloc[0]:.4f} Ha | "
          f"z0={df.proj_z.iloc[0]:.1f} z_final={df.proj_z.iloc[-1]:.2f} vz_final={df.proj_vz.iloc[-1]:.4f}")
print(f"\\nE_total(0) DIRECT − POISSON = {(D.energy_total.iloc[0]-P.energy_total.iloc[0])*HA:.2f} eV")
print("(= the constant charge-dependent background/self-energy offset; it cancels in ΔE_total.)")""")

# ---------------- THE headline: exit transient ----------------
md("""## 3. The exit transient — per-step ΔE_total near the far-face crossing (HEADLINE)

Per-step ΔE_total = E_total(step) − E_total(step−1), in eV, plotted vs proj_z over the
far-face window (proj_z 41..44 = approach to the box edge at +42.5). POISSON is expected to
**spike to ~+40 eV/step** (the G=0-background lurch as the clipped charge crosses the box
face). DIRECT should stay **smooth** (no spike). *Is the abrupt change gone?*""")
code("""def per_step_dE(df):
    z = df.proj_z.values
    dE = np.diff(df.energy_total.values) * HA     # eV/step
    zc = 0.5*(z[1:] + z[:-1])                       # midpoint proj_z
    return zc, dE
zD, dED = per_step_dE(D)
zP, dEP = per_step_dE(P)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4))
# left: full window; right: zoom proj_z 40..45
for ax, (lo, hi, ttl) in zip((a1, a2), [(-32, 45, "full run"), (40, 45, "far-face zoom (box edge 42.5)")]):
    mP = (zP >= lo) & (zP <= hi); mD = (zD >= lo) & (zD <= hi)
    ax.plot(zP[mP], dEP[mP], lw=1.3, c="C3", label="POISSON (pilot)")
    ax.plot(zD[mD], dED[mD], lw=1.3, c="C0", label="DIRECT (pilot_direct)")
    ax.axvline(SLAB_HALF, ls=":", c="0.6"); ax.axvline(BOX_HALF, ls="--", c="0.4")
    ax.axhline(0, lw=0.6, c="0.7")
    ax.set_xlabel("proj_z (Bohr)"); ax.set_ylabel("per-step ΔE_total (eV/step)"); ax.set_title(ttl)
    ax.legend(fontsize=8)
a2.annotate("box face +42.5", (BOX_HALF, a2.get_ylim()[1]*0.9), fontsize=7, ha="right")
plt.tight_layout(); plt.savefig(os.path.join(HYP, "ab_exit_transient.png"), dpi=120); plt.show()

# quantify: peak |per-step ΔE| in the far-face window (proj_z 40..44)
def peak_far(zc, dE):
    m = (zc >= 40) & (zc <= 44)
    if not m.any(): return np.nan, np.nan
    i = np.argmax(np.abs(dE[m]))
    return dE[m][i], zc[m][i]
pP, zpP = peak_far(zP, dEP); pD, zpD = peak_far(zD, dED)
print(f"POISSON: peak per-step ΔE in proj_z 40..44 = {pP:+.2f} eV/step at proj_z={zpP:.2f}")
print(f"DIRECT : peak per-step ΔE in proj_z 40..44 = {pD:+.2f} eV/step at proj_z={zpD:.2f}")
# typical (median |dE|) in that window for DIRECT (smoothness)
mD = (zD >= 40) & (zD <= 44)
print(f"DIRECT : median |per-step ΔE| in that window = {np.median(np.abs(dED[mD])):.3f} eV/step")
print(f"\\nVERDICT: transient {'GONE' if abs(pD) < 5 else 'STILL PRESENT'} for DIRECT "
      f"(|peak|={abs(pD):.2f} eV vs POISSON |peak|={abs(pP):.2f} eV).")""")

# ---------------- raw E_total(t) overlay ----------------
md("""## 4. Raw energy_total(t) — both runs overlaid

INQ's native `energy_total`. The DIRECT curve sits at a higher absolute baseline (the missing
neutralizing-background offset), so we show BOTH the raw curves and each shifted to its own
t=0 so the *shapes* can be compared.""")
code("""fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4))
a1.plot(P.t, P.energy_total, lw=1.3, c="C3", label="POISSON")
a1.plot(D.t, D.energy_total, lw=1.3, c="C0", label="DIRECT")
a1.set_xlabel("t (a.u.)"); a1.set_ylabel("energy_total (Ha)"); a1.set_title("raw E_total(t)"); a1.legend()
a2.plot(P.t, (P.energy_total - P.energy_total.iloc[0])*HA, lw=1.3, c="C3", label="POISSON")
a2.plot(D.t, (D.energy_total - D.energy_total.iloc[0])*HA, lw=1.3, c="C0", label="DIRECT")
a2.set_xlabel("t (a.u.)"); a2.set_ylabel("E_total − E_total(0)  (eV)")
a2.set_title("E_total shifted to own t=0 (shape compare)"); a2.legend()
plt.tight_layout(); plt.savefig(os.path.join(HYP, "ab_energy_total.png"), dpi=120); plt.show()
print("POISSON E_total: 0 ->", f"{P.energy_total.iloc[0]:.4f}", "final", f"{P.energy_total.iloc[-1]:.4f} Ha")
print("DIRECT  E_total: 0 ->", f"{D.energy_total.iloc[0]:.4f}", "final", f"{D.energy_total.iloc[-1]:.4f} Ha")""")

# ---------------- S comparison ----------------
md(r"""## 5. Stopping-power comparison

Two estimators, computed the same way for both runs:

- **S_def2 = [E_total(plateau) − E_total(0)] / 25.** The energy deposited in the slab per Bohr.
  Because the DIRECT run has a *different* constant external-energy offset than POISSON, the
  shared `E_GS`=207.183 Ha is only valid for POISSON. We therefore report S_def2 with **each
  run's OWN E_total(0) baseline** — a run-independent measure — and separately show the
  POISSON-style `plateau − E_GS`/25 for reference. `plateau` = mean E_total over the flat tail
  once the projectile has cleared the box (proj_z > 44.5).
- **S_KEloss = −ΔKE_proj / 25**, ΔKE_proj = KE_proj(exit) − KE_proj(entry) across the slab
  (entry z=−12.5, exit z=+12.5). Box-offset-free by construction.

POISSON reference (measured): S_def2≈1.08 eV/Bohr, S_KEloss≈0.93 eV/Bohr.

**Baseline subtlety (important).** For POISSON, E_total(0) is measured with the projectile at
z=−30 where its charge already imposes the negative neutralizing-background offset, so
E_total(0)=191.8 Ha ≠ the clean-slab GS (207.18 Ha). The plateau (208.17 Ha) is measured after
the charge has fully left (offset gone). Hence the *physical* POISSON deposition uses
`plateau − E_GS` (≈1.08 eV/Bohr), NOT `plateau − E_total(0)` (which spuriously includes the
17-eV/Bohr offset swing). The DIRECT run has **no such time-varying offset** (the erf/r
potential carries no neutralizing background), so for DIRECT the `plateau − E_total(0)` measure
is the physically correct one and should agree with `plateau − E_GS`. Read the S_def2 numbers
per-run accordingly; S_KEloss is offset-free for both and is the cleanest cross-run comparison.""")
code("""def plateau_Etot(df):
    z = df.proj_z.values; E = df.energy_total.values
    exited = z > (BOX_HALF + 2.0)                  # proj_z > 44.5, clear of box
    i = np.argmax(exited) if exited.any() else len(E) - max(10, len(E)//10)
    tail = E[i:]
    return tail.mean(), tail.std(), i

def ke_at_z(df, target):
    z = df.proj_z.values; ke = df.energy_proj_ke.values
    for i in range(1, len(z)):
        if z[i-1] <= target <= z[i]:
            f = (target - z[i-1]) / (z[i] - z[i-1] + 1e-30)
            return ke[i-1] + f*(ke[i]-ke[i-1]), i
    return np.nan, -1

rows = []
for name, df in [("POISSON", P), ("DIRECT", D)]:
    E0 = df.energy_total.iloc[0]
    plat, pstd, ip = plateau_Etot(df)
    S_own_ha  = (plat - E0) / L_SLAB
    S_egs_ha  = (plat - E_GS) / L_SLAB
    ke_in, _  = ke_at_z(df, -SLAB_HALF)
    ke_out, _ = ke_at_z(df, +SLAB_HALF)
    S_ke_ha   = -(ke_out - ke_in) / L_SLAB
    rows.append(dict(run=name, E0_ha=E0, plateau_ha=plat, plat_std_ev=pstd*HA,
                     S_def2_own_ev=S_own_ha*HA, S_def2_egs_ev=S_egs_ha*HA,
                     ke_in=ke_in, ke_out=ke_out, S_KEloss_ev=S_ke_ha*HA))
T = pd.DataFrame(rows)
pd.set_option("display.float_format", lambda v: f"{v:.4f}")
print(T.to_string(index=False))

d = T.set_index("run")
print("\\n=== S comparison (eV/Bohr) ===")
print(f"  S_def2 (own baseline):  DIRECT={d.loc['DIRECT','S_def2_own_ev']:.3f}  "
      f"POISSON={d.loc['POISSON','S_def2_own_ev']:.3f}")
print(f"  S_def2 (vs E_GS 207.18): DIRECT={d.loc['DIRECT','S_def2_egs_ev']:.3f}  "
      f"POISSON={d.loc['POISSON','S_def2_egs_ev']:.3f}  (E_GS valid for POISSON only)")
print(f"  S_KEloss:               DIRECT={d.loc['DIRECT','S_KEloss_ev']:.3f}  "
      f"POISSON={d.loc['POISSON','S_KEloss_ev']:.3f}")
# consistency — compare the PHYSICALLY-appropriate S_def2 per run:
#   DIRECT: own baseline (no offset swing);  POISSON: vs E_GS (offset removed).
sd = d.loc['DIRECT','S_def2_own_ev']; sp = d.loc['POISSON','S_def2_egs_ev']
print(f"\\n  physical S_def2: DIRECT(own)={sd:.3f}  vs  POISSON(E_GS)={sp:.3f} eV/Bohr  "
      f"Δ={sd-sp:+.3f} ({100*(sd-sp)/sp:+.1f}%) — {'CONSISTENT' if abs(sd-sp) < 0.25 else 'DIFFERENT'}")
kd = d.loc['DIRECT','S_KEloss_ev']; kp = d.loc['POISSON','S_KEloss_ev']
print(f"  S_KEloss (offset-free): DIRECT={kd:.3f}  vs  POISSON={kp:.3f} eV/Bohr  "
      f"Δ={kd-kp:+.3f} ({100*(kd-kp)/kp:+.1f}%) — {'CONSISTENT' if abs(kd-kp) < 0.25 else 'DIFFERENT'}")""")

md("""### Bar-chart summary of the S estimators""")
code("""fig, ax = plt.subplots(figsize=(8, 4.2))
labels = ["S_def2 (own base)", "S_def2 (vs E_GS)", "S_KEloss"]
xw = np.arange(len(labels)); w = 0.36
poi = [d.loc['POISSON','S_def2_own_ev'], d.loc['POISSON','S_def2_egs_ev'], d.loc['POISSON','S_KEloss_ev']]
dir_ = [d.loc['DIRECT','S_def2_own_ev'], d.loc['DIRECT','S_def2_egs_ev'], d.loc['DIRECT','S_KEloss_ev']]
ax.bar(xw - w/2, poi, w, label="POISSON", color="C3")
ax.bar(xw + w/2, dir_, w, label="DIRECT", color="C0")
for i,(a,b) in enumerate(zip(poi, dir_)):
    ax.text(i-w/2, a, f"{a:.2f}", ha="center", va="bottom", fontsize=7)
    ax.text(i+w/2, b, f"{b:.2f}", ha="center", va="bottom", fontsize=7)
ax.set_xticks(xw); ax.set_xticklabels(labels); ax.set_ylabel("S (eV/Bohr)")
ax.set_title("Stopping-power estimators: DIRECT vs POISSON"); ax.legend()
plt.tight_layout(); plt.savefig(os.path.join(HYP, "ab_S_comparison.png"), dpi=120); plt.show()""")

# ---------------- trajectory overlay ----------------
md("""## 6. Trajectory overlay — z(t) and vz(t)

If the force is genuinely unchanged (gradient insensitive to the background offset), the
DIRECT and POISSON trajectories should be near-identical.""")
code("""fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4))
a1.plot(P.t, P.proj_z, lw=1.5, c="C3", label="POISSON")
a1.plot(D.t, D.proj_z, lw=1.4, ls="--", c="C0", label="DIRECT")
a1.axhline(SLAB_HALF, ls=":", c="0.6"); a1.axhline(-SLAB_HALF, ls=":", c="0.6")
a1.set_xlabel("t (a.u.)"); a1.set_ylabel("proj_z (Bohr)"); a1.set_title("z(t)"); a1.legend()
a2.plot(P.t, P.proj_vz, lw=1.5, c="C3", label="POISSON")
a2.plot(D.t, D.proj_vz, lw=1.4, ls="--", c="C0", label="DIRECT")
a2.set_xlabel("t (a.u.)"); a2.set_ylabel("proj_vz (a.u.)"); a2.set_title("vz(t)"); a2.legend()
plt.tight_layout(); plt.savefig(os.path.join(HYP, "ab_trajectory.png"), dpi=120); plt.show()
# agreement on common time support
tmax = min(D.t.max(), P.t.max()); ti = np.linspace(0, tmax, 400)
zD = np.interp(ti, D.t, D.proj_z); zP = np.interp(ti, P.t, P.proj_z)
vD = np.interp(ti, D.t, D.proj_vz); vP = np.interp(ti, P.t, P.proj_vz)
print(f"z(t)  DIRECT vs POISSON: max|Δz|={np.abs(zD-zP).max():.4f} Bohr, RMS={np.sqrt(np.mean((zD-zP)**2)):.4f} Bohr")
print(f"vz(t) DIRECT vs POISSON: max|Δvz|={np.abs(vD-vP).max():.5f} a.u., RMS={np.sqrt(np.mean((vD-vP)**2)):.5f} a.u.")
print(f"final: z DIRECT={D.proj_z.iloc[-1]:.3f} POISSON={P.proj_z.iloc[-1]:.3f} | "
      f"vz DIRECT={D.proj_vz.iloc[-1]:.4f} POISSON={P.proj_vz.iloc[-1]:.4f}")""")

md("""---
*Auto-built by `build_ab_notebook.py`. Executed end-to-end (0 errors). Figures saved into this
directory; DIRECT density GIFs base64-embedded above. VTIs read via `inqview.load_vti`
(physical order, never fftshift'd).*""")

nb["cells"] = cells
nb["metadata"]["kernelspec"] = {"display_name": "inqview-venv", "language": "python", "name": "inqview-venv"}

import nbformat
from nbclient import NotebookClient
nb_path = os.path.join(HYP, "pilot_direct_ab_notebook.ipynb")
nbformat.write(nb, nb_path)
print("wrote notebook:", nb_path)
client = NotebookClient(nb, timeout=1800, kernel_name="inqview-venv",
                        resources={"metadata": {"path": HYP}})
client.execute()
nbformat.write(nb, nb_path)
print("EXECUTED (0 errors):", nb_path)
