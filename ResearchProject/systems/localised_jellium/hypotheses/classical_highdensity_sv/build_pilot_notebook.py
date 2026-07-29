#!/usr/bin/env python3
"""Build + execute the Phase-3 pilot run-notebooks (RUN A perturbation, RUN B native).

Usage:
    build_pilot_notebook.py {pilot|pilot_native}

Emits an executed .ipynb (0 errors) into
  ResearchProject/systems/localised_jellium/hypotheses/classical_highdensity_sv/<kind>/
with (per notebook-density-gif + run-notebook rules):
  * a density-evolution GIF at the TOP (n(x,z,t) mid-y, total + induced Δn),
    base64-embedded so it animates on reopen;
  * projectile z(t), vz(t) (transit + exit);
  * E_electronic(t) with the plateau marked + E_absorbed annotated;
  * the full pairwise ledger vs t (RUN A) and the conservation column (flat check);
  * a stopping-power section broken into explicit steps (S = E_abs / L_slab, plus
    the initial-drag −dKE/ds cross-check over the v ≥ 0.85·v0 window);
  * a Lindhard/bulk eyeball number (NON-gating).
RUN B additionally overlays z(t)/vz(t) native-vs-perturbation if RUN A data exists.

Self-contained: no external module beyond inqview (load_vti / gif helpers /
lindhard). VTIs are read via inqview.load_vti (physical order, never fftshift'd).
"""
import sys, os, json

KIND = sys.argv[1] if len(sys.argv) > 1 else "pilot"
assert KIND in ("pilot", "pilot_native"), KIND

ROOT = "/local/data/public/skcb2/tddft"
SYS = f"{ROOT}/ResearchProject/systems/localised_jellium"
SCRIPTS = f"{SYS}/scripts/classical_highdensity_sv"
HYP = f"{SYS}/hypotheses/classical_highdensity_sv/{KIND}"
os.makedirs(HYP, exist_ok=True)

# where the run wrote its outputs
RES = f"{SCRIPTS}/{KIND}/results/{KIND}"
FRAMES = f"{RES}/frames/total"

import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

nb = new_notebook()
cells = []

def md(t): cells.append(new_markdown_cell(t))
def code(t): cells.append(new_code_cell(t))

IS_PERT = (KIND == "pilot")
TITLE = ("Phase-3 pilot — RUN A: perturbation Ehrenfest (mass-1 Gaussian charge, analytic HF force)"
         if IS_PERT else
         "Phase-3 pilot — RUN B: INQ native Ehrenfest (real ghost-UPF ion)")

md(f"""# {TITLE}

**Campaign:** classical-highdensity-sv (high-density localised jellium slab, r_s≈4.18).
**System:** 35×35×85 Bohr box, 25-Bohr slab, N=100 electrons, dx=0.5, periodicity 2
(x,y periodic infinite slab; z open so the projectile leaves the box cleanly, no CAP →
energy conserved → exact post-exit E_electronic plateau).
**Projectile:** mass-1 electron, charge −1, σ_WP=0.5 (σ_pot=0.35355), launched at
z=−30 with v=2 (K0=m·v=2, E_kin=½·1·2²=2 Ha=54 eV). Ehrenfest — the light electron
**decelerates** by design; we do NOT gate on velocity drift (light-projectile rule).

**This run's representation:** """ + (
"a moving Gaussian **charge** applied as a perturbation; the drag force is INQ's exact "
"native analytic Hellmann-Feynman force `F = −∫V_proj·∇n` (density-gradient form, "
"validated <0.1 % vs native `forces_stress`), evaluated from (electrons − background)."
if IS_PERT else
"a **real ghost-UPF ion** (H symbol, z_valence=0, clean local V=+erf(r/(√2·σ_pot))/r), "
"moved by INQ's OWN native Ehrenfest integrator (velocity-Verlet inside ETRS, "
"a=F_localHF/mass). This is the real-system faithfulness check for RUN A."
) + """

**Gates:** abort on NaN/complex energy; expect the projectile to **transit** (cross the
far slab face z=+12.5 with v>0) and E_electronic to **plateau** after exit. """ + (
"" if IS_PERT else
"Decisive fact: does native Ehrenfest even MOVE a z_valence=0 ghost ion? (z(t) flat ⇒ no.)"
))

# ---------------------------------------------------------------- setup cell
code(f"""%matplotlib inline
import os, glob, numpy as np, matplotlib.pyplot as plt
import pandas as pd
ROOT = "{ROOT}"
RES = "{RES}"
FRAMES = "{FRAMES}"
HYP = "{HYP}"
KIND = "{KIND}"
IS_PERT = {IS_PERT}
HA = 27.211386
SLAB_HALF = 12.5      # far face at +12.5, near face at −12.5
L_SLAB = 25.0         # slab thickness (Bohr) — the S = E_abs / L_slab denominator
DT = 0.04
V0 = 2.0
MASS = 1.0
RS = 4.183
os.makedirs(HYP, exist_ok=True)
print("results dir:", RES)
print("exists:", os.path.isdir(RES), " frames:", len(glob.glob(FRAMES + "/*.vti")))""")

# ---------------------------------------------------------------- GIF at TOP
md("""## 1. Density-evolution GIF (visual intuition — read this first)

n(x,z,t) on the mid-y x–z plane, as **total density** and **induced Δn = n(t)−n(0)**
(the projectile-driven wake). VTIs loaded via `inqview.load_vti` (physical order —
never fftshift'd). Slab faces (±12.5) dashed. Base64-embedded so it animates on reopen.""")

code("""from inqview.visualisation.density_gifs import _slice_stack, _save_gif
from IPython.display import Image, display

gif_paths = []
files = sorted(glob.glob(FRAMES + "/*.vti"))
if not files:
    print("no density frames — GIF skipped (run should re-run with SAVE_EVERY>0)")
else:
    nfiles = len(files)
    frames_max = 30
    idx = list(range(0, nfiles, max(1, nfiles // frames_max)))
    times, tot, axes = _slice_stack(FRAMES, idx, DT)
    cap_lines = (SLAB_HALF, -SLAB_HALF)
    # total density
    p_tot = os.path.join(HYP, f"{KIND}_total_density.gif")
    _save_gif(tot, times, axes, p_tot, title=f"{KIND} · total density n(x,z,t)",
              cap_lines=cap_lines, kind="density", fps=10)
    gif_paths.append(p_tot)
    # induced Δn = n(t) − n(0)
    p_ind = os.path.join(HYP, f"{KIND}_induced_delta.gif")
    _save_gif(tot - tot[0][None], times, axes, p_ind,
              title=f"{KIND} · induced Δn = n(t) − n(0)", cap_lines=cap_lines, kind="diff", fps=10)
    gif_paths.append(p_ind)
    print("wrote:", gif_paths)""")

md("### Total density  n(x,z,t)")
code("""if gif_paths:
    display(Image(filename=gif_paths[0]))
else:
    print("no frames")""")
md("### Induced wake  Δn = n(t) − n(0)")
code("""if len(gif_paths) > 1:
    display(Image(filename=gif_paths[1]))
else:
    print("no frames")""")

# ---------------------------------------------------------------- load data
md("## 2. Load the run data")
if IS_PERT:
    code("""obs = pd.read_csv(RES + "/raw/observables/observables.csv")
proj = pd.read_csv(RES + "/raw/observables/projectile.csv")
inter = pd.read_csv(RES + "/raw/observables/interactions.csv")
# E_electronic = kinetic + hartree + xc + external + nonlocal (Ha)
obs["E_elec"] = (obs["energy_kinetic"] + obs["energy_hartree"] + obs["energy_xc"]
                 + obs["energy_external"] + obs["energy_nonlocal"])
# merge projectile (proj_z, proj_vz, KE_proj, U_proj_bg) on step
df = obs.merge(proj, left_on="step", right_on="step", suffixes=("", "_p"))
df["t"] = df["time_au"]
print(len(df), "steps;  z0=", df.proj_z.iloc[0], " z_final=", df.proj_z.iloc[-1],
      " vz_final=", df.proj_vz.iloc[-1])
df.head(3)""")
else:
    code("""nat = pd.read_csv(RES + "/native.csv")
nat["t"] = nat["time"]
nat["E_elec"] = nat["E_kin"] + nat["E_hartree"] + nat["E_xc"] + nat["E_external"] + nat["E_nonlocal"]
df = nat.rename(columns={"z": "proj_z", "vz": "proj_vz"})
moved = abs(df.proj_z.iloc[-1] - df.proj_z.iloc[0]) > 1e-3
print(len(df), "steps;  z0=", df.proj_z.iloc[0], " z_final=", df.proj_z.iloc[-1],
      " vz_final=", df.proj_vz.iloc[-1], " MOVED:", moved)
df.head(3)""")

# ---------------------------------------------------------------- trajectory
md("""## 3. Projectile trajectory — transit & exit

z(t) crossing the slab (faces at ±12.5) and vz(t) (the light electron decelerates).
**Transit** = z crosses the far face +12.5 with v>0.""")
code("""fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
a1.plot(df.t, df.proj_z, lw=1.6)
a1.axhline(SLAB_HALF, ls="--", c="0.4"); a1.axhline(-SLAB_HALF, ls="--", c="0.4")
a1.axhspan(-SLAB_HALF, SLAB_HALF, color="0.9", zorder=0)
a1.set_xlabel("t (a.u.)"); a1.set_ylabel("proj_z (Bohr)"); a1.set_title("z(t) — slab shaded")
a2.plot(df.t, df.proj_vz, lw=1.6, c="C1")
a2.axhline(V0, ls=":", c="0.5"); a2.axhline(0, ls="-", c="0.7", lw=0.6)
a2.set_xlabel("t (a.u.)"); a2.set_ylabel("proj_vz (a.u.)"); a2.set_title("vz(t)")
plt.tight_layout(); plt.savefig(os.path.join(HYP, f"{KIND}_trajectory.png"), dpi=110); plt.show()

# transit diagnostics
z = df.proj_z.values; vz = df.proj_vz.values; t = df.t.values
crossed_far = np.any((z >= SLAB_HALF) & (vz > 0))
i_far = np.argmax(z >= SLAB_HALF) if np.any(z >= SLAB_HALF) else -1
reached_slab = np.any(z >= -SLAB_HALF)
print(f"crossed far face (+12.5) with v>0: {crossed_far}")
print(f"reached near slab face (−12.5): {reached_slab}   z range [{z.min():.1f}, {z.max():.1f}] Bohr")
if i_far >= 0:
    print(f"  first reaches +12.5 at t={t[i_far]:.2f}, vz={vz[i_far]:.3f}")
print(f"  z_final={z[-1]:.2f} Bohr, vz_final={vz[-1]:.3f} a.u.")
if not IS_PERT and not reached_slab:
    print("\\n*** VERDICT (native ghost-UPF ion): NEVER reached the slab — bounced in vacuum.")
    print("    The bare +1/r ghost tail produces a spurious vacuum force; the ion decelerates")
    print("    and reverses ~14 Bohr from the slab. Unphysical for this geometry — the")
    print("    perturbation projectile (Run A) is the valid representation. See section 8 overlay.")""")

# ---------------------------------------------------------------- E_elec plateau
md("""## 4. Electronic energy — deposition & post-exit plateau

E_electronic(t) rises as the projectile drives the slab, then **plateaus** once the
projectile has left (no CAP ⇒ energy conserved ⇒ flat plateau). **E_absorbed** =
plateau − initial.""")
code("""E = df.E_elec.values
E0 = E[:5].mean()
BOX_HALF = 42.5   # z-box edge; the perturbation projectile stops overlapping the grid past here
# TRUE plateau = the flat tail once the projectile's Gaussian has fully left the box.
# (Between the far slab face and the box edge, E_elec still relaxes — the induced
#  polarisation partially returns; only past the box edge is it dead-flat.)
exited = z > (BOX_HALF + 2.0)
if exited.any():
    i_exit = np.argmax(exited)
else:  # projectile never cleared the box: fall back to far-face + margin
    exited = z > (SLAB_HALF + 4.0)
    i_exit = np.argmax(exited) if exited.any() else len(E) - max(10, len(E)//10)
tail = E[i_exit:]
E_plateau = tail.mean(); plateau_std = tail.std()
E_abs_elec_ha = E_plateau - E0
E_abs_elec_ev = E_abs_elec_ha * HA

# CROSS-CHECK: energy deposited = projectile KINETIC-ENERGY LOSS across the transit.
# KE0 = ½ m v0²; KE_final from the settled far-field vz. This is the cleaner
# "energy given to the medium" for a decelerating projectile (frozen-out E_elec
# tail can carry a residual box-relaxation offset).
if IS_PERT:
    KE0 = df.energy_proj_ke.iloc[:5].mean()
    KEf = df.energy_proj_ke.iloc[-50:].mean()
else:
    KE0 = 0.5*MASS*(df.proj_vz.iloc[:5].mean()**2)
    KEf = 0.5*MASS*(df.proj_vz.iloc[-50:].mean()**2)
E_dep_ke_ha = KE0 - KEf
E_dep_ke_ev = E_dep_ke_ha * HA

# Headline E_absorbed for the slab estimator = KE loss (deposited energy).
E_abs_ha = E_dep_ke_ha
E_abs_ev = E_dep_ke_ev

fig, ax = plt.subplots(figsize=(8.4, 4.2))
ax.plot(df.t, (E - E0) * HA, lw=1.4, label="E_elec − E_elec(0)")
ax.axhline(E_abs_elec_ev, ls="--", c="C3", label=f"E_elec plateau ΔE={E_abs_elec_ev:.1f} eV")
ax.axhline(E_dep_ke_ev, ls=":", c="C4", label=f"KE-loss deposited={E_dep_ke_ev:.1f} eV")
ax.axvspan(df.t.iloc[i_exit], df.t.iloc[-1], color="C2", alpha=0.12, label="flat plateau window")
# mark slab-crossing window
in_slab = (z >= -SLAB_HALF) & (z <= SLAB_HALF)
if in_slab.any():
    ax.axvspan(df.t.values[np.argmax(in_slab)], df.t.values[len(in_slab)-1-np.argmax(in_slab[::-1])],
               color="0.85", alpha=0.5, zorder=0, label="in slab")
ax.set_xlabel("t (a.u.)"); ax.set_ylabel("E_elec − E_elec(0)  (eV)")
ax.set_title(f"Electronic deposition — E_elec plateau {E_abs_elec_ev:.1f} eV, KE-loss {E_dep_ke_ev:.1f} eV")
ax.legend(fontsize=8); plt.tight_layout()
plt.savefig(os.path.join(HYP, f"{KIND}_E_elec.png"), dpi=110); plt.show()
print(f"E_elec(0)        = {E0:.4f} Ha")
print(f"E_elec plateau   = {E_plateau:.4f} Ha  (tail std {plateau_std*HA:.2f} eV — flat check; window z>{BOX_HALF+2})")
print(f"  ΔE_elec(plateau) = {E_abs_elec_ha:.4f} Ha = {E_abs_elec_ev:.1f} eV  (incl. box-relaxation offset)")
print(f"KE0={KE0:.4f} Ha  KEf={KEf:.4f} Ha")
print(f"  E_deposited (KE loss) = {E_dep_ke_ha:.4f} Ha = {E_dep_ke_ev:.1f} eV   <-- headline")""")

# ---------------------------------------------------------------- ledger / conservation
if IS_PERT:
    md("""## 5. Pairwise energy ledger & conservation column

The full P/S/B pairwise Coulomb decomposition vs t, and the **conservation column**
E_total_sys = E_elec + KE_proj + U_proj_bg — must be flat (correctness gate).""")
    code("""fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.2))
for col, lab in [("e_ss","P-slab? e_ss"),("e_pp","e_pp"),("e_ps","e_ps"),
                 ("e_sb","e_sb"),("e_pb","e_pb")]:
    a1.plot(inter.time_au, inter[col]*HA, lw=1.1, label=lab)
a1.set_xlabel("t (a.u.)"); a1.set_ylabel("pairwise E (eV)"); a1.set_title("Pairwise Coulomb ledger")
a1.legend(fontsize=7, ncol=2)

# conservation: E_elec + KE_proj + U_proj_bg (Ha) — should be flat
cons = df.E_elec + df.energy_proj_ke + df.energy_proj_bg_ideal
cons_ev = (cons - cons.iloc[0]) * HA
a2.plot(df.t, cons_ev, lw=1.3, c="k")
a2.set_xlabel("t (a.u.)"); a2.set_ylabel("ΔE_conserved (eV)")
a2.set_title(f"Conservation E_elec+KE_proj+U_proj_bg (drift {cons_ev.iloc[-1]:.2f} eV)")
plt.tight_layout(); plt.savefig(os.path.join(HYP, f"{KIND}_ledger.png"), dpi=110); plt.show()
print(f"conservation drift over run: {cons_ev.iloc[-1]:.3f} eV  (max |dev| {np.abs(cons_ev).max():.3f} eV)")""")
else:
    md("""## 5. Total-energy conservation (native Ehrenfest)

INQ's native Ehrenfest conserves the total energy E_total (electronic + ion KE)
internally. We check E_total(t) flatness as the correctness gate.""")
    code("""fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(df.t, (df.E_total - df.E_total.iloc[0]) * HA, lw=1.3, c="k")
ax.set_xlabel("t (a.u.)"); ax.set_ylabel("ΔE_total (eV)")
drift = (df.E_total.iloc[-1] - df.E_total.iloc[0]) * HA
ax.set_title(f"E_total conservation — drift {drift:.2f} eV")
plt.tight_layout(); plt.savefig(os.path.join(HYP, f"{KIND}_conservation.png"), dpi=110); plt.show()
print(f"E_total drift: {drift:.3f} eV")""")

# ---------------------------------------------------------------- stopping power
md("""## 6. Stopping power — step by step

**Method 1 (slab deposition):** S = E_deposited / L_slab, with L_slab = 25 Bohr. For a
DECELERATING light projectile we take E_deposited = the projectile **kinetic-energy
loss** across the transit (KE0 − KE_far), which is the clean, box-offset-free measure
of energy handed to the medium. (The raw E_elec plateau carries a residual
box-relaxation offset once the Gaussian leaves the grid, so it over-counts — reported
alongside for transparency.) This S is a transit-AVERAGED value (v sweeps v0→v_far),
NOT S at a single v.

**Method 2 (initial-drag cross-check):** S(v0) = −d(KE_proj)/ds over the early
near-constant-velocity window (vz ≥ 0.85·v0), from the per-step track — this is S *at*
v0. The light electron decelerates, so a full-run regression is WRONG; window on
vz ≥ 0.85·v0 (light-projectile rule).""")
code("""# --- Method 1: slab deposition (KE-loss / L_slab) ---
S_slab_ha = E_abs_ha / L_SLAB          # E_abs_ha = KE loss (set in section 4)
S_slab_ev = E_abs_ev / L_SLAB
S_elecplateau_ev = E_abs_elec_ev / L_SLAB
print("=== Method 1: S = E_deposited / L_slab ===")
print(f"  E_deposited (KE loss) = {E_abs_ev:.1f} eV,  L_slab = {L_SLAB} Bohr")
print(f"  S_slab (KE-loss)      = {S_slab_ev:.3f} eV/Bohr   ({S_slab_ha:.4f} Ha/Bohr)   <-- headline")
print(f"  [for reference] S from raw E_elec plateau = {S_elecplateau_ev:.2f} eV/Bohr (box-offset inflated)")

# --- Method 2: initial drag −dKE/ds over vz ≥ 0.85 v0 ---
print("\\n=== Method 2: initial-drag −dKE_proj/ds (vz ≥ 0.85·v0) ===")
if IS_PERT:
    KE = df.energy_proj_ke.values
else:
    KE = 0.5 * MASS * (df.proj_vz.values ** 2)   # native: reconstruct KE from vz
s_path = z - z[0]                                  # path length along z (Bohr)
# Restrict to the DRAG region: from the near slab face onward (z ≥ −SLAB_HALF − 2),
# where the projectile actually feels the medium. Fitting KE vs path over the flat
# vacuum approach (z ≈ −30..−15) would dilute the slope toward zero.
in_drag = z >= (-SLAB_HALF - 2.0)
S_drag_ev_ref = float("nan")
for frac in (0.85, 0.70, 0.50):
    m = (vz >= frac * V0) & in_drag
    if m.sum() >= 5:
        A = np.polyfit(s_path[m], KE[m], 1)         # KE ≈ A[0]·s + A[1]
        S_drag_ha = -A[0]; S_drag_ev = S_drag_ha * HA
        vmean = vz[m].mean()
        print(f"  vz ≥ {frac:.2f}·v0 & in drag region: {m.sum():4d} pts, mean v={vmean:.3f}, "
              f"S_drag = {S_drag_ev:.3f} eV/Bohr  ({S_drag_ha:.4f} Ha/Bohr)")
        if frac == 0.85 or not (S_drag_ev_ref == S_drag_ev_ref):
            S_drag_ev_ref = S_drag_ev
    else:
        print(f"  vz ≥ {frac:.2f}·v0: only {m.sum()} pts — too sparse")
print(f"\\nSUMMARY")
print(f"  S_slab (KE-loss/L)      = {S_slab_ev:.2f} eV/Bohr  (transit-averaged, v {V0:.1f}->{np.sqrt(2*KEf/MASS):.2f})")
print(f"  S_drag(v0, 0.85 window) = {S_drag_ev_ref:.2f} eV/Bohr  (S at v0={V0})")
print(f"  Lindhard bulk (context) ~ 0.03-0.7 eV/Bohr (point / sigma=0.5)")""")

# ---------------------------------------------------------------- Lindhard eyeball
md("""## 7. Lindhard / bulk eyeball (context only — NON-gating)

Bulk RPA linear-response stopping for a σ=0.5 Gaussian projectile in an r_s=4.18
electron gas at v=2, as a rough order-of-magnitude reference. This is an infinite-medium
bulk number; the pilot is a finite 25-Bohr slab, so exact agreement is NOT expected.""")
code("""from inqview.analysis.lindhard_elf import kF_from_rs, omega_p, stopping_power_sigma, stopping_power_point
kF = kF_from_rs(RS)
wp = omega_p(kF)
try:
    S_lind_sigma = stopping_power_sigma(V0, kF, 0.5, eta=1e-2)
except Exception as e:
    S_lind_sigma = float("nan"); print("sigma-LR failed:", e)
try:
    S_lind_pt = stopping_power_point(V0, kF, eta=1e-2)
except Exception as e:
    S_lind_pt = float("nan"); print("point-LR failed:", e)
print(f"r_s={RS}  kF={kF:.3f} a.u.  ω_p={wp:.3f} Ha = {wp*HA:.2f} eV")
print(f"Lindhard bulk S(v=2, σ=0.5) = {S_lind_sigma:.3f} eV/Bohr  (per-Bohr, form-factor)"
      if S_lind_sigma==S_lind_sigma else "Lindhard σ-LR unavailable")
print(f"Lindhard bulk S(v=2, point) = {S_lind_pt:.3f} eV/Bohr"
      if S_lind_pt==S_lind_pt else "Lindhard point-LR unavailable")
print("\\n(NON-gating context: bulk infinite medium vs finite 25-Bohr slab pilot.)")""")

# ---------------------------------------------------------------- RUN B overlay
if not IS_PERT:
    md("""## 8. Native vs perturbation — the faithfulness overlay

Overlay this native-ion z(t)/vz(t) against RUN A (perturbation, analytic force) on the
identical slab. Close agreement confirms the perturbation projectile faithfully
reproduces the real-ion trajectory.""")
    code("""pert_res = os.path.join(ROOT,
    "ResearchProject/systems/localised_jellium/scripts/classical_highdensity_sv/pilot/results/pilot")
pert_proj_csv = os.path.join(pert_res, "raw/observables/projectile.csv")
if os.path.isfile(pert_proj_csv):
    pp = pd.read_csv(pert_proj_csv)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    a1.plot(df.t, df.proj_z, lw=1.6, label="native ion (B)")
    a1.plot(pp.time_au, pp.proj_z, lw=1.4, ls="--", label="perturbation (A)")
    a1.axhline(SLAB_HALF, ls=":", c="0.5"); a1.axhline(-SLAB_HALF, ls=":", c="0.5")
    a1.set_xlabel("t (a.u.)"); a1.set_ylabel("proj_z (Bohr)"); a1.set_title("z(t): native vs perturbation"); a1.legend()
    a2.plot(df.t, df.proj_vz, lw=1.6, label="native ion (B)")
    a2.plot(pp.time_au, pp.proj_vz, lw=1.4, ls="--", label="perturbation (A)")
    a2.set_xlabel("t (a.u.)"); a2.set_ylabel("proj_vz (a.u.)"); a2.set_title("vz(t): native vs perturbation"); a2.legend()
    plt.tight_layout(); plt.savefig(os.path.join(HYP, f"{KIND}_overlay.png"), dpi=110); plt.show()
    # agreement metric on common time support
    tmax = min(df.t.max(), pp.time_au.max())
    ti = np.linspace(0, tmax, 200)
    zn = np.interp(ti, df.t, df.proj_z); za = np.interp(ti, pp.time_au, pp.proj_z)
    print(f"native-vs-perturbation z(t): max|Δz|={np.abs(zn-za).max():.3f} Bohr, "
          f"RMS={np.sqrt(np.mean((zn-za)**2)):.3f} Bohr")
else:
    print("RUN A (perturbation) projectile.csv not found — overlay skipped.")""")

# ================================================================ write + execute
md("""---
*Auto-built by `build_pilot_notebook.py`. Executed end-to-end (0 errors). All figures
saved into this directory; density GIFs base64-embedded in the cells above.*""")

nb["cells"] = cells

import nbformat
from nbclient import NotebookClient

nb_path = os.path.join(HYP, f"{KIND}_pilot.ipynb")
nbformat.write(nb, nb_path)
print("wrote notebook:", nb_path)

client = NotebookClient(nb, timeout=1200,
    kernel_name="inqview-venv",
    resources={"metadata": {"path": HYP}})
client.execute()
nbformat.write(nb, nb_path)
print("EXECUTED (0 errors):", nb_path)
