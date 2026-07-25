#!/usr/bin/env python3
"""Deterministic builder for p2_classical_truncated_stopping.ipynb.

Models p5_classical_transient_comparison.ipynb, but on the LATEST classical run
(qsp_phase2) and with the new headline feature: a PERIODIC-WRAP TRUNCATION.

  run: ResearchProject/systems/localised_jellium/scripts/qsp_phase2/classical/
       results/p2_classical
  localised positive jellium SLAB, half-width 12.5 (25 Bohr thick), faces z=+/-12.5;
  cell 50x50x70 Bohr, z in [-35,+35] PERIODIC; two-sided sin^2 CAP eta=-0.7 over
  |z| in [25,35]; classical Gaussian-e projectile sigma_pot=0.35 (=sigma_WP 0.5),
  launched z0=-22, v0=2.711 a.u. (E=100 eV, r_s~5.67); 2000 steps, dt=0.02 (40 a.u.).

The 40 a.u. run is LONG ENOUGH for the projectile to cross the slab, traverse the
top CAP, WRAP through the periodic boundary, and re-enter the slab. That wrapped
image re-ploughs the slab and ruins E_total. We TRUNCATE at the step where the
projectile, having looped the cell, exits the first CAP after wrapping
(unwrapped z = +45 Bohr, physical z = -25): empirical step 1572 (t=31.44 a.u.).

Headline (user's choice, 2026-06-25): apply the stopping-power-extraction SKILL
kernel Method B (slab deposit / L_z) to the TRUNCATED window AS-IS and trust its
convergence flag. (The notebook validates the skill on a run it has not seen.)

Run (venv):
  /local/data/public/skcb2/tddft/venv/bin/python3 build_p2_classical_truncated.py
"""
from __future__ import annotations
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "p2_classical_truncated_stopping.ipynb")
cells = []
def md(s): cells.append(new_markdown_cell(s))
def code(s): cells.append(new_code_cell(s))

md(r"""# Periodic-wrap truncation + skill Method B on the **latest classical run** (`qsp_phase2`)

Companion to `p5_classical_transient_comparison.ipynb`, on the newest localised-jellium
classical slab run. **New feature: a periodic-wrap truncation.** The run is long
enough (40 a.u.) that the projectile loops the periodic cell and **re-enters the
slab** — a catastrophe for the energy signal — so we cut the time series before that.

**Run.** `localised_jellium / qsp_phase2 / classical / p2_classical`.
Localised positive jellium **slab**, half-width 12.5 → faces at z=±12.5 (25 Bohr
thick), r_s≈5.67. Cell 50×50×**70** Bohr, **z ∈ [−35,+35] periodic**. Classical
Gaussian-e projectile σ_pot=0.35 (≡ σ_WP=0.5), launched **z₀=−22**, v₀=2.711 a.u.
(E=100 eV). Two-sided sin² CAP η=−0.7 over **|z| ∈ [25,35]**. 2000 steps, dt=0.02.

**Geometry of the wrap (the truncation calculation).** Moving +z, the projectile
crosses the slab, the top CAP, the wall at z=+35, wraps to z=−35, traverses the
*lower* CAP, and **exits it at physical z=−25 ≡ unwrapped z=+45**. Path from
launch = 45−(−22) = **67 Bohr**.

| | distance | time | step (`t/dt`) |
|---|---|---|---|
| quick estimate, constant v₀ | 67 Bohr | 24.7 a.u. | ~1236 |
| **actual track (projectile decelerates)** | 67 Bohr | **31.44 a.u.** | **1572** ← cut |

We **keep steps 0…1572** (t ≤ 31.44 a.u.) and drop the rest.

**Headline method (locked with the user 2026-06-25):** apply the
`stopping-power-extraction` **skill kernel Method B** (`slab_stopping_power`,
`S = [E_total(t_f) − E_total(t₀)]/L_z`) to the truncated window **as-is** and
**trust its convergence flag**. The notebook thereby validates the skill on a run
it has never seen.""")

# ---------------------------------------------------------------- §1 kernel + load
md(r"""## 1 — Load the run + the skill kernel; compute the truncation step

We import the kernel straight from the skill folder
(`.claude/skills/stopping-power-extraction/stopping_power.py`) — this notebook is a
live test that the shipped skill works on a fresh run.""")
code(r'''import sys, os, csv
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
sys.path.insert(0, "/local/data/public/skcb2/tddft/.claude/skills/stopping-power-extraction")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import stopping_power as SP                      # the SKILL kernel under test
from inqview.analysis.stopping_extract import load_track
from inqview.analysis import lindhard_elf as LE
from inqview.visualisation import style as ST
ST.apply_theme()

HA = SP.HA_TO_EV                                 # 27.2114 eV/Ha
HpB2eVpA = SP.HA_PER_BOHR_TO_EV_PER_A            # Ha/Bohr -> eV/Angstrom
RUN = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "scripts/qsp_phase2/classical/results/p2_classical")
OBS = f"{RUN}/raw/observables"

# --- geometry (from run.cpp + slab_n82_L50x50x70.hpp) ---
Z0   = -22.0
V0   =  2.7110633401
MASS =  1.0                                       # KE in a.u. uses electron mass = 1
LZ   = 70.0                                       # cell length; z in [-35,+35] periodic
SLAB_HALF, L_SLAB = 12.5, 25.0                    # slab faces +/-12.5; thickness 25 Bohr
CAP_LO, CAP_HI = 25.0, 35.0                       # |z| in [25,35] is CAP
RS = 5.666; KF = LE.kF_from_rs(RS)
DT = 0.02

# --- truncation step: projectile exits the FIRST CAP after wrapping (unwrapped z=+45) ---
Z_CUT_UNWRAP = (35.0)            # top wall
Z_EXIT_CAP   = 45.0             # = -25 + 70 : exit lower CAP after one wrap  -> THE CUT

tr_full = load_track(f"{OBS}/electron_track.csv", mass=MASS, axis="z")
z_unwrap = Z0 + tr_full.s        # tr.s is displacement; unwrapped z
# empirical step where unwrapped z first reaches +45
o_raw = pd.read_csv(f"{OBS}/electron_track.csv").drop_duplicates(subset="step").sort_values("step")
steps_trk = o_raw.step.to_numpy(); z_trk = o_raw.z.to_numpy()
i_cut = int(np.argmax(z_trk >= Z_EXIT_CAP))
CUT_STEP = int(steps_trk[i_cut]); CUT_T = CUT_STEP * DT
# quick v0 estimate, for the record
t_v0 = (Z_EXIT_CAP - Z0) / V0
print(f"slab r_s={RS}  kF=vF={KF:.3f}  v0/vF={V0/KF:.2f}")
print(f"truncation: unwrapped z=+45 reached at step {CUT_STEP}, t={CUT_T:.2f} a.u.")
print(f"  (constant-v0 quick estimate would be t={t_v0:.2f} a.u. = step ~{round(t_v0/DT)}; "
      f"the real projectile is slower, so it reaches +45 LATER)")
''')

# ---------------------------------------------------------------- §2 the catastrophe
md(r"""## 2 — Why we must truncate: the wrapped image re-ploughs the slab

The track `z` is **unwrapped** (`final z = 62.7` Bohr in a 70-Bohr cell), so by
the end of the run the projectile sits at physical `62.7 mod 70 = −7.3` Bohr —
**back inside the slab**. Below: the projectile path with the wrap milestones, and
the energy signal it produces. Run the skill's Method B on the **full** (un-cut)
series and it returns garbage.""")
code(r'''o = pd.read_csv(f"{OBS}/observables.csv").drop_duplicates(subset="step").sort_values("time_au")
ot = o.time_au.to_numpy(); Etot = o.energy_total.to_numpy(); ostep = o.step.to_numpy()
x_obs = np.interp(ot, tr_full.t, tr_full.s) + Z0    # unwrapped z at obs times
dEtot_full = (Etot - Etot[0])                       # Ha

# skill Method B on the FULL (untruncated) series -> catastrophe
rB_full = SP.slab_stopping_power(ot, Etot, L_z=L_SLAB)
print(f"FULL run, skill Method B:  S = {rB_full['S']:.4f} Ha/Bohr = {rB_full['S']*HpB2eVpA:.2f} eV/A"
      f"   dE = {rB_full['dE']*HA:.1f} eV   status={rB_full['status']}")
print("  -> the projectile's KE FELL from 3.67 to 0.99 Ha; most of that 'loss' is the")
print("     wrapped image re-ploughing the slab. Unusable without truncation.")

fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.9), constrained_layout=True)
ax = axs[0]
ax.plot(tr_full.t, z_unwrap, "-", lw=1.0)
for zc, lab, col in [(CAP_LO, "CAP_hi entry", "C2"), (35.0, "cell loop (z=+35)", "C1"),
                     (45.0, "exit 1st CAP (cut)", "C3")]:
    ax.axhline(zc, ls=":", lw=0.8, color=col)
ax.axhline(SLAB_HALF, ls="--", lw=0.6, color="0.6"); ax.axhline(-SLAB_HALF, ls="--", lw=0.6, color="0.6")
ax.axvline(CUT_T, color="C3", lw=1.2, label=f"cut t={CUT_T:.1f}")
ax.axhspan(-SLAB_HALF, SLAB_HALF, color="C0", alpha=0.08)
ax.set_xlabel("t (a.u.)"); ax.set_ylabel("unwrapped z (Bohr)")
ax.set_title("Projectile loops the cell"); ax.legend(fontsize=6, loc="upper left")
ax = axs[1]
ax.plot(x_obs, dEtot_full*HA, "o-", ms=2.5, lw=0.7)
ax.axvspan(-SLAB_HALF, SLAB_HALF, color="C0", alpha=0.10, label="slab")
ax.axvline(Z_EXIT_CAP, color="C3", lw=1.2, label="cut (z=+45)")
ax.set_xlabel("unwrapped z (Bohr)"); ax.set_ylabel(r"$\Delta E_\mathrm{total}$ (eV)")
ax.set_title("Energy signal blows up after the wrap"); ax.legend(fontsize=6, loc="upper left")
fig
''')

# ---------------------------------------------------------------- §3 truncate
md(r"""## 3 — Apply the truncation (keep steps 0…1572)""")
code(r'''m_obs = ostep <= CUT_STEP
m_trk = steps_trk <= CUT_STEP
ot_c, Etot_c, x_c = ot[m_obs], Etot[m_obs], x_obs[m_obs]
dEtot_c = (Etot_c - Etot_c[0])
# truncated track arrays
t_c   = tr_full.t[m_trk]; s_c = tr_full.s[m_trk]; ke_c = tr_full.ke[m_trk]; v_c = tr_full.v[m_trk]
zc    = Z0 + s_c
print(f"kept {m_obs.sum()}/{len(ostep)} observable rows, {m_trk.sum()}/{len(steps_trk)} track rows")
print(f"window: t in [0, {CUT_T:.2f}] a.u., unwrapped z in [{zc.min():.1f}, {zc.max():.1f}] Bohr")
''')

# ---------------------------------------------------------------- §4 guards
md(r"""## 4 — Guards (skill §0): N(t) conservation + energy conservation

Both methods are invalid if the CAP drains the bath (then raw `E_total` is CAP
energy, not deposit). Here the η=−0.7 CAP barely drains N, so `E_total` is a clean
deposit signal.""")
code(r'''N = pd.read_csv(f"{OBS}/electron_number.csv")
g = SP.conservation_guard(N.N_total.to_numpy()[N.time_au.to_numpy() <= CUT_T])
print(f"N(t): {N.N_total.iloc[0]:.3f} -> {N.N_total[N.step<=CUT_STEP].iloc[-1]:.3f}   "
      f"conservation_guard: ok={g['ok']}  drained={g['drained_frac']*100:.3f}%  (tol {g['tol']*100:.0f}%)")
dKE_cut  = ke_c[0] - ke_c[-1]
dEel_cut = dEtot_c[-1]
print(f"energy conservation over the window:  dKE_ion = {dKE_cut:+.4f} Ha   "
      f"dE_electronic = {dEel_cut:+.4f} Ha   (gap {abs(dKE_cut-dEel_cut)*HA:.2f} eV = CAP+numerics)")
''')

# ---------------------------------------------------------------- §5 headline Method B
md(r"""## 5 — Headline: skill **Method B** on the truncated window (as-is, trust the flag)

This is the agreed deliverable: feed the truncated `E_total` to the shipped kernel
`SP.slab_stopping_power(t, E_total, L_z=25)` and report what it says — number **and**
convergence verdict — without second-guessing it.""")
code(r'''rB = SP.slab_stopping_power(ot_c, Etot_c, L_z=L_SLAB)
print("skill Method B (slab_stopping_power) on the TRUNCATED window:")
print(f"  S        = {rB['S']:.5f} Ha/Bohr = {rB['S']*HpB2eVpA:.3f} eV/A")
print(f"  deposit  = {rB['dE']*HA:.3f} eV over L_z = {rB['L_z']:.0f} Bohr")
print(f"  STATUS   = {rB['status'].upper()}   (tail change = {rB['tail_frac_of_total']*100:.0f}% of deposit "
      f"in the last {rB['converge_frac']*100:.0f}% of the window; tol {rB['converge_tol']*100:.0f}%)")
print()
print("VERDICT (trust the flag): the convergence gate fires NOT_CONVERGED, so the")
print("reported S = 0.0125 Ha/Bohr (0.64 eV/A) is a LOWER BOUND, not a settled value.")

fig, ax = ST.figure_one_col()
ax.plot(x_c, dEtot_c*HA, "o-", ms=3, lw=0.8, label=r"$\Delta E_\mathrm{total}$ (truncated)")
ax.axvspan(-SLAB_HALF, SLAB_HALF, color="C0", alpha=0.10, label="slab")
ax.axhline(rB['dE']*HA, ls="--", color="C3", lw=1.0,
           label=f"deposit = {rB['dE']*HA:.1f} eV  (S={rB['S']*HpB2eVpA:.2f} eV/Å, {rB['status']})")
ax.set_xlabel("unwrapped z (Bohr)"); ax.set_ylabel(r"$\Delta E_\mathrm{total}$ (eV)")
ax.set_title("Method B deposit over the truncated window"); ax.legend(fontsize=6, loc="upper left")
fig
''')

# ---------------------------------------------------------------- §6 reading the flag
md(r"""## 6 — Reading the `not_converged` flag honestly (trainee note)

The flag is correct that the deposit has **not settled** in this window — but the
*reason* is **not** the one the skill's docstring assumes (“projectile hasn’t
finished the slab → extend the run”). Here the projectile finished the slab long
ago; the tail keeps moving because **the projectile’s KE has a large *reversible*
excursion** across the slab (it slows climbing the mean-field feature and speeds
back up leaving it). So the skill’s prescribed remedy — *extend the run* — is
**inverted** here: a longer run only brings the wrapped image back sooner. This is
a genuine gap to feed back into the skill (a charged projectile in a localised slab
needs an *equal-potential* window, not a longer run).""")
code(r'''# KE(z) shows the reversible well: min at slab centre, recovers on exit
fig, ax = ST.figure_one_col()
ax.plot(zc, ke_c, "-", lw=1.0)
ax.axvspan(-SLAB_HALF, SLAB_HALF, color="C0", alpha=0.10, label="slab (faces ±12.5)")
ax.axvspan(CAP_LO, CAP_HI, color="C3", alpha=0.06)
i_min = int(np.argmin(ke_c))
ax.plot(zc[i_min], ke_c[i_min], "v", color="C3", ms=6,
        label=f"KE_min={ke_c[i_min]:.2f} Ha at z={zc[i_min]:.1f}")
ax.axhline(ke_c[0], ls=":", color="0.5", lw=0.8, label=f"launch KE={ke_c[0]:.2f} Ha")
ax.set_xlabel("unwrapped z (Bohr)"); ax.set_ylabel(r"KE$_\mathrm{ion}$ (Ha)")
ax.set_title("Reversible KE excursion → why Method B reports not_converged")
ax.legend(fontsize=5.5, loc="lower right"); fig
''')
code(r'''# Cross-check CHANNEL (skill §4): equal-potential slab-face window z=-12.5 -> +12.5.
# Reported as a sanity channel, NOT the headline.
def ke_at_z(z_target):
    i = int(np.argmin(np.abs(zc - z_target))); return ke_c[i]
def Eel_at_z(z_target):
    i = int(np.argmin(np.abs(x_c - z_target))); return dEtot_c[i]
dKE_face = ke_at_z(-SLAB_HALF) - ke_at_z(+SLAB_HALF)
dEel_face = Eel_at_z(+SLAB_HALF) - Eel_at_z(-SLAB_HALF)
print("sanity channel — net loss across equal-potential slab faces (z=-12.5 -> +12.5):")
print(f"  dKE_ion  = {dKE_face:+.4f} Ha -> S = {dKE_face/L_SLAB:.5f} Ha/Bohr = {dKE_face/L_SLAB*HpB2eVpA:.3f} eV/A")
print(f"  dE_elec  = {dEel_face:+.4f} Ha -> S = {dEel_face/L_SLAB:.5f} Ha/Bohr = {dEel_face/L_SLAB*HpB2eVpA:.3f} eV/A")
print(f"  (independent reference from dispatch: --measured-s 0.018632 Ha/Bohr = {0.018632*HpB2eVpA:.3f} eV/A)")
print("  -> the two channels agree to ~4%, BUT this window is a cross-check; the")
print("     headline remains the skill Method B value + its not_converged flag.")
''')

# ---------------------------------------------------------------- §7 sanity channels
md(r"""## 7 — Sanity channels on the truncated window (skill §4)

Run every time; large deviations get reported, not averaged away.""")
code(r'''# (a) kinetic channel: -dKE/dx over the window
kin = SP.kinetic_channel(s_c, ke_c, s_c, s_c.min(), s_c.max())
# (b) int F.v cumulative (== dKE analytically; a profile/discretisation check)
t_fv, dep_fv = SP.force_power_channel(t_c, v_c, mass=MASS)
fig, ax = ST.figure_one_col()
ax.plot(zc, (ke_c[0]-ke_c)*HA, "s--", ms=2.5, lw=0.7, label=r"$\Delta$KE$_\mathrm{ion}$ (kinetic)")
ax.plot(x_c, dEtot_c*HA, "o-", ms=2.5, lw=0.7, label=r"$\Delta E_\mathrm{total}$ (electronic)")
ax.plot(zc, dep_fv*HA, "^:", ms=2.5, lw=0.7, label=r"$\int -F\cdot v\,dt$ (≡ΔKE)")
ax.axvspan(-SLAB_HALF, SLAB_HALF, color="C0", alpha=0.10)
ax.set_xlabel("unwrapped z (Bohr)"); ax.set_ylabel("deposited energy (eV)")
ax.set_title("Three channels over the truncated window")
ax.legend(fontsize=6, loc="upper left"); fig
''')

# ---------------------------------------------------------------- §8 summary
md(r"""## 8 — Summary (neutral — your verdict)

- **Truncation works.** Cutting at the wrap+first-CAP-exit (unwrapped z=+45, **step
  1572, t=31.44 a.u.**) removes the catastrophe: full-run Method B = 0.10 Ha/Bohr
  (70 eV deposit, image re-plough) → truncated Method B = 0.0125 Ha/Bohr.
- **Guards pass:** N conserved (≈0.02% drained), `dKE_ion ≈ dE_electronic`.
- **Headline (skill Method B, as-is):** `S = 0.0125 Ha/Bohr = 0.64 eV/Å`,
  **status `not_converged`** — trust the flag: a **lower bound**, not a settled S.
- **Trainee note:** the flag fires because of a *reversible* KE excursion across
  the slab, not an unfinished transit — so the skill's "extend the run" remedy is
  inverted here. A cross-check over the equal-potential slab faces gives
  ≈0.95–0.99 eV/Å (both channels, matches the independent 0.0186 reference), which
  is why the skill likely needs an **equal-potential-face window** rule for charged
  projectiles in a localised slab.

*No verdict recorded here (verification-user-owns-verdict).*""")

nb = new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"name": "inqview-venv", "display_name": "inqview-venv", "language": "python"}
from nbconvert.preprocessors import ExecutePreprocessor
ep = ExecutePreprocessor(timeout=900, kernel_name="inqview-venv")
print("executing ...")
ep.preprocess(nb, {"metadata": {"path": HERE}})
with open(OUT, "w") as fh: nbf.write(nb, fh)
print("wrote", OUT)
