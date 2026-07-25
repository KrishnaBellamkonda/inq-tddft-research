#!/usr/bin/env python3
"""Classical stopping baseline — S extraction + headline comparison figure.

Campaign: docs/campaigns/localised_jellium/classical-stopping-baseline.md
Two matched localised-slab classical runs (twin of WP p5_wp_v1p3, sigma_WP=0.5, v=1.3):
  P1 (Ehrenfest, light electron)  -> S(v0) = initial-drag slope -dKE/ds over the
       in-slab v>=0.85*v0 window (light-projectile-stopping rule; E_total is
       contaminated by the moving-charge coupling term so is NOT the clean deposit).
  P2 (prescribed const-v)         -> S = dE_deposited/L_slab (stopping-power-extraction
       Method B); deposit isolated by subtracting the projectile-coupling (e_ps+e_pb)
       from E_total, read at the post-slab plateau.

References overlaid: S_WP(p5_wp_v1p3)=2.37 eV/Bohr [UPPER BOUND], bulk sigma=0.5
classical=0.94, point-charge Lindhard, sigma-Lindhard.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/local/data/public/skcb2/tddft"
sys.path.insert(0, f"{ROOT}/inq-stack/python")
sys.path.insert(0, f"{ROOT}/.claude/skills/stopping-power-extraction")
from inqview.analysis import lindhard_elf as L

RES = f"{ROOT}/ResearchProject/systems/localised_jellium/scripts/classical_slab_stopping/results"
OUT = f"{ROOT}/ResearchProject/systems/localised_jellium/hypotheses/classical_slab_stopping"
HA = 27.211386
V0 = 1.3
HALF = 12.5
L_SLAB = 2 * HALF          # 25 Bohr traversal length
N0 = 82 / (50 * 50 * L_SLAB)
RS = (3 / (4 * np.pi * N0)) ** (1 / 3)

def load(run):
    d = f"{RES}/{run}/raw/observables"
    p = np.genfromtxt(f"{d}/projectile.csv", delimiter=",", names=True)
    o = np.genfromtxt(f"{d}/observables.csv", delimiter=",", names=True)
    ix = np.genfromtxt(f"{d}/interactions.csv", delimiter=",", names=True)
    return p, o, ix

def linfit(x, y):
    """free-intercept least squares; returns slope, intercept, slope_stderr."""
    A = np.vstack([x, np.ones_like(x)]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    slope, icpt = coef
    yhat = A @ coef
    dof = max(len(x) - 2, 1)
    s2 = np.sum((y - yhat) ** 2) / dof
    sxx = np.sum((x - x.mean()) ** 2)
    se = np.sqrt(s2 / sxx) if sxx > 0 else np.nan
    return slope, icpt, se

# ---------- P1: Ehrenfest initial-drag slope ----------
p1, o1, ix1 = load("p1_ehrenfest_v1p3")
z1, vz1, ke1 = p1["proj_z"], p1["proj_vz"], p1["energy_proj_ke"] * HA   # eV
# in-slab AND v>=0.85*v0 window (deposition only starts inside the slab)
vth = 0.85 * V0
mask1 = (z1 >= -HALF) & (vz1 >= vth)
s_path = z1[mask1] - (-HALF)                    # path INTO the slab (Bohr)
ke_win = ke1[mask1]
# S(v0) = -dKE/ds  (light-projectile-stopping rule; energy the projectile loses per Bohr)
slope_ke, icpt_ke, se_ke = linfit(s_path, ke_win)
S_P1 = -slope_ke
S_P1_err = se_ke
mean_v_P1 = vz1[mask1].mean()
n_P1 = int(mask1.sum())
# N-conservation guard (slab norm should be ~const)
nslab1 = ix1["norm_slab"]
Ndrift1 = (nslab1.max() - nslab1.min()) / nslab1[0]

# ---------- P2: const-v deposit / L_slab ----------
p2, o2, ix2 = load("p2_constv_v1p3")
z2 = p2["proj_z"]
Etot2 = o2["energy_total"] * HA
# projectile-coupling terms (Hartree -> eV); subtract to isolate slab excitation
coupl2 = (ix2["e_ps"] + ix2["e_pb"]) * HA
deposit2 = (Etot2 - coupl2) - (Etot2[0] - coupl2[0])   # eV, relative to t=0
# plateau: projectile past far slab face (z>=+HALF), average the settled tail
past = z2 >= HALF
dep_plateau = deposit2[past]
# converged deposit = mean over the final settled portion (last 40% of the past-slab window)
tail = dep_plateau[int(0.6 * len(dep_plateau)):]
E_dep_P2 = tail.mean()
E_dep_P2_err = tail.std()
S_P2 = E_dep_P2 / L_SLAB
S_P2_err = E_dep_P2_err / L_SLAB
# raw E_total/L_slab (the note's naive metric) for comparison
raw_dep_P2 = (Etot2[past] - Etot2[0])
S_P2_raw = raw_dep_P2[int(0.6*len(raw_dep_P2)):].mean() / L_SLAB
nslab2 = ix2["norm_slab"]; Ndrift2 = (nslab2.max()-nslab2.min())/nslab2[0]

# ---------- references ----------
S_WP = 2.374            # p5_wp_v1p3 [UPPER BOUND]
S_bulk = 0.937          # classical sigma=0.5 bulk, v0=1.3
S_lind_pt = L.stopping_power_point(V0, RS)
try:
    S_lind_sig = L.stopping_power_sigma(V0, RS, 0.5 / np.sqrt(2))
except Exception:
    S_lind_sig = np.nan

print(f"r_s={RS:.3f}  L_slab={L_SLAB} Bohr")
print(f"[P1 Ehrenfest] S(v0={V0}) = {S_P1:.3f} +/- {S_P1_err:.3f} eV/Bohr "
      f"(initial-drag -dKE/ds, in-slab v>=1.1, n={n_P1}, mean_v={mean_v_P1:.2f})  "
      f"N_drift={Ndrift1*100:.2f}%")
print(f"[P2 const-v ] S = {S_P2:.3f} +/- {S_P2_err:.3f} eV/Bohr "
      f"(deposit/L_slab, coupling-subtracted; raw E_total/L={S_P2_raw:.3f})  N_drift={Ndrift2*100:.2f}%")
print(f"[refs] S_WP={S_WP} [UB]  bulk_sigma0.5={S_bulk}  Lindhard_point={S_lind_pt:.3f}  Lindhard_sigma={S_lind_sig:.3f}")

# ---------- figures ----------
# (a) P1 initial-drag fit
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(s_path, ke_win, ".", ms=3, label="KE_proj (in-slab, v>=1.1)")
ax[0].plot(s_path, icpt_ke + slope_ke * s_path, "r-",
           label=f"fit: S(v0)={S_P1:.2f}±{S_P1_err:.2f} eV/Bohr")
ax[0].set_xlabel("path into slab s (Bohr)"); ax[0].set_ylabel("KE_proj (eV)")
ax[0].set_title("P1 Ehrenfest — initial-drag slope"); ax[0].legend(fontsize=8)
# (b) P2 deposit plateau
ax[1].plot(z2, deposit2, "-", label="deposit = E_total − (e_ps+e_pb)")
ax[1].axhline(E_dep_P2, color="r", ls="--", label=f"plateau {E_dep_P2:.1f} eV → S={S_P2:.2f}")
ax[1].axvline(-HALF, color="k", ls=":"); ax[1].axvline(HALF, color="k", ls=":")
ax[1].set_xlabel("proj_z (Bohr)"); ax[1].set_ylabel("deposited energy (eV)")
ax[1].set_title("P2 const-v — deposit/L_slab (slab dotted)"); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{OUT}/stopping_extraction.png", dpi=130); plt.close(fig)

# (c) headline comparison bar
fig, ax = plt.subplots(figsize=(7, 4))
labels = ["S_classical\nP1 Ehrenfest\n(initial-drag)", "S_classical\nP2 const-v\n(ΔE/L)",
          "S_WP\np5_wp_v1p3\n[UPPER BOUND]", "bulk σ=0.5\nclassical", "Lindhard\npoint"]
vals = [S_P1, S_P2, S_WP, S_bulk, S_lind_pt]
errs = [S_P1_err, S_P2_err, 0, 0, 0]
colors = ["#1f77b4", "#2ca02c", "#d62728", "#7f7f7f", "#9467bd"]
ax.bar(range(len(vals)), vals, yerr=errs, color=colors, capsize=4)
for i, v in enumerate(vals):
    ax.text(i, v + 0.05, f"{v:.2f}", ha="center", fontsize=9)
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("S (eV/Bohr)")
ax.set_title(f"Localised-slab classical stopping baseline @ v=1.3 (σ_WP=0.5, r_s={RS:.1f})")
fig.tight_layout(); fig.savefig(f"{OUT}/S_comparison.png", dpi=130); plt.close(fig)
print(f"wrote {OUT}/stopping_extraction.png and S_comparison.png")

# summary CSV
import csv
with open(f"{OUT}/S_summary.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["quantity", "S_eVbohr", "err", "note"])
    w.writerow(["P1_ehrenfest_initial_drag", f"{S_P1:.4f}", f"{S_P1_err:.4f}", f"v0=1.3 in-slab v>=1.1 n={n_P1}"])
    w.writerow(["P2_constv_deposit_over_L", f"{S_P2:.4f}", f"{S_P2_err:.4f}", "coupling-subtracted"])
    w.writerow(["P2_constv_raw_Etot_over_L", f"{S_P2_raw:.4f}", "", "note's naive metric"])
    w.writerow(["S_WP_p5_wp_v1p3", f"{S_WP}", "", "UPPER BOUND"])
    w.writerow(["bulk_sigma0p5_classical", f"{S_bulk}", "", "v0=1.3"])
    w.writerow(["lindhard_point", f"{S_lind_pt:.4f}", "", f"rs={RS:.3f}"])
    w.writerow(["lindhard_sigma", f"{S_lind_sig:.4f}", "", "sigma_charge=0.354"])
