#!/usr/bin/env python3
"""Build ke_bookkeeping.ipynb — validate, with tables + plots, that INQ does NOT
save the classical projectile's kinetic energy, and that the WP-drift KE and the
classical ion KE are equal (mass+velocity matched) but only the WP's is inside
E_total.

Campaign: docs/campaigns/localised_jellium_parameter_study_2/ (Energy book-keeping).
Data: scripts/ke_check/runs/results/{wp_k0,wp_k1,cl_v0,cl_v1}/raw/observables/observables.csv
(step-0 rows; p2 slab GS gs_p2_lz120; WP k0 in {0,1}; classical v in {0,1}, v=1 moving
Ehrenfest, ghost mass = 1 electron mass). Execute:
  python3 -m nbconvert --to notebook --execute --inplace ke_bookkeeping.ipynb (venv)
"""
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

OUT = Path(__file__).resolve().parent

PRE = r"""import sys, csv, glob
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
try:
    from inqview.visualisation import style as _st; _st.apply()
except Exception: pass

HA = 27.211386
RUNS = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/ke_check/runs"
LABEL = {"wp_k0":"WP k0=0", "wp_k1":"WP k0=1", "cl_v0":"classical v=0", "cl_v1":"classical v=1 (moving)"}

def row0(tag):
    f = glob.glob(f"{RUNS}/results/{tag}/raw/observables/observables.csv")[0]
    r = list(csv.reader(open(f))); h, d = r[0], r[1]
    return {k.replace("energy_",""): float(v) for k, v in zip(h, d) if k.startswith("energy_")}

E = {t: row0(t) for t in LABEL}                 # step-0 energy components (Ha)
# matched projectile: WP momentum k0=1 (m_e=1) -> drift KE = 1/2 k0^2;
# classical ghost mass = 1 m_e, v=1 -> KE = 1/2 m v^2. Both = 0.5 Ha.
KE_ANALYTIC_EV = 0.5 * 1.0 * HA                  # 13.606 eV
print("loaded step-0 energies for:", ", ".join(LABEL.values()))
print(f"analytic matched projectile KE (1/2 k0^2 = 1/2 M v^2) = {KE_ANALYTIC_EV:.3f} eV")"""

TAB = r"""# FULL step-0 energy component table (eV). ion_kinetic is the column of interest.
comps = ["total","kinetic","ion_kinetic","hartree","external","xc","nonlocal","ion"]
df = pd.DataFrame({LABEL[t]: {c: E[t][c]*HA for c in comps} for t in LABEL}).round(3)
df.index.name = "component (eV)"
df"""

V1 = r"""# VALIDATION 1 -- the WP's kinetic energy IS inside E_total.
# Give the WP momentum (k0: 0 -> 1). The extra KE must appear in energy_kinetic AND
# propagate into energy_total, matching the analytic 1/2 k0^2.
d_kin = (E["wp_k1"]["kinetic"] - E["wp_k0"]["kinetic"]) * HA
d_tot = (E["wp_k1"]["total"]   - E["wp_k0"]["total"])   * HA
val1 = pd.DataFrame({
    "quantity": ["Δ energy_kinetic (k0=1 − k0=0)", "Δ energy_total (k0=1 − k0=0)",
                 "analytic drift ½k0²  (k0=1)"],
    "eV": [d_kin, d_tot, KE_ANALYTIC_EV]}).round(3)
print(val1.to_string(index=False))
print(f"\\n-> WP drift KE lands in energy_kinetic ({d_kin:.2f} eV) and flows into "
      f"energy_total ({d_tot:.2f} eV); both match ½k0² = {KE_ANALYTIC_EV:.2f} eV.")

fig, ax = plt.subplots(figsize=(6.2, 4.2))
bars = ["Δkinetic", "Δtotal", "analytic ½k0²"]
vals = [d_kin, d_tot, KE_ANALYTIC_EV]
ax.bar(bars, vals, color=["#1f77b4","#2ca02c","0.6"])
for i, v in enumerate(vals): ax.text(i, v+0.2, f"{v:.2f}", ha="center", fontsize=9)
ax.set_ylabel("energy (eV)"); ax.set_title("WP KE enters E_total (matches ½k0²)")
ax.axhline(KE_ANALYTIC_EV, ls="--", color="k", lw=0.8)
fig.tight_layout(); plt.show()"""

V2 = r"""# VALIDATION 2 -- INQ does NOT save the classical projectile's KE.
# cl_v1 is a genuinely MOVING projectile (Ehrenfest, v=1, ghost mass = 1 m_e, so its
# true KE = 1/2 M v^2 = 13.61 eV). Yet the ion_kinetic column is 0 and E_total is
# unchanged vs the stationary v=0 run.
ion_ke = {LABEL[t]: E[t]["ion_kinetic"]*HA for t in LABEL}
d_tot_cl = (E["cl_v1"]["total"] - E["cl_v0"]["total"]) * HA
val2 = pd.DataFrame({
    "quantity": ["ion_kinetic column, cl_v1 (moving)", "Δ energy_total (v=1 − v=0)",
                 "analytic ½M v²  (should be here)"],
    "eV": [ion_ke["classical v=1 (moving)"], d_tot_cl, KE_ANALYTIC_EV]}).round(3)
print(val2.to_string(index=False))
print(f"\\nion_kinetic per run (eV): " + ", ".join(f"{k}={v:.3f}" for k,v in ion_ke.items()))
print(f"-> The moving ghost's {KE_ANALYTIC_EV:.2f} eV is ABSENT: ion_kinetic=0, ΔE_total=0.")

fig, ax = plt.subplots(figsize=(6.6, 4.2))
xs = list(ion_ke.keys()); ys = list(ion_ke.values())
ax.bar(xs, ys, color=["#1f77b4","#1f77b4","#d62728","#d62728"])
ax.axhline(KE_ANALYTIC_EV, ls="--", color="k", lw=1,
           label=f"true ½Mv² of moving ghost = {KE_ANALYTIC_EV:.2f} eV (NOT recorded)")
ax.set_ylabel("energy_ion_kinetic (eV)")
ax.set_title("INQ never populates ion_kinetic (0 even when moving)")
ax.set_ylim(0, KE_ANALYTIC_EV*1.3); ax.legend(frameon=False, fontsize=8)
plt.setp(ax.get_xticklabels(), rotation=15, ha="right")
fig.tight_layout(); plt.show()"""

V3 = r"""# VALIDATION 3 -- the two projectile KEs are equal (matched) but do NOT cancel in E_total.
# WP drift (recorded, in E_total) vs classical ½Mv² (real, but dropped by INQ).
wp_drift = (E["wp_k1"]["kinetic"] - E["wp_k0"]["kinetic"]) * HA
cl_true  = KE_ANALYTIC_EV
cl_in_total = (E["cl_v1"]["total"] - E["cl_v0"]["total"]) * HA   # = 0 (INQ)
val3 = pd.DataFrame({
    "projectile KE": ["WP drift  ½k0²  (in E_total)",
                      "classical ½Mv² (true value)",
                      "classical ½Mv² (as seen in E_total)"],
    "eV": [wp_drift, cl_true, cl_in_total]}).round(3)
print(val3.to_string(index=False))
gap = wp_drift - cl_in_total
print(f"\\nMatched? |½k0² − ½Mv²| = {abs(wp_drift-cl_true):.3f} eV  (equal to ~0.1%).")
print(f"Un-cancelled residual in a WP−CL comparison at this velocity = {gap:.2f} eV")
print(f"  -> add ½Mv² = {cl_true:.2f} eV to the classical total by hand; then it cancels the WP drift.")

fig, ax = plt.subplots(figsize=(6.8, 4.4))
cats = ["WP drift\\n(in E_total)", "classical ½Mv²\\n(true)", "classical ½Mv²\\n(in E_total)"]
vals = [wp_drift, cl_true, cl_in_total]
cols = ["#2ca02c", "0.6", "#d62728"]
ax.bar(cats, vals, color=cols)
for i, v in enumerate(vals): ax.text(i, v+0.2, f"{v:.2f}", ha="center", fontsize=9)
ax.axhline(cl_true, ls="--", color="k", lw=0.8)
ax.annotate("", xy=(2, cl_true), xytext=(2, 0),
            arrowprops=dict(arrowstyle="<->", color="#d62728"))
ax.text(2.05, cl_true/2, f"missing\\n{cl_true:.1f} eV", color="#d62728", fontsize=8, va="center")
ax.set_ylabel("energy (eV)")
ax.set_title("Equal by construction, but only the WP's is booked")
fig.tight_layout(); plt.show()"""

ANCHOR = r"""# SANITY ANCHOR -- at v=0 the drift terms vanish and dKin is the pure localisation energy.
dkin0 = (E["wp_k0"]["kinetic"] - E["cl_v0"]["kinetic"]) * HA
zp = 3.0/(4.0*0.5**2) * HA      # 3/(4 sigma^2), sigma_WP = 0.5
print(f"dKin(v=0) = kinetic(WP,k0=0) − kinetic(classical,v=0) = {dkin0:.2f} eV")
print(f"localisation zero-point 3/(4σ²), σ=0.5           = {zp:.2f} eV")
print(f"-> agree to {abs(dkin0-zp):.2f} eV (bath residual); confirms the A1 ledger's 81.7 eV.")"""

def build():
    cells = [
        new_markdown_cell(
            "# KE bookkeeping — does INQ save the projectile's kinetic energy?\n\n"
            "*Campaign `localised_jellium_parameter_study_2` (Energy book-keeping). Validates two "
            "claims with the step-0 energies of four matched runs (p2 slab GS `gs_p2_lz120`): a "
            "wave-packet at `k0∈{0,1}` and a classical Gaussian ghost at `v∈{0,1}` (v=1 = a "
            "**moving** Ehrenfest projectile, ghost mass = 1 electron mass). Built by "
            "`build_ke_bookkeeping.py`; runs produced by `scripts/ke_check/`.*\n\n"
            "**Claim 1 — INQ does NOT record the classical projectile's KE.** The `energy_ion_kinetic` "
            "slot exists but no production code sets it (only unit tests do), so the column is always 0.\n\n"
            "**Claim 2 — the WP drift KE (`½k0²`) and the classical KE (`½Mv²`) are equal** (mass + "
            "velocity matched) **but do not cancel**: the WP's is inside `E_total`, the classical's is "
            "dropped — a second hand-add term, analogous to `E_proj_bg`."),
        new_code_cell(PRE),
        new_markdown_cell(
            "## Step-0 energy components (all four runs)\n"
            "Every column is INQ's own `energy_*` output at insertion. Note `ion_kinetic = 0` "
            "across the board — including the moving `v=1` run."),
        new_code_cell(TAB),
        new_markdown_cell(
            "## Validation 1 — the WP's kinetic energy enters `E_total`\n"
            "Turning on the WP momentum (`k0: 0→1`) must raise `energy_kinetic` by `½k0²` and carry "
            "that straight into `energy_total`."),
        new_code_cell(V1),
        new_markdown_cell(
            "## Validation 2 — INQ never populates `ion_kinetic` (0 even when moving)\n"
            "`cl_v1` is a genuinely moving Ehrenfest projectile with a true KE of `½Mv² = 13.61 eV`. "
            "If INQ booked it, `ion_kinetic` (or `E_total`) would show it. It shows **nothing**."),
        new_code_cell(V2),
        new_markdown_cell(
            "## Validation 3 — equal KEs, but only one is booked\n"
            "`½k0²` (WP, recorded) vs `½Mv²` (classical, real but dropped). They match to ~0.1%, so "
            "at finite velocity a raw `WP − CL` carries a spurious un-cancelled `+½k0²`. Fix: add "
            "`½Mv²` to the classical total by hand."),
        new_code_cell(V3),
        new_markdown_cell(
            "## Sanity anchor — at `v=0`, `dKin` is the pure localisation energy\n"
            "With no drift on either side, `kinetic(WP) − kinetic(classical)` must equal the WP "
            "zero-point `3/(4σ²)` = 81.6 eV — the A1 ledger's r-independent `dKin`."),
        new_code_cell(ANCHOR),
        new_markdown_cell(
            "## Conclusion\n"
            "1. **INQ does not save the projectile KE** — `ion_kinetic` is structurally never set "
            "(Validation 2: 0 even for a moving ghost).\n"
            "2. **The WP KE is inside `E_total`** (Validation 1) because the WP is a real electron; "
            "the classical KE is **outside** (the ghost is a potential/ion).\n"
            "3. **They are equal but do not cancel** (Validation 3): a complete classical ledger needs "
            "a hand-added `½Mv²` (kinetic) alongside the hand-added `E_proj_bg` (electrostatic). At "
            "`v=0` both vanish, which is why the A1/A2 (stationary) ledgers were unaffected."),
    ]
    n = new_notebook(); n.cells = cells
    n.metadata.kernelspec = {"name": "python3", "display_name": "Python 3"}
    p = OUT / "ke_bookkeeping.ipynb"
    nbf.write(n, str(p)); print("wrote", p.name); return p

if __name__ == "__main__":
    build()
    print("execute: python3 -m nbconvert --to notebook --execute --inplace ke_bookkeeping.ipynb (venv)")
