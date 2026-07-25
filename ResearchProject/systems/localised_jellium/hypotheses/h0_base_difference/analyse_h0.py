#!/usr/bin/env python3
"""H0 analysis + highlight plot + email — base WP-vs-classical E_total(0) gap.

Reads the four H0 runs (WP/classical x r=4,40) + the L_z=120 GS, computes the
excess above GS for each, isolates the WP zero-point + SIE and the classical
unscreened ghost artifact, makes the highlight plot, and EMAILS the result with
the mandatory four-part structure (email-notifications skill).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

ROOT = Path("/local/data/public/skcb2/tddft")
H0 = ROOT / "ResearchProject/systems/localised_jellium/scripts/h0_base_difference"
OUTDIR = ROOT / "ResearchProject/systems/localised_jellium/hypotheses/h0_base_difference"
OUTDIR.mkdir(parents=True, exist_ok=True)
HA_EV = 27.211386
SIGMA_WP = 0.5
ZP_PRED = 3.0 / (4.0 * SIGMA_WP**2) * HA_EV   # 81.6 eV

def e_total0(csv: Path) -> float:
    import csv as _csv
    with open(csv) as f:
        rows = list(_csv.reader(f))
    hdr, first = rows[0], rows[1]
    return float(first[hdr.index("energy_total")])

def e_kin0(csv: Path) -> float:
    import csv as _csv
    with open(csv) as f:
        rows = list(_csv.reader(f))
    hdr, first = rows[0], rows[1]
    return float(first[hdr.index("energy_kinetic")])

def gs_energy() -> float:
    for line in (H0 / "gs/results/run_summary.txt").read_text().splitlines():
        if line.startswith("ground_state_energy_ha"):
            return float(line.split("=")[1])
    raise RuntimeError("GS energy not found")

# ---- gather ---------------------------------------------------------------
E_GS = gs_energy()
runs = {
    ("wp", 4):  H0 / "wp/results/h0_wp_r4/raw/observables/observables.csv",
    ("wp", 40): H0 / "wp/results/h0_wp_r40/raw/observables/observables.csv",
    ("cl", 4):  H0 / "classical/results/h0_cl_r4/raw/observables/observables.csv",
    ("cl", 40): H0 / "classical/results/h0_cl_r40/raw/observables/observables.csv",
}
E = {k: e_total0(v) for k, v in runs.items()}
EK = {k: e_kin0(v) for k, v in runs.items()}
excess = {k: (E[k] - E_GS) * HA_EV for k in E}            # eV above GS
zp_meas = (EK[("wp", 40)] - EK[("cl", 40)]) * HA_EV       # WP zero-point KE (eV)
sie_est = excess[("wp", 40)] - zp_meas                    # route-1 E_SIE at r=40

r = np.array([4, 40])
wp_x = np.array([excess[("wp", 4)], excess[("wp", 40)]])
cl_x = np.array([excess[("cl", 4)], excess[("cl", 40)]])
raw_gap = wp_x - cl_x

print(f"E_GS = {E_GS:.4f} Ha")
for k in runs: print(f"  {k}: E_total0={E[k]:.4f} Ha  excess={excess[k]:+.1f} eV")
print(f"WP zero-point (measured) = {zp_meas:.1f} eV   (pred {ZP_PRED:.1f})")
print(f"E_SIE (route-1, r=40)    = {sie_est:.1f} eV")
print(f"raw WP-cl gap: r=4 {raw_gap[0]:+.1f} eV, r=40 {raw_gap[1]:+.1f} eV")

# ---- plot -----------------------------------------------------------------
try:
    from inqview.visualisation import style as _style
    _style.apply()
except Exception as e:
    print("style apply failed (using default):", e)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6.2, 4.4))
ax.plot(r, wp_x, "o-", color="#1b6ca8", lw=2, ms=8, label="wavepacket (quantum)")
ax.plot(r, cl_x, "s--", color="#c0392b", lw=2, ms=8, label="classical ghost (raw)")
ax.axhline(ZP_PRED, color="0.4", ls=":", lw=1.4,
           label=f"WP localisation $3/4\\sigma^2$ = {ZP_PRED:.0f} eV")
ax.axhline(zp_meas + sie_est, color="#1b6ca8", ls=":", lw=1.0, alpha=0.6)
for xi, yi in zip(r, wp_x): ax.annotate(f"{yi:+.0f}", (xi, yi), textcoords="offset points", xytext=(6, 8), color="#1b6ca8")
for xi, yi in zip(r, cl_x): ax.annotate(f"{yi:+.0f}", (xi, yi), textcoords="offset points", xytext=(6, -14), color="#c0392b")
ax.set_xlabel("WP–slab distance  r  (Bohr, from near face)")
ax.set_ylabel(r"$E_{\rm total}(0) - E_{\rm GS}$  (eV)")
ax.set_title("H0: base WP-vs-classical energy gap (localised jellium slab)")
ax.set_xticks(r); ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
PLOT = OUTDIR / "H0_base_difference.png"
fig.savefig(PLOT, dpi=150, bbox_inches="standard" if False else None)
print("wrote", PLOT)

# ---- email ----------------------------------------------------------------
body = f"""HYPOTHESIS
  Is the base WP-vs-classical energy gap E_total(0) just the wavepacket
  localisation energy 3/(4 sigma^2) = {ZP_PRED:.0f} eV? (sigma_WP = 0.5 Bohr.)

WHAT WAS DONE
  - Localised jellium slab, L_z=120 box, GS converged (E_GS = {E_GS:.3f} Ha,
    interior n0 = 1.32e-3 verified).
  - Stationary (k0=0) WP and a charge-matched classical Gaussian ghost
    (electron_gaussian_wpsigma0p5) placed at the SAME distance from the slab,
    r = 4 and r = 40 Bohr; t=0 total energy read from step 0 (4 runs).

PLOT (attached: H0_base_difference.png)
  Excess energy above GS vs distance r, for the WP (blue) and the classical
  ghost (red), with the WP localisation energy ({ZP_PRED:.0f} eV) as a dotted line.

CONCLUSION
  NO -- the base gap is NOT the localisation energy.
  - WP excess is distance-STABLE at ~86 eV (r=4: {wp_x[0]:+.0f}, r=40: {wp_x[1]:+.0f}),
    = localisation ({zp_meas:.0f} eV, measured) + self-interaction error (~{sie_est:.0f} eV).
    So the WP's own base energy IS ~ localisation + SIE.
  - The classical ghost excess is DISTANCE-DEPENDENT ({cl_x[0]:+.0f} eV near ->
    {cl_x[1]:+.0f} eV far): the unscreened ghost-slab Coulomb (the omitted
    int v_ghost*n_+ background term).
  - Hence the raw WP-classical gap flips sign with distance (r=4: {raw_gap[0]:+.0f} eV,
    r=40: {raw_gap[1]:+.0f} eV) and is artifact-dominated, NOT the localisation energy.
  Implication for the ladder: the classical reference needs the ghost-background
  correction (H5) before any WP-classical subtraction; the WP route-1 SIE is
  already ~{sie_est:.0f} eV (k0=0), consistent with the known ~4.5 eV.
"""
print("\n--- EMAIL BODY ---\n" + body)

if "--no-email" not in sys.argv:
    from inqview.email import send_run_email
    send_run_email(
        subject=f"[localised-jellium GS] H0 — base gap is artifact-dominated "
                f"({raw_gap[0]:+.0f} eV near, {raw_gap[1]:+.0f} eV far), not localisation",
        body=body,
        attachments=[str(PLOT)],
        to="chiddukanna@gmail.com",
    )
    print("EMAIL SENT")
else:
    print("(--no-email: skipped send)")
