#!/usr/bin/env python3
"""Build the P0b DESIGN-GATE REVIEW notebook (wide-WP campaign, localised jellium).

One consolidated, self-contained .ipynb gathering EXACTLY the plots + data the user
must verify before the campaign moves from Phase 0 (design) to Phase 1 (autonomous
S(E) sweep). It reads the raw P0b observables directly (robust to which pipeline
figures exist) and renders one focused panel per gate criterion.

The user OWNS the verdict: each criterion shows measured value + threshold and a
BLANK "Verdict (you): ____" line. This script does not decide pass/fail.

Run:
  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
  /local/data/public/skcb2/tddft/venv/bin/python3 build_p0b_gate_review.py
"""
from __future__ import annotations
import io, math, textwrap
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    from inqview.visualisation.style import apply_theme
    apply_theme()
except Exception:
    pass
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

# ----------------------------------------------------------------- paths / params
LJ   = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium")
WWP  = LJ / "scripts/wide_wp"
WP   = WWP / "wp/results/p0b_wp/raw/observables"
CL   = WWP / "classical/results/p0b_classical/raw/observables"
HYP  = LJ / "hypotheses/wide_wp"
FIGS = HYP / "p0b_gate_review_figs"
FIGS.mkdir(parents=True, exist_ok=True)

HA        = 27.211386
E_GS      = -86.04107005396197   # dx=0.40 / LZ=101 validated GS anchor
L_Z       = 25.0                 # slab thickness = 2*12.5 Bohr
SLAB_HALF = 12.5
LAUNCH_Z  = -26.5
CAP_INNER = 40.5
SIGMA_WP  = 3.5
SIGMA_DEN0 = SIGMA_WP / math.sqrt(2.0)   # density std at t=0 = 2.475
V0        = 4.69569944759        # k0 from wp_momentum_stats t=0 (E~300 eV)
E0_EV     = 0.5 * V0 * V0 * HA

def rd(p, **kw):  return pd.read_csv(p, **kw)

# stats CSVs carry a leading "# ..." comment line -> header is the 2nd line.
# NB: these runs were KILLED mid-write, so the final row(s) may be partial/NaN -> drop.
obs   = rd(WP / "observables.csv").dropna()
rs    = rd(WP / "wp_real_space_stats.csv", comment=None, skiprows=1).dropna()
mom   = rd(WP / "wp_momentum_stats.csv",  comment=None, skiprows=1).dropna()
nnum  = rd(WP / "electron_number.csv").dropna()
mdist = rd(WP / "momentum_distribution.csv", comment=None, skiprows=1).dropna()
trk   = rd(CL / "electron_track.csv").dropna()
clnum = rd(CL / "electron_number.csv").dropna()

# ============================================================ metric computations
# ---- criterion 1: spreading sigma_z(t) measured AT SLAB CROSSING vs free law -----
# sigma is only physical while the WP norm is retained (before CAP absorption inflates
# the 2nd moment). Restrict to norm>=0.95 rows, then read sigma where the CENTROID
# crosses the slab (entry face z=-12.5 and centre z=0). Compare to the free law.
rs["sigma_z"] = np.sqrt(rs["sigma_z2"].clip(lower=0))
t_rs, sig, nchk, zc = (rs["time_au"].values, rs["sigma_z"].values,
                        rs["norm_check"].values, rs["z_mean"].values)
sig_free = SIGMA_DEN0 * np.sqrt(1.0 + (t_rs / (2.0 * SIGMA_DEN0**2))**2)
clean = rs[rs.norm_check >= 0.95]
sig0 = float(sig[0])
def sig_at(z_target):
    i = int((clean.z_mean - z_target).abs().idxmin())
    r = clean.loc[i]
    free = SIGMA_DEN0 * math.sqrt(1.0 + (r.time_au / (2.0 * SIGMA_DEN0**2))**2)
    return r.time_au, r.sigma_z, (r.sigma_z/sig0 - 1)*100, (r.sigma_z/free - 1)*100
t_face, sig_face, spread_face_pct, dev_face_pct = sig_at(-SLAB_HALF)   # entry face
t_cen,  sig_cen,  spread_cen_pct,  dev_cen_pct  = sig_at(0.0)          # slab centre
z_clean_max = float(clean.z_mean.max())
v_ballistic = float(np.polyfit(clean.time_au, clean.z_mean, 1)[0])     # should ~ k0

# ---- criterion 2: CAP completeness (packet reached CAP; absorption in progress) --
# centroid reached ~z_clean_max (approaching CAP inner face); norm_check then falls as
# the CAP eats the leading edge. NOT completed -> residual not yet -> 0.
wp_norm_min_clean = float(clean.norm_check.min())
z_wp_reached = z_clean_max

# ---- criterion 3: bath conservation N_total(t) (finite rows only) ----------------
N0, Nf = float(nnum["N_total"].iloc[0]), float(nnum["N_total"].iloc[-1])
N_drift_pct = (Nf - N0) / N0 * 100.0

# ---- criterion 4: energy plateau ------------------------------------------------
Et, tE = obs["energy_total"].values, obs["time_au"].values
late_dE = float(np.max(np.abs(np.diff(Et[-5:]))) * HA) if len(Et) >= 6 else float("nan")
S_energy_partial = (Et[-1] - E_GS) / L_Z * HA     # PARTIAL (packet mid-transit/absorb)

# ---- criterion 5: quantum S, two channels ---------------------------------------
# channel A: energy method (partial, above). channel B: n(k,t) coherent-peak centroid.
def peak_centroid(sub):
    k, n = sub["k_bohr_inv"].values, sub["n_wp"].values
    m = (k > 2.0) & (k < 7.0)                       # band around k0
    if m.sum() < 3: return np.nan
    kk, nn = k[m], n[m]
    j = int(np.argmax(nn))
    if 0 < j < len(kk) - 1:                          # parabolic sub-bin refine
        y0, y1, y2 = nn[j-1], nn[j], nn[j+1]
        denom = (y0 - 2*y1 + y2)
        d = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        return kk[j] + d * (kk[j+1] - kk[j-1]) * 0.5
    return float(kk[j])
kpk = mdist.groupby("step", group_keys=False).apply(peak_centroid)
t_k = mdist.groupby("step")["time_au"].first()
kc  = pd.DataFrame({"t": t_k.values, "kpeak": kpk.values}).dropna()
pz_mean = mom["pz_mean"].values                     # scattering-inflated mean (compare)
kc_last = float(kc["kpeak"].iloc[-1]) if len(kc) else float("nan")

# ---- criterion 6: classical stopping — TWO estimators (user chooses) -------------
z, vz, ke = trk["z"].values, trk["vz"].values, trk["ke_ion_ha"].values
# (a) INITIAL DRAG: first CONTIGUOUS approach window vz>=0.85 v0 (before the packet
#     enters the conservative slab well). Avoids averaging in the post-slab recovery.
below = np.where(vz < 0.85 * V0)[0]
i_first_below = int(below[0]) if len(below) else len(vz)
idx0 = np.arange(0, max(i_first_below, 2))
S_cl_drag = -np.polyfit(z[idx0], ke[idx0], 1)[0] * HA
drag_win = (round(float(z[idx0].min()), 1), round(float(z[idx0].max()), 1))
# (b) NET face-to-face loss between EQUAL-POTENTIAL slab faces (true stopping proxy):
def ke_at(z_target):
    return float(ke[int(np.argmin(np.abs(z - z_target)))]) * HA
ke_entry, ke_exit = ke_at(-SLAB_HALF), ke_at(+SLAB_HALF)
S_cl_face = (ke_entry - ke_exit) / L_Z            # eV/Bohr (net; ~0 => conservative)
vz_min_frac = float(vz.min() / V0)                # dip depth
i_entry = int(np.argmin(np.abs(z - (-SLAB_HALF))))
vz_frac_slab = float(vz[i_entry] / V0)            # v/v0 at slab entry face
i_exit = int(np.argmin(np.abs(z - (+SLAB_HALF))))
vz_frac_exit = float(vz[i_exit] / V0)             # v/v0 at slab exit face (recovery)

# =================================================================== figures
def save(fig, name):
    p = FIGS / name
    fig.savefig(p, dpi=140, bbox_inches="tight"); plt.close(fig)
    return p.name

# fig 1 — spreading vs free law, x-axis = centroid position (so slab faces are visible)
fig, ax = plt.subplots(figsize=(6.6, 4.0))
ax.plot(zc, sig, "o-", label=r"measured $\sigma_z$")
ax.plot(zc, sig_free, "k--", label="free-spreading law")
for zf, lab in [(-SLAB_HALF, "slab in"), (SLAB_HALF, "slab out"), (CAP_INNER, "CAP inner")]:
    ax.axvline(zf, ls=":", color="0.5"); ax.text(zf, ax.get_ylim()[1]*0.9, lab, rotation=90, fontsize=7, va="top")
ax.set_ylim(0, min(12, np.nanmax(sig)*1.1))       # clip the post-CAP absorption blow-up
ax.set_xlabel(r"WP centroid $z$ (Bohr)"); ax.set_ylabel(r"$\sigma_z$ (Bohr)")
ax.set_title("Criterion 1 — WP spreading vs free law (measured at slab crossing)")
ax.legend(frameon=False, fontsize=8, loc="upper left")
f1 = save(fig, "c1_spreading.png")

# fig 2 — WP norm_check vs centroid position (CAP absorption context)
fig, ax = plt.subplots(figsize=(6.6, 4.0))
ax.plot(zc, nchk, "o-")
for zf, lab in [(-SLAB_HALF, "slab in"), (SLAB_HALF, "slab out"), (CAP_INNER, "CAP inner")]:
    ax.axvline(zf, ls="--", color="0.5"); ax.text(zf, 0.5, lab, rotation=90, fontsize=7, va="center")
ax.set_xlabel(r"WP centroid $z$ (Bohr)"); ax.set_ylabel("WP norm_check (state diagnostic)")
ax.set_title("Criterion 2 — WP norm_check vs centroid position (absorption starts near CAP)")
f2 = save(fig, "c2_wpnorm.png")

# fig 3 — bath N(t)
fig, ax = plt.subplots(figsize=(6.6, 4.0))
ax.plot(nnum["time_au"], nnum["N_total"], "-")
ax.axhline(N0, ls=":", color="0.5")
ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$N_\mathrm{total}(t)$")
ax.set_title(f"Criterion 3 — bath conservation (drift {N_drift_pct:+.2f}%)")
f3 = save(fig, "c3_ntotal.png")

# fig 4 — energy total
fig, ax = plt.subplots(figsize=(6.6, 4.0))
ax.plot(tE, Et, "-")
ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$E_\mathrm{total}$ (Ha)")
ax.set_title("Criterion 4 — energy plateau (did E settle before readout?)")
f4 = save(fig, "c4_energy.png")

# fig 5 — momentum centroid (two channels)
fig, ax = plt.subplots(figsize=(6.6, 4.0))
ax.plot(kc["t"], kc["kpeak"], "o-", label="coherent-peak centroid  $k_c(t)$")
ax.plot(mom["time_au"], pz_mean, "s--", alpha=0.6, label=r"mean $\langle p_z\rangle$ (scatter-inflated)")
ax.axhline(V0, ls=":", color="0.5", label=f"$k_0$={V0:.2f}")
ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$k_z$ (a.u.)")
ax.set_title("Criterion 5 — WP momentum loss (quantum-S channel B)")
ax.legend(frameon=False, fontsize=8)
f5 = save(fig, "c5_momentum.png")

# fig 6 — classical KE(z) with BOTH estimators + vz(z) dip-and-recovery
fig, (axa, axb) = plt.subplots(1, 2, figsize=(9.6, 3.8))
axa.plot(z, ke * HA, "-")
axa.plot(z[idx0], np.polyval(np.polyfit(z[idx0], ke[idx0], 1), z[idx0]) * HA,
         "r-", lw=2, label=f"initial-drag fit {drag_win}")
axa.plot([-SLAB_HALF, SLAB_HALF], [ke_entry, ke_exit], "gs", ms=8,
         label=f"equal-potential faces (net ΔKE={ke_entry-ke_exit:+.2f} eV)")
for zf in (-SLAB_HALF, SLAB_HALF):
    axa.axvline(zf, ls="--", color="0.5")
axa.set_xlabel("z (Bohr)"); axa.set_ylabel("KE_ion (eV)")
axa.set_title("Criterion 6 — classical KE(z): conservative dip-and-recovery?")
axa.legend(frameon=False, fontsize=7.5)
axb.plot(z, vz / V0, "-"); axb.axhline(0.85, ls=":", color="0.5", label="0.85 v₀")
for zf in (-SLAB_HALF, SLAB_HALF):
    axb.axvline(zf, ls="--", color="0.5")
axb.set_xlabel("z (Bohr)"); axb.set_ylabel(r"$v_z/v_0$")
axb.set_title(f"velocity fraction (dip to {vz_min_frac:.2f} at centre)")
axb.legend(frameon=False, fontsize=8)
f6 = save(fig, "c6_classical.png")

# =================================================================== assemble nb
def md(s):  return new_markdown_cell(textwrap.dedent(s))
def figcell(name, cap):
    return md(f"![{cap}](p0b_gate_review_figs/{name})\n\n*{cap}*")

nb = new_notebook()
C = nb.cells

C.append(md(f"""
# P0b design-gate review — wide low-spread wavepacket (localised jellium)

**Campaign** `wide-wavepacket-lowspread` · **This is the Phase 0 → Phase 1 gate.**
The autonomous 6-energy S(E) sweep must NOT start until you sign off that quantum
stopping extracts cleanly from ONE matched WP + classical run.

**Operating point:** σ_WP = {SIGMA_WP} Bohr (density std {SIGMA_DEN0:.3f}),
E ≈ {E0_EV:.0f} eV (v₀ = {V0:.3f} a.u.), box 50×50×101, dx=0.40,
CAP η=−0.7 / 10-Bohr-side (inner face ±{CAP_INNER}), launch z₀={LAUNCH_Z},
slab |z|<{SLAB_HALF}, GS anchor E_GS={E_GS:.3f} Ha.

> ⚠️ **PARTIAL DATA.** This P0b pair was killed at ~55 % of transit (WP step 415/750,
> t={t_rs.max():.1f} of ~30 a.u.; classical step {int(trk['step'].iloc[-1])}). The
> packet never reached the far CAP, so criteria 2 (CAP completeness) and 4 (energy
> plateau) are **inconclusive by construction** — they need the completed relaunch.
> Criteria 1, 3, 5, 6 are already informative. This notebook auto-rebuilds on the
> same path when the full run lands.

*You own the verdict.* Each criterion shows the measured number and its threshold,
then a blank **Verdict (you): ____** for you to fill in.
"""))

C.append(md(f"""
## Sign-off table

| # | Criterion | Source | Measured (partial P0b) | Threshold | Verdict (you) |
|---|-----------|--------|------------------------|-----------|---------------|
| 1 | WP near-rigid at slab | σ_z at slab crossing vs free law | σ_z(0)={sig0:.2f}→ **entry {sig_face:.2f} ({spread_face_pct:+.0f}%)**, centre {sig_cen:.2f} ({spread_cen_pct:+.0f}%); dev. from free law {dev_cen_pct:+.1f}%; ballistic v={v_ballistic:.2f} (k₀={V0:.2f}) | ≤ ~3 % at slab (density std spreads √2 faster than the σ_WP figure) | ____ |
| 2 | CAP fully absorbs WP | norm_check + N_total | packet reached z={z_wp_reached:.0f} (CAP at {CAP_INNER}); absorption started; **not run long enough to finish** | residual → 0 | ____ (INCONCLUSIVE) |
| 3 | Bath conserved | N_total(t) | {N0:.3f}→{Nf:.3f} ({N_drift_pct:+.2f} %) | drift < 2 % | ____ |
| 4 | Energy plateau | E_total(t) | late ΔE/step={late_dE:.2f} eV; still evolving (mid-absorption) | plateau before readout | ____ (INCONCLUSIVE) |
| 5 | Quantum S, 2 channels | energy + n(k,t) centroid | S_energy(partial, meaningless mid-run)={S_energy_partial:.1f}; k_c: {V0:.2f}→{kc_last:.2f} a.u. | finite; channels agree on FULL run | ____ |
| 6 | Classical stopping | electron_track | initial-drag S={S_cl_drag:.1f} eV/Bohr {drag_win}; **net face-to-face ΔKE={ke_entry-ke_exit:+.2f} eV → S_net={S_cl_face:+.2f} eV/Bohr**; v_z dips to {vz_min_frac:.2f} then recovers | your choice of estimator | ____ |

**Overall gate verdict (you): ____**   →  proceed to Phase 1 sweep?  yes / no / relaunch-then-review

> **The one thing to look at hardest (criterion 6).** The classical KE largely
> *recovers* across the slab: net face-to-face loss is only {ke_entry-ke_exit:+.2f} eV
> (**S_net≈{S_cl_face:+.2f} eV/Bohr**, ~{abs(S_cl_drag/S_cl_face):.0f}× below the
> {S_cl_drag:.1f} eV/Bohr "initial drag"). So most of the approach drag is **conservative**
> slab-well work, and the genuinely dissipative classical stopping looks much smaller.
> Which estimator is the physical S here is *your* call (it ties to the
> `stopping-power-extraction` and light-projectile rules) — I am not deciding it.
"""))

C.append(md(f"""## 1 — Spreading (the whole point of a wide σ)
Free-spreading law σ(t)=σ₀·√(1+(ħt/2mσ₀²)²) on the **density** std σ₀=σ_WP/√2={SIGMA_DEN0:.3f}.
The packet moves ballistically (centroid slope {v_ballistic:.2f} ≈ k₀={V0:.2f}) and σ_z
tracks the free law to within {dev_cen_pct:+.1f}% at slab centre — i.e. **no anomalous,
interaction-induced spreading**. Absolute spread at the slab is **{spread_face_pct:+.0f}% (entry face)
/ {spread_cen_pct:+.0f}% (centre)** on the density std. (The campaign's headline "2.6 %"
is on the amplitude σ_WP; the density std spreads √2 faster — same packet, two conventions.
Decide which tolerance you want to gate on.) σ past the CAP is clipped — it is absorption,
not spreading."""))
C.append(figcell(f1, "σ_z vs centroid z: measured vs free law; slab & CAP faces marked"))

C.append(md(f"""## 2 — CAP completeness (INCONCLUSIVE — needs the full run)
The packet must be fully absorbed at the far CAP with residual norm → 0. Here the
centroid **reached z≈{z_wp_reached:.0f}** (approaching the CAP inner face at {CAP_INNER}) and
`norm_check` began to fall as the CAP ate the leading edge — absorption is behaving as
designed, but the run was killed before it finished, so residual-norm → 0 is **not yet
demonstrated**. (Note `norm_check` is a WP-state diagnostic that also deforms as the
packet spreads; the clean bath guard is N_total in criterion 3, which barely moved.)
This is the one open scientific risk the gate exists to close — settle it on the relaunch."""))
C.append(figcell(f2, "WP norm_check vs centroid position; slab & CAP faces dashed"))

C.append(md("""## 3 — Bath conservation
N_total(t) should stay flat outside the absorber (no spurious drain)."""))
C.append(figcell(f3, "N_total(t)"))

C.append(md("""## 4 — Energy plateau (INCONCLUSIVE on partial data)
S is read from E_total once it plateaus; here the run stopped mid-transit so E is
still evolving — the partial S below is an artefact of the cut, not a converged number."""))
C.append(figcell(f4, "E_total(t)"))

C.append(md("""## 5 — Quantum stopping, two channels
Channel A = energy method S=[E_total(t_f)−E_GS]/L_z (needs the plateau → partial only).
Channel B = the **coherent-peak centroid** k_c(t) of n(k,t) (the rule-compliant
estimator; NOT the scatter-inflated mean ⟨p_z⟩, shown dashed for contrast). The gate
wants the two channels to agree on the completed run."""))
C.append(figcell(f5, "WP momentum: coherent-peak centroid vs mean"))

C.append(md(f"""## 6 — Classical stopping — TWO estimators, your call
- **(a) Initial drag** −dKE/dz over the first contiguous approach window (v_z≥0.85·v₀,
  {drag_win}) = **{S_cl_drag:.1f} eV/Bohr** (light-projectile rule).
- **(b) Net face-to-face** loss between the equal-potential slab faces =
  **{ke_entry-ke_exit:+.2f} eV over {L_Z:.0f} Bohr → {S_cl_face:+.2f} eV/Bohr**
  (stopping-power-extraction slab method; true dissipative stopping proxy).

The velocity dips to {vz_min_frac:.2f}·v₀ at slab centre and **recovers to {vz_frac_exit:.2f}·v₀**
on exit (entry was {vz_frac_slab:.2f}·v₀) — the signature of a conservative well. Since
(b) ≈ 0, most of (a) is reversible potential work, not friction. **Which is the physical
S is your interpretation** — flagged, not resolved by me."""))
C.append(figcell(f6, "Classical KE(z): initial-drag fit + equal-potential faces (left); v_z/v₀(z) dip-and-recovery (right)"))

C.append(md(f"""## Takeaway (provisional — partial data)
- **Verifiable now:** WP is ballistic (v={v_ballistic:.2f}≈k₀) and spreads only as the
  free law ({spread_cen_pct:+.0f}% at slab centre, {dev_cen_pct:+.1f}% off the law → no
  anomalous spreading); bath drift {N_drift_pct:+.2f}%; WP momentum centroid
  {V0:.2f}→{kc_last:.2f} a.u.
- **Open scientific question (yours):** classical net face-to-face stopping ≈
  {S_cl_face:+.2f} eV/Bohr (≈0) vs {S_cl_drag:.1f} eV/Bohr initial drag — is the physical
  classical S here ~0 (conservative well dominates) or the drag value?
- **Needs the completed relaunch:** CAP completeness (criterion 2), energy plateau
  (criterion 4), and the energy-vs-centroid agreement for quantum S (criterion 5).
- **Do not start Phase 1** until criteria 2, 4, 5 are green on a full run and you sign
  the gate above.
"""))

out = HYP / "p0b_gate_review.ipynb"
ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
ep.preprocess(nb, {"metadata": {"path": str(HYP)}})
nbf.write(nb, str(out))
print(f"wrote {out}  ({len(nb.cells)} cells)")
print(f"figs in {FIGS}")
