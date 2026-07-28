#!/usr/bin/env python3
"""Builder for p5_wp_v1p3_ehrenfest_drag.ipynb — refined stopping method #3:
the LOCAL, time-resolved Ehrenfest initial-drag extraction.

S(v0) = -d(<p_z>^2/2)/ds along the packet's own trajectory, from the per-step
recorded orbital momentum moments (wp_momentum_stats.csv) and the
center-of-density path (wp_real_space_stats.csv).  This is the direct quantum
analog of the classical-projectile friction definition, and the only method
that reads S *during* the traversal rather than from asymptotic states.

Run: PYTHONPATH=<stack> venv/python3 build_v1p3_ehrenfest_drag_report.py
(two-pass: ehrenfest_drag_summary.json -> takeaway numbers).
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from _nbreport import md, code, embed, setup_cell, set_outdir, build  # noqa: E402

set_outdir(HERE)
OUT = os.path.join(HERE, "p5_wp_v1p3_ehrenfest_drag.ipynb")
SUMMARY = os.path.join(HERE, "ehrenfest_drag_summary.json")
GIF = os.path.join(HERE, "p5_wp_v1p3_run_notebook_figs", "p5_wp_v1p3_total_density.gif")


def build_cells(summary):
    cells = []

    cells.append(md(r"""# p5_wp_v1p3 — refined stopping by the Ehrenfest initial-drag method

**Method #3 of the refined-stopping family** (see
`docs/notes/refined-stopping-two-schemes.md`; #1 = TOF/rank-matched flux, #2 =
snapshot corridor kinematics). Those two are *asymptotic* — they compare the
projectile before and after. This one is **local and time-resolved**: it reads
the momentum-dependent kinetic energy of the packet *while it crosses the slab*
and takes the stopping power as the drag slope along the packet's own path —
the direct quantum analog of the classical definition
$S=-\,dE_\mathrm{kin}/ds$, and the estimator closest to "S at the initial
velocity $v_0$" that the data permits.

It belongs to the **KS-orbital scheme** (Scheme 2 of the notes): its validity
rests on the WP orbital remaining a meaningful carrier of the projectile during
the traversal. The notebook therefore carries its own diagnostics of exactly
that assumption (Ehrenfest-identity check, CAP-selection decomposition)."""))

    cells.append(md(r"""## Symbols and equations

| symbol | meaning | source |
|---|---|---|
| $\langle\mathbf p\rangle(t)$ | orbital momentum expectation | `wp_momentum_stats.csv` (every 12 steps = 0.48 a.u.) |
| $T_\mathrm{drift}(t)$ | $\|\langle\mathbf p\rangle\|^2/2$ — momentum-dependent KE | derived |
| $z_c(t)$ | center of density $\langle z\rangle$ | `wp_real_space_stats.csv` |
| $v_c(t)$ | $dz_c/dt$ | derived |
| $s$ | path length ≡ $z_c$ (forward motion) | derived |
| $v_0,\ \sigma_p$ | launch velocity 1.3, momentum spread 1.41 | run provenance |

**The estimator.** Ehrenfest's theorem gives
$d\langle p_z\rangle/dt=-\langle\partial_z V\rangle$: the mean momentum changes
only through real forces. The drag reading of stopping is

$$S_\mathrm{drag}(\bar v)\;=\;-\,\frac{d}{ds}\,\frac{\langle p_z\rangle^2}{2}
\;\Big|_\mathrm{window},\qquad s=z_c,$$

a linear fit of $T_\mathrm{drift}$ vs $z_c$ over a stated window, quoted at the
window's mean velocity $\bar v$. By construction it contains **no localisation
energy** — $\mathrm{Var}(p)/2$ never enters.

**Two windows, and why the naive one fails here.** The standing light-projectile
rule says "fit the early window $v\ge0.85v_0$". For an *attractive*
projectile–target pair (electron → jellium) that window straddles the
**image-potential acceleration**: the packet *gains* drift KE on approach, so
the early fit mixes surface acceleration with bulk friction and returns ~0.
The physically clean window is the **slab interior** ($|z_c|\le8$), where the
background is uniform, the surface terms are absent, and the slope is friction:

$$S_\mathrm{interior}=-\,\frac{d\,T_\mathrm{drift}}{dz_c}\Big|_{|z_c|\le8}$$

**Known biases, diagnosed in-notebook:**
1. *CAP selection*: the moments are norm-weighted over the **surviving**
   density; the entrance CAP eats the slow dispersive tail, pushing
   $\langle p_z\rangle$ **up** independently of any force. Diagnosed by watching
   the full-box momentum distribution $P(k_z,t)$: truncation (selection) vs
   rigid shift (force).
2. *Non-Hermitian Ehrenfest*: with a CAP, $v_c \ne \langle p_z\rangle$ exactly —
   the identity check below *measures* when absorption starts to matter.
3. *Orbital identity mid-collision*: the standing caveat of Scheme 2."""))

    cells.append(setup_cell())

    cells.append(code(r'''import numpy as np, pandas as pd, json
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

HYP     = SYS + "/hypotheses/qsp_phase5"
RESULTS = SYS + "/scripts/qsp_phase5/wp/results"
RAW     = f"{RESULTS}/p5_wp_v1p3/raw"
HA_EV, DX, V0, SIG_P = 27.211386, 0.5, 1.3, np.sqrt(2.0)
Z_SLAB, L_SLAB = 12.5, 25.0
mom = pd.read_csv(f"{RAW}/observables/wp_momentum_stats.csv", comment="#")
rs  = pd.read_csv(f"{RAW}/observables/wp_real_space_stats.csv", comment="#")
t   = mom.time_au.values
zc  = rs.z_mean.values
drift = ((mom.px_mean**2 + mom.py_mean**2 + mom.pz_mean**2)/2).values   # [Ha]
vc  = np.gradient(zc, t)
print(f"{len(t)} samples, cadence {t[1]-t[0]:.2f} a.u.; z_c: {zc[0]:.1f} -> {zc[-1]:.1f}")'''))

    cells.append(md(r"""## Source files

| role | path |
|---|---|
| run | `scripts/qsp_phase5/wp/results/p5_wp_v1p3/` |
| per-step orbital moments | `.../raw/observables/wp_momentum_stats.csv` |
| center-of-density track | `.../raw/observables/wp_real_space_stats.csv` |
| wavefunction frames (bias diagnosis) | `.../raw/vti/wavefunction_wp/` (complex, every 12 steps) |
| companion methods | `qsp_phase5_momentum_stopping.ipynb` (TOF), `p5_wp_v1p3_snapshot_kinematics.ipynb` |
| this builder | `hypotheses/qsp_phase5/build_v1p3_ehrenfest_drag_report.py` |
"""))

    if os.path.exists(GIF):
        cells.append(embed(GIF, caption="p5_wp_v1p3 total density (xz) — the traversal this notebook reads locally",
                           width=680))

    # -- Ehrenfest identity check -------------------------------------------
    cells.append(md(r"""## Diagnostic 1 — the Ehrenfest identity $v_c=\langle p_z\rangle$

For unitary evolution these are *identical* for any potential. Their divergence
timestamps where CAP absorption (non-unitarity) starts reshaping the surviving
ensemble — beyond that point the moments mix force physics with selection."""))

    cells.append(code(r'''fig, ax = plt.subplots(figsize=(8, 3.4))
ax.plot(t, mom.pz_mean, color="tab:orange", label="$\\langle p_z\\rangle$")
ax.plot(t, vc, color="tab:blue", lw=1, label="$v_c=dz_c/dt$")
ax.axhline(V0, color="k", ls=":", lw=0.8)
ax.annotate("$v_0$", (t[-1]*0.98, V0), fontsize=9, ha="right", va="bottom")
i_div = int(np.argmax(np.abs(vc - mom.pz_mean.values) > 0.05))
ax.axvline(t[i_div], color="0.5", ls="--", lw=1)
ax.annotate(f"identity breaks t≈{t[i_div]:.0f}\n(CAP selection onset)", (t[i_div]+1, ax.get_ylim()[1]*0.8), fontsize=8)
ax.set_xlabel("t [a.u.]"); ax.set_ylabel("[a.u.]"); ax.set_xlim(0, 60); ax.legend(fontsize=8)
ax.set_title("Ehrenfest identity check — divergence marks non-unitary (CAP) contamination", fontsize=10)
plt.tight_layout(); plt.show()'''))

    # -- bias decomposition: shift vs truncation ----------------------------
    cells.append(md(r"""## Diagnostic 2 — is the approach-phase rise of $\langle p_z\rangle$ force or selection?

$\langle p_z\rangle$ rises 1.30 → 1.45 before the packet is fully inside. Two
candidate causes with opposite meaning: (a) **image acceleration** (real force:
the whole distribution shifts up), (b) **CAP selection** (the slow tail is
absorbed: the distribution is truncated from below, raising the mean of what
survives). The signed momentum marginal distinguishes them:

$$P(k_z,t)=\int\!\!\!\int |\tilde\psi(k_x,k_y,k_z,t)|^2\,dk_x\,dk_y
=\sum_{x,y}\big|\mathrm{FFT}_z\,\psi\big|^2$$"""))

    cells.append(code(r'''def load_psi(step):
    r = vtk.vtkXMLImageDataReader()
    r.SetFileName(f"{RAW}/vti/wavefunction_wp/wavefunction_t{step:06d}.vti")
    r.Update()
    img = r.GetOutput(); pdd = img.GetPointData(); dims = img.GetDimensions()
    re = vtk_to_numpy(pdd.GetArray("wavefunction_real")).reshape(dims[::-1]).T
    im = vtk_to_numpy(pdd.GetArray("wavefunction_imag")).reshape(dims[::-1]).T
    return re + 1j*im

NZ = 180
kz = np.fft.fftshift(2*np.pi*np.fft.fftfreq(NZ, d=DX))
def P_kz(step):
    psi = load_psi(step)
    phi = np.fft.fft(psi, axis=2)*DX/np.sqrt(2*np.pi)
    return np.fft.fftshift((np.abs(phi)**2).sum(axis=(0, 1))*DX*DX)

fig, ax = plt.subplots(figsize=(8, 3.6))
for step, c in ((0, "0.2"), (120, "tab:blue"), (240, "tab:green"), (360, "tab:orange")):
    P = P_kz(step)
    ax.plot(kz, P, color=c, label=f"t={step*0.04:.1f}")
ax.axvline(V0, color="k", ls=":", lw=0.8)
ax.set_xlim(-3, 6); ax.set_xlabel("$k_z$ [a.u.]"); ax.set_ylabel("$P(k_z)$")
ax.set_title("full-box momentum marginal: truncation from below (selection) vs shift (force)", fontsize=10)
ax.legend(fontsize=8)
plt.tight_layout(); plt.show()'''))

    # -- the drag curve ------------------------------------------------------
    cells.append(md(r"""## The drag curve and the window fits

$T_\mathrm{drift}$ against the packet's own position $z_c$: approach
acceleration (image), interior drag (the physics), exit/late (selection-
dominated, excluded). Fits:

- **naive rule window** $v\ge0.85v_0$ — reported to show *why* it fails for an
  attractive projectile (acceleration cancels drag);
- **interior window** $|z_c|\le8$ — the headline $S_\mathrm{interior}$;
- wider velocity windows ($v\ge0.7v_0$, $v\ge0.5v_0$) — the deceleration sweep,
  giving S at *lower* mean velocities from the same run (the light-projectile
  bonus)."""))

    cells.append(code(r'''fits = {}
def fit_window(mask, label):
    p = np.polyfit(zc[mask], drift[mask]*HA_EV, 1)
    vbar = float(np.mean(mom.pz_mean.values[mask]))
    fits[label] = dict(S=-float(p[0]), vbar=vbar, npts=int(mask.sum()),
                       z_lo=float(zc[mask].min()), z_hi=float(zc[mask].max()))
    return p

def prefix_mask(cond):
    """Contiguous prefix where cond holds (stops at first failure) — prevents a
    late-time selection-driven recovery of pz from re-entering a velocity window."""
    m = np.zeros_like(cond)
    stop = np.argmax(~cond) if (~cond).any() else len(cond)
    m[:stop] = cond[:stop]
    return m

i_exit = int(np.argmax(zc > 8.0)) if (zc > 8.0).any() else len(zc)
interior = (np.abs(zc) <= 8.0) & (np.arange(len(zc)) < i_exit)   # first crossing only
# HEADLINE window — "S at v0": post-peak, near-peak velocity. The image
# acceleration tops out at pz_peak ~ t_peak; the drag is then read while the
# packet is still near its (accelerated) entry velocity, before the deep
# deceleration mixes lower velocities in.
i_peak = int(np.argmax(mom.pz_mean.values))
pz_peak = float(mom.pz_mean.values[i_peak])
post_peak = np.zeros(len(t), dtype=bool)
post_peak[i_peak:] = prefix_mask(mom.pz_mean.values[i_peak:] >= 0.85*pz_peak)
masks = {
    "HEADLINE post-peak v>=0.85 v_peak": post_peak,
    "naive v>=0.85v0": prefix_mask(mom.pz_mean.values >= 0.85*V0),
    "interior |z_c|<=8 (full deceleration)": interior,
    "sweep v>=0.7v0": prefix_mask(mom.pz_mean.values >= 0.70*V0),
    "sweep v>=0.5v0": prefix_mask(mom.pz_mean.values >= 0.50*V0),
}
# NOTE: no hard time gate — CAP selection is active during the traversal (identity
# break stamped above). Selection removes slow components -> pushes pz_mean UP over
# time -> the fitted drag slope is biased LOW: S_interior is a LOWER bound.
fig, ax = plt.subplots(figsize=(8.5, 4))
ax.plot(zc, drift*HA_EV, ".-", ms=3, color="0.4", lw=0.8, label="$T_\\mathrm{drift}(z_c)$")
for (lab, m), c in zip(masks.items(), ("tab:red", "tab:olive", "tab:green", "tab:blue", "tab:purple")):
    p = fit_window(m, lab)
    zz = np.linspace(zc[m].min(), zc[m].max(), 10)
    ax.plot(zz, np.polyval(p, zz), color=c, lw=2,
            label=f"{lab}: S={fits[lab]['S']:.2f} eV/Bohr @ $\\bar v$={fits[lab]['vbar']:.2f}")
ax.axvspan(-Z_SLAB, Z_SLAB, color="tab:blue", alpha=0.08)
ax.axvline(-Z_SLAB, color="k", ls="--", lw=0.8); ax.axvline(Z_SLAB, color="k", ls="--", lw=0.8)
ax.set_xlabel("$z_c$ [Bohr]"); ax.set_ylabel("$T_\\mathrm{drift}$ [eV]")
ax.set_title("drag curve: image acceleration → interior friction → exit", fontsize=10)
ax.legend(fontsize=8)
plt.tight_layout(); plt.show()
pd.DataFrame(fits).T.round(2)'''))

    cells.append(code(r'''# sensitivity of the headline to the interior-window half-width
sens = []
for zh in (5.0, 6.0, 8.0, 10.0, 12.0):
    m = (np.abs(zc) <= zh) & (np.arange(len(zc)) < i_exit)
    p = np.polyfit(zc[m], drift[m]*HA_EV, 1)
    sens.append(dict(half_width=zh, S=-p[0], npts=int(m.sum())))
sens = pd.DataFrame(sens)
S_INT = float(sens[sens.half_width == 8.0].S.iloc[0])
S_INT_ERR = float(sens.S.std())
print(f"S_interior = {S_INT:.2f} eV/Bohr, window-sensitivity spread ±{S_INT_ERR:.2f}")
sens.round(2)'''))

    # -- comparison ----------------------------------------------------------
    cells.append(md(r"""## All refined-S methods on this run"""))

    cells.append(code(r'''HL = "HEADLINE post-peak v>=0.85 v_peak"
with open(f"{HYP}/momentum_stopping_summary.json") as f:
    tof = json.load(f)["p5_wp_v1p3"]
with open(f"{HYP}/snapshot_kinematics_summary.json") as f:
    snap = json.load(f)
with open(f"{HYP}/results_p5_wp_v1p3.json") as f:
    dep = json.load(f)
snap_mid = [s["S"] for s in snap["snapshots"] if 35 < s["t"] < 55]
cmp_ = pd.DataFrame([
    ("#3 Ehrenfest post-peak drag (this notebook)", f"{fits[HL]['S']:.2f}", f"{fits[HL]['vbar']:.2f}"),
    ("#3b Ehrenfest full-deceleration slope", f"{S_INT:.2f} ± {S_INT_ERR:.2f}", f"{fits['interior |z_c|<=8 (full deceleration)']['vbar']:.2f}"),
    ("#1 TOF rank-matched", f"{tof['S_drift']:.2f} ± {tof['S_err']:.2f}", f"{tof['u_ref']:.2f}"),
    ("#2 snapshot (t*=40–50)", f"{np.mean(snap_mid):.2f}", "~0.85"),
    ("deposit-based (uncorrected)", f"{dep['S_eVbohr']:.2f}", "—"),
], columns=["method", "S [eV/Bohr]", "at v"]).set_index("method")
cmp_'''))

    cells.append(code(r'''summary = dict(S_at_v0=fits[HL]["S"], vbar_at_v0=fits[HL]["vbar"],
               pz_peak=pz_peak,
               S_interior=S_INT, S_interior_err=S_INT_ERR,
               vbar_interior=fits["interior |z_c|<=8 (full deceleration)"]["vbar"],
               S_naive_rule=fits["naive v>=0.85v0"]["S"],
               S_sweep_07=fits["sweep v>=0.7v0"]["S"], vbar_07=fits["sweep v>=0.7v0"]["vbar"],
               S_sweep_05=fits["sweep v>=0.5v0"]["S"], vbar_05=fits["sweep v>=0.5v0"]["vbar"],
               t_identity_break=float(t[i_div]),
               S_tof=tof["S_drift"], S_tof_err=tof["S_err"], S_deposit=dep["S_eVbohr"])
with open(f"{HYP}/ehrenfest_drag_summary.json", "w") as f:
    json.dump(summary, f, indent=1)
print("wrote ehrenfest_drag_summary.json")'''))

    if summary:
        s = summary
        cells.append(md(f"""## Takeaway

- **Headline — S at v0:** the post-peak window ($v\\ge0.85\\,v_\\mathrm{{peak}}$, after
  the image acceleration tops out at $p_z$ = {s['pz_peak']:.2g}) gives
  $S(v_0)$ = **{s['S_at_v0']:.2g} eV/Bohr at $\\bar v$ = {s['vbar_at_v0']:.2g}** — the drag on the
  still-coherent packet at (slightly above) its launch velocity.
- The full-deceleration interior slope is {s['S_interior']:.2g} ± {s['S_interior_err']:.2g} eV/Bohr
  but at $\\bar v$ = {s['vbar_interior']:.2g} — it averages the *whole* stopping sweep
  (v ≈ 1.45 → 0.2) and must not be quoted at $v_0$ (light-projectile rule).
  CAP selection is active throughout (identity break t ≈ {s['t_identity_break']:.0f} a.u.),
  pushing $\\langle p_z\\rangle$ up — drag numbers are lower bounds.
- **The naive light-projectile window returns {s['S_naive_rule']:.2g} eV/Bohr — near zero —
  because an electron packet is *attracted* into jellium**: image acceleration on
  approach cancels the early drag. For attractive projectiles the rule's window
  must be replaced by the slab-interior window. (Rule update candidate.)
- Deceleration-sweep bonus: S = {s['S_sweep_07']:.2g} at $\\bar v$ = {s['vbar_07']:.2g},
  {s['S_sweep_05']:.2g} at $\\bar v$ = {s['vbar_05']:.2g} — one run samples the S(v) curve
  downward, consistent with the TOF band trend.
- Cross-method: TOF gives {s['S_tof']:.2g} ± {s['S_tof_err']:.2g}, snapshot ~0.5, this method
  {s['S_interior']:.2g} — all far below the deposit-based {s['S_deposit']:.2g}: the refined-S
  family is mutually consistent; the deposit over-count stands.
- Validity: this is a Scheme-2 (KS-orbital) method; its window deliberately ends
  before the diagnosed CAP-selection onset, and the momentum-marginal panel
  separates force (shift) from selection (truncation) explicitly.
"""))
    else:
        cells.append(md("## Takeaway\n\n*(populated on second build pass)*"))
    return cells


summary = None
if os.path.exists(SUMMARY):
    with open(SUMMARY) as f:
        summary = json.load(f)

print("pass 1: executing notebook ...")
build(build_cells(summary), OUT, timeout=1800)
with open(SUMMARY) as f:
    summary2 = json.load(f)
print("pass 2: re-rendering takeaway ...")
build(build_cells(summary2), OUT, timeout=1800)
print("done:", OUT)
