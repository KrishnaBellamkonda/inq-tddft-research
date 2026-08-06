#!/usr/bin/env python3
"""
Build + execute the run notebooks for the WP twin of the high-density classical
S(v) benchmark, and the synthesis notebook that overlays the quantum curve on the
classical one.

    python3 build_run_notebooks.py 2.0          # one velocity
    python3 build_run_notebooks.py all          # all four + synthesis
    python3 build_run_notebooks.py synthesis    # synthesis only

Each run notebook follows `.claude/rules/notebook-density-gif.md` (density-matrix
GIFs first, displayed inline) and shows the stopping-power derivation step by
step, per the user's request.

Plan: docs/plans/wavepacket-highdensity-sv-twin.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import nbformat as nbf
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import wp_hd_stopping as W  # noqa: E402


def md(t: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(t)


def code(t: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(t)


# --------------------------------------------------------------------------- #
def build_run(v: float, sigma: float = 0.5,
              launch_z: float = W.LAUNCH_Z_FAR) -> nbf.NotebookNode:
    W.set_campaign(sigma)
    W.set_launch(launch_z)
    name = W.name_for(v)
    nsteps = W.N_STEPS[v]
    alias = W.aliasing_bias_pct(v)
    aliaserr = f"{alias['sigma_pz2_err_pct']:.2f}" if alias["sigma_pz2_err_pct"] > 1e-3 \
        else f"{alias['sigma_pz2_err_pct']:.1e}"
    is_legacy = abs(sigma - 0.5) < 1e-9
    cells: list[nbf.NotebookNode] = []

    _intro = ("""Quantum twin of campaign `classical-highdensity-sv`, whose classical
Gaussian-charge projectile is replaced by an electron **wavepacket** with every
other physical parameter held fixed.""" if is_legacy else f"""Member of the **σ_WP = {sigma:g} Bohr** campaign (user instruction 2026-07-31).
Identical in every respect to the σ_WP = 0.5 twin of `classical-highdensity-sv`
— same ground state, dx, launch position, CAPs, dt, step count and velocity grid
— **except the wavepacket width**. That single change is the point: at σ = 0.5
the packet disperses at 1/(√2σ) = 1.414 Bohr/a.u. and is already 4.7–8.1 Bohr
wide when it reaches the slab face, so it is no longer a localised projectile
when it does the physics. At σ = {sigma:g} the rate is
{1/(2**0.5*sigma):.3f} Bohr/a.u. and the packet arrives at
{W.sigma_d((abs(W.LAUNCH_Z)-W.SLAB_HALF)/v):.2f} Bohr.""")

    cells.append(md(f"""# Wavepacket stopping power in a localised jellium slab — **σ = {sigma:g}, v = {v}**

{_intro}

| | |
|---|---|
| Cell | 35 × 35 × 85 Bohr, orthorhombic, `periodicity(2)` (x,y periodic; z open electrostatically) |
| Grid | dx = {W.DX_PRODUCTION} Bohr |
| Slab | 25 Bohr thick (half-width {W.SLAB_HALF}), erfc edge 1.0, N = {W.N_ELECTRONS}, r_s = {W.R_S} |
| Wavepacket | σ_WP = {W.SIGMA_WP} Bohr, k₀ = v = {v}, launched at z = {W.LAUNCH_Z} |
| CAP | two sin² bands, {W.CAP_L} Bohr per z face, η = {W.CAP_ETA} Ha, over z ∈ ±[{W.CAP_INNER}, {W.LZ/2}] |
| Propagation | dt = {W.DT}, {nsteps} steps, t = {nsteps*W.DT:.1f} a.u. |
| Fit window | t ∈ [{W.FIT_T0}, {W.FIT_T1}] a.u. |

## The width convention

σ always denotes **σ_WP**, the ψ-width (`.claude/rules/sigma-wp-convention.md`),
so the wavepacket's **density** standard deviation at t = 0 is σ_WP/√2 =
{sigma/2**0.5:.5f} Bohr. {"The classical campaign used a Gaussian *charge* of that same standard deviation (σ_pot = σ_WP/√2 = 0.35355), which is what makes this run its like-for-like quantum twin. Both halves are labelled σ = 0.5; the √2 lives inside each binary, never at the call site." if is_legacy else f"There is **no classical twin at this σ** — the classical benchmark exists only at σ_WP = 0.5. Comparisons in this notebook are therefore against the σ = 0.5 *quantum* campaign, not against classical."}

## Three limitations, stated up front

1. **The fit window is t ≤ {W.FIT_T1} a.u.** The packet disperses at 1/(√2σ) =
   {1/(2**0.5*sigma):.3f} Bohr per a.u., so its transverse periodic images overlap
   (6σ_d = L_xy = 35 Bohr) at t = {W.T_TRANSVERSE:.2f} a.u. Beyond this the packet
   is not a localised projectile and the slope is not stopping power.
   {"" if is_legacy else f"(At σ = 0.5 this window was only 4.12 a.u.; the {W.T_TRANSVERSE/4.117:.1f}× widening here is a direct consequence of the slower dispersion.)"}
2. **The CAP itself decelerates the packet.** In an empty box with no bath at
   all, the CAP alone drags ⟨p_z⟩ from 2.00 to 0.61 over 48 a.u., because the
   packet's leading edge reaches the absorbing band first. This notebook
   therefore differences the run against a **vacuum control** at the same
   velocity and the same σ.
3. **Momentum aliasing.** σ_p = 1/(√2σ_WP) = {1/(2**0.5*sigma):.3f} Bohr⁻¹ is fixed
   by σ, so the k-distribution folds at k_Nyq = π/dx = {3.14159265/W.DX_PRODUCTION:.2f}.
   At v = {v} on dx = {W.DX_PRODUCTION} this biases σ_pz² by **{aliaserr} %**.
   {"(v = 4.0 and 4.5 were excluded from the sweep for this reason: +17.9 % and +55.1 %.)" if is_legacy else "At this σ the bias is zero to machine precision at every velocity up to 4.5 — nothing reaches the fold. The four-point grid is kept only for comparability with the σ = 0.5 campaign."}

**Energy is not conserved in this run** — the CAP is non-Hermitian by
construction. The correctness gates are the ledger closure residuals and the
norm/absorption trace, not ΔE_total.
"""))

    cells.append(code(f"""import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from IPython.display import Image, display

sys.path.insert(0, "{HERE}")
import wp_hd_stopping as W

# Select the sigma campaign. This sets the run-directory prefix, the fit window
# (which is sigma-dependent via the transverse-image-overlap time) and every
# sigma-derived label downstream. Must come before anything else that uses W.
W.set_campaign({sigma})
# Select the LAUNCH DISTANCE. Orthogonal to the sigma campaign and composes with
# it (nl_ + s<sigma>_ + v<velocity>). Omitting this would silently resolve the
# FAR-launch run of the same sigma, because the near-launch tag defaults to "".
W.set_launch({launch_z})

V = {v}
NAME = W.name_for(V)
RUN = W.WP_RESULTS / NAME
print("run dir:", RUN)

summ = W.run_summary(RUN)
for k in ("run_completed", "wp_state_index", "norm_after", "max_overlap",
          "spacing_bohr", "wp_sigma_bohr", "wp_sigma_density", "wp_k0_bohr_inv",
          "cap", "cap_eta_ha", "cap_band_hi_bohr", "rt_num_steps", "wall_time_s"):
    if k in summ:
        print(f"  {{k:22s}} = {{summ[k]}}")"""))

    # ---------------------------------------------------------------- GIFs
    cells.append(md(f"""---
## 1. Density-matrix animations

Per `.claude/rules/notebook-density-gif.md` these come **first**: the animated
real-space density is the most direct picture of what the wavepacket does, and a
static carpet compresses time onto one axis and hides it.

Mid-y **xz** slices of the total density n(x,z,t), the wavepacket orbital density
|ψ_WP|², and the induced Δn = n(t) − n(0), in linear and log scaling. Slab faces
(±{W.SLAB_HALF} Bohr) and the CAP inner edges (±{W.CAP_INNER} Bohr) are marked.

VTIs are written in **physical order** by inqkit
(`.claude/rules/vti-coordinate-mapping.md`), so they are loaded through
`inqview.load_vti` and never `fftshift`ed."""))

    cells.append(code(f"""from inqview.visualisation import make_density_gif_battery

OUT = RUN / "report"
OUT.mkdir(parents=True, exist_ok=True)

save_every = int(summ.get("save_every", 12))
gifs, vmax = make_density_gif_battery(
    str(RUN), str(OUT),
    run_label=NAME, dt=W.DT * save_every,
    slab_face={W.SLAB_HALF}, cap_inner={W.CAP_INNER},
    cap_lines=({-W.LZ/2}, {-W.CAP_INNER}, {W.CAP_INNER}, {W.LZ/2}),
    frames_max=40, fps=10,
    run_title=f"WP projectile, v={{V}} — localised jellium slab",
)
print(f"built {{len(gifs)}} gifs, shared density scale vmax = {{vmax:.3e}}")"""))

    cells.append(code("""for cat, kind, path, ttl in gifs:
    print(f"=== {ttl} ===")
    display(Image(filename=path))"""))

    # ------------------------------------------------------------- loading
    cells.append(md("""---
## 2. Loading the kinematics

`wp_momentum_stats.csv` and `wp_real_space_stats.csv` are written **every step** —
they are the measurement, not a diagnostic. Segment-suffixed files from any
resume are concatenated in step order."""))

    cells.append(code("""run = W.load_run(V)
vac = W.load_vacuum(V)
print(f"slab run : {len(run.t)} steps, t = 0 .. {run.t[-1]:.1f} a.u.")
print(f"vacuum   : {'loaded' if vac is not None else 'MISSING — CAP baseline unavailable'}")

i0 = 0
print("\\nt = 0 values (analytic expectations in brackets):")
print(f"  <p_z>      = {run.pz[i0]:.4f}   [{V}]")
print(f"  T1         = {run.T1[i0]:.4f} Ha  [{0.5*(V**2 + 3/(2*W.SIGMA_WP**2)):.4f}]")
print(f"  T2         = {run.T2[i0]:.4f} Ha  [{0.5*V**2:.4f}]")
print(f"  T1-T2      = {(run.T1[i0]-run.T2[i0])*W.HA_TO_EV:.2f} eV  "
      f"[{3/(4*W.SIGMA_WP**2)*W.HA_TO_EV:.2f}]")
print(f"  s3, s4     = {run.s3[i0]:.3f}, {run.s4[i0]:.3f} Bohr  [{W.LAUNCH_Z}]")
print(f"  norm       = {run.norm[i0]:.5f}")"""))

    # ---------------------------------------------------- how S is computed
    cells.append(md(r"""---
## 3. How the stopping power is computed, step by step

The stopping power is the energy the projectile loses per unit path length,

$$S \;=\; -\frac{dT}{ds}$$

For a wavepacket both $T$ and $s$ admit **two** definitions, and the contrast
between them is the physics this study is after
(`docs/plans/bulk-jellium-ks-stopping.md` §4).

**Kinetic energy.**

$$T_1 \;=\; \frac{\langle p^2\rangle}{2m}
\qquad\qquad
T_2 \;=\; \frac{\langle p\rangle^2}{2m}$$

$T_1$ is the full orbital kinetic energy (INQ's native measure, column
`e_kin_ha`). $T_2$ keeps only the **drift**, discarding the momentum spread.
Their difference is the localisation + scattering energy,

$$T_1 - T_2 \;=\; \tfrac{1}{2}\sum_{d}\sigma_{p_d}^2 \;=\; \frac{3}{4\sigma_{WP}^2}
\;=\; """ + f"{3/(4*sigma**2):.4g}" + r"""\ \mathrm{Ha} \;=\; """ + f"{3/(4*sigma**2)*W.HA_TO_EV:.3g}" + r"""\ \mathrm{eV\ at}\ \sigma_{WP}=""" + f"{sigma:g}" + r""" .$$

A drop in $T_1$ can therefore mean *either* real drift loss *or* momentum-width
broadening (angular scattering); comparing $S_{1j}$ with $S_{2j}$ separates them.

**Position.**

$$s_3 \;=\; \langle z\rangle_{\rm circ}
\;=\; \frac{L_z}{2\pi}\arg\big\langle \psi\big|e^{i2\pi z/L_z}\big|\psi\big\rangle
\qquad\qquad
s_4 \;=\; z_0 + \int_0^t \langle p_z\rangle\,dt'$$

$s_3$ uses the **circular** (Resta phase) estimator, not the naive
$\int z|\psi|^2$: the naive centroid is discontinuous across a periodic face, and
this box wraps. $s_4$ integrates the mean momentum instead.

**Ehrenfest cross-check.** With no ions the Kohn–Sham Hamiltonian is purely local,
so $d\langle z\rangle/dt=\langle p_z\rangle/m$ *exactly* and $s_3$ and $s_4$ must
agree. Any divergence localises to CAP non-unitarity or orbital norm leakage — it
is a **validation**, not a second physics channel.

Four combinations follow, $S_{ij}=-dT_i/ds_j$, fitted by OLS over the window."""))

    # ------------------------------------------------------------ KE terms
    cells.append(md("""---
## 4. Each kinetic-energy term against time"""))
    cells.append(code("""fig, ax = plt.subplots(1, 3, figsize=(15, 4))
WIN = dict(color="0.85", zorder=0)
for a in ax: a.axvspan(W.FIT_T0, W.FIT_T1, **WIN)

ax[0].plot(run.t, run.T1*W.HA_TO_EV, "C0")
ax[0].set(xlabel="t (a.u.)", ylabel=r"$T_1=\\langle p^2\\rangle/2m$ (eV)",
          title="$T_1$ — full orbital KE")
ax[1].plot(run.t, run.T2*W.HA_TO_EV, "C1")
ax[1].set(xlabel="t (a.u.)", ylabel=r"$T_2=\\langle p\\rangle^2/2m$ (eV)",
          title="$T_2$ — drift only")
ax[2].plot(run.t, run.localisation_energy, "C2")
ax[2].axhline(3/(4*W.SIGMA_WP**2)*W.HA_TO_EV, ls=":", c="k",
              label=r"$3/4\\sigma^2$ = 81.6 eV")
ax[2].set(xlabel="t (a.u.)", ylabel=r"$T_1-T_2$ (eV)",
          title="localisation + scattering energy")
ax[2].legend(fontsize=8)
for a in ax: a.grid(alpha=.3)
fig.suptitle(f"v = {V}: kinetic-energy definitions (shaded = fit window)")
fig.tight_layout(); plt.show()"""))

    # ----------------------------------------------------- position terms
    cells.append(md("""---
## 5. Each position term against time, and the Ehrenfest residual

`s3_naive` is shown only to demonstrate *why* the circular estimator is
required — it breaks the moment the packet straddles a cell face."""))
    cells.append(code("""fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for a in ax: a.axvspan(W.FIT_T0, W.FIT_T1, **WIN)

ax[0].plot(run.t, run.s3, "C0", label=r"$s_3$ circular (unwrapped)")
ax[0].plot(run.t, run.s3_naive, "C3--", lw=1, label=r"$s_3$ naive $\\int z|\\psi|^2$")
for y in (-W.SLAB_HALF, W.SLAB_HALF): ax[0].axhline(y, ls=":", c="k", lw=.8)
for y in (-W.CAP_INNER, W.CAP_INNER): ax[0].axhline(y, ls="--", c="C3", lw=.8)
ax[0].set(xlabel="t (a.u.)", ylabel="z (Bohr)", title="centroid definitions")
ax[0].legend(fontsize=8)

ax[1].plot(run.t, run.s3, "C0", label=r"$s_3$ centroid")
ax[1].plot(run.t, run.s4, "C1--", label=r"$s_4=\\int\\langle p_z\\rangle dt$")
ax[1].set(xlabel="t (a.u.)", ylabel="z (Bohr)", title="$s_3$ vs $s_4$")
ax[1].legend(fontsize=8)

ax[2].plot(run.t, run.ehrenfest_residual, "C2")
ax[2].axhline(0, ls=":", c="k")
ax[2].set(xlabel="t (a.u.)", ylabel=r"$s_3-s_4$ (Bohr)",
          title="Ehrenfest residual (must be ~0)")
for a in ax: a.grid(alpha=.3)
fig.tight_layout(); plt.show()

m = (run.t >= W.FIT_T0) & (run.t <= W.FIT_T1)
print(f"max |s3-s4| inside the fit window = {np.abs(run.ehrenfest_residual[m]).max():.3e} Bohr")"""))

    # ------------------------------------------------------ CAP diagnostics
    cells.append(md("""---
## 6. CAP diagnostics — separating absorption from stopping

This is the section that decides whether the fitted slope means anything. The
vacuum control is the same wavepacket, same grid, same CAP, same number of steps,
in an **empty box**: whatever it shows is pure CAP attrition, because there is
nothing there to stop it."""))
    cells.append(code("""fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for a in ax: a.axvspan(W.FIT_T0, W.FIT_T1, **WIN)

ax[0].plot(run.t, run.norm/run.norm[0], "C0", label="slab run")
if vac is not None:
    ax[0].plot(vac.t, vac.norm/vac.norm[0], "C3--", label="vacuum control")
ax[0].set(xlabel="t (a.u.)", ylabel="norm / norm(0)", yscale="log",
          title="A. WP norm — CAP absorption")
ax[0].legend(fontsize=8)

ax[1].plot(run.t, run.pz, "C0", label="slab run")
if vac is not None:
    ax[1].plot(vac.t, vac.pz, "C3--", label="vacuum (CAP only)")
ax[1].axhline(V, ls=":", c="k", label=f"launch k0 = {V}")
ax[1].set(xlabel="t (a.u.)", ylabel=r"$\\langle p_z\\rangle$", title=r"B. $\\langle p_z\\rangle$")
ax[1].legend(fontsize=8)

if vac is not None:
    cc = W.cap_corrected(run, vac)
    ax[2].plot(cc.t, cc.dpz, "C2")
    ax[2].axhline(0, ls=":", c="k")
    ax[2].set(xlabel="t (a.u.)",
              ylabel=r"$\\langle p_z\\rangle_{slab}-\\langle p_z\\rangle_{vac}$",
              title="C. CAP-corrected momentum loss\\n(this is the bath's doing)")
else:
    ax[2].text(.5, .5, "vacuum control missing", ha="center", transform=ax[2].transAxes)
for a in ax: a.grid(alpha=.3)
fig.tight_layout(); plt.show()

print(f"norm remaining at the end of the fit window: "
      f"{run.norm[m][-1]/run.norm[0]:.4f}")
if vac is not None:
    print(f"vacuum control norm there:                   "
          f"{np.interp(W.FIT_T1, vac.t, vac.norm)/vac.norm[0]:.4f}")"""))

    # -------------------------------------------------------------- the fits
    cells.append(md("""---
## 7. The four stopping powers

$S_{ij} = -dT_i/ds_j$ by OLS over the fit window. $S_{13}\\approx S_{14}$ and
$S_{23}\\approx S_{24}$ by the Ehrenfest identity, so the physics contrast is
**$S_{1j}$ against $S_{2j}$**: how much of the apparent stopping is drift loss
versus momentum-width broadening."""))
    cells.append(md(f"""**Two windows are reported.** The *localised* window
t ∈ [{W.FIT_T0}, {W.FIT_T1}] a.u. is where the packet is still compact and CAP-free — but with
the campaign-matched launch at z = −24 there is 11.5 Bohr of vacuum standoff, so
the centroid has not reached the slab face (−12.5) yet. A slope fitted there
measures the packet **accelerating down the slab's attractive gradient in vacuum**
and comes out negative. The *in-slab transit* window is where the centroid is
actually inside the medium, which is where a stopping power is defined. Neither
is silently preferred — read them together, and against the vacuum control."""))

    cells.append(code(f"""t_in, t_out = W.slab_window(V)
print(f"localised window : t in [{{W.FIT_T0}}, {{W.FIT_T1}}] a.u.")
print(f"in-slab transit  : t in [{{t_in:.2f}}, {{t_out:.2f}}] a.u.  (centroid inside +/-12.5)")
fits_slab = W.fit_all_in_slab(run, V)
rows = [dict(combination=k, S_eV_per_Bohr=f.S_ev_per_bohr) for k, f in fits_slab.items()]
print()
print("S over the IN-SLAB TRANSIT:")
display(pd.DataFrame(rows).round(4))"""))

    cells.append(code("""fits = W.fit_all(run)
rows = []
for k, f in fits.items():
    rows.append(dict(combination=k, S_eV_per_Bohr=f.S_ev_per_bohr,
                     stderr=getattr(f, "S_stderr_ev_per_bohr", float("nan")),
                     n_points=getattr(f, "n_points", len(run.t[m]))))
tab = pd.DataFrame(rows)
display(tab.round(4))

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for (k, f), c in zip(fits.items(), ["C0", "C1", "C2", "C3"]):
    s = run.s3 if k.endswith("3") else run.s4
    T = run.T1 if k[2] == "1" else run.T2
    ax[0].plot(s, T*W.HA_TO_EV, c, label=k, lw=1)
ax[0].set(xlabel="s (Bohr)", ylabel="T (eV)", title="T vs s (all four)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)

ax[1].bar(tab.combination, tab.S_eV_per_Bohr, color=["C0","C1","C2","C3"])
ax[1].set(ylabel="S (eV/Bohr)", title=f"v = {V}: the four definitions")
ax[1].grid(alpha=.3, axis="y")
fig.tight_layout(); plt.show()"""))

    # ------------------------------------- WP-kinetic normalisation correction
    cells.append(md(r"""---
## 8. Energy ledger — with INQ's norm-divided kinetic term corrected

**Read this before reading any energy number below.** INQ reports each orbital's
kinetic energy as $\mathrm{occ}\cdot\langle\psi|T|\psi\rangle/\langle\psi|\psi\rangle$
(`inq/src/hamiltonian/energy.hpp:50-55`, used only at `:83`). Every other term is
density-based and already extensive. Under a CAP the wavepacket's norm decays, so
its kinetic contribution keeps reporting the **per-particle mean** instead of
leaving the ledger — `energy_kinetic`, and hence `energy_total`, are inflated.
This is the known CAP energy-oscillation artefact.

Verified on this machine 2026-07-30: `inq-study/src/hamiltonian/energy.hpp` is
**byte-identical** to stock `inq`'s, so the engine-level fix is *not* applied to
these runs. The correction is therefore applied here, in post-processing:

$$E_{\rm total}^{\rm corr}(t) \;=\; E_{\rm total}^{\rm rep}(t)
\;-\; \mathrm{occ}\cdot\langle T\rangle(t)\Big(\tfrac{1}{\mathrm{norm}(t)}-1\Big)
\;=\; E_{\rm total}^{\rm rep}(t) \;-\; \mathrm{occ}\cdot T_1(t)\,\big(1-\mathrm{norm}(t)\big)$$

since $T_1$ (`e_kin_ha`) *is* the norm-divided $\langle T\rangle/\mathrm{norm}$, and
occ = 1. All three inputs are written **every step**, so this is exact at full
cadence — unlike the original `wp_kinetic_normalization_fix.py`, which had to
reconstruct $\langle T\rangle$ and the norm from ~100 sparse wavefunction VTIs.

Two checks below: the correction must be **exactly zero at t = 0** (norm = 1), and
the bare WP kinetic content $T_1\cdot\mathrm{norm}$ must fall by the same amount
the corrected total drifts, up to what the bath actually absorbs."""))

    cells.append(code("""corr = W.wp_kinetic_norm_correction(RUN)

fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for a in ax: a.axvspan(W.FIT_T0, W.FIT_T1, **WIN)

ax[0].plot(corr.time_au, corr.e_total_raw_ev, "C3--", lw=1.2,
           label="reported (norm-divided kinetic)")
ax[0].plot(corr.time_au, corr.e_total_corrected_ev, "C0", lw=1.6,
           label="corrected (extensive)")
ax[0].set(xlabel="t (a.u.)", ylabel=r"$E_{total}$ (eV)",
          title="A. Energy ledger, raw vs corrected")
ax[0].legend(fontsize=8)

ax[1].plot(corr.time_au, corr.correction_ev, "C2")
ax[1].set(xlabel="t (a.u.)", ylabel="correction (eV)",
          title=r"B. $T_1(1-\\mathrm{norm})$ — the artefact removed")

ax[2].plot(corr.time_au, corr.norm_wp, "C1", label="WP norm")
ax[2].set(xlabel="t (a.u.)", ylabel="norm", title="C. WP norm (drives the artefact)")
ax[2].legend(fontsize=8)
for a in ax: a.grid(alpha=.3)
fig.tight_layout(); plt.show()

print(f"correction at t=0            : {corr.correction_ev.iloc[0]:.6f} eV  (must be 0)")
lost = corr.wp_kinetic_bare_ev.iloc[0] - corr.wp_kinetic_bare_ev.iloc[-1]
draw = corr.e_total_raw_ev.iloc[0] - corr.e_total_raw_ev.iloc[-1]
dcor = corr.e_total_raw_ev.iloc[0] - corr.e_total_corrected_ev.iloc[-1]
print(f"bare WP kinetic lost         : {lost:8.2f} eV")
print(f"CORRECTED E_total drift      : {dcor:8.2f} eV   <- tracks the line above")
print(f"raw (uncorrected) drift      : {draw:8.2f} eV   <- hides "
      f"{dcor-draw:.1f} eV of absorbed kinetic energy")
print(f"residual (not carried off by the packet) : {dcor-lost:6.2f} eV")"""))

    # ------------------------------------------------------ energy ledger
    cells.append(md("""---
## 8b. Full energy decomposition and the pairwise Coulomb ledger

Every component is written each step. The ledger uses the WP closure of
`inqkit::jellium::interaction_energies`,

$$E_{\\rm hartree}=E_{SS}+E_{PS}+E_{PP}\\ ,\\qquad E_{\\rm external}=E_{SB}+E_{PB}$$

with P = the wavepacket orbital density and S = bath. The two `*_check` columns
are computed independently and must equal INQ's own scalars — that residual is
the correctness gate that replaces energy conservation here."""))
    cells.append(code("""try:
    ix = W.load_interactions(RUN)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    for a in ax: a.axvspan(W.FIT_T0, W.FIT_T1, **WIN)
    obs_cols = [c for c in ("energy_total","energy_kinetic","energy_hartree",
                            "energy_xc","energy_external") if c in ix]
    for c in obs_cols:
        ax[0].plot(ix.time_au, ix[c]*W.HA_TO_EV, lw=1, label=c.replace("energy_",""))
    ax[0].set(xlabel="t (a.u.)", ylabel="E (eV)", title="energy decomposition")
    ax[0].legend(fontsize=7)

    for c, lab in [("e_ss","$E_{SS}$"),("e_pp","$E_{PP}$"),("e_ps","$E_{PS}$"),
                   ("e_sb","$E_{SB}$"),("e_pb","$E_{PB}$")]:
        if c in ix: ax[1].plot(ix.time_au, ix[c]*W.HA_TO_EV, lw=1, label=lab)
    ax[1].set(xlabel="t (a.u.)", ylabel="E (eV)", title="pairwise Coulomb ledger")
    ax[1].legend(fontsize=7)

    for c, lab in [("hartree_residual","Hartree"),("external_residual","external")]:
        if c in ix: ax[2].plot(ix.time_au, ix[c]*W.HA_TO_EV, lw=1, label=lab)
    ax[2].axhline(0, ls=":", c="k")
    ax[2].set(xlabel="t (a.u.)", ylabel="residual (eV)",
              title="closure residuals (gate: ~0)")
    ax[2].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); plt.show()

    for c in ("hartree_residual","external_residual"):
        if c in ix:
            print(f"max |{c}| = {ix[c].abs().max()*W.HA_TO_EV:.3e} eV")
except Exception as e:
    print("interactions.csv unavailable:", e)"""))

    # ---------------------------------------------------------- comparison
    cells.append(md("""---
## 9. Against the classical benchmark at this velocity"""))
    cells.append(code("""cl = W.classical_reference()
row = cl.loc[(cl.v - V).abs().idxmin()]
print(f"classical at v = {row.v}:  S = {row.S_eV_per_Bohr:.3f} eV/Bohr "
      f"(E_absorbed = {row.E_absorbed_eV:.1f} eV, v_final = {row.v_final:.2f})")
print("quantum:")
for k, f in fits.items():
    print(f"  {k}: {f.S_ev_per_bohr:8.3f} eV/Bohr")
print("\\nNOTE: the classical S is a DEPOSIT measure (E_total plateau - E_GS)/L_slab.")
print("The quantum S is -dT/ds of the projectile orbital. They answer the same")
print("physical question but are not the same estimator — see the plan, section 4.")"""))

    cells.append(md(f"""---
## 10. Takeaway and what is NOT established

Fill in the numbers above, then read them against these caveats:

- The fit window is **t ≤ {W.FIT_T1} a.u.** ({int(W.FIT_T1/W.DT)} steps of {nsteps}).
  Beyond it the packet is transversely overlapping its own periodic images and the
  CAP is removing it; the later steps are recorded for the animations and the
  absorption physics, not for the slope.
- σ_pz² carries a **+{aliaserr} %** aliasing bias at this velocity and dx = {W.DX_PRODUCTION}.
- `energy_total` is not conserved (non-Hermitian CAP). The closure residuals in
  section 8 are the correctness evidence instead.
- The classical comparison is against **published summary values only** — the raw
  classical per-step data did not survive the machine migration, so no per-step
  overlay or WP−classical difference GIF is possible.
- The classical and quantum stopping powers are **different estimators** of the
  same physics (energy deposit vs projectile KE loss), not the same number
  computed two ways.
"""))

    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python",
                                 "name": "python3"}
    return nb


# --------------------------------------------------------------------------- #
def build_synthesis(sigma: float = 0.5,
                    launch_z: float = W.LAUNCH_Z_FAR) -> nbf.NotebookNode:
    W.set_campaign(sigma)
    W.set_launch(launch_z)
    is_legacy = abs(sigma - 0.5) < 1e-9 and abs(launch_z - W.LAUNCH_Z_FAR) < 1e-9
    sfx = "" if is_legacy else "_" + _campaign_tag(sigma, launch_z)
    _why4 = ("""The quantum curve has **four** points (v = 2.0, 2.5, 3.0, 3.5). v = 4.0 and 4.5
were excluded: σ_WP = 0.5 fixes the packet's momentum width at σ_p = 1.414 Bohr⁻¹,
whose k-distribution folds at k_Nyq = π/dx, biasing σ_pz² by +17.9 % and +55.1 %
at dx = 0.4. Recoverable at dx = 0.30 (≤0.11 % everywhere).""" if is_legacy else
f"""The curve has four points (v = 2.0, 2.5, 3.0, 3.5), matching the σ = 0.5
campaign's grid. Unlike that campaign, the exclusion of v = 4.0/4.5 does **not**
apply here: σ_p = 1/(√2σ) = {1/(2**0.5*sigma):.3f} Bohr⁻¹ against k_Nyq = 7.85, so
aliasing is zero to machine precision at every velocity. The grid is held fixed
purely for point-for-point comparability.""")
    cells = [md(f"""# Wavepacket S(v) — **synthesis, σ_WP = {sigma:g} Bohr**

Quantum S(v) across the four production velocities{", overlaid on the six-point classical benchmark curve of campaign `classical-highdensity-sv`" if is_legacy else " (the classical benchmark is shown for orientation only — it exists at σ_WP = 0.5, not at this σ)"}.

{_why4}
"""),
        code(f"""import sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, "{HERE}")
import wp_hd_stopping as W
W.set_campaign({sigma})
W.set_launch({launch_z})
print("campaign:", W.campaign_label(), "launch z =", W.current_launch_z(),
      " fit window t <=", W.FIT_T1, "a.u.")

rows = []
for v in W.VELOCITIES:
    try:
        run = W.load_run(v)
    except FileNotFoundError as e:
        print(f"v={{v}}: MISSING ({{e}})"); continue
    t0, t1 = W.slab_window(v)
    slab = W.fit_all_in_slab(run, v)      # centroid INSIDE the medium — the physical S
    loc  = W.fit_all(run)                 # localised window — diagnostic only
    i_ex = int(np.argmin(np.abs(run.t - t1)))
    d = dict(v=v, t_in=round(t0, 2), t_out=round(t1, 2),
             norm_at_slab_exit=run.norm[i_ex], norm_end=run.norm[-1]/run.norm[0],
             aliasing_sigma_pz2_pct=W.aliasing_bias_pct(v)["sigma_pz2_err_pct"])
    for k, f in slab.items():
        d[k] = f.S_ev_per_bohr            # S_13 .. S_24, IN-SLAB
    d["broadening_S13_minus_S23"] = slab["S_13"].S_ev_per_bohr - slab["S_23"].S_ev_per_bohr
    dep = W.deposit_stopping(v)           # classical Definition 2 applied to the WP
    d.update(dep)
    for k, f in loc.items():
        d[k + "_localised"] = f.S_ev_per_bohr
    rows.append(d)
q = pd.DataFrame(rows)
display(q.round(4))
q.to_csv("wp_S_summary{sfx}.csv", index=False)
print("wrote wp_S_summary{sfx}.csv")
print()
print("S_1j uses the FULL orbital KE <p^2>/2m; S_2j the DRIFT only <p>^2/2m.")
print(f"The *_localised columns are the [{{W.FIT_T0}}, {{W.FIT_T1}}] a.u. window, which sits BEFORE")
print("the packet reaches the slab (11.5 Bohr standoff) and is therefore negative")
print("— it measures acceleration down the slab's attractive gradient, not stopping.")"""),
        code(f"""cl = W.classical_reference()

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(cl.v, cl.S_eV_per_Bohr, "ko-", lw=2,
        label="classical benchmark (sigma=0.5, 6 pts)")
ax.plot(q.v, q.S_13, "C0s--", label="quantum S_13: full <p^2>/2m")
ax.plot(q.v, q.S_23, "C1^--", label="quantum S_23: drift <p>^2/2m")
ax.fill_between(q.v, q.S_23, q.S_13, color="C2", alpha=0.18,
                label="momentum-width broadening")
ax.plot(q.v, q.S_deposit_corrected, "C3v-.",
        label="quantum S_deposit: (E_total(t_f) - E_GS)/L_slab")
ax.set(xlabel="v (a.u.)", ylabel="S (eV/Bohr)",
       title="Stopping power, r$_s$=4.18 slab, $\\\\sigma_{{WP}}$ = {sigma:g} Bohr")
ax.grid(alpha=.3); ax.legend()
fig.tight_layout(); fig.savefig("wp_vs_classical_Sv{sfx}.png", dpi=140); plt.show()
print("wrote wp_vs_classical_Sv{sfx}.png")"""),
        md("""## The third quantum curve — the classical estimator applied to the wavepacket

`S_deposit = (E_total(t_final) - E_GS)/L_slab` is the SAME formula the classical
benchmark uses, now evaluated on the wavepacket runs (E_GS = the dx = 0.40
production ground state, L_slab = 25 Bohr). It is plotted from the **norm-corrected**
ledger.

It comes out far BELOW the classical curve, and that is expected rather than a
discrepancy. Classically the projectile is an *external perturbation*: it is never
in the electronic ledger, and the CAP-free z-open box lets it leave without taking
ledger energy with it, so `plateau - E_GS` is purely the slab's gain. Here the
wavepacket **is** part of the system and the CAP **removes** it, so
`E_total(t_final) - E_GS` is only what is still in the box once everything the CAP
absorbed has already been subtracted. It is a **lower bound** on the deposit.

The raw-vs-corrected contrast is worth noting on its own: on INQ's uncorrected
ledger this estimator is velocity-INDEPENDENT at ~2.44 eV/Bohr — unphysical, since
S must fall with v here — because the norm-divided kinetic term dominates the
residual. After the correction it falls monotonically, as it must. Both are in the
CSV (`S_deposit_raw`, `S_deposit_corrected`).

## Reading this comparison honestly

The classical S is an **energy-deposit** measure, `(E_total(plateau) − E_GS)/L_slab`,
read after the projectile has fully left the box. The quantum S is
`−dT/ds` of the projectile orbital over a short early window. They address the
same physics but are not the same estimator, and the quantum one is additionally
limited by the transverse-overlap window, CAP attrition, and momentum aliasing
documented in each run notebook. Differences between the two curves should not be
attributed to quantum effects without first accounting for those three.
""")]
    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python",
                                 "name": "python3"}
    return nb


# --------------------------------------------------------------------------- #
def build_sigma_sweep() -> nbf.NotebookNode:
    """
    The cross-campaign figure the sigma sweep was run for (user instruction
    2026-07-31): S(v) by the localised-jellium DEPOSIT definition,

        S = (E_total(t_final) - E_GS) / L_slab_z,

    with one trace per sigma_WP in {0.5, 2.0, 3.0}.
    """
    cells = [md(r"""# Stopping power vs projectile localisation — **σ sweep**

$$S \;=\; \frac{E_{\rm total}(t_{\rm final}) - E_{\rm GS}}{L_{{\rm slab},z}}
\qquad L_{{\rm slab},z} = 25\ \mathrm{Bohr}$$

The localised-jellium **deposit** definition, evaluated on three wavepacket
campaigns that differ in exactly one parameter: the wavepacket width σ_WP.
Everything else — ground state, dx = 0.40, launch z = −24, the two 12.5 Bohr CAPs
with η = −1, dt = 0.04, step counts, and the four velocities — is held fixed.

## Why σ is the interesting axis

σ_WP sets how fast the packet stops being a localised object. A free Gaussian's
density width grows as

$$s(t) \;=\; \sqrt{\tfrac{\sigma^2}{2} + \tfrac{t^2}{2\sigma^2}},$$

so the asymptotic spreading rate is $1/(\sqrt2\,\sigma)$: **1.414** Bohr/a.u. at
σ = 0.5, **0.354** at σ = 2, **0.236** at σ = 3. Launched at z = −24 with the slab
face at −12.5, the σ = 0.5 packet is already 4.7–8.1 Bohr wide when it arrives and
15–26 Bohr wide when it leaves — wider than the 25 Bohr slab. The σ = 2 and 3
packets arrive at 1.8–2.5 Bohr and leave at 3.2–6.6 Bohr.

σ also fixes the momentum width, $\sigma_p = 1/(\sqrt2\,\sigma)$, hence the
zero-point kinetic content $T_1 - T_2 = 3/(4\sigma^2)$ = **81.6 / 5.10 / 2.27 eV**
for σ = 0.5 / 2 / 3. A broad packet is closer to a plane wave of definite momentum,
which is the limit in which "stopping power at velocity v" is even well posed.

## What this estimator does and does not measure

It is the same formula the classical benchmark uses, but the two are **not** the
same measurement, and the difference matters more here than anywhere else in the
campaign. Classically the projectile is an external perturbation, never in the
electronic ledger, so `plateau − E_GS` is purely the slab's energy gain. Here the
wavepacket **is** part of the system and the CAP removes it, so
`E_total(t_final) − E_GS` is *what is left in the box* once everything the CAP
absorbed has been subtracted: a **lower bound** on the deposit.

Both the raw and the norm-corrected ledgers are plotted. The correction removes
INQ's norm-divided kinetic term (`energy.hpp`), which otherwise inflates the
residual as the CAP eats the packet — on the raw ledger at σ = 0.5 this estimator
comes out velocity-*independent* at ~2.44 eV/Bohr, which no stopping power can be.
"""),
             code(f"""import sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, "{HERE}")
import wp_hd_stopping as W

rows = []
for s in W.SIGMAS:
    if not W.has_campaign(s):
        print(f"sigma={{s}}: no results on disk — skipped")
        continue
    W.set_campaign(s)
    for v in W.VELOCITIES:
        try:
            dep = W.deposit_stopping(v, sigma=s)
        except (FileNotFoundError, KeyError, IndexError) as e:
            print(f"sigma={{s}} v={{v}}: MISSING ({{type(e).__name__}})"); continue
        rows.append(dict(sigma=s, v=v,
                         s0_bohr=s/np.sqrt(2),
                         s_at_entry=float(W.sigma_d((abs(W.LAUNCH_Z)-W.SLAB_HALF)/v, s)),
                         s_at_exit=float(W.sigma_d((abs(W.LAUNCH_Z)+W.SLAB_HALF)/v, s)),
                         **dep))
d_all = pd.DataFrame(rows)

# COMPLETENESS GATE. A run that is still propagating (or was killed) yields a
# mid-flight E_total that looks like a perfectly plausible deposit but is not one:
# the packet has not deposited and the CAP has not removed it. Incomplete points
# are reported, then EXCLUDED from the figure.
incomplete = d_all[~d_all.complete]
if len(incomplete):
    print("EXCLUDED — run not finished:")
    display(incomplete[["sigma", "v", "steps_done", "steps_target",
                        "t_final_au", "norm_final", "S_deposit_corrected"]].round(3))
d = d_all[d_all.complete].copy()
print(f"{{len(d)}} of {{len(d_all)}} points complete")
display(d.round(4))
d_all.to_csv("sigma_sweep_S_deposit.csv", index=False)
print("wrote sigma_sweep_S_deposit.csv (all points, with a `complete` column)")"""),
             md("""## The figure

Solid = norm-corrected ledger (the physical trace). Dotted = raw INQ ledger, kept
visible because its flatness is the diagnostic that motivated the correction."""),
             code("""MK = {0.5: ("C0", "o"), 2.0: ("C1", "s"), 3.0: ("C2", "^")}

fig, ax = plt.subplots(figsize=(8.5, 5.8))
for s, g in d.groupby("sigma"):
    c, m = MK.get(s, ("C7", "x"))
    g = g.sort_values("v")
    ax.plot(g.v, g.S_deposit_corrected, color=c, marker=m, lw=2,
            label=f"$\\\\sigma_{{WP}}$ = {s:g} Bohr (corrected)")
    ax.plot(g.v, g.S_deposit_raw, color=c, marker=m, lw=1, ls=":", alpha=.55,
            label=f"$\\\\sigma_{{WP}}$ = {s:g} Bohr (raw ledger)")

try:
    cl = W.classical_reference()
    ax.plot(cl.v, cl.S_eV_per_Bohr, "k--", lw=1.5, alpha=.7,
            label="classical benchmark ($\\\\sigma$=0.5, deposit)")
except Exception as e:
    print("classical reference unavailable:", e)

ax.set(xlabel="v (a.u.)", ylabel="S (eV/Bohr)",
       title="Deposit stopping power $(E_{tot}(t_f)-E_{GS})/L_{slab}$ vs projectile width")
ax.grid(alpha=.3); ax.legend(fontsize=8, ncol=2)
fig.tight_layout(); fig.savefig("sigma_sweep_S_deposit.png", dpi=140); plt.show()
print("wrote sigma_sweep_S_deposit.png")"""),
             md("""## Companion: how wide the packet actually is when it does the physics

The deposit curve above is only interpretable alongside this — it is the reason
the σ campaigns were run. `s_at_entry` is the density standard deviation when the
centroid reaches the slab face at z = −12.5; `s_at_exit` when it reaches +12.5."""),
             code("""fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for s, g in d.groupby("sigma"):
    c, m = MK.get(s, ("C7", "x"))
    g = g.sort_values("v")
    ax[0].plot(g.v, g.s_at_entry, color=c, marker=m, label=f"$\\\\sigma$={s:g}")
    ax[1].plot(g.v, g.s_at_exit, color=c, marker=m, label=f"$\\\\sigma$={s:g}")
ax[0].set(xlabel="v (a.u.)", ylabel="density std at slab entry (Bohr)",
          title="Width on arrival")
ax[1].set(xlabel="v (a.u.)", ylabel="density std at slab exit (Bohr)",
          title="Width on departure")
ax[1].axhline(2*W.SLAB_HALF, ls="--", c="k", lw=1)
ax[1].annotate("slab thickness (25 Bohr)", (2.05, 2*W.SLAB_HALF*1.02), fontsize=8)
for a in ax: a.grid(alpha=.3); a.legend()
fig.tight_layout(); fig.savefig("sigma_sweep_widths.png", dpi=140); plt.show()
print("wrote sigma_sweep_widths.png")"""),
             md("""## Reading this honestly — three caveats

1. **The σ = 3 packet starts inside the −z CAP.** Launch is z = −24 and the CAP
   inner edge is z = −30, so the 6 Bohr clearance is 4.2 density-std at σ = 2 but
   only 2.8 at σ = 3, where ≈0.23 % of the packet is absorbed immediately. That is
   a one-off norm loss, not a stopping effect; it is reproduced by the matching
   vacuum control (`vac_s3p0_*`) and should be checked there before any σ = 3 vs
   σ = 2 difference is called physics.
2. **This is a lower bound, and the bound is σ-dependent.** How much of the packet
   the CAP has removed by t_final differs between campaigns, so part of any σ
   trend is CAP bookkeeping rather than deposition. The vacuum controls at each σ
   are what separate the two.
3. **No classical twin exists at σ = 2 or 3.** The classical benchmark was run at
   σ_pot = 0.354 = 0.5/√2 only. The dashed classical curve is orientation, not a
   like-for-like comparison, for the two wider campaigns.

A deposit measure genuinely comparable to the classical one would be the bath's
own internal energy change, ΔE_SS + ΔE_SB from the pairwise ledger — available
every step in `interactions.csv` for all three campaigns, and not yet built.
""")]
    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python",
                                 "name": "python3"}
    return nb


def execute(nb: nbf.NotebookNode, out: Path) -> int:
    nbf.write(nb, str(out))
    print(f"wrote {out} ({len(nb.cells)} cells); executing...")
    NotebookClient(nb, timeout=3600, kernel_name="python3",
                   resources={"metadata": {"path": str(HERE)}},
                   allow_errors=True).execute()
    nbf.write(nb, str(out))
    n_err = sum(1 for c in nb.cells if c.cell_type == "code"
                for o in c.get("outputs", []) if o.get("output_type") == "error")
    print(f"  executed with {n_err} error(s) -> {out}")
    return n_err


def _campaign_tag(sigma: float, launch_z: float = W.LAUNCH_Z_FAR) -> str:
    """'' far+0.5, 's2p0' far+2.0, 'nl' near+0.5, 'nl_s2p0' near+2.0."""
    return (W.launch_tag(launch_z) + W.sigma_tag(sigma)).rstrip("_")


def _synth_name(sigma: float, launch_z: float = W.LAUNCH_Z_FAR) -> str:
    tag = _campaign_tag(sigma, launch_z)
    return "synthesis.ipynb" if not tag else f"synthesis_{tag}.ipynb"


def main() -> int:
    """
    Targets:
      all         every campaign present on disk: run notebooks + per-campaign
                  synthesis + the cross-sigma sweep figure
      synthesis   syntheses + sweep only (no per-run notebooks)
      sweep       the cross-sigma S_deposit figure only
      <v>         one velocity of the sigma = 0.5 campaign, e.g. 3.0
      <sigma>:<v> one velocity of one campaign,       e.g. 2.0:3.0
    """
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    arg, errs = sys.argv[1], 0

    launch_z = W.LAUNCH_Z_FAR
    if arg.startswith("nl:"):                        # near-launch campaign
        launch_z, arg = W.LAUNCH_Z_NEAR, arg[3:]

    if ":" in arg:                                   # <sigma>:<v>
        s_str, v_str = arg.split(":", 1)
        sigma, v = float(s_str), float(v_str)
        W.set_campaign(sigma); W.set_launch(launch_z)
        return 1 if execute(build_run(v, sigma, launch_z),
                            HERE / f"run_{W.name_for(v)}.ipynb") else 0

    if arg == "sweep":
        return 1 if execute(build_sigma_sweep(), HERE / "sigma_sweep.ipynb") else 0

    if arg not in ("all", "synthesis"):              # bare velocity -> sigma 0.5
        v = float(arg)
        W.set_campaign(0.5); W.set_launch(launch_z)
        return 1 if execute(build_run(v, 0.5, launch_z),
                            HERE / f"run_{W.name_for(v)}.ipynb") else 0

    # Campaign-wide targets. Only campaigns with results on disk are built, so a
    # partially-finished sweep still yields everything it can (the notebook stage
    # is chained afterANY for exactly this reason).
    # A campaign is a (launch distance, sigma) pair. The far-launch sigma sweep
    # and the near-launch effective-sigma test are both enumerated here, so the
    # chained notebook stage picks up whichever have results without being told.
    campaigns = [(lz, s) for lz in (W.LAUNCH_Z_FAR, W.LAUNCH_Z_NEAR)
                 for s in W.SIGMAS if W.has_campaign(s, lz)]
    print("campaigns with results on disk: "
          + (", ".join(f"z={lz:g},sigma={s:g}" for lz, s in campaigns) or "NONE"))
    for launch_z, sigma in campaigns:
        if arg == "all":
            W.set_campaign(sigma)
            W.set_launch(launch_z)
            for v in W.VELOCITIES:
                errs += execute(build_run(v, sigma, launch_z),
                                HERE / f"run_{W.name_for(v)}.ipynb")
        errs += execute(build_synthesis(sigma, launch_z),
                        HERE / _synth_name(sigma, launch_z))

    # The cross-sigma figure is the point of the sweep — always attempted last, so
    # it sees every campaign the run above just refreshed.
    errs += execute(build_sigma_sweep(), HERE / "sigma_sweep.ipynb")
    return 1 if errs else 0


if __name__ == "__main__":
    raise SystemExit(main())
