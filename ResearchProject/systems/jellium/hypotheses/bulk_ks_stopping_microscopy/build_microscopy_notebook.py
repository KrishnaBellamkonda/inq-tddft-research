#!/usr/bin/env python3
"""Build (and execute) the classical-vs-wavepacket MICROSCOPY notebook.

    venv/bin/python build_microscopy_notebook.py [--stride 8]

One notebook, one projectile velocity (100 eV, k0 = 2.7111), covering EVERY
bulk-jellium density run made so far. Each density gets its own self-contained
section; a final section compares them.

Plan: docs/plans/bulk-jellium-ks-stopping.md
Handover: docs/handovers/bulk-jellium-ks-stopping.md

The notebook narrates; `microscopy.py` and `ks_stopping.py` do the arithmetic,
so this notebook and the per-run notebooks cannot disagree about a number.
"""
from __future__ import annotations

import sys
from pathlib import Path

import nbformat as nbf
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
SYSTEM = REPO / "ResearchProject/systems/jellium"
KS_HYP = SYSTEM / "hypotheses/bulk_ks_stopping"

# Every density pair made at 100 eV, sigma_WP = 2. Ordered dilute -> dense.
DENSITIES = [
    dict(key="rs5p7", label="r_s = 5.702", base="scripts/bulk_ks_stopping",
         LX=46.0, LZ=80.0, N_E=218, DX=0.40, RS=5.702, N0=1.287807e-3,
         HW_P=3.462, VF=0.3366, TP=49.39, Z0=-32.0, FIT_T1=18.4270),
    dict(key="rs4", label="r_s = 3.987", base="scripts/bulk_ks_stopping_rs4",
         LX=40.0, LZ=80.0, N_E=482, DX=0.50, RS=3.987, N0=3.765625e-3,
         HW_P=5.919, VF=0.4813, TP=28.88, Z0=-32.0, FIT_T1=18.4270),
]
COMMON_T0, COMMON_T1 = 4.0, 9.37   # the cross-sigma / cross-density window

# Every (density, sigma) pair. z0 follows the boundary rule -L_z/2 + 4 sigma, so
# it DIFFERS per sigma — which is exactly why the comparison must be made on a
# common TIME window rather than on raw path.
SIGMA_RUNS = [
    dict(bath="r_s=5.702", sigma=1.0, base="scripts/bulk_ks_stopping_sigma1",
         z0=-36.0, n0=1.287807e-3, LZ=80.0),
    dict(bath="r_s=5.702", sigma=2.0, base="scripts/bulk_ks_stopping",
         z0=-32.0, n0=1.287807e-3, LZ=80.0),
    dict(bath="r_s=5.702", sigma=3.0, base="scripts/bulk_ks_stopping_sigma3",
         z0=-28.0, n0=1.287807e-3, LZ=80.0),
    dict(bath="r_s=3.987", sigma=1.0, base="scripts/bulk_ks_stopping_rs4_sigma1",
         z0=-36.0, n0=3.765625e-3, LZ=80.0),
    dict(bath="r_s=3.987", sigma=2.0, base="scripts/bulk_ks_stopping_rs4",
         z0=-32.0, n0=3.765625e-3, LZ=80.0),
    dict(bath="r_s=3.987", sigma=3.0, base="scripts/bulk_ks_stopping_rs4_sigma3",
         z0=-28.0, n0=3.765625e-3, LZ=80.0),
]


def md(t): return nbf.v4.new_markdown_cell(t)
def code(s): return nbf.v4.new_code_cell(s)


def build(stride: int) -> nbf.NotebookNode:
    c: list = []

    c.append(md(r"""# Why does the classical projectile slow down so much faster?

**One velocity (100 eV, $k_0 = 2.711$ a.u.), every bulk-jellium density, both
projectile representations.**

Two completed twin pairs agree on the *size* of the effect and say nothing about
its cause:

| bath | $S_{\rm classical}$ | $S_2$ (WP drift) | ratio |
|---|---|---|---|
| $r_s = 5.702$ | 0.375 | 0.058 | **6.5×** |
| $r_s = 3.987$ | 0.890 | 0.157 | **5.7×** |

A 2.92× density increase moved that ratio by only 13 %, which ruled the density
lever out as the explanation. This notebook stops fitting slopes and looks at
what the electron gas actually *does* around each projectile.

**What drag is, microscopically.** A charged projectile polarises the gas. The
induced density does not sit symmetrically about it — the gas responds at a
finite rate (set by $\omega_p$) while the projectile keeps moving, so the
polarisation cloud **lags**. More induced charge ends up behind than in front,
and the induced field of that imbalance pulls backwards. Stopping power *is* that
asymmetry.

So "why is the classical projectile stopped harder?" splits into two measurables:

1. **Coupling strength** — does it induce a *larger* $\Delta n$?
2. **Lag** — is its $\Delta n$ more *asymmetric* front-to-back?

Both are computed below, in the projectile's own moving frame."""))

    c.append(code(f"""import sys, math
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

HERE   = Path("{HERE}")
SYSTEM = Path("{SYSTEM}")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path("{KS_HYP}")))
sys.path.insert(0, str(Path("{REPO}") / "inq-stack/python"))

import microscopy as M
import ks_stopping as K
from ks_stopping import HA_TO_EV
from inqview import load_vti

DENSITIES = {DENSITIES!r}
STRIDE = {stride}
COMMON_T0, COMMON_T1 = {COMMON_T0}, {COMMON_T1}
DT = 0.04

plt.rcParams.update({{"figure.dpi": 120, "axes.grid": True, "grid.alpha": 0.25,
                     "axes.spines.top": False, "axes.spines.right": False}})
(HERE / "figures").mkdir(exist_ok=True)

# Load everything once. VTI reads dominate, so stride the frames.
DATA = {{}}
for d in DENSITIES:
    base = SYSTEM / d["base"]
    e = {{}}
    e["cl"]   = K.load_classical_run(base / "classical", box_length_z=d["LZ"])
    e["wp"]   = K.load_wp_run(base / "wp", box_length_z=d["LZ"], z0=d["Z0"])
    e["en_cl"] = K.load_energies(base / "classical")
    e["en_wp"] = K.load_energies(base / "wp")
    e["base"] = base
    DATA[d["key"]] = e
    print(f"{{d['label']}}: classical {{len(e['cl'].t)}} steps, wp {{len(e['wp'].t)}} steps")
print("\\nscalar observables loaded; density frames are read per-section below")"""))

    c.append(md(r"""---
## Conventions, and what each panel is computed from

| symbol | meaning | source |
|---|---|---|
| $n_0$ | uniform bath density | $N/V$ |
| $\Delta n$ | induced density, $n(t)-n(0)$ | derived from `density_total` |
| $\zeta$ | distance from the projectile, $z - z_{\rm proj}$, wrapped | derived |
| $T_1,\;T_2$ | $\langle p^2\rangle/2m$, $\langle p\rangle^2/2m$ | `wp_momentum_stats.csv` |
| $s_3$ | circular (periodic-safe) centroid | `wp_real_space_stats.csv` |
| $z_{\rm proj}$ | classical track / WP circular centroid | `electron_track.csv` |

**One decomposition matters more than any plot here.** The induced density is
$\Delta n = n_{\rm total}(t)-n_{\rm total}(0)$, formed here by subtracting the
$t=0$ frame of `density_total`. (Runs used to also write this as a separate
`density_delta` field; it was a byte-for-byte duplicate at ~10 GB per run and is
no longer produced.) For the **classical** run that *is* the bath
response, because the projectile is an external potential and never enters $n$.

For the **wavepacket** run it is not — the WP is an occupied orbital, so
`density_delta` contains the packet itself moving (a large positive blob where it
now is, negative where it launched) on top of the bath response. Comparing that
with the classical `density_delta` would show the wavepacket inducing a far
*larger* response, which is an artefact of comparing a projectile against a
polarisation cloud.

So every WP panel below uses the **explicit bath**:

$$n_{\rm bath}(t) = n_{\rm total}(t) - n_{\rm WP}(t), \qquad
\Delta n_{\rm bath}(t) = n_{\rm bath}(t) - n_{\rm bath}(0)$$

**Transverse reduction.** Profiles are averaged over a cylinder of radius 3 Bohr
about the trajectory, not taken on a single grid line — a one-line cut through a
3-D density is noisy and grid-orientation dependent, whereas the charge in a tube
around the path is what the drag integrates over.

**Periodicity.** All $\zeta$ profiles wrap into $(-L_z/2, +L_z/2]$, and the WP
position is the *circular* centroid throughout — the naive $\int z|\psi|^2$
slides to a wrong value near a cell face."""))

    for d in DENSITIES:
        k, lab = d["key"], d["label"]
        c.append(md(f"""---
---
# Section: **{lab}** bath

`{d['base']}` — {d['LX']:.0f} × {d['LX']:.0f} × {d['LZ']:.0f} Bohr, N = {d['N_E']},
dx = {d['DX']}, $n_0$ = {d['N0']:.4e} e/Bohr³.
$\\hbar\\omega_p$ = {d['HW_P']} eV, $v_F$ = {d['VF']}, plasma period {d['TP']} a.u.

Projectile: 100 eV electron, $k_0 = 2.711$, $\\sigma_{{\\rm WP}} = 2$, launched at
$z_0 = {d['Z0']:.0f}$. The classical twin presents the identical charge cloud
$\\exp(-r^2/\\sigma^2)$ via its Gaussian UPF."""))

        # ---- energies ---------------------------------------------------
        c.append(md(f"""## {lab} — 1. Decomposed energies

Every component INQ tracks, both halves. Absolute values first, then the
**change** from $t=0$, which is what says where the deposited energy went.

Remember the asymmetry in the bookkeeping: the WP *is* inside `energy_total`
(closed system ⇒ total must be constant), while for the classical run
`energy_ion_kinetic` is left at zero by INQ, so `energy_total` is the
**electronic** energy alone and is supposed to **rise** by the projectile's loss."""))

        c.append(code(f"""d = [x for x in DENSITIES if x["key"] == "{k}"][0]
E = DATA["{k}"]
fig, ax = plt.subplots(2, 2, figsize=(13, 8))
for col, (half, en) in enumerate((("classical", E["en_cl"]), ("wavepacket", E["en_wp"]))):
    t = en["time_au"].to_numpy()
    for cname in M.ENERGY_COLS:
        if cname not in en.columns: continue
        v = en[cname].to_numpy() * HA_TO_EV
        if np.all(v == 0.0): continue
        ax[0, col].plot(t, v, lw=1.2, label=cname.replace("energy_", ""))
    ax[0, col].set_title(f"{{half}} — components (absolute)")
    ax[0, col].set_ylabel("E [eV]"); ax[0, col].legend(fontsize=7, frameon=False, ncol=2)

    ch = M.energy_channels(en)
    for cname, v in ch.items():
        if cname == "time_au": continue
        ax[1, col].plot(ch["time_au"], v, lw=1.4,
                        ls="--" if cname == "energy_total" else "-",
                        label=cname.replace("energy_", ""))
    ax[1, col].axhline(0, color="k", lw=0.8)
    ax[1, col].axvspan(COMMON_T0, COMMON_T1, color="0.85", zorder=0)
    ax[1, col].set_title(f"{{half}} — CHANGE from t=0")
    ax[1, col].set_xlabel("time [a.u.]"); ax[1, col].set_ylabel(r"$\\Delta E$ [eV]")
    ax[1, col].legend(fontsize=7, frameon=False, ncol=2)
fig.suptitle(f"{{d['label']}} — energy decomposition (shaded = common window)", y=1.00)
fig.tight_layout(); fig.savefig(HERE/"figures"/f"energies_{k}.png", dpi=150, bbox_inches="tight")
plt.show()

ch = M.energy_channels(E["en_cl"])
print("CLASSICAL — where the energy goes, over the whole run (eV):")
for cname, v in sorted(ch.items(), key=lambda kv: -abs(kv[1][-1]) if kv[0] != "time_au" else 0):
    if cname == "time_au": continue
    print(f"   {{cname:26s}} {{v[-1]:+9.3f}}")
print(f"\\n   projectile KE loss        {{(E['cl'].T[0]-E['cl'].T[-1])*HA_TO_EV:+9.3f}}")
chw = M.energy_channels(E["en_wp"])
print("\\nWAVEPACKET — total must stay constant (closed system):")
print(f"   energy_total drift        {{chw['energy_total'][-1]:+9.2e}} eV")"""))

        # ---- KS orbital / projectile profile ----------------------------
        c.append(md(f"""## {lab} — 2. The projectile itself: KS orbital vs classical charge

The classical projectile is a **fixed** Gaussian ($\\sigma_{{\\rm charge}}=1.414$
Bohr, unchanging). The wavepacket KS orbital **spreads**. This panel is the
clearest statement of what physically differs between the two runs."""))

        c.append(code(f"""E = DATA["{k}"]; d = [x for x in DENSITIES if x["key"] == "{k}"][0]
wp = E["wp"]
fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.8))
m = (wp.t >= COMMON_T0) & (wp.t <= COMMON_T1)
ax[0].plot(wp.t, wp.sigma_z, lw=1.6, color="C0", label=r"WP $\\sigma_z$ (circular)")
sd_free = np.sqrt(2.0**2/2 + wp.t**2/(2*2.0**2))
ax[0].plot(wp.t, sd_free, ls=":", color="k", lw=1.2, label=r"free Gaussian $\\sigma_d(t)$")
ax[0].axhline(2.0/np.sqrt(2), ls="--", color="C3", lw=1.2, label="classical (fixed) 1.414")
ax[0].axvspan(COMMON_T0, COMMON_T1, color="0.85", zorder=0)
ax[0].set_xlabel("time [a.u.]"); ax[0].set_ylabel(r"width [Bohr]")
ax[0].set_title("projectile width"); ax[0].legend(fontsize=7, frameon=False)

ax[1].plot(wp.t, wp.norm, lw=1.5, color="C2")
ax[1].set_xlabel("time [a.u.]"); ax[1].set_ylabel("WP orbital norm")
ax[1].set_title("norm (leakage into the bath)")

ax[2].plot(wp.t, wp.localisation_energy, lw=1.5, color="C3")
ax[2].axhline(3/(4*2.0**2)*HA_TO_EV, ls=":", color="k", lw=1)
ax[2].axvspan(COMMON_T0, COMMON_T1, color="0.85", zorder=0)
ax[2].set_xlabel("time [a.u.]"); ax[2].set_ylabel(r"$T_1-T_2$ [eV]")
ax[2].set_title("localisation energy")
fig.suptitle(f"{{d['label']}} — the projectile's own state", y=1.03)
fig.tight_layout(); fig.savefig(HERE/"figures"/f"projectile_{k}.png", dpi=150, bbox_inches="tight")
plt.show()
print(f"WP sigma_z: {{wp.sigma_z[0]:.3f}} -> {{wp.sigma_z[-1]:.3f}} Bohr "
      f"(x{{wp.sigma_z[-1]/wp.sigma_z[0]:.2f}});  classical stays 1.414 throughout")"""))

        # ---- kinematics -------------------------------------------------
        c.append(md(f"""## {lab} — 3. Position and kinetic energy

Classical: $z(t)$, $v_z(t)$, $\\tfrac12 m v^2$ from the Ehrenfest track.
Wavepacket: the circular centroid $s_3$, the integrated-momentum position $s_4$,
and both KE definitions."""))

        c.append(code(f"""E = DATA["{k}"]; d = [x for x in DENSITIES if x["key"] == "{k}"][0]
cl, wp = E["cl"], E["wp"]
fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.8))
ax[0].plot(cl.t, cl.z, lw=1.6, color="C3", label="classical $z$")
ax[0].plot(wp.t, wp.s3, lw=1.6, color="C0", label="WP $s_3$ (circular)")
ax[0].plot(wp.t, wp.s4, lw=1.0, ls="--", color="C1", label=r"WP $s_4=\\int\\langle p\\rangle dt$")
ax[0].axvspan(COMMON_T0, COMMON_T1, color="0.85", zorder=0)
ax[0].set_xlabel("time [a.u.]"); ax[0].set_ylabel("z [Bohr]")
ax[0].set_title("position"); ax[0].legend(fontsize=7, frameon=False)

ax[1].plot(cl.t, cl.vz, lw=1.6, color="C3", label="classical $v_z$")
ax[1].plot(wp.t, wp.pz, lw=1.6, color="C0", label=r"WP $\\langle p_z\\rangle$")
ax[1].axvspan(COMMON_T0, COMMON_T1, color="0.85", zorder=0)
ax[1].set_xlabel("time [a.u.]"); ax[1].set_ylabel("velocity [a.u.]")
ax[1].set_title("velocity (m = 1 so p = v)"); ax[1].legend(fontsize=7, frameon=False)

ax[2].plot(cl.t, (cl.T-cl.T[0])*HA_TO_EV, lw=1.6, color="C3", label="classical")
ax[2].plot(wp.t, (wp.T1-wp.T1[0])*HA_TO_EV, lw=1.6, color="C0", label="WP $T_1$")
ax[2].plot(wp.t, (wp.T2-wp.T2[0])*HA_TO_EV, lw=1.6, color="C1", label="WP $T_2$")
ax[2].axvspan(COMMON_T0, COMMON_T1, color="0.85", zorder=0)
ax[2].axhline(0, color="k", lw=0.8)
ax[2].set_xlabel("time [a.u.]"); ax[2].set_ylabel(r"$\\Delta$KE [eV]")
ax[2].set_title("kinetic-energy loss"); ax[2].legend(fontsize=7, frameon=False)
fig.suptitle(f"{{d['label']}} — kinematics", y=1.03)
fig.tight_layout(); fig.savefig(HERE/"figures"/f"kinematics_{k}.png", dpi=150, bbox_inches="tight")
plt.show()

mc = (cl.t >= COMMON_T0) & (cl.t <= COMMON_T1); mw = (wp.t >= COMMON_T0) & (wp.t <= COMMON_T1)
print(f"over the common window [{{COMMON_T0}}, {{COMMON_T1}}] a.u.:")
print(f"   classical  KE loss {{(cl.T[mc][0]-cl.T[mc][-1])*HA_TO_EV:7.3f}} eV over "
      f"{{cl.z[mc][-1]-cl.z[mc][0]:6.2f}} Bohr")
print(f"   WP  T2     loss    {{(wp.T2[mw][0]-wp.T2[mw][-1])*HA_TO_EV:7.3f}} eV over "
      f"{{wp.s3[mw][-1]-wp.s3[mw][0]:6.2f}} Bohr")"""))

        # ---- the wake ---------------------------------------------------
        c.append(md(f"""## {lab} — 4. The induced density: the wake, in the projectile frame

This is the microscopic answer. Both panels show $\\Delta n$ against
$\\zeta = z - z_{{\\rm proj}}$, so the projectile sits at $\\zeta = 0$ and the
polarisation cloud is seen in its own frame. **Left of 0 is behind.**

The classical panel uses $\Delta n$ derived from `density_total`. The wavepacket panel uses
$\\Delta n_{{\\rm bath}} = (n_{{\\rm total}}-n_{{\\rm WP}}) -$ its $t=0$ value,
for the reason given in the conventions section."""))

        c.append(code(f"""d = [x for x in DENSITIES if x["key"] == "{k}"][0]
E = DATA["{k}"]; base = E["base"]
print("reading density frames (this is the slow cell)...", flush=True)
cl_s = M.induced_series(base/"classical", dt_au=DT,
                        stride=M.stride_for(base/"classical","density_total",40),
                        radius_bohr=3.0, loader=load_vti)
wp_s = M.bath_series(base/"wp", dt_au=DT,
                     stride=M.stride_for(base/"wp","density_total",40),
                     radius_bohr=3.0, loader=load_vti)
zp_cl = np.interp(cl_s.t, E["cl"].t, E["cl"].z)
zp_wp = np.interp(wp_s.t, E["wp"].t, E["wp"].s3)
zeta_c, V_c = M.to_projectile_frame(cl_s, zp_cl, d["LZ"])
zeta_w, V_w = M.to_projectile_frame(wp_s, zp_wp, d["LZ"])
A_c = M.wake_asymmetry(zeta_c, V_c); A_w = M.wake_asymmetry(zeta_w, V_w)
print(f"   classical {{V_c.shape[0]}} frames, wavepacket {{V_w.shape[0]}} frames")

vmax = max(np.abs(V_c).max(), np.abs(V_w).max())
fig, ax = plt.subplots(2, 2, figsize=(13, 8))
for col, (ttl, zt, V) in enumerate((("classical", zeta_c, V_c),
                                    ("wavepacket (bath only)", zeta_w, V_w))):
    im = ax[0, col].pcolormesh(zt, np.arange(V.shape[0]), V, cmap="RdBu_r",
                               vmin=-vmax, vmax=vmax, shading="auto")
    ax[0, col].set_yticks(np.linspace(0, V.shape[0]-1, 6))
    ax[0, col].set_yticklabels([f"{{x:.0f}}" for x in np.linspace(
        (cl_s.t if col == 0 else wp_s.t)[0], (cl_s.t if col == 0 else wp_s.t)[-1], 6)])
    ax[0, col].axvline(0, color="k", lw=1.0, ls="--")
    ax[0, col].set_xlim(-20, 20); ax[0, col].set_xlabel(r"$\\zeta = z - z_{{\\rm proj}}$ [Bohr]")
    ax[0, col].set_ylabel("time [a.u.]"); ax[0, col].set_title(f"{{ttl}}: $\\Delta n(\\zeta, t)$")
    fig.colorbar(im, ax=ax[0, col], label=r"$\\Delta n$")

for col, (ttl, zt, V, tt) in enumerate((("classical", zeta_c, V_c, cl_s.t),
                                        ("wavepacket (bath)", zeta_w, V_w, wp_s.t))):
    sel = [i for i in range(len(tt)) if COMMON_T0 <= tt[i] <= COMMON_T1]
    sel = sel[::max(1, len(sel)//4)] if sel else [len(tt)//2]
    for i in sel:
        ax[1, col].plot(zt, V[i], lw=1.3, label=f"t={{tt[i]:.1f}}")
    ax[1, col].axvline(0, color="k", lw=1.0, ls="--")
    ax[1, col].axhline(0, color="k", lw=0.8)
    ax[1, col].set_xlim(-20, 20); ax[1, col].set_ylim(-vmax, vmax)
    ax[1, col].set_xlabel(r"$\\zeta$ [Bohr]"); ax[1, col].set_ylabel(r"$\\Delta n$")
    ax[1, col].set_title(f"{{ttl}}: profiles in the window")
    ax[1, col].legend(fontsize=7, frameon=False)
fig.suptitle(f"{{d['label']}} — induced density in the projectile frame "
             f"(shared colour scale; behind = $\\zeta<0$)", y=1.00)
fig.tight_layout(); fig.savefig(HERE/"figures"/f"wake_{k}.png", dpi=150, bbox_inches="tight")
plt.show()

np.savez(HERE/f"wake_{k}.npz", zeta_c=zeta_c, V_c=V_c, t_c=cl_s.t,
         zeta_w=zeta_w, V_w=V_w, t_w=wp_s.t,
         asym_c=A_c["asymmetry"], asym_w=A_w["asymmetry"],
         atproj_c=A_c["at_projectile"], atproj_w=A_w["at_projectile"],
         depl_c=A_c["peak_depletion"], depl_w=A_w["peak_depletion"])
print(f"n0 = {{d['N0']:.3e}}")
print(f"   classical : peak |dn| = {{np.abs(V_c).max():.3e}} = {{np.abs(V_c).max()/d['N0']*100:5.2f}} % of n0")
print(f"   wavepacket: peak |dn| = {{np.abs(V_w).max():.3e}} = {{np.abs(V_w).max()/d['N0']*100:5.2f}} % of n0")"""))

        c.append(code(f"""d = [x for x in DENSITIES if x["key"] == "{k}"][0]
Z = np.load(HERE/f"wake_{k}.npz")
fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.8))
ax[0].plot(Z["t_c"], Z["atproj_c"], lw=1.6, color="C3", label="classical")
ax[0].plot(Z["t_w"], Z["atproj_w"], lw=1.6, color="C0", label="wavepacket (bath)")
ax[0].axhline(0, color="k", lw=0.8); ax[0].axvspan(COMMON_T0, COMMON_T1, color="0.85", zorder=0)
ax[0].set_xlabel("time [a.u.]"); ax[0].set_ylabel(r"$\\Delta n$ at $\\zeta=0$")
ax[0].set_title("density AT the projectile"); ax[0].legend(fontsize=7, frameon=False)

ax[1].plot(Z["t_c"], Z["depl_c"], lw=1.6, color="C3", label="classical")
ax[1].plot(Z["t_w"], Z["depl_w"], lw=1.6, color="C0", label="wavepacket (bath)")
ax[1].axvspan(COMMON_T0, COMMON_T1, color="0.85", zorder=0)
ax[1].set_xlabel("time [a.u.]"); ax[1].set_ylabel(r"min $\\Delta n$")
ax[1].set_title("peak depletion (coupling strength)"); ax[1].legend(fontsize=7, frameon=False)

ax[2].plot(Z["t_c"], Z["asym_c"], lw=1.6, color="C3", label="classical")
ax[2].plot(Z["t_w"], Z["asym_w"], lw=1.6, color="C0", label="wavepacket (bath)")
ax[2].axhline(0, color="k", lw=0.8); ax[2].axvspan(COMMON_T0, COMMON_T1, color="0.85", zorder=0)
ax[2].set_xlabel("time [a.u.]"); ax[2].set_ylabel(r"$\\int_{{-w}}^{{0}}\\Delta n - \\int_0^{{+w}}\\Delta n$")
ax[2].set_title("front/back ASYMMETRY (the drag)"); ax[2].legend(fontsize=7, frameon=False)
fig.suptitle(f"{{d['label']}} — the two drag ingredients", y=1.03)
fig.tight_layout(); fig.savefig(HERE/"figures"/f"drag_{k}.png", dpi=150, bbox_inches="tight")
plt.show()

mc = (Z["t_c"] >= COMMON_T0) & (Z["t_c"] <= COMMON_T1)
mw = (Z["t_w"] >= COMMON_T0) & (Z["t_w"] <= COMMON_T1)
print("mean over the common window:")
print(f"   peak depletion : classical {{Z['depl_c'][mc].mean():+.3e}}   "
      f"wavepacket {{Z['depl_w'][mw].mean():+.3e}}   "
      f"ratio {{Z['depl_c'][mc].mean()/Z['depl_w'][mw].mean():5.2f}}x")
print(f"   asymmetry      : classical {{Z['asym_c'][mc].mean():+.3e}}   "
      f"wavepacket {{Z['asym_w'][mw].mean():+.3e}}   "
      f"ratio {{Z['asym_c'][mc].mean()/Z['asym_w'][mw].mean():5.2f}}x")
print("\\ncompare with the stopping-power ratio for this bath "
      "(6.5x at r_s=5.702, 5.7x at r_s=3.987 on the full window)")"""))

    # ---------------- cross-density synthesis ------------------------------
    c.append(md(r"""---
---
# Cross-density synthesis

The two sections above are self-contained. This one puts the microscopic numbers
from both baths side by side against the stopping powers they are supposed to
explain."""))

    c.append(code("""rows = []
for d in DENSITIES:
    Z = np.load(HERE/f"wake_{d['key']}.npz")
    mc = (Z["t_c"] >= COMMON_T0) & (Z["t_c"] <= COMMON_T1)
    mw = (Z["t_w"] >= COMMON_T0) & (Z["t_w"] <= COMMON_T1)
    E = DATA[d["key"]]
    fw = K.fit_all_wp(E["wp"], COMMON_T0, COMMON_T1)["S_23"]
    cl = E["cl"]; mk = (cl.t >= COMMON_T0) & (cl.t <= COMMON_T1)
    xk = cl.z[mk] - cl.z[0]
    dKE = (cl.T[mk][0] - cl.T[mk]) * HA_TO_EV
    Scl = np.polyfit(xk, dKE, 1)[0]
    rows.append(dict(bath=d["label"], n0=d["N0"],
                     S_cl=Scl, S_wp=fw.S_ev_per_bohr, S_ratio=Scl/fw.S_ev_per_bohr,
                     depl_cl=Z["depl_c"][mc].mean(), depl_wp=Z["depl_w"][mw].mean(),
                     depl_ratio=Z["depl_c"][mc].mean()/Z["depl_w"][mw].mean(),
                     asym_cl=Z["asym_c"][mc].mean(), asym_wp=Z["asym_w"][mw].mean(),
                     asym_ratio=Z["asym_c"][mc].mean()/Z["asym_w"][mw].mean()))
tab = pd.DataFrame(rows)
display(tab.style.format({c: "{:.4g}" for c in tab.columns if c != "bath"}))

print("\\nDoes the microscopic asymmetry ratio explain the stopping ratio?")
for r in rows:
    print(f"  {r['bath']}: S ratio {r['S_ratio']:5.2f}x   "
          f"asymmetry ratio {r['asym_ratio']:5.2f}x   "
          f"depletion ratio {r['depl_ratio']:5.2f}x")"""))

    c.append(md(r"""## How to read the synthesis table

Three ratios, each classical ÷ wavepacket, all on the same window:

- **S ratio** — the thing to be explained.
- **depletion ratio** — how much more strongly the classical projectile
  *couples*. A compact charge reaches shorter wavelengths than a smeared one, so
  a ratio $\gg 1$ says coupling strength is doing the work.
- **asymmetry ratio** — how much more *lagged* the classical cloud is. A ratio
  $\gg 1$ with a depletion ratio near 1 would mean both projectiles polarise the
  gas equally but only the classical one leaves it behind.

If the **asymmetry** ratio tracks the **S** ratio across both baths, the drag is
being set by the lag and the microscopic picture is closed. If it does not, the
missing factor is elsewhere — most likely in how the induced field couples back
to an *extended* projectile, since the wavepacket feels a force averaged over its
own width while the classical point charge feels it at one place.

*This notebook reports the measurement. Where the ratios do not close, that is
recorded as open, not explained away.*"""))

    # ---------------- sigma dependence of the microscopic ratios -----------
    c.append(md(r"""---
---
# Width dependence: does the residual move with $\sigma$?

The cross-density synthesis above found that the stopping ratio **factorises**:

$$S_{\rm ratio} \;=\; \underbrace{\text{lag asymmetry}}_{\text{carries the density dependence}}
\;\times\; \underbrace{\approx 2.15}_{\text{density-independent residual}}$$

That residual was attributed — as an *inference*, not a measurement — to
**back-action averaging**: the classical point charge samples the induced field at
one point, while the wavepacket samples it averaged over its own width.

**If that attribution is right, the residual is geometric and must move with
$\sigma$.** Specifically it should be *larger* for the wider in-flight packet
($\sigma=1$, which disperses to 7.2 Bohr) and *smaller* for the narrower one
($\sigma=3$, 3.5 Bohr). If instead it sits at ~2.15 across all $\sigma$, the
back-action explanation is **wrong** and something else is doing the work.

This section runs that test over all six (density, $\sigma$) pairs.

**Note on $z_0$:** the launch point follows the boundary rule
$z_0 = -L_z/2 + 4\sigma$, so it differs per $\sigma$ ($-36, -32, -28$). Path
lengths therefore differ, which is why every comparison here is made on the
common **time** window $[4, 9.37]$ a.u. and never on raw path."""))

    c.append(code("""import stopping_power as SP
SIGMA_RUNS = """ + repr(SIGMA_RUNS) + """
SSTRIDE = 16          # coarser than the per-density sections: these are
                      # window-AVERAGED ratios, not shapes, so ~20 frames suffice

rows = []
for r in SIGMA_RUNS:
    base = SYSTEM / r["base"]
    tag = f"{r['bath']} s={r['sigma']:.0f}"
    try:
        wp = K.load_wp_run(base/"wp", box_length_z=r["LZ"], z0=r["z0"])
        cl = K.load_classical_run(base/"classical", box_length_z=r["LZ"])
        en = K.load_energies(base/"classical")
    except Exception as exc:
        print(f"  {tag:20s} SKIPPED — {type(exc).__name__}: {exc}")
        continue

    # stopping powers on the common window
    fw = K.fit_all_wp(wp, COMMON_T0, COMMON_T1)["S_23"]
    te = en["time_au"].to_numpy()
    dE = (en["energy_total"] - en["energy_total"].iloc[0]).to_numpy() * HA_TO_EV
    x  = np.interp(te, cl.t, cl.z) - cl.z[0]
    mE = (te >= COMMON_T0) & (te <= COMMON_T1)
    fc = SP.free_fit(x[mE], dE[mE], x[mE].min(), x[mE].max())

    # microscopic ingredients, projectile frame
    # TIME-matched sampling, not a fixed frame stride: runs written before the
    # 2026-07-31 cadence change have 301-347 frames, after it ~87. A blanket
    # stride samples them 4x differently and manufactured a spurious residual
    # of 3.32 for this very pair (true value 2.34). See M.stride_for().
    sc = M.stride_for(base/"classical", "density_total", 40)
    sw = M.stride_for(base/"wp", "density_total", 40)
    cl_s = M.induced_series(base/"classical", dt_au=DT, stride=sc,
                            radius_bohr=3.0, loader=load_vti)
    wp_s = M.bath_series(base/"wp", dt_au=DT, stride=sw,
                         radius_bohr=3.0, loader=load_vti)
    zc, Vc = M.to_projectile_frame(cl_s, np.interp(cl_s.t, cl.t, cl.z), r["LZ"])
    zw, Vw = M.to_projectile_frame(wp_s, np.interp(wp_s.t, wp.t, wp.s3), r["LZ"])
    Ac, Aw = M.wake_asymmetry(zc, Vc), M.wake_asymmetry(zw, Vw)
    mc = (cl_s.t >= COMMON_T0) & (cl_s.t <= COMMON_T1)
    mw = (wp_s.t >= COMMON_T0) & (wp_s.t <= COMMON_T1)
    mwp = (wp.t >= COMMON_T0) & (wp.t <= COMMON_T1)

    asym_r = Ac["asymmetry"][mc].mean() / Aw["asymmetry"][mw].mean()
    S_r    = fc["S"] / fw.S_ev_per_bohr
    rows.append(dict(pair=tag, bath=r["bath"], sigma=r["sigma"],
                     w_cl=r["sigma"]/np.sqrt(2), w_wp=wp.sigma_z[mwp][-1],
                     S_cl=fc["S"], S_wp=fw.S_ev_per_bohr, S_ratio=S_r,
                     depl_ratio=Ac["peak_depletion"][mc].mean()/Aw["peak_depletion"][mw].mean(),
                     asym_ratio=asym_r, residual=S_r/asym_r))
    print(f"  {tag:20s} S_ratio {S_r:5.2f}  asym {asym_r:5.2f}  "
          f"residual {S_r/asym_r:5.2f}")

tab = pd.DataFrame(rows)
display(tab.style.format({c: "{:.4g}" for c in tab.columns if c not in ("pair","bath")}))"""))

    c.append(code(r"""if len(tab) >= 2:
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0))
    for bath, mk in (("r_s=5.702", "o-"), ("r_s=3.987", "s--")):
        sub = tab[tab.bath == bath].sort_values("w_wp")
        if sub.empty: continue
        ax[0].plot(sub.w_wp, sub.S_ratio,    mk, label=bath)
        ax[1].plot(sub.w_wp, sub.asym_ratio, mk, label=bath)
        ax[2].plot(sub.w_wp, sub.residual,   mk, label=bath)
    for a, t, yl in zip(ax,
            ("stopping ratio $S_{cl}/S_{wp}$", "lag asymmetry ratio",
             "RESIDUAL = $S$ ratio / asymmetry ratio"),
            ("ratio", "ratio", "residual")):
        a.set_xlabel(r"WP in-flight width at $t=9.37$  [Bohr]")
        a.set_ylabel(yl); a.set_title(t, fontsize=9)
        a.legend(fontsize=7, frameon=False)
    ax[2].axhline(2.15, ls=":", color="k", lw=1.2)
    ax[2].annotate("2.15 (density-sweep value)", xy=(0.04, 0.06),
                   xycoords="axes fraction", fontsize=8)
    fig.suptitle("Does the back-action residual depend on packet width?", y=1.02)
    fig.tight_layout()
    fig.savefig(HERE/"figures"/"sigma_residual.png", dpi=150, bbox_inches="tight")
    plt.show()

    print("VERDICT")
    sp = tab.residual.max() / tab.residual.min() if tab.residual.min() > 0 else float("nan")
    print(f"  residual spans {tab.residual.min():.2f} - {tab.residual.max():.2f} "
          f"(factor {sp:.2f}) over WP in-flight widths "
          f"{tab.w_wp.min():.2f} - {tab.w_wp.max():.2f} Bohr")
    mono = all(
        (lambda v: all(v[i] <= v[i+1] for i in range(len(v)-1))
                or all(v[i] >= v[i+1] for i in range(len(v)-1)))(
            list(tab[tab.bath == b].sort_values("w_wp").residual))
        for b in tab.bath.unique())
    print(f"  monotonic in width within each bath? {'YES' if mono else 'NO'}")
    if sp > 1.5 and mono:
        print("  => the residual DOES move with width, monotonically: consistent")
        print("     with a geometric back-action form factor.")
    elif sp > 1.5:
        print("  => spread exceeds 1.5x but is NOT monotonic — that is scatter, not")
        print("     a form factor. Check sampling before believing a width trend.")
    else:
        print("  => the residual is FLAT in width: the back-action form-factor")
        print("     explanation is NOT supported. Something density- AND")
        print("     width-independent is suppressing the wavepacket. Report as open.")
else:
    print("not enough pairs available yet (a run may still be executing)")"""))

    c.append(md(r"""## Reading the verdict

The middle panel is the control: the lag asymmetry is a property of *how the gas
responds*, and it should track the stopping ratio if the microscopic picture is
right.

The right panel is the test. Two outcomes, both informative:

- **Residual varies with width** → it is geometric, the back-action
  interpretation stands, and the classical/WP gap is fully accounted for as
  (lag) × (extended-object averaging).
- **Residual is flat** → the inference is wrong. Something suppresses the
  wavepacket that depends on *neither* density *nor* width — which would point at
  the representation itself (self-interaction of the KS orbital, or Pauli
  exclusion against the occupied manifold), not at geometry at all.

A flat residual would be the more interesting result, and it is the one this
notebook is built to be able to report honestly."""))

    nb = nbf.v4.new_notebook(cells=c)
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python",
                                 "name": "python3"}
    return nb


def main() -> int:
    stride = 8
    if "--stride" in sys.argv:
        stride = int(sys.argv[sys.argv.index("--stride") + 1])
    nb = build(stride)
    out = HERE / "classical_vs_wavepacket_microscopy.ipynb"
    nbf.write(nb, str(out))
    print(f"wrote {out} ({len(nb.cells)} cells, stride={stride})")

    print("executing...")
    client = NotebookClient(nb, timeout=7200, kernel_name="python3",
                            resources={"metadata": {"path": str(HERE)}},
                            allow_errors=True)
    client.execute()
    nbf.write(nb, str(out))
    n_err = sum(1 for c in nb.cells if c.cell_type == "code"
                for o in c.get("outputs", []) if o.get("output_type") == "error")
    print(f"executed with {n_err} error(s) -> {out}")
    return 1 if n_err else 0


if __name__ == "__main__":
    raise SystemExit(main())
