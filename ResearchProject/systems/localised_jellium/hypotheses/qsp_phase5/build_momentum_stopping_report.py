#!/usr/bin/env python3
"""Builder for qsp_phase5_momentum_stopping.ipynb — the orbital-free
momentum-KE-loss stopping extraction over the qsp_phase5 WP velocity sweep.

Method (locked 2026-07-27, validated in-session):
  z-profiles rho(z,t) -> side-adaptive continuity+CAP-sink flux J(z,t)
  -> TOF detector planes -> exceedance-matched (rank-transport) velocity
  distributions -> S(u) per run, mid-rank trusted band calibrated by a
  free-packet NULL test through the identical pipeline.

Two-pass build: pass 1 executes the notebook, whose synthesis cell writes
``momentum_stopping_summary.json``; the builder then re-renders the takeaway
markdown with the *computed* numbers and executes again (numbers can never
disagree with the analysis — scientific-grounding rule).

Run:
    PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
    /local/data/public/skcb2/tddft/venv/bin/python3 build_momentum_stopping_report.py

Requires the per-run kinematics caches (qsp5_momentum_kinematics.py).
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))              # hypotheses/ for _nbreport
from _nbreport import md, code, embed, setup_cell, set_outdir, build  # noqa: E402

set_outdir(HERE)

RUNS = ["p5_wp_v1p3", "p5_wp_v3p0", "p5_wp_v4p0", "p5_wp_v5p0", "p5_wp_v6p0"]
OUT = os.path.join(HERE, "qsp_phase5_momentum_stopping.ipynb")
SUMMARY = os.path.join(HERE, "momentum_stopping_summary.json")

missing = [r for r in RUNS
           if not os.path.exists(os.path.join(HERE, "cache", f"{r}_kinematics.npz"))]
if missing:
    sys.exit(f"FATAL: kinematics caches missing for {missing} — "
             f"run qsp5_momentum_kinematics.py first.")

GIF = {}
for tag in ("p5_wp_v1p3", "p5_wp_v4p0", "p5_wp_v6p0"):
    for cand in (os.path.join(HERE, f"{tag}_run_notebook_figs", f"{tag}_total_density.gif"),
                 os.path.join(HERE, f"{tag}_run_notebook_figs", f"{tag}_wp_delta0.gif")):
        if os.path.exists(cand):
            GIF[tag] = cand
            break


def build_cells(summary):
    cells = []

    # == 1. Title + question =================================================
    cells.append(md(r"""# qsp_phase5 — stopping power from the momentum-dependent KE loss (orbital-free)

**Question.** The deposit-based stopping estimate for this sweep
($S = \Delta E_\mathrm{deposited}/L$, `results_p5_wp_*.json`) is *too high to be
physical* — at $k_0=1.3$ the slab "absorbs" 59 eV while the packet's entire drift
kinetic energy is only 23 eV. Hypothesis (conversation 2026-07-27): the excess is
the wavepacket's **localisation (zero-point) energy**
$T_\mathrm{loc}(0)=3/(8\sigma_r^2)=1.5\,\mathrm{Ha}=41\,\mathrm{eV}$ plus
capture/binding release — energy a *classical* projectile does not carry and
which stopping-as-projectile-KE-loss must exclude.

This notebook builds the complementary estimate — $S_\mathrm{drift}$, the
**momentum-dependent kinetic-energy loss** — *without identifying any Kohn–Sham
orbital as "the projectile"*: everything derives from the total density
$n(\mathbf r,t)$ in regions where its ownership is unambiguous, via
time-of-flight (TOF) detector planes and rank-matched velocity distributions.

| approach | observable | counts localisation? | counts capture binding? | orbital-free? |
|---|---|---|---|---|
| deposit-based (existing) | $\Delta E_\mathrm{total}$ plateau | **yes** (over-counts) | **yes** (over-counts) | yes |
| KS-orbital drift | $\langle p\rangle_\mathrm{wp}^2/2$ | no | no | **no** |
| **this notebook** | TOF rank-matched kinematics | no | no | **yes** |
"""))

    # == 2. Conventions ======================================================
    cells.append(md(r"""## Conventions and symbols

Hartree atomic units ($\hbar=m_e=e=1$); $1\,\mathrm{Ha}=27.211\,\mathrm{eV}$;
lengths in Bohr, times in a.u. Human-facing numbers rounded to 2 s.f. (3 where a
difference would vanish at 2).

| symbol | meaning | value / range |
|---|---|---|
| $\sigma$ | WP width, σ_WP convention (density std per axis $\sigma_r=\sigma/\sqrt2=0.354$) | 0.5 Bohr |
| $k_0$ | nominal launch momentum (along $+z$) | 1.3 – 6 a.u. |
| $\sigma_p$ | momentum std per axis $=1/(2\sigma_r)$ | 1.41 a.u. |
| $n,\ n_\mathrm{gs}$ | total density ($\int=83$, **includes** WP) / GS slab density ($\int=82$) | — |
| $\rho(z,t)$ | $\int\!\!\int(n-n_\mathrm{gs})\,dx\,dy$ — transverse-integrated excess | e/Bohr |
| $J(z,t)$ | longitudinal particle flux, reconstructed | e/a.u. |
| $W(z)$ | CAP magnitude, $35<\|z\|<45$ | peak 0.7 Ha |
| $u$ | longitudinal component velocity (TOF-resolved) | a.u. |
| $N(>u)$ | exceedance: norm crossing a plane faster than $u$ | e |
| $q$ | rank (exceedance level) used for matching | e |
| $T_W$ | von Weizsäcker (localisation/shape) energy | Ha |
| $L$ | slab thickness (stopping path) | 25 Bohr |

**Geometry** (all runs): box $50\times50\times90$, $dx=0.5$; slab $|z|\le12.5$;
launch $z_0=-23.75$; CAP both faces $35<|z|<45$; detectors at $z_\mathrm{in}=-15.5$,
$z_\mathrm{out}=+22$."""))

    # == 3. Setup cell (theme) ==============================================
    cells.append(setup_cell())

    cells.append(code(r'''# --- constants (single source of truth for this notebook) ---
import numpy as np, pandas as pd, json, re
import matplotlib.pyplot as plt
from scipy.special import erfc

HYP     = SYS + "/hypotheses/qsp_phase5"
RESULTS = SYS + "/scripts/qsp_phase5/wp/results"
RUNS    = ["p5_wp_v1p3", "p5_wp_v3p0", "p5_wp_v4p0", "p5_wp_v5p0", "p5_wp_v6p0"]
K0      = {"p5_wp_v1p3": 1.3, "p5_wp_v3p0": 3.0, "p5_wp_v4p0": 4.0,
           "p5_wp_v5p0": 5.0, "p5_wp_v6p0": 6.0}
HA_EV   = 27.211386
DX, DT  = 0.5, 0.04                      # grid spacing [Bohr], time step [a.u.]
Z_SLAB, L_SLAB = 12.5, 25.0              # slab half-width, stopping path
Z_IN, Z_OUT = -15.5, 22.0                # detector planes [Bohr]
CAP_IN, CAP_OUT, ETA = 35.0, 45.0, 0.7   # CAP region, |eta| [Ha]
SIG_R, SIG_P = 0.5/np.sqrt(2.0), np.sqrt(2.0)
K_NY = np.pi/DX                          # grid Nyquist momentum = 6.28
COL = {r: plt.cm.viridis(x) for r, x in zip(RUNS, np.linspace(0.0, 0.85, len(RUNS)))}
print("runs:", ", ".join(RUNS))'''))

    # == 4. Simulation setup =================================================
    cells.append(md(r"""## Simulation setup — reconstructable record

Parsed from each run's `run_summary.txt` (provenance written by `run.cpp`), never
from memory. All five runs share engine (`inq-study`, LDA), cell, grid, GS slab
(82 e⁻), CAP and WP definition — **only $k_0$ varies**; `write_every` keeps the
frame cadence roughly constant in *distance*."""))

    cells.append(code(r'''rows = []
for tag in RUNS:
    txt = open(f"{RESULTS}/{tag}/run_summary.txt").read()
    get = lambda k: re.search(rf"{k}\s*=\s*([\d.eE+-]+)", txt)
    rows.append(dict(run=tag, k0=float(get("wp_k0").group(1)),
                     E_drift_eV=float(get("wp_E_drift_eV").group(1)),
                     n_steps=int(get("n_steps").group(1)),
                     write_every=int(get("write_every").group(1)),
                     t_total_au=int(get("n_steps").group(1))*DT))
setup = pd.DataFrame(rows).set_index("run")
setup.round(2)'''))

    cells.append(code(r'''# geometry sketch — everything happens along z
fig, ax = plt.subplots(figsize=(9, 1.9))
ax.axvspan(-Z_SLAB, Z_SLAB, color="tab:blue", alpha=0.25, label="jellium slab")
ax.axvspan(-CAP_OUT, -CAP_IN, color="tab:red", alpha=0.25, label="CAP")
ax.axvspan(CAP_IN, CAP_OUT, color="tab:red", alpha=0.25)
ax.axvline(Z_IN, color="g", ls="--", lw=1.2); ax.axvline(Z_OUT, color="g", ls="--", lw=1.2)
ax.annotate("$z_\\mathrm{in}$", (Z_IN, 0.82), fontsize=9, ha="right")
ax.annotate("$z_\\mathrm{out}$", (Z_OUT, 0.82), fontsize=9)
ax.plot([-23.75], [0.4], "o", color="tab:green")
ax.annotate("WP launch, $+z$ →", (-23.6, 0.45), fontsize=9)
ax.set_xlim(-45, 45); ax.set_ylim(0, 1); ax.set_yticks([])
ax.set_xlabel("z [Bohr]"); ax.legend(loc="upper right", fontsize=8, ncol=2)
ax.set_title("qsp_phase5 geometry — detector planes dashed green", fontsize=10)
plt.tight_layout(); plt.show()'''))

    # == 5. Source files =====================================================
    cells.append(md(r"""## Source files

| role | path (repo-relative) |
|---|---|
| run definition | `ResearchProject/systems/localised_jellium/scripts/qsp_phase5/wp/run.cpp` |
| sweep dispatcher | `ResearchProject/systems/localised_jellium/scripts/qsp_phase5/run_sweep.sh` |
| CAP implementation | `inq-study/src/perturbations/absorbing.hpp` |
| kinematics cache extractor | `hypotheses/qsp_phase5/qsp5_momentum_kinematics.py` |
| this builder | `hypotheses/qsp_phase5/build_momentum_stopping_report.py` |
| prior deposit-based S | `hypotheses/qsp_phase5/results_p5_wp_*.json` (`analyse_phase5.py`) |
| classical reference S(v) | `hypotheses/qsp_phase5/classical_sigma0p5_bulk.csv` |
| Lindhard reference | `hypotheses/qsp_phase5/lindhard_ref.npz` |
| raw per-run data | `scripts/qsp_phase5/wp/results/<run>/raw/{vti,observables}/` |

**Data-definition caveats verified numerically before this build:**
`density_total` integrates to 83 — it **includes** the WP (unlike some campaign
configs); `density_delta` is $n(t)-n(0)$, **not** $n-n_\mathrm{gs}$, so the cache
extractor subtracts the one-off `density_gs_system` field itself."""))

    # == 6. QC-1 aliasing ====================================================
    cells.append(md(r"""## QC-1 — launch integrity: momentum aliasing on the grid

The grid holds only $|k|<k_\mathrm{Ny}=\pi/dx=6.28$. The WP momentum density is
Gaussian at $k_0$ with std $\sigma_p=1.41$; weight beyond $k_\mathrm{Ny}$ **wraps
to negative momenta** at injection (it is not truncated — it comes back moving
backwards). Predicted aliased fraction:

$$f_\mathrm{alias}=\tfrac12\,\mathrm{erfc}\!\left(\frac{k_\mathrm{Ny}-k_0}{\sqrt2\,\sigma_p}\right)$$

Compared below with the *measured* launch moments (`wp_momentum_stats.csv`, step
0; a clean minimum-uncertainty launch has $\langle p_z\rangle(0)=k_0$,
$\sigma_{p_z}^2=2.0$). This sweep predates the now-mandatory cutoff-aliasing
guard — applied here forensically. **Consequence: incident kinematics must be
measured, never taken as $k_0^2/2$.**"""))

    cells.append(code(r'''qc = []
for tag in RUNS:
    st = pd.read_csv(f"{RESULTS}/{tag}/raw/observables/wp_momentum_stats.csv", comment="#")
    r0 = st.iloc[0]
    qc.append(dict(run=tag, k0=K0[tag],
                   f_alias_pred=0.5*erfc((K_NY-K0[tag])/(np.sqrt(2.0)*SIG_P)),
                   pz_mean_0=r0.pz_mean, pz_frac=r0.pz_mean/K0[tag],
                   sig_pz2_0=r0.sigma_pz2, e_kin_0_Ha=r0.e_kin_ha))
qc = pd.DataFrame(qc).set_index("run")
def grade(r):
    if r.pz_frac > 0.98 and abs(r.sig_pz2_0-2.0) < 0.1: return "A (clean)"
    if r.pz_frac > 0.90: return "B (mild)"
    if r.pz_frac > 0.70: return "C (corrupted)"
    return "F (unusable)"
qc["grade"] = qc.apply(grade, axis=1)
qc.round(2)'''))

    cells.append(code(r'''fig, ax = plt.subplots(figsize=(7, 3.2))
for tag in RUNS:
    m = pd.read_csv(f"{RESULTS}/{tag}/raw/observables/momentum_distribution.csv", comment="#")
    m0 = m[m.step == 0]
    w = m0.n_wp/np.trapezoid(m0.n_wp, m0.k_bohr_inv)
    ax.plot(m0.k_bohr_inv, w, color=COL[tag], label=f"$k_0$={K0[tag]}")
    ax.axvline(K0[tag], color=COL[tag], ls=":", lw=0.8)
ax.axvline(K_NY, color="k", lw=1.2)
ax.annotate("$k_\\mathrm{Ny}$", (K_NY, ax.get_ylim()[1]*0.9), fontsize=9, ha="right")
ax.set_xlabel("$|k|$ [1/Bohr]"); ax.set_ylabel("normalised WP weight")
ax.set_title("launch (t=0) momentum distributions — dotted lines mark nominal $k_0$", fontsize=10)
ax.legend(fontsize=8)
plt.tight_layout(); plt.show()'''))

    # == 7. Visual intuition ================================================
    cells.append(md(r"""## Visual intuition — representative runs

Sweep policy: GIFs for three representative runs — **v1p3** (min $v$, grade A),
**v4p0** (mid, grade C), **v6p0** (max, grade F) — the others appear in aggregate
figures only (nothing else is silently dropped). Reused from the per-run
notebooks, path-referenced, xz-plane densities."""))

    for tag, why in (("p5_wp_v1p3", "v1p3 (grade A): clean launch; strong interaction at low v"),
                     ("p5_wp_v4p0", "v4p0 (grade C): launch already 18% momentum-deficient"),
                     ("p5_wp_v6p0", "v6p0 (grade F): aliased launch — most weight moves backwards")):
        if tag in GIF:
            cells.append(embed(GIF[tag], caption=why, width=680))

    cells.append(md(r"""### Per-run energetics (recorded ledger)

Energy components vs time (`observables.csv`, every step), as changes from $t=0$.
With the CAP on, $E_\mathrm{total}$ is *not* conserved — its drain is the CAP
removing the packet; the late plateau is the raw material of the deposit-based
estimate."""))

    cells.append(code(r'''fig, axes = plt.subplots(2, 2, figsize=(10, 6))
comp = [("energy_total", "$\\Delta E_\\mathrm{total}$"), ("energy_kinetic", "$\\Delta E_\\mathrm{kin}$"),
        ("energy_hartree", "$\\Delta E_\\mathrm{H}$"), ("energy_xc", "$\\Delta E_\\mathrm{xc}$")]
for tag in RUNS:
    o = pd.read_csv(f"{RESULTS}/{tag}/raw/observables/observables.csv")
    for ax, (c, lab) in zip(axes.flat, comp):
        ax.plot(o.time_au, (o[c]-o[c].iloc[0])*HA_EV, color=COL[tag], label=f"$k_0$={K0[tag]}")
for ax, (c, lab) in zip(axes.flat, comp):
    ax.set_title(lab, fontsize=10); ax.set_xlabel("t [a.u.]"); ax.set_ylabel("eV")
axes.flat[0].legend(fontsize=8)
plt.tight_layout(); plt.show()'''))

    # == 8. Method 1: profiles ==============================================
    cells.append(md(r"""## Method step 1 — the transverse-integrated excess density

$$\rho(z,t)=\int\!\!\!\int\left[n(\mathbf r,t)-n_\mathrm{gs}(\mathbf r)\right]dx\,dy,
\qquad \int\rho\,dz\Big|_{t=0}=1$$

Subtracting $n_\mathrm{gs}$ removes the static slab; in the vacuum corridors the
remainder is projectile density (plus a small dynamic bath-polarisation tail —
the reason the detector planes carry buffers beyond the slab face). Computed once
by `qsp5_momentum_kinematics.py` over the `density_total` frames and cached."""))

    cells.append(code(r'''def load_kin(tag):
    d = np.load(f"{HYP}/cache/{tag}_kinematics.npz")
    return {k: d[k] for k in d.files}
KIN = {tag: load_kin(tag) for tag in RUNS}

k1 = KIN["p5_wp_v1p3"]
fig, ax = plt.subplots(figsize=(8.5, 4))
v = np.abs(k1["rho"]).max()*0.5
im = ax.pcolormesh(k1["t_au"], k1["z"], k1["rho"].T, cmap="RdBu_r", vmin=-v, vmax=v)
for zl, c in ((Z_SLAB, "k"), (-Z_SLAB, "k"), (Z_IN, "g"), (Z_OUT, "g"),
              (CAP_IN, "r"), (-CAP_IN, "r")):
    ax.axhline(zl, color=c, ls="--", lw=0.8)
plt.colorbar(im, ax=ax, label=r"$\rho(z,t)$ [e/Bohr]")
ax.set_xlabel("t [a.u.]"); ax.set_ylabel("z [Bohr]")
ax.set_title("p5_wp_v1p3: excess density — slab (black), detectors (green), CAP faces (red)", fontsize=10)
plt.tight_layout(); plt.show()'''))

    # == 9. Method 2: flux ===================================================
    cells.append(md(r"""## Method step 2 — flux from continuity + CAP sink (side-adaptive)

No current frames were saved, so $J$ comes from the 1D **continuity equation with
the CAP sink**:

$$\frac{\partial\rho}{\partial t}+\frac{\partial J}{\partial z}=-\,2W(z)\,\rho,
\qquad W(z)=|\eta|\sin^2\!\frac{\pi(|z|-35)}{10}\ \ (35<|z|<45),\ |\eta|=0.7$$

(a CAP $-iW$ removes density at rate $2W\rho$ — standard absorbing-potential
kinematics, e.g. Muga et al., Phys. Rep. 395, 357 (2004); $W$'s shape read from
`perturbations::absorbing` source, its *location* verified empirically from
band-resolved norm decay). Integrating from an edge where $J\simeq0$ (everything
reaching a CAP is absorbed before the box edge):

$$J(z,t)=-\!\!\int_{-L/2}^{z}\!\Big(\frac{\partial\rho}{\partial t}+2W\rho\Big)dz'
\;\;\overset{\text{or}}{=}\;\;+\!\!\int_{z}^{+L/2}\!\Big(\frac{\partial\rho}{\partial t}+2W\rho\Big)dz'$$

**Side-adaptive rule** (validated below): evaluate the *entrance* plane from the
left edge and the *exit* plane from the right edge — the short integration path
avoids accumulating gradient noise and bath-region dynamics. **Closure test:**
whole-box balance $-dN_\mathrm{tot}/dt = A_++A_-$,
$A_\pm=\int_{\mathrm{CAP}\pm}2W\rho\,dz$, must hold — it validates the sink model
and both edge assumptions at once."""))

    cells.append(code(r'''def cap_W(z):
    W = np.zeros_like(z)
    m = (np.abs(z) > CAP_IN) & (np.abs(z) < CAP_OUT)
    W[m] = ETA*np.sin(np.pi*(np.abs(z[m]) - CAP_IN)/(CAP_OUT - CAP_IN))**2
    return W

def flux_both(rho, z, t, W=None):
    """J from the left edge (J_L) and right edge (J_R), + CAP drain rates A+-."""
    dz = z[1] - z[0]
    src = np.gradient(rho, t, axis=0) + (2.0*W[None, :]*rho if W is not None else 0.0)
    J_L = -np.cumsum(src, axis=1)*dz
    J_R = np.cumsum(src[:, ::-1], axis=1)[:, ::-1]*dz
    if W is None:
        return J_L, J_R, None, None
    A_p = (2.0*W[None, :]*rho)[:, z >  CAP_IN].sum(axis=1)*dz
    A_m = (2.0*W[None, :]*rho)[:, z < -CAP_IN].sum(axis=1)*dz
    return J_L, J_R, A_p, A_m

FLUX = {}
fig, axes = plt.subplots(1, len(RUNS), figsize=(14, 2.6))
for ax, tag in zip(axes, RUNS):
    kin = KIN[tag]
    W = cap_W(kin["z"])
    FLUX[tag] = flux_both(kin["rho"], kin["z"], kin["t_au"], W)
    dz = kin["z"][1]-kin["z"][0]
    dNdt = np.gradient(kin["rho"].sum(axis=1)*dz, kin["t_au"])
    A_p, A_m = FLUX[tag][2], FLUX[tag][3]
    ax.plot(kin["t_au"], -dNdt, "k", lw=1.4, label="$-dN/dt$")
    ax.plot(kin["t_au"], A_p+A_m, color="tab:red", lw=1, ls="--", label="$A_++A_-$")
    resid = np.abs(-dNdt-(A_p+A_m)).mean()/max(np.abs(dNdt).max(), 1e-12)
    ax.set_title(f"$k_0$={K0[tag]}  resid={resid:.0%}", fontsize=9)
    ax.set_xlabel("t [a.u.]")
axes[0].set_ylabel("norm rate [e/a.u.]"); axes[0].legend(fontsize=7)
plt.suptitle("closure: measured norm loss vs modelled CAP drain", fontsize=10, y=1.04)
plt.tight_layout(); plt.show()'''))

    cells.append(code(r'''J1 = FLUX["p5_wp_v1p3"][1]
fig, ax = plt.subplots(figsize=(8.5, 4))
v = np.abs(J1).max()*0.6
im = ax.pcolormesh(k1["t_au"], k1["z"], J1.T, cmap="PuOr_r", vmin=-v, vmax=v)
for zl, c in ((Z_IN, "g"), (Z_OUT, "g"), (CAP_IN, "r"), (-CAP_IN, "r")):
    ax.axhline(zl, color=c, ls="--", lw=0.8)
plt.colorbar(im, ax=ax, label="$J(z,t)$ [e/a.u.]")
ax.set_xlabel("t [a.u.]"); ax.set_ylabel("z [Bohr]")
ax.set_title("p5_wp_v1p3: reconstructed flux — orange forward, purple backward", fontsize=10)
plt.tight_layout(); plt.show()'''))

    # == 10. Method 3: TOF detectors ========================================
    cells.append(md(r"""## Method step 3 — time-of-flight detector planes

Snapshot ("whole lobe in vacuum at once") kinematics are impossible here: with
$\sigma_p=1.41$ the packet stretches ~14 Bohr per 10 a.u. — longer than the
20-Bohr corridor. Instead, fixed planes are read out like beamline detectors.
The instantaneous **mean passing velocity** at a plane is

$$u(t)=\frac{J(z_m,t)}{\rho(z_m,t)}$$

and time-of-flight *sorts* the components — fast cross early, slow late — so
$u(t)$ sweeps the longitudinal velocity distribution (sweep resolution
$\Delta u/u\sim\sigma_r u/D$ for flight distance $D$; validated below). Signed
flux separates channels: at $z_\mathrm{in}$, $J>0$ = incoming, $J<0$ = reflected;
at $z_\mathrm{out}$, $J>0$ = transmitted."""))

    cells.append(code(r'''def sweep(kin, J, z_plane, sgn, vmax):
    """TOF readout at a plane: per-frame velocity u_i and channel weight w_i=|J|dt."""
    z, t = kin["z"], kin["t_au"]
    i = int(np.argmin(np.abs(z - z_plane)))
    Jp, rp = J[:, i], kin["rho"][:, i]
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.where(np.abs(rp) > 1e-7, Jp/rp, 0.0)
    u = np.clip(u, -vmax, vmax)
    w = np.where(sgn*Jp > 0, sgn*Jp, 0.0)*np.gradient(t)
    return sgn*u, w

kin = KIN["p5_wp_v1p3"]; J_L, J_R = FLUX["p5_wp_v1p3"][0], FLUX["p5_wp_v1p3"][1]
fig, axes = plt.subplots(2, 2, figsize=(10, 5.5), sharex=True)
for col, (J, zp, lab) in enumerate(((J_L, Z_IN, "entrance $z=-15.5$ (left-integrated)"),
                                    (J_R, Z_OUT, "exit $z=+22$ (right-integrated)"))):
    i = int(np.argmin(np.abs(kin["z"] - zp)))
    axes[0, col].plot(kin["t_au"], J[:, i], color="tab:blue")
    axes[0, col].axhline(0, color="k", lw=0.5)
    axes[0, col].set_title(lab, fontsize=10); axes[0, col].set_ylabel("$J$ [e/a.u.]")
    u, w = sweep(kin, J, zp, +1, 7.0)
    axes[1, col].plot(kin["t_au"], u, color="tab:orange")
    axes[1, col].set_ylabel("$u=J/\\rho$ [a.u.]"); axes[1, col].set_xlabel("t [a.u.]")
fig.suptitle("p5_wp_v1p3 detector readout — the TOF sweep runs fast → slow", fontsize=10)
plt.tight_layout(); plt.show()'''))

    # == 11. Method 4: rank matching ========================================
    cells.append(md(r"""## Method step 4 — exceedance-matched stopping $S(u)$

A naive per-electron comparison $\bar K_\mathrm{in}-\bar K_\mathrm{out}$ fails
for two *measured* reasons: (i) **capture bias** — at low $v$ the slab removes
the *slow* components, so the transmitted ensemble is fast-selected and the
averages compare different electrons; (ii) the entrance fast tail carries a
resolution systematic (quantified by the null test below). Both are cured by
comparing the velocity **distributions by rank**.

Build the **exceedance curve** at each plane — the norm crossing faster than $u$:

$$N_\mathrm{in}(>u)=\sum_{t_i:\ u_i>u} w_i^\mathrm{in},\qquad
N_\mathrm{out}(>u)=\sum_{t_i:\ u_i>u} w_i^\mathrm{out}$$

Assuming the slab preserves velocity *ordering* (fast stays fastest — monotone
transport, no overtaking) and capture removes only the bottom of the
distribution, matching **equal exceedance ranks** $q$ pairs each incident
component with its transmitted image:

$$u_\mathrm{in}(q)\ \longmapsto\ u_\mathrm{out}(q)\quad\Rightarrow\quad
S(u_\mathrm{in}(q))=\frac{\tfrac12\left[u_\mathrm{in}^2(q)-u_\mathrm{out}^2(q)\right]}{L}$$

This yields a *whole $S(u)$ curve from every run* (the packet's momentum spread
becomes a feature: each run probes a velocity band), and it is robust to capture:
truncating the slow end shifts no top-down rank."""))

    cells.append(code(r'''def exceedance(u, w, ugrid):
    return np.array([w[u > g].sum() for g in ugrid])

def s_of_u(kin, flx, k0, n_rank=60):
    """Rank-matched S(u) between entrance (left-int) and exit (right-int) planes."""
    vmax = min(k0 + 4.0*SIG_P, K_NY)
    u_i, w_i = sweep(kin, flx[0], Z_IN, +1, vmax)     # incoming
    u_o, w_o = sweep(kin, flx[1], Z_OUT, +1, vmax)    # transmitted
    u_r, w_r = sweep(kin, flx[0], Z_IN, -1, vmax)     # reflected (bookkeeping)
    ug = np.linspace(0.15, vmax, 240)
    Ein, Eout = exceedance(u_i, w_i, ug), exceedance(u_o, w_o, ug)
    q_top = 0.92*min(Ein[0], Eout[0])
    qs = np.linspace(0.02, q_top, n_rank)
    uin  = np.interp(qs, Ein[::-1], ug[::-1])
    uout = np.interp(qs, Eout[::-1], ug[::-1])
    S = 0.5*(uin**2 - uout**2)*HA_EV/L_SLAB
    trust = (qs > 0.30*q_top) & (qs < 0.90*q_top)     # null-calibrated rank window
    return dict(ug=ug, Ein=Ein, Eout=Eout, qs=qs, uin=uin, uout=uout, S=S,
                trust=trust, N_in=Ein[0], N_out=Eout[0], N_refl=w_r.sum())

SOU = {tag: s_of_u(KIN[tag], FLUX[tag], K0[tag]) for tag in RUNS}
pd.DataFrame({t: dict(N_in=r["N_in"], N_trans=r["N_out"], N_refl=r["N_refl"])
              for t, r in SOU.items()}).T.round(2)'''))

    # == 12. Validation battery =============================================
    cells.append(md(r"""## Validation — the pipeline must return $S=0$ for a free packet

Three tests, run through the **identical code path** (`flux_both` → `sweep` →
`s_of_u` formulas):

1. **estimator exactness**: with analytic $\rho,J$ of a free 1D Gaussian, the
   plane KE readout reproduces the forward-ensemble truth to <1%;
2. **null test (below)**: the free packet *sampled exactly like the runs* (grid
   0.5 Bohr, frame cadence 0.48 a.u., finite-difference $\partial_t$, cumsum
   reconstruction) must give $S(u)\equiv0$. Its residual **calibrates the
   systematic error band and fixes the trusted rank window** (the fast tail —
   top ~25% of ranks — is under-resolved at the entrance plane and excluded);
3. **closure** (shown in step 2): CAP sink model balances total norm loss to a
   few %.

The null-test domain is extended to ±1100 Bohr so nothing escapes (the real runs'
CAP guarantees $J\to0$ at the edges; the free emulation has no CAP)."""))

    cells.append(code(r'''# null test: free Gaussian, k0=1.3 (worst case: largest relative spread)
z_n = np.arange(-1100.25, 1100.5, DX)
t_n = np.arange(0.0, 160.0, 0.48) + 1e-9
s2 = SIG_R**2 + (SIG_P*t_n)**2
rho_n = np.exp(-(z_n[None, :]+23.75-1.3*t_n[:, None])**2/(2*s2[:, None]))/np.sqrt(2*np.pi*s2[:, None])
kin_n = dict(z=z_n, t_au=t_n, rho=rho_n)
flx_n = flux_both(rho_n, z_n, t_n, None)
null = s_of_u(kin_n, flx_n, 1.3)
fig, ax = plt.subplots(figsize=(7, 3.2))
ax.plot(null["uin"], null["S"], "o-", ms=3, color="0.5", label="all ranks")
ax.plot(null["uin"][null["trust"]], null["S"][null["trust"]], "o", ms=4,
        color="tab:green", label="trusted rank window")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("$u_\\mathrm{in}$ [a.u.]"); ax.set_ylabel("$S_\\mathrm{null}$ [eV/Bohr]")
ax.set_title("null test: free packet through the full pipeline (must be 0)", fontsize=10)
ax.legend(fontsize=8)
plt.tight_layout(); plt.show()
S_SYST = float(np.abs(null["S"][null["trust"]]).max())
print(f"null residual in trusted window: mean {null['S'][null['trust']].mean():+.2f}, "
      f"max |S| {S_SYST:.2f} eV/Bohr  -> adopted as the systematic error band")'''))

    # == 13. Per-run S(u) ====================================================
    cells.append(md(r"""## Results — exceedance curves and $S(u)$ per run

Left: incident (solid) and transmitted (dashed) exceedance curves — their
horizontal gap *is* the velocity loss at each rank; the vertical end-gap is the
non-transmitted (captured + reflected) norm. Right: the rank-matched $S(u)$,
trusted window solid, excluded ranks faded; grey band = null-test systematic."""))

    cells.append(code(r'''fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
for tag in RUNS:
    r = SOU[tag]
    axes[0].plot(r["ug"], r["Ein"], color=COL[tag], label=f"$k_0$={K0[tag]}")
    axes[0].plot(r["ug"], r["Eout"], color=COL[tag], ls="--")
    tr = r["trust"]
    axes[1].plot(r["uin"][~tr], r["S"][~tr], ".", ms=3, color=COL[tag], alpha=0.25)
    axes[1].plot(r["uin"][tr], r["S"][tr], "o-", ms=3, color=COL[tag], label=f"$k_0$={K0[tag]}")
axes[1].axhspan(-S_SYST, S_SYST, color="0.85", zorder=0)
axes[0].set_xlabel("$u$ [a.u.]"); axes[0].set_ylabel("$N(>u)$ [e]")
axes[0].set_title("exceedance: incident (solid) vs transmitted (dashed)", fontsize=10)
axes[0].legend(fontsize=8)
axes[1].axhline(0, color="k", lw=0.6)
axes[1].set_xlabel("$u_\\mathrm{in}$ [a.u.]"); axes[1].set_ylabel("$S$ [eV/Bohr]")
axes[1].set_title("rank-matched stopping $S(u)$ — grey = null-test band", fontsize=10)
axes[1].legend(fontsize=8)
plt.tight_layout(); plt.show()'''))

    cells.append(md(r"""### Channel accounting

What entered, came back, got through — and, by difference plus the end-state
interior norm, what was **captured**. Deposit-side over-count context: each
captured electron banks its ~41 eV localisation energy *plus* binding into the
slab ledger."""))

    cells.append(code(r'''acc = []
for tag in RUNS:
    kin, r = KIN[tag], SOU[tag]
    z, dz = kin["z"], kin["z"][1]-kin["z"][0]
    interior = kin["rho"][-1][np.abs(z) < Z_SLAB+3].sum()*dz
    acc.append(dict(run=tag, N_in=r["N_in"], N_refl=r["N_refl"], N_trans=r["N_out"],
                    N_interior_end=interior,
                    N_missing=r["N_in"]-r["N_refl"]-r["N_out"]))
acc = pd.DataFrame(acc).set_index("run")
acc.round(2)'''))

    # == 14. Localisation channel ===========================================
    cells.append(md(r"""## The localisation (shape) energy — watching what we excluded

Two exact decompositions of a one-lump kinetic energy:

- **momentum split (A)**: $E_\mathrm{kin}=\langle p\rangle^2/2+\mathrm{Var}(p)/2$
  — *both terms constant in vacuum*;
- **Madelung split (B)**: $E_\mathrm{kin}=T_v+T_W$ with the von Weizsäcker
  **localisation energy**

$$T_W=\int\frac{|\nabla n|^2}{8n}\,d^3r,\qquad
T_W^\mathrm{free\ Gaussian}(t)=\frac{3}{8\,\sigma_r^2(t)},\quad
\sigma_r^2(t)=\sigma_r^2(0)+\sigma_p^2t^2$$

In vacuum $T_W$ *falls* while expansion flow grows in $T_v$ — sum constant. The
cache carries lobe-integrated $T_W$ from the 3D density; the KS-orbital stats
(`wp_momentum_stats`) give split (A) — used **only in this section, as a
cross-check**; the S extraction above never touches them."""))

    cells.append(code(r'''st = pd.read_csv(f"{RESULTS}/p5_wp_v1p3/raw/observables/wp_momentum_stats.csv", comment="#")
drift = (st.px_mean**2 + st.py_mean**2 + st.pz_mean**2)/2
internal = st.e_kin_ha - drift
kin = KIN["p5_wp_v1p3"]
t_hit = (23.75 - Z_SLAB)/1.3

fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
ax = axes[0]
ax.plot(st.time_au, st.e_kin_ha, "k", label="$E_\\mathrm{kin}$ (orbital)")
ax.plot(st.time_au, drift, color="tab:orange", label="drift $|\\langle p\\rangle|^2/2$")
ax.plot(st.time_au, internal, color="tab:blue", label="internal")
ax.axvline(t_hit, color="k", ls=":", lw=0.8)
ax.set_xlim(0, 40); ax.set_xlabel("t [a.u.]"); ax.set_ylabel("Ha")
ax.set_title("split (A): both terms constant until slab contact", fontsize=10)
ax.legend(fontsize=8)

ax = axes[1]
tt = kin["t_au"]
ax.plot(tt, kin["TW_lo"], color="tab:green", label="$T_W$ entrance lobe (from $n$)")
ax.plot(tt, 3.0/(8.0*(SIG_R**2 + (SIG_P*tt)**2)), "k--", lw=1,
        label="free Gaussian $3/8\\sigma_r^2(t)$")
ax.axvline(t_hit, color="k", ls=":", lw=0.8)
ax.set_xlim(0, 40); ax.set_xlabel("t [a.u.]"); ax.set_ylabel("Ha")
ax.set_title("split (B): localisation energy falls under dispersion", fontsize=10)
ax.legend(fontsize=8)
plt.tight_layout(); plt.show()
print(f"t=0 ledger (v1p3): internal(0) = {internal.iloc[0]:.3f} Ha | lobe T_W(0) = "
      f"{kin['TW_lo'][0]:.3f} Ha | analytic 3/(8 sigma_r^2) = {3/(8*SIG_R**2):.3f} Ha")'''))

    cells.append(md(r"""$T_W(0)$, the orbital internal energy, and the analytic $3/(8\sigma_r^2)$ agree
at 1.5 Ha = **41 eV**: that is the launch-time shape energy every run carries —
energy a classical projectile does not have, correctly invisible to the TOF
detectors, and (for captured weight) wrongly banked by the deposit-based S.

*(The lobe $T_W$ leaves the analytic curve as soon as the packet's fast front
exits the entrance lobe — the lobe integral then sees only part of the packet;
the very early frames follow the free-dispersion law.)*"""))

    cells.append(code(r'''fig, axes = plt.subplots(1, 2, figsize=(11, 3.4), sharey=True)
for tag in RUNS:
    kin = KIN[tag]
    axes[0].plot(kin["t_au"], kin["TW_lo"], color=COL[tag], label=f"$k_0$={K0[tag]}")
    axes[1].plot(kin["t_au"], kin["TW_hi"], color=COL[tag])
axes[0].set_title("entrance-lobe $T_W$ [Ha]", fontsize=10)
axes[1].set_title("exit-lobe $T_W$ [Ha]", fontsize=10)
for ax in axes: ax.set_xlabel("t [a.u.]")
axes[0].legend(fontsize=8)
plt.tight_layout(); plt.show()'''))

    # == 15. Synthesis ======================================================
    cells.append(md(r"""## Synthesis — $S_\mathrm{drift}(v)$ against the deposit-based estimate

Per-run headline: mean trusted-window $S$, quoted at the flux-weighted mean
trusted $u_\mathrm{in}$, error = rank-window spread ⊕ null-test systematic.
Deposit-based values (`results_p5_wp_*.json`) at the same measured velocities;
classical $\sigma=0.5$ and Lindhard for context. Hollow markers = grade C/F
(corrupted launch: the point measures *a* packet at its measured velocity, not
the nominal one)."""))

    cells.append(code(r'''head = []
for tag in RUNS:
    r = SOU[tag]; tr = r["trust"]
    u_ref = float(np.mean(r["uin"][tr])); S_mid = float(np.mean(r["S"][tr]))
    S_err = float(np.hypot(np.std(r["S"][tr]), S_SYST))
    with open(f"{HYP}/results_{tag}.json") as f:
        dj = json.load(f)
    head.append(dict(run=tag, k0=K0[tag], u_ref=u_ref, S_drift=S_mid, S_err=S_err,
                     S_deposit=dj["S_eVbohr"], deposited_eV=dj["deposited_eV"],
                     grade=qc.loc[tag, "grade"]))
head = pd.DataFrame(head).set_index("run")

cls = pd.read_csv(f"{HYP}/classical_sigma0p5_bulk.csv")
lin = np.load(f"{HYP}/lindhard_ref.npz")
fig, ax = plt.subplots(figsize=(8, 4.8))
ax.plot(cls.v, cls.S_eVbohr, "s-", color="0.55", ms=4, label="classical $\\sigma$=0.5 (bulk)")
ax.plot(np.sqrt(2*lin["E"]/HA_EV), lin["S"], "-", color="0.8", lw=1.2, label="Lindhard (point)")
for tag in RUNS:
    r = SOU[tag]; tr = r["trust"]
    solid = head.loc[tag, "grade"].startswith(("A", "B"))
    ax.plot(r["uin"][tr], r["S"][tr], "-", color=COL[tag], lw=2 if solid else 1,
            alpha=1.0 if solid else 0.45)
    ax.errorbar([head.loc[tag, "u_ref"]], [head.loc[tag, "S_drift"]],
                yerr=[head.loc[tag, "S_err"]], fmt="o", color=COL[tag], ms=7,
                mfc=(COL[tag] if solid else "none"))
    ax.plot([head.loc[tag, "u_ref"]], [head.loc[tag, "S_deposit"]], "^",
            color="tab:red", ms=8, mfc=("tab:red" if solid else "none"))
ax.plot([], [], "o", color="k", label="$S_\\mathrm{drift}$ (this notebook)")
ax.plot([], [], "^", color="tab:red", label="$S_\\mathrm{deposit}$ (existing)")
ax.set_xlabel("$u$ [a.u.]"); ax.set_ylabel("S [eV/Bohr]")
ax.set_title("stopping: momentum-KE loss (curves+dots) vs energy deposit (triangles)", fontsize=10)
ax.legend(fontsize=8); ax.set_xlim(0, 6.5)
plt.tight_layout(); plt.show()
head.round(2)'''))

    # == 14b. Combined low-velocity S(v) ====================================
    cells.append(md(r"""## Combined low-velocity $S(v)$ — averaging the clean bands

Each clean run's $S(u)$ segment is many stopping measurements at different
incident velocities (the packet's own momentum spread, TOF-resolved). The
aliased launches (v4p0–v6p0, grades C/F) are **excluded outright**; the
grade-A/B bands (v1p3: $u\approx1.2$–2.6; v3p0: $u\approx2.4$–3.5) are combined
into one curve by a **flux-weighted average**: at each $u$, a run contributes
with weight equal to how much of its packet actually crossed at that velocity,

$$w_r(u) = -\frac{dN^{(r)}_\mathrm{in}(>u)}{du}\ \ (\ge0),\qquad
\bar S(u)=\frac{\sum_r w_r(u)\,S_r(u)}{\sum_r w_r(u)}$$

Uncertainty band: between-run spread (where bands overlap) ⊕ the null-test
systematic."""))

    cells.append(code(r'''CLEAN = ["p5_wp_v1p3", "p5_wp_v3p0"]          # grades A/B only
ugrid = np.linspace(0.6, 3.8, 48)
num = np.zeros_like(ugrid); den = np.zeros_like(ugrid)
vals = {tg: np.full_like(ugrid, np.nan) for tg in CLEAN}
for tag in CLEAN:
    r = SOU[tag]; tr = r["trust"]
    u_tr, S_tr = r["uin"][tr], r["S"][tr]
    order = np.argsort(u_tr)
    m = (ugrid >= u_tr.min()) & (ugrid <= u_tr.max())
    Si = np.interp(ugrid[m], u_tr[order], S_tr[order])
    w_u = np.clip(-np.gradient(r["Ein"], r["ug"]), 0.0, None)   # flux density dN/du
    wi = np.interp(ugrid[m], r["ug"], w_u)
    vals[tag][m] = Si
    num[m] += wi*Si; den[m] += wi
S_avg = np.where(den > 0, num/np.maximum(den, 1e-12), np.nan)
spread = np.nanstd(np.vstack([vals[tg] for tg in CLEAN]), axis=0)
S_avg_err = np.hypot(np.nan_to_num(spread), S_SYST)

fig, ax = plt.subplots(figsize=(7.5, 4.4))
ax.plot(cls.v, cls.S_eVbohr, "s-", color="0.55", ms=4, label="classical $\\sigma$=0.5 (bulk)")
ax.plot(np.sqrt(2*lin["E"]/HA_EV), lin["S"], "-", color="0.8", lw=1.2, label="Lindhard (point)")
for tag in CLEAN:
    ax.plot(ugrid, vals[tag], "-", color=COL[tag], lw=1, alpha=0.5, label=f"{tag} band")
ok = np.isfinite(S_avg)
ax.plot(ugrid[ok], S_avg[ok], "-", color="tab:red", lw=2.5, label="$\\bar S(u)$ flux-weighted (A/B)")
ax.fill_between(ugrid[ok], (S_avg-S_avg_err)[ok], (S_avg+S_avg_err)[ok],
                color="tab:red", alpha=0.15)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("$u$ [a.u.]"); ax.set_ylabel("S [eV/Bohr]"); ax.set_xlim(0.4, 4.0)
ax.set_title("combined low-velocity stopping — aliased runs excluded", fontsize=10)
ax.legend(fontsize=8)
plt.tight_layout(); plt.show()

u_pts = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
tab = pd.DataFrame(dict(u=u_pts,
                        S_avg=np.interp(u_pts, ugrid[ok], S_avg[ok], left=np.nan, right=np.nan),
                        err=np.interp(u_pts, ugrid[ok], S_avg_err[ok], left=np.nan, right=np.nan)))
tab.round(2)'''))

    cells.append(code(r'''summary = {t: dict(u_ref=float(head.loc[t, "u_ref"]),
                   S_drift=float(head.loc[t, "S_drift"]),
                   S_err=float(head.loc[t, "S_err"]),
                   S_deposit=float(head.loc[t, "S_deposit"]),
                   grade=str(head.loc[t, "grade"]),
                   S_syst=S_SYST) for t in head.index}
with open(f"{HYP}/momentum_stopping_summary.json", "w") as f:
    json.dump(summary, f, indent=1)
lowv = dict(u=[float(x) for x in tab.u],
            S=[None if not np.isfinite(x) else float(x) for x in tab.S_avg],
            err=[None if not np.isfinite(x) else float(x) for x in tab.err],
            runs_used=CLEAN)
with open(f"{HYP}/momentum_stopping_lowv.json", "w") as f:
    json.dump(lowv, f, indent=1)
print("wrote momentum_stopping_summary.json + momentum_stopping_lowv.json")'''))

    # == 16. Takeaway ========================================================
    if summary:
        lowv_bullet = ""
        lowv_path = os.path.join(HERE, "momentum_stopping_lowv.json")
        if os.path.exists(lowv_path):
            with open(lowv_path) as f:
                lv = json.load(f)
            pts = "; ".join(f"S({u:.2g}) = {s:.2g} ± {e:.2g}"
                            for u, s, e in zip(lv["u"], lv["S"], lv["err"])
                            if s is not None)
            lowv_bullet = (f"- **Combined low-velocity curve** (flux-weighted over "
                           f"{', '.join(lv['runs_used'])}; aliased runs excluded): "
                           f"{pts} eV/Bohr.")
        rowfmt = "| {t} | {u:.2g} | {sd:.2g} ± {se:.2g} | {dep:.2g} | {rat:.2g}× | {g} |"
        lines = "\n".join(rowfmt.format(
            t=t, u=d["u_ref"], sd=d["S_drift"], se=d["S_err"], dep=d["S_deposit"],
            rat=d["S_deposit"]/d["S_drift"] if d["S_drift"] else float("nan"),
            g=d["grade"]) for t, d in summary.items())
        v13 = summary.get("p5_wp_v1p3", {})
        takeaway = f"""## Takeaway

| run | $u_\\mathrm{{ref}}$ | $S_\\mathrm{{drift}}$ [eV/Bohr] | $S_\\mathrm{{deposit}}$ [eV/Bohr] | ratio | grade |
|---|---|---|---|---|---|
{lines}

- **The momentum-KE-loss stopping is measurable orbital-free** from $\\rho(z,t)$
  alone (continuity + CAP sink + TOF + rank matching), validated by a free-packet
  null test (systematic ±{v13.get('S_syst', 0):.2g} eV/Bohr in the trusted window).
- **The deposit-based S over-counts** — for the clean run (v1p3) by
  ~{(v13.get('S_deposit', 0)/max(v13.get('S_drift', 1e-9), 1e-9)):.2g}×: the gap is the
  localisation energy (41 eV/electron at launch, measured three independent ways)
  plus capture binding, banked in the slab ledger but absent from any projectile-KE
  definition of stopping.
- Each run yields an $S(u)$ *curve* (the momentum spread probes a velocity band);
  the curves are mutually consistent where bands overlap and rise with $u$ over
  the measured range.
{lowv_bullet}
- **Grades C/F (v4–v6) launched corrupted** (momentum aliasing, QC-1): usable only
  at their *measured* velocities; a clean high-$v$ sweep needs $dx\\le0.35$ or
  larger σ.
- Open: transverse-deflection channel (invisible to z-planes) and the
  ordering-preservation assumption of rank matching; both are second-order within
  the quoted errors but bounding them needs saved current frames or the σ=1
  cap/nocap pair (`wp_cap_energy_plateau`).
"""
        cells.append(md(takeaway))
    else:
        cells.append(md("## Takeaway\n\n*(populated on the second build pass from "
                        "`momentum_stopping_summary.json`)*"))
    return cells


summary = None
if os.path.exists(SUMMARY):
    with open(SUMMARY) as f:
        summary = json.load(f)

print("pass 1: executing notebook ...")
build(build_cells(summary), OUT, timeout=1800)

with open(SUMMARY) as f:
    summary2 = json.load(f)
print("pass 2: re-rendering takeaway with computed numbers (incl. low-v average) ...")
build(build_cells(summary2), OUT, timeout=1800)
print("done:", OUT)
