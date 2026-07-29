#!/usr/bin/env python3
"""Builder for wp_cap_energy_plateau_momentum_stopping.ipynb — plateau
dissection + orbital-free momentum-KE-loss stopping on the sigma=1, E=100 eV
cap/nocap pair (the clean-launch plateau runs).

Reuses the TOF/rank-matching method validated in
hypotheses/qsp_phase5/build_momentum_stopping_report.py (same inline pipeline,
this pair's geometry), and adds what qsp_phase5 could not do: a term-by-term
energy-ledger dissection of the CAP plateau (pairwise interactions.csv, in eV),
including a direct test of the naive estimator
"S from (E(0) - E_plateau) minus the initial localisation energy".

Two-pass build (summary json -> takeaway numbers). Run:
    PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
    /local/data/public/skcb2/tddft/venv/bin/python3 build_plateau_momentum_report.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from _nbreport import md, code, embed, setup_cell, set_outdir, build  # noqa: E402

set_outdir(HERE)

OUT = os.path.join(HERE, "wp_cap_energy_plateau_momentum_stopping.ipynb")
SUMMARY = os.path.join(HERE, "plateau_dissection_summary.json")
RESULTS = os.path.normpath(os.path.join(
    HERE, "..", "..", "scripts", "wp_cap_energy_plateau", "wp", "results"))

for r in ("cap", "nocap"):
    if not os.path.exists(os.path.join(HERE, "cache", f"{r}_kinematics.npz")):
        sys.exit(f"FATAL: cache missing for '{r}' — run plateau_kinematics.py first.")

GIF = {}
for r in ("cap", "nocap"):
    p = os.path.join(RESULTS, r, "report", "wp_total_density.gif")
    if os.path.exists(p):
        GIF[r] = p


def build_cells(summary):
    cells = []

    # == 1. Title + question =================================================
    cells.append(md(r"""# wp_cap_energy_plateau — dissecting the energy plateau; stopping from momentum-KE loss (σ = 1)

**Question.** In CAP runs the total energy *plateaus* once the wavepacket has
been absorbed, and $\Delta E = E(0)-E_\mathrm{plateau}$ has been read as "energy
absorbed" — giving stopping powers too large to be physical. The 2026-07-27
analysis (see `hypotheses/qsp_phase5/qsp_phase5_momentum_stopping.ipynb`)
identified the culprits: the packet's **localisation (zero-point) energy** and
capture/binding release, neither of which a classical projectile carries. This
pair — **clean launch** (no momentum aliasing, unlike qsp_phase5's fast runs),
matched **cap/nocap** twins, and the **pairwise energy decomposition**
(`interactions.csv`) — allows what qsp_phase5 could not: a *term-by-term
dissection of the plateau*, including a direct test of the naive estimator

$$S_\mathrm{naive}\;\;\text{from}\;\;\Delta E_\mathrm{plateau}-T_\mathrm{loc}(0)
\quad\text{(assumes the internal energy is conserved through the collision)}$$

against the measured orbital-free momentum-KE-loss stopping $S_\mathrm{drift}$.

| run | CAP | what it provides |
|---|---|---|
| `cap` | η = −0.7 Ha, $60<\|z\|<70$ | the plateau; TOF channels; drain bookkeeping |
| `nocap` | off | energy-conservation control; TOF cross-check (pre-wrap window) |
"""))

    # == 2. Conventions ======================================================
    cells.append(md(r"""## Conventions and symbols

Hartree atomic units; $1\,\mathrm{Ha}=27.211\,\mathrm{eV}$. Numbers to 2 s.f.
(3 where a difference would vanish). Same method symbols as the qsp_phase5
notebook ($\rho$, $J$, $W$, $u$, $N(>u)$, $q$, $T_W$); the σ=1 packet changes
the numerology:

| symbol | meaning | value |
|---|---|---|
| $\sigma$ | WP width (σ_WP convention) | 1.0 Bohr |
| $\sigma_r$ | density std per axis $=\sigma/\sqrt2$ | 0.707 Bohr |
| $\sigma_p$ | momentum std per axis $=1/(2\sigma_r)$ | 0.707 a.u. |
| $k_0$ | launch momentum | 2.711 a.u. ($E_\mathrm{drift}=100$ eV) |
| $T_\mathrm{loc}(0)$ | localisation energy $3/(8\sigma_r^2)$ | 0.75 Ha = **20.4 eV** |
| $\sigma_{p_z}^2$(min-unc) | clean-launch variance benchmark | 0.50 |
| $L$ | slab thickness (stopping path) | 25 Bohr |
| $z_\mathrm{in},z_\mathrm{out}$ | detector planes | −16.5, +30 Bohr |

**Geometry** (both runs): box $25\times25\times140$, $dx=0.5$; slab
$|z|\le12.5$; launch $z_0=-20.5$; grid Nyquist $k_\mathrm{Ny}=6.28$;
$k_0+3\sigma_p=4.8<k_\mathrm{Ny}$ → **clean by construction**. dt = 0.02,
5000 steps = 100 a.u.; density frames every 0.4 a.u. (251 frames)."""))

    # == 3. Setup ============================================================
    cells.append(setup_cell())

    cells.append(code(r'''import numpy as np, pandas as pd, json, re
import matplotlib.pyplot as plt
from scipy.special import erfc

HYP     = SYS + "/hypotheses/wp_cap_energy_plateau"
RESULTS = SYS + "/scripts/wp_cap_energy_plateau/wp/results"
QSP5    = SYS + "/hypotheses/qsp_phase5"
RUNS    = ["cap", "nocap"]
HA_EV   = 27.211386
DX, DT  = 0.5, 0.02
Z_SLAB, L_SLAB = 12.5, 25.0
Z_IN, Z_OUT = -16.5, 30.0
CAP_IN, CAP_OUT, ETA = 60.0, 70.0, 0.7      # cap run only
K0, SIG_R, SIG_P = 2.7110633403, 1.0/np.sqrt(2.0), np.sqrt(2.0)/2.0
K_NY = np.pi/DX
T_LOC0 = 3.0/(8.0*SIG_R**2)                  # 0.75 Ha
E_GS = -830.0242258                          # slab-only GS [Ha] (shared_gs run_summary)
COL = {"cap": "tab:red", "nocap": "tab:blue"}
print(f"T_loc(0) = {T_LOC0:.3f} Ha = {T_LOC0*HA_EV:.1f} eV;  "
      f"E_drift = {0.5*K0**2*HA_EV:.1f} eV;  k0+3sigma_p = {K0+3*SIG_P:.2f} < k_Ny = {K_NY:.2f}")'''))

    cells.append(md(r"""## Simulation setup — reconstructable record

Both runs are identical except the CAP (`WP_CAP_ETA` env switch in `run.cpp`):
engine `inq-study`, LDA, slab GS with 102 e⁻ (`shared_gs/slab_n102_L25x25x140_w0p5_h0p5`,
$E_\mathrm{GS}=-830.02$ Ha), WP injected as the 75th state
(`wp_state_index=74`, $\int n_\mathrm{wp}=1$, total $\int n=103$ — verified:
$n(0)=n_\mathrm{gs}+n_\mathrm{wp}(0)$ to machine precision)."""))

    cells.append(code(r'''rows = []
for tag in RUNS:
    txt = open(f"{RESULTS}/{tag}/run_summary.txt").read()
    get = lambda k, d=None: (re.search(rf"{k}\s*=\s*([\d.eE+-]+)", txt) or d)
    rows.append(dict(run=tag,
                     cap_eta=float(get("cap_eta_ha").group(1)),
                     n_steps=int(get("n_steps").group(1)),
                     dt=float(get("dt_au").group(1)),
                     dens_every=int(get("dens_every").group(1)),
                     launch_z=float(get("launch_z").group(1)),
                     wall_h=float(get("wall_time_s").group(1))/3600.0))
pd.DataFrame(rows).set_index("run").round(2)'''))

    cells.append(code(r'''fig, ax = plt.subplots(figsize=(9.5, 1.9))
ax.axvspan(-Z_SLAB, Z_SLAB, color="tab:blue", alpha=0.25, label="jellium slab")
ax.axvspan(-CAP_OUT, -CAP_IN, color="tab:red", alpha=0.25, label="CAP (cap run)")
ax.axvspan(CAP_IN, CAP_OUT, color="tab:red", alpha=0.25)
ax.axvline(Z_IN, color="g", ls="--", lw=1.2); ax.axvline(Z_OUT, color="g", ls="--", lw=1.2)
ax.annotate("$z_\\mathrm{in}$", (Z_IN, 0.82), fontsize=9, ha="right")
ax.annotate("$z_\\mathrm{out}$", (Z_OUT, 0.82), fontsize=9)
ax.plot([-20.5], [0.4], "o", color="tab:green")
ax.annotate("launch, $+z$ →", (-20.2, 0.45), fontsize=9)
ax.set_xlim(-70, 70); ax.set_ylim(0, 1); ax.set_yticks([])
ax.set_xlabel("z [Bohr]"); ax.legend(loc="upper right", fontsize=8, ncol=2)
ax.set_title("wp_cap_energy_plateau geometry — 47.5 Bohr vacuum corridor each side", fontsize=10)
plt.tight_layout(); plt.show()'''))

    # == 4. Source files =====================================================
    cells.append(md(r"""## Source files

| role | path (repo-relative) |
|---|---|
| run definition (cap/nocap via env) | `ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau/wp/run.cpp` |
| orchestrator | `scripts/wp_cap_energy_plateau/orchestrate.sh` |
| shared GS | `shared_gs/slab_n102_L25x25x140_w0p5_h0p5/` |
| kinematics cache extractor | `hypotheses/wp_cap_energy_plateau/plateau_kinematics.py` |
| this builder | `hypotheses/wp_cap_energy_plateau/build_plateau_momentum_report.py` |
| method-mother notebook | `hypotheses/qsp_phase5/qsp_phase5_momentum_stopping.ipynb` |
| qsp_phase5 S(u) headlines (σ=0.5) | `hypotheses/qsp_phase5/momentum_stopping_summary.json` |
| per-run data | `scripts/wp_cap_energy_plateau/wp/results/{cap,nocap}/raw/` |

**Ledger units:** `energies.csv` is in **Ha**; `interactions.csv` is in **eV**
(verified: `e_ss+e_pp+e_ps` = `hartree_inq` = `energies.hartree`×27.211).
Pairwise labels: p = WP, s = slab electrons, b = jellium background."""))

    # == 5. QC-1 launch ======================================================
    cells.append(md(r"""## QC-1 — launch integrity

$f_\mathrm{alias}=\tfrac12\mathrm{erfc}[(k_\mathrm{Ny}-k_0)/(\sqrt2\sigma_p)]$
should be negligible here, and the measured launch moments should sit at
$\langle p_z\rangle(0)=k_0=2.711$, $\sigma_{p_z}^2=0.50$ (min-uncertainty for
σ=1). This is the clean-launch condition the qsp_phase5 fast runs violated."""))

    cells.append(code(r'''qcr = []
for tag in RUNS:
    st = pd.read_csv(f"{RESULTS}/{tag}/raw/observables/wp_momentum_stats.csv", comment="#")
    r0 = st.iloc[0]
    qcr.append(dict(run=tag, f_alias_pred=0.5*erfc((K_NY-K0)/(np.sqrt(2)*SIG_P)),
                    pz_mean_0=r0.pz_mean, pz_frac=r0.pz_mean/K0,
                    sig_pz2_0=r0.sigma_pz2, e_kin_0_Ha=r0.e_kin_ha))
qc = pd.DataFrame(qcr).set_index("run")
qc["grade"] = np.where((qc.pz_frac > 0.98) & (np.abs(qc.sig_pz2_0-0.5) < 0.05),
                       "A (clean)", "check!")
qc.round(3)'''))

    # == 6. Visual intuition ================================================
    cells.append(md(r"""## Visual intuition — both runs

Both runs are significant (each answers a distinct question), so both get their
density GIF (reused from `results/<run>/report/`, path-referenced) and their
energetics. First the packets themselves, then **the plateau figure this run-set
exists for**."""))

    if "cap" in GIF:
        cells.append(embed(GIF["cap"], caption="cap run: total density (xz) — the packet is absorbed at the far CAP", width=680))
    if "nocap" in GIF:
        cells.append(embed(GIF["nocap"], caption="nocap run: total density (xz) — closed box, energy must be conserved", width=680))

    cells.append(md(r"""### The energy plateau (and its control)

$\Delta E_\mathrm{total}(t)=E(t)-E(0)$ for both runs, with the WP norm
$N_\mathrm{wp}(t)$ (from `interactions.csv`). The nocap curve is the
**conservation control**: any drift there bounds the numerics. The cap curve
drains as the CAP removes the packet and then **plateaus** — that plateau value
is the raw material of every deposit-based S estimate."""))

    cells.append(code(r'''EN, INTER = {}, {}
for tag in RUNS:
    EN[tag] = pd.read_csv(f"{RESULTS}/{tag}/raw/observables/energies.csv")
    INTER[tag] = pd.read_csv(f"{RESULTS}/{tag}/raw/observables/interactions.csv")

fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
for tag in RUNS:
    e = EN[tag]
    axes[0].plot(e.time_au, (e.total-e.total.iloc[0])*HA_EV, color=COL[tag], label=tag)
    axes[1].plot(INTER[tag].time_au, INTER[tag].N_wp, color=COL[tag], label=tag)
axes[0].axhline(0, color="k", lw=0.5)
axes[0].set_xlabel("t [a.u.]"); axes[0].set_ylabel("$\\Delta E_\\mathrm{total}$ [eV]")
axes[0].set_title("total energy: plateau (cap) vs conservation control (nocap)", fontsize=10)
axes[1].set_xlabel("t [a.u.]"); axes[1].set_ylabel("$N_\\mathrm{wp}$ [e]")
axes[1].set_title("WP norm: CAP drain bookkeeping", fontsize=10)
for ax in axes: ax.legend(fontsize=8)
plt.tight_layout(); plt.show()

e = EN["cap"]; late = e.time_au > 80.0
E_plat = float(e.total[late].mean())
slope = float(np.polyfit(e.time_au[late], e.total[late]*HA_EV, 1)[0])
econ = EN["nocap"]
drift_nocap = float((econ.total.iloc[-1]-econ.total.iloc[0])*HA_EV)
print(f"E(0) = {e.total.iloc[0]:.4f} Ha | E_plateau (t>80) = {E_plat:.4f} Ha | "
      f"late slope {slope:.2e} eV/a.u.")
print(f"nocap conservation drift over 100 a.u.: {drift_nocap:+.2f} eV  (numerics bound)")
N_wp_end = float(INTER["cap"].N_wp.iloc[-1])
print(f"cap: N_wp(end) = {N_wp_end:.3f} -> absorbed fraction {1-N_wp_end:.3f}")'''))

    # == 7. TOF extraction ==================================================
    cells.append(md(r"""## Orbital-free TOF extraction (method of the qsp_phase5 notebook)

Same validated pipeline, this geometry: $\rho(z,t)$ from `density_total` minus
the GS field (cached by `plateau_kinematics.py`); side-adaptive continuity flux

$$J(z,t)=-\!\!\int_{-L/2}^{z}\!\Big(\partial_t\rho+2W\rho\Big)dz'
\quad\text{(entrance)},\qquad
J(z,t)=+\!\!\int_{z}^{+L/2}\!\Big(\partial_t\rho+2W\rho\Big)dz'
\quad\text{(exit)}$$

with $W(z)=0.7\sin^2[\pi(|z|-60)/10]$ on $60<|z|<70$ for the cap run and
$W\equiv0$ for nocap. **nocap validity window:** with no CAP the edge condition
$J(\pm L/2)=0$ fails once the packet's front reaches the box edge — the window
is detected from the data and later times are excluded. Closure test for the
cap run as before."""))

    cells.append(code(r'''KIN = {tag: {k: v for k, v in np.load(f"{HYP}/cache/{tag}_kinematics.npz").items()}
       for tag in RUNS}

def cap_W(z, on):
    W = np.zeros_like(z)
    if on:
        m = (np.abs(z) > CAP_IN) & (np.abs(z) < CAP_OUT)
        W[m] = ETA*np.sin(np.pi*(np.abs(z[m]) - CAP_IN)/(CAP_OUT - CAP_IN))**2
    return W

def flux_both(rho, z, t, W):
    dz = z[1] - z[0]
    src = np.gradient(rho, t, axis=0) + 2.0*W[None, :]*rho
    J_L = -np.cumsum(src, axis=1)*dz
    J_R = np.cumsum(src[:, ::-1], axis=1)[:, ::-1]*dz
    A_p = (2.0*W[None, :]*rho)[:, z >  CAP_IN].sum(axis=1)*dz
    A_m = (2.0*W[None, :]*rho)[:, z < -CAP_IN].sum(axis=1)*dz
    return J_L, J_R, A_p, A_m

FLUX, T_VALID = {}, {}
for tag in RUNS:
    kin = KIN[tag]
    FLUX[tag] = flux_both(kin["rho"], kin["z"], kin["t_au"], cap_W(kin["z"], tag == "cap"))
    edge = np.abs(kin["rho"][:, np.abs(kin["z"]) > 68.0]).sum(axis=1)*DX
    bad = np.where(edge > 2e-3)[0]
    T_VALID[tag] = float(kin["t_au"][bad[0]]) if len(bad) else float(kin["t_au"][-1])
print("validity windows [a.u.]:", {k: round(v, 1) for k, v in T_VALID.items()})

# closure QC (cap run)
kin = KIN["cap"]; _, _, A_p, A_m = FLUX["cap"]
dz = kin["z"][1]-kin["z"][0]
dNdt = np.gradient(kin["rho"].sum(axis=1)*dz, kin["t_au"])
fig, ax = plt.subplots(figsize=(6.5, 2.8))
ax.plot(kin["t_au"], -dNdt, "k", lw=1.4, label="$-dN/dt$ (measured)")
ax.plot(kin["t_au"], A_p+A_m, color="tab:red", ls="--", label="$A_++A_-$ (sink model)")
resid = np.abs(-dNdt-(A_p+A_m)).mean()/max(np.abs(dNdt).max(), 1e-12)
ax.set_xlabel("t [a.u.]"); ax.set_ylabel("norm rate [e/a.u.]")
ax.set_title(f"cap-run closure: residual {resid:.0%}", fontsize=10); ax.legend(fontsize=8)
plt.tight_layout(); plt.show()'''))

    cells.append(code(r'''# rho(z,t) carpet, cap run, with geometry annotated
kin = KIN["cap"]
fig, ax = plt.subplots(figsize=(8.5, 4))
v = np.abs(kin["rho"]).max()*0.5
im = ax.pcolormesh(kin["t_au"], kin["z"], kin["rho"].T, cmap="RdBu_r", vmin=-v, vmax=v)
for zl, c in ((Z_SLAB, "k"), (-Z_SLAB, "k"), (Z_IN, "g"), (Z_OUT, "g"),
              (CAP_IN, "r"), (-CAP_IN, "r")):
    ax.axhline(zl, color=c, ls="--", lw=0.8)
plt.colorbar(im, ax=ax, label=r"$\rho(z,t)$ [e/Bohr]")
ax.set_xlabel("t [a.u.]"); ax.set_ylabel("z [Bohr]")
ax.set_title("cap: excess density — slab (black), detectors (green), CAP faces (red)", fontsize=10)
plt.tight_layout(); plt.show()'''))

    # == 8. Rank-matched S(u) + null ========================================
    cells.append(md(r"""### Rank-matched $S(u)$ and the geometry-specific null test

TOF readout $u(t)=J/\rho$ at the planes; signed-flux channels; exceedance curves
$N(>u)$; equal-rank matching $u_\mathrm{in}(q)\mapsto u_\mathrm{out}(q)$;
$S=\tfrac12[u_\mathrm{in}^2-u_\mathrm{out}^2]/L$ (full derivation and failure
analysis of the naive per-electron difference: qsp_phase5 notebook, step 4).
The **null test** (free σ=1 packet, this launch/planes/cadence, extended domain,
identical code path) recalibrates the trusted rank window and systematic for
*this* geometry — the entrance flight here is short (4 Bohr), which the null
must price."""))

    cells.append(code(r'''def sweep(kin, J, z_plane, sgn, vmax, t_max=None):
    z, t = kin["z"], kin["t_au"]
    i = int(np.argmin(np.abs(z - z_plane)))
    Jp, rp = J[:, i], kin["rho"][:, i]
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.where(np.abs(rp) > 1e-7, Jp/rp, 0.0)
    u = np.clip(u, -vmax, vmax)
    w = np.where(sgn*Jp > 0, sgn*Jp, 0.0)*np.gradient(t)
    if t_max is not None:
        w = np.where(t <= t_max, w, 0.0)
    return sgn*u, w

def s_of_u(kin, flx, t_max=None, n_rank=60):
    vmax = min(K0 + 4.0*SIG_P, K_NY)
    u_i, w_i = sweep(kin, flx[0], Z_IN, +1, vmax, t_max)
    u_o, w_o = sweep(kin, flx[1], Z_OUT, +1, vmax, t_max)
    u_r, w_r = sweep(kin, flx[0], Z_IN, -1, vmax, t_max)
    ug = np.linspace(0.15, vmax, 240)
    Ein = np.array([w_i[u_i > g].sum() for g in ug])
    Eout = np.array([w_o[u_o > g].sum() for g in ug])
    q_top = 0.92*min(Ein[0], Eout[0])
    qs = np.linspace(0.02, q_top, n_rank)
    uin = np.interp(qs, Ein[::-1], ug[::-1])
    uout = np.interp(qs, Eout[::-1], ug[::-1])
    S = 0.5*(uin**2 - uout**2)*HA_EV/L_SLAB
    trust = (qs > 0.30*q_top) & (qs < 0.90*q_top)
    return dict(ug=ug, Ein=Ein, Eout=Eout, qs=qs, uin=uin, uout=uout, S=S,
                trust=trust, N_in=Ein[0], N_out=Eout[0], N_refl=w_r.sum())

# null test in THIS geometry (free sigma=1 packet, extended domain)
z_n = np.arange(-1100.25, 1100.5, DX)
t_n = np.arange(0.0, 100.0, 0.4) + 1e-9
s2 = SIG_R**2 + (SIG_P*t_n)**2
rho_n = np.exp(-(z_n[None, :]+20.5-K0*t_n[:, None])**2/(2*s2[:, None]))/np.sqrt(2*np.pi*s2[:, None])
kin_n = dict(z=z_n, t_au=t_n, rho=rho_n)
flx_n = flux_both(rho_n, z_n, t_n, np.zeros_like(z_n))
null = s_of_u(kin_n, flx_n)
S_SYST = float(np.abs(null["S"][null["trust"]]).max())

SOU = {tag: s_of_u(KIN[tag], FLUX[tag], t_max=T_VALID[tag]) for tag in RUNS}

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(null["uin"], null["S"], "o-", ms=3, color="0.6", label="null, all ranks")
axes[0].plot(null["uin"][null["trust"]], null["S"][null["trust"]], "o", ms=4,
             color="tab:green", label="trusted window")
axes[0].axhline(0, color="k", lw=0.8)
axes[0].set_xlabel("$u_\\mathrm{in}$ [a.u.]"); axes[0].set_ylabel("$S_\\mathrm{null}$ [eV/Bohr]")
axes[0].set_title(f"null test (σ=1 geometry): syst ±{S_SYST:.2f} eV/Bohr", fontsize=10)
axes[0].legend(fontsize=8)
for tag in RUNS:
    r = SOU[tag]; tr = r["trust"]
    axes[1].plot(r["uin"][~tr], r["S"][~tr], ".", ms=3, color=COL[tag], alpha=0.25)
    axes[1].plot(r["uin"][tr], r["S"][tr], "o-", ms=3, color=COL[tag], label=tag)
axes[1].axhspan(-S_SYST, S_SYST, color="0.85", zorder=0)
axes[1].axhline(0, color="k", lw=0.6)
axes[1].set_xlabel("$u_\\mathrm{in}$ [a.u.]"); axes[1].set_ylabel("$S$ [eV/Bohr]")
axes[1].set_title("rank-matched S(u): cap vs nocap must agree (independent absorbers)", fontsize=10)
axes[1].legend(fontsize=8)
plt.tight_layout(); plt.show()

head = {}
for tag in RUNS:
    r = SOU[tag]; tr = r["trust"]
    head[tag] = dict(u_ref=float(np.mean(r["uin"][tr])),
                     S_drift=float(np.mean(r["S"][tr])),
                     S_err=float(np.hypot(np.std(r["S"][tr]), S_SYST)),
                     N_in=float(r["N_in"]), N_trans=float(r["N_out"]),
                     N_refl=float(r["N_refl"]))
pd.DataFrame(head).T.round(2)'''))

    # == 9. Plateau dissection ==============================================
    cells.append(md(r"""## The centrepiece — term-by-term dissection of the plateau

Every energy in one ledger (eV per injected electron; ledger identities from the
conversation of 2026-07-27):

$$E_\mathrm{arrival} = \underbrace{\tfrac12 k_0^2}_{\text{drift, 100 eV}}
+\underbrace{3/(8\sigma_r^2)}_{T_\mathrm{loc},\ 20.4\ \mathrm{eV}}
+\underbrace{\text{self/interaction terms}}_{\text{from the DFT ledger}}$$

$$R \equiv E(0)-E_\mathrm{plateau}\ \text{(CAP-removed)},\qquad
D \equiv E_\mathrm{plateau}-E_\mathrm{arrival\ ledger\ start}\ \text{(kept by the slab)}$$

The **naive estimator** under test — "the localisation energy that existed
initially also exists finally, so subtract it once":

$$S_\mathrm{naive}=\frac{E_\mathrm{drift}(0)-\big[R-T_\mathrm{loc}(0)\big]}{L}
\qquad\text{vs the measured}\qquad
S_\mathrm{drift}\ \text{(TOF, above)}$$

The DFT arrival ledger is checkable here: $E(0)-E_\mathrm{GS}$ should equal
drift + $T_\mathrm{loc}$ + (WP self-Hartree `e_pp(0)` + xc self-interaction +
WP–slab/background interaction at launch) — the pairwise table supplies each
piece."""))

    cells.append(code(r'''e = EN["cap"]; icap = INTER["cap"]
E0, e0 = float(e.total.iloc[0]), icap.iloc[0]
arrival_ledger = (E0 - E_GS)*HA_EV                    # what the WP brought, per DFT
R = (E0 - E_plat)*HA_EV                               # CAP-removed
D = arrival_ledger - R                                # kept by the slab
drift0, loc0 = 0.5*K0**2*HA_EV, T_LOC0*HA_EV

ledger = pd.DataFrame([
    ("E_drift(0) = k0^2/2",                 drift0),
    ("T_loc(0) = 3/(8 sigma_r^2)",          loc0),
    ("e_pp(0)  (WP self-Hartree)",          float(e0.e_pp)),
    ("e_ps(0)+e_pb(0) (WP-slab/bg at launch)", float(e0.e_ps + e0.e_pb)),
    ("sum of the above",                    drift0 + loc0 + float(e0.e_pp) + float(e0.e_ps + e0.e_pb)),
    ("E(0) - E_GS  (DFT arrival ledger)",   arrival_ledger),
    ("R = E(0) - E_plateau (CAP-removed)",  R),
    ("D = arrival - R (kept by slab)",      D),
], columns=["term", "eV"]).set_index("term")
ledger.round(1)'''))

    cells.append(md(r"""*(The 'sum of the above' vs 'DFT arrival ledger' row-pair is the closure of the
arrival budget; their residual is the xc self-interaction + charged-cell
convention terms — the parts of the one-electron ledger LDA gets wrong, worth
knowing but not part of any stopping definition.)*

Now the estimator shoot-out. Note $R-T_\mathrm{loc}(0)$ is *supposed* to be the
surviving drift KE under the naive assumption; the TOF exit channel *measures*
the surviving longitudinal KE per transmitted electron — so the comparison
directly prices the "internal energy is conserved" assumption, and the pairwise
`e_pp(t)` (WP self-Hartree, tracking $N_\mathrm{wp}$ and the packet's shape)
shows *where* the internal energy actually went."""))

    cells.append(code(r'''r = SOU["cap"]; tr = r["trust"]
S_drift = head["cap"]["S_drift"]
K_out_per_e = float((0.5*np.interp(np.mean(r["qs"][tr]), r["qs"], r["uout"]**2))*HA_EV)  # at mean trusted rank
S_naive = (drift0 - (R - loc0))/L_SLAB
S_deposit = D/L_SLAB
res = pd.DataFrame([
    ("S_drift (TOF rank-matched, measured)", S_drift, head["cap"]["S_err"]),
    ("S_naive = [E_drift0 - (R - T_loc0)]/L", S_naive, np.nan),
    ("S_deposit = D/L (no correction)",       S_deposit, np.nan),
], columns=["estimator", "S [eV/Bohr]", "err"]).set_index("estimator")
print(f"R - T_loc(0) = {R - loc0:.1f} eV  vs  TOF surviving KE/e (trusted mean) = {K_out_per_e:.1f} eV")
res.round(2)'''))

    cells.append(code(r'''# where the internal energy went: pairwise terms + T_W lobes vs time
fig, axes = plt.subplots(1, 2, figsize=(12, 3.8))
for tag in RUNS:
    i = INTER[tag]
    axes[0].plot(i.time_au, i.e_pp - i.e_pp.iloc[0], color=COL[tag], label=f"{tag}: $\\Delta e_{{pp}}$")
    axes[0].plot(i.time_au, i.e_ps - i.e_ps.iloc[0], color=COL[tag], ls="--", label=f"{tag}: $\\Delta e_{{ps}}$")
axes[0].set_xlabel("t [a.u.]"); axes[0].set_ylabel("eV")
axes[0].set_title("pairwise: WP self-Hartree ($e_{pp}$) and WP–slab ($e_{ps}$)", fontsize=10)
axes[0].legend(fontsize=7)
for tag in RUNS:
    kin = KIN[tag]
    axes[1].plot(kin["t_au"], kin["TW_lo"]*HA_EV, color=COL[tag], label=f"{tag}: entrance lobe")
    axes[1].plot(kin["t_au"], kin["TW_hi"]*HA_EV, color=COL[tag], ls="--", label=f"{tag}: exit lobe")
tt = KIN["cap"]["t_au"]
axes[1].plot(tt, 3.0/(8.0*(SIG_R**2 + (SIG_P*tt)**2))*HA_EV, "k:", lw=1.2,
             label="free Gaussian $3/8\\sigma_r^2(t)$")
axes[1].set_xlabel("t [a.u.]"); axes[1].set_ylabel("$T_W$ [eV]")
axes[1].set_title("localisation (shape) energy in the vacuum lobes", fontsize=10)
axes[1].legend(fontsize=7)
plt.tight_layout(); plt.show()
print(f"T_W(0) entrance lobe = {KIN['cap']['TW_lo'][0]*HA_EV:.1f} eV  vs  analytic {loc0:.1f} eV")'''))

    cells.append(md(r"""*(Charged-cell caveat: individual Hartree/external pairwise values are
$G=0$-convention dependent, and the cell's net charge* **changes in time** *as
the CAP drains the WP — so pairwise* differences *are physically meaningful only
while $N_\mathrm{wp}$ is roughly constant, and trend-level once draining starts.
The total-energy plateau itself is unaffected by the split.)*"""))

    # == 10. Synthesis =======================================================
    cells.append(md(r"""## Synthesis — this pair on the S(u) map, against qsp_phase5 (σ = 0.5)

The σ=1 point lands on the same axes as the qsp_phase5 σ=0.5 curves. They need
not coincide — stopping of an extended charge is σ-dependent (the packet spans
several screening lengths) — but cap/nocap must agree with each other, and the
deposit-based points must sit above the drift-based ones by the localisation +
capture margin."""))

    cells.append(code(r'''with open(f"{QSP5}/momentum_stopping_summary.json") as f:
    q5 = json.load(f)
cls = pd.read_csv(f"{QSP5}/classical_sigma0p5_bulk.csv")
lin = np.load(f"{QSP5}/lindhard_ref.npz")

fig, ax = plt.subplots(figsize=(8, 4.6))
ax.plot(cls.v, cls.S_eVbohr, "s-", color="0.55", ms=4, label="classical σ=0.5 (bulk)")
ax.plot(np.sqrt(2*lin["E"]/HA_EV), lin["S"], "-", color="0.8", lw=1.2, label="Lindhard (point)")
for t5, d in q5.items():
    solid = d["grade"].startswith(("A", "B"))
    ax.errorbar([d["u_ref"]], [d["S_drift"]], yerr=[d["S_err"]], fmt="D",
                color="0.35", ms=5, mfc=("0.35" if solid else "none"), lw=1)
ax.plot([], [], "D", color="0.35", label="qsp_phase5 $S_\\mathrm{drift}$ (σ=0.5)")
for tag in RUNS:
    r = SOU[tag]; tr = r["trust"]
    ax.plot(r["uin"][tr], r["S"][tr], "-", color=COL[tag], lw=2)
    ax.errorbar([head[tag]["u_ref"]], [head[tag]["S_drift"]], yerr=[head[tag]["S_err"]],
                fmt="o", color=COL[tag], ms=8, label=f"{tag} $S_\\mathrm{{drift}}$ (σ=1)")
ax.plot([head["cap"]["u_ref"]], [S_naive], "*", color="tab:purple", ms=14,
        label="$S_\\mathrm{naive}$ ($R-T_\\mathrm{loc}$ corrected)")
ax.plot([head["cap"]["u_ref"]], [S_deposit], "^", color="tab:red", ms=9, mfc="none",
        label="$S_\\mathrm{deposit}$ (uncorrected)")
ax.set_xlabel("$u$ [a.u.]"); ax.set_ylabel("S [eV/Bohr]"); ax.set_xlim(0, 6.5)
ax.set_title("stopping estimates: σ=1 plateau pair on the σ=0.5 sweep map", fontsize=10)
ax.legend(fontsize=7)
plt.tight_layout(); plt.show()'''))

    cells.append(code(r'''summary = dict(
    E0_Ha=E0, E_plateau_Ha=E_plat, E_GS_Ha=E_GS,
    arrival_ledger_eV=arrival_ledger, R_eV=R, D_eV=D,
    drift0_eV=drift0, T_loc0_eV=loc0, e_pp0_eV=float(e0.e_pp),
    R_minus_Tloc_eV=R-loc0, K_out_per_e_eV=K_out_per_e,
    S_drift=head["cap"]["S_drift"], S_err=head["cap"]["S_err"],
    S_drift_nocap=head["nocap"]["S_drift"], S_err_nocap=head["nocap"]["S_err"],
    S_naive=S_naive, S_deposit=S_deposit, S_syst=S_SYST,
    u_ref=head["cap"]["u_ref"], N_wp_end=N_wp_end,
    nocap_drift_eV=drift_nocap, t_valid_nocap=T_VALID["nocap"])
with open(f"{HYP}/plateau_dissection_summary.json", "w") as f:
    json.dump(summary, f, indent=1)
print("wrote plateau_dissection_summary.json")'''))

    # == 11. Takeaway ========================================================
    if summary:
        s = summary
        agree = abs(s["S_drift"]-s["S_drift_nocap"]) <= (s["S_err"]+s["S_err_nocap"])
        takeaway = f"""## Takeaway

| estimator | S [eV/Bohr] |
|---|---|
| $S_\\mathrm{{drift}}$ (TOF, cap) | **{s['S_drift']:.2g} ± {s['S_err']:.2g}** |
| $S_\\mathrm{{drift}}$ (TOF, nocap control) | {s['S_drift_nocap']:.2g} ± {s['S_err_nocap']:.2g} |
| $S_\\mathrm{{naive}}$ ($R-T_\\mathrm{{loc}}$ corrected) | {s['S_naive']:.2g} |
| $S_\\mathrm{{deposit}}$ (uncorrected plateau) | {s['S_deposit']:.2g} |

- **Plateau ledger** (per injected electron): arrival {s['arrival_ledger_eV']:.0f} eV
  (= 100 drift + {s['T_loc0_eV']:.0f} localisation + {s['e_pp0_eV']:.0f} self-Hartree
  + launch interaction); CAP removed R = {s['R_eV']:.0f} eV; kept by slab
  D = {s['D_eV']:.0f} eV.
- **The naive correction is testable here and {'lands within errors of' if abs(s['S_naive']-s['S_drift'])<=s['S_err'] else 'differs from'} the measured
  $S_\\mathrm{{drift}}$** — $R-T_\\mathrm{{loc}}(0) = {s['R_minus_Tloc_eV']:.0f}$ eV vs TOF surviving
  KE/e ≈ {s['K_out_per_e_eV']:.0f} eV: the "internal energy is conserved" assumption is
  priced quantitatively for a clean traversal at u ≈ {s['u_ref']:.2g}.
- cap and nocap TOF estimates {'agree within errors' if agree else 'DISAGREE — investigate'} —
  independent absorbing boundaries, same physics window.
- nocap conservation control: |ΔE_total| = {abs(s['nocap_drift_eV']):.2g} eV over 100 a.u.
  bounds the ledger numerics; nocap TOF valid to t = {s['t_valid_nocap']:.0f} a.u. (edge wrap).
- Systematic (null test, this geometry): ±{s['S_syst']:.2g} eV/Bohr.
- σ-dependence caveat: these are σ = 1 numbers on a σ = 0.5 map — same-u
  differences vs qsp_phase5 conflate physics (extended-charge screening) with σ.
"""
        cells.append(md(takeaway))
    else:
        cells.append(md("## Takeaway\n\n*(populated on the second build pass from "
                        "`plateau_dissection_summary.json`)*"))
    return cells


summary = None
if os.path.exists(SUMMARY):
    with open(SUMMARY) as f:
        summary = json.load(f)

print("pass 1: executing notebook ...")
build(build_cells(summary), OUT, timeout=1800)

with open(SUMMARY) as f:
    summary2 = json.load(f)
if summary2 != summary:
    print("pass 2: takeaway numbers changed — re-rendering ...")
    build(build_cells(summary2), OUT, timeout=1800)
print("done:", OUT)
