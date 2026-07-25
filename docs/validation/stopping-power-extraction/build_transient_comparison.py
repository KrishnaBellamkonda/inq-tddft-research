#!/usr/bin/env python3
"""Deterministic builder for transient_method_comparison.ipynb.

Compares two ways of choosing the transient cutoff x0 when extracting the
electronic stopping power S = d(E_total)/dx of a classical projectile:

  METHOD-1  fixed-fraction       : discard the first 20% of the run, free-intercept fit.
  METHOD-2  slope-plateau agent  : detect_x0_and_stopping_power (user-supplied spec):
                                    sweep x0, find where S(x0) plateaus within tolerance,
                                    40% hard gate, status flags.

Both use a FREE-INTERCEPT fit (dE = S*x + E0). The comparison therefore isolates
ONLY the transient-cut choice. The decision (which method) is the USER's; this
notebook only presents the evidence neutrally.

Run (venv):
  /local/data/public/skcb2/tddft/venv/bin/python3 build_transient_comparison.py

Writes + EXECUTES transient_method_comparison.ipynb (kernel inqview-venv).
"""
from __future__ import annotations
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "transient_method_comparison.ipynb")

cells = []
def md(s): cells.append(new_markdown_cell(s))
def code(s): cells.append(new_code_cell(s))

# ----------------------------------------------------------------------------- 0
md(r"""# Transient-cutoff method comparison — stopping power $S = dE_{\rm total}/dx$

**Decision to be made (by the user):** which way of choosing the transient cutoff
$x_0$ is most applicable for this project —

| | method | rule |
|---|---|---|
| **M1** | fixed-fraction | discard the first **20 %** of the run, then fit |
| **M2** | slope-plateau "agent" | detect $x_0$ where the fitted slope $S(x_0)$ stops depending on $x_0$; **40 %** hard gate; status flags |

Both fit the *same* free-intercept line $\Delta E_{\rm total}=S\,x+E_0$ over
$[x_0,x_T]$, so the **only** thing being compared is the transient-cut choice.
This notebook presents the evidence; it does **not** pick a winner.

**Physics (grounded).** For a classical Ehrenfest projectile the electronic
energy gained equals the projectile's kinetic-energy loss, and its steady rate
per unit path is the electronic stopping power
($S=\langle dE/dt\rangle/v$, **Correa 2018 Eq. 10**; transient excluded per
**Correa Fig. 8**, `docs/sources/correa-2018-electronic-stopping-power.md`).
The independent variable is the projectile **displacement** $x(t)$ from the
track (velocity is not exactly constant under Ehrenfest).

**Data.** Classical erf-Gaussian projectiles, $r_s=5.69$ jellium ($N=162$,
$L=50$, $dx=0.40$), $\sigma$-sweep runs. $\sigma$ in the run name = projectile
**charge** std $\sigma_q$ ($\sigma_{\rm WP}=\sqrt2\,\sigma_q$, so $\sigma_q{=}0.35\Leftrightarrow\sigma_{\rm WP}{=}0.5$).
Reference: point-charge Lindhard `lindhard_elf.stopping_power_point`.""")

# ----------------------------------------------------------------------------- 1 kernel
md(r"""## 1 — The extraction kernel (both methods + sanity channels)

Everything below is computed from these functions. They are the candidate
content of the future `stopping-power-extraction` skill, so they are shown in
full (no black boxes).""")

code(r'''import sys, os, csv
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
import numpy as np
import matplotlib.pyplot as plt
from inqview.analysis.stopping_extract import load_track
from inqview.analysis import lindhard_elf as LE
from inqview.visualisation import style as ST
ST.apply_theme()

JELL = "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium"
RS   = 5.69
KF   = LE.kF_from_rs(RS)           # Fermi wavevector = vF
MASS = 1.0                         # projectile = m_e
HA_PER_BOHR_TO_EV_PER_A = 27.211386245988 / 0.529177210903   # 1 Ha/Bohr -> eV/Angstrom

# --- run inventory: (sigma_q, run_dir, layout). "nested" = results/<v>/raw/observables/;
#     "flat" = observables directly under results/<v>/ (the sigma=0.5 anchor set). ---
SIGMAS = [
    (0.15, "run_classical_n162_L50_sv_sigma0p15", "nested"),
    (0.25, "run_classical_n162_L50_sv_sigma0p25", "nested"),
    (0.35, "run_classical_n162_L50_sv_sigma0p35", "nested"),
    (0.50, "run_sv_sigma0p5",                     "flat"),
    (3.00, "run_classical_n162_L50_sv_sigma3p0",  "nested"),
]
VTAGS = [("v0p2", 0.2), ("v0p6", 0.6), ("v0p8", 0.8), ("v1p0", 1.0),
         ("v1p3", 1.3), ("v2p0", 2.0), ("v3p0", 3.0)]

def _obsdir(run_dir, vtag, layout):
    base = os.path.join(JELL, run_dir, "results", vtag)
    return base if layout == "flat" else os.path.join(base, "raw", "observables")

def load_obs(obsdir):
    """(t_au, E_total) from observables.csv, deduped by step, time-sorted."""
    t, Et, seen = [], [], set()
    with open(os.path.join(obsdir, "observables.csv")) as fh:
        for row in csv.DictReader(fh):
            if row.get("energy_total") in (None, ""):
                continue
            if row["step"] in seen:
                continue
            seen.add(row["step"])
            t.append(float(row["time_au"])); Et.append(float(row["energy_total"]))
    t = np.array(t); Et = np.array(Et); o = np.argsort(t)
    return t[o], Et[o]

def load_run(run_dir, vtag, layout):
    """Return dict(x, Ed, t_o, tr) or None. x = projectile displacement at the
    observable times; Ed = E_total(t)-E_total(0); tr = full fine-cadence Track."""
    obsdir = _obsdir(run_dir, vtag, layout)
    if not os.path.exists(os.path.join(obsdir, "observables.csv")):
        return None
    t_o, Et = load_obs(obsdir)
    tr = load_track(os.path.join(obsdir, "electron_track.csv"), mass=MASS, axis="z")
    if t_o.size < 6 or tr.t.size < 6:
        return None
    x = np.interp(t_o, tr.t, tr.s)          # displacement at the energy samples
    return dict(x=x, Ed=Et - Et[0], t_o=t_o, tr=tr)

# --- free-intercept least squares: dE = S*x + E0 ---------------------------------
def free_fit(x, E, x0, xT):
    m = (x >= x0) & (x <= xT); n = int(m.sum())
    if n < 5:
        return None
    A = np.vstack([x[m], np.ones(n)]).T
    (S, c), *_ = np.linalg.lstsq(A, E[m], rcond=None)
    r = E[m] - (S * x[m] + c)
    dof = max(n - 2, 1)
    sxx = np.sum((x[m] - x[m].mean()) ** 2)
    se = np.sqrt(np.sum(r ** 2) / dof / sxx) if sxx > 0 else np.inf
    ybar = E[m].mean(); sst = np.sum((E[m] - ybar) ** 2)
    r2 = 1.0 - np.sum(r ** 2) / sst if sst > 0 else np.nan
    return dict(S=float(S), E0=float(c), se=float(se), r2=float(r2), n=n,
                x0=float(x0), xT=float(xT), resid=r, xfit=x[m])

# --- METHOD 1: fixed fraction ----------------------------------------------------
def fixed_fraction(x, E, frac=0.20, xT=None):
    xT = x.max() if xT is None else xT
    x0 = x.min() + frac * (x.max() - x.min())
    f = free_fit(x, E, x0, xT)
    if f is None:
        return dict(status="range_too_short")
    f["status"] = "ok"; f["frac"] = frac
    return f

# --- METHOD 2: slope-plateau agent (user spec, verbatim logic) -------------------
def detect_x0_and_stopping_power(x, E, xT=None, remain_min=0.30, rel_tol=0.02,
                                 k_sigma=2.0, resid_ratio_max=2.0, gate=0.40,
                                 grid_step=None):
    x = np.asarray(x, float); E = np.asarray(E, float)
    o = np.argsort(x); x, E = x[o], E[o]
    x_min, x_max = x.min(), x.max(); L = x_max - x_min
    if xT is None: xT = x_max
    if grid_step is None: grid_step = max(np.median(np.diff(x)), L / 400)
    status = "ok"
    # endpoint (image re-entry) check
    end = free_fit(x, E, 0.7 * (xT - x_min) + x_min, xT)
    mid = free_fit(x, E, x_min + 0.4 * (xT - x_min), x_min + 0.7 * (xT - x_min))
    if end and mid and (end["S"] - mid["S"]) > rel_tol * abs(mid["S"]) + 3 * mid["se"]:
        status = "endpoint_contaminated"
    # sweep
    x0cap = x_min + (1 - remain_min) * (xT - x_min)
    grid = np.arange(x_min, x0cap, grid_step)
    fits = [free_fit(x, E, g, xT) for g in grid]
    keep = [f is not None for f in fits]
    grid = grid[keep]; fits = [f for f, k in zip(fits, keep) if k]
    if len(grid) < 3:
        return dict(status="range_too_short", x0=None, grid=grid)
    S = np.array([f["S"] for f in fits]); serr = np.array([f["se"] for f in fits])
    rem = (xT - grid) / (xT - x_min); ref = (rem >= 0.30) & (rem <= 0.55)
    if ref.sum() < 3: ref = np.ones_like(rem, bool)
    Spl = float(np.median(S[ref]))
    tol = max(rel_tol * abs(Spl), k_sigma * float(np.median(serr[ref])))
    ok = np.abs(S - Spl) <= tol
    x0 = None
    for i in range(len(grid)):
        if ok[i] and ok[i:].all():
            x0 = float(grid[i]); break
    if x0 is None:
        return dict(status="no_plateau", x0=None, Spl=Spl, tol=float(tol),
                    grid=grid, Sgrid=S, segrid=serr)
    f = free_fit(x, E, x0, xT)
    lead = abs(f["resid"][0]) / np.sqrt(np.mean(f["resid"] ** 2))
    tf = (x0 - x_min) / L
    gate_pass = tf <= gate
    if not gate_pass: status = "range_too_short"
    f.update(status=status, x0=x0, Spl=Spl, tol=float(tol), transient_fraction=float(tf),
             steady_fraction=float((xT - x0) / L), lead_resid_ratio=float(lead),
             gate_pass=bool(gate_pass), grid=grid, Sgrid=S, segrid=serr)
    return f

# --- sanity channels: kinetic (-dKE/dx) and cumulative power int F.v -------------
def sanity_channels(run, x0, xT):
    """Return dicts for the kinetic and force-power channels, fitted on [x0,xT]
    in displacement. F = m*dv/dt from the track (so int F.v dt == Delta KE
    analytically: a consistency/profile check, NOT an independent physics channel)."""
    tr = run["tr"]
    dKE_loss = tr.ke[0] - tr.ke                          # energy lost by projectile (Ha)
    a = np.gradient(tr.v, tr.t)                          # dv/dt
    P = MASS * a * tr.v                                  # power on projectile = dKE/dt (<0)
    work_on_proj = np.concatenate([[0.0], np.cumsum(0.5 * (P[1:] + P[:-1]) * np.diff(tr.t))])
    dep_from_force = -work_on_proj                       # energy deposited to medium
    # interpolate both onto the energy-sample displacements, fit on the window
    x_obs = run["x"]
    dKE_obs = np.interp(x_obs, tr.s, dKE_loss)
    dF_obs  = np.interp(x_obs, tr.s, dep_from_force)
    fKE = free_fit(x_obs, dKE_obs, x0, xT)
    fF  = free_fit(x_obs, dF_obs, x0, xT)
    return dict(x_obs=x_obs, dKE_obs=dKE_obs, dF_obs=dF_obs, fKE=fKE, fF=fF,
                tr_s=tr.s, dKE_loss=dKE_loss, dep_from_force=dep_from_force)

print("kernel ready.  kF = vF =", round(KF, 4), " omega_p =", round(LE.omega_p(KF), 4),
      "a.u.  tau_p =", round(2*np.pi/LE.omega_p(KF), 1), "a.u.")
''')

# ----------------------------------------------------------------------------- 2 detailed
md(r"""## 2 — Detailed walk: $\sigma_q=0.35$, $v_0=3.0$ (fast, clean)

We open one run and look at the raw signal, then run both methods on it. $v_0=3.0$
is the cleanest case: $v$ falls only ~1.4 %, and energy conservation
$\Delta E_{\rm total}\approx\Delta\mathrm{KE}_{\rm loss}$ holds to <1 %.""")

code(r'''run = load_run("run_classical_n162_L50_sv_sigma0p35", "v3p0", "nested")
x, Ed = run["x"], run["Ed"]
L = x.max() - x.min()
print(f"displacement span L = {L:.2f} Bohr over {len(x)} energy samples")
print(f"dE_total total      = {Ed[-1]:.4f} Ha")
print(f"dKE_loss  total     = {run['tr'].ke[0]-run['tr'].ke[-1]:.4f} Ha   "
      f"(energy-conservation check: agree to "
      f"{100*abs(Ed[-1]-(run['tr'].ke[0]-run['tr'].ke[-1]))/Ed[-1]:.1f} %)")

fig, ax = ST.figure_one_col()
ax.plot(x, Ed, "o-", ms=3, lw=0.8)
ax.set_xlabel("projectile displacement  x  (Bohr)")
ax.set_ylabel(r"$\Delta E_\mathrm{total}$  (Ha)")
ax.set_title(r"The raw signal $\Delta E_\mathrm{total}(x)$ — note the upward curvature")
fig
''')

md(r"""### 2.1 — The central diagnostic: the slope curve $S(x_0)$

If a steady state existed, $S(x_0)$ would **fall through the transient then go
flat** (a plateau). What we actually see decides the comparison.""")

code(r'''ag = detect_x0_and_stopping_power(x, Ed)
ff = fixed_fraction(x, Ed, frac=0.20)
grid, Sg, seg = ag["grid"], ag["Sgrid"], ag["segrid"]
frac_grid = (grid - x.min()) / L
SLR = LE.stopping_power_point(3.0, KF)

fig, ax = ST.figure_one_col()
ax.errorbar(frac_grid, Sg, yerr=seg, fmt="o-", ms=3, lw=0.7, capsize=1.5,
            label=r"$S(x_0)$ free-intercept fit")
ax.axhline(SLR, ls="--", color="k", lw=1.0, label=f"Lindhard (point) = {SLR:.4f}")
if "Spl" in ag:
    ax.axhspan(ag["Spl"]-ag["tol"], ag["Spl"]+ag["tol"], color="green", alpha=0.12,
               label=f"agent tolerance band (±{ag['tol']:.4f})")
ax.axvline(0.20, ls=":", color="C3", lw=1.2, label="M1 fixed 20% cut")
ax.axvline(0.40, ls="-.", color="gray", lw=1.0, label="40% hard gate")
if ag.get("x0") is not None:
    ax.axvline(ag["transient_fraction"], color="C2", lw=1.2,
               label=f"M2 detected x0 ({ag['transient_fraction']:.2f})")
ax.set_xlabel(r"transient fraction  $(x_0-x_\mathrm{min})/L$")
ax.set_ylabel("fitted slope  S  (Ha/Bohr)")
ax.set_title(r"$S(x_0)$: does the slope ever stop moving?")
ax.legend(fontsize=5.5, loc="upper left")
print("M2 agent status :", ag["status"], " x0:", ag.get("x0"),
      " S:", round(ag["S"],4) if ag.get("S") else None)
print(f"M1 fixed-20%     : S = {ff['S']:.4f} +/- {ff['se']:.4f} Ha/Bohr  (r2={ff['r2']:.4f})")
print(f"Lindhard point   : S = {SLR:.4f}   ->  M1/LR = {ff['S']/SLR:.2f}")
fig
''')

md(r"""**Read this plot.** If $S(x_0)$ is flat over a range, that range is the
plateau and M2 reports its left edge. If $S(x_0)$ keeps **drifting** (no flat
region within the green tolerance band), M2 returns `no_plateau` and **refuses**
to report $S$ — by design. M1 instead reads the slope at the dotted 20 % line
regardless. The gap between the dotted (M1) line's height and the dashed Lindhard
line is M1's error against the reference.""")

md(r"""### 2.2 — Both fits drawn on the signal""")

code(r'''fig, ax = ST.figure_one_col()
ax.plot(x, Ed, "o", ms=3, color="0.4", label="data")
xs = np.linspace(x.min(), x.max(), 100)
ax.plot(xs, ff["S"]*xs + ff["E0"], "-", color="C3", lw=1.4,
        label=f"M1 fixed-20%:  S={ff['S']:.4f}")
ax.axvline(x.min()+0.20*L, ls=":", color="C3", lw=1.0)
if ag.get("x0") is not None:
    ax.plot(xs, ag["S"]*xs + ag["E0"], "-", color="C2", lw=1.4,
            label=f"M2 agent:  S={ag['S']:.4f}")
    ax.axvline(ag["x0"], ls=":", color="C2", lw=1.0)
else:
    ax.plot([], [], " ", label=f"M2 agent: {ag['status']} (no S reported)")
ax.set_xlabel("x  (Bohr)"); ax.set_ylabel(r"$\Delta E_\mathrm{total}$  (Ha)")
ax.set_title("Free-intercept fits from each method")
ax.legend(fontsize=6, loc="upper left")
fig
''')

md(r"""### 2.3 — Sanity channels (run every time; large deviations get flagged)

Three independent-ish readouts of the deposited energy on the **same** fit window:
the primary $\Delta E_{\rm total}$, the projectile kinetic loss $\Delta$KE, and the
cumulative force-power $\int(-F\!\cdot\!v)\,dt$.

> **Honest caveat.** The track stores only $x,v$ (no force column), so
> $F=m\,dv/dt$ makes $\int(-F\!\cdot\!v)\,dt\equiv\Delta$KE *analytically* — it is a
> deposition-**profile** / discretisation check, **not** a third independent
> physics channel. The genuinely independent comparison is $\Delta E_{\rm total}$
> (electronic) vs $\Delta$KE (mechanical) — their agreement is energy
> conservation.""")

code(r'''x0_use = ff["x0"]; xT_use = x.max()
sc = sanity_channels(run, x0_use, xT_use)
fig, ax = ST.figure_one_col()
ax.plot(x, Ed, "o-", ms=3, lw=0.7, label=r"$\Delta E_\mathrm{total}$ (primary)")
ax.plot(sc["x_obs"], sc["dKE_obs"], "s--", ms=2.5, lw=0.7, label=r"$\Delta$KE loss (sanity a)")
ax.plot(sc["x_obs"], sc["dF_obs"], "^:", ms=2.5, lw=0.7, label=r"$\int(-F\cdot v)dt$ (sanity b)")
ax.set_xlabel("x  (Bohr)"); ax.set_ylabel("deposited energy  (Ha)")
ax.set_title("Three energy channels (should overlay if energy is conserved)")
ax.legend(fontsize=6, loc="upper left")
S_prim = ff["S"]; S_ke = sc["fKE"]["S"]; S_f = sc["fF"]["S"]
print(f"slope on the 20% window:")
print(f"  primary  dE_total : S = {S_prim:.4f} Ha/Bohr")
print(f"  sanity-a dKE      : S = {S_ke:.4f}   (dev {100*(S_ke-S_prim)/S_prim:+.1f} %)")
print(f"  sanity-b int F.v  : S = {S_f:.4f}   (dev {100*(S_f-S_prim)/S_prim:+.1f} %)")
DEV_FLAG = 10.0
for nm, Sx in [("dKE", S_ke), ("intFv", S_f)]:
    if abs(100*(Sx-S_prim)/S_prim) > DEV_FLAG:
        print(f"  ** FLAG: {nm} deviates >{DEV_FLAG:.0f}% from primary -> investigate **")
fig
''')

# ----------------------------------------------------------------------------- 3 sweep
md(r"""## 3 — Full sweep: both methods on every $(\sigma_q, v_0)$ run

`M2 status` ∈ {`ok`, `no_plateau`, `range_too_short`, `endpoint_contaminated`}.
Where M2 is not `ok`, it reports **no** $S$ — that refusal is the comparison's
main signal. `M1/LR` and `M2/LR` are the ratios to point-charge Lindhard.""")

code(r'''rows = []
print(f"{'sigma_q':>7} {'v0':>4} {'np':>3} | {'M1_S':>7} {'+-':>6} {'M1/LR':>5} | "
      f"{'M2_S':>7} {'x0frac':>6} {'M2/LR':>5} {'status':>17} | {'S_LR':>7}")
for sigma, rd, layout in SIGMAS:
    for vtag, v0 in VTAGS:
        run = load_run(rd, vtag, layout)
        if run is None:
            continue
        x, Ed = run["x"], run["Ed"]; L = x.max()-x.min()
        ff = fixed_fraction(x, Ed, 0.20)
        ag = detect_x0_and_stopping_power(x, Ed)
        slr = LE.stopping_power_point(v0, KF)
        m2S = ag.get("S"); tf = ag.get("transient_fraction")
        rows.append(dict(sigma=sigma, v0=v0, M1=ff["S"], M1e=ff["se"],
                         M2=m2S, tf=tf, status=ag["status"], slr=slr))
        print(f"{sigma:>7.2f} {v0:>4.1f} {len(x):>3} | "
              f"{ff['S']:>7.4f} {ff['se']:>6.4f} {ff['S']/slr:>5.2f} | "
              f"{(f'{m2S:.4f}' if m2S else '   --'):>7} "
              f"{(f'{tf:.2f}' if tf is not None else '--'):>6} "
              f"{(f'{m2S/slr:.2f}' if m2S else '--'):>5} {ag['status']:>17} | {slr:>7.4f}")
n_ok = sum(1 for r in rows if r['status']=='ok')
print(f"\nM2 reported S on {n_ok}/{len(rows)} runs; M1 reported S on all {len(rows)}.")
''')

# ----------------------------------------------------------------------------- 4 S(v)
md(r"""## 4 — $S(v)$ vs Lindhard — the evaluation plot

The test you asked for: each method's $S(v_0)$ over the point-charge Lindhard
reference, one panel. M1 points (filled) exist for every run; M2 points (ringed)
appear only where the detector returned `ok`.""")

code(r'''fig, ax = ST.figure_one_col()
vg = np.linspace(0.12, 3.25, 80)
slr_curve = np.array([LE.stopping_power_point(float(v), KF) for v in vg])
ax.plot(vg, slr_curve, "-", color="k", lw=1.6, label="Lindhard (point charge)")
cmap = ST.plt.get_cmap(ST.cmap_for("sequential"))
sig_list = [s for s, _, _ in SIGMAS]
for i, sigma in enumerate(sig_list):
    col = cmap(0.82 - 0.64 * i / max(len(sig_list)-1, 1))
    m1 = [(r["v0"], r["M1"], r["M1e"]) for r in rows if r["sigma"]==sigma]
    m2 = [(r["v0"], r["M2"]) for r in rows if r["sigma"]==sigma and r["status"]=="ok"]
    if m1:
        v_, s_, e_ = zip(*sorted(m1))
        ax.errorbar(v_, s_, yerr=e_, fmt="o", color=col, ms=4.5, capsize=2,
                    elinewidth=0.7, mec="k", mew=0.3, zorder=4)
    if m2:
        v2_, s2_ = zip(*sorted(m2))
        ax.plot(v2_, s2_, "o", mfc="none", mec=col, mew=1.6, ms=9, zorder=5)
    ax.plot([], [], "o", color=col, mec="k", mew=0.3, label=rf"$\sigma_q$={sigma}")
ax.plot([], [], "o", mfc="none", mec="0.3", mew=1.6, ms=9, label="M2 (ringed): detector ok")
ax.plot([], [], "o", color="0.3", label="M1 (filled): fixed-20%")
ax.axvline(KF, ls="--", color="gray", lw=0.8); ax.text(KF, ax.get_ylim()[1]*0.96, " k_F",
        fontsize=6, color="gray", va="top")
ax.set_xlabel("v  (a.u.)"); ax.set_ylabel("S(v)  (Ha/Bohr)"); ax.set_ylim(bottom=0)
ax.set_title(f"S(v): M1 vs M2 vs Lindhard  (r_s={RS}, k_F={KF:.3f})", fontsize=7.5)
ax.legend(fontsize=5.5, loc="upper right")
fig
''')

# ----------------------------------------------------------------------------- 5 findings
md(r"""## 5 — What the comparison shows (neutral — the verdict is yours)

- **M1 (fixed-20 %)** always returns a number; on the fast, clean runs it lands
  within a few % of point-charge Lindhard.
- **M2 (slope-plateau agent)** returns $S$ **only** when $S(x_0)$ genuinely
  flattens within its 2 % tolerance; otherwise it flags `no_plateau` /
  `range_too_short` and reports nothing.
- The deciding physical fact is in §2.1: $\Delta E_{\rm total}(x)$ is **convex**
  (slope rises monotonically with the cut), because the wake-formation time
  ($\sim1/\omega_p\approx8$ a.u.) is comparable to / longer than these runs
  (6–20 a.u.) — steady state is **approached but not reached**. M2 is strict
  about that; M1 reports the bulk-average slope anyway.

**Your decision:** is M1's pragmatic bulk-average the right operating choice for
this project, or is M2's refusal-on-no-plateau the safer gate (at the cost of
discarding most runs / demanding longer simulations)? Once you decide, we encode
the winner — or a hybrid (M1 default + M2 as a flag) — into the
`stopping-power-extraction` skill.

*No verdict is recorded here by design (verification-user-owns-verdict).*""")

# ----------------------------------------------------------------------------- execute
nb = new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"name": "inqview-venv", "display_name": "inqview-venv",
                             "language": "python"}
from nbconvert.preprocessors import ExecutePreprocessor
ep = ExecutePreprocessor(timeout=600, kernel_name="inqview-venv")
print("executing notebook ...")
ep.preprocess(nb, {"metadata": {"path": HERE}})
with open(OUT, "w") as fh:
    nbf.write(nb, fh)
print("wrote", OUT)
