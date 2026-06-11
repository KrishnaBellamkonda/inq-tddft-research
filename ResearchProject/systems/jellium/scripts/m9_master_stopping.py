#!/usr/bin/env python3
"""M9 + M10 — master stopping-power comparison: loss-function S(v) vs classical
vs wave-packet, in ONE figure, with a classical-limit ratio panel.

M9 ask: compute the loss-function stopping power and visualise it against the
other results (classical and wave-packets) so it is clear whether L-based S(v)
is a usable metric.  M10 ask: is classical an upper bound / do we recover the
classical stopping in the classical limit (high v)?

Curves / points (all eV/Bohr):
  * loss-function S(v)  : medium L(q,ω) from E15 long run (draft5 method),
                          area-normalised to the box-q analytic Lindhard curve
                          (arb units — power spectrum, NOT Im[-1/ε]; absolute
                          scale needs FDT calibration).  ONE curve, all v.
  * analytic Lindhard   : full kinematic q, and box-matched q window (absolute).
  * classical Ehrenfest : L50 classical runs (absolute, momentum/force method).
  * WP S_drift          : −Δ(⟨p_z⟩²/2)/Δz   (genuine drift deceleration).
  * WP S_ztot           : −Δ(⟨p_z²⟩/2)/Δz   (drift + momentum broadening).

Panel B: ratio WP/classical and (loss-fn)/classical vs v — the classical-limit
diagnostic (does the quantum stopping approach classical as v grows?).

Known-case (printed): loss-fn m=1 peak near ω_p; classical & WP points finite.
Output: batch2_figures/m9_master_stopping.png
"""
from __future__ import annotations
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt
from inqview.pipeline import lindhard
from inqview.report1 import stopping_power_data as spd

HA = 27.211386245988
ROOT = Path("/local/data/public/skcb2/tddft")
JB = ROOT / "ResearchProject/systems/jellium"
OUT = ROOT / "docs/presentations/storyline/tasks/batch2_figures"
NQ = JB / "run_plasmon_n162_L50_E15/results/analysis/observables/n_q_vs_time.csv"
L_BOHR, N = 50.0, 162
N_DENS = N / L_BOHR**3
kF = (3 * np.pi**2 * N_DENS) ** (1 / 3)
omega_p = np.sqrt(4 * np.pi * N_DENS); vF = kF
OMEGA_CAP = 16.0 / HA
v0_E15 = 1.04999


# ---- loss-function L(q,ω) and S(v) (medium property) ----
def loss_qw(csv):
    df = pd.read_csv(csv); modes = sorted(df["m"].unique())
    qs, om0, Ls = [], None, []
    for m in modes:
        sub = df[df["m"] == m].sort_values("time_au")
        t = sub["time_au"].values
        nq = sub["re_n_q"].values + 1j * sub["im_n_q"].values
        nq = nq - nq.mean(); Nn = len(t); q = sub["q_au"].values[0]
        fft = np.fft.fft(nq * np.hanning(Nn))
        fr = np.fft.fftfreq(Nn, d=t[1] - t[0]); pos = fr >= 0
        om = fr[pos] * 2 * np.pi; L = (np.abs(fft[pos]) ** 2) / q**2
        if om0 is None: om0 = om
        qs.append(q); Ls.append(L)
    return np.array(qs), om0, np.vstack(Ls)


def S_from_L(v, q, om, Lmat):
    g = np.zeros_like(q)
    for i, qi in enumerate(q):
        wmax = min(qi * v, OMEGA_CAP); sel = (om > 0) & (om <= wmax)
        if sel.sum() < 2: continue
        g[i] = np.trapz(om[sel] * Lmat[i, sel], om[sel]) / qi
    return (2.0 / (np.pi * v**2)) * np.trapz(g, q)


q_m, om_au, Lqw = loss_qw(NQ)
ipk = int(np.argmax(Lqw[0, 1:]) + 1)
print(f"[KC] loss-fn m=1 peak ω={om_au[ipk]*HA:.2f} eV (ω_p={omega_p*HA:.2f})")
v_grid = np.linspace(0.2, 6.8, 80)
S_LF = np.array([S_from_L(v, q_m, om_au, Lqw) for v in v_grid])
S_box = np.array([lindhard.stopping_power(v, kF, qmin=q_m[0], qmax=q_m[-1]) for v in v_grid]) * HA
S_full = np.array([lindhard.stopping_power(v, kF) for v in v_grid]) * HA
mask = S_box > 0
scale = np.trapz(S_box[mask], v_grid[mask]) / max(np.trapz(S_LF[mask], v_grid[mask]), 1e-30)
S_LF *= scale


# ---- classical points ----
cls_v, cls_S = [], []
for sp in spd.get_L50_classical_runs():
    try:
        p = ROOT / sp.run_dir; pr = spd.parse_run_summary(p); win = spd.compute_time_window(pr)
        S, _ = spd.compute_classical_S(p, win); v = np.sqrt(2 * sp.energy_eV / HA)
        if np.isfinite(S): cls_v.append(v); cls_S.append(abs(S))
    except Exception as e:
        print(f"  classical {sp.run_dir}: skip ({e})")
cls_v, cls_S = np.array(cls_v), np.array(cls_S)


# ---- WP drift / total points ----
def wp_metrics(run_dir):
    pr = spd.parse_run_summary(run_dir); win = spd.compute_time_window(pr)
    csv = run_dir / "results/raw/observables/wp_momentum_stats.csv"
    if not csv.exists(): return None
    df = pd.read_csv(csv, comment="#"); df = df[df["time_au"] <= win.t_end]
    if len(df) < 2: return None
    pz, pz2 = df["pz_mean"].values, df["pz2_mean"].values
    dz = pr.wp_k0_z * win.t_end
    if abs(dz) < 0.1: return None
    return (abs(pr.wp_k0_z),
            -(pz[-1]**2/2 - pz[0]**2/2) * HA / dz,
            -(pz2[-1]/2 - pz2[0]/2) * HA / dz)


wp_v, wp_Sd, wp_Sz = [], [], []
for sp in spd.get_L50_wp_sigma1_runs():
    rd = ROOT / sp.run_dir if not str(sp.run_dir).startswith("/") else Path(sp.run_dir)
    m = wp_metrics(rd)
    if m: wp_v.append(m[0]); wp_Sd.append(m[1]); wp_Sz.append(m[2])
wp_v, wp_Sd, wp_Sz = np.array(wp_v), np.array(wp_Sd), np.array(wp_Sz)
print(f"[KC] classical pts={len(cls_v)}  WP pts={len(wp_v)}")

# ---- figure ----
fig, (ax, axr) = plt.subplots(1, 2, figsize=(13.5, 5))
ax.plot(v_grid, S_full, "-", color="#1f77b4", lw=2, label="analytic Lindhard (full q)")
ax.plot(v_grid, S_box, "--", color="#1f77b4", lw=1.4, label=f"Lindhard (box q {q_m[0]:.2f}–{q_m[-1]:.2f})")
ax.plot(v_grid, S_LF, "-", color="#d62728", lw=2, label="loss-function S(v) (E15, arb→box-norm)")
if len(cls_v): ax.plot(cls_v, cls_S, "k^", ms=8, label="classical Ehrenfest")
if len(wp_v):
    ax.plot(wp_v, wp_Sd, "o", color="#2ca02c", ms=7, label=r"WP $S_{\rm drift}$ (σ=1)")
    ax.plot(wp_v, wp_Sz, "s", color="#9467bd", ms=6, label=r"WP $S_{z\rm tot}$ (σ=1)")
ax.axvline(v0_E15, color="green", ls=":", lw=1.0)
ax.set_xlabel("v (a.u.)"); ax.set_ylabel("S(v) (eV/Bohr)")
ax.set_title("M9 — loss-function S(v) vs analytic, classical, WP")
ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.set_xlim(0, 6.8); ax.set_ylim(bottom=0)

# ratio panel (classical-limit diagnostic)
if len(cls_v):
    order = np.argsort(cls_v); cvo, cso = cls_v[order], cls_S[order]
    S_LF_at_cls = np.interp(cvo, v_grid, S_LF)
    axr.plot(cvo, S_LF_at_cls / cso, "rd-", lw=1.6, ms=6, label="loss-fn / classical")
    if len(wp_v):
        wp_cls = np.interp(wp_v, cvo, cso)
        good = wp_cls > 0
        axr.plot(wp_v[good], wp_Sd[good] / wp_cls[good], "o-", color="#2ca02c", lw=1.6,
                 label=r"WP $S_{\rm drift}$ / classical")
    axr.axhline(1.0, color="0.4", ls="--", lw=1.0, label="classical-limit (ratio=1)")
    axr.set_xlabel("v (a.u.)"); axr.set_ylabel("S / S_classical")
    axr.set_title("M10 — approach to the classical limit")
    axr.legend(fontsize=8); axr.grid(alpha=0.3); axr.set_ylim(bottom=0)
fig.tight_layout()
fp = OUT / "m9_master_stopping.png"
fig.savefig(fp, dpi=150); plt.close(fig)
print(f"wrote {fp}")

# numeric summary
print("\n=== summary ===")
for e, v in [(20, 1.21), (100, 2.71), (300, 4.70)]:
    print(f"  v={v:.2f}: loss-fn S={np.interp(v, v_grid, S_LF):.3f} eV/Bohr")
if len(cls_v) and len(wp_v):
    print(f"  WP S_drift / classical ranges {np.nanmin(wp_Sd/np.interp(wp_v,np.sort(cls_v),cls_S[np.argsort(cls_v)])):.3f}"
          f"–{np.nanmax(wp_Sd/np.interp(wp_v,np.sort(cls_v),cls_S[np.argsort(cls_v)])):.3f} (→1 = classical limit)")
