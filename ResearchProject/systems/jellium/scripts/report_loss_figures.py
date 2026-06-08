#!/usr/bin/env python3
"""Report-standard INDIVIDUAL panels for the loss-function family (report1 style).

Remakes m8 / m9 / per-run-stopping / projectile_loss as publication panels:
apply_style() (usetex), NO titles/suptitles, one PNG per panel, 600 DPI,
inferno for L(q,ω) magnitude maps, report palettes + references for line plots.

Panels written to .../batch2_figures/report_standard/:
  fig_loss_map_<run>.png     L(q,ω) map (inferno LogNorm) + plasmon (cyan) + wrap (green)
  fig_loss_cut_m{1..4}.png   per-mode ω·L overlay across runs (E3.4/E15[/E25])
  fig_stopping_sv.png        S(v): analytic + loss-fn + classical + WP drift/total
  fig_stopping_ratio.png     classical-limit ratio
  fig_stopping_perrun_norm.png  per-run loss-fn S(v) shape consistency (unit-area)
"""
from __future__ import annotations
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from inqview.postprocess import lindhard
from inqview.report1 import apply_style, panel_label, palette_sweep3, palette_regime3, references
from inqview.report1._shared_style import ONE_COL_IN, STYLE_CONFIG, fix_one_col_axes

apply_style()
HA = 27.211386245988
ROOT = Path("/local/data/public/skcb2/tddft")
JB = ROOT / "ResearchProject/systems/jellium"
OUT = ROOT / "docs/presentations/storyline/tasks/batch2_figures/report_standard"
OUT.mkdir(parents=True, exist_ok=True)
DPI = STYLE_CONFIG["save_dpi"]
L_BOHR, N = 50.0, 162
N_DENS = N / L_BOHR**3
kF = (3 * np.pi**2 * N_DENS) ** (1 / 3); omega_p = np.sqrt(4 * np.pi * N_DENS); vF = kF
q1 = 2 * np.pi / L_BOHR; OMEGA_CAP = 16.0 / HA


def bohm_gross(q):
    return np.sqrt(omega_p**2 + 0.6 * vF**2 * q**2 + 0.25 * q**4)


def loss_qw(csv):
    df = pd.read_csv(csv); modes = sorted(df["m"].unique())
    qs, om0, Ls = [], None, []
    for m in modes:
        sub = df[df["m"] == m].sort_values("time_au"); t = sub["time_au"].values
        nq = sub["re_n_q"].values + 1j * sub["im_n_q"].values
        nq = nq - nq.mean(); Nn = len(t); q = sub["q_au"].values[0]
        fft = np.fft.fft(nq * np.hanning(Nn)); fr = np.fft.fftfreq(Nn, d=t[1] - t[0]); pos = fr >= 0
        om = fr[pos] * 2 * np.pi; L = (np.abs(fft[pos]) ** 2) / q**2
        if om0 is None: om0 = om
        qs.append(q); Ls.append(L)
    return np.array(qs), om0, np.vstack(Ls)


def S_from_L(v, q, om, Lmat):
    g = np.zeros_like(q)
    for i, qi in enumerate(q):
        wmax = min(qi * v, OMEGA_CAP); sel = (om > 0) & (om <= wmax)
        if sel.sum() < 2: continue
        g[i] = np.trapezoid(om[sel] * Lmat[i, sel], om[sel]) / qi
    return (2.0 / (np.pi * v**2)) * np.trapezoid(g, q)


RUNS = [("E3p4", 0.500, JB / "run_plasmon_n162_L50_E3p4_varyv"),
        ("E15", 1.050, JB / "run_plasmon_n162_L50_E15")]
e25 = JB / "run_plasmon_n162_L50_E25"
if (e25 / "results/analysis/observables/n_q_vs_time.csv").exists():
    RUNS.append(("E25", 1.356, e25))
loss = {tag: (v, *loss_qw(r / "results/analysis/observables/n_q_vs_time.csv")) for tag, v, r in RUNS}

# ---- loss map per run (individual) ----
for tag, (v, q, om, L) in loss.items():
    fig, ax = plt.subplots(figsize=ONE_COL_IN)
    sel = (om * HA <= 9) & (om > 0)
    Q, W = np.meshgrid(q, om[sel] * HA, indexing="ij"); Lp = L[:, sel]
    vmax = np.percentile(Lp[Lp > 0], 99.5)
    pc = ax.pcolormesh(Q, W, Lp, shading="auto", cmap="inferno", norm=LogNorm(vmax * 1e-4, vmax))
    qq = np.linspace(q.min(), q.max(), 80)
    ax.plot(qq, bohm_gross(qq) * HA, color="#00d0d0", lw=1.3)
    ax.plot(q, np.arange(1, len(q) + 1) * v * q1 * HA, color="#33cc33", ls="--", lw=1.2)
    ax.set_xlabel(r"$q$ (Bohr$^{-1}$)"); ax.set_ylabel(r"$\omega$ (eV)"); ax.set_ylim(0, 9)
    cb = fig.colorbar(pc, ax=ax, pad=0.02, fraction=0.046); cb.set_label(r"$L(q,\omega)$ (arb.)")
    fig.subplots_adjust(left=0.165, right=0.88, bottom=0.16, top=0.97)
    fig.savefig(OUT / f"fig_loss_map_{tag}.png", dpi=DPI); plt.close(fig)
    print(f"wrote fig_loss_map_{tag}.png")

# ---- per-mode cuts (overlay runs), individual ----
runtags = list(loss); nq = min(len(loss[t][1]) for t in runtags)
for mi in range(min(4, nq)):
    fig, ax = plt.subplots(figsize=ONE_COL_IN)
    for k, tag in enumerate(runtags):
        v, q, om, L = loss[tag]; s = (om * HA <= 9) & (om > 0)
        wl = om[s] * L[mi, s]; ax.plot(om[s] * HA, wl / max(np.trapezoid(wl, om[s] * HA), 1e-30),
                                       color=palette_sweep3[k % 3], label=tag.replace("E3p4", "E3.4"))
    qm = loss[runtags[0]][1][mi]
    ax.axvline(bohm_gross(qm) * HA, color="#00a0a0", lw=1.0)
    ax.axvline(omega_p * HA, **references["theory"])
    ax.set_xlabel(r"$\omega$ (eV)"); ax.set_ylabel(r"$\omega\,L(\omega)$ (norm.)"); ax.set_xlim(0, 9)
    ax.legend(loc="upper right"); panel_label(ax, f"$m={mi+1}$")
    fix_one_col_axes(fig)
    fig.savefig(OUT / f"fig_loss_cut_m{mi+1}.png", dpi=DPI); plt.close(fig)
    print(f"wrote fig_loss_cut_m{mi+1}.png")

# ---- stopping S(v) + ratio (reuse m9 logic) ----
from inqview.report1 import apply_style as _a  # noqa
try:
    from inqview.report1 import stopping_power_data as spd
except Exception:
    import importlib; spd = importlib.import_module("inqview.report1.stopping_power_data")

qE, omE, LE = loss["E15"][1], loss["E15"][2], loss["E15"][3]
v_grid = np.linspace(0.2, 6.8, 80)
S_LF = np.array([S_from_L(v, qE, omE, LE) for v in v_grid])
S_box = np.array([lindhard.stopping_power(v, kF, qmin=qE[0], qmax=qE[-1]) for v in v_grid]) * HA
S_full = np.array([lindhard.stopping_power(v, kF) for v in v_grid]) * HA
m = S_box > 0
S_LF *= np.trapezoid(S_box[m], v_grid[m]) / max(np.trapezoid(S_LF[m], v_grid[m]), 1e-30)
cls_v, cls_S = [], []
for sp in spd.get_L50_classical_runs():
    try:
        p = ROOT / sp.run_dir; pr = spd.parse_run_summary(p); win = spd.compute_time_window(pr)
        S, _ = spd.compute_classical_S(p, win); vv = np.sqrt(2 * sp.energy_eV / HA)
        if np.isfinite(S): cls_v.append(vv); cls_S.append(abs(S))
    except Exception: pass
cls_v, cls_S = np.array(cls_v), np.array(cls_S)
wp_v, wp_Sd = [], []
for sp in spd.get_L50_wp_sigma1_runs():
    rd = ROOT / sp.run_dir if not str(sp.run_dir).startswith("/") else Path(sp.run_dir)
    pr = spd.parse_run_summary(rd); win = spd.compute_time_window(pr)
    csv = rd / "results/raw/observables/wp_momentum_stats.csv"
    if not csv.exists(): continue
    df = pd.read_csv(csv, comment="#"); df = df[df["time_au"] <= win.t_end]
    if len(df) < 2: continue
    pz = df["pz_mean"].values; dz = pr.wp_k0_z * win.t_end
    if abs(dz) < 0.1: continue
    wp_v.append(abs(pr.wp_k0_z)); wp_Sd.append(-(pz[-1]**2/2 - pz[0]**2/2) * HA / dz)
wp_v, wp_Sd = np.array(wp_v), np.array(wp_Sd)

fig, ax = plt.subplots(figsize=ONE_COL_IN)
ax.plot(v_grid, S_full, color="#185070", lw=1.4, label="Lindhard")
ax.plot(v_grid, S_box, color="#185070", ls="--", lw=1.0, label="Lindhard (box $q$)")
ax.plot(v_grid, S_LF, color="#881818", lw=1.4, label="loss-fn $S(v)$")
if len(cls_v): ax.plot(cls_v, cls_S, "^", color="#000000", ms=5, label="classical")
if len(wp_v): ax.plot(wp_v, wp_Sd, "o", color="#188048", ms=5, label=r"WP $S_{\rm drift}$")
ax.set_xlabel(r"$v$ (a.u.)"); ax.set_ylabel(r"$S(v)$ (eV/Bohr)"); ax.set_xlim(0, 6.8); ax.set_ylim(bottom=0)
ax.legend(loc="upper right"); fix_one_col_axes(fig)
fig.savefig(OUT / "fig_stopping_sv.png", dpi=DPI); plt.close(fig); print("wrote fig_stopping_sv.png")

if len(cls_v):
    fig, ax = plt.subplots(figsize=ONE_COL_IN)
    o = np.argsort(cls_v); cvo, cso = cls_v[o], cls_S[o]
    ax.plot(cvo, np.interp(cvo, v_grid, S_LF) / cso, "d-", color="#881818", ms=5, label="loss-fn / cl.")
    if len(wp_v):
        wc = np.interp(wp_v, cvo, cso); g = wc > 0
        ax.plot(wp_v[g], wp_Sd[g] / wc[g], "o-", color="#188048", ms=5, label=r"WP $S_{\rm drift}$ / cl.")
    ax.axhline(1.0, **references["asymptote"])
    ax.set_xlabel(r"$v$ (a.u.)"); ax.set_ylabel(r"$S/S_{\rm classical}$"); ax.set_ylim(bottom=0)
    ax.legend(loc="upper right"); fix_one_col_axes(fig)
    fig.savefig(OUT / "fig_stopping_ratio.png", dpi=DPI); plt.close(fig); print("wrote fig_stopping_ratio.png")

# ---- per-run S(v) shape consistency (unit-area) ----
fig, ax = plt.subplots(figsize=ONE_COL_IN)
for k, tag in enumerate(runtags):
    v, q, om, L = loss[tag]; S = np.array([S_from_L(vv, q, om, L) for vv in v_grid])
    ax.plot(v_grid, S / max(np.trapezoid(S, v_grid), 1e-30), color=palette_sweep3[k % 3],
            label=tag.replace("E3p4", "E3.4"))
ax.plot(v_grid, S_full / np.trapezoid(S_full, v_grid), **references["theory"], label="Lindhard")
ax.set_xlabel(r"$v$ (a.u.)"); ax.set_ylabel(r"$S(v)$ (unit area)"); ax.set_xlim(0, 5); ax.set_ylim(bottom=0)
ax.legend(loc="upper right"); fix_one_col_axes(fig)
fig.savefig(OUT / "fig_stopping_perrun_norm.png", dpi=DPI); plt.close(fig)
print("wrote fig_stopping_perrun_norm.png")
