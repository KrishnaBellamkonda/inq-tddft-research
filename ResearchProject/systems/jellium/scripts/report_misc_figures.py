#!/usr/bin/env python3
"""Report-standard INDIVIDUAL panels for m5 (diffraction), m6 (r_s-metal),
classical-confidence. report1 style (apply_style, usetex, 600 DPI, no titles,
one PNG per panel). Output -> batch2_figures/report_standard/.
"""
from __future__ import annotations
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import RegularGridInterpolator
from applications.report1 import apply_style, panel_label, palette_sweep3, references
from applications.report1._shared_style import ONE_COL_IN, STYLE_CONFIG, fix_one_col_axes

apply_style()
HA = 27.211386245988
ROOT = Path("/local/data/public/skcb2/tddft")
JB = ROOT / "ResearchProject/systems/jellium"
OUT = ROOT / "docs/presentations/storyline/tasks/batch2_figures/report_standard"
OUT.mkdir(parents=True, exist_ok=True)
DPI = STYLE_CONFIG["save_dpi"]

# ============================== M6: r_s vs metals ==============================
BOHR2ANG = 0.529177210903
USER = [("Li", 3.26, 1.727), ("Na", 3.99, 2.110), ("Cs", 5.75, 3.042)]
AM = {"Al": 2.07, "Mg": 2.66, "Li": 3.25, "Na": 3.93, "K": 4.86, "Rb": 5.20, "Cs": 5.62}
METALS = sorted([("Al", 2.07), ("Mg", 2.66), ("Li", 3.26), ("Na", 3.99), ("K", 4.86),
                 ("Rb", 5.20), ("Cs", 5.75)], key=lambda t: t[1])
fig, ax = plt.subplots(figsize=(STYLE_CONFIG.get("twocol", 7.0), 2.4))
for name, rs in METALS:
    ax.plot(rs, 0, "o", color="#404040", ms=6, zorder=3)
    ax.annotate(name, (rs, 0), xytext=(0, 8), textcoords="offset points", ha="center", fontsize=9)
    ax.annotate(f"{rs:.2f}", (rs, 0), xytext=(0, -14), textcoords="offset points", ha="center",
                fontsize=7, color="#666666")
for rs, lab in [(5.69, r"$r_s=5.69$"), (3.41, r"$r_s=3.41$")]:
    ax.plot(rs, 0, "*", color="#881818", ms=15, zorder=4)
    ax.annotate(lab, (rs, 0), xytext=(0, 26), textcoords="offset points", ha="center",
                fontsize=8.5, color="#881818",
                arrowprops=dict(arrowstyle="->", color="#881818", lw=1.0))
ax.axhline(0, color="#808080", lw=0.9); ax.set_xlim(1.8, 6.1); ax.set_ylim(-0.6, 1.0)
ax.set_yticks([]); ax.set_xlabel(r"$r_s$ (Bohr)")
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
fig.subplots_adjust(left=0.04, right=0.99, bottom=0.28, top=0.80)
fig.savefig(OUT / "fig_m6_rs_numberline.png", dpi=DPI); plt.close(fig); print("wrote fig_m6_rs_numberline.png")

fig, ax = plt.subplots(figsize=ONE_COL_IN); ax.axis("off")
rows = [[n, f"{rs:.2f}", f"{ang:.3f}", f"{AM[n]:.2f}", f"{rs-AM[n]:+.2f}"] for n, rs, ang in USER]
tbl = ax.table(cellText=rows, colLabels=["alkali", r"$r_s$ (a$_0$)", r"$r_s$ (\AA)", "A\\&M", r"$\Delta$"],
               loc="center", cellLoc="center")
tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1, 1.5)
for j in range(5): tbl[0, j].set_facecolor("#dddddd")
fig.subplots_adjust(left=0.02, right=0.98, bottom=0.05, top=0.95)
fig.savefig(OUT / "fig_m6_rs_table.png", dpi=DPI); plt.close(fig); print("wrote fig_m6_rs_table.png")

# ============================== M5: diffraction ===============================
d = np.load(JB.parent / "coronene/run_broadening_35x35x80/results/analysis/momentum/momentum_scatter_arrays.npz")
kz, kx, ky, dP3 = d["kz"], d["kx"], d["ky"], d["dP3"].astype(float); k0 = float(d["k0"])
A_G = 2.46 / BOHR2ANG; G1 = 4 * np.pi / (np.sqrt(3) * A_G)
G_ORD = [(r"$|G_1|$", G1), (r"$\sqrt{3}|G_1|$", np.sqrt(3) * G1), (r"$2|G_1|$", 2 * G1)]
interp = RegularGridInterpolator((kz, kx, ky), dP3, bounds_error=False, fill_value=0.0)
nth, nph = 90, 180
th = np.linspace(0, np.pi / 2, nth); ph = np.linspace(-np.pi, np.pi, nph)
TH, PH = np.meshgrid(th, ph, indexing="ij")
shell = interp(np.stack([(k0*np.cos(TH)).ravel(), (k0*np.sin(TH)*np.cos(PH)).ravel(),
                         (k0*np.sin(TH)*np.sin(PH)).ravel()], axis=1)).reshape(TH.shape)
# polar map
fig = plt.figure(figsize=(ONE_COL_IN[0], ONE_COL_IN[0])); ax = fig.add_subplot(111, projection="polar")
vmax = np.percentile(np.abs(shell), 99.5)
ax.pcolormesh(ph, np.degrees(th), shell, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
for _, G in G_ORD:
    if G < k0: ax.plot(ph, np.full_like(ph, np.degrees(np.arcsin(G / k0))), "k--", lw=0.7)
ax.set_rlabel_position(135)
fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
fig.savefig(OUT / "fig_m5_shell_polar.png", dpi=DPI); plt.close(fig); print("wrote fig_m5_shell_polar.png")
# I(phi), I(kperp)
solid = np.sin(TH); kperp = k0 * np.sin(th)
I_phi = (np.clip(shell, 0, None) * solid).sum(axis=0); I_th = (np.clip(shell, 0, None) * solid).sum(axis=1)
fig, ax = plt.subplots(figsize=ONE_COL_IN)
ax.plot(np.degrees(ph), I_phi, color=palette_sweep3[2])
for a in (45, 135, -45, -135): ax.axvline(a, color="#cc8800", ls=":", lw=0.8)
ax.set_xlabel(r"azimuth $\varphi$ (deg)"); ax.set_ylabel(r"$I(\varphi)$ (gain)")
panel_label(ax, "(a)"); fix_one_col_axes(fig)
fig.savefig(OUT / "fig_m5_I_phi.png", dpi=DPI); plt.close(fig); print("wrote fig_m5_I_phi.png")
fig, ax = plt.subplots(figsize=ONE_COL_IN)
ax.plot(kperp, I_th, color=palette_sweep3[0])
for lbl, G in G_ORD:
    if G < k0:
        ax.axvline(G, **references["theory"]); ax.text(G, ax.get_ylim()[1]*0.9, lbl, rotation=90, va="top", fontsize=7)
ax.set_xlabel(r"$k_\perp = k_0\sin\theta$ (Bohr$^{-1}$)"); ax.set_ylabel(r"$I(k_\perp)$ (gain)")
panel_label(ax, "(b)"); fix_one_col_axes(fig)
fig.savefig(OUT / "fig_m5_I_kperp.png", dpi=DPI); plt.close(fig); print("wrote fig_m5_I_kperp.png")

# ========================= classical confidence ==============================
def loss_qw(csv):
    df = pd.read_csv(csv); modes = sorted(df["m"].unique()); qs, om0, Ls = [], None, []
    for m in modes:
        sub = df[df["m"] == m].sort_values("time_au"); t = sub["time_au"].values
        nq = sub["re_n_q"].values + 1j * sub["im_n_q"].values; nq -= nq.mean()
        Nn = len(t); q = sub["q_au"].values[0]; fft = np.fft.fft(nq * np.hanning(Nn))
        fr = np.fft.fftfreq(Nn, d=t[1]-t[0]); pos = fr >= 0
        if om0 is None: om0 = fr[pos]*2*np.pi
        qs.append(q); Ls.append((np.abs(fft[pos])**2)/q**2)
    return np.array(qs), om0, np.vstack(Ls)
N_DENS = 162/50**3; omega_p = np.sqrt(4*np.pi*N_DENS); q1 = 2*np.pi/50
CL = JB / "run_classical_n162_L50_E100_v2"
qc, omc, Lc = loss_qw(CL / "results/analysis/observables/n_q_vs_time.csv")
# A: classical loss map
fig, ax = plt.subplots(figsize=ONE_COL_IN); sel = (omc*HA <= 12) & (omc > 0)
Q, W = np.meshgrid(qc, omc[sel]*HA, indexing="ij"); Lp = Lc[:, sel]; vmax = np.percentile(Lp[Lp>0], 99.5)
pc = ax.pcolormesh(Q, W, Lp, shading="auto", cmap="inferno", norm=LogNorm(vmax*1e-4, vmax))
ax.axhline(omega_p*HA, color="#00d0d0", ls=":", lw=1.0)
ax.set_xlabel(r"$q$ (Bohr$^{-1}$)"); ax.set_ylabel(r"$\omega$ (eV)"); ax.set_ylim(0, 12)
cb = fig.colorbar(pc, ax=ax, pad=0.02, fraction=0.046); cb.set_label(r"$L(q,\omega)$ (arb.)")
panel_label(ax, "(a)"); fig.subplots_adjust(left=0.165, right=0.88, bottom=0.16, top=0.97)
fig.savefig(OUT / "fig_classical_lossmap.png", dpi=DPI); plt.close(fig); print("wrote fig_classical_lossmap.png")
# C: projectile v_z
tr = pd.read_csv(CL / "results/raw/observables/electron_track.csv")
fig, ax = plt.subplots(figsize=ONE_COL_IN)
ax.plot(tr["time_au"], tr["vz"], color="#881818")
ax.set_xlabel(r"$t$ (a.u.)"); ax.set_ylabel(r"projectile $v_z$ (a.u.)")
panel_label(ax, "(b)"); fix_one_col_axes(fig)
fig.savefig(OUT / "fig_classical_vz.png", dpi=DPI); plt.close(fig); print("wrote fig_classical_vz.png")
# D: WP momentum dist before/after
WP = JB / "run_wp_n162_L50_E100_sigma1_v2"
md = pd.read_csv(WP / "results/raw/observables/momentum_distribution.csv", comment="#")
t0, tN = md["time_au"].min(), md["time_au"].max()
b = md[np.isclose(md["time_au"], t0)]; a = md[np.isclose(md["time_au"], tN)]
fig, ax = plt.subplots(figsize=ONE_COL_IN)
ax.plot(b["k_bohr_inv"], b["n_wp"], color=palette_sweep3[0], label=f"$t={t0:.0f}$")
ax.plot(a["k_bohr_inv"], a["n_wp"], color=palette_sweep3[2], label=f"$t={tN:.0f}$")
ax.axvline(np.sqrt(2*100/HA), **references["asymptote"])
ax.set_xlabel(r"$|k|$ (Bohr$^{-1}$)"); ax.set_ylabel(r"$n_{\rm wp}(|k|)$"); ax.set_xlim(0, 6)
ax.legend(loc="upper right"); panel_label(ax, "(c)"); fix_one_col_axes(fig)
fig.savefig(OUT / "fig_classical_wp_momentum.png", dpi=DPI); plt.close(fig); print("wrote fig_classical_wp_momentum.png")
