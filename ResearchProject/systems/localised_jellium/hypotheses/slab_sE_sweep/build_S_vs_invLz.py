#!/usr/bin/env python3
"""S vs 1/L_z for ONE projectile energy — classical & wavepacket slab traces with
a quadratic (empirical) extrapolation to the bulk (L_z→∞) limit.

Usage:  build_S_vs_invLz.py <E_eV>

Reads the campaign's sweep_results.csv (rows for this energy), plots S=E_abs/L_z
(slab-stopping-power skill) vs 1/L_z with y-error bars = |E_abs(0.9 t_f) −
E_abs(t_f)|/L_z (late-time drift), fits S = a + b/L_z + c/L_z², reads off the
1/L_z=0 intercept for classical and WP, and APPENDS that bulk point to
sE_results.csv (consumed by build_S_vs_E.py). No separate periodic-bulk run
(user 2026-09-11: bulk limit from extrapolation only).

Output: hypotheses/slab_sE_sweep/S_vs_invLz_E<E_eV>.png
"""
import csv, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import inqview.visualisation.style as st

st.apply_theme()
HA = 27.211386245988
LJ = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium")
BASE = LJ / "scripts/slab_sE_sweep"
RUNS = BASE / "runs"
RESULTS = BASE / "sweep_results.csv"
SE_CSV = BASE / "sE_results.csv"
HYP = Path(__file__).resolve().parent

sys.path.insert(0, "/home/raid/skcb2/skcb2/tddft/.claude/skills/slab-stopping-power")
import slab_stopping as SS  # noqa: E402

EV = int(sys.argv[1])
OUT = HYP / f"S_vs_invLz_E{EV}.png"


def _safe(fn, *a):
    """Return fn(*a) or 0.0 if it raises (a missing drift error just drops the bar)."""
    try:
        return fn(*a)
    except Exception as ex:
        print(f"  y-error {a}: {ex}")
        return 0.0

df = pd.read_csv(RESULTS)
df = df[df["E_eV"] == EV].sort_values("L_z")
if df.empty:
    print(f"no rows for E={EV} eV in {RESULTS}"); sys.exit(1)
v0 = float(df["v0"].iloc[0])
lz = df["L_z"].to_numpy(float)
S_cl = df["S_classical"].to_numpy(float)
S_wp = df["S_wp"].to_numpy(float)
egs = {int(r["L_z"]): float(r["E_GS_ha"]) for _, r in df.iterrows()}

inv = 1.0 / lz
mcl, mwp = np.isfinite(S_cl), np.isfinite(S_wp)
inv_cl, S_cl_pts, lz_cl = inv[mcl], S_cl[mcl], lz[mcl]
inv_wp, S_wp_pts, lz_wp = inv[mwp], S_wp[mwp], lz[mwp]
oc, ow = np.argsort(inv_cl), np.argsort(inv_wp)
inv_cl, S_cl_pts, lz_cl = inv_cl[oc], S_cl_pts[oc], lz_cl[oc]
inv_wp, S_wp_pts, lz_wp = inv_wp[ow], S_wp_pts[ow], lz_wp[ow]

print(f"E={EV} eV (v0={v0}):")
for L, c, w in zip(lz, S_cl, S_wp):
    print(f"  L_z={L:.0f}: S_cl={c}  S_wp={w}")


def s_error(L, half):
    """|E_abs(0.9 t_f) − E_abs(t_f)|/L_z, using the skill's corrected total."""
    od = RUNS / f"E{EV}_Lz{int(L)}/results" / half / "raw/observables"
    e_gs = egs[int(L)]
    corr = SS.extensive_kinetic_total(od)
    if corr is not None:
        steps = corr["step"].to_numpy(); etot = corr["e_total_corr_ha"].to_numpy()
    else:
        o = SS.concat_obs(od, "observables")
        steps = o["step"].to_numpy(); etot = o["energy_total"].to_numpy()
    ix = SS.concat_obs(od, "interactions") if half == "classical" else None

    def e_abs_at(ts):
        i = int(np.abs(steps - ts).argmin()); e = etot[i]
        if ix is not None:
            j = int(np.abs(ix["step"].to_numpy() - steps[i]).argmin())
            e = e + float(ix["e_pb"].to_numpy()[j])
        return (e - e_gs) * HA

    sf = steps.max()
    return abs(e_abs_at(sf) - e_abs_at(0.9 * sf)) / L


yerr_cl = np.array([_safe(s_error, L, "classical") for L in lz_cl]) if len(lz_cl) else np.array([])
yerr_wp = np.array([_safe(s_error, L, "wp") for L in lz_wp]) if len(lz_wp) else np.array([])


def bulk_extrap(inv_x, S_y, ax, color):
    """Quadratic in 1/L_z if ≥3 pts, else linear; return the 1/L_z=0 intercept."""
    deg = 2 if len(inv_x) >= 3 else 1
    c = np.polyfit(inv_x, S_y, deg)
    b0 = float(np.polyval(c, 0.0))
    xs = np.linspace(0.0, inv_x.min(), 40)
    ax.plot(xs, np.polyval(c, xs), ls=":", color=color, lw=1.6, zorder=2)
    ax.plot(0.0, b0, marker="D", mfc="white", mec=color, mew=1.6, ms=9, zorder=5)
    ax.annotate(f"{b0:.2f}", (0.0, b0), textcoords="offset points", xytext=(7, 0),
                va="center", fontsize=12, color=color)
    return b0


C_CL, C_WP = st.CLASSICAL, st.WP
fig, ax = plt.subplots(figsize=(5.2, 4.4), layout="constrained")

common = sorted(set(lz_cl) & set(lz_wp))
if common:
    cdict, wdict = dict(zip(lz_cl, S_cl_pts)), dict(zip(lz_wp, S_wp_pts))
    inv_sh = np.array([1.0 / L for L in common])
    Sc_sh = np.array([cdict[L] for L in common]); Sw_sh = np.array([wdict[L] for L in common])
    sh = np.argsort(inv_sh)
    ax.fill_between(inv_sh[sh], Sc_sh[sh], Sw_sh[sh], color="0.5", alpha=0.18, zorder=1)

if len(lz_cl):
    ax.errorbar(inv_cl, S_cl_pts, yerr=yerr_cl, fmt="-s", color=C_CL, lw=1.6, ms=7,
                capsize=5, capthick=1.6, elinewidth=1.6, ecolor=C_CL, zorder=5)
if len(lz_wp):
    ax.errorbar(inv_wp, S_wp_pts, yerr=yerr_wp, fmt="-o", color=C_WP, lw=1.6, ms=7,
                capsize=5, capthick=1.6, elinewidth=1.6, ecolor=C_WP, zorder=5)

b_cl = bulk_extrap(inv_cl, S_cl_pts, ax, C_CL) if len(inv_cl) >= 2 else float("nan")
b_wp = bulk_extrap(inv_wp, S_wp_pts, ax, C_WP) if len(inv_wp) >= 2 else float("nan")
print(f"bulk extrapolation: classical={b_cl:.4f}  wavepacket={b_wp:.4f}  "
      f"(quantum {b_wp-b_cl:+.4f} eV/Bohr)")

ax.set_xlabel(r"$1/L_z$  (Bohr$^{-1}$)", fontsize=14)
ax.set_ylabel(r"$S$  (eV/Bohr)", fontsize=14)
ax.set_title(f"E = {EV} eV  ($v_0$ = {v0} a.u.)", fontsize=13)
ax.axvline(0.0, color="k", lw=2.0, ls="--", zorder=0)
ax.tick_params(direction="in", which="both", top=True, right=True, length=4)
handles = [
    Line2D([], [], color=C_CL, marker="s", ls="-", ms=7, label="classical"),
    Line2D([], [], color=C_WP, marker="o", ls="-", ms=7, label="wavepacket"),
    Patch(facecolor="0.5", alpha=0.18, label="quantum correction"),
    Line2D([], [], color="0.4", marker="D", mfc="white", mew=1.4, ls=":", ms=8,
           label="bulk extrapolation"),
]
ax.legend(handles=handles, loc="best", frameon=True, fancybox=False,
          edgecolor="0.7", framealpha=1.0, fontsize=11)
fig.savefig(OUT, dpi=200)
print(f"wrote {OUT}")

# ---- append the bulk point to sE_results.csv (upsert this energy) ----
cols = ["E_eV", "v0", "S_cl_bulk", "S_wp_bulk", "quantum_correction",
        "n_pts_cl", "n_pts_wp"]
rows = []
if SE_CSV.exists():
    rows = [r for r in csv.DictReader(open(SE_CSV)) if int(r["E_eV"]) != EV]
rows.append({"E_eV": EV, "v0": v0, "S_cl_bulk": f"{b_cl:.6f}", "S_wp_bulk": f"{b_wp:.6f}",
             "quantum_correction": f"{b_wp-b_cl:.6f}", "n_pts_cl": len(inv_cl),
             "n_pts_wp": len(inv_wp)})
rows.sort(key=lambda r: float(r["E_eV"]))
with open(SE_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
print(f"updated {SE_CSV}")
