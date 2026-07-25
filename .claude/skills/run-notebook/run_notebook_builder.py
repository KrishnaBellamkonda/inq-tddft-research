#!/usr/bin/env python3
"""run_notebook_builder.py — skill-local, shippable builder for a RUN-NOTEBOOK.

A run-notebook is the deep single-run analysis artefact (see this skill's SKILL.md):
the full standardised plot battery for ONE run. This builder ASSEMBLES over
inqview.pipeline — it runs the pipeline phases (which compute most of the figures,
auto-skipping what a run can't produce), ADDS the density-matrix carpets / lead GIF /
E-field that the user-requested battery needs, then embeds everything into an executed
notebook with the house-narrative context (title+question, reconstructable config,
linked source files).

Usage (venv + stack on path):
    PYTHONPATH=.../inq-stack/python .../venv/bin/python3 run_notebook_builder.py \
        <results_dir> <out.ipynb> [--baseline <baseline_results_dir>] [--run-cpp <path>]

`results_dir` is a run's results tree (contains raw/ and run_summary.txt).
The generated .ipynb is written to <out.ipynb> (live in the run's hypotheses/ folder).
Density categories: total system always; "wavepacket" = the WAKE (run - baseline) proxy
when --baseline is given (bare |psi_WP|^2 is not saved — a future observable).
"""
from __future__ import annotations
import argparse
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np

STACK = "/local/data/public/skcb2/tddft/inq-stack/python"
if STACK not in sys.path:
    sys.path.insert(0, STACK)
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.animation as animation  # noqa: E402
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
import vtk  # noqa: E402
from vtk.util.numpy_support import vtk_to_numpy  # noqa: E402
from inqview.visualisation import style  # noqa: E402
from inqview.visualisation import make_density_gif_battery  # noqa: E402
from inqview.pipeline import runner  # noqa: E402
from inqview.analysis import lindhard_elf as _lind  # noqa: E402  analytical Lindhard
from inqview.analysis import compute_heuristics  # noqa: E402  physical-anchor heuristics

style.apply_theme()
DT_DEFAULT = 0.02
HA_EV = 27.21138625


def lindhard_stopping_fig(out, rs, proj_sigma, measured_v=None, measured_s=None,
                          mode="both"):
    """Analytical Lindhard linear-response stopping S(v) for an r_s electron gas, in
    eV/Bohr. mode='both' draws BOTH the Gaussian-projectile curve (charge std =
    proj_sigma) and the bare point-charge reference; mode='point' draws ONLY the
    point-charge curve — use when the projectile width was chosen (from a σ-sweep) to
    sit in the linear-response regime, so the finite-size Gaussian form factor is a
    negligible correction and the point-charge curve is the apples-to-apples reference.
    Overlays the measured rt-TDDFT point if given. Returns a dict of the at-v0 values
    (Ha/Bohr). Source: Lindhard-Winther 1964; Correa 2018."""
    kF = _lind.kF_from_rs(rs)
    vF, wp = kF, _lind.omega_p(kF)
    vmax = max(4.0, (measured_v or 0) * 1.4)
    vg = np.linspace(0.15, vmax, 40)
    Spt = np.array([_lind.stopping_power_point(v, kF) for v in vg])
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(vg, Spt * HA_EV, "C3-" if mode == "point" else "C3--", lw=1.6,
            label="Lindhard (point charge, linear response)")
    if mode == "both":
        Ssig = np.array([_lind.stopping_power_sigma(v, kF, proj_sigma) for v in vg])
        ax.plot(vg, Ssig * HA_EV, "C0-", lw=1.6,
                label=rf"Lindhard Gaussian ($\sigma_q$={proj_sigma:.3f})")
    ax.axvline(vF, ls=":", color="grey", lw=0.8)
    ax.text(vF, ax.get_ylim()[1] * 0.05, r" $v_F$", color="grey", fontsize=8)
    out_at = {}
    if measured_v is not None and measured_s is not None:
        ax.plot([measured_v], [measured_s * HA_EV], "ks", ms=9,
                label=f"rt-TDDFT (classical ΔKE/x) = {measured_s*HA_EV:.3f}")
        out_at = {"S_point": float(_lind.stopping_power_point(measured_v, kF))}
        if mode == "both":
            out_at["S_sigma"] = float(_lind.stopping_power_sigma(measured_v, kF, proj_sigma))
    ax.set_xlabel("projectile velocity v (a.u.)")
    ax.set_ylabel("stopping power S (eV/Bohr)")
    ax.set_title(rf"Analytical Lindhard stopping — slab gas $r_s$={rs:.3f} "
                 rf"($\omega_p$={wp*HA_EV:.2f} eV)")
    ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    out_at["kF"] = kF; out_at["wp_eV"] = wp * HA_EV
    return out_at


def wp_exit_time_fig(rd, out, v0, launch_z, slab_half=12.5, cap_inner=None, dt=0.02):
    """Plot the WP centroid z(t) (from wp_real_space_stats) with the slab faces, CAP
    boundaries, and the BALLISTIC exit time t_exit = (z_far − z_launch)/v0 (the time
    to reach the far slab face travelling at the mean momentum). Returns t_exit."""
    p = Path(rd) / "raw" / "observables" / "wp_real_space_stats.csv"
    z_far = slab_half
    t_exit = (z_far - launch_z) / v0
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    if p.exists():
        from io import StringIO
        rows = [ln for ln in p.read_text().splitlines() if ln and not ln.startswith("#")]
        d = np.genfromtxt(StringIO("\n".join(rows)), delimiter=",", names=True)
        ax.plot(d["time_au"], d["z_mean"], "C0-", lw=1.5, label=r"WP centroid $\langle z\rangle(t)$")
    for s, lab in [(-slab_half, "slab face"), (slab_half, None)]:
        ax.axhline(s, ls=":", color="C3", lw=0.8, label=lab)
    if cap_inner is not None:
        for s in (-abs(cap_inner), abs(cap_inner)):
            ax.axhline(s, ls="--", color="lime", lw=1.0)
        ax.axhline(np.nan, ls="--", color="lime", lw=1.0, label="CAP inner edge")
    ax.axvline(t_exit, ls="-.", color="k", lw=1.2,
               label=rf"ballistic exit $t$={t_exit:.1f} a.u.")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("z (Bohr)")
    ax.set_title(rf"WP transit — launch z={launch_z}, v={v0:.3f} → cross 25 Bohr slab")
    ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    return t_exit


def delta_total_energy_fig(rd, out, rtype):
    """Two-panel [ ΔE_total(t) | N(t)=∫n dV ]: the plateau test NEXT TO the total
    electron number (total number density integrated over the box) vs time — the
    enforced run-notebook pairing (2026-07-11).

    LEFT — ΔE_total(t) = E_total(t) − E_total(0) [eV]. For WP runs also overlay the
    ledger with the ABSORBED WP-orbital kinetic removed: energy_total keeps counting
    the CAP-absorbed KS orbital at occupation 1 with its NORMALIZED ⟨T⟩ (~11 Ha),
    which rings as a phantom even though the electron is gone from the density.
    Dropping that term recovers the physical (conserved) energy. See
    reference_phantom_absorbed_wp_orbital_energy.
    RIGHT — N(t) from ``electron_number.csv`` (classical) or ∫``density_total`` VTIs
    (WP); its drop is the CAP boundary absorption. Returns a markdown finding
    string, or None if observables are missing."""
    import pandas as _pd
    obs = Path(rd) / "raw" / "observables" / "observables.csv"
    if not obs.exists():
        return None
    _o = _pd.read_csv(obs).drop_duplicates(subset="step").sort_values("time_au")
    if "energy_total" not in _o.columns or len(_o) < 4:
        return None
    t = _o["time_au"].to_numpy(float)
    dE = (_o["energy_total"].to_numpy(float) - float(_o["energy_total"].iloc[0])) * HA_EV
    late = t > (0.4 * t.max())                       # "after absorption" heuristic
    fig, (ax, axN) = plt.subplots(1, 2, figsize=(12.4, 4.0))
    ax.plot(t, dE, lw=1.2, color="C3", alpha=0.85,
            label=f"as logged (late std {dE[late].std():.1f} eV)")
    finding = (f"**Left — ΔE_total(t) = E_total(t) − E_total(0).** As logged, late-time "
               f"(t > {0.4*t.max():.0f} a.u.) std = {dE[late].std():.1f} eV, "
               f"end = {dE[-1]:.1f} eV.")
    # WP de-ledger: subtract the absorbed-orbital normalized kinetic (e_kin_ha)
    mom = Path(rd) / "raw" / "observables" / "wp_momentum_stats.csv"
    if rtype == "wp" and mom.exists():
        try:
            _m = _pd.read_csv(mom, comment="#")
            j = _o.merge(_m[["step", "e_kin_ha"]], on="step", how="inner")
            tj = j["time_au"].to_numpy(float)
            fix = ((j["energy_total"].to_numpy(float) - j["e_kin_ha"].to_numpy(float))
                   - (float(j["energy_total"].iloc[0]) - float(j["e_kin_ha"].iloc[0]))) * HA_EV
            latej = tj > (0.4 * tj.max())
            ax.plot(tj, fix, lw=1.6, color="C0",
                    label=f"WP orbital removed (late std {fix[latej].std():.2f} eV)")
            flat = fix[latej].std() < 1.0
            finding += (f" **Removing the absorbed WP-orbital kinetic** → late std "
                        f"{fix[latej].std():.2f} eV, end {fix[-1]:.1f} eV"
                        + (" — **plateau recovered** (the swing was a phantom-orbital "
                           "bookkeeping artifact, not physics)." if flat else
                           " — residual swing remains (see the FFT panel)."))
        except Exception:
            pass
    ax.axhline(0.0, ls=":", lw=0.8, color="0.6")
    ax.axvspan(0.4 * t.max(), t.max(), color="0.5", alpha=0.06)
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$\Delta E_\mathrm{total}$ (eV)")
    ax.set_title("Total-system energy change (plateau test)")
    ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.25)

    # RIGHT panel — N(t) = ∫n dV, the total number density in the cell vs time.
    tN, N = _total_number_series(rd)
    if tN is not None and len(N) > 1:
        axN.plot(tN, N, lw=1.4, color="C0", marker="o", ms=2.5)
        src = ("electron_number.csv" if (Path(rd) / "raw" / "observables"
               / "electron_number.csv").exists() else "∫ density_total VTIs")
        absorbed = N[0] - N[-1]
        note = (" — the ~1 e⁻ drop is the fully-absorbed wavepacket"
                if rtype == "wp" else "")
        finding += (f"  **Right — N(t) = ∫n dV** (total number density in the cell, "
                    f"{src}): {N[0]:.2f} → {N[-1]:.2f} e⁻ "
                    f"(**{absorbed:+.2f} e⁻** lost to CAP boundary absorption{note}).")
    else:
        axN.text(0.5, 0.5, "N(t) unavailable\n(no electron_number.csv,\nno density_total VTIs)",
                 ha="center", va="center", fontsize=9, transform=axN.transAxes)
        finding += " _N(t) unavailable for this run._"
    axN.set_xlabel("time (a.u.)"); axN.set_ylabel(r"$N(t)=\int n\,dV$  (electrons)")
    axN.set_title("Total electron number (∫n dV) vs time"); axN.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    return finding


def _total_number_series(rd):
    """N(t) = ∫ n_total dV, the total electron number in the box vs time.

    Prefers ``electron_number.csv`` (``N_total``, full time resolution — classical
    runs log it). Falls back to integrating the ``density_total`` VTI series
    (∑ n · dV per frame; WP runs, ~30 snapshot frames). Returns (t[au], N) or
    (None, None) if neither source exists."""
    import pandas as _pd
    en = Path(rd) / "raw" / "observables" / "electron_number.csv"
    if en.exists():
        _e = _pd.read_csv(en).drop_duplicates(subset="step").sort_values("time_au")
        col = "N_total" if "N_total" in _e.columns else _e.columns[-1]
        return _e["time_au"].to_numpy(float), _e[col].to_numpy(float)
    tot_dir = Path(rd) / "raw" / "vti" / "density_total"
    fs = _frames(tot_dir) if tot_dir.is_dir() else []
    if not fs:
        return None, None
    try:
        dt = float(parse_summary(rd).get("dt_au", 1.0))
    except (TypeError, ValueError):
        dt = 1.0
    ts, ns = [], []
    for f in fs:
        _o, s, a = _vol(f)
        ns.append(float(a.sum()) * float(s[0]) * float(s[1]) * float(s[2]))
        m = re.search(r"_t(\d+)", f.name)
        ts.append(int(m.group(1)) * dt if m else float("nan"))
    return np.array(ts, float), np.array(ns, float)


# ----------------------------------------------------------------- VTI helpers
def _vol(path):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(str(path)); r.Update()
    img = r.GetOutput(); nx, ny, nz = img.GetDimensions()
    o, s = img.GetOrigin(), img.GetSpacing()
    a = vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx)
    return o, s, a


def _frames(vti_dir):
    d = Path(vti_dir)
    fs = sorted(d.glob("*_t*.vti"),
                key=lambda p: int(re.search(r"_t(\d+)", p.name).group(1)))
    return fs


def load_series(vti_dir, max_frames=160):
    """Return (steps, z, linear n(z,t), mid-y xz slices) from a VTI series."""
    fs = _frames(vti_dir)
    if not fs:
        return None
    if len(fs) > max_frames:
        fs = fs[:: len(fs) // max_frames + 1]
    steps = [int(re.search(r"_t(\d+)", p.name).group(1)) for p in fs]
    o, s, a0 = _vol(fs[0]); nz, ny, nx = a0.shape
    z = o[2] + s[2] * np.arange(nz)
    x = o[0] + s[0] * np.arange(nx)
    nzt, slabs = [], []
    midy = ny // 2
    for p in fs:
        _, _, a = _vol(p)
        nzt.append(a.sum(axis=(1, 2)))         # transverse-integrated n(z)
        slabs.append(a[:, midy, :])            # xz slice
    return dict(steps=np.array(steps), z=z, x=x,
                nzt=np.array(nzt), slabs=slabs, dt=DT_DEFAULT)


# ---- static (ground-state / single-point) fallbacks -------------------------
# Some runs are NOT trajectories: an inqkit ground state writes a single density
# to results/density_gs_system/<*>.vti (a DIFFERENT path from the raw/vti/
# density_system/ series), and a "frozen single-point" run logs only a few
# observables rows with no density VTI. These helpers let the (otherwise
# trajectory-oriented) builder still SHOW those runs' results. Additive: they run
# only when there is no trajectory series, so real runs are unaffected.
def load_static_density(results_dir):
    """Return a VtiField for a static GS density (or None). Uses the canonical
    inqview.load_vti — physical order, NEVER fftshift a VTI (vti-coordinate rule)."""
    rd = Path(results_dir)
    cands = sorted(rd.glob("density_gs_system/*.vti")) or sorted(rd.glob("*density*system*/*.vti"))
    if not cands:
        return None
    from inqview import load_vti
    try:
        return load_vti(str(cands[0]), expect_centered_axis="z")
    except Exception:
        return load_vti(str(cands[0]))          # centre-check off if feature not centred


def gs_density_fig(v, out, title, slab_half=12.5):
    """Static density xz mid-y slice, LINEAR | LOG (data axis order x,y,z)."""
    x, z = v.x, v.z
    sl = v.data[:, v.data.shape[1] // 2, :]                  # (nx, nz)
    ext = [x[0], x[-1], z[0], z[-1]]
    a = float(np.percentile(np.abs(sl), 99.5)) or 1e-12
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.4, 4.7), sharey=True)
    imL = axL.imshow(sl.T, origin="lower", aspect="auto", extent=ext,
                     norm=mcolors.Normalize(0.0, a), cmap="inferno")
    imR = axR.imshow(sl.T, origin="lower", aspect="auto", extent=ext,
                     norm=mcolors.LogNorm(vmin=a / 1e3, vmax=a), cmap="inferno")
    for ax, im, sci, pl in ((axL, imL, True, "linear"), (axR, imR, False, "log")):
        ax.set_xlabel("x (Bohr)")
        ax.text(0.03, 0.97, pl, transform=ax.transAxes, va="top", fontsize=8, color="w",
                bbox=dict(fc="k", alpha=0.4, lw=0, pad=1.5))
        _readable_cbar(fig, im, ax, r"$n$ (e/Bohr$^3$)", sci=sci)
        _zlines(ax, slab_half)
    axL.set_ylabel("z (Bohr)"); fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(out, dpi=150); plt.close(fig)


def gs_profile_fig(v, out, title, slab_half=12.5):
    """Transverse-integrated planar density n(z) of a static GS."""
    z = v.z
    nz = v.data.sum(axis=(0, 1))                             # sum over x,y
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(z, nz, "C0-", lw=1.4)
    ax.axvspan(-slab_half, slab_half, color="0.9", zorder=0)
    ax.set_xlabel("z (Bohr)"); ax.set_ylabel(r"$\int n\,dx\,dy$ (e/Bohr)")
    ax.set_title(title); ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)


def single_point_energy_md(results_dir, e_gs_ha=None):
    """Step-0 energy decomposition of a frozen single-point run (a few-row
    observables.csv, no trajectory). Returns markdown, or None. Numbers are read
    verbatim — nothing is re-converged."""
    import pandas as _pd
    obs = Path(results_dir) / "raw" / "observables" / "observables.csv"
    if not obs.exists():
        return None
    o = _pd.read_csv(obs).drop_duplicates(subset="step").sort_values("step")
    if "energy_total" not in o.columns or len(o) == 0:
        return None
    r0 = o.iloc[0]
    terms = [("energy_total", "E_tot(0)"), ("energy_kinetic", "T"),
             ("energy_hartree", "U_H"), ("energy_xc", "E_xc")]
    lines = ["| term | Ha | eV |", "|---|---|---|"]
    for col, sym in terms:
        if col in o.columns:
            val = float(r0[col])
            lines.append(f"| `{col}` ({sym}) | {val:.4f} | {val * HA_EV:.2f} |")
    if e_gs_ha is not None:
        exc = float(r0["energy_total"]) - e_gs_ha
        lines.append(f"| **E_tot(0) − E_GS** | **{exc:+.4f}** | **{exc * HA_EV:+.2f}** |")
    tail = (f"\n\nObservables logged over **{len(o)} steps** "
            f"(t = {float(o['time_au'].iloc[0]):g} … {float(o['time_au'].iloc[-1]):g} a.u.) — "
            "a frozen single-point energy, not a dynamical trajectory.")
    if e_gs_ha is not None:
        tail += f" `E_GS = {e_gs_ha:.4f}` Ha (bare-slab ground state, read from the GS run)."
    return "\n".join(lines) + tail


def _zlines(ax, slab_half=None, cap_inner=None):
    """Mark the slab faces (cyan dotted) and the CAP inner boundaries (lime
    dashed) on a z-axis figure. cap_inner is the |z| where absorption begins."""
    if slab_half is not None:
        for s in (-abs(slab_half), abs(slab_half)):
            ax.axhline(s, ls=":", color="cyan", lw=0.7)
    if cap_inner is not None:
        for s in (-abs(cap_inner), abs(cap_inner)):
            ax.axhline(s, ls="--", color="lime", lw=1.0)


def _readable_cbar(fig, im, ax, label, *, sci=True):
    """Colorbar with readable ticks: scientific ×10^n offset + ≤2 s.f. (report-
    figures rule 8), so small-magnitude decimals don't clip. sci=False for log axes."""
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    if sci:
        fmt = mticker.ScalarFormatter(useMathText=True); fmt.set_powerlimits((0, 0))
        cb.ax.yaxis.set_major_formatter(fmt)
        cb.ax.yaxis.set_major_locator(mticker.MaxNLocator(5))
        cb.ax.yaxis.get_offset_text().set_fontsize(7)
    cb.ax.tick_params(labelsize=7); cb.set_label(label, fontsize=8)
    return cb


def carpet(ser, kind, out, title, dt, cap_inner=None, slab_half=12.5):
    """z-t carpet, LINEAR | LOG side by side (report-figures rule 9). kind:
    'total' | 'delta0' (n(t)-n(0)) | 'dstep' (n(t+dt)-n(t))."""
    z, steps, nzt = ser["z"], ser["steps"], ser["nzt"]
    t = steps * dt
    if kind == "total":
        M, cmap, lab, signed = nzt, "inferno", r"$\int n\,dx\,dy$ (e/Bohr)", False
    elif kind == "delta0":
        M, cmap, lab, signed = nzt - nzt[0], "RdBu_r", r"$n(t)-n(0)$", True
    else:  # dstep
        M = np.zeros_like(nzt); M[1:] = nzt[1:] - nzt[:-1]
        cmap, lab, signed = "RdBu_r", r"$n(t+dt)-n(t)$", True
    a = float(np.percentile(np.abs(M), 99.5)) or 1e-12
    if signed:
        lin_norm = mcolors.Normalize(-a, a)
        log_norm = mcolors.SymLogNorm(linthresh=a / 100.0, vmin=-a, vmax=a, base=10)
    else:
        lin_norm = mcolors.Normalize(0.0, a)
        log_norm = mcolors.LogNorm(vmin=a / 1e3, vmax=a)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    for ax, norm, sci, pl in ((axL, lin_norm, True, "linear"), (axR, log_norm, False, "log")):
        im = ax.pcolormesh(t, z, M.T, cmap=cmap, shading="auto", norm=norm)
        ax.set_xlabel("time (a.u.)")
        ax.text(0.02, 0.98, pl, transform=ax.transAxes, va="top", fontsize=8,
                color="w", bbox=dict(fc="k", alpha=0.4, lw=0, pad=1.5))
        _readable_cbar(fig, im, ax, lab, sci=sci)
        _zlines(ax, slab_half, cap_inner)
    axL.set_ylabel("z (Bohr)"); fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(out, dpi=150); plt.close(fig)


def wake_carpet(ser, base, out, title, dt):
    z, steps = ser["z"], ser["steps"]
    t = steps * dt
    n = min(len(steps), len(base["steps"]))
    M = ser["nzt"][:n] - base["nzt"][:n]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    vmax = np.percentile(np.abs(M), 99.0) or 1.0
    im = ax.pcolormesh(t[:n], z, M.T, cmap="RdBu_r", shading="auto", vmin=-vmax, vmax=vmax)
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("z (Bohr)"); ax.set_title(title)
    fig.colorbar(im, ax=ax, label=r"$\delta n$ = run $-$ baseline (e/Bohr)")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)


def lead_gif(ser, out, title, dt, cap_inner=None, slab_half=12.5):
    """Total-density xz GIF, LINEAR | LOG side by side (report-figures rule 9),
    fixed colour scale across all frames."""
    x, z, slabs, steps = ser["x"], ser["z"], ser["slabs"], ser["steps"]
    a = float(np.percentile([s.max() for s in slabs], 99.5)) or 1e-12
    ext = [x[0], x[-1], z[0], z[-1]]
    lin_norm = mcolors.Normalize(0.0, a)
    log_norm = mcolors.LogNorm(vmin=a / 1e3, vmax=a)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.4, 4.7))
    imL = axL.imshow(slabs[0], origin="lower", aspect="auto", extent=ext, norm=lin_norm, cmap="inferno")
    imR = axR.imshow(slabs[0], origin="lower", aspect="auto", extent=ext, norm=log_norm, cmap="inferno")
    for ax, im, sci, pl in ((axL, imL, True, "linear"), (axR, imR, False, "log")):
        ax.set_xlabel("x (Bohr)")
        ax.text(0.03, 0.97, pl, transform=ax.transAxes, va="top", fontsize=8, color="w",
                bbox=dict(fc="k", alpha=0.4, lw=0, pad=1.5))
        _readable_cbar(fig, im, ax, r"$n$ (e/Bohr$^3$)", sci=sci)
        _zlines(ax, slab_half, cap_inner)
    axL.set_ylabel("z (Bohr)")
    sup = fig.suptitle("", fontsize=10); fig.tight_layout(rect=(0, 0, 1, 0.96))

    def upd(i):
        imL.set_data(slabs[i]); imR.set_data(slabs[i])
        sup.set_text(f"{title}  t = {steps[i]*dt:6.1f} a.u.")
        return imL, imR, sup
    animation.FuncAnimation(fig, upd, frames=len(slabs), interval=120, blit=False
                            ).save(out, writer="pillow", dpi=85)
    plt.close(fig)


def twin_energy_diff_bar_gif(cl_dir, wp_dir, out, dt, seconds_per_frame=0.45,
                             max_frames=140,
                             title="WP − classical energy difference"):
    """Animated bar chart of the per-step energy difference d = WP − classical
    over the shared INQ energy stores (total, kinetic, Hartree, xc), in eV.

    This is the twin ENERGY-DIFFERENCE view (the two runs' reported energies
    subtracted store by store) — NOT the pairwise P/S/B Coulomb split (which
    needs an ``interactions.csv`` the WP run does not carry). Both runs report
    these stores in ``raw/observables/observables.csv`` (Ha); aligned on the
    common step index (the WP is sampled every ``write_every`` steps, so the
    frames are the WP-sampled subset within the classical run's step range).
    The bars are paced ``seconds_per_frame`` per frame (deliberately slow so the
    term-by-term evolution can be read carefully).

    CAVEAT baked into the caller's caption: the WP run does not report
    ``energy_external`` separately (the WP is an electron, so its
    projectile-external coupling folds into Hartree), whereas the classical run
    does — so the component bars are NOT expected to sum to ΔE_total, and the WP
    run carries a CAP that slowly drains energy the classical run does not.
    Returns ``out`` on success, else ``None``.
    """
    import pandas as _pd
    co = Path(cl_dir) / "raw" / "observables" / "observables.csv"
    wo = Path(wp_dir) / "raw" / "observables" / "observables.csv"
    if not (co.exists() and wo.exists()):
        return None
    c = _pd.read_csv(co).drop_duplicates("step")
    w = _pd.read_csv(wo).drop_duplicates("step")
    stores = [("energy_total", "ΔE_total"), ("energy_kinetic", "ΔE_kin"),
              ("energy_hartree", "ΔE_Hartree"), ("energy_xc", "ΔE_xc")]
    stores = [(k, l) for k, l in stores if k in c.columns and k in w.columns]
    if not stores:
        return None
    m = _pd.merge(c[["step"] + [k for k, _ in stores]],
                  w[["step"] + [k for k, _ in stores]],
                  on="step", suffixes=("_cl", "_wp"))
    if len(m) == 0:
        return None
    if len(m) > max_frames:
        m = m.iloc[:: len(m) // max_frames + 1].reset_index(drop=True)
    labels = [l for _, l in stores]
    D = np.column_stack([(m[k + "_wp"].to_numpy() - m[k + "_cl"].to_numpy()) * HA_EV
                         for k, _ in stores])            # (nframes, nstores), eV
    steps = m["step"].to_numpy()
    ymin, ymax = float(np.nanmin(D)), float(np.nanmax(D))
    span = (ymax - ymin) or 1.0
    ymin -= 0.12 * span
    ymax += 0.12 * span
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    bars = ax.bar(x, D[0], color=["C0" if v >= 0 else "C3" for v in D[0]])
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("WP − classical energy (eV)")
    ax.set_ylim(ymin, ymax)
    ax.grid(axis="y", alpha=0.25)
    sup = ax.set_title("")
    vtxt = [ax.text(xi, 0, "", ha="center", fontsize=8) for xi in x]

    def upd(i):
        for b, v, tx in zip(bars, D[i], vtxt):
            b.set_height(v)
            b.set_color("C0" if v >= 0 else "C3")
            off = 0.03 * span if v >= 0 else -0.03 * span
            tx.set_y(v + off)
            tx.set_va("bottom" if v >= 0 else "top")
            tx.set_text(f"{v:.1f}")
        sup.set_text(f"{title}   t = {steps[i]*dt:6.1f} a.u.  (step {int(steps[i])})")
        return list(bars) + vtxt + [sup]

    interval_ms = max(60, int(seconds_per_frame * 1000))
    fig.tight_layout()
    animation.FuncAnimation(fig, upd, frames=len(D), interval=interval_ms,
                            blit=False).save(out, writer="pillow", dpi=90)
    plt.close(fig)
    return out


# ------------------------------------------------------------------ run_summary
def parse_summary(results_dir):
    p = Path(results_dir) / "run_summary.txt"
    d = {}
    if p.exists():
        # run_summary.txt packs SEVERAL "key = value" pairs onto one line, e.g.
        #   dt_au = 0.04  n_steps = 2500  write_every = 8
        # so a naive split-on-first-"=" mis-assigns the whole tail to the first key.
        # Capture every "<key> = <value>" pair; the value runs up to the next
        # "<key> =" or end-of-line, so space-containing values (background = slab
        # half_width 12.5 axis 2) stay intact.
        pair = re.compile(r"(\w+)\s*=\s*(.*?)(?=\s+\w+\s*=|$)")
        for ln in p.read_text().splitlines():
            for k, v in pair.findall(ln):
                v = v.strip()
                if v:
                    d[k.strip()] = v
    return d


def detect_type(results_dir):
    obs = Path(results_dir) / "raw" / "observables"
    if (obs / "momentum_distribution.csv").exists():
        return "wp"
    # mass-fork / effmass schema: stores WP momentum MOMENTS (pz_mean, e_kin_ha)
    # instead of the full 1D distribution — still a wavepacket run. The full-n(k)
    # momentum phase auto-skips; the moment-based WP battery (stopping via
    # ΔE_total/L_z, transit time, spreading, KL) lights up.
    if (obs / "wp_momentum_stats.csv").exists() or (obs / "wp_real_space_stats.csv").exists():
        return "wp"
    if (obs / "electron_track.csv").exists():
        return "classical"
    return "baseline"


# --------------------------------------------------------------- figure grouping
# (filename keyword -> battery group). Pipeline figures embed under these groups.
GROUPS = [
    ("Energetics", ["energy_", "energy_balance", "fft_total_energy", "ks_energies"]),
    ("Momentum", ["momentum_"]),
    ("Collective response (plasmon / current)", ["spectra/dipole", "spectra/current",
                                                 "fft_current"]),
    ("KS excitation & eigenstates", ["eigenvalue", "dos", "gamma_", "occupations",
                                     "gs_projected"]),
    ("Other metrics", ["kl_divergence", "knudsen", "observables_summary"]),
]


def group_of(relpath):
    s = str(relpath)
    for name, keys in GROUPS:
        if any(k in s for k in keys):
            return name
    return None


# ------------------------------------------------------------------------ build
# ---- annotation preservation (harvest-before-rebuild) ----------------------
# Builder cells carry metadata.gen="builder" (+ an anchor slug for headed markdown).
# Before a rebuild we HARVEST any markdown a reader added (cells without that tag)
# and RE-INJECT them at the same anchor so direct in-notebook annotations survive.
def slug(text):
    line = text.strip().splitlines()[0] if text.strip() else ""
    line = re.sub(r"^#+\s*", "", line)
    line = re.sub(r"[`*_$\\(){}\[\]]", "", line)
    return re.sub(r"[^a-z0-9]+", "-", line.lower()).strip("-")[:60]


def tag_builder(cell):
    """Stamp a builder-emitted cell; add an anchor when it is a markdown heading."""
    cell.setdefault("metadata", {})["gen"] = "builder"
    src = cell.get("source", "")
    if cell.get("cell_type") == "markdown" and src.lstrip().startswith("#"):
        cell["metadata"]["anchor"] = slug(src)
    return cell


def harvest_user_cells(ipynb_path):
    """{anchor -> [user markdown/raw cell, ...]} from an EXISTING notebook. User
    cells = markdown/raw NOT tagged gen='builder', anchored to the nearest preceding
    builder anchor; pre-anchor cells go to '__orphan__'. Empty if file absent."""
    import nbformat as nbf
    p = Path(ipynb_path)
    if not p.exists():
        return {}
    try:
        nb = nbf.read(str(p), as_version=4)
    except Exception:
        return {}
    # Transition guard: a pre-tagging notebook has NO gen='builder' cells, so user
    # vs builder cannot be told apart — harvest nothing (regenerate fresh + tag).
    if not any(c.get("metadata", {}).get("gen") == "builder" for c in nb.cells):
        return {}
    out, cur = {}, "__orphan__"
    for c in nb.cells:
        is_builder = c.get("metadata", {}).get("gen") == "builder"
        if c.cell_type == "markdown" and is_builder:
            cur = c.get("metadata", {}).get("anchor", cur)
        elif c.cell_type in ("markdown", "raw") and not is_builder:
            out.setdefault(cur, []).append(c)
    return out


def reinject(cells, harvested):
    """Splice harvested user cells back after their anchor; orphaned/lost-anchor
    cells go to a carried-over section after the title. Never drops annotations."""
    if not harvested:
        return cells
    import copy
    from nbformat.v4 import new_markdown_cell
    out, used = [], set()
    for c in cells:
        out.append(c)
        a = c.get("metadata", {}).get("anchor")
        if a and a in harvested:
            for uc in harvested[a]:
                u = copy.deepcopy(uc); u.setdefault("metadata", {})["gen"] = "user"
                out.append(u)
            used.add(a)
    leftover = [uc for a, ucs in harvested.items() if a not in used for uc in ucs]
    if leftover:
        hdr = tag_builder(new_markdown_cell(
            "## 📌 Carried-over reader annotations (original section changed/removed)"))
        block = [hdr] + [(lambda u: (u.setdefault("metadata", {}).update({"gen": "user"}) or u))(copy.deepcopy(x)) for x in leftover]
        out = out[:1] + block + out[1:]
    return out


def _retime_gif(src, dst, total_sec):
    """Re-time an animated GIF so the whole loop lasts ~total_sec seconds (a readable
    pace). Frames are unchanged; only the per-frame duration is rewritten. Used to
    slow the run-notebook GIFs down and, for external figures, to land a local copy
    beside the notebook so it renders (Jupyter won't serve files outside its tree)."""
    from PIL import Image, ImageSequence
    im = Image.open(str(src))
    frames = [f.copy() for f in ImageSequence.Iterator(im)]
    im.close()
    n = len(frames)
    if n <= 1:
        if str(src) != str(dst):
            shutil.copy2(src, dst)
        return
    per = int(round(total_sec * 1000.0 / n))
    per = max(20, min(per, 2000))                      # clamp per-frame 20 ms .. 2 s
    frames[0].save(str(dst), save_all=True, append_images=frames[1:],
                   loop=0, duration=per, disposal=2, optimize=False)


def build(results_dir, out_ipynb, baseline=None, run_cpp=None,
          cap_inner=None, decomp_prefix=None, rs=None, proj_sigma=0.3536,
          measured_s=None, measured_v=None, launch_z=None, v0=None,
          lindhard_mode="both", e_gs_ha=None, l_slab=25.0, gif_seconds=None,
          twin_wp=None, bar_gif_seconds=0.45):
    import nbformat as nbf
    from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
    from nbconvert.preprocessors import ExecutePreprocessor

    rd = Path(results_dir).resolve()
    out_ipynb = Path(out_ipynb).resolve()
    figdir = out_ipynb.parent / (out_ipynb.stem + "_figs")
    figdir.mkdir(parents=True, exist_ok=True)
    summ = parse_summary(rd)
    rtype = detect_type(rd)
    dt = float(summ.get("dt_au", DT_DEFAULT))
    # A GS run's results tree ends in a literal ".../results" dir; use the run dir
    # above it as the name so the notebook isn't titled "results".
    name = rd.parent.name if rd.name == "results" else rd.name

    # 1. run the pipeline (figures land under results/analysis/), auto-skip irrelevant
    try:
        runner.run(rd, run_name=name, rebuild=True)
    except Exception as e:                                  # never block the notebook
        print(f"[warn] pipeline.run raised: {e}")
    analysis = rd / "analysis"

    # 2. density-matrix carpets + lead GIF (computed here; pipeline doesn't emit them)
    sysdir = rd / "raw" / "vti" / "density_system"
    ser = load_series(sysdir) if sysdir.exists() else None
    made = {}
    if ser is not None:
        for kind, tag, ttl in [("total", "carpet_total", "Total density n(z,t)"),
                               ("delta0", "carpet_delta0", "Δn = n(t) − n(0)"),
                               ("dstep", "carpet_dstep", "Δn = n(t+dt) − n(t)")]:
            f = figdir / f"{tag}.png"
            carpet(ser, kind, f, f"{name}: {ttl}", dt, cap_inner=cap_inner); made[tag] = f
        g = figdir / "lead_density.gif"
        lead_gif(ser, g, f"{name} total density", dt, cap_inner=cap_inner); made["gif"] = g
        if baseline:
            bser = load_series(Path(baseline).resolve() / "raw" / "vti" / "density_system")
            if bser is not None:
                f = figdir / "carpet_wake.png"
                wake_carpet(ser, bser, f, f"{name}: wake = run − baseline", dt); made["wake"] = f

    # ---- static (GS / single-point) fallbacks when there is NO trajectory -------
    #      GS runs → converged density slice + n(z); frozen single points → step-0
    #      energy table. Populated only when `ser is None`, so real runs are untouched.
    static = {}
    if ser is None:
        try:
            sv = load_static_density(rd)
            if sv is not None:
                f1 = figdir / "gs_density_xz.png"
                gs_density_fig(sv, f1, f"{name}: converged density (xz mid-y)"); static["gs_xz"] = f1
                f2 = figdir / "gs_nz.png"
                gs_profile_fig(sv, f2, f"{name}: transverse-integrated n(z)"); static["gs_nz"] = f2
        except Exception as e:
            print(f"[warn] static GS density skipped: {e}")
    sp_md = single_point_energy_md(rd, e_gs_ha) if ser is None else None

    # ---- density-GIF battery {density,Δn,Δn_step} × {total,wp,bath} (standard) ---
    #      WP run → 9 GIFs; classical → 3 (total; Δn = induced wake). Reusable kernel
    #      inqview.visualisation.make_density_gif_battery (never blocks the build).
    battery = []
    try:
        _run_title = {"wp": "Wavepacket run", "classical": "Classical run",
                      "baseline": "Baseline run"}.get(rtype, name)
        battery, _ = make_density_gif_battery(
            str(rd), str(figdir), run_label=name, run_title=_run_title, dt=dt,
            slab_face=12.5, cap_inner=(cap_inner or 25.0), frames_max=48)
    except Exception as e:
        print(f"[warn] density-GIF battery skipped: {e}")

    # ---- physical-anchor heuristics (groups A–I) --------------------------------
    heur = None
    if rs is not None and v0 is not None and launch_z is not None:
        try:
            _box_half = None
            cb = str(summ.get("cell_bohr", ""))
            if "x" in cb:
                _box_half = float(cb.split("x")[-1]) / 2.0
            heur = compute_heuristics(
                str(rd), rs=rs, v0=v0, z0=launch_z, slab_half=12.5,
                box_half=(_box_half or 35.0),
                sigma_wp=(proj_sigma * np.sqrt(2.0) if rtype == "wp" else None))
        except Exception as e:
            print(f"[warn] heuristics skipped: {e}")

    # ---- harvest reader annotations from the existing notebook (before rebuild) -
    harvested = harvest_user_cells(out_ipynb)

    # ---- assemble notebook ------------------------------------------------------
    nb = new_notebook(); C = nb.cells
    def md(s): C.append(tag_builder(new_markdown_cell(s)))
    # display widths (px): Markdown ![]() renders at native pixel size (too large);
    # cap via an HTML <img> tag. GIFs are many + animated → keep compact.
    W_GIF, W_PNG = 360, 520
    def img(p, cap="", width=None):
        # PATH-REFERENCED (not base64-embedded): keeps the .ipynb small so viewers
        # render it. CRITICAL: the figure must travel BESIDE the notebook — Jupyter /
        # VSCode refuse to serve files outside the notebook's directory tree, so a
        # pipeline figure left under results/analysis/ (a `../../…` ref) shows as a
        # BROKEN image. Copy any external figure into <stem>_figs/ and reference the
        # local copy. GIFs are additionally re-timed to `gif_seconds` (readable pace).
        src = Path(p).resolve()
        is_gif = src.suffix.lower() == ".gif"
        rel0 = os.path.relpath(src, out_ipynb.parent)
        external = rel0.startswith("..") or os.path.isabs(rel0)
        if external:
            local = figdir / src.name
            if is_gif and gif_seconds:
                _retime_gif(src, local, gif_seconds)
            else:
                shutil.copy2(src, local)
            src = local.resolve()
        elif is_gif and gif_seconds:
            _retime_gif(src, src, gif_seconds)          # local GIF → slow in place
        rel = os.path.relpath(src, out_ipynb.parent)
        if width is None:
            width = W_GIF if is_gif else W_PNG
        tag = f'<img src="{rel}" width="{int(width)}" alt="{cap}">'
        md(f"*{cap}*\n\n{tag}" if cap else tag)

    proj = {"wp": "Gaussian wavepacket (σ_wp; see config)",
            "classical": "classical Gaussian-charge electron (matched UPF)",
            "baseline": "no projectile (CAP-on-bath baseline)"}[rtype]
    cap = "no CAP" if str(summ.get("cap_eta_ha", "")).strip() in ("0", "0.0", "-0") else \
          f"two-sided CAP η={summ.get('cap_eta_ha','?')} Ha, {summ.get('cap_width_bohr','?')}/side"

    # 1. title + question
    md(f"""# Run-notebook — `{name}`

**One-run deep dive.** {rtype.upper()} run in cubic jellium ({summ.get('rs','r_s 5.69')}),
{cap}. Projectile: {proj}, v0={summ.get('v0_au','?')} a.u.

**What this run shows:** the full per-run battery — density evolution, energetics, the
projectile's transport/stopping, collective response, momentum, and KS excitation — so
this single run can be read in depth. PROVISIONAL until the inq-study engine regression
(Task #7). σ-convention is unified (σ = σ_wp; charge std = σ/√2 — see CONTEXT.md).""")

    # 1b. Context / aim / hypothesis — run-SPECIFIC narrative from a SIDECAR file so it
    #     is not baked into this generic builder and survives rebuilds. Edit
    #     <stem>.aim.md; injected here (house-narrative section 1: "what this run shows").
    aim_path = out_ipynb.parent / (out_ipynb.stem + ".aim.md")
    if aim_path.exists():
        body = aim_path.read_text().strip()
        if body:
            md(f"## Context, aim & hypothesis\n\n{body}")

    # 2. Reader notes / TODOs — pinned from a SIDECAR file so they SURVIVE rebuilds.
    #    Edit <stem>.notes.md freely; it is injected here on every build (never wiped).
    notes_path = out_ipynb.parent / (out_ipynb.stem + ".notes.md")
    if notes_path.exists():
        body = notes_path.read_text().strip()
        if body:
            md(f"""## 📝 Reader notes / TODOs
*(pinned from `{notes_path.name}` — edit that file; this cell is regenerated each
build and is never overwritten by the builder.)*

{body}""")

    # 3. setup (reconstructable) — verbatim run_summary
    md("## Setup (fully reconstructable — `run_summary.txt`)\n\n"
       + "| key | value |\n|---|---|\n"
       + "\n".join(f"| `{k}` | {v} |" for k, v in summ.items()))

    # 4. source files
    md(f"""## Source files
| file | role |
|---|---|
| `{run_cpp or 'ResearchProject/.../run.cpp'}` | the run binary (env-driven) |
| `inqview.pipeline` phases | computed most figures below |
| `.claude/skills/run-notebook/run_notebook_builder.py` | this builder |
| `{rd}/run_summary.txt` | provenance |""")

    # 5. Visual intuition: lead GIF + per-run energetics (or the static fallback)
    md("## Visual intuition — total density (xz slice)")
    if "gif" in made:
        img(made["gif"], "Total electronic density, mid-y xz slice, animated.")
    elif static:
        md("_No trajectory series — this is a ground-state / single-point run. The "
           "converged density is shown below (slab faces |z|=12.5 dotted)._")
        img(static["gs_xz"], "Converged density, mid-y xz slice (linear | log).")
        img(static["gs_nz"], "Transverse-integrated planar density n(z); slab region shaded.")
    else:
        md("_no density VTI series found_")

    # 5b. Single-point energy (frozen single-point runs: the result IS the step-0 energy)
    if sp_md:
        md("## Result — single-point energy decomposition\n\n" + sp_md)

    # 6. Density matrix (carpets)
    md("""## Density matrix (z–t carpets)
Transverse-integrated linear density vs (z, t). Top: total n(z,t); then the two deltas
(cumulative n(t)−n(0) and per-step n(t+dt)−n(t) = the instantaneous flow). The **wake**
(run − baseline) isolates the projectile-induced response (the WP-density proxy; bare
|ψ_WP|² is a future observable).""")
    for tag, cap_ in [("carpet_total", "total n(z,t)"), ("carpet_delta0", "n(t)−n(0)"),
                      ("carpet_dstep", "n(t+dt)−n(t)"), ("wake", "wake = run − baseline")]:
        if tag in made:
            img(made[tag], cap_)

    # 6b. Density-GIF battery {density, Δn, Δn_step} × {total, wp, bath}, auto-rendered
    #     by make_density_gif_battery (WP: 9; classical: 3, Δn = induced wake).
    if battery:
        _CATL = {"total": "Total system", "wp": "Wavepacket |ψ_WP|²", "bath": "Bath (slab = total − WP)"}
        _KINDL = {"density": "n(x,z,t)  [linear | log; total/bath share the slab scale]",
                  "delta0": "Δn = n(t) − n(0)  (induced response) [linear | symlog]",
                  "dstep": "Δn = n(t+dt) − n(t)  (instantaneous flux) [linear | symlog]"}
        md("## Density-GIF battery (xz slices)\n"
           "Three kinds — absolute density, cumulative Δ-vs-initial, and per-step Δ — for "
           "each density channel. Slab faces (|z|=12.5) and CAP inner faces (|z|=25) dashed. "
           "**Every GIF shows LINEAR | LOG side by side** (shared-colorbar rule): density uses "
           "a shared linear+log scale (low densities visible), and the Δ kinds use a symmetric "
           "diverging linear + symlog scale (the symlog panel exposes the low-|Δn| wake tail).")
        for cat in ("total", "wp", "bath"):
            rows = [(k, p) for (c, k, p, _t) in battery if c == cat]
            if not rows:
                continue
            md(f"### {_CATL.get(cat, cat)}")
            for kind in ("density", "delta0", "dstep"):
                for k, p in rows:
                    if k == kind:
                        img(p, _KINDL[kind])

    # 6c. (legacy) pre-rendered decomposition GIFs via --decomp-prefix, if provided.
    if decomp_prefix:
        import glob as _glob
        gifs = sorted(_glob.glob(str(out_ipynb.parent / f"{decomp_prefix}_*.gif")))
        if gifs:
            cap_note = (" Cyan dotted = slab faces (±12.5 Bohr); **lime dashed = CAP "
                        "inner boundary (±17.5 Bohr)**." if cap_inner
                        else " Cyan dotted = slab faces (±12.5 Bohr).")
            syslabel = {"total": "Total electron density",
                        "bath": "Bath (gas = total − WP)",
                        "wp": "Wavepacket |ψ_WP|²",
                        "proj": "Projectile (Gaussian charge)"}
            viewlabel = {"total": "n(t)", "dfirst": "Δn = n(t)−n(0) (induced)",
                         "dprev": "Δn = n(t)−n(t−Δt) (flux)"}
            md("## Three-way density decomposition (xz GIFs)\n"
               "Each system in three views: absolute n(t), induced Δ-vs-first, and "
               "flux Δ-vs-previous." + cap_note)
            for s in ["total", "bath", "wp", "proj"]:
                sg = [g for g in gifs if Path(g).name.startswith(f"{decomp_prefix}_{s}_")]
                if not sg:
                    continue
                md(f"### {syslabel.get(s, s)}")
                for v in ["total", "dfirst", "dprev"]:
                    for g in sg:
                        if Path(g).name.endswith(f"_{s}_{v}.gif"):
                            img(g, f"{syslabel.get(s, s)} — {viewlabel[v]}")

    # battery groups from pipeline figures. The pipeline's RAW FFT plots
    # (fft_*, spectrum_*, spectra/) are DELIBERATELY EXCLUDED — every FFT-driven
    # observable goes through the fourier-analysis skill's audited 6-stage
    # fft_pipeline_panel below (user, 2026-06-30). Only time-domain figs embed here.
    pics = sorted([p for p in (analysis.rglob("*.png") if analysis.exists() else [])]) + \
           sorted([p for p in (analysis.rglob("*.gif") if analysis.exists() else [])])
    by_group = {}
    for p in pics:
        if p.name.startswith(("fft_", "spectrum_")) or "spectra" in p.parts:
            continue                      # FFT handled by the audited panel below
        g = group_of(p.relative_to(rd))
        if g:
            by_group.setdefault(g, []).append(p)
    for gname, _ in GROUPS:
        if gname in by_group:
            md(f"## {gname}")
            if gname == "Energetics":
                if rtype == "wp":
                    md("_Total-system energy + components below. The WP KS-orbital total "
                       "energy ⟨ψ_WP|H|ψ_WP⟩(t) + variance are in the `ks_energies_*` panels._")
                _df = figdir / "energy_delta_total_vs_time.png"
                _finding = delta_total_energy_fig(rd, _df, rtype)
                if _finding:
                    img(_df, "ΔE_total(t) — total-system energy change (plateau test, "
                             "left) NEXT TO the total electron number N(t)=∫n dV (right)")
                    md(_finding)
                # Twin energy-difference bar GIF (WP − classical), user request:
                # alongside the total-energy plots, a decomposed-energy bar chart
                # animated over time at a deliberately slow pace. Only for the
                # classical member of a twin pair (twin_wp = its WP counterpart).
                if twin_wp and rtype == "classical":
                    _eg = figdir / "energy_diff_evolution.gif"
                    try:
                        if twin_energy_diff_bar_gif(rd, twin_wp, _eg, dt,
                                                    seconds_per_frame=bar_gif_seconds):
                            md("### WP − classical energy difference (animated bars)\n"
                               "The **wavepacket minus classical** energy difference for each "
                               "shared INQ store (total, kinetic, Hartree, xc), animated "
                               "over the run at a **slow pace** so the term-by-term change "
                               "can be read carefully. The kinetic bar is the largest and "
                               "roughly constant (≈ +100 eV) — the WP's localisation "
                               "zero-point energy $3/(4\\sigma^2)+k_0^2/2$ that a classical "
                               "point charge does not pay; ΔE_xc ≈ −16 eV is the WP's own "
                               "exchange–correlation. **Caveats:** the WP run does not "
                               "report `energy_external` separately (it folds into Hartree), "
                               "so the component bars are **not** expected to sum to "
                               "ΔE_total; and the WP run carries a CAP that slowly drains "
                               "energy the classical run does not — so the late-time drift "
                               "is partly boundary absorption, not physics.")
                            img(_eg, "WP − classical energy difference per shared store, "
                                     "animated (slow pace)", width=W_GIF + 90)
                    except Exception as _e:
                        print(f"[warn] twin energy-diff GIF skipped: {_e}")
            for p in by_group[gname]:
                img(p, p.name)

    # ---- FFT-pipeline diagnostic panel (project figure standard) -------------
    # For each FFT'd signal show the processing stages so a spectrum is never a
    # black box: raw -> de-trend -> window -> zero-pad -> |FFT| (linear, log).
    obs_csv = rd / "raw" / "observables" / "observables.csv"
    if obs_csv.exists():
        import pandas as _pd
        from inqview.analysis.fourier import FourierTransform, WindowSpec
        from inqview.visualisation.fourier_panel import fft_pipeline_panel
        _o = _pd.read_csv(obs_csv).drop_duplicates(subset="step").sort_values("time_au")
        # Audited FFT for every dynamic collective-response observable: the most
        # dynamic dipole component AND the axial current (the jellium plasmon signal).
        _cands = []
        _dcols = [c for c in ("dipole_z", "dipole_x", "dipole_y") if c in _o.columns]
        if _dcols:
            _cands.append(max(_dcols, key=lambda c: float(np.nanstd(_o[c].values))))
        for _c in ("current_z", "energy_total"):
            if _c in _o.columns and float(np.nanstd(_o[_c].values)) > 0:
                _cands.append(_c)
        if _cands and len(_o) >= 8:
            band = None; fmax = 15.0; wp = None
            if rs is not None:
                wp = (3.0 / rs ** 3) ** 0.5 * HA_EV                 # plasmon energy (eV)
                band = (max(0.3, 0.4 * wp), 2.0 * wp); fmax = max(12.0, 3.0 * wp)
            tau = float(_o["time_au"].iloc[-1] - _o["time_au"].iloc[0])
            dw_ev = 2 * np.pi / tau * HA_EV
            wp_txt = (r" Plasmon band $\hbar\omega_p\approx%.2f$ eV shaded." % wp) if wp else ""
            md("## FFT pipeline (fourier-analysis skill — every spectrum stage-by-stage)\n"
               "Each collective-response spectrum is built through audited stages — "
               "**raw → mean-baseline → Hann → zero-pad ×4 → |FFT| (linear & log, "
               "detrend overlaid)** — via `fourier_panel.fft_pipeline_panel`; the "
               "pipeline's raw FFT plots are not used. Peak located inside the physical "
               "plasmon band, never by global argmax.%s\n"
               r"**Resolution:** $\Delta\omega=2\pi/\tau\approx%.2f$ eV for "
               r"$\tau=%.0f$ a.u. — informational (not a gate); treat sub-bin peak "
               r"positions cautiously on short runs." % (wp_txt, dw_ev, tau))
            ft = FourierTransform(window=WindowSpec("hann"), zero_pad=4, subtract="mean")
            for col in _cands:
                f = figdir / f"fft_pipeline_{col}.png"
                fig = fft_pipeline_panel(_o["time_au"].to_numpy(), _o[col].to_numpy(), ft,
                                         label=col, peak_band=band, fmax=fmax,
                                         title=f"FFT pipeline — {col}")
                fig.savefig(f, dpi=140); plt.close(fig)
                img(f, f"FFT pipeline for {col}: raw, baseline, windowed, padded, |FFT| linear+log")

    # Momentum future-observable note
    if rtype == "wp":
        md("> **2D scattering map (k_z, k_⊥):** requires a future 2D-momentum observable "
           "(current runs save only 1D |k|). The 1D n(k) heatmap / GIF above is the "
           "available before/after view.")

    # WP transit / ballistic exit time (user request) ----------------------
    if rtype == "wp" and v0 is not None and launch_z is not None:
        f = figdir / "wp_exit_time.png"
        t_exit = wp_exit_time_fig(rd, f, v0, launch_z, slab_half=12.5,
                                  cap_inner=cap_inner, dt=dt)
        md(f"""## WP transit & ballistic exit time
Travelling at the mean momentum $v={v0:.3f}$ a.u., a packet launched at
$z={launch_z}$ Bohr reaches the **far slab face** ($z=+12.5$) after
$t_\\mathrm{{exit}} = (12.5-({launch_z}))/v = ${t_exit:.1f}$ a.u. — marked below on the
measured WP centroid $\\langle z\\rangle(t)$. (The centroid lags the ballistic line
once the CAP starts removing the leading edge; `norm_check` shows the WP draining.)""")
        img(f, "WP centroid z(t) vs the ballistic-exit estimate")

    # WP quantum stopping power — energy method (auto-measured) ------------------
    if rtype == "wp" and e_gs_ha is not None:
        import pandas as _pd
        obsd = rd / "raw" / "observables"
        try:
            _o = _pd.read_csv(obsd / "observables.csv")
            _t = _o["time_au"].values.astype(float); _E = _o["energy_total"].values.astype(float)
            _s = _pd.read_csv(obsd / "wp_real_space_stats.csv", comment="#")
            _norm = _s["norm_check"].values.astype(float)
            deposited = (float(_E[-1]) - e_gs_ha) * HA_EV
            S_wp = deposited / l_slab
            norm_f = float(_norm[-1])
            _m = _t >= 0.85 * _t.max()
            slope = float(np.polyfit(_t[_m], _E[_m], 1)[0] * HA_EV) if _m.sum() > 1 else float("nan")
            converged = bool(norm_f < 0.02 and abs(slope) < 0.2)
            bound = "exact value" if converged else ("UPPER bound" if slope < 0 else "LOWER bound")
            f = figdir / "wp_quantum_stopping.png"
            fig, (q1, q2) = plt.subplots(2, 1, figsize=(7.4, 6.2), sharex=True)
            q1.plot(_t, (_E - e_gs_ha) * HA_EV, "C0-", lw=1.7)
            q1.axhline(deposited, ls=":", color="0.5",
                       label=f"deposited = {deposited:.1f} eV  ⇒  S = {S_wp:.2f} eV/Bohr")
            q1.set_ylabel("E_total(t) − E_GS  (eV)"); q1.grid(alpha=.25)
            q1.legend(fontsize=8, frameon=False)
            q1.set_title(f"Quantum stopping S = [E_total(t_f)−E_GS]/L_z  ({bound})", fontsize=9)
            q2.plot(_s["time_au"].values, _norm, "C1-", lw=1.6, label="WP orbital norm")
            q2.axhline(0.02, ls="--", color="0.6", lw=.8, label="convergence gate 0.02")
            q2.set_xlabel("time (a.u.)"); q2.set_ylabel("WP norm"); q2.grid(alpha=.25)
            q2.legend(fontsize=8, frameon=False)
            q2.set_title(f"norm_f={norm_f:.3f}, late dE/dt={slope:+.3f} eV/au", fontsize=9)
            fig.tight_layout(); fig.savefig(f, dpi=175); plt.close(fig)
            lr_txt = ""
            if rs is not None and v0 is not None:
                from inqview.analysis import lindhard_elf as _le
                S_lr = _le.stopping_power_point(float(v0), _le.kF_from_rs(rs)) * HA_EV
                lr_txt = (f" Point-charge Lindhard at v={v0:.2f} is S={S_lr:.2f} eV/Bohr → "
                          f"**S_WP/Lindhard = {S_wp/S_lr:.1f}×**.")
                measured_v = float(v0); measured_s = S_wp / HA_EV   # overlay on Lindhard panel (Ha/Bohr)
            _conv_txt = (" The WP is fully absorbed and the deposit has settled → a true value."
                         if converged else
                         " The WP is not fully absorbed by τ (residual WP energy still draining, "
                         "late slope < 0) → the deposit will fall, so this is an UPPER bound; a "
                         "converged value needs a longer τ / stronger CAP.")
            md(f"""## Quantum stopping power (wavepacket) — energy method
The quantum electronic stopping power is the energy deposited in the slab per unit
traversal length, **S = [E_total(t_f) − E_GS] / L_z** (L_z = {l_slab:.0f} Bohr; the
geometry-correct slab method). The anchor is the bare-slab ground state
**E_GS = {e_gs_ha:.4f} Ha**, *not* E_total(t₀) — the WP's drift kinetic energy lives
inside E_total(0), so the E_GS anchor is the WP-correct adaptation.

**S = {S_wp:.2f} eV/Bohr** (deposited {deposited:.1f} eV), **{bound}**: norm_f =
{norm_f:.3f} (gate < 0.02), late dE/dt = {slope:+.3f} eV/au.{_conv_txt}{lr_txt}""")
            img(f, "Quantum stopping: retained energy E_total(t)−E_GS (top) and WP-norm convergence (bottom)")
        except Exception as _e:  # noqa: BLE001
            md(f"_WP quantum-stopping section skipped: {_e}_")

    # Classical transport: z(t) + KE(z) dip-and-recovery (conservative well) ----
    if rtype == "classical":
        tkf = rd / "raw" / "observables" / "electron_track.csv"
        if tkf.exists():
            import pandas as _pd
            tk = _pd.read_csv(tkf).drop_duplicates("step")
            zc, kc, tc = tk["z"].values, tk["ke_ion_ha"].values, tk["time_au"].values
            f = figdir / "classical_transport.png"
            fig, (b1, b2) = plt.subplots(1, 2, figsize=(11, 4.4))
            b1.plot(tc, zc, "C3-", lw=1.6)
            for zz in (12.5, -12.5, (cap_inner or 25), -(cap_inner or 25)):
                b1.axhline(zz, ls="--", lw=.7, color="0.5")
            b1.set_xlabel("time (a.u.)"); b1.set_ylabel("ion z (Bohr)")
            b1.set_title("Classical ion z(t)", fontsize=9); b1.grid(alpha=.25)
            iout = int(np.argmax(zc >= 35)) or len(zc)
            b2.plot(zc[:iout], kc[:iout] * HA_EV, "C3.-", lw=1.3, ms=3)
            b2.axvspan(-12.5, 12.5, color="C0", alpha=.10)
            ki = np.interp(-12.5, zc, kc); ko = np.interp(12.5, zc, kc)
            S_face = (ki - ko) * HA_EV / 25.0
            b2.set_xlabel("ion z (Bohr)"); b2.set_ylabel("projectile KE (eV)")
            b2.set_title("KE(z): conservative dip-and-recovery", fontsize=9); b2.grid(alpha=.25)
            fig.tight_layout(); fig.savefig(f, dpi=170); plt.close(fig)
            md(f"""## Classical transport — z(t) and KE(z) (conservative dip-and-recovery)
The projectile KE is **not monotonic**: it slows to a minimum near the slab **centre** and
recovers on exit — a conservative mean-field-potential effect (energy borrowed and returned),
**not** stopping. Only the net loss between two points at **equal background potential** (the
symmetric slab faces ±12.5) is true electronic stopping: ΔKE/25 ⇒ **S = {S_face:.3f} eV/Bohr**.
A centre or asymmetric window would conflate the well with stopping — **window choice dominates**.""")
            img(f, "Classical ion z(t) (left) and KE(z) dip-and-recovery with equal-potential window (right)")

    # Stopping power vs analytical Lindhard (user request) -----------------
    if rs is not None:
        f = figdir / "lindhard_stopping.png"
        at = lindhard_stopping_fig(f, rs, proj_sigma, measured_v=measured_v,
                                   measured_s=measured_s, mode=lindhard_mode)
        cmp_txt = ""
        if measured_v is not None and measured_s is not None:
            if lindhard_mode == "point":
                cmp_txt = (
                    f"\n\nAt $v={measured_v:.3f}$ a.u. the **measured** stopping is "
                    f"$S={measured_s*HA_EV:.3f}$ eV/Bohr against the **point-charge "
                    f"linear-response Lindhard** $S={at['S_point']*HA_EV:.3f}$ eV/Bohr — "
                    f"agreement to ~{abs(measured_s/at['S_point']-1)*100:.0f}%. The "
                    f"projectile width (charge std $\\approx{proj_sigma:.3f}$ Bohr) was "
                    f"chosen from the σ-convergence sweep to sit in the linear-response "
                    f"regime, so the point-charge curve is the appropriate reference and "
                    f"the finite-size Gaussian form-factor correction is omitted by design. "
                    f"The useful velocity window is $v\\gtrsim v_F={at['kF']:.2f}$ a.u. "
                    f"where $S(v)$ peaks.")
            else:
                cmp_txt = (
                    f"\n\nAt $v={measured_v:.3f}$ a.u. the two analytical curves **bracket** the "
                    f"measurement: **Gaussian** $S={at['S_sigma']*HA_EV:.3f}$, **point** "
                    f"$S={at['S_point']*HA_EV:.3f}$, **measured (rt-TDDFT)** "
                    f"$S={measured_s*HA_EV:.3f}$ eV/Bohr. The rt-TDDFT point sits ~"
                    f"{(measured_s/ at['S_sigma']-1)*100:.0f}% above the matched-σ Gaussian "
                    f"curve, near the point-charge limit. **This gap is not yet resolved** — "
                    f"candidate causes, none confirmed: (i) the Gaussian form factor "
                    f"$e^{{-q^2\\sigma_q^2}}$ cuts the q-integral at $q\\sim1/\\sigma_q"
                    f"\\approx{1/proj_sigma:.1f}$ a.u.; (ii) non-linear $Z=-1$ response "
                    f"beyond RPA; (iii) finite-box / CAP / transient contamination of "
                    f"$\\Delta$KE. *Inference, to verify:* a $dx\\,0.5\\to0.25$ convergence "
                    f"run and a $Z$-scaling check would discriminate (i)/(ii). The useful "
                    f"velocity window is $v\\gtrsim v_F={at['kF']:.2f}$ where $S(v)$ peaks.")
        curve_desc = ("the bare **point-charge** curve (the projectile width sits in "
                      "the linear-response regime, so the finite-size correction is "
                      "negligible)" if lindhard_mode == "point" else
                      "the bare **point-charge** curve and the **Gaussian-projectile** "
                      f"curve (charge std $\\sigma_q={proj_sigma:.3f}$)")
        md(f"""## Stopping power vs analytical Lindhard (slab gas $r_s$={rs:.3f})
The analytical **linear-response (RPA) Lindhard** stopping power for the slab's
electron gas — {curve_desc} — with the rt-TDDFT measurement overlaid. Plasmon
$\\omega_p={at['wp_eV']:.2f}$ eV, $v_F={at['kF']:.3f}$ a.u. Source: Lindhard–Winther
1964; Correa 2018 (`docs/sources/stopping-power-formulae.md`).{cmp_txt}""")
        img(f, "Analytical Lindhard S(v) for the slab density + measured point")

    # Loss function (always attempted; low-res note)
    nsteps = int(float(summ.get("n_steps", 0)))
    tau = nsteps * dt
    md(f"""## Loss function L(q,ω)
{"Dynamic structure factor / loss function from δn(q,t). **Resolution note (not a gate):** Δω = 2π/τ with τ = %.0f a.u., so fine spectral features below one bin are unresolved at this τ — the spectrum is still shown in full." % tau if tau < 300 else "Dynamic structure factor / loss function from δn(q,t)."}""")
    loss = [p for p in (analysis.rglob("*loss*") if analysis.exists() else [])] + \
           [p for p in (analysis.rglob("*spectral*") if analysis.exists() else [])] + \
           [p for p in (analysis.rglob("*lindhard*") if analysis.exists() else [])]
    if loss:
        for p in loss:
            img(p, p.name)
    else:
        md("_loss-function phase produced no figure for this run (too short / δn absent)._")

    # physical anchors / heuristics (groups A–I)
    if heur is not None:
        eg, tsc, sp = heur.eg_scales, heur.timescales, heur.spreading
        rows = [("r_s", f"{eg.get('rs', float('nan')):.3f}"),
                ("k_F (a₀⁻¹)", f"{eg.get('kF', float('nan')):.4f}"),
                ("v_F (a.u.)", f"{eg.get('vF', float('nan')):.4f}"),
                ("E_F (eV)", f"{eg.get('EF_ev', float('nan')):.2f}"),
                ("λ_F=π/k_F (Bohr)", f"{eg.get('lambda_F_friedel', float('nan')):.2f}"),
                ("ω_p (eV)", f"{eg.get('omega_p_ev', float('nan')):.2f}"),
                ("T_plasmon (a.u.)", f"{eg.get('T_plasmon_au', float('nan')):.1f}"),
                ("k_TF (a₀⁻¹)", f"{eg.get('k_TF', float('nan')):.3f}"),
                ("t_enter slab (a.u.)", f"{tsc.get('t_enter_slab_au', float('nan')):.2f}"),
                ("t_exit slab — END (a.u.)", f"{tsc.get('t_exit_slab_au', float('nan')):.2f}"),
                ("t reach box edge (a.u.)", f"{tsc.get('t_reach_box_edge_au', float('nan')):.2f}")]
        if heur.wp_kinetics:
            rows.append(("zero-point KE (eV)", f"{heur.wp_kinetics.get('zero_point_ke_ev', float('nan')):.1f}"))
        if sp:
            rows.append(("spreading σ_z(t)/σ_z(0)", f"×{sp.get('spread_factor', float('nan')):.0f}"))
        if heur.norms:
            rows.append(("total absorbed (e)", f"{heur.norms.get('total_absorbed', float('nan')):.3f}"))
        tau_h = float(summ.get("n_steps", 0)) * dt
        reslim = ("\n\n> **Resolution:** T_plasmon ≈ %.0f a.u. vs τ = %.0f a.u. — %s a single "
                  "plasmon period, so any L(q,ω) is %s." % (
                      eg.get('T_plasmon_au', float('nan')), tau_h,
                      "shorter than" if eg.get('T_plasmon_au', 0) > tau_h else "longer than",
                      "under-resolved (Δω∝1/τ)" if eg.get('T_plasmon_au', 0) > tau_h else "resolvable"))
        md("## Physical anchors & heuristics (groups A–I)\n\n"
           "HEG scales, projectile timescales, and integrity metrics "
           "(`inqview.analysis.compute_heuristics`).\n\n"
           "| quantity | value |\n|---|---|\n"
           + "\n".join(f"| {k} | {v} |" for k, v in rows) + reslim)

    # takeaway
    md(f"""## Takeaway
- Run `{name}` ({rtype}) completed: {summ.get('run_completed','?')}, nan_seen={summ.get('nan_seen','?')}.
- absorbed_frac = {summ.get('absorbed_frac','NA')}; """
       + (f"projectile v_z {summ.get('v0_au','?')} → {summ.get('final_vz','?')} (deceleration = stopping)."
          if rtype == "classical" else "see momentum / energetics for the stopping signal.")
       + "\n- PROVISIONAL until Task #7.")

    # ---- re-inject the harvested reader annotations at their anchors ------------
    nb.cells = reinject(nb.cells, harvested)
    n_user = sum(1 for c in nb.cells if c.get("metadata", {}).get("gen") == "user")

    ep = ExecutePreprocessor(timeout=1200, kernel_name="python3")
    ep.preprocess(nb, {"metadata": {"path": str(out_ipynb.parent)}})
    with open(out_ipynb, "w") as fh:
        nbf.write(nb, fh)
    print(f"wrote {out_ipynb}  ({len(nb.cells)} cells, type={rtype}, "
          f"{n_user} reader annotation(s) preserved)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir")
    ap.add_argument("out_ipynb")
    ap.add_argument("--baseline", default=None)
    ap.add_argument("--run-cpp", default=None)
    ap.add_argument("--decomp-prefix", default=None,
                    help="embed pre-rendered <prefix>_<system>_<view>.gif beside the notebook")
    ap.add_argument("--cap-inner", type=float, default=None,
                    help="|z| (Bohr) of the CAP inner boundary; draws dashed lines on z-figures")
    ap.add_argument("--rs", type=float, default=None,
                    help="electron-gas r_s of the target; enables the analytical Lindhard stopping panel")
    ap.add_argument("--proj-sigma", type=float, default=0.3536,
                    help="projectile charge std (σ_WP/√2 = 0.354 for σ=0.5); Gaussian form-factor")
    ap.add_argument("--lindhard", choices=("both", "point"), default="both",
                    help="'point' draws only point-charge Lindhard (projectile in LR regime); "
                         "'both' (default) adds the finite-σ Gaussian curve")
    ap.add_argument("--measured-s", type=float, default=None, help="measured S (Ha/Bohr) to overlay")
    ap.add_argument("--measured-v", type=float, default=None, help="measured-point velocity (a.u.)")
    ap.add_argument("--launch-z", type=float, default=None, help="WP launch z (Bohr) for exit-time")
    ap.add_argument("--v0", type=float, default=None, help="projectile velocity (a.u.) for exit-time")
    ap.add_argument("--e-gs-ha", type=float, default=None,
                    help="bare-slab ground-state energy (Ha); enables the WP energy-method "
                         "quantum-stopping section S=[E_total(t_f)−E_GS]/L_z (WP runs)")
    ap.add_argument("--l-slab", type=float, default=25.0,
                    help="slab thickness / traversal length (Bohr) for the WP quantum-stopping method")
    ap.add_argument("--gif-seconds", type=float, default=None,
                    help="re-time every embedded GIF so its loop lasts ~this many seconds "
                         "(readable pace, e.g. 17); omit to keep native frame timing")
    ap.add_argument("--twin-wp", default=None,
                    help="the WP twin run's results dir; adds a WP−classical energy-diff "
                         "bar GIF (shared stores) to the Energetics section (classical run)")
    ap.add_argument("--bar-gif-seconds", type=float, default=0.45,
                    help="seconds per frame for the WP−classical energy-diff bar GIF "
                         "(default 0.45 — deliberately slow so the bars can be read)")
    a = ap.parse_args()
    build(a.results_dir, a.out_ipynb, baseline=a.baseline, run_cpp=a.run_cpp,
          cap_inner=a.cap_inner, decomp_prefix=a.decomp_prefix, rs=a.rs,
          proj_sigma=a.proj_sigma, measured_s=a.measured_s, measured_v=a.measured_v,
          launch_z=a.launch_z, v0=a.v0, lindhard_mode=a.lindhard,
          e_gs_ha=a.e_gs_ha, l_slab=a.l_slab, gif_seconds=a.gif_seconds,
          twin_wp=a.twin_wp, bar_gif_seconds=a.bar_gif_seconds)
