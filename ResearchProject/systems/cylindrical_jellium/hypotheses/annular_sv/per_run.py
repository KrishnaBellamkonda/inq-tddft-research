#!/usr/bin/env python3
"""Per-run figure battery for the annular-tube S(v) sweep — TUBE-aware.

For one projectile run (classical or WP) this generates, into
``per_run_figs/<label>/``:

  1. **Density matrix** GIFs {density, Δn=n(t)−n(0), Δn=n(t+dt)−n(t)} × {total[,wp,bath]}
     — classical: 3 (total only); WP: 9 (bath = total − wp). xz mid-y slices with
     **vertical wall-radius markers x=±R_in,±R_out** (the tube geometry; NOT slab
     faces — there are none, and no CAP). The **moving projectile** (classical ion,
     or WP centroid) is overlaid from THAT run's own track.
  2. **z–t carpets** (total n, Δ-vs-0, per-step Δ) with the projectile z(t) overlaid.
  3. **Stopping power** via the ``stopping-power-extraction`` skill kernels (Correa
     2018): PRIMARY = electronic deposit slope ``dE_total(x)`` (Method A, continuous
     traversal); cross-checks = ``−dKE_ion/dx`` (energy conservation) and the
     ``N(t)≈const`` guard. The window is the **early v≥0.85·v0 segment** (the light
     free-Ehrenfest electron decelerates and stops — `.claude/rules/light-projectile-
     stopping.md`), NOT the skill's default 20%-time cut. Divergent channels / poor
     r² are FLAGGED, not averaged.

High-value scalar/spectral observables (energy decomposition, current+FFT, dipole,
momentum incl/excl WP, KL metric) come from ``inqview.pipeline`` into
``results/analysis`` via :func:`run_pipeline`; the report embeds those by path.

Geometry-agnostic loading uses the canonical ``inqview.load_vti`` (physical order,
NEVER fftshift — vti-coordinate-mapping rule). Run-tied analysis, so this lives in
``hypotheses/annular_sv/`` (file-placement rule), not in the shared slab library.
"""
from __future__ import annotations

import glob
import os
import re
import sys
from pathlib import Path

import numpy as np

STACK = "/local/data/public/skcb2/tddft/inq-stack/python"
SKILL_SP = "/home/raid/skcb2/skcb2/tddft/.claude/skills/stopping-power-extraction"
for p in (STACK, SKILL_SP):
    if p not in sys.path:
        sys.path.insert(0, p)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import animation  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402

from inqview import load_vti  # noqa: E402  canonical physical-order loader
from inqview.visualisation import style  # noqa: E402
import stopping_power as sp  # noqa: E402  the stopping-power-extraction skill kernels

try:
    style.apply()
except Exception:  # pragma: no cover - theme is cosmetic
    pass

HA_EV = 27.21138625
R_IN, R_OUT = 5.0, 13.0           # tube wall radii (Bohr) — locked geometry
_TPAT = re.compile(r"_t(\d+)\.vti$")


# --------------------------------------------------------------- VTI slice load
def _frame_time(path: str, dt: float) -> float:
    m = _TPAT.search(os.path.basename(path))
    return (int(m.group(1)) * dt) if m else float("nan")


def _slice_stack(vti_dir: str, idx, dt: float):
    """(times[T], slices[T, nz, nx], (x, z)) for the xz mid-y plane at frames idx.

    load_vti.data is (nx, ny, nz) in physical order; the xz slice is data[:, iy, :].T
    → (nz, nx) so imshow(origin='lower', extent=[x0,x1,z0,z1]) shows x horizontal,
    z vertical (the projectile glides up the z axis through the bore)."""
    files = sorted(glob.glob(os.path.join(vti_dir, "*.vti")))
    if not files:
        return None, None, None
    files = [files[k] for k in idx]
    first = load_vti(files[0])
    iy = first.data.shape[1] // 2
    times = np.array([_frame_time(f, dt) for f in files])
    sl = np.empty((len(files), first.data.shape[2], first.data.shape[0]))  # (T,nz,nx)
    for t, f in enumerate(files):
        sl[t] = load_vti(f).data[:, iy, :].T
    return times, sl, (np.asarray(first.x), np.asarray(first.z))


def _wall_lines(ax):
    """Vertical wall-radius markers at x=±R_in, ±R_out (the tube wall)."""
    for xx in (R_IN, -R_IN, R_OUT, -R_OUT):
        ax.axvline(xx, ls="--", lw=0.7, color="0.4")


def _save_gif(slices, times, axes, out_path, *, title, kind, vmax=None, vmin=None,
              fps=10, ion_xz=None, ion_label="projectile"):
    """Animate xz slices; if ion_xz (T,2) of (x,z) is given, overlay the moving
    projectile (its own-run trajectory interpolated to each frame time)."""
    x, z = axes
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    if kind == "density":
        vmax = vmax if vmax is not None else float(np.percentile(slices[len(slices) // 2], 99.5)) or 1e-12
        vmin = vmin if vmin is not None else vmax * 1e-3
        data0 = np.clip(slices[0], vmin, None)
        im = ax.imshow(data0, origin="lower", aspect="auto",
                       extent=[x[0], x[-1], z[0], z[-1]],
                       cmap="viridis", norm=LogNorm(vmin=vmin, vmax=vmax))
        cbar_label = "n (a₀⁻³, log)"
    else:  # diverging difference
        vmax = vmax if vmax is not None else float(np.percentile(np.abs(slices), 99.0)) or 1e-12
        vmin = -vmax
        im = ax.imshow(slices[0], origin="lower", aspect="auto",
                       extent=[x[0], x[-1], z[0], z[-1]],
                       cmap="RdBu_r", vmin=vmin, vmax=vmax)
        cbar_label = "Δn (a₀⁻³)"
    _wall_lines(ax)
    marker = trail = None
    if ion_xz is not None:
        ax.plot(ion_xz[:, 0], ion_xz[:, 1], "-", color="1.0", lw=0.8, alpha=0.45)
        trail, = ax.plot([], [], "-", color="cyan", lw=1.4, alpha=0.8)
        marker, = ax.plot([], [], "o", mfc="cyan", mec="k", mew=0.8, ms=9,
                          label=ion_label)
        ax.legend(loc="upper right", fontsize=7, framealpha=0.5)
    ax.set_xlabel("x (Bohr)"); ax.set_ylabel("z — tube axis (Bohr)")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(cbar_label, fontsize=8)
    ttl = ax.set_title("", fontsize=9)

    def upd(k):
        d = np.clip(slices[k], vmin, None) if kind == "density" else slices[k]
        im.set_data(d)
        ttl.set_text(f"{title} — t = {times[k]:.1f} a.u.")
        arts = [im, ttl]
        if marker is not None:
            marker.set_data([ion_xz[k, 0]], [ion_xz[k, 1]])
            trail.set_data(ion_xz[:k + 1, 0], ion_xz[:k + 1, 1])
            arts += [marker, trail]
        return arts

    an = animation.FuncAnimation(fig, upd, frames=len(slices), blit=False)
    an.save(out_path, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)
    return vmax


def _ion_at(times, ion_traj):
    """Interpolate the projectile (x,z) onto the GIF frame times. ion_traj =
    (t, x, z) arrays; returns (T,2) or None."""
    if ion_traj is None:
        return None
    t_i, x_i, z_i = ion_traj
    return np.column_stack([np.interp(times, t_i, x_i), np.interp(times, t_i, z_i)])


def density_matrix(results_dir, out_dir, label, dt, frames_max=40, fps=10,
                   ion_traj=None, ion_label="projectile"):
    """The {density, delta0, dstep} × {total[, wp, bath]} GIF matrix for one run,
    with the moving projectile overlaid. Returns [(category, kind, path, caption)]."""
    raw = os.path.join(results_dir, "raw", "vti")
    os.makedirs(out_dir, exist_ok=True)
    tot_dir = os.path.join(raw, "density_system")
    wp_dir = os.path.join(raw, "density_wp")
    nfiles = len(glob.glob(os.path.join(tot_dir, "*.vti")))
    if nfiles == 0:
        return []
    idx = list(range(0, nfiles, max(1, nfiles // frames_max)))
    times, tot, axes = _slice_stack(tot_dir, idx, dt)
    ion_xz = _ion_at(times, ion_traj)
    cats = {"total": tot}
    if os.path.isdir(wp_dir) and glob.glob(os.path.join(wp_dir, "*.vti")):
        _, wp, _ = _slice_stack(wp_dir, idx, dt)
        cats["wp"] = wp
        cats["bath"] = tot - wp
    has_wp = "wp" in cats
    base = cats.get("bath", tot)
    dens_vmax = float(np.percentile(base[len(base) // 2], 99.7)) or 1e-12

    KIND = {"density": "density n(x,z,t)", "delta0": "Δn = n(t) − n(0)",
            "dstep": "Δn = n(t+dt) − n(t)"}
    CAT = {"total": "Total system", "wp": "Wavepacket |ψ|²", "bath": "Bath (wall, total − WP)"}
    out = []
    for cat, stack in cats.items():
        for kind in ("density", "delta0", "dstep"):
            if kind == "density":
                series, vmax = stack, (dens_vmax if cat in ("total", "bath") else None)
                ixz = ion_xz
            elif kind == "delta0":
                series, vmax = stack - stack[0][None], None
                ixz = ion_xz
            else:
                series, vmax = np.diff(stack, axis=0), None
                ixz = ion_xz[1:] if ion_xz is not None else None
            klabel = KIND[kind]
            if cat == "total" and kind == "delta0" and not has_wp:
                klabel = "Δn = n(t) − n(0)  (induced wake)"
            tt = times[1:] if kind == "dstep" else times
            f = os.path.join(out_dir, f"{label}_{cat}_{kind}.gif")
            _save_gif(series, tt, axes, f,
                      title=f"{label} · {CAT[cat]} · {klabel}",
                      kind=("density" if kind == "density" else "diff"),
                      vmax=vmax, fps=fps, ion_xz=ixz, ion_label=ion_label)
            out.append((cat, kind, f, f"{CAT[cat]} — {klabel}"))
    return out


# ------------------------------------------------------------------- z-t carpets
def carpets(results_dir, out_dir, label, dt, frames_max=160, ion_traj=None,
            ion_label="projectile"):
    """z–t carpets (total n(z,t), Δ-vs-0, per-step Δ) with the projectile z(t)."""
    tot_dir = os.path.join(results_dir, "raw", "vti", "density_system")
    files = sorted(glob.glob(os.path.join(tot_dir, "*.vti")))
    if not files:
        return []
    if len(files) > frames_max:
        files = files[:: len(files) // frames_max + 1]
    first = load_vti(files[0])
    z = np.asarray(first.z)
    steps = np.array([int(_TPAT.search(os.path.basename(f)).group(1)) for f in files])
    nzt = np.array([load_vti(f).data.sum(axis=(0, 1)) for f in files])  # ∫dx dy n → n(z)
    t = steps * dt
    os.makedirs(out_dir, exist_ok=True)
    z_ion = (np.interp(t, ion_traj[0], ion_traj[2]) if ion_traj is not None else None)
    specs = [("total", nzt, "inferno", r"$\int n\,dx\,dy$ (e/Bohr)", False, "total n(z,t)"),
             ("delta0", nzt - nzt[0], "RdBu_r", r"$n(t)-n(0)$", True, "Δn = n(t) − n(0)"),
             ("dstep", np.vstack([np.zeros_like(nzt[0]), np.diff(nzt, axis=0)]),
              "RdBu_r", r"$n(t+dt)-n(t)$", True, "Δn = n(t+dt) − n(t)")]
    out = []
    for tag, M, cmap, lab, signed, cap in specs:
        a = float(np.percentile(np.abs(M), 99.5)) or 1e-12
        norm = Normalize(-a, a) if signed else Normalize(0.0, a)
        fig, ax = plt.subplots(figsize=(6.6, 4.2))
        im = ax.pcolormesh(t, z, M.T, cmap=cmap, shading="auto", norm=norm)
        if z_ion is not None:
            ax.plot(t, z_ion, "-", color="cyan", lw=1.4, label=f"{ion_label} z(t)")
            ax.legend(loc="upper left", fontsize=7, framealpha=0.5)
        ax.set_xlabel("time (a.u.)"); ax.set_ylabel("z — tube axis (Bohr)")
        ax.set_title(f"{label}: {cap}", fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label=lab)
        f = os.path.join(out_dir, f"{label}_carpet_{tag}.png")
        fig.tight_layout(); fig.savefig(f, dpi=150); plt.close(fig)
        out.append((cap, f))
    return out


# --------------------------------------------------- stopping (skill-compliant)
def _find(results_dir, name):
    return next(Path(results_dir).glob(f"**/{name}"), None)


def stopping_analysis(results_dir, out_dir, label, rs, v0, vfrac=0.85):
    """Skill-compliant classical stopping (stopping-power-extraction; Correa 2018).

    Continuous-traversal geometry (periodic tube) → Method A: PRIMARY = free-intercept
    slope of the electronic deposit dE_total(x). The window is the EARLY v≥vfrac·v0
    segment (light free-Ehrenfest electron decelerates; S is read AT v0), NOT the
    skill's default 20%-time cut. Guards: N(t)≈const; dE_total ≈ −dKE_ion. Sanity:
    kinetic channel −dKE_ion/dx. Flags channel divergence (>10%) and poor r²."""
    import pandas as pd
    obs, trk, num = (_find(results_dir, "observables.csv"),
                     _find(results_dir, "electron_track.csv"),
                     _find(results_dir, "electron_number.csv"))
    if not (obs and trk):
        return None
    O = pd.read_csv(obs); T = pd.read_csv(trk).drop_duplicates("step")
    tE = O["time_au"].to_numpy()
    E = (O["energy_total"] - O["energy_total"].iloc[0]).to_numpy()  # electronic deposit
    z0 = T["z"].to_numpy()[0]
    s_tr = np.abs(T["z"].to_numpy() - z0)
    ke_tr = T["ke_ion_ha"].to_numpy()
    if len(tE) < 10 or not np.all(np.isfinite(E)):
        return None
    x = np.interp(tE, T["time_au"], s_tr)
    vz = np.interp(tE, T["time_au"], T["vz"])
    # early near-constant-velocity window (widen if sparse)
    used_vf = vfrac
    for vf in (vfrac, 0.70, 0.50):
        m = vz >= vf * v0
        m[:max(2, int(0.03 * len(x)))] = False
        used_vf = vf
        if m.sum() >= 20:
            break
    if m.sum() < 6:
        return None
    x0, xT = float(x[m].min()), float(x[m].max())
    prim = sp.free_fit(x, E, x0, xT)                 # PRIMARY: dE_total slope
    kin = sp.kinetic_channel(s_tr, ke_tr, x, x0, xT)  # sanity: −dKE_ion/dx
    nguard = (sp.conservation_guard(pd.read_csv(num)["N_total"].to_numpy())
              if num else {"ok": None, "drained_frac": float("nan")})
    if prim is None or kin is None:
        return None
    S, Sk = prim["S"], kin["S"]
    ratio = S / Sk if Sk else float("nan")
    # energy-conservation over the window: deposit vs ion KE loss
    dE_win = float(np.interp(xT, x, E) - np.interp(x0, x, E))
    dKE_win = float(np.interp(x0, s_tr, ke_tr) - np.interp(xT, s_tr, ke_tr))
    econs = dE_win / dKE_win if dKE_win else float("nan")
    flags = []
    if abs(ratio - 1) > 0.10:
        flags.append(f"channel divergence {abs(ratio-1)*100:.0f}% (dE_total vs −dKE_ion)")
    if prim["r2"] < 0.80:
        flags.append(f"poor linear fit r²={prim['r2']:.2f}")
    if nguard.get("ok") is False:
        flags.append(f"N drained {nguard['drained_frac']*100:.1f}%")

    # plot: vz(t) window | dE_total(s)+fit | both deposit channels overlaid
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.2))
    a1.plot(tE, vz, "C0-", lw=1.4)
    a1.axhline(used_vf * v0, ls="--", color="0.5", lw=0.8,
               label=f"{used_vf:.2f}·v0 = {used_vf*v0:.3f}")
    a1.axhspan(used_vf * v0, max(vz) * 1.02, color="C0", alpha=0.08)
    a1.set_xlabel("time (a.u.)"); a1.set_ylabel(r"$v_z$ (a.u.)")
    a1.set_title(f"{label}: decelerates (v0={v0})", fontsize=9)
    a1.legend(fontsize=8, frameon=False); a1.grid(alpha=.25)

    a2.plot(x, E * HA_EV, ".", color="0.75", ms=3, label="all steps")
    a2.plot(x[m], E[m] * HA_EV, "C0.", ms=4, label="early window")
    xs = np.linspace(x0, xT, 50)
    a2.plot(xs, (prim["S"] * xs + prim["E0"]) * HA_EV, "k-", lw=1.4,
            label=f"S = {S:.4f} Ha/Bohr (r²={prim['r2']:.2f})")
    a2.set_xlabel("path s = |z − z₀| (Bohr)")
    a2.set_ylabel(r"$\Delta E_\mathrm{total}$ electronic deposit (eV)")
    a2.set_title("PRIMARY — defined method: S = dE_total/dx", fontsize=9, weight="bold")
    a2.legend(fontsize=8, frameon=False); a2.grid(alpha=.25)

    dKE_loss = ke_tr[0] - np.interp(x, s_tr, ke_tr)
    a3.plot(x[m], E[m] * HA_EV, "C0.", ms=4, label=f"PRIMARY ΔE_total → S={S:.4f}")
    a3.plot(x[m], dKE_loss[m] * HA_EV, "C3.", ms=4, label=f"sanity −ΔKE_ion → S={Sk:.4f}")
    a3.set_xlabel("path s (Bohr)"); a3.set_ylabel("deposited energy (eV)")
    a3.set_title(f"SANITY CHECK: KE method vs primary (ratio {ratio:.2f})", fontsize=9)
    a3.legend(fontsize=8, frameon=False); a3.grid(alpha=.25)

    if flags:
        fig.text(0.5, 0.005, "⚠ " + "; ".join(flags), ha="center", color="C3",
                 fontsize=9)
    os.makedirs(out_dir, exist_ok=True)
    f = os.path.join(out_dir, f"{label}_stopping.png")
    fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(f, dpi=160); plt.close(fig)
    return dict(path=f, S=float(S), S_err=float(prim["se"]), r2=float(prim["r2"]),
                S_kinetic=float(Sk), ratio=float(ratio),
                N_drained=float(nguard.get("drained_frac", float("nan"))),
                econs_ratio=float(econs), v_mean=float(vz[m].mean()),
                npts=int(m.sum()), window=[x0, xT], vfrac=float(used_vf),
                flags=flags)


# ------------------------------------------------------------------ the pipeline
def run_pipeline(results_dir, run_name):
    """Run the high-value-observable pipeline phases; figures land in
    results/analysis. Returns the PipelineResult (never raises — best effort)."""
    from inqview.pipeline import runner
    phases = ["observables", "momentum", "kl_divergence"]
    try:
        return runner.run(results_dir, run_name=run_name, phases=phases, rebuild=True)
    except Exception as e:  # noqa: BLE001
        print(f"[warn] pipeline on {run_name}: {e}")
        return None


def collect_pipeline_figs(results_dir):
    """Group the pipeline figures under results/analysis into report buckets.

    NOTE: the pipeline's own FFT/spectrum PNGs (`fft_*`, `spectra/`) are
    DELIBERATELY EXCLUDED — every FFT-driven observable is re-done through the
    `fourier-analysis` skill (`fft_panels`, the audited 6-stage panel). Only
    TIME-DOMAIN pipeline figures are embedded here."""
    an = Path(results_dir) / "analysis"
    if not an.exists():
        return {}
    def has(*keys):
        out = []
        for p in sorted(an.rglob("*.png")) + sorted(an.rglob("*.gif")):
            if "fft_" in p.name or "spectrum_" in p.name or "/spectra/" in str(p):
                continue                      # FFT handled by the fourier-analysis skill
            if any(k in p.name for k in keys):
                out.append((p.name, str(p)))
        return out
    groups = {
        "Energy decomposition": has("all_energies_vs_time", "energy_total_vs_time",
                                    "energy_kinetic", "energy_hartree", "energy_xc"),
        "Induced current & dipole (time domain)": has("current_components",
                                                      "dipole_components"),
        "Momentum (1D n(k), incl/excl WP)": has("momentum_heatmap", "momentum_distribution"),
        "KL metric (WP momentum drift)": has("kl_divergence_vs_t"),
    }
    return {k: v for k, v in groups.items() if v}


def fft_panels(results_dir, out_dir, label, rs):
    """Well-made FFT via the `fourier-analysis` skill standard — the audited 3×2
    6-stage `fft_pipeline_panel` (Hann + mean baseline + ×4 zero-pad + coherent
    gain + ANGULAR ħω axis; detrend overlaid for audit). Replaces the pipeline's
    raw FFT plots. Plasmon band ħω_p=√(3/r_s³)·27.211 eV shaded; Δω=2π/τ annotated."""
    import pandas as pd
    from inqview.visualisation.fourier_panel import fft_pipeline_panel
    obs = _find(results_dir, "observables.csv")
    if obs is None:
        return []
    O = pd.read_csv(obs).drop_duplicates("step").sort_values("time_au")
    if len(O) < 8:
        return []
    t = O["time_au"].to_numpy()
    wp = (3.0 / rs ** 3) ** 0.5 * HA_EV                # plasmon energy (eV)
    band = (0.4 * wp, 2.0 * wp)
    fmax = max(12.0, 3.0 * wp)
    tau = float(t[-1] - t[0])
    dw = 2 * np.pi / tau * HA_EV
    os.makedirs(out_dir, exist_ok=True)
    out = []
    for col, desc in [("current_z", "induced axial wall current (hydrovoltaic signal)"),
                      ("energy_total", "electronic total energy")]:
        if col not in O.columns or float(np.nanstd(O[col].to_numpy())) == 0.0:
            continue
        fig = fft_pipeline_panel(t, O[col].to_numpy(), label=col, peak_band=band,
                                 fmax=fmax,
                                 title=f"{label} — FFT pipeline: {col}  (ħω_p≈{wp:.2f} eV)")
        f = os.path.join(out_dir, f"{label}_fft_{col}.png")
        fig.savefig(f, dpi=140)
        plt.close(fig)
        out.append((f"FFT pipeline — {col} ({desc}). Plasmon band ħω_p≈{wp:.2f} eV "
                    f"shaded; resolution Δω=2π/τ≈{dw:.1f} eV (τ={tau:.0f} a.u., coarse — "
                    f"informational, not a gate). fourier-analysis skill standard.", f))
    return out


# -------------------------------------------------------- projectile trajectory
def ion_trajectory(results_dir, rtype):
    """(t, x, z, label) for the moving projectile — classical ion from the track,
    WP centroid from wp_real_space_stats. None if unavailable."""
    import pandas as pd
    if rtype == "wp":
        p = _find(results_dir, "wp_real_space_stats.csv")
        if p is None:
            return None
        R = pd.read_csv(p, comment="#")
        return (R["time_au"].to_numpy(), R["x_mean"].to_numpy(),
                R["z_mean"].to_numpy(), "WP centroid")
    p = _find(results_dir, "electron_track.csv")
    if p is None:
        return None
    T = pd.read_csv(p).drop_duplicates("step")
    return (T["time_au"].to_numpy(), T["x"].to_numpy(), T["z"].to_numpy(),
            "classical ion")


# --------------------------------------------------------------------- assemble
def generate(label, rs, v0, run_dir, results_dir, figroot, rtype):
    """Generate the full per-run battery; return {group: [(caption, path), ...]}."""
    dt = 0.02
    out_dir = os.path.join(figroot, label)
    os.makedirs(out_dir, exist_ok=True)
    traj = ion_trajectory(results_dir, rtype)
    ion_traj = (traj[0], traj[1], traj[2]) if traj else None
    ion_label = traj[3] if traj else "projectile"
    print(f"[{label}] pipeline …", flush=True)
    run_pipeline(results_dir, label)
    print(f"[{label}] density matrix (+ion overlay) …", flush=True)
    matrix = density_matrix(results_dir, out_dir, label, dt, ion_traj=ion_traj,
                            ion_label=ion_label)
    print(f"[{label}] carpets …", flush=True)
    cps = carpets(results_dir, out_dir, label, dt, ion_traj=ion_traj,
                  ion_label=ion_label)
    sp_res = None
    if rtype == "classical":
        print(f"[{label}] stopping (skill) …", flush=True)
        sp_res = stopping_analysis(results_dir, out_dir, label, rs, v0)
    print(f"[{label}] FFT panels (fourier-analysis skill) …", flush=True)
    fft = fft_panels(results_dir, out_dir, label, rs)
    pipe = collect_pipeline_figs(results_dir)
    return dict(matrix=matrix, carpets=cps, stopping=sp_res, fft=fft, pipeline=pipe)
