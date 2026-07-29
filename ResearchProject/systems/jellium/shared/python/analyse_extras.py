"""analyse_extras.py — custom post-pipeline analyses shared by every
per-run analyse.py in the jellium project.

These plots are *not* in inqview's pipeline. Each analyse.py imports the
helpers below and invokes them on the run's results/ tree. Keeping them
here avoids drift between the per-run analyse.py templates.

Plots produced (all 2D / 1D, no 3D rendering):

  classical_only:
    classical_force_fixed.png        F_z = m·dv/dt  in eV/Bohr, vs z
    delta_E_total_vs_time.png        bath energy gain vs time (eV) +
                                     window annotation
    delta_E_total_vs_z.png           bath energy gain vs projectile z,
                                     with windowed linear fit per plan
                                     §6.1 (Δz ∈ [3, 28] Bohr) → S(v)
    running_slope_vs_z.png           box-deficit diagnostic per plan
                                     §7 deliverable 4 (rolling 0.4 a.u.
                                     window slope dE_bath/dt vs z)

  both:
    overlap_heatmap_t_end.png        |O_ij(t_end)|² heatmap, fixed [0,1]
    effective_gs_occupations_t_end.png   n_i^GS(t_end) = Σ_j f_j|O_ij|²

All plots follow docs/visualisation-instructions-v1.md:
  - run name in title
  - axes labelled with units
  - 3 significant figures on annotations
  - legends where multiple series

References:
  - F_z = dv/dt fix: classical electron_track.csv emits fz=0.0 placeholders
    because INQ's StepContext doesn't expose ehrenfest forces. We recover
    F_z from the velocity column.
  - Plan §6.1 windowed stopping power: fit delta_E_bath vs time over
    t ∈ [3/v, 28/v] a.u. (equivalent to Δz ∈ [3, 28] Bohr), report
    slope/v in eV/Bohr with regression-covariance standard error.
  - Effective GS-projected occupation:
      n_i^GS(t) = Σ_j f_j(0) |⟨ψ_i^GS|ψ_j(t)⟩|² = Σ_j f_j(0) |O_ij(t)|²
    With f_j(0) = (initial occupations of the evolved orbitals; for
    closed-shell N-electron systems all occupied f_j=2, unoccupied f_j=0).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HA_TO_EV = 27.211386245988


# ---------------------------------------------------------------------------
#  Classical: force from dv/dt (in eV/Bohr)
# ---------------------------------------------------------------------------
def plot_classical_force_fixed(results_dir: Path, run_name: str) -> Path | None:
    """F_z(z) = m·dv/dt recovered from electron_track.csv (fz column is 0).

    Output: results/analysis/observables/classical_force_fixed.png
    """
    track_csv = results_dir / "raw" / "observables" / "electron_track.csv"
    if not track_csv.exists():
        return None
    df = pd.read_csv(track_csv).sort_values("step").reset_index(drop=True)
    t_au = df["time_au"].to_numpy()
    z    = df["z"].to_numpy()
    vz   = df["vz"].to_numpy()

    # Central differences: F_z = m·dv/dt with m=1 a.u.
    dt_au = np.gradient(t_au)
    fz_au = np.gradient(vz, t_au)                       # Ha/Bohr
    fz_eV_per_Bohr = fz_au * HA_TO_EV

    out_dir = results_dir / "analysis" / "observables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "classical_force_fixed.png"

    mean_fz = float(np.median(fz_eV_per_Bohr))          # robust central value

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(z, fz_eV_per_Bohr, "C0-", lw=1.3)
    ax.axhline(mean_fz, color="C3", lw=1.0, ls="--",
               label=f"median F_z = {mean_fz:.3g} eV/Bohr")
    ax.axhline(0.0, color="0.5", lw=0.5)
    ax.set_xlabel("Projectile z / Bohr")
    ax.set_ylabel("F_z on projectile / (eV / Bohr)")
    ax.set_title(f"{run_name}: stopping force F_z = m·dv/dt vs z")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
#  Classical: ΔE_total vs time, and ΔE_total vs projectile z with windowed
#  slope fit (plan §6.1 rule: Δz ∈ [3, 28] Bohr)
# ---------------------------------------------------------------------------
def _read_v_initial(results_dir: Path) -> float | None:
    """Parse classical projectile's initial v_z (a.u.) from run_summary.txt.

    Two formats supported:
      velocity_atu    = 0 0 2.7110633...   (run.cpp emits 3-tuple)
      projectile_v    = 0 0 2.7110633...   (stub format)
    """
    rs = results_dir / "run_summary.txt"
    if not rs.exists():
        return None
    import re
    text = rs.read_text()
    for key in ("velocity_atu", "projectile_v"):
        m = re.search(rf"^\s*{key}\s*=\s*\S+\s+\S+\s+(\S+)", text,
                       flags=re.MULTILINE)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                continue
    return None


def plot_classical_excess_energy(results_dir: Path, run_name: str
                                  ) -> tuple[Path | None, Path | None,
                                              float, float]:
    """Plots ΔE_total of the electronic system vs (a) time, (b) projectile z,
    with a windowed linear fit per plan §6.1.

    Window: Δz ∈ [3, 28] Bohr, equivalent to t ∈ [3/v, 28/v] a.u.,
    where v is parsed from run_summary.txt (PROJ_VEL_Z). Outside this
    window the slope is contaminated by (early) bath-relaxation transient
    or (late) periodic-image-wake overlap.

    Outputs:
      delta_E_total_vs_time.png   (with vertical window markers)
      delta_E_total_vs_z.png      (with vertical window markers at
                                   z=3, z=28 and the windowed linear
                                   fit annotated with slope ± SE)

    Returns (path_time, path_z, S_eV_per_Bohr, S_uncertainty_eV_per_Bohr).
    S is the stopping power = slope of ΔE_bath vs z (or equivalently
    slope_vs_time / v_initial). For a perfectly-uniform-velocity classical
    projectile under .ehrenfest(), these two routes agree.
    """
    obs_csv   = results_dir / "raw" / "observables" / "observables.csv"
    track_csv = results_dir / "raw" / "observables" / "electron_track.csv"
    if not obs_csv.exists() or not track_csv.exists():
        return (None, None, float("nan"), float("nan"))

    obs   = pd.read_csv(obs_csv)
    track = pd.read_csv(track_csv).sort_values("step").reset_index(drop=True)

    # Inner-join on step (observables every WRITE_EVERY; track every step).
    merged = pd.merge(obs[["step", "time_au", "energy_total"]],
                      track[["step", "z"]],
                      on="step", how="inner")
    t_au  = merged["time_au"].to_numpy()
    E     = merged["energy_total"].to_numpy()
    z     = merged["z"].to_numpy()
    dE_eV = (E - E[0]) * HA_TO_EV

    # Δz from launch (launch z is z[0]; for our setup z[0] = -10 Bohr).
    z0       = float(z[0])
    delta_z  = z - z0

    # Plan §6.1 windowing: Δz ∈ [3, 28] Bohr. Fall back to 20-80% in time
    # if v_initial is not parseable or the window collapses.
    v_init = _read_v_initial(results_dir)
    fit_label = ""
    if v_init is not None and v_init > 1e-6:
        t_start = 3.0 / v_init
        t_end   = 28.0 / v_init
        mask    = (t_au >= t_start) & (t_au <= t_end)
        if mask.sum() >= 5:
            z_fit   = z[mask]
            dE_fit  = dE_eV[mask]
            t_fit   = t_au[mask]
            fit_label = (f"window Δz ∈ [3, 28] Bohr  "
                         f"(t ∈ [{t_start:.3g}, {t_end:.3g}] a.u.)")
        else:
            # Window collapsed; fall back to mid-60%.
            n = len(z); i0 = int(0.2 * n); i1 = int(0.8 * n)
            z_fit, dE_fit, t_fit = z[i0:i1], dE_eV[i0:i1], t_au[i0:i1]
            fit_label = "fallback window (mid 60% of trajectory)"
    else:
        n = len(z); i0 = int(0.2 * n); i1 = int(0.8 * n)
        z_fit, dE_fit, t_fit = z[i0:i1], dE_eV[i0:i1], t_au[i0:i1]
        fit_label = "fallback window (no v_initial; mid 60% of trajectory)"

    # Linear regression with covariance (for uncertainty estimate).
    # polyfit gives the slope and intercept; we recompute residuals to
    # get the standard error on the slope.
    coeffs, cov = np.polyfit(z_fit, dE_fit, 1, cov=True)
    slope_eVperBohr   = float(coeffs[0])
    intercept         = float(coeffs[1])
    se_slope          = float(np.sqrt(cov[0, 0]))

    out_dir = results_dir / "analysis" / "observables"
    out_dir.mkdir(parents=True, exist_ok=True)

    # (a) ΔE vs time
    out_t = out_dir / "delta_E_total_vs_time.png"
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(t_au, dE_eV, "C0-", lw=1.5)
    ax.axhline(0.0, color="0.5", lw=0.5)
    if v_init is not None and v_init > 1e-6:
        ax.axvline(3.0 / v_init,  color="C3", ls=":", lw=1.0,
                   label=f"Δz=3 Bohr (t={3.0/v_init:.3g} a.u.)")
        ax.axvline(28.0 / v_init, color="C3", ls=":", lw=1.0,
                   label=f"Δz=28 Bohr (t={28.0/v_init:.3g} a.u.)")
        ax.legend(loc="best", fontsize=9)
    ax.set_xlabel("Time / a.u.")
    ax.set_ylabel("Δ E_total of bath / eV")
    ax.set_title(f"{run_name}: bath excess energy vs time")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_t, dpi=150)
    plt.close(fig)

    # (b) ΔE vs projectile z, with windowed linear fit
    out_z = out_dir / "delta_E_total_vs_z.png"
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(z, dE_eV, "C0-", lw=1.5, label="bath energy gain")
    z_line = np.linspace(z_fit.min(), z_fit.max(), 100)
    ax.plot(z_line, slope_eVperBohr * z_line + intercept, "C3--", lw=1.6,
            label=(f"linear fit  S = {slope_eVperBohr:.3g} ± "
                   f"{se_slope:.2g} eV/Bohr"))
    # Vertical window markers at z=z0+3 and z=z0+28.
    ax.axvline(z0 + 3.0,  color="C2", ls=":", lw=1.0,
               label=f"Δz=3 Bohr (z={z0+3.0:.1f})")
    ax.axvline(z0 + 28.0, color="C2", ls=":", lw=1.0,
               label=f"Δz=28 Bohr (z={z0+28.0:.1f})")
    ax.axhline(0.0, color="0.5", lw=0.5)
    ax.set_xlabel("Projectile z / Bohr")
    ax.set_ylabel("Δ E_total of bath / eV")
    ax.set_title(f"{run_name}: bath excess energy vs projectile z  "
                 f"[{fit_label}]")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_z, dpi=150)
    plt.close(fig)

    return (out_t, out_z, slope_eVperBohr, se_slope)


# ---------------------------------------------------------------------------
#  Classical: running slope dE_bath/dz smoothed over 0.4 a.u. windows
#             (plan §7 deliverable 4: box-deficit diagnostic)
# ---------------------------------------------------------------------------
def plot_classical_running_slope(results_dir: Path, run_name: str
                                  ) -> Path | None:
    """Box-deficit diagnostic: rolling slope of dE_bath/dt vs z.

    For each time t, compute the local slope of energy_total over a
    ±0.4 a.u. window centred on t (so a 0.8 a.u. window total). The
    running slope is dE/dt at that time, in Ha/a.u. Plot it as a function
    of projectile z = z(t) and in eV/Bohr (divide by v_initial to convert
    from dE/dt to dE/dz).

    The plot reveals when periodic-image-wake overlap starts contaminating
    the stopping-power estimate (the slope drops at late z), and any
    early-time transient at small z.
    """
    obs_csv   = results_dir / "raw" / "observables" / "observables.csv"
    track_csv = results_dir / "raw" / "observables" / "electron_track.csv"
    if not obs_csv.exists() or not track_csv.exists():
        return None

    obs   = pd.read_csv(obs_csv)
    track = pd.read_csv(track_csv).sort_values("step").reset_index(drop=True)
    merged = pd.merge(obs[["step", "time_au", "energy_total"]],
                      track[["step", "z"]],
                      on="step", how="inner")
    t_au = merged["time_au"].to_numpy()
    E    = merged["energy_total"].to_numpy()
    z    = merged["z"].to_numpy()

    v_init = _read_v_initial(results_dir)
    if v_init is None or v_init < 1e-6:
        return None

    # Rolling slope over half-width = 0.4 a.u.
    half_window_au = 0.4
    n = len(t_au)
    slope_eVperBohr = np.full(n, np.nan)
    for i in range(n):
        t_i = t_au[i]
        mask = (t_au >= t_i - half_window_au) & (t_au <= t_i + half_window_au)
        if mask.sum() < 3:
            continue
        coeffs = np.polyfit(t_au[mask], E[mask], 1)
        # slope in Ha / a.u. = dE/dt; stopping S = (dE/dt) / v
        slope_eVperBohr[i] = float(coeffs[0]) * HA_TO_EV / v_init

    out_dir = results_dir / "analysis" / "observables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "running_slope_vs_z.png"

    z0 = float(z[0])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(z, slope_eVperBohr, "C0-", lw=1.3,
            label="running slope (±0.4 a.u. window)")
    ax.axvline(z0 + 3.0,  color="C2", ls=":", lw=1.0,
               label=f"Δz=3 Bohr (clean window start)")
    ax.axvline(z0 + 28.0, color="C2", ls=":", lw=1.0,
               label=f"Δz=28 Bohr (clean window end)")
    ax.axhline(0.0, color="0.5", lw=0.5)
    ax.set_xlabel("Projectile z / Bohr")
    ax.set_ylabel("Local dE_bath/dz / (eV/Bohr)")
    ax.set_title(f"{run_name}: box-deficit diagnostic — "
                 f"running stopping power vs z")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
#  Both: |O_ij(t_end)|² heatmap and effective GS-projected occupations
# ---------------------------------------------------------------------------
def _read_overlap_full_at(results_dir: Path, step_query: int | None = None
                          ) -> tuple[np.ndarray, int, float] | None:
    """Read the overlap_full row at the requested step (or last if None).

    Returns (O_squared, step, time_au) or None if no data.

    Format of raw/observables/overlap_full/:
      index.csv     step,time_au,file
      overlap_NNNNNN.csv   header line, then n_ref rows of n_evolved cols
                            each entry = |⟨ψ_i^GS | ψ_j(t)⟩|² already
                            (per inqkit/observables/orbital_overlap.hpp).
    """
    full_dir = results_dir / "raw" / "observables" / "overlap_full"
    idx = full_dir / "index.csv"
    if not idx.exists():
        return None
    df = pd.read_csv(idx)
    if df.empty:
        return None
    if step_query is None:
        row = df.iloc[-1]
    else:
        row = df.iloc[(df["step"] - step_query).abs().argmin()]
    csv = full_dir / row["file"]
    if not csv.exists():
        return None
    arr = pd.read_csv(csv, comment="#", header=None).to_numpy()
    return arr, int(row["step"]), float(row["time_au"])


def plot_overlap_heatmap_at_end(results_dir: Path, run_name: str) -> Path | None:
    """Plot |O_ij(t_end)|² as a heatmap with fixed [0, 1] colour scale.

    rows    = GS KS orbital index i
    columns = evolved KS orbital index j
    colour  = overlap magnitude |⟨ψ_i^GS | ψ_j(t_end)⟩|²
    """
    data = _read_overlap_full_at(results_dir, step_query=None)
    if data is None:
        return None
    O2, step, t_au = data

    out_dir = results_dir / "analysis" / "overlap"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "overlap_heatmap_t_end.png"

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    im = ax.imshow(O2, origin="lower", aspect="auto",
                   vmin=0.0, vmax=1.0, cmap="viridis")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(r"$|\langle\psi_i^{GS}|\psi_j(t)\rangle|^2$")
    ax.set_xlabel("Evolved KS orbital index j")
    ax.set_ylabel("Ground-state KS orbital index i")
    ax.set_title(f"{run_name}: overlap matrix, step {step}, "
                 f"t = {t_au:.3g} a.u.")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return out_png


def plot_effective_gs_occupations_at_end(
    results_dir: Path,
    run_name: str,
    occupations_t0: np.ndarray | None = None,
) -> Path | None:
    """Plot n_i^GS(t_end) = Σ_j f_j(0) |O_ij(t_end)|².

    If occupations_t0 not supplied, read from
    raw/observables/eigenvalues/occupations.csv (the GS occupations).
    For a closed-shell jellium run these are all 2 for the occupied
    shells, 0 for the unoccupied / WP-injected state.
    """
    data = _read_overlap_full_at(results_dir, step_query=None)
    if data is None:
        return None
    O2, step, t_au = data  # shape (n_ref, n_evolved)
    n_ref, n_evolved = O2.shape

    # Default occupations: read GS occupations CSV.
    if occupations_t0 is None:
        occ_csv = (results_dir / "raw" / "observables" / "eigenvalues"
                   / "occupations.csv")
        if not occ_csv.exists():
            return None
        occ_df = pd.read_csv(occ_csv).sort_values("state_index")
        f_evolved = occ_df["occupation"].to_numpy()
        if f_evolved.size < n_evolved:
            # Pad with zeros (e.g. WP slot doesn't appear in GS occupations).
            pad = np.zeros(n_evolved - f_evolved.size)
            f_evolved = np.concatenate([f_evolved, pad])
        elif f_evolved.size > n_evolved:
            f_evolved = f_evolved[:n_evolved]
    else:
        f_evolved = occupations_t0[:n_evolved]

    # n_i^GS = Σ_j f_j(0) |O_ij|²
    n_gs = O2 @ f_evolved
    # Reference: GS occupations (i = j ordering aligned for n_ref ≤ n_evolved)
    n_gs_ref = f_evolved[:n_ref]

    # Conservation check.
    n_total = float(n_gs.sum())
    n_ref_total = float(n_gs_ref.sum())

    out_dir = results_dir / "analysis" / "observables" / "gs_projected_occupations"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "effective_gs_occupations_t_end.png"

    idx = np.arange(n_ref)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    width = 0.4
    ax.bar(idx - width/2, n_gs_ref, width=width, color="0.75",
           label=f"GS occupations f_i  (Σ = {n_ref_total:.3g} e⁻)")
    ax.bar(idx + width/2, n_gs,     width=width, color="C0",
           label=fr"n_i^GS(t_end) = Σ_j f_j |O_ij|²  (Σ = {n_total:.3g} e⁻)")
    ax.set_xlabel("Ground-state KS orbital index i")
    ax.set_ylabel("Effective occupation")
    ax.set_title(f"{run_name}: GS-projected occupations at "
                 f"t = {t_au:.3g} a.u.  (step {step})")
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

    # Also write a CSV so the user can table these numbers.
    csv_out = out_dir / "effective_gs_occupations_t_end.csv"
    pd.DataFrame({
        "gs_orbital_index": idx,
        "f_i_gs":           n_gs_ref,
        "n_i_gs_t_end":     n_gs,
        "delta":            n_gs - n_gs_ref,
    }).to_csv(csv_out, index=False)
    return out_png


# ---------------------------------------------------------------------------
#  Driver: call from analyse.py main(), after pipeline.run().
# ---------------------------------------------------------------------------
def run_extras(results_dir: Path, run_name: str, *, is_classical: bool
               ) -> dict[str, str]:
    """Run all custom analyses and return a {label: status} log."""
    log: dict[str, str] = {}

    # Heatmap (both classical and WP).
    p = plot_overlap_heatmap_at_end(results_dir, run_name)
    log["overlap_heatmap_t_end"] = "[ok] " + str(p) if p else "[skip] no overlap_full data"

    # Effective GS-projected occupations at t_end.
    p = plot_effective_gs_occupations_at_end(results_dir, run_name)
    log["effective_gs_occupations_t_end"] = "[ok] " + str(p) if p else "[skip] no overlap_full data"

    if is_classical:
        p = plot_classical_force_fixed(results_dir, run_name)
        log["classical_force_fixed"] = "[ok] " + str(p) if p else "[skip] no electron_track.csv"

        pt, pz, slope, se = plot_classical_excess_energy(results_dir, run_name)
        if pt and pz:
            log["delta_E_total_vs_time"] = "[ok] " + str(pt)
            log["delta_E_total_vs_z"]    = (
                f"[ok] {pz}  S = {slope:.4g} ± {se:.2g} eV/Bohr (plan §6.1 windowed)"
            )
        else:
            log["delta_E_total_vs_time"] = "[skip] missing inputs"
            log["delta_E_total_vs_z"]    = "[skip] missing inputs"

        # Box-deficit diagnostic (plan §7 deliverable 4).
        p = plot_classical_running_slope(results_dir, run_name)
        log["running_slope_vs_z"] = "[ok] " + str(p) if p else "[skip] no electron_track.csv or v_initial"

    return log
