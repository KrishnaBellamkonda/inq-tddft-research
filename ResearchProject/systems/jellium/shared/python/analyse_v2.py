"""analyse_v2.py — reusable post-processing plots for jellium TDDFT runs.

Second-generation analysis plots covering energy decomposition, KS
eigenenergy evolution, GS-basis decomposition, and overlap heatmaps.
Designed to be imported from per-run analyse.py scripts.

Each public function takes (results_dir, run_name, ...) and returns
the output PNG path (or None if data is missing).

Plots produced:
  1. plot_energy_decomposition        6-panel energy components vs Dz
  2. plot_energy_bookkeeping_bar      cross-run bar chart at t_IFW
  3. plot_ks_eigenenergy_evolution    per-orbital <H>(t) evolution
  4. plot_gs_basis_decomposition      Delta n_i^GS bar chart at t_end
  5. plot_overlap_heatmap             |O_ij(t_end)|^2 heatmap (log)
  6. plot_overlap_heatmap_diff        WP - classical overlap difference

References:
  - Energy decomposition pattern: build_figures.py (meeting 2026-05-14)
  - KS eigenenergy pattern: case_study_E100eV.py figure 10
  - GS-basis decomposition pattern: case_study_E100eV.py figure 11
  - Overlap heatmap pattern: analyse_extras.py
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HA_TO_EV = 27.211386245988


# =============================================================================
# Internal helpers
# =============================================================================

def _parse_run_summary(results_dir: Path) -> dict[str, str]:
    """Parse run_summary.txt into a flat key=value dict."""
    rs = results_dir / "run_summary.txt"
    if not rs.exists():
        return {}
    out: dict[str, str] = {}
    for line in rs.read_text().splitlines():
        m = re.match(r"^\s*(\w[\w_]*)\s*=\s*(.+)$", line)
        if m:
            out[m.group(1).strip()] = m.group(2).strip()
    return out


def _get_wp_state_index(results_dir: Path) -> int | None:
    """Read wp_state_index from run_summary.txt."""
    params = _parse_run_summary(results_dir)
    val = params.get("wp_state_index")
    if val is not None:
        try:
            return int(val)
        except ValueError:
            pass
    return None


def _get_n_occupied(results_dir: Path) -> int | None:
    """Read n_occupied from run_summary.txt."""
    params = _parse_run_summary(results_dir)
    val = params.get("n_occupied")
    if val is not None:
        try:
            return int(val)
        except ValueError:
            pass
    return None


def _get_wp_sigma(results_dir: Path) -> float:
    """Read wp_sigma_bohr from run_summary.txt; default 5.0."""
    params = _parse_run_summary(results_dir)
    val = params.get("wp_sigma_bohr")
    if val is not None:
        try:
            return float(val)
        except ValueError:
            pass
    return 5.0


def _get_wp_energy_ev(results_dir: Path) -> float | None:
    """Read wp_energy_ev from run_summary.txt."""
    params = _parse_run_summary(results_dir)
    val = params.get("wp_energy_ev")
    if val is not None:
        try:
            return float(val)
        except ValueError:
            pass
    return None


def _get_cell_bohr(results_dir: Path) -> float:
    """Read cell side length from run_summary.txt; default 50.0."""
    params = _parse_run_summary(results_dir)
    val = params.get("cell_bohr")
    if val is not None:
        # Parse "50^3 (cubic, periodic)" or "30.0"
        m = re.match(r"([\d.]+)", val)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                pass
    return 50.0


def _get_launch_z(results_dir: Path) -> float:
    """Read WP launch z coordinate from run_summary.txt."""
    params = _parse_run_summary(results_dir)
    val = params.get("wp_center_bohr")
    if val is not None:
        parts = val.split()
        if len(parts) >= 3:
            try:
                return float(parts[2])
            except ValueError:
                pass
    return -21.0  # common default


def _v0_from_ekin(E_eV: float) -> float:
    """Projectile velocity in a.u. from kinetic energy in eV."""
    return np.sqrt(2.0 * E_eV / HA_TO_EV)


def _compute_ifw_window(results_dir: Path) -> tuple[float, float]:
    """Compute the interference-free window [z_min, z_max] in Dz (Bohr).

    Uses the boundary rule: launch clearance = 3 Bohr,
    stop = L/2 - 3*sigma - 2 Bohr (or 28 Bohr if L=50 and sigma >= 3).
    """
    sigma = _get_wp_sigma(results_dir)
    L = _get_cell_bohr(results_dir)
    z_min = 3.0
    z_max = L / 2.0 - 3.0 * sigma - 2.0
    if z_max <= z_min:
        # Fallback for very large sigma: use 80% of half-box
        z_max = 0.8 * L / 2.0
    return z_min, z_max


def _read_overlap_at(overlap_dir: Path, step_query: int | None = None
                     ) -> tuple[np.ndarray, int, float] | None:
    """Read an overlap matrix CSV at the given step (last if None).

    Returns (O_squared[n_ref, n_evolved], step, time_au) or None.
    """
    idx_csv = overlap_dir / "index.csv"
    if not idx_csv.exists():
        return None
    df = pd.read_csv(idx_csv)
    if df.empty:
        return None
    if step_query is None:
        row = df.iloc[-1]
    else:
        row = df.iloc[(df["step"] - step_query).abs().argmin()]
    csv_path = overlap_dir / row["file"]
    if not csv_path.exists():
        return None
    arr = pd.read_csv(csv_path, comment="#", header=None).to_numpy()
    return arr, int(row["step"]), float(row["time_au"])


def _read_occupations(results_dir: Path) -> np.ndarray | None:
    """Read GS occupations from eigenvalues/occupations.csv."""
    occ_csv = results_dir / "raw" / "observables" / "eigenvalues" / "occupations.csv"
    if not occ_csv.exists():
        return None
    df = pd.read_csv(occ_csv).sort_values("state_index")
    return df["occupation"].to_numpy()


def _read_v_initial(results_dir: Path) -> float | None:
    """Parse classical projectile initial v_z (a.u.) from run_summary.txt."""
    rs = results_dir / "run_summary.txt"
    if not rs.exists():
        return None
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


# =============================================================================
# 1. Energy decomposition — 6-panel plot
# =============================================================================

def _energy_series(results_dir: Path, *, is_wp: bool) -> pd.DataFrame | None:
    """Build unified energy-component DataFrame for a single run.

    Columns: step, time_au, dz_au, dE_total_ev, dE_kinetic_ev,
    dE_hartree_ev, dE_xc_ev, dE_bath_ev, dE_wp_slot_ev.
    """
    obs_csv = results_dir / "raw" / "observables" / "observables.csv"
    se_csv = results_dir / "raw" / "observables" / "state_energies.csv"
    if not obs_csv.exists() or not se_csv.exists():
        return None

    obs = pd.read_csv(obs_csv)
    se = pd.read_csv(se_csv)
    wp_idx = _get_wp_state_index(results_dir) if is_wp else None

    # Compute v0 from WP energy or classical velocity
    E_eV = _get_wp_energy_ev(results_dir)
    v_init = _read_v_initial(results_dir)
    if E_eV is not None:
        v0 = _v0_from_ekin(E_eV)
    elif v_init is not None:
        v0 = v_init
    else:
        v0 = 1.0  # fallback

    obs["dz_au"] = v0 * obs["time_au"]

    for col, base in [
        ("dE_total_ev",   "energy_total"),
        ("dE_kinetic_ev", "energy_kinetic"),
        ("dE_hartree_ev", "energy_hartree"),
        ("dE_xc_ev",      "energy_xc"),
    ]:
        obs[col] = (obs[base] - obs[base].iloc[0]) * HA_TO_EV

    # Bath orbital state sum (exclude WP slot)
    if wp_idx is not None:
        bath = se[se["state_index"] != wp_idx].copy()
    else:
        bath = se.copy()
    bath["contrib"] = bath["weight"] * bath["occupation"] * bath["E_expect_ha"]
    bath_g = bath.groupby(["step", "time_au"])["contrib"].sum().reset_index()
    bath_g["dE_bath_ev"] = (bath_g["contrib"] - bath_g["contrib"].iloc[0]) * HA_TO_EV

    # WP slot energy
    if wp_idx is not None:
        wp_se = se[se["state_index"] == wp_idx].sort_values("step").copy()
        wp_se["dE_wp_slot_ev"] = (wp_se["E_expect_ha"] - wp_se["E_expect_ha"].iloc[0]) * HA_TO_EV
        df = obs.merge(bath_g[["step", "dE_bath_ev"]], on="step", how="inner")
        df = df.merge(wp_se[["step", "dE_wp_slot_ev"]], on="step", how="inner")
    else:
        df = obs.merge(bath_g[["step", "dE_bath_ev"]], on="step", how="inner")
        df["dE_wp_slot_ev"] = np.nan

    return df


def plot_energy_decomposition(
    results_dir: Path,
    run_name: str,
    classical_results_dir: Optional[Path] = None,
) -> Path | None:
    """6-panel energy decomposition: total, kinetic, Hartree, xc, bath sum,
    WP slot vs Dz (projectile position from launch).

    If classical_results_dir is provided, overlay classical curves.

    Output: results/analysis/observables/energy_decomposition.png
    """
    print(f"  [analyse_v2] plot_energy_decomposition for {run_name}...")

    # Determine whether this run is WP or classical
    wp_idx = _get_wp_state_index(results_dir)
    is_wp = wp_idx is not None

    df = _energy_series(results_dir, is_wp=is_wp)
    if df is None:
        print("    [skip] missing observables or state_energies")
        return None

    df_cls = None
    if classical_results_dir is not None:
        df_cls = _energy_series(classical_results_dir, is_wp=False)

    # IFW window
    z_min, z_max = _compute_ifw_window(results_dir)

    out_dir = results_dir / "analysis" / "observables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "energy_decomposition.png"

    fig, axes = plt.subplots(3, 2, figsize=(13, 11), sharex=True)

    panels = [
        ("dE_total_ev",   r"$\Delta E_{\mathrm{total}}$ [eV]",      "Total energy"),
        ("dE_kinetic_ev", r"$\Delta E_{\mathrm{kinetic}}$ [eV]",    "Kinetic energy"),
        ("dE_hartree_ev", r"$\Delta E_{\mathrm{Hartree}}$ [eV]",    "Hartree energy"),
        ("dE_xc_ev",      r"$\Delta E_{\mathrm{xc}}$ [eV]",         "Exchange-correlation"),
        ("dE_bath_ev",    r"$\Delta E_{\mathrm{bath}}$ [eV]",        "Bath orbital sum"),
        ("dE_wp_slot_ev", r"$\Delta E_{\mathrm{WP\,slot}}$ [eV]",   "WP slot energy"),
    ]

    label_primary = "WP" if is_wp else "Classical"
    color_primary = "firebrick" if is_wp else "navy"

    for ax, (col, ylabel, subtitle) in zip(axes.flat, panels):
        # IFW shading
        ax.axvspan(z_min, z_max, alpha=0.10, color="grey", zorder=0)

        # Primary run
        if col == "dE_wp_slot_ev" and not is_wp:
            ax.text(0.5, 0.5, "N/A (classical run)",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=10, color="0.5", style="italic")
        else:
            ax.plot(df["dz_au"], df[col], "-", color=color_primary, lw=1.8,
                    label=label_primary)

        # Classical overlay
        if df_cls is not None and col != "dE_wp_slot_ev":
            ax.plot(df_cls["dz_au"], df_cls[col], "--", color="navy", lw=1.5,
                    alpha=0.7, label="Classical")

        ax.axhline(0, color="black", lw=0.5, alpha=0.5)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(subtitle, fontsize=10)
        ax.grid(True, ls=":", lw=0.4, alpha=0.6)
        ax.tick_params(direction="in", which="both")
        ax.legend(loc="best", fontsize=8, framealpha=0.92)

    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\Delta z = v_0 \cdot t$ [Bohr]", fontsize=10.5)

    fig.suptitle(f"{run_name}: energy decomposition vs projectile displacement\n"
                 f"IFW: $\\Delta z \\in [{z_min:.0f},\\, {z_max:.0f}]$ Bohr "
                 f"(grey shaded)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {out_png}")
    return out_png


# =============================================================================
# 2. Energy bookkeeping bar chart — cross-run comparison at t_IFW
# =============================================================================

def plot_energy_bookkeeping_bar(
    results_dirs: list[Path],
    run_names: list[str],
    output_path: Optional[Path] = None,
) -> Path | None:
    """Cross-run grouped bar chart of energy components at t_IFW.

    For each run, shows DeltaE_kinetic, DeltaE_hartree, DeltaE_xc
    at the end of the IFW window.

    output_path: if None, writes to the first run's analysis dir.
    """
    print(f"  [analyse_v2] plot_energy_bookkeeping_bar for {len(results_dirs)} runs...")

    if len(results_dirs) == 0:
        return None

    # Collect data for each run
    run_data = []
    for rdir, rname in zip(results_dirs, run_names):
        wp_idx = _get_wp_state_index(rdir)
        is_wp = wp_idx is not None
        df = _energy_series(rdir, is_wp=is_wp)
        if df is None:
            print(f"    [skip] {rname}: missing data")
            continue

        z_min, z_max = _compute_ifw_window(rdir)

        # Find row closest to z_max (end of IFW)
        idx_end = (df["dz_au"] - z_max).abs().idxmin()
        row = df.loc[idx_end]

        run_data.append({
            "name": rname,
            "dE_kinetic": row["dE_kinetic_ev"],
            "dE_hartree": row["dE_hartree_ev"],
            "dE_xc": row["dE_xc_ev"],
            "z_ifw": z_max,
        })

    if len(run_data) == 0:
        return None

    if output_path is None:
        output_path = results_dirs[0] / "analysis" / "observables" / "energy_bookkeeping_bar.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    names = [d["name"] for d in run_data]
    dE_kin = [d["dE_kinetic"] for d in run_data]
    dE_har = [d["dE_hartree"] for d in run_data]
    dE_xc = [d["dE_xc"] for d in run_data]

    x = np.arange(len(names))
    width = 0.25

    fig, ax = plt.subplots(figsize=(max(8, 3 * len(names)), 6))
    ax.bar(x - width, dE_kin, width, color="steelblue", edgecolor="black",
           lw=0.6, label=r"$\Delta E_{\mathrm{kinetic}}$", zorder=3)
    ax.bar(x, dE_har, width, color="coral", edgecolor="black",
           lw=0.6, label=r"$\Delta E_{\mathrm{Hartree}}$", zorder=3)
    ax.bar(x + width, dE_xc, width, color="mediumpurple", edgecolor="black",
           lw=0.6, label=r"$\Delta E_{\mathrm{xc}}$", zorder=3)

    # Value labels
    for bars, vals in [(x - width, dE_kin), (x, dE_har), (x + width, dE_xc)]:
        for bx, val in zip(bars, vals):
            y_off = 0.02 if val >= 0 else -0.02
            ax.text(bx, val + y_off, f"{val:+.3f}",
                    ha="center", va="bottom" if val >= 0 else "top",
                    fontsize=8, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=15, ha="right", fontsize=10)
    ax.set_ylabel(r"$\Delta E$ at IFW end [eV]", fontsize=11)
    ax.set_title("Energy bookkeeping at end of IFW window", fontsize=11.5)
    ax.axhline(0, color="black", lw=0.6)
    ax.grid(True, axis="y", ls=":", lw=0.4, alpha=0.6)
    ax.legend(loc="best", fontsize=10)
    ax.tick_params(direction="in", which="both")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {output_path}")
    return output_path


# =============================================================================
# 3. KS eigenenergy evolution
# =============================================================================

def plot_ks_eigenenergy_evolution(
    results_dir: Path,
    run_name: str,
) -> Path | None:
    """Plot all KS orbital energies E_i(t) vs time, coloured by initial
    occupation (blue=occupied, red=virtual), with HOMO and LUMO highlighted.

    Output: results/analysis/observables/ks_eigenenergy_evolution.png
    """
    print(f"  [analyse_v2] plot_ks_eigenenergy_evolution for {run_name}...")

    se_csv = results_dir / "raw" / "observables" / "state_energies.csv"
    if not se_csv.exists():
        print("    [skip] no state_energies.csv")
        return None

    se = pd.read_csv(se_csv)
    if se.empty:
        print("    [skip] state_energies.csv is empty")
        return None

    wp_idx = _get_wp_state_index(results_dir)
    n_occ = _get_n_occupied(results_dir)

    # Get initial occupations for colour coding
    t0 = se[se["step"] == se["step"].min()]
    occ_map = dict(zip(t0["state_index"], t0["occupation"]))

    state_indices = sorted(se["state_index"].unique())

    # If n_occ not available from run_summary, infer from occupations
    if n_occ is None:
        n_occ = sum(1 for si in state_indices if occ_map.get(si, 0) > 0.5)

    # Determine x-axis: if we have WP energy, use Dz = v0*t; else time
    E_eV = _get_wp_energy_ev(results_dir)
    v_init = _read_v_initial(results_dir)
    use_dz = False
    v0 = None
    if E_eV is not None:
        v0 = _v0_from_ekin(E_eV)
        use_dz = True
    elif v_init is not None:
        v0 = v_init
        use_dz = True

    # IFW window
    z_min, z_max = _compute_ifw_window(results_dir)

    out_dir = results_dir / "analysis" / "observables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "ks_eigenenergy_evolution.png"

    fig, ax = plt.subplots(figsize=(11, 6.5))

    homo_idx = n_occ - 1 if n_occ > 0 else None
    lumo_idx = n_occ if n_occ < len(state_indices) else None

    for si in state_indices:
        if si == wp_idx:
            continue  # plot WP separately

        sub = se[se["state_index"] == si].sort_values("step")
        if len(sub) < 2:
            continue

        # Delta E in meV
        dE_mev = (sub["E_expect_ha"].to_numpy()
                  - sub["E_expect_ha"].iloc[0]) * HA_TO_EV * 1000

        if use_dz:
            x_vals = v0 * sub["time_au"].to_numpy()
        else:
            x_vals = sub["time_au"].to_numpy()

        occ_i = occ_map.get(si, 0)
        is_occupied = occ_i > 0.5

        # Colour: blue for occupied, red for virtual
        if si == homo_idx:
            ax.plot(x_vals, dE_mev, "-", color="royalblue", lw=2.0,
                    zorder=4, label=f"HOMO (state {si})")
        elif si == lumo_idx:
            ax.plot(x_vals, dE_mev, "-", color="orangered", lw=2.0,
                    zorder=4, label=f"LUMO (state {si})")
        elif is_occupied:
            ax.plot(x_vals, dE_mev, "-", color="steelblue", lw=0.5,
                    alpha=0.4, zorder=2)
        else:
            ax.plot(x_vals, dE_mev, "-", color="salmon", lw=0.5,
                    alpha=0.4, zorder=2)

    # WP slot in bright red
    if wp_idx is not None:
        sub = se[se["state_index"] == wp_idx].sort_values("step")
        if len(sub) >= 2:
            dE_mev = (sub["E_expect_ha"].to_numpy()
                      - sub["E_expect_ha"].iloc[0]) * HA_TO_EV * 1000
            if use_dz:
                x_vals = v0 * sub["time_au"].to_numpy()
            else:
                x_vals = sub["time_au"].to_numpy()
            ax.plot(x_vals, dE_mev, "-", color="red", lw=2.5,
                    zorder=5, label=f"WP slot (state {wp_idx})")

    # IFW shading
    if use_dz:
        ax.axvspan(z_min, z_max, alpha=0.08, color="grey", zorder=0,
                   label=f"IFW [{z_min:.0f}, {z_max:.0f}] Bohr")

    # Dummy entries for occupied/virtual legend
    ax.plot([], [], "-", color="steelblue", lw=1.5, alpha=0.6,
            label=f"Occupied (states 0-{n_occ-1})")
    ax.plot([], [], "-", color="salmon", lw=1.5, alpha=0.6,
            label=f"Virtual (states {n_occ}+)")

    if use_dz:
        ax.set_xlabel(r"$\Delta z = v_0 \cdot t$ [Bohr]", fontsize=11)
    else:
        ax.set_xlabel("Time [a.u.]", fontsize=11)
    ax.set_ylabel(r"$\langle H \rangle_i(t) - \langle H \rangle_i(0)$ [meV]",
                  fontsize=11)
    ax.set_title(f"{run_name}: KS eigenenergy evolution", fontsize=11.5)
    ax.grid(True, ls=":", lw=0.4, alpha=0.6)
    ax.legend(loc="best", fontsize=9, framealpha=0.92)
    ax.tick_params(direction="in", which="both")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {out_png}")
    return out_png


# =============================================================================
# 4. GS-basis decomposition
# =============================================================================

def plot_gs_basis_decomposition(
    results_dir: Path,
    run_name: str,
) -> Path | None:
    """Bar chart of Delta n_i^GS = change in GS-projected occupation at t_end.

    Occupied orbitals (i <= n_occ) show depletion (negative, red bars).
    Virtual orbitals (i > n_occ) show excitation (positive, blue bars).
    Fermi level marked. Inset: charge-conservation check.

    Reads from overlap_full (or overlap_proxies as fallback).

    Output: results/analysis/observables/gs_basis_decomposition.png
    """
    print(f"  [analyse_v2] plot_gs_basis_decomposition for {run_name}...")

    # Try overlap_full first, then overlap_proxies
    full_dir = results_dir / "raw" / "observables" / "overlap_full"
    proxy_dir = results_dir / "raw" / "observables" / "overlap_proxies"

    data = _read_overlap_at(full_dir)
    overlap_source = "overlap_full"
    if data is None:
        data = _read_overlap_at(proxy_dir)
        overlap_source = "overlap_proxies"
    if data is None:
        print("    [skip] no overlap_full or overlap_proxies data")
        return None

    O2, step, t_au = data  # shape (n_ref, n_evolved)
    n_ref, n_evolved = O2.shape

    # Read GS occupations
    occupations = _read_occupations(results_dir)
    if occupations is None:
        print("    [skip] no occupations.csv")
        return None

    # Pad/trim occupations to match n_evolved
    f_evolved = occupations.copy()
    if f_evolved.size < n_evolved:
        f_evolved = np.concatenate([f_evolved, np.zeros(n_evolved - f_evolved.size)])
    elif f_evolved.size > n_evolved:
        f_evolved = f_evolved[:n_evolved]

    # n_i^GS(t) = sum_j f_j(0) |<psi_i^GS | psi_j(t)>|^2
    n_gs = O2 @ f_evolved
    # Reference: GS occupations for i < n_ref
    n_gs_ref = f_evolved[:n_ref]

    delta = n_gs - n_gs_ref

    n_occ = _get_n_occupied(results_dir)
    if n_occ is None:
        # Infer from occupations (count entries > 0.5)
        n_occ = int(np.sum(f_evolved[:n_ref] > 0.5))

    # Charge-conservation check
    sum_depletion = float(delta[:n_occ].sum())
    sum_excitation = float(delta[n_occ:].sum()) if n_occ < n_ref else 0.0
    total_delta = float(delta.sum())
    unaccounted = -total_delta  # leakage into unsaved virtuals

    out_dir = results_dir / "analysis" / "observables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "gs_basis_decomposition.png"

    fig, (ax_main, ax_inset) = plt.subplots(
        1, 2, figsize=(13, 5.5),
        gridspec_kw={"width_ratios": [10, 1.5]})

    idx = np.arange(n_ref)

    # Colour each bar: red for depleted (occupied), blue for excited (virtual)
    colors = []
    for i in range(n_ref):
        if i < n_occ:
            colors.append("indianred")
        else:
            colors.append("steelblue")

    ax_main.bar(idx, delta * 1000, color=colors, edgecolor="black",
                lw=0.3, width=0.85, zorder=3)

    # Fermi level
    ax_main.axvline(n_occ - 0.5, color="orange", lw=1.5, ls="--", alpha=0.8,
                    label=f"Fermi level (n_occ={n_occ})")
    ax_main.axhline(0, color="black", lw=0.5)

    # Annotate top-3 depletions and excitations
    delta_occ = delta[:n_occ]
    if len(delta_occ) > 0:
        top_dep_idx = np.argsort(delta_occ)[:3]
        for i in top_dep_idx:
            if abs(delta[i]) > 1e-6:
                ax_main.annotate(f"#{i}", xy=(i, delta[i] * 1000),
                                 xytext=(0, -10), textcoords="offset points",
                                 fontsize=7, color="darkred", ha="center")

    delta_virt = delta[n_occ:]
    if len(delta_virt) > 0:
        top_exc_idx = np.argsort(-delta_virt)[:3] + n_occ
        for i in top_exc_idx:
            if i < n_ref and delta[i] > 1e-6:
                ax_main.annotate(f"#{i}", xy=(i, delta[i] * 1000),
                                 xytext=(0, 8), textcoords="offset points",
                                 fontsize=7, color="darkblue", ha="center")

    ax_main.set_xlabel("GS orbital index $i$", fontsize=11)
    ax_main.set_ylabel(r"$\delta n_i^{\mathrm{GS}}(t_{\mathrm{end}})$ [$\times 10^{-3}$]",
                       fontsize=11)
    ax_main.set_title(f"{run_name}: GS-basis decomposition at t = {t_au:.3g} a.u. "
                      f"(step {step}, {overlap_source})", fontsize=10.5)
    ax_main.grid(True, axis="y", ls=":", lw=0.4, alpha=0.6)
    ax_main.legend(loc="best", fontsize=9, framealpha=0.92)
    ax_main.tick_params(direction="in", which="both")

    # Conservation check annotation
    ax_main.text(0.99, 0.05,
                 f"$\\Sigma_{{occ}}\\delta n_i$ = {sum_depletion*1000:+.2f}e-3   "
                 f"$\\Sigma_{{virt}}\\delta n_i$ = {sum_excitation*1000:+.2f}e-3   "
                 f"total = {total_delta*1000:+.2f}e-3   "
                 f"unaccounted = {unaccounted*1000:+.2f}e-3",
                 transform=ax_main.transAxes, ha="right", va="bottom",
                 fontsize=8, bbox=dict(boxstyle="round,pad=0.25",
                                       fc="lightyellow", ec="black",
                                       lw=0.5, alpha=0.9))

    # Inset: charge conservation bars
    bar_labels = ["depletion\n(occ)", "excitation\n(virt)", "unaccounted"]
    bar_vals = [sum_depletion * 1000, sum_excitation * 1000, unaccounted * 1000]
    bar_colors = ["indianred", "steelblue", "grey"]
    ax_inset.bar(range(3), bar_vals, color=bar_colors, edgecolor="black",
                 lw=0.8, width=0.6, zorder=3)
    ax_inset.axhline(0, color="black", lw=0.5)
    ax_inset.set_xticks(range(3))
    ax_inset.set_xticklabels(bar_labels, fontsize=7.5)
    ax_inset.set_ylabel(r"$\Sigma\,\delta n_i$ [$\times 10^{-3}$]", fontsize=9)
    ax_inset.set_title("Conservation", fontsize=9)
    ax_inset.grid(True, axis="y", ls=":", lw=0.4, alpha=0.6)
    ax_inset.tick_params(direction="in", which="both", labelsize=8)

    # Value labels on inset bars
    for i, val in enumerate(bar_vals):
        y_off = 0.02 if val >= 0 else -0.02
        ax_inset.text(i, val + y_off, f"{val:+.1f}",
                      ha="center", va="bottom" if val >= 0 else "top",
                      fontsize=7.5, fontweight="bold")

    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {out_png}")
    return out_png


# =============================================================================
# 5. Overlap heatmap (log scale)
# =============================================================================

def plot_overlap_heatmap(
    results_dir: Path,
    run_name: str,
    log_scale: bool = True,
) -> Path | None:
    """Heatmap of |<psi_i^GS | psi_j(t_end)>|^2 with diagonal masked,
    log colour scale, and occupied/virtual boundary.

    Output: results/analysis/observables/overlap_heatmap_log.png
    """
    print(f"  [analyse_v2] plot_overlap_heatmap for {run_name}...")

    full_dir = results_dir / "raw" / "observables" / "overlap_full"
    data = _read_overlap_at(full_dir)
    if data is None:
        print("    [skip] no overlap_full data")
        return None

    O2, step, t_au = data
    n_ref, n_evolved = O2.shape

    n_occ = _get_n_occupied(results_dir)
    if n_occ is None:
        occupations = _read_occupations(results_dir)
        if occupations is not None:
            n_occ = int(np.sum(occupations > 0.5))
        else:
            n_occ = n_ref // 2  # rough fallback

    out_dir = results_dir / "analysis" / "observables"
    out_dir.mkdir(parents=True, exist_ok=True)
    suffix = "log" if log_scale else "linear"
    out_png = out_dir / f"overlap_heatmap_{suffix}.png"

    # Mask diagonal for off-diagonal focus
    O2_masked = O2.copy()
    diag_size = min(n_ref, n_evolved)
    for i in range(diag_size):
        O2_masked[i, i] = np.nan

    fig, ax = plt.subplots(figsize=(8, 7))

    if log_scale:
        from matplotlib.colors import LogNorm
        # Replace zeros/nans for log scale
        O2_plot = np.where((O2_masked > 0) & np.isfinite(O2_masked),
                           O2_masked, np.nan)
        vmin = np.nanmin(O2_plot[O2_plot > 0]) if np.any(O2_plot > 0) else 1e-8
        vmax = np.nanmax(O2_plot) if np.any(np.isfinite(O2_plot)) else 1.0
        im = ax.imshow(O2_plot, origin="lower", aspect="auto",
                       norm=LogNorm(vmin=max(vmin, 1e-8), vmax=vmax),
                       cmap="viridis")
    else:
        im = ax.imshow(O2_masked, origin="lower", aspect="auto",
                       vmin=0, vmax=1.0, cmap="viridis")

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(r"$|\langle\psi_i^{\mathrm{GS}}|\psi_j(t)\rangle|^2$"
                 + (" (diagonal masked)" if True else ""))

    # Occupied/virtual boundary lines
    if n_occ < n_ref:
        ax.axhline(n_occ - 0.5, color="white", lw=1.5, ls="--", alpha=0.8)
    if n_occ < n_evolved:
        ax.axvline(n_occ - 0.5, color="white", lw=1.5, ls="--", alpha=0.8)

    ax.set_xlabel("Evolved KS orbital index $j$", fontsize=11)
    ax.set_ylabel("Ground-state KS orbital index $i$", fontsize=11)
    ax.set_title(f"{run_name}: overlap matrix at step {step}, "
                 f"t = {t_au:.3g} a.u.\n"
                 f"Diagonal masked; {'log' if log_scale else 'linear'} scale; "
                 f"dashed lines = Fermi level (n_occ={n_occ})",
                 fontsize=10.5)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {out_png}")
    return out_png


# =============================================================================
# 6. Overlap heatmap difference (WP - classical)
# =============================================================================

def plot_overlap_heatmap_diff(
    wp_results_dir: Path,
    classical_results_dir: Path,
    run_name: str,
    output_path: Optional[Path] = None,
) -> Path | None:
    """Difference heatmap: WP overlap - classical overlap at t_end.

    Uses a diverging (RdBu) colourmap.
    """
    print(f"  [analyse_v2] plot_overlap_heatmap_diff for {run_name}...")

    wp_full = wp_results_dir / "raw" / "observables" / "overlap_full"
    cls_full = classical_results_dir / "raw" / "observables" / "overlap_full"

    data_wp = _read_overlap_at(wp_full)
    data_cls = _read_overlap_at(cls_full)

    if data_wp is None or data_cls is None:
        print("    [skip] missing overlap_full in one or both runs")
        return None

    O2_wp, step_wp, t_wp = data_wp
    O2_cls, step_cls, t_cls = data_cls

    # Align to common shape
    n_ref = min(O2_wp.shape[0], O2_cls.shape[0])
    n_evolved = min(O2_wp.shape[1], O2_cls.shape[1])
    diff = O2_wp[:n_ref, :n_evolved] - O2_cls[:n_ref, :n_evolved]

    n_occ = _get_n_occupied(wp_results_dir)
    if n_occ is None:
        n_occ = n_ref // 2

    if output_path is None:
        output_path = (wp_results_dir / "analysis" / "observables"
                       / "overlap_heatmap_diff.png")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    vmax = max(abs(np.nanmin(diff)), abs(np.nanmax(diff)))
    if vmax < 1e-10:
        vmax = 1.0  # avoid degenerate scale

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(diff, origin="lower", aspect="auto",
                   vmin=-vmax, vmax=vmax, cmap="RdBu_r")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(r"$\Delta |O_{ij}|^2$ (WP $-$ classical)")

    # Occupied/virtual boundary
    if n_occ < n_ref:
        ax.axhline(n_occ - 0.5, color="black", lw=1.5, ls="--", alpha=0.7)
    if n_occ < n_evolved:
        ax.axvline(n_occ - 0.5, color="black", lw=1.5, ls="--", alpha=0.7)

    ax.set_xlabel("Evolved KS orbital index $j$", fontsize=11)
    ax.set_ylabel("Ground-state KS orbital index $i$", fontsize=11)
    ax.set_title(f"{run_name}: overlap difference (WP - classical)\n"
                 f"WP: step {step_wp}, t={t_wp:.3g} a.u.  |  "
                 f"Classical: step {step_cls}, t={t_cls:.3g} a.u.\n"
                 f"Dashed = Fermi level (n_occ={n_occ}); "
                 f"RdBu diverging colourmap",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {output_path}")
    return output_path


# =============================================================================
# Convenience: run all single-run analyses
# =============================================================================

def run_all(results_dir: Path, run_name: str, *,
            classical_results_dir: Optional[Path] = None,
            ) -> dict[str, str]:
    """Run all single-run analyse_v2 plots and return a {label: status} log."""
    log: dict[str, str] = {}

    p = plot_energy_decomposition(results_dir, run_name,
                                  classical_results_dir=classical_results_dir)
    log["energy_decomposition"] = f"[ok] {p}" if p else "[skip] missing data"

    p = plot_ks_eigenenergy_evolution(results_dir, run_name)
    log["ks_eigenenergy_evolution"] = f"[ok] {p}" if p else "[skip] missing data"

    p = plot_gs_basis_decomposition(results_dir, run_name)
    log["gs_basis_decomposition"] = f"[ok] {p}" if p else "[skip] missing data"

    p = plot_overlap_heatmap(results_dir, run_name, log_scale=True)
    log["overlap_heatmap_log"] = f"[ok] {p}" if p else "[skip] missing data"

    return log


# =============================================================================
# CLI entry point
# =============================================================================

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python analyse_v2.py <results_dir> <run_name> [classical_results_dir]")
        sys.exit(1)
    rdir = Path(sys.argv[1])
    rname = sys.argv[2]
    cdir = Path(sys.argv[3]) if len(sys.argv) > 3 else None
    log = run_all(rdir, rname, classical_results_dir=cdir)
    print("\n=== analyse_v2 summary ===")
    for k, v in log.items():
        print(f"  {k}: {v}")
