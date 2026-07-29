"""stopping_power_data — Compute stopping power S(E) from jellium run data.

Reads raw CSV observables from WP and classical Ehrenfest runs, computes
interference-free time windows, and extracts stopping power via two
definitions (momentum S₁, KS orbital energy S₂) plus classical linear
regression.

Run:
    python -m applications.report1.stopping_power_data   # prints all data points
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import linregress

HA_TO_EV = 27.2114
AU_TO_FS = 0.02418884  # 1 a.u. of time = 0.02419 fs

JBASE = Path("ResearchProject/systems/jellium")


@dataclass
class RunParams:
    """Parameters extracted from run_summary.txt."""
    cell_bohr: float
    wp_state_index: int
    wp_center_z: float
    wp_sigma: float
    wp_k0_z: float
    dt_au: float
    rt_num_steps: int
    is_classical: bool = False


@dataclass
class TimeWindow:
    """Interference-free propagation window."""
    t_start: float = 0.0  # a.u.
    t_end: float = 0.0    # a.u.
    n_steps_end: int = 0
    boundary_clear: bool = True


@dataclass
class StoppingResult:
    """Result of a stopping power computation."""
    run_dir: str
    energy_eV: float
    sigma: float
    S1_eV_per_bohr: float = np.nan  # momentum definition
    S2_eV_per_bohr: float = np.nan  # KS orbital definition
    S_classical: float = np.nan
    S_classical_err: float = np.nan
    window: Optional[TimeWindow] = None
    is_v2: bool = False
    is_boundary_compromised: bool = False
    dE1_eV: float = np.nan
    dE2_eV: float = np.nan
    dz_bohr: float = np.nan
    E_kin_0_eV: float = np.nan
    E_kin_end_eV: float = np.nan


def parse_run_summary(path: Path) -> RunParams:
    """Parse run_summary.txt for simulation parameters."""
    text = (path / "results" / "run_summary.txt").read_text()

    cell_match = re.search(r"cell_bohr\s*=\s*(\d+(?:\.\d+)?)", text)
    cell = float(cell_match.group(1)) if cell_match else 50.0

    idx_match = re.search(r"wp_state_index\s*=\s*(\d+)", text)
    wp_idx = int(idx_match.group(1)) if idx_match else 100

    center_match = re.search(r"wp_center_bohr\s*=\s*([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)", text)
    center_z = float(center_match.group(3)) if center_match else 0.0

    sigma_match = re.search(r"wp_sigma_bohr\s*=\s*([\d.e+-]+)", text)
    sigma = float(sigma_match.group(1)) if sigma_match else 1.0

    k0_match = re.search(r"wp_k0_bohr_inv\s*=\s*([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)", text)
    k0_z = float(k0_match.group(3)) if k0_match else 0.0

    dt_match = re.search(r"dt_au\s*=\s*([\d.e+-]+)", text)
    dt = float(dt_match.group(1)) if dt_match else 0.01

    nsteps_match = re.search(r"rt_num_steps\s*=\s*(\d+)", text)
    nsteps = int(nsteps_match.group(1)) if nsteps_match else 0

    is_classical = "wp_sigma_bohr" not in text

    return RunParams(
        cell_bohr=cell,
        wp_state_index=wp_idx,
        wp_center_z=center_z,
        wp_sigma=sigma,
        wp_k0_z=k0_z,
        dt_au=dt,
        rt_num_steps=nsteps,
        is_classical=is_classical,
    )


def compute_time_window(params: RunParams, n_sigma: float = 4.0) -> TimeWindow:
    """Compute interference-free window: t_end when leading edge hits L/2.

    σ_r(t) = sqrt(σ²/2 + t²/(2σ²))   — density envelope width
    z_leading(t) = z₀ + k₀*t + n_sigma * σ_r(t)
    Solve z_leading(t_end) = L/2
    """
    sigma = params.wp_sigma
    z0 = params.wp_center_z
    k0 = params.wp_k0_z
    L = params.cell_bohr
    boundary = L / 2.0

    def sigma_r(t):
        return np.sqrt(sigma**2 / 2.0 + t**2 / (2.0 * sigma**2))

    def z_leading(t):
        return z0 + k0 * t + n_sigma * sigma_r(t)

    t_max = params.dt_au * params.rt_num_steps

    if z_leading(0) >= boundary:
        return TimeWindow(t_end=0, boundary_clear=False)

    if z_leading(t_max) < boundary:
        n_end = params.rt_num_steps
        return TimeWindow(t_end=t_max, n_steps_end=n_end, boundary_clear=True)

    t_end = brentq(lambda t: z_leading(t) - boundary, 0, t_max)
    n_end = int(t_end / params.dt_au)

    rear_clearance = z0 - n_sigma * sigma_r(0) + L / 2.0
    boundary_clear = rear_clearance > 0

    return TimeWindow(t_end=t_end, n_steps_end=n_end, boundary_clear=boundary_clear)


def compute_wp_S1(run_path: Path, window: TimeWindow, params: RunParams) -> tuple[float, float, float, float, float]:
    """Compute S₁ = -ΔE_kin / Δz (momentum definition).

    E_kin = ⟨p_z⟩²/2 + σ_pz²/2  (atomic units, m_e=1)

    Returns: (S1_eV_per_bohr, dE_eV, dz_bohr, E_kin_0_eV, E_kin_end_eV)
    """
    csv_path = run_path / "results" / "raw" / "observables" / "wp_momentum_stats.csv"
    if not csv_path.exists():
        return (np.nan,) * 5

    df = pd.read_csv(csv_path, comment="#")

    mask = df["time_au"] <= window.t_end
    df_win = df[mask]
    if len(df_win) < 2:
        return (np.nan,) * 5

    pz_0 = df_win["pz_mean"].iloc[0]
    sig_pz2_0 = df_win["sigma_pz2"].iloc[0]
    E_kin_0 = pz_0**2 / 2.0 + sig_pz2_0 / 2.0  # Hartree

    pz_end = df_win["pz_mean"].iloc[-1]
    sig_pz2_end = df_win["sigma_pz2"].iloc[-1]
    E_kin_end = pz_end**2 / 2.0 + sig_pz2_end / 2.0

    dE_ha = E_kin_end - E_kin_0
    dE_eV = dE_ha * HA_TO_EV

    dz = params.wp_k0_z * window.t_end
    S1 = -dE_eV / dz if abs(dz) > 0.1 else np.nan

    return (S1, dE_eV, dz, E_kin_0 * HA_TO_EV, E_kin_end * HA_TO_EV)


def compute_wp_S2(run_path: Path, window: TimeWindow, params: RunParams) -> tuple[float, float]:
    """Compute S₂ = -ΔE_KS / Δz (KS orbital energy definition).

    Returns: (S2_eV_per_bohr, dE_KS_eV)
    """
    csv_path = run_path / "results" / "raw" / "observables" / "state_energies.csv"
    if not csv_path.exists():
        return (np.nan, np.nan)

    df = pd.read_csv(csv_path)
    wp_state = df[df["state_index"] == params.wp_state_index].copy()
    if wp_state.empty:
        return (np.nan, np.nan)

    wp_state = wp_state.sort_values("time_au")
    mask = wp_state["time_au"] <= window.t_end
    wp_win = wp_state[mask]
    if len(wp_win) < 2:
        return (np.nan, np.nan)

    E0 = wp_win["E_expect_ha"].iloc[0]
    E_end = wp_win["E_expect_ha"].iloc[-1]

    dE_ha = E_end - E0
    dE_eV = dE_ha * HA_TO_EV

    dz = params.wp_k0_z * window.t_end
    S2 = -dE_eV / dz if abs(dz) > 0.1 else np.nan

    return (S2, dE_eV)


def compute_classical_S(run_path: Path, window: TimeWindow) -> tuple[float, float]:
    """Compute classical S from linear regression of KE(z).

    Returns: (S_eV_per_bohr, stderr_eV_per_bohr)
    """
    csv_path = run_path / "results" / "raw" / "observables" / "electron_track.csv"
    if not csv_path.exists():
        return (np.nan, np.nan)

    df = pd.read_csv(csv_path)
    df = df.drop_duplicates(subset=["step"], keep="last")
    df = df.dropna(subset=["vx", "vy", "vz"])

    mask = df["time_au"] <= window.t_end
    df_win = df[mask]
    if len(df_win) < 5:
        return (np.nan, np.nan)

    KE_eV = 0.5 * (df_win["vx"]**2 + df_win["vy"]**2 + df_win["vz"]**2) * HA_TO_EV
    z = df_win["z"].values

    result = linregress(z, KE_eV)
    S = -result.slope
    S_err = result.stderr

    return (S, S_err)


# ─── Run registries ───────────────────────────────────────────────────────────

@dataclass
class WPRunSpec:
    energy_eV: float
    sigma: float
    run_dir: str
    is_v2: bool
    L: float = 50.0


@dataclass
class ClassicalRunSpec:
    energy_eV: float
    run_dir: str
    is_v2: bool
    L: float = 50.0


def get_L50_wp_sigma1_runs() -> list[WPRunSpec]:
    return [
        WPRunSpec(20, 1.0, f"{JBASE}/run_wp_n162_L50_E20_sigma1_v2", True),
        WPRunSpec(25, 1.0, f"{JBASE}/run_wp_n162_L50_E25_sigma1_v2", True),
        WPRunSpec(50, 1.0, f"{JBASE}/run_wp_n162_L50_E50_sigma1_v2", True),
        WPRunSpec(100, 1.0, f"{JBASE}/run_wp_n162_L50_E100_sigma1_v2", True),
        WPRunSpec(200, 1.0, f"{JBASE}/run_wp_n162_L50_E200_sigma1_v2", True),
        WPRunSpec(300, 1.0, f"{JBASE}/run_wp_n162_L50_E300_sigma1_v2", True),
    ]


def get_L50_wp_sigma5_runs() -> list[WPRunSpec]:
    return [
        WPRunSpec(20, 5.0, f"{JBASE}/run_wp_n162_L50_E20", False),
        WPRunSpec(25, 5.0, f"{JBASE}/run_wp_n162_L50_E25", False),
        WPRunSpec(50, 5.0, f"{JBASE}/run_wp_n162_L50_E50_v2", True),
        WPRunSpec(100, 5.0, f"{JBASE}/run_wp_n162_L50_E100_v2", True),
        WPRunSpec(300, 5.0, f"{JBASE}/run_wp_n162_L50_E300_v2", True),
        WPRunSpec(600, 5.0, f"{JBASE}/run_wp_n162_L50_E600_v2", True),
    ]


def get_L50_wp_supplementary() -> list[WPRunSpec]:
    return [
        WPRunSpec(100, 0.5, f"{JBASE}/run_wp_n162_L50_E100_sigma0p5", False),
        WPRunSpec(100, 3.0, f"{JBASE}/run_wp_n162_L50_E100_sigma3", False),
        WPRunSpec(100, 8.0, f"{JBASE}/run_wp_n162_L50_E100_sigma8", False),
    ]


def get_L50_classical_runs() -> list[ClassicalRunSpec]:
    return [
        ClassicalRunSpec(20, f"{JBASE}/run_classical_n162_L50_E20", False),
        ClassicalRunSpec(25, f"{JBASE}/run_classical_n162_L50_E25", False),
        ClassicalRunSpec(50, f"{JBASE}/run_classical_n162_L50_E50_v2", True),
        ClassicalRunSpec(100, f"{JBASE}/run_classical_n162_L50_E100_v2", True),
        ClassicalRunSpec(600, f"{JBASE}/run_classical_n162_L50_E600_v2", True),
    ]


def get_L30_wp_sigma1_runs() -> list[WPRunSpec]:
    return [
        WPRunSpec(50, 1.0, f"{JBASE}/run_wp_n162_L30_E50_highdens_sigma1_v2", True, L=30),
        WPRunSpec(100, 1.0, f"{JBASE}/run_wp_n162_L30_E100_highdens_sigma1_v2", True, L=30),
        WPRunSpec(200, 1.0, f"{JBASE}/run_wp_n162_L30_E200_highdens_sigma1_v2", True, L=30),
        WPRunSpec(300, 1.0, f"{JBASE}/run_wp_n162_L30_E300_highdens_sigma1_v2", True, L=30),
    ]


def get_L30_wp_supplementary() -> list[WPRunSpec]:
    return [
        WPRunSpec(100, 0.5, f"{JBASE}/run_wp_n162_L30_E100_highdens", False, L=30),
    ]


def get_L30_classical_runs() -> list[ClassicalRunSpec]:
    return [
        ClassicalRunSpec(50, f"{JBASE}/run_classical_n162_L30_E50_highdens", False, L=30),
        ClassicalRunSpec(100, f"{JBASE}/run_classical_n162_L30_E100_highdens", False, L=30),
        ClassicalRunSpec(200, f"{JBASE}/run_classical_n162_L30_E200_highdens", False, L=30),
        ClassicalRunSpec(300, f"{JBASE}/run_classical_n162_L30_E300_highdens", False, L=30),
    ]


# ─── Master extraction ────────────────────────────────────────────────────────

def extract_wp_stopping(spec: WPRunSpec) -> Optional[StoppingResult]:
    """Extract stopping power from a single WP run."""
    path = Path(spec.run_dir)
    if not path.exists():
        print(f"  SKIP (not found): {spec.run_dir}")
        return None

    params = parse_run_summary(path)
    window = compute_time_window(params)

    S1, dE1, dz, E0, Ef = compute_wp_S1(path, window, params)
    S2, dE2 = compute_wp_S2(path, window, params)

    is_compromised = (spec.sigma >= 8.0)

    result = StoppingResult(
        run_dir=spec.run_dir,
        energy_eV=spec.energy_eV,
        sigma=spec.sigma,
        S1_eV_per_bohr=S1,
        S2_eV_per_bohr=S2,
        window=window,
        is_v2=spec.is_v2,
        is_boundary_compromised=is_compromised,
        dE1_eV=dE1,
        dE2_eV=dE2,
        dz_bohr=dz,
        E_kin_0_eV=E0,
        E_kin_end_eV=Ef,
    )

    _print_wp_provenance(result, window, params)
    return result


def extract_classical_stopping(
    spec: ClassicalRunSpec,
    wp_window: Optional[TimeWindow] = None,
) -> Optional[StoppingResult]:
    """Extract stopping power from a classical run.

    Uses wp_window if provided (matched-window protocol), otherwise
    uses the full run duration.
    """
    path = Path(spec.run_dir)
    if not path.exists():
        print(f"  SKIP (not found): {spec.run_dir}")
        return None

    track_path = path / "results" / "raw" / "observables" / "electron_track.csv"
    if not track_path.exists():
        print(f"  SKIP (no electron_track): {spec.run_dir}")
        return None

    df = pd.read_csv(track_path)
    df = df.drop_duplicates(subset=["step"], keep="last")
    df = df.dropna(subset=["vx", "vy", "vz"])
    t_max_run = df["time_au"].max()

    if wp_window is not None:
        t_end = min(wp_window.t_end, t_max_run)
    else:
        t_end = t_max_run

    window = TimeWindow(t_end=t_end, n_steps_end=int(t_end / 0.01))

    S, S_err = compute_classical_S(path, window)

    z0 = df["z"].iloc[1]  # skip duplicate row 0
    df_win = df[df["time_au"] <= t_end]
    z_end = df_win["z"].iloc[-1]
    dz = z_end - z0

    result = StoppingResult(
        run_dir=spec.run_dir,
        energy_eV=spec.energy_eV,
        sigma=0.0,
        S_classical=S,
        S_classical_err=S_err,
        window=window,
        is_v2=spec.is_v2,
        dz_bohr=dz,
    )

    _print_classical_provenance(result, window, t_max_run)
    return result


def _print_wp_provenance(r: StoppingResult, w: TimeWindow, p: RunParams):
    """Print provenance line for a WP data point."""
    v_tag = "v2" if r.is_v2 else "v1"
    flag = " [BOUNDARY-COMPROMISED]" if r.is_boundary_compromised else ""
    print(
        f"  WP σ={r.sigma:.1f} E={r.energy_eV:.0f}eV ({v_tag}){flag}\n"
        f"    dir: {r.run_dir}\n"
        f"    window: [0, {w.t_end:.2f}] a.u. = [0, {w.t_end * AU_TO_FS:.3f}] fs\n"
        f"    4σ clear at launch: {w.boundary_clear}\n"
        f"    E_kin(0)={r.E_kin_0_eV:.3f} eV  E_kin(end)={r.E_kin_end_eV:.3f} eV\n"
        f"    ΔE₁={r.dE1_eV:.4f} eV  ΔE₂={r.dE2_eV:.4f} eV  Δz={r.dz_bohr:.2f} Bohr\n"
        f"    S₁={r.S1_eV_per_bohr:.5f} eV/Bohr  S₂={r.S2_eV_per_bohr:.5f} eV/Bohr"
    )


def _print_classical_provenance(r: StoppingResult, w: TimeWindow, t_max: float):
    """Print provenance line for a classical data point."""
    v_tag = "v2" if r.is_v2 else "v1"
    trunc = " [TRUNCATED]" if w.t_end < t_max * 0.99 else ""
    print(
        f"  Classical E={r.energy_eV:.0f}eV ({v_tag}){trunc}\n"
        f"    dir: {r.run_dir}\n"
        f"    window: [0, {w.t_end:.2f}] a.u. = [0, {w.t_end * AU_TO_FS:.3f}] fs\n"
        f"    Δz={r.dz_bohr:.2f} Bohr\n"
        f"    S_cl={r.S_classical:.5f} ± {r.S_classical_err:.5f} eV/Bohr"
    )


# ─── Collect all data ─────────────────────────────────────────────────────────

@dataclass
class MasterData:
    """All stopping power data for one density (L=50 or L=30)."""
    wp_sigma1: list[StoppingResult] = field(default_factory=list)
    wp_sigma5: list[StoppingResult] = field(default_factory=list)
    wp_supplementary: list[StoppingResult] = field(default_factory=list)
    classical: list[StoppingResult] = field(default_factory=list)


def collect_L50_data() -> MasterData:
    """Collect all L=50 stopping power data."""
    print("=" * 60)
    print("L=50 standard density (r_s ≈ 5.69)")
    print("=" * 60)

    data = MasterData()

    print("\n── WP σ=1 energy sweep ──")
    wp_s1_results = []
    for spec in get_L50_wp_sigma1_runs():
        r = extract_wp_stopping(spec)
        if r is not None:
            wp_s1_results.append(r)
    data.wp_sigma1 = wp_s1_results

    print("\n── WP σ=5 energy sweep ──")
    wp_s5_results = []
    for spec in get_L50_wp_sigma5_runs():
        r = extract_wp_stopping(spec)
        if r is not None:
            wp_s5_results.append(r)
    data.wp_sigma5 = wp_s5_results

    print("\n── WP supplementary (E=100, various σ) ──")
    wp_supp = []
    for spec in get_L50_wp_supplementary():
        r = extract_wp_stopping(spec)
        if r is not None:
            wp_supp.append(r)
    data.wp_supplementary = wp_supp

    print("\n── Classical Ehrenfest ──")
    # Build window lookup from σ=1 runs (preferred) then σ=5
    window_by_energy: dict[float, TimeWindow] = {}
    for r in data.wp_sigma1:
        if r.window is not None:
            window_by_energy[r.energy_eV] = r.window
    for r in data.wp_sigma5:
        if r.window is not None and r.energy_eV not in window_by_energy:
            window_by_energy[r.energy_eV] = r.window

    cl_results = []
    for spec in get_L50_classical_runs():
        wp_win = window_by_energy.get(spec.energy_eV)
        r = extract_classical_stopping(spec, wp_win)
        if r is not None:
            cl_results.append(r)
    data.classical = cl_results

    return data


def collect_L30_data() -> MasterData:
    """Collect all L=30 stopping power data."""
    print("\n" + "=" * 60)
    print("L=30 high density (r_s ≈ 3.41)")
    print("=" * 60)

    data = MasterData()

    print("\n── WP σ=1 energy sweep ──")
    wp_s1_results = []
    for spec in get_L30_wp_sigma1_runs():
        r = extract_wp_stopping(spec)
        if r is not None:
            wp_s1_results.append(r)
    data.wp_sigma1 = wp_s1_results

    print("\n── WP supplementary (σ=0.5) ──")
    wp_supp = []
    for spec in get_L30_wp_supplementary():
        r = extract_wp_stopping(spec)
        if r is not None:
            wp_supp.append(r)
    data.wp_supplementary = wp_supp

    print("\n── Classical Ehrenfest ──")
    window_by_energy: dict[float, TimeWindow] = {}
    for r in data.wp_sigma1:
        if r.window is not None:
            window_by_energy[r.energy_eV] = r.window

    cl_results = []
    for spec in get_L30_classical_runs():
        wp_win = window_by_energy.get(spec.energy_eV)
        r = extract_classical_stopping(spec, wp_win)
        if r is not None:
            cl_results.append(r)
    data.classical = cl_results

    return data


# ─── Bethe–Bloch reference ────────────────────────────────────────────────────

def bethe_curve(E_eV: np.ndarray, n_el: float) -> np.ndarray:
    """Bethe stopping: S = (4π n / v²) ln(2v² / ω_p) in eV/Bohr."""
    v = np.sqrt(2 * E_eV / HA_TO_EV)
    omega_p = np.sqrt(4 * np.pi * n_el)
    L_arg = 2 * v**2 / omega_p
    S = 4 * np.pi * n_el / v**2 * np.log(L_arg) * HA_TO_EV
    return np.maximum(S, 0)


if __name__ == "__main__":
    data_L50 = collect_L50_data()
    data_L30 = collect_L30_data()

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"L=50: {len(data_L50.wp_sigma1)} σ=1 pts, "
          f"{len(data_L50.wp_sigma5)} σ=5 pts, "
          f"{len(data_L50.wp_supplementary)} supp pts, "
          f"{len(data_L50.classical)} classical pts")
    print(f"L=30: {len(data_L30.wp_sigma1)} σ=1 pts, "
          f"{len(data_L30.wp_supplementary)} supp pts, "
          f"{len(data_L30.classical)} classical pts")
