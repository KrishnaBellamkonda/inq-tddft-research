#!/usr/bin/env python3
"""Plot stopping power S(v) vs projectile energy for ALL jellium runs.

Compares classical Ehrenfest and WP runs across densities and sigma values.
For classical: S = KE_loss / distance (from REPORT.md).
For WP: S_eff = |ΔE_kinetic| / traversal (effective stopping power).

Output: stopping_power_vs_energy_all.png
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

HA_TO_EV = 27.211386245988

JELLIUM_DIR = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium")


def read_summary(rs_path: Path) -> dict:
    if not rs_path.exists():
        return {}
    text = rs_path.read_text()
    d = {}
    for line in text.splitlines():
        m = re.match(r"\s*(\w+)\s*=\s*(.+)", line)
        if m:
            key, val = m.group(1).strip(), m.group(2).strip()
            try:
                d[key] = float(val)
            except ValueError:
                d[key] = val
    return d


def extract_dE_kinetic(obs_csv: Path) -> float | None:
    if not obs_csv.exists():
        return None
    obs = pd.read_csv(obs_csv)
    if "energy_kinetic" not in obs.columns:
        return None
    dE = (obs["energy_kinetic"].iloc[-1] - obs["energy_kinetic"].iloc[0]) * HA_TO_EV
    return float(dE)


def extract_classical_S(report_md: Path) -> float | None:
    if not report_md.exists():
        return None
    text = report_md.read_text()
    m = re.search(r"Stopping power S\(v\)\*\*\s*\|\s*\*\*([\d.]+)", text)
    if m:
        return float(m.group(1))
    return None


def extract_classical_KE_loss_and_distance(report_md: Path) -> tuple[float, float] | None:
    if not report_md.exists():
        return None
    text = report_md.read_text()
    m_ke = re.search(r"KE loss\*\*\s*\|\s*\*\*[+-]?([\d.]+)", text)
    m_dist = re.search(r"Distance traveled\s*\|\s*([\d.]+)", text)
    if m_ke and m_dist:
        return float(m_ke.group(1)), float(m_dist.group(1))
    return None


def get_traversal_distance(sigma: float, L: float, v: float, dt: float = 0.02) -> float:
    """Compute effective traversal from self-spread cap or boundary rule."""
    launch_z = -L / 2 + 4 * sigma
    stop_z = L / 2 - sigma

    if sigma <= 2.0:
        import scipy.optimize as opt
        def spread_eq(t):
            sigma_dens = (sigma / np.sqrt(2)) * np.sqrt(1 + (t / sigma**2)**2)
            return launch_z + v * t + 3 * sigma_dens - L / 2
        try:
            t_star = opt.brentq(spread_eq, 0.1, 100.0)
        except ValueError:
            t_star = (stop_z - launch_z) / v
        z_end = launch_z + v * t_star
        return z_end - launch_z
    else:
        return stop_z - launch_z


def collect_data() -> pd.DataFrame:
    rows = []

    for d in sorted(JELLIUM_DIR.glob("run_*")):
        if not d.is_dir():
            continue
        rs_path = d / "results" / "run_summary.txt"
        rs = read_summary(rs_path)
        if not rs:
            continue
        if str(rs.get("run_completed", "")).lower() != "true":
            continue

        name = d.name
        is_classical = "classical" in name
        is_free = "free" in name
        is_plasmon = "plasmon" in name
        is_base_legacy = name in ("run_base_n138", "run_base_n162_L50_E1p5",
                                   "run_base_n138_L30_E5")
        if is_free or is_plasmon or is_base_legacy:
            continue
        # Skip known-bad runs (GPU hang, incomplete propagation)
        if "E700" in name or "E800" in name or "E900" in name or \
           "E1000" in name or "E1100" in name:
            continue

        # Energy
        E = rs.get("wp_energy_ev") or rs.get("wp_ekin_ev") or rs.get("projectile_KE_eV")
        if E is None:
            E_match = re.search(r"[Ee](\d+)", name)
            if E_match:
                E = float(E_match.group(1))
        if E is None:
            continue
        E = float(E)

        # Sigma
        sigma = rs.get("wp_sigma_bohr")
        if sigma is not None:
            sigma = float(sigma)

        # Cell size
        L_str = str(rs.get("cell_bohr", ""))
        L_match = re.search(r"([\d.]+)", L_str)
        L = float(L_match.group(1)) if L_match else None

        if L is None:
            continue

        density_class = "high" if L < 40 else "standard"

        obs_csv = d / "results" / "raw" / "observables" / "observables.csv"
        report_md = d / "results" / "analysis" / "REPORT.md"

        if is_classical:
            S = extract_classical_S(report_md)
            if S is None:
                result = extract_classical_KE_loss_and_distance(report_md)
                if result:
                    ke_loss, dist = result
                    S = ke_loss / dist if dist > 0 else None
            if S is not None:
                rows.append({
                    "name": name, "type": "classical", "E_eV": E,
                    "sigma": None, "L": L, "density": density_class,
                    "S_eV_per_Bohr": S,
                })
        else:
            dE_kin = extract_dE_kinetic(obs_csv)
            if dE_kin is None:
                continue

            if sigma is None:
                continue

            v = np.sqrt(2 * E / HA_TO_EV)
            traversal = get_traversal_distance(sigma, L, v)
            S_eff = abs(dE_kin) / traversal if traversal > 0 else None

            if S_eff is not None:
                rows.append({
                    "name": name, "type": "wp", "E_eV": E,
                    "sigma": sigma, "L": L, "density": density_class,
                    "S_eV_per_Bohr": S_eff,
                    "dE_kin_eV": dE_kin,
                })

    return pd.DataFrame(rows)


def deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    """When multiple runs exist at the same (type, E, sigma, density), keep
    the one with the latest directory name (v2 > v1 > base)."""
    if df.empty:
        return df
    key_cols = ["type", "E_eV", "density"]
    if "sigma" in df.columns:
        key_cols.append("sigma")
    return df.sort_values("name").groupby(key_cols, dropna=False).last().reset_index()


def main():
    df = collect_data()
    if df.empty:
        print("No data found!")
        return 1

    df = deduplicate(df)

    fig, ax = plt.subplots(figsize=(12, 7))

    # --- Classical series ---
    cl_std = df[(df["type"] == "classical") & (df["density"] == "standard")]
    cl_hd = df[(df["type"] == "classical") & (df["density"] == "high")]

    if not cl_std.empty:
        cl_std = cl_std.sort_values("E_eV")
        ax.plot(cl_std["E_eV"], cl_std["S_eV_per_Bohr"],
                "k-o", ms=7, lw=2, zorder=10,
                label=r"Classical, $r_s \approx 5.69$ (L=50)")

    if not cl_hd.empty:
        cl_hd = cl_hd.sort_values("E_eV")
        ax.plot(cl_hd["E_eV"], cl_hd["S_eV_per_Bohr"],
                "k--s", ms=7, lw=2, zorder=10,
                label=r"Classical, $r_s \approx 3.41$ (L=30)")

    # --- WP series: group by (sigma, density) ---
    wp = df[df["type"] == "wp"].copy()

    colors = {0.5: "C0", 1.0: "C1", 3.0: "C2", 5.0: "C3", 8.0: "C4"}
    markers_std = {0.5: "v", 1.0: "D", 3.0: "^", 5.0: "o", 8.0: "p"}
    markers_hd = {0.5: "v", 1.0: "D", 3.0: "^", 5.0: "o", 8.0: "p"}

    for (sigma, density), group in wp.groupby(["sigma", "density"]):
        group = group.sort_values("E_eV")
        color = colors.get(sigma, "C5")
        filled = density == "standard"
        marker = markers_std.get(sigma, "x") if filled else markers_hd.get(sigma, "x")
        facecolor = color if filled else "none"
        density_label = r"$r_s \approx 5.69$" if density == "standard" else r"$r_s \approx 3.41$"

        label = rf"WP $\sigma={sigma}$ Bohr, {density_label}"

        ax.plot(group["E_eV"], group["S_eV_per_Bohr"],
                color=color, marker=marker, ms=8, lw=1.5,
                markerfacecolor=facecolor, markeredgecolor=color,
                markeredgewidth=1.5, linestyle="-" if filled else "--",
                label=label, zorder=5)

    # --- Bethe-Bloch reference line (v^-2 scaling) ---
    E_ref = np.logspace(np.log10(15), np.log10(2000), 100)
    S_bethe = cl_std["S_eV_per_Bohr"].iloc[len(cl_std) // 2] * \
              (cl_std["E_eV"].iloc[len(cl_std) // 2] / E_ref) if not cl_std.empty else None
    if S_bethe is not None:
        ax.plot(E_ref, S_bethe, ":", color="gray", lw=1, alpha=0.5,
                label=r"$\propto E^{-1}$ (Bethe guide)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Projectile kinetic energy / eV", fontsize=13)
    ax.set_ylabel(r"Stopping power $S$ / (eV / Bohr)", fontsize=13)
    ax.set_title("Stopping power vs projectile energy — all jellium runs", fontsize=14)

    ax.set_xlim(15, 2000)
    y_vals = df["S_eV_per_Bohr"].dropna()
    if not y_vals.empty:
        ax.set_ylim(y_vals.min() * 0.5, y_vals.max() * 2)

    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xticks([20, 50, 100, 200, 300, 500, 1000, 1500])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())

    ax.legend(fontsize=9, loc="upper right", framealpha=0.9,
              ncol=1, handlelength=2.5)
    ax.grid(True, which="both", alpha=0.2)
    fig.tight_layout()

    out = JELLIUM_DIR / "stopping_power_vs_energy_all.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved: {out}")

    # Print data table
    print("\n=== DATA TABLE ===")
    print(f"{'Name':<55} {'Type':<10} {'E(eV)':<8} {'σ':<6} {'L':<5} {'ρ':<10} {'S(eV/Bohr)':<12}")
    print("-" * 110)
    for _, r in df.sort_values(["type", "density", "sigma", "E_eV"]).iterrows():
        sigma_str = f"{r['sigma']:.1f}" if pd.notna(r.get("sigma")) else "n/a"
        print(f"{r['name']:<55} {r['type']:<10} {r['E_eV']:<8.0f} {sigma_str:<6} "
              f"{r['L']:<5.0f} {r['density']:<10} {r['S_eV_per_Bohr']:<12.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
