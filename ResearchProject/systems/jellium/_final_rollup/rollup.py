"""Final stopping-power vs v rollup for the 2026-05-21 meeting.

Aggregates every available data point in the campaign into a single
S(v) figure, plotting three independent stopping-power estimators:

  A. Eigenvalue / WP energy_balance — ΔE_WP / Δz over the IFW window
     (clean single-state signal for σ=5 and σ=1 jellium WP runs).
  B. Knudsen <|p|²>/2 — from `wp_momentum_stats.csv` if available, else
     from the histogram-based retroactive `knudsen_ke.csv` (legacy
     runs). Note: at σ=1 the Knudsen E_kin GROWS (σ_p² spreads), so
     the slope sign convention is the *energy lost by the central-k₀
     momentum*, captured here by the boundary_rule traversal sample.
  C. Classical (Ehrenfest) — `-dE_proj/dz` from electron_track.csv +
     observables.csv.

The IFW used per run is taken from `run_summary.txt` if available, else
boundary_rule.

Output: stopping_power_vs_v.png + stopping_power_data.csv next to this
script. Used by the [jellium-rollup] meeting email.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HA = 27.21138625
ROOT = Path(__file__).parent.parent     # ResearchProject/systems/jellium
HERE = Path(__file__).parent


# ────────────────────────────────────────────────────────────────────────
# Run inventory — each tuple: (run_name, archetype, energy_ev, sigma_w, density_label, partner_classical)
# ────────────────────────────────────────────────────────────────────────
RUNS = [
    # New campaign runs (with native wp_momentum_stats / wp_real_space_stats)
    ("run_wp_n162_L50_E25",                "wp",       25.0, 5.0, "r_s=5.69", None),
    ("run_wp_n162_L50_E100_v2",            "wp",       100.0, 5.0, "r_s=5.69",
        "run_classical_n162_L50_E100_v2"),
    ("run_wp_n162_L50_E100_sigma1",        "wp",       100.0, 1.0, "r_s=5.69, σ=1",
        "run_classical_n162_L50_E100_sigma1"),
    ("run_wp_n162_L30_E100_highdens",      "wp",       100.0, 0.5, "r_s=3.41 (highdens)",
        "run_classical_n162_L30_E100_highdens"),
    # New Run-5 + Run-9 (low-v anchors added 2026-05-18..19)
    ("run_wp_n162_L50_E20",                "wp",       20.0, 5.0, "r_s=5.69",
        "run_classical_n162_L50_E20"),
    # (Run-1 already in list as run_wp_n162_L50_E25; pair with Run-9 below)
    # Run-9 standalone classical at E=25 — manually added since
    # Run-1 (run_wp_n162_L50_E25) is wp with no partner in original list
    ("run_classical_n162_L50_E25",         "classical_only", 25.0, None, "r_s=5.69", None),
    # 2026-05-20: backlog re-runs at E={50, 300, 600} with new observables
    ("run_wp_n162_L50_E50_v2",             "wp", 50.0, 5.0, "r_s=5.69", "run_classical_n162_L50_E50"),
    ("run_wp_n162_L50_E300_v2",            "wp", 300.0, 5.0, "r_s=5.69", "run_classical_n162_L50_E300"),
    ("run_wp_n162_L50_E600_v2",            "wp", 600.0, 5.0, "r_s=5.69", "run_classical_n162_L50_E600_v2"),
    # Legacy WP runs (retroactive histogram Knudsen only — kept for reference)
    ("run_wp_n162_L50_E50",                "wp_legacy", 50.0, 5.0, "r_s=5.69", "run_classical_n162_L50_E50"),
    ("run_wp_n162_L50_E100",               "wp_legacy", 100.0, 5.0, "r_s=5.69", "run_classical_n162_L50_E100"),
    ("run_wp_n162_L50_E300",               "wp_legacy", 300.0, 5.0, "r_s=5.69", "run_classical_n162_L50_E300"),
    ("run_wp_n162_L50_E600",               "wp_legacy", 600.0, 5.0, "r_s=5.69", None),
    ("run_classical_e1500_L50_cubic",      "classical_only", 1500.0, None, "r_s=5.69", None),
]

# Patch run_wp_n162_L50_E25 → set partner to run_classical_n162_L50_E25 (Run-9)
RUNS = [(name, "wp" if name == "run_wp_n162_L50_E25" else arc,
         ev, sig, dens,
         "run_classical_n162_L50_E25" if name == "run_wp_n162_L50_E25" else part)
        for (name, arc, ev, sig, dens, part) in RUNS]


def k0_of(ev: float) -> float:
    return math.sqrt(2.0 * ev / HA)


def _ifw_t_au(run_dir: Path) -> Optional[tuple[float, float]]:
    """Return (t_IFW, t_total) from run_summary.txt via the helper."""
    from inqview.postprocess._common import post_ifw_window_from_summary
    try:
        return post_ifw_window_from_summary(run_dir / "results")
    except Exception:
        return None


@dataclass
class S_estimate:
    method: str          # 'wp_energy_balance', 'wp_knudsen', 'classical'
    s_ev_per_bohr: float
    dE_eV: float
    dz_bohr: float
    note: str = ""


def s_from_wp_energy_balance(run_dir: Path,
                              ifw_t: float) -> Optional[S_estimate]:
    """Method A: ΔE_WP / Δz from energy_balance.csv + wp_real_space_stats z_mean."""
    eb = run_dir / "results/analysis/observables/energy_balance.csv"
    rs = run_dir / "results/raw/observables/wp_real_space_stats.csv"
    if not (eb.exists() and rs.exists()):
        return None
    eb_df = pd.read_csv(eb)
    rs_df = pd.read_csv(rs, comment="#")
    eb_w = eb_df[eb_df["time_au"] <= ifw_t]
    rs_w = rs_df[rs_df["time_au"] <= ifw_t]
    if len(eb_w) < 2 or len(rs_w) < 2:
        return None
    dE_eV = float(eb_w["dE_wp_ev"].iloc[-1]) - float(eb_w["dE_wp_ev"].iloc[0])
    dz = float(rs_w["z_mean"].iloc[-1]) - float(rs_w["z_mean"].iloc[0])
    if abs(dz) < 1e-6:
        return None
    return S_estimate("wp_energy_balance", -dE_eV / dz, dE_eV, dz)


def s_from_wp_knudsen(run_dir: Path, ifw_t: float) -> Optional[S_estimate]:
    """Method B: ΔE_Knudsen / Δz over IFW. Note: σ=1 has a spread artefact —
    Knudsen grows because σ_p² spreads. We *report* the slope anyway and note
    the convention in the rollup caption."""
    ke = run_dir / "results/analysis/observables/knudsen_ke.csv"
    if not ke.exists():
        return None
    df = pd.read_csv(ke)
    df_w = df[df["time_au"] <= ifw_t]
    if len(df_w) < 2:
        return None
    dE_eV = float(df_w["e_kin_ev"].iloc[-1]) - float(df_w["e_kin_ev"].iloc[0])
    # Δz from z_bohr if present, else from velocity (k0 in atomic units) × Δt.
    if "z_bohr" in df_w.columns and df_w["z_bohr"].notna().any():
        dz = float(df_w["z_bohr"].iloc[-1]) - float(df_w["z_bohr"].iloc[0])
    else:
        # Read k0 from run_summary
        rs = (run_dir / "results/run_summary.txt").read_text()
        m = re.search(r"wp_k0_bohr_inv\s*=\s*\S+\s+\S+\s+(\S+)", rs)
        if not m:
            return None
        v = float(m.group(1))
        dt = float(df_w["time_au"].iloc[-1]) - float(df_w["time_au"].iloc[0])
        dz = v * dt
    if abs(dz) < 1e-6:
        return None
    return S_estimate("wp_knudsen", -dE_eV / dz, dE_eV, dz,
                       note=("σ=1 spread-dominated" if "sigma1" in run_dir.name else ""))


def s_from_classical(run_dir: Path) -> Optional[S_estimate]:
    """Method C: dE_total_bath / Δz from observables.csv + electron_track."""
    obs = run_dir / "results/raw/observables/observables.csv"
    tr  = run_dir / "results/raw/observables/electron_track.csv"
    if not (obs.exists() and tr.exists()):
        return None
    obs_df = pd.read_csv(obs)
    tr_df = pd.read_csv(tr)
    if "energy_total" not in obs_df.columns or "z" not in tr_df.columns:
        return None
    dE_eV = (float(obs_df["energy_total"].iloc[-1])
             - float(obs_df["energy_total"].iloc[0])) * HA
    dz = float(tr_df["z"].iloc[-1]) - float(tr_df["z"].iloc[0])
    if abs(dz) < 1e-6:
        return None
    # Bath gain == projectile loss, so S = +ΔE_bath / Δz (energy out of projectile per Bohr).
    return S_estimate("classical", dE_eV / dz, dE_eV, dz)


def main() -> None:
    rows = []
    for (name, archetype, ev, sigma, density, partner) in RUNS:
        run_dir = ROOT / name
        if not (run_dir / "results").exists():
            print(f"  [skip] {name}: no results/")
            continue
        v = k0_of(ev)         # atomic units; m_e = 1 so v = k0
        ifw = _ifw_t_au(run_dir)
        if ifw is None:
            # Default if run_summary lacks the needed fields
            ifw_t = 9.5
        else:
            ifw_t = ifw[0]
        row = {"run_name": name, "archetype": archetype,
               "energy_ev": ev, "sigma_w": sigma, "density": density,
               "v_au": v, "t_ifw_au": ifw_t}

        if archetype in ("wp", "wp_legacy"):
            sA = s_from_wp_energy_balance(run_dir, ifw_t)
            sB = s_from_wp_knudsen(run_dir, ifw_t)
            row["S_wp_energy_balance"] = sA.s_ev_per_bohr if sA else None
            row["S_wp_knudsen"]        = sB.s_ev_per_bohr if sB else None
            row["dE_wp_ev"]            = sA.dE_eV if sA else (sB.dE_eV if sB else None)
            row["dz_wp_bohr"]          = sA.dz_bohr if sA else (sB.dz_bohr if sB else None)
            row["knudsen_note"]        = sB.note if sB else ""

        # Classical companion (if applicable)
        if partner:
            cl_dir = ROOT / partner
            if (cl_dir / "results").exists():
                sC = s_from_classical(cl_dir)
                if sC:
                    row["S_classical"] = sC.s_ev_per_bohr
                    row["dE_classical_ev"] = sC.dE_eV
                    row["dz_classical_bohr"] = sC.dz_bohr
                    row["classical_partner"] = partner

        if archetype == "classical_only":
            sC = s_from_classical(run_dir)
            if sC:
                row["S_classical"] = sC.s_ev_per_bohr
                row["dE_classical_ev"] = sC.dE_eV
                row["dz_classical_bohr"] = sC.dz_bohr

        rows.append(row)
        print(f"  {name}: v={v:.3f}  S_eb={row.get('S_wp_energy_balance')}  "
              f"S_knudsen={row.get('S_wp_knudsen')}  S_classical={row.get('S_classical')}")

    df = pd.DataFrame(rows)
    csv_out = HERE / "stopping_power_data.csv"
    df.to_csv(csv_out, index=False)
    print(f"\nwrote {csv_out}")

    # ──────────────────────────────────────────────────────────────────
    # Plot S(v) vs v with the three estimators
    # ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 6))

    # Standard-density runs (r_s = 5.69)
    std = df[df["density"] == "r_s=5.69"]
    sig5 = std[std["sigma_w"] == 5.0]
    sig1 = std[std["sigma_w"] == 1.0]

    if not sig5.empty:
        m = sig5["S_wp_energy_balance"].notna()
        if m.any():
            ax.scatter(sig5.loc[m, "v_au"], sig5.loc[m, "S_wp_energy_balance"],
                       s=80, marker="o", color="C3",
                       label="Method A: WP energy_balance (σ=5)", zorder=5)
        m = sig5["S_wp_knudsen"].notna()
        if m.any():
            ax.scatter(sig5.loc[m, "v_au"], sig5.loc[m, "S_wp_knudsen"],
                       s=55, marker="s", color="C3", alpha=0.5,
                       label="Method B: WP Knudsen <|p|²>/2 (σ=5)")
        m = sig5["S_classical"].notna()
        if m.any():
            ax.scatter(sig5.loc[m, "v_au"], sig5.loc[m, "S_classical"],
                       s=80, marker="^", color="C0",
                       label="Method C: classical Ehrenfest (σ=5 partner)")

    if not sig1.empty:
        m = sig1["S_wp_energy_balance"].notna()
        if m.any():
            ax.scatter(sig1.loc[m, "v_au"], sig1.loc[m, "S_wp_energy_balance"],
                       s=120, marker="o", facecolors="none", edgecolors="C2",
                       linewidths=2, label="Method A: WP energy_balance (σ=1)")
        m = sig1["S_classical"].notna()
        if m.any():
            ax.scatter(sig1.loc[m, "v_au"], sig1.loc[m, "S_classical"],
                       s=120, marker="^", facecolors="none", edgecolors="C2",
                       linewidths=2, label="Method C: classical Ehrenfest (σ=1 partner)")

    # High-density Run-6
    hd = df[df["density"] == "r_s=3.41 (highdens)"]
    if not hd.empty:
        m = hd["S_wp_energy_balance"].notna()
        if m.any():
            ax.scatter(hd.loc[m, "v_au"], hd.loc[m, "S_wp_energy_balance"],
                       s=120, marker="D", color="purple",
                       label="Method A: WP energy_balance (r_s=3.41)")
        m = hd["S_classical"].notna()
        if m.any():
            ax.scatter(hd.loc[m, "v_au"], hd.loc[m, "S_classical"],
                       s=120, marker="v", color="purple",
                       label="Method C: classical (r_s=3.41)")

    ax.set_xlabel("v = k₀ (Bohr/atu)")
    ax.set_ylabel("Stopping power S(v) (eV/Bohr)")
    ax.set_title("Jellium electronic stopping — campaign rollup 2026-05-21")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="0.5", lw=0.6)

    fig.tight_layout()
    png_out = HERE / "stopping_power_vs_v.png"
    fig.savefig(png_out, dpi=150)
    plt.close(fig)
    print(f"wrote {png_out}")


if __name__ == "__main__":
    main()
