"""Phase: ``knudsen_ke`` — projectile kinetic-energy stopping-power (Method B).

Computes the **Knudsen-style** stopping-power estimator

    E_kin_WP(t) = <|p|^2>_WP(t) / 2          (Hartree)

i.e. the kinetic energy of the wave-packet orbital extracted from its
momentum-space distribution. This is the second of the two independent
stopping-power estimators used in the 2026-05-21 meeting campaign (the
first being the eigenvalue ``epsilon_WP(t)`` trajectory in
``stopping.py``); see ``docs/plans/jellium-meeting-2026-05-21.md``
§"Stopping-power estimators" and Knudsen et al., arXiv 2605.12854,
*Ultrafast electron dynamics of electron-irradiated graphene*.

Two input paths, in priority order:

1. **Native (new runs):** ``raw/observables/wp_momentum_stats.csv``
   written by ``inqkit::observables::WPMomentumStats``. Column
   ``e_kin_ha`` is already ``<|p|^2>/2`` for the WP orbital. Direct.

2. **Retroactive (existing runs):** ``raw/observables/momentum_distribution.csv``
   written by ``inqkit::observables::MomentumDistribution``. Column
   ``n_wp`` is a |k|-binned histogram for the WP only; we approximate

       <|p|^2>_WP(t) ~= sum_bins k_bin^2 * n_wp(bin, t) / sum_bins n_wp(bin, t)

   This is degraded by the 1D radial binning but gives a usable curve
   for the E={50, 300, 600, 1500} eV runs that pre-date Infra-4.

The trajectory z(t) is read from ``observables.csv``
(``cod_z_bohr`` column for WP runs, matching the column name written by
the jellium run-template — see ``wp_trajectory.py``). Outputs:

* ``analysis/observables/knudsen_ke.csv``
   — long-format ``step,time_au,z_bohr,e_kin_ha,e_kin_ev`` with the
   value of ``e_kin_ha`` derived per the priority above. Two extra
   columns ``de_kin_ev`` and ``stopping_power_ev_per_bohr`` give the
   shifted-to-zero kinetic energy and the running ``-dE_kin/dz`` slope
   (one-sided finite difference, sign chosen so positive means the WP
   slowed down).
* ``analysis/observables/knudsen_ke_vs_t.png``
* ``analysis/observables/knudsen_ke_vs_z.png`` (only if z is available)
* ``analysis/observables/knudsen_stopping_power_vs_z.png``

Phase is skipped cleanly if both source CSVs are missing.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import _common


HA_TO_EV = 27.21138625


def _from_wp_momentum_stats(csv_path: Path) -> pd.DataFrame:
    """Native path: e_kin_ha column is already Method B."""
    df = pd.read_csv(csv_path, comment="#")
    out = pd.DataFrame({
        "step":      df["step"].astype(int),
        "time_au":   df["time_au"].astype(float),
        "e_kin_ha":  df["e_kin_ha"].astype(float),
    })
    out["source"] = "wp_momentum_stats"
    return out


def _from_momentum_distribution(csv_path: Path) -> pd.DataFrame:
    """Retroactive degraded path: histogram-weighted <|k|^2>/2.

    The CSV is long-format with one row per (step, |k|-bin). We compute
    per step:
        N(t)        = sum_bin n_wp(bin, t)
        <|k|^2>(t)  = sum_bin k_bin^2 * n_wp(bin, t) / N(t)
        e_kin_ha    = <|k|^2>(t) / 2
    """
    df = pd.read_csv(csv_path, comment="#")
    # n_wp column may be 0 at every bin if the WP was not tagged with the
    # right state index in the run — guard against that.
    if "n_wp" not in df.columns or df["n_wp"].abs().sum() == 0:
        raise ValueError("momentum_distribution.csv has no n_wp data "
                         "(check wp_idx in the MomentumDistribution writer)")
    grp = df.groupby(["step", "time_au"], sort=True)
    rows = []
    for (step, t), g in grp:
        k = g["k_bohr_inv"].to_numpy(dtype=float)
        w = g["n_wp"].to_numpy(dtype=float)
        n = w.sum()
        if n <= 0:
            continue
        k2_mean = float((k * k * w).sum() / n)
        rows.append((int(step), float(t), 0.5 * k2_mean))
    out = pd.DataFrame(rows, columns=["step", "time_au", "e_kin_ha"])
    out["source"] = "momentum_distribution_histogram"
    return out


def _merge_trajectory(ke_df: pd.DataFrame, results_dir: Path) -> pd.DataFrame:
    """Attach z(t) = cod_z_bohr from observables.csv when available."""
    obs = results_dir / "raw" / "observables" / "observables.csv"
    if not obs.exists():
        ke_df["z_bohr"] = np.nan
        return ke_df
    o = pd.read_csv(obs)
    # cod_z_bohr is the canonical column name written by the jellium
    # run-template; older runs may use cod_z. Accept either.
    z_col = next((c for c in ("cod_z_bohr", "cod_z") if c in o.columns), None)
    if z_col is None:
        ke_df["z_bohr"] = np.nan
        return ke_df
    merged = pd.merge(ke_df, o[["step", z_col]], on="step", how="left")
    merged = merged.rename(columns={z_col: "z_bohr"})
    return merged


def _compute_stopping_power(df: pd.DataFrame) -> pd.DataFrame:
    """One-sided ``-d(E_kin) / dz`` slope. Positive = WP slowing down.

    Uses forward differences in z (skipping NaN); ill-defined where dz
    is small or non-monotonic, but adequate for an unsmoothed first cut.
    """
    out = df.copy()
    out["de_kin_ev"] = (out["e_kin_ha"] - out["e_kin_ha"].iloc[0]) * HA_TO_EV
    z  = out["z_bohr"].to_numpy()
    e  = out["e_kin_ha"].to_numpy() * HA_TO_EV
    sp = np.full(len(out), np.nan, dtype=float)
    for i in range(len(out) - 1):
        dz = z[i+1] - z[i]
        if np.isfinite(dz) and abs(dz) > 1e-6:
            sp[i] = -(e[i+1] - e[i]) / dz   # eV / Bohr; +ve = slowdown
    out["stopping_power_ev_per_bohr"] = sp
    return out


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    obs_dir = results_dir / "raw" / "observables"
    wpms = obs_dir / "wp_momentum_stats.csv"
    md   = obs_dir / "momentum_distribution.csv"

    if wpms.exists():
        ke = _from_wp_momentum_stats(wpms)
        source = "wp_momentum_stats.csv (native)"
    elif md.exists():
        ke = _from_momentum_distribution(md)
        source = "momentum_distribution.csv (retroactive histogram)"
    else:
        return {"skipped": "neither wp_momentum_stats.csv nor "
                           "momentum_distribution.csv present"}

    ke = _merge_trajectory(ke, results_dir)
    ke = _compute_stopping_power(ke)
    ke["e_kin_ev"] = ke["e_kin_ha"] * HA_TO_EV

    out_dir = _common.ensure_dir(results_dir / "analysis" / "observables")
    artefacts: list[str] = []

    # Post-IFW shading on time-axis plots (campaign rule §4).
    ifw = _common.post_ifw_window_from_summary(results_dir)

    csv_out = out_dir / "knudsen_ke.csv"
    if _common.need_rebuild(csv_out, rebuild):
        ke.to_csv(csv_out, index=False)
    artefacts.append(str(csv_out))

    # ---- E_kin(t) -------------------------------------------------------
    p_t = out_dir / "knudsen_ke_vs_t.png"
    if _common.need_rebuild(p_t, rebuild):
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(ke["time_au"], ke["e_kin_ev"], "C1-", lw=1.4,
                label="$E_{\\rm kin,WP}$ Method-B")
        if ifw is not None:
            _common.ifw_highlight(ax, ifw[0])
            ax.legend(loc="best", fontsize=9)
        ax.set_xlabel("time (a.u.)")
        ax.set_ylabel("E_kin_WP = <|p|²>/2  (eV)")
        ax.set_title(_common.title(run_name,
            f"Knudsen E_kin_WP(t) — Method B [{source}]"))
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(p_t, dpi=150)
        plt.close(fig)
        artefacts.append(str(p_t))

    # ---- E_kin(z), only if z is finite ----------------------------------
    has_z = ke["z_bohr"].notna().any()
    if has_z:
        p_z = out_dir / "knudsen_ke_vs_z.png"
        if _common.need_rebuild(p_z, rebuild):
            fig, ax = plt.subplots(figsize=(8, 4.5))
            ax.plot(ke["z_bohr"], ke["e_kin_ev"], "C1-", lw=1.4)
            ax.set_xlabel("WP cod_z (Bohr)")
            ax.set_ylabel("E_kin_WP = <|p|²>/2  (eV)")
            ax.set_title(_common.title(run_name,
                f"Knudsen E_kin_WP(z) — Method B [{source}]"))
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            fig.savefig(p_z, dpi=150)
            plt.close(fig)
            artefacts.append(str(p_z))

        p_sp = out_dir / "knudsen_stopping_power_vs_z.png"
        if _common.need_rebuild(p_sp, rebuild):
            fig, ax = plt.subplots(figsize=(8, 4.5))
            ax.plot(ke["z_bohr"], ke["stopping_power_ev_per_bohr"],
                    "C3.", ms=3, lw=0)
            ax.axhline(0.0, color="0.5", lw=0.6)
            ax.set_xlabel("WP cod_z (Bohr)")
            ax.set_ylabel("-dE_kin_WP / dz  (eV / Bohr)")
            ax.set_title(_common.title(run_name,
                f"Knudsen stopping-power proxy — Method B [{source}]"))
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            fig.savefig(p_sp, dpi=150)
            plt.close(fig)
            artefacts.append(str(p_sp))

    return {"source": source, "n_steps": int(len(ke)),
            "artefacts": artefacts}
