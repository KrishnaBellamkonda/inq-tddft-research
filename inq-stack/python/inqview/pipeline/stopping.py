"""Phase: ``stopping`` — projectile stopping-power diagnostics.

Designed for the WP-vs-classical electron comparison (see
``docs/plans/the-objective-in-this-dapper-moon.md``). Produces:

* ``analysis/observables/dE_kinetic_vs_z.png``
  — change in kinetic energy of the *system* (``observables.csv``,
  ``energy_kinetic`` column) plotted against the projectile's z. For WP
  runs the projectile-z is taken from ``cod_z`` in observables.csv; for
  classical runs (no cod_z, no WP slot) it is taken from
  ``electron_track.csv``. The slope is a direct stopping-power proxy.

* ``analysis/observables/stopping_force_vs_z.png`` (classical run only)
  — F_z on the projectile (from ``electron_track.csv``) vs trajectory z.
  This is the integrand of the stopping power: integrating -F_z dz over
  the trajectory gives the energy lost by the projectile to the bath.

* ``analysis/observables/sigma_z_vs_time.png`` (WP run only)
  — measured WP packet width sigma_z(t) (from the projectile-only density
  if a WP density VTI series is available; otherwise SKIPPED with a
  message — we don't reload VTI here just for this overlay). Compared
  against the analytic free-particle expansion
  sigma(t) = sigma_0 sqrt(1 + (t/(2 m sigma_0^2))^2)
  with m=1 (electron mass), sigma_0 read from ``run_summary.txt``.
  Departure indicates many-body slowdown beyond free-particle quantum
  spreading.

Pure Python on existing CSVs (and optionally VTI directory listing for
the sigma_z overlay; reading VTI for sigma_z is left to a downstream
analysis script and not done here).
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import _common


HA_TO_EV = 27.21138625


def _read_run_summary_value(results_dir: Path, key: str) -> float | None:
    rs = results_dir / "run_summary.txt"
    if not rs.exists():
        return None
    m = re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", rs.read_text(),
                  flags=re.MULTILINE)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def _is_classical_run(results_dir: Path) -> bool:
    """Detect classical-projectile run by presence of electron_track.csv."""
    return (results_dir / "raw" / "observables" / "electron_track.csv").exists()


def _wp_sigma0_bohr(results_dir: Path) -> float | None:
    return _read_run_summary_value(results_dir, "wp_sigma_bohr")


def _free_particle_sigma(t_au: np.ndarray, sigma0: float, m_au: float = 1.0
                         ) -> np.ndarray:
    """sigma(t) = sigma0 * sqrt(1 + (t / (2 m sigma0^2))^2) — Gaussian
    free-particle quantum spreading in atomic units."""
    return sigma0 * np.sqrt(1.0 + (t_au / (2.0 * m_au * sigma0 ** 2)) ** 2)


def run(results_dir: Path, *, run_name: str, rebuild: bool, **opts) -> dict:
    obs_csv = results_dir / "raw" / "observables" / "observables.csv"
    if not obs_csv.exists():
        return {"skipped": f"missing: {obs_csv}"}

    obs = pd.read_csv(obs_csv)
    if obs.empty:
        return {"skipped": "observables.csv is empty"}

    out_dir = _common.ensure_dir(results_dir / "analysis" / "observables")
    classical = _is_classical_run(results_dir)
    out: dict = {"classical": classical, "outputs": []}

    # ----- 1. dE_kinetic vs trajectory z --------------------------------
    if classical:
        track_csv = results_dir / "raw" / "observables" / "electron_track.csv"
        track = pd.read_csv(track_csv)
        # Align on step (electron_track has every step; observables has every
        # WRITE_EVERY). Inner-join on step.
        #
        # fz is OPTIONAL (2026-07-30). The force columns are not part of the
        # minimum electron_track schema: a run may legitimately record only the
        # kinematics (step, t, x/y/z, vx/vy/vz) plus its own ke_ion_ha. The
        # kinetic-energy-vs-z panel below needs only `z`, so a missing fz must
        # degrade to skipping the force panel, NOT crash the whole phase — which
        # is what a hard `track[["step","z","fz"]]` selection did (KeyError
        # "['fz'] not in index", bulk_ks_stopping/classical).
        # Deliberately NOT defaulted to zeros: forces that were never recorded
        # are absent, not zero, and a flat F_z(z)=0 curve would read as physics.
        has_fz = "fz" in track.columns
        cols = ["step", "z"] + (["fz"] if has_fz else [])
        merged = pd.merge(obs, track[cols], on="step",
                          how="inner", suffixes=("", "_ion"))
        z_col   = "z"
        z_label = "projectile z (Bohr) [from electron_track.csv]"
    else:
        if "cod_z" not in obs.columns:
            return {"skipped": "WP run but no cod_z column in observables.csv"}
        merged   = obs.copy()
        z_col    = "cod_z"
        z_label  = "WP center-of-density z (Bohr) [from cod_z]"

    if "energy_kinetic" not in merged.columns:
        return {"skipped": "no energy_kinetic column in observables.csv"}

    e_kin = merged["energy_kinetic"].to_numpy()
    de_kin_ev = (e_kin - e_kin[0]) * HA_TO_EV

    out_png = out_dir / "dE_kinetic_vs_z.png"
    if _common.need_rebuild(out_png, rebuild):
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(merged[z_col], de_kin_ev, "k-", lw=1.5)
        ax.set_xlabel(z_label)
        ax.set_ylabel("dE_kinetic of system (eV)")
        ax.set_title(_common.title(run_name,
            "system kinetic-energy change vs. projectile z"))
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(out_png, dpi=150)
        plt.close(fig)
        out["outputs"].append(out_png)

    # ----- 2. Stopping force vs z (classical only, needs fz) ------------
    if classical and not has_fz:
        out.setdefault("notes", []).append(
            "stopping_force_vs_z skipped: electron_track.csv has no fz column "
            "(this run recorded projectile kinematics + ke_ion_ha but not forces)")
    if classical and has_fz:
        out_force = out_dir / "stopping_force_vs_z.png"
        if _common.need_rebuild(out_force, rebuild):
            fig, ax = plt.subplots(figsize=(8, 4.5))
            ax.plot(merged["z"], merged["fz"], "C0-", lw=1.0)
            ax.axhline(0.0, color="0.5", lw=0.5)
            ax.set_xlabel("projectile z (Bohr)")
            ax.set_ylabel("F_z on projectile (Ha/Bohr)")
            ax.set_title(_common.title(run_name,
                "stopping force on classical electron"))
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            fig.savefig(out_force, dpi=150)
            plt.close(fig)
            out["outputs"].append(out_force)

    # ----- 3. WP sigma_z(t) vs free-particle expansion (WP run only) ----
    if not classical:
        sigma0 = _wp_sigma0_bohr(results_dir)
        if sigma0 is None:
            out["skipped_sigma_z"] = "no wp_sigma_bohr in run_summary.txt"
        else:
            out_sigma = out_dir / "sigma_z_analytic_vs_time.png"
            if _common.need_rebuild(out_sigma, rebuild):
                t = obs["time_au"].to_numpy()
                sig_an = _free_particle_sigma(t, sigma0, m_au=1.0)
                fig, ax = plt.subplots(figsize=(8, 4.5))
                ax.plot(t, sig_an, "C2-", lw=1.5,
                        label="analytic free particle")
                ax.axhline(sigma0, color="0.6", lw=0.5,
                           label=f"sigma_0 = {sigma0:.2f} Bohr")
                ax.set_xlabel("time (a.u.)")
                ax.set_ylabel("sigma_z (Bohr)")
                ax.set_title(_common.title(run_name,
                    "WP analytic free-particle width vs. time"))
                ax.legend(loc="best", fontsize=9)
                ax.grid(True, alpha=0.3)
                fig.tight_layout()
                fig.savefig(out_sigma, dpi=150)
                plt.close(fig)
                out["outputs"].append(out_sigma)
                # Note: we DO NOT measure sigma_z from the WP density VTI
                # series here. That requires loading every density_wp_*.vti
                # and computing second moments on a real grid — better suited
                # to a dedicated downstream script. The analytic curve alone
                # is informative as a baseline (~4% growth at t=14 a.u. for
                # sigma_0=5).

    return out
