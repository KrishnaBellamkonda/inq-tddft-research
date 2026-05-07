"""Phase: ``wp_trajectory`` — WP centre-of-density tracking.

Reads ``results/raw/observables/observables.csv`` (which now carries
``cod_x_bohr``, ``cod_y_bohr``, ``cod_z_bohr``, ``density_l2`` columns
from the jellium run-template) and produces:

* ``analysis/observables/wp_position_vs_time.png``
   — one panel per axis showing <r>(t).
* ``analysis/observables/wp_velocity_vs_time.png``
   — finite-difference velocity v_z(t) (analytic check: at t=0,
   v_z should equal the launch wavenumber k_0 in atomic units).
* ``analysis/observables/density_fluctuation_l2.png``
   — sigma^2_n(t) = integral |dn|^2 dV from the dn observable.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import _common


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    csv_path = results_dir / "raw" / "observables" / "observables.csv"
    if not csv_path.exists():
        return {"skipped": f"missing: {csv_path}"}
    df = pd.read_csv(csv_path)

    out_dir = _common.ensure_dir(results_dir / "analysis" / "observables")
    artefacts: list[str] = []

    if {"cod_x_bohr", "cod_y_bohr", "cod_z_bohr"}.issubset(df.columns):
        fig, axes = plt.subplots(3, 1, figsize=(7, 7), sharex=True)
        for ax, col, label in zip(
            axes,
            ("cod_x_bohr", "cod_y_bohr", "cod_z_bohr"),
            ("<x>", "<y>", "<z>"),
        ):
            ax.plot(df["time_au"], df[col], lw=1.5)
            ax.set_ylabel(f"{label} (Bohr)")
            ax.grid(alpha=0.3)
        axes[-1].set_xlabel("time (a.u.)")
        fig.suptitle(f"{run_name}: WP centre of density vs time")
        fig.tight_layout()
        out = out_dir / "wp_position_vs_time.png"
        if rebuild or not out.exists():
            fig.savefig(out, dpi=130)
        plt.close(fig)
        artefacts.append(str(out))

        t = df["time_au"].to_numpy()
        z = df["cod_z_bohr"].to_numpy()
        vz = np.gradient(z, t)
        fig, ax = plt.subplots(figsize=(7, 3))
        ax.plot(t, vz, lw=1.2, color="C2")
        ax.set_xlabel("time (a.u.)")
        ax.set_ylabel("v_z = d<z>/dt (Bohr / a.u.)")
        ax.set_title(f"{run_name}: WP velocity along z")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        out = out_dir / "wp_velocity_vs_time.png"
        if rebuild or not out.exists():
            fig.savefig(out, dpi=130)
        plt.close(fig)
        artefacts.append(str(out))

    if "density_l2" in df.columns:
        fig, ax = plt.subplots(figsize=(7, 3))
        ax.plot(df["time_au"], df["density_l2"], lw=1.2, color="C3")
        ax.set_xlabel("time (a.u.)")
        ax.set_ylabel(r"$\sigma^2_n(t) = \int |\Delta n|^2 \, dV$")
        ax.set_title(f"{run_name}: integrated local density fluctuation")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        out = out_dir / "density_fluctuation_l2.png"
        if rebuild or not out.exists():
            fig.savefig(out, dpi=130)
        plt.close(fig)
        artefacts.append(str(out))

    return {"artefacts": artefacts}
