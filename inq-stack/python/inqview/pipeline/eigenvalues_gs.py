"""Phase: ``eigenvalues_gs`` — ground-state band structure visualisations.

Reads ``results/raw/observables/eigenvalues/eigenvalues.csv``. The full
schema written by ``inqkit::observables::dump_eigenvalues`` is
``kpoint_index, kx, ky, kz, weight, state_index, eigenvalue_ha,
eigenvalue_ev, occupation``. The simpler schema used by the older
``jellium::eigenvalues::dump`` is just ``state_index, eigenvalue_ha,
eigenvalue_ev`` plus a sister ``occupations.csv`` (``state_index,
occupation``); the loader synthesises defaults for the missing columns
(Γ-only single-k-point, weight=1, kx=ky=kz=0) so the same plots can
still be produced.

Outputs:

* ``analysis/observables/eigenvalues/eigenvalues_levels.png``
   — A3: horizontal level diagram, one column per k-point, eV scale,
   colour-coded by occupation; Fermi level dashed.
* ``analysis/observables/eigenvalues/dos.png`` (+ ``dos.csv``)
   — A4: Gaussian-broadened weighted density of states, eV axis.
* ``analysis/observables/eigenvalues/eigenvalue_table.txt``
   — A5: plain-text dump of (k, state, ε_eV, occ) for quick reference.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import _common

HA_TO_EV = 27.21138625


def _band_fraction(df: pd.DataFrame) -> pd.Series:
    """INQ stores ``occupations[ik][i] = f_i * weight_k * spin_factor``, so a
    fully-occupied band at a 1/8-weighted k-point has the value 0.25. Divide
    out the (weight × spin) factor so the column is back in [0, 1] for
    threshold filtering.

    Spin factor is 2 for the spin-restricted runs we use. For collinear runs
    INQ does not multiply by 2 (each spin channel contributes separately) — we
    auto-detect by checking the maximum of ``occupation/weight``.
    """
    raw = df["occupation"] / df["weight"]
    spin_factor = 2.0 if raw.max() > 1.5 else 1.0
    return (raw / spin_factor).clip(lower=0.0, upper=1.0)


def _fermi_level_ev(df: pd.DataFrame, occ_threshold: float = 0.5) -> float:
    bf = _band_fraction(df)
    occ = df[bf > occ_threshold]
    unocc = df[bf <= occ_threshold]
    if occ.empty or unocc.empty:
        return float(df["eigenvalue_ev"].median())
    e_homo = float(occ["eigenvalue_ev"].max())
    e_lumo = float(unocc["eigenvalue_ev"].min())
    return 0.5 * (e_homo + e_lumo)


def _plot_levels(df: pd.DataFrame, out_path: Path, *, run_name: str,
                 e_window_ev: tuple[float, float] | None = None) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    kpts = sorted(df["kpoint_index"].unique())
    e_fermi = _fermi_level_ev(df)

    bf = _band_fraction(df)
    for ik in kpts:
        sub_idx = df.index[df["kpoint_index"] == ik]
        x = ik
        for i in sub_idx:
            row = df.loc[i]
            colour = plt.cm.RdBu_r(0.85) if bf.loc[i] > 0.5 else plt.cm.RdBu_r(0.15)
            ax.hlines(row["eigenvalue_ev"] - e_fermi,
                      xmin=x - 0.35, xmax=x + 0.35,
                      colors=colour, lw=1.2, alpha=0.85)

    ax.axhline(0.0, color="black", lw=0.8, ls="--", alpha=0.6,
               label=f"E_F = {e_fermi:.3f} eV")
    if e_window_ev is not None:
        ax.set_ylim(*e_window_ev)
    ax.set_xticks(kpts)
    ax.set_xticklabels([f"k{ik}" for ik in kpts])
    ax.set_xlabel("k-point index")
    ax.set_ylabel(r"$\varepsilon - E_F$ (eV)")
    ax.set_title(f"{run_name}: KS eigenvalues (red=occupied, blue=empty)")
    ax.grid(alpha=0.3, axis="y")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def _plot_dos(df: pd.DataFrame, out_png: Path, out_csv: Path,
              *, run_name: str, sigma_ev: float = 0.1) -> None:
    eps = df["eigenvalue_ev"].to_numpy()
    w   = df["weight"].to_numpy()
    e_min, e_max = float(eps.min()) - 1.0, float(eps.max()) + 1.0
    grid = np.linspace(e_min, e_max, 2000)
    dos = np.zeros_like(grid)
    norm = 1.0 / (sigma_ev * np.sqrt(2.0 * np.pi))
    for e_i, w_i in zip(eps, w):
        dos += w_i * norm * np.exp(-0.5 * ((grid - e_i) / sigma_ev) ** 2)

    pd.DataFrame({"energy_ev": grid, "dos": dos}).to_csv(out_csv, index=False)

    e_fermi = _fermi_level_ev(df)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(grid - e_fermi, dos, color="black", lw=1.4)
    ax.fill_between(grid - e_fermi, 0.0, dos,
                    where=(grid <= e_fermi),
                    color=plt.cm.RdBu_r(0.85), alpha=0.35,
                    label="occupied")
    ax.axvline(0.0, color="black", lw=0.8, ls="--",
               label=f"E_F = {e_fermi:.3f} eV")
    ax.set_xlabel(r"$\varepsilon - E_F$ (eV)")
    ax.set_ylabel(rf"DOS (Gaussian σ={sigma_ev} eV)")
    ax.set_title(f"{run_name}: weighted density of states")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def _dump_table(df: pd.DataFrame, out_path: Path) -> None:
    cols = ["kpoint_index", "kx", "ky", "kz", "weight",
            "state_index", "eigenvalue_ev", "occupation"]
    sub = df[cols].copy()
    sub["eigenvalue_ev"] = sub["eigenvalue_ev"].map(lambda x: f"{x: .6f}")
    sub["occupation"]   = sub["occupation"].map(lambda x: f"{x: .4f}")
    sub["weight"]       = sub["weight"].map(lambda x: f"{x: .4f}")
    for c in ("kx", "ky", "kz"):
        sub[c] = sub[c].map(lambda x: f"{x: .4f}")
    with open(out_path, "w") as f:
        f.write(f"# Eigenvalue table\n")
        f.write(f"# columns: {' '.join(cols)}\n")
        f.write(sub.to_string(index=False))
        f.write("\n")


def _load_eigenvalues_with_defaults(csv_path: Path) -> pd.DataFrame:
    """Load eigenvalues.csv; if columns are missing (older Γ-only dumps),
    backfill with sensible defaults so the rest of the module can run.
    Joins a sister ``occupations.csv`` when ``occupation`` is missing."""
    df = pd.read_csv(csv_path)

    if "kpoint_index" not in df.columns:
        df["kpoint_index"] = 0
    for c in ("kx", "ky", "kz"):
        if c not in df.columns:
            df[c] = 0.0
    if "weight" not in df.columns:
        df["weight"] = 1.0

    if "occupation" not in df.columns:
        occ_path = csv_path.parent / "occupations.csv"
        if occ_path.exists():
            occ = pd.read_csv(occ_path)
            df = df.merge(occ, on="state_index", how="left")
        if "occupation" not in df.columns or df["occupation"].isna().any():
            df["occupation"] = df.get("occupation", 0.0)
            df["occupation"] = df["occupation"].fillna(0.0)
    return df


def run(results_dir: Path, *, run_name: str, rebuild: bool, **opts) -> dict:
    csv_path = results_dir / "raw" / "observables" / "eigenvalues" / "eigenvalues.csv"
    if not csv_path.exists():
        return {"skipped": f"missing: {csv_path}"}

    df = _load_eigenvalues_with_defaults(csv_path)
    out_dir = _common.ensure_dir(
        results_dir / "analysis" / "observables" / "eigenvalues")

    artefacts: list[str] = []

    levels_png = out_dir / "eigenvalues_levels.png"
    if rebuild or not levels_png.exists():
        _plot_levels(df, levels_png, run_name=run_name,
                     e_window_ev=opts.get("e_window_ev"))
    artefacts.append(str(levels_png))

    dos_png = out_dir / "dos.png"
    dos_csv = out_dir / "dos.csv"
    if rebuild or not dos_png.exists() or not dos_csv.exists():
        _plot_dos(df, dos_png, dos_csv, run_name=run_name,
                  sigma_ev=float(opts.get("dos_sigma_ev", 0.1)))
    artefacts.extend([str(dos_png), str(dos_csv)])

    tbl = out_dir / "eigenvalue_table.txt"
    if rebuild or not tbl.exists():
        _dump_table(df, tbl)
    artefacts.append(str(tbl))

    return {"artefacts": artefacts}
