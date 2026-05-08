"""Phase: ``gs_projected_occupations`` — KS-orbital occupations
projected onto the *initial* GS basis.

Given the time-evolved KS orbitals psi_j(t), the overlap matrix
``O_{ij}(t) = |<psi_i^GS | psi_j(t)>|^2`` lets us redistribute the
initial occupations f_j(0) onto the frozen GS basis indexed by i:

    n_i^GS(t)  =  sum_j  f_j(0)  *  |<psi_i^GS | psi_j(t)>|^2

At t=0 this equals f_i(0) exactly (overlap is the identity). At later t
it spreads the initial occupations over the GS basis according to where
each evolved orbital has acquired a component. The contribution of
initially-occupied evolved orbitals to *unoccupied* GS orbitals is
exactly the single-particle-excitation occupation transfer — the
cleanest TDDFT definition of "excitation into orbital i".

Inputs (read from results/raw/observables/):

* overlap_full/  — full O matrix (n_ref x n_evolved) at one or more
  snapshots, with index.csv listing (step, time_au, file). This is the
  expensive but exact input.
* overlap_proxies/ — proxy O matrix (n_ref x n_proxies) at many
  snapshots, with shells.csv giving the shell structure (shell_id,
  proxy_indices, degeneracy, evolved_indices_in_shell). Cheaper input
  with shell-averaging weights to estimate the full sum.
* eigenvalues/occupations.csv — initial occupations f_j(0), generated
  by the inqkit::observables::dump_eigenvalues() in the GS save script.

Outputs in analysis/observables/gs_projected_occupations/:

* n_i_gs_vs_time.csv — time series of n_i^GS(t) per i, long format:
  step,time_au,gs_state_index,n_gs.
* excitation_total_vs_time.png — total excitation occupation
  ``sum_{i above E_F} n_i^GS(t)`` over time; the projectile-driven
  excitation amplitude.
* gs_projection_t_<NN>.png — bar chart at each available snapshot;
  blue = initially-occupied i, red = initially-empty i (excitation
  occupation).
* heatmap_overlap_t_<NN>.png — full |O_{ij}(t)|^2 as a heatmap (when
  full matrix is available; skipped for proxy-only input).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import _common


def _load_overlap_csv(path: Path) -> tuple[np.ndarray, dict]:
    """Read one overlap_NNNNNN.csv into an array and metadata dict.

    The first line is a comment with metadata of the form
        # step=N time_au=T n_ref=R n_evolved=E [mode=...]
    Subsequent lines are R rows of E comma-separated |O|^2 values.
    """
    with open(path) as fh:
        header = fh.readline().lstrip("# ").strip()
    meta = {}
    for tok in header.split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            try:
                meta[k] = int(v)
            except ValueError:
                try:
                    meta[k] = float(v)
                except ValueError:
                    meta[k] = v
    data = pd.read_csv(path, comment="#", header=None).values
    return data, meta


def _load_index(idx_path: Path) -> pd.DataFrame:
    return pd.read_csv(idx_path)


def _load_initial_occupations(results_dir: Path) -> Optional[np.ndarray]:
    """Read raw/observables/eigenvalues/occupations.csv into a 1D array
    indexed by state_index. Returns None if missing."""
    p = results_dir / "raw" / "observables" / "eigenvalues" / "occupations.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    df = df.sort_values("state_index").reset_index(drop=True)
    return df["occupation"].to_numpy(dtype=float)


def _load_shells(results_dir: Path) -> Optional[pd.DataFrame]:
    """Read overlap_proxies/shells.csv if present.

    Expected columns: shell_id, degeneracy, n_proxies, proxy_indices
    (comma-separated within a quoted field), evolved_indices
    (comma-separated)."""
    p = results_dir / "raw" / "observables" / "overlap_proxies" / "shells.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def _project_full(O_sq: np.ndarray, f_init: np.ndarray) -> np.ndarray:
    """Compute n_i^GS = sum_j f_j(0) * O_{ij} from a full matrix.

    O_sq has shape (n_ref, n_evolved). f_init has length >= n_evolved.
    Returns array of length n_ref.
    """
    n_evolved = O_sq.shape[1]
    return O_sq @ f_init[:n_evolved]


def _project_proxy(O_sq: np.ndarray,
                   f_init: np.ndarray,
                   shells: pd.DataFrame) -> np.ndarray:
    """Compute n_i^GS using shell-averaged proxy overlaps.

    For each shell s with degeneracy g_s and proxy indices P_s, the
    contribution from shell s is

        contrib_s(i) = (g_s / |P_s|) * sum_{j in P_s} f_j(0) * O_{i,j}^proxy

    where the columns of ``O_sq`` correspond to the proxies in some
    order given by ``shells['proxy_indices']``. The total is summed
    over shells.
    """
    n_ref = O_sq.shape[0]
    n_gs = np.zeros(n_ref)
    # Build map from proxy column index to (shell_id, proxy_state_index)
    col = 0
    for _, row in shells.iterrows():
        degeneracy = int(row["degeneracy"])
        proxy_idxs = [int(x) for x in str(row["proxy_indices"]).split(",")]
        n_proxies = len(proxy_idxs)
        for j_state in proxy_idxs:
            # f_init is indexed by absolute orbital index
            weight = (degeneracy / n_proxies) * f_init[j_state]
            n_gs += weight * O_sq[:, col]
            col += 1
    if col != O_sq.shape[1]:
        raise ValueError(
            f"shell book-keeping mismatch: "
            f"{col} proxy cols accounted vs {O_sq.shape[1]} in matrix"
        )
    return n_gs


def _detect_homo_index(f_init: np.ndarray, occ_threshold: float = 0.5
                       ) -> int:
    """Highest state index with occupation >= threshold.

    With INQ's spin-paired storage f_max = 2; threshold 0.5 covers
    smeared band-edge states."""
    occupied = np.where(f_init >= occ_threshold)[0]
    return int(occupied.max()) if len(occupied) else -1


def run(results_dir: Path, *, run_name: str, rebuild: bool, **opts) -> dict:
    obs_root = results_dir / "raw" / "observables"
    full_dir = obs_root / "overlap_full"
    proxy_dir = obs_root / "overlap_proxies"

    f_init = _load_initial_occupations(results_dir)
    if f_init is None:
        return {"skipped": "no eigenvalues/occupations.csv at GS"}

    homo = _detect_homo_index(f_init)
    n_states = len(f_init)

    # Pick data source: prefer proxy if present (typically denser in
    # time), fall back to full.
    snapshots = []
    if proxy_dir.exists() and (proxy_dir / "index.csv").exists():
        idx = _load_index(proxy_dir / "index.csv")
        shells = _load_shells(results_dir)
        if shells is None:
            return {"skipped": "overlap_proxies/ present but shells.csv missing"}
        for _, row in idx.iterrows():
            csv_path = proxy_dir / row["file"]
            if not csv_path.exists():
                continue
            O_sq, meta = _load_overlap_csv(csv_path)
            n_gs_t = _project_proxy(O_sq, f_init, shells)
            snapshots.append({
                "step": int(row["step"]),
                "time_au": float(row["time_au"]),
                "n_gs": n_gs_t,
                "O_sq": None,         # proxy is not square; skip heatmap
                "source": "proxy",
            })
    elif full_dir.exists() and (full_dir / "index.csv").exists():
        idx = _load_index(full_dir / "index.csv")
        for _, row in idx.iterrows():
            csv_path = full_dir / row["file"]
            if not csv_path.exists():
                continue
            O_sq, meta = _load_overlap_csv(csv_path)
            n_gs_t = _project_full(O_sq, f_init)
            snapshots.append({
                "step": int(row["step"]),
                "time_au": float(row["time_au"]),
                "n_gs": n_gs_t,
                "O_sq": O_sq,
                "source": "full",
            })
    else:
        return {"skipped": "no overlap_full/ or overlap_proxies/ data"}

    if not snapshots:
        return {"skipped": "overlap dirs exist but no snapshots loaded"}

    out_dir = _common.ensure_dir(
        results_dir / "analysis" / "observables" / "gs_projected_occupations"
    )

    # ----- CSV: long-format n_i^GS(t) ----------------------------------
    rows = []
    for s in snapshots:
        for i, n in enumerate(s["n_gs"]):
            rows.append({
                "step": s["step"],
                "time_au": s["time_au"],
                "gs_state_index": i,
                "n_gs": float(n),
                "source": s["source"],
            })
    df = pd.DataFrame(rows)
    csv_path = out_dir / "n_i_gs_vs_time.csv"
    df.to_csv(csv_path, index=False)

    # ----- Total excitation occupation vs time -------------------------
    excitation_rows = []
    for s in snapshots:
        n_gs = s["n_gs"]
        # Total occupation projected onto initially-empty (i > homo) and
        # initially-filled (i <= homo) GS orbitals.
        unocc = float(n_gs[homo + 1 :].sum()) if homo + 1 < len(n_gs) else 0.0
        occ = float(n_gs[: homo + 1].sum())
        # "Excitation occupation" = mass that left occupied orbitals.
        # This is sum_{j occupied} f_j(0) - n_j^GS(t), summed:
        sum_initial_occ = float(f_init[: homo + 1].sum())
        loss_from_occ = sum_initial_occ - occ
        excitation_rows.append({
            "step": s["step"],
            "time_au": s["time_au"],
            "n_gs_sum_unoccupied": unocc,
            "n_gs_sum_occupied":   occ,
            "loss_from_occupied":  loss_from_occ,
            "total_initial_occ":   sum_initial_occ,
        })
    exc_df = pd.DataFrame(excitation_rows)
    exc_csv = out_dir / "excitation_total_vs_time.csv"
    exc_df.to_csv(exc_csv, index=False)

    # ----- Plot total excitation -------------------------------------
    out_exc = out_dir / "excitation_total_vs_time.png"
    if _common.need_rebuild(out_exc, rebuild):
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(exc_df["time_au"], exc_df["loss_from_occupied"],
                "o-", color="C3", label="loss from occupied (= excitation)")
        ax.plot(exc_df["time_au"], exc_df["n_gs_sum_unoccupied"],
                "s-", color="C0", label="sum n_i^GS over unoccupied i")
        ax.set_xlabel("time (a.u.)")
        ax.set_ylabel("excitation occupation")
        ax.set_title(_common.title(run_name,
            "GS-projected excitation: occupation transferred out of "
            f"occupied subspace (HOMO=index {homo})"))
        ax.legend(loc="best", fontsize=9)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(out_exc, dpi=150)
        plt.close(fig)

    # ----- Per-snapshot bar charts -----------------------------------
    bar_pngs = []
    for s in snapshots:
        out_bar = out_dir / f"gs_projection_t{s['step']:06d}.png"
        if _common.need_rebuild(out_bar, rebuild):
            fig, ax = plt.subplots(figsize=(10, 4.5))
            i_arr = np.arange(len(s["n_gs"]))
            colors = ["C0" if i <= homo else "C3" for i in i_arr]
            ax.bar(i_arr, s["n_gs"], color=colors, width=1.0)
            ax.axvline(homo + 0.5, color="k", linestyle="--", linewidth=0.7,
                       label=f"HOMO at i={homo}")
            ax.set_xlabel("GS orbital index i")
            ax.set_ylabel("n_i^GS(t)")
            ax.set_title(_common.title(run_name,
                f"GS-projected occupation at step {s['step']} "
                f"(t={s['time_au']:.3f} a.u., source={s['source']})"))
            ax.legend(loc="upper right", fontsize=9)
            ax.grid(True, axis="y", alpha=0.3)
            fig.tight_layout()
            fig.savefig(out_bar, dpi=150)
            plt.close(fig)
            bar_pngs.append(out_bar)

    # ----- Heatmaps for full-matrix snapshots ------------------------
    heatmap_pngs = []
    for s in snapshots:
        if s["O_sq"] is None:
            continue
        out_hm = out_dir / f"heatmap_overlap_t{s['step']:06d}.png"
        if _common.need_rebuild(out_hm, rebuild):
            fig, ax = plt.subplots(figsize=(7, 6))
            # log-color so the off-diagonal is visible
            data = np.log10(np.maximum(s["O_sq"], 1e-8))
            im = ax.imshow(data, aspect="auto", origin="lower",
                           cmap="viridis", vmin=-6, vmax=0)
            ax.set_xlabel("evolved orbital index j")
            ax.set_ylabel("GS reference orbital index i")
            ax.set_title(_common.title(run_name,
                f"|O_ij(t)|^2 (log10) at step {s['step']} "
                f"(t={s['time_au']:.3f} a.u.)"))
            fig.colorbar(im, ax=ax, label="log10 |O_ij|^2")
            fig.tight_layout()
            fig.savefig(out_hm, dpi=150)
            plt.close(fig)
            heatmap_pngs.append(out_hm)

    return {
        "n_snapshots": len(snapshots),
        "source":      snapshots[0]["source"] if snapshots else None,
        "homo_index":  homo,
        "n_states":    n_states,
        "csv":         csv_path,
        "excitation_csv": exc_csv,
        "excitation_png": out_exc,
        "n_bar_pngs":  len(bar_pngs),
        "n_heatmap_pngs": len(heatmap_pngs),
    }
