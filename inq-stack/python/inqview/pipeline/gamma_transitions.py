"""Phase: ``gamma_transitions`` — D5 single-k-point inter-band transition
histogram.

Mirrors the iter-5 published analysis
(``QuantumKickExtension/iter 5/image6.png``): for one k-point (Γ or its
closest discrete approximant on the chosen MP grid), enumerate all
(occupied → unoccupied) Kohn-Sham single-particle pairs at the ground
state, histogram their transition energies ε_n' − ε_n, shade the gap
region, and overlay the dominant FFT peak energies of the excess-energy
time-series so the kinematics of the absorbed energy are visible at a
glance.

Inputs:
* ``raw/observables/eigenvalues/eigenvalues.csv``
* (optional) ``raw/observables/fft_excess_energy.csv`` — peaks overlaid
  if present.

Outputs:
* ``raw/observables/eigenvalues/gamma_transitions.csv``
   columns: n_initial, n_final, eigenvalue_initial_ev,
            eigenvalue_final_ev, transition_energy_ev
* ``analysis/observables/eigenvalues/gamma_gamma_transitions.png``

Reference: ``iter 5/image6.png`` from the published QBall analysis.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import _common


def _select_gamma_kpoint(df: pd.DataFrame) -> tuple[int, np.ndarray, str]:
    """Return (kpoint_index, k-coord, label) for the k-point closest to Γ.

    For a 2x2x2 shifted MP grid no point sits exactly at Γ — we pick the
    one with the smallest |k| and label it explicitly. For an unshifted
    grid Γ is in the set and the label reads "Γ" exactly.
    """
    by_k = (df.drop_duplicates("kpoint_index")
              .set_index("kpoint_index")[["kx", "ky", "kz"]])
    norms = by_k.apply(lambda r: float(np.linalg.norm(r.values)), axis=1)
    ik = int(norms.idxmin())
    kvec = by_k.loc[ik].to_numpy(dtype=float)
    label = ("Γ" if float(norms.loc[ik]) < 1e-6
             else f"k closest to Γ "
                  f"({kvec[0]:+.3f}, {kvec[1]:+.3f}, {kvec[2]:+.3f}) "
                  "(2π/a units)")
    return ik, kvec, label


def _band_fraction(df: pd.DataFrame) -> pd.Series:
    """INQ stores ``occupations[ik][i] = f * weight * spin_factor``. Divide
    out (weight × spin) so the result is the per-band Fermi-Dirac occupation
    in [0, 1]. See :mod:`eigenvalues_gs._band_fraction` for the full note.
    """
    raw = df["occupation"] / df["weight"]
    spin_factor = 2.0 if raw.max() > 1.5 else 1.0
    return (raw / spin_factor).clip(lower=0.0, upper=1.0)


def _build_transitions(df_k: pd.DataFrame,
                       occ_threshold: float = 0.99,
                       unocc_threshold: float = 0.01) -> pd.DataFrame:
    bf = _band_fraction(df_k)
    occ = df_k[bf >= occ_threshold]
    unocc = df_k[bf <= unocc_threshold]
    rows = []
    for _, n in occ.iterrows():
        for _, np_ in unocc.iterrows():
            rows.append({
                "n_initial": int(n["state_index"]),
                "n_final":   int(np_["state_index"]),
                "eigenvalue_initial_ev": float(n["eigenvalue_ev"]),
                "eigenvalue_final_ev":   float(np_["eigenvalue_ev"]),
                "transition_energy_ev":  float(np_["eigenvalue_ev"] -
                                                n["eigenvalue_ev"]),
            })
    return pd.DataFrame(rows)


def _gap_ev(df_k: pd.DataFrame,
            occ_threshold: float = 0.99,
            unocc_threshold: float = 0.01) -> float | None:
    bf = _band_fraction(df_k)
    occ = df_k[bf >= occ_threshold]
    unocc = df_k[bf <= unocc_threshold]
    if occ.empty or unocc.empty:
        return None
    return float(unocc["eigenvalue_ev"].min() - occ["eigenvalue_ev"].max())


def _load_fft_peaks(results_dir: Path,
                    *, max_peaks: int = 3,
                    e_max_ev: float = 12.0) -> list[tuple[float, float]]:
    """Return up to ``max_peaks`` strongest peaks (ω_eV, amplitude) from
    ``raw/observables/fft_excess_energy.csv`` if present, restricted to the
    [0, e_max_ev] window. Empty if file missing.
    """
    fft_csv = results_dir / "raw" / "observables" / "fft_excess_energy.csv"
    if not fft_csv.exists():
        return []
    df = pd.read_csv(fft_csv)
    # Column-name detection: qball / inqview / analyse_inq.py have all
    # used slightly different names; fall back to positional if no match.
    omega_candidates = ["omega_eV", "omega_ev", "energy_ev"]
    amp_candidates   = ["sp_norm", "power_norm_0_20ev",
                        "power_norm", "amplitude", "power_unnormalised"]
    omega_col = next((c for c in omega_candidates if c in df.columns),
                     df.columns[0])
    amp_col   = next((c for c in amp_candidates if c in df.columns),
                     df.columns[1])
    sub = df[(df[omega_col] >= 0.0) & (df[omega_col] <= e_max_ev)].copy()
    if sub.empty:
        return []
    omega = sub[omega_col].to_numpy()
    amp   = sub[amp_col].to_numpy()
    # Cheap local-maxima detection
    peaks: list[tuple[float, float]] = []
    for i in range(1, len(amp) - 1):
        if amp[i] > amp[i-1] and amp[i] > amp[i+1] and amp[i] > 0.05:
            peaks.append((float(omega[i]), float(amp[i])))
    peaks.sort(key=lambda p: -p[1])
    return peaks[:max_peaks]


def _plot_histogram(transitions: pd.DataFrame,
                    out_png: Path, *,
                    run_name: str,
                    k_label: str,
                    gap_ev: float | None,
                    fft_peaks: list[tuple[float, float]],
                    bin_width_ev: float = 0.1,
                    e_max_ev: float = 8.0) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    if not transitions.empty:
        bins = np.arange(0.0, e_max_ev + bin_width_ev, bin_width_ev)
        ax.hist(transitions["transition_energy_ev"], bins=bins,
                color="firebrick", alpha=0.85, edgecolor="firebrick")

    if gap_ev is not None and gap_ev > 0.0:
        ax.axvspan(0.0, gap_ev, color="grey", alpha=0.18)
        ax.text(gap_ev * 0.5, ax.get_ylim()[1] * 0.92,
                f"gap = {gap_ev:.2f} eV",
                ha="center", va="top", color="black", fontsize=9)

    palette = ["royalblue", "seagreen", "darkorange"]
    for (omega, amp), colour in zip(fft_peaks, palette):
        ax.axvline(omega, ls="--", color=colour, lw=1.5,
                   label=f"FFT peak {omega:.2f} eV")

    ax.set_xlim(0.0, e_max_ev)
    ax.set_xlabel(r"Transition energy $\varepsilon_{n'} - \varepsilon_n$ (eV)")
    ax.set_ylabel("Number of single-particle transitions")
    ax.set_title(f"{run_name}: occupied → unoccupied transitions at {k_label}")
    if fft_peaks:
        ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def run(results_dir: Path, *, run_name: str, rebuild: bool, **opts) -> dict:
    eig_csv = results_dir / "raw" / "observables" / "eigenvalues" / "eigenvalues.csv"
    if not eig_csv.exists():
        return {"skipped": f"missing: {eig_csv}"}

    df = pd.read_csv(eig_csv)

    # Older Γ-only dumps (jellium::eigenvalues::dump) lack
    # ``kpoint_index, kx, ky, kz, weight, occupation`` — synthesise them so
    # the rest of this phase can run unmodified. The legacy split puts
    # occupation in a sister occupations.csv.
    if "kpoint_index" not in df.columns:
        df["kpoint_index"] = 0
    for c in ("kx", "ky", "kz"):
        if c not in df.columns:
            df[c] = 0.0
    if "weight" not in df.columns:
        df["weight"] = 1.0
    if "occupation" not in df.columns:
        occ_path = eig_csv.parent / "occupations.csv"
        if occ_path.exists():
            occ = pd.read_csv(occ_path)
            df = df.merge(occ, on="state_index", how="left")
        if "occupation" not in df.columns:
            df["occupation"] = 0.0
        df["occupation"] = df["occupation"].fillna(0.0)

    ik_gamma, kvec, k_label = _select_gamma_kpoint(df)
    df_k = df[df["kpoint_index"] == ik_gamma].copy()

    occ_thr   = float(opts.get("occ_threshold", 0.99))
    unocc_thr = float(opts.get("unocc_threshold", 0.01))
    transitions = _build_transitions(df_k, occ_thr, unocc_thr)
    gap = _gap_ev(df_k, occ_thr, unocc_thr)

    out_csv = (results_dir / "raw" / "observables" / "eigenvalues"
               / "gamma_transitions.csv")
    transitions.to_csv(out_csv, index=False)

    out_png = (_common.ensure_dir(
        results_dir / "analysis" / "observables" / "eigenvalues")
        / "gamma_gamma_transitions.png")

    fft_peaks = _load_fft_peaks(
        results_dir,
        max_peaks=int(opts.get("fft_max_peaks", 3)),
        e_max_ev=float(opts.get("fft_e_max_ev", 12.0)),
    )

    e_max = float(opts.get("hist_e_max_ev", 8.0))
    bin_w = float(opts.get("hist_bin_width_ev", 0.1))
    if rebuild or not out_png.exists():
        _plot_histogram(transitions, out_png, run_name=run_name,
                        k_label=k_label, gap_ev=gap, fft_peaks=fft_peaks,
                        bin_width_ev=bin_w, e_max_ev=e_max)

    return {
        "artefacts": [str(out_csv), str(out_png)],
        "n_transitions": int(len(transitions)),
        "gap_ev": gap,
        "k_label": k_label,
        "kpoint_index": ik_gamma,
        "kvec_2pia": kvec.tolist(),
    }
