"""Wave-packet integrity metrics (IV-M05).

Quantifies how much of the WP is preserved during propagation, via two
complementary, translation-aware measures:

- **momentum_kl** — KL(P_t‖P_0) of the WP momentum distribution. Translation-
  INVARIANT (a rigidly moving WP keeps its momentum spectrum), so it rises only
  when the bath genuinely scatters/redistributes momentum.
- **real-space spread σ_r(t) + ipr** — spatial localisation. A free WP spreads
  as σ_r(t)=σ₀√(1+(t/τ)²); the inverse participation ratio falls as it
  delocalises.

Pure numpy (no matplotlib / VTK). KL is reimplemented here rather than imported
from ``postprocess.kl_divergence`` (which pulls matplotlib) to keep
``inqview.analysis`` deps-clean.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import numpy as np

_EPS = 1e-300


def _normalise(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    s = float(p.sum())
    return p / s if s > 0 else np.zeros_like(p)


def momentum_kl(p_t: np.ndarray, p_0: np.ndarray) -> float:
    """KL(P_t‖P_0) in nats; bins with P_0≤0 are dropped. ≥0; 0 iff P_t==P_0."""
    pt, q0 = _normalise(p_t), _normalise(p_0)
    mask = q0 > 0
    if not mask.any():
        return float("nan")
    pp = np.where(mask, pt, 0.0) + _EPS
    qq = np.where(mask, q0, 1.0)
    return float(np.where(mask, pp * np.log(pp / qq), 0.0).sum())


def kl_series(P_series: np.ndarray, reference: str = "initial") -> np.ndarray:
    """KL(P_t‖P_ref) time series for a stack of distributions (n_t, n_bins).

    ``reference='initial'`` → drift from launch KL(P_t‖P_0) (rises monotonically
    as the WP scatters). ``reference='previous'`` → frame-to-frame drift RATE
    KL(P_t‖P_{t−1}) (instantaneous change; ~0 for a steady WP). Both start at 0.
    """
    if reference not in ("initial", "previous"):
        raise ValueError("reference must be 'initial' or 'previous'")
    P = np.asarray(P_series, dtype=float)
    nt = P.shape[0]
    out = np.zeros(nt)
    for i in range(nt):
        ref = P[0] if reference == "initial" else P[max(i - 1, 0)]
        out[i] = momentum_kl(P[i], ref)
    return out


def ipr(rho: np.ndarray, dV: float = 1.0) -> float:
    """Inverse participation ratio (∫ρ² dV)/(∫ρ dV)² — localisation (larger =
    more localised). Units 1/volume; for normalised ρ it is ∫ρ²dV."""
    rho = np.asarray(rho, dtype=float)
    num = float((rho * rho).sum()) * dV
    den = (float(rho.sum()) * dV) ** 2
    return num / den if den > 0 else 0.0


def real_space_variance(rho: np.ndarray, coord: np.ndarray) -> float:
    """Density-weighted variance ⟨(x−⟨x⟩)²⟩ along one axis (Bohr²)."""
    rho = np.clip(np.asarray(rho, dtype=float), 0.0, None)
    s = float(rho.sum())
    if s <= 0:
        return 0.0
    x = np.asarray(coord, dtype=float)
    mean = float((rho * x).sum() / s)
    return float((rho * (x - mean) ** 2).sum() / s)


@dataclass(frozen=True)
class WPIntegrity:
    """WP-integrity time series (assembled per-frame from a run's outputs)."""

    time_au: np.ndarray
    kl_mom: np.ndarray          # momentum KL(P_t‖P_0) (nats)
    sigma_r: np.ndarray         # real-space spread (Bohr)
    ipr: np.ndarray             # localisation (1/Bohr³); NaN if WP density absent


# --- from-run assembly ------------------------------------------------------
# Default locations within a run directory (results/raw/observables/...).
_MOM_REL = "results/raw/observables/momentum_distribution.csv"
_RS_REL = "results/raw/observables/wp_real_space_stats.csv"


def assemble_from_run(
    run_dir: Union[str, Path],
    *,
    momentum_csv: Optional[Union[str, Path]] = None,
    real_space_csv: Optional[Union[str, Path]] = None,
    reference: str = "initial",
) -> WPIntegrity:
    """Assemble a :class:`WPIntegrity` time series from a run's CSV outputs.

    Reads ``momentum_distribution.csv`` (long format
    ``step,time_au,k_bohr_inv,n_total,n_wp`` — the WP momentum spectrum is the
    ``n_wp`` column per step) and ``wp_real_space_stats.csv`` (per-step
    ``sigma_x2/y2/z2``). Computes:

    - ``kl_mom`` = ``kl_series`` of the per-step WP distributions vs ``reference``
      (``'initial'`` → drift from launch; ``'previous'`` → frame-to-frame rate),
    - ``sigma_r`` = ``sqrt(sigma_x2 + sigma_y2 + sigma_z2)`` (total spread, Bohr).

    ``ipr`` is returned as all-NaN: the per-step WP-only density is not saved as a
    standalone field by the current pipeline (only total/system/delta VTIs), and
    deriving it is run-vintage dependent (system-vs-bath convention — see memory
    ``reference_canonical_bath_density``). Compute ``ipr`` separately with
    :func:`ipr` once a WP density frame series is available.

    Pure numpy + pandas (deps-clean).
    """
    import pandas as pd  # local import keeps module-load free of pandas

    run_dir = Path(run_dir)
    mom_path = Path(momentum_csv) if momentum_csv is not None else run_dir / _MOM_REL
    rs_path = Path(real_space_csv) if real_space_csv is not None else run_dir / _RS_REL

    mom = pd.read_csv(mom_path, comment="#")
    rs = pd.read_csv(rs_path, comment="#")

    # Per-step WP momentum distribution P_t(|k|): n_wp ordered by k, stacked by step.
    steps = sorted(mom["step"].unique())
    P = np.stack([
        mom[mom["step"] == s].sort_values("k_bohr_inv")["n_wp"].to_numpy(dtype=float)
        for s in steps
    ])
    t_mom = np.array([float(mom[mom["step"] == s]["time_au"].iloc[0]) for s in steps])

    # Real-space spread on the SAME steps (intersection, to keep a common base).
    rs_by_step = rs.set_index("step")
    common = [s for s in steps if s in rs_by_step.index]
    keep = np.array([s in rs_by_step.index for s in steps])
    P, t_mom = P[keep], t_mom[keep]
    sigma_r = np.array([
        float(np.sqrt(rs_by_step.loc[s, "sigma_x2"]
                      + rs_by_step.loc[s, "sigma_y2"]
                      + rs_by_step.loc[s, "sigma_z2"]))
        for s in common
    ])

    kl_mom = kl_series(P, reference=reference)
    ipr_series = np.full(len(common), np.nan)
    return WPIntegrity(time_au=t_mom, kl_mom=kl_mom, sigma_r=sigma_r, ipr=ipr_series)
