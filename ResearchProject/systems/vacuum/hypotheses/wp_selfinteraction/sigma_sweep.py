"""How the wavepacket self-interaction scales with sigma_WP.

Runs: `.../wp_selfinteraction/results/sweep_s{TAG}_{theory}`, produced by
`shared/bin/run-wp-si-sweep.slurm` (array 32625669).
Companion to `selfinteraction.py`, which analyses the single sigma = 4 case.

THE PROTOCOL, AND WHY IT MAKES THE SWEEP CLEAN
----------------------------------------------
Free evolution is invariant under `r -> lambda r`, `t -> lambda^2 t`. The sweep
scales EVERYTHING off sigma:

    L = 18 sigma      h = 0.125 sigma      dt = 0.00125 sigma^2      1500 steps

so every sigma is the SAME discrete problem: always a 144^3 grid, always
`sigma_dens(0)/h = 5.657`, always a box half-width of `5.99 sigma_dens(t_end)`,
always `t_end = 1.875 sigma^2` (i.e. the same DIMENSIONLESS time
`tau = t/sigma^2`). Discretisation error is therefore identical across sigma
rather than varying with it, and sigma = 4 reduces exactly to the completed
fixed-box run.

What is left varying is a single DIMENSIONLESS COUPLING. Kinetic energy scales
as 1/sigma^2, while both the Hartree self-energy (`E_PP ~ 1/sigma`) and the LDA
exchange potential (`v_x ~ n^(1/3) ~ 1/sigma`) scale as 1/sigma. Their ratio to
the kinetic term is therefore ~ sigma: a WIDER packet is MORE strongly
self-interacting per unit kinetic energy, at fixed tau.

E_PP(0) IS EXACTLY 1/sigma HERE — INCLUDING THE BOX TERM
--------------------------------------------------------
A Gaussian of density std `a = sigma/sqrt2` has free-space self-energy
`1/(2 a sqrt(pi)) = 1/(sigma sqrt(2 pi))`. INQ drops G=0 in a charged periodic
cell, which subtracts the Madelung term `xi/(2L)` with `xi = 2.8373`. Because
`L = 18 sigma` scales too, that term is ALSO ~ 1/sigma:

    E_PP(0) = [ 1/sqrt(2 pi) - xi/36 ] / sigma = 0.32013 / sigma  Ha

so `E_PP(0) * sigma` is a CONSTANT here — a sharp gate on the whole protocol.
(In a FIXED box it would instead be `A/sigma + C`; do not carry this formula
across to a fixed-box sweep.)

TWO WAYS TO READ THE SWEEP, AND THEY ANSWER DIFFERENT QUESTIONS
---------------------------------------------------------------
- **fixed tau** (`sigma_table`) — every run compared at the same DIMENSIONLESS
  spreading. This isolates the coupling and is the scientifically clean trend.
- **fixed physical time** (`at_fixed_physical_time`) — what a real run of a given
  duration actually suffers. Only available where `T/sigma^2 <= 1.875`, i.e.
  `sigma >= sqrt(T/1.875)`; for T = 30 a.u. that is sigma >= 4. Smaller sigma at
  fixed T is OUTSIDE these runs and would need longer ones — the function reports
  which sigma it had to drop rather than silently truncating.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]

HA_TO_EV = 27.211386245988
XI_MADELUNG = 2.837297          # simple-cubic (Wigner) Madelung constant
TAU_END = 1.875                 # t_end / sigma^2, fixed by the protocol
L_OVER_SIGMA = 18.0

RESULTS = Path(os.environ.get(
    "WPSI_RESULTS",
    REPO / "ResearchProject/systems/vacuum/scripts/wp_selfinteraction/results"))

SIGMAS = (1.0, 2.0, 3.0, 4.0, 6.0, 8.0)
THEORIES = ("noninteracting", "hartree", "lda", "sic_pzrun")


def tag_of(sigma: float) -> str:
    """`4.0 -> '4p0'` — the run-name encoding used by the sweep driver."""
    return f"{sigma:.1f}".replace(".", "p")


def epp0_predicted(sigma: float) -> float:
    """Analytic E_PP(0) in Ha for the SCALED box. See module docstring."""
    return (1.0 / math.sqrt(2.0 * math.pi) - XI_MADELUNG / (2.0 * L_OVER_SIGMA)) / sigma


def sigma_dens_free(t, sigma: float) -> np.ndarray:
    return np.sqrt(sigma**2 / 2.0 + np.asarray(t, float) ** 2 / (2.0 * sigma**2))


@dataclass
class SweepRun:
    sigma: float
    theory: str
    t: np.ndarray
    tau: np.ndarray             # t / sigma^2 — the protocol's natural clock
    sigma_iso: np.ndarray       # 3-D geometric-mean density width, Bohr
    var_p3d: np.ndarray
    e_pp: np.ndarray            # Ha
    complete: bool

    @property
    def sigma_scaled(self) -> np.ndarray:
        """Width in units of sigma. Identical across sigma if evolution is free."""
        return self.sigma_iso / self.sigma


def _read(obs: Path, stem: str) -> pd.DataFrame:
    files = sorted(obs.glob(f"{stem}*.csv"))
    if not files:
        raise FileNotFoundError(f"no {stem}*.csv under {obs}")
    df = pd.concat([pd.read_csv(f, comment="#") for f in files], ignore_index=True)
    return df.sort_values("step").drop_duplicates("step", keep="last").reset_index(drop=True)


def load(sigma: float, theory: str) -> SweepRun:
    run_dir = RESULTS / f"sweep_s{tag_of(sigma)}_{theory}"
    obs = run_dir / "raw" / "observables"
    if not obs.is_dir():
        raise FileNotFoundError(f"no observables under {run_dir}")
    pos, mom, ie = _read(obs, "wp_real_space_stats"), _read(obs, "wp_momentum_stats"), _read(obs, "interactions")
    t = pos["time_au"].to_numpy()
    smry = (run_dir / "run_summary.txt")
    done = smry.is_file() and "run_completed = true" in smry.read_text()
    return SweepRun(
        sigma=sigma, theory=theory, t=t, tau=t / sigma**2,
        sigma_iso=(pos["sigma_x_circ"] * pos["sigma_y_circ"]
                   * pos["sigma_z_circ"]).to_numpy() ** (1.0 / 3.0),
        var_p3d=(mom["sigma_px2"] + mom["sigma_py2"] + mom["sigma_pz2"]).to_numpy(),
        e_pp=ie["e_pp"].to_numpy(), complete=done)


def load_all() -> dict[tuple[float, str], SweepRun]:
    out = {}
    for s in SIGMAS:
        for th in THEORIES:
            try:
                out[(s, th)] = load(s, th)
            except FileNotFoundError:
                continue
    return out


# ---------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------

def protocol_gate(runs: dict) -> pd.DataFrame:
    """Per sigma: does the scaled protocol hold, and is the reference free?

    Three independent checks:
      `epp0_x_sigma`      must be the SAME number at every sigma (the scaling)
      `max_rel_sigma_err` reference width vs the closed form (the grid)
      `max_var_p_drift`   var(p) is EXACTLY conserved under free evolution
    """
    rows = []
    for s in SIGMAS:
        ref = runs.get((s, "noninteracting"))
        if ref is None:
            continue
        rel = np.abs(ref.sigma_iso / sigma_dens_free(ref.t, s) - 1.0)
        vp0 = 3.0 / (2.0 * s**2)                # 3-D: 3 * 1/(2 sigma^2)
        rows.append({
            "sigma": s,
            "epp0_ha": float(ref.e_pp[0]),
            "epp0_x_sigma": float(ref.e_pp[0] * s),
            "epp0_vs_predicted": float(ref.e_pp[0] / epp0_predicted(s)),
            "max_rel_sigma_err": float(rel.max()),
            "max_var_p_drift": float(np.abs(ref.var_p3d / vp0 - 1.0).max()),
            "complete": ref.complete,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the trend
# ---------------------------------------------------------------------------

def _excess_at(run: SweepRun, ref: SweepRun, tau: float) -> float:
    """Width relative to the non-interacting run at dimensionless time `tau`.

    Taken against the REFERENCE RUN rather than the analytic formula so that
    discretisation error — identical between them by construction — cancels.

    Returns NaN if EITHER run stops short of `tau`. This guard exists because its
    absence produced a wrong answer that looked right: analysing the sweep while
    the `lda` runs were still in flight, `argmin` silently returned their last
    written row (a SMALLER tau) and compared it against a completed `hartree` run
    at full tau. The resulting "LDA excess" was too low by ~40 %, was smooth and
    monotonic across sigma, and was mistaken for a cross-binary artefact.
    A partial run must therefore be a REFUSAL, never a nearest-neighbour match.
    """
    tol = 1e-9
    if run.tau[-1] < tau - tol or ref.tau[-1] < tau - tol:
        return float("nan")
    n = min(run.tau.size, ref.tau.size)
    i = int(np.argmin(np.abs(run.tau[:n] - tau)))
    return float(run.sigma_iso[i] / ref.sigma_iso[i])


def sigma_table(runs: dict, tau: float = TAU_END) -> pd.DataFrame:
    """The headline trend: excess spreading vs sigma, at fixed DIMENSIONLESS time.

    `xc_cancellation` is the fraction of the Hartree self-repulsion that LDA xc
    removes: `1 - (lda-1)/(hartree-1)`. If it is sigma-independent, the
    cancellation is a property of the functional; if it drifts, LDA correlation
    (which does not follow the 1/sigma similarity) is responsible.
    """
    rows = []
    for s in SIGMAS:
        ref = runs.get((s, "noninteracting"))
        if ref is None:
            continue
        row = {"sigma": s, "tau": tau}
        for th in ("hartree", "lda", "sic_pzrun"):
            r = runs.get((s, th))
            row[f"excess_{th}"] = _excess_at(r, ref, tau) if r is not None else np.nan
        h, l = row.get("excess_hartree", np.nan), row.get("excess_lda", np.nan)
        row["xc_cancellation"] = (1.0 - (l - 1.0) / (h - 1.0)
                                  if np.isfinite(h) and np.isfinite(l) and h > 1.0
                                  else np.nan)
        # the SIC run must land back on the reference: |excess - 1| is the residual
        row["sic_residual"] = (abs(row["excess_sic_pzrun"] - 1.0)
                               if np.isfinite(row.get("excess_sic_pzrun", np.nan))
                               else np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def at_fixed_physical_time(runs: dict, t_au: float = 30.0) -> tuple[pd.DataFrame, list]:
    """Excess spreading after `t_au` of REAL time — the practically useful view.

    Returns `(table, dropped)`. A sigma is DROPPED when `t_au` lies beyond that
    run's end (`t_au > 1.875 sigma^2`, i.e. sigma < sqrt(t_au/1.875)); extending
    to it would need longer runs. Reporting the dropped list is the point — a
    silent truncation would read as "covered every sigma".
    """
    rows, dropped = [], []
    for s in SIGMAS:
        ref = runs.get((s, "noninteracting"))
        if ref is None:
            continue
        if t_au > ref.t[-1] + 1e-9:
            dropped.append(s)
            continue
        tau = t_au / s**2
        row = {"sigma": s, "t_au": t_au, "tau": tau}
        for th in ("hartree", "lda", "sic_pzrun"):
            r = runs.get((s, th))
            row[f"excess_{th}"] = _excess_at(r, ref, tau) if r is not None else np.nan
        rows.append(row)
    return pd.DataFrame(rows), dropped
