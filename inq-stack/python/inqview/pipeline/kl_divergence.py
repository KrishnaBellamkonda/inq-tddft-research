"""Phase: ``kl_divergence`` — WP momentum-distribution drift from t=0.

Computes the Kullback-Leibler divergence of the WP momentum distribution
P(|k|, t) relative to its initial-state reference P(|k|, 0):

    KL(P_t || P_0) = sum_k P_t(k) * log( P_t(k) / P_0(k) )

KL = 0 at t=0 by construction; rises monotonically (or near-monotonically)
as the jellium bath redistributes the WP's momentum content. Used in the
final stopping-power rollup to compare how aggressively each WP run
spreads relative to its launch state, and to contextualise the (more
quantitative) Knudsen `<|p|^2>/2` curve from ``knudsen_ke``.

Notes:

* The 1D radial histogram in ``momentum_distribution.csv`` makes this
  an angle-averaged ("isotropised") KL distance, which is the natural
  invariant under WP rotation and matches what the rollup needs.
* Bins where ``P_0(k) <= 0`` are dropped from the sum (cannot define
  log(P_t / 0)). A tiny floor (``EPS = 1e-300``) is also added to P_t
  to keep ``P_t log P_t`` finite when a bin transiently drops to zero.
* If only the native ``wp_momentum_stats.csv`` is present (no histogram),
  KL cannot be computed — the file lacks the per-bin distribution, only
  the moments. Phase skips cleanly in that case.

Outputs:

* ``analysis/observables/kl_divergence.csv`` — ``step,time_au,kl_div``.
* ``analysis/observables/kl_divergence_vs_t.png``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# TODO: Explain to me, what reallt is this import
from . import _common


EPS = 1e-300

# TODO: Can the KL divergence be calculated as a time series? How expensive would
# it be?

# TODO: KL value could actually be a metric that can quantify how much of the wavepacket
# is preserved. A good sanity check on the shape of the orbital. this is because
# the idea is taken from information theory and might be interesting. 

# TODO: Runfeng a while ago suggested, I make contour visualisations for the wavepacket. 
# Need to brainstorm how this can be doe. 

def _normalise(p: np.ndarray) -> np.ndarray:
    """Return p / sum(p) (safe; returns zeros if sum is zero)."""
    s = float(p.sum())
    if s <= 0:
        return np.zeros_like(p)
    return p / s


def _kl(p: np.ndarray, q: np.ndarray) -> float:
    """sum_k p_k log(p_k / q_k), with q_k > 0 only."""
    mask = q > 0
    if not mask.any():
        return float("nan")
    pp = np.where(mask, p, 0.0) + EPS
    qq = np.where(mask, q, 1.0)
    contrib = np.where(mask, pp * np.log(pp / qq), 0.0)
    return float(contrib.sum())


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    md = results_dir / "raw" / "observables" / "momentum_distribution.csv"
    if not md.exists():
        return {"skipped": f"missing: {md}"}

    df = pd.read_csv(md, comment="#")
    if "n_wp" not in df.columns or df["n_wp"].abs().sum() == 0:
        return {"skipped": "momentum_distribution.csv has no n_wp data"}

    # Pivot to (time x k-bin); rows = ordered by time_au, cols = k_bohr_inv.
    times = sorted(df["time_au"].unique())
    k_vals = sorted(df["k_bohr_inv"].unique())
    grid = (df.pivot_table(index="time_au", columns="k_bohr_inv",
                           values="n_wp", aggfunc="sum")
              .reindex(index=times, columns=k_vals).to_numpy())
    if grid.shape[0] < 2:
        return {"skipped": "need at least 2 time samples for KL(t)"}

    # Reference: P_0 = WP distribution at first recorded time.
    ref = _normalise(grid[0])
    if ref.sum() == 0:
        return {"skipped": "P_0 is zero everywhere (WP not visible at t=0)"}

    steps_per_time = (df.groupby("time_au")["step"].first()
                        .reindex(times).to_numpy().astype(int))

    rows = []
    for i, t in enumerate(times):
        pt = _normalise(grid[i])
        rows.append((int(steps_per_time[i]), float(t), _kl(pt, ref)))
    kl_df = pd.DataFrame(rows, columns=["step", "time_au", "kl_div"])

    out_dir = _common.ensure_dir(results_dir / "analysis" / "observables")
    artefacts: list[str] = []

    csv_out = out_dir / "kl_divergence.csv"
    if _common.need_rebuild(csv_out, rebuild):
        kl_df.to_csv(csv_out, index=False)
    artefacts.append(str(csv_out))

    ifw = _common.post_ifw_window_from_summary(results_dir)

    p_t = out_dir / "kl_divergence_vs_t.png"
    if _common.need_rebuild(p_t, rebuild):
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(kl_df["time_au"], kl_df["kl_div"], "C4-", lw=1.4,
                label="KL$(P_t \\| P_0)$")
        if ifw is not None:
            _common.ifw_highlight(ax, ifw[0])
            ax.legend(loc="best", fontsize=9)
        ax.set_xlabel("time (a.u.)")
        ax.set_ylabel(r"KL$(P_t \,||\, P_0)$  (nats)")
        ax.set_title(_common.title(run_name,
            "WP momentum-distribution KL divergence from t=0"))
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(p_t, dpi=150)
        plt.close(fig)
        artefacts.append(str(p_t))

    return {"n_steps": int(len(kl_df)),
            "kl_max": float(kl_df["kl_div"].max()),
            "artefacts": artefacts}
