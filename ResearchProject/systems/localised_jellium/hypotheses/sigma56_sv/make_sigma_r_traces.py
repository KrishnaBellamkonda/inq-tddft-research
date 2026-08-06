#!/usr/bin/env python
"""sigma_r(t) for every sigma_WP = 6 wavepacket run — and where it stops meaning anything.

Writes s6_sigma_r_traces.png (2x2, one panel per velocity) and
s6_sigma_r_traces_overlay.png (all four on one axis).

WHY THIS PLOT EXISTS. The effective-width label <sigma_r> is a time average of
    sigma_r(t) = sqrt(sigma_x^2 + sigma_y^2 + sigma_z^2).
Averaged over the WHOLE run it came out ~20 Bohr for a packet that launches at
7.35 and is 8-9 when it leaves the slab. These traces show why: once the CAP has
absorbed the packet (norm -> 1e-9) the second moment is being computed on a
numerically negligible remnant spread across the cell, and sigma_r runs away to
24-30 Bohr. Those steps carry equal weight in a plain mean.

Each panel therefore shows, together:
  * sigma_r(t) measured                          (solid)
  * sqrt(3) * sigma_d(t), the FREE-Gaussian law  (dashed) -- what an unabsorbed,
    undisturbed packet would do: sigma_d = sqrt(sigma^2/2 + t^2/(2 sigma^2))
  * the surviving norm on a log right-hand axis  (grey)
  * the in-slab transit window                   (shaded)
  * the point where norm drops below 1e-3        (vertical rule): to the right of
    it sigma_r is measuring numerical residue, not a projectile.
  * the three candidate averages as horizontal rules.

The free-Gaussian curve is the control: where the solid trace tracks the dashed
one the measurement is trustworthy, and where it departs upward while the norm
collapses it is not.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import s56_stopping as S                                       # noqa: E402

SIGMA = 6.0
VS = (2.0, 2.5, 3.0, 3.5)
COL = {2.0: "#4C72B0", 2.5: "#DD8452", 3.0: "#55A868", 3.5: "#C44E52"}
NORM_FLOOR = 1e-3          # below this, sigma_r is residue, not packet


def trace(v: float):
    rs = S._concat(S.run_dir(SIGMA, v, "wp") / "raw" / "observables",
                   "wp_real_space_stats")
    t = rs.time_au.to_numpy()
    sr = np.sqrt((rs.sigma_x2 + rs.sigma_y2 + rs.sigma_z2).to_numpy())
    nm = rs.norm_check.to_numpy() if "norm_check" in rs else np.ones_like(t)
    return t, sr, nm


def free_sigma_r(t: np.ndarray, sigma: float) -> np.ndarray:
    """sqrt(3) * sigma_d(t) for an isotropic free Gaussian."""
    return np.sqrt(3.0) * np.sqrt(sigma**2 / 2.0 + t**2 / (2.0 * sigma**2))


def main() -> int:
    # ---------- 2x2, one panel per velocity -------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.0), sharex=False)
    for ax, v in zip(axes.ravel(), VS):
        t, sr, nm = trace(v)
        ti, to = S.transit_window(v)
        m = (t >= ti) & (t <= to)
        w = nm / nm.sum()

        ax.plot(t, sr, color=COL[v], lw=1.6, label=r"$\sigma_r(t)$ measured")
        ax.plot(t, free_sigma_r(t, SIGMA), color="0.35", lw=1.1, ls="--",
                label=r"$\sqrt{3}\,\sigma_d(t)$ free Gaussian")
        ax.axvspan(ti, to, color="0.85", zorder=0, label="in-slab transit")

        dead = np.where(nm < NORM_FLOOR)[0]
        if dead.size:
            ax.axvline(t[dead[0]], color="#d62728", lw=1.0, ls=":")
            ax.axvspan(t[dead[0]], t[-1], color="#d62728", alpha=0.06, zorder=0)

        ax.axhline(sr.mean(), color=COL[v], lw=0.9, ls="-.", alpha=0.9)
        ax.axhline(sr[m].mean() if m.any() else np.nan, color="#2ca02c",
                   lw=0.9, ls=":")
        ax.axhline((sr * w).sum(), color="#7b3ba8", lw=0.9, ls=":")

        ax.set_title(rf"$\sigma_{{\mathrm{{WP}}}}=6$, $v={v}$   "
                     rf"($\langle\sigma_r\rangle$ full {sr.mean():.1f}, "
                     rf"transit {sr[m].mean():.1f}, "
                     rf"norm-wtd {(sr*w).sum():.1f})", fontsize=9)
        ax.set_xlabel("t (a.u.)")
        ax.set_ylabel(r"$\sigma_r$ (Bohr)")

        axn = ax.twinx()
        axn.semilogy(t, np.maximum(nm, 1e-12), color="0.55", lw=0.9, alpha=0.8)
        axn.axhline(NORM_FLOOR, color="#d62728", lw=0.7, ls=":")
        axn.set_ylabel(r"$\|\psi_{\rm WP}\|^2$", color="0.45", fontsize=8)
        axn.tick_params(axis="y", labelsize=7, colors="0.45")

    axes[0, 0].legend(fontsize=7, frameon=False, loc="upper left")
    fig.suptitle(r"$\sigma_r(t)$ for the $\sigma_{\mathrm{WP}}=6$ wavepackets — "
                 r"dash-dot = full-run mean, dotted green = transit, "
                 r"dotted purple = norm-weighted", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = HERE / "s6_sigma_r_traces.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # ---------- all four on one axis --------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    for v in VS:
        t, sr, nm = trace(v)
        ok = nm >= NORM_FLOOR
        ax.plot(t[ok], sr[ok], color=COL[v], lw=1.7, label=f"v={v}  (norm > 1e-3)")
        ax.plot(t[~ok], sr[~ok], color=COL[v], lw=1.0, ls=":", alpha=0.5)
    tt = np.linspace(0, 180, 400)
    ax.plot(tt, free_sigma_r(tt, SIGMA), color="0.3", lw=1.2, ls="--",
            label=r"$\sqrt{3}\,\sigma_d(t)$ free Gaussian")
    ti, to = S.transit_window(VS[0])
    ax.axvspan(ti, to, color="0.85", zorder=0)
    ax.set_xlabel("t (a.u.)")
    ax.set_ylabel(r"$\sigma_r$ (Bohr)")
    ax.set_title(r"$\sigma_{\mathrm{WP}}=6$: solid where the packet still exists, "
                 r"dotted once the CAP has taken it", fontsize=9)
    ax.legend(fontsize=7.5, frameon=False)
    fig.tight_layout()
    out2 = HERE / "s6_sigma_r_traces_overlay.png"
    fig.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out2}")

    # ---------- provenance -------------------------------------------------
    print(f"\n{'v':>4} {'t_f':>7} {'sr(0)':>7} {'sr@exit':>8} "
          f"{'t(norm<1e-3)':>13} {'sr there':>9} {'sr(t_f)':>8}")
    for v in VS:
        t, sr, nm = trace(v)
        ti, to = S.transit_window(v)
        m = (t >= ti) & (t <= to)
        dead = np.where(nm < NORM_FLOOR)[0]
        td = t[dead[0]] if dead.size else np.nan
        sd = sr[dead[0]] if dead.size else np.nan
        print(f"{v:>4.1f} {t[-1]:>7.1f} {sr[0]:>7.2f} "
              f"{sr[m][-1] if m.any() else np.nan:>8.2f} {td:>13.1f} "
              f"{sd:>9.2f} {sr[-1]:>8.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
