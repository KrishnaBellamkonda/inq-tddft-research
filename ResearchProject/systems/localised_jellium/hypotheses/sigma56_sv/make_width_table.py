#!/usr/bin/env python
"""Per-run effective-width table — where <sigma_r> actually comes from.

Writes s6_effective_width_table.md / .csv.

THE POINT OF THIS TABLE. The upstream label <sigma_r> is a FULL-RUN time average
of sigma_r(t) = sqrt(sigma_x^2 + sigma_y^2 + sigma_z^2). For a packet that the CAP
removes early, most of that average is taken over steps where there is no packet
left: norm falls to ~1e-9, and the second moment is then measured on a numerically
negligible residue smeared across the whole cell. sigma_r(t_f) reaches 24-30 Bohr
for sigma_WP = 6 purely for that reason, and it drags the full-run mean to ~20 --
which is why the sigma = 6 curve was labelled <sigma_r> ~ 20 despite the packet
being ~7.3 Bohr wide at launch and only ~8-9 Bohr when it leaves the slab.

Columns, so the choice of window is visible rather than buried:

  sigma_r(0)      analytic at launch: sqrt(3)*sigma_WP/sqrt(2)
  sigma_r@exit    when the centroid clears the far slab face
  sigma_r(t_f)    last step -- dominated by the absorbed remnant, NOT physical
  <sigma_r> full  upstream definition: plain mean over every step
  <sigma_r> transit   mean over the in-slab window only
  <sigma_r> norm-wtd  mean weighted by the surviving packet norm

The last two are the ones that describe the object the slab actually saw.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import s56_stopping as S                                      # noqa: E402

VS = (2.0, 2.5, 3.0, 3.5)


def row(sigma: float, v: float) -> dict | None:
    try:
        rs = S._concat(S.run_dir(sigma, v, "wp") / "raw" / "observables",
                       "wp_real_space_stats")
    except Exception:
        return None
    if rs.empty or not {"sigma_x2", "sigma_y2", "sigma_z2"} <= set(rs.columns):
        return None
    t = rs.time_au.to_numpy()
    sr = np.sqrt((rs.sigma_x2 + rs.sigma_y2 + rs.sigma_z2).to_numpy())
    nm = rs.norm_check.to_numpy() if "norm_check" in rs else np.ones_like(t)
    ti, to = S.transit_window(v)
    m = (t >= ti) & (t <= to)
    try:
        p = S.measure(sigma, v, "wp")
        steps, done = p.steps_target, p.steps_done
        complete = p.complete
    except Exception:
        steps = done = -1
        complete = False
    w = nm / nm.sum() if nm.sum() > 0 else np.full_like(nm, 1.0 / len(nm))
    m1, t1, n1, hit1 = S.sigma_r_window(rs)      # the adopted 1 %-norm-loss window
    return {
        "sigma_WP": sigma, "v": v,
        "run": S.run_name(sigma, v, "wp"),
        "steps": f"{done}/{steps}", "complete": complete,
        "t_final_au": round(float(t[-1]), 1),
        "sigma_r_0": round(float(sr[0]), 2),
        "sigma_r_exit": round(float(sr[m][-1]), 2) if m.any() else np.nan,
        "sigma_r_tfinal": round(float(sr[-1]), 2),
        "mean_full": round(float(sr.mean()), 2),
        "mean_transit": round(float(sr[m].mean()), 2) if m.any() else np.nan,
        "mean_normwtd": round(float((sr * w).sum()), 2),
        "mean_1pct": round(m1, 2), "t_1pct": round(t1, 1),
        "steps_1pct": n1, "reached_1pct": hit1,
        "t_slab_exit": round(float(to), 1),
        "norm_final": float(nm[-1]),
    }


def main() -> int:
    rows = [r for s in (5.0, 6.0) for v in VS if (r := row(s, v))]
    d = pd.DataFrame(rows)
    d.to_csv(HERE / "s6_effective_width_table.csv", index=False)

    lines = ["# Effective width per run — sigma56_sv wavepackets\n",
             "`sigma_r(t) = sqrt(sigma_x^2+sigma_y^2+sigma_z^2)`; all widths in Bohr.\n",
             "**⟨σ_r⟩ 1%-window** is the adopted label: the mean of σ_r(t) from",
             "t = 0 until the packet norm has fallen by 1 % (user decision,",
             "2026-08-03). `t(1%)` is where that window closes; compare it with",
             "`t_exit`, the in-slab transit end — the window should cover the",
             "crossing. `⟨σ_r⟩ full` is the superseded full-run mean, which keeps",
             "averaging long after the CAP has destroyed the packet (see",
             "`norm(t_f)`) and so measures a smeared remnant, not the projectile.\n",
             "| run | v | steps | t_f | σ_r(0) | σ_r@exit | σ_r(t_f) | **⟨σ_r⟩ 1%-window** | t(1%) | t_exit | ⟨σ_r⟩ full | norm(t_f) |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(
            f"| `{r['run']}` | {r['v']} | {r['steps']} | {r['t_final_au']} | "
            f"{r['sigma_r_0']} | {r['sigma_r_exit']} | {r['sigma_r_tfinal']} | "
            f"**{r['mean_1pct']}** | {r['t_1pct']} | {r['t_slab_exit']} | "
            f"{r['mean_full']} | {r['norm_final']:.1e} |")

    lines.append("\n## Per-sigma means (what a legend label would read)\n")
    lines.append("| σ_WP | **⟨σ_r⟩ 1%-window** | ⟨σ_r⟩ full (superseded) | n runs |")
    lines.append("|---|---|---|---|")
    for s, g in d.groupby("sigma_WP"):
        lines.append(f"| {s:g} | **{g.mean_1pct.mean():.2f}** | "
                     f"{g.mean_full.mean():.2f} | {len(g)} |")
    (HERE / "s6_effective_width_table.md").write_text("\n".join(lines) + "\n")

    print("\n".join(lines))
    print(f"\nwrote {HERE/'s6_effective_width_table.md'}")
    print(f"wrote {HERE/'s6_effective_width_table.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
