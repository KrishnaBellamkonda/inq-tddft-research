#!/usr/bin/env python3
"""Figures + notebooks for the CYLINDRICAL PROXIMITY LADDER.

Plan: docs/plans/cylindrical-proximity-ladder.md

REUSE, NOT REIMPLEMENTATION
---------------------------
The channeling twin's `build_report_figures.py` already draws every panel the
user asked for — interaction energies, trajectory, T1/T2/var(p), classical
1/2 m v^2, the loss-definition overlay, momentum-loss maps, and S from all three
estimators with uncertainties. Its run sets are just a dict of
{tag: {wp_results, cl_results, wp_name, cl_name}}, so a LADDER RUNG is simply
another tag. This module therefore imports that builder and drives it, rather
than forking 1300 lines that would then drift out of sync.

The two-pass driver matters here more than it did for the twin: `collect_limits`
is computed across ALL rungs at once, so every rung's figures share axis limits
and are directly comparable by eye. A per-rung autoscale would silently rescale
the very effect the ladder exists to show.

DEGRADES, DOES NOT DIE
----------------------
The stage is chained with `afterany`, so it must cope with a rung whose runs
failed or are still missing: such a rung is REPORTED and skipped, and the rest
are written up. A ladder analysis that refuses to run because one of five rungs
is absent would waste the four that worked.

THE COUPLING COORDINATE
-----------------------
Rungs are compared against the MEASURED occupancy the packet sees, not against
the nominal R_in, because the bore is not empty: rung r10's ground state already
has 16 of its 160 electrons inside it, so "distance to the background edge" is
not "distance to the electrons". `f_wall(t)` from radial_occupancy is the measured
proxy, and it exists at every rung INCLUDING the filled one (where `f_bore` is
identically zero and meaningless).

WHAT IS AND IS NOT A CONFOUND HERE (corrected 2026-08-03 against measured data)
------------------------------------------------------------------------------
It is tempting to argue that a fixed-TIME fit window cannot compare rungs because
the packet spreads and "the rungs merge in time". That argument is WRONG for this
ladder, and acting on it would have discarded a valid comparison:

  * every rung shares sigma_WP = 4, dt, n_steps and v0, so the packet spreads
    IDENTICALLY across the ladder. At a fixed t, all rungs have the same packet
    size and differ ONLY in where the wall sits — the coupling difference is the
    independent variable, not a contaminant;
  * velocity is not the issue either. Measured v/v0 runs 0.99 -> 0.93 over
    t = 9-25 at r08, comfortably inside the light-projectile v >= 0.85 v0
    criterion (.claude/rules/light-projectile-stopping.md).

The REAL limitation is coupling drift WITHIN the window: f_wall moved 4.4x at r10
and 2.0x at r08 over t = 9-25. So each fitted S is an average over a RANGE of
couplings, and the range differs by rung. Constraining that range to 1.5x would
require t < 2.4 a.u. (~120 steps), too short for a stable fit — the tension is
intrinsic to a spreading wavepacket, not something a better window can remove.

Hence `rung_summary_row` reports mean AND span of f_wall over the window, and
`L04_ratio_vs_measured_coupling` draws S(coupling) with horizontal error bars.
The trend is a curve over coupling ranges, not a set of points at single values.
"""
from __future__ import annotations

import argparse
import importlib
import json
import os
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
SYS_DIR = REPO / "ResearchProject/systems/cylindrical_jellium"
SCRIPTS = SYS_DIR / "scripts"
TWIN_HYP = SYS_DIR / "hypotheses/channeling_twin"

# the twin builder is the drawing engine; its directory also carries refined.py,
# channeling_stopping.py and ks_stopping.py which Ctx re-imports per run set.
if str(TWIN_HYP) not in sys.path:
    sys.path.insert(0, str(TWIN_HYP))

OUT_ROOT = HERE / "figures"

# R_in in Bohr, for the secondary (readability) axis. sigma_WP = 4.
RUNG_RIN = {"r10": 10.0, "r08": 8.0, "r06": 6.0, "r04": 4.0, "r00": 0.0,
            "r04n160": 4.0}
RUNG_ORDER = ["r10", "r08", "r06", "r04", "r00"]
# The primary S fit window, inherited from the channeling twin's "T1 9-25".
# Used here only to report the COUPLING sampled over it (see rung_summary_row).
FIT_T0, FIT_T1 = 9.0, 25.0
SIGMA_WP = 4.0


def rung_run_sets(rungs: list[str]) -> dict:
    """{tag: cfg} in the shape `build_report_figures.Ctx` expects.

    r10 is the COMPLETED channeling twin and lives in a different script folder
    under different run names; every ladder rung shares one pair of folders and
    is distinguished by run name (= the rung tag).
    """
    out = {}
    for r in rungs:
        if r == "r10":
            out[r] = dict(wp_results=SCRIPTS / "channeling_twin/wp/results",
                          cl_results=SCRIPTS / "channeling_twin/classical/results",
                          wp_name="wp", cl_name="classical")
        else:
            out[r] = dict(wp_results=SCRIPTS / "proximity_ladder/wp/results",
                          cl_results=SCRIPTS / "proximity_ladder/classical/results",
                          wp_name=r, cl_name=r)
    return out


def rung_is_present(tag: str, cfg: dict) -> tuple[bool, str]:
    """Both halves must carry a completed run_summary before we try to load them."""
    for half, root, name in (("wp", cfg["wp_results"], cfg["wp_name"]),
                             ("classical", cfg["cl_results"], cfg["cl_name"])):
        s = Path(root) / name / "run_summary.txt"
        if not s.is_file():
            return False, f"{half}: no run_summary.txt at {s}"
        if "run_completed = true" not in s.read_text():
            return False, f"{half}: run at {s.parent} did not complete"
    return True, "ok"


# ---------------------------------------------------------------------------
# cross-rung summary
# ---------------------------------------------------------------------------

def rung_summary_row(tag: str, C) -> dict:
    """One row of the ladder's headline table, from an already-built Ctx.

    `f_wall` is the measured fraction of the wavepacket inside the jellium — the
    coupling coordinate. It is read at t=0 and at the END of the fit window, so
    the table shows BOTH the nominal starting coupling and how far the packet had
    already spread by the time S was measured.
    """
    wp, cl = C.wp, C.cl
    row = {"rung": tag, "R_in": RUNG_RIN.get(tag, np.nan),
           "R_in_over_sigma": RUNG_RIN.get(tag, np.nan) / SIGMA_WP,
           "n_steps": len(wp), "t_end_au": float(wp.t.iloc[-1])}

    # measured coupling, if radial occupancy was written
    for col, key in (("f_wall", "f_wall_t0"), ("f_bore", "f_bore_t0")):
        if col in wp.columns:
            row[key] = float(wp[col].iloc[0])
            row[key.replace("_t0", "_end")] = float(wp[col].iloc[-1])

    # COUPLING OVER THE FIT WINDOW — the honest x-axis for S.
    #
    # Every rung shares sigma_WP, dt, n_steps and v0, so the packet spreads
    # IDENTICALLY across the ladder: a fixed-time window compares the same packet
    # stage at different wall positions, which is a controlled comparison and not
    # a confound. Velocity is likewise not the issue (0.99 -> 0.93 v0 across the
    # window at r08, so the light-projectile v >= 0.85 v0 criterion is satisfied).
    #
    # What IS true is that the coupling DRIFTS WITHIN the window — measured 4.4x
    # at r10 and 2.0x at r08 over t = 9-25 — because the packet spreads into the
    # wall as it flies. So a fitted S is an AVERAGE over a range of couplings, and
    # the range differs by rung. Holding the coupling to within 1.5x of its initial
    # value would need t < 2.4 a.u. (~120 steps), far too short for a stable fit:
    # the tension is intrinsic to a spreading wavepacket, not an analysis defect.
    #
    # Therefore report the mean AND the span, so S(coupling) is drawn with
    # horizontal error bars rather than as points at a single coupling.
    if "f_wall" in wp.columns:
        t = wp.t.to_numpy()
        fw = wp.f_wall.to_numpy()
        m = (t >= FIT_T0) & (t <= FIT_T1)
        if m.any():
            row["fw_fit_mean"] = float(fw[m].mean())
            row["fw_fit_lo"] = float(fw[m].min())
            row["fw_fit_hi"] = float(fw[m].max())
            row["fw_fit_drift"] = float(fw[m].max() / max(fw[m].min(), 1e-12))

    # energy loss over the run, both representations
    row["dE_wp_T1_ev"] = float(wp.T1_drift_ev.iloc[0] - wp.T1_drift_ev.iloc[-1])
    row["dE_wp_T2_ev"] = float(wp.T2_total_ev.iloc[0] - wp.T2_total_ev.iloc[-1])
    row["dE_cl_ke_ev"] = float(cl.ke_ev.iloc[0] - cl.ke_ev.iloc[-1])
    row["frac_loss_cl"] = row["dE_cl_ke_ev"] / float(cl.ke_ev.iloc[0])

    # stopping powers from the shared fit table Ctx already computed
    if getattr(C, "fits", None) is not None and len(C.fits):
        f = C.fits
        for _, r in f.iterrows():
            lbl = str(r.get("window", "")) + "|" + str(r.get("estimator", ""))
            for src in ("S_wp", "S_cl", "S_wp_err", "S_cl_err"):
                if src in r.index:
                    row[f"{src}[{lbl}]"] = float(r[src])
    return row


def draw_ladder_comparison(rows: pd.DataFrame, outdir: Path, dpi: int) -> list[str]:
    """The cross-rung headline: how the classical/quantum agreement moves.

    Plotted against the MEASURED coupling where available, with R_in/sigma as a
    secondary readability axis (see module docstring for why nominal R_in alone
    is not the honest x-axis).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir.mkdir(parents=True, exist_ok=True)
    names: list[str] = []
    if rows.empty:
        return names

    r = rows.sort_values("R_in", ascending=False)
    x = r["R_in_over_sigma"].to_numpy()

    def save(fig, stem):
        p = outdir / f"{stem}.png"
        fig.savefig(p, dpi=dpi)
        plt.close(fig)
        names.append(p.name)

    # 1. fractional energy loss — the "weak -> strong" claim, measured
    fig, ax = plt.subplots(figsize=(3.6, 3.2), constrained_layout=True)
    ax.plot(x, 100 * r["frac_loss_cl"], "o-", color="tab:blue", label="classical")
    if "dE_wp_T1_ev" in r:
        ke0 = r["dE_cl_ke_ev"] / r["frac_loss_cl"]
        ax.plot(x, 100 * r["dE_wp_T1_ev"] / ke0, "s--", color="tab:red",
                label=r"wavepacket $T_1$")
    ax.set_xlabel(r"$R_\mathrm{in}/\sigma_\mathrm{WP}$")
    ax.set_ylabel("energy lost over the run  (%)")
    ax.invert_xaxis()                      # weak (left) -> strong (right)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.3)
    save(fig, "L01_fractional_energy_loss")

    # 2. measured coupling: what the packet actually sits in
    if "f_wall_t0" in r.columns:
        fig, ax = plt.subplots(figsize=(3.6, 3.2), constrained_layout=True)
        ax.plot(x, 100 * r["f_wall_t0"], "o-", color="tab:green", label="t = 0")
        ax.plot(x, 100 * r["f_wall_end"], "s--", color="tab:olive", label="end of run")
        ax.set_xlabel(r"$R_\mathrm{in}/\sigma_\mathrm{WP}$")
        ax.set_ylabel(r"WP charge inside the jellium, $f_\mathrm{wall}$  (%)")
        ax.invert_xaxis()
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=0.3)
        # the rungs MERGE in time — this is the panel that shows it
        ax.set_title("rungs separate early, merge late", fontsize=8)
        save(fig, "L02_measured_coupling")

    # 3b. S ratio against the MEASURED coupling, with horizontal spans
    if "fw_fit_mean" in r.columns:
        sc = next((c for c in r.columns if c.startswith("S_wp[") and "T1" in c), None)
        cc = sc.replace("S_wp[", "S_cl[") if sc else None
        if sc and cc in r.columns:
            fig, ax = plt.subplots(figsize=(4.0, 3.2), constrained_layout=True)
            ax.errorbar(r["fw_fit_mean"], r[sc] / r[cc],
                        xerr=[r["fw_fit_mean"] - r["fw_fit_lo"],
                              r["fw_fit_hi"] - r["fw_fit_mean"]],
                        fmt="o-", color="tab:purple", lw=1.2, ms=5, capsize=3)
            for _, q in r.iterrows():
                ax.annotate(q["rung"], (q["fw_fit_mean"], q[sc] / q[cc]),
                            textcoords="offset points", xytext=(6, 5), fontsize=7)
            ax.axhline(1.0, color="k", lw=0.8, ls=":")
            ax.set_xscale("log")
            ax.set_xlabel(r"mean $f_\mathrm{wall}$ over the fit window")
            ax.set_ylabel(r"$S_\mathrm{WP}/S_\mathrm{classical}$  ($T_1$)")
            ax.set_title("bars = coupling drift within the window", fontsize=8)
            ax.grid(alpha=0.3)
            save(fig, "L04_ratio_vs_measured_coupling")

    # 3. the headline ratio, per available fit window
    scols = [c for c in r.columns if c.startswith("S_wp[")]
    if scols:
        fig, ax = plt.subplots(figsize=(4.0, 3.2), constrained_layout=True)
        for sc in scols:
            cc = sc.replace("S_wp[", "S_cl[")
            if cc not in r.columns:
                continue
            lbl = sc[len("S_wp["):-1]
            ax.plot(x, r[sc] / r[cc], "o-", label=lbl, lw=1.2, ms=4)
        ax.axhline(1.0, color="k", lw=0.8, ls=":")
        ax.set_xlabel(r"$R_\mathrm{in}/\sigma_\mathrm{WP}$")
        ax.set_ylabel(r"$S_\mathrm{WP}/S_\mathrm{classical}$")
        ax.invert_xaxis()
        ax.legend(frameon=False, fontsize=6)
        ax.grid(alpha=0.3)
        save(fig, "L03_stopping_ratio_vs_rung")

    return names


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rungs", default="r10,r08,r06,r04,r00")
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--skip-per-rung", action="store_true",
                    help="only rebuild the cross-rung comparison")
    a = ap.parse_args()

    import build_report_figures as B          # the drawing engine
    from inqview.visualisation import style

    style.apply_theme()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    wanted = [r.strip() for r in a.rungs.split(",") if r.strip()]
    cfgs = rung_run_sets(wanted)

    present, missing = {}, {}
    for tag, cfg in cfgs.items():
        ok, why = rung_is_present(tag, cfg)
        (present if ok else missing).__setitem__(tag, cfg if ok else why)

    print(f"rungs requested : {wanted}")
    print(f"rungs present   : {sorted(present)}")
    for tag, why in missing.items():
        print(f"  [SKIP] {tag}: {why}")
    if not present:
        print("\nNo rung has both halves complete — nothing to draw.")
        (OUT_ROOT / "manifest.json").write_text(
            json.dumps({"present": [], "missing": missing}, indent=1))
        return 0

    # ---- load every present rung, then share axis limits across all of them --
    ctxs, failed = {}, {}
    for tag in [t for t in RUNG_ORDER if t in present]:
        try:
            print(f"loading {tag} ...")
            ctxs[tag] = B.Ctx(tag, present[tag])
            w = ctxs[tag].wp
            print(f"  {len(w)} steps, t = {w.t.iloc[0]:.2f}..{w.t.iloc[-1]:.2f} a.u.")
        except Exception as exc:                       # noqa: BLE001
            failed[tag] = f"{type(exc).__name__}: {exc}"
            print(f"  [FAIL] {tag}: {failed[tag]}")
            traceback.print_exc()

    manifest: dict[str, object] = {"missing": missing, "load_failed": failed}
    if not ctxs:
        (OUT_ROOT / "manifest.json").write_text(json.dumps(manifest, indent=1))
        print("\nEvery present rung failed to load.")
        return 1

    if not a.skip_per_rung:
        print(f"\npass 1: shared axis limits across {len(ctxs)} rung(s) — so the "
              "rungs are comparable BY EYE, not just by number")
        lims = B.collect_limits(ctxs)
        clim = B.momentum_map_clims(ctxs)

        for tag, C in ctxs.items():
            d = OUT_ROOT / tag
            print(f"pass 2: writing {tag}")
            try:
                names = B.draw_set(C, lims, d, a.dpi)
                names += B.draw_momentum_maps(C, clim, d, a.dpi)
                manifest[tag] = names
                print(f"  wrote {len(names)} figures to {d}")
            except Exception as exc:                   # noqa: BLE001
                manifest[tag] = f"draw failed: {type(exc).__name__}: {exc}"
                print(f"  [FAIL] drawing {tag}: {exc}")
                traceback.print_exc()

        (OUT_ROOT / "shared_limits.json").write_text(
            json.dumps({"axis_limits": lims, "momentum_clim": clim}, indent=1))

    # ---- the cross-rung comparison ------------------------------------------
    print("\ncross-rung comparison")
    rows = pd.DataFrame([rung_summary_row(t, C) for t, C in ctxs.items()])
    rows.to_csv(OUT_ROOT / "ladder_summary.csv", index=False)
    manifest["comparison"] = draw_ladder_comparison(rows, OUT_ROOT / "comparison",
                                                    a.dpi)
    print(rows.to_string(index=False))

    (OUT_ROOT / "manifest.json").write_text(json.dumps(manifest, indent=1, default=str))
    n = sum(len(v) for v in manifest.values() if isinstance(v, list))
    print(f"\n{n} figures under {OUT_ROOT}")
    print(f"summary table: {OUT_ROOT / 'ladder_summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
