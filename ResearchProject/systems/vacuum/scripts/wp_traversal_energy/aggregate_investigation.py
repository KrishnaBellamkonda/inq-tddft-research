#!/usr/bin/env python3
"""Cross-run aggregation for the CAP energy-normalization investigation.

For every completed vacuum WP run, extract the decisive quantities that test the
hypothesis (docs/plans/cap-energy-normalization-validation.md):

  E_reported(0), E_reported(T), dE_reported   [Ha, eV]  -- the "shoots up" signal
  norm(0)=1, norm(T)                                    -- physical orbital norm
  E_ext(T) = E_reported(T)*norm(T)                      -- extensive energy
  reflection_residual                                    -- inner-region norm (Phase 1)

Physical norm per-step is taken from the momentum-stats norm_check (write_every=1),
normalized to its t=0 value (k-space integral is proportional to <psi|psi>).
e_kin_ha == energies.csv:kinetic is the per-particle mean KE (INQ's /norm quantity).

Writes: investigation_summary.csv (one row per run) to the results/ root.
"""
from __future__ import annotations
import sys, csv, glob
from pathlib import Path
import numpy as np
import pandas as pd

HA_EV = 27.211386

# run_tag -> (phase, group, param_label)
RUNS = {
    "nocap":            ("0", "baseline",   "no CAP"),
    "cap":              ("0", "baseline",   "1-sided CAP eta=-3.5"),
    "nocap_long":       ("0", "baseline",   "no CAP (1600 steps)"),
    "exp1a_eta-0.3":    ("1", "eta_sweep",  "eta=-0.3"),
    "exp1a_eta-0.7":    ("1", "eta_sweep",  "eta=-0.7"),
    "exp1a_eta-1.0":    ("1", "eta_sweep",  "eta=-1.0"),
    "exp1a_eta-2.0":    ("1", "eta_sweep",  "eta=-2.0"),
    "exp1a_eta-3.5":    ("1", "eta_sweep",  "eta=-3.5"),
    "exp2_N0p1":        ("2", "partial_abs","target N_abs~0.1"),
    "exp2_N0p3":        ("2", "partial_abs","target N_abs~0.3"),
    "exp2_N0p5":        ("2", "partial_abs","target N_abs~0.5"),
    "exp3a_mask_etrs":  ("3", "decisive",   "mask + ETRS (norm-losing)"),
    "exp3b_mask_cn":    ("3", "decisive",   "mask + CN (norm-preserving)"),
}


def _read_obs(run_dir: Path, name: str) -> pd.DataFrame | None:
    # comment lines start with '#'; header is the first non-# line
    p = run_dir / "raw" / "observables" / name
    if not p.exists():
        return None
    return pd.read_csv(p, comment="#")


def summarise(results_root: Path, tag: str, phase: str, group: str, param: str) -> dict | None:
    run_dir = results_root / tag
    if not (run_dir / "raw" / "observables").is_dir():
        return None
    en = _read_obs(run_dir, "energies.csv")
    mom = _read_obs(run_dir, "wp_momentum_stats.csv")
    if en is None or mom is None or len(en) < 2:
        return None
    E0 = float(en["kinetic"].iloc[0])
    ET = float(en["kinetic"].iloc[-1])
    # physical norm(t) = norm_check_mom(t) / norm_check_mom(0)
    nc = mom["norm_check"].to_numpy(dtype=float)
    normT = float(nc[-1] / nc[0]) if nc[0] != 0 else np.nan
    Eext_T = ET * normT
    row = {
        "run": tag, "phase": phase, "group": group, "param": param,
        "E_rep0_Ha": E0, "E_repT_Ha": ET, "dE_rep_Ha": ET - E0,
        "dE_rep_eV": (ET - E0) * HA_EV,
        "dE_rep_pct": 100.0 * (ET - E0) / E0 if E0 else np.nan,
        "norm_T": normT, "frac_absorbed": 1.0 - normT,
        "E_ext_T_Ha": Eext_T, "E_ext_T_eV": Eext_T * HA_EV,
        "E_ext_frac": Eext_T / E0 if E0 else np.nan,
        "n_steps": int(en["step"].iloc[-1]),
    }
    return row


def main(argv):
    results_root = Path(argv[1]) if len(argv) > 1 else Path("results")
    rows = []
    for tag, (phase, group, param) in RUNS.items():
        r = summarise(results_root, tag, phase, group, param)
        if r is None:
            print(f"[agg] skip {tag} (no observables yet)")
            continue
        rows.append(r)
        print(f"[agg] {tag:18s} dE_rep={r['dE_rep_eV']:+7.2f} eV "
              f"({r['dE_rep_pct']:+5.1f}%)  norm_T={r['norm_T']:.3f}  "
              f"E_ext/E0={r['E_ext_frac']:.3f}")
    if not rows:
        print("[agg] no completed runs found")
        return 1
    df = pd.DataFrame(rows).sort_values(["phase", "run"])
    out = results_root / "investigation_summary.csv"
    df.to_csv(out, index=False)
    print(f"[agg] wrote {out}  ({len(df)} runs)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
