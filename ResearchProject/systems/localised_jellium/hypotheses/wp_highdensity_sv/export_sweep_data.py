"""
export_sweep_data — snapshot the non-VTI observables of the WP S(v) sweep into
the tracked hypotheses tree, with the corrected total-energy columns attached.

The production runs live under scripts/wp_highdensity_sv/wp/results/, which the
repo gitignores wholesale (**/results/) because it also holds VTIs, checkpoints
and GIF reports. This script copies ONLY the light per-step data of the twelve
production points (sigma_WP in {0.5, 2, 3} x v in {2.0, 2.5, 3.0, 3.5}) into
sweep_data/<run_name>/ so the sweep is versioned and readable from any machine:

  - verbatim: every .csv/.txt in raw/observables (including resume segments
    observables.from<N>.csv etc.), run_summary.txt, rt_state.txt;
  - derived:  observables_corrected.csv — the segment-concatenated observables
    ledger merged with wp_hd_stopping.wp_kinetic_norm_correction(), i.e. the
    CAP norm-corrected total energy

        E_total_corrected = E_total_raw - occ * T1 * (1 - norm_WP)

    (INQ reports the WP orbital's kinetic energy norm-divided, so under a CAP
    the raw ledger is inflated by the surviving fraction's per-particle mean;
    see wp_hd_stopping.py for the derivation and the 2026-07-30 verification).

Columns added to observables_corrected.csv on top of the raw ledger:
  norm_wp, correction_ev, e_total_raw_ev, e_total_corrected_ev,
  wp_kinetic_bare_ev, energy_total_corrected (Ha, same convention as
  energy_total).

Self-check: the final-step corrected energy must reproduce the published
S_deposit_corrected = (E_total_corrected(t_f) - E_GS) / 25 Bohr in
sigma_sweep_S_deposit.csv for every run. The script fails loudly if not.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import wp_hd_stopping as W  # noqa: E402

K = W.K
OUT = HERE / "sweep_data"
L_SLAB = 2.0 * W.SLAB_HALF                      # 25 Bohr
E_GS_EV = W.E_GS_HA_DX040 * W.HA_TO_EV


def export_run(sigma: float, v: float) -> dict:
    name = W.name_for(v, sigma)
    src = W.WP_RESULTS / name
    obs_dir = src / "raw" / "observables"
    if not obs_dir.is_dir():
        raise FileNotFoundError(obs_dir)
    dst = OUT / name
    dst.mkdir(parents=True, exist_ok=True)

    # 1. verbatim non-VTI copies (segments included)
    for f in sorted(obs_dir.iterdir()):
        if f.suffix in (".csv", ".txt"):
            shutil.copy2(f, dst / f.name)
    for extra in ("run_summary.txt", "rt_state.txt"):
        p = src / extra
        if p.exists():
            shutil.copy2(p, dst / extra)

    # 2. corrected total energy, at full cadence, over all segments
    corr = W.wp_kinetic_norm_correction(src)
    obs = K._concat_segments(obs_dir, "observables")
    df = pd.merge(obs, corr[["step", "norm_wp", "correction_ev",
                             "e_total_raw_ev", "e_total_corrected_ev",
                             "wp_kinetic_bare_ev"]], on="step")
    df["energy_total_corrected"] = df["e_total_corrected_ev"] / W.HA_TO_EV
    df.to_csv(dst / "observables_corrected.csv", index=False)

    ef_cor = float(df["e_total_corrected_ev"].iloc[-1]) - E_GS_EV
    return {
        "sigma": sigma, "v": v, "name": name,
        "steps_done": int(df["step"].iloc[-1]),
        "S_deposit_corrected": ef_cor / L_SLAB,
    }


def main() -> None:
    ref = pd.read_csv(HERE / "sigma_sweep_S_deposit.csv")
    rows, failures = [], []
    for sigma in W.SIGMAS:
        for v in W.VELOCITIES:
            r = export_run(sigma, v)
            rows.append(r)
            m = ref[(np.isclose(ref.sigma, sigma)) & (np.isclose(ref.v, v))]
            if len(m) != 1:
                failures.append(f"{r['name']}: no reference row")
                continue
            want = float(m.S_deposit_corrected.iloc[0])
            got = r["S_deposit_corrected"]
            ok = np.isclose(got, want, rtol=0, atol=5e-4)
            print(f"{r['name']:>10}  steps={r['steps_done']:>5}  "
                  f"S_dep_corr={got:.4f}  ref={want:.4f}  {'OK' if ok else 'MISMATCH'}")
            if not ok:
                failures.append(f"{r['name']}: {got:.6f} != {want:.6f}")
    if failures:
        raise SystemExit("self-check FAILED:\n  " + "\n  ".join(failures))
    print(f"\nexported {len(rows)} runs to {OUT} — all match sigma_sweep_S_deposit.csv")


if __name__ == "__main__":
    main()
