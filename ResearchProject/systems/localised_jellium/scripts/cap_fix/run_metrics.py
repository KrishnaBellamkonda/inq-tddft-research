#!/usr/bin/env python3
"""LOCKED EVAL HARNESS (metric extractor) — cap_fix campaign.

Parses one completed cap_fix run and emits `METRIC name=value` lines for the
autoresearch loop. Per the autoresearch skill this file is OFF LIMITS to
experiments: experiments vary the CAP config, never how the metric is measured.

Metrics (all from raw per-step data, HA_EV = 27.211386):
  artifact_rise_eV  PRIMARY, lower is better. E_total(final) - E_total(min):
                    the post-minimum rise of the reported total energy. The
                    diagnosed artifact IS this rise (drain-then-rise); a clean
                    absorbing run is monotone -> rise ~ 0 (propagator floor).
  excursion_eV      max(0, max_t[E_total(t) - E_total(0)]): how far the ledger
                    climbs ABOVE the t=0 reference (the "crosses zero" form).
  drain_eV          E_total(min) - E_total(0): absorption actually happening
                    (sanity context; a CAP that absorbs nothing has rise 0 AND
                    drain 0 — that is a failed absorber, not a fix).
  absorbed_e        N(0) - N(final) from charge.csv (electrons removed).
  t_min_au          time of the E_total minimum.
"""
import sys
import pandas as pd

HA_EV = 27.211386


def main(out_dir: str) -> int:
    obs = pd.read_csv(f"{out_dir}/raw/observables/observables.csv")
    dE = (obs["energy_total"] - obs["energy_total"].iloc[0]) * HA_EV
    i_min = int(dE.idxmin())

    rise = float(dE.iloc[-1] - dE[i_min])
    excursion = float(max(0.0, dE.max()))
    drain = float(dE[i_min])
    t_min = float(obs["time_au"][i_min])

    absorbed = float("nan")
    try:
        ch = pd.read_csv(f"{out_dir}/raw/observables/charge.csv")
        # step 0 appears twice: a stale pre-propagator snapshot (WP not yet in
        # the density) and the in-propagator value. Keep the LAST of each step.
        ch = ch.drop_duplicates(subset="step", keep="last")
        absorbed = float(ch["n_total"].iloc[0] - ch["n_total"].iloc[-1])
    except Exception as exc:  # charge.csv is new; tolerate absence loudly
        print(f"WARNING charge.csv unreadable: {exc}", file=sys.stderr)

    print(f"METRIC artifact_rise_eV={rise:.6f}")
    print(f"METRIC excursion_eV={excursion:.6f}")
    print(f"METRIC drain_eV={drain:.6f}")
    print(f"METRIC absorbed_e={absorbed:.6f}")
    print(f"METRIC t_min_au={t_min:.3f}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: run_metrics.py <results_out_dir>", file=sys.stderr)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
