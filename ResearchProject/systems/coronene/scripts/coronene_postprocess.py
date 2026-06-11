#!/usr/bin/env python3
"""coronene_postprocess.py — thin CLI wrapper around inqview.pipeline.

Two subcommands:

  * ``run`` — process a single run's ``results/`` tree, producing every
    spec-required figure / GIF / summary.

  * ``hypothesis`` — collate the ``results/`` trees of several runs into a
    single ``hypotheses/<NN>_*/`` folder of comparison artefacts.

Coronene-specific defaults (cmap, overlap-axis range, fixed dt for the
density GIF time-axis labels) are baked in; everything else is delegated to
``inqview.pipeline``.

Usage examples
--------------

  # Single-run postprocess (auto-detects run name from parent dir)
  python coronene_postprocess.py run --results <abs_path_to_run_dir>/results

  # Force-rebuild every figure
  python coronene_postprocess.py run --results <...>/results --rebuild

  # Run only the screens phase
  python coronene_postprocess.py run --results <...>/results --phases screens

  # Hypothesis: collate three runs into hypotheses/01_wp_energy_spread/
  python coronene_postprocess.py hypothesis \\
      --hypothesis-dir .../hypotheses/01_wp_energy_spread \\
      --runs run_E30=<.../run_E30> run_base=<.../run_base> run_E800=<.../run_E800>
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Tsubonoya base parameters expressed in atomic units. Used for default plot
# axis labels (so we don't have to read run_summary.txt to label time-stepped
# plots). The C++ template also uses these.
_DT_AU_DEFAULT = 0.020
_WRITE_EVERY_DEFAULT = 10


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="coronene_postprocess",
                                description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser("run", help="postprocess a single run's results/")
    p_run.add_argument("--results", type=Path, required=True,
                       help="path to <run_dir>/results")
    p_run.add_argument("--run-name", type=str, default=None,
                       help="run name shown in plot titles (defaults to "
                            "parent directory name)")
    p_run.add_argument("--phases", type=str, default=None,
                       help="comma-separated subset of phases "
                            "(default: all). Valid: summary,gs,observables,"
                            "density,screens,overlap,orbitals,paraview")
    p_run.add_argument("--rebuild", action="store_true",
                       help="regenerate even if outputs exist")
    p_run.add_argument("--with-paraview", action="store_true",
                       help="enable the paraview phase (slow; off by default)")
    p_run.add_argument("--dt-au", type=float, default=_DT_AU_DEFAULT,
                       help="propagation dt in atomic units (used for time "
                            "labels on density frames; default 0.020)")
    p_run.add_argument("--write-every", type=int, default=_WRITE_EVERY_DEFAULT,
                       help="density write_every (default 10)")
    p_run.add_argument("--percentile", type=float, default=99.0,
                       help="percentile clip for density GIF colour scale "
                            "(default 99)")

    p_h = sub.add_parser("hypothesis",
                         help="collate several runs into a hypotheses folder")
    p_h.add_argument("--hypothesis-dir", type=Path, required=True,
                     help="path to hypotheses/<NN>_*/")
    p_h.add_argument("--runs", nargs="+", type=str, required=True,
                     help="list of label=<abs_path_to_run_dir> entries")
    p_h.add_argument("--rebuild", action="store_true")

    return p


def cmd_run(args: argparse.Namespace) -> int:
    from inqview.pipeline import run as pipeline_run
    log = logging.getLogger("coronene_postprocess")
    if not log.handlers:
        log.setLevel(logging.INFO)
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter("%(message)s"))
        log.addHandler(h)

    phases = (args.phases.split(",")
              if args.phases else None)
    res = pipeline_run(
        args.results,
        run_name=args.run_name,
        phases=phases,
        rebuild=args.rebuild,
        skip_paraview=not args.with_paraview,
        dt_au=args.dt_au,
        write_every=args.write_every,
        percentile=args.percentile,
        logger=log,
    )
    print()
    print(f"=== {res.run_name}: postprocess summary ===")
    for ph in res.ok:
        print(f"  ok    {ph} ({res.durations_s.get(ph, 0):.1f}s)")
    for ph, why in res.skipped.items():
        print(f"  skip  {ph} — {why}")
    for ph, why in res.failed.items():
        print(f"  FAIL  {ph} — {why}")
    return 0 if res.succeeded() else 1


def cmd_hypothesis(args: argparse.Namespace) -> int:
    from inqview.pipeline.compare import run_hypothesis
    runs: list[tuple[str, Path]] = []
    for spec in args.runs:
        if "=" not in spec:
            print(f"bad --runs entry {spec!r}: must be label=path",
                  file=sys.stderr)
            return 2
        label, path = spec.split("=", 1)
        runs.append((label, Path(path)))
    out = run_hypothesis(args.hypothesis_dir, runs, rebuild=args.rebuild)
    print()
    print(f"=== hypothesis: {args.hypothesis_dir.name} ===")
    for k, v in out.items():
        print(f"  {k}: {v}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.cmd == "run":
        return cmd_run(args)
    if args.cmd == "hypothesis":
        return cmd_hypothesis(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())
