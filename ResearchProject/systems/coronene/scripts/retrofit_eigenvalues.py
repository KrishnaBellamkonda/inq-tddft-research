#!/usr/bin/env python3
"""retrofit_eigenvalues.py — copy eigenvalues.csv + occupations.csv into every
existing run's results/raw/observables/eigenvalues/.

Designed for the one-time retrofit after Phase-4: the eigenvalue writer
landed after most runs were produced, so every run loaded from a
checkpoint that didn't yet carry the CSVs. Once the three save_gs/<sig>/
scripts have been rebuilt and re-run (cheap, ~10 min each), this script
walks every run_*/results/run_summary.txt, parses the checkpoint_dir,
and copies the eigenvalue CSVs from there into the run's results tree.

Usage::

    python3 scripts/retrofit_eigenvalues.py [--dry-run]

No args = process every run_*/ sibling under
ResearchProject/systems/coronene/. ``--dry-run`` prints the planned
copies without writing.
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

CKPT_RE = re.compile(r"^checkpoint_dir\s*=\s*(\S+)\s*$")


def find_checkpoint_dir(summary_path: Path) -> Path | None:
    """Parse ``checkpoint_dir = ...`` from results/run_summary.txt."""
    if not summary_path.exists():
        return None
    for line in summary_path.read_text().splitlines():
        m = CKPT_RE.match(line.strip())
        if m:
            return Path(m.group(1))
    return None


def retrofit_one(run_dir: Path, *, dry_run: bool) -> str:
    summary = run_dir / "results" / "run_summary.txt"
    ckpt_dir = find_checkpoint_dir(summary)
    if ckpt_dir is None:
        return f"{run_dir.name}: skipped (no run_summary.txt or no checkpoint_dir)"
    if not ckpt_dir.exists():
        return f"{run_dir.name}: checkpoint missing on disk ({ckpt_dir})"

    out_dir = run_dir / "results" / "raw" / "observables" / "eigenvalues"
    copied: list[str] = []
    missing: list[str] = []
    for fname in ("eigenvalues.csv", "occupations.csv"):
        src = ckpt_dir / fname
        dst = out_dir / fname
        if not src.exists():
            missing.append(fname)
            continue
        if dry_run:
            copied.append(f"would copy {src} -> {dst}")
        else:
            out_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            copied.append(fname)
    msg = f"{run_dir.name}: copied={copied}"
    if missing:
        msg += f" missing={missing}"
    return msg


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args(argv)

    here = Path(__file__).resolve().parent.parent  # systems/coronene/
    run_dirs = sorted(d for d in here.iterdir()
                      if d.is_dir() and d.name.startswith("run_"))
    for run in run_dirs:
        # Skip the legacy run_propagate_paper_replica / run_save_gs_paper_replica
        # which are not part of the new framework.
        if run.name in ("run_propagate_paper_replica", "run_save_gs_paper_replica"):
            continue
        print(retrofit_one(run, dry_run=args.dry_run))
    return 0


if __name__ == "__main__":
    sys.exit(main())
