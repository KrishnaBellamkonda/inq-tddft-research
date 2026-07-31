#!/usr/bin/env python3
"""
Rewind a wavepacket run's OUTPUTS to its last retained checkpoint, so LJ_RESUME=1
can restart cleanly.

WHY THIS IS NEEDED. Checkpoints are written every CKPT_EVERY steps, but
observables and VTI frames are written continuously. A run that dies between
checkpoints therefore has outputs AHEAD of the state it can resume from:

    s3p0_v2p0   died at step 2461,  last checkpoint 2172   -> 289 orphan steps

On resume the engine recomputes 2172 -> 2461 and hits the frames already on disk:

    VTIImageDataWriter: file already exists and overwrite=false:
      results/s3p0_v2p0/raw/vti/density_total/density_t002184.vti

overwrite=false is deliberate on resume (final-timestep-checkpoint.md: segments
must not clobber each other), so the run aborts with SIGABRT. Worse, had it not
aborted, steps 2172-2460 would appear TWICE in the concatenated CSVs -- once in
the base segment and once in the .from2172 segment -- silently double-counting.

WHAT THIS DOES. For each run, reads last_step from rt_state.txt, then:
  * deletes VTI frames with step > last_step in every field directory;
  * truncates the base CSVs to rows with step <= last_step (comments/header kept);
  * deletes stale .from<N>.csv segments left by the failed resume attempt.

Everything removed is exactly what the resume recomputes. The checkpoint itself
and density_gs_system (the t=0 reference) are never touched.

Usage:  rewind_to_checkpoint.py <results_root> <run_name> [...] [--apply]
        (default is a DRY RUN)
"""
import re
import sys
from pathlib import Path

STEP_RE = re.compile(r"_t(\d+)\.vti$")
SEG_RE = re.compile(r"\.from\d+\.csv$")
CSVS = ("observables.csv", "interactions.csv",
        "wp_momentum_stats.csv", "wp_real_space_stats.csv")
# density_gs_system is the t=0 bath reference, written once; never rewind it.
SKIP_VTI_DIRS = {"density_gs_system"}


def last_step_of(run: Path) -> int | None:
    rt = run / "rt_state.txt"
    if not rt.exists():
        return None
    m = re.search(r"last_step\s*=\s*(\d+)", rt.read_text())
    return int(m.group(1)) if m else None


def rewind(run: Path, apply: bool) -> None:
    last = last_step_of(run)
    if last is None:
        print(f"{run.name}: SKIP — no readable rt_state.txt")
        return
    print(f"\n{run.name}: rewinding to last_step={last}")

    vti_root = run / "raw" / "vti"
    for d in sorted(p for p in vti_root.iterdir() if p.is_dir()):
        if d.name in SKIP_VTI_DIRS:
            continue
        doomed = []
        for f in d.iterdir():
            m = STEP_RE.search(f.name)
            if m and int(m.group(1)) > last:
                doomed.append(f)
        if doomed and apply:
            for f in doomed:
                f.unlink()
        print(f"  vti/{d.name:22s} remove {len(doomed):5d} frames past {last}")

    # Base CSVs AND any .from<N> segments. A segment written by an EARLIER resume
    # can hold rows that are still valid (steps <= the new last_step) as well as
    # orphans past it, so segments are TRUNCATED like the base file, not deleted
    # wholesale. Only a segment that begins entirely after last_step is dropped.
    # (Deleting them unconditionally was wrong once the rolling checkpoint had
    # advanced past the segment's start -- it would have thrown away real data.)
    obs = run / "raw" / "observables"
    targets = list(CSVS) + sorted(p.name for p in obs.iterdir()
                                  if SEG_RE.search(p.name))
    for name in targets:
        f = obs / name
        if not f.exists():
            continue
        kept, dropped, out = 0, 0, []
        for line in f.read_text().splitlines(keepends=True):
            s = line.lstrip()
            if s.startswith("#") or s.startswith("step,"):
                out.append(line)
                continue
            head = s.split(",", 1)[0]
            try:
                step = int(head)
            except ValueError:
                out.append(line)
                continue
            if step <= last:
                out.append(line)
                kept += 1
            else:
                dropped += 1
        if kept == 0 and SEG_RE.search(name):
            # Segment lies entirely past the checkpoint — nothing of it survives.
            print(f"  {name:26s} DROP (starts after {last})")
            if apply:
                f.unlink()
            continue
        if apply and dropped:
            f.write_text("".join(out))
        print(f"  {name:26s} keep {kept:5d} rows, drop {dropped:5d}")


def main() -> int:
    args = [a for a in sys.argv[1:] if a != "--apply"]
    apply = "--apply" in sys.argv
    root = Path(args[0])
    for name in args[1:]:
        rewind(root / name, apply)
    print("\n" + ("APPLIED" if apply else "DRY RUN — pass --apply to execute"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
