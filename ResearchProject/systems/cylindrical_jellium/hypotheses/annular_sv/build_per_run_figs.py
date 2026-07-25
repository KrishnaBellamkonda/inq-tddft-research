#!/usr/bin/env python3
"""Driver: generate the per-run figure battery for every projectile run and write
a manifest (per_run_manifest.json) the report notebook embeds.

Runs `per_run.generate` for each of the 10 runs (9 classical + 1 WP). Figure paths
in the manifest are RELATIVE to this folder so the notebook embeds them portably.

    venv/bin/python3 build_per_run_figs.py        # generate all (heavy: GIFs)
"""
from __future__ import annotations
import json
import os
import sys
from pathlib import Path

HYP = Path(__file__).resolve().parent
SWEEP = HYP.parent.parent / "annular_sv"
sys.path.insert(0, str(HYP))
import per_run  # noqa: E402

FIGROOT = HYP / "per_run_figs"
MANIFEST = HYP / "per_run_manifest.json"

VELS = {"v0p15": 0.15, "v0p30": 0.30, "v0p45": 0.45}
RUNS = [(f"rs{rs}_{vk}", rs, vv, "classical")
        for rs in (6, 4, 2) for vk, vv in VELS.items()]
RUNS.append(("wp_rs6_v0p30", 6, 0.30, "wp"))


def rel(p):
    return os.path.relpath(str(p), str(HYP)) if p else None


def _launch_z(summ_path):
    import re
    m = re.search(r"launch_z\s*=\s*(-?\d+(?:\.\d+)?)", summ_path.read_text())
    return float(m.group(1)) if m else None


def main():
    manifest = {}
    for label, rs, v0, rtype in RUNS:
        run_dir = SWEEP / label
        summ = next(run_dir.glob("**/run_summary.txt"), None)
        if summ is None:
            print(f"[skip] {label}: no run_summary.txt")
            continue
        results = str(summ.parent)
        print(f"==== {label} (rs={rs}, v0={v0}, {rtype}) ====", flush=True)
        out = per_run.generate(label, rs, v0, str(run_dir), results,
                               str(FIGROOT), rtype)
        st = out["stopping"]
        if st:
            st = {**st, "path": rel(st["path"])}
        manifest[label] = {
            "rs": rs, "v0": v0, "rtype": rtype, "launch_z": _launch_z(summ),
            "matrix": [[c, k, rel(p), cap] for c, k, p, cap in out["matrix"]],
            "carpets": [[cap, rel(p)] for cap, p in out["carpets"]],
            "stopping": st,
            "fft": [[cap, rel(p)] for cap, p in out.get("fft", [])],
            "pipeline": {g: [[n, rel(p)] for n, p in v]
                         for g, v in out["pipeline"].items()},
        }
        MANIFEST.write_text(json.dumps(manifest, indent=1))  # checkpoint each run
        print(f"[done] {label}: {len(manifest[label]['matrix'])} gifs, "
              f"{sum(len(v) for v in manifest[label]['pipeline'].values())} pipeline figs",
              flush=True)
    print(f"\nWROTE {MANIFEST}  ({len(manifest)} runs)")


if __name__ == "__main__":
    main()
