#!/usr/bin/env python3
"""check_twin — the twin-run-generation gate.

Validates that two runs form a legal *twin pair* and carry the full energy
decomposition, then writes ``twin_manifest.json``. This is the guarantee the
twin-run-analysis skill relies on; run it immediately after a paired dispatch.

Checks (all must pass):
  1. Both runs completed         (run_summary: run_completed = true)
  2. Config parity               (periodicity, Lz, spacing, N, sigma_WP, launch_z,
                                  gs_dir agree within tolerance)
  3. Projectile actually differs (the twins are NOT identical runs)
  4. Full decomposition present  (observables.csv has energy_total/kinetic/hartree/
                                  xc/external in BOTH twins)
  5. U_proj_bg available         (classical: energy_proj_bg_ideal column OR
                                  U_proj_bg in run_summary)

Self-contained: only the stdlib (own loose summary parser); no dependency on the
analysis skill, so the generation skill ships alone.

Usage:
  check_twin.py <pair_dir>                       # expects <pair_dir>/wp, <pair_dir>/classical
  check_twin.py --wp DIR --classical DIR [--pair-dir DIR-for-manifest]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

OBS_REL = "raw/observables/observables.csv"
SUMMARY = "run_summary.txt"

PARITY_FIELDS = ("periodicity", "lz", "spacing", "n", "sigma_wp", "launch_z", "gs_dir")
PARITY_ATOL = 1e-4
REQUIRED_ENERGY_COLS = (
    "energy_total", "energy_kinetic", "energy_hartree", "energy_xc", "energy_external",
)
# Additional per-step columns a DYNAMIC (Rung-2) pair must carry.
DYNAMIC_CLASSICAL_COLS = ("energy_proj_bg_ideal", "proj_z")
DYNAMIC_WP_COLS = ("wp_centroid_z",)


def infer_representation(d: dict) -> str | None:
    if "representation" in d:
        return d["representation"].lower()
    blob = " ".join(str(v) for v in d.values()).lower()
    if "perturbation" in blob:
        return "perturbation"
    if "ghost" in blob or "pseudopotential" in blob or ".upf" in blob or "z_valence" in blob:
        return "pseudopotential"
    if "wavepacket" in blob:
        return "wavepacket"
    return None


def parse_summary(path: Path) -> dict:
    """Flat lowercased dict from a free-form run_summary.txt (key=value or key:value)."""
    out: dict[str, str] = {}
    for k, v in re.findall(r"([A-Za-z_]\w*)\s*[=:]\s*(\S+)", path.read_text()):
        out[k.lower()] = v
    return out


def _num(d: dict, *keys):
    for k in keys:
        if k in d:
            try:
                return float(d[k])
            except ValueError:
                return d[k]
    return None


def header_cols(path: Path) -> list[str]:
    with path.open() as f:
        return next(csv.reader(f))


def check(wp_dir: Path, cl_dir: Path, dynamic: bool = False) -> tuple[bool, list[str], dict]:
    errs: list[str] = []
    wp_sum = parse_summary(wp_dir / SUMMARY)
    cl_sum = parse_summary(cl_dir / SUMMARY)

    # 1. completion
    for tag, s in (("wp", wp_sum), ("classical", cl_sum)):
        if s.get("run_completed", "false").lower() != "true":
            errs.append(f"[{tag}] run_completed != true")

    # 2. parity
    checked = []
    for f in PARITY_FIELDS:
        a, b = _num(wp_sum, f), _num(cl_sum, f)
        if a is None or b is None:
            continue
        checked.append(f)
        same = (math.isclose(a, b, abs_tol=PARITY_ATOL)
                if isinstance(a, (int, float)) and isinstance(b, (int, float)) else a == b)
        if not same:
            errs.append(f"parity {f}: wp={a!r} classical={b!r}")

    # 3. projectile must differ (guard against dispatching two identical runs)
    wp_proj = wp_sum.get("projectile") or wp_sum.get("mode")
    cl_proj = cl_sum.get("projectile") or cl_sum.get("mode") or cl_sum.get("run")
    if wp_proj is not None and wp_proj == cl_proj:
        errs.append("projectile field identical in both twins (not a real pair)")

    # 4. full decomposition columns present in both
    for tag, d in (("wp", wp_dir), ("classical", cl_dir)):
        obs = d / OBS_REL
        if not obs.exists():
            errs.append(f"[{tag}] missing {OBS_REL}")
            continue
        cols = set(header_cols(obs))
        missing = [c for c in REQUIRED_ENERGY_COLS if c not in cols]
        if missing:
            errs.append(f"[{tag}] observables missing columns: {missing}")

    # 5. U_proj_bg available on the classical side (observables.csv column, the
    #    auxiliary projectile.csv column, or the run_summary constant)
    cl_obs = cl_dir / OBS_REL
    cl_cols = set(header_cols(cl_obs)) if cl_obs.exists() else set()
    wp_cols = set(header_cols(wp_dir / OBS_REL)) if (wp_dir / OBS_REL).exists() else set()
    proj_csv = cl_dir / OBS_REL.replace("observables.csv", "projectile.csv")
    proj_cols = set(header_cols(proj_csv)) if proj_csv.exists() else set()
    has_col = ("energy_proj_bg_ideal" in cl_cols) or ("energy_proj_bg_ideal" in proj_cols)
    has_sum = ("u_proj_bg_ev" in cl_sum) or ("u_proj_bg_ha" in cl_sum)
    if not (has_col or has_sum):
        errs.append("classical U_proj_bg unavailable (no energy_proj_bg_ideal column in "
                    "observables.csv/projectile.csv, no U_proj_bg in run_summary)")

    # 6. dynamic (Rung-2) trajectory: the classical twin must emit projectile.csv
    #    (step,proj_z,proj_vz,energy_proj_ke,energy_proj_bg_ideal), an auxiliary CSV
    #    alongside observables.csv.
    if dynamic:
        proj_csv = cl_dir / OBS_REL.replace("observables.csv", "projectile.csv")
        if not proj_csv.exists():
            errs.append(f"[classical] dynamic run missing {proj_csv.name} (trajectory/proj-KE)")
        else:
            pc = set(header_cols(proj_csv))
            for c in ("proj_z", "energy_proj_ke", "energy_proj_bg_ideal"):
                if c not in pc:
                    errs.append(f"[classical] projectile.csv missing column: {c}")

    manifest = {
        "wp_dir": str(wp_dir), "classical_dir": str(cl_dir),
        "representation": infer_representation(cl_sum),
        "dynamic": dynamic,
        "parity_fields_checked": checked,
        "shared_config": {f: _num(cl_sum, f) for f in PARITY_FIELDS if _num(cl_sum, f) is not None},
        "wp_projectile": wp_proj, "classical_projectile": cl_proj,
        "u_proj_bg_source": "column" if has_col else ("summary" if has_sum else None),
        "valid": not errs,
    }
    return (not errs), errs, manifest


def main(argv=None):
    ap = argparse.ArgumentParser(description="Validate a twin run pair; write twin_manifest.json.")
    ap.add_argument("pair", nargs="?", help="pair dir containing wp/ and classical/")
    ap.add_argument("--wp"); ap.add_argument("--classical")
    ap.add_argument("--pair-dir", help="where to write twin_manifest.json (default: pair dir or CWD)")
    ap.add_argument("--dynamic", action="store_true",
                    help="also require the Rung-2 per-step trajectory/proj-KE columns")
    args = ap.parse_args(argv)

    if args.pair:
        pair = Path(args.pair); wp, cl = pair / "wp", pair / "classical"
    elif args.wp and args.classical:
        wp, cl = Path(args.wp), Path(args.classical)
        pair = Path(args.pair_dir) if args.pair_dir else Path.cwd()
    else:
        ap.error("give <pair_dir> or --wp DIR --classical DIR")

    ok, errs, manifest = check(wp, cl, dynamic=args.dynamic)
    out = (Path(args.pair_dir) if args.pair_dir else pair if args.pair else Path.cwd()) / "twin_manifest.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2))

    if ok:
        print(f"TWIN OK  representation={manifest['representation']}  "
              f"{'dynamic' if args.dynamic else 'static'}  "
              f"parity={manifest['parity_fields_checked']}  "
              f"U_proj_bg={manifest['u_proj_bg_source']}\nwrote {out}")
        return 0
    print("TWIN FAIL:", file=sys.stderr)
    for e in errs:
        print("  -", e, file=sys.stderr)
    print(f"wrote {out}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
