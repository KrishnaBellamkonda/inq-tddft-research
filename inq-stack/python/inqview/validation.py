"""Minimum-observable-set post-run validator (ADR 0006).

Reads a run's ``results/observables_manifest.json`` (written by the C++ run at
startup) and checks every declared observable at four tiers:

  1. existence — the file (and CSV column) is present
  2. schema    — declared columns present; non-empty
  3. finite    — no NaN/Inf in numeric data
  4. invariant — the manifest-declared physical check (drift, norm band, zero-at-t0…)

A run is PASS iff every REQUIRED observable passes tiers 1–3 and any declared
tier-4 invariant holds. Optional observables are validated only if present and
never fail the run. Deps-clean (json + numpy + pandas) so a headless node can
audit a run; VTI-backed observables get tier-1 (existence via glob) only — their
schema/finite/invariant tiers need VTK and are reported as skipped.

CLI:  python -m inqview.validation <run_dir>   (exit 0 PASS, 1 FAIL/no-manifest)
"""
from __future__ import annotations

import glob
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import pandas as pd

HA_TO_MHA = 1000.0


@dataclass(frozen=True)
class TierResult:
    tier: str            # existence | schema | finite | invariant
    passed: bool
    detail: str = ""


@dataclass(frozen=True)
class ObservableResult:
    name: str
    required: bool
    tiers: list[TierResult] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return all(t.passed for t in self.tiers)


@dataclass(frozen=True)
class ValidationReport:
    run_dir: str
    run_type: str
    observables: list[ObservableResult]
    note: str = ""

    @property
    def passed(self) -> bool:
        # Required observables must pass all their tiers; optional never fail.
        return all(o.passed for o in self.observables if o.required)

    def summary(self) -> str:
        head = f"{'PASS' if self.passed else 'FAIL'}  {self.run_type}  ({self.run_dir})"
        lines = [head]
        for o in self.observables:
            bad = [t for t in o.tiers if not t.passed]
            if bad:
                tag = "FAIL" if o.required else "warn"   # optional misses don't fail the run
                opt = "" if o.required else " (optional)"
                for t in bad:
                    lines.append(f"  {tag} [{t.tier:9s}] {o.name}{opt}: {t.detail}")
            else:
                lines.append(f"  ok        {o.name}")
        if self.note:
            lines.append(f"  note: {self.note}")
        return "\n".join(lines)


# ── tier-4 invariant registry ──────────────────────────────────────────────
def _inv_drift_max(df: pd.DataFrame, col: str, spec: dict) -> tuple[bool, str]:
    v = df[col].to_numpy(dtype=float)
    drift_mha = float(np.nanmax(np.abs(v - v[0]))) * HA_TO_MHA
    lim = float(spec["value_mHa"])
    return drift_mha <= lim, f"drift {drift_mha:.3g} mHa vs limit {lim}"


def _inv_norm_band(df: pd.DataFrame, col: str, spec: dict) -> tuple[bool, str]:
    v = df[col].to_numpy(dtype=float)
    lo, hi = float(spec["lo"]), float(spec["hi"])
    ok = bool(np.all((v >= lo) & (v <= hi)))
    return ok, f"{col} range [{np.nanmin(v):.4g},{np.nanmax(v):.4g}] vs [{lo},{hi}]"


def _inv_zero_at_t0(df: pd.DataFrame, col: str, spec: dict) -> tuple[bool, str]:
    v0 = float(df[col].to_numpy(dtype=float)[0])
    tol = float(spec.get("atol", 1e-9))
    return abs(v0) <= tol, f"{col}(t0)={v0:.3g} vs atol {tol}"


def _inv_value_band(df: pd.DataFrame, col: str, spec: dict) -> tuple[bool, str]:
    return _inv_norm_band(df, col, spec)


def _inv_monotone(df: pd.DataFrame, col: str, spec: dict) -> tuple[bool, str]:
    v = df[col].to_numpy(dtype=float)
    ok = bool(np.all(np.diff(v) >= -float(spec.get("atol", 0.0))))
    return ok, f"{col} monotone non-decreasing"


_INVARIANTS: dict[str, Callable[[pd.DataFrame, str, dict], tuple[bool, str]]] = {
    "drift_max": _inv_drift_max,
    "norm_band": _inv_norm_band,
    "zero_at_t0": _inv_zero_at_t0,
    "value_band": _inv_value_band,
    "monotone_increasing": _inv_monotone,
}


def _check_csv(obs: dict, run_dir: Path) -> list[TierResult]:
    out: list[TierResult] = []
    path = run_dir / "results" / obs["file"]
    col = obs.get("column")
    # tier 1: existence
    if not path.exists():
        out.append(TierResult("existence", False, f"missing file {obs['file']}"))
        return out
    df = pd.read_csv(path, comment="#")
    if col is not None and col not in df.columns:
        out.append(TierResult("existence", False, f"missing column {col}"))
        return out
    out.append(TierResult("existence", True))
    # tier 2: schema (declared columns present) + non-empty
    schema = obs.get("schema")
    if schema is not None:
        missing = [c for c in schema if c not in df.columns]
        out.append(TierResult("schema", not missing,
                              "" if not missing else f"missing columns {missing}"))
    if len(df) == 0:
        out.append(TierResult("schema", False, "empty (0 rows)"))
        return out
    out.append(TierResult("schema", True, f"{len(df)} rows"))
    # tier 3: finite
    num = df.select_dtypes("number")
    nonfinite = int(np.sum(~np.isfinite(num.to_numpy(dtype=float))))
    out.append(TierResult("finite", nonfinite == 0,
                          "" if nonfinite == 0 else f"{nonfinite} NaN/Inf values"))
    if nonfinite:
        return out
    # tier 4: declared invariant
    inv = obs.get("invariant")
    if inv is not None:
        kind = inv.get("kind")
        icol = inv.get("col", col)
        fn = _INVARIANTS.get(kind)
        if fn is None:
            out.append(TierResult("invariant", True, f"skipped: unknown kind {kind!r}"))
        elif icol is None or icol not in df.columns:
            out.append(TierResult("invariant", False, f"invariant col {icol!r} absent"))
        else:
            ok, detail = fn(df, icol, inv)
            out.append(TierResult("invariant", ok, detail))
    return out


def _check_vti(obs: dict, run_dir: Path) -> list[TierResult]:
    # existence via glob; deeper tiers need VTK -> reported as skipped.
    frames = glob.glob(str(run_dir / "results" / obs["path"]))
    if not frames:
        return [TierResult("existence", False, f"no frames at {obs['path']}")]
    return [TierResult("existence", True, f"{len(frames)} frames"),
            TierResult("invariant", True, "skipped: VTI schema/finite/invariant need VTK")]


def _check_text(obs: dict, run_dir: Path) -> list[TierResult]:
    path = run_dir / "results" / obs.get("file", obs.get("path", ""))
    ok = path.exists() and path.stat().st_size > 0
    return [TierResult("existence", ok, "" if ok else f"missing/empty {path.name}")]


def validate_run(run_dir) -> ValidationReport:
    """Validate a run against its observable manifest. Returns a ValidationReport."""
    run_dir = Path(run_dir)
    man_path = run_dir / "results" / "observables_manifest.json"
    if not man_path.exists():
        return ValidationReport(str(run_dir), "unknown", [],
                                note="no observables_manifest.json (run predates ADR-0006)")
    manifest = json.loads(man_path.read_text())
    results: list[ObservableResult] = []
    for obs in manifest.get("observables", []):
        fmt = obs.get("format", "csv")
        if fmt == "csv":
            tiers = _check_csv(obs, run_dir)
        elif fmt == "vti":
            tiers = _check_vti(obs, run_dir)
        else:
            tiers = _check_text(obs, run_dir)
        results.append(ObservableResult(obs["name"], bool(obs.get("required", True)), tiers))
    return ValidationReport(str(run_dir), manifest.get("run_type", "unknown"), results)


def main(argv: Optional[list[str]] = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if not argv:
        print("usage: python -m inqview.validation <run_dir>")
        return 2
    report = validate_run(argv[0])
    print(report.summary())
    return 0 if report.passed and not report.note.startswith("no observables") else 1


if __name__ == "__main__":
    raise SystemExit(main())
