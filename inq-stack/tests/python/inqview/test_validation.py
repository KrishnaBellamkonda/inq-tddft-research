"""Tests for the minimum-observable-set validator (ADR 0006).

Builds synthetic run dirs (manifest + results CSVs) that each break exactly one
tier, with the expected verdict defined up front:
  complete      -> PASS
  dropped column-> FAIL tier 1 (existence)
  NaN injected  -> FAIL tier 3 (finite)
  norm = 2.0    -> FAIL tier 4 (norm_band invariant)
  big drift     -> FAIL tier 4 (drift_max invariant)
  no manifest   -> report.note set, not PASS
Pure json/numpy/pandas, portable.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from inqview.validation import validate_run

pytestmark = pytest.mark.io

_MANIFEST = {
    "run_type": "test-wp", "schema_version": 1, "write_every": 2, "n_steps": 4,
    "observables": [
        {"name": "energy_total", "required": True,
         "file": "raw/observables/observables.csv", "column": "energy_total",
         "format": "csv", "cadence": "step",
         "invariant": {"kind": "drift_max", "value_mHa": 1.0}},
        {"name": "wp_momentum_stats", "required": True,
         "file": "raw/observables/wp_momentum_stats.csv", "format": "csv",
         "schema": ["step", "time_au", "norm_check"],
         "invariant": {"kind": "norm_band", "col": "norm_check", "lo": 0.97, "hi": 1.03}},
    ],
}


def _write_run(tmp_path, *, energy=None, norm=None, drop_energy_col=False,
               nan=False, manifest=True):
    obs = tmp_path / "results" / "raw" / "observables"
    obs.mkdir(parents=True)
    if manifest:
        (tmp_path / "results" / "observables_manifest.json").write_text(json.dumps(_MANIFEST))
    e = energy if energy is not None else [-0.6156, -0.61559, -0.6156, -0.61561]
    if nan:
        e = [-0.6156, np.nan, -0.6156, -0.6156]
    edf = pd.DataFrame({"step": [0, 2, 4, 6], "time_au": [0.0, 0.04, 0.08, 0.12],
                        "energy_total": e})
    if drop_energy_col:
        edf = edf.drop(columns=["energy_total"])
    edf.to_csv(obs / "observables.csv", index=False)
    n = norm if norm is not None else [1.0, 0.999, 1.001, 1.0]
    pd.DataFrame({"step": [0, 2, 4, 6], "time_au": [0.0, 0.04, 0.08, 0.12],
                  "norm_check": n}).to_csv(obs / "wp_momentum_stats.csv", index=False)
    return tmp_path


def test_complete_run_passes(tmp_path):
    r = validate_run(_write_run(tmp_path))
    assert r.passed, r.summary()


def test_dropped_required_column_fails_tier1(tmp_path):
    r = validate_run(_write_run(tmp_path, drop_energy_col=True))
    assert not r.passed
    et = next(o for o in r.observables if o.name == "energy_total")
    assert any(t.tier == "existence" and not t.passed for t in et.tiers)


def test_nan_fails_tier3(tmp_path):
    r = validate_run(_write_run(tmp_path, nan=True))
    assert not r.passed
    et = next(o for o in r.observables if o.name == "energy_total")
    assert any(t.tier == "finite" and not t.passed for t in et.tiers)


def test_bad_norm_fails_tier4(tmp_path):
    r = validate_run(_write_run(tmp_path, norm=[1.0, 2.0, 1.0, 1.0]))
    assert not r.passed
    wp = next(o for o in r.observables if o.name == "wp_momentum_stats")
    assert any(t.tier == "invariant" and not t.passed for t in wp.tiers)


def test_energy_drift_fails_tier4(tmp_path):
    # drift 2 mHa > 1 mHa limit
    r = validate_run(_write_run(tmp_path, energy=[-0.6156, -0.6156, -0.6156, -0.6136]))
    assert not r.passed
    et = next(o for o in r.observables if o.name == "energy_total")
    assert any(t.tier == "invariant" and not t.passed for t in et.tiers)


def test_no_manifest_is_noted_not_passed(tmp_path):
    r = validate_run(_write_run(tmp_path, manifest=False))
    assert "no observables_manifest" in r.note
    assert r.observables == []
