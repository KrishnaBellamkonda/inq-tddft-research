"""Shared fixtures for the inqview test suite.

The whole suite is **pure-Python / numpy**, runs on any machine, needs no GPU
and no INQ engine (ADR 0005). Expected values are derived **analytically up
front**, never captured from current code output (anti-circularity, IV-M10).

Synthetic-signal helpers live in ``_signals.py`` (importable by the tests too).
"""
from __future__ import annotations

import os
import sys

# Tests now live in mirrored subfolders (tests/python/inqview/<layer>/...); make
# the shared `_signals` helper importable from any of them by putting this
# directory (tests/python/) on sys.path regardless of pytest's import mode.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pytest

from _signals import Tone, make_tone


@pytest.fixture
def tone() -> Tone:
    """A clean zero-mean on-bin tone (A=2.5, f0=0.15625 /a.u.)."""
    return make_tone()
