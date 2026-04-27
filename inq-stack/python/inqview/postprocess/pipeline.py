"""Postprocess pipeline dispatcher.

Each phase is a callable ``run(results_dir, run_name, **opts) -> dict``. The
dispatcher logs a one-line ``[ok] phase`` / ``[skip] phase reason`` /
``[fail] phase msg`` per phase and aggregates results into a PipelineResult.
"""

from __future__ import annotations

import logging
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable

from . import (
    density,
    ground_state,
    layout,
    observables,
    orbitals,
    overlap,
    run_summary,
    screens,
)

PathLike = str | Path

PHASES: tuple[str, ...] = (
    "summary",
    "gs",
    "layout",
    "observables",
    "density",
    "screens",
    "overlap",
    "orbitals",
    "paraview",
    "paraview_3d",
)

# Phase -> module entry point
_PHASE_FUNCS: dict[str, Callable] = {
    "summary":     run_summary.run,
    "gs":          ground_state.run,
    "layout":      layout.run,
    "observables": observables.run,
    "density":     density.run,
    "screens":     screens.run,
    "overlap":     overlap.run,
    "orbitals":    orbitals.run,
}


@dataclass
class PipelineResult:
    """Aggregate result from a pipeline run."""

    results_dir: Path
    run_name: str
    ok: dict[str, dict] = field(default_factory=dict)
    skipped: dict[str, str] = field(default_factory=dict)
    failed: dict[str, str] = field(default_factory=dict)
    durations_s: dict[str, float] = field(default_factory=dict)

    def succeeded(self) -> bool:
        return not self.failed


def _norm_phases(phases: Iterable[str] | None) -> list[str]:
    if phases is None:
        return list(PHASES)
    out = []
    for p in phases:
        if p not in PHASES:
            raise ValueError(f"unknown phase {p!r}; valid: {PHASES}")
        out.append(p)
    return out


def run(
    results_dir: PathLike,
    *,
    run_name: str | None = None,
    phases: Iterable[str] | None = None,
    rebuild: bool = False,
    skip_paraview: bool = True,
    paraview_runner: Callable | None = None,
    logger: logging.Logger | None = None,
    **opts,
) -> PipelineResult:
    """Run the postprocess pipeline against ``results_dir``.

    Parameters
    ----------
    results_dir : path to a run's ``results/`` tree.
    run_name : run name shown in plot titles. Inferred from the parent
        directory name if not given.
    phases : iterable of phase names. ``None`` runs every phase.
    rebuild : if True, regenerate artefacts even when output files exist.
    skip_paraview : if True (default), the ``paraview`` phase is skipped
        because pvbatch is heavy. Set False once a renderer is plumbed in.
    paraview_runner : callable plugged in when ``skip_paraview=False``;
        signature ``(results_dir, run_name, rebuild) -> dict``.
    logger : optional logger; defaults to a stderr logger named
        ``inqview.postprocess``.
    **opts : forwarded verbatim to every phase function (so phases can pick
        up things like ``cmap``, ``percentile`` overrides).
    """
    log = logger or _default_logger()
    results_dir = Path(results_dir).resolve()
    if run_name is None:
        run_name = results_dir.parent.name or "run"
    requested = _norm_phases(phases)

    res = PipelineResult(results_dir=results_dir, run_name=run_name)
    for phase in requested:
        t0 = time.monotonic()
        try:
            if phase == "paraview":
                if skip_paraview or paraview_runner is None:
                    res.skipped[phase] = (
                        "paraview skipped (skip_paraview=True or no runner)"
                    )
                    log.info("[skip] %s — %s", phase, res.skipped[phase])
                    continue
                out = paraview_runner(results_dir, run_name, rebuild)
            elif phase == "paraview_3d":
                if skip_paraview:
                    res.skipped[phase] = (
                        "paraview_3d skipped (skip_paraview=True; use --with-paraview)"
                    )
                    log.info("[skip] %s — %s", phase, res.skipped[phase])
                    continue
                from . import paraview_3d
                out = paraview_3d.run(
                    results_dir, run_name=run_name, rebuild=rebuild, **opts)
            else:
                fn = _PHASE_FUNCS[phase]
                out = fn(
                    results_dir,
                    run_name=run_name,
                    rebuild=rebuild,
                    **opts,
                )
            if out is None:
                out = {}
            res.ok[phase] = out
            log.info("[ok] %s (%.1fs)", phase, time.monotonic() - t0)
        except _SkipPhase as exc:
            res.skipped[phase] = str(exc)
            log.info("[skip] %s — %s", phase, exc)
        except Exception as exc:
            res.failed[phase] = f"{type(exc).__name__}: {exc}"
            log.error(
                "[fail] %s — %s\n%s", phase, exc, traceback.format_exc()
            )
        finally:
            res.durations_s[phase] = time.monotonic() - t0
    return res


class _SkipPhase(Exception):
    """Raised inside a phase to declare a clean skip rather than a failure."""


def skip(reason: str):
    raise _SkipPhase(reason)


def _default_logger() -> logging.Logger:
    log = logging.getLogger("inqview.postprocess")
    if not log.handlers:
        log.setLevel(logging.INFO)
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter("%(message)s"))
        log.addHandler(h)
    return log
