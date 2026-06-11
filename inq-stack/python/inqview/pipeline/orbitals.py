"""Phase: ``orbitals`` — RT orbital plots.

Currently a soft no-op placeholder. The C++ run template does not export
RT orbital VTIs by default (only GS orbital VTIs). When that changes, this
phase will produce ``results/analysis/orbitals/`` artefacts including the
WP orbital density GIF and a HOMO/LUMO panel.

For now this phase reports skipped with a clear reason.
"""

from __future__ import annotations

from pathlib import Path

from . import runner as _pipeline


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    rt_orb_dir = results_dir / "raw" / "vti" / "orbitals"
    if not rt_orb_dir.exists() or not any(rt_orb_dir.glob("*.vti")):
        _pipeline.skip(
            "no RT orbital VTIs (raw/vti/orbitals/ is empty); "
            "GS orbitals are handled by the gs phase."
        )
    return {"out_dir": str(results_dir / "analysis" / "orbitals"),
            "note": "RT orbitals present but no plotter implemented yet"}
