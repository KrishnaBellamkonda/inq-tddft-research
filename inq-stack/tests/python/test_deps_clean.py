"""Executable ADR-0003 invariant: importing an analysis kernel must pull in NO
matplotlib and NO VTK (so a headless cluster node can compute observables with
numpy/scipy only).

Checked in a clean subprocess (sys.modules in the test process is already
polluted by pytest/matplotlib). Now ENFORCED: the kernel lives at
``inqview.analysis.fourier`` and ``inqview/__init__`` is lazy (PEP 562), so
neither importing the package nor the kernel drags in plotting deps. We probe
both the whole ``inqview.analysis`` package and the kernel module, and also
assert that merely ``import inqview`` stays matplotlib/VTK-free.
"""
from __future__ import annotations

import subprocess
import sys
import textwrap

import pytest

pytestmark = pytest.mark.analysis

# Modules that MUST import cleanly (no matplotlib, no VTK), per ADR 0003.
_CLEAN_TARGETS = ("inqview", "inqview.analysis", "inqview.analysis.fourier")


def _probe(module: str) -> str:
    return textwrap.dedent(
        f"""
        import sys, importlib
        importlib.import_module({module!r})
        bad = [m for m in ('matplotlib', 'vtk') if m in sys.modules]
        if bad:
            print('deps-clean violated by {module}:', ','.join(bad))
            sys.exit(1)
        sys.exit(0)
        """
    )


@pytest.mark.parametrize("module", _CLEAN_TARGETS)
def test_analysis_import_is_matplotlib_and_vtk_free(module):
    result = subprocess.run(
        [sys.executable, "-c", _probe(module)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr


if __name__ == "__main__":
    import subprocess as sp

    sys.exit(sp.call(["pytest", "-v", __file__]))
