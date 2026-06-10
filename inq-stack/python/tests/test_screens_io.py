"""I/O test for the LEED-pattern loader (``inqview.io.leed.load_leed_pattern``).

Writes a tiny hand-built ``.dat`` (deterministic, so exact comparison is valid —
ADR 0005) and checks header parsing, the fftshift centring that fixes the
"four-corner-split" failure mode, and the physical extent.

The loader now lives at ``inqview.io.leed`` (ADR 0003 restructure — resolving
the loader-vs-phase name collision with ``postprocess.screens``); only the
import line moved.
"""
from __future__ import annotations

import numpy as np
import pytest

from inqview.io.leed import load_leed_pattern

pytestmark = pytest.mark.io

# A 4×4 screen whose single nonzero sample sits at FFT-natural index (0,0) =
# the physical origin. After load_leed_pattern's np.fft.fftshift, a length-4
# axis maps index 0 → index 2, so the peak must land at array centre [2,2].
_DAT = """\
# label=test_screen z=5.0 total_time=2.0 n_accum=4
# nx=4 ny=4 dx=1.0 dy=1.0 origin_x=0 origin_y=0
9 0 0 0
0 0 0 0
0 0 0 0
0 0 0 0
"""


@pytest.fixture
def dat_file(tmp_path):
    p = tmp_path / "screen_00.dat"
    p.write_text(_DAT)
    return p


def test_header_fields_parsed(dat_file):
    pat = load_leed_pattern(dat_file)
    assert pat.label == "test_screen"
    assert pat.z_bohr == pytest.approx(5.0)
    assert pat.total_time_au == pytest.approx(2.0)
    assert pat.n_accum == 4
    assert (pat.nx, pat.ny) == (4, 4)
    assert pat.dx_bohr == pytest.approx(1.0)
    assert pat.dy_bohr == pytest.approx(1.0)


def test_fftshift_moves_origin_peak_to_centre(dat_file):
    """The corner peak (FFT-natural origin) must be centred after loading."""
    pat = load_leed_pattern(dat_file)
    assert pat.data.shape == (4, 4)
    assert pat.data.sum() == pytest.approx(9.0)          # mass conserved
    assert pat.data[2, 2] == pytest.approx(9.0)          # centred, not in a corner
    # all four corners are now empty (the failure mode this guards against)
    for iy, ix in [(0, 0), (0, 3), (3, 0), (3, 3)]:
        assert pat.data[iy, ix] == pytest.approx(0.0)


def test_origin_and_extent_are_centred(dat_file):
    """Loader overrides origin to −L/2 so extent spans [−Lx/2, +Lx/2, …]."""
    pat = load_leed_pattern(dat_file)
    assert pat.origin_x_bohr == pytest.approx(-2.0)      # -0.5 * nx * dx
    assert pat.origin_y_bohr == pytest.approx(-2.0)
    assert pat.x_axis[0] == pytest.approx(-2.0)
    assert pat.extent_bohr == pytest.approx((-2.0, 2.0, -2.0, 2.0))


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_leed_pattern(tmp_path / "nope.dat")


if __name__ == "__main__":
    import subprocess
    import sys

    sys.exit(subprocess.call(["pytest", "-v", __file__]))
