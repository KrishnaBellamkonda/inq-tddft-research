"""From-run assembly test for WP-integrity (IV-M05 glue).

Builds a tiny synthetic run directory with the real CSV layout
(``results/raw/observables/{momentum_distribution,wp_real_space_stats}.csv``)
and asserts ``assemble_from_run`` reconstructs the WPIntegrity series with the
analytically known values defined up front (NOT captured from code output):

- ``sigma_r = sqrt(sigma_x2+sigma_y2+sigma_z2)`` per step,
- ``kl_mom`` via KL of the per-step ``n_wp`` distributions vs the initial one,
- ``ipr`` is NaN (no WP-only density frames are saved by the pipeline).

Pure-numpy/pandas, portable (ADR 0005); no real run data needed.
"""
from __future__ import annotations

import numpy as np
import pytest

from inqview.analysis.wp_integrity import assemble_from_run

pytestmark = pytest.mark.analysis

# WP momentum distributions (n_wp), 3 k-bins, 2 steps:
#   P0 = [1,2,1] -> normalised [0.25, 0.50, 0.25]
#   P1 = [2,2,0] -> normalised [0.50, 0.50, 0.00]
# KL(P1||P0) = 0.5 ln(0.5/0.25) + 0.5 ln(0.5/0.5) + 0  = 0.5 ln 2 = 0.3465735903
_KL_01 = 0.5 * np.log(2.0)

_MOM_CSV = """\
# l_bohr=10  n_bins=3  wp_idx=100
step,time_au,k_bohr_inv,n_total,n_wp
0,0.0,0.1,9,1.0
0,0.0,0.2,9,2.0
0,0.0,0.3,9,1.0
2,0.04,0.1,9,2.0
2,0.04,0.2,9,2.0
2,0.04,0.3,9,0.0
"""

# sigma_r(step0) = sqrt(0.5+0.5+0.5)=sqrt(1.5); sigma_r(step2)=sqrt(1+2+1)=2
_RS_CSV = """\
# wp_state_index=100  write_every=2
step,time_au,x_mean,y_mean,z_mean,x2_mean,y2_mean,z2_mean,sigma_x2,sigma_y2,sigma_z2,norm_check
0,0.0,0,0,0,0.5,0.5,0.5,0.5,0.5,0.5,1
2,0.04,0,0,0,1,2,1,1.0,2.0,1.0,1
"""


@pytest.fixture()
def run_dir(tmp_path):
    obs = tmp_path / "results" / "raw" / "observables"
    obs.mkdir(parents=True)
    (obs / "momentum_distribution.csv").write_text(_MOM_CSV)
    (obs / "wp_real_space_stats.csv").write_text(_RS_CSV)
    return tmp_path


def test_assemble_time_and_sigma(run_dir):
    wi = assemble_from_run(run_dir)
    assert np.allclose(wi.time_au, [0.0, 0.04])
    assert np.allclose(wi.sigma_r, [np.sqrt(1.5), 2.0])


def test_assemble_kl_known(run_dir):
    wi = assemble_from_run(run_dir)
    assert wi.kl_mom[0] == pytest.approx(0.0, abs=1e-12)   # KL(P0||P0)
    assert wi.kl_mom[1] == pytest.approx(_KL_01, abs=1e-9)


def test_assemble_ipr_is_nan(run_dir):
    wi = assemble_from_run(run_dir)
    assert wi.ipr.shape == (2,)
    assert np.isnan(wi.ipr).all()


def test_assemble_reference_previous_first_is_zero(run_dir):
    wi = assemble_from_run(run_dir, reference="previous")
    # frame-to-frame: first frame compares to itself -> 0; second is KL(P1||P0).
    assert wi.kl_mom[0] == pytest.approx(0.0, abs=1e-12)
    assert wi.kl_mom[1] == pytest.approx(_KL_01, abs=1e-9)
