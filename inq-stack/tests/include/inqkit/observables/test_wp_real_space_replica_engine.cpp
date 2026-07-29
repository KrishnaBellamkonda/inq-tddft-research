// Engine-tier physics validation for wp_real_space_stats' moments: compute the
// real-space ⟨r⟩ and per-axis variance directly from the injected WP field and
// check ⟨r⟩ ≈ the injected centre and Var ≈ σ²/2 (|ψ|² ∝ exp(-r²/σ²) for a WP
// envelope exp(-r²/2σ²)). Independent in-test replica of the GPU reduction the
// class performs; no source change.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/detail/grid_layout.hpp>

#include <complex>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;
using inqkit::detail::grid_layout::flatten_index;

TEST_CASE("wp_real_space moments: ⟨r⟩ = injected centre, Var = σ²/2", "[observables][wp_real_space][engine]") {
  const double CX = 1.0, CY = -1.0, CZ = 0.5, SIGMA = 1.5;

  systems::ions ions(systems::cell::cubic(16.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(2));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  auto report = inqkit::WavePacket{}
                    .center(CX, CY, CZ)
                    .sigma(SIGMA)
                    .k0(0.0, 0.0, 0.0)
                    .inject_into_last_extra_state(electrons, 1.0);

  auto psi = inqkit::fields::orbital::wavefunction(electrons, report.state_index);

  // NODE convention (matches INQ rvector = symmetric_coord·dx): physical
  // coordinate of inqkit field index ix is origin + ix·dx, NOT (ix+0.5)·dx.
  // (center_of_density.hpp uses the +0.5 form → +dx/2 offset; see E04.)
  long double norm = 0.0L, mx = 0.0L, my = 0.0L, mz = 0.0L;
  long double mx2 = 0.0L, my2 = 0.0L, mz2 = 0.0L;
  for (int ix = 0; ix < psi.nx; ++ix) {
    const double x = psi.origin_x_bohr + ix * psi.dx_bohr;
    for (int iy = 0; iy < psi.ny; ++iy) {
      const double y = psi.origin_y_bohr + iy * psi.dy_bohr;
      for (int iz = 0; iz < psi.nz; ++iz) {
        const double z = psi.origin_z_bohr + iz * psi.dz_bohr;
        const double w = std::norm(psi.values[flatten_index(ix, iy, iz, psi.ny, psi.nz)]);
        norm += w;
        mx += w * x; my += w * y; mz += w * z;
        mx2 += w * x * x; my2 += w * y * y; mz2 += w * z * z;
      }
    }
  }
  const double cx = static_cast<double>(mx / norm);
  const double cy = static_cast<double>(my / norm);
  const double cz = static_cast<double>(mz / norm);
  const double vx = static_cast<double>(mx2 / norm) - cx * cx;

  CHECK(cx == Approx(CX).margin(0.05));
  CHECK(cy == Approx(CY).margin(0.05));
  CHECK(cz == Approx(CZ).margin(0.05));
  // |ψ|² is a Gaussian of variance σ²/2 per axis.
  CHECK(vx == Approx(SIGMA * SIGMA / 2.0).margin(0.15));
}
