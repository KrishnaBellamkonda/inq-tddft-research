// Pure-tier test of inqkit::absorbers mask shapes (mask_shape.hpp): the
// single-sided sin^2 mask (Eq. 13) and the NEW two-sided symmetric mask used by
// the two-sided CAP-vs-mask study (docs/plans/twosided-cap-vs-mask.md).
//
// Reference: De Giovannini, Larsen & Rubio, arXiv:1409.1689 (2014), Eq. 13.
// The two-sided mask must be (a) symmetric in z, (b) =1 in the inner region,
// (c) =0 at/beyond both walls, (d) =1/2 at the per-end midpoint (sin²(π/4)=1/2),
// and (e) equal to the single-sided ramp evaluated at |z|.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inqkit/absorbers/mask_shape.hpp>

using inqkit::absorbers::sin2_mask_value;
using inqkit::absorbers::sin2_mask_value_twosided;
using Catch::Approx;

TEST_CASE("single-sided sin^2 mask: boundary and midpoint values", "[absorbers][mask][pure]") {
  const double z0 = 20.0, L = 10.0;       // ramp over [20, 30]
  CHECK(sin2_mask_value(0.0, z0, L)  == Approx(1.0));   // inner region
  CHECK(sin2_mask_value(20.0, z0, L) == Approx(1.0));   // ramp start
  CHECK(sin2_mask_value(30.0, z0, L) == Approx(0.0));   // wall
  CHECK(sin2_mask_value(35.0, z0, L) == Approx(0.0));   // beyond wall
  CHECK(sin2_mask_value(25.0, z0, L) == Approx(0.5));   // midpoint: 1-sin²(π/4)=0.5
  // monotonic non-increasing across the ramp
  CHECK(sin2_mask_value(22.0, z0, L) > sin2_mask_value(28.0, z0, L));
}

TEST_CASE("two-sided sin^2 mask: symmetry, inner=1, walls=0, midpoint=1/2", "[absorbers][mask][pure]") {
  const double z_in = 20.0, Lh = 5.0;     // ramp over |z| in [20, 25]
  // inner region
  CHECK(sin2_mask_value_twosided(0.0,  z_in, Lh) == Approx(1.0));
  CHECK(sin2_mask_value_twosided(19.0, z_in, Lh) == Approx(1.0));
  CHECK(sin2_mask_value_twosided(20.0, z_in, Lh) == Approx(1.0));   // ramp start
  // walls (both signs)
  CHECK(sin2_mask_value_twosided( 25.0, z_in, Lh) == Approx(0.0));
  CHECK(sin2_mask_value_twosided(-25.0, z_in, Lh) == Approx(0.0));
  CHECK(sin2_mask_value_twosided( 99.0, z_in, Lh) == Approx(0.0));
  // per-end midpoint |z| = 22.5 -> 0.5
  CHECK(sin2_mask_value_twosided( 22.5, z_in, Lh) == Approx(0.5));
  CHECK(sin2_mask_value_twosided(-22.5, z_in, Lh) == Approx(0.5));
}

TEST_CASE("two-sided mask is symmetric and matches single-sided at |z|", "[absorbers][mask][pure]") {
  const double z_in = 20.0, Lh = 5.0;
  for (double z = -30.0; z <= 30.0; z += 1.0) {
    // symmetry M(-z) == M(z)
    CHECK(sin2_mask_value_twosided(-z, z_in, Lh) ==
          Approx(sin2_mask_value_twosided(z, z_in, Lh)));
    // equals the single-sided ramp (start z_in, width Lh) evaluated at |z|
    CHECK(sin2_mask_value_twosided(z, z_in, Lh) ==
          Approx(sin2_mask_value(std::fabs(z), z_in, Lh)));
  }
}
