// Pure-tier unit test for inqkit::detail::Vec3 (supports T07/T09).

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inqkit/detail/vec3.hpp>

using inqkit::detail::Vec3;
using Catch::Approx;

TEST_CASE("Vec3: default is the zero vector", "[vec3][pure]") {
  Vec3 v;
  CHECK(v.x == 0.0);
  CHECK(v.y == 0.0);
  CHECK(v.z == 0.0);
}

TEST_CASE("Vec3: operator[] reads and writes x/y/z (D1)", "[vec3][pure]") {
  Vec3 v{1.0, 2.0, 3.0};
  CHECK(v[0] == 1.0);
  CHECK(v[1] == 2.0);
  CHECK(v[2] == 3.0);
  v[0] = 10.0; v[1] = 20.0; v[2] = 30.0;   // non-const write
  CHECK(v.x == 10.0);
  CHECK(v.y == 20.0);
  CHECK(v.z == 30.0);
  Vec3 const c{4.0, 5.0, 6.0};             // const read
  CHECK(c[2] == 6.0);
}

TEST_CASE("Vec3: dot, norm2, norm", "[vec3][pure]") {
  Vec3 a{1.0, 2.0, 2.0};
  CHECK(a.norm2() == Approx(9.0));   // 1 + 4 + 4
  CHECK(a.norm() == Approx(3.0));

  Vec3 b{3.0, 0.0, -1.0};
  CHECK(a.dot(b) == Approx(1.0));    // 3 + 0 - 2

  // Orthogonal vectors have zero dot product.
  Vec3 ex{1.0, 0.0, 0.0}, ey{0.0, 1.0, 0.0};
  CHECK(ex.dot(ey) == Approx(0.0));
}

TEST_CASE("Vec3: arithmetic operators", "[vec3][pure]") {
  Vec3 a{1.0, 2.0, 3.0}, b{4.0, 5.0, 6.0};

  Vec3 sum = a + b;
  CHECK((sum == Vec3{5.0, 7.0, 9.0}));

  Vec3 diff = b - a;
  CHECK((diff == Vec3{3.0, 3.0, 3.0}));

  CHECK(((2.0 * a) == Vec3{2.0, 4.0, 6.0}));
  CHECK(((a * 2.0) == Vec3{2.0, 4.0, 6.0}));

  Vec3 c = a;
  c += b;
  CHECK((c == Vec3{5.0, 7.0, 9.0}));

  CHECK(a != b);
}

TEST_CASE("Vec3: project a displacement onto a unit direction (k-hat use case)", "[vec3][pure]") {
  // The COD docstring mentions projecting a centroid onto k-hat. Sanity-check
  // that dot with a unit vector recovers the along-axis component.
  Vec3 disp{3.0, 4.0, 0.0};
  Vec3 khat{0.0, 1.0, 0.0};
  CHECK(disp.dot(khat) == Approx(4.0));
}
