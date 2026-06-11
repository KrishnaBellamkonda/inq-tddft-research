// Pure-tier test of inqkit::jellium::shells: the free-electron-gas shell
// enumeration in a cubic box. Asserts the documented (|G|^2, degeneracy) table
// and the cumulative orbital count at the 162-electron closed shell (|G|^2=6),
// plus partial-shell truncation.
//
// Reference: header docstring + docs/sources/free-electron-gas-magic-numbers.md.
// Orbital fill to shells 0..6 = 1+6+12+8+6+24+24 = 81 orbitals = 162 electrons.

#include <catch2/catch_test_macros.hpp>

#include <inqkit/jellium/shells.hpp>

#include <vector>

using namespace inqkit::jellium::shells;

TEST_CASE("jellium shells: documented shell table up to the 162-electron closure", "[jellium][shells][pure]") {
  // 81 orbitals == shells |G|^2 = 0..6 fully filled.
  auto shells = enumerate_for_n_states(81);

  const std::vector<int> gsq_expected  = {0, 1, 2, 3, 4, 5, 6};
  const std::vector<int> deg_expected  = {1, 6, 12, 8, 6, 24, 24};

  REQUIRE(shells.size() == gsq_expected.size());

  int cumulative = 0;
  for (std::size_t k = 0; k < shells.size(); ++k) {
    CHECK(shells[k].shell_id == static_cast<int>(k));
    CHECK(shells[k].gsq == gsq_expected[k]);
    CHECK(shells[k].degeneracy == deg_expected[k]);
    CHECK(shells[k].first_index == cumulative);
    CHECK(shells[k].members.size() == static_cast<std::size_t>(deg_expected[k]));
    if (!shells[k].members.empty()) {
      CHECK(shells[k].members.front() == cumulative);
      CHECK(shells[k].members.back() == cumulative + deg_expected[k] - 1);
    }
    cumulative += shells[k].degeneracy;
  }
  CHECK(cumulative == 81);  // 81 orbitals -> 162 electrons (closed shell)
}

TEST_CASE("jellium shells: partial last shell is truncated to n_states", "[jellium][shells][pure]") {
  // 5 orbitals: shell 0 (deg 1) full, shell 1 (deg 6) clipped to 4.
  auto shells = enumerate_for_n_states(5);
  REQUIRE(shells.size() == 2);
  CHECK(shells[0].degeneracy == 1);
  CHECK(shells[1].gsq == 1);
  CHECK(shells[1].degeneracy == 4);                 // clipped from 6
  CHECK(shells[1].members.size() == 4u);
  int total = 0;
  for (auto const &s : shells) total += s.degeneracy;
  CHECK(total == 5);
}
