// ============================================================================
// syntax_check.cpp — verifies that the new inqkit headers
//   inqkit/jellium/shells.hpp
//   inqkit/observables/orbital_overlap.hpp::snapshot_proxies
// compile and link cleanly without running anything.
//
// We construct a small finite-cell setup, an OrbitalOverlapMatrix, and call
// snapshot_proxies on a trivial proxy list. Cell is 5 Bohr cubic with H
// pseudopotential (a real INQ pseudopotential, just for instantiation), so
// no custom UPF games. We DO call ground_state::initial_guess to populate
// the orbital_set so snapshot_proxies has something to extract; but we do
// NOT run any SCF iterations.
// ============================================================================

#include <inq/inq.hpp>

#include <inqkit/jellium/shells.hpp>
#include <inqkit/observables/orbital_overlap.hpp>

#include <iostream>

using namespace inq;
using namespace inq::magnitude;

int main() {
    std::cout << "=== syntax_check_proxies ===\n";

    // 1. Shell enumeration.
    auto shells = inqkit::jellium::shells::enumerate_for_n_states(101);
    std::cout << "  enumerate_for_n_states(101) -> " << shells.size() << " shells\n";
    for (auto const &s : shells)
        std::cout << "    shell_id=" << s.shell_id
                  << " gsq=" << s.gsq
                  << " degeneracy=" << s.degeneracy
                  << " first_index=" << s.first_index << "\n";

    auto proxies = inqkit::jellium::shells::pick_proxies(shells, 2);
    std::cout << "  pick_proxies -> " << proxies.size() << " indices: ";
    for (int p : proxies) std::cout << p << " ";
    std::cout << "\n";

    // Test write_shells_csv.
    inqkit::jellium::shells::write_shells_csv(shells, proxies, "shells_test_out");
    std::cout << "  wrote shells_test_out/shells.csv\n";

    // 2. Tiny INQ system to verify snapshot_proxies compiles + links.
    auto cell = systems::cell::cubic(10.0_b).finite();
    systems::ions ions(cell);
    ions.insert("H", {0.0_b, 0.0_b, 0.0_b});

    auto electrons = systems::electrons(
        ions,
        options::electrons{}.cutoff(15.0_Ha).extra_states(2).temperature(300.0_K));
    ground_state::initial_guess(ions, electrons);

    inqkit::observables::OrbitalOverlapMatrix obs(electrons, /*n_ref=*/2,
                                                  "overlap_test_out");
    std::vector<int> small_proxies = {0, 1};
    obs.snapshot_proxies(electrons, small_proxies, /*time=*/0.0, /*step=*/0);
    std::cout << "  snapshot_proxies completed\n";

    std::cout << "[syntax_check_proxies] PASS\n";
    return 0;
}
