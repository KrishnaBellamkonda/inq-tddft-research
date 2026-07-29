// ============================================================================
// smoke_C1.cpp — Step C1 of the classical-electron rollout.
//
// Goal: verify that INQ can parse the new electron-ONCV-1.2.upf and that
// `ionic::species("H").pseudo_file(path)` returns a usable species.
//
// This is the cheapest possible test: a tiny finite cell, a single ion, no
// SCF, no propagation. Just construct the species + ions container and print
// what INQ thinks it is.
//
// Pass criteria (per docs/plans/the-objective-in-this-dapper-moon.md, Step C1):
//   - INQ does not throw or print a parser error on the UPF.
//   - `ions << std::cout` prints the ion at (0,0,0) with symbol "H".
//   - `ions.species(0).mass()` prints the *default* (proton) mass — we have
//     not called `.mass(...)` here. Default mass for "H" via INQ is the
//     periodic-table value × 1822.8885 ≈ 1837 atomic units.
//   - File path printed via has_file()/file_path() matches the UPF we copied.
// ============================================================================

#include <inq/inq.hpp>
#include <iostream>
#include <iomanip>

int main() {
    using namespace inq;
    using namespace inq::magnitude;

    auto env = input::environment{};

    std::cout << std::setprecision(12);
    std::cout << "=== smoke_C1: load electron-ONCV-1.2.upf ===\n\n";

    auto cell = systems::cell::cubic(10.0_b).finite();   // tiny finite box
    systems::ions ions(cell);

    auto sp = ionic::species("H").pseudo_file(
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "shared/pseudopotentials/electron-ONCV-1.2.upf");

    ions.insert(sp, {0.0_b, 0.0_b, 0.0_b});

    std::cout << "n_ions          = " << ions.size() << "\n";
    std::cout << "species(0) symbol      = " << ions.species(0).symbol() << "\n";
    std::cout << "species(0) atomic_num  = " << ions.species(0).atomic_number() << "\n";
    std::cout << "species(0) mass [a.u.] = " << ions.species(0).mass() << "  (expected default: ~1837 = proton mass)\n";
    std::cout << "species(0) has_file    = " << std::boolalpha << ions.species(0).has_file() << "\n";
    if (ions.species(0).has_file()) {
        std::cout << "species(0) file_path   = " << ions.species(0).file_path() << "\n";
    }

    std::cout << "\n--- ions container ---\n" << ions << "\n";
    std::cout << "[smoke_C1] DONE\n";
    return 0;
}
