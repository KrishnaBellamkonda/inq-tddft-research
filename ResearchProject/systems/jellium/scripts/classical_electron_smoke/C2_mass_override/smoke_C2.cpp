// ============================================================================
// smoke_C2.cpp — Step C2 of the classical-electron rollout.
//
// Goal: confirm that `.mass(amu)` actually changes ions.species(0).mass() to
// 1.0 atomic units (i.e. m_e), and that the kinetic_energy() returned by
// systems::ions matches ½·m·v² at v = 8.5732 bohr/atu (the WP velocity).
//
// Pass criteria (per docs/plans/the-objective-in-this-dapper-moon.md, C2):
//   - mass_au ∈ [0.999, 1.001]  (electron mass in atomic units to ≥3 digits)
//   - KE_au  ≈ 36.75 Ha  (1 keV in Ha)
//   - KE_eV  ∈ [999, 1001]
// ============================================================================

#include <inq/inq.hpp>
#include <iostream>
#include <iomanip>

int main() {
    using namespace inq;
    using namespace inq::magnitude;

    auto env = input::environment{};

    std::cout << std::setprecision(12);
    std::cout << "=== smoke_C2: mass override = m_e ===\n\n";

    // 1.0 / 1822.8885 ≈ 5.485799082e-4 amu — exactly converts back to 1.0 a.u.
    constexpr double m_e_amu = 1.0 / inq::ionic::species::amu_to_atomic_units;
    constexpr double v_z     = 8.5732;   // bohr/atu, matches WP k_0 at 1000 eV
    constexpr double Ha_to_eV = 27.211386245988;

    auto cell = systems::cell::cubic(10.0_b).finite();
    systems::ions ions(cell);

    auto sp = ionic::species("H")
        .pseudo_file("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
                     "shared/pseudopotentials/electron-ONCV-1.2.upf")
        .mass(m_e_amu);

    ions.insert(sp, {0.0_b, 0.0_b, 0.0_b});
    ions.velocities()[0] = vector3<double>{0.0, 0.0, v_z};

    double mass_au   = ions.species(0).mass();
    double KE_au     = ions.kinetic_energy();
    double KE_expect = 0.5 * 1.0 * v_z * v_z;
    double KE_eV     = KE_au * Ha_to_eV;

    std::cout << "m_e_amu (literal)     = " << m_e_amu     << "  (= 1.0 / 1822.8885)\n";
    std::cout << "mass_au               = " << mass_au     << "  (expect ≈ 1.0)\n";
    std::cout << "v_z [bohr/atu]        = " << v_z         << "\n";
    std::cout << "KE_au   (½·m·v²)      = " << KE_au       << "\n";
    std::cout << "KE_au expected (m=1)  = " << KE_expect   << "  (= 36.75)\n";
    std::cout << "KE_eV                 = " << KE_eV       << "  (expect ≈ 1000)\n\n";

    bool mass_ok = std::abs(mass_au - 1.0) < 1e-3;
    bool ke_ok   = std::abs(KE_eV - 1000.0) < 1.0;

    std::cout << "[smoke_C2] mass_au within tol     : " << std::boolalpha << mass_ok << "\n";
    std::cout << "[smoke_C2] KE_eV within tol       : " << std::boolalpha << ke_ok   << "\n";
    std::cout << "[smoke_C2] " << ((mass_ok && ke_ok) ? "PASS" : "FAIL") << "\n";
    return (mass_ok && ke_ok) ? 0 : 1;
}
