// ============================================================================
// smoke_C3.cpp — Step C3 of the classical-electron rollout.
//
// Goal: confirm the impulsive propagator advances ion positions exactly as
//   r(t+dt) = r(t) + dt * v
// for our electron-mass species. Pure kinematic test — bypasses real_time::
// propagate and the electrons machinery (which requires num_electrons > 0
// per inq/src/systems/electrons.hpp:191).
//
// We call inq::ionic::propagator::impulsive::propagate_positions() directly,
// 20 times at dt = 0.005 a.u., and verify the final z = z0 + 20*dt*v_z to
// machine precision.
//
// Pass criteria (per docs/plans/the-objective-in-this-dapper-moon.md, C3):
//   - |Δz/N_steps - dt*v_z| < 1e-12  (impulsive update is exact)
//   - |x| = |y| = 0 to machine precision (no transverse drift)
// ============================================================================

#include <inq/inq.hpp>
#include <ionic/propagator.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

int main() {
    using namespace inq;
    using namespace inq::magnitude;

    auto env = input::environment{};

    std::cout << std::setprecision(16);
    std::cout << "=== smoke_C3: free-flight impulsive propagator ===\n\n";

    constexpr double m_e_amu = 1.0 / inq::ionic::species::amu_to_atomic_units;
    constexpr double v_z     = 8.5732;
    constexpr double dt      = 0.005;
    constexpr int    N_STEPS = 20;
    constexpr double z0      = 5.0;

    auto cell = systems::cell::cubic(10.0_b).finite();
    systems::ions ions(cell);
    auto sp = ionic::species("H")
        .pseudo_file("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
                     "shared/pseudopotentials/electron-ONCV-1.2.upf")
        .mass(m_e_amu);

    ions.insert(sp, {0.0_b, 0.0_b, z0 * 1.0_b});
    ions.velocities()[0] = vector3<double>{0.0, 0.0, v_z};

    inq::ionic::propagator::impulsive prop;
    std::vector<vector3<double>> dummy_forces(ions.size(), vector3<double>{0.0, 0.0, 0.0});

    std::ofstream csv("trajectory.csv");
    csv << std::setprecision(16);
    csv << "step,t_au,x,y,z\n";
    csv << 0 << "," << 0.0 << ","
        << ions.positions()[0][0] << ","
        << ions.positions()[0][1] << ","
        << ions.positions()[0][2] << "\n";

    for (int n = 1; n <= N_STEPS; ++n) {
        prop.propagate_positions(dt, ions, dummy_forces);
        csv << n << "," << (n * dt) << ","
            << ions.positions()[0][0] << ","
            << ions.positions()[0][1] << ","
            << ions.positions()[0][2] << "\n";
    }

    auto const& pos = ions.positions()[0];
    double z_final  = pos[2];
    double z_expect = z0 + N_STEPS * dt * v_z;
    double dz_err   = std::abs(z_final - z_expect);
    double xy_err   = std::abs(pos[0]) + std::abs(pos[1]);

    std::cout << "v_z [bohr/atu]        = " << v_z      << "\n";
    std::cout << "dt [au]               = " << dt       << "\n";
    std::cout << "N_STEPS               = " << N_STEPS  << "\n";
    std::cout << "z0                    = " << z0       << "\n";
    std::cout << "z_final               = " << z_final  << "\n";
    std::cout << "z_expect (= z0+N·dt·v)= " << z_expect << "\n";
    std::cout << "|z_err|               = " << dz_err   << "  (expect < 1e-12)\n";
    std::cout << "|x|+|y|               = " << xy_err   << "  (expect = 0)\n\n";

    bool z_ok  = dz_err < 1e-12;
    bool xy_ok = xy_err < 1e-15;

    std::cout << "[smoke_C3] z within tol           : " << std::boolalpha << z_ok  << "\n";
    std::cout << "[smoke_C3] xy within tol          : " << std::boolalpha << xy_ok << "\n";
    std::cout << "[smoke_C3] " << ((z_ok && xy_ok) ? "PASS" : "FAIL") << "\n";
    return (z_ok && xy_ok) ? 0 : 1;
}
