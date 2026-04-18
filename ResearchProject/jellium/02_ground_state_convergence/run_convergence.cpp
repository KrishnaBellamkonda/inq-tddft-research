// 02_ground_state_convergence/run_convergence.cpp
//
// WHAT THIS DOES
//   Runs two systematic convergence tests for the jellium system
//   (N electrons, cubic cell of side L, r_s fixed = 7.26 bohr).
//
// TEST A — Grid-spacing (E_cut) convergence
//   Fix N=40, L=40 bohr.  Vary the real-space grid spacing h from 0.70 to
//   0.30 bohr.  The kinetic energy cutoff E_cut = π²/(2h²) spans ~10–55 Ha.
//   For plane-wave jellium, convergence is rapid (exponential in E_cut).
//   Goal: identify the coarsest h that keeps ΔE_total < 0.01 Ha.
//
// TEST B — Shell-closure (finite-N) convergence
//   Use the optimal h from Test A.  Run at "magic" electron counts N where
//   the Gamma-point shell structure is exactly closed:
//     N = 2  (|n|²=0 full), 14 (|n|²≤1 full), 38 (|n|²≤2 full),
//     54 (|n|²≤3 full), 66 (|n|²≤4 full).
//   For each N, L is set so that r_s = 7.26 bohr is preserved.
//   T_s/N should oscillate around the Thomas-Fermi value (3/5)E_F and converge
//   to it from above/below as N → ∞.
//
// OUTPUT
//   Writes results/convergence_results.csv with two CSV blocks:
//     Block A header: # TEST_A spacing_bohr E_cut_Ha E_total_Ha T_s_Ha E_xc_Ha n_iter
//     Block B header: # TEST_B N L_bohr k0_inv_bohr Ts_Ha Ts_per_N T_TF_per_N n_iter
//   Progress messages go to stderr.  Plot with plot_convergence.py.
//
// Build and run:
//   inq-run

#include <inq/inq.hpp>
#include "../01_ground_state/jellium_utils.hpp"

#include <iostream>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <cmath>
#include <vector>
#include <tuple>

int main() {
    using namespace inq;
    using namespace inq::magnitude;
    namespace fs = std::filesystem;

    // ── Physical constants ───────────────────────────────────────────────────

    constexpr double PI       = M_PI;
    constexpr double SMEAR_EV = 0.00862;  // 100 K Fermi smearing

    // ── System parameters (shared r_s for both tests) ────────────────────────

    // r_s is set by the N=40, L=40 bohr reference system.
    // All shell-closure runs use the same r_s = 7.26 bohr.
    constexpr double RS_BOHR = 7.2557;   // Wigner-Seitz radius (bohr)

    // Helper: L for N electrons at the reference r_s.
    auto L_from_N = [&](int N) -> double {
        return RS_BOHR * std::cbrt(4.0 * PI * N / 3.0);
    };

    // ── Output file ──────────────────────────────────────────────────────────

    fs::create_directories("results");
    std::ofstream out("results/convergence_results.csv");
    out << std::fixed << std::setprecision(8);

    // ── Output header ─────────────────────────────────────────────────────────

    out << "# Jellium convergence tests\n";
    out << "# r_s = " << RS_BOHR << " bohr  (N=40, L=40 bohr reference)\n";
    out << "# Smearing = " << SMEAR_EV << " eV  (100 K Fermi-Dirac)\n";
    out << "#\n";

    // ════════════════════════════════════════════════════════════════════════
    // TEST A — Grid-spacing convergence
    // ════════════════════════════════════════════════════════════════════════

    constexpr int    N_A = 40;
    constexpr double L_A = 40.0;   // bohr

    // Different grid spacings would mean having different
    // E_cut values for jellium
    std::vector<double> spacings = {0.70, 0.60, 0.50, 0.45, 0.40, 0.35, 0.30};

    out << "# TEST_A spacing_bohr,E_cut_Ha,E_total_Ha,T_s_Ha,E_xc_Ha,n_iter\n";

    for (double h : spacings) {
        double E_cut = PI * PI / (2.0 * h * h);   // Ha

        systems::ions ions(systems::cell::cubic(L_A * 1.0_b).periodic());

        systems::electrons electrons(ions,
            options::electrons{}
                .spacing(h * 1.0_b)
                .extra_electrons(N_A)
                .extra_states(8)
                .temperature(SMEAR_EV * 1.0_eV),
            input::kpoints::gamma()
        );

        ground_state::initial_guess(ions, electrons);

        auto gs = ground_state::calculate(ions, electrons,
            options::theory{}.lda(),
            options::ground_state{}.energy_tolerance(1e-6_Ha)  // slightly looser for sweep
        );

        out << "# TEST_A "
            << h                    << ","
            << E_cut                << ","
            << gs.energy.total()    << ","
            << gs.energy.kinetic()  << ","
            << gs.energy.xc()       << ","
            << gs.total_iter        << "\n";
        out.flush();

        std::cerr << "[Test A] h=" << h << " bohr  E_cut=" << E_cut
                  << " Ha  E_total=" << gs.energy.total() << " Ha\n";
    }

    out << "#\n";

    // ════════════════════════════════════════════════════════════════════════
    // TEST B — Shell-closure convergence at fixed r_s
    // ════════════════════════════════════════════════════════════════════════
    //
    // Magic closed-shell electron counts for the cubic Gamma-point:
    //   N = 2  : |n|²=0 full (1 state,  2 electrons)
    //   N = 14 : |n|²≤1 full (7 states, 14 electrons)
    //   N = 38 : |n|²≤2 full (19 states, 38 electrons)
    //   N = 54 : |n|²≤3 full (27 states, 54 electrons)
    //   N = 66 : |n|²≤4 full (33 states, 66 electrons)
    //
    // For each N, L is chosen to preserve r_s = 7.26 bohr.
    // Use h = 0.50 bohr (E_cut ≈ 20 Ha) — converged for smooth jellium.

    constexpr double H_B = 0.50;   // grid spacing for shell-closure test (bohr)

    // Thomas-Fermi kinetic energy density: T_TF/N = (3/5) E_F
    double kF_ref = fermi_wavevector(N_A, L_A);
    double EF_ref = fermi_energy(N_A, L_A);
    double T_TF_per_N = (3.0 / 5.0) * EF_ref;   // Ha per electron (bulk limit)

    // Closed-shell systems to test
    //
    struct ClosedSystem { int N; double L; int extra_states; };
    // Contains a list of different number of electrons
    // to use
    std::vector<ClosedSystem> systems_B = {
        {2,  L_from_N(2),  12},   // next shell has degeneracy 6 → 12 states
        {14, L_from_N(14),  8},
        {38, L_from_N(38),  8},
        {54, L_from_N(54),  8},
        {66, L_from_N(66),  8},
    };

    out << "# TEST_B N,L_bohr,k0_inv_bohr,Ts_Ha,Ts_per_N,T_TF_per_N,n_iter\n";
    out << "# T_TF_per_N = " << T_TF_per_N << " Ha  (Thomas-Fermi bulk limit)\n";
	
    // For differnt number of electrons (for the shell test)
    // obtaining the ground state energy 
    for (auto const & sys : systems_B) {
        int    N = sys.N;
        double L = sys.L;
        double k0 = 2.0 * PI / L;

        systems::ions ions(systems::cell::cubic(L * 1.0_b).periodic());

        systems::electrons electrons(ions,
            options::electrons{}
                .spacing(H_B * 1.0_b)
                .extra_electrons(N)
                .extra_states(sys.extra_states)
                .temperature(SMEAR_EV * 1.0_eV),
            input::kpoints::gamma()
        );

        ground_state::initial_guess(ions, electrons);

        auto gs = ground_state::calculate(ions, electrons,
            options::theory{}.lda(),
            options::ground_state{}.energy_tolerance(1e-6_Ha)
        );

        double Ts      = gs.energy.kinetic();
        double Ts_per_N = Ts / N; // Measured average K.E per electron

        out << "# TEST_B "
            << N            << ","
            << L            << ","
            << k0           << ","
            << Ts           << ","
            << Ts_per_N     << ","
            << T_TF_per_N   << ","
            << gs.total_iter << "\n";
        out.flush();

        std::cerr << "[Test B] N=" << N << "  L=" << L << "  Ts/N=" << Ts_per_N
                  << "  T_TF/N=" << T_TF_per_N << "\n";
    }

    out << "#\n# Done.\n";
    return 0;
}
