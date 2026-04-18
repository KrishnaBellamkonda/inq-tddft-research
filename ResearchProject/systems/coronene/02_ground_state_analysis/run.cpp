// ============================================================================
// 02_ground_state_analysis/run.cpp
//
// Purpose: Compute the LDA ground state of coronene (C24H12) in an isolated
//          finite cell.  Produces:
//   results/gs_summary.txt  — total energy, forces, orbital eigenvalues
//
// This run serves two purposes:
//   1. Verify the molecular geometry (forces near zero = correct equilibrium)
//   2. Establish the ground-state KS orbitals that will be used as the target
//      in the TDDFT WP scattering simulation.
//
// Run with: inq-run (GPU) or inq-run --cpu (CPU)
// ============================================================================

#include "config.hpp"
#include <fstream>
#include <iomanip>

int main(){
    using namespace inq;
    using namespace inq::magnitude;

    // ── 1. Ions: load from XYZ ───────────────────────────────────────────────
    // parse::xyz reads positions in Angstrom by default.
    // Cell is provided explicitly (finite, isolated).
    auto cell = cfg::make_cell();
    auto ions = systems::ions::parse(cfg::CORONENE_XYZ, cell);

    // ── 2. Electrons ─────────────────────────────────────────────────────────
    // extra_states(3): 1 for WP orbital + 2 SCF buffer
    auto electrons = systems::electrons(ions,
        options::electrons{}.cutoff(cfg::ECUT_HA * 1.0_Ha)
                            .extra_states(cfg::EXTRA_STATES));

    // ── 3. Ground state SCF ──────────────────────────────────────────────────
    ground_state::initial_guess(ions, electrons);

    auto result = ground_state::calculate(ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(cfg::SCF_TOL * 1.0_Ha)
            .max_steps(300)
            .calculate_forces());

    // ── 4. Output ─────────────────────────────────────────────────────────────
    if(electrons.root()){
        std::cout << "\n=== Coronene C24H12 ground state ===\n";
        std::cout << "  Atoms loaded : " << ions.size() << "  (expect 36)\n";
        std::cout << "  Cell (bohr)  : " << cfg::LX_BOHR << " x "
                  << cfg::LY_BOHR << " x " << cfg::LZ_BOHR << "\n";
        std::cout << "  E_cut        : " << cfg::ECUT_HA  << " Ha\n";
        std::cout << "  Extra states : " << cfg::EXTRA_STATES << "\n\n";

        // ── Energy breakdown ─────────────────────────────────────────────────
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "=== Energy breakdown (Ha) ===\n";
        std::cout << "  Total        : " << result.energy.total()    << "\n";
        std::cout << "  Kinetic      : " << result.energy.kinetic()  << "\n";
        std::cout << "  Hartree      : " << result.energy.hartree()  << "\n";
        std::cout << "  XC (LDA)     : " << result.energy.xc()      << "\n";
        std::cout << "  External     : " << result.energy.external() << "\n";
        std::cout << "  Non-local PP : " << result.energy.non_local()<< "\n";
        std::cout << "  Ion-ion      : " << result.energy.ion()      << "\n";
        std::cout << "  SCF steps    : " << result.total_iter        << "\n";
        std::cout << "  Total (eV)   : " << result.energy.total() * cfg::HA_TO_EV << "\n\n";

        // ── Forces ───────────────────────────────────────────────────────────
        std::cout << "=== Forces on atoms (Ha/bohr) — expect < 1e-2 for idealized geometry ===\n";
        double fmax = 0.0;
        for(int i = 0; i < ions.size(); i++){
            auto f = result.forces[i];
            double fn = std::sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
            fmax = std::max(fmax, fn);
            if(i < 10 || fn > 0.01)
                std::cout << "  atom[" << std::setw(2) << i << "] "
                          << (i < 24 ? "C" : "H")
                          << "  |F|=" << std::setw(8) << fn
                          << "  (" << f[0] << ", " << f[1] << ", " << f[2] << ")\n";
        }
        std::cout << "  Max |force|: " << fmax << " Ha/bohr\n\n";

        // ── Orbital occupations and eigenvalues ──────────────────────────────
        std::cout << "=== KS orbital summary ===\n";
        std::cout << "  Total states : " << electrons.kpin()[0].spinor_set_size() << "\n";
        auto const & occ = electrons.occupations()[0];
        for(int ist = 0; ist < electrons.kpin()[0].set_part().local_size(); ist++){
            auto istg = electrons.kpin()[0].set_part().local_to_global(ist).value();
            std::cout << "  State " << std::setw(3) << istg
                      << "  occ=" << std::setw(5) << occ[ist] << "\n";
        }

        // ── Validation: total electrons from density ─────────────────────────
        // (Visual sanity check — should sum to 108 for coronene ground state)
        std::cout << "\n=== Validation ===\n";
        double total_occ = 0.0;
        for(int ist = 0; ist < (int)occ.size(); ist++) total_occ += occ[ist];
        std::cout << "  Sum of occupations : " << total_occ
                  << "  (expect 108.0 for coronene GS)\n";
        std::cout << "  SCF converged?     : " << (result.total_iter < 300 ? "YES" : "CHECK") << "\n";

        // ── Save summary file ─────────────────────────────────────────────────
        std::ofstream out("results/gs_summary.txt");
        out << std::fixed << std::setprecision(8);
        out << "# Coronene C24H12 — LDA ground state\n";
        out << "# E_cut=" << cfg::ECUT_HA << " Ha\n";
        out << "E_total_Ha   " << result.energy.total()   << "\n";
        out << "E_kinetic_Ha " << result.energy.kinetic() << "\n";
        out << "E_xc_Ha      " << result.energy.xc()     << "\n";
        out << "E_ion_Ha     " << result.energy.ion()    << "\n";
        out << "scf_steps    " << result.total_iter       << "\n";
        out << "max_force_Ha_per_bohr  " << fmax          << "\n";
        out.close();
        std::cout << "\nWrote results/gs_summary.txt\n";

        // ── Save forces CSV (all atoms, atomic + SI units) ────────────────────
        // Conversion: 1 Ha/bohr = 8.23872e-8 N = 8.23872e-2 nN
        // Reference: CODATA 2018 values
        //   1 Ha = 4.3597447222071e-18 J
        //   1 bohr = 5.29177210903e-11 m
        //   → 1 Ha/bohr = 4.3597447222071e-18 / 5.29177210903e-11 = 8.2387234983e-8 N
        const double HA_BOHR_TO_N  = 8.2387234983e-8;   // [N]
        const double HA_BOHR_TO_NN = HA_BOHR_TO_N * 1e9; // [nN]

        std::ofstream fcsv("results/forces.csv");
        fcsv << std::scientific << std::setprecision(5);
        fcsv << "# Forces on coronene C24H12 atoms — LDA ground state (E_cut=" << cfg::ECUT_HA << " Ha)\n";
        fcsv << "# Atomic units: Ha/bohr.  SI: nN (nanonewtons; 1 Ha/bohr = 82.387 nN).\n";
        fcsv << "# atom_idx, element, "
             << "Fx_Ha_per_bohr, Fy_Ha_per_bohr, Fz_Ha_per_bohr, |F|_Ha_per_bohr, "
             << "Fx_nN, Fy_nN, Fz_nN, |F|_nN\n";
        for(int i = 0; i < ions.size(); i++){
            auto f  = result.forces[i];
            double fn = std::sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
            std::string el = (i < 24 ? "C" : "H");
            fcsv << i << ", " << el << ", "
                 << f[0] << ", " << f[1] << ", " << f[2] << ", " << fn << ", "
                 << f[0]*HA_BOHR_TO_NN << ", " << f[1]*HA_BOHR_TO_NN
                 << ", " << f[2]*HA_BOHR_TO_NN << ", " << fn*HA_BOHR_TO_NN << "\n";
        }
        fcsv.close();
        std::cout << "Wrote results/forces.csv  (all " << ions.size()
                  << " atoms, Ha/bohr + nN)\n";
    }

    return 0;
}
