// ============================================================================
// 03_ecut_convergence/run.cpp
//
// Purpose: Establish the converged energy cutoff for coronene.
//
// The paper uses a grid spacing of 0.16 Å ≈ 0.302 bohr, corresponding to
// E_cut ≈ (π/0.302)² / 2 ≈ 54 Ha.  We sweep over a range of cutoffs and
// look for where the total energy and HOMO eigenvalue plateau (< 1 meV change).
//
// Output: results/ecut_convergence.csv  (E_cut, E_total, n_grid_points)
// ============================================================================

#include "config.hpp"
#include <fstream>
#include <iomanip>
#include <vector>

int main(){
    using namespace inq;
    using namespace inq::magnitude;

    // Cutoff sweep (Ha)
    // Paper target: 54 Ha (0.16 Å grid).  Start from 20 Ha to locate threshold.
    std::vector<double> ecut_list = {20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0};

    auto cell = cfg::make_cell();

    std::ofstream csv;
    if(/* root check below */ true){
        csv.open("results/ecut_convergence.csv");
        csv << "# E_cut(Ha), E_total(Ha), E_total(eV), grid_points, scf_steps\n";
    }

    for(double ecut : ecut_list){

        auto ions = systems::ions::parse(cfg::CORONENE_XYZ, cell);
        auto electrons = systems::electrons(ions,
            options::electrons{}.cutoff(ecut * 1.0_Ha)
                                .extra_states(cfg::EXTRA_STATES));

        ground_state::initial_guess(ions, electrons);
        auto result = ground_state::calculate(ions, electrons,
            options::theory{}.lda(),
            options::ground_state{}
                .energy_tolerance(cfg::SCF_TOL * 1.0_Ha)
                .max_steps(300));

        if(electrons.root()){
            int ngp = electrons.states_basis().size();
            double etot = result.energy.total();

            std::cout << std::fixed << std::setprecision(6);
            std::cout << "E_cut=" << std::setw(6) << ecut << " Ha"
                      << "  E_total=" << etot << " Ha"
                      << "  (" << etot * cfg::HA_TO_EV << " eV)"
                      << "  grid=" << ngp << " pts"
                      << "  steps=" << result.total_iter << "\n";

            csv << std::fixed << std::setprecision(8)
                << ecut << ", " << etot << ", " << etot*cfg::HA_TO_EV
                << ", " << ngp << ", " << result.total_iter << "\n";
            csv.flush();
        }
    }

    if(/* root */ true) csv.close();

    return 0;
}
