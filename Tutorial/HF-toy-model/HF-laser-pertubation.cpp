/*
 * The aim of this script is to visualise simple objects and start building an intuition
 * for - 
 * 1. orbitals/wavefunctions in real space
 * 2. harmonic oscillation of stretched molecules
 * 3. oscillation of molecules upon being incident by light (INQ laser pertubation)
 * 4. density cloud visulisation
 * */



#include <inq/inq.hpp>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <filesystem>

using namespace inq;
using namespace inq::magnitude;


/* Utility functions for data processing */

// 1. Density function
void write_density_snapshot(systems::electrons const& electrons, int step) {
    auto density = electrons.density();
    auto basis   = density.basis();

    std::ostringstream name;

    // This is the name of the file that is being created
    name << "results/density/density_step_" << std::setw(6) << std::setfill('0') << step << ".dat";
    std::ofstream out(name.str());

    out << "# ix iy iz x y z rho\n";

    for (int ix = 0; ix < basis.sizes()[0]; ix++) {
        for (int iy = 0; iy < basis.sizes()[1]; iy++) {
            for (int iz = 0; iz < basis.sizes()[2]; iz++) {
                // rr is the cartesian coordinate stored at this index
                auto rr  = basis.point_op().rvector_cartesian(ix, iy, iz);
                // density value stored at this position
                auto rho = density.cubic()[ix][iy][iz];

                out << ix << ' '
                    << iy << ' '
                    << iz << ' '
                    << rr[0] << ' '
                    << rr[1] << ' '
                    << rr[2] << ' '
                    << rho   << '\n';
            }
        }
    }
}


// 2. Positions Snapshot
void write_positions_snapshot(const std::vector<vector3<double>>& positions, int step, double time) {
    std::ostringstream name;
    name << "results/ion-positions/positions_step_" << std::setw(6) << std::setfill('0') << step << ".dat";
    std::ofstream out(name.str());

    out << "# step " << step << "\n";
    out << "# time " << time << "\n";
    out << "# atom x y z\n";

    for (std::size_t i = 0; i < positions.size(); i++) {
        out << i << ' '
            << positions[i][0] << ' '
            << positions[i][1] << ' '
            << positions[i][2] << '\n';
    }
}


/* Task 1: Monitoring charge distribution of an extended bond length over time
 *
 * */

int main(){
        std::filesystem::create_directories("results/density");
        std::filesystem::create_directories("results/ion-positions");


        // Step 1: Calculate the ground state of the HF molecule
        systems::cell cell = systems::cell::cubic(8.0_bohr).finite();
        systems::ions ions(cell);

        auto bond_length = 0.917_angstrom;
        auto zero = 0.0_angstrom;

        ions.insert("H", {zero, zero, -bond_length/2});
        ions.insert("F", {zero, zero,  bond_length/2});

        systems::electrons electrons(ions,
                                     options::electrons{}
                                     .cutoff(40.0_Ry)
                                     );

        ground_state::initial_guess(ions, electrons);
        auto gs = ground_state::calculate(ions, electrons,
                                options::theory{}.pbe(),
                                options::ground_state{}
                                .energy_tolerance(1e-6_Ha)
                                .mixing(0.1)
                                .max_steps(1000)
                                );

        auto gs_energy = gs.energy.total();
        std::cout << "Ground state energy (Ha): " << gs_energy << std::endl;

        // HF fundamental IR stretch ~ 3961 cm^-1 = 0.4911 eV = 118.75 THz
        auto laser_frequency = 0.4911_eV;
        auto laser = perturbations::laser({0.0, 0.0, 1.0e-3},
                                          laser_frequency,
                                          perturbations::gauge::length);

        // Step 2: propagate the system in real time and observe this new observable
        int n_timesteps_per_write = 50;
        int n_timesteps_per_log = 1;
        auto data_processor = [&](auto const& data){
                if (data.every(n_timesteps_per_log)){
                    auto step = data.iter();
                    auto time = data.time();

                     // optional screen output
                    std::cout << "step = " << step
                              << " time = " << time
                              << " E = " << data.energy().total()
                              << '\n';
                }

                if (data.every(n_timesteps_per_write)) {
                    auto step = data.iter();
                    auto time = data.time();

                    write_positions_snapshot(data.positions(), step, time);
                    write_density_snapshot(electrons, step);
                }
        };

        real_time::propagate(ions, electrons, data_processor,
                             options::theory{}.pbe(),
                             options::real_time{}
                             .ehrenfest()
                             .observables_dipole()
                             .num_steps(20000)
                             .dt(0.0326_atomictime),
                             laser
                        );


        return 0;

}
