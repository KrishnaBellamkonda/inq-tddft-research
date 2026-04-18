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
#include <complex>

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

// 3. Calculate charge
std::vector<double> calculate_relative_charge(systems::electrons const& electrons){

        // Creating a new observable
        /* Density is an instance of the INQs inq::basis::field class
           The total electron density is represented as a scalar 3d field sampled 
           on a real space grid. (I believe the inq::electrons holds spin_density
           within it. When the internal function .density() is called, then it converts
           the spin_density into a total density. 
         */ 
        auto density = electrons.density();
           


        /* The basis function contains the grid dimensions via sizes() function. 
         * It also contains the mapping from grid indices (as the density stored in a 1d array)
           , the mappings between the indices to the points are stored in this basis object.
         * They also contain the geometry/grid metadata
         */
        auto basis = density.basis();


        auto h_charge = 0.0;
        auto f_charge = 0.0;

        for (int ix=0; ix<basis.sizes()[0]; ix++){
            for (int iy=0; iy<basis.sizes()[1]; iy++){
                for (int iz=0; iz<basis.sizes()[2]; iz++){
                        // The following function converts the given indices ix, iy and iz 
                        // for the density grid into cartesian coordinates
                        auto rr = basis.point_op().rvector_cartesian(ix, iy, iz);
                        // If the cartesian coordinates reveal that the point is closer to H
                        // or F, this density is added to the corresponding atom

                        // density.cubic() turns a 1D object into 3D with each point storing
                        // the corresponding charge density
                        if (rr[2] < 0.0) h_charge += density.cubic()[ix][iy][iz];
                        if (rr[2] >= 0.0) f_charge += density.cubic()[ix][iy][iz];
                }
            }
        }

        // To get the charge from charge density, we multiply by the volume element of the grid
        h_charge *= basis.volume_element();
        f_charge *= basis.volume_element();


        std::cout << "H charge = " << h_charge - 1.0 << '\n'; // delta charge

        std::cout << "F charge = " << f_charge - 7.0 << '\n';

        std::vector<double> relative_charges {h_charge, f_charge}; 

        return relative_charges;

}


// Task 4: Write Ground-State Orbital snapshots in a human readable format
void write_gs_orbital_metadata(systems::electrons const& electrons) {
    auto const& basis = electrons.states_basis();

    if (!basis.comm().root()) return;

    std::filesystem::create_directories("results/gs-orbitals/human-readable-orbitals");

    std::ofstream out("results/gs-orbitals/human-readable-orbitals/metadata.dat");

    out << "# Ground-state human-readable orbital export metadata\n";
    out << "# global_grid nx ny nz\n";
    out << basis.sizes()[0] << ' '
        << basis.sizes()[1] << ' '
        << basis.sizes()[2] << '\n';

    out << "# volume_element\n";
    out << basis.volume_element() << '\n';

    out << "# number_of_kpin_blocks\n";
    out << electrons.kpin().size() << '\n';

    out << "# columns in orbital text files:\n";
    out << "# ix_global iy_global iz_global x y z re im abs2\n";
}

void write_gs_orbital_snapshot(systems::electrons const& electrons,
                               bool occupied_only = true) {
    auto const& basis = electrons.states_basis();
    int rank = basis.comm().rank();

    std::string basedir = "results/gs-orbitals/human-readable-orbitals";

    std::ostringstream rankdir;
    rankdir << basedir
            << "/rank_" << std::setw(3) << std::setfill('0') << rank;

    std::filesystem::create_directories(rankdir.str());

    if (basis.comm().root()) {
        std::ofstream manifest(basedir + "/manifest.dat");
        manifest << "# Ground-state orbital export\n";
        manifest << "# one subdirectory per real-space MPI rank\n";
        manifest << "# each orbital file columns are:\n";
        manifest << "# ix_global iy_global iz_global x y z re im abs2\n";
    }

    for (std::size_t ilot = 0; ilot < electrons.kpin().size(); ++ilot) {
        auto const& phi = electrons.kpin()[ilot];

        int kindex = electrons.kpoint_index(phi);
        int sindex = phi.spin_index();

        for (int ist = 0; ist < phi.spinor_set_part().local_size(); ++ist) {
            int orbital_index = phi.spinor_set_part().local_to_global(ist).value();
            double occ        = electrons.occupations()[ilot][ist];

            if (occupied_only && occ < 1.0e-12) continue;

            std::ostringstream fname;
            fname << rankdir.str()
                  << "/orbital_k"
                  << std::setw(3) << std::setfill('0') << kindex
                  << "_s" << sindex
                  << "_n" << std::setw(4) << std::setfill('0') << orbital_index
                  << ".dat";

            std::ofstream out(fname.str());

            out << "# kpoint_index " << kindex << "\n";
            out << "# spin_index " << sindex << "\n";
            out << "# orbital_index " << orbital_index << "\n";
            out << "# occupation " << occ << "\n";
            out << "# ix_global iy_global iz_global x y z re im abs2\n";

            for (int ix = 0; ix < basis.local_sizes()[0]; ++ix) {
                for (int iy = 0; iy < basis.local_sizes()[1]; ++iy) {
                    for (int iz = 0; iz < basis.local_sizes()[2]; ++iz) {
                        auto ixg = basis.cubic_part(0).local_to_global(ix).value();
                        auto iyg = basis.cubic_part(1).local_to_global(iy).value();
                        auto izg = basis.cubic_part(2).local_to_global(iz).value();

                        auto rr  = basis.point_op().rvector_cartesian(ix, iy, iz);
                        auto psi = phi.hypercubic()[ix][iy][iz][ist];

                        out << ixg << ' '
                            << iyg << ' '
                            << izg << ' '
                            << rr[0] << ' '
                            << rr[1] << ' '
                            << rr[2] << ' '
                            << real(psi) << ' '
                            << imag(psi) << ' '
                            << norm(psi) << '\n';
                    }
                }
            }
        }
    }
}



/* Task 1: Monitoring charge distribution of an extended bond length over time
 *
 * */

int main(){
        std::filesystem::create_directories("results/density");
        std::filesystem::create_directories("results/ion-positions");
        std::filesystem::create_directories("results/gs-orbitals/human-readable-orbitals");


        // Step 1: Calculate the ground state of the HF molecule
        systems::cell cell = systems::cell::cubic(8.0_bohr).finite();
        systems::ions ions(cell);

        auto original_bond_length = 0.917_angstrom;
        auto bond_length = original_bond_length + 0.1_angstrom; // Don't know if this is enough to cause oscillations
        auto zero = 0.0_angstrom;

        ions.insert("H", {zero, zero, -bond_length/2});
        ions.insert("F", {zero, zero, bond_length/2});

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

        // Task 2: Saving the ground state orbitals to a human readable file that can be read and visualised
        write_gs_orbital_metadata(electrons);
        write_gs_orbital_snapshot(electrons, true);


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
                             .num_steps(11000)
                             .dt(0.0326_atomictime)
                        );


        return 0;

}
