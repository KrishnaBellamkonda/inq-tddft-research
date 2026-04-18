/* The aim of this script is to understand the ions class in INQ
 * Specifically, to understand how fractional and cartesian coordinates work
 * Cartesian: Present the exact coordinates in which the ions (the nuclei) are to be placed
 * Fractional: Uses the lattice vectors (explicitly defined in the cell by using 
 *	       systems::cell::lattice(vec1, vec2, vec3), or implicitly calculated when using
 *	       systems::cell::cubic or systems:cell:orthorhombic. Fractional coordinates
 *	       spawn the ions in coordinates (f1, f2, f3) as defined in the basis defined by
 *	       (vec1, vec2, vec3) for the lattice vectors. 
 * */


#include <inq/inq.hpp>
#include <vector>
#include <iostream>

#include <fstream>
#include <iomanip>
#include <string>

using namespace inq;
using namespace inq::magnitude;

// A good exercise to run would be a quantum kick here. This way, I will be able to define
// the lattice, first and foremost using cartesian coordinates and then by using fractional
// coordinates. Then, I would be able to quantum kick them to get a result. 

auto make_tddft_output(double gs_energy_ha,
                       std::ofstream& energy_out,
                       std::ofstream& current_out){

    energy_out  << "# iter time_au total_energy_Ha excess_energy_Ha\n";
    current_out << "# iter time_au Jx Jy Jz\n";

    return [&](auto data){
        const int iter = data.iter();
        const double time_au = data.time();
        const double etot = data.energy().total();
        const double de = etot - gs_energy_ha;
        const auto J = data.current();

        energy_out
            << iter << " "
            << std::setprecision(16) << time_au << " "
            << etot << " "
            << de << "\n";

        current_out
            << iter << " "
            << std::setprecision(16) << time_au << " "
            << J[0] << " "
            << J[1] << " "
            << J[2] << "\n";
    };
}

void run_impulsive_tddft(systems::ions& ions,
                         systems::electrons& electrons,
                         double gs_energy_ha,
                         int num_steps,
                         auto dt){

    std::ofstream energy_out("energy_vs_time.dat");
    std::ofstream current_out("current_vs_time.dat");

    auto output = make_tddft_output(gs_energy_ha, energy_out, current_out);

    real_time::propagate(
        ions,
        electrons,
        output,
        options::theory{}.pbe(),
        options::real_time{}
            .num_steps(num_steps)
            .dt(dt)
            .impulsive()
            .observables_current()
    );
}

void print_coordinates(const auto& input_vector){
	auto it = input_vector.begin();
    	std::cout << it[0] << ", " << it[1] << ", " << it[2] << std::endl;
}


void set_uniform_ionic_kick_x(systems::ions& ions, double v_kick){
    for(int i = 0; i < ions.size(); ++i){
        ions.velocities()[i] = vector3<double>{v_kick, 0.0, 0.0};
    }
}

void add_ions_by_cartesian(systems::ions& ions){
	// We are making a 3x3x3 supercell of Li BCC. We need to find the coordinates
	// for such a task. 


	// Lattice parameter
	auto lattice_parameter = 3.51_angstrom;
	auto half_lattice_parameter = lattice_parameter/2;
	

	// Insert the ions into this class using a for loop
	for (int k=0; k<3; k++){
		for (int j=0; j<3; j++){
			for (int i=0; i<3; i++){
				// Each combination of i, j and k defines a unit cell within
				// the supercell. Each unit cell has two atoms at the basis
				// functions defined above.
				
				std::cout << "Unit cell " << i << "," << j << "," << k << std::endl;
				
				auto atom_1 = vector3{i*lattice_parameter, 
					       j*lattice_parameter,
					       k*lattice_parameter};

				auto atom_2 = vector3{i*lattice_parameter + half_lattice_parameter,
					       j*lattice_parameter + half_lattice_parameter,
					       k*lattice_parameter + half_lattice_parameter};
				
				print_coordinates(atom_1);
				print_coordinates(atom_2);

				ions.insert("Li", atom_1);
				ions.insert("Li", atom_2);

				std::cout << "-------------------";

			}
		
		}
	}

	return;
}

int main(){
	

	// System
	systems::cell cell = systems::cell::cubic(10.53_angstrom).periodic();
	
	// Ions
	systems::ions ions = systems::ions(cell);
	add_ions_by_cartesian(ions);

	// Electronic system
	// Options for the electronic system include extra states and temperature
	systems::electrons electrons(ions, 
				     options::electrons{}
				     .cutoff(74.0_Ry)
				     .extra_states(20)
				     .temperature(0.086_eV));

	// Ground state
	ground_state::initial_guess(ions, electrons);
	auto gs = ground_state::calculate(ions, electrons,
					options::theory{}.pbe(),
		options::ground_state{}.energy_tolerance(1e-6_Ha).mixing(0.5));

		
	// Report the value
	std::cout << "Number of ions: " << ions.size() << '\n';
	std::cout << "GS Total energy: " << gs.energy.total() << std::endl;

	// Ground-state total energy
	double gs_energy = gs.energy.total();

	// Quantum kick and TDDFT
	double v_kick = 0.05028; // bohr / atomic time
	set_uniform_ionic_kick_x(ions, v_kick);

	// Real-time propagation
	run_impulsive_tddft(
    	ions,
    	electrons,
    	gs_energy,
    	500,
    	0.04_atomictime
	);

	return 0;
}
