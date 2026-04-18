#include <inq/inq.hpp>
#include <vector>

int main(){
	
	using namespace inq;
	using namespace inq::magnitude;

	// Consider a cubic lattice
	auto lattice_param = 4.05_angstrom;
	auto zero = 0.0_angstrom;

	// The different types of cells that can be initialised are - 
	// I believe, cell in INQ is slightly different to how qball defines a cell
	// In QBall, cell defines the space in which the simuation happens. 
	// However, in INQ, cell contains both the simulation box and also the
	// lattice structure. 
	
	// In all of the cases listed below, I have used periodic boundary conditions
	// However, the other options available are finite()
	// periodicity() 
	
	// periodic() is the same as periodicity(3).Essentiall, the given cell 
	// in the simulation is repeated in the x, y and z directions. This is useful
	// in simulating bulk materials. 
	//
	// finite() is an individual cell that is surrounded by vacuum. This is useful in 
	// simulating isolated materials
	//
	// periodicity(2) for example repeats the unit cell in two dimensions (not the third). 
	// This is useful in simulating, for example, a 2d slab. 
		
	// 1. cubic (I don't know if this is simple cubic or FCC or BCC, etc. )
	systems::cell cubic_cell = systems::cell::cubic(lattice_param).periodic();

	// 2. orthorhombic cell
	systems::cell orthorhombic_cell = systems::cell::orthorhombic(lattice_param,
								      lattice_param,
								      lattice_param).periodic();

	// 3. general unit cell (with lattice parameters)
	systems::cell lattice_cell = systems::cell::lattice(
			vector3<quantity<magnitude::length>>{lattice_param, zero, zero},
			vector3<quantity<magnitude::length>>{zero, lattice_param, zero},
			vector3<quantity<magnitude::length>>{zero, zero, lattice_param}
			).periodic();

	// Run a ground state enegy calculation for the three different lattices and compare the answers
	//std::vector<systems::cell> cell_list = {cubic_cell, orthorhombic_cell, lattice_cell};

	// Coordinates of ions
	std::vector<systems::cell> cell_list = {cubic_cell, orthorhombic_cell, lattice_cell};
    	std::vector<std::string> cell_names = {"Cubic", "Orthorhombic", "General Lattice"};
    
    // 1. Create a vector to store the results
    std::vector<double> results;
	
	// Fix the atom coordinate types
    auto atom_1 = vector3{0.0_angstrom, 0.0_angstrom, 0.0_angstrom};
    auto atom_2 = vector3{0.5_angstrom, 0.5_angstrom, 0.0_angstrom};
    auto atom_3 = vector3{0.5_angstrom, 0.0_angstrom, 0.5_angstrom};
    auto atom_4 = vector3{0.0_angstrom, 0.5_angstrom, 0.5_angstrom};

    for (auto& cell : cell_list) {
        systems::ions ions(cell);
        ions.insert("Al", atom_1);
        ions.insert("Al", atom_2);
        ions.insert("Al", atom_3);
        ions.insert("Al", atom_4);

        systems::electrons electrons(ions, options::electrons{}.cutoff(30.0_Ry));
        ground_state::initial_guess(ions, electrons);

        auto gs = ground_state::calculate(ions,
                                         electrons,
                                         options::theory{}.lda(),
                                         options::ground_state{}.energy_tolerance(1e-6_Ha)
                                                                .mixing(0.1)
                                                                .max_steps(1000));

        // 2. Store the total energy (converted to a double)
        results.push_back(gs.energy.total());
        
        // Optional: Print a progress marker so you know the loop finished one iteration
        std::cout << "Finished one calculation..." << std::endl;
    }

    // 3. Print all results at the end for easy comparison
    std::cout << "\n--- Final Comparison ---\n";
    for (size_t i = 0; i < results.size(); ++i) {
        std::cout << cell_names[i] << " Total Energy: " << results[i] << " Ha" << std::endl;
    }

	return 0;
}
