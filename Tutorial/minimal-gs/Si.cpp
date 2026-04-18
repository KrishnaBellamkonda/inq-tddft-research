// The idea of this file is to learn how to write a cpp file for INQ. 
// I understand that INQ is a header only package. Meaning, the only cpp 
// file is the final configuration file that the user write (for example this file). 
// So, the user would have to compile the entire code (that includes the entire library)
// everytime they want to run a program. 

// Need to add the package
#include <inq/inq.hpp>

int main(){
	
	// These two commands make to so
	using namespace inq; // this gives access to the inq objects without having to write inq:: so on....
	using namespace inq::magnitude; // this gives access to the units in inq

	// systems is the top most layer that contains cell, electrons and ions objects
	// 1. We start by defining the cell object
	
	//auto cell =systems::cell::cubic(5.0_angstrom).periodic();
	systems::cell periodic_cell = systems::cell::cubic(5.0_angstrom).periodic();
	// I suspect, we could have initialised the cell variable by writing
	// system::cell cell(arguments);


	// 2. Then we define the ions which take in as the argument the cell
	systems::ions ions(periodic_cell);
	// Then, we insert the ions in their respective locations. I believe this can
	// be done in a multiplitude of ways. For example providing the lattice supercell
	// size and the lattice vectors, fractional coordinates etc. 
	// I am going to use cartesian coordinates to start with
	ions.insert("Si", {0.0_b, 0.0_b, 0.0_b});

	// 3. We then define the electronic system using this command
	systems::electrons electrons(ions, options::electrons{}.cutoff(30.0_Ry));

	// 4. We then make an initial guess for ground state calculation. 
	// ground_state is a first layer packages on the same layer as systems.
	// It contains a function initial guess
	// Initial guess (takes as argument ions, electrons) for the wave function
	ground_state::initial_guess(ions, electrons);


	// 5. We then calculate the ground state using a SCF iteration. This is done
	// using the ground state calculate function.
	// the calculate function takes in additional arguments. The functional to use
	// and the energy tolerance for the scf loop to terminate ( I think 
	auto gs = ground_state::calculate(ions, electrons,
		       	options::theory{}.lda(),
			options::ground_state{}.energy_tolerance(1e-6_Ha).mixing(0.1).max_steps(1000));


	return 0;
}
