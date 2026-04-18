// In this file, I am going to write the code for running a ground state INQ simulation for
// N2 molecule. 

#include <inq/inq.hpp>
#include <iostream>
#include <iomanip>

int main(){

	// Namespace inq added for access to ions, electrons, ground_state etc. 
	using namespace inq;
	// Name space inq::magnitude for unit values
	using namespace inq::magnitude;

	auto half_bond = 1.10_angstrom/2; // N-N bond length
	

	// 1. Initialise the system box and ions 
	systems::ions ions(systems::cell::cubic(20.0_bohr).finite());
	ions.insert("N", {0.0_angstrom, 0.0_angstrom, -half_bond});
	ions.insert("N", {0.0_angstrom, 0.0_angstrom, half_bond});

	// 2. Add electronic structure
	// This function takes care of the KS orbitals. The specification for
	// these orbitals is provided in the options
	systems::electrons electrons(ions, options::electrons{}.cutoff(80.0_Ry) );

	// 3. SCF 
	// Initialise ground state wavefunction
	ground_state::initial_guess(ions, electrons);
	// Here we specify the XC functional and an option to calculate forces
	auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(),
						options::ground_state{}.calculate_forces());

	// 4. Results
	std::cout<< std::fixed << std::setprecision(6);

	std::cout << "\n=== Energy breakdown (Ha) ===\n";
        std::cout << "  Total        : " << result.energy.total()       << "\n";
    	std::cout << "  Kinetic      : " << result.energy.kinetic()     << "\n";
    	std::cout << "  Hartree      : " << result.energy.hartree()     << "\n";
    	std::cout << "  XC (PBE)     : " << result.energy.xc()         << "\n";
    	std::cout << "  External     : " << result.energy.external()    << "\n";
    	std::cout << "  Non-local PP : " << result.energy.non_local()   << "\n";
    	std::cout << "  Ion-ion      : " << result.energy.ion()         << "\n";
    	std::cout << "  SCF iters    : " << result.total_iter           << "\n";

	// 5. Convert total energy to eV for reference
	// All the results are given in atomic units. So we convert them to the requried usits
	double total_eV = result.energy.total() / in_atomic_units(1.0_eV);
	std::cout << "\n  Total		: " << total_eV << "  eV\n";



	// 5. Forces
	std::cout << "\n ===  Forces (Ha/bohr) === \n";
	for (int i=0; i< ions.size(); i++){
		std::cout << "	N[" << i << "]  Fx= " << result.forces[i][0]
			<< "   Fy= " << result.forces[i][1]
			<< "   Fz= " << result.forces[i][2] << "\n";
	}
	
	// ── 6. Dipole ───────────────────────────────────────────────────────────
    	// H2 is homonuclear so the dipole should be ~0. Good sanity check.
    	std::cout << "\n=== Dipole moment (a.u.) ===\n";
    	std::cout << "  d = (" << result.dipole[0] << ", "
                           << result.dipole[1] << ", "
                           << result.dipole[2] << ")\n";



	return 0;
}
