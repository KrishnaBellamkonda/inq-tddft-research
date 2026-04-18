/* The aim of this script is to learn about making new observables in INQ. 
 * Specifically, we are going to model HF molecule. This as we know is a polar molecule. 
 * Hence, each atom gains a partial charge. We want to monitor the charge each atom gains.
 * This is done by calculating the charge around each atom. 
 * 
 * */

#include <inq/inq.hpp>
#include <vector>
#include <fstream>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;

void insert_nf_atoms(){



}


int main(){
	
	// Step 1: Calculate the ground state of the HF molecule
	systems::cell cell = systems::cell::cubic(8.0_bohr).finite();
	systems::ions ions(cell);
	
	auto bond_length = 0.917_angstrom;
	auto zero = 0.0_angstrom;

	ions.insert("H", {zero, zero, -bond_length/2});
	ions.insert("F", {zero, zero, bond_length/2});

	systems::electrons electrons(ions,
				     options::electrons{}
				     .cutoff(30.0_Ry), 
				     input::kpoints::gamma()
				     );

	ground_state::initial_guess(ions, electrons);
	auto gs = ground_state::calculate(ions, electrons, 
				options::theory{}.non_interacting(),
				options::ground_state{}
				.energy_tolerance(1e-6_Ha)
				.mixing(0.1)
				.max_steps(1000)
				);

	auto gs_energy = gs.energy.total();
	std::cout << "Ground state energy (Ha): " << gs_energy << std::endl;

	// Step 2: Create a new observable
	
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
	
	
	return 0;

}
