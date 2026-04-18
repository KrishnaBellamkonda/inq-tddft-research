// Test: verify hypercubic() index ordering and rvector_cartesian units in INQ
//
// Answers:
//  1. rvector_cartesian(ix,iy,iz) units, axis mapping, subtraction validity
//  2. hypercubic()[i0][i1][i2][ist] ordering — which index is x, y, z?
//  3. flat matrix[][ist] vs hypercubic: which grid axis is fastest-varying?
//
// Tests 1 and 2 require no converged ground state.
// Test 3 (wavefunction peak) uses a converged N2 system.

#include <inq/inq.hpp>

int main(){
	using namespace inq;
	using namespace inq::magnitude;

	// ── TEST 1 + 2: coordinate and hypercubic ordering ───────────────────────
	// Asymmetric cell: Lx=10, Ly=12, Lz=14 bohr → spacings differ → axes unambiguous
	{
		auto ions = systems::ions(systems::cell::orthorhombic(10.0_b, 12.0_b, 14.0_b));
		ions.insert("He", {0.0_b, 0.0_b, 0.0_b});

		auto electrons = systems::electrons(ions,
			options::electrons{}.cutoff(20.0_Ha).extra_states(3));

		if(electrons.root()){
			auto & basis = electrons.states_basis();
			auto po = basis.point_op();
			int Nx = basis.sizes()[0], Ny = basis.sizes()[1], Nz = basis.sizes()[2];

			// ── TEST 1: coordinate units and axis mapping ─────────────────
			std::cout << "\n=== TEST 1: rvector_cartesian ===\n";
			std::cout << "Cell Lx=10 Ly=12 Lz=14 bohr   Grid " << Nx << "x" << Ny << "x" << Nz << "\n";
			auto r000 = po.rvector_cartesian(0, 0, 0);
			auto r100 = po.rvector_cartesian(1, 0, 0);
			auto r010 = po.rvector_cartesian(0, 1, 0);
			auto r001 = po.rvector_cartesian(0, 0, 1);
			std::cout << "rvector_cartesian(0,0,0) = [" << r000[0] << ", " << r000[1] << ", " << r000[2] << "]\n";
			std::cout << "rvector_cartesian(1,0,0) = [" << r100[0] << ", " << r100[1] << ", " << r100[2] << "]"
			          << "   (expect x≈+" << 10.0/Nx << ", y=0, z=0)\n";
			std::cout << "rvector_cartesian(0,1,0) = [" << r010[0] << ", " << r010[1] << ", " << r010[2] << "]"
			          << "   (expect x=0, y≈+" << 12.0/Ny << ", z=0)\n";
			std::cout << "rvector_cartesian(0,0,1) = [" << r001[0] << ", " << r001[1] << ", " << r001[2] << "]"
			          << "   (expect x=0, y=0, z≈+" << 14.0/Nz << ")\n";
			auto dr = r100 - r000;
			std::cout << "r100 - r000 = [" << dr[0] << ", " << dr[1] << ", " << dr[2]
			          << "]   subtraction valid: yes\n\n";

			// ── TEST 2: hypercubic marker test (no GS needed) ────────────
			// Write marker 111 at flat index 0, marker 222 at flat index 1.
			// Then scan hypercubic to find where 222 lands.
			// This tells us: which hypercubic index is fastest-varying (= z, y, or x).
			std::cout << "=== TEST 2: flat matrix → hypercubic ordering ===\n";

			auto & phi = electrons.kpin()[0];
			int Nx_loc = phi.basis().cubic_part(0).local_size();
			int Ny_loc = phi.basis().cubic_part(1).local_size();
			int Nz_loc = phi.basis().cubic_part(2).local_size();
			int Nst_loc = phi.set_part().local_size();

			std::cout << "Hypercubic local sizes from cubic_part:\n";
			std::cout << "  cubic_part(0) = " << Nx_loc << "\n";
			std::cout << "  cubic_part(1) = " << Ny_loc << "\n";
			std::cout << "  cubic_part(2) = " << Nz_loc << "\n";
			std::cout << "  set_part (states) = " << Nst_loc << "\n\n";

			// Write markers
			phi.matrix()[0][0] = complex{111.0, 0.0};
			phi.matrix()[1][0] = complex{222.0, 0.0};
			phi.matrix()[2][0] = complex{333.0, 0.0};

			auto hc = phi.hypercubic();

			std::cout << "matrix[0][0]=111  → hc[0][0][0][0] = " << hc[0][0][0][0].real() << "\n";

			// Find 222 and 333
			int f2_i0=-1,f2_i1=-1,f2_i2=-1;
			int f3_i0=-1,f3_i1=-1,f3_i2=-1;
			for(int i0=0;i0<Nx_loc;i0++) for(int i1=0;i1<Ny_loc;i1++) for(int i2=0;i2<Nz_loc;i2++){
				double v = hc[i0][i1][i2][0].real();
				if(std::abs(v-222.0)<1.0){ f2_i0=i0; f2_i1=i1; f2_i2=i2; }
				if(std::abs(v-333.0)<1.0){ f3_i0=i0; f3_i1=i1; f3_i2=i2; }
			}
			std::cout << "matrix[1][0]=222  → hc[" << f2_i0 << "][" << f2_i1 << "][" << f2_i2 << "][0]\n";
			std::cout << "matrix[2][0]=333  → hc[" << f3_i0 << "][" << f3_i1 << "][" << f3_i2 << "][0]\n\n";

			// Determine fastest-varying index
			// If flat[1]→hc[0][0][1]: z is fastest (ip=ix*Ny*Nz + iy*Nz + iz)
			// If flat[1]→hc[0][1][0]: y is fastest
			// If flat[1]→hc[1][0][0]: x is fastest
			std::cout << "Conclusion:\n";
			if(f2_i2==1 && f2_i0==0 && f2_i1==0)
				std::cout << "  flat[1] → hc[0][0][1]  →  iz (third index) is fastest-varying\n"
				          << "  storage: ip = i0*Ny_loc*Nz_loc + i1*Nz_loc + i2\n"
				          << "  NOTE: i0 = ix (x-axis), i1 = iy (y-axis), i2 = iz (z-axis)\n"
				          << "  → to use with rvector_cartesian: rvector_cartesian(i0, i1, i2)\n";
			else if(f2_i1==1 && f2_i0==0 && f2_i2==0)
				std::cout << "  flat[1] → hc[0][1][0]  →  i1 (second index) is fastest-varying\n";
			else if(f2_i0==1 && f2_i1==0 && f2_i2==0)
				std::cout << "  flat[1] → hc[1][0][0]  →  i0 (first index) is fastest-varying\n";
			else
				std::cout << "  flat[1] → hc[" << f2_i0 << "][" << f2_i1 << "][" << f2_i2 << "] — unexpected, check manually\n";
		}
	}

	// ── TEST 3: wavefunction peak on converged N2 ground state ───────────────
	// N2 along z-axis: atoms at (0,0,±d/2). The HOMO should peak near atoms.
	// We check that the hypercubic peak is at the expected Cartesian position.
	{
		const double bond = 2.074;   // N≡N bond length in bohr (plain double)
		auto ions = systems::ions(systems::cell::cubic(20.0_b).finite());
		ions.insert("N", {0.0_b, 0.0_b,  (bond/2)*1.0_b});
		ions.insert("N", {0.0_b, 0.0_b, -(bond/2)*1.0_b});

		auto electrons = systems::electrons(ions,
			options::electrons{}.cutoff(40.0_Ha));  // no extra_states → stable SCF

		ground_state::initial_guess(ions, electrons);
		ground_state::calculate(ions, electrons,
			options::theory{}.lda(),
			options::ground_state{}.max_steps(200).energy_tolerance(1e-6_Ha));

		if(electrons.root()){
			auto & phi = electrons.kpin()[0];
			auto & basis = phi.basis();
			int Nx_loc = basis.cubic_part(0).local_size();
			int Ny_loc = basis.cubic_part(1).local_size();
			int Nz_loc = basis.cubic_part(2).local_size();
			int Nst = phi.spinor_set_size();
			auto po = basis.point_op();

			std::cout << "\n=== TEST 3: N2 wavefunction peak locations ===\n";
			std::cout << "N2: atoms at z=+" << bond/2 << " and z=-" << bond/2 << " bohr\n";
			std::cout << "Number of states: " << Nst << "\n\n";

			auto hc = phi.hypercubic();

			// Find peak of each occupied state
			for(int ist = 0; ist < phi.set_part().local_size(); ist++){
				auto istg = phi.set_part().local_to_global(ist).value();
				double max_v = -1.0;
				int mx=-1,my=-1,mz=-1;
				for(int i0=0;i0<Nx_loc;i0++) for(int i1=0;i1<Ny_loc;i1++) for(int i2=0;i2<Nz_loc;i2++){
					double v = norm(hc[i0][i1][i2][ist]);
					if(v>max_v){ max_v=v; mx=i0; my=i1; mz=i2; }
				}
				auto r = po.rvector_cartesian(mx, my, mz);
				std::cout << "State " << istg << " peak: hc[" << mx << "][" << my << "][" << mz << "]"
				          << "  r=[" << r[0] << ", " << r[1] << ", " << r[2] << "] bohr"
				          << "  |psi|²=" << max_v << "\n";
			}
			std::cout << "\n(σ and π states should peak near z=±" << bond/2 << " bohr;\n"
			          << " non-bonding states may peak at z=0)\n";
		}
	}

	return 0;
}
