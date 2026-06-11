/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__REAL_TIME__VIEWABLES
#define INQ__REAL_TIME__VIEWABLES

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <systems/ions.hpp>
#include <operations/overlap_diagonal.hpp>
#include <observables/current.hpp>
#include <observables/dipole.hpp>
#include <observables/forces_stress.hpp>
#include <perturbations/none.hpp>
#include <systems/electrons.hpp>
#include <real_time/crank_nicolson.hpp>
#include <real_time/etrs.hpp>
#include <utils/profiling.hpp>

#include <chrono>

namespace inq {
namespace real_time {

template <class ForcesType, class HamiltonianType, class Perturbation>
class viewables {
	bool last_iter_;
	int iter_;
	double time_;
	systems::ions const & ions_;
	systems::electrons const & electrons_;
	hamiltonian::energy const & energy_;
	ForcesType forces_;
  HamiltonianType const & ham_;
	Perturbation const & pert_;
	
public:

	viewables(bool last_iter, int iter, double time, systems::ions const & ions, systems::electrons const & electrons, hamiltonian::energy const & energy, ForcesType const & forces, HamiltonianType const & ham, Perturbation const & pert)
		:last_iter_(last_iter), iter_(iter), time_(time), ions_(ions), electrons_(electrons), energy_(energy), forces_(forces), ham_(ham), pert_(pert){
	}

	auto iter() const {
		return iter_;
	}

	auto last_iter() const {
		return last_iter_;
	}

	auto every(int every_iter) const {
		if(iter() == 0) return false;
		return (iter()%every_iter == 0) or last_iter(); 
	}

	auto root() const {
		return electrons_.full_comm().root();
	}
	
	auto time() const {
		return time_;
	}
	
	auto positions() const {
		return ions_.positions();
	}
	
	auto velocities() const {
		return ions_.velocities();
	}

	auto forces() {
		if(forces_.size() == 0) forces_ = observables::forces_stress{ions_, electrons_, ham_, energy_}.forces;
		return forces_;
	}

	auto energy() const {
		return energy_;
	}

	auto dipole() const {
		return observables::dipole(ions_, electrons_);
	}

	auto laser_field() const {
		return pert_.uniform_electric_field(time_);
	}

	auto uniform_vector_potential() const{
		return ions_.cell().metric().to_cartesian(ham_.uniform_vector_potential());
	}

	auto num_electrons() const {
		return operations::integral(electrons_.density());
	}

	auto magnetization() const {
		return observables::total_magnetization(electrons_.spin_density());
	}

	auto uniform_magnetic_field() const {
		return pert_.uniform_magnetic_field(time_);
	}

  auto current() const {
    return ions_.cell().metric().to_cartesian(observables::current(ions_, electrons_, ham_));
  }
	auto projected_occupation(const systems::electrons & gs) {
		auto calc = [] (auto occ, auto v) {
			return occ*norm(v);
		};
		gpu::array<double, 2> occ({gs.kpin_part().size(), gs.kpin()[0].set_size()}, 0.0);
		for(int ilot = 0; ilot < gs.kpin_size(); ilot++) {

			auto ortho = matrix::all_gather(operations::overlap(electrons_.kpin()[ilot], gs.kpin()[ilot]));

			for (int it = 0; it < get<0>(sizes(ortho)); it++) {
				auto start = electrons_.kpin()[ilot].set_part().start();
				auto finish = electrons_.kpin()[ilot].set_part().end();
				occ[ilot + gs.kpin_part().start()][it] = operations::sum(electrons_.occupations()[ilot], ortho[it]({start, finish}), calc)/electrons_.kpin_weights()[ilot];
			}
		}

		if(electrons_.kpin_states_comm().size() > 1){
			electrons_.kpin_states_comm().all_reduce_in_place_n(raw_pointer_cast(occ.data_elements()), occ.num_elements(), std::plus<>{});
		}

		return occ;
	}
	
// ============ projected_occupation ============	
auto projected_occupation_array(const systems::electrons & gs) {

    auto n_kpoints = gs.kpin_part().size();
    auto n_td_states = electrons_.kpin()[0].set_size(); 
    auto n_gs_states = gs.kpin()[0].set_size();        


    gpu::array<double, 3> proj_array({n_kpoints, n_td_states, n_gs_states}, 0.0);
    auto k_start_global = gs.kpin_part().start();

    for(int ilot = 0; ilot < gs.kpin_size(); ilot++) {
        // ortho[it][igs] = <TD_i | GS_j>
        auto ortho = matrix::all_gather(operations::overlap(electrons_.kpin()[ilot], gs.kpin()[ilot]));

        double weight = electrons_.kpin_weights()[ilot];
        
        // f_i^{TD}
        const auto& td_occs = electrons_.occupations()[ilot]; 

        auto start = electrons_.kpin()[ilot].set_part().start();
        auto finish = electrons_.kpin()[ilot].set_part().end();

        for (int i_loc = 0; i_loc < (finish - start); i_loc++) {
            int i_global = start + i_loc;
            double f_i_td = td_occs[i_loc]; 

            for (int j_gs = 0; j_gs < n_gs_states; j_gs++) {
                
                // P(i,j) = f_i^{TD} * |<TD_i | GS_j>|^2
                using std::norm;
                double contribution = f_i_td * norm(ortho[i_global][j_gs]);

                proj_array[ilot + k_start_global][i_global][j_gs] = contribution / weight;
            }
        }
    }

    if(electrons_.kpin_states_comm().size() > 1){
        electrons_.kpin_states_comm().all_reduce_in_place_n(
            raw_pointer_cast(proj_array.data_elements()), 
            proj_array.num_elements(), 
            std::plus<>{}
        );
    }

    return proj_array;
}

	auto electrons() const {
		return electrons_;
	}

// ============ energy expectation ============

private:
	mutable std::vector<gpu::array<double, 1>> exp_cache_;
	mutable std::vector<gpu::array<double, 1>> var_cache_;
	mutable std::vector<double> host_print_buffer_; 

public:

	auto const& state_energy_expectations() const {
		CALI_CXX_MARK_FUNCTION;
		
		int nk = electrons_.kpin_size();
		
		if (exp_cache_.size() != nk) {
			exp_cache_.resize(nk);
		}
		
		for(int ik = 0; ik < nk; ik++) {
			auto const& phi = electrons_.kpin()[ik];
			
			auto hphi = ham_(phi);
			
			//  <φ_i|H|φ_i> 
			auto diag_complex = operations::overlap_diagonal(phi, hphi);
			auto nstates = diag_complex.size();
			

			if (exp_cache_[ik].size() != nstates) {
				exp_cache_[ik].reextent({nstates});
			}
			
			for(int ist = 0; ist < nstates; ist++) {
				exp_cache_[ik][ist] = real(diag_complex[ist]);
			}
		}
		
		return exp_cache_;
	}
	

	auto const& simple_gathered_state_energy_expectations() const {

		state_energy_expectations();
		
		if(electrons_.kpin_states_comm().size() > 1) {
			for(int ik = 0; ik < exp_cache_.size(); ik++) {
				auto* data_ptr = raw_pointer_cast(exp_cache_[ik].data_elements());
				auto size = exp_cache_[ik].size();
				
				electrons_.kpin_states_comm().all_reduce_in_place_n(data_ptr, size, std::plus<>{});
			}
		}
		
		return exp_cache_;
	}

	auto gathered_state_energy_expectations() const {
		return simple_gathered_state_energy_expectations();
	}


	void write_state_energies(std::ostream& out) const {
		if(!electrons_.full_comm().root()) return;
		

		auto const& expectations = simple_gathered_state_energy_expectations();
		
		out << "\n=== State Energy Expectations at time " << time_ << " ===\n";
		for(int ik = 0; ik < expectations.size(); ik++) {
			auto size = expectations[ik].size();
			
			if(host_print_buffer_.size() < size) {
				host_print_buffer_.resize(size);
			}
			

			for(int i = 0; i < size; ++i) {
				host_print_buffer_[i] = expectations[ik][i];
			}

			out << "k-point " << ik << " (weight: " << electrons_.kpin_weights()[ik] << "):\n";
			for(int ist = 0; ist < size; ist++) {
				out << "  State " << std::setw(4) << ist 
					<< ": E = " << std::setw(12) << std::setprecision(6) << std::fixed 
					<< host_print_buffer_[ist] << " Ha"
					<< " (occ: " << std::setw(8) << std::setprecision(4) 
					<< electrons_.occupations()[ik][ist] << ")\n";
			}
		}
		out << std::endl;
	}

	auto const& state_energy_variance() const {
		CALI_CXX_MARK_FUNCTION;
		
		int nk = electrons_.kpin_size();
		
		if (var_cache_.size() != nk) {
			var_cache_.resize(nk);
		}
		
		for(int ik = 0; ik < nk; ik++) {
			auto const& phi = electrons_.kpin()[ik];
			
			auto hphi = ham_(phi);
			auto hhphi = ham_(hphi);
			
			auto h_diag = operations::overlap_diagonal(phi, hphi);
			auto h2_diag = operations::overlap_diagonal(phi, hhphi);
			
			auto nstates = h_diag.size();
			if(var_cache_[ik].size() != nstates) {
				var_cache_[ik].reextent({nstates});
			}
			
			for(int ist = 0; ist < nstates; ist++) {
				double h_exp = real(h_diag[ist]);
				double h2_exp = real(h2_diag[ist]);
				var_cache_[ik][ist] = h2_exp - h_exp * h_exp;
			}
		}
		
		return var_cache_;
	}
	
	auto const& simple_gathered_state_energy_variance() const {
		state_energy_variance();
		
		if(electrons_.kpin_states_comm().size() > 1) {
			for(int ik = 0; ik < var_cache_.size(); ik++) {
				auto* data_ptr = raw_pointer_cast(var_cache_[ik].data_elements());
				auto size = var_cache_[ik].size();
				
				electrons_.kpin_states_comm().all_reduce_in_place_n(data_ptr, size, std::plus<>{});
			}
		}
		
		return var_cache_;
	}

	auto gathered_state_energy_variance() const {
		return simple_gathered_state_energy_variance();
	}

	void write_state_energy_variance(std::ostream& out) const {
		if(!electrons_.full_comm().root()) return;
		
		auto const& variances = simple_gathered_state_energy_variance();
		
		out << "\n=== State Energy Variance at time " << time_ << " ===\n";
		for(int ik = 0; ik < variances.size(); ik++) {
			auto size = variances[ik].size();
			
			if(host_print_buffer_.size() < size) host_print_buffer_.resize(size);
			for(int i = 0; i < size; ++i) host_print_buffer_[i] = variances[ik][i];

			out << "k-point " << ik << " (weight: " << electrons_.kpin_weights()[ik] << "):\n";
			for(int ist = 0; ist < size; ist++) {
				out << "  State " << std::setw(4) << ist 
					<< ": ΔE² = " << std::setw(12) << std::setprecision(6) << std::scientific 
					<< host_print_buffer_[ist] << " Ha²"
					<< " (σ = " << std::setw(10) << std::setprecision(4) 
					<< sqrt(std::max(0.0, host_print_buffer_[ist])) << " Ha)\n";
			}
		}
		out << std::endl;
	}

	void write_state_energy_and_variance(std::ostream& out) const {
		if(!electrons_.full_comm().root()) return;
		
		auto const& expectations = simple_gathered_state_energy_expectations();
		auto const& variances = simple_gathered_state_energy_variance();
		
		out << "\n=== State Energy and Variance at time " << time_ << " ===\n";
		
		std::vector<double> host_exp;
		std::vector<double> host_var;

		for(int ik = 0; ik < expectations.size(); ik++) {
			auto size = expectations[ik].size();
			host_exp.resize(size);
			host_var.resize(size);
			
			for(int i = 0; i < size; ++i) host_exp[i] = expectations[ik][i];
			for(int i = 0; i < size; ++i) host_var[i] = variances[ik][i];
			
			out << "k-point " << ik << " (weight: " << electrons_.kpin_weights()[ik] << "):\n";
			for(int ist = 0; ist < size; ist++) {
				out << "  State " << std::setw(4) << ist+1 
					<< ": E = " << std::setw(12) << std::setprecision(6) << std::fixed 
					<< host_exp[ist] << " Ha"
					<< ", ΔE² = " << std::setw(12) << std::setprecision(6) << std::scientific 
					<< host_var[ist] << " Ha²"
					<< " (σ = " << std::setw(10) << std::setprecision(4) 
					<< sqrt(std::max(0.0, host_var[ist])) << " Ha)" << "\n";
			}
		}
		out << std::endl;
	}

	
};

}
}
#endif
