// self_hartree.cpp — EMPIRICAL wavepacket self-Hartree in the ACTUAL run cell.
//
// Computes E_self = 1/2 ∫ n_WP · φ_WP with φ_WP = poisson(n_WP), where n_WP is the
// wavepacket DENSITY (a normalised Gaussian of density-std σ_ρ = σ_WP/√2, ∫=1 —
// identical to the classical projectile's Gaussian charge). The point: INQ's Poisson
// solver picks the boundary-matched kernel from the cell periodicity —
//   periodicity 3 → fully-periodic FFT  (Makov–Payne / p3)
//   periodicity 2 → Rozzi et al. 2D Coulomb-cutoff (open-z / p2)  [inq/src/solvers/poisson.hpp:190]
// so E_self is the OPEN-Z self-Hartree for PER=2 — the empirical reference that
// replaces the free-space analytic 1/(σ_WP√2π) for the p2 production runs.
//
// No GS load, no dynamics: electrons built only to supply the density basis + comm.
// Env: LJ_LX LJ_LY LJ_LZ LJ_N LJ_PERIODICITY LJ_SPACING LJ_SIGMA LJ_LAUNCH_Z.
#include <inq/inq.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>   // gaussian_density

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
static double env_d(const char* k,double d){const char*v=std::getenv(k);return v?std::atof(v):d;}
static int    env_i(const char* k,int d){const char*v=std::getenv(k);return v?std::atoi(v):d;}

int main(){
    const double HA=27.211386;
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120);
    const int N=env_i("LJ_N",82); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5), SIGMA_WP=env_d("LJ_SIGMA",0.5);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",-24.5);
    const double SIGMA_POT=SIGMA_WP/std::sqrt(2.0);   // = σ_ρ, the WP density std

    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(20).temperature(0.00862*1.0_eV),input::kpoints::gamma());

    auto basis = electrons.density().basis();
    auto nwp   = inqkit::jellium::gaussian_density(basis, {0.0,0.0,LAUNCH_Z}, SIGMA_POT);
    const double norm = operations::integral(nwp);
    auto phiwp = solvers::poisson::solve(nwp);                       // boundary-matched kernel
    const double E_self = 0.5*operations::integral_product(nwp, phiwp);

    if(electrons.root()){
        std::cout<<std::setprecision(12)
            <<"SELF_HARTREE per="<<PER<<" Lz="<<LZ<<" dx="<<SPACING<<" sigma_wp="<<SIGMA_WP
            <<" z="<<LAUNCH_Z<<" norm="<<norm
            <<" E_self_ha="<<E_self<<" E_self_ev="<<E_self*HA<<"\n";
    }
    return 0;
}
