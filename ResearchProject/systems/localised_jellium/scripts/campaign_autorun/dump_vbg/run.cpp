// ============================================================================
// localised_jellium / scripts/campaign_autorun / dump_vbg / run.cpp
//
// Dump the STATIC background potential v_bg used by the localised-jellium
// perturbation, as a z-lineout, for the theoretical_slab_model sanity check
// (campaign localised_jellium_parameter_study_2, notebook section
// "v_bg from the Poisson solver vs the infinite plate").
//
// Builds the SAME cell + background as the campaign_autorun runs, then:
//   n_+(r) = make_localised_background(basis, params)   (positive slab charge)
//   phi(r) = inq::solvers::poisson::solve(n_+)          (INQ p2 Rozzi kernel)
//   v_bg(r) = -phi(r)                                   (electron well)
// v_bg is translationally invariant in x,y (slab), so we write only the z
// lineout at x=y=0 to CSV: z_bohr, n_plus, phi_ha, v_bg_ha.
//
// Env (defaults = p2 A1 geometry): LJ_LX(50) LJ_LY(50) LJ_LZ(120) LJ_HALF(12.5)
//      LJ_N(82) LJ_EDGE_W(0) LJ_PERIODICITY(2) LJ_SPACING(0.5) LJ_OUT(vbg).
// Build against INQ_SOURCE=inq-study.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/detail/grid_layout.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <solvers/poisson.hpp>

using namespace inq;
using namespace inq::magnitude;
static double env_d(const char* k,double d){const char*v=std::getenv(k);return v?std::atof(v):d;}
static int    env_i(const char* k,int d){const char*v=std::getenv(k);return v?std::atoi(v):d;}
static std::string env_s(const char* k,const std::string&d){const char*v=std::getenv(k);return v?std::string(v):d;}

int main(){
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",82); const double EDGE_W=env_d("LJ_EDGE_W",0); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5);
    const std::string OUT=env_s("LJ_OUT","vbg");
    const double N0=double(N)/(LX*LY*(2.0*HALF));

    std::cout<<"\n=== dump_vbg per="<<PER<<" Lz="<<LZ<<" half="<<HALF<<" N0="<<N0<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(20).temperature(0.00862*1.0_eV),input::kpoints::gamma());

    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;

    auto basis = electrons.density().basis();
    auto nplus = inqkit::jellium::make_localised_background(basis, bg);
    auto phi   = inq::solvers::poisson::solve(nplus);   // phi = poisson(n_+)

    const int nx=basis.sizes()[0], ny=basis.sizes()[1], nz=basis.sizes()[2];
    auto dz = basis.rspacing()[2];
    auto z0 = basis.symmetric_range_begin(2)*dz;        // origin (-L/2)

    boost::multi::array<double,3> hphi{phi.cubic()};
    boost::multi::array<double,3> hnp{nplus.cubic()};

    // central x,y OUTPUT column (x=y=0 sits at output index n/2 for symmetric grid)
    const int ixo=nx/2, iyo=ny/2;
    const int sx=inqkit::detail::grid_layout::fft_shift_index(ixo,nx);
    const int sy=inqkit::detail::grid_layout::fft_shift_index(iyo,ny);

    std::filesystem::create_directories("results/"+OUT);
    std::ofstream f("results/"+OUT+"/vbg_lineout.csv");
    f<<std::setprecision(10)<<"z_bohr,n_plus,phi_ha,v_bg_ha\n";
    for(int izo=0;izo<nz;++izo){
        int sz=inqkit::detail::grid_layout::fft_shift_index(izo,nz);
        double z=z0+izo*dz;
        double ph=hphi[sx][sy][sz];
        double np=hnp[sx][sy][sz];
        f<<z<<","<<np<<","<<ph<<","<<(-ph)<<"\n";
    }
    std::cout<<"  wrote results/"<<OUT<<"/vbg_lineout.csv  (nz="<<nz<<", dz="<<dz<<", N0="<<N0<<")\n";
    if(electrons.root()){std::ofstream s("results/"+OUT+"/run_summary.txt");
        s<<std::setprecision(12)<<"run = campaign_autorun/dump_vbg\nengine = inq-study\n"
         <<"periodicity = "<<PER<<"  Lz = "<<LZ<<"  half_width = "<<HALF<<"  edge_width = "<<EDGE_W<<"\n"
         <<"N = "<<N<<"  n0 = "<<N0<<"  spacing = "<<SPACING<<"\nrun_completed = true\n";}
    return 0;
}
