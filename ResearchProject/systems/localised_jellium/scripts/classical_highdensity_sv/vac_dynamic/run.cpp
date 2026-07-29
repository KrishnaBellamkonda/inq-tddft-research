// ============================================================================
// classical_highdensity_sv / vac_dynamic / run.cpp
//
// Phase 1b (vacuum DYNAMIC exit test). Unlike vac_exit (static snapshots at
// prescribed centers), this is a REAL propagation: the inqkit::dynamics::Projectile
// is advanced by velocity-Verlet every timestep and the
// moving_gaussian_projectile_perturbation reads proj.R() each step and re-solves
// φ_proj = poisson(n_proj). This validates the DYNAMIC coupling — that during an
// actual real_time::propagate the perturbation tracks the moving projectile and
// leaves the z-open box smoothly.
//
// VACUUM: empty box (no jellium background), 2 spectator electrons only (needed
// for a valid basis + a propagator to drive). Const-velocity drive (zero force),
// so proj_z(t) = z0 + v*t exactly and φ_proj(t) is sampled at the real timestep
// positions. Runs until the projectile is >= Lz BEYOND the far face (+Lz/2 + Lz).
//
// Emits per SAVE_EVERY: φ_proj VTI (physical order) + a CSV row
// (step,time_au,proj_z,proj_vz,phi_peak). periodicity(2), z OPEN.
//
// Env: VD_LX(35) VD_LY(35) VD_LZ(85) VD_SPACING(0.5) VD_SIGMA_WP(0.5)
//      VD_MASS(1) VD_V(2.0) VD_LAUNCH_Z(-30) VD_DT(0.04) VD_SAVE_EVERY(0=auto)
// No inq/ or inq-study/ edit — wrapper-only.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/real_field_3d.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>          // gaussian_density
#include <inqkit/dynamics/projectile.hpp>
#include <inqkit/dynamics/moving_gaussian_projectile_perturbation.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;
static double env_d(const char* k,double d){const char*v=std::getenv(k);return v?std::atof(v):d;}
static int    env_i(const char* k,int d){const char*v=std::getenv(k);return v?std::atoi(v):d;}

// INQ real-space field (FFT-natural) -> physical-order inqkit RealField3D
// (mirrors inqkit::fields::density::total; same fft_shift so the VTI loads via
// inqview.load_vti).
template <class Field>
inqkit::fields::RealField3D field_to_physical(Field const & f) {
	auto const & basis = f.basis();
	auto const nx=basis.sizes()[0], ny=basis.sizes()[1], nz=basis.sizes()[2];
	auto const sp = basis.rspacing();
	inqkit::fields::RealField3D field;
	field.nx=nx; field.ny=ny; field.nz=nz;
	field.origin_x_bohr=basis.symmetric_range_begin(0)*sp[0];
	field.origin_y_bohr=basis.symmetric_range_begin(1)*sp[1];
	field.origin_z_bohr=basis.symmetric_range_begin(2)*sp[2];
	field.dx_bohr=sp[0]; field.dy_bohr=sp[1]; field.dz_bohr=sp[2];
	field.values.resize(static_cast<std::size_t>(nx)*ny*nz);
	boost::multi::array<double,3> hc{f.cubic()};
	for(int ix=0;ix<nx;++ix){ int sx=inqkit::detail::grid_layout::fft_shift_index(ix,nx);
		for(int iy=0;iy<ny;++iy){ int sy=inqkit::detail::grid_layout::fft_shift_index(iy,ny);
			for(int iz=0;iz<nz;++iz){ int sz=inqkit::detail::grid_layout::fft_shift_index(iz,nz);
				auto flat=inqkit::detail::grid_layout::flatten_index(ix,iy,iz,ny,nz);
				field.values[flat]=hc[sx][sy][sz]; }}}
	return field;
}
static std::string tag6(int n){std::ostringstream o;o<<std::setw(6)<<std::setfill('0')<<n;return o.str();}

int main(){
	const double LX=env_d("VD_LX",35),LY=env_d("VD_LY",35),LZ=env_d("VD_LZ",85);
	const double SPACING=env_d("VD_SPACING",0.5), SIGMA_WP=env_d("VD_SIGMA_WP",0.5);
	const double MASS=env_d("VD_MASS",1.0), V=env_d("VD_V",2.0), LAUNCH_Z=env_d("VD_LAUNCH_Z",-30.0);
	const double DT=env_d("VD_DT",0.04);
	const double SIGMA_POT=SIGMA_WP/std::sqrt(2.0);
	const double ZFAR=LZ/2.0;
	// run until proj is >= Lz beyond the far face: end at +ZFAR + LZ
	const double Z_END = ZFAR + LZ;
	const int N_STEPS = (int)std::ceil((Z_END - LAUNCH_Z)/(V*DT));
	const int SAVE_EVERY = env_i("VD_SAVE_EVERY",0)>0 ? env_i("VD_SAVE_EVERY",0)
	                       : std::max(1,(int)std::round(N_STEPS/300.0));
	const std::string OUT="results";

	std::cout<<std::setprecision(10)<<"\n=== vac_dynamic  "<<LX<<"x"<<LY<<"x"<<LZ
	         <<" periodicity(2)  v="<<V<<"  launch_z="<<LAUNCH_Z
	         <<"  far_face=+"<<ZFAR<<"  z_end=+"<<Z_END
	         <<"  N_STEPS="<<N_STEPS<<"  save_every="<<SAVE_EVERY<<" ===\n";

	auto cell=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b).periodicity(2);
	auto ions=systems::ions(cell);
	auto electrons=systems::electrons(ions,
		options::electrons{}.spacing(SPACING*1.0_b).extra_electrons(2),
		input::kpoints::gamma());
	ground_state::initial_guess(ions, electrons);

	inqkit::dynamics::Projectile proj(MASS, -1.0,
		inqkit::detail::Vec3{0.0,0.0,LAUNCH_Z}, inqkit::detail::Vec3{0.0,0.0,V});
	inqkit::dynamics::moving_gaussian_projectile_perturbation proj_pert(proj, SIGMA_POT);
	auto basis=electrons.density().basis();

	std::filesystem::create_directories(OUT+"/frames/phi");
	std::ofstream csv;
	if(electrons.root()){ csv.open(OUT+"/phi_time.csv");
		csv<<std::setprecision(10)<<"step,time_au,proj_z,proj_vz,phi_peak\n"; }

	inqkit::RealTimeSession rt(ions,electrons,1);
	rt.add([&](inqkit::StepContext const& ctx){
		auto Rn=proj.R();
		inq::vector3<double> center{Rn.x,Rn.y,Rn.z};
		if(ctx.step % SAVE_EVERY == 0){
			auto nproj=inqkit::jellium::gaussian_density(basis, center, SIGMA_POT);
			auto phi=solvers::poisson::solve(nproj);
			auto phys=field_to_physical(phi);
			double phi_peak=0.0; for(double v:phys.values) phi_peak=std::max(phi_peak,v);
			inqkit::io::RealField3DWriter wr(OUT+"/frames/phi",
				{.field_name="phi_proj",.include_meta=false,.emit_raw=false,.emit_vti=true,
				 .vti_format=inqkit::io::VTIWriteOptions::Format::binary},{.overwrite=true});
			wr.write(phys,"phi_t"+tag6(ctx.step));
			if(electrons.root()) csv<<ctx.step<<","<<ctx.time_au<<","<<Rn.z<<","<<proj.V().z<<","<<phi_peak<<"\n";
		}
		proj.advance(inqkit::detail::Vec3{0.0,0.0,0.0}, DT);   // const-v: zero force
	});

	auto step_fn=[&](auto const&data){rt.step(data);};
	auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
	real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,proj_pert,0);

	if(electrons.root()){ csv.close(); std::ofstream s(OUT+"/run_summary.txt");
		s<<std::setprecision(10)<<"run = classical_highdensity_sv/vac_dynamic\n"
		 <<"test = vacuum_dynamic_moving_perturbation\nengine = inq\nperiodicity = 2\n"
		 <<"cell_bohr = "<<LX<<"x"<<LY<<"x"<<LZ<<"\nspacing_bohr = "<<SPACING<<"\n"
		 <<"sigma_wp = "<<SIGMA_WP<<"  sigma_pot = "<<SIGMA_POT<<"\n"
		 <<"mass = "<<MASS<<"  v = "<<V<<"  launch_z = "<<LAUNCH_Z<<"\n"
		 <<"far_face = "<<ZFAR<<"  z_end = "<<Z_END<<"  (>= Lz beyond face)\n"
		 <<"dt = "<<DT<<"  n_steps = "<<N_STEPS<<"  save_every = "<<SAVE_EVERY<<"\n"
		 <<"proj_z_final = "<<proj.R().z<<"  proj_vz_final = "<<proj.V().z<<"\nrun_completed = true\n"; }
	std::cout<<"  done  proj_z_final="<<proj.R().z<<"\n"; return 0;
}
