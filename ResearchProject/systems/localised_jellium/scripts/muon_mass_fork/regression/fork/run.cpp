// ============================================================================
// localised_jellium / scripts/muon_mass_fork / regression / run.cpp
//
// PHASE-3 bit-for-bit electron regression for the per-state mass fork.
//
// A small closed-shell electron system (He atom, LDA) is taken through a full
// SCF ground state + a short kicked real-time propagation. The run NEVER touches
// electrons.inverse_mass(), so every orbital keeps mass 1.0: the fork's
// empty-factor guard (deviation Sum|im-1| == 0) MUST route the ORIGINAL scalar
// kinetic path in both the GS (calculator.hpp) and RT (propagate.hpp).
//
// Built TWICE -- against inq-study (the fork) and against pristine inq -- the two
// binaries must produce IDENTICAL E_total/E_kinetic/E_hartree/E_xc (GS + every RT
// step) and density to machine precision. Any difference means the fork perturbs
// the mass-1 path => the fork is BROKEN (hard trust gate; stop, report no muon
// physics). compare_regression.py does the diff.
//
// Env: REG_OUT(reg)  REG_L(10.0)  REG_CUTOFF_RY(30)  REG_DT(0.02)
//      REG_TSTEPS(40)  REG_KICK(0.01)
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

int main() {
    const std::string OUT   = "results/" + env_s("REG_OUT", "reg");
    const double L          = env_d("REG_L", 10.0);
    const double CUTOFF_RY  = env_d("REG_CUTOFF_RY", 30.0);
    const double DT_AU      = env_d("REG_DT", 0.02);
    const int    N_STEPS    = env_i("REG_TSTEPS", 40);
    const double KICK       = env_d("REG_KICK", 0.01);

    std::filesystem::create_directories(OUT + "/raw/observables");

    auto cell = systems::cell::cubic(L * 1.0_b).finite();
    systems::ions ions(cell);
    ions.insert("He", {L/2 * 1.0_b, L/2 * 1.0_b, L/2 * 1.0_b});

    auto electrons = systems::electrons(
        ions, options::electrons{}.cutoff(CUTOFF_RY * 1.0_Ry), input::kpoints::gamma());
    // NOTE: inverse_mass() is deliberately LEFT at its 1.0 default -> scalar path.
    ground_state::initial_guess(ions, electrons);

    auto gs = ground_state::calculate(
        ions, electrons, options::theory{}.lda(),
        options::ground_state{}.energy_tolerance(1e-8_Ha).max_steps(400).broyden_mixing());

    if (electrons.root()) {
        std::ofstream g(OUT + "/raw/observables/gs_energy.csv");
        g << std::setprecision(14) << "quantity,value_ha\n"
          << "total,"    << gs.energy.total()    << "\n"
          << "kinetic,"  << gs.energy.kinetic()  << "\n"
          << "hartree,"  << gs.energy.hartree()  << "\n"
          << "external," << gs.energy.external() << "\n"
          << "xc,"       << gs.energy.xc()       << "\n";
        std::cout << "  GS total = " << std::setprecision(14) << gs.energy.total() << " Ha\n";
    }

    // GS density VTI (physical order; compare_regression.py diffs the raw array).
    {
        inqkit::io::RealField3DLayout lay{};
        lay.field_name = "density"; lay.emit_raw = true; lay.emit_vti = false;
        inqkit::io::RealField3DWriter(OUT + "/raw/observables/gs_density", lay,
            inqkit::io::RealField3DWriteOptions{.overwrite = true})
            .write(inqkit::fields::density::total(electrons), "gs_density");
    }

    // ----- kicked RT: energy trace each step ---------------------------------
    std::ofstream rt;
    if (electrons.root()) {
        rt.open(OUT + "/raw/observables/rt_energy.csv");
        rt << std::setprecision(14) << "step,time_au,total,kinetic,hartree,external,xc\n";
    }
    auto step_fn = [&](auto const& data){
        if (data.root()) {
            auto const& e = data.energy();
            rt << data.iter() << "," << (data.iter()*DT_AU) << ","
               << e.total() << "," << e.kinetic() << "," << e.hartree() << ","
               << e.external() << "," << e.xc() << "\n";
        }
    };

    perturbations::kick kick(cell, {0.0, 0.0, KICK});
    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime);
    real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, kick);

    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << "run = localised_jellium/muon_mass_fork/regression/" << env_s("REG_OUT","reg") << "\n"
          << "system = He atom LDA, L=" << L << " cutoff=" << CUTOFF_RY << " Ry\n"
          << "gs_total_ha = " << std::setprecision(14) << gs.energy.total() << "\n"
          << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  kick_z = " << KICK << "\n"
          << "inverse_mass = untouched (all 1.0 -> scalar path)\n"
          << "run_completed = true\n";
        std::cout << "  regression run done.\n";
    }
    return 0;
}
