// ============================================================================
// localised_jellium / scripts/sigma1_masspair / run.cpp
//
// sigma=1 mass-pair decay runs (plan: docs/plans/sigma1-masspair-decay-runs.md).
// CLEAN qsp_phase3 geometry (50x50x90, dx 0.5, N=82, two-sided sin2 CAP
// eta -0.7 region +/-35..+/-45) with ONLY the wavepacket changed:
//   sigma_WP=1, k0=4.5, launch_z=-16.5 (4 sigma from slab face),
//   mass fork m=2 (E=138 eV) / m=3 (E=92 eV) via per-state inverse_mass.
//
// Grafts onto the qsp_phase3 full observable suite:
//   * checkpoint/resume: interior ckpt every LJ_CKPT_EVERY steps + final;
//     LJ_RESUME=1 reloads, re-applies inverse_mass (NOT persisted by save),
//     writes segment-suffixed CSVs (.from<START>).
//   * pairwise Coulomb decomposition interactions.csv per step via
//     inqkit/jellium/interaction_energies.hpp (P=WP, S=slab, B=background).
//   * ledger extended with energy_external + energy_nonlocal (twin contract).
//
// Env: LJ_OUT(REQ) LJ_INV_MASS(0.5) LJ_K0(4.5) LJ_SIGMA_WP(1.0)
//      LJ_LAUNCH_Z(-16.5) LJ_N_STEPS(2500) LJ_DT(0.04) LJ_WRITE_EVERY(8)
//      LJ_WF_EVERY(8) LJ_CKPT_EVERY(200) LJ_RESUME(0).
// Build vs INQ_SOURCE=inq-study. GPU per dispatcher (CUDA_VISIBLE_DEVICES).
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/observables/state_energy_writer.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/interaction_energies.hpp>

#include "../../../shared/configs/slab_n82_L50x50x90.hpp"
#include "../../../../jellium/shared/cpp/eigenvalues_writer.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN82_L50x50x90;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }
static double read_state_d(const std::string& path, const char* key, double def){
    std::ifstream f(path); std::string line, k = std::string(key) + "=";
    while (std::getline(f, line)) { auto p = line.find(k);
        if (p != std::string::npos) return std::atof(line.substr(p + k.size()).c_str()); }
    return def;
}

int main() {
    auto t0 = std::chrono::steady_clock::now();

    const std::string OUT    = "results/" + env_s("LJ_OUT", "wp_m2_k4p5");
    const double INV_MASS    = env_d("LJ_INV_MASS", 0.5);
    const double K0          = env_d("LJ_K0", 4.5);
    const double SIGMA_WP    = env_d("LJ_SIGMA_WP", 1.0);
    const double LAUNCH_Z    = env_d("LJ_LAUNCH_Z", -16.5);
    const double DT_AU       = env_d("LJ_DT", 0.04);
    const int    N_STEPS     = env_i("LJ_N_STEPS", 2500);
    const int    WRITE_EVERY = env_i("LJ_WRITE_EVERY", 8);
    const int    WF_EVERY    = env_i("LJ_WF_EVERY", 8);
    const int    CKPT_EVERY  = env_i("LJ_CKPT_EVERY", 200);
    const bool   RESUME      = env_i("LJ_RESUME", 0) != 0;
    const double CAP_ETA = -0.7, CAP_MID = 40.0/90.0, CAP_WIDTH = 10.0/90.0;  // p3 two-sided CAP, region [+/-35, +/-45]

    const double MASS = 1.0/INV_MASS, V0 = K0*INV_MASS;
    const double EKIN_EV = 0.5*MASS*V0*V0*27.211386;

    const std::string GS_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "shared_gs/slab_n82_L50x50x90";
    if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2; }

    // ----- checkpoint / resume -------------------------------------------
    const std::string RT_CKPT_DIR = OUT + "/rt_ckpt";
    const std::string RT_STATE    = RT_CKPT_DIR + "/rt_state.txt";
    int START = 0;
    if (RESUME) {
        START = (int)read_state_d(RT_STATE, "last_step", -1);
        if (START < 0) { std::cerr << "FATAL: LJ_RESUME=1 but unreadable " << RT_STATE << "\n"; return 2; }
        if (START >= N_STEPS) { std::cout << "Already at/after target (" << START << ">=" << N_STEPS << "); nothing to do.\n"; return 0; }
    }
    const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

    std::cout << "\n=== sigma1_masspair (out=" << OUT << (RESUME?"  [RESUME]":"  [FRESH]") << ") ===\n"
              << "  sigma_WP=" << SIGMA_WP << " k0=" << K0 << " m=" << MASS << " v0=" << V0
              << " E=" << EKIN_EV << " eV  launch_z=" << LAUNCH_Z << "\n"
              << "  start=" << START << " -> " << N_STEPS << " dt=" << DT_AU
              << " ckpt_every=" << CKPT_EVERY << "\n";

    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, Cfg::LZ_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    electrons.load(RESUME ? RT_CKPT_DIR : GS_DIR);   // resume ckpt already holds the propagated WP
    const int n_states = electrons.states().num_states();

    int wp_idx; double wp_norm_after = 1.0;
    if (RESUME) {
        wp_idx = (int)read_state_d(RT_STATE, "wp_idx", -1);
        if (wp_idx < 0) { std::cerr << "FATAL: no wp_idx in " << RT_STATE << "\n"; return 2; }
        std::cout << "  [RESUME] ckpt at step " << START << "  wp_idx=" << wp_idx << "\n";
    } else {
        jellium::eigenvalues::copy_from_checkpoint(GS_DIR, OUT + "/raw/observables/eigenvalues");
        auto wp = inqkit::WavePacket{}
                      .center(0.0, 0.0, LAUNCH_Z).sigma(SIGMA_WP)
                      .k0(0.0, 0.0, K0).orthogonalise_against_occupied(electrons);
        auto report = wp.inject_into_last_extra_state(electrons, 1.0);
        wp_idx = report.state_index; wp_norm_after = report.norm_after;
        std::cout << "  WP injected: idx=" << wp_idx << " norm_after=" << report.norm_after
                  << " max_overlap=" << report.max_overlap << "\n";
    }
    electrons.inverse_mass()[0][wp_idx] = INV_MASS;   // the mass fork (save/load does NOT persist it)

    // ----- output skeleton ----------------------------------------------
    for (auto sub : {"density_total","density_system","density_gs_system","density_wp",
                     "wavefunction_wp","density_delta","density_delta_coarse"})
        std::filesystem::create_directories(OUT + "/raw/vti/" + std::string(sub));
    std::filesystem::create_directories(OUT + "/raw/observables/overlap");
    std::filesystem::create_directories(OUT + "/raw/observables/overlap_full");
    std::filesystem::create_directories(RT_CKPT_DIR);

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};

    if (!RESUME) {
        inqkit::io::RealField3DWriter gs_wr(OUT + "/raw/vti/density_gs_system", vti_layout, {.overwrite=true});
        gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system");
    }

    // ----- background well + p3 two-sided CAP ----------------------------
    inqkit::jellium::localised_background_params bg;
    bg.shape = inqkit::jellium::background_shape::slab;
    bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
    bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);
    perturbations::absorbing cap_lo(CAP_ETA * 1.0_Ha, -CAP_MID, CAP_WIDTH);
    perturbations::absorbing cap_hi(CAP_ETA * 1.0_Ha,  CAP_MID, CAP_WIDTH);
    auto pert_cap = perturbations::sum(bg_pert, perturbations::sum(cap_lo, cap_hi));

    // ----- density writers ----------------------------------------------
    inqkit::io::RealField3DWriter total_wr (OUT + "/raw/vti/density_total",  vti_layout, {.overwrite = !RESUME});
    inqkit::io::RealField3DWriter system_wr(OUT + "/raw/vti/density_system", vti_layout, {.overwrite = !RESUME});
    if (!RESUME) { auto s0 = inqkit::fields::density::total(electrons); total_wr.write(s0,0.0,0); system_wr.write(s0,0.0,0); }
    inqkit::io::RealField3DWriter wp_density_wr(OUT + "/raw/vti/density_wp", vti_layout, {.overwrite = !RESUME});
    inqkit::io::ComplexField3DWriter wp_wf_wr(
        OUT + "/raw/vti/wavefunction_wp",
        {.field_name="wavefunction", .include_meta=false, .emit_raw=false,
         .emit_vti=true, .vti_format=inqkit::io::VTIWriteOptions::Format::binary},
        {.overwrite = !RESUME});

    // ----- scalar/observable writers (twin contract: external+nonlocal on)
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.energy_external = sel.energy_nonlocal = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables" + SEG + ".csv", sel);
    obs_writer.write_header();
    inqkit::observables::StateEnergyWriter state_energy_wr(OUT + "/raw/observables/state_energies" + SEG + ".csv", true);
    inqkit::observables::OccupationsWriter occupations_wr(OUT + "/raw/observables/occupations_vs_time" + SEG + ".csv");
    inqkit::observables::DensityDelta density_delta(
        OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
        {.emit_raw_vti=true, .emit_coarse_vti=true, .compute_l2=true, .coarse_bin_bohr=3.0});

    // ----- pairwise interaction decomposition (P=WP, S=slab, B=background)
    auto ie_basis   = electrons.density().basis();
    auto ie_nplus   = bg_pert.background_density(ie_basis);
    auto ie_phiplus = solvers::poisson::solve(ie_nplus);
    const double E_BB = inqkit::jellium::background_self_energy(ie_nplus, ie_phiplus);
    std::ofstream ix;   // Hartree; e_ss+e_ps+e_pp == E_hartree, e_sb+e_pb == E_external
    if (electrons.root()) { ix.open(OUT + "/raw/observables/interactions" + SEG + ".csv");
        ix << std::setprecision(12) << "step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,"
              "e_hartree_check,e_external_check,norm_wp,norm_total\n"; }

    // ----- WP-specific observables --------------------------------------
    inqkit::observables::OrbitalOverlapMatrix overlap_obs(electrons, wp_idx, OUT + "/raw/observables/overlap");
    inqkit::observables::OrbitalOverlapMatrix overlap_full_obs(electrons, n_states - 1, OUT + "/raw/observables/overlap_full");
    if (!RESUME) {
        overlap_full_obs.snapshot(electrons, 0.0, 0);
        overlap_obs.snapshot_wp_only(electrons, 0.0, 0);
    }
    inqkit::observables::MomentumDistribution momentum_dist(
        OUT + "/raw/observables/momentum_distribution" + SEG + ".csv", wp_idx, Cfg::LZ_BOHR,
        {.n_bins=64, .k_max_bohr_inv=0.0, .write_every=WRITE_EVERY});
    inqkit::observables::WPMomentumStats wp_momentum_stats(
        OUT + "/raw/observables/wp_momentum_stats" + SEG + ".csv", wp_idx, {.write_every=WRITE_EVERY});
    inqkit::observables::WPRealSpaceStats wp_real_space_stats(
        OUT + "/raw/observables/wp_real_space_stats" + SEG + ".csv", wp_idx, {.write_every=WRITE_EVERY});

    std::ofstream nlog(OUT + "/raw/observables/electron_number" + SEG + ".csv");
    nlog << std::setprecision(12) << "step,time_au,N_total\n";

    // ----- real-time session ---------------------------------------------
    inqkit::RealTimeSession rt_obs(ions, electrons, WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const& ctx) {
        auto sys_f = inqkit::fields::density::total(*ctx.electrons);
        system_wr.write(sys_f, ctx.time_au, ctx.step);
        total_wr.write (sys_f, ctx.time_au, ctx.step);
        const double l2 = density_delta.snapshot(sys_f, ctx.time_au, ctx.step);
        inqkit::StepContext c = ctx; c.density_l2 = l2; obs_writer.append(c);
        if (ctx.step % WF_EVERY == 0) {
            wp_density_wr.write(inqkit::fields::density::orbital(*ctx.electrons, wp_idx), ctx.time_au, ctx.step);
            char nm[64]; std::snprintf(nm, sizeof(nm), "wavefunction_t%06d", ctx.step);
            wp_wf_wr.write(inqkit::fields::orbital::wavefunction(*ctx.electrons, wp_idx), std::string(nm));
        }
        if (ctx.step % 10 == 0) overlap_obs.snapshot_wp_only(*ctx.electrons, ctx.time_au, ctx.step);
    });

    auto write_state = [&](int step){
        std::ofstream st(RT_STATE);
        st << std::setprecision(12) << "last_step=" << step << "\ntime_au=" << (step*DT_AU)
           << "\ndt=" << DT_AU << "\nwp_idx=" << wp_idx << "\ninv_mass=" << INV_MASS << "\n";
    };

    auto step_fn = [&](auto const& data) {
        rt_obs.step(data);
        const int it = data.iter();
        if (it % (5 * WRITE_EVERY) == 0) { state_energy_wr.snapshot(data); occupations_wr.snapshot(data); }
        momentum_dist.maybe_accumulate(data);
        wp_momentum_stats.maybe_accumulate(data);
        wp_real_space_stats.maybe_accumulate(data);
        { // pairwise decomposition every step (2 Poisson solves, ~ms on this grid)
            auto n_wp = inqkit::jellium::orbital_density_field(electrons, wp_idx);
            auto ct   = inqkit::jellium::compute_coulomb_wp(electrons.density(), n_wp, ie_phiplus);
            if (data.root())
                ix << it << "," << (it*DT_AU) << "," << ct.e_ss << "," << ct.e_pp << "," << ct.e_ps << ","
                   << ct.e_sb << "," << ct.e_pb << "," << E_BB << "," << ct.e_hartree_check << ","
                   << ct.e_external_check << "," << ct.norm_wp << "," << ct.norm_total << "\n";
        }
        if (data.root()) nlog << it << "," << (it*DT_AU) << "," << data.num_electrons() << "\n";
        if (it > START && it % CKPT_EVERY == 0 && it < N_STEPS) {   // interior checkpoint
            electrons.save(RT_CKPT_DIR);
            if (data.root()) write_state(it);
        }
    };

    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime)
                       .observables_current().observables_dipole();
    real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, pert_cap, START);

    overlap_full_obs.snapshot(electrons, DT_AU * N_STEPS, N_STEPS);

    // final checkpoint (extendable: LJ_RESUME=1 + larger LJ_N_STEPS continues)
    electrons.save(RT_CKPT_DIR);
    if (electrons.root()) write_state(N_STEPS);
    if (electrons.root()) ix.close();

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = localised_jellium/sigma1_masspair/" << env_s("LJ_OUT","wp_m2_k4p5") << "\n"
          << "engine = inq-study\n"
          << "projectile = wavepacket sigma " << SIGMA_WP << " k0 " << K0
          << "  inverse_mass = " << INV_MASS << "  m_eff = " << MASS << "\n"
          << "velocity_au = " << V0 << "  E_eV = " << EKIN_EV << "\n"
          << "cap = on (two-sided sin2, eta -0.7 Ha, 10 Bohr/side, region +/-35..+/-45)\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR << "  spacing = " << Cfg::SPACING_BOHR << "\n"
          << "background = slab half_width " << Cfg::SLAB_HALF_WIDTH << " axis " << Cfg::SLAB_AXIS << "\n"
          << "n_electrons = " << Cfg::N_ELECTRONS << "  n_states = " << n_states << "  wp_state_index = " << wp_idx << "\n"
          << "wp_norm_after = " << wp_norm_after << "  launch_z = " << LAUNCH_Z << "\n"
          << "start_step = " << START << "  n_steps = " << N_STEPS << "  dt_au = " << DT_AU
          << "  write_every = " << WRITE_EVERY << "  ckpt_every = " << CKPT_EVERY << "\n"
          << "gs_dir = " << GS_DIR << "\n"
          << "wall_time_s = " << wall << "\nrun_completed = true\n";
    }
    std::cout << "  done. wall=" << wall << "s\n";
    return 0;
}
