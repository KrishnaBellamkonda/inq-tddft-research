// ============================================================================
// localised_jellium / scripts/fullsuite_wp / run.cpp
//
// FULL-SUITE wavepacket re-run (D1/D2, 2026-06-22). Mirrors the proven gold WP
// template (run_wp_n162_*) + the localised background well + an env-gated
// two-sided sin² CAP. Build ONCE against inq-study; the GS checkpoint
// (shared_gs/slab_n234_L50, converged WITH the background) is reused and the
// background perturbation re-applied so the well persists.
//
//   LJ_CAP=0 -> Phase 3 baseline (no CAP, closed box). launch_z default -23.
//   LJ_CAP=1 -> Phase 5 WP+CAP.                        launch_z default -15.5.
//
// Emits the canonical suite: density VTIs total/system/gs_system/wp +
// wavefunction_wp, delta(+coarse); observables(energies/current/dipole/L2),
// state_energies, occupations, eigenvalues; WP: momentum_distribution,
// wp_momentum_stats, wp_real_space_stats, overlap_full (excitations) + overlap.
//
// Env: LJ_OUT(p3_wp) LJ_CAP(0) LJ_N_STEPS LJ_WRITE_EVERY(10) LJ_WF_EVERY(10)
//      LJ_LAUNCH_Z LJ_DT(0.02).
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

#include "../../../shared/configs/slab_n82_L50x50x90_E54.hpp"
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
using Cfg = localised_jellium::config::SlabN82_L50x50x90_E54;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

int main() {
    auto t0 = std::chrono::steady_clock::now();

    const bool   USE_CAP     = env_i("LJ_CAP", 1) != 0;
    const std::string OUT    = "results/" + env_s("LJ_OUT", "p3_wp");
    const double DT_AU       = env_d("LJ_DT", 0.02);
    const int    N_STEPS     = env_i("LJ_N_STEPS", 2000);
    const int    WRITE_EVERY = env_i("LJ_WRITE_EVERY", 10);
    const int    WF_EVERY    = env_i("LJ_WF_EVERY", 10);
    const double LAUNCH_Z    = env_d("LJ_LAUNCH_Z", -23.75);
    const double CAP_ETA = -0.7, CAP_MID = 40.0/90.0, CAP_WIDTH = 10.0/90.0;  // two-sided, 10 Bohr/side, region [+/-35, +/-45]

    const std::string GS_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "shared_gs/slab_n82_L50x50x90";
    if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2; }

    std::cout << "\n=== fullsuite_wp (cap=" << (USE_CAP?"on":"off") << ", out=" << OUT << ") ===\n"
              << "  N_STEPS=" << N_STEPS << " dt=" << DT_AU << " WRITE_EVERY=" << WRITE_EVERY
              << " launch_z=" << LAUNCH_Z << "\n";

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
    electrons.load(GS_DIR);
    std::cout << "  Loaded GS (inq-study) from " << GS_DIR << "\n";
    jellium::eigenvalues::copy_from_checkpoint(GS_DIR, OUT + "/raw/observables/eigenvalues");
    const int n_states = electrons.states().num_states();

    // ----- output skeleton ----------------------------------------------
    for (auto sub : {"density_total","density_system","density_gs_system","density_wp",
                     "wavefunction_wp","density_delta","density_delta_coarse"})
        std::filesystem::create_directories(OUT + "/raw/vti/" + sub);
    std::filesystem::create_directories(OUT + "/raw/observables/overlap");
    std::filesystem::create_directories(OUT + "/raw/observables/overlap_full");

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};

    { inqkit::io::RealField3DWriter gs_wr(OUT + "/raw/vti/density_gs_system", vti_layout, {.overwrite=true});
      gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system"); }

    // ----- WP injection --------------------------------------------------
    auto wp = inqkit::WavePacket{}
                  .center(0.0, 0.0, LAUNCH_Z).sigma(Cfg::WP_SIGMA_BOHR)
                  .k0(0.0, 0.0, Cfg::WP_KZ).orthogonalise_against_occupied(electrons);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report.state_index;
    std::cout << "  WP injected: idx=" << wp_idx << " norm_after=" << report.norm_after
              << " max_overlap=" << report.max_overlap << "\n";

    // ----- background well (+ optional two-sided CAP) --------------------
    inqkit::jellium::localised_background_params bg;
    bg.shape = inqkit::jellium::background_shape::slab;
    bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
    bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);
    perturbations::absorbing cap_lo(CAP_ETA * 1.0_Ha, -CAP_MID, CAP_WIDTH);
    perturbations::absorbing cap_hi(CAP_ETA * 1.0_Ha,  CAP_MID, CAP_WIDTH);
    auto pert_bg  = bg_pert;
    auto pert_cap = perturbations::sum(bg_pert, perturbations::sum(cap_lo, cap_hi));

    // ----- density writers ----------------------------------------------
    inqkit::io::RealField3DWriter total_wr (OUT + "/raw/vti/density_total",  vti_layout, {.overwrite=true});
    inqkit::io::RealField3DWriter system_wr(OUT + "/raw/vti/density_system", vti_layout, {.overwrite=true});
    { auto s0 = inqkit::fields::density::total(electrons); total_wr.write(s0,0.0,0); system_wr.write(s0,0.0,0); }
    inqkit::io::RealField3DWriter wp_density_wr(OUT + "/raw/vti/density_wp", vti_layout, {.overwrite=true});
    inqkit::io::ComplexField3DWriter wp_wf_wr(
        OUT + "/raw/vti/wavefunction_wp",
        {.field_name="wavefunction", .include_meta=false, .emit_raw=false,
         .emit_vti=true, .vti_format=inqkit::io::VTIWriteOptions::Format::binary},
        {.overwrite=true});

    // ----- scalar/observable writers ------------------------------------
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables.csv", sel);
    obs_writer.write_header();
    inqkit::observables::StateEnergyWriter state_energy_wr(OUT + "/raw/observables/state_energies.csv", true);
    inqkit::observables::OccupationsWriter occupations_wr(OUT + "/raw/observables/occupations_vs_time.csv");
    inqkit::observables::DensityDelta density_delta(
        OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
        {.emit_raw_vti=true, .emit_coarse_vti=true, .compute_l2=true, .coarse_bin_bohr=3.0});

    // ----- WP-specific observables --------------------------------------
    inqkit::observables::OrbitalOverlapMatrix overlap_obs(electrons, wp_idx, OUT + "/raw/observables/overlap");
    inqkit::observables::OrbitalOverlapMatrix overlap_full_obs(electrons, n_states - 1, OUT + "/raw/observables/overlap_full");
    overlap_full_obs.snapshot(electrons, 0.0, 0);
    overlap_obs.snapshot_wp_only(electrons, 0.0, 0);
    inqkit::observables::MomentumDistribution momentum_dist(
        OUT + "/raw/observables/momentum_distribution.csv", wp_idx, Cfg::LZ_BOHR,
        {.n_bins=64, .k_max_bohr_inv=0.0, .write_every=WRITE_EVERY});
    inqkit::observables::WPMomentumStats wp_momentum_stats(
        OUT + "/raw/observables/wp_momentum_stats.csv", wp_idx, {.write_every=WRITE_EVERY});
    inqkit::observables::WPRealSpaceStats wp_real_space_stats(
        OUT + "/raw/observables/wp_real_space_stats.csv", wp_idx, {.write_every=WRITE_EVERY});

    std::ofstream nlog(OUT + "/raw/observables/electron_number.csv");
    nlog << std::setprecision(12) << "step,time_au,N_total\n";

    // ----- real-time density/observable session -------------------------
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

    auto step_fn = [&](auto const& data) {
        rt_obs.step(data);
        const int it = data.iter();
        if (it % (5 * WRITE_EVERY) == 0) { state_energy_wr.snapshot(data); occupations_wr.snapshot(data); }
        momentum_dist.maybe_accumulate(data);
        wp_momentum_stats.maybe_accumulate(data);
        wp_real_space_stats.maybe_accumulate(data);
        if (data.root()) nlog << it << "," << (it*DT_AU) << "," << data.num_electrons() << "\n";
    };

    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime)
                       .observables_current().observables_dipole();
    if (USE_CAP) real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, pert_cap);
    else         real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, pert_bg);

    overlap_full_obs.snapshot(electrons, DT_AU * N_STEPS, N_STEPS);

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = localised_jellium/fullsuite_wp/" << env_s("LJ_OUT","p3_wp") << "\n"
          << "engine = inq-study\n"
          << "projectile = wavepacket sigma " << Cfg::WP_SIGMA_BOHR << " E " << Cfg::WP_EKIN_EV << " eV k0 " << Cfg::WP_KZ << "\n"
          << "cap = " << (USE_CAP?"on (two-sided sin2, eta -0.7 Ha, 10 Bohr/side, region +/-35..+/-45, inner faces +/-35)":"off") << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR << "  spacing = " << Cfg::SPACING_BOHR << "\n"
          << "background = slab half_width " << Cfg::SLAB_HALF_WIDTH << " axis " << Cfg::SLAB_AXIS << "\n"
          << "n_electrons = " << Cfg::N_ELECTRONS << "  n_states = " << n_states << "  wp_state_index = " << wp_idx << "\n"
          << "wp_norm_after = " << report.norm_after << "  launch_z = " << LAUNCH_Z << "\n"
          << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY << "\n"
          << "wall_time_s = " << wall << "\nrun_completed = true\n";
    }
    std::cout << "  done. wall=" << wall << "s\n";
    return 0;
}
