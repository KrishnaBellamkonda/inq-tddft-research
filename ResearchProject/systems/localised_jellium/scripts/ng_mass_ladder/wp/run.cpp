// ============================================================================
// systems/localised_jellium/scripts/ng_mass_ladder/wp/run.cpp
//
// WAVEPACKET half of the Nazarov-Gross mass ladder.
// Plan: docs/plans/nazarov-gross-slab-mass-ladder.md
//
// A wide Gaussian electron wavepacket, injected as an extra KS orbital with a
// PER-STATE mass (the inq-study inverse_mass fork), fired at fixed VELOCITY
// through a dense (r_s = 2.5) jellium slab with two-sided CAPs.
//
// THE SWEEP INVARIANT IS VELOCITY, NOT ENERGY. Nazarov-Gross compare projectiles
// of the same charge at the same v but different M, so every rung runs at
// v0 = Cfg::V0_AU = 1.40 v_F and carries k0 = M*v0, E = M*v0^2/2.
//
// dt = 0.08 * min(M,1) * h^2. The min() is not a typo: one dt advances ALL 124
// orbitals and the 103 bath states have m = 1, so a heavy projectile cannot
// loosen the ceiling — only a light one can tighten it.
//
// ONLY THE WP'S MASS CHANGES. electrons.inverse_mass()[0][wp_idx] touches the
// single extra state; the 103 bath orbitals keep m = 1. save/load does NOT
// persist it, so it is re-applied on resume (validated in muon_mass_fork).
//
// Emits, every callback: observables.csv (INQ scalars) + interactions.csv (the
// P/S/B pairwise decomposition, rule decomposed-interaction-energies.md) +
// wp_momentum_stats / wp_real_space_stats (the WIDTH sigma(t), which is the
// mechanism this campaign tests) + density frames for the notebook GIF.
//
// Env: NG_MASS(1.0) NG_V(Cfg::V0_AU) NG_SIGMA_WP(Cfg::WP_SIGMA_BOHR)
//      NG_SPACING(0.50) NG_GS_DIR(REQUIRED) NG_OUT(wp_m1) NG_N_STEPS
//      NG_DT(0.02*M) NG_LAUNCH_Z(-25) NG_WRITE_EVERY(10) NG_WF_EVERY(50)
//      NG_SAVE_EVERY(50) NG_CKPT_EVERY(200) NG_RESUME(0) NG_CAP(1)
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/interaction_energies.hpp>

#include "../../../shared/configs/slab_n206_L30x30x120_rs2p5.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN206_L30x30x120_rs2p5;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

static int read_last_step(const std::string& path) {
    std::ifstream f(path); if (!f) return -1;
    std::string line; int last = -1;
    while (std::getline(f, line)) if (line.rfind("last_step=",0)==0) last = std::atoi(line.c_str()+10);
    return last;
}
static int read_kv_int(const std::string& path, const std::string& key, int dflt) {
    std::ifstream f(path); if (!f) return dflt;
    std::string line;
    while (std::getline(f, line)) if (line.rfind(key+"=",0)==0) return std::atoi(line.c_str()+key.size()+1);
    return dflt;
}
static std::string tag6(int n){ std::ostringstream o; o<<std::setw(6)<<std::setfill('0')<<n; return o.str(); }

int main() {
    auto t0 = std::chrono::steady_clock::now();
    const double HA = 27.211386245988;

    const double MASS        = env_d("NG_MASS", 1.0);
    const double VEL         = env_d("NG_V", Cfg::V0_AU);       // fixed-velocity invariant
    const double SIGMA_WP    = env_d("NG_SIGMA_WP", Cfg::WP_SIGMA_BOHR);
    const double KZ          = MASS * VEL;                       // k0 = m*v
    const double INV_MASS    = 1.0 / MASS;
    const double SPACING     = env_d("NG_SPACING", Cfg::SPACING_BOHR);
    const std::string TAG    = env_s("NG_OUT", "wp_m1");
    const std::string OUT    = "results/" + TAG;
    const double DT_AU       = env_d("NG_DT", 0.08 * std::min(MASS, 1.0) * SPACING * SPACING);
    const int    N_STEPS     = env_i("NG_N_STEPS", 2560);
    const int    WRITE_EVERY = env_i("NG_WRITE_EVERY", 10);
    // DISK CADENCES, sized from measured artefact sizes on THIS grid
    // (60 x 60 x 240 = 864k points, 123 states), user instruction 2026-08-05:
    //   density VTI      864k x 8 B   =  6.9 MB   -> ~24 frames = 170 MB per field
    //   wavefunction VTI 864k x 16 B  = 13.8 MB   -> 4 frames   =  55 MB
    //   RT checkpoint    123 x 864k x 16 B = 1.7 GB  (measured: an equivalent
    //     dx=0.5 slab GS with 74 states is 975 MB -> 16 B/state/point)
    // Interior checkpoints ROLL INTO ONE DIRECTORY, so N of them cost 1.7 GB, not
    // N x 1.7 GB. CKPT_EVERY defaults to N_STEPS/3 => 2 interior writes + the
    // mandatory FINAL one (rule final-timestep-checkpoint.md), which is the
    // minimum that still bounds a crash to ~1/3 of a run.
    const int    WF_EVERY    = env_i("NG_WF_EVERY",   std::max(1, N_STEPS/4));
    const int    SAVE_EVERY  = env_i("NG_SAVE_EVERY", std::max(1, N_STEPS/24));
    const double LAUNCH_Z    = env_d("NG_LAUNCH_Z", Cfg::WP_CZ_BOHR);
    const int    CKPT_EVERY  = env_i("NG_CKPT_EVERY", std::max(1, N_STEPS/3));
    const bool   RESUME      = env_i("NG_RESUME", 0) != 0;
    const bool   USE_CAP     = env_i("NG_CAP", 1) != 0;
    const double E_EV        = 0.5 * MASS * VEL * VEL * HA;
    const std::string RT_CKPT_DIR  = OUT + "/rt_ckpt";
    const std::string RT_STATE_TXT = RT_CKPT_DIR + "/rt_state.txt";

    const std::string GS_DIR = env_s("NG_GS_DIR", "");
    if (GS_DIR.empty() || !std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: NG_GS_DIR missing or unset: '" << GS_DIR << "'\n"; return 2; }

    // ---- aliasing guard (plan section 4.1). A packet whose k0 + 3 sigma_k
    //      exceeds the grid Nyquist folds back and the run is silently garbage.
    const double k_need = KZ + 3.0/(2.0*SIGMA_WP);
    const double k_max  = M_PI / SPACING;
    if (k_need > k_max) {
        std::cerr << "FATAL: aliasing. k0+3sigma_k = " << k_need << " > k_max = " << k_max
                  << " (M=" << MASS << ", v=" << VEL << ", sigma_WP=" << SIGMA_WP
                  << ", h=" << SPACING << "). Lower M*v, widen the packet, or refine h.\n";
        return 2;
    }
    // ---- timestep guard (plan section 4.2), calibrated on p3 and sigma1_masspair.
    //
    // THE CEILING IS SET BY THE LIGHTEST STATE IN THE BOX, NOT BY THE PROJECTILE.
    // The propagator advances all 124 orbitals with ONE dt, and the 103 bath
    // orbitals have m = 1. So min(M, 1) — a HEAVY wavepacket buys nothing (its own
    // kinetic operator is gentler, but the bath's is not), while a LIGHT one
    // tightens the ceiling for everybody. Cost therefore scales as 1/min(M,1):
    // flat for M >= 1, and 1/M below it.
    const double dt_ceiling = 0.08 * std::min(MASS, 1.0) * SPACING * SPACING;
    if (DT_AU > dt_ceiling * 1.001) {
        std::cerr << "FATAL: dt = " << DT_AU << " exceeds the calibrated ceiling "
                  << dt_ceiling << " = 0.08*M*h^2 (M=" << MASS << ", h=" << SPACING << ").\n";
        return 2;
    }

    int START = 0;
    if (RESUME) {
        if (!std::filesystem::exists(RT_CKPT_DIR)) { std::cerr << "FATAL: resume but no RT ckpt: " << RT_CKPT_DIR << "\n"; return 2; }
        START = read_last_step(RT_STATE_TXT);
        if (START < 0) { std::cerr << "FATAL: resume but unreadable " << RT_STATE_TXT << "\n"; return 2; }
        if (START >= N_STEPS) { std::cout << "Already at/after target (" << START << ">=" << N_STEPS << "); nothing to do.\n"; return 0; }
    }
    const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

    std::cout << std::setprecision(12)
              << "\n=== ng_mass_ladder wp (out=" << OUT << ") ===\n"
              << "  mass=" << MASS << " v=" << VEL << " (" << VEL/Cfg::KF_AU << " v_F, "
              << VEL/Cfg::V_BRAGG_AU << " of Bragg peak) k0=" << KZ << " E=" << E_EV << " eV\n"
              << "  sigma_WP=" << SIGMA_WP << "  h=" << SPACING << "  dt=" << DT_AU
              << " (ceiling " << dt_ceiling << ")\n"
              << "  k0+3sk=" << k_need << " <= k_max=" << k_max << "  [aliasing OK]\n"
              << "  N_STEPS=" << N_STEPS << " launch_z=" << LAUNCH_Z
              << "  cap=" << (USE_CAP?"on":"off") << "\n";

    auto cell0 = systems::cell::orthorhombic(Cfg::LX_BOHR*1.0_b, Cfg::LY_BOHR*1.0_b, Cfg::LZ_BOHR*1.0_b);
    auto cell  = (Cfg::PERIODICITY == 2) ? cell0.periodicity(2) : cell0.periodic();
    auto ions  = systems::ions(cell);
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());

    int wp_idx_pre = -1;
    if (RESUME) {
        electrons.load(RT_CKPT_DIR);
        wp_idx_pre = read_kv_int(RT_STATE_TXT, "wp_idx", electrons.states().num_states() - 1);
        std::cout << "  [RESUME] loaded RT ckpt at step " << START << "\n";
    } else {
        electrons.load(GS_DIR);
        std::cout << "  Loaded GS from " << GS_DIR << "\n";
    }
    const int n_states = electrons.states().num_states();

    std::filesystem::create_directories(OUT + "/raw/observables");
    for (auto sub : {"density_total","density_wp","wavefunction_wp"})
        std::filesystem::create_directories(OUT + "/raw/vti/" + sub);
    std::filesystem::create_directories(RT_CKPT_DIR);

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};

    // ----- WP injection + per-state mass (the inq-study fork) ------------
    int wp_idx = wp_idx_pre;
    double wp_norm_after = 1.0, wp_max_overlap = 0.0;
    if (!RESUME) {
        auto wp = inqkit::WavePacket{}
                      .center(0.0, 0.0, LAUNCH_Z).sigma(SIGMA_WP)
                      .k0(0.0, 0.0, KZ).orthogonalise_against_occupied(electrons);
        auto report = wp.inject_into_last_extra_state(electrons, 1.0);
        wp_idx = report.state_index;
        wp_norm_after = report.norm_after;
        wp_max_overlap = report.max_overlap;
        std::cout << "  WP injected: idx=" << wp_idx << " norm_after=" << report.norm_after
                  << " max_overlap=" << report.max_overlap << "\n";
    }
    electrons.inverse_mass()[0][wp_idx] = INV_MASS;   // the mass fork (re-applied on resume)
    std::cout << "  inverse_mass[" << wp_idx << "]=" << INV_MASS << "  (bath states keep m=1)\n";

    // ----- background well + CAPs ----------------------------------------
    inqkit::jellium::localised_background_params bg;
    bg.shape = inqkit::jellium::background_shape::slab;
    bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
    bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    perturbations::absorbing cap_lo(Cfg::CAP_ETA_HA * 1.0_Ha, -Cfg::CAP_MID_FRAC, Cfg::CAP_WIDTH_FRAC);
    perturbations::absorbing cap_hi(Cfg::CAP_ETA_HA * 1.0_Ha,  Cfg::CAP_MID_FRAC, Cfg::CAP_WIDTH_FRAC);

    // constant background fields for the pairwise decomposition
    auto basis   = electrons.density().basis();
    auto nplus   = bg_pert.background_density(basis);
    auto phiplus = solvers::poisson::solve(nplus);
    const double E_BB = inqkit::jellium::background_self_energy(nplus, phiplus);

    // ----- writers --------------------------------------------------------
    // FULL energy decomposition — every INQ scalar the writer exposes. These are
    // the GROSS ledger; interactions.csv below carries the representation-
    // independent P/S/B pairwise terms that are actually comparable between the
    // classical and wavepacket halves (rule decomposed-interaction-energies.md:
    // for a WP, energy_external is identically 0 and energy_hartree silently
    // contains E_SS + E_PS + E_PP, so a raw scalar comparison is meaningless).
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.energy_external = sel.energy_nonlocal = sel.energy_ion = true;
    sel.energy_ion_kinetic = sel.energy_exact_exchange = true;
    sel.energy_nvxc = sel.energy_eigenvalues = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
    sel.cod_x = sel.cod_y = sel.cod_z = true;
    inqkit::io::ObservablesWriter obs(OUT + "/raw/observables/observables" + SEG + ".csv", sel);
    obs.write_header();

    std::ofstream ix;
    if (electrons.root()) {
        ix.open(OUT + "/raw/observables/interactions" + SEG + ".csv");
        ix << std::setprecision(12)
           << "step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,e_hartree_check,e_external_check,"
              "norm_wp,norm_total\n";
    }
    std::ofstream nlog;
    if (electrons.root()) {
        nlog.open(OUT + "/raw/observables/electron_number" + SEG + ".csv");
        nlog << std::setprecision(12) << "step,time_au,N_total\n";
    }

    inqkit::observables::WPMomentumStats wp_momentum_stats(
        OUT + "/raw/observables/wp_momentum_stats" + SEG + ".csv", wp_idx, {.write_every=WRITE_EVERY});
    inqkit::observables::WPRealSpaceStats wp_real_space_stats(
        OUT + "/raw/observables/wp_real_space_stats" + SEG + ".csv", wp_idx, {.write_every=WRITE_EVERY});

    inqkit::io::RealField3DWriter total_wr(OUT + "/raw/vti/density_total", vti_layout, {.overwrite = !RESUME});
    inqkit::io::RealField3DWriter wp_wr   (OUT + "/raw/vti/density_wp",    vti_layout, {.overwrite = !RESUME});
    inqkit::io::ComplexField3DWriter wf_wr(
        OUT + "/raw/vti/wavefunction_wp",
        {.field_name="wavefunction", .include_meta=false, .emit_raw=false,
         .emit_vti=true, .vti_format=inqkit::io::VTIWriteOptions::Format::binary},
        {.overwrite = !RESUME});

    // ----- real-time session ---------------------------------------------
    inqkit::RealTimeSession rt(ions, electrons, 1);
    rt.add([&](inqkit::StepContext const& ctx) {
        obs.append(ctx);

        // Pairwise P/S/B decomposition EVERY step (rule decomposed-interaction-
        // energies.md: two Poisson solves are negligible against the propagator).
        auto n_total = ctx.electrons->density();
        auto n_wp    = inqkit::jellium::orbital_density_field(*ctx.electrons, wp_idx);
        auto ct      = inqkit::jellium::compute_coulomb_wp(n_total, n_wp, phiplus);
        if (ctx.electrons->root()) {
            ix << ctx.step << "," << ctx.time_au << "," << ct.e_ss << "," << ct.e_pp << ","
               << ct.e_ps << "," << ct.e_sb << "," << ct.e_pb << "," << E_BB << ","
               << ct.e_hartree_check << "," << ct.e_external_check << ","
               << ct.norm_wp << "," << ct.norm_total << "\n";
            nlog << ctx.step << "," << ctx.time_au << "," << ct.norm_total << "\n";
        }

        if (SAVE_EVERY > 0 && ctx.step % SAVE_EVERY == 0) {
            total_wr.write(inqkit::fields::density::total(*ctx.electrons), ctx.time_au, ctx.step);
            wp_wr.write(inqkit::fields::density::orbital(*ctx.electrons, wp_idx), ctx.time_au, ctx.step);
        }
        if (WF_EVERY > 0 && ctx.step % WF_EVERY == 0)
            wf_wr.write(inqkit::fields::orbital::wavefunction(*ctx.electrons, wp_idx),
                        "wavefunction_t" + tag6(ctx.step));
    });

    auto step_fn = [&](auto const& data) {
        rt.step(data);
        const int it = data.iter();
        wp_momentum_stats.maybe_accumulate(data);
        wp_real_space_stats.maybe_accumulate(data);
        if (it > START && it % CKPT_EVERY == 0 && it < N_STEPS) {
            electrons.save(RT_CKPT_DIR);
            if (data.root()) {
                std::ofstream st(RT_STATE_TXT, std::ios::trunc);
                st << "last_step=" << it << "\ntime_au=" << (it*DT_AU) << "\nwp_idx=" << wp_idx
                   << "\nn_states=" << n_states << "\nn_steps_target=" << N_STEPS
                   << "\ndt=" << DT_AU << "\ninv_mass=" << INV_MASS << "\n";
                std::cout << "  [ckpt] step " << it << " (t=" << it*DT_AU << ")\n";
            }
        }
    };

    auto pert    = perturbations::sum(bg_pert, perturbations::sum(cap_lo, cap_hi));
    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime)
                       .observables_current().observables_dipole();
    if (USE_CAP) real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, pert,    START);
    else         real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, bg_pert, START);

    // FINAL checkpoint (rule final-timestep-checkpoint.md)
    electrons.save(RT_CKPT_DIR);
    if (electrons.root()) {
        std::ofstream st(RT_STATE_TXT, std::ios::trunc);
        st << "last_step=" << N_STEPS << "\ntime_au=" << N_STEPS*DT_AU << "\nwp_idx=" << wp_idx
           << "\nn_states=" << n_states << "\nn_steps_target=" << N_STEPS
           << "\ndt=" << DT_AU << "\ninv_mass=" << INV_MASS << "\n";
    }

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (electrons.root()) {
        ix.close(); nlog.close();
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = localised_jellium/ng_mass_ladder/" << TAG << "\n"
          << "plan = docs/plans/nazarov-gross-slab-mass-ladder.md\n"
          << "engine = inq-study (per-state mass fork)\nxc = LDA (no SIC)\n"
          << "representation = wavepacket_orbital\n"
          << "projectile = wavepacket sigma_WP " << SIGMA_WP << " mass " << MASS
          << " velocity " << VEL << " k0 " << KZ << " E " << E_EV << " eV\n"
          << "inverse_mass = " << INV_MASS << "  wp_state_index = " << wp_idx << "\n"
          << "v_over_vF = " << VEL/Cfg::KF_AU << "  v_over_vBragg = " << VEL/Cfg::V_BRAGG_AU << "\n"
          << "cap = " << (USE_CAP?"on":"off") << "  cap_eta_ha = " << Cfg::CAP_ETA_HA
          << "  cap_region_bohr = +/-[" << Cfg::CAP_INNER_FACE_BOHR << "," << Cfg::LZ_BOHR/2.0 << "]\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR
          << "  periodicity = " << Cfg::PERIODICITY << "  spacing = " << SPACING << "\n"
          << "background = slab half_width " << Cfg::SLAB_HALF_WIDTH << " axis " << Cfg::SLAB_AXIS
          << " edge " << Cfg::EDGE_WIDTH_BOHR << "\n"
          << "n0_a0m3 = " << Cfg::N0 << "  r_s = " << Cfg::RS_BOHR
          << "  kf_au = " << Cfg::KF_AU << "  omega_p_au = " << Cfg::OMEGA_P_AU << "\n"
          << "n_electrons = " << Cfg::N_ELECTRONS << "  n_states = " << n_states << "\n"
          << "wp_norm_after = " << wp_norm_after << "  wp_max_overlap = " << wp_max_overlap
          << "  launch_z = " << LAUNCH_Z << "\n"
          << "gs_dir = " << GS_DIR << "\n"
          << "dt = " << DT_AU << "  dt_ceiling = " << dt_ceiling
          << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY << "\n"
          << "aliasing_k_need = " << k_need << "  k_max = " << k_max << "\n"
          << "e_bb_ha = " << E_BB << "\n"
          << "ckpt_every = " << CKPT_EVERY << "  rt_ckpt_dir = " << RT_CKPT_DIR << "\n"
          << "resume = " << (RESUME?"true":"false") << "  start_step = " << START
          << "  segment_suffix = " << (SEG.empty()?"(none)":SEG) << "\n"
          << "wall_time_s = " << wall << " (this segment)\nrun_completed = true\n";
    }
    std::cout << "  done. wall=" << wall << "s\n";
    return 0;
}
