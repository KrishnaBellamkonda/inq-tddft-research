// ============================================================================
// localised_jellium / scripts/slab_ks_wrap / wp / run.cpp
//
// KS-ORBITAL STOPPING POWER ON THE JELLIUM SLAB, CAP-FREE, WITH WRAP-AROUND.
// Plan: docs/plans/slab-ks-orbital-stopping-wrap.md
//
// Forked from scripts/wp_highdensity_sv/wp/run.cpp. Same slab, same box, same
// ground states, same launch point, same dt. THREE things change:
//
//   1. NO CAP (LJ_CAP_ETA defaults to 0). The packet is no longer absorbed at
//      the z faces: it crosses the slab, leaves at +z, RE-ENTERS at -z, and
//      crosses again — about 14 slab crossings over 362 Bohr of path.
//   2. sigma_WP = 2.0 (was 0.5), matching the bulk KS-stopping study exactly so
//      the slab-vs-bulk comparison is like-for-like.
//   3. The density is an ENV parameter (LJ_N), so ONE binary serves both the
//      r_s = 4.18 (N=100) and r_s = 5.67 (N=40) campaigns. Everything else is
//      compile-time pinned to the shared Cfg, so density is the only variable.
//
// ---------------------------------------------------------------------------
// WHY REMOVING THE CAP IS THE POINT, NOT A SHORTCUT
// The CAP'd campaign could only fit a stopping power over ~4 a.u. (sigma = 0.5)
// or ~16 a.u. (sigma = 2) before the packet reached the absorbing bands. A
// stopping power fitted over such a window is not defensible. Removing the CAP
// turns the run into a multi-pass trajectory whose whole length is fittable.
//
// The wrap needs NO boundary-condition change. periodicity(2) is consulted only
// by the Poisson solver (inq/src/solvers/poisson.hpp:189,206); the wavefunction
// basis and the kinetic operator are a plain 3-D FFT, periodic in ALL THREE
// directions (inq/src/basis/fourier_space.hpp:60-151,
// inq/src/hamiltonian/ks_hamiltonian.hpp:200-204). A KS orbital travelling +z
// ALREADY wraps — that is exactly why the CAP had to be added in the first
// place. Switching it off RESTORES the wrap rather than introducing it, and
// periodicity(2) is kept so the slab still has no spurious z images.
//
// TWO CONSEQUENCES, both good:
//   * The Hamiltonian is time-independent and Hermitian again, so energy_total
//     is CONSERVED. That is a far stronger correctness gate than the norm
//     monitoring a CAP'd run has to settle for. Reported at the end of the run.
//   * The WP norm is conserved too, so <p> and <p^2> are expectation values of
//     the WHOLE packet, with none of the CAP's preferential removal of the slow
//     spread tail (the qsp5 bias, docs/handovers/qsp5-momentum-stopping.md).
//
// ---------------------------------------------------------------------------
// STOPPING POWER (definitions from docs/plans/bulk-jellium-ks-stopping.md sec 4)
//   T1 = <p^2>/2m   -> wp_momentum_stats.csv, e_kin_ha
//   T2 = <p>^2/2m   -> 0.5*(px_mean^2+py_mean^2+pz_mean^2), same file
//   s3 = WP density centroid -> wp_real_space_stats.csv, z_mean_circ (CIRCULAR;
//        post-processing UNWRAPS it into a monotone path across the wraps)
//   s4 = integral of <p_z> dt -> cumulative trapezoid of pz_mean
//   s5 = IN-SLAB path (NEW here) -> integral of f(t) <p_z>/m dt, with f the
//        in-slab occupancy written every step to wp_slab_occupancy.csv.
// S_ij = -dT_i/ds_j.
//
// WHY s5 EXISTS. Stopping power is a force: energy lost per unit path INSIDE the
// medium. The slab is 25 of the 85 Bohr the packet traverses, and by t ~ 35 a.u.
// the packet is WIDER than the slab (sigma_d(t) = sqrt(sigma^2/2 + t^2/(2
// sigma^2)) = sqrt(2 + t^2/8)), so it is inside and outside at the same time.
// Fitting -dT/ds against the centroid path then averages the drag over slab AND
// vacuum. With ds5 = f v dt and dT/dt = -F v f, -dT/ds5 = F exactly, in both the
// localised limit (f -> 1) and the delocalised one (f -> 25/85 = 0.294). f is
// MEASURED on the grid every step, not modelled from a Gaussian ansatz, because
// the packet stops being Gaussian as soon as it scatters.
//
// ACCEPTED LIMITATION, stated so the notebooks repeat it: transverse periodic
// images overlap (6 sigma_d = L_xy = 35) at t = 16 a.u., so past that point the
// object being dragged is a periodic ARRAY of packets, not one packet. This is
// not removable without changing L_xy, which would change r_s and break the
// like-for-like comparison with the classical benchmark.
//
// ---------------------------------------------------------------------------
// CHECKPOINTS (.claude/rules/final-timestep-checkpoint.md + checkpoint-dont-block)
//   * results/<name>/checkpoint        — rolling latest; what LJ_RESUME=1 loads
//   * results/<name>/ckpt_step<N>      — RETAINED numbered snapshots (3 MAX,
//                                     oldest pruned; the FINAL step is always one)
//   * results/<name>/rt_state.txt      — last_step, time_au, dt, wp_idx
//
// ENGINE: inq-study (kept for parity with every other run in this line of work;
// with the CAP off the complexified potential is not exercised).
//
// Env: LJ_K0(2.0) LJ_N_STEPS(4529) LJ_DT(0.04) LJ_SIGMA(2.0) LJ_LAUNCH_Z(-24)
//      LJ_N(100) LJ_SPACING(0.40) LJ_CAP_ETA(0=OFF) LJ_CAP_L(12.5)
//      LJ_SAVE_EVERY(38) LJ_WF_EVERY(228) LJ_CKPT_EVERY(0=auto N/3) LJ_MAX_CKPT(3)
//      LJ_RESUME(0) LJ_OUT(REQUIRED) LJ_GS_DIR(REQUIRED)
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/observables/slab_occupancy.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/injection_report.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/interaction_energies.hpp>
#include <inqkit/jellium/analytics.hpp>

#include "../../../shared/configs/slab_n100_L35x35x85.hpp"

#include <chrono>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN100_L35x35x85;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

// rt_state.txt is written as `key=value` lines; tolerate several per line.
static double read_state_d(const std::string& path, const char* key, double def){
    std::ifstream f(path); std::string line; const std::string k = std::string(key) + "=";
    while (std::getline(f, line)) {
        auto p = line.find(k);
        if (p != std::string::npos) return std::atof(line.substr(p + k.size()).c_str());
    }
    return def;
}
static std::string tag6(int n){ std::ostringstream o; o << std::setw(6) << std::setfill('0') << n; return o.str(); }
static std::string iso_now(){
    auto t = std::time(nullptr); auto tm = *std::localtime(&t);
    char b[64]; std::strftime(b, sizeof(b), "%Y-%m-%dT%H:%M:%S", &tm); return std::string(b);
}

int main() {
    auto t_wall0 = std::chrono::steady_clock::now();
    const double HA = 27.211386245988;

    // ---- parameters (geometry pinned to the classical campaign's Cfg) -------
    const double SIGMA_WP = env_d("LJ_SIGMA", 2.0);                  // sigma_WP
    // PRODUCTION GRID = 0.40, NOT the classical campaign's 0.50 (user decision
    // 2026-07-30). dx is a numerical, not physical, parameter: at dx = 0.5 the
    // sigma_WP = 0.5 packet is ONE grid point per sigma and aliases 10.4 % of its
    // z-momentum weight past k_Nyq = pi/dx at v = 4.5, biasing T1 = <p^2>/2 — the
    // primary stopping observable. At 0.40 that falls to <= 0.9 % everywhere, for
    // 2x the grid. The campaign's own <resolved_decisions> called dx = 0.50
    // PROVISIONAL and pre-authorised 0.35-0.40. The GS_DIR must match this dx.
    const double SPACING  = env_d("LJ_SPACING", 0.40);
    const double K0       = env_d("LJ_K0", 2.0);                     // = v (m = 1)
    const double LAUNCH_Z = env_d("LJ_LAUNCH_Z", -24.0);
    const double DT       = env_d("LJ_DT", 0.04);
    const int    N_STEPS  = env_i("LJ_N_STEPS", 4529);               // v=2.0 point
    const double CAP_ETA  = env_d("LJ_CAP_ETA", 0.0);                // 0 = NO CAP (this study)
    const double CAP_L    = env_d("LJ_CAP_L", 12.5);                 // Bohr per face
    const int    SAVE_EVERY = env_i("LJ_SAVE_EVERY", 15);
    const int    WF_EVERY   = env_i("LJ_WF_EVERY", 45);
    const bool   RESUME     = env_i("LJ_RESUME", 0) != 0;
    const std::string OUT    = "results/" + env_s("LJ_OUT", "wp");
    const std::string GS_DIR = env_s("LJ_GS_DIR", "");

    int CKPT_EVERY = env_i("LJ_CKPT_EVERY", 0);
    if (CKPT_EVERY <= 0) CKPT_EVERY = std::max(1, N_STEPS / 3);      // -> 3 retained

    if (GS_DIR.empty() || !std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: LJ_GS_DIR missing or unset: '" << GS_DIR << "'\n";
        return 2;
    }

    const double CAP_WIDTH_FRAC = CAP_L / Cfg::LZ_BOHR;
    const double CAP_MID_FRAC   = 0.5 - CAP_WIDTH_FRAC / 2.0;
    const double z_cap_in       = Cfg::LZ_BOHR/2.0 - CAP_L;          // +30.0
    const bool   CAP_ON         = (CAP_ETA != 0.0);
    const double sigma_p2       = 1.0 / (2.0 * SIGMA_WP * SIGMA_WP); // = 2.0
    // DENSITY IS THE ONLY CAMPAIGN VARIABLE. N comes from the environment so a
    // single binary serves both r_s = 4.18 (N=100, shared/configs/
    // slab_n100_L35x35x85.hpp) and r_s = 5.67 (N=40, slab_n40_L35x35x85.hpp).
    // The slab volume is compile-time, so n0 follows from N alone and the two
    // campaigns cannot differ in anything else by construction.
    const int    N_ELEC         = env_i("LJ_N", Cfg::N_ELECTRONS);
    const double N0             = double(N_ELEC) / Cfg::V_INSIDE_BOHR3;
    const double rs             = inqkit::jellium::rs_from_n0(N0);

    // ---- resume state -------------------------------------------------------
    const std::string CKPT = OUT + "/checkpoint", RT_STATE = OUT + "/rt_state.txt";
    int START = 0, wp_idx = -1;
    if (RESUME) {
        START  = (int)read_state_d(RT_STATE, "last_step", -1);
        wp_idx = (int)read_state_d(RT_STATE, "wp_idx", -1);
        if (START < 0 || wp_idx < 0) {
            std::cerr << "FATAL: LJ_RESUME=1 but no readable " << RT_STATE << "\n"; return 2;
        }
        if (START >= N_STEPS) {
            std::cout << "Already at/after target (" << START << " >= " << N_STEPS
                      << "); nothing to do.\n"; return 0;
        }
    }
    const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

    std::cout << std::setprecision(10)
              << "\n=== slab_ks_wrap/wp  OUT=" << OUT << " ===\n"
              << "  cell      = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR << " x "
              << Cfg::LZ_BOHR << " Bohr, periodicity(2), dx=" << SPACING << "\n"
              << "  slab      = 25 Bohr (half " << Cfg::SLAB_HALF_WIDTH << ", edge "
              << Cfg::EDGE_WIDTH_BOHR << "), N=" << N_ELEC
              << ", n0=" << N0 << ", r_s=" << rs << "\n"
              << "  WP        = sigma_WP " << SIGMA_WP << " (density std "
              << SIGMA_WP/std::sqrt(2.0) << " == classical sigma_pot)"
              << "  k0=" << K0 << " (v=" << K0 << ")  E_drift="
              << 0.5*K0*K0*HA << " eV\n"
              << "  launch_z  = " << LAUNCH_Z << "\n"
              << "  CAP       = " << (CAP_ON ? "ON" : "OFF") << "  eta=" << CAP_ETA
              << " Ha  width=" << CAP_L << " Bohr/face  bands +/-["
              << z_cap_in << "," << Cfg::LZ_BOHR/2.0 << "]\n"
              << "  dt=" << DT << "  START=" << START << " -> N_STEPS=" << N_STEPS
              << "  t_total=" << (DT*N_STEPS) << " a.u." << (RESUME ? "  [RESUME]" : "") << "\n"
              << "  cadence   : density/" << SAVE_EVERY << "  wavefn/" << WF_EVERY
              << "  stats/1  ckpt/" << CKPT_EVERY << "\n"
              << "  GS        = " << GS_DIR << "\n"
              << "  spreading : rate " << 1.0/(std::sqrt(2.0)*SIGMA_WP)
              << " Bohr/a.u.;  T1-T2 = " << 3.0/(4.0*SIGMA_WP*SIGMA_WP)*HA << " eV\n\n";

    // ---- system -------------------------------------------------------------
    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b,
                                            Cfg::LY_BOHR * 1.0_b,
                                            Cfg::LZ_BOHR * 1.0_b).periodicity(2);
    auto ions = systems::ions(cell);                    // jellium: no nuclei
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(N_ELEC)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    const int n_states = electrons.states().num_states();

    // ---- output skeleton ----------------------------------------------------
    const std::string OBS = OUT + "/raw/observables";
    const std::string VTI = OUT + "/raw/vti";
    for (auto const& d : {OBS, VTI + "/density_total", VTI + "/density_wp",
                          VTI + "/density_delta", VTI + "/density_delta_coarse",
                          VTI + "/density_gs_system", VTI + "/wavefunction_wp"})
        std::filesystem::create_directories(d);

    inqkit::InjectionReport report{};

    if (RESUME) {
        electrons.load(CKPT);
        std::cout << "  RESUMED from step " << START << " (wp_idx=" << wp_idx << ")\n";

        // PURGE STALE FRAMES (step > START). They were written by the aborted
        // tail after the checkpoint we are resuming from, so the run is about to
        // recompute those very steps. Two reasons this must happen:
        //
        //  1. CORRECTNESS. Leaving them would mix frames from a discarded
        //     trajectory into the retained one.
        //  2. IT OTHERWISE CRASHES. The writers use overwrite=false on resume
        //     (segment safety, .claude/rules/final-timestep-checkpoint.md), and
        //     frames are named by STEP. When LJ_SAVE_EVERY changes between
        //     segments, the old and new cadences collide on their common
        //     multiples and the writer aborts with "file already exists". That
        //     is exactly what killed the 2026-07-31 resume after the cadence was
        //     reduced 12 -> 30 (they collide every 60 steps).
        //
        // Frames at steps <= START are the retained segment and are NEVER
        // touched, so overwrite=false still guards what it is meant to guard.
        if (electrons.root()) {
            int purged = 0;
            for (auto const& dir : {VTI + "/density_total", VTI + "/density_wp",
                                    VTI + "/density_delta", VTI + "/density_delta_coarse",
                                    VTI + "/wavefunction_wp"}) {
                if (!std::filesystem::exists(dir)) continue;
                for (auto const& e : std::filesystem::directory_iterator(dir)) {
                    if (!e.is_regular_file()) continue;
                    const std::string stem = e.path().stem().string();
                    // trailing run of digits = the step index (density_t002220)
                    std::size_t i = stem.size();
                    while (i > 0 && std::isdigit(static_cast<unsigned char>(stem[i-1]))) --i;
                    if (i == stem.size()) continue;              // no digits: not a frame
                    const int step_of = std::atoi(stem.c_str() + i);
                    if (step_of > START) {
                        std::error_code ec;
                        std::filesystem::remove(e.path(), ec);
                        if (!ec) ++purged;
                    }
                }
            }
            std::cout << "  purged " << purged << " stale frame(s) with step > "
                      << START << " (they will be recomputed)\n" << std::flush;
        }
    } else {
        electrons.load(GS_DIR);
        std::cout << "  Loaded GS from " << GS_DIR << "\n";

        // t=0 BATH density BEFORE the WP goes in — the baseline every induced
        // density in post-processing is measured against.
        inqkit::io::RealField3DLayout lay{
            .field_name = "density", .include_meta = false, .emit_raw = false,
            .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
        inqkit::io::RealField3DWriter gs_wr(VTI + "/density_gs_system", lay, {.overwrite = true});
        gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system");

        auto wp = inqkit::WavePacket{}
                      .center(0.0, 0.0, LAUNCH_Z)
                      .sigma(SIGMA_WP)
                      .k0(0.0, 0.0, K0)
                      .orthogonalise_against_occupied(electrons);
        report = wp.inject_into_last_extra_state(electrons, 1.0);
        wp_idx = report.state_index;
        std::cout << "  WP injected: state_index=" << wp_idx
                  << "  norm_after=" << report.norm_after
                  << "  max_overlap=" << report.max_overlap << "\n";

        std::ofstream f(OBS + "/wp_config.txt");
        f << std::setprecision(16)
          << "wp_center_bohr = 0 0 " << LAUNCH_Z << "\n"
          << "wp_sigma_bohr  = " << SIGMA_WP << "\n"
          << "wp_sigma_density = " << SIGMA_WP/std::sqrt(2.0) << "\n"
          << "wp_k0_bohr_inv = " << K0 << "\n"
          << "wp_state_index = " << wp_idx << "\n"
          << "norm_after     = " << report.norm_after << "\n"
          << "max_overlap    = " << report.max_overlap << "\n";
    }

    // ---- background well + the two CAPs ------------------------------------
    inqkit::jellium::localised_background_params bg;
    bg.shape      = inqkit::jellium::background_shape::slab;
    bg.n0         = N0;
    bg.half_width = Cfg::SLAB_HALF_WIDTH;
    bg.slab_axis  = Cfg::SLAB_AXIS;
    bg.center     = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR};
    bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    // The absorbing bands are CONSTRUCTED ONLY IF REQUESTED. This study runs
    // with CAP_ETA = 0, and an eta = 0 absorbing perturbation is not the same
    // thing as no perturbation at all: it still complexifies the scalar
    // potential and takes the engine down the non-Hermitian path. Branching on
    // CAP_ON keeps the CAP-free run exactly unitary, which is what makes the
    // energy-conservation gate below meaningful.
    const bool WRAP_EXPECTED = !CAP_ON;

    // constant background fields for the pairwise ledger
    auto basis   = electrons.density().basis();
    auto nplus   = bg_pert.background_density(basis);
    auto phiplus = solvers::poisson::solve(nplus);
    const double E_BB = inqkit::jellium::background_self_energy(nplus, phiplus);

    // ---- writers (segment-suffixed on resume) -------------------------------
    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter total_wr(VTI + "/density_total", vti_layout, {.overwrite = (START == 0)});
    inqkit::io::RealField3DWriter wp_wr   (VTI + "/density_wp",    vti_layout, {.overwrite = (START == 0)});
    inqkit::io::ComplexField3DWriter wf_wr(
        VTI + "/wavefunction_wp",
        {.field_name = "wavefunction", .include_meta = false, .emit_raw = false,
         .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary},
        {.overwrite = (START == 0)});

    // FULL energy decomposition — the user asked for every component.
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.energy_external = sel.energy_nonlocal = sel.energy_ion = true;
    sel.energy_ion_kinetic = sel.energy_exact_exchange = true;
    sel.energy_nvxc = sel.energy_eigenvalues = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs(OBS + "/observables" + SEG + ".csv", sel);
    obs.write_header();

    inqkit::observables::WPMomentumStats  wp_mom(OBS + "/wp_momentum_stats"  + SEG + ".csv", wp_idx, {.write_every = 1});
    inqkit::observables::WPRealSpaceStats wp_pos(OBS + "/wp_real_space_stats" + SEG + ".csv", wp_idx, {.write_every = 1});

    // snapshot() is called EVERY step (the density_delta L2 is an every-step
    // scalar in observables.csv), but the 18 MB delta FIELD is only wanted at the
    // density cadence. Without .emit_every this wrote 3624 frames per run instead
    // of 302 -- 66 GB each, which filled /rds on 2026-07-31 and killed three
    // sigma=3 runs mid-flight. Aligning on SAVE_EVERY also matches delta frames
    // one-to-one with density_total / density_wp, which is what the GIF battery
    // pairs them against.
    inqkit::observables::DensityDelta density_delta(
        VTI + "/density_delta", VTI + "/density_delta_coarse",
        {.emit_raw_vti = true, .emit_coarse_vti = true,
         .compute_l2 = true, .coarse_bin_bohr = 3.0,
         .emit_every = (SAVE_EVERY > 0 ? SAVE_EVERY : 1)});

    // IN-SLAB OCCUPANCY f(t) — the measurement that makes the whole-run fit
    // possible (plan sec 3, Window B). One extra grid reduction per step.
    std::ofstream occ_f;
    if (electrons.root()) {
        occ_f.open(OBS + "/wp_slab_occupancy" + SEG + ".csv");
        occ_f << std::setprecision(12)
              << "# slab band: axis=2 center=" << Cfg::SLAB_CENTER_BOHR
              << " half_width=" << Cfg::SLAB_HALF_WIDTH
              << "  (geometric filling factor "
              << (2.0*Cfg::SLAB_HALF_WIDTH/Cfg::LZ_BOHR) << ")\n"
              << "step,time_au,f_in_slab,norm_in_slab,norm_total\n";
    }
    const inqkit::observables::SlabBand SLAB_BAND{
        .axis = Cfg::SLAB_AXIS, .center = Cfg::SLAB_CENTER_BOHR,
        .half_width = Cfg::SLAB_HALF_WIDTH};

    std::ofstream ix;   // pairwise Coulomb ledger, Hartree units
    if (electrons.root()) {
        ix.open(OBS + "/interactions" + SEG + ".csv");
        ix << std::setprecision(12)
           << "step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,"
              "e_hartree_check,e_external_check,norm_wp,norm_total\n";
    }

    // ---- t = 0 analytic gates (abort before burning GPU time) ---------------
    if (START == 0) {
        auto m0 = wp_mom.compute(electrons);
        auto r0 = wp_pos.compute(electrons);
        const double T1 = m0.ekin;
        const double T2 = 0.5*(m0.px*m0.px + m0.py*m0.py + m0.pz*m0.pz);
        int fails = 0;
        // RELATIVE tolerances: sigma_WP = 0.5 on a dx = 0.5 grid is ONE grid point
        // per sigma, so the momentum moments carry a real O(1%) discretisation
        // error (measured in cap_check: sigma_pz^2 +1.6%, <p_z> -0.5% at k0 = 2).
        // Percent-level bounds accept that while still catching a factor-2 blunder.
        auto gate_rel = [&](char const* nm, double got, double want, double relpc){
            const double rel = (want != 0.0) ? 100.0*(got-want)/std::abs(want) : 0.0;
            const bool ok = std::abs(rel) <= relpc;
            std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << nm << ": " << got
                      << "  (expect " << want << ", dev " << rel << " %, tol +/-"
                      << relpc << " %)\n";
            if (!ok) ++fails;
        };
        auto gate_abs = [&](char const* nm, double got, double want, double tol){
            const bool ok = std::abs(got - want) <= tol;
            std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << nm << ": " << got
                      << "  (expect " << want << " +/- " << tol << ")\n";
            if (!ok) ++fails;
        };
        std::cout << "\n  --- t=0 analytic gates ---\n";
        gate_abs("norm (real space)",         r0.N,  1.0, 0.02);
        gate_rel("<p_z> = k0",                m0.pz, K0,  2.0);
        gate_rel("sigma_pz^2 = 1/(2 s^2)",    m0.sz2, sigma_p2, 10.0);
        gate_rel("T1 = (k0^2+3 sp2)/2 (Ha)",  T1, 0.5*(K0*K0 + 3.0*sigma_p2), 3.0);
        gate_rel("T1-T2 = 3/(4 s^2) (Ha)",    T1-T2, 3.0/(4.0*SIGMA_WP*SIGMA_WP), 5.0);
        gate_abs("centroid z (circular)",     r0.zc, LAUNCH_Z, 0.05);
        gate_rel("density std = s/sqrt2",     std::sqrt(r0.sz2), SIGMA_WP/std::sqrt(2.0), 5.0);
        std::cout << "  [info] max_overlap with occupied manifold = "
                  << report.max_overlap << " (want < 1e-3)\n";

        // Momentum-space aliasing diagnostic. The WP k-distribution is centred at
        // k0 with std sigma_p = 1/(sqrt2 sigma) = 1.414; whatever lies beyond
        // k_Nyq = pi/dx folds back and corrupts <p^2> (T1) worst of all. At
        // dx = 0.5 this is 0.12 % at v = 2.0 but 10.4 % at v = 4.5 — a
        // WP-SPECIFIC problem the classical Gaussian CHARGE never had, and the
        // reason the campaign flagged dx = 0.5 as PROVISIONAL. Reported, not
        // gated, so the number lands in every run log and notebook.
        {
            const double sp   = 1.0/(std::sqrt(2.0)*SIGMA_WP);
            const double knyq = M_PI / SPACING;
            const double zsc  = (knyq - K0)/sp;
            const double tail = 0.5*std::erfc(zsc/std::sqrt(2.0));
            std::cout << "  [info] ALIASING: k_Nyq=" << knyq << " sigma_p=" << sp
                      << " (k_Nyq-k0)/sigma_p=" << zsc << " -> " << 100.0*tail
                      << " % of the z-momentum weight is beyond Nyquist\n";
            if (tail > 0.02)
                std::cout << "  [WARN] > 2 % aliased: T1 = <p^2>/2 is biased at this "
                             "velocity on dx=" << SPACING << ".\n";
        }
        if (fails > 0) {
            std::cerr << "\nFATAL: " << fails << " t=0 gate(s) failed — the injected "
                         "packet is not the one this run claims. Aborting.\n";
            return 4;
        }
        std::cout << "  all t=0 gates PASSED\n\n";

        total_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);
        wp_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0);
    }

    // ---- real-time -----------------------------------------------------------
    // With no CAP the Hamiltonian is time-independent and Hermitian, so
    // energy_total must be constant. That is the headline correctness gate for
    // this study and it is one a CAP'd run cannot use at all.
    double e_tot_first = 0.0, e_tot_last = 0.0;
    bool   e_tot_seen  = false;

    inqkit::RealTimeSession rt(ions, electrons, 1);     // callback EVERY step
    rt.add([&](inqkit::StepContext const& ctx) {
        const int step = ctx.step;

        auto n_tot = inqkit::fields::density::total(*ctx.electrons);
        const double l2 = density_delta.snapshot(n_tot, ctx.time_au, step);
        inqkit::StepContext out = ctx; out.density_l2 = l2;
        obs.append(out);

        // pairwise Coulomb ledger with the WP as charge group P (Poisson-linearity
        // form, so the terms close EXACTLY against INQ's own E_hartree/E_external)
        auto n_wp_f = inqkit::jellium::orbital_density_field(*ctx.electrons, wp_idx);
        auto ct = inqkit::jellium::compute_coulomb_wp(ctx.electrons->density(), n_wp_f, phiplus);
        if (ctx.electrons->root())
            ix << step << ',' << ctx.time_au << ',' << ct.e_ss << ',' << ct.e_pp << ','
               << ct.e_ps << ',' << ct.e_sb << ',' << ct.e_pb << ',' << E_BB << ','
               << ct.e_hartree_check << ',' << ct.e_external_check << ','
               << ct.norm_wp << ',' << ct.norm_total << '\n';

        // in-slab occupancy of the WP orbital, exact on the grid
        auto occ = inqkit::observables::slab_occupancy(*ctx.electrons, wp_idx, SLAB_BAND);
        if (ctx.electrons->root())
            occ_f << step << ',' << ctx.time_au << ',' << occ.fraction << ','
                  << occ.norm_in << ',' << occ.norm_total << '\n';

        // energy-conservation tracking (meaningful ONLY because there is no CAP)
        if (!e_tot_seen) { e_tot_first = ctx.energy_total; e_tot_seen = true; }
        e_tot_last = ctx.energy_total;

        if (SAVE_EVERY > 0 && step % SAVE_EVERY == 0) {
            total_wr.write(n_tot, ctx.time_au, step);
            wp_wr.write(inqkit::fields::density::orbital(*ctx.electrons, wp_idx), ctx.time_au, step);
        }
        if (WF_EVERY > 0 && step % WF_EVERY == 0) {
            auto wf = inqkit::fields::orbital::wavefunction(*ctx.electrons, wp_idx);
            wf_wr.write(wf, "wavefunction_t" + tag6(step));
        }
    });

    // At most MAX_CKPT retained numbered snapshots (user instruction
    // 2026-07-31). Names are zero-padded, so lexicographic order IS step order
    // and the oldest sort first. The rolling `checkpoint` that LJ_RESUME loads
    // is a separate directory and is never pruned.
    const int MAX_CKPT = env_i("LJ_MAX_CKPT", 3);
    auto prune_ckpts = [&](){
        if (!electrons.root()) return;
        std::vector<std::filesystem::path> snaps;
        for (auto const& e : std::filesystem::directory_iterator(OUT))
            if (e.is_directory() && e.path().filename().string().rfind("ckpt_step", 0) == 0)
                snaps.push_back(e.path());
        std::sort(snaps.begin(), snaps.end());
        while ((int)snaps.size() > MAX_CKPT) {
            std::error_code ec;
            std::filesystem::remove_all(snaps.front(), ec);
            std::cout << "  [ckpt] pruned " << snaps.front().filename().string()
                      << " (keeping newest " << MAX_CKPT << ")\n" << std::flush;
            snaps.erase(snaps.begin());
        }
    };

    auto write_rt_state = [&](int last){
        if (!electrons.root()) return;
        std::ofstream st(RT_STATE);
        st << std::setprecision(12)
           << "last_step=" << last << "\ntime_au=" << (last*DT) << "\ndt=" << DT
           << "\nwp_idx=" << wp_idx << "\nk0=" << K0 << "\nsigma_wp=" << SIGMA_WP << "\n";
    };

    // The step callback is identical either way; only the perturbation type
    // differs, so the propagation is wrapped in a generic lambda and called from
    // exactly one of the two branches below.
    auto step_callback = [&](auto const& data) {
        rt.step(data);
        wp_mom.maybe_accumulate(data);
        wp_pos.maybe_accumulate(data);
        // Interior checkpoints: rolling `checkpoint` (what LJ_RESUME loads) plus
        // a RETAINED numbered snapshot, pruned to the newest MAX_CKPT.
        if (data.iter() > 0 && data.iter() % CKPT_EVERY == 0) {
            electrons.save(CKPT);
            write_rt_state(data.iter());
            const std::string snap = OUT + "/ckpt_step" + tag6(data.iter());
            electrons.save(snap);
            prune_ckpts();
            std::cout << "  [ckpt] step " << data.iter() << " -> " << snap << "\n" << std::flush;
        }
    };
    auto do_propagate = [&](auto const& pert) {
        real_time::propagate(
            ions, electrons, step_callback,
            options::theory{}.lda(),
            options::real_time{}.num_steps(N_STEPS).dt(DT * 1.0_atomictime),
            pert, START);
    };

    if (CAP_ON) {
        perturbations::absorbing cap_lo(CAP_ETA * 1.0_Ha, -CAP_MID_FRAC, CAP_WIDTH_FRAC);
        perturbations::absorbing cap_hi(CAP_ETA * 1.0_Ha,  CAP_MID_FRAC, CAP_WIDTH_FRAC);
        do_propagate(perturbations::sum(bg_pert, perturbations::sum(cap_lo, cap_hi)));
    } else {
        do_propagate(bg_pert);          // unitary: the packet wraps, nothing absorbs
    }

    // ---- final checkpoint ---------------------------------------------------
    // The LAST TIMESTEP is saved TWICE on purpose (.claude/rules/final-timestep-
    // checkpoint.md + user instruction 2026-07-31 "ensure the last timestep is
    // timestamped"): once as the rolling `checkpoint` that LJ_RESUME=1 loads,
    // and once as a STEP-STAMPED snapshot ckpt_step<N_STEPS>, so the final state
    // is identifiable by step number instead of being an anonymous directory
    // whose provenance has to be read out of rt_state.txt. Pruning runs AFTER it
    // is written and it sorts last, so the final state can never be pruned away.
    electrons.save(CKPT);
    {
        const std::string final_snap = OUT + "/ckpt_step" + tag6(N_STEPS);
        electrons.save(final_snap);
        prune_ckpts();
        std::cout << "  [ckpt] FINAL step " << N_STEPS << " -> " << final_snap
                  << "  (t = " << (N_STEPS*DT) << " a.u.)\n" << std::flush;
    }
    write_rt_state(N_STEPS);

    const double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_wall0).count();

    // ---- energy-conservation gate (the CAP-free study's headline check) -----
    // A time-independent Hermitian H must conserve energy_total exactly. Any
    // real drift means the propagation itself is untrustworthy and every S below
    // is meaningless, so it is REPORTED LOUDLY rather than buried in a CSV.
    const double e_drift_ev = (e_tot_last - e_tot_first) * HA;
    if (!CAP_ON && e_tot_seen) {
        const bool ok = std::abs(e_drift_ev) < 1.0e-3;
        std::cout << (ok ? "\n  [PASS] " : "\n  [WARN] ")
                  << "energy conservation: E_total drift = " << e_drift_ev
                  << " eV over " << (N_STEPS - START) << " steps"
                  << " (want < 1e-3 eV; no CAP => H is Hermitian and t-independent)\n";
    }

    if (electrons.root()) {
        ix.close();
        occ_f.close();
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(14)
          << "run = localised_jellium/slab_ks_wrap/wp/" << env_s("LJ_OUT","wp") << "\n"
          << "run_type = wavepacket projectile, localised jellium slab TDDFT (ALDA)\n"
          << "campaign = slab-ks-orbital-stopping-wrap (CAP-free, wrap-around)\n"
          << "plan = docs/plans/slab-ks-orbital-stopping-wrap.md\n"
          << "cap_enabled = " << (CAP_ON ? "yes" : "no") << "\n"
          << "wrap_around = " << (WRAP_EXPECTED ? "yes" : "no")
          << "  (KS orbitals wrap on the FFT grid; periodicity(2) is electrostatic only)\n"
          << "energy_total_first_ha = " << e_tot_first << "\n"
          << "energy_total_last_ha = " << e_tot_last << "\n"
          << "energy_total_drift_ev = " << e_drift_ev << "\n"
          << "slab_filling_factor = " << (2.0*Cfg::SLAB_HALF_WIDTH/Cfg::LZ_BOHR) << "\n"
          << "engine = inq-study\nxc = LDA (ALDA in TDDFT)\n"
          << "date_finished = " << iso_now() << "\nwall_time_s = " << wall << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR << "\n"
          << "periodicity = 2\nspacing_bohr = " << SPACING << "\n"
          << "slab_half_width = " << Cfg::SLAB_HALF_WIDTH
          << "  edge_width = " << Cfg::EDGE_WIDTH_BOHR << "\n"
          << "n_electrons = " << N_ELEC << "  n0_a0m3 = " << N0
          << "  r_s = " << rs << "\n"
          << "n_states = " << n_states << "  extra_states = " << Cfg::EXTRA_STATES << "\n"
          << "wp_enabled = yes\nwp_state_index = " << wp_idx << "\n"
          << "wp_sigma_bohr = " << SIGMA_WP << "\n"
          << "wp_sigma_note = wavepacket sigma (psi width); density std = this/sqrt2 "
             "= classical sigma_pot\n"
          << "wp_sigma_density = " << SIGMA_WP/std::sqrt(2.0) << "\n"
          << "wp_k0_bohr_inv = " << K0 << "  wp_velocity = " << K0 << "\n"
          << "wp_drift_energy_ev = " << 0.5*K0*K0*HA << "\n"
          << "wp_localisation_energy_ev = " << 3.0/(4.0*SIGMA_WP*SIGMA_WP)*HA << "\n"
          << "launch_z = " << LAUNCH_Z << "\n"
          << "norm_after = " << report.norm_after << "  max_overlap = " << report.max_overlap << "\n"
          << "cap = " << (CAP_ON ? "on" : "off") << "  cap_eta_ha = " << CAP_ETA
          << "  cap_width_bohr = " << CAP_L << " per face\n"
          << "cap_mid_frac = " << CAP_MID_FRAC << "  cap_width_frac = " << CAP_WIDTH_FRAC << "\n"
          << "cap_band_hi_bohr = [" << z_cap_in << "," << Cfg::LZ_BOHR/2.0 << "]\n"
          << "cap_band_lo_bohr = [" << -Cfg::LZ_BOHR/2.0 << "," << -z_cap_in << "]\n"
          << "cap_note = energy_total NOT conserved (non-Hermitian); gate on norm instead\n"
          << "start_step = " << START << "  rt_num_steps = " << N_STEPS
          << "  dt_au = " << DT << "  total_time_au = " << (DT*N_STEPS) << "\n"
          << "save_every = " << SAVE_EVERY << "  wf_every = " << WF_EVERY
          << "  stats_every = 1  ckpt_every = " << CKPT_EVERY << "\n"
          << "spread_rate_bohr_per_au = " << 1.0/(std::sqrt(2.0)*SIGMA_WP) << "\n"
          << "t_transverse_overlap_au = 4.12  (6 sigma_d = L_xy = 35; step "
          << int(4.12/DT) << ")\n"
          << "analysis_note = fit S over the pre-overlap window; later steps are "
             "recorded but the packet is delocalised\n"
          << "gs_dir = " << GS_DIR << "\nrun_completed = true\n";
    }
    std::cout << "\nDone. Wall time " << wall << " s.\n";
    return 0;
}
