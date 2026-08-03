// ============================================================================
// cylindrical_jellium / scripts/proximity_ladder / wp / run.cpp
//
// QUANTUM HALF of the twin at ONE RUNG of the proximity ladder: an electron
// WAVEPACKET (sigma_WP = 4 Bohr, drift k0 = 1.917 = 50 eV) injected as an occupied
// KS orbital and launched on-axis at z = -28 down a periodic r_s = 3 jellium tube
// whose bore radius R_in is the rung parameter.
// Plan: docs/plans/cylindrical-proximity-ladder.md
//
// Twin: ../classical/run.cpp — same GS, cell, grid, N, sigma, launch point, dt,
// n_steps. The ONLY difference is the projectile representation.
//
// ---------------------------------------------------------------------------
// BUILD ONCE, RUN PER RUNG. CJ_RUNG selects the geometry at RUNTIME from
// shared/configs/proximity_ladder_rs3.hpp; there is no default, because a
// silently-defaulted rung yields a plausible run at the wrong geometry:
//
//   rung   R_in  R_out   N_e  states  WP charge in wall @t=0   shape
//   r10    10.0  14.000  160    104          0.19 %           annulus (DONE)
//   r08     8.0  14.000  220    143          1.83 %           annulus
//   r06     6.0  13.986  266    173         10.54 %           annulus
//   r04     4.0  14.000  300    195         36.79 %           annulus
//   r00     0.0  13.986  326    212          100 %            CYLINDER (filled)
//
// r_s is held at 3.000000 on every rung (R_out is SOLVED, not transcribed), the
// projectile energy at 50 eV, and the cell/grid/dt/n_steps unchanged — so the
// rungs differ ONLY in how much bath the projectile is immersed in.
//
// WHAT THE LADDER CANNOT REACH, stated so the notebooks repeat it: the Gaussian
// form factor exp(-q^2 sigma_pot^2/2) is 0.018 at q = 1 and 3e-26 at q = 2 v0.
// This projectile couples to the COLLECTIVE response and essentially nothing else,
// at EVERY rung. The ladder is weak-collective -> strong-collective; the
// electron-hole pair channel is a sigma_WP axis, not an R_in axis.
//
// ---------------------------------------------------------------------------
// WHAT THIS RUN IS FOR
// Validate a KS-ORBITAL definition of stopping power against the classical
// deltaE/ds one. In BULK the wavepacket's kinetic energy was contaminated: the
// momentum-SPREAD term var(p) grew by +6.8 eV through interaction with the bath
// while the DRIFT term 1/2<p>^2 stayed flat, so -dT1/ds was not a stopping power.
// var(p) is CONSERVED under free evolution, so the growth was interaction and
// nothing else. Channeling suppresses that interaction — the packet flies through
// the vacuum bore and couples to the wall only through the smooth image force —
// so dT/dt should collapse onto the drift channel, which IS the classical
// projectile kinetic energy.
//
// THE CLAIM IS THEREFORE TESTABLE THREE WAYS, all measured here:
//   (a) S_2j (drift-based) agrees with the classical S            <- the aim
//   (b) var(p) = sigma_pz^2(t) stays FROZEN                       <- the mechanism
//   (c) f_bore(t) stays ~1                                        <- the premise
// (b) and (c) are what turn (a) from a coincidence into an explanation.
//
// ---------------------------------------------------------------------------
// STOPPING-POWER DEFINITIONS (S_ij = -dT_i/ds_j)
//   T1 = <p^2>/2m          wp_momentum_stats.csv, e_kin_ha
//   T2 = <p>^2/2m          0.5*(px_mean^2+py_mean^2+pz_mean^2), same file
//   s3 = density centroid  wp_real_space_stats.csv, z_mean_circ (CIRCULAR)
//   s4 = integral <p_z> dt cumulative trapezoid of pz_mean
// T1 - T2 = var(p)/2m is the localisation+scattering term; it is CONSTANT at
// 3/(4 sigma_WP^2) = 1.2755 eV for free evolution, so its drift is the direct
// readout of diagnostic (b).
//
// NO s5 / IN-MEDIUM PATH CORRECTION IS NEEDED, and that is a genuine advantage of
// the tube over the slab. The tube is UNIFORM along z, so the medium fills every
// z the projectile visits and the in-medium path IS the path. The slab study
// needed s5 = integral f v dt only because 25 of its 85 Bohr were vacuum. What
// the tube needs instead is the RADIAL question — is the packet still in the
// bore — which is f_bore(t) from inqkit::observables::radial_occupancy.
//
// ---------------------------------------------------------------------------
// THE CIRCULAR CENTROID IS MANDATORY HERE
// The packet is launched 2 Bohr from the -z face with a density std of 2.83, so
// ~24 % of it is on the far side of the periodic cell AT t = 0. The naive <z> is
// meaningless in that state — it averages the two halves and slides smoothly to a
// wrong answer rather than jumping, so it cannot be repaired afterwards. z_mean_circ
// (Resta phase estimator) is exact and is what s3 uses; the t=0 gates below check
// z_mean_circ, not z_mean.
//
// AND THE PACKET MUST BE INJECTED WITH minimum_image(true). It is tempting to say
// <p_z>, T1 and T2 are unaffected by the straddle because they are momentum-space
// expectation values — that is TRUE of a genuinely wrapped packet and FALSE of a
// CLIPPED one. WavePacket builds its Gaussian from a plain Cartesian displacement
// unless told otherwise, so without the flag the packet is truncated at the face,
// its sharp edge is broadband in momentum, and var(p_z) came out FIFTEEN TIMES too
// large (measured 2026-08-01: six of the nine t=0 gates below failed, all of them
// in z, while x and y were perfect). A KS orbital wraps exactly under PROPAGATION;
// its INJECTION does not wrap unless asked.
//
// ---------------------------------------------------------------------------
// NO CAP, SO ENERGY IS CONSERVED
// The Hamiltonian is time-independent and Hermitian, so energy_total must be
// constant; that is the headline correctness gate and is reported at the end. The
// projectile covers 57.5 Bohr in 30 a.u. and ends at z = +29.5, i.e. ONE traversal
// with no wrap — but a KS orbital wraps exactly on the FFT basis anyway, so a
// resumed longer run needs no boundary change.
//
// ACCEPTED LIMITATIONS, stated so the notebooks repeat them:
//   * The packet SPREADS (sigma_d: 2.83 at t=0 -> 6.01 at t=30), so a rung's
//     nominal R_in/sigma understates the coupling at late time and the rungs
//     MERGE: the wall-overlap fraction spans a factor of 190 across the ladder at
//     t = 0 but only a factor of 3 by t = 30. The rungs are distinct only EARLY.
//     Consequently S must be fitted over a window defined by a COMMON MEASURED
//     COUPLING across rungs, never a common time — a fixed-time fit compares
//     different couplings and would manufacture a trend. f_bore(t) / f_wall(t)
//     from radial_occupancy are the measured versions; gate on those, not on the
//     Rayleigh formula.
//   * 2 sigma_d reaches the bore wall at t = 23.3 (r10), 17.9 (r08), 10.6 (r06)
//     and 0 a.u. (r04 — the packet already overlaps at injection). r00 is
//     immersed from the start.
//   * transverse periodic images overlap (6 sigma_d = L_xy = 40) only at
//     t = 34.1 a.u., AFTER the run ends — unlike the slab study, this run never
//     drags a periodic array.
//
// Env: CH_GS_DIR(REQUIRED) CH_OUT(wp) CH_N(160) CH_V0 CH_SIGMA CH_LAUNCH_Z
//      CH_N_STEPS CH_DT CH_SPACING CH_SAVE_EVERY CH_WF_EVERY
//      CH_CKPT_EVERY(0=auto N/3) CH_MAX_CKPT(3) CH_RESUME(0)
// Build against INQ_SOURCE=inq-study; runtime shares from inq/install/share.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/radial_occupancy.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/injection_report.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/interaction_energies.hpp>
#include <inqkit/jellium/analytics.hpp>

#include "../../../shared/configs/proximity_ladder_rs3.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
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
namespace cfg = cylindrical_jellium::config;
using Cfg = cfg::ChannelingTubeRs3;   // everything SHARED across rungs

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

static double read_state_d(const std::string& path, const char* key, double def){
    std::ifstream f(path); std::string line; const std::string k = std::string(key) + "=";
    while(std::getline(f, line)){ auto p = line.find(k);
        if(p != std::string::npos) return std::atof(line.substr(p + k.size()).c_str()); }
    return def;
}
static std::string tag6(int n){ std::ostringstream o; o << std::setw(6) << std::setfill('0') << n; return o.str(); }
static std::string iso_now(){
    auto t = std::time(nullptr); auto tm = *std::localtime(&t);
    char b[64]; std::strftime(b, sizeof(b), "%Y-%m-%dT%H:%M:%S", &tm); return std::string(b);
}

int main(){
    auto t_wall0 = std::chrono::steady_clock::now();
    const double HA = 27.211386245988;

    // ---- rung selection: required, never defaulted ------------------------
    // A silently-defaulted rung produces a perfectly plausible run at the wrong
    // geometry, which no downstream gate would catch.
    const cfg::Rung* rung = cfg::rung_from_env();
    if(rung == nullptr){
        const char* got = std::getenv("CJ_RUNG");
        std::cerr << "FATAL: CJ_RUNG "
                  << (got ? "= '" + std::string(got) + "' is not a known rung" : "is unset")
                  << ".\n       Valid: ";
        for(int i = 0; i < cfg::LADDER_N; ++i) std::cerr << cfg::LADDER[i].label << " ";
        std::cerr << cfg::CONTROL_SAME_N.label << "\n";
        return 2;
    }
    if(!(cfg::max_density_error() < 1e-12)){
        std::cerr << "FATAL: rung table density error " << cfg::max_density_error()
                  << " >= 1e-12 — R_out is not solving to the r_s = 3 target.\n";
        return 2;
    }
    const double R_IN    = rung->r_in;
    const double R_OUT   = rung->r_out();
    const int    EXTRA_ST= rung->extra_states;
    const bool   FILLED  = rung->filled();

    const int    N_ELEC   = env_i("CH_N",         rung->n_electrons);
    const double SPACING  = env_d("CH_SPACING",   Cfg::SPACING_BOHR);
    const double SIGMA_WP = env_d("CH_SIGMA",     Cfg::SIGMA_WP_BOHR);
    const double K0       = env_d("CH_V0",        Cfg::PROJ_V0);      // = v (m = 1)
    const double LAUNCH_Z = env_d("CH_LAUNCH_Z",  Cfg::LAUNCH_Z_BOHR);
    const double DT       = env_d("CH_DT",        Cfg::DT_AU);
    const int    N_STEPS  = env_i("CH_N_STEPS",   Cfg::N_STEPS);
    const int    SAVE_EVERY = env_i("CH_SAVE_EVERY", Cfg::SAVE_EVERY);
    const int    WF_EVERY   = env_i("CH_WF_EVERY",   Cfg::WF_EVERY);
    const bool   RESUME     = env_i("CH_RESUME", 0) != 0;
    const std::string OUT    = "results/" + env_s("CH_OUT", "wp");
    const std::string GS_DIR = env_s("CH_GS_DIR", "");

    int CKPT_EVERY = env_i("CH_CKPT_EVERY", 0);
    if(CKPT_EVERY <= 0) CKPT_EVERY = std::max(1, N_STEPS / 3);

    if(GS_DIR.empty() || !std::filesystem::exists(GS_DIR)){
        std::cerr << "FATAL: CH_GS_DIR missing or unset: '" << GS_DIR << "'\n"; return 2;
    }

    const double SIGMA_POT = SIGMA_WP / std::sqrt(2.0);   // = WP DENSITY std
    const double sigma_p2  = 1.0 / (2.0 * SIGMA_WP * SIGMA_WP);
    const double V_ANN     = rung->v_jellium();
    const double N0        = double(N_ELEC) / V_ANN;
    const double RS        = inqkit::jellium::rs_from_n0(N0);
    const double OMEGA_P   = std::sqrt(4.0*M_PI*N0);
    const double VF        = inqkit::jellium::k_fermi_n0(N0);

    // ---- resume state -------------------------------------------------------
    const std::string CKPT = OUT + "/checkpoint", RT_STATE = OUT + "/rt_state.txt";
    int START = 0, wp_idx = -1;
    if(RESUME){
        START  = (int)read_state_d(RT_STATE, "last_step", -1);
        wp_idx = (int)read_state_d(RT_STATE, "wp_idx", -1);
        if(START < 0 || wp_idx < 0){
            std::cerr << "FATAL: CH_RESUME=1 but no readable " << RT_STATE << "\n"; return 2; }
        if(START >= N_STEPS){
            std::cout << "Already at/after target (" << START << " >= " << N_STEPS
                      << "); nothing to do.\n"; return 0; }
    }
    const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

    std::cout << std::setprecision(10)
              << "\n=== proximity_ladder/wp  rung=" << rung->label << "  OUT=" << OUT << " ===\n"
              << "  tube      = R_in " << R_IN << "  R_out " << R_OUT
              << "  L_z " << Cfg::LZ_BOHR << "  (fully periodic)\n"
              << "  cell      = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR << " x "
              << Cfg::LZ_BOHR << "  dx = " << SPACING << "\n"
              << "  bath      = " << N_ELEC << " e   n0 = " << N0 << "   r_s = " << RS
              << "   omega_p = " << OMEGA_P*HA << " eV   v_F = " << VF << "\n"
              << "  WP        = sigma_WP " << SIGMA_WP << " (density std " << SIGMA_POT
              << " == classical sigma_pot)  k0 = " << K0
              << "  E_drift = " << 0.5*K0*K0*HA << " eV  v/v_F = " << K0/VF << "\n"
              << "  launch    = (0,0," << LAUNCH_Z << ")\n"
              << "  dt = " << DT << "  START = " << START << " -> N_STEPS = " << N_STEPS
              << "  t_total = " << DT*N_STEPS << " a.u." << (RESUME ? "  [RESUME]" : "") << "\n"
              << "  cadence   : density/" << SAVE_EVERY << "  wavefn/" << WF_EVERY
              << "  stats/1  ckpt/" << CKPT_EVERY << "\n"
              << "  dispersion: sigma_d(0) = " << Cfg::sigma_d(0.0)
              << "  sigma_d(t_end) = " << Cfg::sigma_d(DT*N_STEPS)
              << " ;  2 sigma_d reaches R_in at t = "
              << std::sqrt(2.0*SIGMA_WP*SIGMA_WP*((R_IN/2.0)*(R_IN/2.0)
                                                  - SIGMA_WP*SIGMA_WP/2.0)) << " a.u.\n"
              << "  T1-T2     = " << 3.0/(4.0*SIGMA_WP*SIGMA_WP)*HA
              << " eV (constant if var(p) stays frozen)\n"
              << "  GS        = " << GS_DIR << "\n\n";

    // ---- system -------------------------------------------------------------
    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b,
                                            Cfg::LY_BOHR * 1.0_b,
                                            Cfg::LZ_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);                   // jellium: no nuclei
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(N_ELEC)
            .extra_states(EXTRA_ST)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    const int n_states = electrons.states().num_states();

    // ---- output skeleton ----------------------------------------------------
    const std::string OBS = OUT + "/raw/observables";
    const std::string VTI = OUT + "/raw/vti";
    for(auto const& d : {OBS, VTI + "/density_total", VTI + "/density_wp",
                         VTI + "/density_delta", VTI + "/density_delta_coarse",
                         VTI + "/density_gs_system", VTI + "/wavefunction_wp"})
        std::filesystem::create_directories(d);

    inqkit::InjectionReport report{};

    if(RESUME){
        electrons.load(CKPT);
        std::cout << "  RESUMED from step " << START << " (wp_idx = " << wp_idx << ")\n";
        // Purge frames written by the ABORTED tail (step > START): the run is about
        // to recompute those steps, and the writers use overwrite=false on resume
        // so a name collision would abort the job.
        if(electrons.root()){
            int purged = 0;
            for(auto const& dir : {VTI + "/density_total", VTI + "/density_wp",
                                   VTI + "/density_delta", VTI + "/density_delta_coarse",
                                   VTI + "/wavefunction_wp"}){
                if(!std::filesystem::exists(dir)) continue;
                for(auto const& e : std::filesystem::directory_iterator(dir)){
                    if(!e.is_regular_file()) continue;
                    const std::string stem = e.path().stem().string();
                    std::size_t i = stem.size();
                    while(i > 0 && std::isdigit(static_cast<unsigned char>(stem[i-1]))) --i;
                    if(i == stem.size()) continue;
                    if(std::atoi(stem.c_str() + i) > START){
                        std::error_code ec; std::filesystem::remove(e.path(), ec);
                        if(!ec) ++purged;
                    }
                }
            }
            std::cout << "  purged " << purged << " stale frame(s) with step > " << START << "\n" << std::flush;
        }
    } else {
        electrons.load(GS_DIR);
        std::cout << "  Loaded GS from " << GS_DIR << "\n";

        // t = 0 BATH density BEFORE the WP goes in — the baseline every induced
        // density in post-processing is measured against.
        inqkit::io::RealField3DLayout lay{
            .field_name = "density", .include_meta = false, .emit_raw = false,
            .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
        inqkit::io::RealField3DWriter gs_wr(VTI + "/density_gs_system", lay, {.overwrite = true});
        gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system");

        // minimum_image: MANDATORY here. The launch point is 2 Bohr = 0.71 density
        // sigma from the -z face, so a plain-Cartesian Gaussian is TRUNCATED, not
        // wrapped — and the truncation's sharp edge inflates var(p_z) by ~15x.
        // It also matches the classical twin, which already builds its charge
        // with gaussian_density_minimum_image; a clipped packet against a wrapped
        // charge is not a twin at the boundary this study introduces on purpose.
        auto wp = inqkit::WavePacket{}
                      .center(0.0, 0.0, LAUNCH_Z)
                      .sigma(SIGMA_WP)
                      .k0(0.0, 0.0, K0)
                      .minimum_image(true)
                      .orthogonalise_against_occupied(electrons);
        report = wp.inject_into_last_extra_state(electrons, 1.0);
        wp_idx = report.state_index;
        std::cout << "  WP injected: state_index = " << wp_idx
                  << "  norm_after = " << report.norm_after
                  << "  max_overlap = " << report.max_overlap << "\n";

        std::ofstream f(OBS + "/wp_config.txt");
        f << std::setprecision(16)
          << "wp_center_bohr = 0 0 " << LAUNCH_Z << "\n"
          << "wp_sigma_bohr  = " << SIGMA_WP << "\n"
          << "wp_sigma_density = " << SIGMA_POT << "\n"
          << "wp_k0_bohr_inv = " << K0 << "\n"
          << "wp_state_index = " << wp_idx << "\n"
          << "norm_after     = " << report.norm_after << "\n"
          << "max_overlap    = " << report.max_overlap << "\n";
    }

    // ---- background + constant fields for the ledger ------------------------
    inqkit::jellium::localised_background_params bg;
    // FILLED rungs are a DIFFERENT SHAPE, not annulus with R_in = 0: the erfc step

    // is centred ON its nominal edge, so a degenerate inner edge would put n_+ = n0/2

    // EXACTLY on the tube axis — precisely where this projectile flies.

    bg.shape        = FILLED ? inqkit::jellium::background_shape::cylinder

                             : inqkit::jellium::background_shape::annulus;
    bg.n0           = N0;
    bg.half_width   = R_OUT;
    bg.inner_radius = R_IN;
    bg.slab_axis    = Cfg::TUBE_AXIS;
    bg.center       = {Cfg::CENTER_X, Cfg::CENTER_Y, Cfg::CENTER_Z};
    bg.edge_width   = Cfg::EDGE_W_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

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
    // Full 1-D n(k) at the density cadence: it feeds the momentum battery of the
    // run-notebook (n(k) before/after, momentum carpet, KL divergence) and is the
    // distribution behind the var(p) scalar.
    inqkit::observables::MomentumDistribution momentum_dist(
        OBS + "/momentum_distribution" + SEG + ".csv", wp_idx, Cfg::LZ_BOHR,
        {.n_bins = 128, .k_max_bohr_inv = 0.0,
         .write_every = (SAVE_EVERY > 0 ? SAVE_EVERY : 1)});

    inqkit::observables::DensityDelta density_delta(
        VTI + "/density_delta", VTI + "/density_delta_coarse",
        {.emit_raw_vti = true, .emit_coarse_vti = true,
         .compute_l2 = true, .coarse_bin_bohr = 3.0,
         .emit_every = (SAVE_EVERY > 0 ? SAVE_EVERY : 1)});

    // THE CHANNELING DIAGNOSTIC. f_bore(t) is the premise of the whole study made
    // falsifiable; <r_perp>(t) says whether the packet is still centred on the axis.
    const inqkit::observables::RadialBand TUBE_BAND{
        .axis = Cfg::TUBE_AXIS,
        .center = {Cfg::CENTER_X, Cfg::CENTER_Y, Cfg::CENTER_Z},
        .r_inner = R_IN, .r_outer = R_OUT};

    std::ofstream occ_f, ix;
    if(electrons.root()){
        occ_f.open(OBS + "/wp_radial_occupancy" + SEG + ".csv");
        occ_f << std::setprecision(12)
              << "# tube: axis=" << Cfg::TUBE_AXIS << " R_in=" << R_IN
              << " R_out=" << R_OUT << "\n"
              << "step,time_au,f_bore,f_wall,f_outside,r_mean,r2_mean,sigma_r,norm_total\n";

        ix.open(OBS + "/interactions" + SEG + ".csv");
        ix << std::setprecision(12)
           << "step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,"
              "e_hartree_check,e_external_check,norm_proj,norm_electrons\n";
    }

    // ---- t = 0 analytic gates (abort before burning GPU time) ---------------
    if(START == 0){
        auto m0 = wp_mom.compute(electrons);
        auto r0 = wp_pos.compute(electrons);
        auto o0 = inqkit::observables::radial_occupancy(electrons, wp_idx, TUBE_BAND);
        const double T1 = m0.ekin;
        const double T2 = 0.5*(m0.px*m0.px + m0.py*m0.py + m0.pz*m0.pz);
        int fails = 0;
        auto gate_rel = [&](char const* nm, double got, double want, double relpc){
            const double rel = (want != 0.0) ? 100.0*(got-want)/std::abs(want) : 0.0;
            const bool ok = std::abs(rel) <= relpc;
            std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << nm << ": " << got
                      << "  (expect " << want << ", dev " << rel << " %, tol +/-" << relpc << " %)\n";
            if(!ok) ++fails;
        };
        auto gate_abs = [&](char const* nm, double got, double want, double tol){
            const bool ok = std::abs(got - want) <= tol;
            std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << nm << ": " << got
                      << "  (expect " << want << " +/- " << tol << ")\n";
            if(!ok) ++fails;
        };
        std::cout << "\n  --- t=0 analytic gates ---\n";
        gate_abs("norm (real space)",           r0.N,  1.0, 0.02);
        gate_rel("<p_z> = k0",                  m0.pz, K0,  1.0);
        gate_rel("sigma_pz^2 = 1/(2 s^2)",      m0.sz2, sigma_p2, 5.0);
        gate_rel("T1 = (k0^2 + 3 sp2)/2 (Ha)",  T1, 0.5*(K0*K0 + 3.0*sigma_p2), 2.0);
        gate_rel("T1 - T2 = 3/(4 s^2) (Ha)",    T1-T2, 3.0/(4.0*SIGMA_WP*SIGMA_WP), 5.0);
        // CIRCULAR centroid / spread: the packet straddles the -z face at launch,
        // so the naive <z> is meaningless here BY CONSTRUCTION and is not gated.
        gate_abs("centroid z (CIRCULAR)",       r0.zc,  LAUNCH_Z, 0.05);
        gate_rel("spread sigma_z (CIRCULAR)",   r0.szc, SIGMA_POT, 5.0);
        gate_rel("transverse std = s/sqrt2",    std::sqrt(r0.sx2), SIGMA_POT, 5.0);
        // The coupling premise at t = 0, from the transverse Rayleigh tail:
        //   f_bore = 1 - exp(-R_in^2 / 2 sigma_d^2)
        // = 0.9981 (r10), 0.9817 (r08), 0.8946 (r06), 0.6321 (r04).
        //
        // A FILLED rung has NO bore, so that expression is identically 0 — and
        // gate_rel treats want == 0 as an automatic pass, i.e. the gate would
        // silently stop testing anything on precisely the rung where the geometry
        // is newest. Gate the meaningful statement instead: the packet must lie
        // wholly INSIDE the cylinder, f_wall = 1 - exp(-R_out^2/2 sigma_d^2)
        // ~ 0.999995 for R_out = 13.99, sigma_d = 2.83.
        // THE TOLERANCE MUST BE GRID-AWARE, AND A FIXED 1 % IS NOT.
        // f_bore is the integral of a smooth Gaussian over a SHARP cylinder, so a
        // Cartesian grid resolves that boundary only to a staircase of ~dx/2. The
        // resulting error in f_bore is (df/dR)*dR, and the RELATIVE error is
        // (df/dR)/f * dR — which is NOT constant along the ladder:
        //
        //   rung   R_in   f_bore    (df/dR)/f      1 % means dR <
        //   r10     10    0.9981     0.0024/Bohr      4.1 Bohr   (meaningless)
        //   r08      8    0.9817     0.0187/Bohr      0.54 Bohr
        //   r06      6    0.8946     0.0884/Bohr      0.11 Bohr
        //   r04      4    0.6321     0.2910/Bohr      0.034 Bohr (7 % of a cell!)
        //
        // A factor of 120 in sensitivity. Measured 2026-08-03: r04 came in at
        // -1.80 %, i.e. an effective radius error of 0.062 Bohr = 12 % of ONE grid
        // cell — the SAME dR that reads as 0.01 % at r10 and 0.12 % at r08, both of
        // which passed. So the fixed tolerance was silently tightening as the wall
        // closed in and aborted a perfectly good run.
        //
        // That the packet was fine is provable independently: <r_perp>(0), a SMOOTH
        // moment of the same density, matched the analytic value to 0.016 %. A truly
        // mis-injected packet moves both; only the sharp cut moved.
        //
        // So: floor of 1 % (other numerical effects) in quadrature with the grid
        // geometry term. Gives 1.0 / 1.1 / 2.4 / 7.3 % at r10/r08/r06/r04 — still a
        // real gate everywhere, but one that asks for grid-achievable accuracy.
        auto rayleigh_tol_pc = [&](double R){
            const double a2 = SIGMA_POT*SIGMA_POT;
            const double f  = 1.0 - std::exp(-R*R/(2.0*a2));
            if(f <= 0.0) return 1.0;
            const double dfdR = (R/a2) * std::exp(-R*R/(2.0*a2));
            const double geo  = 100.0 * (dfdR/f) * (SPACING/2.0);
            return std::sqrt(1.0 + geo*geo);
        };
        if(FILLED){
            gate_abs("f_bore(0) = 0 (no bore exists)", o0.f_bore, 0.0, 1e-12);
            gate_rel("f_wall(0) (Rayleigh, filled)", o0.f_wall,
                     1.0 - std::exp(-R_OUT*R_OUT/(2.0*SIGMA_POT*SIGMA_POT)),
                     rayleigh_tol_pc(R_OUT));
        } else {
            gate_rel("f_bore(0) (Rayleigh)", o0.f_bore,
                     1.0 - std::exp(-R_IN*R_IN/(2.0*SIGMA_POT*SIGMA_POT)),
                     rayleigh_tol_pc(R_IN));
        }
        gate_rel("<r_perp>(0) = s_d sqrt(pi/2)", o0.r_mean, SIGMA_POT*std::sqrt(M_PI/2.0), 5.0);
        std::cout << "  [info] max_overlap with the occupied manifold = "
                  << report.max_overlap << " (want < 1e-3)\n";

        // Momentum-space aliasing: at sigma_WP = 4 the packet is 22x narrower in k
        // than at sigma_WP = 0.5, so this is expected to be ~0 %. Reported, not
        // gated, so the number lands in every run log.
        {
            const double sp   = 1.0/(std::sqrt(2.0)*SIGMA_WP);
            const double knyq = M_PI / SPACING;
            const double zsc  = (knyq - K0)/sp;
            const double tail = 0.5*std::erfc(zsc/std::sqrt(2.0));
            std::cout << "  [info] ALIASING: k_Nyq = " << knyq << "  sigma_p = " << sp
                      << "  (k_Nyq-k0)/sigma_p = " << zsc << "  -> " << 100.0*tail
                      << " % of the z-momentum weight beyond Nyquist\n";
            if(tail > 0.02)
                std::cout << "  [WARN] > 2 % aliased: T1 = <p^2>/2 is biased at dx = "
                          << SPACING << ".\n";
        }
        if(fails > 0){
            std::cerr << "\nFATAL: " << fails << " t=0 gate(s) failed — the injected "
                         "packet is not the one this run claims. Aborting.\n";
            return 4;
        }
        std::cout << "  all t=0 gates PASSED\n\n";

        total_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);
        wp_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0);
    }

    // ---- real-time -----------------------------------------------------------
    double e_tot_first = 0.0, e_tot_last = 0.0; bool e_tot_seen = false;
    // var(p) freeze + channeling diagnostics, carried out of the loop for the summary.
    double var_pz_first = 0.0, var_pz_last = 0.0;
    double f_bore_first = 0.0, f_bore_last = 0.0, f_bore_min = 1.0;
    bool   f_bore_seen  = false;
    { auto m_init = wp_mom.compute(electrons); var_pz_first = m_init.sz2; }

    inqkit::RealTimeSession rt(ions, electrons, 1);      // callback EVERY step
    rt.add([&](inqkit::StepContext const& ctx){
        const int step = ctx.step;

        auto n_tot = inqkit::fields::density::total(*ctx.electrons);
        const double l2 = density_delta.snapshot(n_tot, ctx.time_au, step);
        inqkit::StepContext out = ctx; out.density_l2 = l2;
        obs.append(out);

        // pairwise Coulomb ledger with the WP as charge group P (Poisson-linearity
        // form, so the terms close EXACTLY against INQ's own E_hartree/E_external).
        // e_pp is the WP SELF-Hartree — the term with no classical counterpart and
        // the leading suspect for any residual classical/WP discrepancy.
        auto n_wp_f = inqkit::jellium::orbital_density_field(*ctx.electrons, wp_idx);
        auto ct = inqkit::jellium::compute_coulomb_wp(ctx.electrons->density(), n_wp_f, phiplus);

        auto occ = inqkit::observables::radial_occupancy(*ctx.electrons, wp_idx, TUBE_BAND);
        if(!f_bore_seen){ f_bore_first = occ.f_bore; f_bore_seen = true; }
        f_bore_last = occ.f_bore;
        f_bore_min  = std::min(f_bore_min, occ.f_bore);

        if(ctx.electrons->root()){
            ix << step << ',' << ctx.time_au << ',' << ct.e_ss << ',' << ct.e_pp << ','
               << ct.e_ps << ',' << ct.e_sb << ',' << ct.e_pb << ',' << E_BB << ','
               << ct.e_hartree_check << ',' << ct.e_external_check << ','
               << ct.norm_wp << ',' << ct.norm_total << '\n';
            occ_f << step << ',' << ctx.time_au << ',' << occ.f_bore << ',' << occ.f_wall
                  << ',' << occ.f_outside << ',' << occ.r_mean << ',' << occ.r2_mean
                  << ',' << occ.sigma_r << ',' << occ.norm_total << '\n';
        }

        if(!e_tot_seen){ e_tot_first = ctx.energy_total; e_tot_seen = true; }
        e_tot_last = ctx.energy_total;

        if(SAVE_EVERY > 0 && step % SAVE_EVERY == 0){
            total_wr.write(n_tot, ctx.time_au, step);
            wp_wr.write(inqkit::fields::density::orbital(*ctx.electrons, wp_idx), ctx.time_au, step);
        }
        if(WF_EVERY > 0 && step % WF_EVERY == 0){
            auto wf = inqkit::fields::orbital::wavefunction(*ctx.electrons, wp_idx);
            wf_wr.write(wf, "wavefunction_t" + tag6(step));
        }
    });

    const int MAX_CKPT = env_i("CH_MAX_CKPT", 3);
    auto prune_ckpts = [&](){
        if(!electrons.root()) return;
        std::vector<std::filesystem::path> snaps;
        for(auto const& e : std::filesystem::directory_iterator(OUT))
            if(e.is_directory() && e.path().filename().string().rfind("ckpt_step", 0) == 0)
                snaps.push_back(e.path());
        std::sort(snaps.begin(), snaps.end());
        while((int)snaps.size() > MAX_CKPT){
            std::error_code ec; std::filesystem::remove_all(snaps.front(), ec);
            std::cout << "  [ckpt] pruned " << snaps.front().filename().string() << "\n" << std::flush;
            snaps.erase(snaps.begin());
        }
    };
    auto write_rt_state = [&](int last){
        if(!electrons.root()) return;
        std::ofstream st(RT_STATE);
        st << std::setprecision(12)
           << "last_step=" << last << "\ntime_au=" << (last*DT) << "\ndt=" << DT
           << "\nwp_idx=" << wp_idx << "\nk0=" << K0 << "\nsigma_wp=" << SIGMA_WP << "\n";
    };

    auto step_callback = [&](auto const& data){
        rt.step(data);
        wp_mom.maybe_accumulate(data);
        wp_pos.maybe_accumulate(data);
        momentum_dist.maybe_accumulate(data);
        if(data.iter() > 0 && data.iter() % CKPT_EVERY == 0){
            electrons.save(CKPT);
            write_rt_state(data.iter());
            const std::string snap = OUT + "/ckpt_step" + tag6(data.iter());
            electrons.save(snap);
            prune_ckpts();
            std::cout << "  [ckpt] step " << data.iter() << " -> " << snap << "\n" << std::flush;
        }
    };

    // No CAP: the propagation is unitary, so energy_total is conserved and the
    // wavepacket norm is too. That is a far stronger gate than norm monitoring.
    real_time::propagate(ions, electrons, step_callback, options::theory{}.lda(),
                         options::real_time{}.num_steps(N_STEPS).dt(DT * 1.0_atomictime),
                         bg_pert, START);

    // ---- final checkpoint ---------------------------------------------------
    electrons.save(CKPT);
    {
        const std::string final_snap = OUT + "/ckpt_step" + tag6(N_STEPS);
        electrons.save(final_snap);
        prune_ckpts();
        std::cout << "  [ckpt] FINAL step " << N_STEPS << " -> " << final_snap
                  << "  (t = " << (N_STEPS*DT) << " a.u.)\n" << std::flush;
    }
    write_rt_state(N_STEPS);

    // ---- end-of-run diagnostics ---------------------------------------------
    { auto mF = wp_mom.compute(electrons); var_pz_last = mF.sz2; }
    const double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_wall0).count();
    const double e_drift_ev = (e_tot_last - e_tot_first) * HA;
    {
        const bool ok = std::abs(e_drift_ev) < 1.0e-3;
        std::cout << (ok ? "\n  [PASS] " : "\n  [WARN] ")
                  << "energy conservation: E_total drift = " << e_drift_ev
                  << " eV over " << (N_STEPS - START) << " steps"
                  << " (want < 1e-3 eV; no CAP => H is Hermitian and t-independent)\n";
        std::cout << "  [info] var(p_z): " << var_pz_first << " -> " << var_pz_last
                  << "  (" << 100.0*(var_pz_last-var_pz_first)/var_pz_first
                  << " %; free-evolution value " << sigma_p2 << ")\n"
                  << "  [info] f_bore: " << f_bore_first << " -> " << f_bore_last
                  << "  (min over the run " << f_bore_min << ")\n"
                  << "         a FROZEN var(p_z) and f_bore ~ 1 are the two conditions "
                     "under which -dT2/ds is the classical stopping power.\n";
    }

    if(electrons.root()){
        ix.close(); occ_f.close();
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(14)
          << "run = cylindrical_jellium/proximity_ladder/wp/" << rung->label << "/" << env_s("CH_OUT","wp") << "\n"
          << "rung = " << rung->label << "\n"
          << "run_type = wavepacket projectile, annular jellium tube TDDFT (ALDA)\n"
          << "campaign = cylindrical channeling KS-orbital stopping (twin)\n"
          << "plan = docs/plans/cylindrical-proximity-ladder.md\n"
          << "twin_role = wavepacket\ntwin_partner = ../classical/results/"
          << env_s("CH_OUT_CL","classical") << "\n"
          << "engine = inq-study\nxc = LDA (ALDA in TDDFT)\n"
          << "representation = wavepacket\n"
          // Single-token `projectile` value; check_twin.py compares it against the
          // classical twin's to prove the two runs really differ.
          << "projectile = wavepacket_orbital\n"
          << "projectile_detail = occupied KS orbital, Gaussian envelope with drift k0\n"
          << "geometry = annular_tube\n"
          << "r_in_bohr = " << R_IN << "  r_out_bohr = " << R_OUT
          << "  edge_width_bohr = " << Cfg::EDGE_W_BOHR << "  tube_axis = " << Cfg::TUBE_AXIS << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR << "\n"
          << "Lz = " << Cfg::LZ_BOHR << "\nperiodicity = 3\nspacing = " << SPACING
          << "\nspacing_bohr = " << SPACING << "\n"
          << "N = " << N_ELEC << "\nn_electrons = " << N_ELEC << "  (bath; +1 for the WP)\n"
          << "n0_a0m3 = " << N0 << "\nr_s = " << RS << "\n"
          << "omega_p_au = " << OMEGA_P << "  omega_p_ev = " << OMEGA_P*HA
          << "  lambda_p_bohr = " << 2.0*M_PI*K0/OMEGA_P << "\n"
          << "v_fermi = " << VF << "  v_over_vf = " << K0/VF << "\n"
          << "n_states = " << n_states << "  extra_states = " << EXTRA_ST << "\n"
          << "wp_enabled = yes\nwp_state_index = " << wp_idx << "\n"
          << "sigma_WP = " << SIGMA_WP << "\nsigma_wp = " << SIGMA_WP
          << "\nwp_sigma_bohr = " << SIGMA_WP << "\n"
          << "wp_sigma_density = " << SIGMA_POT << "\nsigma_pot = " << SIGMA_POT << "\n"
          << "wp_sigma_note = the run is LABELLED by sigma_WP; the density std "
             "sigma_WP/sqrt2 equals the twin's classical sigma_pot\n"
          << "wp_k0_bohr_inv = " << K0 << "  wp_velocity = " << K0
          << "  k0 = " << K0 << "  v0 = " << K0 << "\n"
          << "wp_drift_energy_ev = " << 0.5*K0*K0*HA << "\n"
          << "projectile_energy_ev = " << 0.5*K0*K0*HA << "\n"
          << "wp_localisation_energy_ev = " << 3.0/(4.0*SIGMA_WP*SIGMA_WP)*HA << "\n"
          << "launch_z = " << LAUNCH_Z << "\n"
          << "norm_after = " << report.norm_after << "  max_overlap = " << report.max_overlap << "\n"
          << "cap = off\ncap_enabled = no\n"
          << "cap_note = no absorbing potential; H is Hermitian and t-independent so "
             "energy_total is CONSERVED and is the correctness gate\n"
          << "start_step = " << START << "  rt_num_steps = " << N_STEPS
          << "  n_steps = " << N_STEPS << "  dt_au = " << DT << "  dt = " << DT
          << "  total_time_au = " << (DT*N_STEPS) << "\n"
          << "save_every = " << SAVE_EVERY << "  wf_every = " << WF_EVERY
          << "  stats_every = 1  ckpt_every = " << CKPT_EVERY << "\n"
          << "energy_total_first_ha = " << e_tot_first << "\n"
          << "energy_total_last_ha = " << e_tot_last << "\n"
          << "energy_total_drift_ev = " << e_drift_ev << "\n"
          << "var_pz_free_value = " << sigma_p2 << "\n"
          << "var_pz_first = " << var_pz_first << "\nvar_pz_final = " << var_pz_last << "\n"
          << "var_pz_growth_pct = " << 100.0*(var_pz_last-var_pz_first)/var_pz_first << "\n"
          << "f_bore_first = " << f_bore_first << "\nf_bore_final = " << f_bore_last
          << "\nf_bore_min = " << f_bore_min << "\n"
          << "spread_rate_bohr_per_au = " << 1.0/(std::sqrt(2.0)*SIGMA_WP) << "\n"
          << "t_2sigma_reaches_bore_au = "
          << std::sqrt(2.0*SIGMA_WP*SIGMA_WP*((R_IN/2.0)*(R_IN/2.0)
                                              - SIGMA_WP*SIGMA_WP/2.0)) << "\n"
          << "t_transverse_overlap_au = "
          << std::sqrt(2.0*SIGMA_WP*SIGMA_WP*((Cfg::LX_BOHR/6.0)*(Cfg::LX_BOHR/6.0)
                                              - SIGMA_WP*SIGMA_WP/2.0))
          << "  (6 sigma_d = L_xy; AFTER the end of the run)\n"
          << "analysis_note = fit S over the window where f_bore stays near 1; later "
             "steps are recorded but the packet has begun to touch the wall\n"
          << "date_finished = " << iso_now() << "\nwall_time_s = " << wall << "\n"
          << "gs_dir = " << GS_DIR << "\nrun_completed = true\n";
    }
    std::cout << "\nDone. Wall time " << wall << " s.\n";
    return 0;
}
