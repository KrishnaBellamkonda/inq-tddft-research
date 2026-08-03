// ============================================================================
// cylindrical_jellium / scripts/channeling_sic / wp / run.cpp
//
// SELF-INTERACTION-CORRECTED re-run of the channeling twin's QUANTUM half.
// Plan: docs/plans/wp-self-interaction-correction.md (reviewed 2026-08-02).
// Identical in EVERY physical parameter to scripts/channeling_twin/wp/run.cpp
// (same shared config header, same GS, same dt/N_STEPS/cadences/observables);
// the ONLY addition is the per-step SIC kick on the wavepacket orbital:
//
//     psi_wp <- N . Q . exp(+i dt_eff (v_H[n_wp] [+ v_xc^unpol[n_wp]])) psi_wp
//
// applied through inqkit::SelfInteractionCorrection (CH_SIC = pzrun | h |
// none; none reproduces the uncorrected twin bit-for-bit, the regression
// control). Strang scheduling: half-kick at the step-0 callback, full kick
// after every interior step, half-kick at the final step.
//
// WHAT THIS RUN MEASURES (plan §5). The completed uncorrected twin showed the
// WP stopping 20 % below its classical twin, with the deficit driven by the
// LDA self-interaction blowing the packet apart (E_PP released 1.64 eV; 1.47x
// excess expansion; impulse ratio tracking f_bore at r = +0.98). Removing the
// self-interaction splits that deficit into its SIE-artefact part (removed)
// and its genuine quantum-kinematic part (remains; Nazarov & Gross 2025,
// arXiv:2510.26222, predict a nonzero one). Three-way comparison: classical /
// WP uncorrected / WP+SIC, all on the same GS.
//
// CONSERVATION BOOKKEEPING (plan §0/D2 — read before gating!). E_total is NO
// LONGER the conserved quantity; the corrected functional is
//     E_corr = E_KS - U[n_wp]                    (CH_SIC=h)
//     E_corr = E_KS - U[n_wp] - E_xc^unpol[n_wp] (CH_SIC=pzrun)
// and even E_corr is conserved only APPROXIMATELY in the jellium: the
// projected (one-sided Lagrange-multiplier) scheme leaves the drift channel
// dE/dt = 2 Im sum_j <wp|h|j><j|v_SIC|wp>, which vanishes in vacuum but not
// where occupied bath orbitals overlap the packet. It is MEASURED (sic.csv)
// and soft-gated, never assumed zero. Orthogonality, by contrast, is EXACT:
// max_j |<psi_j|psi_wp>| after each step's projection is the hard invariant.
//
// TIERS (CH_SIC_TIER): "b" = 200-step bath-integrity check with HARD gates
// (cum_norm_removed < 1e-3, max_overlap_pre < 1e-3, |dE_corr| < 0.02 eV);
// "prod" = the 1500-step production run; gates become WARN-and-report
// (.claude/rules/checkpoint-dont-block.md: a completed-tiers production run
// never self-blocks; it reports).
//
// RESUME (.claude/rules/final-timestep-checkpoint.md). Interior checkpoints
// sit at FULL-kick boundaries (the callback kicks before the ckpt save), so a
// resumed segment continues the Strang pattern with no action. The FINAL
// checkpoint sits after the CLOSING half-kick; rt_state.txt records
// sic_boundary=closed and the resume branch applies the compensating opening
// half-kick after electrons.load().
//
// Everything else (circular centroid, minimum_image injection, t=0 analytic
// gates, no CAP, interactions.csv ledger, observable cadences) is inherited
// verbatim from the twin — see scripts/channeling_twin/wp/run.cpp for the
// full rationale of each.
//
// Env: CH_GS_DIR(REQUIRED) CH_OUT(wp_sic) CH_N(160) CH_V0 CH_SIGMA CH_LAUNCH_Z
//      CH_N_STEPS CH_DT CH_SPACING CH_SAVE_EVERY CH_WF_EVERY
//      CH_CKPT_EVERY(0=auto N/3) CH_MAX_CKPT(3) CH_RESUME(0)
//      CH_SIC(pzrun) CH_SIC_TIER(prod)
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
#include <inqkit/wavepacket/self_interaction_correction.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/interaction_energies.hpp>
#include <inqkit/jellium/analytics.hpp>

#include "../../../shared/configs/channeling_tube_rs3.hpp"

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
using Cfg = cylindrical_jellium::config::ChannelingTubeRs3;
using SIC = inqkit::SelfInteractionCorrection;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

static double read_state_d(const std::string& path, const char* key, double def){
    std::ifstream f(path); std::string line; const std::string k = std::string(key) + "=";
    while(std::getline(f, line)){ auto p = line.find(k);
        if(p != std::string::npos) return std::atof(line.substr(p + k.size()).c_str()); }
    return def;
}
static std::string read_state_s(const std::string& path, const char* key, const std::string& def){
    std::ifstream f(path); std::string line; const std::string k = std::string(key) + "=";
    while(std::getline(f, line)){ auto p = line.find(k);
        if(p != std::string::npos){
            std::string v = line.substr(p + k.size());
            while(!v.empty() && (v.back()=='\n' || v.back()=='\r' || v.back()==' ')) v.pop_back();
            return v; } }
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

    const int    N_ELEC   = env_i("CH_N",         Cfg::N_ELECTRONS);
    const double SPACING  = env_d("CH_SPACING",   Cfg::SPACING_BOHR);
    const double SIGMA_WP = env_d("CH_SIGMA",     Cfg::SIGMA_WP_BOHR);
    const double K0       = env_d("CH_V0",        Cfg::PROJ_V0);
    const double LAUNCH_Z = env_d("CH_LAUNCH_Z",  Cfg::LAUNCH_Z_BOHR);
    const double DT       = env_d("CH_DT",        Cfg::DT_AU);
    const int    N_STEPS  = env_i("CH_N_STEPS",   Cfg::N_STEPS);
    const int    SAVE_EVERY = env_i("CH_SAVE_EVERY", Cfg::SAVE_EVERY);
    const int    WF_EVERY   = env_i("CH_WF_EVERY",   Cfg::WF_EVERY);
    const bool   RESUME     = env_i("CH_RESUME", 0) != 0;
    const std::string OUT    = "results/" + env_s("CH_OUT", "wp_sic");
    const std::string GS_DIR = env_s("CH_GS_DIR", "");
    const std::string SICSTR = env_s("CH_SIC", "pzrun");
    const std::string TIER   = env_s("CH_SIC_TIER", "prod");
    const auto SIC_MODE = SIC::mode_from_string(SICSTR);
    const bool tier_b = (TIER == "b" || TIER == "B");

    int CKPT_EVERY = env_i("CH_CKPT_EVERY", 0);
    if(CKPT_EVERY <= 0) CKPT_EVERY = std::max(1, N_STEPS / 3);

    if(GS_DIR.empty() || !std::filesystem::exists(GS_DIR)){
        std::cerr << "FATAL: CH_GS_DIR missing or unset: '" << GS_DIR << "'\n"; return 2;
    }

    const double SIGMA_POT = SIGMA_WP / std::sqrt(2.0);
    const double sigma_p2  = 1.0 / (2.0 * SIGMA_WP * SIGMA_WP);
    const double V_ANN     = Cfg::v_annulus();
    const double N0        = double(N_ELEC) / V_ANN;
    const double RS        = inqkit::jellium::rs_from_n0(N0);
    const double OMEGA_P   = std::sqrt(4.0*M_PI*N0);
    const double VF        = inqkit::jellium::k_fermi_n0(N0);

    // ---- resume state -------------------------------------------------------
    const std::string CKPT = OUT + "/checkpoint", RT_STATE = OUT + "/rt_state.txt";
    int START = 0, wp_idx = -1;
    std::string resume_boundary = "full";
    if(RESUME){
        START  = (int)read_state_d(RT_STATE, "last_step", -1);
        wp_idx = (int)read_state_d(RT_STATE, "wp_idx", -1);
        resume_boundary = read_state_s(RT_STATE, "sic_boundary", "full");
        const std::string prev_mode = read_state_s(RT_STATE, "sic_mode", "?");
        if(prev_mode != "?" && prev_mode != SIC::mode_name(SIC_MODE)){
            std::cerr << "FATAL: resume with CH_SIC=" << SIC::mode_name(SIC_MODE)
                      << " but the checkpoint was written with sic_mode=" << prev_mode
                      << " — a mixed-mode trajectory is not the run it claims.\n"; return 2; }
        if(START < 0 || wp_idx < 0){
            std::cerr << "FATAL: CH_RESUME=1 but no readable " << RT_STATE << "\n"; return 2; }
        if(START >= N_STEPS){
            std::cout << "Already at/after target (" << START << " >= " << N_STEPS
                      << "); nothing to do.\n"; return 0; }
    }
    const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

    std::cout << std::setprecision(10)
              << "\n=== channeling_sic/wp  OUT=" << OUT << " ===\n"
              << "  SIC       = " << SIC::mode_name(SIC_MODE) << "   tier = " << TIER << "\n"
              << "  tube      = R_in " << Cfg::R_IN_BOHR << "  R_out " << Cfg::R_OUT_BOHR
              << "  L_z " << Cfg::LZ_BOHR << "  (fully periodic)\n"
              << "  cell      = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR << " x "
              << Cfg::LZ_BOHR << "  dx = " << SPACING << "\n"
              << "  bath      = " << N_ELEC << " e   n0 = " << N0 << "   r_s = " << RS
              << "   omega_p = " << OMEGA_P*HA << " eV   v_F = " << VF << "\n"
              << "  WP        = sigma_WP " << SIGMA_WP << " (density std " << SIGMA_POT
              << ")  k0 = " << K0 << "  E_drift = " << 0.5*K0*K0*HA << " eV  v/v_F = " << K0/VF << "\n"
              << "  launch    = (0,0," << LAUNCH_Z << ")\n"
              << "  dt = " << DT << "  START = " << START << " -> N_STEPS = " << N_STEPS
              << "  t_total = " << DT*N_STEPS << " a.u." << (RESUME ? "  [RESUME]" : "") << "\n"
              << "  cadence   : density/" << SAVE_EVERY << "  wavefn/" << WF_EVERY
              << "  stats/1  ckpt/" << CKPT_EVERY << "\n"
              << "  GS        = " << GS_DIR << "\n\n";

    // ---- system -------------------------------------------------------------
    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b,
                                            Cfg::LY_BOHR * 1.0_b,
                                            Cfg::LZ_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
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
    for(auto const& d : {OBS, VTI + "/density_total", VTI + "/density_wp",
                         VTI + "/density_delta", VTI + "/density_delta_coarse",
                         VTI + "/density_gs_system", VTI + "/wavefunction_wp"})
        std::filesystem::create_directories(d);

    inqkit::InjectionReport report{};

    if(RESUME){
        electrons.load(CKPT);
        std::cout << "  RESUMED from step " << START << " (wp_idx = " << wp_idx
                  << ", sic_boundary = " << resume_boundary << ")\n";
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

        inqkit::io::RealField3DLayout lay{
            .field_name = "density", .include_meta = false, .emit_raw = false,
            .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
        inqkit::io::RealField3DWriter gs_wr(VTI + "/density_gs_system", lay, {.overwrite = true});
        gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system");

        // minimum_image mandatory: launch is 0.71 density-sigma from the -z face
        // (see channeling_twin/wp/run.cpp for the 15x-var(p_z) failure this avoids).
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
          << "max_overlap    = " << report.max_overlap << "\n"
          << "sic_mode       = " << SIC::mode_name(SIC_MODE) << "\n";
    }

    // ---- SIC ----------------------------------------------------------------
    const SIC sic(SIC_MODE, wp_idx);
    const bool sic_on = (SIC_MODE != SIC::Mode::none);
    if(sic_on && RESUME && resume_boundary == "closed"){
        // The final checkpoint of a completed segment holds the state AFTER the
        // closing half-kick; re-open the Strang pattern with the compensating
        // opening half-kick before propagating further.
        auto r0 = sic.apply(electrons, 0.5*DT);
        std::cout << "  [sic] applied OPENING half-kick on resume from closed boundary"
                  << " (u_self = " << r0.u_self*HA << " eV)\n";
    }

    // ---- background + constant fields for the ledger ------------------------
    inqkit::jellium::localised_background_params bg;
    bg.shape        = inqkit::jellium::background_shape::annulus;
    bg.n0           = N0;
    bg.half_width   = Cfg::R_OUT_BOHR;
    bg.inner_radius = Cfg::R_IN_BOHR;
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
    inqkit::observables::MomentumDistribution momentum_dist(
        OBS + "/momentum_distribution" + SEG + ".csv", wp_idx, Cfg::LZ_BOHR,
        {.n_bins = 128, .k_max_bohr_inv = 0.0,
         .write_every = (SAVE_EVERY > 0 ? SAVE_EVERY : 1)});

    inqkit::observables::DensityDelta density_delta(
        VTI + "/density_delta", VTI + "/density_delta_coarse",
        {.emit_raw_vti = true, .emit_coarse_vti = true,
         .compute_l2 = true, .coarse_bin_bohr = 3.0,
         .emit_every = (SAVE_EVERY > 0 ? SAVE_EVERY : 1)});

    const inqkit::observables::RadialBand TUBE_BAND{
        .axis = Cfg::TUBE_AXIS,
        .center = {Cfg::CENTER_X, Cfg::CENTER_Y, Cfg::CENTER_Z},
        .r_inner = Cfg::R_IN_BOHR, .r_outer = Cfg::R_OUT_BOHR};

    std::ofstream occ_f, ix, sc;
    if(electrons.root()){
        occ_f.open(OBS + "/wp_radial_occupancy" + SEG + ".csv");
        occ_f << std::setprecision(12)
              << "# tube: axis=" << Cfg::TUBE_AXIS << " R_in=" << Cfg::R_IN_BOHR
              << " R_out=" << Cfg::R_OUT_BOHR << "\n"
              << "step,time_au,f_bore,f_wall,f_outside,r_mean,r2_mean,sigma_r,norm_total\n";

        ix.open(OBS + "/interactions" + SEG + ".csv");
        ix << std::setprecision(12)
           << "step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,"
              "e_hartree_check,e_external_check,norm_proj,norm_electrons\n";

        // The SIC ledger (plan §2). E_PP in interactions.csv stays as the
        // packet-size diagnostic; u_self here is the SAME integral evaluated
        // by the correction itself (they must track each other — a free
        // cross-check between two independent code paths).
        sc.open(OBS + "/sic" + SEG + ".csv");
        sc << std::setprecision(12)
           << "step,time_au,u_self_ha,exc_self_ha,e_corrected_ha,"
              "max_overlap_pre,norm_removed,cum_norm_removed\n";
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
        gate_abs("centroid z (CIRCULAR)",       r0.zc,  LAUNCH_Z, 0.05);
        gate_rel("spread sigma_z (CIRCULAR)",   r0.szc, SIGMA_POT, 5.0);
        gate_rel("transverse std = s/sqrt2",    std::sqrt(r0.sx2), SIGMA_POT, 5.0);
        gate_rel("f_bore(0) (Rayleigh)", o0.f_bore,
                 1.0 - std::exp(-Cfg::R_IN_BOHR*Cfg::R_IN_BOHR/(2.0*SIGMA_POT*SIGMA_POT)), 1.0);
        gate_rel("<r_perp>(0) = s_d sqrt(pi/2)", o0.r_mean, SIGMA_POT*std::sqrt(M_PI/2.0), 5.0);
        std::cout << "  [info] max_overlap with the occupied manifold = "
                  << report.max_overlap << " (want < 1e-3)\n";
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
    double var_pz_first = 0.0, var_pz_last = 0.0;
    double f_bore_first = 0.0, f_bore_last = 0.0, f_bore_min = 1.0;
    bool   f_bore_seen  = false;
    { auto m_init = wp_mom.compute(electrons); var_pz_first = m_init.sz2; }

    // SIC accumulators for the summary + gates
    double e_corr_first = 0.0, e_corr_last = 0.0, e_corr_mid = 0.0;
    double cum_removed = 0.0, max_ov_run = 0.0;
    bool   e_corr_seen = false;
    std::string sic_boundary = "full";   // what the NEXT rt_state write records

    inqkit::RealTimeSession rt(ions, electrons, 1);
    rt.add([&](inqkit::StepContext const& ctx){
        const int step = ctx.step;

        auto n_tot = inqkit::fields::density::total(*ctx.electrons);
        const double l2 = density_delta.snapshot(n_tot, ctx.time_au, step);
        inqkit::StepContext out = ctx; out.density_l2 = l2;
        obs.append(out);

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
           << "\nwp_idx=" << wp_idx << "\nk0=" << K0 << "\nsigma_wp=" << SIGMA_WP
           << "\nsic_mode=" << SIC::mode_name(SIC_MODE)
           << "\nsic_boundary=" << sic_boundary << "\n";
    };

    auto step_callback = [&](auto const& data){
        rt.step(data);
        wp_mom.maybe_accumulate(data);
        wp_pos.maybe_accumulate(data);
        momentum_dist.maybe_accumulate(data);

        // ---- SIC kick, AFTER every observable (uniform PRE-kick convention,
        // plan §0) and BEFORE any checkpoint save (so interior checkpoints sit
        // at full-kick boundaries and resume needs no action).
        const int step = (int)data.iter();
        if(sic_on && (step == 0 || step > START)){
            const double f = (step == 0 || step == N_STEPS) ? 0.5 : 1.0;
            auto srep = sic.apply(electrons, f * DT);
            const double e_corr = sic.corrected_energy(data.energy().total(), srep);
            cum_removed += srep.norm_removed;
            max_ov_run = std::max(max_ov_run, srep.max_overlap_pre);
            if(!e_corr_seen){ e_corr_first = e_corr; e_corr_seen = true; }
            if(step == (START + N_STEPS)/2) e_corr_mid = e_corr;
            e_corr_last = e_corr;
            sic_boundary = (step == N_STEPS) ? "closed" : "full";
            if(electrons.root())
                sc << step << ',' << data.time() << ',' << srep.u_self << ','
                   << srep.exc_self << ',' << e_corr << ',' << srep.max_overlap_pre
                   << ',' << srep.norm_removed << ',' << cum_removed << '\n';
        }

        if(data.iter() > 0 && data.iter() % CKPT_EVERY == 0){
            electrons.save(CKPT);
            write_rt_state(data.iter());
            const std::string snap = OUT + "/ckpt_step" + tag6(data.iter());
            electrons.save(snap);
            prune_ckpts();
            std::cout << "  [ckpt] step " << data.iter() << " -> " << snap << "\n" << std::flush;
        }
    };

    real_time::propagate(ions, electrons, step_callback, options::theory{}.lda(),
                         options::real_time{}.num_steps(N_STEPS).dt(DT * 1.0_atomictime),
                         bg_pert, START);

    // ---- final checkpoint (post closing half-kick => boundary "closed") ------
    electrons.save(CKPT);
    {
        const std::string final_snap = OUT + "/ckpt_step" + tag6(N_STEPS);
        electrons.save(final_snap);
        prune_ckpts();
        std::cout << "  [ckpt] FINAL step " << N_STEPS << " -> " << final_snap
                  << "  (t = " << (N_STEPS*DT) << " a.u.)\n" << std::flush;
    }
    write_rt_state(N_STEPS);

    // ---- end-of-run diagnostics ----------------------------------------------
    { auto mF = wp_mom.compute(electrons); var_pz_last = mF.sz2; }
    const double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_wall0).count();
    const double e_drift_ev = (e_tot_last - e_tot_first) * HA;
    const double e_corr_drift_ev = sic_on ? (e_corr_last - e_corr_first) * HA : 0.0;

    int sic_fails = 0;
    if(sic_on){
        // E_total is NOT conserved by design (the kick does work on the WP);
        // report it as information, not a gate. E_corr is the soft gate.
        std::cout << "\n  [info] E_total drift = " << e_drift_ev
                  << " eV (NOT a gate under SIC — the correction removes energy "
                     "terms from the WP's Hamiltonian)\n";
        auto sic_gate = [&](char const* nm, double got, double tol){
            const bool ok = std::abs(got) <= tol;
            if(ok)           std::cout << "  [PASS] " << nm << ": " << got << "  (tol " << tol << ")\n";
            else if(tier_b){ std::cout << "  [FAIL] " << nm << ": " << got << "  (tol " << tol << ")\n"; ++sic_fails; }
            else             std::cout << "  [WARN] " << nm << ": " << got << "  (tol " << tol
                                       << "; production reports, never self-blocks)\n";
        };
        std::cout << "  --- SIC gates (" << (tier_b ? "TIER B, hard" : "production, soft") << ") ---\n";
        sic_gate("cum_norm_removed", cum_removed, 1e-3);
        sic_gate("max_overlap_pre (run max)", max_ov_run, 1e-3);
        sic_gate("E_corrected drift (eV)", e_corr_drift_ev,
                 0.1 * double(N_STEPS - START) / 1500.0);
        std::cout << "  [info] E_corr first/mid/last (Ha) = " << e_corr_first << " / "
                  << e_corr_mid << " / " << e_corr_last
                  << "   (secularity read from sic.csv in analysis)\n";
    } else {
        const bool ok = std::abs(e_drift_ev) < 1.0e-3;
        std::cout << (ok ? "\n  [PASS] " : "\n  [WARN] ")
                  << "energy conservation: E_total drift = " << e_drift_ev
                  << " eV over " << (N_STEPS - START) << " steps (want < 1e-3 eV)\n";
    }
    std::cout << "  [info] var(p_z): " << var_pz_first << " -> " << var_pz_last
              << "  (" << 100.0*(var_pz_last-var_pz_first)/var_pz_first
              << " %; free-evolution value " << sigma_p2 << ")\n"
              << "  [info] f_bore: " << f_bore_first << " -> " << f_bore_last
              << "  (min over the run " << f_bore_min << ")\n";

    if(electrons.root()){
        ix.close(); occ_f.close(); sc.close();
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(14)
          << "run = cylindrical_jellium/channeling_sic/wp/" << env_s("CH_OUT","wp_sic") << "\n"
          << "run_type = SIC-corrected wavepacket projectile, annular jellium tube TDDFT (ALDA)\n"
          << "campaign = wp self-interaction correction (plan docs/plans/wp-self-interaction-correction.md)\n"
          << "plan = docs/plans/wp-self-interaction-correction.md\n"
          << "twin_role = wavepacket_sic\n"
          << "twin_partner = ../../channeling_twin/classical/results/classical\n"
          << "uncorrected_sibling = ../../channeling_twin/wp/results/wp\n"
          << "engine = inq-study\nxc = LDA (ALDA in TDDFT)\n"
          << "representation = wavepacket\n"
          << "projectile = wavepacket_orbital_sic\n"
          << "projectile_detail = occupied KS orbital, Gaussian envelope with drift k0, "
             "per-step projected SIC kick\n"
          << "sic_mode = " << SIC::mode_name(SIC_MODE) << "\n"
          << "sic_tier = " << TIER << "\n"
          << "geometry = annular_tube\n"
          << "r_in_bohr = " << Cfg::R_IN_BOHR << "  r_out_bohr = " << Cfg::R_OUT_BOHR
          << "  edge_width_bohr = " << Cfg::EDGE_W_BOHR << "  tube_axis = " << Cfg::TUBE_AXIS << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR << "\n"
          << "Lz = " << Cfg::LZ_BOHR << "\nperiodicity = 3\nspacing = " << SPACING
          << "\nspacing_bohr = " << SPACING << "\n"
          << "N = " << N_ELEC << "\nn_electrons = " << N_ELEC << "  (bath; +1 for the WP)\n"
          << "n0_a0m3 = " << N0 << "\nr_s = " << RS << "\n"
          << "omega_p_au = " << OMEGA_P << "  omega_p_ev = " << OMEGA_P*HA
          << "  lambda_p_bohr = " << 2.0*M_PI*K0/OMEGA_P << "\n"
          << "v_fermi = " << VF << "  v_over_vf = " << K0/VF << "\n"
          << "n_states = " << n_states << "  extra_states = " << Cfg::EXTRA_STATES << "\n"
          << "wp_enabled = yes\nwp_state_index = " << wp_idx << "\n"
          << "sigma_WP = " << SIGMA_WP << "\nsigma_wp = " << SIGMA_WP
          << "\nwp_sigma_bohr = " << SIGMA_WP << "\n"
          << "wp_sigma_density = " << SIGMA_POT << "\nsigma_pot = " << SIGMA_POT << "\n"
          << "wp_k0_bohr_inv = " << K0 << "  wp_velocity = " << K0
          << "  k0 = " << K0 << "  v0 = " << K0 << "\n"
          << "wp_drift_energy_ev = " << 0.5*K0*K0*HA << "\n"
          << "projectile_energy_ev = " << 0.5*K0*K0*HA << "\n"
          << "launch_z = " << LAUNCH_Z << "\n"
          << "norm_after = " << report.norm_after << "  max_overlap = " << report.max_overlap << "\n"
          << "cap = off\ncap_enabled = no\n"
          << "start_step = " << START << "  rt_num_steps = " << N_STEPS
          << "  n_steps = " << N_STEPS << "  dt_au = " << DT << "  dt = " << DT
          << "  total_time_au = " << (DT*N_STEPS) << "\n"
          << "save_every = " << SAVE_EVERY << "  wf_every = " << WF_EVERY
          << "  stats_every = 1  ckpt_every = " << CKPT_EVERY << "\n"
          << "energy_total_first_ha = " << e_tot_first << "\n"
          << "energy_total_last_ha = " << e_tot_last << "\n"
          << "energy_total_drift_ev = " << e_drift_ev << "\n"
          << "energy_total_drift_note = NOT conserved by design under SIC; the "
             "conserved quantity is e_corrected (sic.csv), see plan section 0/D2\n"
          << "sic_e_corrected_first_ha = " << e_corr_first << "\n"
          << "sic_e_corrected_mid_ha = " << e_corr_mid << "\n"
          << "sic_e_corrected_last_ha = " << e_corr_last << "\n"
          << "sic_e_corrected_drift_ev = " << e_corr_drift_ev << "\n"
          << "sic_cum_norm_removed = " << cum_removed << "\n"
          << "sic_max_overlap_pre = " << max_ov_run << "\n"
          << "sic_gate_failures = " << sic_fails << "\n"
          << "var_pz_free_value = " << sigma_p2 << "\n"
          << "var_pz_first = " << var_pz_first << "\nvar_pz_final = " << var_pz_last << "\n"
          << "var_pz_growth_pct = " << 100.0*(var_pz_last-var_pz_first)/var_pz_first << "\n"
          << "f_bore_first = " << f_bore_first << "\nf_bore_final = " << f_bore_last
          << "\nf_bore_min = " << f_bore_min << "\n"
          << "analysis_note = compare against BOTH channeling_twin halves: classical "
             "(unchanged reference) and uncorrected wp (isolates the SIE artefact)\n"
          << "date_finished = " << iso_now() << "\nwall_time_s = " << wall << "\n"
          << "gs_dir = " << GS_DIR << "\nrun_completed = true\n";
    }
    if(sic_fails > 0){
        std::cerr << "\nFATAL: " << sic_fails << " Tier-B SIC gate(s) failed — do not "
                     "run production (plan §4).\n";
        return 4;
    }
    std::cout << "\nDone. Wall time " << wall << " s.\n";
    return 0;
}
