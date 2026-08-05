// ============================================================================
// localised_jellium / scripts/lz_bulk_sweep / wp / run.cpp
//
// WAVEPACKET half of the slab->bulk L_slab sweep.
// Plan: docs/plans/jellium-slab-extend-Lz.md
//
// Cloned from scripts/sigma56_sv/wp/run.cpp (the validated sigma = 5/6
// production binary) with exactly one structural change: the geometry moves
// from the compile-time Cfg to the RUNTIME box preset (env LZB_CFG,
// shared/configs/lzb_boxes.hpp), so ONE binary serves all four boxes:
//
//   preset     L_z   L_slab  N_e  launch z   family
//   s0p5_L15    75     15     60    -19.0    sigma = 0.5 (standoff 11.5)
//   s0p5_L35    95     35    140    -29.0    sigma = 0.5
//   s5p0_L15    95     15     60    -22.5    sigma = 5   (standoff 15)
//   s5p0_L35   115     35    140    -32.5    sigma = 5
//
// Every t=0 gate, the orthogonalisation budget, the pairwise ledger, the
// checkpoint/resume block and the CAP wiring are UNCHANGED from sigma56_sv, so
// nothing validated there needs re-validating. The CAP fractions are derived
// from the preset's L_z, so they follow the box automatically; the launch
// standoff is per sigma FAMILY so each family's runs match its existing
// L_slab = 25 anchor (the user's comparability rule, 2026-08-05).
//
// See the sigma56_sv wp/run.cpp header for the full rationale of the CAP
// (orbitals wrap on the FFT grid regardless of periodicity(2)), the width
// mapping sigma_pot = sigma_WP/sqrt(2), the E_absorbed-vs-ledger discussion and
// the checkpoint rules — all of it applies verbatim here.
//
// Env: LZB_CFG(REQUIRED) LJ_K0 LJ_N_STEPS LJ_DT(0.04) LJ_SIGMA(family default)
//      LJ_LAUNCH_Z(preset) LJ_CAP_ETA(-1.0) LJ_CAP_L(12.5) LJ_SAVE_EVERY
//      LJ_WF_EVERY LJ_CKPT_EVERY(0=auto N/5) LJ_RESUME(0) LJ_OUT(REQUIRED)
//      LJ_GS_DIR(REQUIRED)
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
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/injection_report.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/interaction_energies.hpp>
#include <inqkit/jellium/analytics.hpp>

#include "../../../shared/configs/lzb_boxes.hpp"

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

using namespace inq;
using namespace inq::magnitude;
using Sh = localised_jellium::config::LzbShared;

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

    const auto B = localised_jellium::config::lzb_box_from_env();

    // ---- parameters (geometry pinned to the LZB_CFG preset) -----------------
    const double SIGMA_WP = env_d("LJ_SIGMA", B.SIGMA_DEFAULT);
    const double SPACING  = env_d("LJ_SPACING", Sh::SPACING_BOHR);
    const double K0       = env_d("LJ_K0", 3.0);                 // pilot velocity
    // Per-FAMILY launch: both halves of a family and every L_slab share it, so
    // arrival width depends on (sigma, v) only — the comparability rule.
    const double LAUNCH_Z = env_d("LJ_LAUNCH_Z", B.LAUNCH_Z);
    const double DT       = env_d("LJ_DT", 0.04);
    const int    N_STEPS  = env_i("LJ_N_STEPS", 2543);
    const double CAP_ETA  = env_d("LJ_CAP_ETA", Sh::CAP_ETA_HA); // Ha, <0 absorbs
    const double CAP_L    = env_d("LJ_CAP_L", Sh::CAP_L_BOHR);   // Bohr per face
    const int    SAVE_EVERY = env_i("LJ_SAVE_EVERY", 25);        // ~100 frames
    const int    WF_EVERY   = env_i("LJ_WF_EVERY", 64);          // ~40 wavefunctions
    const bool   RESUME     = env_i("LJ_RESUME", 0) != 0;
    const std::string OUT    = "results/" + env_s("LJ_OUT", B.name);
    const std::string GS_DIR = env_s("LJ_GS_DIR", "");

    int CKPT_EVERY = env_i("LJ_CKPT_EVERY", 0);
    if (CKPT_EVERY <= 0) CKPT_EVERY = std::max(1, N_STEPS / 5);  // >= 4 retained

    if (GS_DIR.empty() || !std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: LJ_GS_DIR missing or unset: '" << GS_DIR << "'\n";
        return 2;
    }

    const double CAP_WIDTH_FRAC = CAP_L / B.LZ_BOHR;
    const double CAP_MID_FRAC   = 0.5 - CAP_WIDTH_FRAC / 2.0;
    const double z_cap_in       = B.LZ_BOHR/2.0 - CAP_L;
    const bool   CAP_ON         = (CAP_ETA != 0.0);
    const double sigma_p2       = 1.0 / (2.0 * SIGMA_WP * SIGMA_WP);
    const double rs             = inqkit::jellium::rs_from_n0(B.n0());
    // transverse periodic-image overlap: 6 sigma_d = L_xy
    const double t_ov_arg = 2.0 * std::pow(Sh::LX_BOHR/6.0, 2) - SIGMA_WP*SIGMA_WP;
    const double t_ov     = (t_ov_arg > 0.0) ? SIGMA_WP * std::sqrt(t_ov_arg) : 0.0;

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
              << "\n=== lz_bulk_sweep WP [" << B.name << "]  OUT=" << OUT << " ===\n"
              << "  cell      = " << Sh::LX_BOHR << " x " << Sh::LY_BOHR << " x "
              << B.LZ_BOHR << " Bohr, periodicity(2), dx=" << SPACING << "\n"
              << "  slab      = " << B.l_slab() << " Bohr (half " << B.SLAB_HALF
              << ", edge " << Sh::EDGE_WIDTH_BOHR << "), N=" << B.N_ELECTRONS
              << ", n0=" << B.n0() << ", r_s=" << rs << "\n"
              << "  WP        = sigma_WP " << SIGMA_WP << " (density std "
              << SIGMA_WP/std::sqrt(2.0) << " == classical sigma_pot)"
              << "  k0=" << K0 << " (v=" << K0 << ")  E_drift="
              << 0.5*K0*K0*HA << " eV\n"
              << "  launch_z  = " << LAUNCH_Z
              << "  (standoff " << (-B.SLAB_HALF - LAUNCH_Z) << " Bohr)\n"
              << "  CAP       = " << (CAP_ON ? "ON" : "OFF") << "  eta=" << CAP_ETA
              << " Ha  width=" << CAP_L << " Bohr/face  bands +/-["
              << z_cap_in << "," << B.LZ_BOHR/2.0 << "]\n"
              << "  dt=" << DT << "  START=" << START << " -> N_STEPS=" << N_STEPS
              << "  t_total=" << (DT*N_STEPS) << " a.u." << (RESUME ? "  [RESUME]" : "") << "\n"
              << "  cadence   : density/" << SAVE_EVERY << "  wavefn/" << WF_EVERY
              << "  stats/1  ckpt/" << CKPT_EVERY << "\n"
              << "  GS        = " << GS_DIR << "\n"
              << "  spreading : rate " << 1.0/(std::sqrt(2.0)*SIGMA_WP)
              << " Bohr/a.u.;  T1-T2 = " << 3.0/(4.0*SIGMA_WP*SIGMA_WP)*HA << " eV"
              << ";  t_ov = " << t_ov << " a.u.\n\n";

    // ---- system -------------------------------------------------------------
    auto cell = systems::cell::orthorhombic(Sh::LX_BOHR * 1.0_b,
                                            Sh::LY_BOHR * 1.0_b,
                                            B.LZ_BOHR * 1.0_b).periodicity(2);
    auto ions = systems::ions(cell);                    // jellium: no nuclei
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(B.N_ELECTRONS)
            .extra_states(B.EXTRA_STATES)
            .temperature(Sh::TEMPERATURE_EV * 1.0_eV),
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
                  << "  max_overlap=" << report.max_overlap << "\n"
                  << "  ortho loss : removed_weight=" << report.removed_weight
                  << " (" << 100.0*report.removed_weight << " %)"
                  << "  closure_res=" << report.ortho_closure_residual() << "\n";

        std::ofstream f(OBS + "/wp_config.txt");
        f << std::setprecision(16)
          << "wp_center_bohr = 0 0 " << LAUNCH_Z << "\n"
          << "wp_sigma_bohr  = " << SIGMA_WP << "\n"
          << "wp_sigma_density = " << SIGMA_WP/std::sqrt(2.0) << "\n"
          << "wp_k0_bohr_inv = " << K0 << "\n"
          << "wp_state_index = " << wp_idx << "\n"
          << "norm_after     = " << report.norm_after << "\n"
          << "max_overlap    = " << report.max_overlap << "\n"
          // Orthogonalisation LOSS: norm_after is measured AFTER renormalisation
          // and cannot report how much the Gram-Schmidt projection removed.
          << "norm_pre_ortho = " << report.norm_pre_ortho << "\n"
          << "norm_pre_renorm = " << report.norm_pre_renorm << "\n"
          << "removed_weight = " << report.removed_weight << "\n"
          << "removed_percent = " << 100.0*report.removed_weight << "\n"
          << "sum_overlap_sq = " << report.sum_overlap_sq << "\n"
          << "ortho_closure_residual = " << report.ortho_closure_residual() << "\n";
    }

    // ---- background well + the two CAPs ------------------------------------
    inqkit::jellium::localised_background_params bg;
    bg.shape      = inqkit::jellium::background_shape::slab;
    bg.n0         = B.n0();
    bg.half_width = B.SLAB_HALF;
    bg.slab_axis  = Sh::SLAB_AXIS;
    bg.center     = {0.0, 0.0, Sh::SLAB_CENTER};
    bg.edge_width = Sh::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    perturbations::absorbing cap_lo(CAP_ETA * 1.0_Ha, -CAP_MID_FRAC, CAP_WIDTH_FRAC);
    perturbations::absorbing cap_hi(CAP_ETA * 1.0_Ha,  CAP_MID_FRAC, CAP_WIDTH_FRAC);
    auto pert_with_cap = perturbations::sum(bg_pert, perturbations::sum(cap_lo, cap_hi));

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

    // FULL energy decomposition — every component, every step.
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

    // Delta FIELD only at the density cadence (the every-frame version once
    // filled /rds — see the sigma56 header); the L2 scalar stays every-step.
    inqkit::observables::DensityDelta density_delta(
        VTI + "/density_delta", VTI + "/density_delta_coarse",
        {.emit_raw_vti = true, .emit_coarse_vti = true,
         .compute_l2 = true, .coarse_bin_bohr = 3.0,
         .emit_every = (SAVE_EVERY > 0 ? SAVE_EVERY : 1)});

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
        // RELATIVE tolerances — see sigma56_sv wp/run.cpp for the calibration.
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
        auto pc_dev = [](double got, double want){
            return (want != 0.0) ? 100.0*(got-want)/std::abs(want) : 0.0; };
        std::cout << "\n  --- t=0 analytic gates ---\n";
        gate_abs("norm (real space)",         r0.N,  1.0, 0.02);
        gate_rel("<p_z> = k0",                m0.pz, K0,  2.0);
        gate_rel("sigma_pz^2 = 1/(2 s^2)",    m0.sz2, sigma_p2, 10.0);
        gate_rel("T1 = (k0^2+3 sp2)/2 (Ha)",  T1, 0.5*(K0*K0 + 3.0*sigma_p2), 3.0);
        gate_rel("T1-T2 = 3/(4 s^2) (Ha)",    T1-T2, 3.0/(4.0*SIGMA_WP*SIGMA_WP), 5.0);
        gate_abs("centroid z (circular)",     r0.zc, LAUNCH_Z, 0.05);

        // Orthogonalisation loss — ALWAYS strict (user criterion, 2026-08-01).
        gate_abs("ortho removed weight < 3 %", 100.0*report.removed_weight, 0.0,
                 env_d("LJ_ORTHO_TOL_PC", 3.0));

        // Real-space width: raw second moment for a tail-free far launch,
        // consistency check against the ortho tail otherwise (see sigma56).
        {
            const double var_meas  = r0.sz2;
            const double var_gauss = SIGMA_WP*SIGMA_WP/2.0;
            const bool   tail_free = (report.removed_weight < 1.0e-5);
            if (tail_free) {
                gate_rel("density std = s/sqrt2", std::sqrt(var_meas),
                         SIGMA_WP/std::sqrt(2.0), 5.0);
            } else {
                const double dvar = var_meas - var_gauss;
                const double d_implied = (dvar > 0.0 && report.removed_weight > 0.0)
                    ? std::sqrt(dvar / report.removed_weight) : 0.0;
                std::cout << "  [info] density std = " << std::sqrt(var_meas)
                          << " vs Gaussian " << SIGMA_WP/std::sqrt(2.0)
                          << " (" << pc_dev(std::sqrt(var_meas), SIGMA_WP/std::sqrt(2.0))
                          << " %) — inflated by the PHYSICAL orthogonality tail,"
                             " not a deformed core\n"
                          << "  [info] implied tail RMS distance = " << d_implied
                          << " Bohr (launch " << LAUNCH_Z << ", slab face "
                          << -B.SLAB_HALF << ")\n";
                const bool ok = (d_implied > 3.0*SIGMA_WP/std::sqrt(2.0))
                             && (d_implied < B.LZ_BOHR/2.0);
                std::cout << (ok ? "  [PASS] " : "  [FAIL] ")
                          << "excess variance consistent with the ortho tail"
                          << "  (d_implied " << d_implied << " Bohr, want "
                          << 3.0*SIGMA_WP/std::sqrt(2.0) << " < d < "
                          << B.LZ_BOHR/2.0 << ")\n";
                if (!ok) ++fails;
            }
        }
        std::cout << "  [info] max_overlap with occupied manifold = "
                  << report.max_overlap << " (want < 1e-3)\n";

        // Momentum-space aliasing diagnostic. Matters again at sigma = 0.5
        // (sigma_p = 1.414): the recorded dx = 0.40 table gives sigma_pz^2 bias
        // +0.05/+0.26/+1.24/+5.06 % at v = 2.0-3.5 — reported, not gated.
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
    inqkit::RealTimeSession rt(ions, electrons, 1);     // callback EVERY step
    rt.add([&](inqkit::StepContext const& ctx) {
        const int step = ctx.step;

        auto n_tot = inqkit::fields::density::total(*ctx.electrons);
        const double l2 = density_delta.snapshot(n_tot, ctx.time_au, step);
        inqkit::StepContext out = ctx; out.density_l2 = l2;
        obs.append(out);

        // pairwise Coulomb ledger with the WP as charge group P
        auto n_wp_f = inqkit::jellium::orbital_density_field(*ctx.electrons, wp_idx);
        auto ct = inqkit::jellium::compute_coulomb_wp(ctx.electrons->density(), n_wp_f, phiplus);
        if (ctx.electrons->root())
            ix << step << ',' << ctx.time_au << ',' << ct.e_ss << ',' << ct.e_pp << ','
               << ct.e_ps << ',' << ct.e_sb << ',' << ct.e_pb << ',' << E_BB << ','
               << ct.e_hartree_check << ',' << ct.e_external_check << ','
               << ct.norm_wp << ',' << ct.norm_total << '\n';

        if (SAVE_EVERY > 0 && step % SAVE_EVERY == 0) {
            total_wr.write(n_tot, ctx.time_au, step);
            wp_wr.write(inqkit::fields::density::orbital(*ctx.electrons, wp_idx), ctx.time_au, step);
        }
        if (WF_EVERY > 0 && step % WF_EVERY == 0) {
            auto wf = inqkit::fields::orbital::wavefunction(*ctx.electrons, wp_idx);
            wf_wr.write(wf, "wavefunction_t" + tag6(step));
        }
    });

    auto write_rt_state = [&](int last){
        if (!electrons.root()) return;
        std::ofstream st(RT_STATE);
        st << std::setprecision(12)
           << "last_step=" << last << "\ntime_au=" << (last*DT) << "\ndt=" << DT
           << "\nwp_idx=" << wp_idx << "\nk0=" << K0 << "\nsigma_wp=" << SIGMA_WP << "\n";
    };

    real_time::propagate(
        ions, electrons,
        [&](auto const& data) {
            rt.step(data);
            wp_mom.maybe_accumulate(data);
            wp_pos.maybe_accumulate(data);
            // Interior checkpoints: rolling `checkpoint` plus a RETAINED
            // numbered snapshot, so >= 4 distinct resume points survive.
            if (data.iter() > 0 && data.iter() % CKPT_EVERY == 0) {
                electrons.save(CKPT);
                write_rt_state(data.iter());
                const std::string snap = OUT + "/ckpt_step" + tag6(data.iter());
                electrons.save(snap);
                std::cout << "  [ckpt] step " << data.iter() << " -> " << snap << "\n" << std::flush;
            }
        },
        options::theory{}.lda(),
        options::real_time{}.num_steps(N_STEPS).dt(DT * 1.0_atomictime),
        pert_with_cap, START);

    // ---- final checkpoint ---------------------------------------------------
    electrons.save(CKPT);
    write_rt_state(N_STEPS);

    const double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_wall0).count();

    if (electrons.root()) {
        ix.close();
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(14)
          << "run = localised_jellium/lz_bulk_sweep/wp/" << env_s("LJ_OUT", B.name) << "\n"
          << "run_type = wavepacket projectile, localised jellium slab TDDFT (ALDA)\n"
          << "campaign = lz_bulk_sweep (slab->bulk L_slab extrapolation)\n"
          << "plan = docs/plans/jellium-slab-extend-Lz.md\n"
          << "engine = inq-study\nxc = LDA (ALDA in TDDFT)\n"
          << "date_finished = " << iso_now() << "\nwall_time_s = " << wall << "\n"
          << "box_preset = " << B.name << "\n"
          << "cell_bohr = " << Sh::LX_BOHR << "x" << Sh::LY_BOHR << "x" << B.LZ_BOHR << "\n"
          << "periodicity = 2\nspacing_bohr = " << SPACING << "\n"
          << "slab_half_width = " << B.SLAB_HALF
          << "  slab_thickness = " << B.l_slab()
          << "  edge_width = " << Sh::EDGE_WIDTH_BOHR << "\n"
          << "n_electrons = " << B.N_ELECTRONS << "  n0_a0m3 = " << B.n0()
          << "  r_s = " << rs << "\n"
          << "n_states = " << n_states << "  extra_states = " << B.EXTRA_STATES << "\n"
          << "wp_enabled = yes\nwp_state_index = " << wp_idx << "\n"
          << "wp_sigma_bohr = " << SIGMA_WP << "\n"
          << "wp_sigma_note = wavepacket sigma (psi width); density std = this/sqrt2 "
             "= classical sigma_pot\n"
          << "wp_sigma_density = " << SIGMA_WP/std::sqrt(2.0) << "\n"
          << "wp_k0_bohr_inv = " << K0 << "  wp_velocity = " << K0 << "\n"
          << "wp_drift_energy_ev = " << 0.5*K0*K0*HA << "\n"
          << "wp_localisation_energy_ev = " << 3.0/(4.0*SIGMA_WP*SIGMA_WP)*HA << "\n"
          << "launch_z = " << LAUNCH_Z << "  standoff_bohr = " << (-B.SLAB_HALF - LAUNCH_Z) << "\n"
          << "norm_after = " << report.norm_after << "  max_overlap = " << report.max_overlap << "\n"
          << "cap = " << (CAP_ON ? "on" : "off") << "  cap_eta_ha = " << CAP_ETA
          << "  cap_width_bohr = " << CAP_L << " per face\n"
          << "cap_mid_frac = " << CAP_MID_FRAC << "  cap_width_frac = " << CAP_WIDTH_FRAC << "\n"
          << "cap_band_hi_bohr = [" << z_cap_in << "," << B.LZ_BOHR/2.0 << "]\n"
          << "cap_band_lo_bohr = [" << -B.LZ_BOHR/2.0 << "," << -z_cap_in << "]\n"
          << "cap_note = energy_total NOT conserved (non-Hermitian); gate on norm instead\n"
          << "start_step = " << START << "  rt_num_steps = " << N_STEPS
          << "  dt_au = " << DT << "  total_time_au = " << (DT*N_STEPS) << "\n"
          << "save_every = " << SAVE_EVERY << "  wf_every = " << WF_EVERY
          << "  stats_every = 1  ckpt_every = " << CKPT_EVERY << "\n"
          << "spread_rate_bohr_per_au = " << 1.0/(std::sqrt(2.0)*SIGMA_WP) << "\n"
          << "t_transverse_overlap_au = " << t_ov << "  (6 sigma_d = L_xy = "
          << Sh::LX_BOHR << "; step " << (t_ov > 0.0 ? int(t_ov/DT) : 0) << ")\n"
          << "analysis_note = deposit S = [E_total(t_f) - E_GS - E_PS(t_f)] / "
          << B.l_slab() << " Bohr, norm-corrected; per-box E_GS from the gs run\n"
          << "gs_dir = " << GS_DIR << "\nrun_completed = true\n";
    }
    std::cout << "\nDone. Wall time " << wall << " s.\n";
    return 0;
}
