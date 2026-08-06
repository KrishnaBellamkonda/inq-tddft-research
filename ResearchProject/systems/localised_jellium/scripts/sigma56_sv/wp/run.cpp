// ============================================================================
// localised_jellium / scripts/sigma56_sv / wp / run.cpp
//
// WAVEPACKET half of the sigma_WP = 5 and 6 Bohr twin campaign.
// Plan: docs/plans/sigma56-sv-twin.md
//
// Cloned VERBATIM from scripts/wp_highdensity_sv/wp/run.cpp (the validated
// sigma = 0.5/2/3 production binary) with exactly three changes:
//   (1) Cfg -> SlabN100_L35x35x105  (L_z 85 -> 105; slab, N and r_s unchanged)
//   (2) launch z default -24.0 -> Cfg::LAUNCH_Z_BOHR = -27.5
//   (3) N_STEPS / sigma defaults retuned for the new box (see the table below)
// Every t=0 gate, the orthogonalisation budget, the pairwise ledger, the
// checkpoint/resume block and the CAP wiring are unchanged, so nothing that was
// validated for the 85-Bohr campaigns needs re-validating here. The CAP
// fractions are DERIVED from Cfg::LZ_BOHR, so they follow the box automatically.
//
// ---------------------------------------------------------------------------
// WHY sigma = 5 AND 6, AND WHY THE BOX HAD TO GROW
// sigma_d(t) = sqrt(sigma^2/2 + t^2/(2 sigma^2)): sigma_WP fixes both the initial
// width AND the spreading rate 1/(sqrt2 sigma). At sigma = 5/6 the packet grows
// by only x1.23/x1.12 over the in-slab transit (against x3.2 at sigma = 0.5), so
// it is effectively CONSTANT-WIDTH and its label agrees with its time-average —
// which is what makes a classical twin at a fixed sigma_pot a fair comparison.
// But sigma_d(0) = sigma_WP/sqrt2 = 4.243 Bohr at sigma = 6, and the 85-Bohr box
// has only 17.5 Bohr between the slab face (-12.5) and the CAP inner edge (-30);
// 3 sigma_d of clearance on each side needs 25.5. Hence L_z = 105: the gap
// becomes 27.5 Bohr and launch z = -27.5 leaves 12.5 Bohr to the CAP (2.95
// sigma_d, 0.16 % of the packet inside at t=0) and 15.0 to the slab (3.54
// sigma_d, 0.020 % inside) — at or below the 0.23 % already accepted at sigma=3.
//
// ---------------------------------------------------------------------------
// THE WIDTH MAPPING (.claude/rules/sigma-wp-convention.md)
// sigma always means sigma_WP, the psi-width. The classical twin runs a Gaussian
// erf/r potential of std sigma_pot = sigma_WP/sqrt(2) (3.5355 at sigma = 5,
// 4.2426 at sigma = 6), which equals this packet's DENSITY std at t = 0 — so the
// two projectiles present the identical charge cloud at launch. Both halves are
// labelled sigma_WP. LJ_SIGMA is sigma_WP in BOTH binaries — the sqrt(2) lives
// inside each one, never at the call site.
//
// ---------------------------------------------------------------------------
// WHY THERE IS A CAP HERE WHEN THE CLASSICAL CAMPAIGN IS CAP-FREE
// The classical campaign is CAP-free by design: periodicity(2) makes z open for
// the ELECTROSTATICS, so the moving Gaussian CHARGE leaves the box and energy is
// conserved, giving an exact post-exit plateau. That mechanism does NOT transfer
// to a wavepacket. periodicity(2) is consulted only by the Poisson solver
// (inq/src/solvers/poisson.hpp:189,206), ionic replicas and the kick gauge; the
// wavefunction basis and kinetic operator are a plain 3-D FFT periodic in ALL
// three directions (inq/src/basis/fourier_space.hpp:60-151,
// inq/src/hamiltonian/ks_hamiltonian.hpp:200-204). A KS ORBITAL travelling +z
// therefore WRAPS and re-enters at -z, arriving behind the slab and
// contaminating the bath. (Confirmed in-repo: docs/handovers/pbc-open-z-
// oscillation.md:20 "wavefunction always wraps on the FFT grid (p2 switches
// electrostatics only)".)
//
// So two absorbing bands are added at the z faces (user decision 2026-07-30):
//   width 12.5 Bohr per face, |eta| = 1 Ha
// perturbations::absorbing takes FRACTIONAL cell coordinates (it compares
// point_op.rvector()[2], which uses the contravariant spacing and lies in
// [-0.5,0.5) — see real_space.hpp:105,129; the constructor's
// assert(mid_pos in [-0.5,0.5)) is the tell):
//   CAP_WIDTH_FRAC = 12.5/105 = 0.119047619048
//   CAP_MID_FRAC   = 0.5 - W/2 = 0.440476190476   ( = 46.25 Bohr )
//   +z band z in [+40.0,+52.5];  -z band z in [-52.5,-40.0]
// eta < 0 ABSORBS (exp(-iVt) = exp(eta sin^2 t)); +1.0 would be a gain medium.
// The CAP geometry and strength are validated independently, before production,
// by scripts/sigma56_sv/vac/ (free WP, no slab) — which also supplies the
// per-(sigma, v) CAP-only baselines.
//
// CONSEQUENCES OF THE CAP, which the notebooks must state:
//   * energy_total is NO LONGER CONSERVED (the CAP is non-Hermitian), so the
//     classical correctness gate is replaced by norm/absorption monitoring.
//   * the WP norm decays; WPMomentumStats divides every moment by the CURRENT
//     norm, so T1/T2 stay valid expectation values OF THE SURVIVING PACKET —
//     but the CAP preferentially removes the slow, spread tail, which biases
//     <p_z> upward. This is the qsp5 lesson (docs/handovers/qsp5-momentum-
//     stopping.md) and is why norm_check is recorded every step.
//
// ---------------------------------------------------------------------------
// WHY E_absorbed HAS NO DIRECT WP ANALOGUE
// In the classical run the projectile is an EXTERNAL moving perturbation doing
// work, so energy_total climbs and plateaus and E_absorbed = plateau - E_GS.
// Here the projectile IS part of the system, so (CAP aside) the Hamiltonian is
// time-independent and energy_total is conserved: there is no plateau to read.
// The deposited energy must come from the pairwise Coulomb ledger instead —
// hence interactions.csv is written EVERY step via compute_coulomb_wp, whose
// documented WP closure is
//     E_hartree = E_SS + E_PS + E_PP ,  E_external = E_SB + E_PB
// and which emits e_hartree_check / e_external_check as closure gates against
// INQ's own scalars.
//
// ---------------------------------------------------------------------------
// STOPPING POWER (definitions locked in docs/plans/bulk-jellium-ks-stopping.md §4)
//   T1 = <p^2>/2m   -> wp_momentum_stats.csv, e_kin_ha
//   T2 = <p>^2/2m   -> 0.5*(px_mean^2+py_mean^2+pz_mean^2), same file
//   s3 = WP density centroid -> wp_real_space_stats.csv, z_mean_circ (CIRCULAR:
//        the naive centroid is discontinuous across the periodic z face)
//   s4 = integral of <p_z> dt -> cumulative trapezoid of pz_mean
// S_ij = -dT_i/ds_j. Both stats files are written EVERY step: they ARE the
// measurement, and one extra FFT per step is negligible against ~74 orbitals.
//
// ANALYSIS WINDOW — WIDE OPEN AT THESE WIDTHS, unlike sigma = 0.5.
// The two clocks that bound a usable fit are (a) the packet centroid being inside
// the slab and (b) the packet not yet overlapping its own transverse periodic
// images (6 sigma_d = L_xy = 35, i.e. t_ov = sigma*sqrt(2(L_xy/6)^2 - sigma^2)).
// At sigma = 0.5 those windows did not intersect AT ALL at v = 2.0/2.5
// (t_ov = 4.12 a.u. against a slab entry at 5.75). Here:
//     t_ov = 32.8 a.u. (sigma=5) / 34.0 a.u. (sigma=6)
//     in-slab transit = [7.50, 20.00] a.u. at v = 2.0, shrinking with v
// so the ENTIRE transit is transversely clean at every velocity, and the leading
// tail does not reach the +z CAP (inner edge +40) until t ~ 26-34 a.u. Both the
// deposit estimator and the localised slope estimators are usable.
// T1 - T2 = 3/(4 sigma^2) = 0.030/0.021 Ha = 0.82/0.57 eV of localisation energy,
// against 81.6 eV at sigma = 0.5 — the packet is now overwhelmingly drift energy.
//
// ---------------------------------------------------------------------------
// CHECKPOINTS (.claude/rules/final-timestep-checkpoint.md + checkpoint-dont-block)
// The user asked for AT LEAST 4 resumable checkpoints per run. This writes:
//   * results/checkpoint          — rolling latest; what LJ_RESUME=1 loads
//   * results/ckpt_step<N>        — RETAINED numbered snapshots, one per
//                                   LJ_CKPT_EVERY (default N_STEPS/5 -> 5 of them)
//   * results/rt_state.txt        — last_step, time_au, dt, wp_idx
// so a killed job loses at most one interval and any quarter can be resumed from.
//
// ENGINE: inq-study, NOT stock inq (the CAP's complex scalar potential does not
// compile, let alone propagate, against stock inq).
//
// STEP COUNTS (dispatcher supplies them; defaults below are the v = 2.0 point).
// Same convention as the 85-Bohr campaigns, re-evaluated for the longer path:
//     N_STEPS = round( 4.36 * (|launch_z| + L_z/2) / (v * dt) ),  dt = 0.04
// (that formula reproduces the recorded 3623 at v=2.0/z0=-24/L_z=85 to 1 step).
//   idx  v     N_steps   t (a.u.)   save/  wf/  ckpt/
//     0  2.0     4360      174.4      14   43    872
//     1  2.5     3488      139.5      12   35    698
//     2  3.0     2907      116.3      10   29    581
//     3  3.5     2491       99.6       8   25    498
// The packet clears the slab by t = 20 a.u. and reaches the CAP by t ~ 34 a.u.
// at the slowest velocity, so the remainder is plateau time for the deposit.
//
// Env: LJ_K0(2.0) LJ_N_STEPS(4360) LJ_DT(0.04) LJ_SIGMA(6.0) LJ_LAUNCH_Z(-27.5)
//      LJ_CAP_ETA(-1.0) LJ_CAP_L(12.5) LJ_SAVE_EVERY(14) LJ_WF_EVERY(43)
//      LJ_CKPT_EVERY(0=auto N/5) LJ_RESUME(0) LJ_OUT(REQUIRED) LJ_GS_DIR(REQUIRED)
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

#include "../../../shared/configs/slab_n100_L35x35x105.hpp"

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
using Cfg = localised_jellium::config::SlabN100_L35x35x105;

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

    // ---- parameters (geometry pinned to Cfg = SlabN100_L35x35x105) ---------
    const double SIGMA_WP = env_d("LJ_SIGMA", Cfg::WP_SIGMA_BOHR);   // 6.0; sweep {5,6}
    // dx = 0.40, inherited from the 85-Bohr production runs so the two campaigns
    // share a grid. It is generous here: the aliasing that forced 0.50 -> 0.40 at
    // sigma_WP = 0.5 scales as sigma_p = 1/(sqrt2 sigma) = 0.141/0.118 at sigma =
    // 5/6, against k_Nyq = pi/0.40 = 7.85, i.e. ~1e-14 % moment bias at every
    // velocity. dx is held anyway: it must match the GS_DIR's spacing, and a
    // coarser grid would break comparability with the existing traces.
    const double SPACING  = env_d("LJ_SPACING", Cfg::SPACING_BOHR);
    const double K0       = env_d("LJ_K0", 2.0);                     // = v (m = 1)
    // Common to BOTH halves and BOTH sigma (plan §4): the classical twin must
    // start at the identical z or the comparison is not a twin.
    const double LAUNCH_Z = env_d("LJ_LAUNCH_Z", Cfg::LAUNCH_Z_BOHR);  // -27.5
    const double DT       = env_d("LJ_DT", 0.04);
    const int    N_STEPS  = env_i("LJ_N_STEPS", 4360);               // v = 2.0 point
    const double CAP_ETA  = env_d("LJ_CAP_ETA", Cfg::CAP_ETA_HA);    // Ha, <0 absorbs
    const double CAP_L    = env_d("LJ_CAP_L", Cfg::CAP_L_BOHR);      // Bohr per face
    const int    SAVE_EVERY = env_i("LJ_SAVE_EVERY", 14);            // ~310 frames
    const int    WF_EVERY   = env_i("LJ_WF_EVERY", 43);              // ~100 wavefunctions
    const bool   RESUME     = env_i("LJ_RESUME", 0) != 0;
    const std::string OUT    = "results/" + env_s("LJ_OUT", "wp");
    const std::string GS_DIR = env_s("LJ_GS_DIR", "");

    int CKPT_EVERY = env_i("LJ_CKPT_EVERY", 0);
    if (CKPT_EVERY <= 0) CKPT_EVERY = std::max(1, N_STEPS / 5);      // >= 4 retained

    if (GS_DIR.empty() || !std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: LJ_GS_DIR missing or unset: '" << GS_DIR << "'\n";
        return 2;
    }

    const double CAP_WIDTH_FRAC = CAP_L / Cfg::LZ_BOHR;
    const double CAP_MID_FRAC   = 0.5 - CAP_WIDTH_FRAC / 2.0;
    const double z_cap_in       = Cfg::LZ_BOHR/2.0 - CAP_L;          // +40.0
    const bool   CAP_ON         = (CAP_ETA != 0.0);
    const double sigma_p2       = 1.0 / (2.0 * SIGMA_WP * SIGMA_WP); // 0.02 / 0.0139
    const double rs             = inqkit::jellium::rs_from_n0(Cfg::N0);

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
              << "\n=== sigma56_sv WP (L_z=105)  OUT=" << OUT << " ===\n"
              << "  cell      = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR << " x "
              << Cfg::LZ_BOHR << " Bohr, periodicity(2), dx=" << SPACING << "\n"
              << "  slab      = 25 Bohr (half " << Cfg::SLAB_HALF_WIDTH << ", edge "
              << Cfg::EDGE_WIDTH_BOHR << "), N=" << Cfg::N_ELECTRONS
              << ", n0=" << Cfg::N0 << ", r_s=" << rs << "\n"
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
            .extra_electrons(Cfg::N_ELECTRONS)
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
          // Orthogonalisation LOSS. norm_after is measured AFTER renormalisation
          // and is ~1 by construction, so it cannot report how much of the packet
          // the Gram-Schmidt projection removed. That matters as soon as the
          // launch point moves into the slab's electronic spill-out (the
          // near-launch campaign, docs/plans/effective-sigma-near-launch.md);
          // it is ~1e-6 at the far-launch z = -24.
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
    bg.n0         = Cfg::N0;
    bg.half_width = Cfg::SLAB_HALF_WIDTH;
    bg.slab_axis  = Cfg::SLAB_AXIS;
    bg.center     = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR};
    bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
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
        auto pc_dev = [](double got, double want){
            return (want != 0.0) ? 100.0*(got-want)/std::abs(want) : 0.0; };
        std::cout << "\n  --- t=0 analytic gates ---\n";
        gate_abs("norm (real space)",         r0.N,  1.0, 0.02);
        gate_rel("<p_z> = k0",                m0.pz, K0,  2.0);
        gate_rel("sigma_pz^2 = 1/(2 s^2)",    m0.sz2, sigma_p2, 10.0);
        gate_rel("T1 = (k0^2+3 sp2)/2 (Ha)",  T1, 0.5*(K0*K0 + 3.0*sigma_p2), 3.0);
        gate_rel("T1-T2 = 3/(4 s^2) (Ha)",    T1-T2, 3.0/(4.0*SIGMA_WP*SIGMA_WP), 5.0);
        gate_abs("centroid z (circular)",     r0.zc, LAUNCH_Z, 0.05);

        // ---- orthogonalisation loss (near-launch campaign) ------------------
        // ALWAYS strict: the packet that propagates must still be mostly the
        // Gaussian we asked for. 3 % is the user's criterion (2026-08-01).
        gate_abs("ortho removed weight < 3 %", 100.0*report.removed_weight, 0.0,
                 env_d("LJ_ORTHO_TOL_PC", 3.0));

        // ---- real-space width ------------------------------------------------
        // The raw second moment is the RIGHT gate for a far launch and the WRONG
        // one for a near-slab launch, so it is applied conditionally.
        //
        // WHY. Launched inside the slab's electronic spill-out, the packet is
        // Pauli-orthogonalised against occupied states that extend through the
        // whole slab, so it necessarily acquires a small, PHYSICAL far tail.
        // Variance weights by (z-zc)^2, so ~0.1 % of the weight sitting ~10 Bohr
        // away inflates sqrt(sz2) by tens of percent while the packet CORE is
        // untouched. Measured at z = -14 (job 32528019): full-profile std 0.4684
        // (+32.5 %) but CORE std 0.3559 vs the Gaussian 0.3536 — +0.65 %.
        // Gating the raw second moment there rejects a CORRECT packet (it did:
        // job 32528175 aborted on exactly this).
        //
        // Nothing is lost by relaxing it, because sigma is probed by two
        // TAIL-IMMUNE gates that stay strict in both regimes: sigma_pz^2 and
        // T1-T2 = 3/(4 sigma^2). Momentum moments are not distance-weighted, so
        // a genuine factor-2 sigma blunder is still caught. What replaces it is
        // a CONSISTENCY check: the excess variance must be explainable by the
        // measured loss placing weight at a physically sensible distance,
        //     dvar = var_measured - sigma^2/2  ~=  w_tail * d^2,
        // so d = sqrt(dvar / removed_weight) must land beyond the packet core
        // and inside the box. At z = -14 this gives d = 9.3 Bohr, which is
        // exactly the distance from the launch point to the slab orbitals.
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
                          << -Cfg::SLAB_HALF_WIDTH << ")\n";
                // A sensible tail sits beyond the core (> 3 sigma_d) and inside
                // the box. Outside that range the excess is NOT an orthogonality
                // tail and something really is wrong.
                const bool ok = (d_implied > 3.0*SIGMA_WP/std::sqrt(2.0))
                             && (d_implied < Cfg::LZ_BOHR/2.0);
                std::cout << (ok ? "  [PASS] " : "  [FAIL] ")
                          << "excess variance consistent with the ortho tail"
                          << "  (d_implied " << d_implied << " Bohr, want "
                          << 3.0*SIGMA_WP/std::sqrt(2.0) << " < d < "
                          << Cfg::LZ_BOHR/2.0 << ")\n";
                if (!ok) ++fails;
            }
        }
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
            // Interior checkpoints: rolling `checkpoint` (what LJ_RESUME loads) plus
            // a RETAINED numbered snapshot, so >= 4 distinct resume points survive.
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
          << "run = localised_jellium/sigma56_sv/wp/" << env_s("LJ_OUT","wp") << "\n"
          << "run_type = wavepacket projectile, localised jellium slab TDDFT (ALDA)\n"
          << "campaign = classical-highdensity-sv (WP twin)\n"
          << "plan = docs/plans/sigma56-sv-twin.md\n"
          << "engine = inq-study\nxc = LDA (ALDA in TDDFT)\n"
          << "date_finished = " << iso_now() << "\nwall_time_s = " << wall << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR << "\n"
          << "periodicity = 2\nspacing_bohr = " << SPACING << "\n"
          << "slab_half_width = " << Cfg::SLAB_HALF_WIDTH
          << "  edge_width = " << Cfg::EDGE_WIDTH_BOHR << "\n"
          << "n_electrons = " << Cfg::N_ELECTRONS << "  n0_a0m3 = " << Cfg::N0
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
