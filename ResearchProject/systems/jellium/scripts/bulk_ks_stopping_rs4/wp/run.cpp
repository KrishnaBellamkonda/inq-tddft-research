// ============================================================================
// scripts/bulk_ks_stopping_rs4/wp/run.cpp
//
// WAVEPACKET half of the bulk-jellium KS-orbital stopping-power twin pair.
// Plan: docs/plans/bulk-jellium-ks-stopping.md (user design lock 2026-07-30).
//
// A Gaussian electron wavepacket (sigma_WP = 2 Bohr, E = 100 eV, k0 = 2.7111)
// is injected as an extra KS orbital into a fully periodic 40 x 40 x 80 Bohr
// jellium bath (N = 482, r_s = 3.987, dx = 0.50) and propagated 646 steps at dt = 0.04
// (t = 25.84 a.u.) from z = -32 to the +z face.
//
// WHAT THIS RUN IS FOR. The stopping power S = -dT/ds is extracted from four
// combinations of two KE definitions and two position definitions:
//
//   T1 = <p^2>/2m   -> wp_momentum_stats.csv, column e_kin_ha
//   T2 = <p>^2/2m   -> 0.5*(px_mean^2 + py_mean^2 + pz_mean^2), same file
//   s3 = WP density centroid          -> wp_real_space_stats.csv, z_mean_circ
//   s4 = integral of <p_z> dt         -> cumulative trapezoid of pz_mean
//
// Both stats files are written EVERY step (STATS_EVERY = 1): they are the
// measurement, and one extra FFT per step is negligible against the ~129
// orbital FFTs the propagator already does.
//
// s3 USES THE CIRCULAR COLUMN, NOT z_mean. The naive integral of z|psi|^2 is
// discontinuous across a periodic face; the circular (phase) estimator added to
// WPRealSpaceStats on 2026-07-30 is exact in a periodic cell. Both are written
// so the run can cross-check them: they must agree to numerical precision until
// the packet nears the +z face, and their divergence marks where the naive one
// stopped being usable.
//
// EHRENFEST CROSS-CHECK. This run has NO ions, so the KS Hamiltonian is purely
// local (kinetic + Hartree + ALDA) and there is no CAP. Ehrenfest's theorem then
// gives d<z>/dt = <p_z>/m EXACTLY, so s3 and s4 MUST coincide. Their comparison
// is a validation of the propagation, not a second physics channel. (Contrast the
// qsp5 runs, where CAP non-unitarity broke this identity at t ~ 5 a.u.)
//
// ANALYSIS WINDOW. Fit over t in [4.0, 18.97] a.u. (steps 100..474): the lower
// edge drops the injection/orthogonalisation transient, the upper edge is the
// DISPERSION-AWARE interference-free limit, where the spreading packet's leading
// 3-sigma tail reaches the +z face. Do NOT use the legacy static ifw_end_z rule,
// which claims 24.3 a.u. here (a 22 % over-estimate — it ignores spreading).
// Steps beyond 474 are still recorded: they carry the exit/interference physics
// and the density GIFs, they are simply excluded from the slope fit.
//
// t = 0 ANALYTIC GATES (checked below, run aborts on failure):
//   <p_z>  = 2.7111 Bohr^-1        (= k0)
//   sigma_pz^2 = 1/(2 sigma^2) = 0.125       (Heisenberg-saturating)
//   <p^2>/2 = 0.5*(k0^2 + 3 sigma_p^2) = 3.86242 Ha = 105.10 eV
//   T1 - T2 = 3/(4 sigma^2) = 0.1875 Ha = 5.102 eV   (the localisation energy)
//
// RESUME. Supported (this run has no ions, so INQ's start_step restart path is
// available). BKS_RESUME=1 loads results/checkpoint + rt_state.txt and continues
// to BKS_N_STEPS, writing segment-suffixed CSVs. See
// .claude/rules/final-timestep-checkpoint.md.
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/jellium/interaction_energies.hpp>   // P/S/B pairwise decomposition
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/observables/state_energy_writer.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/injection_report.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include "../../../shared/configs/bulk_ks_stopping_L40x40x80_rs4.hpp"
#include "../../../shared/cpp/eigenvalues_writer.hpp"
#include "../../../shared/cpp/rt_state.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = jellium::config::Bulk_KS_Stopping_L40x40x80_rs4_WP;

static std::string iso_now() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return std::string(buf);
}

static int env_int(char const* k, int dflt) {
    char const* v = std::getenv(k);
    return v ? std::atoi(v) : dflt;
}

int main() {
    auto t_wall0 = std::chrono::steady_clock::now();

    const std::string RUN_NAME = "bulk_ks_stopping_rs4_wp";
    const std::string REPO =
        "/rds/user/skcb2/hpc-work/tddft/inq-tddft-research";
    const std::string GS_DIR = REPO +
        "/ResearchProject/systems/jellium/checkpoints/"
        "gs_L40x40x80_orth_N482_dx0p50";

    const int  N_STEPS = env_int("BKS_N_STEPS", Cfg::N_STEPS);
    const bool RESUME  = env_int("BKS_RESUME", 0) != 0;

    const double HA_TO_EV = 27.211386;
    // sigma_p^2 = 1/(2 sigma^2) for psi ~ exp(-r^2/2 sigma^2): the FT is
    // exp(-sigma^2 (k-k0)^2/2), so |psi~|^2 has variance 1/(2 sigma^2). Together
    // with the real-space density std sigma/sqrt2 this saturates Heisenberg
    // (sigma_d sigma_p = 1/2). NOT 1/(4 sigma^2) — that was a factor-2 error in
    // the wp_momentum_stats.hpp docstring (corrected 2026-07-30) which aborted
    // job 32401321 on a wrong gate.
    const double sigma_p2 = 1.0 / (2.0 * Cfg::WP_SIGMA_BOHR * Cfg::WP_SIGMA_BOHR);

    std::cout << "\n=== " << RUN_NAME << " ===\n"
              << "  cell     = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR
              << " x " << Cfg::LZ_BOHR << " Bohr (orthorhombic, PERIODIC xyz)\n"
              << "  N_e      = " << Cfg::N_ELECTRONS
              << "   spacing = " << Cfg::SPACING_BOHR << " Bohr\n"
              << "  WP       : sigma=" << Cfg::WP_SIGMA_BOHR
              << "  k0=" << Cfg::WP_K0 << "  E=" << Cfg::WP_EKIN_EV << " eV\n"
              << "  launch   = (" << Cfg::WP_CX_BOHR << ", " << Cfg::WP_CY_BOHR
              << ", " << Cfg::WP_CZ_BOHR << ") Bohr\n"
              << "  dt=" << Cfg::DT_AU << "  N_steps=" << N_STEPS
              << "  t_total=" << (Cfg::DT_AU * N_STEPS) << " a.u.\n"
              << "  fit window = [" << Cfg::FIT_T0_AU << ", " << Cfg::FIT_T1_AU
              << "] a.u.  (IFW=" << Cfg::T_IFW_AU
              << ", transverse=" << Cfg::T_TRANSVERSE_AU << ")\n"
              << "  cadence  : density every " << Cfg::WRITE_EVERY
              << ", wavefunction every " << Cfg::WF_WRITE_EVERY
              << ", stats every " << Cfg::STATS_EVERY << "\n"
              << "  GS       = " << GS_DIR << "\n"
              << "  resume   = " << (RESUME ? "yes" : "no") << "\n\n";

    if (!std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: GS checkpoint missing at " << GS_DIR
                  << " — run save_gs/gs_L40x40x80_orth_N482_dx0p50 first.\n";
        return 2;
    }

    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b,
                                            Cfg::LY_BOHR * 1.0_b,
                                            Cfg::LZ_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);          // jellium: no nuclei at all

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());

    const int n_states = electrons.states().num_states();

    // ----- Output skeleton ------------------------------------------------
    const std::string OBS = "results/raw/observables";
    const std::string VTI = "results/raw/vti";
    for (auto const& d : {OBS, VTI + "/density_total", VTI + "/density_wp",
                          VTI + "/density_delta", VTI + "/density_delta_coarse",
                          VTI + "/density_gs_system", VTI + "/wavefunction_wp"})
        std::filesystem::create_directories(d);

    int START = 0;
    int wp_idx = -1;
    inqkit::InjectionReport report{};

    if (RESUME) {
        auto st = jellium::rt_state::read("results/rt_state.txt");
        START  = st.last_step;
        wp_idx = st.wp_idx;
        if (START >= N_STEPS) {
            std::cout << "  Nothing to do: last_step=" << START
                      << " >= N_STEPS=" << N_STEPS << ". Clean exit.\n";
            return 0;
        }
        electrons.load("results/checkpoint");
        std::cout << "  RESUMED from step " << START
                  << " (wp_idx=" << wp_idx << ")\n";
    } else {
        electrons.load(GS_DIR);
        std::cout << "  Loaded GS from " << GS_DIR << "\n";
        jellium::eigenvalues::copy_from_checkpoint(GS_DIR, OBS + "/eigenvalues");

        // t=0 bath density BEFORE the WP goes in (the baseline every induced
        // density in post-processing is measured against).
        inqkit::io::RealField3DLayout lay{
            .field_name = "density", .include_meta = false,
            .emit_raw = false, .emit_vti = true,
            .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
        inqkit::io::RealField3DWriter gs_wr(VTI + "/density_gs_system", lay,
                                            {.overwrite = true});
        gs_wr.write(inqkit::fields::density::total(electrons),
                    "density_gs_system");

        auto wp = inqkit::WavePacket{}
                      .center(Cfg::WP_CX_BOHR, Cfg::WP_CY_BOHR, Cfg::WP_CZ_BOHR)
                      .sigma(Cfg::WP_SIGMA_BOHR)
                      .k0(Cfg::WP_KX, Cfg::WP_KY, Cfg::WP_KZ)
                      .orthogonalise_against_occupied(electrons);
        report = wp.inject_into_last_extra_state(electrons, 1.0);
        wp_idx = report.state_index;
        std::cout << "  WP injected: state_index=" << wp_idx
                  << "  norm_after=" << report.norm_after
                  << "  max_overlap=" << report.max_overlap << "\n";

        std::ofstream f(OBS + "/wp_config.txt");
        f << std::setprecision(16)
          << "wp_center_bohr = " << Cfg::WP_CX_BOHR << " " << Cfg::WP_CY_BOHR
          << " " << Cfg::WP_CZ_BOHR << "\n"
          << "wp_sigma_bohr  = " << Cfg::WP_SIGMA_BOHR << "\n"
          << "wp_k0_bohr_inv = " << Cfg::WP_KZ << "\n"
          << "wp_energy_ev   = " << Cfg::WP_EKIN_EV << "\n"
          << "wp_state_index = " << wp_idx << "\n"
          << "norm_after     = " << report.norm_after << "\n"
          << "max_overlap    = " << report.max_overlap << "\n"
          << "fit_t0_au      = " << Cfg::FIT_T0_AU << "\n"
          << "fit_t1_au      = " << Cfg::FIT_T1_AU << "\n"
          << "t_ifw_au       = " << Cfg::T_IFW_AU << "\n"
          << "t_transverse_au= " << Cfg::T_TRANSVERSE_AU << "\n";
    }

    // ----- Writers (segment-suffixed on resume) ---------------------------
    const std::string seg = (START > 0) ? (".from" + std::to_string(START)) : "";

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false,
        .emit_raw = false, .emit_vti = true,
        .vti_format = inqkit::io::VTIWriteOptions::Format::binary};

    inqkit::io::RealField3DWriter total_wr(VTI + "/density_total", vti_layout,
                                           {.overwrite = (START == 0)});
    inqkit::io::RealField3DWriter wp_density_wr(VTI + "/density_wp", vti_layout,
                                                {.overwrite = (START == 0)});
    inqkit::io::ComplexField3DWriter wp_wf_wr(
        VTI + "/wavefunction_wp",
        {.field_name = "wavefunction", .include_meta = false,
         .emit_raw = false, .emit_vti = true,
         .vti_format = inqkit::io::VTIWriteOptions::Format::binary},
        {.overwrite = (START == 0)});

    // FULL energy decomposition — the user asked for every component, and the
    // eight terms below sum to energy_total (nvxc / eigenvalues are diagnostics
    // outside the total).
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.energy_external = sel.energy_nonlocal = sel.energy_ion = true;
    sel.energy_ion_kinetic = sel.energy_exact_exchange = true;
    sel.energy_nvxc = sel.energy_eigenvalues = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs_writer(
        OBS + "/observables" + seg + ".csv", sel);
    obs_writer.write_header();

    inqkit::observables::WPMomentumStats wp_mom(
        OBS + "/wp_momentum_stats" + seg + ".csv", wp_idx,
        {.write_every = Cfg::STATS_EVERY});
    inqkit::observables::WPRealSpaceStats wp_pos(
        OBS + "/wp_real_space_stats" + seg + ".csv", wp_idx,
        {.write_every = Cfg::STATS_EVERY});

    inqkit::observables::StateEnergyWriter state_energy_wr(
        OBS + "/state_energies" + seg + ".csv", /*emit_variance=*/true);
    inqkit::observables::OccupationsWriter occupations_wr(
        OBS + "/occupations_vs_time" + seg + ".csv");

    // ---- pairwise interaction-energy decomposition (P = wavepacket, S = bath
    // electrons, B = neutralising background) -----------------------------
    // BULK NOTE: the background here is UNIFORM, so its Poisson solution is pure
    // G=0 - which INQ drops. phi_plus is therefore IDENTICALLY ZERO, and with it
    // E_SB = E_PB = E_BB = 0. They are still written (as zeros) so the schema
    // matches the localised-jellium runs, but the physics in bulk lives in
    // E_SS (bath-bath), E_PP (projectile SELF-Hartree) and E_PS (projectile-bath).
    //
    // CLOSURE, checkable offline: E_SS + E_PS + E_PP == INQ energy_hartree.
    inq::basis::field<inq::basis::real_space, double> ie_phiplus(
        electrons.density().basis());
    ie_phiplus.fill(0.0);
    std::ofstream ix;
    if (electrons.root()) {
        ix.open(OBS + "/interactions" + seg + ".csv");
        ix << std::setprecision(12)
           << "step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,"
              "e_hartree_check,e_external_check,norm_wp,norm_total\n";
    }

    inqkit::observables::DensityDelta density_delta(
        VTI + "/density_delta", VTI + "/density_delta_coarse",
        {.emit_raw_vti = false, .emit_coarse_vti = true,   // raw delta = n(t)-n(0), exactly recomputable from density_total
         .compute_l2 = true, .coarse_bin_bohr = 3.0});

    // ----- t = 0 analytic gates ------------------------------------------
    if (START == 0) {
        auto m0 = wp_mom.compute(electrons);
        auto r0 = wp_pos.compute(electrons);
        const double T1_0 = m0.ekin;
        const double T2_0 = 0.5 * (m0.px*m0.px + m0.py*m0.py + m0.pz*m0.pz);
        const double expect_T1 =
            0.5 * (Cfg::WP_K0 * Cfg::WP_K0 + 3.0 * sigma_p2);

        int fails = 0;
        auto gate = [&](char const* name, double got, double want, double tol) {
            const bool ok = std::abs(got - want) <= tol;
            std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << name << ": " << got
                      << " (expect " << want << " +/- " << tol << ")\n";
            if (!ok) ++fails;
        };
        std::cout << "\n  --- t=0 analytic gates ---\n";
        // Real-space norm carries an explicit dV and IS ~1. The momentum-space
        // Parseval constant m0.N does NOT (INQ's FFT prefactor sets its scale,
        // ~5e7 here); every moment is divided by it, so it is reported, not gated.
        gate("norm (real space)", r0.N, 1.0, 0.02);
        std::cout << "  [info] momentum-space Parseval constant N = " << m0.N
                  << " (grid/FFT-prefactor dependent; must stay CONSTANT in time)\n";
        gate("<p_z> = k0", m0.pz, Cfg::WP_K0, 0.01);
        gate("sigma_pz^2 = 1/(2 sigma^2)", m0.sz2, sigma_p2, 0.005);
        gate("T1 = <p^2>/2", T1_0, expect_T1, 0.02);
        gate("T2 = <p>^2/2 (Ha)", T2_0, 0.5 * Cfg::WP_K0 * Cfg::WP_K0, 0.02);
        gate("T1-T2 = 3/(4 sigma^2) (eV)", (T1_0 - T2_0) * HA_TO_EV,
             3.0 / (4.0 * Cfg::WP_SIGMA_BOHR * Cfg::WP_SIGMA_BOHR) * HA_TO_EV,
             0.15);
        gate("centroid z (naive)", r0.z, Cfg::WP_CZ_BOHR, 0.05);
        gate("centroid z (circular)", r0.zc, Cfg::WP_CZ_BOHR, 0.05);
        gate("sigma_z (density std)", std::sqrt(r0.sz2),
             Cfg::WP_SIGMA_BOHR / std::sqrt(2.0), 0.05);
        if (fails > 0) {
            std::cerr << "\nFATAL: " << fails
                      << " t=0 gate(s) failed — the injected packet is not the "
                         "one this run claims. Aborting before burning GPU time.\n";
            return 4;
        }
        std::cout << "  all t=0 gates PASSED\n\n";

        auto sys0 = inqkit::fields::density::total(electrons);
        total_wr.write(sys0, 0.0, 0);
        auto wp0 = inqkit::fields::density::orbital(electrons, wp_idx);
        wp_density_wr.write(wp0, 0.0, 0);
    }

    // ----- Real-time callbacks -------------------------------------------
    // Callback fires at OBS_EVERY (energy series); the VTI writes inside are
    // gated on WRITE_EVERY so storage and the energy cadence are independent.
    inqkit::RealTimeSession rt(ions, electrons, Cfg::OBS_EVERY);
    rt.add([&](inqkit::StepContext const& ctx) {
        const int abs_step = ctx.step;
        const bool write_vti = (abs_step % Cfg::WRITE_EVERY == 0);

        auto sys_f = inqkit::fields::density::total(*ctx.electrons);
        if (write_vti) total_wr.write(sys_f, ctx.time_au, abs_step);

        const double l2 = density_delta.snapshot(sys_f, ctx.time_au, abs_step);
        inqkit::StepContext out = ctx;
        out.density_l2 = l2;
        obs_writer.append(out);

        // WP-only density: the field the density GIFs and the induced
        // decomposition are built from. Per-element GPU->host extraction, so it
        // runs at the density cadence, not every step.
        auto wp_dens = inqkit::fields::density::orbital(*ctx.electrons, wp_idx);
        if (write_vti) wp_density_wr.write(wp_dens, ctx.time_au, abs_step);

        // Pairwise decomposition, EVERY callback (cheap: two Poisson solves).
        auto n_wp_f = inqkit::jellium::orbital_density_field(*ctx.electrons, wp_idx);
        auto ct = inqkit::jellium::compute_coulomb_wp(
            ctx.electrons->density(), n_wp_f, ie_phiplus);
        if (ctx.electrons->root())
            ix << abs_step << "," << ctx.time_au << "," << ct.e_ss << ","
               << ct.e_pp << "," << ct.e_ps << "," << ct.e_sb << ","
               << ct.e_pb << "," << 0.0 << "," << ct.e_hartree_check << ","
               << ct.e_external_check << "," << ct.norm_wp << ","
               << ct.norm_total << "\n";

        if (abs_step % Cfg::WF_WRITE_EVERY == 0) {
            auto wf = inqkit::fields::orbital::wavefunction(*ctx.electrons, wp_idx);
            char nm[64];
            std::snprintf(nm, sizeof(nm), "wavefunction_t%06d", abs_step);
            wp_wf_wr.write(wf, std::string(nm));
        }
    });

    real_time::propagate(
        ions, electrons,
        [&](auto const& data) {
            rt.step(data);
            wp_mom.maybe_accumulate(data);
            wp_pos.maybe_accumulate(data);
            if (data.iter() % (10 * Cfg::WRITE_EVERY) == 0) {
                state_energy_wr.snapshot(data);
                occupations_wr.snapshot(data);
            }
            // Interior checkpoint (.claude/rules/checkpoint-dont-block.md): a
            // killed job loses at most CKPT_EVERY steps.
            if (data.iter() > 0 && data.iter() % Cfg::CKPT_EVERY == 0) {
                electrons.save("results/checkpoint");
                jellium::rt_state::write("results/rt_state.txt",
                                         {.last_step = data.iter(),
                                          .time_au   = data.time(),
                                          .dt        = Cfg::DT_AU,
                                          .wp_idx    = wp_idx});
            }
        },
        options::theory{}.lda(),
        options::real_time{}
            .num_steps(N_STEPS)
            .dt(Cfg::DT_AU * 1.0_atomictime)
            .observables_current()
            .observables_dipole(),
        {}, START);

    // ----- Final checkpoint (.claude/rules/final-timestep-checkpoint.md) ---
    electrons.save("results/checkpoint");
    jellium::rt_state::write("results/rt_state.txt",
                             {.last_step = N_STEPS,
                              .time_au   = Cfg::DT_AU * N_STEPS,
                              .dt        = Cfg::DT_AU,
                              .wp_idx    = wp_idx});

    const double wall =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t_wall0).count();

    if (electrons.root()) {
        std::ofstream s("results/run_summary.txt");
        s << std::setprecision(16);
        s << "RUN SUMMARY\n===========\n\n"
          << "1. Run identity\n---------------\n"
          << "run_name        = " << RUN_NAME << "\n"
          << "run_type        = wave-packet projectile, bulk jellium TDDFT (ALDA)\n"
          << "plan            = docs/plans/bulk-jellium-ks-stopping.md\n"
          << "date_finished   = " << iso_now() << "\n"
          << "wall_time_s     = " << wall << "\n"
          << "start_step      = " << START << "\n\n"
          << "3. System configuration\n-----------------------\n"
          << "cell_bohr       = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR
          << " x " << Cfg::LZ_BOHR << "\n"
          << "cell_geometry   = orthorhombic, periodic xyz\n"
          << "n_electrons     = " << Cfg::N_ELECTRONS << "\n"
          << "n_states        = " << n_states << "\n"
          << "extra_states    = " << Cfg::EXTRA_STATES << "\n"
          << "wp_state_index  = " << wp_idx << "\n"
          << "spacing_bohr    = " << Cfg::SPACING_BOHR << "\n"
          << "cutoff_ha       = " << Cfg::CUTOFF_HA << "\n"
          << "xc_functional   = LDA (ALDA in TDDFT)\n"
          << "gs_checkpoint   = " << GS_DIR << "\n\n"
          << "5. Wavepacket configuration\n---------------------------\n"
          << "wp_enabled      = yes\n"
          << "wp_center_bohr  = " << Cfg::WP_CX_BOHR << " " << Cfg::WP_CY_BOHR
          << " " << Cfg::WP_CZ_BOHR << "\n"
          << "wp_sigma_bohr   = " << Cfg::WP_SIGMA_BOHR << "\n"
          << "wp_sigma_note   = wavepacket sigma (psi width); density std = this/sqrt2\n"
          << "wp_sigma_density= " << Cfg::WP_SIGMA_BOHR / std::sqrt(2.0) << "\n"
          << "wp_k0_bohr_inv  = " << Cfg::WP_KZ << "\n"
          << "wp_energy_ev    = " << Cfg::WP_EKIN_EV << "\n"
          << "norm_after      = " << report.norm_after << "\n"
          << "max_overlap     = " << report.max_overlap << "\n\n"
          << "6. Real-time configuration\n--------------------------\n"
          << "rt_num_steps    = " << N_STEPS << "\n"
          << "dt_au           = " << Cfg::DT_AU << "\n"
          << "total_time_au   = " << (Cfg::DT_AU * N_STEPS) << "\n"
          << "write_every     = " << Cfg::WRITE_EVERY << "\n"
          << "wf_write_every  = " << Cfg::WF_WRITE_EVERY << "\n"
          << "stats_every     = " << Cfg::STATS_EVERY << "\n\n"
          << "7. Analysis window\n------------------\n"
          << "fit_t0_au       = " << Cfg::FIT_T0_AU << "\n"
          << "fit_t1_au       = " << Cfg::FIT_T1_AU << "\n"
          << "t_ifw_au        = " << Cfg::T_IFW_AU << "\n"
          << "t_ifw_note      = dispersion-aware; the static rule would say 24.34\n"
          << "t_transverse_au = " << Cfg::T_TRANSVERSE_AU << "\n"
          << "plasma_period_au= 49.39\n"
          << "plasma_period_note = EXCEEDS the run: S is initial-drag, not "
             "steady-state (see plan section 3.4)\n\n"
          << "9. End-of-run diagnostics\n-------------------------\n"
          << "run_completed   = true\n";
    }

    std::cout << "Done. Wall time " << wall << " s.\n";
    return 0;
}
