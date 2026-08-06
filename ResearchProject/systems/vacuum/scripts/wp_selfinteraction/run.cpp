// ============================================================================
// vacuum / scripts/wp_selfinteraction / run.cpp
//
// HOW MUCH DOES THE SELF-HARTREE TERM SPREAD A WAVEPACKET, AND WHAT DOES IT COST
// IN ENERGY?
//
// Plan: docs/plans/wp-self-interaction-correction.md
// Motivating measurement: docs/handovers/cylindrical-channeling-ks-stopping.md
// (2026-08-02, 5th) — the channeling wavepacket expanded 1.47x faster than free
// dispersion and its stopping power sat ~20 % below its classical twin.
//
// THE IDEA. A single electron has, exactly, NO self-interaction. So for ONE
// electron alone in a vacuum box the exact answer is free-particle dispersion,
// known in closed form. Running the SAME initial Gaussian at three theory levels
// therefore measures the self-interaction directly, by DIFFERENCE, with no
// self-interaction correction implemented at all:
//
//   WP_THEORY=noninteracting   no self-interaction         <- the reference
//   WP_THEORY=hartree          + Hartree self-interaction
//   WP_THEORY=lda              + Hartree AND LDA xc self-interaction
//
//   (lda - noninteracting)   = the TOTAL self-interaction error
//   (hartree - noninteracting) = its HARTREE part alone
//   (lda - hartree)          = its XC part alone
//
// The last split matters for the correction design: it says whether removing the
// Hartree self-term is sufficient or whether the xc self-term must go too.
//
// WHY THE INITIAL STATE IS IDENTICAL ACROSS ALL THREE. The ground state is ALWAYS
// computed non-interacting and is then OVERWRITTEN by inject_into_last_extra_state
// (extra_states(0) + extra_electrons(1.0) => exactly one state, which the packet
// replaces). So the three runs differ ONLY in the theory used to PROPAGATE. Any
// divergence between them is self-interaction and nothing else.
//
// WHY THE PACKET IS STATIONARY (k0 = 0 by default). The question is about
// SPREADING, and for a free particle spreading is frame-independent — the
// self-Hartree depends only on the packet's shape in its own rest frame. A
// stationary packet needs no traversal length, so the box can be small enough to
// afford a fine grid, and it cannot wrap in z. Set WP_K0 to reproduce the moving
// case as a cross-check.
//
// EXACT REFERENCE (free evolution, m_e = 1, atomic units). For
// psi_0 ~ exp(-r^2/(2 sigma^2)) the DENSITY |psi|^2 has per-axis standard
// deviation
//     sigma_dens(t) = sqrt( sigma^2/2 + t^2/(2 sigma^2) )
// and the momentum distribution NEVER changes:
//     <p> = k0   and   var(p_d) = 1/(2 sigma^2)   for all t.
// var(p) is the sharpest gate of the three: it is EXACTLY conserved under free
// evolution, so any growth is interaction, with no discretisation excuse.
//
// GEOMETRY. sigma_WP = 4 Bohr, matched to the channeling run. Cubic box L = 72
// Bohr at h = 0.5 => 144^3. At t = 30 the free density width is 6.01 Bohr, so the
// box half-width is ~6 sigma; if self-interaction widens the packet substantially
// beyond that, the naive and circular second moments in wp_real_space_stats.csv
// diverge, which is the built-in wrap diagnostic (no extra observable needed).
//
// CUTOFF GUARD. sigma_p = 1/(sqrt2 sigma) = 0.177, so the packet occupies
// |k| <~ k0 + 4 sigma_p. k_max = pi/h = 6.28 — enormous headroom.
//
// ---------------------------------------------------------------------------
// TIER V OF THE SIC CAMPAIGN (added 2026-08-02, plan §4 after review). The same
// binary now also RUNS the correction, not only the uncorrected difference:
//
//   WP_SIC=none    (default) the three-theory difference measurement above
//   WP_SIC=h       subtract Q v_H[n_wp] Q            (SIC-H)
//   WP_SIC=pzrun   subtract Q (v_H + v_xc^unpol) Q   (SIC-PZrun)
//
// SIC requires WP_THEORY=lda (it removes terms the LDA Hamiltonian contains).
// Strang scheduling: half-kick at the step-0 callback, full kick after every
// interior step, half-kick at the final step. In vacuum Q = 1 (no occupied
// bath), so this tier validates the KICK + BOOKKEEPING only; the projection is
// validated by test_wp_sic_engine and the jellium Tier B (plan §0/D6).
//
// CLOSED-FORM GATES (evaluated at end of run, exit 4 on hard failure):
//   pzrun: var(p_z) drift < 0.1 %; sigma_dens(t_end) within 0.5 % of analytic;
//          |<p_z> - k0| < 1e-4; E_corrected drift PASS < 1e-5 eV
//          (WARN < 1e-3 eV — reported, chain continues; FAIL above).
//   h:     |<p_z> - k0| < 1e-4 (zero-force) and the E_corrected ladder; sigma
//          is EXPECTED to under-spread (v_xc self-binding remains) — reported,
//          not gated.
// Expected qualitative ordering across runs: lda spreads FASTEST (full SIE),
// sic-h SLOWEST (xc binding uncancelled), sic-pzrun == free (the gate).
//
// Env: WP_THEORY(lda) WP_SIC(none) WP_SIGMA(4.0) WP_K0(0.0) WP_L(72) WP_H(0.5)
//      WP_DT(0.02) WP_NSTEPS(1500) WP_SAVE_EVERY(150) WP_WF_EVERY(150) WP_OUT
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/jellium/interaction_energies.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/wavepacket/self_interaction_correction.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
namespace fs = std::filesystem;
namespace obs_ = inqkit::observables;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }
static std::string tag6(int s){ std::ostringstream o; o<<std::setw(6)<<std::setfill('0')<<s; return o.str(); }

int main() {
    const double HA_TO_EV = 27.211386245988;

    const std::string THEORY = env_s("WP_THEORY", "lda");
    const std::string SICSTR = env_s("WP_SIC", "none");
    const double SIGMA   = env_d("WP_SIGMA", 4.0);          // matched to channeling_twin
    const double K0      = env_d("WP_K0", 0.0);             // stationary by default
    const double L       = env_d("WP_L", 72.0);
    const double H       = env_d("WP_H", 0.5);
    const double DT      = env_d("WP_DT", 0.02);
    const int    NSTEPS  = env_i("WP_NSTEPS", 1500);
    const int    SAVE    = env_i("WP_SAVE_EVERY", 150);
    const int    WF      = env_i("WP_WF_EVERY", 150);
    const std::string OUT = "results/" + env_s("WP_OUT", THEORY);

    if (THEORY != "noninteracting" && THEORY != "hartree" && THEORY != "lda") {
        std::cerr << "FATAL: WP_THEORY must be noninteracting | hartree | lda; got '"
                  << THEORY << "'\n";
        return 2;
    }
    using SIC = inqkit::SelfInteractionCorrection;
    const auto SIC_MODE = SIC::mode_from_string(SICSTR);
    if (SIC_MODE != SIC::Mode::none && THEORY != "lda") {
        std::cerr << "FATAL: WP_SIC=" << SICSTR << " requires WP_THEORY=lda (the "
                     "correction removes terms only the LDA Hamiltonian contains).\n";
        return 2;
    }

    // --- analytic references, computed up front so they can be printed and gated
    const double sigma_dens0 = SIGMA / std::sqrt(2.0);       // density std at t=0
    const double var_p_free  = 1.0 / (2.0 * SIGMA * SIGMA);  // per Cartesian axis
    const double t_end       = NSTEPS * DT;
    const double sigma_end   = std::sqrt(SIGMA*SIGMA/2.0 + t_end*t_end/(2.0*SIGMA*SIGMA));
    // Localisation energy (3/2) sum_d var(p_d)/2 = 3/(4 sigma^2), the t=0 value of
    // T2 - T1 in the channeling analysis. Constant under free evolution.
    const double loc_ev      = 3.0 / (4.0 * SIGMA * SIGMA) * HA_TO_EV;

    fs::create_directories(OUT + "/raw/observables");
    fs::create_directories(OUT + "/raw/vti/density_total");
    fs::create_directories(OUT + "/raw/vti/wavefunction_wp");

    std::cout << std::setprecision(8)
        << "\n=== wp_selfinteraction (out=" << OUT << ") ===\n"
        << "  theory      = " << THEORY << "   <- the ONLY difference between runs\n"
        << "  sigma_WP    = " << SIGMA << " Bohr   k0 = " << K0 << "\n"
        << "  box         = " << L << "^3 Bohr   h = " << H
        << "   (" << int(std::lround(L/H)) << "^3 grid)\n"
        << "  dt          = " << DT << "   n_steps = " << NSTEPS
        << "   t_end = " << t_end << " a.u.\n"
        << "  ANALYTIC FREE REFERENCE:\n"
        << "    sigma_dens(0)   = " << sigma_dens0 << " Bohr\n"
        << "    sigma_dens(end) = " << sigma_end   << " Bohr   (box half-width "
        << L/2.0 << " = " << (L/2.0)/sigma_end << " sigma)\n"
        << "    var(p_d)        = " << var_p_free << "  (CONSTANT under free evolution)\n"
        << "    3/(4 sigma^2)   = " << loc_ev << " eV\n";

    auto cell = systems::cell::cubic(L*1.0_b).periodic();
    auto ions = systems::ions(cell);
    // Exactly ONE electron in ONE state; the packet REPLACES it, so
    // density_total == density_wp and there is no background gas.
    auto electrons = systems::electrons(
        ions, options::electrons{}.spacing(H*1.0_b).extra_states(0).extra_electrons(1.0));

    // The ground state is ALWAYS non-interacting and is discarded by the
    // injection below. Keeping it identical across the three runs is what makes
    // this a controlled comparison.
    ground_state::initial_guess(ions, electrons);
    ground_state::calculate(ions, electrons, options::theory{}.non_interacting(),
                            options::ground_state{}.energy_tolerance(1.0e-8_Ha).max_steps(200));

    auto rep = inqkit::WavePacket{}
                   .center(0.0, 0.0, 0.0)
                   .sigma(SIGMA)
                   .k0(0.0, 0.0, K0)
                   .minimum_image(true)
                   .inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = rep.state_index;
    std::cout << "  WP injected: state_index=" << wp_idx
              << "  norm_after=" << rep.norm_after << "\n";

    // PER-STATE MASS (added 2026-08-05 for the Nazarov-Gross mass ladder,
    // docs/plans/nazarov-gross-slab-mass-ladder.md Phase 1, steps 5-6).
    //
    // WHY IT BELONGS HERE. That campaign's whole signal is width-mediated, and
    // LDA self-interaction ALSO widens the packet. If the SIE-driven excess width
    // is mass-DEPENDENT it forges the campaign's result with the right sign, so
    // it has to be measured, and this vacuum box is the only place the exact
    // answer is known (a lone electron has no self-interaction, so free
    // dispersion is the truth). WP_INV_MASS = 1/M reproduces the DEFAULT scalar
    // path exactly (INQ keeps the scalar branch while every inverse mass is 1.0),
    // so every previously published run in this directory is bit-unchanged.
    //
    // The analytic reference generalises as sigma_dens(t) =
    // sqrt(sigma^2/2 + t^2/(2 M^2 sigma^2)) and var(p) stays EXACTLY 1/(2 sigma^2)
    // for any M — var(p) is a property of the initial state, not of the mass.
    const double INV_MASS = env_d("WP_INV_MASS", 1.0);
    if (INV_MASS != 1.0) {
        electrons.inverse_mass()[0][wp_idx] = INV_MASS;
        std::cout << "  inverse_mass[" << wp_idx << "]=" << INV_MASS
                  << "  (M = " << 1.0/INV_MASS << ")\n";
    }

    const SIC sic(SIC_MODE, wp_idx);
    const bool sic_on = (SIC_MODE != SIC::Mode::none);
    if (sic_on)
        std::cout << "  SIC ENABLED: " << SIC::mode_name(SIC_MODE)
                  << "  (Strang: half kick at step 0 and step " << NSTEPS
                  << ", full kicks between)\n";

    // ---- observables --------------------------------------------------------
    std::ofstream en(OUT + "/raw/observables/energies.csv");
    en << std::setprecision(12)
       << "step,time_au,total,kinetic,hartree,external,non_local,xc,"
          "exact_exchange,ion,ion_kinetic\n";

    // interactions.csv (.claude/rules/decomposed-interaction-energies.md).
    // VACUUM SPECIAL CASE: there is no bath and no background, so
    //   n_S = 0,  phi_+ = 0  =>  E_SS = E_PS = E_SB = E_PB = E_BB = 0
    // and E_PP is the ONLY non-zero term — which is exactly the quantity this
    // run exists to measure. The columns are still written so the schema matches
    // every other system.
    //
    // CLOSURE GATE, and it is a strong one: the wavepacket is the only charge in
    // the box, so INQ's own energy_hartree MUST equal our offline E_PP for the
    // hartree and lda runs (and be identically 0 for noninteracting, while the
    // offline E_PP stays non-zero as a pure diagnostic of the packet's size).
    std::ofstream ie(OUT + "/raw/observables/interactions.csv");
    ie << std::setprecision(12)
       << "step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,"
          "e_hartree_inq,closure_pp_minus_hartree,norm_wp\n";

    // SIC per-step ledger. E_corrected is the variant's conserved functional
    // (exact in vacuum: Q = 1, plan §0/D2); u_self/exc_self are the energies
    // the correction removes; max_overlap_pre/norm_removed are projection
    // diagnostics (identically 0 here — asserted by their column staying 0).
    std::ofstream sc(OUT + "/raw/observables/sic.csv");
    sc << std::setprecision(12)
       << "step,time_au,u_self_ha,exc_self_ha,e_corrected_ha,"
          "max_overlap_pre,norm_removed,cum_norm_removed\n";
    double e_corr_first = 0.0, e_corr_last = 0.0, cum_removed = 0.0;
    double max_ov_run = 0.0; bool e_corr_seen = false;

    obs_::MomentumDistribution mom_dist(
        OUT + "/raw/observables/momentum_distribution.csv", wp_idx, L,
        {.n_bins = 64, .k_max_bohr_inv = 0.0, .write_every = 1});
    obs_::WPMomentumStats  wp_mom(OUT + "/raw/observables/wp_momentum_stats.csv",
                                  wp_idx, {.write_every = 1});
    obs_::WPRealSpaceStats wp_rs (OUT + "/raw/observables/wp_real_space_stats.csv",
                                  wp_idx, {.write_every = 1});

    using RLay = inqkit::io::RealField3DLayout;
    const auto vtifmt = inqkit::io::VTIWriteOptions::Format::binary;
    RLay lay_tot{.field_name="density_total", .include_meta=false, .emit_raw=false,
                 .emit_vti=true, .vti_format=vtifmt};
    inqkit::io::RealField3DWriter tot_wr(OUT + "/raw/vti/density_total", lay_tot,
                                         {.overwrite=true});
    // The wavefunction itself, not just its density — needed for the 2-D
    // (k_z, k_perp) momentum maps (inqview kz_kperp_map).
    inqkit::io::ComplexField3DWriter wf_wr(
        OUT + "/raw/vti/wavefunction_wp",
        {.field_name="wavefunction", .include_meta=false, .emit_raw=false,
         .emit_vti=true, .vti_format=vtifmt},
        {.overwrite=true});

    tot_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);
    wf_wr.write(inqkit::fields::orbital::wavefunction(electrons, wp_idx),
                "wavefunction_t" + tag6(0));

    auto step_fn = [&](auto const& data) {
        const int step = data.iter();
        auto e = data.energy();

        // pairwise decomposition: one Poisson solve on the packet's own density
        auto n_wp   = inqkit::jellium::orbital_density_field(electrons, wp_idx);
        auto phi_wp = inq::solvers::poisson::solve(n_wp);
        const double e_pp   = 0.5 * inq::operations::integral_product(n_wp, phi_wp);
        const double norm_p = inq::operations::integral(n_wp);

        if (data.root()) {
            en << step << ',' << data.time() << ',' << e.total() << ',' << e.kinetic()
               << ',' << e.hartree() << ',' << e.external() << ',' << e.non_local()
               << ',' << e.xc() << ',' << e.exact_exchange() << ',' << e.ion()
               << ',' << e.ion_kinetic() << '\n';
            ie << step << ',' << data.time()
               << ",0," << e_pp << ",0,0,0,0,"
               << e.hartree() << ',' << (e_pp - e.hartree()) << ',' << norm_p << '\n';
        }

        mom_dist.maybe_accumulate(data);
        wp_mom.maybe_accumulate(data);
        wp_rs.maybe_accumulate(data);

        if (SAVE > 0 && step % SAVE == 0)
            tot_wr.write(inqkit::fields::density::total(electrons), data.time(), step);
        if (WF > 0 && step % WF == 0)
            wf_wr.write(inqkit::fields::orbital::wavefunction(electrons, wp_idx),
                        "wavefunction_t" + tag6(step));

        // ---- SIC kick, AFTER every observable (uniform pre-kick convention).
        // Strang factors: 1/2 at the step-0 callback (opening), 1/2 at the
        // final step (closing), 1 between. A real kick leaves the density
        // untouched, so nothing above is invalidated.
        if (sic_on) {
            const double f = (step == 0 || step == NSTEPS) ? 0.5 : 1.0;
            auto srep = sic.apply(electrons, f * DT);
            const double e_corr = sic.corrected_energy(e.total(), srep);
            cum_removed += srep.norm_removed;
            max_ov_run = std::max(max_ov_run, srep.max_overlap_pre);
            if (!e_corr_seen) { e_corr_first = e_corr; e_corr_seen = true; }
            e_corr_last = e_corr;
            if (data.root())
                sc << step << ',' << data.time() << ',' << srep.u_self << ','
                   << srep.exc_self << ',' << e_corr << ',' << srep.max_overlap_pre
                   << ',' << srep.norm_removed << ',' << cum_removed << '\n';
        }
    };

    auto theory = options::theory{}.non_interacting();
    if (THEORY == "hartree") theory = options::theory{}.hartree();
    else if (THEORY == "lda") theory = options::theory{}.lda();

    const auto wall0 = std::chrono::steady_clock::now();
    real_time::propagate(ions, electrons, step_fn, theory,
        options::real_time{}.num_steps(NSTEPS).dt(DT*1.0_atomictime));
    const double wall_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall0).count();
    en.close();
    ie.close();
    sc.close();

    // ---- SIC closed-form gates (plan §4, amended §0/D1-D3) -------------------
    // The final state carries the closing half-kick; the analytic references are
    // insensitive to that O(dt v_SIC) phase at these tolerances.
    int gate_fails = 0;
    std::string gate_lines;
    if (sic_on) {
        auto pF = wp_mom.compute(electrons);
        auto rF = wp_rs.compute(electrons);
        const double HA = HA_TO_EV;
        std::ostringstream g;
        g << std::setprecision(10);
        auto gate = [&](bool hard, char const* nm, double got, double want,
                        double tol, bool relative) {
            const double dev = relative ? (want != 0.0 ? (got - want)/std::abs(want) : got)
                                        : (got - want);
            const bool ok = std::abs(dev) <= tol;
            g << (ok ? "  [PASS] " : (hard ? "  [FAIL] " : "  [WARN] "))
              << nm << ": " << got << "  (expect " << want << ", dev " << dev
              << ", tol " << tol << (relative ? " rel" : " abs") << ")\n";
            if (!ok && hard) ++gate_fails;
        };
        const bool pz = (SIC_MODE == SIC::Mode::pz_run);
        g << "\n  --- SIC end-of-run gates (" << SIC::mode_name(SIC_MODE) << ") ---\n";
        // zero-force: both variants
        gate(true,  "<p_z> - k0 (abs)", pF.pz, K0, 1e-4, false);
        if (pz) {
            gate(true, "var(p_z) vs 1/(2s^2) (rel)", pF.sz2, var_p_free, 1e-3, true);
            gate(true, "sigma_dens(t_end) vs analytic (rel)", rF.szc, sigma_end, 5e-3, true);
        } else {
            g << "  [info] sigma_dens(t_end) = " << rF.szc << " vs free " << sigma_end
              << "  (SIC-H EXPECTED to under-spread: xc self-binding remains)\n";
            g << "  [info] var(p_z) = " << pF.sz2 << " vs free " << var_p_free << "\n";
        }
        // E_corrected ladder: PASS < 1e-5 eV, WARN < 1e-3 eV, FAIL above.
        const double dE = (e_corr_last - e_corr_first) * HA;
        if (std::abs(dE) <= 1e-5)
            g << "  [PASS] E_corrected drift = " << dE << " eV (< 1e-5)\n";
        else if (std::abs(dE) <= 1e-3)
            g << "  [WARN] E_corrected drift = " << dE << " eV (1e-5..1e-3: split-"
                 "operator residual, reported, chain continues)\n";
        else { g << "  [FAIL] E_corrected drift = " << dE << " eV (> 1e-3)\n"; ++gate_fails; }
        g << "  [info] max_overlap_pre over run = " << max_ov_run
          << "   cum_norm_removed = " << cum_removed
          << "   (vacuum: both must be ~0; occupied manifold is empty)\n";
        gate_lines = g.str();
        std::cout << gate_lines;
    }

    // final checkpoint (.claude/rules/final-timestep-checkpoint.md)
    electrons.save(OUT + "/checkpoint");
    {
        std::ofstream rt(OUT + "/rt_state.txt");
        rt << std::setprecision(17)
           << "last_step=" << NSTEPS << "\ntime_au=" << NSTEPS*DT
           << "\ndt=" << DT << "\nwp_idx=" << wp_idx << "\n";
    }

    std::ofstream sum(OUT + "/run_summary.txt");
    sum << std::setprecision(12)
        << "run = vacuum/wp_selfinteraction/" << env_s("WP_OUT", THEORY) << "\n"
        << "purpose = quantify the self-interaction of ONE electron in vacuum by "
           "comparing theory levels on an identical initial packet; WP_SIC != none "
           "additionally REMOVES it and gates against the closed-form free evolution\n"
        << "engine = inq-study\n"
        << "theory = " << THEORY << "\n"
        << "sic_mode = " << SIC::mode_name(SIC_MODE) << "\n"
        << "gs_theory = non_interacting (overwritten by the injected packet)\n"
        << "wp = gaussian sigma " << SIGMA << " k0 " << K0 << " mass "
        << (1.0/INV_MASS) << ", centred at origin\n"
        << "wp_inverse_mass = " << INV_MASS << "\n"
        << "wp_state_index = " << wp_idx << "  norm_after_injection = " << rep.norm_after << "\n"
        << "cell_bohr = " << L << "^3 cubic periodic   spacing = " << H << "\n"
        << "dt_au = " << DT << "  n_steps = " << NSTEPS << "  t_end_au = " << t_end << "\n"
        << "save_every = " << SAVE << "  wf_every = " << WF << "\n"
        << "analytic_sigma_dens_0 = " << sigma_dens0 << "\n"
        << "analytic_sigma_dens_end = " << sigma_end << "\n"
        << "analytic_var_p_per_axis = " << var_p_free << "  (constant under free evolution)\n"
        << "analytic_localisation_ev = " << loc_ev << "\n"
        << "box_half_widths_of_final_sigma = " << (L/2.0)/sigma_end << "\n"
        << "propagate_wall_s = " << wall_s
        << "  per_step_ms = " << 1000.0*wall_s/NSTEPS << "\n";
    if (sic_on) {
        sum << "sic_e_corrected_first_ha = " << e_corr_first << "\n"
            << "sic_e_corrected_last_ha = " << e_corr_last << "\n"
            << "sic_e_corrected_drift_ev = " << (e_corr_last - e_corr_first)*HA_TO_EV << "\n"
            << "sic_max_overlap_pre = " << max_ov_run << "\n"
            << "sic_cum_norm_removed = " << cum_removed << "\n"
            << "sic_gate_failures = " << gate_fails << "\n"
            << "sic_gates:\n" << gate_lines;
    }
    sum << "run_completed = true\n";
    sum.close();

    if (gate_fails > 0) {
        std::cerr << "\nFATAL: " << gate_fails << " SIC gate(s) failed — the "
                     "correction is not doing what the derivation says. Do not "
                     "run production (plan §4 decision rule).\n";
        return 4;
    }
    std::cout << "  done -> " << OUT << "/   (" << wall_s << " s)\n";
    return 0;
}
