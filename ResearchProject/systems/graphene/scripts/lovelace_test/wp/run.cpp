// ============================================================================
// graphene/scripts/lovelace_test/wp/run.cpp
//
// WAVEPACKET half of the lovelace_test graphene stopping campaign.
// Adapted from twodef_sv/wp/run.cpp; key changes:
//   - spacing(DX) instead of cutoff(50 Ha); default DX = 0.5 Bohr
//   - Lz = 90 Bohr, CAP = 10 Bohr/face ([−45,−35] and [+35,+45])
//   - WP launch z = −19 Bohr (midway between CAP inner edge and entry layer)
//   - sigma_WP = 4.0 Bohr, K0 = 3.320 (150 eV), dt = 0.04, N_STEPS = 500
//
// Observables written:
//   observables.csv    : energy_total, energy_kinetic (KS KE), energy_hartree,
//                        energy_xc, energy_external, energy_nonlocal,
//                        energy_ion, energy_ion_kinetic (Ehrenfest ionic KE)
//   wp_momentum_stats  : ekin = <p²>/2m,  pz/px/py = components of <p>,
//                        sz2 = var(pz)  → drift KE = (pz²+px²+py²)/2
//   wp_real_space_stats: WP centroid and spatial spread vs time
//   interactions.csv   : pairwise P/S/B decomposition (E_SS,E_PP,E_PS; E_SB=E_PB=E_BB=0)
//
// Env: GR_VARIANT(bi|mono)  GR_GEOM(path)  GR_LZ_BOHR(90)  GR_DX_BOHR(0.5)
//      LJ_K0(3.320)  LJ_SIGMA(4.0)  LJ_LAUNCH_Z(-19)
//      LJ_IMPACT_X(2.3244)  LJ_IMPACT_Y(1.3420)
//      LJ_DT(0.04)  LJ_N_STEPS(500)  LJ_CAP_ETA(-1.0)  LJ_CAP_L(10.0)
//      LJ_SAVE_EVERY(10)  LJ_WF_EVERY(100)  LJ_CKPT_EVERY(0=auto)
//      LJ_RESUME(0)  LJ_OUT(REQUIRED)  LJ_GS_DIR(REQUIRED)
//
// Requires inq-study (CAP compile; set INQ_SOURCE=<tddft_root>/inq-study).
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
#include <inqkit/jellium/interaction_energies.hpp>

#include "../../../shared/configs/twodef_gs.hpp"

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
namespace Cfg = graphene_twodef;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

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
    char b[64]; std::strftime(b, sizeof(b), "%Y-%m-%dT%H:%M:%S", &tm); return b;
}

int main() {
    auto t_wall0 = std::chrono::steady_clock::now();
    const double HA = 27.211386245988;

    const std::string VARIANT = env_s("GR_VARIANT", "bi");
    const bool  BI      = (VARIANT == "bi");
    const int   N_ATOMS = BI ? Cfg::N_C_BI : Cfg::N_C_MONO;
    const int   N_ELEC  = BI ? Cfg::N_ELEC_BI : Cfg::N_ELEC_MONO;
    const int   EXTRA   = env_i("GR_EXTRA", BI ? Cfg::EXTRA_STATES_BI : Cfg::EXTRA_STATES_MONO);

    const std::string LOCAL_ROOT = std::string(
        std::getenv("TDDFT_ROOT") ? std::getenv("TDDFT_ROOT") : "/local/data/public/skcb2/tddft");
    const std::string DEF_GEOM = LOCAL_ROOT + "/ResearchProject/systems/graphene/shared/geometry/"
        + (BI ? "graphene_3x2_bilayer.xyz" : "graphene_3x2.xyz");
    const std::string GEOM = env_s("GR_GEOM", DEF_GEOM);

    const double LZ       = env_d("GR_LZ_BOHR", 90.0);
    const double DX       = env_d("GR_DX_BOHR",  0.5);
    const double SIGMA_WP = env_d("LJ_SIGMA",    4.0);
    const double K0       = env_d("LJ_K0",       3.3204);  // 150 eV
    const double LAUNCH_Z = env_d("LJ_LAUNCH_Z", -19.0);
    const double IMPACT_X = env_d("LJ_IMPACT_X",  2.3244);
    const double IMPACT_Y = env_d("LJ_IMPACT_Y",  1.3420);
    const double DT       = env_d("LJ_DT",        0.04);
    const double CAP_ETA  = env_d("LJ_CAP_ETA",  -1.0);
    const double CAP_L    = env_d("LJ_CAP_L",     10.0);
    const bool   RESUME   = env_i("LJ_RESUME",    0) != 0;
    const std::string OUT    = "results/" + env_s("LJ_OUT", "");
    const std::string GS_DIR = env_s("LJ_GS_DIR", "");

    const int N_STEPS  = env_i("LJ_N_STEPS", 500);  // 20 a.u. / 0.04 dt
    const int SAVE_EVERY = env_i("LJ_SAVE_EVERY", std::max(1, N_STEPS / 50));
    const int WF_EVERY   = env_i("LJ_WF_EVERY", 100);
    int CKPT_EVERY = env_i("LJ_CKPT_EVERY", 0);
    if (CKPT_EVERY <= 0) CKPT_EVERY = std::max(1, N_STEPS / 3);

    const double sigma_p2 = 1.0 / (2.0 * SIGMA_WP * SIGMA_WP);

    if (env_s("LJ_OUT","").empty()) { std::cerr << "FATAL: LJ_OUT unset\n"; return 2; }
    if (GS_DIR.empty() || !std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: LJ_GS_DIR missing/not found: '" << GS_DIR << "'\n"; return 2; }
    if (CAP_ETA >= 0.0 || CAP_L <= 0.0) {
        std::cerr << "FATAL: CAP required (LJ_CAP_ETA<0, LJ_CAP_L>0)\n"; return 2; }

    const double CAP_WIDTH_FRAC = CAP_L / LZ;
    const double CAP_MID_FRAC   = 0.5 - CAP_WIDTH_FRAC / 2.0;
    const double z_cap_in       = LZ / 2.0 - CAP_L;

    const std::string CKPT = OUT + "/checkpoint", RT_STATE = OUT + "/rt_state.txt";
    int START = 0; int wp_idx = -1;
    if (RESUME) {
        START  = (int)read_state_d(RT_STATE, "last_step", -1);
        wp_idx = (int)read_state_d(RT_STATE, "wp_idx",    -1);
        if (START < 0 || wp_idx < 0) {
            std::cerr << "FATAL: LJ_RESUME=1 but no readable " << RT_STATE << "\n"; return 2; }
        if (START >= N_STEPS) { std::cout << "Already at/after target; nothing to do.\n"; return 0; }
    }
    const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

    std::cout << std::setprecision(10)
              << "\n=== lovelace_test WP [" << VARIANT << "]  OUT=" << OUT << " ===\n"
              << "  cell     = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR << " x " << LZ
              << " Bohr, periodicity(2), spacing=" << DX << " Bohr\n"
              << "  ions     = " << N_ATOMS << " C (Ehrenfest), " << N_ELEC << " valence e-\n"
              << "  WP       = sigma_WP " << SIGMA_WP << "  k0=" << K0
              << "  E_drift=" << 0.5*K0*K0*HA << " eV  launch=("
              << IMPACT_X << "," << IMPACT_Y << "," << LAUNCH_Z << ")\n"
              << "  CAP      = eta=" << CAP_ETA << " Ha  W=" << CAP_L
              << " Bohr/face  inner-edge=+/-" << z_cap_in << " Bohr\n"
              << "  dt=" << DT << "  N_STEPS=" << N_STEPS << "  total_time="
              << DT * N_STEPS << " a.u." << (RESUME ? "  [RESUME]" : "") << "\n"
              << "  cadence  : density/" << SAVE_EVERY << "  wavefn/" << WF_EVERY
              << "  stats/1  ckpt/" << CKPT_EVERY << "\n"
              << "  GS       = " << GS_DIR << "\n\n";

    // ---- cell and ions -------------------------------------------------------
    auto cell = systems::cell::orthorhombic(
        Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, LZ * 1.0_b).periodicity(2);
    auto ions = systems::ions::parse(GEOM, cell);
    if (int(ions.size()) != N_ATOMS) {
        std::cerr << "FATAL: " << ions.size() << " atoms, expected " << N_ATOMS << "\n"; return 2; }

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(DX * 1.0_b)
            .extra_states(EXTRA)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());

    // ---- output skeleton -----------------------------------------------------
    const std::string OBS = OUT + "/raw/observables";
    const std::string VTI = OUT + "/raw/vti";
    for (auto const& d : {OBS, VTI+"/density_total", VTI+"/density_wp",
                          VTI+"/density_delta", VTI+"/density_delta_coarse",
                          VTI+"/density_gs_system", VTI+"/wavefunction_wp"})
        std::filesystem::create_directories(d);

    inqkit::InjectionReport report{};

    if (RESUME) {
        electrons.load(CKPT);
        std::cout << "  RESUMED from step " << START << " (wp_idx=" << wp_idx << ")\n";
    } else {
        electrons.load(GS_DIR);
        inqkit::io::RealField3DLayout lay{.field_name="density",.include_meta=false,
            .emit_raw=false,.emit_vti=true,.vti_format=inqkit::io::VTIWriteOptions::Format::binary};
        inqkit::io::RealField3DWriter gs_wr(VTI+"/density_gs_system", lay, {.overwrite=true});
        gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system");

        auto wp = inqkit::WavePacket{}
                      .center(IMPACT_X, IMPACT_Y, LAUNCH_Z)
                      .sigma(SIGMA_WP)
                      .k0(0.0, 0.0, K0)
                      .orthogonalise_against_occupied(electrons);
        report = wp.inject_into_last_extra_state(electrons, 1.0);
        wp_idx = report.state_index;
        std::cout << "  WP injected: state=" << wp_idx
                  << "  norm_after=" << report.norm_after
                  << "  max_overlap=" << report.max_overlap << "\n";

        std::ofstream f(OBS + "/wp_config.txt");
        f << std::setprecision(16)
          << "wp_center_bohr = " << IMPACT_X << " " << IMPACT_Y << " " << LAUNCH_Z << "\n"
          << "wp_sigma_bohr  = " << SIGMA_WP  << "\n"
          << "wp_sigma_density = " << SIGMA_WP/std::sqrt(2.0) << "\n"
          << "wp_k0_bohr_inv = " << K0 << "\n"
          << "wp_state_index = " << wp_idx << "\n"
          << "norm_after     = " << report.norm_after << "\n"
          << "max_overlap    = " << report.max_overlap << "\n";
    }

    // ---- CAPs ----------------------------------------------------------------
    perturbations::absorbing cap_lo(CAP_ETA * 1.0_Ha, -CAP_MID_FRAC, CAP_WIDTH_FRAC);
    perturbations::absorbing cap_hi(CAP_ETA * 1.0_Ha,  CAP_MID_FRAC, CAP_WIDTH_FRAC);
    auto pert_with_cap = perturbations::sum(cap_lo, cap_hi);

    // phi_plus = 0: real C ions are not a background charge group
    inq::basis::field<inq::basis::real_space, double> ie_phiplus(electrons.density().basis());
    ie_phiplus.fill(0.0);
    const double E_BB = 0.0;

    // ---- writers (segment-suffixed on resume) --------------------------------
    inqkit::io::RealField3DLayout vti_layout{.field_name="density",.include_meta=false,
        .emit_raw=false,.emit_vti=true,.vti_format=inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter total_wr(VTI+"/density_total", vti_layout, {.overwrite=(START==0)});
    inqkit::io::RealField3DWriter wp_wr   (VTI+"/density_wp",    vti_layout, {.overwrite=(START==0)});
    inqkit::io::ComplexField3DWriter wf_wr(VTI+"/wavefunction_wp",
        {.field_name="wavefunction",.include_meta=false,.emit_raw=false,
         .emit_vti=true,.vti_format=inqkit::io::VTIWriteOptions::Format::binary},
        {.overwrite=(START==0)});

    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.energy_external = sel.energy_nonlocal = sel.energy_ion = true;
    sel.energy_ion_kinetic = true;   // Ehrenfest ionic KE
    sel.energy_exact_exchange = sel.energy_nvxc = sel.energy_eigenvalues = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs(OBS + "/observables" + SEG + ".csv", sel);
    obs.write_header();

    // wp_momentum_stats: ekin = <p²>/2m  (columns pz,px,py give <p>; drift KE = (pz²+px²+py²)/2)
    inqkit::observables::WPMomentumStats  wp_mom(OBS+"/wp_momentum_stats"+SEG+".csv",  wp_idx, {.write_every=1});
    inqkit::observables::WPRealSpaceStats wp_pos(OBS+"/wp_real_space_stats"+SEG+".csv", wp_idx, {.write_every=1});

    inqkit::observables::DensityDelta density_delta(
        VTI+"/density_delta", VTI+"/density_delta_coarse",
        {.emit_raw_vti=true,.emit_coarse_vti=true,.compute_l2=true,.coarse_bin_bohr=3.0,
         .emit_every=(SAVE_EVERY>0?SAVE_EVERY:1)});

    std::ofstream ix;
    if (electrons.root()) {
        ix.open(OBS + "/interactions" + SEG + ".csv");
        ix << std::setprecision(12)
           << "step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,"
              "e_hartree_check,e_external_check,norm_wp,norm_total\n";
    }
    std::ofstream ionst;
    if (electrons.root()) {
        ionst.open(OBS + "/ions_track" + SEG + ".csv");
        ionst << "step,time_au";
        for (int i = 0; i < N_ATOMS; ++i) ionst << ",z" << i;
        ionst << "\n";
    }

    // ---- t=0 gates ----------------------------------------------------------
    if (START == 0) {
        auto m0 = wp_mom.compute(electrons);
        auto r0 = wp_pos.compute(electrons);
        int fails = 0;
        auto gate_abs = [&](const char* nm, double got, double want, double tol){
            const bool ok = std::abs(got-want) <= tol;
            std::cout << (ok?"  [PASS] ":"  [FAIL] ") << nm << ": " << got
                      << "  (expect " << want << " +/- " << tol << ")\n";
            if (!ok) ++fails;
        };
        auto gate_rel = [&](const char* nm, double got, double want, double pc){
            const double r = (want!=0.0)?100.0*(got-want)/std::abs(want):0.0;
            const bool ok = std::abs(r) <= pc;
            std::cout << (ok?"  [PASS] ":"  [FAIL] ") << nm << ": " << got
                      << "  (expect " << want << ", dev " << r << " %)\n";
            if (!ok) ++fails;
        };
        std::cout << "\n  --- t=0 analytic gates ---\n";
        gate_abs("norm (real space)",    r0.N,   1.0,  0.02);
        gate_rel("<p_z> = k0",          m0.pz,  K0,   2.0);
        gate_rel("sigma_pz^2=1/(2s^2)", m0.sz2, sigma_p2, 10.0);
        gate_rel("T1=(k0^2+3sp2)/2",    m0.ekin, 0.5*(K0*K0+3.0*sigma_p2), 3.0);
        gate_abs("centroid z",          r0.zc,  LAUNCH_Z, 0.05);
        gate_abs("ortho removed < 3%",  100.0*report.removed_weight, 0.0,
                 env_d("LJ_ORTHO_TOL_PC", 3.0));
        {
            const double sp   = 1.0/(std::sqrt(2.0)*SIGMA_WP);
            const double knyq = std::sqrt(2.0) * M_PI / DX;  // from spacing
            const double tail = 0.5*std::erfc((knyq-K0)/sp/std::sqrt(2.0));
            std::cout << "  [info] ALIASING: k_Nyq=" << knyq << " sp=" << sp
                      << " tail=" << 100.0*tail << " %\n";
            if (tail > 0.02) std::cout << "  [WARN] > 2% aliased.\n";
        }
        if (fails > 0) { std::cerr << "\nFATAL: " << fails << " gate(s) failed.\n"; return 4; }
        std::cout << "  all t=0 gates PASSED\n\n";

        total_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);
        wp_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0);
    }

    // ---- real-time (Ehrenfest ions) ------------------------------------------
    inqkit::RealTimeSession rt(ions, electrons, 1);
    rt.add([&](inqkit::StepContext const& ctx) {
        const int step = ctx.step;
        auto n_tot = inqkit::fields::density::total(*ctx.electrons);
        const double l2 = density_delta.snapshot(n_tot, ctx.time_au, step);
        inqkit::StepContext out = ctx; out.density_l2 = l2;
        obs.append(out);

        auto n_wp_f = inqkit::jellium::orbital_density_field(*ctx.electrons, wp_idx);
        auto ct = inqkit::jellium::compute_coulomb_wp(ctx.electrons->density(), n_wp_f, ie_phiplus);
        if (ctx.electrons->root())
            ix << step << ',' << ctx.time_au << ',' << ct.e_ss << ',' << ct.e_pp << ','
               << ct.e_ps << ',' << ct.e_sb << ',' << ct.e_pb << ',' << E_BB << ','
               << ct.e_hartree_check << ',' << ct.e_external_check << ','
               << ct.norm_wp << ',' << ct.norm_total << '\n';

        if (SAVE_EVERY > 0 && step % SAVE_EVERY == 0) {
            total_wr.write(n_tot, ctx.time_au, step);
            wp_wr.write(inqkit::fields::density::orbital(*ctx.electrons, wp_idx), ctx.time_au, step);
            if (ctx.electrons->root()) {
                ionst << step << ',' << ctx.time_au;
                for (int i = 0; i < N_ATOMS; ++i) ionst << ',' << ions.positions()[i][2];
                ionst << '\n';
            }
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
        std::ofstream ic(OUT + "/ckpt_ions_step" + tag6(last) + ".csv");
        ic << std::setprecision(16) << "i,x,y,z,vx,vy,vz\n";
        for (int i = 0; i < N_ATOMS; ++i) {
            auto p = ions.positions()[i]; auto v = ions.velocities()[i];
            ic << i << ',' << p[0] << ',' << p[1] << ',' << p[2] << ','
               << v[0] << ',' << v[1] << ',' << v[2] << '\n';
        }
    };

    real_time::propagate(
        ions, electrons,
        [&](auto const& data) {
            rt.step(data);
            wp_mom.maybe_accumulate(data);
            wp_pos.maybe_accumulate(data);
            if (data.iter() > 0 && data.iter() % CKPT_EVERY == 0 && data.iter() < N_STEPS) {
                electrons.save(CKPT); write_rt_state(data.iter());
                std::cout << "  [ckpt] step " << data.iter() << "\n" << std::flush;
            }
        },
        options::theory{}.lda(),
        options::real_time{}.num_steps(N_STEPS).dt(DT * 1.0_atomictime).ehrenfest(),
        pert_with_cap, START);

    electrons.save(CKPT);
    write_rt_state(N_STEPS);

    const double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_wall0).count();

    if (electrons.root()) {
        ix.close(); ionst.close();
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(14)
          << "run = graphene/lovelace_test/wp/" << env_s("LJ_OUT","") << "\n"
          << "run_type = wavepacket projectile, graphene " << VARIANT << " TDDFT (ALDA), Ehrenfest ions\n"
          << "date_finished = " << iso_now() << "\nwall_time_s = " << wall << "\n"
          << "variant = " << VARIANT << "\nn_atoms = " << N_ATOMS << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << LZ << "\n"
          << "periodicity = 2\ndx_bohr = " << DX << "\n"
          << "n_electrons = " << N_ELEC << "\nextra_states = " << EXTRA << "\n"
          << "ion_dynamics = ehrenfest\n"
          << "wp_state_index = " << wp_idx << "\n"
          << "wp_sigma_bohr = " << SIGMA_WP << "\nwp_k0_bohr_inv = " << K0 << "\n"
          << "wp_drift_energy_ev = " << 0.5*K0*K0*HA << "\n"
          << "impact_xy_bohr = " << IMPACT_X << " " << IMPACT_Y << "\nlaunch_z = " << LAUNCH_Z << "\n"
          << "cap_eta_ha = " << CAP_ETA << "  cap_width_bohr = " << CAP_L << "\n"
          << "start_step = " << START << "  rt_num_steps = " << N_STEPS
          << "  dt_au = " << DT << "  total_time_au = " << (DT*N_STEPS) << "\n"
          << "save_every = " << SAVE_EVERY << "  ckpt_every = " << CKPT_EVERY << "\n"
          << "gs_dir = " << GS_DIR << "\nrun_completed = true\n";
    }
    std::cout << "\nDone. Wall " << wall << " s.\n";
    return 0;
}
