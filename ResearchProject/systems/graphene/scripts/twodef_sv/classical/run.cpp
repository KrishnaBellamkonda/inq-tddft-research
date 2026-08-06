// ============================================================================
// systems/graphene/scripts/twodef_sv/classical/run.cpp   (v2 — DIRECT potential)
//
// CLASSICAL half of the TWO-DEFINITIONS stopping campaign twins (graphene).
// Plan: docs/plans/real-material-stopping-comparison.md  (Phase 2)
//
// v2 REWORK (2026-08-05). v1 used an INQ ion for the projectile (cap_cl idiom)
// and FAILED: INQ's ion-ion term treats the projectile as a POINT -1 charge
// (the erf smoothing lives only in the electron-side UPF), and the Bernal
// impact site is atom-atop in layer 2 -> 1/r^2 slingshot (E100: vz 2.71->95,
// z->744; results quarantined in results/invalid_pointion_*). This v2 uses the
// corpus-validated DIRECT design (sigma56_sv/classical): the projectile is a
// moving external erf/r potential (moving_gaussian_projectile_potential) with
// its own velocity-Verlet (inqkit::dynamics::Projectile), NEVER an INQ ion.
//
// Forces on the projectile (both erf-smoothed with sigma_pot = sigma_WP/sqrt2):
//   F_el   = projectile_force_direct_z(n_e)      analytic HF force, sign
//            validated by energy conservation in the jellium campaigns;
//   F_core = sum over C cores (+4 each, LIVE Ehrenfest positions, xy minimum
//            image): F = (-q_P*Z_C) * g'(r) * r_hat, g(r) = erf(r/(sigma*sqrt2))/r.
//            g'(0) = 0 -> smooth at contact, NO singularity by construction.
// KNOWN LIMITATION (recorded): the carbons do not feel the projectile's
// reaction force (perturbation potentials act on electrons only) — the NUCLEAR
// stopping channel is not back-reacted. Kinematic bound: 4mM/(m+M)^2 ~ 1.8e-4
// of E, i.e. <= 0.05 eV at 300 eV. Electronic stopping is unaffected.
//
// The carbon lattice itself stays INQ-Ehrenfest (.ehrenfest()), as in the WP
// twin. Projectile is stopped BEFORE the CAP inner edge by run sizing (a CAP
// absorbs wavefunctions, not classical charges); since we own the coordinate,
// the run simply ends there — no wrap, no out-of-box evaluation, ledger clean.
//
// interactions.csv: compute_coulomb_direct with zero nplus/phiplus (real ions
// are not a background charge group). Classical closure: e_ss == energy_hartree
// (gated in analysis). Schema identical to the WP half.
//
// Env: GR_VARIANT(bi|mono) LJ_K0(2.71) LJ_SIGMA(2.0 = sigma_WP) LJ_LAUNCH_Z(-12)
//      LJ_IMPACT_X(2.3244) LJ_IMPACT_Y(1.3420) LJ_DT(0.025) LJ_N_STEPS(auto)
//      LJ_CAP_ETA(REQUIRED<0) LJ_CAP_L(REQUIRED) LJ_SAVE_EVERY(auto ~50 frames)
//      LJ_CKPT_EVERY(0=auto N/3) LJ_RESUME(0) LJ_OUT(REQUIRED) LJ_GS_DIR(REQUIRED)
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/interaction_energies.hpp>       // compute_coulomb_direct
#include <inqkit/jellium/gaussian_potential.hpp>
#include <inqkit/dynamics/projectile.hpp>
#include <inqkit/dynamics/projectile_force.hpp>
#include <inqkit/dynamics/moving_gaussian_projectile_potential.hpp>

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
    char b[64]; std::strftime(b, sizeof(b), "%Y-%m-%dT%H:%M:%S", &tm); return std::string(b);
}
// minimum-image displacement in a periodic direction
static double min_image(double d, double L){
    while (d >  L/2) d -= L;
    while (d < -L/2) d += L;
    return d;
}

int main() {
    auto t_wall0 = std::chrono::steady_clock::now();
    const double HA = 27.211386245988;
    const double Q_P = -1.0;     // projectile charge (electron)
    const double Z_C =  4.0;     // ONCV-C valence core charge

    const std::string VARIANT = env_s("GR_VARIANT", "bi");
    const bool  BI      = (VARIANT == "bi");
    const int   N_ATOMS = BI ? Cfg::N_C_BI : Cfg::N_C_MONO;
    const int   N_ELEC  = BI ? Cfg::N_ELEC_BI : Cfg::N_ELEC_MONO;
    const int   EXTRA   = env_i("GR_EXTRA", BI ? Cfg::EXTRA_STATES_BI : Cfg::EXTRA_STATES_MONO);
    const std::string GEOM = BI ? Cfg::GEOM_BI : Cfg::GEOM_MONO;

    const double CUTOFF   = env_d("GR_CUTOFF_HA", Cfg::CUTOFF_HA);
    const double LZ       = env_d("GR_LZ_BOHR", Cfg::LZ_BOHR);
    const double SIGMA_WP = env_d("LJ_SIGMA", 2.0);
    const double SIGMA_POT = SIGMA_WP / std::sqrt(2.0);
    const double K0       = env_d("LJ_K0", 2.71);
    const double LAUNCH_Z = env_d("LJ_LAUNCH_Z", -12.0);
    const double IMPACT_X = env_d("LJ_IMPACT_X", 2.3244);
    const double IMPACT_Y = env_d("LJ_IMPACT_Y", 1.3420);
    const double DT       = env_d("LJ_DT", 0.025);
    const double CAP_ETA  = env_d("LJ_CAP_ETA", 0.0);
    const double CAP_L    = env_d("LJ_CAP_L", 0.0);
    const bool   RESUME   = env_i("LJ_RESUME", 0) != 0;
    const std::string OUT    = "results/" + env_s("LJ_OUT", "");
    const std::string GS_DIR = env_s("LJ_GS_DIR", "");

    const double z_cap_in = LZ/2.0 - CAP_L;
    const int N_STEPS_AUTO = int(std::ceil((std::abs(LAUNCH_Z) + z_cap_in) / K0 / DT));
    const int N_STEPS      = env_i("LJ_N_STEPS", N_STEPS_AUTO);
    const int SAVE_EVERY   = env_i("LJ_SAVE_EVERY", std::max(1, N_STEPS / 50));
    int CKPT_EVERY = env_i("LJ_CKPT_EVERY", 0);
    if (CKPT_EVERY <= 0) CKPT_EVERY = std::max(1, N_STEPS / 3);

    if (env_s("LJ_OUT", "").empty()) { std::cerr << "FATAL: LJ_OUT unset\n"; return 2; }
    if (GS_DIR.empty() || !std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: LJ_GS_DIR missing/unset\n"; return 2;
    }
    if (CAP_ETA >= 0.0 || CAP_L <= 0.0) {
        std::cerr << "FATAL: CAP required (twin parity with the WP half).\n"; return 2;
    }

    const double CAP_WIDTH_FRAC = CAP_L / LZ;
    const double CAP_MID_FRAC   = 0.5 - CAP_WIDTH_FRAC / 2.0;

    const std::string CKPT = OUT + "/checkpoint", RT_STATE = OUT + "/rt_state.txt";
    int START = 0;
    double R0z = LAUNCH_Z, V0z = K0;
    if (RESUME) {
        START = (int)read_state_d(RT_STATE, "last_step", -1);
        R0z   = read_state_d(RT_STATE, "proj_z", LAUNCH_Z);
        V0z   = read_state_d(RT_STATE, "proj_vz", K0);
        if (START < 0) { std::cerr << "FATAL: LJ_RESUME=1, no " << RT_STATE << "\n"; return 2; }
        if (START >= N_STEPS) { std::cout << "Nothing to do.\n"; return 0; }
    }
    const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

    std::cout << std::setprecision(10)
              << "\n=== twodef_sv CLASSICAL v2 (direct erf/r) [" << VARIANT
              << "]  OUT=" << OUT << " ===\n"
              << "  cell     = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR << " x " << LZ
              << " Bohr, periodicity(2), cutoff=" << CUTOFF << " Ha\n"
              << "  ions     = " << N_ATOMS << " C (EHRENFEST); projectile = moving "
                 "external erf/r potential (NOT an ion)\n"
              << "  proj     = sigma_pot=" << SIGMA_POT << " (sigma_WP " << SIGMA_WP
              << "), q=" << Q_P << ", m=1, v0=" << K0 << "  E=" << 0.5*K0*K0*HA << " eV\n"
              << "  launch   = (" << IMPACT_X << "," << IMPACT_Y << "," << R0z << ")"
              << "  run ends before z=" << z_cap_in << "\n"
              << "  CAP      = eta=" << CAP_ETA << " Ha  W=" << CAP_L
              << " Bohr/face (secondary-electron flux)\n"
              << "  dt=" << DT << "  START=" << START << " -> N_STEPS=" << N_STEPS
              << (RESUME ? "  [RESUME]" : "") << "\n"
              << "  GS       = " << GS_DIR << "\n\n";

    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b,
                                            Cfg::LY_BOHR * 1.0_b,
                                            LZ * 1.0_b).periodicity(2);
    auto ions = systems::ions::parse(GEOM, cell);
    if (int(ions.size()) != N_ATOMS) {
        std::cerr << "FATAL: parsed " << ions.size() << " != " << N_ATOMS << "\n"; return 2;
    }
    // Plain graphene electrons — NO projectile charge in the cell, so the
    // construction matches the GS checkpoint exactly (192 e, 144 states for bi).
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .cutoff(CUTOFF * 1.0_Ha)
            .extra_states(EXTRA)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    const int n_states = electrons.states().num_states();

    const std::string OBS = OUT + "/raw/observables";
    const std::string VTI = OUT + "/raw/vti";
    for (auto const& d : {OBS, VTI + "/density_total", VTI + "/density_delta",
                          VTI + "/density_delta_coarse", VTI + "/density_gs_system"})
        std::filesystem::create_directories(d);

    electrons.load(RESUME ? CKPT : GS_DIR);
    if (!RESUME) {
        inqkit::io::RealField3DLayout lay{
            .field_name = "density", .include_meta = false, .emit_raw = false,
            .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
        inqkit::io::RealField3DWriter gs_wr(VTI + "/density_gs_system", lay, {.overwrite = true});
        gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system");
    }

    inqkit::dynamics::Projectile proj(1.0, Q_P,
        inqkit::detail::Vec3{IMPACT_X, IMPACT_Y, R0z},
        inqkit::detail::Vec3{0.0, 0.0, V0z});
    inqkit::dynamics::moving_gaussian_projectile_potential proj_pert(proj, SIGMA_POT);

    perturbations::absorbing cap_lo(CAP_ETA * 1.0_Ha, -CAP_MID_FRAC, CAP_WIDTH_FRAC);
    perturbations::absorbing cap_hi(CAP_ETA * 1.0_Ha,  CAP_MID_FRAC, CAP_WIDTH_FRAC);
    auto pert = perturbations::sum(proj_pert, perturbations::sum(cap_lo, cap_hi));

    // zero background fields (real ions are not a background charge group)
    inq::basis::field<inq::basis::real_space, double> zero_n(electrons.density().basis());
    inq::basis::field<inq::basis::real_space, double> zero_phi(electrons.density().basis());
    zero_n.fill(0.0); zero_phi.fill(0.0);
    const double E_BB = 0.0;
    auto basis = electrons.density().basis();

    // erf-smoothed core force on the projectile: F = (-Q_P*Z_C) g'(r) r_hat,
    // g(r) = erf(r/(sigma*sqrt2))/r; g'(0)=0 (smooth). xy minimum image; live
    // Ehrenfest carbon positions.
    auto core_force_z = [&](double px, double py, double pz) -> double {
        double Fz = 0.0;
        const double s = SIGMA_POT;
        for (int i = 0; i < N_ATOMS; ++i) {
            auto const& c = ions.positions()[i];
            const double dx = min_image(px - c[0], Cfg::LX_BOHR);
            const double dy = min_image(py - c[1], Cfg::LY_BOHR);
            const double dz = pz - c[2];
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double r  = std::sqrt(r2);
            if (r < 1e-8) continue;                     // g'(0)=0
            const double gp = std::sqrt(2.0/M_PI) * std::exp(-r2/(2.0*s*s)) / (s*r)
                              - std::erf(r/(s*std::sqrt(2.0))) / r2;
            Fz += (-Q_P * Z_C) * gp * (dz / r);
        }
        return Fz;
    };

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter total_wr(VTI + "/density_total", vti_layout, {.overwrite = (START == 0)});

    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.energy_external = sel.energy_nonlocal = sel.energy_ion = true;
    sel.energy_ion_kinetic = sel.energy_exact_exchange = true;
    sel.energy_nvxc = sel.energy_eigenvalues = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs(OBS + "/observables" + SEG + ".csv", sel);
    obs.write_header();

    inqkit::observables::DensityDelta density_delta(
        VTI + "/density_delta", VTI + "/density_delta_coarse",
        {.emit_raw_vti = true, .emit_coarse_vti = true,
         .compute_l2 = true, .coarse_bin_bohr = 3.0,
         .emit_every = (SAVE_EVERY > 0 ? SAVE_EVERY : 1)});

    std::ofstream ix, track, ionst;
    if (electrons.root()) {
        ix.open(OBS + "/interactions" + SEG + ".csv");
        ix << std::setprecision(12)
           << "step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,"
              "e_hartree_check,e_external_check,norm_wp,norm_total\n";
        track.open(OBS + "/electron_track" + SEG + ".csv");
        track << std::setprecision(12)
              << "step,time_au,x,y,z,vz,ke_ion_ha,fz_el,fz_core\n";
        ionst.open(OBS + "/ions_track" + SEG + ".csv");
        ionst << "step,time_au";
        for (int i = 0; i < N_ATOMS; ++i) ionst << ",z" << i;
        ionst << "\n";
    }

    if (START == 0)
        total_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);

    inqkit::RealTimeSession rt(ions, electrons, 1);
    rt.add([&](inqkit::StepContext const& ctx) {
        const int step = ctx.step;

        auto n_tot = inqkit::fields::density::total(*ctx.electrons);
        const double l2 = density_delta.snapshot(n_tot, ctx.time_au, step);
        inqkit::StepContext out = ctx; out.density_l2 = l2;
        obs.append(out);

        auto Rn = proj.R();
        inq::vector3<double> center{Rn.x, Rn.y, Rn.z};
        const double Fel   = inqkit::dynamics::projectile_force_direct_z(
                                 ctx.electrons->density(), cell, center, SIGMA_POT);
        const double Fcore = core_force_z(Rn.x, Rn.y, Rn.z);

        auto vproj = inqkit::jellium::gaussian_potential(basis, center, SIGMA_POT);
        auto ct = inqkit::jellium::compute_coulomb_direct(
            ctx.electrons->density(), vproj, zero_n, zero_phi, SIGMA_POT);

        proj.advance(inqkit::detail::Vec3{0.0, 0.0, Fel + Fcore}, DT);

        if (ctx.electrons->root()) {
            ix << step << ',' << ctx.time_au << ',' << ct.e_ss << ',' << ct.e_pp << ','
               << ct.e_ps << ',' << ct.e_sb << ',' << ct.e_pb << ',' << E_BB << ','
               << ct.e_ss << ',' << 0.0 << ','
               << ct.norm_p << ',' << ct.norm_slab << '\n';
            track << step << ',' << ctx.time_au << ',' << Rn.x << ',' << Rn.y << ','
                  << Rn.z << ',' << proj.V().z << ',' << proj.ke() << ','
                  << Fel << ',' << Fcore << '\n';
        }

        if (SAVE_EVERY > 0 && step % SAVE_EVERY == 0) {
            total_wr.write(n_tot, ctx.time_au, step);
            if (ctx.electrons->root()) {
                ionst << step << ',' << ctx.time_au;
                for (int i = 0; i < N_ATOMS; ++i) ionst << ',' << ions.positions()[i][2];
                ionst << '\n';
            }
        }
    });

    auto write_rt_state = [&](int last){
        if (!electrons.root()) return;
        std::ofstream st(RT_STATE);
        st << std::setprecision(16)
           << "last_step=" << last << "\ntime_au=" << (last*DT) << "\ndt=" << DT
           << "\nproj_z=" << proj.R().z << "\nproj_vz=" << proj.V().z
           << "\nproj_x=" << proj.R().x << "\nproj_y=" << proj.R().y
           << "\nproj_mass=1\nproj_charge=" << Q_P << "\nsigma_wp=" << SIGMA_WP << "\n";
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
            if (data.iter() > 0 && data.iter() % CKPT_EVERY == 0 && data.iter() < N_STEPS) {
                electrons.save(CKPT);
                write_rt_state(data.iter());
                std::cout << "  [ckpt] step " << data.iter() << "\n" << std::flush;
            }
        },
        options::theory{}.lda(),
        options::real_time{}.num_steps(N_STEPS).dt(DT * 1.0_atomictime).ehrenfest(),
        pert, START);

    electrons.save(CKPT);
    write_rt_state(N_STEPS);

    const double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_wall0).count();

    if (electrons.root()) {
        ix.close(); track.close(); ionst.close();
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(14)
          << "run = graphene/twodef_sv/classical/" << env_s("LJ_OUT","") << "\n"
          << "run_type = classical electron projectile (DIRECT erf/r potential, "
             "velocity-Verlet), graphene " << VARIANT << " TDDFT (ALDA), Ehrenfest lattice\n"
          << "representation = direct_potential (v2; v1 point-ion INVALID, quarantined)\n"
          << "campaign = two-definitions stopping (real material)\n"
          << "plan = docs/plans/real-material-stopping-comparison.md\n"
          << "engine = inq-study\nxc = LDA (ALDA)\npseudos = ONCV-C (lattice only)\n"
          << "date_finished = " << iso_now() << "\nwall_time_s = " << wall << "\n"
          << "variant = " << VARIANT << "\nn_atoms = " << N_ATOMS << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << LZ << "\n"
          << "periodicity = 2\ncutoff_ha = " << CUTOFF << "\n"
          << "n_electrons = " << N_ELEC << "\nn_states = " << n_states << "\n"
          << "ion_dynamics = ehrenfest lattice; projectile external (no back-reaction "
             "on cores — nuclear channel <= 1.8e-4 of E, recorded)\n"
          << "proj_sigma_wp_label = " << SIGMA_WP
          << "\nproj_sigma_pot = " << SIGMA_POT << " (= sigma_WP/sqrt2)\n"
          << "proj_mass = 1 (m_e)\nproj_charge = " << Q_P << "\n"
          << "proj_v0 = " << K0 << "  proj_energy_ev = " << 0.5*K0*K0*HA << "\n"
          << "impact_xy_bohr = " << IMPACT_X << " " << IMPACT_Y << "\nlaunch_z = " << LAUNCH_Z << "\n"
          << "proj_final_z = " << proj.R().z << "  proj_final_vz = " << proj.V().z << "\n"
          << "cap_eta_ha = " << CAP_ETA << "  cap_width_bohr = " << CAP_L << "\n"
          << "interactions_note = compute_coulomb_direct, zero background; closure "
             "e_ss == energy_hartree (analysis gate)\n"
          << "start_step = " << START << "  rt_num_steps = " << N_STEPS
          << "  dt_au = " << DT << "  total_time_au = " << (DT*N_STEPS) << "\n"
          << "save_every = " << SAVE_EVERY << "  track_every = 1  ckpt_every = " << CKPT_EVERY << "\n"
          << "gs_dir = " << GS_DIR << "\nrun_completed = true\n";
    }
    std::cout << "\nDone. proj final z=" << proj.R().z << " vz=" << proj.V().z
              << ". Wall " << wall << " s.\n";
    return 0;
}
