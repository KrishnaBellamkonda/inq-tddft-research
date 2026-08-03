// ============================================================================
// cylindrical_jellium / scripts/proximity_ladder / classical / run.cpp
//
// CLASSICAL HALF of the twin at ONE RUNG of the proximity ladder: a rigid
// Gaussian CHARGE of std sigma_pot = sigma_WP/sqrt2 = 2.8284 Bohr, charge -1,
// mass m_e, launched on-axis at z = -28 with v = 1.917 a.u. (50 eV) down a
// periodic r_s = 3 jellium tube whose bore radius R_in is the rung parameter,
// advanced by velocity-Verlet Ehrenfest from its own Hellmann-Feynman force.
//
// Plan: docs/plans/cylindrical-proximity-ladder.md
//
// ---------------------------------------------------------------------------
// BUILD ONCE, RUN PER RUNG. CJ_RUNG selects the geometry at RUNTIME from
// shared/configs/proximity_ladder_rs3.hpp; there is no default, because a
// silently-defaulted rung yields a plausible run at the wrong geometry:
//
//   rung   R_in  R_out   N_e  states   shape
//   r10    10.0  14.000  160    104    annulus (DONE — channeling_twin)
//   r08     8.0  14.000  220    143    annulus
//   r06     6.0  13.986  266    173    annulus
//   r04     4.0  14.000  300    195    annulus
//   r00     0.0  13.986  326    212    CYLINDER (filled)
//
// r_s is held at 3.000000 on every rung (R_out is SOLVED, not transcribed).
//
// THE CLASSICAL TWIN DOES NOT SPREAD, AND THAT MATTERS MORE AS R_in SHRINKS.
// The wavepacket's density std grows 2.83 -> 6.01 Bohr over the run while this
// rigid charge keeps its width forever, so the two halves sample DIFFERENT
// amounts of bath at late time — by construction, not by error. It is also why
// the WP-minus-classical gap at the strong rungs is not purely "dispersion":
// this projectile contributes nothing to E_xc (it is an external potential, not
// a density), whereas the WP adds ~32 % to the on-axis density and LDA xc is
// non-linear in n. Measure that asymmetry; do not attribute it to dispersion.
//
// ---------------------------------------------------------------------------
// Twin: ../wp/run.cpp — IDENTICAL in every physical parameter (same GS, cell,
// grid, N, sigma, launch point, dt, n_steps); the ONLY difference is that the
// projectile here is an external Gaussian charge and there it is an occupied KS
// orbital. That is the whole point: a difference in S between the two is a
// quantum effect and nothing else.
//
// ---------------------------------------------------------------------------
// WHY A MOVING GAUSSIAN CHARGE AND NOT A GHOST UPF
// The ghost-UPF classical projectile carries an r_cut, and outside r_cut the
// pseudopotential is the bare -1/r of a POINT charge rather than the erf/r of the
// intended Gaussian. That mismatch inflates S. A perturbation built as
// v_proj = +poisson(n_proj) with n_proj an exactly normalised Gaussian has no
// r_cut and no aliasing tail, and it sources EXACTLY the potential the
// wavepacket's own charge distribution sources when sigma_pot = sigma_WP/sqrt2 —
// which is the condition that makes the twins comparable at all.
// (.claude/rules/sigma-wp-convention.md; reference_projectile_charge_sheet_inflates_stopping.)
//
// ---------------------------------------------------------------------------
// MINIMUM IMAGE IS MANDATORY HERE, NOT OPTIONAL
// The launch point is 2 Bohr from the -z face = 0.71 sigma_pot. A plain
// Cartesian Gaussian would keep only Phi(2/2.83) = 76 % of its charge at t = 0,
// and — because the +delta and -delta finite-difference blobs are clipped by
// DIFFERENT amounts — the loss does not cancel out of the force. The wavepacket
// twin has no such problem: a KS orbital lives on a plain 3-D FFT basis and
// wraps exactly. So both the perturbation AND the force use the minimum-image
// kernel (inqkit::jellium::gaussian_density_minimum_image, and the
// minimum_image=true option on projectile_force_axis added for this study).
//
// ---------------------------------------------------------------------------
// FULL 3-D EHRENFEST, ON PURPOSE
// The force is computed in all three Cartesian directions and the projectile is
// free to leave the axis. By the 4-fold symmetry of the tube on this grid the
// transverse force vanishes at r_perp = 0, so a correct run keeps x = y = 0 to
// round-off — and that is then a MEASURED statement about channeling stability
// rather than a constraint imposed by only integrating z. F_x is written every
// step for exactly this check.
//
// ---------------------------------------------------------------------------
// TRAJECTORY. 1500 steps x dt 0.02 = 30 a.u. covers 57.5 Bohr, so the projectile
// runs z = -28 -> +29.5: ONE traversal of the 60-Bohr periodic tube, no wrap.
// The wrap machinery (proj_z_unwrapped, n_wraps, minimum image) is nevertheless
// wired so the run can be EXTENDED past one traversal with LJ_RESUME=1
// (.claude/rules/final-timestep-checkpoint.md).
//
// EXPECTED BEHAVIOUR — NOT A BUG. A mass-1 electron at v = 1.917 carries 50 eV
// and DECELERATES (.claude/rules/light-projectile-stopping.md). S is extracted as
// the INITIAL DRAG over the early near-constant-velocity window (v >= 0.85 v0),
// never as a full-run regression.
//
// OUTPUTS (raw/observables/)
//   observables.csv     full INQ energy decomposition, every step
//   electron_track.csv  step,time_au,x,y,z,vx,vy,vz,ke_ion_ha  <- the CANONICAL
//                       classical schema: it is what ks_stopping.load_classical_run
//                       reads and what the run-notebook builder uses to detect a
//                       classical run. Do not rename its columns.
//   projectile.csv      the same trajectory plus U_proj_bg, the unwrapped path,
//                       the wrap count and the three force components
//   interactions.csv    pairwise P/S/B Coulomb ledger (.claude/rules/
//                       decomposed-interaction-energies.md), every step
//   electron_number.csv bath norm conservation
//
// Env: CH_GS_DIR(REQUIRED) CH_OUT(classical) CH_N(160) CH_V0 CH_SIGMA CH_LAUNCH_Z
//      CH_N_STEPS CH_DT CH_SPACING CH_SAVE_EVERY CH_CKPT_EVERY(0=auto N/3)
//      CH_MAX_CKPT(3) CH_RESUME(0) CH_CONST_V(0) CH_FD_DELTA
// Build against INQ_SOURCE=inq-study; runtime shares from inq/install/share.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>
#include <inqkit/jellium/interaction_energies.hpp>
#include <inqkit/jellium/analytics.hpp>
#include <inqkit/dynamics/projectile.hpp>
#include <inqkit/dynamics/projectile_force.hpp>
#include <inqkit/dynamics/moving_gaussian_projectile_perturbation.hpp>

#include "../../../shared/configs/proximity_ladder_rs3.hpp"

#include <algorithm>
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
    const double V0       = env_d("CH_V0",        Cfg::PROJ_V0);
    const double LAUNCH_Z = env_d("CH_LAUNCH_Z",  Cfg::LAUNCH_Z_BOHR);
    const double DT       = env_d("CH_DT",        Cfg::DT_AU);
    const int    N_STEPS  = env_i("CH_N_STEPS",   Cfg::N_STEPS);
    const int    SAVE_EVERY = env_i("CH_SAVE_EVERY", Cfg::SAVE_EVERY);
    const double FD_DELTA = env_d("CH_FD_DELTA",  Cfg::FD_DELTA);
    const double MASS     = Cfg::PROJ_MASS;
    const bool   CONST_V  = env_i("CH_CONST_V", 0) != 0;
    const bool   RESUME   = env_i("CH_RESUME", 0) != 0;
    const std::string OUT    = "results/" + env_s("CH_OUT", "classical");
    const std::string GS_DIR = env_s("CH_GS_DIR", "");

    int CKPT_EVERY = env_i("CH_CKPT_EVERY", 0);
    if(CKPT_EVERY <= 0) CKPT_EVERY = std::max(1, N_STEPS / 3);

    if(GS_DIR.empty() || !std::filesystem::exists(GS_DIR)){
        std::cerr << "FATAL: CH_GS_DIR missing or unset: '" << GS_DIR << "'\n"; return 2;
    }

    const double SIGMA_POT = SIGMA_WP / std::sqrt(2.0);
    const double V_ANN     = rung->v_jellium();
    const double N0        = double(N_ELEC) / V_ANN;
    const double RS        = inqkit::jellium::rs_from_n0(N0);
    const double OMEGA_P   = std::sqrt(4.0*M_PI*N0);
    const double VF        = inqkit::jellium::k_fermi_n0(N0);

    // ---- resume state -------------------------------------------------------
    const std::string CKPT = OUT + "/checkpoint", RT_STATE = OUT + "/rt_state.txt";
    int START = 0, NWRAP0 = 0;
    double R0x = 0.0, R0y = 0.0, R0z = LAUNCH_Z;
    double V0x = 0.0, V0y = 0.0, V0z = V0 / MASS;
    double Z_UNWRAP0 = LAUNCH_Z;
    if(RESUME){
        START = (int)read_state_d(RT_STATE, "last_step", -1);
        if(START < 0){ std::cerr << "FATAL: CH_RESUME=1 but no readable " << RT_STATE << "\n"; return 2; }
        if(START >= N_STEPS){
            std::cout << "Already at/after target (" << START << " >= " << N_STEPS
                      << "); nothing to do.\n"; return 0; }
        R0x = read_state_d(RT_STATE, "proj_x",  0.0);
        R0y = read_state_d(RT_STATE, "proj_y",  0.0);
        R0z = read_state_d(RT_STATE, "proj_z",  LAUNCH_Z);
        V0x = read_state_d(RT_STATE, "proj_vx", 0.0);
        V0y = read_state_d(RT_STATE, "proj_vy", 0.0);
        V0z = read_state_d(RT_STATE, "proj_vz", V0 / MASS);
        Z_UNWRAP0 = read_state_d(RT_STATE, "proj_z_unwrapped", R0z);
        NWRAP0    = (int)read_state_d(RT_STATE, "n_wraps", 0);
    }
    const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

    std::cout << std::setprecision(12)
              << "\n=== proximity_ladder/classical  rung=" << rung->label << "  OUT=" << OUT << " ===\n"
              << "  tube      = R_in " << R_IN << "  R_out " << R_OUT
              << "  L_z " << Cfg::LZ_BOHR << "  (fully periodic)\n"
              << "  cell      = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR << " x "
              << Cfg::LZ_BOHR << "  dx = " << SPACING << "\n"
              << "  bath      = " << N_ELEC << " e   n0 = " << N0 << "   r_s = " << RS
              << "   omega_p = " << OMEGA_P*HA << " eV   v_F = " << VF << "\n"
              << "  projectile= Gaussian CHARGE  sigma_WP " << SIGMA_WP
              << " (sigma_pot " << SIGMA_POT << ")  q -1  m " << MASS << "\n"
              << "  launch    = (0,0," << LAUNCH_Z << ")  v0 = " << V0
              << "  (E = " << 0.5*MASS*V0*V0*HA << " eV, v/v_F = " << V0/VF << ")\n"
              << "  drive     = " << (CONST_V ? "const_velocity" : "ehrenfest (3-D force)") << "\n"
              << "  dt = " << DT << "  START = " << START << " -> N_STEPS = " << N_STEPS
              << "  t_total = " << DT*N_STEPS << " a.u." << (RESUME ? "  [RESUME]" : "") << "\n"
              << "  cadence   : density/" << SAVE_EVERY << "  scalars/1  ckpt/" << CKPT_EVERY << "\n"
              << "  GS        = " << GS_DIR << "\n\n";

    // ---- system -------------------------------------------------------------
    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b,
                                            Cfg::LY_BOHR * 1.0_b,
                                            Cfg::LZ_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);                    // NO projectile ion
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(N_ELEC)
            .extra_states(EXTRA_ST)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    electrons.load(RESUME ? CKPT : GS_DIR);
    const int n_states = electrons.states().num_states();

    // ---- background + live projectile + moving perturbation -----------------
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

    inqkit::dynamics::Projectile proj(MASS, Cfg::PROJ_CHARGE,
        inqkit::detail::Vec3{R0x, R0y, R0z},
        inqkit::detail::Vec3{V0x, V0y, V0z});
    // minimum_image = true: see the header note. Same flag on the force below.
    inqkit::dynamics::moving_gaussian_projectile_perturbation proj_pert(proj, SIGMA_POT, true);
    auto pert = perturbations::sum(bg_pert, proj_pert);

    // constant background fields for the force and the pairwise ledger
    auto basis   = electrons.density().basis();
    auto nplus   = bg_pert.background_density(basis);
    auto phiplus = solvers::poisson::solve(nplus);
    const double E_BB = inqkit::jellium::background_self_energy(nplus, phiplus);

    // ---- output skeleton ----------------------------------------------------
    const std::string OBS = OUT + "/raw/observables";
    const std::string VTI = OUT + "/raw/vti";
    for(auto const& d : {OBS, VTI + "/density_total", VTI + "/density_delta",
                         VTI + "/density_delta_coarse", VTI + "/density_gs_system"})
        std::filesystem::create_directories(d);

    if(START == 0){
        inqkit::io::RealField3DWriter gs_wr(VTI + "/density_gs_system",
            {.field_name = "density", .include_meta = false, .emit_raw = false,
             .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary},
            {.overwrite = true});
        gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system");
    }

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter total_wr(VTI + "/density_total", vti_layout,
                                           {.overwrite = (START == 0)});

    // FULL energy decomposition — the twin contract requires it in BOTH halves.
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

    std::ofstream trk, pj, ix, nlog;
    if(electrons.root()){
        // CANONICAL classical trajectory schema — ks_stopping.load_classical_run
        // and run_notebook_builder.detect_type both key off this file and these
        // column names.
        trk.open(OBS + "/electron_track" + SEG + ".csv");
        trk << std::setprecision(16) << "step,time_au,x,y,z,vx,vy,vz,ke_ion_ha\n";

        pj.open(OBS + "/projectile" + SEG + ".csv");
        pj << std::setprecision(12)
           << "step,time_au,proj_z,proj_vz,energy_proj_ke,energy_proj_bg_ideal,"
              "proj_z_unwrapped,n_wraps,proj_x,proj_y,force_x,force_y,force_z\n";

        // Pairwise Coulomb ledger, Hartree. Same 12 columns as the WP twin so one
        // loader reads both. Closure (classical): e_hartree_check == INQ
        // energy_hartree and e_external_check == INQ energy_external.
        ix.open(OBS + "/interactions" + SEG + ".csv");
        ix << std::setprecision(12)
           << "step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,"
              "e_hartree_check,e_external_check,norm_proj,norm_electrons\n";

        nlog.open(OBS + "/electron_number" + SEG + ".csv");
        nlog << std::setprecision(12) << "step,time_au,N_total\n";
    }

    // continuous-path bookkeeping across wraps
    double z_unwrapped = Z_UNWRAP0, z_prev_wrapped = R0z;
    int    n_wraps     = NWRAP0;

    // energy-conservation tracking: with a rigid classical charge the conserved
    // quantity is E_electronic + KE_proj + U_proj_bg, not energy_total alone.
    double e_cons_first = 0.0, e_cons_last = 0.0; bool e_cons_seen = false;

    inqkit::RealTimeSession rt(ions, electrons, 1);     // callback EVERY step
    rt.add([&](inqkit::StepContext const& ctx){
        const int step = ctx.step;

        auto n_tot = inqkit::fields::density::total(*ctx.electrons);
        const double l2 = density_delta.snapshot(n_tot, ctx.time_au, step);
        inqkit::StepContext out = ctx; out.density_l2 = l2;
        obs.append(out);

        auto Rn = proj.R();
        inq::vector3<double> center{Rn.x, Rn.y, Rn.z};

        // Ehrenfest: HF force from (electrons - background). projectile_force_axis
        // is LINEAR in the field, so F(phi_e - phi_+) = F(phi_e) - F(phi_+); phi_+
        // is constant and solved once. minimum_image=true matches the perturbation.
        double Fx = 0.0, Fy = 0.0, Fz = 0.0;
        if(not CONST_V){
            auto phi_e = solvers::poisson::solve(ctx.electrons->density());
            for(int ax = 0; ax < 3; ++ax){
                const double f =
                    inqkit::dynamics::projectile_force_axis(phi_e,   center, SIGMA_POT, ax, FD_DELTA, true)
                  - inqkit::dynamics::projectile_force_axis(phiplus, center, SIGMA_POT, ax, FD_DELTA, true);
                if(ax == 0) Fx = f; else if(ax == 1) Fy = f; else Fz = f;
            }
        }

        // U_proj_bg(t) at the current centre (ideal, r_cut-free) + the pairwise
        // ledger. n_proj is rebuilt at the CURRENT position every step because the
        // classical projectile is NOT part of electrons.density()
        // (.claude/rules/decomposed-interaction-energies.md).
        auto nproj = inqkit::jellium::gaussian_density_minimum_image(basis, center, SIGMA_POT);
        const double Uprojbg = -operations::integral_product(nproj, phiplus);
        auto ct = inqkit::jellium::compute_coulomb(ctx.electrons->density(), nproj, phiplus);

        const double e_cons = ctx.energy_total + proj.ke() + Uprojbg;
        if(!e_cons_seen){ e_cons_first = e_cons; e_cons_seen = true; }
        e_cons_last = e_cons;

        proj.advance(inqkit::detail::Vec3{Fx, Fy, Fz}, DT);      // V->V_n, R->R_{n+1}

        // PERIODIC WRAP (z only — the projectile never approaches a transverse
        // face). Relabelling across a face is not a dynamical event: velocity and
        // accumulated path are untouched, only the reported coordinate moves by
        // one lattice vector. proj_z_unwrapped is the continuous path the
        // stopping fit uses.
        {
            const double z_before = proj.R().z;
            if(proj.wrap_into_cell(inqkit::detail::Vec3{0.0, 0.0, Cfg::LZ_BOHR})){
                ++n_wraps;
                std::cout << "  [wrap] step " << step << "  z " << z_before
                          << " -> " << proj.R().z << "  (wrap " << n_wraps << ")\n" << std::flush;
            }
        }
        z_unwrapped += proj.R().z - z_prev_wrapped;
        if(proj.R().z - z_prev_wrapped < -0.5*Cfg::LZ_BOHR) z_unwrapped += Cfg::LZ_BOHR;
        z_prev_wrapped = proj.R().z;

        if(ctx.electrons->root()){
            auto const V = proj.V();
            trk << step << ',' << ctx.time_au << ',' << Rn.x << ',' << Rn.y << ',' << Rn.z
                << ',' << V.x << ',' << V.y << ',' << V.z << ',' << proj.ke() << '\n';
            pj  << step << ',' << ctx.time_au << ',' << Rn.z << ',' << V.z << ','
                << proj.ke() << ',' << Uprojbg << ',' << z_unwrapped << ',' << n_wraps << ','
                << Rn.x << ',' << Rn.y << ',' << Fx << ',' << Fy << ',' << Fz << '\n';
            // Closure identities for the classical representation:
            //   E_hartree = E_SS ;  E_external = E_SB + E_PS
            ix  << step << ',' << ctx.time_au << ',' << ct.e_ss << ',' << ct.e_pp << ','
                << ct.e_ps << ',' << ct.e_sb << ',' << ct.e_pb << ',' << E_BB << ','
                << ct.e_ss << ',' << (ct.e_sb + ct.e_ps) << ','
                << ct.norm_p << ',' << ct.norm_slab << '\n';
        }

        if(SAVE_EVERY > 0 && step % SAVE_EVERY == 0) total_wr.write(n_tot, ctx.time_au, step);
    });

    // ---- checkpoint plumbing (rules: final-timestep-checkpoint + dont-block) --
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
           << "\nproj_x=" << proj.R().x << "\nproj_y=" << proj.R().y
           << "\nproj_z=" << proj.R().z
           << "\nproj_vx=" << proj.V().x << "\nproj_vy=" << proj.V().y
           << "\nproj_vz=" << proj.V().z
           << "\nproj_z_unwrapped=" << z_unwrapped << "\nn_wraps=" << n_wraps
           << "\nproj_mass=" << MASS << "\nproj_charge=" << Cfg::PROJ_CHARGE
           << "\nsigma_wp=" << SIGMA_WP << "\n";
    };

    auto step_fn = [&](auto const& data){
        rt.step(data);
        if(electrons.root())
            nlog << data.iter() << ',' << (data.iter()*DT) << ',' << data.num_electrons() << '\n';
        if(data.iter() > 0 && data.iter() % CKPT_EVERY == 0){
            electrons.save(CKPT);
            write_rt_state(data.iter());
            const std::string snap = OUT + "/ckpt_step" + tag6(data.iter());
            electrons.save(snap);
            prune_ckpts();
            std::cout << "  [ckpt] step " << data.iter() << " -> " << snap << "\n" << std::flush;
        }
    };

    real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(),
                         options::real_time{}.num_steps(N_STEPS).dt(DT * 1.0_atomictime),
                         pert, START);

    // ---- final checkpoint (saved TWICE: rolling + step-stamped) --------------
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

    // ---- correctness gate: the Ehrenfest conserved quantity ------------------
    const double e_cons_drift_ev = (e_cons_last - e_cons_first) * HA;
    {
        const bool ok = std::abs(e_cons_drift_ev) < 0.05;
        std::cout << (ok ? "\n  [PASS] " : "\n  [WARN] ")
                  << "Ehrenfest conservation: drift of (E_electronic + KE_proj + U_proj_bg) = "
                  << e_cons_drift_ev << " eV over " << (N_STEPS - START) << " steps"
                  << " (want < 0.05 eV; a large drift means the force and the "
                     "perturbation disagree)\n";
    }
    const double v_frac = (V0 != 0.0) ? proj.V().z / V0 : 0.0;
    std::cout << "  [info] final v_z/v0 = " << v_frac
              << "  (deceleration is EXPECTED for a light projectile; S is the "
                 "INITIAL drag, not a full-run slope)\n";

    if(electrons.root()){
        trk.close(); pj.close(); ix.close(); nlog.close();
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(14)
          << "run = cylindrical_jellium/proximity_ladder/classical/" << rung->label << "/" << env_s("CH_OUT","classical") << "\n"
          << "rung = " << rung->label << "\n"
          << "run_type = classical Gaussian-charge projectile, annular jellium tube TDDFT (ALDA)\n"
          << "campaign = cylindrical channeling KS-orbital stopping (twin)\n"
          << "plan = docs/plans/cylindrical-proximity-ladder.md\n"
          << "twin_role = classical\ntwin_partner = ../wp/results/" << env_s("CH_OUT_WP","wp") << "\n"
          << "engine = inq-study\nxc = LDA (ALDA in TDDFT)\n"
          << "representation = perturbation\n"
          // Single-token `projectile` value: check_twin.py's summary parser reads
          // one \S+ token, and it compares this field between the twins to prove
          // they are not two copies of the same run.
          << "projectile = gaussian_charge_perturbation\n"
          << "projectile_detail = moving Gaussian charge (velocity-Verlet, 3-D Ehrenfest)\n"
          << "drive = " << (CONST_V ? "const_velocity" : "ehrenfest") << "\n"
          << "minimum_image = yes  (charge AND force)\n"
          << "geometry = annular_tube\n"
          << "r_in_bohr = " << R_IN << "  r_out_bohr = " << R_OUT
          << "  edge_width_bohr = " << Cfg::EDGE_W_BOHR << "  tube_axis = " << Cfg::TUBE_AXIS << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR << "\n"
          << "Lz = " << Cfg::LZ_BOHR << "\nperiodicity = 3\nspacing = " << SPACING
          << "\nspacing_bohr = " << SPACING << "\n"
          << "N = " << N_ELEC << "\nn_electrons = " << N_ELEC
          << "\nn0_a0m3 = " << N0 << "\nr_s = " << RS << "\n"
          << "omega_p_au = " << OMEGA_P << "  omega_p_ev = " << OMEGA_P*HA
          << "  lambda_p_bohr = " << 2.0*M_PI*V0/OMEGA_P << "\n"
          << "v_fermi = " << VF << "  v_over_vf = " << V0/VF << "\n"
          << "n_states = " << n_states << "  extra_states = " << EXTRA_ST << "\n"
          << "sigma_WP = " << SIGMA_WP << "\nsigma_wp = " << SIGMA_WP
          << "\nsigma_pot = " << SIGMA_POT << "\n"
          << "sigma_note = the run is LABELLED by sigma_WP; sigma_pot = sigma_WP/sqrt2 "
             "is the internal charge std that matches the WP density std\n"
          << "launch_z = " << LAUNCH_Z << "\nk0 = " << V0 << "  v0 = " << V0
          << "  mass = " << MASS << "  charge = " << Cfg::PROJ_CHARGE << "\n"
          << "projectile_energy_ev = " << 0.5*MASS*V0*V0*HA << "\n"
          << "delta_fd = " << FD_DELTA << "\n"
          << "start_step = " << START << "  n_steps = " << N_STEPS << "  dt = " << DT
          << "  dt_au = " << DT << "  total_time_au = " << (DT*N_STEPS) << "\n"
          << "save_every = " << SAVE_EVERY << "  ckpt_every = " << CKPT_EVERY << "\n"
          << "proj_x_final = " << proj.R().x << "  proj_y_final = " << proj.R().y << "\n"
          << "proj_z_final = " << proj.R().z << "  proj_vz_final = " << proj.V().z << "\n"
          << "proj_z_unwrapped_final = " << z_unwrapped << "  n_wraps = " << n_wraps << "\n"
          << "v_final_over_v0 = " << v_frac << "\n"
          << "energy_conserved_first_ha = " << e_cons_first << "\n"
          << "energy_conserved_last_ha = " << e_cons_last << "\n"
          << "energy_conserved_drift_ev = " << e_cons_drift_ev << "\n"
          << "energy_conserved_note = E_electronic + KE_proj + U_proj_bg\n"
          << "date_finished = " << iso_now() << "\nwall_time_s = " << wall << "\n"
          << "gs_dir = " << GS_DIR << "\nrun_completed = true\n";
    }
    std::cout << "\nDone. Wall time " << wall << " s.  proj_z_final = "
              << proj.R().z << "  v_z_final = " << proj.V().z << "\n";
    return 0;
}
