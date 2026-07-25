// ============================================================================
// localised_jellium / scripts/fullsuite_classical / run.cpp
//
// FULL-SUITE classical-projectile re-run (Phase 5 headline, D2 2026-06-22).
// Mirrors the proven cap_baselines b2 path (Ehrenfest Gaussian-electron ion) +
// the localised background well + two-sided sin² CAP. Build ONCE against
// inq-study; reuse shared_gs/slab_n234_L50 with the background re-applied.
//
// The projectile is a classical Gaussian electron (electron-Gaussian UPF,
// sigma_pot = sigma_WP/sqrt2 = 0.35; mass = m_e), launched +z. Its kinetic-energy
// loss across the slab faces gives the stopping power S = ΔKE_ion / x.
//
// Emits: density VTIs total/system/gs_system + delta(+coarse); observables
// (energies/current/dipole/L2), state_energies, occupations, eigenvalues,
// overlap_full (bath excitations); electron_track.csv (z,v,KE every step);
// electron_number.csv (bath norm -> over-drain check).
//
// Env: LJ_OUT(p5_classical) LJ_CAP(1) LJ_N_STEPS(900) LJ_WRITE_EVERY(10)
//      LJ_LAUNCH_Z(-15.5) LJ_DT(0.02).
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/observables/state_energy_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

#include "../../../../shared/configs/slab_n82_L50x50x90_E54.hpp"
#include "../../../../../jellium/shared/cpp/eigenvalues_writer.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN82_L50x50x90_E54;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

// Gaussian-charge electron UPF at sigma_pot = sigma_WP/sqrt2 = 1.414 (the sigma=2
// companion of the quantum WP; generated + data-verified 2026-07-07).
static const char* PROJ_PSEUDO =
    "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
    "shared/pseudopotentials/electron_gaussian_wpsigma2p0.upf";

int main() {
    auto t0 = std::chrono::steady_clock::now();
    // Effective-mass projectile: m_eff = 3.0852 m_e (paired with the quantum WP's
    // inverse_mass). INQ .mass() takes amu -> divide by 1822.8885.
    const double M_PROJ = env_d("EM_MASS_ME", 3.0852) / 1822.8885;

    const bool   USE_CAP     = env_i("EM_CAP", 1) != 0;
    const std::string OUT    = "results/" + env_s("EM_OUT", "classical");
    const double SPACING     = env_d("EM_SPACING", 0.33333);   // 1.5x finer grid (matches quantum)
    const double DT_AU       = env_d("EM_DT", 0.04);
    const int    N_STEPS     = env_i("EM_N_STEPS", 2000);
    const int    WRITE_EVERY = env_i("EM_WRITE_EVERY", 20);
    const double LAUNCH_Z    = env_d("EM_LAUNCH_Z", -16.743);
    const double V0          = env_d("EM_V0", 2.7111);          // velocity (a.u.); m*v=k0=8.36 (matched momentum)
    const double CAP_ETA = -0.7, CAP_MID = 40.0/90.0, CAP_WIDTH = 10.0/90.0;

    const std::string GS_DIR = env_s("EM_GS_DIR",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "shared_gs/slab_n82_L50x50x90_dx0p333");
    if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2; }

    std::cout << "\n=== fullsuite_classical (cap=" << (USE_CAP?"on":"off") << ", out=" << OUT << ") ===\n"
              << "  N_STEPS=" << N_STEPS << " dt=" << DT_AU << " launch_z=" << LAUNCH_Z << " v0=" << V0 << "\n";

    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, Cfg::LZ_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto sp = ionic::species("H").pseudo_file(PROJ_PSEUDO).mass(M_PROJ);
    ions.insert(sp, {0.0 * 1.0_b, 0.0 * 1.0_b, LAUNCH_Z * 1.0_b});

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    electrons.load(GS_DIR);
    std::cout << "  Loaded GS (inq-study) from " << GS_DIR << "\n";
    jellium::eigenvalues::copy_from_checkpoint(GS_DIR, OUT + "/raw/observables/eigenvalues");
    const int n_states = electrons.states().num_states();

    ions.velocities()[0] = vector3<double>{0.0, 0.0, V0};

    // ----- background well (+ optional CAP) ------------------------------
    inqkit::jellium::localised_background_params bg;
    bg.shape = inqkit::jellium::background_shape::slab;
    bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
    bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);
    perturbations::absorbing cap_lo(CAP_ETA * 1.0_Ha, -CAP_MID, CAP_WIDTH);
    perturbations::absorbing cap_hi(CAP_ETA * 1.0_Ha,  CAP_MID, CAP_WIDTH);
    auto pert_bg  = bg_pert;
    auto pert_cap = perturbations::sum(bg_pert, perturbations::sum(cap_lo, cap_hi));

    // ----- output skeleton + writers ------------------------------------
    for (auto sub : {"density_total","density_system","density_gs_system",
                     "density_delta","density_delta_coarse"})
        std::filesystem::create_directories(OUT + "/raw/vti/" + sub);
    std::filesystem::create_directories(OUT + "/raw/observables/overlap_full");

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
    { inqkit::io::RealField3DWriter gs_wr(OUT + "/raw/vti/density_gs_system", vti_layout, {.overwrite=true});
      gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system"); }
    inqkit::io::RealField3DWriter total_wr (OUT + "/raw/vti/density_total",  vti_layout, {.overwrite=true});
    inqkit::io::RealField3DWriter system_wr(OUT + "/raw/vti/density_system", vti_layout, {.overwrite=true});
    { auto s0 = inqkit::fields::density::total(electrons); total_wr.write(s0,0.0,0); system_wr.write(s0,0.0,0); }

    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables.csv", sel);
    obs_writer.write_header();
    inqkit::observables::StateEnergyWriter state_energy_wr(OUT + "/raw/observables/state_energies.csv", true);
    inqkit::observables::OccupationsWriter occupations_wr(OUT + "/raw/observables/occupations_vs_time.csv");
    inqkit::observables::DensityDelta density_delta(
        OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
        {.emit_raw_vti=true, .emit_coarse_vti=true, .compute_l2=true, .coarse_bin_bohr=3.0});
    inqkit::observables::OrbitalOverlapMatrix overlap_full_obs(electrons, n_states - 1, OUT + "/raw/observables/overlap_full");
    overlap_full_obs.snapshot(electrons, 0.0, 0);

    // KE in atomic units uses the electron mass = 1 a.u. (NOT the amu value
    // M_PROJ passed to INQ's .mass(), which is 1/1822.8885 amu). Logging with
    // M_PROJ would under-scale KE by ~1822x.
    auto ke_ion = [&](){ auto const& v = ions.velocities()[0];
        return 0.5*1.0*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); };
    std::ofstream trk(OUT + "/raw/observables/electron_track.csv");
    trk << std::setprecision(16) << "step,time_au,x,y,z,vx,vy,vz,ke_ion_ha\n";
    { auto const& p = ions.positions()[0]; auto const& v = ions.velocities()[0];
      trk << 0 << ",0," << p[0] << "," << p[1] << "," << p[2] << "," << v[0] << "," << v[1] << "," << v[2] << "," << ke_ion() << "\n"; }
    std::ofstream nlog(OUT + "/raw/observables/electron_number.csv");
    nlog << std::setprecision(12) << "step,time_au,N_total\n";

    // ----- two-segment park+neutralise propagation --------------------------
    // Segment 1: Ehrenfest WITH the projectile, issued in fresh start_step=0
    // CHUNKS. (INQ restart with moving ions is unsupported — propagate.hpp:15 —
    // so each chunk is a *fresh* propagate; the persistent ions/electrons objects
    // carry the state forward.) After each chunk we test the projectile z; once
    // |z| >= the CAP inner face (35 Bohr) the projectile is PARKED and removed via
    // ions.remove(0), so its Gaussian-Coulomb potential becomes EXACTLY zero.
    // Segment 2 then continues projectile-free to N_STEPS so the slab relaxes and
    // the CAP absorbs the induced flux. If the projectile never reaches |z|>=35
    // within N_STEPS, the whole run is plain chunked Ehrenfest (park never fires).
    //
    // Numbering: each propagate fires the callback at iter=0 (the carried-over
    // state) then iter=1..n (propagate.hpp:81,123). We keep a GLOBAL step `g`,
    // SKIP iter==0 on every chunk after the first (it duplicates the previous
    // chunk's last state), and stamp physical time as g*DT_AU — NOT data.time(),
    // which is chunk-local and would reset each chunk.
    const double CAP_INNER_BOHR  = env_d("LJ_PARK_Z", 35.0);   // |z| park trigger (CAP inner face)
    const int    CHUNK           = env_i("LJ_CHUNK", 25);      // steps/chunk -> boundary granularity
    const int    PARK_STEP_FORCE = env_i("LJ_PARK_STEP", 0);   // SMOKE hook: force park at this global step

    int    g               = 0;            // global unique step index
    bool   first_chunk     = true;
    bool   projectile_live = true;         // false once parked + removed
    double parked_x = 0.0, parked_y = 0.0, parked_z = LAUNCH_Z;

    auto func = [&](auto const& data) {
        if (data.iter() == 0 && !first_chunk) return;     // skip carried-over duplicate state
        const double t = g * DT_AU;
        if (g % WRITE_EVERY == 0) {
            auto sys_f = inqkit::fields::density::total(electrons);
            system_wr.write(sys_f, t, g);
            total_wr.write (sys_f, t, g);
            const double l2 = density_delta.snapshot(sys_f, t, g);
            inqkit::StepContext c;
            c.step = g; c.time_au = t;
            c.energy_total   = data.energy().total();
            c.energy_kinetic = data.energy().kinetic();
            c.energy_hartree = data.energy().hartree();
            c.energy_xc      = data.energy().xc();
            try { auto cu = data.current(); c.current = {cu[0], cu[1], cu[2]}; } catch (...) {}
            try { auto dp = data.dipole();  c.dipole  = {dp[0], dp[1], dp[2]}; } catch (...) {}
            c.density_l2 = l2;
            obs_writer.append(c);
        }
        if (g % (5 * WRITE_EVERY) == 0) { state_energy_wr.snapshot(data); occupations_wr.snapshot(data); }
        if (data.root()) {
            nlog << g << "," << t << "," << data.num_electrons() << "\n";
            double px, py, pz, vx, vy, vz, ke;
            if (projectile_live) {
                auto const& p = ions.positions()[0]; auto const& v = ions.velocities()[0];
                px=p[0]; py=p[1]; pz=p[2]; vx=v[0]; vy=v[1]; vz=v[2]; ke=ke_ion();
                parked_x=px; parked_y=py; parked_z=pz;     // remember last live state
            } else {
                px=parked_x; py=parked_y; pz=parked_z; vx=0; vy=0; vz=0; ke=0;  // inert: zero potential
            }
            trk << g << "," << t << "," << px << "," << py << "," << pz << ","
                << vx << "," << vy << "," << vz << "," << ke << "\n";
        }
        ++g;
    };

    int done = 0, park_step = -1;
    // SEGMENT 1 — projectile present, chunked until |z| >= CAP inner face
    while (done < N_STEPS && projectile_live) {
        const int n = (CHUNK < N_STEPS - done) ? CHUNK : (N_STEPS - done);
        auto opts = options::real_time{}.num_steps(n).dt(DT_AU * 1.0_atomictime)
                        .observables_current().observables_dipole().ehrenfest();
        if (USE_CAP) real_time::propagate(ions, electrons, func, options::theory{}.lda(), opts, pert_cap);
        else         real_time::propagate(ions, electrons, func, options::theory{}.lda(), opts, pert_bg);
        first_chunk = false;
        done += n;
        const double z  = ions.positions()[0][2];
        const double az = (z < 0.0 ? -z : z);
        if (az >= CAP_INNER_BOHR || (PARK_STEP_FORCE > 0 && g >= PARK_STEP_FORCE)) {
            parked_x = ions.positions()[0][0]; parked_y = ions.positions()[0][1]; parked_z = z;
            projectile_live = false; park_step = g;
            std::cout << "  [park] projectile reached |z|=" << az << " at step " << g
                      << " (t=" << g*DT_AU << " au) -> ions.remove(0): zero potential\n";
        }
    }
    // NEUTRALISE — remove the projectile so its radial potential is exactly zero
    const bool parked = !projectile_live;
    if (parked) ions.remove(0);
    // SEGMENT 2 — projectile-free, remaining steps to N_STEPS
    const int rem = N_STEPS - done;
    if (parked && rem > 0) {
        auto opts2 = options::real_time{}.num_steps(rem).dt(DT_AU * 1.0_atomictime)
                        .observables_current().observables_dipole();   // no ions to move
        if (USE_CAP) real_time::propagate(ions, electrons, func, options::theory{}.lda(), opts2, pert_cap);
        else         real_time::propagate(ions, electrons, func, options::theory{}.lda(), opts2, pert_bg);
        done += rem;
    }

    overlap_full_obs.snapshot(electrons, DT_AU * g, g);

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        const double fz  = parked ? parked_z : ions.positions()[0][2];
        const double fvz = parked ? 0.0      : ions.velocities()[0][2];
        s << std::setprecision(12)
          << "run = localised_jellium/qsp_phase4_classical/" << env_s("LJ_OUT","p4_classical") << "\n"
          << "engine = inq-study\n"
          << "projectile = classical Gaussian-e ion (sigma_pot 0.35, mass m_e); Ehrenfest then park+remove at |z|>=35 (zero potential)\n"
          << "cap = " << (USE_CAP?"on (two-sided sin2, eta -0.7 Ha, 10 Bohr/side, region +/-35..+/-45, inner faces +/-35)":"off") << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR << "  spacing = " << SPACING << "\n"
          << "background = slab half_width " << Cfg::SLAB_HALF_WIDTH << " axis " << Cfg::SLAB_AXIS << "\n"
          << "n_electrons = " << Cfg::N_ELECTRONS << "  n_states = " << n_states << "\n"
          << "launch_z = " << LAUNCH_Z << "  v0 = " << V0 << "  ke_ion_initial_ha = " << (0.5*V0*V0) << "\n"
          << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY << "  chunk = " << CHUNK << "\n"
          << "park_triggered = " << (parked?"true":"false") << "  park_step = " << park_step
          << "  park_t_au = " << (park_step >= 0 ? park_step*DT_AU : -1.0) << "  park_z = " << parked_z << "\n"
          << "steps_run = " << g << "  final_z = " << fz << "  final_vz = " << fvz << "\n"
          << "wall_time_s = " << wall << "\nrun_completed = true\n";
    }
    std::cout << "  done. parked=" << (parked?"yes":"no")
              << " final_z=" << (parked ? parked_z : ions.positions()[0][2]) << " wall=" << wall << "s\n";
    return 0;
}
