// ============================================================================
// localised_jellium / muon_mass_fork / effmass_sigma1 / classical / run.cpp
//
// CLASSICAL TWIN of the sigma_WP=1 chirped-focus WP run (LEAN re-run; plan:
// docs/plans/effmass-sigma1-lean-rerun.md). Rigid Gaussian-charge electron ion
// (electron_gaussian_wpsigma1p0.upf, sigma_pot = 1/sqrt2 = 0.707, generated +
// V(r)-verified 2026-07-09), effective mass m = 2.10 m_e (matched to the WP's
// inverse_mass), v0 = 2.7111, free Ehrenfest through the density-matched
// 40x40x80 (N=52, r_s=5.68) slab at dx=0.33333, dt=0.04, N=900 steps.
//
// Mirrors the proven fullsuite_classical path: chunked Ehrenfest (INQ restart
// with moving ions is unsupported -> fresh start_step=0 chunks over persistent
// ions/electrons), park+remove once |z| >= CAP inner face (25 Bohr for L80).
//
// CHECKPOINTS (user req: 3, to extend later): every EM_CKPT_EVERY (300) global
// steps, at a chunk boundary: electrons.save(rt_ckpt) + rt_state.txt holding
// last_step, projectile_live, ion position+velocity. Extension = relaunch with
// EM_RESUME=1 (loads ckpt, restores ion state, CSVs get .from<START> suffix;
// density_delta reference re-anchored to the GS density).
//
// FIX vs effmass_pair template: ke_ion uses the ACTUAL projectile mass
// (0.5*M_ME*v^2), not 0.5*1.0*v^2 — with m=2.10 the template would under-log
// KE by 2.1x. (INQ .mass() itself takes amu = M_ME/1822.8885.)
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/observables/state_energy_writer.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

#include "../../../../shared/configs/slab_n52_L40x40x80.hpp"
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
using Cfg = localised_jellium::config::SlabN52_L40x40x80;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

// Gaussian-charge electron UPF at sigma_pot = sigma_WP/sqrt2 = 0.707 (the sigma=1
// companion of the quantum WP; generated + data-verified by V(r) 2026-07-09).
static const char* PROJ_PSEUDO =
    "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
    "shared/pseudopotentials/electron_gaussian_wpsigma1p0.upf";

static int read_kv_int(const std::string& path, const std::string& key, int dflt) {
    std::ifstream f(path); if (!f) return dflt;
    std::string line;
    while (std::getline(f, line))
        if (line.rfind(key + "=", 0) == 0) return std::atoi(line.c_str() + key.size() + 1);
    return dflt;
}
static double read_kv_d(const std::string& path, const std::string& key, double dflt) {
    std::ifstream f(path); if (!f) return dflt;
    std::string line;
    while (std::getline(f, line))
        if (line.rfind(key + "=", 0) == 0) return std::atof(line.c_str() + key.size() + 1);
    return dflt;
}

int main() {
    auto t0 = std::chrono::steady_clock::now();
    // Effective-mass projectile: m_eff = 2.10 m_e (matched to the quantum WP's
    // inverse_mass = 0.476190). INQ .mass() takes amu -> divide by 1822.8885.
    const double M_ME   = env_d("EM_MASS_ME", 2.10);
    const double M_PROJ = M_ME / 1822.8885;

    const bool   USE_CAP     = env_i("EM_CAP", 1) != 0;
    const std::string OUT    = "results/" + env_s("EM_OUT", "classical");
    const double SPACING     = env_d("EM_SPACING", 0.33333);   // matches quantum twin
    const double DT_AU       = env_d("EM_DT", 0.04);
    const int    N_STEPS     = env_i("EM_N_STEPS", 900);
    const int    WRITE_EVERY = env_i("EM_WRITE_EVERY", 20);
    const int    CKPT_EVERY  = env_i("EM_CKPT_EVERY", 300);    // 3 checkpoints over 900
    const bool   RESUME      = env_i("EM_RESUME", 0) != 0;
    const double LAUNCH_Z    = env_d("EM_LAUNCH_Z", -16.5);    // same standoff as the WP
    const double V0          = env_d("EM_V0", 2.7111);         // m*v = k0 = 5.693 (matched momentum)
    // CAP for Lz=80: region [25,40] each side — centre |z|=32.5, full width 15 Bohr
    // (reflectivity-tuned width kept), inner faces +/-25.
    const double CAP_ETA = -1.0, CAP_MID = 32.5/80.0, CAP_WIDTH = 15.0/80.0;
    const double CAP_INNER_BOHR  = env_d("LJ_PARK_Z", 25.0);   // |z| park trigger (CAP inner face)
    const int    CHUNK           = env_i("LJ_CHUNK", 25);      // steps/chunk -> boundary granularity
    const int    PARK_STEP_FORCE = env_i("LJ_PARK_STEP", 0);   // SMOKE hook: force park at this global step

    const std::string GS_DIR = env_s("EM_GS_DIR",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "shared_gs/slab_n52_L40x40x80_dx0p333");
    const std::string RT_CKPT_DIR  = env_s("EM_RT_CKPT_DIR", OUT + "/rt_ckpt");
    const std::string RT_STATE_TXT = RT_CKPT_DIR + "/rt_state.txt";
    if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2; }

    // ----- resume bookkeeping ---------------------------------------------
    int START = 0; bool resume_live = true;
    double r_x=0, r_y=0, r_z=LAUNCH_Z, r_vx=0, r_vy=0, r_vz=V0;
    if (RESUME) {
        if (!std::filesystem::exists(RT_CKPT_DIR)) { std::cerr << "FATAL: resume but no RT ckpt: " << RT_CKPT_DIR << "\n"; return 2; }
        START = read_kv_int(RT_STATE_TXT, "last_step", -1);
        if (START < 0) { std::cerr << "FATAL: resume but unreadable " << RT_STATE_TXT << "\n"; return 2; }
        if (START >= N_STEPS) { std::cout << "Already at/after target (" << START << ">=" << N_STEPS << "); nothing to do.\n"; return 0; }
        resume_live = read_kv_int(RT_STATE_TXT, "projectile_live", 1) != 0;
        r_x  = read_kv_d(RT_STATE_TXT, "x",  0.0);  r_y  = read_kv_d(RT_STATE_TXT, "y",  0.0);
        r_z  = read_kv_d(RT_STATE_TXT, "z",  LAUNCH_Z);
        r_vx = read_kv_d(RT_STATE_TXT, "vx", 0.0);  r_vy = read_kv_d(RT_STATE_TXT, "vy", 0.0);
        r_vz = read_kv_d(RT_STATE_TXT, "vz", V0);
    }
    const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

    std::cout << "\n=== effmass_sigma1 CLASSICAL twin (cap=" << (USE_CAP?"on":"off")
              << ", out=" << OUT << (RESUME?"  [RESUME]":"  [FRESH]") << ") ===\n"
              << "  m=" << M_ME << " m_e  N_STEPS=" << N_STEPS << " dt=" << DT_AU
              << " launch_z=" << LAUNCH_Z << " v0=" << V0 << " start=" << START << "\n";

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

    // ----- observables plumbing that must know the reference density -------
    inqkit::observables::DensityDelta density_delta(
        OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
        {.emit_raw_vti=true, .emit_coarse_vti=true, .compute_l2=true, .coarse_bin_bohr=3.0});
    std::filesystem::create_directories(OUT + "/raw/vti/density_delta");
    std::filesystem::create_directories(OUT + "/raw/vti/density_delta_coarse");

    if (RESUME) {
        electrons.load(GS_DIR);                                   // (1) anchor delta ref to the GS
        density_delta.set_reference(inqkit::fields::density::total(electrons));
        electrons.load(RT_CKPT_DIR);                              // (2) then the propagated state
        ions.positions()[0]  = vector3<double>{r_x, r_y, r_z};
        ions.velocities()[0] = vector3<double>{r_vx, r_vy, r_vz};
        std::cout << "  Resumed RT ckpt @step " << START << " (live=" << resume_live
                  << " z=" << r_z << " vz=" << r_vz << ")\n";
    } else {
        electrons.load(GS_DIR);
        std::cout << "  Loaded GS from " << GS_DIR << "\n";
        ions.velocities()[0] = vector3<double>{0.0, 0.0, V0};
    }
    jellium::eigenvalues::copy_from_checkpoint(GS_DIR, OUT + "/raw/observables/eigenvalues");
    const int n_states = electrons.states().num_states();

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
    for (auto sub : {"density_total","density_system","density_gs_system"})
        std::filesystem::create_directories(OUT + std::string("/raw/vti/") + sub);
    std::filesystem::create_directories(OUT + "/raw/observables/overlap_full");

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
    if (!RESUME) {
        inqkit::io::RealField3DWriter gs_wr(OUT + "/raw/vti/density_gs_system", vti_layout, {.overwrite=true});
        gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system");
    }
    inqkit::io::RealField3DWriter total_wr (OUT + "/raw/vti/density_total",  vti_layout, {.overwrite=true});
    inqkit::io::RealField3DWriter system_wr(OUT + "/raw/vti/density_system", vti_layout, {.overwrite=true});
    if (!RESUME) { auto s0 = inqkit::fields::density::total(electrons); total_wr.write(s0,0.0,0); system_wr.write(s0,0.0,0); }

    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables" + SEG + ".csv", sel);
    obs_writer.write_header();
    inqkit::observables::StateEnergyWriter state_energy_wr(OUT + "/raw/observables/state_energies" + SEG + ".csv", true);
    inqkit::observables::OccupationsWriter occupations_wr(OUT + "/raw/observables/occupations_vs_time" + SEG + ".csv");
    inqkit::observables::OrbitalOverlapMatrix overlap_full_obs(electrons, n_states - 1, OUT + "/raw/observables/overlap_full");
    if (!RESUME) overlap_full_obs.snapshot(electrons, 0.0, 0);

    // KE ledger uses the ACTUAL projectile mass in a.u. (m_e units): 0.5*M_ME*v^2.
    // (NOT the amu value M_PROJ handed to .mass(); and NOT the template's 0.5*1.0*v^2,
    // which under-logs KE by M_ME x for an effective-mass projectile.)
    auto ke_ion = [&](){ auto const& v = ions.velocities()[0];
        return 0.5*M_ME*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); };
    std::ofstream trk(OUT + "/raw/observables/electron_track" + SEG + ".csv");
    trk << std::setprecision(16) << "step,time_au,x,y,z,vx,vy,vz,ke_ion_ha\n";
    std::ofstream nlog(OUT + "/raw/observables/electron_number" + SEG + ".csv");
    nlog << std::setprecision(12) << "step,time_au,N_total\n";

    // ----- chunked park+neutralise propagation (see fullsuite_classical) ----
    // Chunked fresh-start propagates (moving-ion restart unsupported in INQ);
    // global step g; iter==0 skipped on every chunk after the first (duplicate).
    int    g               = START;        // global unique step index
    bool   first_chunk     = !RESUME;      // on resume the ckpt state is already logged
    bool   projectile_live = RESUME ? resume_live : true;
    double parked_x = r_x, parked_y = r_y, parked_z = r_z;

    auto save_ckpt = [&](int at_step){
        std::filesystem::create_directories(RT_CKPT_DIR);
        electrons.save(RT_CKPT_DIR);
        if (electrons.root()) {
            double px=parked_x, py=parked_y, pz=parked_z, vx=0, vy=0, vz=0;
            if (projectile_live) {
                auto const& p = ions.positions()[0]; auto const& v = ions.velocities()[0];
                px=p[0]; py=p[1]; pz=p[2]; vx=v[0]; vy=v[1]; vz=v[2];
            }
            std::ofstream s(RT_STATE_TXT);
            s << std::setprecision(17)
              << "last_step=" << at_step << "\nprojectile_live=" << (projectile_live?1:0) << "\n"
              << "x=" << px << "\ny=" << py << "\nz=" << pz << "\n"
              << "vx=" << vx << "\nvy=" << vy << "\nvz=" << vz << "\n";
            std::cout << "  [ckpt] step " << at_step << " saved -> " << RT_CKPT_DIR << "\n";
        }
    };

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

    int done = START, park_step = -1;
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
        if (done % CKPT_EVERY == 0 && done < N_STEPS) save_ckpt(done);
    }
    // NEUTRALISE — remove the projectile so its radial potential is exactly zero
    const bool parked = !projectile_live;
    if (parked && (!RESUME || resume_live)) ions.remove(0);
    else if (RESUME && !resume_live) ions.remove(0);
    // SEGMENT 2 — projectile-free, chunked so checkpoints keep firing
    while (done < N_STEPS) {
        const int n = (CHUNK < N_STEPS - done) ? CHUNK : (N_STEPS - done);
        auto opts2 = options::real_time{}.num_steps(n).dt(DT_AU * 1.0_atomictime)
                        .observables_current().observables_dipole();   // no ions to move
        if (USE_CAP) real_time::propagate(ions, electrons, func, options::theory{}.lda(), opts2, pert_cap);
        else         real_time::propagate(ions, electrons, func, options::theory{}.lda(), opts2, pert_bg);
        first_chunk = false;
        done += n;
        if (done % CKPT_EVERY == 0 && done < N_STEPS) save_ckpt(done);
    }
    save_ckpt(done);   // final checkpoint (enables extension beyond N_STEPS)

    overlap_full_obs.snapshot(electrons, DT_AU * g, g);

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        const double fz  = parked ? parked_z : ions.positions()[0][2];
        const double fvz = parked ? 0.0      : ions.velocities()[0][2];
        s << std::setprecision(12)
          << "run = localised_jellium/muon_mass_fork/effmass_sigma1/" << env_s("EM_OUT","classical") << "\n"
          << "engine = inq-study\n"
          << "projectile = classical Gaussian-e ion (sigma_pot 0.7071 = sigma_WP 1/sqrt2, mass " << M_ME
          << " m_e); Ehrenfest then park+remove at |z|>=" << CAP_INNER_BOHR << "\n"
          << "cap = " << (USE_CAP?"on (two-sided sin2, eta -1.0 Ha, 15 Bohr/side, region +/-25..+/-40, inner faces +/-25)":"off") << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR << "  spacing = " << SPACING << "\n"
          << "background = slab half_width " << Cfg::SLAB_HALF_WIDTH << " axis " << Cfg::SLAB_AXIS << "\n"
          << "n_electrons = " << Cfg::N_ELECTRONS << "  n_states = " << n_states << "\n"
          << "launch_z = " << LAUNCH_Z << "  v0 = " << V0 << "  mass_me = " << M_ME
          << "  ke_ion_initial_ha = " << (0.5*M_ME*V0*V0) << "\n"
          << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY
          << "  chunk = " << CHUNK << "  ckpt_every = " << CKPT_EVERY << "\n"
          << "park_triggered = " << (parked?"true":"false") << "  park_step = " << park_step
          << "  park_t_au = " << (park_step >= 0 ? park_step*DT_AU : -1.0) << "  park_z = " << parked_z << "\n"
          << "steps_run = " << g << "  final_z = " << fz << "  final_vz = " << fvz << "\n"
          << "wall_time_s = " << wall << "\nrun_completed = true\n";
    }
    std::cout << "  done. parked=" << (parked?"yes":"no")
              << " final_z=" << (parked ? parked_z : ions.positions()[0][2]) << " wall=" << wall << "s\n";
    return 0;
}
