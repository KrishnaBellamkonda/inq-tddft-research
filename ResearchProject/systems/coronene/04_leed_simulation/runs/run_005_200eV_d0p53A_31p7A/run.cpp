// ============================================================================
// run_005 — Coronene TDDFT LEED simulation (paper parameters)
//
// Reproduces Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014), Fig. 2(a):
//   200 eV electron WP scattering off coronene C24H12, reflection LEED.
//
// Key differences from run_004:
//   - Lz = 31.7 Å  (paper value)
//   - D  = 6.35 Å  (paper value; WP spread at flake ≈ 1.7 Å, near C-C scale)
//   - T2 = 0.25 fs  (paper value; ends when transmitted WP reaches z=0)
//   - LEED = ∫_{t1}^{t2} n_total(z_obs) dt  (paper Eq. 5; total density)
//   - t1  = T1 = D/k0 = 0.077 fs  (start: WP arrives at flake)
//   - One screen at z_obs
//
// Removed from run_004:
//   - 3D GS orbital saves (large and slow, not needed for LEED)
//   - 3D density saves (large, not needed)
//   - WP orbital 3D saves (large)
//   - 3 LEED screens → 1 screen
//   - Background subtraction from LEED accumulator → matches paper exactly
//
// Source: Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014)
// ============================================================================

#include "config.hpp"
#include "utils.hpp"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#define GPU_SYNC() cudaDeviceSynchronize()
#else
#define GPU_SYNC() ((void)0)
#endif

int main(){
    using namespace inq;
    using namespace inq::magnitude;

    auto wall_start = std::chrono::steady_clock::now();

    auto pad6 = [](int n){ std::ostringstream s; s << std::setfill('0') << std::setw(6) << n; return s.str(); };

    // ── STEP 1: Geometry ──────────────────────────────────────────────────────
    auto cell = cfg::make_cell();
    auto ions = systems::ions::parse(cfg::CORONENE_XYZ, cell);

    std::cout << "\n=== run_005: Coronene TDDFT LEED (paper parameters) ===\n";
    std::cout << "  Cell   : " << cfg::LX_BOHR << " × " << cfg::LY_BOHR
              << " × " << cfg::LZ_BOHR << " bohr  [finite]\n";
    std::cout << "  Lz     : " << cfg::LZ_ANG  << " Å  (paper: 31.7 Å)\n";
    std::cout << "  d      : " << cfg::WP_D_ANG << " Å  (paper: 0.53 Å)\n";
    std::cout << "  D      : " << cfg::WP_D_IMPACT_ANG << " Å  (paper: 6.35 Å)\n";
    std::cout << "  T1     : " << cfg::T1_FS << " fs  (paper: 0.077 fs)\n";
    std::cout << "  T2     : " << cfg::T2_FS << " fs  (paper: 0.25 fs)\n";
    std::cout << "  Atoms  : " << ions.size() << "  (expect 36)\n";

    // ── STEP 2: Output directories ────────────────────────────────────────────
    for(auto d : {"results", "results/energy", "results/momentum",
                  "results/grid", "results/leed_pattern",
                  "results/density_snapshots", "results/density_obs_snapshots",
                  "results/wp_trajectory"})
        leed_utils::mkdir_p(d);

    // ── STEP 3: Ground state SCF ──────────────────────────────────────────────
    auto electrons = systems::electrons(ions,
        options::electrons{}.cutoff(cfg::ECUT_HA * 1.0_Ha)
                            .extra_states(cfg::EXTRA_STATES));

    std::cout << "\n=== STEP 2: Ground state SCF ===\n";
    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(cfg::SCF_TOL * 1.0_Ha)
            .max_steps(cfg::SCF_MAX_STEPS)
            .broyden_mixing()
            .mixing_ndim(cfg::SCF_MIXING_NDIM)
            .mixing(cfg::SCF_MIXING));

    double E_gs = gs.energy.total();
    std::cout << "  GS energy : " << E_gs << " Ha  SCF steps: " << gs.total_iter << "\n";
    if(gs.total_iter >= cfg::SCF_MAX_STEPS)
        std::cout << "  WARNING: SCF may not be fully converged.\n";

    // Save grid coords and GS energy
    leed_utils::save_grid_coords(electrons, "results/grid");
    {
        std::ofstream f("results/energy/gs_energy.txt");
        f << "GS_energy_Ha  " << std::fixed << std::setprecision(10) << E_gs << "\n"
          << "SCF_steps     " << gs.total_iter << "\n";
    }

    // Capture GS baseline at z_obs before WP injection.
    // Saved for post-processing background subtraction in analysis.py.
    // The in-run LEED accumulator uses total density (paper Eq.5).
    GPU_SYNC();
    auto gs_slice_obs = leed_utils::extract_density_slice(electrons, cfg::Z_SCREEN_BOHR());
    {
        std::ofstream f("results/leed_pattern/gs_baseline_z_obs.txt");
        f << "# GS density at z_obs=" << cfg::Z_SCREEN_BOHR() << " bohr\n";
        f << "# Rows: iy=0.." << gs_slice_obs.size()-1
          << "  Cols: ix=0.." << gs_slice_obs[0].size()-1 << "\n";
        for(auto& row : gs_slice_obs){
            for(int ix = 0; ix < (int)row.size(); ix++){
                if(ix) f << " ";
                f << std::scientific << std::setprecision(8) << row[ix];
            }
            f << "\n";
        }
    }
    std::cout << "  Saved GS baseline at z_obs for post-processing.\n";

    // ── STEP 4: WP injection ──────────────────────────────────────────────────
    std::cout << "\n=== STEP 3: WP injection ===\n";
    const double k0     = cfg::wp_k0();
    const int    ist_wp = electrons.kpin()[0].set_part().local_size() - 1;

    leed_utils::inject_wp(electrons,
        cfg::WP_BX(), cfg::WP_BY(), cfg::WP_BZ(),
        0.0, 0.0, -k0);
    electrons.occupations()[0][ist_wp] = cfg::WP_OCCUPATION;

    auto [wp_norm, wp_ke] = leed_utils::validate_wp(electrons);
    std::cout << "  WP centre: (" << cfg::WP_BX() << ", " << cfg::WP_BY()
              << ", " << cfg::WP_BZ() << ") bohr\n";
    std::cout << "  z_flake  : " << cfg::Z_FLAKE_BOHR() << " bohr = "
              << cfg::Z_FLAKE_BOHR()*cfg::BOHR_TO_ANG << " Å\n";
    std::cout << "  z_obs    : " << cfg::Z_SCREEN_BOHR() << " bohr = "
              << cfg::Z_SCREEN_BOHR()*cfg::BOHR_TO_ANG << " Å\n";
    std::cout << "  k0       : " << k0 << " bohr⁻¹  (kz = -k0)\n";
    std::cout << "  WP norm  : " << wp_norm << "\n";
    if(std::abs(wp_norm - 1.0) > 0.05)
        std::cout << "  WARNING: norm deviates >5%\n";

    // ── STEP 5: TDDFT propagation ─────────────────────────────────────────────
    std::cout << "\n=== STEP 4: TDDFT propagation ===\n";
    std::cout << "  Δt   : " << cfg::DT_AU << " a.u.  (" << cfg::DT_FS << " fs)\n";
    std::cout << "  t1   : " << cfg::T1_AU << " a.u.  (" << cfg::T1_FS << " fs)  <- LEED start\n";
    std::cout << "  t2   : " << cfg::T2_AU << " a.u.  (" << cfg::T2_FS << " fs)  <- LEED end\n";
    std::cout << "  N    : " << cfg::N_STEPS << " steps\n";
    std::cout << "  LEED : ∫_{t1}^{t2} n_total(z_obs) dt  (paper Eq. 5)\n\n";

    auto const & basis = electrons.states_basis();
    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    int Nz_g = basis.sizes()[2];

    // LEED accumulator: paper Eq. 5 — total density, no background subtraction
    // I(x,y) = ∫_{t1}^{t2} n_total(x,y,z_obs,t) dt
    std::vector<std::vector<double>> leed(Ny_g, std::vector<double>(Nx_g, 0.0));
    double t1 = cfg::T1_AU;
    double t2 = cfg::T2_AU;

    std::ofstream energy_csv("results/energy/energy_vs_time.csv");
    energy_csv << "# step,t_au,E_total_Ha\n";

    std::ofstream momentum_csv("results/momentum/momentum_vs_time.csv");
    momentum_csv << "# step,t_au,Jx_au,Jy_au,Jz_au\n";

    std::ofstream zprofile_csv("results/wp_trajectory/density_z_profile_vs_time.csv");
    zprofile_csv << "# step,t_au";
    for(int iz = 0; iz < Nz_g; iz++) zprofile_csv << ",n_iz" << iz;
    zprofile_csv << "\n";

    auto obs_callback = [&](auto && obs){
        GPU_SYNC();

        double t    = obs.time();
        int    iter = obs.iter();
        auto   e    = obs.energy();
        auto   J    = obs.current();

        // Every step: energy and momentum
        energy_csv << iter << "," << std::fixed << std::setprecision(8)
                   << t << "," << e.total() << "\n";
        momentum_csv << iter << "," << std::fixed << std::setprecision(8)
                     << t << ","
                     << std::scientific << std::setprecision(8)
                     << J[0] << "," << J[1] << "," << J[2] << "\n";

        // LEED: paper Eq. 5 — accumulate total density during [t1, t2]
        // t1 = WP arrives at flake; t2 = WP reaches bottom boundary.
        // During this window, z_obs sees only the backscattered electrons.
        if(t >= t1 && t <= t2){
            auto s = leed_utils::extract_density_slice(electrons, cfg::Z_SCREEN_BOHR());
            for(int iy = 0; iy < Ny_g; iy++)
                for(int ix = 0; ix < Nx_g; ix++)
                    leed[iy][ix] += s[iy][ix] * cfg::DT_AU;
        }

        // Every SNAPSHOT_INTERVAL steps: lightweight 2D slices + z-profile
        if(iter % cfg::SNAPSHOT_INTERVAL != 0) return;

        std::cout << "  step=" << iter << "  t=" << std::fixed << std::setprecision(4)
                  << t << " a.u.  E=" << std::setprecision(6) << e.total()
                  << " Ha  Jz=" << std::scientific << std::setprecision(3) << J[2]
                  << (t >= t1 && t <= t2 ? "  [LEED]" : "") << "\n";

        // 2D density at z_flake (molecular plane)
        {
            auto s = leed_utils::extract_density_slice(electrons, cfg::Z_FLAKE_BOHR());
            std::ofstream f("results/density_snapshots/snapshot_t" + pad6(iter) + ".txt");
            f << "# t=" << std::fixed << std::setprecision(6) << t
              << " z=" << cfg::Z_FLAKE_BOHR()*cfg::BOHR_TO_ANG << "\n";
            for(auto& row : s){
                for(int ix = 0; ix < (int)row.size(); ix++){
                    if(ix) f << " ";
                    f << std::scientific << std::setprecision(6) << row[ix];
                }
                f << "\n";
            }
        }

        // 2D density at z_obs
        {
            auto s = leed_utils::extract_density_slice(electrons, cfg::Z_SCREEN_BOHR());
            std::ofstream f("results/density_obs_snapshots/snapshot_t" + pad6(iter) + ".txt");
            f << "# t=" << std::fixed << std::setprecision(6) << t
              << " z=" << cfg::Z_SCREEN_BOHR()*cfg::BOHR_TO_ANG << "\n";
            for(auto& row : s){
                for(int ix = 0; ix < (int)row.size(); ix++){
                    if(ix) f << " ";
                    f << std::scientific << std::setprecision(6) << row[ix];
                }
                f << "\n";
            }
        }

        // 1D z-profile
        {
            auto const & phi = electrons.kpin()[0];
            auto hc = begin(phi.hypercubic());
            int n_st = phi.set_part().local_size();
            double dV = basis.volume_element();
            int ix_c = Nx_g / 2;
            int iy_c = Ny_g / 2;
            zprofile_csv << iter << "," << std::fixed << std::setprecision(8) << t;
            for(int iz = 0; iz < Nz_g; iz++){
                double n = 0.0;
                for(int ist = 0; ist < n_st; ist++){
                    auto v = hc[ix_c][iy_c][iz][ist];
                    n += v.real()*v.real() + v.imag()*v.imag();
                }
                zprofile_csv << "," << std::scientific << std::setprecision(6) << n;
            }
            zprofile_csv << "\n";
        }
    };

    real_time::propagate(ions, electrons,
        obs_callback,
        options::theory{}.lda(),
        options::real_time{}
            .num_steps(cfg::N_STEPS)
            .dt(cfg::DT_AU * 1.0_atomictime)
            .observables_current());

    // ── STEP 6: Write LEED pattern ────────────────────────────────────────────
    std::cout << "\n=== STEP 5: Writing outputs ===\n";
    {
        std::ofstream f("results/leed_pattern/leed_screen.txt");
        f << "# LEED pattern I(x,y) = integral_{t1}^{t2} n_total(x,y,z_obs,t) dt\n";
        f << "# Paper Eq.5: Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014)\n";
        f << "# z_obs=" << std::fixed << std::setprecision(4)
          << cfg::Z_SCREEN_BOHR() << " bohr"
          << "  t1=" << cfg::T1_AU << " a.u.  t2=" << cfg::T2_AU << " a.u.\n";
        f << "# Rows: iy=0.." << Ny_g-1 << "  Cols: ix=0.." << Nx_g-1 << "\n";
        for(int iy = 0; iy < Ny_g; iy++){
            for(int ix = 0; ix < Nx_g; ix++){
                if(ix) f << " ";
                f << std::scientific << std::setprecision(8) << leed[iy][ix];
            }
            f << "\n";
        }
        std::cout << "  Wrote results/leed_pattern/leed_screen.txt\n";
    }

    // Simulation summary
    auto wall_end  = std::chrono::steady_clock::now();
    double wall_s  = std::chrono::duration<double>(wall_end - wall_start).count();
    {
        std::ofstream f("results/sim_summary.txt");
        f << "# run_005 simulation summary\n";
        f << std::fixed << std::setprecision(10);
        f << "GS_energy_Ha         " << E_gs                       << "\n";
        f << "SCF_steps            " << gs.total_iter              << "\n";
        f << "WP_d_ang             " << cfg::WP_D_ANG              << "\n";
        f << "WP_d_bohr            " << cfg::WP_D_BOHR             << "\n";
        f << "WP_D_ang             " << cfg::WP_D_IMPACT_ANG       << "\n";
        f << "WP_D_bohr            " << cfg::WP_D_IMPACT_BOHR      << "\n";
        f << "WP_bz_bohr           " << cfg::WP_BZ()               << "\n";
        f << "WP_k0_bohr_inv       " << cfg::wp_k0()               << "\n";
        f << "WP_norm              " << wp_norm                    << "\n";
        f << "Lz_ang               " << cfg::LZ_ANG                << "\n";
        f << "z_flake_bohr         " << cfg::Z_FLAKE_BOHR()        << "\n";
        f << "z_screen_bohr        " << cfg::Z_SCREEN_BOHR()       << "\n";
        f << "t1_au                " << cfg::T1_AU                 << "\n";
        f << "t2_au                " << cfg::T2_AU                 << "\n";
        f << "dt_au                " << cfg::DT_AU                 << "\n";
        f << "n_steps              " << cfg::N_STEPS               << "\n";
        f << "wall_time_sec        " << wall_s                     << "\n";
        std::cout << "  Wrote results/sim_summary.txt\n";
    }

    std::cout << "  Wall time: " << std::fixed << std::setprecision(1)
              << wall_s << " s  (" << wall_s/3600.0 << " h)\n";
    std::cout << "\n=== run_005 complete. ===\n";
    return 0;
}
