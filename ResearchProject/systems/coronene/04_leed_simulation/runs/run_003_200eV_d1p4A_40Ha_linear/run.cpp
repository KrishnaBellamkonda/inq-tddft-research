// ============================================================================
// run_003 — Coronene TDDFT LEED simulation
//
// Reproduces Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014):
//   200 eV electron wavepacket scattering off coronene C24H12.
//
// Fixes vs run_002:
//   (a) GPU sync:  cudaDeviceSynchronize() at top of callback — prevents stale
//       orbital reads that corrupted WP data in run_002.
//   (b) GPU write: inject_wp uses gpu::run kernel (not CPU loop) to write WP
//       orbital into UVM-managed device memory.
//   (c) Coord fix: WP/centroid computation uses (ix+x0)*dx_sp instead of
//       point_operator capture (parallel::partition GPU-capture bug).
//   (d) Tighter SCF: 1e-6 Ha, linear mixing α=0.05.
//   (e) Full observables: GS orbitals, WP orbital, 3D density, 57×57 overlap
//       matrix, momentum — all saved every SNAPSHOT_INTERVAL=10 steps.
//   (f) Correct WP centre: (Lx/2, Ly/2, Lz/2 + D) — cell corner convention.
//
// Output layout (results/):
//   grid/                   grid coordinates (once)
//   gs_orbitals/            all 57 GS KS orbitals (once, before WP injection)
//   energy/                 gs_energy.txt + energy_vs_time.csv (every step)
//   momentum/               momentum_vs_time.csv  (every step)
//   wp_orbital/             WP orbital 3D text     (every 10 steps)
//   density/                3D total density       (every 10 steps)
//   density_snapshots/      2D slice at z_flake    (every 10 steps)
//   density_mid_snapshots/  2D slice at z_mid      (every 10 steps)
//   density_obs_snapshots/  2D slice at z_obs      (every 10 steps)
//   wp_trajectory/          1D z-profile CSV       (every 10 steps)
//   ks_overlaps/            diagonal S_ii CSV      (every 10 steps)
//   overlap_matrix/         full 57×57 S_ij        (every 10 steps)
//   leed_pattern/           I(x,y) = ∫n dt at z_obs over [t1, t2]
//   sim_summary.txt
//
// Source: Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014)
// All values in atomic units unless annotated.
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

    // ── Helpers ───────────────────────────────────────────────────────────────
    auto pad4 = [](int n){ std::ostringstream s; s << std::setfill('0') << std::setw(4) << n; return s.str(); };
    auto pad6 = [](int n){ std::ostringstream s; s << std::setfill('0') << std::setw(6) << n; return s.str(); };

    // ── STEP 1: Load coronene geometry ────────────────────────────────────────
    auto cell = cfg::make_cell();
    auto ions = systems::ions::parse(cfg::CORONENE_XYZ, cell);

    std::cout << "\n=== run_003: Coronene TDDFT LEED ===\n";
    std::cout << "  Cell   : " << cfg::LX_BOHR << " × " << cfg::LY_BOHR
              << " × " << cfg::LZ_BOHR << " bohr  [finite, corner at origin]\n";
    std::cout << "  E_cut  : " << cfg::ECUT_HA  << " Ha\n";
    std::cout << "  Atoms  : " << ions.size()   << "  (expect 36)\n";

    // ── STEP 2: Ground state SCF ──────────────────────────────────────────────
    auto electrons = systems::electrons(ions,
        options::electrons{}.cutoff(cfg::ECUT_HA * 1.0_Ha)
                            .extra_states(cfg::EXTRA_STATES));

    std::cout << "  States : " << electrons.kpin()[0].spinor_set_size()
              << "  (54 occ + " << cfg::EXTRA_STATES << " extra)\n\n";
    std::cout << "=== STEP 2: Ground state SCF ===\n";
    std::cout << "  tol=" << cfg::SCF_TOL << " Ha  mixing=" << cfg::SCF_MIXING
              << "  max_steps=" << cfg::SCF_MAX_STEPS << "\n";

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
    std::cout << "  GS energy : " << E_gs << " Ha"
              << "  SCF steps : " << gs.total_iter << "\n";

    if(gs.total_iter >= cfg::SCF_MAX_STEPS)
        std::cout << "  WARNING: SCF reached max_steps — may not be fully converged.\n";

    // Save ground state energy
    {
        std::ofstream f("results/energy/gs_energy.txt");
        f << "# Ground state LDA energy — run_003\n";
        f << "# Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014)\n";
        f << std::fixed << std::setprecision(10)
          << "GS_energy_Ha  " << E_gs << "\n"
          << "SCF_steps     " << gs.total_iter << "\n"
          << "SCF_tol_Ha    " << cfg::SCF_TOL << "\n"
          << "E_cut_Ha      " << cfg::ECUT_HA << "\n";
    }
    std::cout << "  Wrote results/energy/gs_energy.txt\n";

    // Save grid coordinates (once, before WP injection)
    leed_utils::save_grid_coords(electrons, "results/grid");
    std::cout << "  Wrote results/grid/\n";

    // ── STEP 3: Save GS orbitals (before WP injection) ────────────────────────
    std::cout << "\n=== STEP 3: Saving GS orbitals ===\n";
    {
        int n_st = electrons.kpin()[0].set_part().local_size();
        for(int ist = 0; ist < n_st; ist++){
            std::string dir = "results/gs_orbitals/orbital_" + pad4(ist);
            leed_utils::mkdir_p(dir);
            leed_utils::save_orbital_3d(electrons, ist, 0.0, -1, dir + "/orbital.txt");
        }
        std::cout << "  Saved " << n_st << " GS orbitals to results/gs_orbitals/\n";
    }

    // Keep a copy of the GS electrons for overlap matrix computation.
    // gs_electrons has occupations of the ground state (all extra states = 0).
    auto gs_electrons = electrons;

    // ── STEP 4: Inject WP and set occupation ─────────────────────────────────
    std::cout << "\n=== STEP 4: WP injection ===\n";

    const double k0    = cfg::wp_k0();
    const int    ist_wp = electrons.kpin()[0].set_part().local_size() - 1;

    leed_utils::inject_wp(electrons,
        cfg::WP_BX(), cfg::WP_BY(), cfg::WP_BZ(),
        0.0, 0.0, -k0);

    electrons.occupations()[0][ist_wp] = cfg::WP_OCCUPATION;

    auto [wp_norm, wp_ke] = leed_utils::validate_wp(electrons);
    std::cout << "  WP centre : (" << cfg::WP_BX() << ", " << cfg::WP_BY()
              << ", " << cfg::WP_BZ() << ") bohr\n";
    std::cout << "  k₀        : " << k0 << " bohr⁻¹  (kz = -k₀)\n";
    std::cout << "  E_kin     : " << cfg::WP_EKIN_EV << " eV  ("
              << cfg::WP_EKIN_HA << " Ha)\n";
    std::cout << "  ⟨ψ|ψ⟩    : " << wp_norm;
    if(std::abs(wp_norm - 1.0) > 0.05)
        std::cout << "  WARNING: norm deviates > 5% — consider increasing E_cut";
    std::cout << "\n";

    // Save WP orbital at step 0 (before propagation)
    {
        std::string dir = "results/wp_orbital/step_000000/kpt_0/orbital_" + pad4(ist_wp);
        leed_utils::mkdir_p(dir);
        leed_utils::save_orbital_3d(electrons, ist_wp, 0.0, 0, dir + "/orbital.txt");
    }
    // Save 3D density at step 0
    leed_utils::save_density_3d(electrons, 0.0, 0, "results/density/density_t000000.txt");
    std::cout << "  Saved step-0 WP orbital and density.\n";

    // ── STEP 5: TDDFT propagation ─────────────────────────────────────────────
    std::cout << "\n=== STEP 5: TDDFT propagation ===\n";
    std::cout << "  Δt          : " << cfg::DT_AU << " a.u.  (" << cfg::DT_FS << " fs)\n";
    std::cout << "  t₁ (arrival): " << cfg::T1_AU << " a.u.  (" << cfg::T1_FS << " fs)\n";
    std::cout << "  t₂ (end)    : " << cfg::T2_AU << " a.u.  (" << cfg::T2_FS << " fs)\n";
    std::cout << "  N_steps     : " << cfg::N_STEPS << "\n";
    std::cout << "  Snapshot cadence: every " << cfg::SNAPSHOT_INTERVAL
              << " steps (~" << cfg::MAX_SNAPSHOTS << " saves)\n\n";

    // Allocate outputs
    auto const & basis = electrons.states_basis();
    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    int Nz_g = basis.sizes()[2];

    // LEED accumulator: I(x,y) = ∫_{t1}^{t2} n(x,y,z_obs,t) dt
    std::vector<std::vector<double>> leed_accum(Ny_g, std::vector<double>(Nx_g, 0.0));

    // Open persistent output files (appended to throughout callback)
    std::ofstream energy_csv("results/energy/energy_vs_time.csv");
    energy_csv << "# step,t_au,E_total_Ha\n";

    std::ofstream momentum_csv("results/momentum/momentum_vs_time.csv");
    momentum_csv << "# step,t_au,Jx_au,Jy_au,Jz_au\n";

    std::ofstream zprofile_csv("results/wp_trajectory/density_z_profile_vs_time.csv");
    zprofile_csv << "# step,t_au";
    for(int iz = 0; iz < Nz_g; iz++) zprofile_csv << ",n_iz" << iz;
    zprofile_csv << "\n";

    std::ofstream proj_occ_csv("results/ks_overlaps/projected_occ_vs_time.csv");
    proj_occ_csv << "# step,t_au";
    for(int i = 0; i < electrons.kpin()[0].set_part().local_size(); i++)
        proj_occ_csv << ",|S_" << i << i << "|^2";
    proj_occ_csv << "\n";

    std::ofstream overlap_f("results/overlap_matrix/overlap_matrix.txt");
    overlap_f << "# Overlap matrix S_ij(t) = <phi_i^GS | phi_j(t)>\n";
    overlap_f << "# run_003 coronene TDDFT LEED — every " << cfg::SNAPSHOT_INTERVAL << " steps\n";
    overlap_f << "# Ref: Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014)\n\n";

    double t1 = cfg::T1_AU;
    double t2 = cfg::T2_AU;

    // Propagation callback
    auto obs_callback = [&](auto && obs){
        GPU_SYNC();  // ensure all GPU kernels complete before any CPU read

        double t    = obs.time();
        int    iter = obs.iter();
        auto   e    = obs.energy();
        auto   J    = obs.current();  // vector3<double>: Jx, Jy, Jz (a.u.)

        // ── Every step: energy and momentum (lightweight) ──────────────────
        energy_csv << iter << "," << std::fixed << std::setprecision(8)
                   << t << "," << e.total() << "\n";

        momentum_csv << iter << "," << std::fixed << std::setprecision(8)
                     << t << ","
                     << std::scientific << std::setprecision(8)
                     << J[0] << "," << J[1] << "," << J[2] << "\n";

        // ── LEED accumulation: t1 ≤ t ≤ t2 (every step, at z_obs) ────────
        if(t >= t1 && t <= t2){
            auto slice_obs = leed_utils::extract_density_slice(electrons, cfg::Z_OBS_BOHR());
            for(int iy = 0; iy < Ny_g; iy++)
                for(int ix = 0; ix < Nx_g; ix++)
                    leed_accum[iy][ix] += slice_obs[iy][ix] * cfg::DT_AU;
        }

        // ── Every SNAPSHOT_INTERVAL steps: all heavy saves ─────────────────
        if(iter % cfg::SNAPSHOT_INTERVAL != 0) return;

        std::cout << "  step=" << iter << "  t=" << std::fixed << std::setprecision(4)
                  << t << " a.u.  E=" << std::setprecision(6) << e.total()
                  << " Ha  Jz=" << std::scientific << std::setprecision(3) << J[2] << "\n";

        // WP orbital 3D save
        {
            std::string dir = "results/wp_orbital/step_" + pad6(iter)
                              + "/kpt_0/orbital_" + pad4(ist_wp);
            leed_utils::mkdir_p(dir);
            leed_utils::save_orbital_3d(electrons, ist_wp, t, iter, dir + "/orbital.txt");
        }

        // 3D total density save
        leed_utils::save_density_3d(electrons, t, iter,
            "results/density/density_t" + pad6(iter) + ".txt");

        // 2D density slices at z_flake, z_mid, z_obs
        auto slice_flake = leed_utils::extract_density_slice(electrons, cfg::Z_FLAKE_BOHR());
        auto slice_mid   = leed_utils::extract_density_slice(electrons, cfg::Z_MID_BOHR());
        auto slice_obs   = leed_utils::extract_density_slice(electrons, cfg::Z_OBS_BOHR());

        leed_utils::save_density_slice(slice_flake, t, cfg::Z_FLAKE_BOHR(),
            "results/density_snapshots/snapshot_t" + pad6(iter) + ".txt");
        leed_utils::save_density_slice(slice_mid, t, cfg::Z_MID_BOHR(),
            "results/density_mid_snapshots/snapshot_t" + pad6(iter) + ".txt");
        leed_utils::save_density_slice(slice_obs, t, cfg::Z_OBS_BOHR(),
            "results/density_obs_snapshots/snapshot_t" + pad6(iter) + ".txt");

        // 1D z-profile at cell centre (ix=Nx/2, iy=Ny/2)
        auto zprof = leed_utils::extract_z_profile(electrons);
        zprofile_csv << iter << "," << std::fixed << std::setprecision(8) << t;
        for(double v : zprof) zprofile_csv << "," << std::scientific << std::setprecision(6) << v;
        zprofile_csv << "\n";

        // Full 57×57 overlap matrix S_ij and diagonal (projected occupation)
        auto S = leed_utils::compute_overlap_matrix(gs_electrons, electrons);
        leed_utils::save_overlap_matrix(S, t, iter, overlap_f);

        int n_st = static_cast<int>(S.size());
        proj_occ_csv << iter << "," << std::fixed << std::setprecision(8) << t;
        for(int i = 0; i < n_st; i++){
            double diag = S[i][i].real()*S[i][i].real() + S[i][i].imag()*S[i][i].imag();
            proj_occ_csv << "," << std::scientific << std::setprecision(6) << diag;
        }
        proj_occ_csv << "\n";
    };

    // Run propagation
    real_time::propagate(ions, electrons,
        obs_callback,
        options::theory{}.lda(),
        options::real_time{}
            .dt(cfg::DT_AU * 1.0_atomictime)
            .num_steps(cfg::N_STEPS)
            .observables_current());

    // ── STEP 6: Write LEED pattern and summary ────────────────────────────────
    std::cout << "\n=== STEP 6: Writing final outputs ===\n";

    // LEED pattern: I(x,y) = ∫_{t1}^{t2} n(x,y,z_obs,t) dt
    {
        std::ofstream f("results/leed_pattern/leed_pattern.txt");
        f << "# LEED pattern I(x,y) = integral_{t1}^{t2} n(x,y,z_obs,t) dt\n";
        f << "# Ref: Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014), Fig. 2\n";
        f << "# z_obs=" << cfg::Z_OBS_BOHR() << " bohr"
          << "  t1=" << t1 << " a.u.  t2=" << t2 << " a.u.\n";
        f << "# Rows: iy=0.." << Ny_g-1 << "  Cols: ix=0.." << Nx_g-1 << "\n";
        f << std::scientific << std::setprecision(8);
        for(auto const & row : leed_accum){
            for(size_t ix = 0; ix < row.size(); ix++){
                f << row[ix];
                if(ix+1 < row.size()) f << " ";
            }
            f << "\n";
        }
        std::cout << "  Wrote results/leed_pattern/leed_pattern.txt\n";
    }

    // Simulation summary
    {
        auto wall_end = std::chrono::steady_clock::now();
        double wall_sec = std::chrono::duration<double>(wall_end - wall_start).count();

        std::ofstream f("results/sim_summary.txt");
        f << "# Coronene WP scattering — TDDFT LEED simulation (run_003)\n";
        f << "# Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014)\n";
        f << std::fixed << std::setprecision(10);
        f << "GS_energy_Ha       " << E_gs                    << "\n";
        f << "SCF_steps          " << gs.total_iter            << "\n";
        f << "SCF_tol_Ha         " << cfg::SCF_TOL             << "\n";
        f << "SCF_mixing         " << cfg::SCF_MIXING          << "\n";
        f << "E_cut_Ha           " << cfg::ECUT_HA             << "\n";
        f << "WP_d_bohr          " << cfg::WP_D_BOHR           << "\n";
        f << "WP_D_bohr          " << cfg::WP_D_IMPACT_BOHR    << "\n";
        f << "WP_bx_bohr         " << cfg::WP_BX()             << "\n";
        f << "WP_by_bohr         " << cfg::WP_BY()             << "\n";
        f << "WP_bz_bohr         " << cfg::WP_BZ()             << "\n";
        f << "WP_k0_bohr_inv     " << cfg::wp_k0()             << "\n";
        f << "WP_Ekin_eV         " << cfg::WP_EKIN_EV          << "\n";
        f << "WP_norm            " << wp_norm                   << "\n";
        f << "WP_occ             " << cfg::WP_OCCUPATION        << "\n";
        f << "dt_au              " << cfg::DT_AU                << "\n";
        f << "t1_au              " << cfg::T1_AU                << "\n";
        f << "t2_au              " << cfg::T2_AU                << "\n";
        f << "n_steps            " << cfg::N_STEPS              << "\n";
        f << "snapshot_interval  " << cfg::SNAPSHOT_INTERVAL    << "\n";
        f << "n_snapshots        " << cfg::MAX_SNAPSHOTS        << "\n";
        f << "wall_time_sec      " << wall_sec                  << "\n";

        std::cout << "  Wrote results/sim_summary.txt\n";
        std::cout << "  Wall time: " << wall_sec << " s\n";
    }

    std::cout << "\n=== Simulation complete. ===\n\n";

    // Post-run validation hints
    std::cout << "Post-run checks:\n";
    std::cout << "  1. Energy drift: diff first/last line of results/energy/energy_vs_time.csv\n";
    std::cout << "  2. WP trajectory: plot density_z_profile_vs_time.csv iz-centroid vs step\n";
    std::cout << "  3. Momentum: results/momentum/momentum_vs_time.csv Jz[step=1] ≈ ±3.834\n";
    std::cout << "  4. Overlap diagonal: |S_ii(step=0)|² ≈ 1.0 in ks_overlaps/\n";
    std::cout << "  5. LEED: python3 analysis.py\n\n";

    return 0;
}
