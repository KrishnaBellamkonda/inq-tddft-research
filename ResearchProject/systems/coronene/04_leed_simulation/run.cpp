// ============================================================================
// 04_leed_simulation/run.cpp
//
// Coronene C24H12 — TDDFT electron wavepacket scattering simulation.
// Replicates Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014), Fig. 1 and Fig. 2.
//
// Simulation flow
// ───────────────
//   Step 1 │ Load coronene geometry from XYZ file
//   Step 2 │ LDA ground state SCF (extra_states=3 for WP + buffer)
//   Step 3 │ Inject Gaussian WP into last extra-state orbital; set occ = 1
//   Step 4 │ Validate WP: check ⟨ψ|ψ⟩ ≈ 1, total electrons ≈ 109
//   Step 5 │ TDDFT propagation (ALDA, 4th-order Taylor, Δt = 0.02 a.u.)
//           │   — every SNAPSHOT_INTERVAL steps: save 2D density slice (Fig. 1)
//           │   — from t₁ to t₂: accumulate density on observation plane (Fig. 2)
//   Step 6 │ Write LEED pattern and density snapshots to results/
//
// Three output files:
//   results/leed_pattern.txt       — I(x,y) = ∫n(x,y,z=D,t)dt  (Fig. 2 data)
//   results/snapshot_t<NNN>.txt    — 2D density at z=0, time t   (Fig. 1 frames)
//   results/sim_summary.txt        — run parameters and timing
//
// Usage: inq-run (GPU) or inq-run --cpu (CPU)
// ============================================================================

#include "config.hpp"
#include "utils.hpp"

#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

int main(){
    using namespace inq;
    using namespace inq::magnitude;

    const bool root = true; // will be replaced by electrons.root() check after init

    // ── STEP 1: Load coronene ions ────────────────────────────────────────────
    auto cell  = cfg::make_cell();
    auto ions  = systems::ions::parse(cfg::CORONENE_XYZ, cell);

    // ── STEP 2: Ground state SCF ──────────────────────────────────────────────
    // Use the converged E_cut from 03_ecut_convergence.
    // extra_states(3): slot 54 = WP orbital; slots 55-56 = SCF buffer (occ = 0)
    auto electrons = systems::electrons(ions,
        options::electrons{}.cutoff(cfg::ECUT_HA_LEED * 1.0_Ha)
                            .extra_states(cfg::EXTRA_STATES));

    if(electrons.root()){
        std::cout << "\n=== STEP 1-2: Coronene ground state (LDA) ===\n";
        std::cout << "  Cell   : " << cfg::LX_BOHR << " × " << cfg::LY_BOHR
                  << " × " << cfg::LZ_BOHR << " bohr  [finite]\n";
        std::cout << "  E_cut  : " << cfg::ECUT_HA_LEED  << " Ha\n";
        std::cout << "  Atoms  : " << ions.size()         << "  (expect 36)\n";
        std::cout << "  States : " << electrons.kpin()[0].spinor_set_size()
                  << "  (54 occ + " << cfg::EXTRA_STATES << " extra)\n";
    }

    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(cfg::SCF_TOL * 1.0_Ha)
            .max_steps(300));

    if(electrons.root()){
        std::cout << "  GS energy : " << gs.energy.total() << " Ha"
                  << "  SCF steps : " << gs.total_iter << "\n";
    }

    // ── STEP 3: Inject WP into last extra-state orbital ───────────────────────
    //
    // WP centre: (0, 0, +D) — D = 12 bohr above the coronene plane.
    // WP momentum: k = (0, 0, −k₀) — directed toward the flake (−z direction).
    //
    // Paper Eq. 1: ψ^WP(r) = (1/(πd²))^{3/4} · exp(−|r−b|²/(2d²)) · exp(ik·r)
    {
        const double k0 = cfg::wp_k0();
        leed_utils::inject_wp(electrons,
            /* bx */ 0.0, /* by */ 0.0, /* bz */ cfg::WP_D_IMPACT_BOHR,
            /* kx */ 0.0, /* ky */ 0.0, /* kz */ -k0);
    }

    // Set WP orbital occupation = 1.0  (singly occupied: one incident electron)
    // occupations()[kpt_idx][local_state_idx]
    {
        int ist_wp = electrons.kpin()[0].set_part().local_size() - 1;
        electrons.occupations()[0][ist_wp] = cfg::WP_OCCUPATION;
    }

    // ── STEP 4: Validate WP ───────────────────────────────────────────────────
    if(electrons.root()){
        auto [wp_norm, wp_ke] = leed_utils::validate_wp(electrons);
        std::cout << "\n=== STEP 3-4: WP injection ===\n";
        std::cout << "  d      = " << cfg::WP_D_BOHR  << " bohr  ("
                  << cfg::WP_D_ANG  << " Å)\n";
        std::cout << "  D      = " << cfg::WP_D_IMPACT_BOHR << " bohr  ("
                  << cfg::WP_D_IMPACT_ANG  << " Å)\n";
        std::cout << "  k₀     = " << cfg::wp_k0() << " bohr⁻¹\n";
        std::cout << "  E_kin  = " << cfg::WP_EKIN_EV << " eV  ("
                  << cfg::WP_EKIN_HA << " Ha)\n";
        std::cout << "  ⟨ψ|ψ⟩ = " << wp_norm
                  << "  (expect 1.0 — deviation indicates discretisation error)\n";
        std::cout << "  occ(WP orbital) = " << cfg::WP_OCCUPATION << "\n\n";

        if(std::abs(wp_norm - 1.0) > 0.05){
            std::cout << "  WARNING: WP norm deviates > 5% from 1.0."
                      << " Consider increasing E_cut.\n";
        }
    }

    // ── STEP 5: TDDFT propagation ─────────────────────────────────────────────
    //
    // Propagation times (from paper):
    //   t₁ = WP arrival time ≈ D/k₀ = 3.13 a.u.
    //   t₂ = 0.25 fs = 10.33 a.u.
    //
    // Observables:
    //   - LEED pattern: I(x,y) = ∫_{t₁}^{t₂} n(x,y,z=D,t) dt
    //     accumulated on observation plane S at z = D (incident position)
    //   - Density snapshots at z=0 (flake plane) for Fig. 1 frames

    // Allocate LEED accumulator on root (2D grid at observation plane)
    auto const & basis = electrons.states_basis();
    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    std::vector<std::vector<double>> leed_accum(Ny_g, std::vector<double>(Nx_g, 0.0));

    double t1   = cfg::t1_au();     // WP arrival time
    double t2   = cfg::T2_AU;       // end time
    int snap_id = 0;
    int step    = 0;

    if(electrons.root()){
        std::cout << "=== STEP 5: TDDFT propagation ===\n";
        std::cout << "  Δt     = " << cfg::DT_AU << " a.u.  (" << cfg::DT_FS << " fs)\n";
        std::cout << "  t₁     = " << t1 << " a.u.\n";
        std::cout << "  t₂     = " << t2 << " a.u.  (" << cfg::T2_FS << " fs)\n";
        std::cout << "  steps  = " << cfg::N_STEPS << "\n\n";
    }

    // Propagation callback
    auto obs_callback = [&](auto const & obs){
        double t = obs.time();

        // ── Save density snapshot every SNAPSHOT_INTERVAL steps (Fig. 1) ──────
        if(obs.iter() % cfg::SNAPSHOT_INTERVAL == 0 && snap_id < cfg::MAX_SNAPSHOTS){
            if(electrons.root()){
                auto slice_flake = leed_utils::extract_density_slice(electrons,
                                        cfg::Z_FLAKE_BOHR);
                std::ostringstream fname;
                fname << "results/snapshot_t" << std::setfill('0')
                      << std::setw(4) << snap_id << ".txt";
                leed_utils::save_density_slice(slice_flake, t, cfg::Z_FLAKE_BOHR, fname.str());
                snap_id++;
                std::cout << "  t=" << std::fixed << std::setprecision(3) << t
                          << " a.u. — saved snapshot " << (snap_id) << "\n";
            }
        }

        // ── Accumulate on observation plane (t₁ to t₂) for LEED (Fig. 2) ─────
        if(t >= t1 && t <= t2){
            if(electrons.root()){
                auto slice_obs = leed_utils::extract_density_slice(electrons,
                                      cfg::Z_OBS_BOHR);
                for(int iy = 0; iy < Ny_g; iy++)
                    for(int ix = 0; ix < Nx_g; ix++)
                        leed_accum[iy][ix] += slice_obs[iy][ix] * cfg::DT_AU;
            }
        }

        step++;
    };

    // Run propagation
    real_time::propagate(ions, electrons,
        obs_callback,
        options::theory{}.lda(),
        options::real_time{}.dt(cfg::DT_AU * 1.0_atomictime).num_steps(cfg::N_STEPS));

    // ── STEP 6: Write results ─────────────────────────────────────────────────
    if(electrons.root()){

        // Write LEED pattern (Fig. 2 data)
        std::ofstream leed_file("results/leed_pattern.txt");
        leed_file << "# LEED pattern I(x,y) = integral_{t1}^{t2} n(x,y,z=D,t) dt\n";
        leed_file << "# z_obs=" << cfg::Z_OBS_BOHR << " bohr"
                  << "  t1=" << t1 << " a.u."
                  << "  t2=" << t2 << " a.u.\n";
        leed_file << "# Rows: iy = 0.." << Ny_g-1
                  << "  Cols: ix = 0.." << Nx_g-1 << "\n";
        for(auto const & row : leed_accum){
            for(size_t ix = 0; ix < row.size(); ix++){
                leed_file << std::scientific << std::setprecision(6) << row[ix];
                if(ix+1 < row.size()) leed_file << " ";
            }
            leed_file << "\n";
        }
        leed_file.close();
        std::cout << "\nWrote results/leed_pattern.txt\n";

        // Write summary
        std::ofstream sum("results/sim_summary.txt");
        sum << "# Coronene WP scattering — TDDFT LEED simulation\n";
        sum << "# Tsubonoya, Hu, Watanabe PRB 90, 035416 (2014)\n";
        sum << "GS_energy_Ha     " << gs.energy.total()     << "\n";
        sum << "E_cut_Ha         " << cfg::ECUT_HA_LEED     << "\n";
        sum << "WP_d_bohr        " << cfg::WP_D_BOHR        << "\n";
        sum << "WP_D_bohr        " << cfg::WP_D_IMPACT_BOHR << "\n";
        sum << "WP_k0_bohr_inv   " << cfg::wp_k0()          << "\n";
        sum << "WP_Ekin_eV       " << cfg::WP_EKIN_EV       << "\n";
        sum << "WP_occ           " << cfg::WP_OCCUPATION    << "\n";
        sum << "dt_au            " << cfg::DT_AU             << "\n";
        sum << "t1_au            " << t1                     << "\n";
        sum << "t2_au            " << t2                     << "\n";
        sum << "n_steps          " << cfg::N_STEPS           << "\n";
        sum << "n_snapshots      " << snap_id                << "\n";
        sum.close();
        std::cout << "Wrote results/sim_summary.txt\n";
        std::cout << "\nSimulation complete.\n";
    }

    return 0;
}
