// ============================================================================
// run_002 — Coronene TDDFT LEED simulation
//
// Key changes from run_001:
//   - E_cut = 54 Ha (paper grid spacing 0.16 Å)
//   - WP width d = 1.4 Å  (illuminates central ring: 35.7% at r=1.421 Å)
//   - Geometry centred at (Lx/2, Ly/2, Lz/2) — no 4-copy artefact
//   - GS electrons copied for projected_occupation analysis
//   - Energy components saved at every step
//   - Projected KS occupation vs GS saved every OVERLAP_INTERVAL steps
//   - WP-only orbital slices saved at snapshot intervals
//   - 1D z-profile at cell centre (WP trajectory density) at snapshots
//
// Reference: Tsubonoya, Hu, Watanabe PRB 90, 035416 (2014)
// ============================================================================

#include "config.hpp"
#include "utils.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>

int main(){
    using namespace inq;
    using namespace inq::magnitude;

    // ── STEP 1: Load centred coronene geometry ────────────────────────────────
    auto cell = cfg::make_cell();
    auto ions = systems::ions::parse(cfg::CORONENE_XYZ, cell);

    // ── STEP 2: Ground state SCF ──────────────────────────────────────────────
    auto electrons = systems::electrons(ions,
        options::electrons{}.cutoff(cfg::ECUT_HA_LEED * 1.0_Ha)
                            .extra_states(cfg::EXTRA_STATES));

    if(electrons.root()){
        std::cout << "\n=== STEP 2: Ground state SCF ===\n";
        std::cout << "  E_cut  : " << cfg::ECUT_HA_LEED << " Ha  (h=0.186 Ang)\n";
        std::cout << "  Atoms  : " << ions.size()         << "  (expect 36)\n";
        std::cout << "  States : " << electrons.kpin()[0].spinor_set_size()
                  << "  (54 occ + " << cfg::EXTRA_STATES << " extra)\n";
    }

    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}.energy_tolerance(cfg::SCF_TOL * 1.0_Ha).max_steps(300)
                               .mixing(0.1).mixing_ndim(8));

    if(electrons.root()){
        std::cout << "  GS energy : " << gs.energy.total()    << " Ha\n";
        std::cout << "  Kinetic   : " << gs.energy.kinetic()  << " Ha\n";
        std::cout << "  Hartree   : " << gs.energy.hartree()  << " Ha\n";
        std::cout << "  XC        : " << gs.energy.xc()       << " Ha\n";
        std::cout << "  External  : " << gs.energy.external() << " Ha\n";
        std::cout << "  Non-local : " << gs.energy.non_local()<< " Ha\n";
        std::cout << "  Ion-ion   : " << gs.energy.ion()      << " Ha\n";
        std::cout << "  SCF steps : " << gs.total_iter        << "\n";

        std::ofstream gsf("results/energy/gs_energy.txt");
        gsf << "# Ground state energy components (Ha)\n";
        gsf << "E_total    " << gs.energy.total()    << "\n";
        gsf << "E_kinetic  " << gs.energy.kinetic()  << "\n";
        gsf << "E_hartree  " << gs.energy.hartree()  << "\n";
        gsf << "E_xc       " << gs.energy.xc()       << "\n";
        gsf << "E_external " << gs.energy.external() << "\n";
        gsf << "E_nonlocal " << gs.energy.non_local()<< "\n";
        gsf << "E_ion      " << gs.energy.ion()      << "\n";
        gsf << "SCF_steps  " << gs.total_iter        << "\n";
    }

    // ── STEP 3: Copy GS electrons (before WP injection) ──────────────────────
    // Used to compute projected occupations during propagation.
    // GPU memory note: this doubles the orbital storage (~2.4 GB at 54 Ha).
    auto gs_electrons = electrons;
    int n_gs_states = gs_electrons.kpin()[0].set_size();  // 57 (54 + EXTRA_STATES)

    // ── STEP 4: Inject WP into last extra-state orbital ───────────────────────
    {
        const double k0 = cfg::wp_k0();
        leed_utils::inject_wp(electrons,
            cfg::WP_BX(), cfg::WP_BY(), cfg::WP_BZ(),
            0.0, 0.0, -k0);
    }

    int ist_wp = electrons.kpin()[0].set_part().local_size() - 1;
    electrons.occupations()[0][ist_wp] = cfg::WP_OCCUPATION;

    // ── STEP 5: Validate WP ───────────────────────────────────────────────────
    if(electrons.root()){
        auto [wp_norm, wp_ke] = leed_utils::validate_wp(electrons);
        std::cout << "\n=== STEP 3-5: WP injection ===\n";
        std::cout << "  d   = " << cfg::WP_D_BOHR  << " bohr  (" << cfg::WP_D_ANG << " Ang)\n";
        std::cout << "  D   = " << cfg::WP_D_IMPACT_BOHR << " bohr\n";
        std::cout << "  bz  = " << cfg::WP_BZ() << " bohr  (Lz/2 + D)\n";
        std::cout << "  k0  = " << cfg::wp_k0() << " bohr^-1\n";
        std::cout << "  Ekin= " << cfg::WP_EKIN_EV << " eV\n";
        std::cout << "  <psi|psi> = " << wp_norm << "  (expect 1.0)\n";
        if(std::abs(wp_norm - 1.0) > 0.05)
            std::cout << "  WARNING: WP norm deviates > 5%\n";
    }

    // ── STEP 6: Open output files ─────────────────────────────────────────────
    std::ofstream energy_csv, overlap_csv, ztraj_csv;
    if(electrons.root()){
        energy_csv.open("results/energy/energy_vs_time.csv");
        energy_csv << "# step, t_au, E_total_Ha, E_kinetic_Ha, E_hartree_Ha,"
                   << " E_xc_Ha, E_external_Ha, E_nonlocal_Ha, E_ion_Ha\n";

        overlap_csv.open("results/ks_overlaps/projected_occ_vs_time.csv");
        overlap_csv << "# step, t_au";
        for(int i = 0; i < n_gs_states; i++) overlap_csv << ", occ_" << i;
        overlap_csv << "\n";

        ztraj_csv.open("results/wp_trajectory/density_z_profile_vs_time.csv");
        ztraj_csv << "# Electron density at cell centre (ix=Nx/2, iy=Ny/2) along z\n";
        ztraj_csv << "# step, t_au, n(z_0), n(z_1), ...\n";
    }

    // LEED accumulator
    auto const & basis = electrons.states_basis();
    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    std::vector<std::vector<double>> leed_accum(Ny_g, std::vector<double>(Nx_g, 0.0));

    const double t1 = cfg::T1_AU;
    const double t2 = cfg::T2_AU;
    int snap_id = 0;

    // ── STEP 7: TDDFT propagation ─────────────────────────────────────────────
    if(electrons.root()){
        std::cout << "\n=== STEP 7: TDDFT propagation ===\n";
        std::cout << "  dt = " << cfg::DT_AU << " a.u.\n";
        std::cout << "  t1 = " << t1 << " a.u.  t2 = " << t2 << " a.u.\n";
        std::cout << "  N  = " << cfg::N_STEPS << " steps\n\n";
    }

    auto obs_callback = [&](auto && obs){
        double t    = obs.time();
        int    iter = obs.iter();
        auto   e    = obs.energy();

        // Save energy every step
        if(electrons.root()){
            energy_csv << iter
                << ", " << std::fixed << std::setprecision(8) << t
                << ", " << e.total()    << ", " << e.kinetic()
                << ", " << e.hartree()  << ", " << e.xc()
                << ", " << e.external() << ", " << e.non_local()
                << ", " << e.ion()      << "\n";
        }

        // Save projected occupation every OVERLAP_INTERVAL steps
        if(iter % cfg::OVERLAP_INTERVAL == 0){
            auto pocc = obs.projected_occupation(gs_electrons);
            if(electrons.root()){
                overlap_csv << iter << ", " << t;
                for(int i = 0; i < n_gs_states; i++)
                    overlap_csv << ", " << pocc[0][i];
                overlap_csv << "\n";
            }
        }

        // Snapshots every SNAPSHOT_INTERVAL steps
        if(iter % cfg::SNAPSHOT_INTERVAL == 0 && snap_id < cfg::MAX_SNAPSHOTS){
            if(electrons.root()){
                // Total density at coronene plane (z_flake)
                {
                    auto sl = leed_utils::extract_density_slice(electrons, cfg::Z_FLAKE_BOHR());
                    std::ostringstream fn;
                    fn << "results/density_snapshots/snapshot_t"
                       << std::setfill('0') << std::setw(4) << snap_id << ".txt";
                    leed_utils::save_density_slice(sl, t, cfg::Z_FLAKE_BOHR(), fn.str());
                }

                // WP orbital density at coronene plane
                {
                    auto sl = leed_utils::extract_wp_slice(electrons, cfg::Z_FLAKE_BOHR());
                    std::ostringstream fn;
                    fn << "results/wp_orbital/wp_slice_t"
                       << std::setfill('0') << std::setw(4) << snap_id << ".txt";
                    leed_utils::save_density_slice(sl, t, cfg::Z_FLAKE_BOHR(), fn.str());
                }

                // 1D z-profile at cell centre
                {
                    auto zp = leed_utils::extract_z_profile(electrons);
                    ztraj_csv << iter << ", " << t;
                    for(double v : zp) ztraj_csv << ", " << v;
                    ztraj_csv << "\n";
                }

                std::cout << "  t=" << std::fixed << std::setprecision(3) << t
                          << " a.u.  snapshot " << snap_id+1 << "\n";
            }
            snap_id++;
        }

        // Accumulate LEED (t1 to t2)
        if(t >= t1 && t <= t2){
            if(electrons.root()){
                auto sl = leed_utils::extract_density_slice(electrons, cfg::Z_OBS_BOHR());
                for(int iy = 0; iy < Ny_g; iy++)
                    for(int ix = 0; ix < Nx_g; ix++)
                        leed_accum[iy][ix] += sl[iy][ix] * cfg::DT_AU;
            }
        }
    };

    real_time::propagate(ions, electrons,
        obs_callback,
        options::theory{}.lda(),
        options::real_time{}.dt(cfg::DT_AU * 1.0_atomictime).num_steps(cfg::N_STEPS));

    // ── STEP 8: Write final results ───────────────────────────────────────────
    if(electrons.root()){
        energy_csv.close();
        overlap_csv.close();
        ztraj_csv.close();

        std::ofstream lf("results/leed_pattern/leed_pattern.txt");
        lf << "# LEED I(x,y) = int_{t1}^{t2} n(x,y,z=z_obs,t) dt\n";
        lf << "# z_obs=" << cfg::Z_OBS_BOHR() << " bohr"
           << "  t1=" << t1 << " t2=" << t2 << " a.u.\n";
        for(auto const & row : leed_accum){
            for(size_t ix = 0; ix < row.size(); ix++){
                lf << std::scientific << std::setprecision(6) << row[ix];
                if(ix+1 < row.size()) lf << " ";
            }
            lf << "\n";
        }
        lf.close();

        std::ofstream sf("results/sim_summary.txt");
        sf << "# run_002 Coronene TDDFT LEED\n";
        sf << "GS_energy_Ha    " << gs.energy.total()      << "\n";
        sf << "E_cut_Ha        " << cfg::ECUT_HA_LEED       << "\n";
        sf << "WP_d_bohr       " << cfg::WP_D_BOHR          << "\n";
        sf << "WP_D_bohr       " << cfg::WP_D_IMPACT_BOHR   << "\n";
        sf << "WP_BZ_bohr      " << cfg::WP_BZ()            << "\n";
        sf << "WP_k0_bohr_inv  " << cfg::wp_k0()            << "\n";
        sf << "WP_Ekin_eV      " << cfg::WP_EKIN_EV         << "\n";
        sf << "dt_au           " << cfg::DT_AU               << "\n";
        sf << "t1_au           " << t1                       << "\n";
        sf << "t2_au           " << t2                       << "\n";
        sf << "n_steps         " << cfg::N_STEPS             << "\n";
        sf << "n_snapshots     " << snap_id                  << "\n";
        sf.close();

        std::cout << "\nWrote all results. Simulation complete.\n";
    }

    return 0;
}
