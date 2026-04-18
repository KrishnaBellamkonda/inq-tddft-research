// 01_ground_state/run.cpp — Jellium ground-state benchmark
//
// PHYSICAL SYSTEM
//   N = 40 electrons in a cubic periodic cell of side L = 40 a₀.
//   Wigner-Seitz radius r_s = 7.26 a₀ (low-density metallic regime).
//   No ionic cores — INQ treats this as uniform jellium:
//     - The positive background exactly neutralises the electrons.
//     - The Kohn-Sham potential is spatially constant (Hartree + external cancel).
//     - KS eigenvalues are exactly ε_i = k_i²/2 + V_xc(n₀).
//
// NOTE ON GRID SPACING
//   SPACING is set to 0.50 bohr as a starting point.
//   Run 02_ground_state_convergence/ first to find the optimal value,
//   then update SPACING and re-run this benchmark.
//
// HOW JELLIUM IS SET UP IN INQ
//   A periodic cell is created but NO ions are inserted.
//   extra_electrons(N) declares N electrons with a compensating uniform
//   positive background — the standard PBC neutrality condition.
//
// OUTPUTS
//   - Energy decomposition and validation tests (stdout)
//   - results/orbitals/grid_slice.txt             — 2D grid coordinates (z=L/2 slice)
//   - results/orbitals/orbital_N_n2_M_real.txt    — Re[ψ_k] for shell M
//   - results/orbitals/orbital_N_n2_M_imag.txt    — Im[ψ_k] for shell M
//   - results/eigenvalues.txt          — KS eigenvalues for XC offset verification
//
// BUILD AND RUN
//   inq-run          (GPU)
//   inq-run --cpu    (CPU)

#include <inq/inq.hpp>
#include "jellium_utils.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <string>
#include <vector>

int main() {
    using namespace inq;
    using namespace inq::magnitude;
    namespace fs = std::filesystem;

    // ── System parameters ────────────────────────────────────────────────────

    constexpr int    N_ELECTRONS = 40;
    constexpr double L_BOHR     = 40.0;   // cell side in bohr; r_s ≈ 7.26 a₀
    constexpr double SPACING     = 0.50;  // bohr — update after running convergence test
                                           // Current E_cut ≈ π²/(2×0.50²) ≈ 19.7 Ha
    // Fermi smearing at 100 K.
    // The |n|²=3 shell is partially filled at T=0 for N=40, so smearing
    // is needed for SCF convergence regardless of the density regime.
    constexpr double SMEAR_EV   = 0.00862;
    constexpr double HA_TO_EV   = 27.211386245988;

    // ── Derived analytical quantities ────────────────────────────────────────

    double rs  = wigner_seitz_radius(N_ELECTRONS, L_BOHR);
    double n0  = mean_density(N_ELECTRONS, L_BOHR);
    double EF  = fermi_energy(N_ELECTRONS, L_BOHR);
    double wp  = plasmon_frequency(N_ELECTRONS, L_BOHR);

    // ── Print run header ─────────────────────────────────────────────────────

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n========================================\n";
    std::cout << " Jellium ground-state benchmark\n";
    std::cout << " N = " << N_ELECTRONS << " electrons,  L = " << L_BOHR << " Bohr\n";
    std::cout << " r_s = " << rs << " a₀,  n₀ = " << n0 << " bohr⁻³\n";
    std::cout << " E_F (free e⁻) = " << EF * HA_TO_EV << " eV\n";
    std::cout << " ω_p (Drude)   = " << wp * HA_TO_EV << " eV\n";
    std::cout << " Grid spacing  = " << SPACING << " bohr"
              << "  (E_cut ≈ " << (M_PI*M_PI)/(2.0*SPACING*SPACING) << " Ha)\n";
    std::cout << "========================================\n\n";

    // ── 1. Cell — cubic, periodic ────────────────────────────────────────────
    // No ions.insert() calls → INQ adds a uniform positive background → jellium.

    systems::ions ions(systems::cell::cubic(L_BOHR * 1.0_b).periodic());

    // ── 2. Electrons ─────────────────────────────────────────────────────────
    // extra_electrons(N): N electrons + uniform compensating positive background.
    // extra_states(8)  : 8 empty states above E_F for smearing + shell inspection.
    //   Total KS states = ceil(40/2) + 8 = 28.
    //
    // Gamma-point only: jellium is translation-invariant; the discrete shell
    // structure |n|²=0,1,2,3,... arises from PBC alone.

    systems::electrons electrons(ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(N_ELECTRONS)
            .extra_states(8)
            .temperature(SMEAR_EV * 1.0_eV),
        input::kpoints::gamma()
    );

    // ── 3. Ground-state SCF ───────────────────────────────────────────────────
    // LDA (Perdew-Zunger parametrisation of Ceperley-Alder data).
    // For the uniform electron gas LDA is exact by construction (it is fitted
    // to HEG data). Tight tolerance (1e-8 Ha) for meaningful TDDFT restarts.

    std::cout << "Running SCF...\n\n";
    ground_state::initial_guess(ions, electrons);

    auto gs = ground_state::calculate(ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}.energy_tolerance(1e-8_Ha)
    );

    // ── 4. Print all energy components ───────────────────────────────────────

    std::cout << "\n=== Energy decomposition (Ha) ===\n";
    std::cout << "  Total           : " << gs.energy.total()     << "\n";
    std::cout << "  Kinetic (T_s)   : " << gs.energy.kinetic()   << "\n";
    std::cout << "  Hartree         : " << gs.energy.hartree()   << "\n";
    std::cout << "  XC              : " << gs.energy.xc()        << "\n";
    std::cout << "  n·V_xc          : " << gs.energy.nvxc()      << "\n";
    std::cout << "  External        : " << gs.energy.external()  << "\n";
    std::cout << "  Non-local PP    : " << gs.energy.non_local() << "\n";
    std::cout << "  Ion-ion         : " << gs.energy.ion()       << "\n";
    std::cout << "  SCF iterations  : " << gs.total_iter         << "\n";
    std::cout << "\n  Total           : " << gs.energy.total() * HA_TO_EV << " eV\n";

    // ── 5. KS eigenvalue shell structure ─────────────────────────────────────
    // For jellium: ε_i = k_i²/2 + V_xc(n₀).
    // Print predicted values alongside V_xc.

    double vxc_pred = vxc_pz81(rs);
    std::cout << "\n=== Free-electron shell structure (predicted) ===\n";
    std::cout << "  V_xc (PZ81)   = " << vxc_pred * HA_TO_EV << " eV\n";
    std::cout << "  Shell |n|²    ε_k (eV)    ε_k + V_xc (eV)    degen\n";

    auto shells = free_electron_shells(L_BOHR, 5);
    for (auto const & sh : shells) {
        std::cout << "        " << std::setw(3) << sh.n2
                  << "      " << std::setw(8) << sh.energy_Ha * HA_TO_EV
                  << "      " << std::setw(8)
                  << (sh.energy_Ha + vxc_pred) * HA_TO_EV
                  << "          " << sh.degeneracy << "\n";
    }

    // ── 6. Rigorous validation tests ─────────────────────────────────────────

    std::cout << "\n=== Validation tests ===\n";
    bool all_passed = true;

    auto check = [&](std::string const & name, double value, double tolerance) {
        bool ok = std::fabs(value) < tolerance;
        std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << name
                  << " = " << std::setw(14) << value
                  << "  (tol = " << tolerance << ")\n";
        if (!ok) all_passed = false;
        return ok;
    };

    // Test 1: Hartree energy ≈ 0
    // Uniform ρ → all G≠0 Fourier components vanish → E_H = (Ω/2)Σ_{G≠0}|ρ_G|²4π/G² = 0.
    // Residual ~10⁻⁴ Ha comes from finite-grid numerical noise in the FFT.
    check("E_Hartree  (Ha)", gs.energy.hartree(), 1.0e-4);

    // Test 2: External energy = 0 exactly (no ionic pseudopotentials)
    check("E_external (Ha)", gs.energy.external(), 1.0e-8);

    // Test 3: Non-local PP energy = 0 (no Kleinman-Bylander projectors)
    check("E_non_local(Ha)", gs.energy.non_local(), 1.0e-8);

    // Test 4: Ion-ion (Ewald) energy = 0 (no nuclei → sum over empty set)
    check("E_ion-ion  (Ha)", gs.energy.ion(), 1.0e-8);

    // Test 5: T_s matches free-electron shell sum (KS orbitals = plane waves)
    double T_s_predicted = kinetic_energy_shells(N_ELECTRONS, L_BOHR);
    double T_s_numerical = gs.energy.kinetic();
    std::cout << "\n  T_s (numerical) : " << T_s_numerical << " Ha\n";
    std::cout << "  T_s (shells)    : " << T_s_predicted  << " Ha  (T=0 discrete sum)\n";
    check("ΔT_s       (Ha)", T_s_numerical - T_s_predicted, 0.5);

    // Test 6: E_total ≈ T_s + N·ε_xc  (all other terms ≈ 0 for jellium)
    double E_total_predicted = predicted_total_energy(N_ELECTRONS, L_BOHR);
    double E_total_numerical = gs.energy.total();
    std::cout << "\n  E_total (INQ)   : " << E_total_numerical << " Ha\n";
    std::cout << "  E_total (HEG)   : " << E_total_predicted  << " Ha  (T_s + N·ε_xc)\n";
    std::cout << "  ε_xc/electron   : " << exc_pz81(rs) * HA_TO_EV << " eV\n";
    std::cout << "  V_xc/electron   : " << vxc_pz81(rs) * HA_TO_EV << " eV\n";
    check("ΔE_total   (Ha)", E_total_numerical - E_total_predicted, 0.5);

    // ── 7. Write KS orbital data to files ────────────────────────────────────
    //
    // For jellium the KS orbitals are analytically exact plane waves:
    //   ψ_k(r) = exp(i k·r) / √Ω,   k = (2π/L)(nx, ny, nz),   Ω = L³
    //
    // The orbital DENSITY |ψ_k|² = 1/Ω is uniform for ALL plane waves —
    // this is expected and confirms the homogeneity of the ground state.
    //
    // What contains spatial information is Re[ψ_k] = cos(k·r)/√Ω, which
    // oscillates with the wavevector k of the shell.
    //
    // We write a 2D slice at z = L/2 for one representative k-vector per shell.
    // These files are readable by plot_results.py for visualisation.

    fs::create_directories("results/orbitals");

    int    N_g   = static_cast<int>(std::round(L_BOHR / SPACING));
    double h_eff = L_BOHR / N_g;   // actual spacing after integer rounding
    int    iz    = N_g / 2;
    double z_sl  = iz * h_eff;
    double Omega = L_BOHR * L_BOHR * L_BOHR;
    double k0    = 2.0 * M_PI / L_BOHR;
    double norm  = 1.0 / std::sqrt(Omega);

    // Write the shared grid file first (common to all orbitals)
    {
        std::ofstream gf("results/orbitals/grid_slice.txt");
        gf << "# 2D slice coordinates at z = " << z_sl << " bohr\n";
        gf << "# Grid: " << N_g << "x" << N_g << " points,  h = " << h_eff << " bohr\n";
        gf << "# Columns: ix  iy  x_bohr  y_bohr\n";
        gf << std::fixed << std::setprecision(6);
        for (int ix = 0; ix < N_g; ++ix) {
            for (int iy = 0; iy < N_g; ++iy) {
                gf << ix << "  " << iy << "  "
                   << ix * h_eff << "  " << iy * h_eff << "\n";
            }
        }
    }
    std::cout << "\n  Written: results/orbitals/grid_slice.txt\n";

    // Representative k-vectors: one per shell (simplest vector in that shell)
    struct OrbDef {
        int    n2;
        double kx, ky, kz;   // bohr⁻¹
        double ev_pred;       // predicted KS eigenvalue (Ha)
        std::string desc;
    };

    std::vector<OrbDef> orb_defs = {
        {0,    0,    0,    0,    0.0*k0*k0*0.5 + vxc_pred, "k=(0,0,0)"},
        {1,   k0,    0,    0,    1.0*k0*k0*0.5 + vxc_pred, "k=(1,0,0)×k₀"},
        {2,   k0,   k0,    0,    2.0*k0*k0*0.5 + vxc_pred, "k=(1,1,0)×k₀"},
        {3,   k0,   k0,   k0,   3.0*k0*k0*0.5 + vxc_pred, "k=(1,1,1)×k₀"},
        {4, 2*k0,    0,    0,    4.0*k0*k0*0.5 + vxc_pred, "k=(2,0,0)×k₀"},
    };

    for (int oi = 0; oi < static_cast<int>(orb_defs.size()); ++oi) {
        auto const & o = orb_defs[oi];
        std::string base   = "results/orbitals/orbital_" + std::to_string(oi)
                           + "_n2_" + std::to_string(o.n2);
        std::string fname_re = base + "_real.txt";
        std::string fname_im = base + "_imag.txt";

        // Common header lines shared by both files
        auto write_header = [&](std::ofstream & f, std::string const & part) {
            f << "# Jellium KS orbital — analytical plane wave  ψ_k(r) = exp(ik·r)/√Ω\n";
            f << "# Part: " << part << "\n";
            f << "# Shell |n|² = " << o.n2 << "  " << o.desc << "\n";
            f << "# k = (" << o.kx << ", " << o.ky << ", " << o.kz << ") bohr⁻¹\n";
            f << "# Predicted KS eigenvalue  ε = k²/2 + V_xc = " << o.ev_pred
              << " Ha  (" << o.ev_pred * HA_TO_EV << " eV)\n";
            f << "# Cell Ω = " << Omega << " bohr³,   norm = 1/√Ω = " << norm << " bohr^{-3/2}\n";
            f << "# Slice  z = " << z_sl << " bohr  (ix=0.." << N_g-1 << ", iy=0.." << N_g-1 << ")\n";
            if (part == "real")
                f << "# Re[ψ_k(r)] = cos(k·r) / √Ω\n";
            else
                f << "# Im[ψ_k(r)] = sin(k·r) / √Ω\n";
            f << "# Columns: ix  iy  x_bohr  y_bohr  psi_value\n";
        };

        std::ofstream rf(fname_re);
        std::ofstream imf(fname_im);
        write_header(rf,  "real");
        write_header(imf, "imag");

        rf  << std::scientific << std::setprecision(8);
        imf << std::scientific << std::setprecision(8);

        for (int ix = 0; ix < N_g; ++ix) {
            double x = ix * h_eff;
            for (int iy = 0; iy < N_g; ++iy) {
                double y  = iy * h_eff;
                double kr = o.kx * x + o.ky * y + o.kz * z_sl;
                std::string coords = std::to_string(ix) + "  " + std::to_string(iy) + "  ";
                // coordinates in fixed, value in scientific
                rf  << ix << "  " << iy << "  "
                    << std::fixed << std::setprecision(6) << x << "  " << y << "  "
                    << std::scientific << std::setprecision(8)
                    << norm * std::cos(kr) << "\n";
                imf << ix << "  " << iy << "  "
                    << std::fixed << std::setprecision(6) << x << "  " << y << "  "
                    << std::scientific << std::setprecision(8)
                    << norm * std::sin(kr) << "\n";
            }
        }
        std::cout << "  Written: " << fname_re << "\n";
        std::cout << "  Written: " << fname_im << "\n";
    }

    // ── 8. Write eigenvalues for XC offset verification ───────────────────────
    //
    // For jellium: ε_i = k_i²/2 + V_xc(n₀).
    // A scatter plot of ε_i (numerical) vs k_i²/2 (analytical) should give
    // a straight line with slope 1 and intercept V_xc(n₀).
    // The file below contains the data for that plot.

    fs::create_directory("results");
    {
        std::ofstream ef("results/eigenvalues.txt");
        ef << "# KS eigenvalues — jellium ground state\n";
        ef << "# V_xc (PZ81) = " << vxc_pred << " Ha  (" << vxc_pred*HA_TO_EV << " eV)\n";
        ef << "# Slope-1 line:  eigenvalue_Ha = k2_over_2_Ha + V_xc_Ha\n";
        ef << "# Columns: state_idx  shell_n2  k2_over_2_Ha  eigenvalue_Ha  "
              "predicted_Ha  residual_Ha\n";

        // Build state→shell mapping (ascending energy order, same as INQ fills)
        struct StateInfo { int n2; double Ek; };
        std::vector<StateInfo> state_info;
        for (auto const & sh : free_electron_shells(L_BOHR, 8)) {
            for (int s = 0; s < sh.degeneracy; ++s) {
                state_info.push_back({sh.n2, sh.energy_Ha});
            }
        }

        auto evals    = electrons.eigenvalues();   // gpu::array<double,2> [kpin][state]
        int  n_states = evals[0].size();

        ef << std::fixed << std::setprecision(8);
        for (int i = 0; i < n_states; ++i) {
            double ev_num  = evals[0][i];
            int    n2      = (i < (int)state_info.size()) ? state_info[i].n2  : -1;
            double Ek      = (i < (int)state_info.size()) ? state_info[i].Ek  : 0.0;
            double ev_pred = Ek + vxc_pred;
            double resid   = ev_num - ev_pred;
            ef << i << "  " << n2 << "  " << Ek << "  "
               << ev_num << "  " << ev_pred << "  " << resid << "\n";
        }
    }
    std::cout << "\n  Written: results/eigenvalues.txt\n";

    // ── 9. Summary ────────────────────────────────────────────────────────────

    std::cout << "\n========================================\n";
    if (all_passed) {
        std::cout << " All tests PASSED.\n";
        std::cout << " INQ is treating this system as jellium (uniform positive background).\n";
    } else {
        std::cout << " Some tests FAILED — check output above.\n";
    }
    std::cout << "========================================\n\n";

    return all_passed ? 0 : 1;
}
