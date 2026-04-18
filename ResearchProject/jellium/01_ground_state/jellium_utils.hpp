#pragma once

// jellium_utils.hpp — Analytical results for the 3D homogeneous electron gas
//
// Provides the analytical quantities needed to validate an INQ jellium
// ground-state calculation: Wigner-Seitz radius, Fermi energy, LDA
// exchange-correlation energy and potential (Perdew-Zunger parametrisation),
// and the non-interacting kinetic energy from the discrete Gamma-point shell
// structure.
//
// All quantities are in Hartree atomic units (ℏ = m_e = e = 1).
//
// Reference for LDA parametrisation:
//   Perdew & Zunger, Phys. Rev. B 23, 5048 (1981) — hereafter PZ81.
//
// NOTE: these utilities are currently local to 01_ground_state. As they
// stabilise and are reused, they will migrate to ResearchProject/jellium/lib/.

#include <cmath>
#include <vector>
#include <algorithm>

// ── Physical constants ───────────────────────────────────────────────────────

constexpr double PI = M_PI;

// ── Basic HEG geometry ───────────────────────────────────────────────────────

// Wigner-Seitz radius r_s (bohr) for N electrons in a cubic cell of side L.
// Definition: (4π/3) r_s³ = 1/n₀  where n₀ = N/L³.
inline double wigner_seitz_radius(int N, double L) {
    double n0 = N / (L * L * L);
    return std::cbrt(3.0 / (4.0 * PI * n0));
}

// Mean electron density n₀ = N/L³ (electrons per bohr³).
inline double mean_density(int N, double L) {
    return N / (L * L * L);
}

// Fermi wavevector k_F = (3π² n₀)^(1/3) in bohr⁻¹.
inline double fermi_wavevector(int N, double L) {
    return std::cbrt(3.0 * PI * PI * mean_density(N, L));
}

// Free-electron Fermi energy E_F = k_F²/2 in Hartree.
inline double fermi_energy(int N, double L) {
    double kF = fermi_wavevector(N, L);
    return 0.5 * kF * kF;
}

// Bulk plasmon frequency ω_p = sqrt(4π n₀) in Hartree (free-electron model).
inline double plasmon_frequency(int N, double L) {
    return std::sqrt(4.0 * PI * mean_density(N, L));
}

// ── LDA exchange-correlation (Perdew-Zunger 1981) ────────────────────────────

// Exchange energy per electron ε_x(r_s) in Hartree.
// From the uniform-gas Dirac exchange: ε_x = -(3/4π)(3π²n)^(1/3)
//   = -3/(4π) × (9π/4)^(1/3) / r_s = -0.4582 / r_s.
// (The numerical coefficient 0.4582 = (3/4)(3/π)^(1/3).)
inline double exchange_energy_pz81(double rs) {
    return -0.4582 / rs;
}

// Exchange potential V_x(r_s) = (4/3) ε_x (from d(n ε_x)/dn).
inline double exchange_potential_pz81(double rs) {
    return (4.0 / 3.0) * exchange_energy_pz81(rs);
}

// Correlation energy per electron ε_c(r_s) in Hartree — PZ81 parametrisation.
// Two regimes matching at r_s = 1:
//   r_s < 1: high-density (perturbative) fit from Ceperley-Alder data.
//   r_s ≥ 1: low-density fit (metallic range, r_s ≈ 2–6 for most metals).
inline double correlation_energy_pz81(double rs) {
    if (rs < 1.0) {
        // High-density regime — PZ81 Eq.(4)
        return 0.0311 * std::log(rs) - 0.048
             + 0.002  * rs * std::log(rs) - 0.0116 * rs;
    } else {
        // Metallic regime — PZ81 Eq.(5)
        // γ, β₁, β₂ from Ceperley-Alder Monte Carlo data (Table I of PZ81)
        constexpr double gamma = -0.1423;
        constexpr double beta1 =  1.0529;
        constexpr double beta2 =  0.3334;
        double sqrtrs = std::sqrt(rs);
        return gamma / (1.0 + beta1 * sqrtrs + beta2 * rs);
    }
}

// Correlation potential V_c = ε_c − (r_s/3) dε_c/dr_s — PZ81 Eq.(6).
inline double correlation_potential_pz81(double rs) {
    if (rs < 1.0) {
        // dε_c/dr_s for high-density regime
        double dec_drs = 0.0311 / rs + 0.002 * (std::log(rs) + 1.0) - 0.0116;
        return correlation_energy_pz81(rs) - (rs / 3.0) * dec_drs;
    } else {
        constexpr double gamma = -0.1423;
        constexpr double beta1 =  1.0529;
        constexpr double beta2 =  0.3334;
        double sqrtrs = std::sqrt(rs);
        double denom  = 1.0 + beta1 * sqrtrs + beta2 * rs;
        double ec = gamma / denom;
        // V_c = ε_c × (1 + (7/6)β₁√r_s + (4/3)β₂r_s) / denom — PZ81 Eq.(6)
        return ec * (1.0 + (7.0/6.0) * beta1 * sqrtrs + (4.0/3.0) * beta2 * rs)
               / denom;
    }
}

// Total LDA xc energy per electron ε_xc = ε_x + ε_c (Hartree).
inline double exc_pz81(double rs) {
    return exchange_energy_pz81(rs) + correlation_energy_pz81(rs);
}

// Total LDA xc potential V_xc = V_x + V_c (Hartree).
// For uniform jellium this is a spatially constant shift of all KS eigenvalues.
inline double vxc_pz81(double rs) {
    return exchange_potential_pz81(rs) + correlation_potential_pz81(rs);
}

// ── Discrete shell structure (Gamma-point) ───────────────────────────────────
//
// For a cubic box with periodic BCs, allowed wavevectors are k = (2π/L)n
// where n ∈ ℤ³. States group into degenerate shells labelled by |n|².
// This is the exact quantum number for the Gamma-point free-electron spectrum.

struct Shell {
    int    n2;          // |n|² = nx² + ny² + nz²
    int    degeneracy;  // number of (nx,ny,nz) triplets with this |n|²
    double energy_Ha;   // k²/2 = 2π²|n|²/L² in Hartree
};

// Enumerate free-electron shells up to |n|²_max for a cubic cell of side L.
inline std::vector<Shell> free_electron_shells(double L, int n2_max = 6) {
    double k0 = 2.0 * PI / L;   // smallest reciprocal wavevector
    std::vector<Shell> shells;

    // Count degeneracy per |n|² by brute-force loop over a small box.
    int nmax = (int)std::sqrt((double)n2_max) + 1;
    std::vector<int> degen(n2_max + 1, 0);
    for (int nx = -nmax; nx <= nmax; ++nx) {
        for (int ny = -nmax; ny <= nmax; ++ny) {
            for (int nz = -nmax; nz <= nmax; ++nz) {
                int n2 = nx*nx + ny*ny + nz*nz;
                if (n2 <= n2_max) degen[n2]++;
            }
        }
    }

    for (int n2 = 0; n2 <= n2_max; ++n2) {
        if (degen[n2] > 0) {
            Shell s;
            s.n2        = n2;
            s.degeneracy = degen[n2];
            s.energy_Ha = 0.5 * n2 * k0 * k0;
            shells.push_back(s);
        }
    }
    return shells;
}

// Non-interacting kinetic energy T_s (Ha) from the Gamma-point shell sum
// for N_electrons in a cubic cell of side L at T = 0.
//
// Fills shells in order of increasing energy until N_electrons is reached.
// The final shell may be partially occupied if N_electrons doesn't exactly
// fill it (as is the case for N=40 at the |n|²=3 shell).
//
// At finite temperature (~100 K for metals), the correction to T_s is
// negligible (kT/E_F ≈ 0.001), so this T=0 estimate is adequate.
inline double kinetic_energy_shells(int N_electrons, double L) {
    auto shells = free_electron_shells(L, 10);

    // Each spatial orbital holds 2 electrons (spin up + spin down).
    int remaining = N_electrons;
    double T_s = 0.0;

    for (auto const & sh : shells) {
        if (remaining <= 0) break;
        int max_electrons = 2 * sh.degeneracy;   // spin factor of 2
        int fill = std::min(remaining, max_electrons);
        // Number of electrons in this shell × ε_k
        T_s += fill * sh.energy_Ha;
        remaining -= fill;
    }
    return T_s;
}

// Predicted total DFT energy for jellium in Hartree:
//   E_total ≈ T_s + N × ε_xc(r_s)
//
// This holds because for uniform jellium in periodic BC:
//   E_Hartree ≈ 0 (G≠0 density components vanish for uniform ρ;
//                   G=0 divergence is cancelled by the positive background)
//   E_external = 0 (no ionic pseudopotentials)
//   E_ion      = 0 (no nuclei)
// So the only terms surviving are T_s and E_xc.
inline double predicted_total_energy(int N_electrons, double L) {
    double rs  = wigner_seitz_radius(N_electrons, L);
    double T_s = kinetic_energy_shells(N_electrons, L);
    double E_xc = N_electrons * exc_pz81(rs);
    return T_s + E_xc;
}
