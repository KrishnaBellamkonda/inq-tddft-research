/*
 * inqkit::jellium::analytics — closed-form jellium quantities for setup,
 * energy bookkeeping, and validation of the localised jellium.
 *
 * Host-only (no INQ engine): plain doubles, usable in run.cpp setup, analyse.py
 * cross-checks, and the validation tests. Atomic units.
 *
 * References: docs/notes/localised-jellium-theory.md (Parts 1–2). HEG limit
 * formulas: Parr & Yang; Gell-Mann–Brueckner kinetic/exchange terms.
 */

#ifndef INQKIT__JELLIUM__ANALYTICS
#define INQKIT__JELLIUM__ANALYTICS

#include <cmath>

namespace inqkit {
namespace jellium {

// Interior density from the Wigner–Seitz radius: n₀ = 3/(4π r_s³).
inline double n0_from_rs(double rs) { return 3.0 / (4.0 * M_PI * rs * rs * rs); }

// And the inverse: r_s from a density.
inline double rs_from_n0(double n0) { return std::cbrt(3.0 / (4.0 * M_PI * n0)); }

// Fermi wavevector of the interior HEG: k_F = (9π/4)^{1/3}/r_s = (3π²n)^{1/3}.
inline double k_fermi_rs(double rs) { return std::cbrt(9.0 * M_PI / 4.0) / rs; }
inline double k_fermi_n0(double n0) { return std::cbrt(3.0 * M_PI * M_PI * n0); }

// Fermi energy E_F = k_F²/2 (Hartree).
inline double e_fermi_rs(double rs) { double k = k_fermi_rs(rs); return 0.5 * k * k; }

// HEG energy-per-electron pieces (Hartree). Kinetic (3/5)E_F and exchange
// −(3/4)(3/π)^{1/3} n^{1/3}. Correlation NOT included (use libxc / a PZ81 fit).
inline double e_kinetic_rs(double rs)  { return 1.104950 / (rs * rs); }
inline double e_exchange_rs(double rs) { return -0.458165 / rs; }

// Sphere radius enclosing N electrons at density r_s: R_cl = r_s · N^{1/3}.
inline double sphere_radius(int n_electrons, double rs) {
	return rs * std::cbrt(static_cast<double>(n_electrons));
}

// Electrostatic SELF-energy of a uniform charged sphere of charge N, radius R:
// E_self = (3/5) N²/R (Hartree). INQ never sees n₊ as a charge, so this must be
// added by hand for the cluster-energy → HEG-limit benchmark (theory Part 3.4).
inline double e_self_sphere(int n_electrons, double radius) {
	double N = static_cast<double>(n_electrons);
	return 0.6 * N * N / radius;
}

} // namespace jellium
} // namespace inqkit

#endif
