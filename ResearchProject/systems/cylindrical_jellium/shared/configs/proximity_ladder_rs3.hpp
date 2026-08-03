// shared/configs/proximity_ladder_rs3.hpp
// ----------------------------------------------------------------------------
// LOCKED rung table for the CYLINDRICAL PROXIMITY LADDER: the channeling twin
// walked from grazing (R_in = 2.5 sigma_WP) to fully immersed (filled cylinder)
// at FIXED r_s = 3, fixed projectile energy, fixed cell and grid.
//
// Plan: docs/plans/cylindrical-proximity-ladder.md
// Rung 1 IS the completed channeling twin (channeling_tube_rs3.hpp) — not re-run.
//
// AIM. The channeling twin validated a KS-orbital definition of stopping power
// against the classical dE/ds one in the WEAK-coupling limit. This ladder brings
// the wall inward to find where — and how — that agreement breaks.
//
// WHAT VARIES, AND WHAT DOES NOT. Only five numbers change between rungs:
// R_in, R_out, N_ELECTRONS, EXTRA_STATES, and whether the tube is hollow or
// filled. EVERYTHING else — cell, spacing, edge width, projectile energy,
// sigma_WP, launch point, dt, N_STEPS, save cadence, SCF settings — is inherited
// verbatim from ChannelingTubeRs3 so the rungs are comparable by construction.
//
// HOW R_out IS DERIVED (not assumed). The requirement is n0 EXACTLY equal to the
// r_s = 3 value, because the run binaries set n0 = N/V and the G=0 cancellation
// needs integral n_+ = N exactly. So: pick R_in; pick the nearest EVEN N (even so
// the closed-shell occupied count N/2 is an integer); then SOLVE for R_out:
//
//     V_needed = N / n0_target
//     R_out    = sqrt( R_in^2 + V_needed / (pi * L_z) )
//
// R_out moves by at most 0.014 Bohr across the whole ladder — far inside the
// 0.5 Bohr erfc edge, i.e. physically the same wall. Every rung below lands on
// r_s = 3.000000000 (verified numerically 2026-08-02).
//
// THE COUPLING IS EXPONENTIAL IN R_in, NOT LINEAR. The wavepacket charge outside
// radius R_in is exp(-R_in^2 / 2 sigma_d(t)^2), so the rung labels 2.5/2.0/1.5/1.0
// sigma read linear but the coupling spans a factor of 190 at t = 0:
//
//   R_in     t=0      t=13     t=30       <- sigma_d(t) = sqrt(s^2/2 + t^2/2s^2)
//    10     0.19 %    2.3 %     25 %
//     8     1.83 %    9.0 %     41 %
//     6    10.54 %   25.8 %     61 %
//     4    36.79 %   54.8 %     80 %
//     0      100 %     100 %    100 %
//
// The rungs are DISTINCT ONLY EARLY — which is where S is fitted, so this is
// workable, but the S fit window must be defined by a common measured coupling
// (<n_bath>_WP), never by a common time. See the plan, section 5.
//
// WHAT THIS LADDER CANNOT REACH. The projectile's Gaussian form factor
// exp(-q^2 sigma_pot^2/2) is 0.37 at q = 0.5, 0.018 at q = 1, and 3e-26 at
// q = 2 v0 = 3.83. The plasmon pole sits at q_min = omega_p/v = 0.174 and the
// electron-hole continuum runs to q = 2v. This projectile couples to the
// COLLECTIVE response and essentially nothing else, at EVERY rung. Shrinking R_in
// scales how much medium responds; it does not harden the projectile. The ladder
// is weak-collective -> strong-collective. Reaching the pair channel is a
// sigma_WP axis, not an R_in axis.
// ----------------------------------------------------------------------------
#pragma once

#include "channeling_tube_rs3.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>

namespace cylindrical_jellium::config {

// The r_s = 3 target density, defined by rung 1 exactly as it was built.
// 160 / (pi (14^2 - 10^2) * 60) = 8.8419412828830736e-3 a0^-3.
inline double n0_target() {
	return double(ChannelingTubeRs3::N_ELECTRONS) / ChannelingTubeRs3::v_annulus();
}

struct Rung {
	const char* label;          // run-name token, e.g. "r08"
	double      r_in;           // Bohr; 0 => FILLED cylinder
	int         n_electrons;    // even, so occupied = N/2 is an integer
	int         extra_states;   // ~30 % of occupied, matching rung 1's 24/80

	// A filled tube is a DIFFERENT SHAPE, not an annulus with R_in = 0:
	// inqkit's erfc step is centred ON its nominal edge, so a degenerate inner
	// edge yields n0/2 exactly on the tube axis — precisely where the projectile
	// flies. run.cpp maps this to background_shape::cylinder vs ::annulus.
	bool filled() const { return r_in <= 0.0; }

	// R_out is DERIVED, never transcribed: it is whatever radius makes n0 land
	// exactly on the r_s = 3 target for this R_in and this N.
	//     V_needed = N / n0_target,  R_out = sqrt(R_in^2 + V_needed/(pi L_z))
	// Hard-coding it to a few decimals leaves a ~1e-10 density error and invites
	// a silent transcription slip; deriving it makes n0() exact by construction.
	double r_out() const {
		return std::sqrt(r_in*r_in
		                 + double(n_electrons)
		                   / (n0_target() * M_PI * ChannelingTubeRs3::LZ_BOHR));
	}

	// V = pi (R_out^2 - R_in^2) L_z; correct for the filled case too (R_in = 0).
	double v_jellium() const {
		const double ro = r_out();
		return M_PI * (ro*ro - r_in*r_in) * ChannelingTubeRs3::LZ_BOHR;
	}
	double n0()      const { return double(n_electrons) / v_jellium(); }
	double r_s()     const { return std::cbrt(3.0 / (4.0*M_PI*n0())); }
	double omega_p() const { return std::sqrt(4.0*M_PI*n0()); }
	double k_fermi() const { return std::cbrt(3.0*M_PI*M_PI*n0()); }
	int    n_states()const { return n_electrons/2 + extra_states; }

	// Fraction of the wavepacket's charge OUTSIDE R_in at time t — the
	// transverse Rayleigh tail exp(-R_in^2 / 2 sigma_d(t)^2). This is the
	// nominal coupling; the run's MEASURED <n_bath>_WP supersedes it.
	double wall_overlap(double t_au) const {
		if(filled()) return 1.0;
		const double sd = ChannelingTubeRs3::sigma_d(t_au);
		return std::exp(-r_in*r_in / (2.0*sd*sd));
	}
};

// ---- the ladder -------------------------------------------------------------
// r10 is the COMPLETED channeling twin, listed so the table is the whole story
// and so analysis can iterate rungs uniformly. Do not re-run it.
// R_out (shown for reading only — r_out() derives it) and n_states:
//   r10  R_out 14.0000  104 states      r04  R_out 14.0000  195 states
//   r08  R_out 14.0000  143 states      r00  R_out 13.9857  212 states
//   r06  R_out 13.9857  173 states
inline constexpr Rung LADDER[] = {
	{"r10", 10.0, 160, 24},   // 2.50 sigma  — done (channeling_twin), do not re-run
	{"r08",  8.0, 220, 33},   // 2.00 sigma
	{"r06",  6.0, 266, 40},   // 1.50 sigma
	{"r04",  4.0, 300, 45},   // 1.00 sigma
	{"r00",  0.0, 326, 49},   // filled cylinder
};
inline constexpr int LADDER_N = int(sizeof(LADDER)/sizeof(LADDER[0]));

// ---- the same-N control -----------------------------------------------------
// The ladder moves THREE things at once: proximity, N_e (160 -> 326), and the
// target's mode spectrum (thin annulus with two coupled surfaces -> solid
// nanowire). In this geometry they are inseparable, because the electrons added
// ARE the close ones.
//
// This control sits at R_in = 4 (rung r04) but keeps the ANNULUS VOLUME of rung
// r10, hence N = 160 and n0 unchanged: R_out = sqrt(4^2 + 96) = sqrt(112).
// Comparing its S against r04 (same R_in, N = 300) separates "wall is closer"
// from "there is more wall". Agreement => the ladder is a proximity sweep and
// the far material is inert; disagreement => S must be reported against
// <n_bath>_WP, never against R_in.
// R_out derives to sqrt(4^2 + 96) = 10.5830 Bohr.
inline constexpr Rung CONTROL_SAME_N = {"r04n160", 4.0, 160, 24};

// ---- lookup -----------------------------------------------------------------
inline const Rung* find_rung(const char* label) {
	if(label == nullptr) return nullptr;
	if(std::strcmp(label, CONTROL_SAME_N.label) == 0) return &CONTROL_SAME_N;
	for(int i = 0; i < LADDER_N; ++i)
		if(std::strcmp(label, LADDER[i].label) == 0) return &LADDER[i];
	return nullptr;
}

// Reads CJ_RUNG. Returns nullptr if unset or unknown — the CALLER must fail
// loudly rather than silently defaulting to a rung, because a wrong rung
// produces a perfectly plausible run at the wrong geometry.
inline const Rung* rung_from_env(const char* var = "CJ_RUNG") {
	return find_rung(std::getenv(var));
}

// Self-check: every rung must land on the target density. Returns the worst
// relative error over the whole table (ladder + control). Runs should assert
// this is < 1e-12 at startup — it catches a mistyped R_out immediately, whereas
// a 1 % density error would just look like slightly different physics.
inline double max_density_error() {
	double worst = 0.0;
	const double tgt = n0_target();
	for(int i = 0; i < LADDER_N; ++i)
		worst = std::fmax(worst, std::fabs(LADDER[i].n0()/tgt - 1.0));
	return std::fmax(worst, std::fabs(CONTROL_SAME_N.n0()/tgt - 1.0));
}

}  // namespace cylindrical_jellium::config
