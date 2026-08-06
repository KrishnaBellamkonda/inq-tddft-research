// Known-case test for WavePacket::minimum_image(), added 2026-08-01.
//
// THE BUG THIS PINS. It is commonly said that a wavepacket needs no boundary
// treatment because a KS orbital lives on a plain 3-D FFT basis and wraps
// exactly. That is true of the PROPAGATION and FALSE of the INJECTION: the
// builder forms the Gaussian from a plain Cartesian displacement, so a packet
// launched within a couple of sigma of a cell face is TRUNCATED.
//
// AND THE NORM WILL NOT TELL YOU. inject_into_last_extra_state renormalises the
// packet ONLY inside its orthogonalisation branch — which every production run
// takes. The missing weight is therefore rescaled away and the norm reads 1.000
// while every momentum observable stays corrupted. Without orthogonalisation the
// same packet keeps a visibly depressed norm (0.856 here). Both behaviours are
// pinned below, because the second is what makes the first dangerous.
//
// The damage is not subtle, and it is worst in exactly the observable a stopping
// study depends on. A truncated Gaussian has a step discontinuity in real space,
// which is broadband in momentum: the production channeling run measured
// var(p_z) FIFTEEN TIMES too large (0.473 vs 0.0313), T1-T2 +471 %, the centroid
// pulled 1.03 Bohr toward the box interior and sigma_z 25 % too small — while x
// and y, far from any face, were perfect to 0.006 %. Six of nine t=0 gates failed.
//
// EXPECTATIONS ARE ANALYTIC, for psi = exp(-|r-b|^2/(2 sigma^2)) exp(i k.r):
//   <p_z>        = k0                       (drift momentum)
//   var(p_z)     = 1/(2 sigma^2)            (minimum-uncertainty packet)
//   <z>_circ     = b_z                      (Resta phase estimator, periodic-exact)
//   sigma_z,circ = sigma/sqrt2              (density std)
// The circular estimators are used because the packet DOES straddle the face —
// that is the point — and the naive moments are meaningless there by construction.
//
// The contrast case (minimum_image OFF, same geometry) is asserted to FAIL these,
// because a test that only checks the fixed path would pass just as happily if
// the flag silently did nothing.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <cmath>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

namespace {

constexpr double L_BOX   = 24.0;   // Bohr, cubic periodic
constexpr double SPACING = 0.5;    // -> 48 points per side
constexpr double SIGMA   = 3.0;    // sigma_WP; density std = SIGMA/sqrt2 = 2.1213
constexpr double K0      = 1.5;    // k0*L = 36 rad = 5.73 x 2pi — deliberately
                                   // NOT a multiple of 2pi, so a wrapped
                                   // amplitude with an unwrapped phase would
                                   // carry a visible seam discontinuity.
constexpr double LAUNCH  = -10.0;  // 2 Bohr from the -z face = 0.94 density sigma

systems::electrons make_electrons() {
	systems::ions ions(systems::cell::cubic(L_BOX * 1.0_bohr).periodic());
	systems::electrons electrons(
		ions, options::electrons{}.spacing(SPACING * 1.0_bohr)
		                          .extra_electrons(2)
		                          .extra_states(4));
	ground_state::initial_guess(ions, electrons);
	return electrons;
}

struct Moments { double pz, var_pz, z_circ, sz_circ, norm; };

// `ortho` mirrors the production call chain. It matters MORE than it looks:
// inject_into_last_extra_state RENORMALISES the packet only inside the
// orthogonalisation branch (wavepacket.hpp: "No orthogonalisation requested:
// ... the packet is left exactly as constructed (not even renormalised)").
// So a truncated packet keeps a visibly depressed norm without ortho, and has
// that evidence SCALED AWAY with it — which is exactly why the production run,
// which always orthogonalises, reported norm = 1.000 while carrying a packet
// missing 17 % of its weight.
Moments inject_and_measure(bool minimum_image, bool ortho = false) {
	auto electrons = make_electrons();
	auto wp = inqkit::WavePacket{}
	              .center(0.0, 0.0, LAUNCH)
	              .sigma(SIGMA)
	              .k0(0.0, 0.0, K0)
	              .minimum_image(minimum_image);
	if(ortho) wp.orthogonalise_against_occupied(electrons);
	auto report = wp.inject_into_last_extra_state(electrons, 1.0);
	const int idx = report.state_index;

	inqkit::observables::WPMomentumStats  mom("/dev/null", idx, {.write_every = 0});
	inqkit::observables::WPRealSpaceStats pos("/dev/null", idx, {.write_every = 0});
	auto m = mom.compute(electrons);
	auto r = pos.compute(electrons);
	return {m.pz, m.sz2, r.zc, r.szc, r.N};
}

const double VAR_P_FREE = 1.0 / (2.0 * SIGMA * SIGMA);   // 0.055556
const double SIGMA_DENS = SIGMA / std::sqrt(2.0);        // 2.12132

} // namespace

TEST_CASE("WavePacket minimum_image: a packet launched near a face is correct",
          "[wavepacket][minimum_image][engine]") {
	auto w = inject_and_measure(true);

	CHECK(w.norm    == Approx(1.0).margin(0.02));
	CHECK(w.pz      == Approx(K0).epsilon(0.01));
	CHECK(w.var_pz  == Approx(VAR_P_FREE).epsilon(0.05));
	CHECK(w.z_circ  == Approx(LAUNCH).margin(0.05));
	CHECK(w.sz_circ == Approx(SIGMA_DENS).epsilon(0.05));
}

TEST_CASE("WavePacket minimum_image: OFF, the same packet is measurably truncated",
          "[wavepacket][minimum_image][engine]") {
	auto off = inject_and_measure(false);
	auto on  = inject_and_measure(true);

	// WITHOUT renormalisation the lost weight is plainly visible: a packet
	// launched 0.94 density sigma from the face keeps only Phi(0.94) = 0.83 of
	// itself. (Measured 0.856 — the residual difference is the wrapped lobe the
	// minimum-image version would have kept.)
	CHECK(off.norm < 0.92);
	CHECK(on.norm  == Approx(1.0).margin(0.02));

	// The sharp real-space edge is broadband in momentum. This is the decisive
	// signature and the one that ruins any observable built on var(p).
	CHECK(off.var_pz > 2.0 * VAR_P_FREE);
	CHECK(on.var_pz  < 1.5 * VAR_P_FREE);

	// Cutting the tail nearest the face pulls the centroid INTO the box ...
	CHECK(off.z_circ > LAUNCH + 0.2);
	// ... and leaves a narrower packet than it claims to be.
	CHECK(off.sz_circ < 0.95 * SIGMA_DENS);

	// The drift momentum is dragged down too, though far less dramatically.
	CHECK(off.pz < K0);
}

TEST_CASE("WavePacket minimum_image: orthogonalisation HIDES the truncation in the norm",
          "[wavepacket][minimum_image][engine]") {
	// This is the production call chain, and the reason the bug survived to a GPU
	// run: orthogonalise_against_occupied renormalises, so the missing 17 % is
	// rescaled away and the norm reads 1.000 — while every momentum observable
	// stays corrupted. The norm is therefore NOT a usable check for this failure;
	// var(p_z) is.
	auto off = inject_and_measure(false, /*ortho=*/true);
	auto on  = inject_and_measure(true,  /*ortho=*/true);

	CHECK(off.norm == Approx(1.0).margin(0.02));   // looks perfectly healthy ...
	CHECK(off.var_pz > 2.0 * VAR_P_FREE);          // ... and is not
	CHECK(on.var_pz  == Approx(VAR_P_FREE).epsilon(0.10));
	CHECK(on.z_circ  == Approx(LAUNCH).margin(0.10));
}

TEST_CASE("WavePacket minimum_image: no effect on a packet far from every face",
          "[wavepacket][minimum_image][engine]") {
	// Backwards-compatibility guarantee: at the box centre the flag must be a
	// no-op, so switching it on cannot perturb any previously published run.
	auto centred = [](bool mi) {
		auto electrons = make_electrons();
		auto wp = inqkit::WavePacket{}
		              .center(0.0, 0.0, 0.0)
		              .sigma(SIGMA)
		              .k0(0.0, 0.0, K0)
		              .minimum_image(mi);
		auto rep = wp.inject_into_last_extra_state(electrons, 1.0);
		inqkit::observables::WPMomentumStats mom("/dev/null", rep.state_index, {.write_every = 0});
		auto m = mom.compute(electrons);
		return std::pair<double,double>{m.pz, m.sz2};
	};
	auto a = centred(false);
	auto b = centred(true);

	// 5.7 density sigma from either face: the wrapped tail is ~1e-8, so the two
	// constructions agree far more tightly than the physics tolerances above.
	CHECK(a.first  == Approx(b.first).epsilon(1e-6));
	CHECK(a.second == Approx(b.second).epsilon(1e-6));
	CHECK(b.second == Approx(VAR_P_FREE).epsilon(0.05));
}
