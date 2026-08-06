// Engine test for inqkit::SelfInteractionCorrection (plan
// docs/plans/wp-self-interaction-correction.md, §4 "engine test").
//
// Three known-case checks, each pinning one failure mode the vacuum/jellium
// tiers cannot separate on their own:
//
//   1. KICK SEMANTICS (1 electron, vacuum). A real multiplicative kick must
//      leave the density and <p> untouched (zero-force: integral n grad
//      v_H[n] = 0) at ANY dt_eff, while a large dt_eff must visibly inflate
//      var(p_z) (the phase gradient is real momentum). Also norm exactly
//      conserved and no projection in vacuum (n_projected == 0).
//
//   2. Q PROJECTION (synthetic occupied manifold). After an exaggerated kick
//      the WP has leaked into the occupied states; apply() must return it to
//      orthogonality at the 1e-10 level, restore norm 1, report the leak it
//      removed (max_overlap_pre > 0, norm_removed > 0), and leave every bath
//      column bit-identical (the correction acts on ONE column only).
//
//   3. RUN-CONSISTENCY OF THE SUBTRACTED ENERGIES (plan §0/D1). For a
//      1-electron system the total density IS n_wp, so INQ's own broadcast
//      scalars are closed-form references: measure().u_self must equal
//      energy_hartree and measure().exc_self must equal energy_xc. This is
//      the check that the xc self-term is evaluated with the run's spin
//      treatment (unpolarised) through the run's own code path — canonical
//      polarised PZ would FAIL it by ~2^(1/3) on the exchange part.
//
// Reference for the SIC scheme: Perdew & Zunger, PRB 23, 5048 (1981);
// Messud, Dinh, Reinhard, Suraud, PRL 101, 096404 (2008).

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/wavepacket/self_interaction_correction.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>

#include <cmath>
#include <filesystem>
#include <vector>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;
namespace fs = std::filesystem;
using SIC = inqkit::SelfInteractionCorrection;

namespace {

// <psi_a|psi_b> on the single-rank gamma set (re, im).
template <class Electrons>
std::pair<double,double> overlap(Electrons const & el, int a, int b) {
	auto & phi = el.kpin()[0];
	const int n_pts = phi.basis().local_size();
	const double dV = phi.basis().volume_element();
	auto mat_ = begin(phi.matrix());
	auto re = gpu::run(1, gpu::reduce(n_pts), 0.0,
		[dV, mat_, a, b] GPU_LAMBDA (auto, auto ip) {
			auto va = mat_[ip][a]; auto vb = mat_[ip][b];
			return dV*(inq::real(va)*inq::real(vb) + inq::imag(va)*inq::imag(vb));
		});
	auto im = gpu::run(1, gpu::reduce(n_pts), 0.0,
		[dV, mat_, a, b] GPU_LAMBDA (auto, auto ip) {
			auto va = mat_[ip][a]; auto vb = mat_[ip][b];
			return dV*(inq::real(va)*inq::imag(vb) - inq::imag(va)*inq::real(vb));
		});
	return {re[0], im[0]};
}

template <class Electrons>
double overlap_abs(Electrons const & el, int a, int b) {
	auto [re, im] = overlap(el, a, b);
	return std::sqrt(re*re + im*im);
}

} // namespace

TEST_CASE("SIC kick: density/<p> invariant, var(p) responds, vacuum has no projection",
          "[wavepacket][sic][engine]") {
	const double L = 20.0, SPACING = 0.5, SIGMA = 2.0, K0Z = 0.8;
	auto cell = systems::cell::orthorhombic(L*1.0_b, L*1.0_b, L*1.0_b).periodic();
	auto ions = systems::ions(cell);
	const double ec_ha = 0.5*std::pow(M_PI/SPACING, 2.0);
	auto electrons = systems::electrons(
		ions, options::electrons{}.cutoff(ec_ha*1.0_Ha).extra_electrons(1.0));
	ground_state::initial_guess(ions, electrons);

	auto report = inqkit::WavePacket{}
		.center(0.0, 0.0, 0.0).sigma(SIGMA).k0(0.0, 0.0, K0Z)
		.inject_into_last_extra_state(electrons, 1.0);
	const int wp_idx = report.state_index;
	REQUIRE(wp_idx == 0);   // 1 electron -> single state

	fs::path dir = fs::temp_directory_path() / "inqkit_test_sic";
	fs::remove_all(dir); fs::create_directories(dir);
	inqkit::observables::WPRealSpaceStats wp_rs((dir/"rs.csv").string(), wp_idx);
	inqkit::observables::WPMomentumStats  wp_mom((dir/"mom.csv").string(), wp_idx);

	auto m0 = wp_rs.compute(electrons);
	auto p0 = wp_mom.compute(electrons);

	SIC sic(SIC::Mode::pz_run, wp_idx);

	// small kick: EVERYTHING invariant to tight tolerance
	auto r1 = sic.apply(electrons, 0.02);
	CHECK(r1.kicked);
	CHECK(r1.n_projected == 0);                 // vacuum: Q = 1
	CHECK(r1.u_self > 0.0);                     // self-Hartree of a positive density
	CHECK(r1.exc_self < 0.0);                   // LDA xc energy is negative
	auto m1 = wp_rs.compute(electrons);
	auto p1 = wp_mom.compute(electrons);
	CHECK(m1.N   == Approx(m0.N).epsilon(1e-10));       // norm untouched
	CHECK(m1.sz2 == Approx(m0.sz2).epsilon(1e-8));      // density untouched
	CHECK(m1.z   == Approx(m0.z).margin(1e-8));
	CHECK(p1.pz  == Approx(p0.pz).margin(1e-5));        // zero-force on itself
	CHECK(p1.sz2 == Approx(p0.sz2).epsilon(2e-3));      // O(dt^2) only

	// exaggerated kick: density still invariant, but var(p_z) must inflate --
	// the phase gradient is genuine momentum content, proof the kick acted.
	auto r2 = sic.apply(electrons, 50.0);
	CHECK(r2.kicked);
	auto m2 = wp_rs.compute(electrons);
	auto p2 = wp_mom.compute(electrons);
	CHECK(m2.N   == Approx(m0.N).epsilon(1e-10));
	CHECK(m2.sz2 == Approx(m0.sz2).epsilon(1e-6));      // still a pure phase
	CHECK(p2.pz  == Approx(p0.pz).margin(5e-3));        // zero-force holds large
	CHECK(p2.sz2 > 1.5*p0.sz2);                         // momentum spread inflated
}

TEST_CASE("SIC projection: exact re-orthogonalisation, leak reported, bath untouched",
          "[wavepacket][sic][engine]") {
	const double L = 16.0, SPACING = 0.5, SIGMA = 2.0;
	auto cell = systems::cell::orthorhombic(L*1.0_b, L*1.0_b, L*1.0_b).periodic();
	auto ions = systems::ions(cell);
	const double ec_ha = 0.5*std::pow(M_PI/SPACING, 2.0);
	// 4 bath electrons -> 2 occupied states + 1 extra slot for the WP.
	auto electrons = systems::electrons(
		ions, options::electrons{}.cutoff(ec_ha*1.0_Ha).extra_electrons(4.0).extra_states(1));
	ground_state::initial_guess(ions, electrons);
	// The guess is not orthonormal; the projector's exactness statement is for
	// an ORTHONORMAL occupied manifold (as in production, where the GS provides
	// it), so manufacture one.
	operations::orthogonalize(electrons.kpin()[0]);

	auto report = inqkit::WavePacket{}
		.center(0.0, 0.0, 0.0).sigma(SIGMA).k0(0.0, 0.0, 0.8)
		.orthogonalise_against_occupied(electrons)
		.inject_into_last_extra_state(electrons, 1.0);
	const int wp_idx = report.state_index;
	REQUIRE(wp_idx == 2);

	// bath fingerprints before
	const double n0_before  = overlap_abs(electrons, 0, 0);
	const double n1_before  = overlap_abs(electrons, 1, 1);
	const double x01_before = overlap_abs(electrons, 0, 1);
	REQUIRE(overlap_abs(electrons, 0, wp_idx) < 1e-8);
	REQUIRE(overlap_abs(electrons, 1, wp_idx) < 1e-8);

	SIC sic(SIC::Mode::hartree, wp_idx);
	// exaggerated dt so the kick leaks measurably into the occupied manifold
	auto r = sic.apply(electrons, 5.0);

	CHECK(r.kicked);
	CHECK(r.n_projected == 2);
	CHECK(r.max_overlap_pre > 1e-8);            // there WAS a leak...
	CHECK(r.norm_removed > 0.0);                // ...and Q removed weight
	CHECK(r.norm_removed < 0.5);                // sanity: not annihilated

	// ...and afterwards orthogonality is restored at numerical precision
	CHECK(overlap_abs(electrons, 0, wp_idx) < 1e-10);
	CHECK(overlap_abs(electrons, 1, wp_idx) < 1e-10);
	CHECK(overlap_abs(electrons, wp_idx, wp_idx) == Approx(1.0).epsilon(1e-9));

	// bath columns untouched by construction
	CHECK(overlap_abs(electrons, 0, 0) == Approx(n0_before).epsilon(1e-12));
	CHECK(overlap_abs(electrons, 1, 1) == Approx(n1_before).epsilon(1e-12));
	CHECK(overlap_abs(electrons, 0, 1) == Approx(x01_before).margin(1e-12));
}

TEST_CASE("SIC energies match INQ's own scalars for a 1-electron system (D1 run-consistency)",
          "[wavepacket][sic][engine]") {
	const double L = 20.0, SPACING = 0.5, SIGMA = 2.0;
	auto cell = systems::cell::orthorhombic(L*1.0_b, L*1.0_b, L*1.0_b).periodic();
	auto ions = systems::ions(cell);
	const double ec_ha = 0.5*std::pow(M_PI/SPACING, 2.0);
	auto electrons = systems::electrons(
		ions, options::electrons{}.cutoff(ec_ha*1.0_Ha).extra_electrons(1.0));
	ground_state::initial_guess(ions, electrons);

	auto report = inqkit::WavePacket{}
		.center(0.0, 0.0, 0.0).sigma(SIGMA).k0(0.0, 0.0, 0.8)
		.inject_into_last_extra_state(electrons, 1.0);
	const int wp_idx = report.state_index;

	SIC sic(SIC::Mode::pz_run, wp_idx);
	auto mine = sic.measure(electrons);

	// One propagate call, 1 step; the step-0 callback fires BEFORE any step, so
	// its energies are those of the injected state — the total density of which
	// IS n_wp (occ 1, single state). INQ's own scalars are then the reference:
	//   energy_hartree = 1/2 integral n phi[n]     == u_self
	//   energy_xc      = E_xc^unpol[n]             == exc_self
	double e_hartree_inq = 0.0, e_xc_inq = 0.0;
	bool got = false;
	real_time::propagate(ions, electrons,
		[&](auto const & data) {
			if(!got) { e_hartree_inq = data.energy().hartree();
			           e_xc_inq = data.energy().xc(); got = true; }
		},
		options::theory{}.lda(),
		options::real_time{}.num_steps(1).dt(0.01*1.0_atomictime));

	REQUIRE(got);
	CHECK(mine.u_self   == Approx(e_hartree_inq).epsilon(1e-9));
	CHECK(mine.exc_self == Approx(e_xc_inq).epsilon(1e-9));
	// and the polarised-PZ exchange would differ by 2^(1/3): assert we are NOT
	// accidentally polarised (the D1 defect this test exists to catch). LDA
	// exchange dominates exc for this density, so a 26 % shift is far outside
	// the tolerance above; this CHECK documents the discrimination explicitly.
	CHECK(std::abs(std::cbrt(2.0)*mine.exc_self - e_xc_inq)
	      > 0.05*std::abs(e_xc_inq));
}
