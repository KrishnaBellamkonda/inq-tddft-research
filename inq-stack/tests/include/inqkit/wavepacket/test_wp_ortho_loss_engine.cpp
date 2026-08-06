// Engine-tier test: the orthogonalisation-LOSS bookkeeping added to
// InjectionReport (norm_pre_ortho, norm_pre_renorm, removed_weight,
// sum_overlap_sq, ortho_closure_residual).
//
// WHY THIS EXISTS. `norm_after` is measured AFTER the post-Gram-Schmidt
// renormalisation, so it is ~1 by construction and cannot express how much of
// the packet the projection carved away; `max_overlap` reports only the single
// largest overlap. Launching a wavepacket close to a metal surface puts it
// inside the electronic spill-out, where the projection removes a real fraction
// of the Gaussian and the renormalise then hides the loss. `removed_weight` is
// the honest measure, and this test pins it to closed-form values.
//
// THE ANALYTIC CASE. Take a single occupied state equal to the CONSTANT
// phi_0 = 1/sqrt(V) — the exact ground state of a bare kinetic Hamiltonian in
// an empty periodic cell. The injector writes
//     psi(r) = (pi sigma^2)^{-3/4} exp(-|r-b|^2 / (2 sigma^2)) exp(i k.r),
// which is continuum-normalised. Then
//     <phi_0|psi> = V^{-1/2} (pi sigma^2)^{-3/4} \int e^{-r^2/(2 sigma^2)} e^{i k0 z} d^3r
//                 = V^{-1/2} (pi sigma^2)^{-3/4} (sqrt(2 pi) sigma)^3 e^{-sigma^2 k0^2 / 2}
//                 = V^{-1/2} 2^{3/2} pi^{3/4} sigma^{3/2} e^{-sigma^2 k0^2 / 2}
// so the weight removed by projecting it out is
//
//     removed_weight = 8 pi^{3/2} sigma^3 exp(-sigma^2 k0^2) / V.               (*)
//
// The k0 dependence is the physically important part: a fast packet is
// EXPONENTIALLY orthogonal to smooth low-k occupied states. That is exactly why
// the slowest velocity is the worst case when launching near a slab.
//
// State 0 is overwritten by hand rather than obtained from the SCF, so the test
// depends on NOTHING but inqkit's own bookkeeping.
//
// Standard Gaussian integrals; see e.g. Sakurai & Napolitano, "Modern Quantum
// Mechanics", free Gaussian wave-packet section.
//
// Plan: docs/plans/effective-sigma-near-launch.md

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <cmath>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

namespace {

constexpr double L       = 20.0;   // Bohr, cubic periodic  -> V = 8000
constexpr double SPACING = 0.5;    // Bohr
constexpr double SIGMA   = 2.0;    // injector sigma (4 grid points per sigma)
constexpr double VOL     = L * L * L;

// Closed-form (*) above.
double analytic_removed_weight(double sigma, double k0, double vol) {
    return 8.0 * std::pow(M_PI, 1.5) * sigma * sigma * sigma
         * std::exp(-sigma * sigma * k0 * k0) / vol;
}

// Overwrite state 0 with the constant 1/sqrt(V). Real, positive, norm 1:
//   \int |phi_0|^2 dV = V * (1/sqrt(V))^2 = 1.
//
// This MUST be a free function, not a constructor body: nvcc rejects an
// extended __device__ lambda inside a constructor ("the enclosing parent
// function must allow its address to be taken").
void impose_constant_state0(systems::electrons& electrons) {
    auto& phi   = electrons.kpin()[0];
    auto& basis = phi.basis();
    const double amp = 1.0 / std::sqrt(VOL);
    auto phicub = begin(phi.hypercubic());
    gpu::run(basis.local_sizes()[2], basis.local_sizes()[1], basis.local_sizes()[0],
             [=] GPU_LAMBDA(auto iz, auto iy, auto ix) {
                 phicub[ix][iy][iz][0] = inq::complex(amp, 0.0);
             });
    INQKIT_GPU_SYNC();
}

// A 2-state system whose state 0 is EXACTLY the constant 1/sqrt(V).
// extra_electrons(2.0) gives INQ a non-empty electron count; extra_states(1)
// supplies the slot the WP is injected into (so ist_wp = 1 and the
// Gram-Schmidt loop runs over state 0 alone).
struct ConstantBathSystem {
    systems::ions      ions;
    systems::electrons electrons;

    ConstantBathSystem()
        : ions(systems::cell::orthorhombic(L * 1.0_b, L * 1.0_b, L * 1.0_b).periodic()),
          electrons(ions, options::electrons{}
                              .spacing(SPACING * 1.0_b)
                              .extra_states(1)
                              .extra_electrons(2.0)) {
        ground_state::initial_guess(ions, electrons);
        impose_constant_state0(electrons);
    }
};

} // namespace

TEST_CASE("WP ortho loss: removed_weight matches the closed form at k0 = 0",
          "[wavepacket][ortho][injection][engine]") {
    ConstantBathSystem sys;

    auto report = inqkit::WavePacket{}
                      .center(0.0, 0.0, 0.0).sigma(SIGMA).k0(0.0, 0.0, 0.0)
                      .orthogonalise_against_occupied(sys.electrons)
                      .inject_into_last_extra_state(sys.electrons, 1.0);

    REQUIRE(report.orthogonalised);
    REQUIRE(report.state_index == 1);

    const double expected = analytic_removed_weight(SIGMA, 0.0, VOL);  // 0.0445464
    INFO("expected removed_weight = " << expected
         << "  got = " << report.removed_weight);
    // 2 % covers the finite-grid quadrature and the e^{-12.5} ~ 4e-6 periodic
    // truncation of the Gaussian tail at |z| = L/2 = 5 sigma.
    CHECK(report.removed_weight == Approx(expected).epsilon(0.02));

    // A non-trivial amount really was removed (guards a silent no-op).
    CHECK(report.removed_weight > 0.01);
    // ...and the raw Gaussian was essentially normalised to begin with, so the
    // loss is orthogonalisation and not a normalisation artefact.
    CHECK(report.norm_pre_ortho == Approx(1.0).epsilon(0.01));
    CHECK(report.norm_pre_renorm < report.norm_pre_ortho);
}

TEST_CASE("WP ortho loss: k0 suppresses the loss as exp(-sigma^2 k0^2)",
          "[wavepacket][ortho][injection][engine]") {
    ConstantBathSystem sys;

    const double K0 = 1.0;
    auto report = inqkit::WavePacket{}
                      .center(0.0, 0.0, 0.0).sigma(SIGMA).k0(0.0, 0.0, K0)
                      .orthogonalise_against_occupied(sys.electrons)
                      .inject_into_last_extra_state(sys.electrons, 1.0);

    // exp(-sigma^2 k0^2) = exp(-4) = 0.0183 -> 8.16e-4, a 55x suppression.
    const double expected = analytic_removed_weight(SIGMA, K0, VOL);
    INFO("expected removed_weight = " << expected
         << "  got = " << report.removed_weight);
    CHECK(report.removed_weight == Approx(expected).epsilon(0.05));

    // The suppression itself, stated as a ratio — this is the property that
    // makes the SLOWEST velocity the worst case for a near-slab launch.
    const double at_k0_zero = analytic_removed_weight(SIGMA, 0.0, VOL);
    CHECK(report.removed_weight < 0.05 * at_k0_zero);
}

TEST_CASE("WP ortho loss: the two independent routes close on each other",
          "[wavepacket][ortho][injection][engine]") {
    ConstantBathSystem sys;

    auto report = inqkit::WavePacket{}
                      .center(0.0, 0.0, 0.0).sigma(SIGMA).k0(0.0, 0.0, 0.0)
                      .orthogonalise_against_occupied(sys.electrons)
                      .inject_into_last_extra_state(sys.electrons, 1.0);

    // sum_i |<psi_i|psi_wp>|^2  ==  ||psi||^2_pre - ||psi||^2_post, exactly,
    // because the KS states are mutually orthonormal (no cross-terms). The two
    // sides come from completely different reductions: a per-state overlap loop
    // versus a single global norm.
    const double lhs = report.sum_overlap_sq;
    const double rhs = report.norm_pre_ortho * report.norm_pre_ortho
                     - report.norm_pre_renorm * report.norm_pre_renorm;
    INFO("sum_overlap_sq = " << lhs << "  norm^2 difference = " << rhs
         << "  residual = " << report.ortho_closure_residual());
    CHECK(lhs > 0.0);
    CHECK(report.ortho_closure_residual() < 1.0e-10 * lhs);

    // max_overlap is the single largest |<psi_i|psi_wp>|; with exactly one
    // state below the WP slot it must equal sqrt(sum_overlap_sq).
    CHECK(report.max_overlap == Approx(std::sqrt(lhs)).epsilon(1.0e-10));
}

TEST_CASE("WP ortho loss: back-compat contract of the existing fields",
          "[wavepacket][ortho][injection][engine]") {
    SECTION("with orthogonalisation the packet is still renormalised to 1") {
        ConstantBathSystem sys;
        auto report = inqkit::WavePacket{}
                          .center(0.0, 0.0, 0.0).sigma(SIGMA).k0(0.0, 0.0, 0.0)
                          .orthogonalise_against_occupied(sys.electrons)
                          .inject_into_last_extra_state(sys.electrons, 1.0);
        // This is precisely why norm_after cannot report the loss.
        CHECK(report.norm_after == Approx(1.0).epsilon(1.0e-9));
        CHECK(report.removed_weight > 0.01);
    }

    SECTION("without orthogonalisation nothing is removed and nothing rescaled") {
        ConstantBathSystem sys;
        auto report = inqkit::WavePacket{}
                          .center(0.0, 0.0, 0.0).sigma(SIGMA).k0(0.0, 0.0, 0.0)
                          .inject_into_last_extra_state(sys.electrons, 1.0);
        CHECK_FALSE(report.orthogonalised);
        CHECK(report.removed_weight == 0.0);
        CHECK(report.sum_overlap_sq == 0.0);
        CHECK(report.norm_pre_renorm == Approx(report.norm_pre_ortho));
        // No renormalisation happens on this path, so norm_after is the raw
        // discrete Gaussian norm — close to 1, but not forced to it.
        CHECK(report.norm_after == Approx(report.norm_pre_ortho));
    }
}
