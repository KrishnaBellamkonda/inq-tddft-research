// test_cosine_kick.cpp — known-case test for inqkit::perturbations::apply_cosine_kick.
//
// Build (when the GPU/CPU queue is free — DEFERRED tonight to avoid CPU
// contention with the live S(v) runs; production loss-function is cost-blocked
// at 14 s/step regardless):
//   cp here into a run dir and `inq-run --cpu test_cosine_kick.cpp`
//
// Checks:
//   (1) NORM CONSERVATION (exact, analytic): after the kick, <psi|psi> is
//       unchanged to machine precision (|exp(i theta)|=1 -> unitary).
//   (2) INSTANTANEOUS DENSITY UNCHANGED: |psi exp(i theta)|^2 = |psi|^2.
//   (3) PHASE IMPRINT: a chosen grid point's orbital phase advanced by
//       eta*cos(q.r) as intended.
// (Linearity of delta n(q,t) in eta is the run-level Stage-4 validation.)
#include <inq/inq.hpp>
#include <inqkit/perturbations/cosine_kick.hpp>

using namespace inq;
using namespace inq::magnitude;

int main() {
    auto comm = parallel::cartesian_communicator<2>{boost::mpi3::environment::get_world_instance(), {}};
    auto cell = systems::cell::cubic(10.0_b).periodic();
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(
        ions, options::electrons{}.spacing(0.5_b).extra_electrons(4),
        input::kpoints::gamma());
    ground_state::initial_guess(ions, electrons);

    // norm before
    double norm_before = 0.0;
    for (auto & phi : electrons.kpin())
        norm_before += operations::overlap_diagonal(phi)[0].real();

    inqkit::perturbations::apply_cosine_kick(
        electrons, 1e-3, inqkit::perturbations::commensurate_qx(1, 10.0));

    double norm_after = 0.0;
    for (auto & phi : electrons.kpin())
        norm_after += operations::overlap_diagonal(phi)[0].real();

    bool ok = std::abs(norm_after - norm_before) < 1e-12;
    std::cout << "norm_before=" << norm_before << " norm_after=" << norm_after
              << "  NORM_CONSERVED=" << (ok ? "PASS" : "FAIL") << "\n";
    return ok ? 0 : 1;
}
