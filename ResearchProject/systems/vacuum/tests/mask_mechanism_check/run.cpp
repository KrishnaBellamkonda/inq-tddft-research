// ============================================================================
// mask_mechanism_check: verify that applying a mask M(z) to the WP orbital
// INSIDE the per-step propagate callback (a) is a faithful no-op when M==1 and
// (b) feeds through into the subsequent ETRS step when M absorbs.
//
// This de-risks the Task-2 mechanism (Eq. 12, psi(t+dt) = M.U.psi(t)) implemented
// entirely in the inq-stack wrapper: the callback mutates the captured non-const
// `electrons` (the same object propagate() holds by reference), so the masked
// orbital is what the next step propagates. inq/ and inq-study stay untouched.
//
// Free particle: non_interacting theory, empty ions, ghost occupied via
// extra_electrons(2.0). WP injected into the single extra state, masked along z.
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>

#include <cmath>
#include <cstdio>
#include <filesystem>

using namespace inq;
using namespace inq::magnitude;
namespace fs = std::filesystem;

// Apply M(z) to the WP orbital in place on the GPU. Mirrors the injection loop
// in WavePacket::inject_into_last_extra_state (rvector_cartesian for correct
// centred-cell coordinates). z0 = absorber start, Labs = absorber width.
// kind: 0 -> M==1 everywhere (no-op fidelity), 1 -> sin^2 absorber (Eq. 13).
static void apply_mask_z(systems::electrons &electrons, double z0, double Labs,
                         int kind) {
  auto &phi = electrons.kpin()[0];
  auto &basis = phi.basis();
  int ist = phi.set_part().local_size() - 1; // WP slot
  auto phicub_ = begin(phi.hypercubic());
  auto point_op_ = basis.point_op();
  gpu::run(basis.local_sizes()[2], basis.local_sizes()[1],
           basis.local_sizes()[0],
           [=] GPU_LAMBDA(auto iz, auto iy, auto ix) {
             auto rvec = point_op_.rvector_cartesian(ix, iy, iz);
             double z = rvec[2];
             double M = 1.0;
             if (kind == 1) {
               if (z <= z0) {
                 M = 1.0;
               } else if (z < z0 + Labs) {
                 double s = std::sin(M_PI * (z - z0) / (2.0 * Labs));
                 M = 1.0 - s * s;
               } else {
                 M = 0.0;
               }
             }
             phicub_[ix][iy][iz][ist] = phicub_[ix][iy][iz][ist] * M;
           });
  INQKIT_GPU_SYNC();
}

// One propagation from a freshly injected WP; returns final (norm, z-centroid).
struct Result { double N, z; };

static Result run_case(systems::ions &ions, systems::electrons &electrons,
                       double sigma, double k0z, double z0_start, int n_steps,
                       double dt, const char *label, bool mask, int mask_kind,
                       double abs_z0, double abs_L) {
  // (re)inject a fresh WP into the extra state
  auto rep = inqkit::WavePacket{}
                 .center(0.0, 0.0, z0_start)
                 .sigma(sigma)
                 .k0(0.0, 0.0, k0z)
                 .inject_into_last_extra_state(electrons, 1.0);
  const int wp_idx = rep.state_index;

  fs::path dir = fs::temp_directory_path() / "mask_check";
  fs::create_directories(dir);
  inqkit::observables::WPRealSpaceStats wp_rs((dir / "rs.csv").string(), wp_idx);

  inqkit::observables::WPRealSpaceMoments mT{};
  real_time::propagate(
      ions, electrons,
      [&](auto const &data) {
        // mask the WP orbital AFTER each ETRS step -> Eq. 12 (M.U.psi)
        if (mask)
          apply_mask_z(electrons, abs_z0, abs_L, mask_kind);
        mT = wp_rs.compute(data.electrons()); // read-only stats (last step kept)
      },
      options::theory{}.non_interacting(),
      options::real_time{}.num_steps(n_steps).dt(dt * 1.0_atomictime));

  std::printf("  [%-10s] final norm = %.8f   z-centroid = %+.4f\n", label, mT.N,
              mT.z);
  return {mT.N, mT.z};
}

int main() {
  std::printf("\n=== mask_mechanism_check ===\n");

  const double L = 40.0;     // Bohr cubic
  const double SPACING = 0.5;
  const double SIGMA = 2.0;
  const double K0Z = 2.0;    // +z
  const double Z0 = -10.0;   // start
  const int N_STEPS = 200;
  const double DT = 0.05;    // T = 10 a.u. -> ballistic centroid Z0 + K0Z*T = +10

  const double ABS_Z0 = 0.0; // absorber occupies z in [0, 15]
  const double ABS_L = 15.0;

  auto cell =
      systems::cell::orthorhombic(L * 1.0_b, L * 1.0_b, L * 1.0_b).periodic();
  auto ions = systems::ions(cell);

  const double ec_ha = 0.5 * std::pow(M_PI / SPACING, 2.0);
  auto electrons = systems::electrons(
      ions, options::electrons{}.cutoff(ec_ha * 1.0_Ha).extra_states(1).extra_electrons(2.0));
  ground_state::initial_guess(ions, electrons);

  // A: no mask (baseline free propagation)
  auto A = run_case(ions, electrons, SIGMA, K0Z, Z0, N_STEPS, DT, "baseline",
                    false, 0, ABS_Z0, ABS_L);
  // B: M==1 every step (no-op fidelity: must equal baseline)
  auto B = run_case(ions, electrons, SIGMA, K0Z, Z0, N_STEPS, DT, "M==1", true,
                    0, ABS_Z0, ABS_L);
  // C: sin^2 absorber every step (feedthrough: norm must drop)
  auto C = run_case(ions, electrons, SIGMA, K0Z, Z0, N_STEPS, DT, "sin2-abs",
                    true, 1, ABS_Z0, ABS_L);

  std::printf("\n--- verdict ---\n");
  std::printf("no-op fidelity  |N_B - N_A| = %.2e   |z_B - z_A| = %.2e\n",
              std::abs(B.N - A.N), std::abs(B.z - A.z));
  std::printf("feedthrough     N_A = %.6f  ->  N_C = %.6f   (drop = %.4f)\n",
              A.N, C.N, A.N - C.N);
  bool fidelity = std::abs(B.N - A.N) < 1e-9 && std::abs(B.z - A.z) < 1e-9;
  bool feedthrough = (A.N - C.N) > 0.05;
  std::printf("FIDELITY %s   FEEDTHROUGH %s\n", fidelity ? "PASS" : "FAIL",
              feedthrough ? "PASS" : "FAIL");
  return (fidelity && feedthrough) ? 0 : 1;
}
