// Engine-tier T28/T04 (USER-LOCKED): cross-validate the WP mean momentum ⟨p⟩
// computed two INDEPENDENT ways on the same injected WP, and check both ≈ k₀.
//
//   Route 1 (real space): ⟨p_d⟩ = Im Σ conj(ψ) ∂_d ψ / Σ|ψ|²  via a CENTRAL
//            FINITE DIFFERENCE on the extracted real-space WP field (no FFT).
//   Route 2 (reciprocal): ⟨p_d⟩ = Σ k_d |ψ̃|² / Σ|ψ̃|²  via to_fourier +
//            point_op.gvector_cartesian — exactly wp_momentum_stats' method.
//
// Agreement of two independent routes validates: k is in Bohr⁻¹ (T28), the full
// complex ψ is used (T04), and the reciprocal momentum machinery is correct.
// WP injected WITHOUT orthogonalisation so ⟨p⟩ is a clean k₀ (the before/after-
// ortho shift is a separate experiment, parked with T26).

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/detail/grid_layout.hpp>

#include <complex>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;
using inqkit::detail::grid_layout::flatten_index;

TEST_CASE("T28/T04: ⟨p⟩ real-space vs reciprocal routes agree and equal k₀", "[observables][wp_momentum][engine]") {
  const double K0 = 1.0;  // Bohr^-1, along x

  systems::ions ions(systems::cell::cubic(12.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(2));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  auto report = inqkit::WavePacket{}
                    .center(0.0, 0.0, 0.0)
                    .sigma(1.5)
                    .k0(K0, 0.0, 0.0)
                    .inject_into_last_extra_state(electrons, 1.0);
  const int wp = report.state_index;

  // ---- Route 1: real-space central finite difference on the WP field ----
  auto psi = inqkit::fields::orbital::wavefunction(electrons, wp);
  long double num_x = 0.0L, norm = 0.0L;
  for (int ix = 1; ix < psi.nx - 1; ++ix)
    for (int iy = 0; iy < psi.ny; ++iy)
      for (int iz = 0; iz < psi.nz; ++iz) {
        auto c  = psi.values[flatten_index(ix, iy, iz, psi.ny, psi.nz)];
        auto cp = psi.values[flatten_index(ix + 1, iy, iz, psi.ny, psi.nz)];
        auto cm = psi.values[flatten_index(ix - 1, iy, iz, psi.ny, psi.nz)];
        std::complex<double> dpsi = (cp - cm) / (2.0 * psi.dx_bohr);
        num_x += std::imag(std::conj(c) * dpsi);   // Im(ψ* ∂ψ) = Re(ψ*(-i∂)ψ)
        norm  += std::norm(c);
      }
  const double px_real = static_cast<double>(num_x / norm);

  // ---- Route 2: reciprocal space via to_fourier (wp_momentum_stats method) ----
  auto fphi = inq::operations::transform::to_fourier(electrons.kpin()[0]);
  auto const &fbasis = fphi.basis();
  auto const sizes = fbasis.local_sizes();
  auto fhc = begin(fphi.hypercubic());
  auto point_op = fbasis.point_op();
  const int ist_l = wp;

  double sum_n = gpu::run(
      gpu::reduce(sizes[2]), gpu::reduce(sizes[1]), gpu::reduce(sizes[0]), 0.0,
      [fhc, ist_l] GPU_LAMBDA(auto iz, auto iy, auto ix) {
        auto v = fhc[ix][iy][iz][ist_l];
        return inq::real(v) * inq::real(v) + inq::imag(v) * inq::imag(v);
      });
  inq::vector3<double> sum_p = gpu::run(
      gpu::reduce(sizes[2]), gpu::reduce(sizes[1]), gpu::reduce(sizes[0]),
      inq::vector3<double>{0.0, 0.0, 0.0},
      [fhc, ist_l, point_op] GPU_LAMBDA(auto iz, auto iy, auto ix) {
        auto v = fhc[ix][iy][iz][ist_l];
        double w = inq::real(v) * inq::real(v) + inq::imag(v) * inq::imag(v);
        auto k = point_op.gvector_cartesian(ix, iy, iz);
        return inq::vector3<double>{k[0] * w, k[1] * w, k[2] * w};
      });
  const double px_recip = sum_p[0] / sum_n;
  const double py_recip = sum_p[1] / sum_n;
  const double pz_recip = sum_p[2] / sum_n;

  // ---- The cross-validation: two independent routes must agree, both ≈ k₀ ----
  // The reciprocal route is spectrally EXACT (px_recip ≈ k₀ to ~1e-8). The
  // real-space CENTRAL DIFFERENCE underestimates by sin(k·dx)/(k·dx) ≈ 0.96 at
  // k·dx = 0.5, so it lands a few % low — that residual is the FD discretization
  // error, not a units/correctness disagreement. Tolerance set accordingly.
  CHECK(px_recip == Approx(K0).margin(0.02));         // reciprocal: tight (exact)
  CHECK(px_real  == Approx(K0).margin(0.12));         // FD: O((k·dx)²) low bias
  CHECK(px_real  == Approx(px_recip).margin(0.12));   // routes agree within FD error
  // Transverse momenta vanish.
  CHECK(py_recip == Approx(0.0).margin(0.05));
  CHECK(pz_recip == Approx(0.0).margin(0.05));
}
