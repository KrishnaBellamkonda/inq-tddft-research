// inqkit/perturbations/cosine_kick.hpp
// ----------------------------------------------------------------------------
// Finite-q delta kick: at t=0 multiply every occupied Kohn-Sham orbital by
//
//     exp( i * eta * cos(q . r) )
//
// the impulse response (delta(t) of strength eta) of the external potential
// eta * cos(q . r). For a box-commensurate q = (2*pi*n/L, 0, 0) this seeds a
// finite-q density response delta n(q, t) whose Fourier transform gives the
// energy-loss function Im[-1/eps(q, omega)] (Stage 4 of the overnight plan).
//
// INQ's perturbations::kick only applies a UNIFORM phase exp(i k . r) (q=0
// dipole boost). This header reuses INQ's exact in-place orbital-mutation idiom
// (perturbations/kick.hpp::zero_step: gpu::run over the real-space grid,
// point_op().rvector, multiply phi.hypercubic()[iz][iy][ix][ist]) but with the
// finite-q cosine phase. No INQ core change required (the orbital_set is
// mutable), so the "duplicate INQ" fallback is not triggered.
//
// Correctness:
//   - Norm conservation is EXACT and analytic: |exp(i theta)| = 1, so the map
//     is unitary on each grid point and preserves <psi|psi> identically.
//   - The instantaneous density is unchanged (|psi exp(i theta)|^2 = |psi|^2);
//     the response delta n(q, t) develops under subsequent propagation.
//   - Linearity in eta (delta n(q,t) proportional to eta for small eta) is the
//     run-level validation in Stage 4.
//
// Usage (apply once, before real_time::propagate):
//   inqkit::perturbations::apply_cosine_kick(electrons, eta, {q, 0.0, 0.0});
// ----------------------------------------------------------------------------
#pragma once

#include <inq/inq.hpp>

namespace inqkit::perturbations {

// Apply exp(i*eta*cos(q.r)) in place to every KS orbital (all k-points).
template <typename Electrons>
void apply_cosine_kick(Electrons & electrons, double eta,
                       inq::vector3<double> qvec) {
    using inq::complex;
    for (auto & phi : electrons.kpin()) {
        inq::gpu::run(
            phi.basis().local_sizes()[2],
            phi.basis().local_sizes()[1],
            phi.basis().local_sizes()[0],
            [pop = phi.basis().point_op(),
             ph  = begin(phi.hypercubic()),
             qvec, eta,
             nst = phi.set_part().local_size()] GPU_LAMBDA(auto iz, auto iy, auto ix) {
                auto rr     = pop.rvector(ix, iy, iz);
                auto factor = exp(complex(0.0, eta * cos(dot(qvec, rr))));
                for (int ist = 0; ist < nst; ist++)
                    ph[ix][iy][iz][ist] *= factor;
            });
    }
}

// Box-commensurate axial wavevector q_n = 2*pi*n/L along +x.
inline inq::vector3<double> commensurate_qx(int n, double L_bohr) {
    return {2.0 * M_PI * n / L_bohr, 0.0, 0.0};
}

}  // namespace inqkit::perturbations
