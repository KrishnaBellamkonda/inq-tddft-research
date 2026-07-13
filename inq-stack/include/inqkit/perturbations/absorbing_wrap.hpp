/* -*- indent-tabs-mode: t -*- */

#ifndef INQKIT__PERTURBATIONS__ABSORBING_WRAP
#define INQKIT__PERTURBATIONS__ABSORBING_WRAP

// inqkit: unified wrap-around complex absorbing potential (CAP).
//
// A SINGLE smooth imaginary bump centred on the periodic z-boundary plane
// (z = ±L/2, which is one plane under periodic boundary conditions), in
// contrast to inq::perturbations::absorbing, whose sin² bump peaks at its
// mid-point and falls to ZERO at the cell edge — i.e. the standard two-sided
// arrangement (two bumps at ±z_mid) has W = 0 exactly where density crosses
// the periodic boundary.
//
// Profile (fractional cell coordinates, matching inq's absorbing.hpp which
// compares point_op.rvector()[2] ∈ [-0.5, 0.5) against fractional bounds):
//
//     d(z)  = 0.5 − |z_frac|                (periodic distance to the boundary)
//     W(z)  = |amplitude| · cos²(π·d/w)      for d < w/2, else 0
//
// so W peaks at the boundary plane (d = 0) and falls smoothly to zero a
// half-width inside the cell on both sides. By construction W is continuous
// and smooth ACROSS the wrap (it is a function of the periodic distance d).
//
// Equal-strength comparison with the two-sided CAP used in the localised
// jellium runs (two sin² bumps of full width w2 = 15 Bohr peaking at
// ±32.5 Bohr, each integrating to |η|·w2/2): a wrap CAP of full width
// w = 2·w2 = 30 Bohr has the SAME footprint (|z| > 25 Bohr in an 80 Bohr
// cell) and the SAME integral ∫W dz = |η|·w2 — the two differ only in
// topology (peak at the boundary vs zero at the boundary).
//
// Adapted from inq/src/perturbations/absorbing.hpp (Mozilla Public License
// 2.0, LLNL); credit: Andrade, Correa, Yao. New code lives in inqkit per the
// inq-immutable rule.

#include <inq_config.h>

#include <math/vector3.hpp>
#include <magnitude/energy.hpp>
#include <perturbations/none.hpp>

#include <cassert>
#include <cmath>

namespace inqkit {
namespace perturbations {

class absorbing_wrap : public inq::perturbations::none {

public:
	absorbing_wrap(inq::quantity<inq::magnitude::energy> amplitude, double width_frac):
		amplitude_(amplitude.in_atomic_units()),
		width_(width_frac)
	{
		assert(width_ > 0.0 and width_ <= 1.0);
	}

	auto has_potential() const {
		return true;
	}

	template<typename PotentialType>
	void potential(const double /*time*/, PotentialType & potential) const {

		using inq::complex;

		gpu::run(potential.basis().local_sizes()[2], potential.basis().local_sizes()[1], potential.basis().local_sizes()[0],
						 [point_op = potential.basis().point_op(), vk = begin(potential.cubic()), width = width_, amplitude = amplitude_] GPU_LAMBDA (auto iz, auto iy, auto ix) {
							 auto rr = point_op.rvector(ix, iy, iz);
							 auto dd = 0.5 - fabs(rr[2]);        // periodic distance to the z boundary plane
							 if (dd < width/2) {
								 vk[ix][iy][iz] += complex(0.0, amplitude*pow(cos(dd*M_PI/width), 2));
							 }
						 });
	}

private:
	double amplitude_;
	double width_;
};

}
}
#endif
