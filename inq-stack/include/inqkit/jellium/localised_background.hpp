/*
 * inqkit::jellium::localised_background — build a confined positive jellium
 * background charge density n₊(r) on an INQ real-space grid.
 *
 * Physics
 * -------
 * A localised jellium replaces INQ's implicit whole-cell uniform compensating
 * background (the dropped G=0 term of the periodic Poisson solve) with an
 * EXPLICIT positive charge confined to a region:
 *
 *     slab:   n₊(r) = n₀ · mask(|z − z₀| ; a)
 *     sphere: n₊(r) = n₀ · mask(|r − r₀| ; R_cl)
 *     box:    n₊(r) = n₀ · Π_i mask(|x_i − x₀_i| ; h_i)
 *     annulus:n₊(r) = n₀ · mask(d ; R_out) · [1 − mask(d ; R_in)], with d the
 *             radial distance ⟂ the tube axis — a hollow periodic tube
 *             (positive background only for R_in < d < R_out, uniform along axis).
 *
 * with n₀ = 3/(4π r_s³) the interior density. Charge neutrality is the CALLER's
 * responsibility: set electrons = round(∫n₊) and set n₀ = N/V_inside so that
 * ∫n₊ = N exactly (this is what makes the G=0 cancellation in the Hartree +
 * background Poisson solves exact — see docs/notes/localised-jellium-theory.md
 * Part 2.4).
 *
 * This header only BUILDS n₊. The electrostatic well v_bg = −poisson(n₊) and its
 * injection into the Kohn–Sham potential live in
 * inqkit/jellium/background_perturbation.hpp.
 *
 * Edge profile
 * ------------
 * edge_width == 0  → sharp Heaviside Θ (may cause Gibbs ringing in v_bg on a
 *                    finite grid; see worksheet Part 8).
 * edge_width  > 0  → complementary-error-function softening of width w:
 *                    mask(d; R) = ½ erfc((d − R)/w), a smooth 1→0 step centred
 *                    on the nominal edge R.
 *
 * Conventions: atomic units / Bohr. INQ centres the cell at the origin, so
 * cartesian coordinates run over [−L/2, +L/2]; rvector_cartesian(ix,iy,iz)
 * returns them directly.
 *
 * Status: UNVERIFIED — pending T0 host unit test
 * (inq-stack/tests/include/inqkit/jellium/test_localised_background.cpp):
 *   ∫n₊ = N, interior n₊ = n₀, edge profile correct.
 */

#ifndef INQKIT__JELLIUM__LOCALISED_BACKGROUND
#define INQKIT__JELLIUM__LOCALISED_BACKGROUND

#include <inq/inq.hpp>

#include <cmath>

namespace inqkit {
namespace jellium {

enum class background_shape { slab, sphere, box, annulus };

// All lengths in Bohr, density in a₀⁻³. center is the cartesian centre r₀.
// For a slab, only half_width (along the slab axis) and slab_axis are used.
// For a sphere, only half_width (= R_cl) is used. For a box, box_half is used.
struct localised_background_params {
	background_shape shape   = background_shape::slab;
	double           n0      = 0.0;             // interior density n₀
	inq::vector3<double> center{0.0, 0.0, 0.0}; // r₀ (Bohr, cartesian)
	double           half_width = 0.0;          // slab half-thickness a / sphere R_cl / annulus R_out
	double           inner_radius = 0.0;        // annulus R_in (hollow bore radius); unused by other shapes
	int              slab_axis  = 2;            // 0=x,1=y,2=z (slab/box normal, or annulus TUBE axis); slab confines this axis
	inq::vector3<double> box_half{0.0, 0.0, 0.0};// box half-extents (Bohr)
	double           edge_width = 0.0;          // 0 = sharp Θ; >0 = erfc softening width w
};

// ½ erfc((d−R)/w) for w>0, else hard Θ(R−d). Returns mask ∈ [0,1].
// TODO: Explain this function to me. Use a concrete case for d, R and w
// Need to use better variable names, d, R and w are very bad names. 
inline GPU_FUNCTION double background_mask(double d, double R, double w) {
	if(w <= 0.0) return (d < R) ? 1.0 : 0.0;
	return 0.5 * erfc((d - R) / w);
}

// Build n₊ on the given real-space basis. Templated on the basis type so it
// works for whatever real-space basis the KS potential field carries.
template <class Basis>
inq::basis::field<Basis, double>
make_localised_background(Basis const & basis, localised_background_params const & p) {

	inq::basis::field<Basis, double> nplus(basis);
	nplus.fill(0.0);

	auto point_op_ = basis.point_op();
	auto cub_      = begin(nplus.cubic());

	// Copy POD params into locals for the device lambda.
	const auto   shape = p.shape;
	const double n0    = p.n0;
	const auto   c     = p.center;
	const double a     = p.half_width;
	const double rin   = p.inner_radius;
	const int    sax   = p.slab_axis;
	const auto   bh    = p.box_half;
	const double w     = p.edge_width;

	gpu::run(basis.local_sizes()[2], basis.local_sizes()[1], basis.local_sizes()[0],
		[=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
			auto r = point_op_.rvector_cartesian(ix, iy, iz);
			double mask = 0.0;
			if(shape == background_shape::slab) {
				double d = fabs(r[sax] - c[sax]);
				mask = background_mask(d, a, w);
			} else if(shape == background_shape::sphere) {
				double dx = r[0]-c[0], dy = r[1]-c[1], dz = r[2]-c[2];
				double d  = sqrt(dx*dx + dy*dy + dz*dz);
				mask = background_mask(d, a, w);
			} else if(shape == background_shape::annulus) {
					// Hollow periodic tube: positive background between two
					// concentric cylinders R_in < d < R_out, axis = sax, uniform
					// along sax. d = radial distance in the plane PERPENDICULAR to
					// the tube axis (sum the two non-axis components).
					double drad2 = 0.0;
					for(int ax = 0; ax < 3; ++ax) {
						if(ax == sax) continue;
						double t = r[ax] - c[ax];
						drad2 += t*t;
					}
					double d = sqrt(drad2);
					// ½erfc((d−R_out)/w) · [1 − ½erfc((d−R_in)/w)]:
					// outer step (1 for d<R_out) times inner complement (0 for d<R_in)
					// ⇒ 1 inside the annulus, erfc-softened at both radial edges.
					mask = background_mask(d, a, w) * (1.0 - background_mask(d, rin, w));
				} else { // box: product of per-axis masks
				double mx = background_mask(fabs(r[0]-c[0]), bh[0], w);
				double my = background_mask(fabs(r[1]-c[1]), bh[1], w);
				double mz = background_mask(fabs(r[2]-c[2]), bh[2], w);
				mask = mx * my * mz;
			}
			cub_[ix][iy][iz] = n0 * mask;
		});

	return nplus;
}

} // namespace jellium
} // namespace inqkit

#endif
