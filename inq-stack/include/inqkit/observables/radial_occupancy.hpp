/*
 * inqkit::observables::radial_occupancy — how a single KS orbital is distributed
 * in the TRANSVERSE radius about a tube axis: the fraction of |ψ|² inside the
 * hollow bore, inside the jellium wall, and outside the tube, plus ⟨r⊥⟩ and
 * σ_r⊥, all measured exactly on the grid.
 *
 *     r⊥(r) = |minimum-image displacement of r from the tube axis|
 *
 *     f_bore(t)    = ∫_{r⊥ <  R_in }  |ψ|² dV  /  ∫ |ψ|² dV
 *     f_wall(t)    = ∫_{R_in ≤ r⊥ < R_out} |ψ|² dV  /  ∫ |ψ|² dV
 *     f_outside(t) = 1 − f_bore − f_wall
 *     ⟨r⊥⟩(t)      = ∫ r⊥ |ψ|² dV / ∫ |ψ|² dV
 *     σ_r⊥(t)      = √(⟨r⊥²⟩ − ⟨r⊥⟩²)
 *
 * WHY THIS EXISTS (docs/plans/cylindrical-channeling-ks-stopping.md §3).
 * ---------------------------------------------------------------------
 * The channeling study rests on ONE physical premise: the wavepacket flies down
 * the hollow bore and couples to the jellium wall only through the smooth image
 * force, so the interaction-driven growth of var(p) that contaminated the bulk
 * KS-stopping measurement is suppressed. That premise is a CLAIM ABOUT WHERE THE
 * PACKET IS, and it is testable only by measuring where the packet is.
 *
 * A Gaussian of ψ-width σ_WP disperses as σ_d(t) = √(σ_WP²/2 + t²/(2σ_WP²)), so
 * a packet that starts comfortably inside a 10-Bohr bore ends up leaking into
 * the wall at large t. f_wall(t) is the leak, measured rather than modelled —
 * which matters because the packet stops being Gaussian the moment it scatters.
 * The KS-stopping fit window must be restricted to the interval over which
 * f_wall stays small; without this observable that window is a guess.
 *
 * It is also the honest replacement for slab_occupancy in a tube geometry. On a
 * z-UNIFORM tube the medium fills every z, so there is no in-medium path
 * fraction to correct for and the stopping power is −dT/ds against the full path
 * (unlike the slab, where s5 = ∫f·v dt was needed). What the tube needs instead
 * is the RADIAL question: is the projectile still channeling?
 *
 * Cost: two reductions over the local grid per call — negligible against a
 * propagation step that touches every orbital.
 *
 * PERIODICITY. The radius is built from the MINIMUM-IMAGE transverse
 * displacement (fractional coordinates wrapped into [−0.5, 0.5), the same window
 * as systems::cell::position_in_cell and slab_occupancy), so a packet that has
 * spread across a transverse face is measured correctly rather than being
 * assigned a spuriously large radius. NOTE that once the packet genuinely
 * overlaps its own transverse periodic image, r⊥ saturates at ~L⊥/2 and the
 * quantity stops meaning "distance from the axis of one packet" — that is a
 * property of the physical setup, not of this estimator, and the run should be
 * analysed only up to that time.
 *
 * AXES. `axis` is the tube axis (0/1/2); the two transverse directions are
 * (axis+1)%3 and (axis+2)%3, and `center` is the full 3-vector of the tube axis
 * position (its component along `axis` is ignored).
 *
 * No inq/ or inq-study/ edit — wrapper-only, same reduction pattern as
 * slab_occupancy.hpp / wp_real_space_stats.hpp.
 */
#ifndef INQKIT__OBSERVABLES__RADIAL_OCCUPANCY
#define INQKIT__OBSERVABLES__RADIAL_OCCUPANCY

#include <inq/inq.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

namespace inqkit {
namespace observables {

// Geometry of the tube, in Bohr. Defaults match the r_s = 3 channeling tube
// (bore R_in = 10, wall out to R_out = 14, axis ∥ z through the cell centre).
struct RadialBand {
	int                  axis     = 2;                  // tube axis: 0=x, 1=y, 2=z
	inq::vector3<double> center   = {0.0, 0.0, 0.0};    // point on the tube axis
	double               r_inner  = 10.0;               // bore radius R_in
	double               r_outer  = 14.0;               // wall outer radius R_out
};

struct RadialOccupancy {
	double norm_total   = 0.0;   // ∫ |ψ|² dV over the cell (the orbital norm)
	double norm_bore    = 0.0;   // ∫ over r⊥ < R_in
	double norm_wall    = 0.0;   // ∫ over R_in ≤ r⊥ < R_out
	double norm_outside = 0.0;   // ∫ over r⊥ ≥ R_out
	double f_bore       = 0.0;   // norm_bore    / norm_total
	double f_wall       = 0.0;   // norm_wall    / norm_total
	double f_outside    = 0.0;   // norm_outside / norm_total
	double r_mean       = 0.0;   // ⟨r⊥⟩   (Bohr)
	double r2_mean      = 0.0;   // ⟨r⊥²⟩  (Bohr²)
	double sigma_r      = 0.0;   // √(⟨r⊥²⟩ − ⟨r⊥⟩²)  (Bohr)
};

// Radial distribution of orbital `state_index` of a gamma-only electrons object.
inline RadialOccupancy
radial_occupancy(inq::systems::electrons const & electrons, int state_index,
                 RadialBand band = {}) {
	using namespace inq;

	if(electrons.kpin_size() != 1)
		throw std::runtime_error("radial_occupancy: only single-kpoint (gamma-only) runs are supported.");
	if(band.axis < 0 or band.axis > 2)
		throw std::runtime_error("radial_occupancy: axis must be 0, 1 or 2.");
	if(not (band.r_inner >= 0.0))
		throw std::runtime_error("radial_occupancy: r_inner must be non-negative.");
	if(not (band.r_outer > band.r_inner))
		throw std::runtime_error("radial_occupancy: r_outer must exceed r_inner.");

	auto const & phi   = electrons.kpin()[0];
	auto const & basis = phi.basis();
	const double dV    = basis.volume_element();

	// Is the requested orbital local to this rank's set partition?
	const long st_start = phi.set_part().start();
	const long st_size  = phi.set_part().local_size();
	const bool local    = (state_index >= st_start and state_index < st_start + st_size);
	const int  ist_l    = local ? static_cast<int>(state_index - st_start) : 0;

	double w_bore = 0.0, w_wall = 0.0, w_all = 0.0, w_r = 0.0, w_r2 = 0.0;

	if(local) {
		auto const sizes = basis.local_sizes();
		auto phic        = begin(phi.hypercubic());
		auto point_op    = basis.point_op();
		const int    ax  = band.axis;
		const int    ia  = (band.axis + 1) % 3;
		const int    ib  = (band.axis + 2) % 3;
		const auto   c   = band.center;
		const double rin = band.r_inner;
		const double rout = band.r_outer;

		// Pass 1: (in-bore weight, in-wall weight, total weight).
		auto shells = gpu::run(
			gpu::reduce(sizes[2]), gpu::reduce(sizes[1]), gpu::reduce(sizes[0]),
			inq::vector3<double>{0.0, 0.0, 0.0},
			[phic, ist_l, dV, point_op, ax, ia, ib, c, rin, rout] GPU_LAMBDA (auto iz, auto iy, auto ix) {
				auto v = phic[ix][iy][iz][ist_l];
				const double w = dV * (inq::real(v) * inq::real(v) + inq::imag(v) * inq::imag(v));

				// minimum-image displacement from the tube axis
				auto r = point_op.rvector_cartesian(ix, iy, iz);
				inq::vector3<double> d{r[0] - c[0], r[1] - c[1], r[2] - c[2]};
				d[ax] = 0.0;                                   // project onto the transverse plane
				auto f = point_op.cell().to_contravariant(d);
				for(int idir = 0; idir < 3; idir++) {
					f[idir] -= floor(f[idir]);
					if(f[idir] >= 0.5) f[idir] -= 1.0;
				}
				auto dc = point_op.cell().to_cartesian(f);
				const double rperp = sqrt(dc[ia] * dc[ia] + dc[ib] * dc[ib]);

				const double in_bore = (rperp <  rin)                       ? w : 0.0;
				const double in_wall = (rperp >= rin and rperp < rout)      ? w : 0.0;
				return inq::vector3<double>{in_bore, in_wall, w};
			});

		// Pass 2: (Σ r⊥ w, Σ r⊥² w, unused).
		auto moments = gpu::run(
			gpu::reduce(sizes[2]), gpu::reduce(sizes[1]), gpu::reduce(sizes[0]),
			inq::vector3<double>{0.0, 0.0, 0.0},
			[phic, ist_l, dV, point_op, ax, ia, ib, c] GPU_LAMBDA (auto iz, auto iy, auto ix) {
				auto v = phic[ix][iy][iz][ist_l];
				const double w = dV * (inq::real(v) * inq::real(v) + inq::imag(v) * inq::imag(v));

				auto r = point_op.rvector_cartesian(ix, iy, iz);
				inq::vector3<double> d{r[0] - c[0], r[1] - c[1], r[2] - c[2]};
				d[ax] = 0.0;
				auto f = point_op.cell().to_contravariant(d);
				for(int idir = 0; idir < 3; idir++) {
					f[idir] -= floor(f[idir]);
					if(f[idir] >= 0.5) f[idir] -= 1.0;
				}
				auto dc = point_op.cell().to_cartesian(f);
				const double rperp = sqrt(dc[ia] * dc[ia] + dc[ib] * dc[ib]);

				return inq::vector3<double>{rperp * w, rperp * rperp * w, 0.0};
			});

		w_bore = shells[0];
		w_wall = shells[1];
		w_all  = shells[2];
		w_r    = moments[0];
		w_r2   = moments[1];
	}

	double host_buf[5] = {w_bore, w_wall, w_all, w_r, w_r2};
	if(basis.comm().size() > 1)
		basis.comm().all_reduce_in_place_n(host_buf, 5, std::plus<>{});
	if(phi.set_comm().size() > 1)
		phi.set_comm().all_reduce_in_place_n(host_buf, 5, std::plus<>{});

	if(not (host_buf[2] > 0.0))
		throw std::runtime_error("radial_occupancy: non-positive norm for orbital "
		                         + std::to_string(state_index) + ".");

	RadialOccupancy out;
	out.norm_bore    = host_buf[0];
	out.norm_wall    = host_buf[1];
	out.norm_total   = host_buf[2];
	out.norm_outside = host_buf[2] - host_buf[0] - host_buf[1];
	out.f_bore       = out.norm_bore    / out.norm_total;
	out.f_wall       = out.norm_wall    / out.norm_total;
	out.f_outside    = out.norm_outside / out.norm_total;
	out.r_mean       = host_buf[3] / out.norm_total;
	out.r2_mean      = host_buf[4] / out.norm_total;
	// Clamp: ⟨r²⟩−⟨r⟩² is non-negative analytically but can go ~−1e−16 on the grid.
	const double var = out.r2_mean - out.r_mean * out.r_mean;
	out.sigma_r      = (var > 0.0) ? std::sqrt(var) : 0.0;
	return out;
}

} // namespace observables
} // namespace inqkit

#endif
