/*
 * inqkit::observables::slab_occupancy — the fraction of a single orbital's
 * density that lies inside a slab band, measured exactly on the grid.
 *
 *     f(t) = ∫_{|x_a − c| ≤ h} |ψ(r,t)|² dV   /   ∫ |ψ(r,t)|² dV
 *
 * with `a` the slab axis (0/1/2), `c` its centre and `h` its half-width, all in
 * Bohr. The denominator is the orbital norm at the same instant, so f is a pure
 * fraction in [0, 1] and is unaffected by any norm drift.
 *
 * WHY THIS EXISTS (docs/plans/slab-ks-orbital-stopping-wrap.md §3, Window B).
 * ------------------------------------------------------------------------
 * Electronic stopping power is a force: energy lost per unit path travelled
 * INSIDE the stopping medium. For a compact classical projectile the path inside
 * the medium is just the trajectory clipped to the slab, and the distinction is
 * invisible. For a WAVEPACKET it is not: a Gaussian of width σ spreads as
 * σ_d(t) = √(σ²/2 + t²/(2σ²)), so over a long CAP-free run the orbital ends up
 * wider than the slab — simultaneously inside and outside it. A stopping power
 * fitted against the centroid path then silently averages the drag over slab AND
 * vacuum and under-reports it by the filling factor.
 *
 * f(t) repairs this. Only the in-slab fraction of the orbital feels drag, so
 *
 *     dT/dt = −F · v · f     and     ds_slab/dt = f · v
 *     ⇒  −dT/ds_slab = F = S_slab      with     s_slab(t) = ∫₀ᵗ f ⟨p_z⟩/m dt'
 *
 * i.e. −dT/ds_slab is the stopping power per Bohr of path travelled inside the
 * slab, in BOTH the localised regime (f → 1, reducing to the ordinary −dT/ds)
 * and the fully delocalised one (f → 2h/L, applying the geometric filling factor
 * automatically). Measuring f on the grid rather than modelling it from a
 * Gaussian ansatz matters because the packet stops being Gaussian as soon as it
 * scatters.
 *
 * Cost: ONE extra reduction over the local grid per call — negligible against a
 * propagation step that touches every orbital.
 *
 * PERIODICITY. The band test is applied to the MINIMUM-IMAGE displacement from
 * the slab centre (fractional coordinate wrapped into [−0.5, 0.5), the same
 * window as systems::cell::position_in_cell), so a slab centred at the cell
 * centre is tested correctly even for density that has wrapped around the far
 * face — which, on a plain 3-D FFT basis, it always eventually has.
 *
 * No inq/ or inq-study/ edit — wrapper-only, same reduction pattern as
 * wp_real_space_stats.hpp.
 */
#ifndef INQKIT__OBSERVABLES__SLAB_OCCUPANCY
#define INQKIT__OBSERVABLES__SLAB_OCCUPANCY

#include <inq/inq.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

namespace inqkit {
namespace observables {

struct SlabOccupancy {
	double norm_in;    // ∫ over the band  |ψ|² dV
	double norm_total; // ∫ over the cell  |ψ|² dV   (the orbital norm)
	double fraction;   // norm_in / norm_total
};

// Geometry of the band, in Bohr. Defaults match the localised-jellium slab
// (25 Bohr thick, centred on z = 0).
struct SlabBand {
	int    axis       = 2;
	double center     = 0.0;
	double half_width = 12.5;
};

// Occupancy of `band` by orbital `state_index` of a gamma-only electrons object.
inline SlabOccupancy
slab_occupancy(inq::systems::electrons const & electrons, int state_index,
               SlabBand band = {}) {
	using namespace inq;

	if(electrons.kpin_size() != 1)
		throw std::runtime_error("slab_occupancy: only single-kpoint (gamma-only) runs are supported.");
	if(band.axis < 0 or band.axis > 2)
		throw std::runtime_error("slab_occupancy: axis must be 0, 1 or 2.");
	if(not (band.half_width > 0.0))
		throw std::runtime_error("slab_occupancy: half_width must be positive.");

	auto const & phi   = electrons.kpin()[0];
	auto const & basis = phi.basis();
	const double dV    = basis.volume_element();

	// Is the requested orbital local to this rank's set partition?
	const long st_start = phi.set_part().start();
	const long st_size  = phi.set_part().local_size();
	const bool local    = (state_index >= st_start and state_index < st_start + st_size);
	const int  ist_l    = local ? static_cast<int>(state_index - st_start) : 0;

	double sum_in = 0.0, sum_all = 0.0;

	if(local) {
		auto const sizes = basis.local_sizes();
		auto phic        = begin(phi.hypercubic());
		auto point_op    = basis.point_op();
		const int    ax  = band.axis;
		const double c   = band.center;
		const double h   = band.half_width;

		// Two sums in one pass: (in-band weight, total weight).
		auto pair = gpu::run(
			gpu::reduce(sizes[2]), gpu::reduce(sizes[1]), gpu::reduce(sizes[0]),
			inq::vector3<double>{0.0, 0.0, 0.0},
			[phic, ist_l, dV, point_op, ax, c, h] GPU_LAMBDA (auto iz, auto iy, auto ix) {
				auto v = phic[ix][iy][iz][ist_l];
				double w = dV * (inq::real(v) * inq::real(v) + inq::imag(v) * inq::imag(v));

				// minimum-image displacement from the band centre along `ax`
				auto r = point_op.rvector_cartesian(ix, iy, iz);
				inq::vector3<double> d{r[0], r[1], r[2]};
				d[ax] -= c;
				auto f = point_op.cell().to_contravariant(d);
				for(int idir = 0; idir < 3; idir++) {
					f[idir] -= floor(f[idir]);
					if(f[idir] >= 0.5) f[idir] -= 1.0;
				}
				auto dc = point_op.cell().to_cartesian(f);

				const double inband = (fabs(dc[ax]) <= h) ? w : 0.0;
				return inq::vector3<double>{inband, w, 0.0};
			});
		sum_in  = pair[0];
		sum_all = pair[1];
	}

	double host_buf[2] = {sum_in, sum_all};
	if(basis.comm().size() > 1)
		basis.comm().all_reduce_in_place_n(host_buf, 2, std::plus<>{});
	if(phi.set_comm().size() > 1)
		phi.set_comm().all_reduce_in_place_n(host_buf, 2, std::plus<>{});

	if(not (host_buf[1] > 0.0))
		throw std::runtime_error("slab_occupancy: non-positive norm for orbital "
		                         + std::to_string(state_index) + ".");

	SlabOccupancy out;
	out.norm_in    = host_buf[0];
	out.norm_total = host_buf[1];
	out.fraction   = host_buf[0] / host_buf[1];
	return out;
}

} // namespace observables
} // namespace inqkit

#endif
