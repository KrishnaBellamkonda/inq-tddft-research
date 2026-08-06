/*
 * inqkit::fields::from_inq_field — export an ARBITRARY INQ real-space field as a
 * RealField3D, so anything the engine can build can be written as a VTI.
 *
 * WHY THIS EXISTS. Until now the only routes into RealField3D were
 * fields::density::total(electrons) and fields::orbital(electrons, i) — both of
 * which take an `electrons` object. Every other real-space field the wrapper
 * builds (the positive jellium background n₊, the Poisson potential φ₊, a
 * projectile charge blob n_proj, a drag field φ_e − φ₊) is a plain
 * inq::basis::field<real_space, double> and had NO way to be visualised at all.
 * That gap was found the hard way: the channeling-twin ground state tried to dump
 * its own annular background beside the electron density — so a reader could SEE
 * the bore/wall boundary instead of taking R_in/R_out on trust — and did not
 * compile (2026-08-01).
 *
 * CONVENTION — this is the whole point of routing it through one function.
 * INQ stores real-space fields in FFT-NATURAL order (index 0 = the origin, the
 * negative half wrapped to the top). RealField3DWriter stamps
 * Origin = −L/2 and expects PHYSICAL order (index 0 = the left edge). The
 * conversion is the fft_shift below, and getting it wrong produces a picture that
 * looks entirely plausible with the centre and the edges swapped — the recurring
 * bug .claude/rules/vti-coordinate-mapping.md exists to prevent. Doing it in one
 * place means a run script can never re-derive it wrongly.
 *
 * The shift, the origin and the spacing are byte-identical to what
 * fields::density::total does; the two SHOULD eventually share this
 * implementation, and density::total is deliberately left untouched for now so
 * this addition cannot perturb any existing run.
 *
 * Single-rank only, like the rest of the export layer (throws otherwise).
 *
 * No inq/ or inq-study/ edit — wrapper-only.
 */
#ifndef INQKIT__FIELDS__INQ_FIELD
#define INQKIT__FIELDS__INQ_FIELD

#include <inq/inq.hpp>
#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/fields/real_field_3d.hpp>

#include <cstddef>
#include <stdexcept>

namespace inqkit {
namespace fields {

// Copy an INQ real-space scalar field to a host RealField3D in PHYSICAL order.
template <class Field>
RealField3D from_inq_field(Field const & f) {
	auto const & basis = f.basis();

	if(basis.comm().size() != 1)
		throw std::runtime_error("inqkit::fields::from_inq_field: multi-rank basis "
		                         "export is not implemented yet.");

	auto const nx = basis.sizes()[0];
	auto const ny = basis.sizes()[1];
	auto const nz = basis.sizes()[2];
	auto const spacing = basis.rspacing();

	RealField3D out;
	out.nx = nx; out.ny = ny; out.nz = nz;
	out.origin_x_bohr = basis.symmetric_range_begin(0) * spacing[0];
	out.origin_y_bohr = basis.symmetric_range_begin(1) * spacing[1];
	out.origin_z_bohr = basis.symmetric_range_begin(2) * spacing[2];
	out.dx_bohr = spacing[0];
	out.dy_bohr = spacing[1];
	out.dz_bohr = spacing[2];
	out.values.resize(static_cast<std::size_t>(nx) * ny * nz);

	// ONE bulk GPU->host copy. Element-wise access on a device-resident cubic
	// view costs a separate fetch per cell (~30 min on a 4.7M-point grid — the
	// lesson already paid for in density::total).
	boost::multi::array<double, 3> host{f.cubic()};

	// ix runs left-to-right in PHYSICAL space (ix = 0 at −L/2); INQ's storage is
	// FFT-natural, so read host[fft_shift_index(ix)].
	for(int ix = 0; ix < nx; ++ix) {
		const int sx = inqkit::detail::grid_layout::fft_shift_index(ix, nx);
		for(int iy = 0; iy < ny; ++iy) {
			const int sy = inqkit::detail::grid_layout::fft_shift_index(iy, ny);
			for(int iz = 0; iz < nz; ++iz) {
				const int sz = inqkit::detail::grid_layout::fft_shift_index(iz, nz);
				const auto flat =
					inqkit::detail::grid_layout::flatten_index(ix, iy, iz, ny, nz);
				out.values[flat] = host[sx][sy][sz];
			}
		}
	}
	return out;
}

} // namespace fields
} // namespace inqkit

#endif
