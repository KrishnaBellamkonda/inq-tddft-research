// ============================================================================
// inqkit/observables/center_of_density.hpp
//
// Centre of density of a real-valued density field f(r) >= 0:
//
//     <r> = (integral r f(r) dV) / (integral f(r) dV)
//
// For a wave-packet orbital this gives the WP centroid <r>_wp(t), which we
// track over time to derive the WP trajectory along the propagation axis
// (z by default; project onto k_hat for tilted launches in postprocess).
//
// Units: bohr (matches RealField3D origin/spacing convention).
// ============================================================================
#pragma once

#include <inqkit/fields/real_field_3d.hpp>
#include <inqkit/detail/grid_layout.hpp>

#include <array>
#include <cstddef>
#include <stdexcept>

namespace inqkit::observables {

struct CenterOfDensityResult {
    double x_bohr = 0.0;
    double y_bohr = 0.0;
    double z_bohr = 0.0;
    double total_weight = 0.0;  // integral of f dV (normalisation diagnostic)
};

inline CenterOfDensityResult
center_of_density(inqkit::fields::RealField3D const& f) {
    if (f.empty()) {
        throw std::runtime_error(
            "inqkit::observables::center_of_density: empty field");
    }

    const double dV = f.dx_bohr * f.dy_bohr * f.dz_bohr;
    long double w = 0.0L, mx = 0.0L, my = 0.0L, mz = 0.0L;

    for (int ix = 0; ix < f.nx; ++ix) {
        const double x = f.origin_x_bohr + (ix + 0.5) * f.dx_bohr;
        for (int iy = 0; iy < f.ny; ++iy) {
            const double y = f.origin_y_bohr + (iy + 0.5) * f.dy_bohr;
            for (int iz = 0; iz < f.nz; ++iz) {
                const double z = f.origin_z_bohr + (iz + 0.5) * f.dz_bohr;
                const auto flat =
                    inqkit::detail::grid_layout::flatten_index(
                        ix, iy, iz, f.ny, f.nz);
                const double v = f.values[flat];
                w  += v;
                mx += v * x;
                my += v * y;
                mz += v * z;
            }
        }
    }

    CenterOfDensityResult r;
    r.total_weight = static_cast<double>(w * dV);
    if (w > 0.0L) {
        r.x_bohr = static_cast<double>(mx / w);
        r.y_bohr = static_cast<double>(my / w);
        r.z_bohr = static_cast<double>(mz / w);
    }
    return r;
}

} // namespace inqkit::observables
