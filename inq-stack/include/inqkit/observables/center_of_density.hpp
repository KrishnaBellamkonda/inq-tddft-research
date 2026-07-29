/*
 * Computes the centre of density (centroid) of a real-valued scalar field
 * f(r) ≥ 0, defined as the density-weighted mean position:
 *
 *   <r> = ( ∫ r f(r) dV ) / ( ∫ f(r) dV )
 *
 * For a wave-packet orbital this yields the WP centroid <r>_wp(t), which
 * can be tracked over time to reconstruct the WP trajectory along the
 * propagation axis. For tilted launches, project the returned centroid
 * onto k̂ in post-processing to extract the along-beam displacement.
 *
 * Units
 * -----
 * All coordinates are in Bohr, matching the origin and spacing convention
 * of RealField3D.
 *
 * Usage
 * -----
 *   auto cod = inqkit::observables::center_of_density(density);
 *   auto centre = cod.center_bohr;  // inqkit::detail::Vec3 {x, y, z} in Bohr
 */
#pragma once

#include <inqkit/fields/real_field_3d.hpp>
#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/detail/vec3.hpp>

#include <array>
#include <cstddef>
#include <stdexcept>

namespace inqkit::observables {


// Centroid of a non-negative scalar field, stored as a Vec3 unit (T07/T09)
// rather than loose x/y/z scalars, plus the integral used as the normaliser.
struct CenterOfDensityResult {
    inqkit::detail::Vec3 center_bohr;  // density-weighted centroid <r> (Bohr)
    double total_weight = 0.0;         // integral of f dV (normalisation diagnostic)
};

/* TODO: In this function, the f variable should be renamed to field for better
   readability. For similar reasons w might be renamed as weight.  
   What does m mean in mx, my and mz? Is there a better way to name this? is this
   moment?
*/
inline CenterOfDensityResult
center_of_density(inqkit::fields::RealField3D const& f) {
    if (f.empty()) {
        throw std::runtime_error(
            "inqkit::observables::center_of_density: empty field");
    }

    const double dV = f.dx_bohr * f.dy_bohr * f.dz_bohr;
    // TODO: What unit does the L suffix represent? Shouldn't this be in Bohr?
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
                /*
                * Calculating the weight and the total center of density. 
                */ 

                w  += v;
                mx += v * x;
                my += v * y;
                mz += v * z;
            }
        }
    }

    CenterOfDensityResult r;
    r.total_weight = static_cast<double>(w * dV);
    // TODO: Explain this condition, and why this was put in. Then, can write this
    // as a succinct comment. 
    if (w > 0.0L) {
        r.center_bohr.x = static_cast<double>(mx / w);
        r.center_bohr.y = static_cast<double>(my / w);
        r.center_bohr.z = static_cast<double>(mz / w);
    }
    return r;
}

} // namespace inqkit::observables
