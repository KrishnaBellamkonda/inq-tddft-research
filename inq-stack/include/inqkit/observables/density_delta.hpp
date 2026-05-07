// ============================================================================
// inqkit/observables/density_delta.hpp
//
// Local density fluctuation observable
//
//     dn(r, t) = n(r, t) - n(r, t0)
//
// where t0 is the user-supplied reference time (typically the first
// post-injection step). The header provides three views on dn:
//
//   1. Raw 3D field          -> VTI series (diverging colourmap downstream)
//   2. Coarse-grained 3D field -> VTI series binned into bin_size_bohr cubes
//                                (suppresses Friedel-oscillation noise; the
//                                macroscopic redistribution view).
//   3. Integrated L2 metric  -> sigma2_n(t) = integral |dn|^2 dV.
//
// All work is done on the host (the underlying RealField3D is host memory
// already, populated via inqkit::fields::density::total). Cost per step at
// 60^3 grid is sub-millisecond plus VTI write time.
// ============================================================================
#pragma once

#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/fields/real_field_3d.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>

namespace inqkit::observables {

struct DensityDeltaConfig {
    bool   emit_raw_vti     = true;
    bool   emit_coarse_vti  = true;
    bool   compute_l2       = true;
    double coarse_bin_bohr  = 3.0;
};

class DensityDelta {
public:
    DensityDelta(std::string raw_dir,
                 std::string coarse_dir,
                 DensityDeltaConfig cfg = {})
        : raw_dir_(std::move(raw_dir)),
          coarse_dir_(std::move(coarse_dir)),
          cfg_(cfg),
          raw_writer_(raw_dir_, make_layout_("density_delta"),
                      {.overwrite = true}),
          coarse_writer_(coarse_dir_, make_layout_("density_delta_coarse"),
                         {.overwrite = true}) {}

    // Capture t0 reference. Must be called once before snapshot().
    void set_reference(inqkit::fields::RealField3D const& ref) {
        ref_ = ref;
        have_ref_ = true;
    }

    bool has_reference() const { return have_ref_; }

    // Compute and write dn for the current density field. Returns the
    // scalar L2 metric (= 0 if compute_l2 disabled).
    //
    // If no reference has been set, the first call captures `current` as
    // the reference and emits a zero delta (it is the t = t_0 frame by
    // construction). This is the recommended way to use the class from
    // inside a real-time callback, since pre-capturing the density
    // outside the propagator can return a stale snapshot — the
    // propagator's first iteration may rebuild the density before the
    // user-supplied callback is reached.
    double snapshot(inqkit::fields::RealField3D const& current,
                    double time_au, int step) {
        if (!have_ref_) {
            set_reference(current);
            // Emit a zero-delta frame for completeness so the GIF cadence
            // doesn't have a missing first frame.
            inqkit::fields::RealField3D zero = current;
            std::fill(zero.values.begin(), zero.values.end(), 0.0);
            if (cfg_.emit_raw_vti)    raw_writer_.write(zero, time_au, step);
            if (cfg_.emit_coarse_vti) {
                auto coarse = coarse_grain_(zero, cfg_.coarse_bin_bohr);
                coarse_writer_.write(coarse, time_au, step);
            }
            return 0.0;
        }
        if (current.values.size() != ref_.values.size()
            || current.nx != ref_.nx || current.ny != ref_.ny
            || current.nz != ref_.nz) {
            throw std::runtime_error(
                "inqkit::observables::DensityDelta::snapshot: "
                "grid mismatch between current density and reference.");
        }

        inqkit::fields::RealField3D delta = current;  // copies metadata
        for (std::size_t i = 0; i < delta.values.size(); ++i) {
            delta.values[i] = current.values[i] - ref_.values[i];
        }

        if (cfg_.emit_raw_vti) {
            raw_writer_.write(delta, time_au, step);
        }

        if (cfg_.emit_coarse_vti) {
            auto coarse = coarse_grain_(delta, cfg_.coarse_bin_bohr);
            coarse_writer_.write(coarse, time_au, step);
        }

        if (!cfg_.compute_l2) return 0.0;

        const double dV = delta.dx_bohr * delta.dy_bohr * delta.dz_bohr;
        long double s = 0.0L;
        for (auto v : delta.values) s += static_cast<long double>(v) * v;
        return static_cast<double>(s * dV);
    }

private:
    static inqkit::io::RealField3DLayout make_layout_(std::string name) {
        return {
            .field_name  = std::move(name),
            .include_meta = false,
            .emit_raw    = false,
            .emit_vti    = true,
            .vti_format  = inqkit::io::VTIWriteOptions::Format::binary,
        };
    }

    // Cubic-bin coarse-grain: each output voxel averages all input voxels
    // whose centre falls inside the bin. Output preserves the input cell
    // extent; output spacing = bin_bohr (rounded to integer multiples of
    // input spacing). For non-divisor bin sizes the residue at the high
    // end is folded into the last bin.
    static inqkit::fields::RealField3D
    coarse_grain_(inqkit::fields::RealField3D const& f, double bin_bohr) {
        const int sx = std::max(1, static_cast<int>(std::round(bin_bohr / f.dx_bohr)));
        const int sy = std::max(1, static_cast<int>(std::round(bin_bohr / f.dy_bohr)));
        const int sz = std::max(1, static_cast<int>(std::round(bin_bohr / f.dz_bohr)));

        const int Nx = std::max(1, f.nx / sx);
        const int Ny = std::max(1, f.ny / sy);
        const int Nz = std::max(1, f.nz / sz);

        inqkit::fields::RealField3D out;
        out.nx = Nx; out.ny = Ny; out.nz = Nz;
        out.dx_bohr = sx * f.dx_bohr;
        out.dy_bohr = sy * f.dy_bohr;
        out.dz_bohr = sz * f.dz_bohr;
        out.origin_x_bohr = f.origin_x_bohr;
        out.origin_y_bohr = f.origin_y_bohr;
        out.origin_z_bohr = f.origin_z_bohr;
        out.values.assign(static_cast<std::size_t>(Nx) * Ny * Nz, 0.0);

        for (int Ix = 0; Ix < Nx; ++Ix) {
            const int ix0 = Ix * sx;
            const int ix1 = (Ix == Nx - 1) ? f.nx : ix0 + sx;
            for (int Iy = 0; Iy < Ny; ++Iy) {
                const int iy0 = Iy * sy;
                const int iy1 = (Iy == Ny - 1) ? f.ny : iy0 + sy;
                for (int Iz = 0; Iz < Nz; ++Iz) {
                    const int iz0 = Iz * sz;
                    const int iz1 = (Iz == Nz - 1) ? f.nz : iz0 + sz;
                    long double sum = 0.0L;
                    long double n   = 0.0L;
                    for (int ix = ix0; ix < ix1; ++ix)
                    for (int iy = iy0; iy < iy1; ++iy)
                    for (int iz = iz0; iz < iz1; ++iz) {
                        auto flat = inqkit::detail::grid_layout::flatten_index(
                            ix, iy, iz, f.ny, f.nz);
                        sum += f.values[flat];
                        n   += 1.0L;
                    }
                    auto out_flat = inqkit::detail::grid_layout::flatten_index(
                        Ix, Iy, Iz, Ny, Nz);
                    out.values[out_flat] =
                        n > 0 ? static_cast<double>(sum / n) : 0.0;
                }
            }
        }
        return out;
    }

    std::string raw_dir_, coarse_dir_;
    DensityDeltaConfig cfg_;
    inqkit::io::RealField3DWriter raw_writer_;
    inqkit::io::RealField3DWriter coarse_writer_;
    inqkit::fields::RealField3D ref_;
    bool have_ref_ = false;
};

} // namespace inqkit::observables
