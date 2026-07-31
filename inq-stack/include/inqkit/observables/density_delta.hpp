/*
 * This file tracks the local density fluctuation field at each time step of a TDDFT
 * run, defined relative to a user-supplied reference density n(r, t₀):
 *
 *   δn(r, t) = n(r, t) − n(r, t₀)
 *
 * Three views are computed and optionally written per call to snapshot():
 *
 *   Raw δn          Full-resolution 3D field written as a VTI series.
 *                   Intended for visualisation with a diverging colourmap.
 *
 *   Coarse δn       δn averaged into cubic bins of coarse_bin_bohr, written
 *                   as a separate VTI series. Suppresses Friedel-oscillation
 *                   noise (which are unwanted oscillations that might exist
 *                   due to impurities in the materials)
 *                   and exposes the macroscopic density redistribution.
 *                   Bin size is rounded to the nearest integer multiple of
 *                   the input spacing; residual voxels at the high edge are
 *                   folded into the last bin.
 *
 *   L2 metric       σ²(t) = ∫ |δn|² dV, returned as a scalar from
 *                   snapshot() and suitable for logging via ObservablesWriter.
 *
 * Reference capture
 * -----------------
 * set_reference() must be called before the first snapshot(). If it is not,
 * the first snapshot() call captures the supplied field as the reference
 * automatically and emits a zero-delta frame so the output series has no
 * missing first frame. In practice, the SCF captured density is different from
 * the density at the first timestep. We are trying to capture the density at
 * the first timestep. This lazy-capture mode is recommended when calling
 * from inside a real-time callback, since the propagator may rebuild the
 * density before the callback is reached, making any pre-captured snapshot
 * stale.
 *
 * Usage
 * -----
 *   DensityDeltaConfig cfg;
 *   cfg.coarse_bin_bohr = 3.0;
 *
 *   DensityDelta dd("/output/dn_raw", "/output/dn_coarse", cfg);
 *   dd.set_reference(n_t0);
 *
 *   // inside the time loop:
 *   double l2 = dd.snapshot(n_current, time_au, step);
 *
 * All computation runs on the host. Cost per step on a 60³ grid is
 * sub-millisecond, plus VTI write time.
 */

// TODO: Verify that the first snapshot taken is that of t=0 iteration during the
// real time propagation phase. Or should I take it from t=dt (1st timestep) after
// a step of real time propagation has happened. Would this remove the artifact of
// the deep hole behind? Run a simulation where both these situations are run and
// visualised. I can then make a decision as to which approach is better.
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

/*
 * Configurations for computing the different densities from the simulation.
 */
struct DensityDeltaConfig {
    bool   emit_raw_vti     = true;
    bool   emit_coarse_vti  = true;
    /* σ²(t) can be thought of as the total statistical fluctuation in the
     * delta density — analogous to variance, it measures the mean-square
     * deviation of n(r,t) from the reference n(r,t₀), integrated over the
     * cell volume:
     *
     *   σ²(t) = ∫ |δn(r,t)|² dV  =  ∫ |n(r,t) − n(r,t₀)|² dV
     *
     * This is distinct from the total-density L2 (∫ |n(r,t)|² dV), which
     * reflects the absolute charge distribution rather than its change.
     *
     * Because the integrand is squared, spatially localised (spiky) changes
     * contribute disproportionately: two systems with the same net displaced
     * charge can have very different σ²(t) values. A higher σ²(t) indicates
     * the density redistribution is more concentrated in space; a lower σ²(t)
     * means the same charge has moved more smoothly across the cell.
     * Comparing σ²(t) across systems or time steps is therefore a measure of
     * how localised the response is, not how large the total charge transfer is.
     */
    bool   compute_l2       = true;
    double coarse_bin_bohr  = 3.0;

    /* VTI write cadence: emit a delta FIELD only when step % emit_every == 0.
     *
     * The scalar L2 above is cheap and is normally wanted EVERY step, but the
     * field is a full grid (18 MB at 35x35x85 / dx=0.4) and is not. Callers that
     * call snapshot() every step to keep the L2 series dense would otherwise also
     * write a field every step: at 3623 steps that is 66 GB per run.
     *
     * That is not hypothetical. It filled /rds (1.0 TB, 100 %) on 2026-07-31 and
     * killed three of the four sigma=3 wavepacket runs mid-flight, their vacuum
     * controls, and the notebook job -- density_delta held 3624 frames where
     * density_total and density_wp held 302 at the same SAVE_EVERY=12 cadence.
     *
     * DEFAULT 1 = write on every call, i.e. the historical behaviour, so none of
     * the ~130 existing run.cpp callers change. Set it to the run's SAVE_EVERY to
     * align delta frames with density_total / density_wp (which is the pairing the
     * GIF builders assume anyway). <= 0 is treated as 1.
     */
    int    emit_every       = 1;
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
    // If no reference has been set, the first call captures `current` as
    // the reference and emits a zero delta (it is the t = t_0 frame by
    // construction). This is the recommended way to use the class from
    // inside a real-time callback, since pre-capturing the density
    // outside the propagator can return a stale snapshot — the
    // propagator's first iteration may rebuild the density before the
    // user-supplied callback is reached.
    double snapshot(inqkit::fields::RealField3D const& current,
                    double time_au, int step) {

        // If the reference density is not explicitly set,
        // treat the first supplied density as the reference (lazy capture).
        // TODO: check — does this mean density at t+1 uses t as base,
        // or is t=0 the base for all timesteps?
        // TODO: It would also be interesting to view the step-by-step changes,
        // i.e. changes induced only within each individual timestep.
        // Field-emission cadence (see DensityDeltaConfig::emit_every). The scalar
        // L2 below is still computed on EVERY call; only the grid write is gated.
        const int  every = (cfg_.emit_every > 0) ? cfg_.emit_every : 1;
        const bool emit  = (step % every == 0);

        if (!have_ref_) {
            set_reference(current);
            // Emit a zero-delta frame so the output series has no missing
            // first frame. (step is 0 here in normal use, so `emit` is true;
            // the guard matters only for a resumed run whose first call lands
            // on an off-cadence step.)
            inqkit::fields::RealField3D zero = current;
            std::fill(zero.values.begin(), zero.values.end(), 0.0);
            if (cfg_.emit_raw_vti && emit)    raw_writer_.write(zero, time_au, step);
            if (cfg_.emit_coarse_vti && emit) {
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

        // delta copies the grid metadata from current (dimensions, spacing,
        // origin). The 3D field is stored as a flat contiguous 1D array
        // (row-major: index = ix*ny*nz + iy*nz + iz), which maps directly
        // to the VTI binary data block with no reshaping needed.
        inqkit::fields::RealField3D delta = current;
        for (std::size_t i = 0; i < delta.values.size(); ++i) {
            delta.values[i] = current.values[i] - ref_.values[i];
        }

        if (cfg_.emit_raw_vti && emit) {
            raw_writer_.write(delta, time_au, step);
        }

        if (cfg_.emit_coarse_vti && emit) {
            auto coarse = coarse_grain_(delta, cfg_.coarse_bin_bohr);
            coarse_writer_.write(coarse, time_au, step);
        }

        if (!cfg_.compute_l2) return 0.0;

        // Discrete approximation of σ²(t) = ∫ |δn|² dV,
        // accumulated in long double to reduce floating-point error
        // when summing many small squared values.
        const double dV = delta.dx_bohr * delta.dy_bohr * delta.dz_bohr;
        long double s = 0.0L;
        for (auto v : delta.values) s += static_cast<long double>(v) * v;
        return static_cast<double>(s * dV);
    }

// TODO: Is using "function_" a convention consistently followed in the wrapper
// library? i.e. is this found throughout the inq-stack C++ and Python codebase?

private:
    static inqkit::io::RealField3DLayout make_layout_(std::string name) {
        return {
            .field_name   = std::move(name),
            .include_meta = false,
            .emit_raw     = false,
            .emit_vti     = true,
            .vti_format   = inqkit::io::VTIWriteOptions::Format::binary,
        };
    }

    // voxel: the 3D equivalent of a 2D pixel.
    //
    // Cubic-bin coarse-grain: each output voxel averages all input voxels
    // whose centre falls inside the bin. Output preserves the input cell
    // extent; output spacing = bin_bohr (rounded to the nearest integer
    // multiple of the input spacing).
    //
    // Boundary handling: if the grid does not divide evenly, the remainder
    // voxels at the high edge are folded into the last bin, which therefore
    // averages more input points than interior bins. The output metadata
    // reports a uniform spacing (sx * dx_bohr) for all bins including the
    // last, so the last bin's physical extent is slightly underreported.
    // This is acceptable for visualisation but should be noted if the coarse
    // field is used for quantitative spatial analysis (e.g. dipole moments).
    static inqkit::fields::RealField3D
    coarse_grain_(inqkit::fields::RealField3D const& f, double bin_bohr) {
        // Number of fine grid points per coarse bin along each axis,
        // rounded to the nearest integer.
        const int sx = std::max(1, static_cast<int>(std::round(bin_bohr / f.dx_bohr)));
        const int sy = std::max(1, static_cast<int>(std::round(bin_bohr / f.dy_bohr)));
        const int sz = std::max(1, static_cast<int>(std::round(bin_bohr / f.dz_bohr)));

        // Number of coarse bins along each axis (floor division).
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
        // Flat 1D storage (row-major), same layout as the input field.
        out.values.assign(static_cast<std::size_t>(Nx) * Ny * Nz, 0.0);

        for (int Ix = 0; Ix < Nx; ++Ix) {
            const int ix0 = Ix * sx;
            // Last bin absorbs any remainder fine points beyond the regular stride.
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
                        // Convert 3D index to flat 1D: ix*ny*nz + iy*nz + iz
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