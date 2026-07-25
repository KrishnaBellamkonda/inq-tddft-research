/*
 * The file that controls the primary observables being measured in the 
 * simulations. It streams scalar observables from an inqkit real-time
 * simulation to a CSV file, one row per time step.
 *
 * Columns are selected at construction time via ObservableSelection. Only
 * enabled columns appear in the file; the header row reflects exactly the
 * same selection so the output is self-describing.
 *
 * Available observables
 * ---------------------
 *   step                 Integer step counter.
 *   time_au              Simulation time in atomic units.
 *   energy_total         Total energy.
 *   energy_kinetic       Kinetic energy.
 *   energy_hartree       Hartree (Coulomb) energy.
 *   energy_xc            Exchange-correlation energy.
 *   energy_external      External-potential energy, E_ext = int n*v_ext.
 *   energy_nonlocal      Non-local pseudopotential energy.
 *   energy_ion           Ion-ion (Ewald / background) energy.
 *   energy_ion_kinetic   Classical ionic kinetic energy.
 *   energy_exact_exchange Exact-exchange energy (0 for pure LDA/GGA).
 *   energy_nvxc          int n*v_xc  (diagnostic; NOT part of the total).
 *   energy_eigenvalues   Sum of occupied eigenvalues (diagnostic; NOT in total).
 *   energy_proj_bg_ideal Classical projectile <-> jellium background Coulomb energy,
 *                        ideal formulation int n_proj*v_bg (DIAGNOSTIC; NOT in total).
 *                        A per-run constant (static ghost) — set via set_proj_bg().
 *   energy_proj_bg_impl  Same, as-implemented formulation -int n+*v_ion (r_cut-dependent).
 *   current_{x,y,z}      Current density vector components.
 *   dipole_{x,y,z}       Dipole moment vector components.
 *   cod_{x,y,z}_bohr     Wave-packet centre-of-density (Bohr).
 *   density_l2           Integrated squared density norm.
 *
 * Typical usage
 * -------------
 *   ObservableSelection sel;
 *   sel.energy_hartree = true;
 *   sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
 *
 *   ObservablesWriter writer("/output/observables.csv", sel);
 *   writer.write_header();               // call once before the time loop
 *
 *   for (auto const& ctx : steps) {
 *       // ... propagate ...
 *       writer.append(ctx);              // one row per step
 *   }
 *   writer.finish();                     // explicit flush + close
 *
 * Flush strategy
 * --------------
 * append() flushes to disk after every row. Without this, an aborted or
 * killed run leaves observables.csv at 0 bytes because the stream destructor
 * never executes. write_header() also flushes immediately so that even a
 * zero-step run produces a readable schema on disk.
 * 
 * Note: single-rank only, consistent with the existing inqkit writers.
 */
#pragma once

#include <inqkit/real_time/step_context.hpp>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace inqkit::io {


// TODO: The current vector and dipole vector must be tracked as a unit
// and not as x, y and z coordinates. This must be taken care of. 
// This is a structure that holds all the observables the simulation wants
// to track
struct ObservableSelection {
    bool step           = true;
    bool time_au        = true;
    bool energy_total   = true;
    bool energy_kinetic = true;
    bool energy_hartree = false;
    bool energy_xc      = false;
    // Full energy decomposition (default off — existing runs' schema unchanged).
    // The first six below + kinetic + hartree + xc sum to energy_total; nvxc and
    // eigenvalues are diagnostics outside the total.
    bool energy_external       = false;
    bool energy_nonlocal       = false;
    bool energy_ion            = false;
    bool energy_ion_kinetic    = false;
    bool energy_exact_exchange = false;
    bool energy_nvxc           = false;
    bool energy_eigenvalues    = false;
    // Projectile<->background Coulomb diagnostics (per-run constants, NOT in the
    // total). Off by default → existing schema unchanged. Values supplied by the
    // run via ObservablesWriter::set_proj_bg(); see projectile_background_energy.hpp.
    bool energy_proj_bg_ideal  = false;
    bool energy_proj_bg_impl   = false;
    bool current_x      = true;
    bool current_y      = true;
    bool current_z      = true;
    bool dipole_x       = false;
    bool dipole_y       = false;
    bool dipole_z       = false;
    // WP centre-of-density and integrated dn^2; populated by the run-template
    // before each call to append() (left zero if disabled).
    bool cod_x          = false;
    bool cod_y          = false;
    bool cod_z          = false;
    bool density_l2     = false;
};

// Streams selected scalar and vector observables to a CSV file.
// One row per call to append(); call write_header() before the first append().
class ObservablesWriter {
public:
    ObservablesWriter(std::string const& path,
                      ObservableSelection sel = {},
                      char separator = ',')
        : sel_(sel), sep_(separator)
    {
        namespace fs = std::filesystem;
        if (auto parent = fs::path(path).parent_path(); !parent.empty())
            fs::create_directories(parent);

        file_.open(path);
        if (!file_) throw std::runtime_error("ObservablesWriter: cannot open '" + path + "'");
    }

    // Supply the per-run projectile<->background Coulomb energies (Hartree).
    // Call once before the time loop; the values are emitted every row.
    void set_proj_bg(double ideal, double impl) {
        proj_bg_ideal_ = ideal;
        proj_bg_impl_  = impl;
    }

    void write_header() {
        bool first = true;
        auto col = [&](char const* name) {
            if (!first) file_ << sep_;
            file_ << name;
            first = false;
        };
        // Flush after the header so that even an early-abort run produces
        // a readable schema on disk.
        struct _flush_at_end {
            std::ofstream& f; ~_flush_at_end() { f.flush(); }
        } _fae{file_};
        if (sel_.step)           col("step");
        if (sel_.time_au)        col("time_au");
        if (sel_.energy_total)   col("energy_total");
        if (sel_.energy_kinetic) col("energy_kinetic");
        if (sel_.energy_hartree) col("energy_hartree");
        if (sel_.energy_xc)      col("energy_xc");
        if (sel_.energy_external)       col("energy_external");
        if (sel_.energy_nonlocal)       col("energy_nonlocal");
        if (sel_.energy_ion)            col("energy_ion");
        if (sel_.energy_ion_kinetic)    col("energy_ion_kinetic");
        if (sel_.energy_exact_exchange) col("energy_exact_exchange");
        if (sel_.energy_nvxc)           col("energy_nvxc");
        if (sel_.energy_eigenvalues)    col("energy_eigenvalues");
        if (sel_.energy_proj_bg_ideal)  col("energy_proj_bg_ideal");
        if (sel_.energy_proj_bg_impl)   col("energy_proj_bg_impl");
        if (sel_.current_x)      col("current_x");
        if (sel_.current_y)      col("current_y");
        if (sel_.current_z)      col("current_z");
        if (sel_.dipole_x)       col("dipole_x");
        if (sel_.dipole_y)       col("dipole_y");
        if (sel_.dipole_z)       col("dipole_z");
        if (sel_.cod_x)          col("cod_x_bohr");
        if (sel_.cod_y)          col("cod_y_bohr");
        if (sel_.cod_z)          col("cod_z_bohr");
        if (sel_.density_l2)     col("density_l2");
        file_ << '\n';
    }

    void append(inqkit::StepContext const& ctx) {
        file_ << std::setprecision(15);
        bool first = true;
        auto val = [&](auto v) {
            if (!first) file_ << sep_;
            file_ << v;
            first = false;
        };
        if (sel_.step)           val(ctx.step);
        if (sel_.time_au)        val(ctx.time_au);
        if (sel_.energy_total)   val(ctx.energy_total);
        if (sel_.energy_kinetic) val(ctx.energy_kinetic);
        if (sel_.energy_hartree) val(ctx.energy_hartree);
        if (sel_.energy_xc)      val(ctx.energy_xc);
        if (sel_.energy_external)       val(ctx.energy_external);
        if (sel_.energy_nonlocal)       val(ctx.energy_nonlocal);
        if (sel_.energy_ion)            val(ctx.energy_ion);
        if (sel_.energy_ion_kinetic)    val(ctx.energy_ion_kinetic);
        if (sel_.energy_exact_exchange) val(ctx.energy_exact_exchange);
        if (sel_.energy_nvxc)           val(ctx.energy_nvxc);
        if (sel_.energy_eigenvalues)    val(ctx.energy_eigenvalues);
        if (sel_.energy_proj_bg_ideal)  val(proj_bg_ideal_);
        if (sel_.energy_proj_bg_impl)   val(proj_bg_impl_);
        if (sel_.current_x)      val(ctx.current[0]);
        if (sel_.current_y)      val(ctx.current[1]);
        if (sel_.current_z)      val(ctx.current[2]);
        if (sel_.dipole_x)       val(ctx.dipole[0]);
        if (sel_.dipole_y)       val(ctx.dipole[1]);
        if (sel_.dipole_z)       val(ctx.dipole[2]);
        if (sel_.cod_x)          val(ctx.wp_center[0]);
        if (sel_.cod_y)          val(ctx.wp_center[1]);
        if (sel_.cod_z)          val(ctx.wp_center[2]);
        if (sel_.density_l2)     val(ctx.density_l2);
        file_ << '\n';
        file_.flush();   // Important: without this an aborted run leaves
                         // observables.csv at 0 bytes, since the destructor
                         // never executes when the process is killed.
    }

    void finish() {
        file_.flush();
        file_.close();
    }

    ~ObservablesWriter() {
        if (file_.is_open()) finish();
    }

private:
    std::ofstream       file_;
    ObservableSelection sel_;
    char                sep_;
    double proj_bg_ideal_ = 0.0;   // per-run projectile<->background Coulomb (set_proj_bg)
    double proj_bg_impl_  = 0.0;
};

} // namespace inqkit::io
