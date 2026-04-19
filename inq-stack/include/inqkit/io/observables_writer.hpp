#pragma once

#include <inqkit/real_time/step_context.hpp>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace inqkit::io {

struct ObservableSelection {
    bool step           = true;
    bool time_au        = true;
    bool energy_total   = true;
    bool energy_kinetic = true;
    bool energy_hartree = false;
    bool energy_xc      = false;
    bool current_x      = true;
    bool current_y      = true;
    bool current_z      = true;
    bool dipole_x       = false;
    bool dipole_y       = false;
    bool dipole_z       = false;
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

    void write_header() {
        bool first = true;
        auto col = [&](char const* name) {
            if (!first) file_ << sep_;
            file_ << name;
            first = false;
        };
        if (sel_.step)           col("step");
        if (sel_.time_au)        col("time_au");
        if (sel_.energy_total)   col("energy_total");
        if (sel_.energy_kinetic) col("energy_kinetic");
        if (sel_.energy_hartree) col("energy_hartree");
        if (sel_.energy_xc)      col("energy_xc");
        if (sel_.current_x)      col("current_x");
        if (sel_.current_y)      col("current_y");
        if (sel_.current_z)      col("current_z");
        if (sel_.dipole_x)       col("dipole_x");
        if (sel_.dipole_y)       col("dipole_y");
        if (sel_.dipole_z)       col("dipole_z");
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
        if (sel_.current_x)      val(ctx.current[0]);
        if (sel_.current_y)      val(ctx.current[1]);
        if (sel_.current_z)      val(ctx.current[2]);
        if (sel_.dipole_x)       val(ctx.dipole[0]);
        if (sel_.dipole_y)       val(ctx.dipole[1]);
        if (sel_.dipole_z)       val(ctx.dipole[2]);
        file_ << '\n';
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
};

} // namespace inqkit::io
