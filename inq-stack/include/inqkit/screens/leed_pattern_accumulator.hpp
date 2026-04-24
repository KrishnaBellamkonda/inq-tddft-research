/*
 * This class uses the PlaneScreen API to collect the LEED pattern for a given setup.
 * .accumulate() function is key and will be written in the .cpp files that we
 * write using this library. 
 *
 *
 * */

#pragma once

// inqkit::screens::LeedPatternAccumulator
//
// Time-integrates the 2D electron density on a PlaneScreen:
//   pattern(x,y) += ρ(x, y, z_screen, t) × dt
//
// Call accumulate() at every RT step (or selectively). Call save() once after
// propagation to write the accumulated pattern to disk.
//
// File format written by save():
//   # label=LABEL z=Z total_time=T n_accum=N
//   # nx=NX ny=NY dx=DX dy=DY origin_x=OX origin_y=OY
//   v00 v01 ... v0(NX-1)
//   v10 ...
//   ...

#include "plane_screen.hpp"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace inqkit {
namespace screens {

class LeedPatternAccumulator {
    PlaneScreen                      screen_;
    std::vector<std::vector<double>> pattern_;   // [Ny][Nx]
    int    nx_           = 0;
    int    ny_           = 0;
    double dx_           = 0.0;
    double dy_           = 0.0;
    bool   initialised_  = false;
    double total_time_au_= 0.0;
    int    n_accum_      = 0;

    void initialise(inq::systems::electrons const & electrons) {
        nx_ = screen_.nx(electrons);
        ny_ = screen_.ny(electrons);
        dx_ = screen_.dx(electrons);
        dy_ = screen_.dy(electrons);
        pattern_.assign(ny_, std::vector<double>(nx_, 0.0));
        initialised_ = true;
    }

public:
    LeedPatternAccumulator() = default;

    explicit LeedPatternAccumulator(PlaneScreen screen)
        : screen_(std::move(screen)) {}

    PlaneScreen const& screen() const { return screen_; }

    // Call at each RT time step.  dt_au is the propagation time step.
    // Occupation weights are read from electrons.occupations() — after
    // WavePacket::inject_into_last_extra_state() the WP has occ=1.0 and is
    // included automatically. No wp_state_global argument needed.
    void accumulate(inq::systems::electrons const & electrons, double dt_au) {
        if (!initialised_) initialise(electrons);

        auto slice = screen_.extract(electrons);
        assert(static_cast<int>(slice.size())    == ny_);
        assert(static_cast<int>(slice[0].size()) == nx_);

        for (int iy = 0; iy < ny_; iy++)
            for (int ix = 0; ix < nx_; ix++)
                pattern_[iy][ix] += slice[iy][ix] * dt_au;

        total_time_au_ += dt_au;
        n_accum_++;
    }

    // Write accumulated pattern.  Creates parent directories automatically.
    void save(std::string const & filename) const {
        if (!initialised_) return;  // nothing accumulated
        std::filesystem::create_directories(
            std::filesystem::path(filename).parent_path());

        std::ofstream f(filename);
        if (!f) throw std::runtime_error("LeedPatternAccumulator: cannot open " + filename);

        f << std::fixed << std::setprecision(6);
        f << "# label=" << screen_.label()
          << " z="          << screen_.z_bohr()
          << " total_time=" << total_time_au_
          << " n_accum="    << n_accum_ << "\n";
        f << "# nx=" << nx_
          << " ny=" << ny_
          << " dx=" << dx_
          << " dy=" << dy_
          << " origin_x=0.000000"
          << " origin_y=0.000000\n";

        for (int iy = 0; iy < ny_; iy++) {
            for (int ix = 0; ix < nx_; ix++) {
                f << std::scientific << std::setprecision(6) << pattern_[iy][ix];
                if (ix + 1 < nx_) f << " ";
            }
            f << "\n";
        }
    }

    std::vector<std::vector<double>> const& pattern()       const { return pattern_;       }
    double                                  total_time_au() const { return total_time_au_; }
    int                                     n_accum()       const { return n_accum_;       }
    bool                                    initialised()   const { return initialised_;   }
};

} // namespace screens
} // namespace inqkit
