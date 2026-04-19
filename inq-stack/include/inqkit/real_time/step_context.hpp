#pragma once

#include <math/vector3.hpp>
#include <systems/ions.hpp>
#include <systems/electrons.hpp>

namespace inqkit {

struct StepContext {
    int    step    = 0;
    double time_au = 0.0;
    inq::systems::ions     const* ions      = nullptr;
    inq::systems::electrons const* electrons = nullptr;

    double energy_total   = 0.0;
    double energy_kinetic = 0.0;
    double energy_hartree = 0.0;
    double energy_xc      = 0.0;
    inq::vector3<double> current = {0.0, 0.0, 0.0};
    inq::vector3<double> dipole  = {0.0, 0.0, 0.0};
};

} // namespace inqkit
