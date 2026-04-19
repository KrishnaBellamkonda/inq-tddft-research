#pragma once

#include <systems/ions.hpp>
#include <systems/electrons.hpp>

namespace inqkit {

struct StepContext {
    int    step    = 0;
    double time_au = 0.0;
    inq::systems::ions     const* ions      = nullptr;
    inq::systems::electrons const* electrons = nullptr;
};

} // namespace inqkit
