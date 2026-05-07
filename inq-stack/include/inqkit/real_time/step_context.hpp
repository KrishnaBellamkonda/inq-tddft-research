/* The idea of this module is simple. As the real_time::propagate() runs, 
 * we are allowed to give this function a callback function which runs
 * every step. We define such callbacks using the Pipeline and Tasks APIs. 
 * However, we need a class to be able to hold the current state of the 
 * real_time propagation. This StepContext class does this by holding 
 * 1. step number
 * 2. time (in atomic units)
 * 3. current ions (handed over from data in the obs_callback)
 * 4. current electrons (handed over from data)
 * 5. Observables we are keeping a track of in the simulation
 *
 * */

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

    // Slots populated by jellium-side per-step callbacks (centre of WP
    // density, integrated dn^2). Left zero when not computed.
    inq::vector3<double> wp_center = {0.0, 0.0, 0.0};
    double density_l2 = 0.0;
};

} // namespace inqkit
