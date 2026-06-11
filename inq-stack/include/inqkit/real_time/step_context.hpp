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

#include <inqkit/detail/vec3.hpp>

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
    // current/dipole as a Vec3 unit (consistency with center_of_density which
    // already returns inqkit::detail::Vec3). CSV columns current_{x,y,z} /
    // dipole_{x,y,z} are unchanged — observables_writer reads ctx.current[0..2]
    // via Vec3::operator[]. Conversion from INQ's vector3 happens at the one
    // callback site (real_time_session.hpp).
    inqkit::detail::Vec3 current = {0.0, 0.0, 0.0};
    inqkit::detail::Vec3 dipole  = {0.0, 0.0, 0.0};

    // Slots populated by jellium-side per-step callbacks (centre of WP
    // density, integrated dn^2). Left zero when not computed.
    inq::vector3<double> wp_center = {0.0, 0.0, 0.0};
    double density_l2 = 0.0;
};

} // namespace inqkit
