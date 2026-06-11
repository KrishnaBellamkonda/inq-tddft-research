/* RealTimeSession is used to handle tasks during real_time::propagate runs. 
 * To an initialised RealTimeSession rt(), rt.add() adds tasks to a stack. 
 * Each task is a function that takes in the StepContext object initialised
 * on data that returned in each timestep of real time session. When the timestep
 * that is include in the WRITE_EVERY steps is included, the tasks are run
 * in sequence. 
 *
 * */
#pragma once

#include <functional> // helps sending functions as arguments to functions
#include <vector>

#include <inqkit/real_time/step_context.hpp>
#include <systems/electrons.hpp>
#include <systems/ions.hpp>

namespace inqkit {

// Minimal real-time session: owns a list of tasks and dispatches every
// write_every steps. Each task is a callable void(StepContext const&). Usage in
// a propagate callback:
//
//   RealTimeSession rt(ions, electrons, /*write_every=*/100);
//   rt.add([&](StepContext const& ctx) { ... });
//   real_time::propagate(ions, electrons,
//       [&](auto const& data) { rt.step(data); }, ...);
class RealTimeSession {
public:
  // Initialisation function takes in ions (ground state), electrons (ground
  // state) and write_every argument (which defines per how many timesteps
  // this session must run)
  RealTimeSession(inq::systems::ions &ions, inq::systems::electrons &electrons,
                  int write_every = 1)
      : ions_(ions), electrons_(electrons), write_every_(write_every) {}

  // Each task is a function that takes in the step context
  // and performs an action using the step context
  void add(std::function<void(StepContext const &)> task) {
    tasks_.push_back(std::move(task));
  }

  // Called inside the INQ propagate lambda with the viewables object.
  // Dispatches tasks every write_every steps (including step 0).
  // 
  template <typename RTData> void step(RTData const &data) {
    if (data.iter() % write_every_ != 0)
      return;

    // Building the step context class using file the data passed
    // to the optional callback funtion in the step iterations
    StepContext ctx;
    ctx.step = data.iter();
    ctx.time_au = data.time();
    ctx.ions = &ions_;
    ctx.electrons = &electrons_;
    ctx.energy_total = data.energy().total();
    ctx.energy_kinetic = data.energy().kinetic();
    ctx.energy_hartree = data.energy().hartree();
    ctx.energy_xc = data.energy().xc();
    // current/dipole computed on demand; guard in case observables aren't
    // active
    // Convert INQ's vector3 -> inqkit::detail::Vec3 (StepContext's type) here,
    // the single callback site, so the pure Vec3 header stays INQ-free.
    try {
      auto c = data.current();
      ctx.current = {c[0], c[1], c[2]};
    } catch (...) {
    }
    try {
      auto d = data.dipole();
      ctx.dipole = {d[0], d[1], d[2]};
    } catch (...) {
    }
    // All of the tasks are executed at the correct timestep
    for (auto &task : tasks_) {
      task(ctx);
    }
  }

private:
  inq::systems::ions &ions_;
  inq::systems::electrons &electrons_;
  int write_every_;
  std::vector<std::function<void(StepContext const &)>> tasks_;
};

} // namespace inqkit
