#pragma once

#include <functional>
#include <vector>

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <inqkit/real_time/step_context.hpp>

namespace inqkit {

// Minimal real-time session: owns a list of tasks and dispatches every write_every steps.
// Each task is a callable void(StepContext const&).
// Usage in a propagate callback:
//
//   RealTimeSession rt(ions, electrons, /*write_every=*/100);
//   rt.add([&](StepContext const& ctx) { ... });
//   real_time::propagate(ions, electrons,
//       [&](auto const& data) { rt.step(data); }, ...);
class RealTimeSession {
public:
    RealTimeSession(inq::systems::ions& ions,
                    inq::systems::electrons& electrons,
                    int write_every = 1)
        : ions_(ions), electrons_(electrons), write_every_(write_every) {}

    void add(std::function<void(StepContext const&)> task) {
        tasks_.push_back(std::move(task));
    }

    // Called inside the INQ propagate lambda with the viewables object.
    // Dispatches tasks every write_every steps (including step 0).
    template <typename RTData>
    void step(RTData const& data) {
        if (data.iter() % write_every_ != 0) return;

        StepContext ctx;
        ctx.step      = data.iter();
        ctx.time_au   = data.time();
        ctx.ions      = &ions_;
        ctx.electrons = &electrons_;

        for (auto& task : tasks_) {
            task(ctx);
        }
    }

private:
    inq::systems::ions&      ions_;
    inq::systems::electrons& electrons_;
    int write_every_;
    std::vector<std::function<void(StepContext const&)>> tasks_;
};

} // namespace inqkit
