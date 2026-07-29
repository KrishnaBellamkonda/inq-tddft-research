// ============================================================================
// run_base — Tsubonoya 2014 base configuration with corrected geometry.
//
// b = 12 Bohr, sigma = 1.0 Bohr, E = 200 eV, cell 35 x 35 x 60 Bohr,
// 600 steps at dt = 0.020 a.u., 20 LEED screens.
//
// GS used (compiled and saved by save_gs/gs_35x35x60_cut40/run.cpp):
/*
    auto cell = systems::cell::orthorhombic(
        Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, Cfg::LZ_BOHR * 1.0_b
    ).finite();
    auto ions = systems::ions::parse(SHARED_GEOMETRY_XYZ, cell);
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .cutoff(Cfg::CUTOFF_HA * 1.0_Ha)        // 40 Ha
            .extra_states(Cfg::EXTRA_STATES));      // 8
    ground_state::initial_guess(ions, electrons);
    ground_state::calculate(ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(1e-6_Ha)
            .max_steps(1000)
            .broyden_mixing()
            .mixing_ndim(8)
            .mixing(0.1));
    electrons.save(GS_CHECKPOINT_DIR);
*/
// ============================================================================

#include "../shared/configs/tsubonoya_2014_base.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return coronene::run_template::run_propagation<coronene::config::Base>(
        "run_base",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/shared/geometry/coronene.xyz",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/checkpoints/gs_35x35x60_cut40");
}
