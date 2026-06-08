// ============================================================================
// run_base_n138_L50_E1p5 - L=50 jellium with low-energy, wide WP.
// Aim: observe inelastic e-h coupling without the periodic-box revival
// contamination of the L=30 base run; confirm WP slowdown is real.
//
// N=138 closed shell |G|^2<=6, dx=1.0 bohr, sigma=5.0 bohr, E_WP=1.5 eV,
// 1500 steps at dt=0.020 -> t_final = 30.0 a.u. (= 0.726 fs).
//
// Cfg = jellium::config::Base_N138_L50_E1p5.
// GS loaded from save_gs/gs_L50_cubic_N138_dx1p0.
// ============================================================================

#include "../shared/configs/base_n138_L50_E1p5.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::Base_N138_L50_E1p5>(
        "run_base_n138_L50_E1p5",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L50_cubic_N138_dx1p0");
}
