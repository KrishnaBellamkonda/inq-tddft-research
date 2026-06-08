// ============================================================================
// run_base — N=128 closed-shell jellium at preserved r_s ≈ 7.38 a₀,
// σ=0.53 Å, E=200 eV, +z. Cell 60^3 bohr (cubic, periodic),
// 320 steps at dt=0.020 a.u. (single-pass), 20 LEED screens.
// GS loaded from save_gs/gs_L60_cubic_N128.
// ============================================================================

#include "../shared/configs/base.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::Base>(
        "run_base",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L60_cubic_N128");
}
