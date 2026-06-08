// ============================================================================
// run_E200_s0p265 — N=38 closed-shell jellium, σ=0.265 Å, E=200 eV, +z (narrow WP)
// ============================================================================

#include "../shared/configs/E200_s0p265.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::E200_s0p265>(
        "run_E200_s0p265",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L40_cubic_N38");
}
