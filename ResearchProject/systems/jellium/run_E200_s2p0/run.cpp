// ============================================================================
// run_E200_s2p0 — N=38 closed-shell jellium, σ=2.0 Å, E=200 eV, +z (wide WP)
// ============================================================================

#include "../shared/configs/E200_s2p0.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::E200_s2p0>(
        "run_E200_s2p0",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L40_cubic_N38");
}
