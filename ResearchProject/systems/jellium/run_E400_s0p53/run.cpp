// ============================================================================
// run_E400_s0p53 — N=38 closed-shell jellium, σ=0.53 Å, E=400 eV, +z
// High-energy WP variant of run_base.
// ============================================================================

#include "../shared/configs/E400_s0p53.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::E400_s0p53>(
        "run_E400_s0p53",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L40_cubic_N38");
}
