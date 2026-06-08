// ============================================================================
// run_E200_s0p53_tilt45 — N=38 closed-shell jellium, σ=0.53 Å, E=200 eV,
// 45° tilt in xz-plane (k_x = k_z = k₀/√2).
// ============================================================================

#include "../shared/configs/E200_s0p53_tilt45.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::E200_s0p53_tilt45>(
        "run_E200_s0p53_tilt45",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L40_cubic_N38");
}
