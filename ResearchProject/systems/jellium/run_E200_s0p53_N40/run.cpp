// ============================================================================
// run_E200_s0p53_N40 — N=40 open-shell jellium, σ=0.53 Å, E=200 eV, +z
// Uses checkpoint gs_L40_cubic_N40 (its own SCF run).
// ============================================================================

#include "../shared/configs/E200_s0p53_N40.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::E200_s0p53_N40>(
        "run_E200_s0p53_N40",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L40_cubic_N40");
}
