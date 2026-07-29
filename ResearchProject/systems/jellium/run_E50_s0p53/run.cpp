// ============================================================================
// run_E50_s0p53 — N=38 closed-shell jellium, σ=0.53 Å, E=50 eV, +z
// Low-energy WP variant of run_base (½ k₀ ⇒ 2× steps for same z-range).
// ============================================================================

#include "../shared/configs/E50_s0p53.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::E50_s0p53>(
        "run_E50_s0p53",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L40_cubic_N38");
}
