// ============================================================================
// run_E800 — high-energy WP probe (E = 800 eV). Everything else at base.
//
// GS used: save_gs/gs_35x35x60_cut40/ (cell 35x35x60, cutoff 40 Ha, LDA).
// See run_base/run.cpp for the commented GS construction code.
// ============================================================================
#include "../shared/configs/E800.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return coronene::run_template::run_propagation<coronene::config::E800>(
        "run_E800",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/shared/geometry/coronene.xyz",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/checkpoints/gs_35x35x60_cut40");
}
