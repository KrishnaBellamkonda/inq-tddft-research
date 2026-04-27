// ============================================================================
// run_35x35x40 — smaller box (Lz = 40 Bohr = 2/3 x base). b at base.
//
// GS used: save_gs/gs_35x35x40_cut40/ (cell 35x35x40, cutoff 40 Ha, LDA).
// See save_gs/gs_35x35x40_cut40/run.cpp for the GS construction code.
// ============================================================================
#include "../shared/configs/cell_35x35x40.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return coronene::run_template::run_propagation<coronene::config::cell_35x35x40>(
        "run_35x35x40",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/shared/geometry/coronene.xyz",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/checkpoints/gs_35x35x40_cut40");
}
