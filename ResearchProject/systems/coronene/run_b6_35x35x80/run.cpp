// ============================================================================
// run_b6_35x35x80 — large box (Lz = 80 Bohr) + b = 6 Bohr (= 0.5 x base).
//
// GS used: save_gs/gs_35x35x80_cut40/ (cell 35x35x80, cutoff 40 Ha, LDA).
// See save_gs/gs_35x35x80_cut40/run.cpp for the GS construction code.
// ============================================================================
#include "../shared/configs/b6_35x35x80.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return coronene::run_template::run_propagation<coronene::config::b6_35x35x80>(
        "run_b6_35x35x80",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/shared/geometry/coronene.xyz",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/checkpoints/gs_35x35x80_cut40");
}
