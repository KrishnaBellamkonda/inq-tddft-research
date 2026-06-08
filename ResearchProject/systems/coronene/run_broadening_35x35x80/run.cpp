// ============================================================================
// run_broadening_35x35x80 — long pre-collision flight (launch z = +30 Bohr) for
// Gaussian wave-packet spreading, centroid kept inside the box (no boundary
// crossing). WRITE_EVERY = 1 for maximum temporal resolution.
//
// GS used: checkpoints/gs_35x35x80_cut40/ (cell 35x35x80, cutoff 40 Ha, LDA).
// Primary deliverable: pre-collision sigma_z(t) Gaussian-broadening plot.
// ============================================================================
#include "../shared/configs/broadening_35x35x80.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return coronene::run_template::run_propagation<coronene::config::broadening_35x35x80>(
        "run_broadening_35x35x80",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/shared/geometry/coronene.xyz",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/checkpoints/gs_35x35x80_cut40");
}
