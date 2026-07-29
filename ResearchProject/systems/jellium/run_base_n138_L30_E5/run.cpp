// ============================================================================
// run_base_n138_L30_E5 - canonical low-energy WP-jellium scattering run
// (project base, set 2026-05-04). N=138 closed shell, L=30 cubic,
// dx=0.85 bohr, E_WP = 5 eV, t_final = 19.80 a.u. (990 steps at dt=0.020).
//
// See `.claude/rules/jellium-base-run-spec.md` for the canonical spec.
// Cfg = jellium::config::Base_N138_L30_E5.
// GS loaded from save_gs/gs_L30_cubic_N138_dx0p85.
// ============================================================================

#include "../shared/configs/base_n138_L30_E5.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::Base_N138_L30_E5>(
        "run_base_n138_L30_E5",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L30_cubic_N138_dx0p85");
}
