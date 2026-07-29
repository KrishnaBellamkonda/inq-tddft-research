// ============================================================================
// run_plasmon_n162_L50_E3p4_varyv — Run D, variable-velocity discriminator
// for the Run B plasmon detection. Same box / GS as run_plasmon_n162_L50_E15;
// only WP_EKIN_EV (15 -> 3.4) and WP_SIGMA_BOHR (5 -> 3) change.
//
// Predicted m=1 channel peak: 1.71 eV (kinematic / wrap) vs 3.59 eV (bath
// plasmon). Separated by 22 x FFT resolution. See
// docs/plans/jellium_plasmon_detection.md and
// docs/reports/plasmon-detection-verdict.md for the protocol.
// ============================================================================

#include "../shared/configs/plasmon_n162_L50_E3p4_varyv.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::Plasmon_N162_L50_E3p4_VaryV>(
        "run_plasmon_n162_L50_E3p4_varyv",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L50_cubic_N162_dx1p0");
}
