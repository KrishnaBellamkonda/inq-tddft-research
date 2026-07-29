// ============================================================================
// run_base_n514 - closed-shell high-density jellium (N=514, |G|^2 <= 16,
// density 2.380e-3 e/bohr^3, r_s ~= 4.64 bohr).
//
// Cell 60^3 bohr (cubic, periodic), sigma=0.53 A, E=200 eV WP launched +z
// from origin, 320 steps at dt=0.020 a.u. (single-pass), 20 LEED screens.
// GS loaded from save_gs/gs_L60_cubic_N514. Cfg = Base_HighN (overrides
// N_ELECTRONS, EXTRA_STATES, SCF_TOL_HA from Base; everything else
// inherited).
// ============================================================================

#include "../shared/configs/base_highN.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::Base_HighN>(
        "run_base_n514",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L60_cubic_N514");
}
