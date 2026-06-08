// ============================================================================
// run_plasmon_n162_L50_E15 — high-velocity WP designed to resonantly excite
// the m=1 axial plasmon mode of the L=50, N=162 closed-shell jellium box.
//
// E_WP = 15.0 eV  ⇒  k0 ≈ 1.05 Bohr⁻¹  ⇒  v ≈ 1.05 a.u.
//
// The bath is the same as run_base_n162_L50_E1p5; only the WP kinetic
// energy and propagation length change. Reuses the existing closed-shell
// GS at gs_L50_cubic_N162_dx1p0.
//
// Predicted m=1 plasmon at ℏω(q1) = 3.59 eV, period T_p = 1.15 fs.
// See `docs/plans/jellium_plasmon_detection.md` for the validated numbers.
// ============================================================================

#include "../shared/configs/plasmon_n162_L50_E15.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::Plasmon_N162_L50_E15>(
        "run_plasmon_n162_L50_E15",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L50_cubic_N162_dx1p0");
}
