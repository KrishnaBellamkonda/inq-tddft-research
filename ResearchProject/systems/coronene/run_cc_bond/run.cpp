// ============================================================================
// run_cc_bond — Off-center WP impact at a Type B (radial spoke) C-C bond.
//
// All parameters identical to run_base (σ=1.0 Bohr, E=200 eV, b=12 Bohr,
// cell 35×35×60 Bohr, dt=0.020 a.u.) except WP_CX_BOHR = 4.028 Bohr
// (the midpoint of the C₁-C₂ bond along the x-axis). Single-variable
// comparison with run_base isolates the effect of beam position on the
// LEED back-scattering pattern.
//
// GS used: checkpoints/gs_35x35x60_cut40 (same as run_base).
// ============================================================================

#include "../shared/configs/cc_bond_35x35x60.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return coronene::run_template::run_propagation<coronene::config::CC_bond_35x35x60>(
        "run_cc_bond",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/shared/geometry/coronene.xyz",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/checkpoints/gs_35x35x60_cut40");
}
