// ============================================================================
// run_base_n138 - closed-shell low-energy jellium (N=138, |G|^2 <= 6,
// density 6.39e-4 e/bohr^3, r_s ~= 7.21 bohr).
//
// Cell 60^3 bohr (cubic, periodic), sigma=0.53 A, E=100 eV WP launched +z
// from origin, 320 steps at dt=0.020 a.u., 20 LEED screens.
// Spacing 0.55 bohr (Nyquist k_max = pi/dx = 5.71 Bohr^-1; resolves WP up
// to k_0 + 3 sigma_k = 2.71 + 2.99 = 5.71 Bohr^-1 - exactly at the edge).
// GS loaded from save_gs/gs_L60_cubic_N138_dx0p55. Cfg = Base_N138.
// ============================================================================

#include "../shared/configs/base_n138.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::Base_N138>(
        "run_base_n138",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L60_cubic_N138_dx0p55");
}
