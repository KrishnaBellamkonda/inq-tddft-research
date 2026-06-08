// ============================================================================
// run_plasmon_n162_L50_E25 — Run E, 25 eV frequency-resolved loss-function
// companion to run_plasmon_n162_L50_E15 (15 eV) and
// run_plasmon_n162_L50_E3p4_varyv (3.4 eV). Same box / GS / propagator length
// (T_sim = 2000 a.u., dE = 0.086 eV); WP_EKIN_EV = 25, WP_SIGMA_BOHR = 3,
// WRITE_EVERY = 10 (high-cadence scalars). Used to test whether the extracted
// loss function L(q,w) is projectile-independent (medium property).
// See shared/configs/plasmon_n162_L50_E25.hpp for the full derivation.
// ============================================================================

#include "../shared/configs/plasmon_n162_L50_E25.hpp"
#include "../shared/cpp/run_template.hpp"

int main() {
    return jellium::run_template::run_propagation<jellium::config::Plasmon_N162_L50_E25>(
        "run_plasmon_n162_L50_E25",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L50_cubic_N162_dx1p0");
}
