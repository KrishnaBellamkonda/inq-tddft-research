// ============================================================================
// shared/cpp/leed_screen_layout.hpp
//
// Deterministic 20-screen z-layout, parameterised by Lz, used by every
// propagation run. Layout matches the existing run_07_paper_replica /
// run_propagate_paper_replica convention so the post-processor and the
// coordinate-check helper find screens at predictable positions.
//
//   * Two boundary screens at +/- (Lz/2 - 1) Bohr.
//   * Eighteen interior screens evenly spaced over [Lz/2 - 2.5, +Lz/2 - 2.5]
//     with a small parity-dependent jitter so multiple screens never coincide
//     on the same FFT grid plane (which would alias).
//
// The layout is symmetric in z so screen 0 is at -Lz/2 + 1 (same side as the
// outgoing/transmitted WP) and screen N-1 is at +Lz/2 - 1 (incident side).
// ============================================================================
#pragma once

#include <array>
#include <cstddef>
#include <string>

namespace coronene::layout {

constexpr int N_SCREENS = 20;

// Per-screen physics window (Phase-3 logic — corrected for LEED purity):
//
// The aim of windowing is to suppress the unscattered Gaussian-WP contribution
// at the screen plane and accumulate only the diffracted / scattered density.
//
// Forward screens (z_screen < 0, transmission side):
//   * t_start = max(0, (b + sigma - z_screen) / |k|)
//                — when the WP trailing edge has cleared the screen plane on
//                  its way down toward the molecule. Any density at z_screen
//                  after this time is the transmitted+diffracted contribution
//                  (the molecule has spread the WP into a non-Gaussian).
//   * t_end   = total_time_au (= N_steps * dt by construction of
//                compute_n_steps), i.e. integrate to end of run.
//
// Backscattering screens (z_screen >= 0, reflection side):
//   * t_start = max(0, (b + sigma - z_screen) / |k|)
//                — same expression. For screens above the initial WP
//                  envelope (z_screen > b + sigma) this is <= 0 and clamps
//                  to 0 (the screen was never under the incoming WP).
//   * t_end   = (b + Lz/2 - sigma) / |k|
//                — when the back-scattered forward leading edge
//                  (rebound centroid + sigma in +z direction) reaches the
//                  +Lz/2 box face, just before periodic-boundary wrap-around
//                  contaminates the signal.
//
// Demarcation: z_screen < 0 vs z_screen >= 0 (the molecule plane).
//
// Convention: WP starts at z = +b moving in -z at speed |k|; molecule at z=0.
// A scattered packet rebounding from the molecule at t = b/|k| moves in +z
// at the same speed; the model treats backscattering screens as receiving
// that rebounded centroid plus its sigma envelope.
struct ScreenWindow {
    double t_start_au;
    double t_end_au;
    bool   is_back;     // true = backscattering side (z_screen >= 0)
};

inline constexpr ScreenWindow compute_screen_window(double z_screen, double b,
                                                    double sigma, double k0,
                                                    double lz_bohr,
                                                    double total_time_au) {
    // Constexpr-friendly max(0, x).
    const double t_start_raw = (b + sigma - z_screen) / k0;
    const double t_start = (t_start_raw > 0.0) ? t_start_raw : 0.0;
    if (z_screen < 0.0) {
        // Forward / transmission side: integrate to end of simulation.
        return ScreenWindow{ t_start, total_time_au, false };
    }
    // Backscattering side: integrate until the rebound leading edge reaches
    // the +Lz/2 boundary.
    const double t_end_back = (b + 0.5 * lz_bohr - sigma) / k0;
    return ScreenWindow{ t_start, t_end_back, true };
}

inline std::array<double, N_SCREENS> screen_z_positions(double lz_bohr) {
    std::array<double, N_SCREENS> z{};
    const double half = 0.5 * lz_bohr;
    z[0]              = -half + 1.0;
    z[N_SCREENS - 1]  = +half - 1.0;
    const int n_interior = N_SCREENS - 2;
    const double z_lo = -half + 2.5;
    const double z_hi = +half - 2.5;
    const double dz = (z_hi - z_lo) / static_cast<double>(n_interior - 1);
    for (int k = 0; k < n_interior; ++k) {
        const double jitter = ((k % 2) == 0 ? +0.07 : -0.13);
        z[k + 1] = z_lo + k * dz + jitter;
    }
    return z;
}

inline std::string screen_label(int k) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "screen_%02d", k);
    return std::string(buf);
}

inline std::string zero_pad6(int n) {
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%06d", n);
    return std::string(buf);
}

inline std::string zero_pad2(int n) {
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%02d", n);
    return std::string(buf);
}

}  // namespace coronene::layout
