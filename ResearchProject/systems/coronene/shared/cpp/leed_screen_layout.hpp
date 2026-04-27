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

// Per-screen physics window:
//   * start_au = time when the WP centroid arrives at the screen.
//   * end_au   = time when the trailing edge has cleared the screen
//                (centroid + sigma for forward screens; centroid - sigma
//                 [back-scattered] for backscattering screens).
//
// Forward screen  (z_screen <  b): centroid arrives at (b - z)/|k|;
//                                  centroid + sigma at (b + sigma - z)/|k|.
// Backscatter     (z_screen >= b): centroid arrives at (b + z)/|k|;
//                                  centroid + sigma at (b + sigma + z)/|k|.
//
// Convention: WP starts at z = +b moving in -z at speed |k|. The molecule
// is at z = 0. A scattered packet rebounding from the molecule at t = b/|k|
// moves in +z at the same speed; the model treats backscattering screens
// as receiving that rebounded centroid plus its sigma envelope.
struct ScreenWindow {
    double t_start_au;
    double t_end_au;
    bool   is_back;     // true = backscattering side (z_screen >= b)
};

inline constexpr ScreenWindow compute_screen_window(double z_screen, double b,
                                                    double sigma, double k0) {
    if (z_screen < b) {
        return ScreenWindow{
            (b - z_screen) / k0,
            (b + sigma - z_screen) / k0,
            false,
        };
    }
    return ScreenWindow{
        (b + z_screen) / k0,
        (b + sigma + z_screen) / k0,
        true,
    };
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
