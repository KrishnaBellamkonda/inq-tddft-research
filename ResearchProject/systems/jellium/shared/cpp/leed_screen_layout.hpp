// ============================================================================
// shared/cpp/leed_screen_layout.hpp  (jellium)
//
// Deterministic 20-screen z-layout for cubic-periodic jellium cells of edge
// L_BOHR, parameterised by L. Adapted from
// ResearchProject/systems/coronene/shared/cpp/leed_screen_layout.hpp.
//
// Coronene used a centred orthorhombic cell [-Lz/2, +Lz/2] with the molecule
// at z=0 and the WP launched from +z moving in -z. Jellium uses a periodic
// cubic cell [0, L] with the WP launched from z = wp_cz (= 5*sigma by default,
// near the z=0 face) moving in +z.
//
// The screen layout below preserves the coronene layout *shape* (two boundary
// screens just inside the box faces, 18 interior screens evenly spaced with
// a small parity-dependent jitter) but maps it to the [0, L] coordinate
// system. Screen 0 sits near z = +1 (the WP-launch side) and screen
// N-1 sits near z = L - 1 (the transmission side).
//
// `compute_screen_window` is documented in the coronene equivalent. It is
// preserved here in the same form (with the WP centroid at offset `b` and
// a forward / back demarcation), but the coronene meaning of `b` and the
// rebound model are physically tuned for a centred cell with a localised
// scatterer at z=0. For jellium (no localised scatterer; periodic cell)
// the window is informational; until a jellium-specific physical model is
// adopted, callers should use `WP_ENVELOPE_SIGMAS = full_propagation` or
// pass `T1_AU = 0, T2_AU = N_STEPS*DT_AU` to get the full-time accumulator
// behaviour.
//
// FOLLOW-UP: re-derive `compute_screen_window` for cubic-periodic jellium.
// Tracked in docs/plans/jellium_reorg.md §14.
// ============================================================================
#pragma once

#include <array>
#include <cstddef>
#include <cstdio>
#include <string>

namespace jellium::layout {

constexpr int N_SCREENS = 20;

// Per-screen physics window (placeholder, see header comment above).
//
// Coronene convention (preserved in port):
//   * Forward side: z_screen "ahead" of WP centroid in launch direction.
//   * Back side:    z_screen "behind" WP centroid (rebound model).
//
// In jellium with WP at z = b moving in +z direction:
//   forward  ↔ z_screen >  b
//   back     ↔ z_screen <= b
struct ScreenWindow {
    double t_start_au;
    double t_end_au;
    bool   is_back;     // true = back-side (WP must rebound to reach)
};

inline constexpr ScreenWindow compute_screen_window(double z_screen, double b,
                                                    double sigma, double k0,
                                                    double l_bohr,
                                                    double total_time_au,
                                                    double n_sigmas = 1.0) {
    const double envelope = n_sigmas * sigma;
    // Time after which the trailing edge of the +z-moving WP has cleared
    // z_screen (forward case) or after which the rebound packet's leading
    // edge would arrive (back case). For jellium with no localised scatterer
    // this is informational — the back model assumes a hard rebound at the
    // far face, which is unphysical for a periodic cell.
    const double t_start_raw = (z_screen - b - envelope) / k0;
    const double t_start = (t_start_raw > 0.0) ? t_start_raw : 0.0;
    if (z_screen > b) {
        // Forward / transmission: integrate to end of run.
        return ScreenWindow{ t_start, total_time_au, false };
    }
    // Back side (z_screen <= b): bound the window by the time at which the
    // rebound centroid (off the z = L face) plus its envelope reaches
    // z_screen on its way back. Placeholder, see header comment.
    const double t_end_back = (l_bohr - b - envelope + (b - z_screen)) / k0;
    return ScreenWindow{ t_start, t_end_back, true };
}

// Parameterised 20-screen layout for a cubic-periodic [0, L] cell.
//   * Two boundary screens at z = +1 (just inside the WP-launch face) and
//     z = L - 1 (just inside the transmission face).
//   * Eighteen interior screens evenly spaced over [2.5, L - 2.5] with the
//     same parity-dependent jitter used by coronene to avoid coinciding
//     with FFT grid planes.
inline std::array<double, N_SCREENS> screen_z_positions(double l_bohr) {
    std::array<double, N_SCREENS> z{};
    z[0]              = 1.0;
    z[N_SCREENS - 1]  = l_bohr - 1.0;
    const int n_interior = N_SCREENS - 2;
    const double z_lo = 2.5;
    const double z_hi = l_bohr - 2.5;
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

}  // namespace jellium::layout
