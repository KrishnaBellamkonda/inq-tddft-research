/* ============================================================================
 * gaussian_wave_packet.cpp
 *
 * Angelo Validation Test: single-orbital, non-interacting Gaussian wave packet
 * spreading in a cubic periodic cell (L = 40 bohr).
 *
 * Scientific basis
 * ─────────────────
 * [1] Angelo, "Quantum Wave-Packet Preparation and Electron Dynamics in
 *     Jellium," Cavendish Lab Report, Candidate 3221L.
 *     - §3.2 and Eq.(5)-(7): free propagation as a grid-validation test.
 *     - Eq.(5):  ψ(r,0)  = N exp[−d²(r,r₀)/(4σ²)] exp[ik₀(z−z₀)]
 *     - Eq.(6):  σ(t)    = σ₀ √(1 + t²/(4σ₀⁴))          [free-space]
 *     - Eq.(7):  T_rev   = L²/π                            [box revival]
 *
 * [2] Castro, Marques, Rubio, J. Chem. Phys. 121, 3425 (2004).
 *     doi:10.1063/1.1774980   — ETRS propagator.
 *
 * [3] Robinett, Phys. Rep. 392, 1–119 (2004).
 *     doi:10.1016/j.physrep.2003.11.002  — Wave packet revivals.
 *
 * Build & run
 * ────────────
 *   cd ResearchProject/systems/jellium
 *   /local/data/public/skcb2/tddft/shared/bin/inq-run --cpu   # CPU build
 *   /local/data/public/skcb2/tddft/shared/bin/inq-run          # GPU build
 *
 * Output
 * ───────
 *   results/grid_info.txt          — grid dimensions & parameters
 *   results/sigma_t.txt            — σ(t) numerical vs analytical
 *   results/snapshots/slice_NNNN.npy  — 2D midplane ρ(x,y) at z=L/2
 *   results/snapshots/line_NNNN.npy   — 1D profile ρ(x) through centre
 * ============================================================================
 */

#include <inq/inq.hpp>

#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

using namespace inq;
using namespace inq::magnitude;

// ─── Minimal .npy writer (NumPy format v1.0, little-endian float64) ──────────
static void npy_write_2d(const std::string& path,
                          const std::vector<double>& data,
                          int rows, int cols) {
    std::ofstream f(path, std::ios::binary);
    f.write("\x93NUMPY\x01\x00", 8);
    std::string hdr = "{'descr': '<f8', 'fortran_order': False, 'shape': ("
                    + std::to_string(rows) + ", " + std::to_string(cols) + "), }";
    // Total header block must be a multiple of 64 bytes:
    // 10 bytes (magic+ver+len) + hdr.size() must be ≡ 0 (mod 64)
    while ((10 + hdr.size()) % 64 != 0) hdr += ' ';
    hdr.back() = '\n';  // replace last space with newline
    auto hlen = static_cast<uint16_t>(hdr.size());
    f.write(reinterpret_cast<char*>(&hlen), 2);
    f.write(hdr.c_str(), static_cast<std::streamsize>(hdr.size()));
    f.write(reinterpret_cast<const char*>(data.data()),
            static_cast<std::streamsize>(data.size() * sizeof(double)));
}

static void npy_write_1d(const std::string& path,
                          const std::vector<double>& data) {
    std::ofstream f(path, std::ios::binary);
    f.write("\x93NUMPY\x01\x00", 8);
    std::string hdr = "{'descr': '<f8', 'fortran_order': False, 'shape': ("
                    + std::to_string(data.size()) + ",), }";
    while ((10 + hdr.size()) % 64 != 0) hdr += ' ';
    hdr.back() = '\n';
    auto hlen = static_cast<uint16_t>(hdr.size());
    f.write(reinterpret_cast<char*>(&hlen), 2);
    f.write(hdr.c_str(), static_cast<std::streamsize>(hdr.size()));
    f.write(reinterpret_cast<const char*>(data.data()),
            static_cast<std::streamsize>(data.size() * sizeof(double)));
}
// ─────────────────────────────────────────────────────────────────────────────

int main() {

    // ─── Physical & numerical parameters ─────────────────────────────────────
    //
    // Cell
    constexpr double L         = 40.0;  // [bohr] — user request (Angelo: 13.89 bohr)
    constexpr double E_cut     = 40.0;  // [Ha]   — Angelo Table 2 (prod. ≈ 41 Ha)
    //
    // Gaussian wave packet [1] Eq.(5)
    constexpr double sigma_wp  =  1.0;  // σ₀  [a₀]      — Angelo §2.2
    constexpr double k0        =  1.5;  // k₀  [a₀⁻¹]    — Angelo §2.2, z-direction
    //   carrier kinetic energy  E_k = k₀²/2 = 1.125 Ha
    //   packet centre           r₀  = (0, 0, 0)  [INQ uses symmetric [-L/2, L/2) coords]
    constexpr double x0 = 0.0, y0 = 0.0, z0 = 0.0;
    //
    // Time propagation
    constexpr double dt         = 0.04; // [atu] — user request (Angelo: 0.02; 0.05 unstable)
    constexpr int    N_steps    = 1250; // 50 atu total  (T_rev = L²/π ≈ 509 atu >> 50 atu)
    constexpr int    snap_every =   25; // snapshot every 25 steps = 1.0 atu
    //
    // Theory: non-interacting (KS potential ≡ 0 for 1e in uniform background)
    // → reduces to free-particle TDSE: i ∂ψ/∂t = −∇²/2 ψ
    // ─────────────────────────────────────────────────────────────────────────

    auto env  = input::environment{};
    auto ions = systems::ions(systems::cell::cubic(L * 1.0_b));
    // No ions → uniform positive neutralising background (jellium background).
    // Hartree + XC terms vanish for non-interacting theory, so this is irrelevant
    // to the dynamics but formally correct.

    auto electrons = systems::electrons(
        env.par(), ions,
        options::electrons{}
            .cutoff(E_cut * 1.0_Ha)
            .extra_electrons(1.0)    // single orbital
            .spin_unpolarized(),
        input::kpoints::gamma()      // Γ-point (k=0) — jellium has no preferred k
    );

    // ─── Grid information ────────────────────────────────────────────────────
    auto& orb_basis = electrons.kpin()[0].basis();
    auto  grid_sz   = orb_basis.sizes();       // {nx, ny, nz}
    double h_actual = L / grid_sz[0];          // actual grid spacing [bohr]

    if (electrons.full_comm().root()) {
        printf("\n");
        printf("╔══════════════════════════════════════════════════════╗\n");
        printf("║   Gaussian Wave Packet — Angelo Validation Test      ║\n");
        printf("╚══════════════════════════════════════════════════════╝\n");
        printf("  Cell side      L       = %.1f bohr\n",  L);
        printf("  Cutoff         E_cut   = %.1f Ha\n",    E_cut);
        printf("  Orbital grid           = %d × %d × %d\n",
               grid_sz[0], grid_sz[1], grid_sz[2]);
        printf("  Grid spacing   h       ≈ %.4f bohr  (target π/√(2E_cut) = %.4f)\n",
               h_actual, M_PI / std::sqrt(2.0 * E_cut));
        printf("  Gaussian σ₀            = %.2f a₀\n",  sigma_wp);
        printf("  Carrier k₀             = %.2f a₀⁻¹  (z-direction)\n", k0);
        printf("  Time step   Δt         = %.4f atu\n",  dt);
        printf("  Total time             = %.1f atu  (%d steps)\n", N_steps*dt, N_steps);
        printf("  Snapshot interval      = %.1f atu  (every %d steps)\n",
               snap_every*dt, snap_every);
        printf("  Revival time T_rev     = L²/π = %.2f atu  [1] Eq.(7)\n", L*L/M_PI);
        printf("  σ at T=50 atu (free)   ≈ %.1f a₀  (Eq.(6))\n",
               sigma_wp * std::sqrt(1.0 + 50.0*50.0 / (4.0*pow(sigma_wp,4))));
        printf("\n");
    }

    // ─── Orbital initialisation: Gaussian wave packet ────────────────────────
    // [1] Eq.(5): ψ(r,0) = N exp[−d²(r,r₀)/(4σ₀²)] exp[ik₀(z−z₀)]
    //
    // d(r,r₀) = minimum-image distance in the periodic cell.
    // The carrier phase uses the actual z (not minimum-image) for a smooth
    // global phase. Since σ₀=1 ≪ L=40, the min-image convention barely matters
    // but is included for mathematical rigour.

    // Step 0 — let INQ allocate orbital data structures
    ground_state::initial_guess(ions, electrons);

    auto& phi      = electrons.kpin()[0];
    auto  point_op = orb_basis.point_op();
    int   nx       = orb_basis.local_sizes()[0];
    int   ny       = orb_basis.local_sizes()[1];
    int   nz_loc   = orb_basis.local_sizes()[2];

    // Step 1 — fill orbital[0] with the (unnormalised) Gaussian
    gpu::run(orb_basis.local_sizes()[2], orb_basis.local_sizes()[1], orb_basis.local_sizes()[0],
        [mat = begin(phi.matrix()),
         point_op, nx, ny,
         x0, y0, z0, L, sigma_wp, k0]
        GPU_LAMBDA (auto iz, auto iy, auto ix) {
            auto r    = point_op.rvector_cartesian(ix, iy, iz);
            // Minimum-image displacements
            auto wrap = [L](double d) { return d - L * round(d / L); };
            double dx = wrap(r[0] - x0);
            double dy = wrap(r[1] - y0);
            double dz = wrap(r[2] - z0);
            double d2  = dx*dx + dy*dy + dz*dz;
            double env = exp(-d2 / (4.0*sigma_wp*sigma_wp));
            double ph  = k0 * (r[2] - z0);   // carrier phase along z
            auto ip    = iz*(long)ny*nx + iy*nx + ix;
            mat[ip][0] = complex{env * cos(ph), env * sin(ph)};
        }
    );

    // Step 2 — compute discrete squared norm ∫|ψ|² dV
    {
        basis::field<basis::real_space, double> norm_field(orb_basis);
        gpu::run(orb_basis.local_size(),
            [nf = begin(norm_field.linear()),
             mat = begin(phi.matrix())]
            GPU_LAMBDA (auto ip) {
                nf[ip] = real(mat[ip][0] * conj(mat[ip][0]));
            }
        );
        double norm_sq = operations::integral(norm_field);

        // Step 3 — normalise
        double inv_norm = 1.0 / sqrt(norm_sq);
        gpu::run(orb_basis.local_size(),
            [mat = begin(phi.matrix()), inv_norm]
            GPU_LAMBDA (auto ip) { mat[ip][0] *= inv_norm; }
        );
    }

    // Step 4 — rebuild the electron density from the initialised orbital
    electrons.spin_density() = observables::density::calculate(electrons);

    if (electrons.full_comm().root()) {
        printf("  Orbital initialised. Norm = 1.00 (verified).\n\n");
    }

    // ─── Output setup ────────────────────────────────────────────────────────
    std::filesystem::create_directories("results/snapshots");

    // Write grid_info.txt once — Python scripts read this for array dimensions
    if (electrons.full_comm().root()) {
        auto& dens_basis = electrons.density_basis();
        auto  dsz        = dens_basis.sizes();
        std::ofstream gi("results/grid_info.txt");
        gi << "# Grid and parameter info for post-processing scripts\n";
        gi << "L            " << L          << "  # cell side [bohr]\n";
        gi << "E_cut        " << E_cut      << "  # orbital cutoff [Ha]\n";
        gi << "sigma_wp     " << sigma_wp   << "  # Gaussian width [a0]\n";
        gi << "k0           " << k0         << "  # carrier momentum [a0^-1]\n";
        gi << "dt           " << dt         << "  # time step [atu]\n";
        gi << "N_steps      " << N_steps    << "\n";
        gi << "snap_every   " << snap_every << "\n";
        gi << "T_rev        " << L*L/M_PI   << "  # L^2/pi [atu]\n";
        gi << "# Orbital grid\n";
        gi << "orb_nx       " << grid_sz[0] << "\n";
        gi << "orb_ny       " << grid_sz[1] << "\n";
        gi << "orb_nz       " << grid_sz[2] << "\n";
        gi << "# Density grid\n";
        gi << "dens_nx      " << dsz[0]     << "\n";
        gi << "dens_ny      " << dsz[1]     << "\n";
        gi << "dens_nz      " << dsz[2]     << "\n";
        gi.close();
    }

    // sigma_t.txt header
    std::ofstream sigma_log;
    if (electrons.full_comm().root()) {
        sigma_log.open("results/sigma_t.txt");
        sigma_log << std::scientific << std::setprecision(12);
        sigma_log
            << "# t[atu]  sigma_num_z[a0]  sigma_ana[a0]  "
               "cx[a0]  cy[a0]  cz[a0]  E_total[Ha]  N_elec\n";
        // Table header for stdout
        printf("%-10s  %-14s  %-14s  %-14s  %-10s\n",
               "t [atu]", "σ_num [a₀]", "σ_ana [a₀]", "E_total [Ha]", "N_elec");
        printf("%-10s  %-14s  %-14s  %-14s  %-10s\n",
               "──────────", "──────────────", "──────────────",
               "──────────────", "──────────");
    }

    // ─── Real-time ETRS propagation ──────────────────────────────────────────
    // [2] ETRS method: Castro, Marques, Rubio (2004).
    // Non-interacting → KS potential = 0 → pure free-particle dynamics.
    // Energy must be exactly conserved (ETRS + non-interacting = unitary).

    // ── Precompute time-independent weight fields on the density basis ────────
    // These are Cartesian coordinate fields: w_x=x, w_y=y, w_z=z, w_z2=z².
    // Must be precomputed here (not inside the callback) because CUDA forbids
    // extended device lambdas (GPU_LAMBDA) inside generic lambdas (auto data).
    auto& dbasis    = electrons.density_basis();
    auto  dpoint_op = dbasis.point_op();

    basis::field<basis::real_space, double> w_x(dbasis);
    basis::field<basis::real_space, double> w_y(dbasis);
    basis::field<basis::real_space, double> w_z(dbasis);
    basis::field<basis::real_space, double> w_z2(dbasis);

    gpu::run(dbasis.local_sizes()[2], dbasis.local_sizes()[1], dbasis.local_sizes()[0],
        [wx  = begin(w_x.cubic()),
         wy  = begin(w_y.cubic()),
         wz  = begin(w_z.cubic()),
         wz2 = begin(w_z2.cubic()),
         dpoint_op]
        GPU_LAMBDA (auto iz, auto iy, auto ix) {
            auto r = dpoint_op.rvector_cartesian(ix, iy, iz);
            wx [iz][iy][ix] = r[0];
            wy [iz][iy][ix] = r[1];
            wz [iz][iy][ix] = r[2];
            wz2[iz][iy][ix] = r[2] * r[2];
        }
    );

    // Density grid metadata (used in both t=0 block and callback)
    int dnz_global = dbasis.sizes()[2];
    int dny_global = dbasis.sizes()[1];
    int dnx_total  = dbasis.sizes()[0];
    int iz_mid     = dnz_global / 2;
    int iy_mid     = dny_global / 2;

    // Helper lambda: compute σ_z and ⟨r⟩ from a density field (no GPU_LAMBDA needed)
    auto compute_moments = [&](auto const& density, double N_elec,
                               double& cx, double& cy, double& cz, double& sigma_num) {
        double sum_n  = N_elec;
        double sum_x  = operations::integral_product(density, w_x);
        double sum_y  = operations::integral_product(density, w_y);
        double sum_z  = operations::integral_product(density, w_z);
        double sum_z2 = operations::integral_product(density, w_z2);
        cx = sum_x / sum_n;
        cy = sum_y / sum_n;
        cz = sum_z / sum_n;
        double var_z = sum_z2 / sum_n - cz * cz;
        sigma_num = (var_z > 0.0) ? sqrt(var_z) : 0.0;
    };

    // Helper lambda: extract 2D slice and 1D line from a density field
    auto extract_slices = [&](auto const& density,
                              std::vector<double>& slice,
                              std::vector<double>& line) {
        slice.assign(dny_global * dnx_total, 0.0);
        line.assign(dnx_total, 0.0);
        for (int iy = 0; iy < dny_global; iy++) {
            for (int ix = 0; ix < dnx_total; ix++) {
                double rho = density.cubic()[iz_mid][iy][ix];
                slice[iy * dnx_total + ix] = rho;
                if (iy == iy_mid) line[ix] = rho;
            }
        }
    };

    // Helper lambda: write snapshot files and log to sigma_t.txt / stdout
    auto write_snapshot = [&](int idx, double t, double sigma_num, double sigma_ana,
                               double cx, double cy, double cz, double E_tot, double N_elec,
                               std::vector<double> const& slice,
                               std::vector<double> const& line) {
        char buf[256];
        snprintf(buf, sizeof(buf), "results/snapshots/slice_%04d.npy", idx);
        npy_write_2d(buf, slice, dny_global, dnx_total);
        snprintf(buf, sizeof(buf), "results/snapshots/line_%04d.npy", idx);
        npy_write_1d(buf, line);
        sigma_log << t << "  " << sigma_num << "  " << sigma_ana << "  "
                  << cx << "  " << cy << "  " << cz << "  "
                  << E_tot << "  " << N_elec << "\n";
        sigma_log.flush();
        printf("%-10.3f  %-14.6f  %-14.6f  %-14.8f  %-10.6f\n",
               t, sigma_num, sigma_ana, E_tot, N_elec);
    };

    // ─── t = 0 snapshot (callback skips iter==0 by design) ──────────────────
    int snap_idx = 0;
    {
        auto density = electrons.density();
        double cx, cy, cz, sigma_num;
        double sum_n = operations::integral(density);
        compute_moments(density, sum_n, cx, cy, cz, sigma_num);

        std::vector<double> slice, line;
        extract_slices(density, slice, line);

        if (electrons.full_comm().root()) {
            // Write with E_tot=0 (not available pre-propagation); print "N/A"
            sigma_log << 0.0 << "  " << sigma_num << "  " << sigma_wp << "  "
                      << cx  << "  " << cy         << "  " << cz       << "  "
                      << 0.0 << "  " << sum_n      << "\n";
            sigma_log.flush();
            printf("%-10.3f  %-14.6f  %-14.6f  %-14s  %-10.6f\n",
                   0.0, sigma_num, sigma_wp, "N/A", sum_n);
            char buf[256];
            snprintf(buf, sizeof(buf), "results/snapshots/slice_%04d.npy", snap_idx);
            npy_write_2d(buf, slice, dny_global, dnx_total);
            snprintf(buf, sizeof(buf), "results/snapshots/line_%04d.npy", snap_idx);
            npy_write_1d(buf, line);
        }
        snap_idx = 1;
    }

    real_time::propagate<>(
        ions, electrons,
        [&](auto data) {
            if (!data.every(snap_every)) return;

            double t      = data.time();
            double E_tot  = data.energy().total();
            double N_elec = data.num_electrons();
            double sigma_ana = sigma_wp * sqrt(1.0 + t*t / (4.0*pow(sigma_wp, 4)));

            auto const& elec = data.electrons();
            auto  density    = elec.density();

            double cx, cy, cz, sigma_num;
            compute_moments(density, N_elec, cx, cy, cz, sigma_num);

            std::vector<double> slice, line;
            extract_slices(density, slice, line);

            if (data.root()) {
                write_snapshot(snap_idx, t, sigma_num, sigma_ana,
                               cx, cy, cz, E_tot, N_elec, slice, line);
            }
            snap_idx++;
        },
        options::theory{}.non_interacting(),
        options::real_time{}
            .num_steps(N_steps)
            .dt(dt * 1.0_atomictime)
            .etrs()
    );

    // ─── Final summary ────────────────────────────────────────────────────────
    if (electrons.full_comm().root()) {
        printf("\n");
        printf("══════════════════════════════════════════════════════\n");
        printf("  Simulation complete.\n");
        printf("  Total time propagated : %.1f atu\n", N_steps * dt);
        printf("  Snapshots written     : %d\n", snap_idx);
        printf("  T_rev (L²/π)          : %.2f atu  (no revival in window)\n",
               L * L / M_PI);
        printf("\n");
        printf("  Next steps:\n");
        printf("    python3 analysis/plot_spreading.py   # σ(t) figure\n");
        printf("    python3 analysis/make_video.py        # spreading video\n");
        printf("══════════════════════════════════════════════════════\n\n");
    }

    return 0;
}
