// Free-particle 3D Gaussian wavepacket propagation via FFT
//
// Grid:       L = 40 bohr, N = 161, dx ≈ 0.2484 bohr (E_cut = 80 Ha)
// Wavepacket: sigma = 1.0 bohr, k0 = [0,0,0], r0 = [L/2, L/2, L/2]
//             ψ(r) = N * exp(-|r-r0|²/(4σ²)) * exp(+i k0·(r-r0))
//             N = (2πσ²)^(-3/4)  (from wavepacket.py, line 35)
// Evolution:  FFT → multiply each k-mode by exp(-i|k|²dt/2) → IFFT
// Reference:  Tutorial/angelo-jellium/jellium/wavepacket.py

#include <complex>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <fftw3.h>

// ── Grid ─────────────────────────────────────────────────────────────────────
static constexpr int    N    = 161;
static constexpr double L    = 40.0;          // bohr
static constexpr double dx   = L / N;
static constexpr long   Ntot = (long)N*N*N;

// ── Wavepacket defaults (from wavepacket.py) ──────────────────────────────
static constexpr double sigma = 1.0;
static constexpr double k0x   = 0.0, k0y = 0.0, k0z = 0.0;
static constexpr double r0x   = L/2, r0y = L/2, r0z = L/2;

// ── Time parameters ───────────────────────────────────────────────────────
static constexpr double DT         = 0.02;    // a.u.
static constexpr double T_MAX      = 12.0;    // a.u.
static constexpr int    SAVE_EVERY = 5;       // frames saved every 5 steps → dt_save = 0.1
static constexpr int    N_STEPS    = static_cast<int>(T_MAX / DT + 0.5);  // 600

// ── Coarse 3D output grid ─────────────────────────────────────────────────
static constexpr int STRIDE = 4;              // sample every 4th point → 41³

// ── Helpers ───────────────────────────────────────────────────────────────
inline long idx(int i, int j, int k) {
    return (long)i*N*N + (long)j*N + k;
}

// FFTW-ordered frequency index → integer wavenumber
inline int kfreq_int(int j) {
    return (j <= N/2) ? j : j - N;
}

// ─────────────────────────────────────────────────────────────────────────────
int main() {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Grid: N=" << N << "  dx=" << dx << " bohr  Ntot=" << Ntot << "\n";
    std::cout << "Propagating t=0 to " << T_MAX << " a.u. in " << N_STEPS << " steps\n\n";

    // ── Allocate ─────────────────────────────────────────────────────────
    fftw_complex* psi  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ntot);
    fftw_complex* psik = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ntot);

    // Plans (MEASURE for better performance; ESTIMATE also fine)
    std::cout << "Planning FFTs...\n";
    fftw_plan fwd = fftw_plan_dft_3d(N, N, N, psi,  psik, FFTW_FORWARD,  FFTW_ESTIMATE);
    fftw_plan bwd = fftw_plan_dft_3d(N, N, N, psik, psi,  FFTW_BACKWARD, FFTW_ESTIMATE);
    std::cout << "Done.\n\n";

    // ── Initialise wavepacket ────────────────────────────────────────────
    // N_norm = (2π σ²)^(-3/4)
    double N_norm = std::pow(2.0 * M_PI * sigma * sigma, -0.75);
    std::cout << "N_norm = " << N_norm << "\n";

    for (int i = 0; i < N; i++) {
        double xi = i * dx - r0x;
        for (int j = 0; j < N; j++) {
            double yj = j * dx - r0y;
            for (int k = 0; k < N; k++) {
                double zk  = k * dx - r0z;
                double r2  = xi*xi + yj*yj + zk*zk;
                double amp = N_norm * std::exp(-r2 / (4.0 * sigma * sigma));
                double phi = k0x*xi + k0y*yj + k0z*zk;
                long   id  = idx(i, j, k);
                psi[id][0] = amp * std::cos(phi);
                psi[id][1] = amp * std::sin(phi);
            }
        }
    }

    // ── Precompute |k|² for each grid point ──────────────────────────────
    double dk = 2.0 * M_PI / L;
    std::vector<double> k2(Ntot);
    for (int i = 0; i < N; i++) {
        double ki = kfreq_int(i) * dk;
        for (int j = 0; j < N; j++) {
            double kj = kfreq_int(j) * dk;
            for (int k = 0; k < N; k++) {
                double kk = kfreq_int(k) * dk;
                k2[idx(i,j,k)] = ki*ki + kj*kj + kk*kk;
            }
        }
    }

    // ── Open CSV ──────────────────────────────────────────────────────────
    std::ofstream csv("results/width_vs_time.csv");
    csv << "t,sigma_x,sigma_y,sigma_z,sigma_analytical\n";
    csv << std::fixed << std::setprecision(8);

    // ── Lambda: save one frame ────────────────────────────────────────────
    auto write_frame = [&](int frame, double t) {
        // --- Compute density moments ---
        double sum_rho = 0, sx = 0, sy = 0, sz = 0;
        double sx2 = 0, sy2 = 0, sz2 = 0;
        for (int i = 0; i < N; i++) {
            double xi = i * dx;
            for (int j = 0; j < N; j++) {
                double yj = j * dx;
                for (int k = 0; k < N; k++) {
                    double zk  = k * dx;
                    long   id  = idx(i, j, k);
                    double rho = psi[id][0]*psi[id][0] + psi[id][1]*psi[id][1];
                    sum_rho += rho;
                    sx  += xi * rho;   sy  += yj * rho;   sz  += zk * rho;
                    sx2 += xi*xi * rho; sy2 += yj*yj * rho; sz2 += zk*zk * rho;
                }
            }
        }
        // dV factors cancel in ratios; use sum_rho as denominator
        double cx = sx / sum_rho, cy = sy / sum_rho, cz = sz / sum_rho;
        double sw_x = std::sqrt(std::max(0.0, sx2/sum_rho - cx*cx));
        double sw_y = std::sqrt(std::max(0.0, sy2/sum_rho - cy*cy));
        double sw_z = std::sqrt(std::max(0.0, sz2/sum_rho - cz*cz));
        double sigma_a = sigma * std::sqrt(1.0 + t*t / (4.0 * std::pow(sigma, 4)));

        csv << t << "," << sw_x << "," << sw_y << "," << sw_z << "," << sigma_a << "\n";

        // --- Write 2D slice at z = N/2 ---
        char fname[128];
        std::snprintf(fname, sizeof(fname), "results/slice_t%03d.txt", frame);
        std::ofstream sf(fname);
        sf << "# t=" << t << " z_index=" << (N/2) << " N=" << N << " dx=" << dx << "\n";
        sf << std::scientific << std::setprecision(6);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                long   id  = idx(i, j, N/2);
                double rho = psi[id][0]*psi[id][0] + psi[id][1]*psi[id][1];
                sf << rho;
                if (j < N-1) sf << ' ';
            }
            sf << '\n';
        }

        // --- Write coarse 3D density (stride=4 → 41³) ---
        std::snprintf(fname, sizeof(fname), "results/density3d_t%03d.txt", frame);
        std::ofstream df(fname);
        int NC = 0;
        for (int i = 0; i < N; i += STRIDE) NC++;  // count actual samples per dim
        df << "# t=" << t << " N=" << N << " dx=" << dx
           << " stride=" << STRIDE << " NC=" << NC << "\n";
        df << std::scientific << std::setprecision(6);
        for (int i = 0; i < N; i += STRIDE)
            for (int j = 0; j < N; j += STRIDE)
                for (int k = 0; k < N; k += STRIDE) {
                    long   id  = idx(i, j, k);
                    double rho = psi[id][0]*psi[id][0] + psi[id][1]*psi[id][1];
                    df << rho << '\n';
                }

        std::cout << "  frame " << frame << "  t=" << t
                  << "  sigma_x=" << sw_x << "\n";
    };

    // ── Time loop ─────────────────────────────────────────────────────────
    int frame = 0;
    write_frame(frame, 0.0);

    double inv_Ntot = 1.0 / (double)Ntot;

    for (int step = 1; step <= N_STEPS; step++) {
        // Forward FFT: ψ(r) → ψ̂(k)
        fftw_execute(fwd);

        // Multiply each k-mode by exp(-i |k|² dt/2)
        double half_dt = 0.5 * DT;
        for (long id = 0; id < Ntot; id++) {
            double phase = -k2[id] * half_dt;
            double cp = std::cos(phase), sp = std::sin(phase);
            double re = psik[id][0], im = psik[id][1];
            psik[id][0] = re*cp - im*sp;
            psik[id][1] = re*sp + im*cp;
        }

        // Backward FFT + FFTW normalisation (IFFT = FFT† / N³)
        fftw_execute(bwd);
        for (long id = 0; id < Ntot; id++) {
            psi[id][0] *= inv_Ntot;
            psi[id][1] *= inv_Ntot;
        }

        if (step % SAVE_EVERY == 0) {
            write_frame(++frame, step * DT);
        }
    }

    csv.close();
    fftw_destroy_plan(fwd);
    fftw_destroy_plan(bwd);
    fftw_free(psi);
    fftw_free(psik);
    fftw_cleanup();

    std::cout << "\nDone. Wrote " << (frame+1) << " frames to results/\n";
    return 0;
}
