// ============================================================================
// coronene-gs-with-inqkit: Ground state of coronene C24H12, field I/O via inqkit
//
// System: same cell and geometry as run_004 (34.771 × 34.771 × 89.856 bohr,
//         finite orthorhombic, LDA, 40 Ha cutoff, no extra states).
//         C24H12: 108 valence electrons → 54 occupied KS orbitals, HOMO = 53.
//
// Writes via inqkit:
//   results/density/           total electron density
//   results/orbital_density/   HOMO |ψ|² density (orbital index 53)
//   results/orbitals/          HOMO complex wavefunction (orbital index 53)
//
// Purpose: validates the inqkit write → inqview read → visualise pipeline.
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;

int main() {
    // Cell dimensions from run_004 config.hpp
    // LX = LY = 18.4 Å = 34.771 bohr, LZ = 47.55 Å = 89.856 bohr
    auto cell = systems::cell::orthorhombic(
        34.771_bohr, 34.771_bohr, 89.856_bohr).finite();

    // Coronene geometry: centred at (Lx/2, Ly/2, Lz/2) — all coords in [0, L]
    auto ions = systems::ions::parse("coronene_centered.xyz", cell);

    std::cout << "\n=== coronene-gs-with-inqkit ===\n";
    std::cout << "  Atoms: " << ions.size() << "  (expect 36: 24 C + 12 H)\n";

    // extra_states(3) matches run_004 — improves Broyden mixing stability
    // 54 occupied + 3 buffer = 57 total; HOMO = index 53
    auto electrons = systems::electrons(
        ions, options::electrons{}.cutoff(40.0_Ha).extra_states(3));

    ground_state::initial_guess(ions, electrons);

    auto gs = ground_state::calculate(
        ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(1e-4_Ha)
            .max_steps(300)
            .broyden_mixing()
            .mixing_ndim(8)
            .mixing(0.1));

    std::cout << "  GS energy = " << gs.energy.total() << " Ha\n";

    // ── Write total electron density ─────────────────────────────────────────
    std::cout << "Writing total density...\n";
    auto rho_total = inqkit::fields::density::total(electrons);
    inqkit::io::RealField3DWriter(
        "results/density",
        {.field_name = "total_density", .include_meta = true},
        {.overwrite = true})
        .write(rho_total, "density_total");

    // ── Write HOMO orbital density (|ψ_53|²) ────────────────────────────────
    std::cout << "Writing HOMO orbital density (index 53)...\n";
    auto rho_homo = inqkit::fields::density::orbital(electrons, 53);
    inqkit::io::RealField3DWriter(
        "results/orbital_density",
        {.field_name = "homo_density", .include_meta = true},
        {.overwrite = true})
        .write(rho_homo, "orbital_0053_density");

    // ── Write HOMO complex wavefunction ─────────────────────────────────────
    std::cout << "Writing HOMO wavefunction (index 53)...\n";
    auto psi_homo = inqkit::fields::orbital::wavefunction(electrons, 53);
    inqkit::io::ComplexField3DWriter(
        "results/orbitals",
        {.field_name = "homo_wavefunction", .include_meta = true},
        {.overwrite = true})
        .write(psi_homo, "orbital_0053");

    std::cout << "\nDone. Output in results/\n";
    return 0;
}
