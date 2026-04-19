// ============================================================================
// hf-gs-with-inqkit: Ground state of HF molecule, field I/O via inqkit
//
// System: hydrogen fluoride in a cubic 8 bohr finite cell, LDA, 30 Ry cutoff.
//         H(1e) + F(7e pseudopotential valence) = 8 electrons → 4 occupied,
//         HOMO = index 3.  Bond length 0.917 Å.
//
// Writes via inqkit:
//   results/density/           total electron density
//   results/orbital_density/   HOMO |ψ|² density (orbital index 3)
//   results/orbitals/          HOMO complex wavefunction (orbital index 3)
//
// Purpose: validates the inqkit write → inqview read → visualise pipeline
//          on a fast small-molecule system (~1 min GPU).
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
    //auto cell = systems::cell::cubic(8.0_bohr).finite();

    //systems::ions ions(cell);
    // HF bond length 0.917 Å; atoms centred at cell origin (same as Tutorial/HF/HF.cpp)
    //auto bond_length = 0.917_angstrom;
    //auto zero = 0.0_angstrom;

    //ions.insert("H", {zero, zero, -bond_length / 2});
    //ions.insert("F", {zero, zero,  bond_length / 2});

	auto L = 16.0_bohr;  // or 18–20 bohr if you want more vacuum
	auto cell = systems::cell::cubic(L).finite();

	systems::ions ions(cell);

	auto bond_length = 0.917_angstrom;
	auto cx = L / 2;
	auto cy = L / 2;
	auto cz = L / 2;

	auto x_center = 4.234_angstrom;
	auto y_center = 4.234_angstrom;
	auto z_center = 4.234_angstrom;

	ions.insert("H", {x_center, y_center, z_center - bond_length / 2});
	ions.insert("F", {x_center, y_center, z_center + bond_length / 2});


    std::cout << "\n=== hf-gs-with-inqkit ===\n";

    // 8 valence electrons → 4 occupied, HOMO = index 3
    auto electrons = systems::electrons(
        ions,
        options::electrons{}.cutoff(30.0_Ry),
        input::kpoints::gamma());

    ground_state::initial_guess(ions, electrons);

    auto gs = ground_state::calculate(
        ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(1e-6_Ha)
            .max_steps(1000)
            .broyden_mixing());

    std::cout << "  GS energy = " << gs.energy.total() << " Ha\n";

    // ── Write total electron density ─────────────────────────────────────────
    std::cout << "Writing total density...\n";
    auto rho_total = inqkit::fields::density::total(electrons);
    inqkit::io::RealField3DWriter(
        "results/density",
        {.field_name = "total_density", .include_meta = true},
        {.overwrite = true})
        .write(rho_total, "density_total");

    // ── Write HOMO orbital density (|ψ_3|²) ─────────────────────────────────
    std::cout << "Writing HOMO orbital density (index 3)...\n";
    auto rho_homo = inqkit::fields::density::orbital(electrons, 3);
    inqkit::io::RealField3DWriter(
        "results/orbital_density",
        {.field_name = "homo_density", .include_meta = true},
        {.overwrite = true})
        .write(rho_homo, "orbital_0003_density");

    // ── Write HOMO complex wavefunction ─────────────────────────────────────
    std::cout << "Writing HOMO wavefunction (index 3)...\n";
    auto psi_homo = inqkit::fields::orbital::wavefunction(electrons, 3);
    inqkit::io::ComplexField3DWriter(
        "results/orbitals",
        {.field_name = "homo_wavefunction", .include_meta = true},
        {.overwrite = true})
        .write(psi_homo, "orbital_0003");

    std::cout << "\nDone. Output in results/\n";
    return 0;
}
