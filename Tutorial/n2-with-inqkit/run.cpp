#include <inq/inq.hpp>

#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

using namespace inq;
using namespace inq::magnitude;

int main() {
  // Simple N2 simulation for testing
  auto half_bond = 1.10_angstrom / 2;
  auto zero = 0.0_angstrom;

  systems::cell cell = systems::cell::cubic(20.0_bohr).finite();

  systems::ions ions(cell);
  ions.insert("N", {zero, zero, half_bond});
  ions.insert("N", {zero, zero, -half_bond});

  systems::electrons electrons(ions, options::electrons{}.cutoff(80.0_Ry));

  ground_state::initial_guess(ions, electrons);
  auto gs = ground_state::calculate(ions, electrons, options::theory{}.pbe());

  // Build orbital 0 as a generic complex field and write it
  auto psi_orb0 = inqkit::fields::orbital::wavefunction(electrons, 0);

  inqkit::io::ComplexField3DWriter orbital_writer(
      "results/orbitals", {.field_name = "ks_orbital", .include_meta = true},
      {.overwrite = true});
  orbital_writer.write(psi_orb0, "orbital_0000");

  // Total electronic density
  auto rho_total = inqkit::fields::density::total(electrons);

  inqkit::io::RealField3DWriter density_writer(
      "results/density", {.field_name = "total_density", .include_meta = true},
      {.overwrite = true});
  density_writer.write(rho_total, "density_total");

  // Density of orbital 0
  auto rho_orb0 = inqkit::fields::density::orbital(electrons, 0);

  inqkit::io::RealField3DWriter orbital_density_writer(
      "results/orbital_density",
      {.field_name = "orbital_density", .include_meta = true},
      {.overwrite = true});
  orbital_density_writer.write(rho_orb0, "orbital_0000_density");

  return 0;
}
