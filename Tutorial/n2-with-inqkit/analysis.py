from pathlib import Path
import numpy as np

from inqview.data import load_complex_field


def main() -> None:
    meta_path = Path("results/orbitals/orbital_0000.meta.txt")

    orbital = load_complex_field(meta_path)

    psi = orbital.values
    density = np.abs(psi) ** 2

    print("Loaded orbital successfully")
    print(f"field_name      : {orbital.field_name}")
    print(f"shape           : {orbital.shape}")
    print(f"spacing_bohr    : {orbital.spacing_bohr}")
    print(f"kpoint          : {orbital.kpoint}")
    print(f"orbital_index   : {orbital.orbital_index}")
    print(f"spin_index      : {orbital.spin_index}")
    print(f"dtype           : {psi.dtype}")
    print(f"norm-like sum   : {density.sum()}")

    center = tuple(s // 2 for s in orbital.shape)
    print(f"psi(center)     : {psi[center]}")
    print(f"|psi|^2(center) : {density[center]}")


if __name__ == "__main__":
    main()
