"""
Logic test for the FFT-shift fix in inqkit::fields::density::total/orbital.

The C++ writer used to publish INQ's FFT-natural array (cell centre at array
index 0) with a metadata Origin of -L/2, which made the .vti reader place
the cell centre at the cell corner. The fix is to read

    hc[fft_shift_index(ix, nx)]
        ...

when filling the output array. This test mirrors the index transformation
in pure NumPy, builds a synthetic FFT-natural array whose only non-zero
voxel is at the cell centre (array index 0), applies the fix, and asserts
that the resulting array places the peak at the metadata centre
(index nx/2, ny/2, nz/2).

If this test fails, the C++ fft_shift_index formula in
inq-stack/include/inqkit/fields/density.hpp is wrong.
"""

from __future__ import annotations

import numpy as np


def fft_shift_index(idx: int, size: int) -> int:
    return (idx + (size + 1) // 2) % size


def apply_writer_shift(arr: np.ndarray) -> np.ndarray:
    nx, ny, nz = arr.shape
    out = np.empty_like(arr)
    for ix in range(nx):
        sx = fft_shift_index(ix, nx)
        for iy in range(ny):
            sy = fft_shift_index(iy, ny)
            for iz in range(nz):
                sz = fft_shift_index(iz, nz)
                out[ix, iy, iz] = arr[sx, sy, sz]
    return out


def test_index_formula_at_centre_and_corner():
    for n in (4, 5, 8, 9, 120, 121, 200, 201):
        # The metadata centre is at output index n/2 (or (n-1)/2 for odd).
        # Under FFT-natural layout the cell centre lives at array index 0.
        # So the writer must map output_idx == n/2 -> array_idx == 0.
        assert fft_shift_index(n // 2, n) == 0, n
        # Output index 0 (= -L/2 in metadata) must map to array_idx == n/2
        # (which is the FFT-natural -L/2 position).
        # For odd n this is (n+1)/2 which is fine.
        assert fft_shift_index(0, n) == (n + 1) // 2 % n, n


def test_writer_shift_centres_a_centre_peak():
    # A synthetic INQ array where only the cell centre (FFT-natural index 0)
    # is non-zero.
    for shape in [(8, 8, 8), (120, 120, 200), (5, 7, 9)]:
        nx, ny, nz = shape
        inq_arr = np.zeros(shape)
        inq_arr[0, 0, 0] = 1.0  # density at cell centre
        out = apply_writer_shift(inq_arr)

        # After the shift, the only non-zero voxel must sit at the
        # metadata centre, which is (nx/2, ny/2, nz/2).
        peak = np.unravel_index(int(np.argmax(out)), out.shape)
        assert peak == (nx // 2, ny // 2, nz // 2), (shape, peak)
        assert out[nx // 2, ny // 2, nz // 2] == 1.0
        # And the corner of the output array must be zero.
        assert out[0, 0, 0] == 0.0


def test_writer_shift_matches_numpy_fftshift():
    rng = np.random.default_rng(0)
    for shape in [(4, 4, 4), (120, 120, 200), (5, 7, 9), (4, 5, 6)]:
        inq_arr = rng.standard_normal(shape)
        manual = apply_writer_shift(inq_arr)
        np.testing.assert_array_equal(manual, np.fft.fftshift(inq_arr))


if __name__ == "__main__":
    test_index_formula_at_centre_and_corner()
    test_writer_shift_centres_a_centre_peak()
    test_writer_shift_matches_numpy_fftshift()
    print("All FFT-shift logic tests passed.")
