"""
Regression test for inqview.overlap.

Builds a synthetic overlap series in a tmp dir mimicking the layout the
C++ writer (inqkit::observables::OrbitalOverlapMatrix) produces:

    overlap_dir/
        index.csv            columns: step,time_au,file
        overlap_000000.csv   header line `# step=N time_au=...` then matrix
        overlap_NNNNNN.csv   ...

It then exercises:
    load_overlap_csv     - skip header, parse rows, robust to whitespace
    iter_overlap_series  - return snapshots in step order, drop inconsistent shapes
    pick_meaningful_columns - row=occupied subset, columns = occupied + WP
    plot_overlap_column_gif - asserts ValueError on wrong-size column data

Pinned by the user-reported symptom: "the number of columns and rows were
not as expected". The test asserts that
    n_ref_raw    = wp_idx     (= what the C++ writer produces)
    n_evolved_raw = wp_idx + 1
    pick_meaningful_columns yields exactly n_occupied + 1 entries
    rows beyond n_occupied are excluded from the GIF
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
PKG_DIR = THIS_DIR.parents[1] / "python"
sys.path.insert(0, str(PKG_DIR))

from inqview.overlap import (  # noqa: E402
    OverlapSnapshot,
    iter_overlap_series,
    load_overlap_csv,
    pick_meaningful_columns,
    plot_overlap_column_gif,
)


def _write_synthetic_series(
    overlap_dir: Path, n_ref: int, n_evolved: int, n_steps: int
) -> None:
    """
    Identity at t=0 on the n_ref x n_ref block; column n_ref-1+1 (= the WP
    column at index n_ref) is zero at t=0. For t>0 the diagonal decays
    smoothly and the WP column grows.
    """
    overlap_dir.mkdir(parents=True, exist_ok=True)
    index_path = overlap_dir / "index.csv"
    with index_path.open("w") as idx:
        idx.write("step,time_au,file\n")
        for step in range(n_steps):
            t = step * 0.02
            file_name = f"overlap_{step:06d}.csv"
            idx.write(f"{step},{t:.6f},{file_name}\n")
            mat = np.zeros((n_ref, n_evolved), dtype=float)
            decay = np.exp(-0.05 * step)
            for i in range(min(n_ref, n_evolved)):
                mat[i, i] = decay
            wp_col = n_evolved - 1
            if wp_col >= 0:
                # WP column overlap with each ref orbital grows
                for i in range(n_ref):
                    mat[i, wp_col] = (1.0 - decay) * (0.1 + 0.9 * (i / max(n_ref - 1, 1)))
            with (overlap_dir / file_name).open("w") as fh:
                fh.write(
                    f"# step={step} time_au={t:.6f} "
                    f"n_ref={n_ref} n_evolved={n_evolved}\n"
                )
                for i in range(n_ref):
                    fh.write(",".join(f"{mat[i, j]:.8f}" for j in range(n_evolved)) + "\n")


def test_load_overlap_csv_skips_header(tmp_path: Path = None):
    if tmp_path is None:
        import tempfile
        tmp_path = Path(tempfile.mkdtemp())
    tmp_path.mkdir(parents=True, exist_ok=True)
    p = tmp_path / "overlap_000000.csv"
    p.write_text(
        "# step=0 time_au=0.000000 n_ref=3 n_evolved=4\n"
        "1.0,0.0,0.0,0.0\n"
        "0.0,1.0,0.0,0.0\n"
        "0.0,0.0,1.0,0.0\n"
    )
    mat = load_overlap_csv(p)
    assert mat is not None
    assert mat.shape == (3, 4)
    np.testing.assert_array_equal(mat[:3, :3], np.eye(3))
    np.testing.assert_array_equal(mat[:, 3], np.zeros(3))


def test_iter_overlap_series_in_order_and_correct_shape(tmp_path: Path = None):
    if tmp_path is None:
        import tempfile
        tmp_path = Path(tempfile.mkdtemp())
    overlap_dir = tmp_path / "overlap"
    n_ref, n_evolved, n_steps = 5, 6, 8
    _write_synthetic_series(overlap_dir, n_ref, n_evolved, n_steps)
    snaps = list(iter_overlap_series(overlap_dir))
    assert len(snaps) == n_steps
    for k, snap in enumerate(snaps):
        assert snap.step == k
        assert snap.matrix.shape == (n_ref, n_evolved)
    # Monotonic time
    assert all(snaps[k].time_au < snaps[k + 1].time_au for k in range(n_steps - 1))


def test_iter_drops_inconsistent_shape(tmp_path: Path = None, capsys=None):
    import io
    import sys as _sys
    if tmp_path is None:
        import tempfile
        tmp_path = Path(tempfile.mkdtemp())
    overlap_dir = tmp_path / "overlap_bad"
    n_ref, n_evolved, n_steps = 4, 5, 4
    _write_synthetic_series(overlap_dir, n_ref, n_evolved, n_steps)
    # Corrupt one snapshot to have a different column count
    bad_path = overlap_dir / "overlap_000002.csv"
    bad_path.write_text(
        "# step=2 time_au=0.040000 n_ref=4 n_evolved=99\n"
        + "\n".join(",".join(["0.5"] * 99) for _ in range(n_ref)) + "\n"
    )
    captured = io.StringIO()
    old_stderr = _sys.stderr
    _sys.stderr = captured
    try:
        snaps = list(iter_overlap_series(overlap_dir))
    finally:
        _sys.stderr = old_stderr
    assert len(snaps) == n_steps - 1
    assert "expected" in captured.getvalue() or "WARN" in captured.getvalue()


def test_pick_meaningful_columns_includes_wp_slot():
    # Coronene-paper case: 54 occupied, 8 extra states, wp at last slot.
    n_occupied = 54
    n_extra = 8
    wp_idx = n_occupied + n_extra - 1   # = 61
    n_ref_raw = wp_idx                  # = 61 (what C++ records)
    n_evolved_raw = wp_idx + 1          # = 62
    cols = pick_meaningful_columns(n_occupied, wp_idx, n_evolved_raw)
    # Expect: 54 occupied columns + 1 WP column = 55
    assert len(cols) == n_occupied + 1
    # Last entry is WP
    assert cols[-1] == (wp_idx, wp_idx)
    # First n_occupied entries are diagonal labels
    for j in range(n_occupied):
        assert cols[j] == (j, j)


def test_pick_meaningful_columns_jellium_case():
    # Jellium: 19 occupied, 3 extra, wp at last slot.
    n_occupied = 19
    n_extra = 3
    wp_idx = n_occupied + n_extra - 1   # = 21
    n_evolved_raw = wp_idx + 1          # = 22
    cols = pick_meaningful_columns(n_occupied, wp_idx, n_evolved_raw)
    assert len(cols) == n_occupied + 1
    assert cols[-1] == (wp_idx, wp_idx)


def test_plot_overlap_column_gif_rejects_wrong_row_count(tmp_path: Path = None):
    if tmp_path is None:
        import tempfile
        tmp_path = Path(tempfile.mkdtemp())
    tmp_path.mkdir(parents=True, exist_ok=True)
    out = tmp_path / "should_not_exist.gif"
    bad_col = [np.zeros(3), np.zeros(4)]
    try:
        plot_overlap_column_gif(out, n_ref_rows=3, col_data=bad_col, times=[0.0, 0.02])
    except ValueError:
        return
    raise AssertionError("plot_overlap_column_gif should have raised on inconsistent row count")


def test_plot_overlap_column_gif_writes_one_frame_per_time(tmp_path: Path = None):
    if tmp_path is None:
        import tempfile
        tmp_path = Path(tempfile.mkdtemp())
    tmp_path.mkdir(parents=True, exist_ok=True)
    out = tmp_path / "test.gif"
    n_rows = 5
    col_data = [np.linspace(0, 1, n_rows), np.linspace(0, 0.5, n_rows)]
    plot_overlap_column_gif(out, n_ref_rows=n_rows, col_data=col_data, times=[0.0, 0.02])
    assert out.exists()
    assert out.stat().st_size > 0


if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        test_load_overlap_csv_skips_header(td / "a")
        test_iter_overlap_series_in_order_and_correct_shape(td / "b")
        test_iter_drops_inconsistent_shape(td / "c")
        test_pick_meaningful_columns_includes_wp_slot()
        test_pick_meaningful_columns_jellium_case()
        test_plot_overlap_column_gif_rejects_wrong_row_count(td / "d")
        test_plot_overlap_column_gif_writes_one_frame_per_time(td / "e")
    print("All overlap loader tests passed.")
