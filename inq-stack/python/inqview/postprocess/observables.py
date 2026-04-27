"""Phase: ``observables`` — time-domain + FFT plots from observables CSV.

Reads ``results/raw/observables/observables.csv`` and produces under
``results/analysis/observables/``:

* ``observables_summary.png``        — energy / current / dipole panels
* ``total_energy_vs_time.png``
* ``current_components_vs_time.png``
* ``dipole_components_vs_time.png``
* ``fft_total_energy.png``
* ``fft_current_x.png`` / ``_y.png`` / ``_z.png``
* ``dipole_spectrum_x.png`` / ``_y.png`` / ``_z.png``

Numerical FFT outputs go in ``results/raw/observables/`` per the spec
(``fft_total_energy.csv`` etc.).
"""

from __future__ import annotations

from pathlib import Path

from . import _common
from . import pipeline as _pipeline


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    csv = results_dir / "raw" / "observables" / "observables.csv"
    if not csv.exists():
        _pipeline.skip(f"observables.csv missing at {csv}")

    out_dir = _common.ensure_dir(results_dir / "analysis" / "observables")
    raw_dir = _common.ensure_dir(results_dir / "raw" / "observables")

    from .. import (
        FourierTransform,
        load_observables,
        plot_all_energy_components_vs_time,
        plot_current_vs_time,
        plot_dipole_vs_time,
        plot_observables_summary,
        plot_spectrum,
        plot_total_energy_vs_time,
    )

    df = load_observables(csv)
    notes: dict = {"out_dir": str(out_dir), "n_rows": len(df)}

    # Summary panel (single 3-row figure)
    out = out_dir / "observables_summary.png"
    if _common.need_rebuild(out, rebuild):
        plot_observables_summary(csv, out)
    notes["summary"] = str(out)

    # Per-quantity plots. total_energy_vs_time.png contains ONLY E_total
    # (TODO 1e); all_energies_vs_time.png contains every component.
    for fn, name in [
        (plot_total_energy_vs_time,         "total_energy_vs_time.png"),
        (plot_all_energy_components_vs_time, "all_energies_vs_time.png"),
        (plot_current_vs_time,              "current_components_vs_time.png"),
        (plot_dipole_vs_time,               "dipole_components_vs_time.png"),
    ]:
        out = out_dir / name
        if _common.need_rebuild(out, rebuild):
            import matplotlib.pyplot as plt
            fig = fn(csv)
            fig.savefig(out)
            plt.close(fig)

    # FFT spectra
    ft = FourierTransform()
    spectra: dict[str, str] = {}

    def _fft_one(col: str, raw_csv_name: str, png_name: str):
        if col not in df.columns:
            return
        result = ft.transform_column(df, col)
        # Save numerical FFT to raw/observables/
        out_csv = raw_dir / raw_csv_name
        if _common.need_rebuild(out_csv, rebuild):
            import numpy as np
            arr = np.column_stack([result.frequency_au, result.amplitude])
            header = f"frequency_au,amplitude  ({col})"
            np.savetxt(out_csv, arr, delimiter=",", header=header)
        # Plot
        out_png = out_dir / png_name
        if _common.need_rebuild(out_png, rebuild):
            plot_spectrum(result, out_png)
        spectra[col] = str(out_png)

    _fft_one("energy_total", "fft_total_energy.csv",  "fft_total_energy.png")
    _fft_one("current_x",    "fft_current_x.csv",     "fft_current_x.png")
    _fft_one("current_y",    "fft_current_y.csv",     "fft_current_y.png")
    _fft_one("current_z",    "fft_current_z.csv",     "fft_current_z.png")
    _fft_one("dipole_x",     "dipole_spectrum_x.csv", "dipole_spectrum_x.png")
    _fft_one("dipole_y",     "dipole_spectrum_y.csv", "dipole_spectrum_y.png")
    _fft_one("dipole_z",     "dipole_spectrum_z.csv", "dipole_spectrum_z.png")

    notes["spectra"] = spectra
    return notes
