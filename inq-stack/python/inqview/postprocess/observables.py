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

    # ── Extended (preprocessed) spectra: 3 variants per quantity ──────────
    # For dipole_z, current_z, energy_total we build:
    #    A. raw-subtracted        s - s(0)
    #    B. mean-subtracted       s - <s>
    #    C. linearly detrended    scipy.signal.detrend(s, type='linear')
    # window each with a Hann window, FFT, convert to omega (a.u.) and eV,
    # save numerical CSVs alongside per-variant PNGs and a 3-curve overlay.
    notes["extended_spectra"] = _extended_spectra(
        df, out_dir, raw_dir, run_name, rebuild)
    return notes


# ──────────────────────────────────────────────────────────────────────────
# Extended-spectrum helpers
# ──────────────────────────────────────────────────────────────────────────

# E[eV] = 2π · f[a.u.] · 27.2114 (since ω[Ha] = 2π·f and 1 Ha = 27.21138625 eV)
_HA_TO_EV = 27.21138625


def _build_variants(signal):
    """Return dict {variant_name: 1-D processed signal (numpy array)}."""
    import numpy as np
    from scipy.signal import detrend
    s = np.asarray(signal, dtype=np.float64)
    return {
        "raw_subtracted":  s - s[0],
        "mean_subtracted": s - s.mean(),
        "detrended":       detrend(s, type="linear"),
    }


def _hann_fft(signal_processed, dt_au: float):
    """Apply Hann window and FFT. Returns (freq_au, omega_au, energy_ev,
    amplitude). amplitude is normalised by N (Hann sum is N/2 in the limit;
    we keep |X|/N to make spectra comparable across variants of the same
    length)."""
    import numpy as np
    s = np.asarray(signal_processed, dtype=np.float64)
    N = s.size
    if N < 4:
        return None
    win = np.hanning(N)
    sw = s * win
    spec = np.fft.rfft(sw)
    freq_au = np.fft.rfftfreq(N, d=dt_au)         # cycles / a.u.-time
    omega_au = 2.0 * np.pi * freq_au              # angular frequency [Ha]
    energy_ev = _HA_TO_EV * omega_au              # photon-energy axis
    amplitude = np.abs(spec) / N
    return freq_au, omega_au, energy_ev, amplitude


def _save_spectrum_csv(out_csv: Path, freq_au, omega_au, energy_ev,
                       amplitude, *, column: str, variant: str) -> None:
    import numpy as np
    arr = np.column_stack([freq_au, omega_au, energy_ev, amplitude])
    header = (f"freq_au,omega_au,energy_ev,amplitude  "
              f"(column={column}, variant={variant})")
    np.savetxt(out_csv, arr, delimiter=",", header=header)


def _plot_one_spectrum(out_png: Path, freq_au, energy_ev, amplitude, *,
                       run_name: str, column: str, variant: str,
                       energy_max_ev: float | None = None) -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
    if energy_max_ev is not None:
        mask = energy_ev <= energy_max_ev
        x = energy_ev[mask]; y = amplitude[mask]
    else:
        x = energy_ev; y = amplitude
    ax.plot(x, y, linewidth=1.0)
    ax.set_xlabel("energy (eV)")
    ax.set_ylabel(f"|FFT({column})|")
    ax.set_title(f"{run_name}: spectrum {column} ({variant})")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def _plot_compare(out_png: Path, results: dict, *,
                  run_name: str, column: str,
                  energy_max_ev: float | None = None) -> None:
    """Overlay the three variants on a single axes."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
    for variant, (_freq, _omega, energy_ev, amplitude) in results.items():
        if energy_max_ev is not None:
            mask = energy_ev <= energy_max_ev
            x = energy_ev[mask]; y = amplitude[mask]
        else:
            x = energy_ev; y = amplitude
        ax.plot(x, y, label=variant, linewidth=1.0)
    ax.set_xlabel("energy (eV)")
    ax.set_ylabel(f"|FFT({column})|")
    ax.set_title(f"{run_name}: spectrum {column} — variant comparison")
    ax.legend(fontsize="x-small")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def _quantity_subfolder(column: str) -> str:
    """Map a CSV column to its compartmentalised spectra subfolder.

    All current_* columns -> 'current/', dipole_* -> 'dipole/',
    energy_* -> 'energy/'.
    """
    if column.startswith("current"):
        return "current"
    if column.startswith("dipole"):
        return "dipole"
    if column.startswith("energy"):
        return "energy"
    return "other"


def _extended_spectra(df, out_dir: Path, raw_dir: Path,
                      run_name: str, rebuild: bool) -> dict:
    """Run the 3-variant spectrum pipeline for dipole_z, current_z,
    energy_total. Outputs are compartmentalised into per-quantity
    subfolders so dipole / current / energy spectra don't crowd a single
    directory:

        results/analysis/observables/spectra/<dipole|current|energy>/
            spectrum_<col>_<variant>.png
            spectrum_<col>_compare.png
        results/raw/observables/spectra/<dipole|current|energy>/
            spectrum_<col>_<variant>.csv
    """
    base_out_specs = _common.ensure_dir(out_dir / "spectra")
    base_raw_specs = _common.ensure_dir(raw_dir / "spectra")

    # Time step from observables.csv (assumed uniform; the C++ writer
    # appends every step at dt_au, so dt = time_au[1] - time_au[0]).
    if "time_au" not in df.columns or len(df) < 4:
        return {"skipped": "time_au column missing or too few rows"}
    t = df["time_au"].to_numpy()
    if t.size < 4:
        return {"skipped": "fewer than 4 time samples"}
    dt_au = float(t[1] - t[0])

    # Cap the displayed energy range. 200 eV easily covers all coronene
    # KS-orbital eigenvalue differences; the spectra at higher energies
    # are dominated by FFT noise from the finite Hann window.
    energy_max_ev = 200.0

    columns = ("dipole_z", "current_z", "energy_total")
    notes: dict = {"dt_au": dt_au, "n": int(t.size), "columns": []}

    for col in columns:
        if col not in df.columns:
            continue
        notes["columns"].append(col)
        sub = _quantity_subfolder(col)
        out_specs = _common.ensure_dir(base_out_specs / sub)
        raw_specs = _common.ensure_dir(base_raw_specs / sub)
        signal = df[col].to_numpy()
        variants = _build_variants(signal)
        per_variant_results: dict = {}
        for variant, processed in variants.items():
            res = _hann_fft(processed, dt_au)
            if res is None:
                continue
            per_variant_results[variant] = res
            freq_au, omega_au, energy_ev, amplitude = res

            out_csv = raw_specs / f"spectrum_{col}_{variant}.csv"
            if _common.need_rebuild(out_csv, rebuild):
                _save_spectrum_csv(out_csv, freq_au, omega_au, energy_ev,
                                   amplitude, column=col, variant=variant)

            out_png = out_specs / f"spectrum_{col}_{variant}.png"
            if _common.need_rebuild(out_png, rebuild):
                _plot_one_spectrum(out_png, freq_au, energy_ev, amplitude,
                                   run_name=run_name, column=col,
                                   variant=variant,
                                   energy_max_ev=energy_max_ev)

        # Three-curve overlay so the variants can be compared at a glance.
        if per_variant_results:
            out_cmp = out_specs / f"spectrum_{col}_compare.png"
            if _common.need_rebuild(out_cmp, rebuild):
                _plot_compare(out_cmp, per_variant_results,
                              run_name=run_name, column=col,
                              energy_max_ev=energy_max_ev)

    return notes
