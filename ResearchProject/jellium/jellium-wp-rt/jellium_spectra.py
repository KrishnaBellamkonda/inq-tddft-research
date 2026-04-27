#!/usr/bin/env python3
"""jellium_spectra.py — extended-spectrum postprocess for the jellium WP-RT runs.

Same pipeline as the coronene framework
(`inq-stack/python/inqview/postprocess/observables.py::_extended_spectra`),
applied to each jellium run's `results/observables.csv`. For dipole_z,
current_z, and energy_total, builds three preprocessed signals
(raw-subtracted, mean-subtracted, linearly detrended), Hann-windows them,
zero-pads by ``pad_factor=4``, FFTs, and emits PNGs + CSVs in
compartmentalised current/, dipole/, energy/ subfolders under each run's
``results/spectra/`` (analysis) and ``results/raw_spectra/`` (CSVs).

Usage::

    python3 jellium_spectra.py [run_dir ...]

If no run_dir arguments are supplied, every ``run_*/`` directory beside
this file is processed.

This script intentionally lives here (not in inqview) because the jellium
runs use a flat ``results/observables.csv`` layout rather than the
coronene ``results/raw/observables/observables.csv`` layout — the
inqview observables phase expects the latter. We re-use the actual
spectrum-building functions from inqview directly.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve()
                       .parents[3] / "inq-stack" / "python"))

import pandas as pd

from inqview.postprocess._common import ensure_dir, need_rebuild
from inqview.postprocess.observables import (
    _build_variants,
    _hann_fft,
    _plot_compare,
    _plot_one_spectrum,
    _quantity_subfolder,
    _save_spectrum_csv,
)


def run_one(run_dir: Path, *, rebuild: bool = True,
            pad_factor: int = 4) -> None:
    csv = run_dir / "results" / "observables.csv"
    if not csv.exists():
        print(f"[skip] {run_dir.name}: no results/observables.csv")
        return

    out_dir = ensure_dir(run_dir / "results" / "spectra")
    raw_dir = ensure_dir(run_dir / "results" / "raw_spectra")

    df = pd.read_csv(csv)
    df.columns = df.columns.str.strip()

    if "time_au" not in df.columns or len(df) < 4:
        print(f"[skip] {run_dir.name}: time_au missing or too few rows")
        return
    t = df["time_au"].to_numpy()
    dt_au = float(t[1] - t[0])
    n = int(t.size)
    energy_max_ev = 200.0
    run_name = run_dir.name

    print(f"[ok] {run_name}: dt_au={dt_au:.6f}, N={n}")

    columns = ("dipole_z", "current_z", "energy_total")
    for col in columns:
        if col not in df.columns:
            continue
        sub = _quantity_subfolder(col)
        out_specs = ensure_dir(out_dir / sub)
        raw_specs = ensure_dir(raw_dir / sub)
        signal = df[col].to_numpy()
        variants = _build_variants(signal)
        per_variant_results: dict = {}
        for variant, processed in variants.items():
            res = _hann_fft(processed, dt_au, pad_factor=pad_factor)
            if res is None:
                continue
            per_variant_results[variant] = res
            freq_au, omega_au, energy_ev, amplitude = res

            out_csv = raw_specs / f"spectrum_{col}_{variant}.csv"
            if need_rebuild(out_csv, rebuild):
                _save_spectrum_csv(out_csv, freq_au, omega_au, energy_ev,
                                   amplitude, column=col, variant=variant)

            out_png = out_specs / f"spectrum_{col}_{variant}.png"
            if need_rebuild(out_png, rebuild):
                _plot_one_spectrum(out_png, freq_au, energy_ev, amplitude,
                                   run_name=run_name, column=col,
                                   variant=variant,
                                   energy_max_ev=energy_max_ev)

        if per_variant_results:
            out_cmp = out_specs / f"spectrum_{col}_compare.png"
            if need_rebuild(out_cmp, rebuild):
                _plot_compare(out_cmp, per_variant_results,
                              run_name=run_name, column=col,
                              energy_max_ev=energy_max_ev)


def main(argv: list[str] | None = None) -> int:
    here = Path(__file__).resolve().parent
    args = sys.argv[1:] if argv is None else argv
    if args:
        targets = [Path(a).resolve() for a in args]
    else:
        targets = sorted(p for p in here.iterdir()
                         if p.is_dir() and p.name.startswith("run_"))
    for t in targets:
        run_one(t, rebuild=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
