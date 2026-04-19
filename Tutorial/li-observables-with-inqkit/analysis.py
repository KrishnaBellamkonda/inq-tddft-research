"""
Li BCC 2x2x2 ionic-kick TDDFT — observable analysis.

Produces time-domain plots and frequency-domain spectra from
results/observables.csv.

Usage:
    python analysis.py
    python analysis.py --csv path/to/observables.csv
"""

import argparse
from pathlib import Path

import inqview
from inqview import FourierTransform, WindowSpec
from inqview.data import load_real_field, load_meta, infer_meta_path
from inqview.plots import plot_spectrum, plot_spectrum_summary, plot_density_slice

CSV_DEFAULT = Path("results/observables.csv")

# Frequency cutoff for spectrum plots: 1.0 Ha/hbar ≈ 27 eV — covers all
# Li phonon and plasmon frequencies while discarding high-frequency noise.
FREQ_MAX_AU = 1.0


def main(csv_path: Path) -> None:
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Observables CSV not found: {csv_path}\n"
            "Run 'inq-run' first to generate the simulation output."
        )

    out = csv_path.parent

    # ------------------------------------------------------------------
    # Ground-state density slices
    # ------------------------------------------------------------------
    gs_raw  = out / "gs_density" / "gs_density.raw"
    gs_meta = out / "gs_density" / "gs_density.meta.txt"
    if gs_raw.exists() and gs_meta.exists():
        rho_gs = load_real_field(str(gs_raw), str(gs_meta))
        for ax, label in enumerate(("x", "y", "z")):
            plot_density_slice(rho_gs, out / f"gs_density_slice_{label}.png", axis=ax)
            print(f"  Saved gs_density_slice_{label}.png")
    else:
        print(f"  GS density not found at {gs_raw} — skipping density plots.")

    print(f"Loading observables from {csv_path}")
    df = inqview.load_observables(csv_path)
    print(f"  {len(df)} rows, columns: {list(df.columns)}")

    # ------------------------------------------------------------------
    # Time-domain plots
    # ------------------------------------------------------------------
    fig = inqview.plot_energy_vs_time(csv_path)
    fig.savefig(out / "energy_vs_time.png")
    print("  Saved energy_vs_time.png")

    fig = inqview.plot_current_vs_time(csv_path)
    fig.savefig(out / "current_vs_time.png")
    print("  Saved current_vs_time.png")

    fig = inqview.plot_dipole_vs_time(csv_path)
    fig.savefig(out / "dipole_vs_time.png")
    print("  Saved dipole_vs_time.png")

    fig = inqview.plot_observables_summary(csv_path)
    fig.savefig(out / "observables_summary.png")
    print("  Saved observables_summary.png")

    # ------------------------------------------------------------------
    # Frequency-domain spectra (Hann window, linear detrend)
    # ------------------------------------------------------------------
    ft = FourierTransform(window=WindowSpec("hann"), detrend=True)

    # Dipole spectra — reveals ionic resonance frequencies excited by the kick
    dipole_spectra = []
    for comp in ("x", "y", "z"):
        col = f"dipole_{comp}"
        if col in df.columns:
            result = ft.transform_dipole(df, component=comp)
            plot_spectrum(result, out / f"spectrum_{col}.png", x_max_au=FREQ_MAX_AU)
            print(f"  Saved spectrum_{col}.png")
            dipole_spectra.append(result)

    if dipole_spectra:
        plot_spectrum_summary(
            dipole_spectra,
            out / "spectrum_dipole_summary.png",
            x_max_au=FREQ_MAX_AU,
        )
        print("  Saved spectrum_dipole_summary.png")

    # Energy spectrum — tracks energy oscillation frequencies
    try:
        e_result = ft.transform_energy(df)
        plot_spectrum(e_result, out / "spectrum_energy.png", x_max_au=FREQ_MAX_AU)
        print("  Saved spectrum_energy.png")
    except ValueError as exc:
        print(f"  Skipping energy spectrum: {exc}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Li ionic-kick TDDFT observables and spectra.")
    parser.add_argument("--csv", type=Path, default=CSV_DEFAULT)
    args = parser.parse_args()
    main(args.csv)
