"""
Li BCC 2x2x2 ionic-kick TDDFT — observable analysis.

Run after the simulation to produce plots from results/observables.csv.

Usage:
    python analysis.py
    python analysis.py --csv path/to/observables.csv
"""

import argparse
from pathlib import Path

import inqview

CSV_DEFAULT = Path("results/observables.csv")


def main(csv_path: Path) -> None:
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Observables CSV not found: {csv_path}\n"
            "Run 'inq-run' first to generate the simulation output."
        )

    out = csv_path.parent
    print(f"Loading observables from {csv_path}")
    df = inqview.load_observables(csv_path)
    print(f"  {len(df)} rows, columns: {list(df.columns)}")

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Li ionic-kick TDDFT observables.")
    parser.add_argument("--csv", type=Path, default=CSV_DEFAULT)
    args = parser.parse_args()
    main(args.csv)
