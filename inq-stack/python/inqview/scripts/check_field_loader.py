#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from inqview.data import load_complex_field, load_real_field


def print_real_summary(field) -> None:
    print("Loaded real field successfully")
    print(f"  field_name   : {field.meta.field_name}")
    print(f"  field_type   : {field.meta.field_type}")
    print(f"  shape        : {field.shape}")
    print(f"  dtype        : {field.dtype}")
    print(f"  origin_bohr  : {field.meta.origin_bohr}")
    print(f"  spacing_bohr : {field.meta.spacing_bohr}")
    print(f"  layout       : {field.meta.layout}")
    print(f"  units        : {field.meta.units}")
    print(f"  time_au      : {field.meta.time_au}")
    print(f"  min          : {field.min:.12e}")
    print(f"  max          : {field.max:.12e}")
    print(f"  mean         : {field.mean:.12e}")


def print_complex_summary(field) -> None:
    mag = field.magnitude
    print("Loaded complex field successfully")
    print(f"  field_name   : {field.meta.field_name}")
    print(f"  field_type   : {field.meta.field_type}")
    print(f"  shape        : {field.shape}")
    print(f"  scalar dtype : {field.real.dtype}")
    print(f"  origin_bohr  : {field.meta.origin_bohr}")
    print(f"  spacing_bohr : {field.meta.spacing_bohr}")
    print(f"  layout       : {field.meta.layout}")
    print(f"  units        : {field.meta.units}")
    print(f"  time_au      : {field.meta.time_au}")
    print(f"  |psi| min    : {float(np.min(mag)):.12e}")
    print(f"  |psi| max    : {float(np.max(mag)):.12e}")
    print(f"  |psi| mean   : {float(np.mean(mag)):.12e}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Sanity-check loading of INQ/inqkit raw field outputs."
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--real",
        type=Path,
        help="Path to a real raw file, e.g. total_density_t000000.raw",
    )
    mode.add_argument(
        "--complex-real",
        type=Path,
        help="Path to a complex real-part file, e.g. homo_t000000_real.raw",
    )

    parser.add_argument(
        "--imag",
        type=Path,
        default=None,
        help="Imaginary raw file for complex fields. Optional if the real file ends with '_real.raw'.",
    )
    parser.add_argument(
        "--meta",
        type=Path,
        default=None,
        help="Metadata sidecar path. Optional if it can be inferred from the raw filename.",
    )

    args = parser.parse_args()

    if args.real is not None:
        field = load_real_field(raw_path=args.real, meta_path=args.meta)
        print_real_summary(field)
        return

    field = load_complex_field(
        real_raw_path=args.complex_real,
        imag_raw_path=args.imag,
        meta_path=args.meta,
    )
    print_complex_summary(field)


if __name__ == "__main__":
    main()
