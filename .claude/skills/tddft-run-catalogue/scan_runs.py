#!/usr/bin/env python3
"""Scan TDDFT run directories and upsert a runs x observables catalogue CSV.

Pure standard library (no VTK/numpy) so it runs anywhere. Detects observables
by parsing run_summary.txt for metadata and by checking for characteristic
files/directories under each run's results/ tree (VTI frame contents are NOT
walked — only the series directory names are recorded, so it stays fast).

Usage:
    python scan_runs.py --all                 # rebuild whole catalogue
    python scan_runs.py --run <run_dir>       # upsert a single run (post-run hook)
    python scan_runs.py --all --root <dir> --out <csv>

Default root: ResearchProject/systems   Default out: docs/runs_catalogue.csv
Run identity key for upsert: run_name (falls back to directory name).
"""
from __future__ import annotations
import argparse, csv, math, os, re, sys
from pathlib import Path

# ---- observable detection tables -------------------------------------------
# file-based: True if the basename appears anywhere under results/ (minus VTI frames)
FILE_OBS = {
    "observables_csv":        lambda f, d: "observables.csv" in f,
    "eigenvalues":            lambda f, d: "eigenvalues.csv" in f,
    "wp_momentum_stats":      lambda f, d: "wp_momentum_stats.csv" in f,
    "wp_realspace_stats":     lambda f, d: "wp_real_space_stats.csv" in f,
    "state_energies":         lambda f, d: "state_energies.csv" in f,
    "occupations":            lambda f, d: "occupations_vs_time.csv" in f,
    "momentum_distribution":  lambda f, d: "momentum_distribution.csv" in f,
    "gamma_transitions":      lambda f, d: "gamma_transitions.csv" in f,
    "electron_track":         lambda f, d: "electron_track.csv" in f,
    "report_md":              lambda f, d: "REPORT.md" in f,
    "loss_function":          lambda f, d: "loss_function.png" in f,
    "energy_decomposition":   lambda f, d: "energy_decomposition_classical_vs_wp.png" in f,
    "gs_basis_decomposition": lambda f, d: "gs_basis_decomposition.png" in f,
    "knudsen_ke":             lambda f, d: any(x.startswith("knudsen_ke_vs_t") for x in f),
    "kl_divergence":          lambda f, d: any(x.startswith("kl_divergence_vs_t") for x in f),
    "bath_energy":            lambda f, d: any(x.startswith("bath_energy_vs_time") for x in f),
    "stopping_curve":         lambda f, d: ("delta_E_total_vs_z.png" in f) or ("stopping_force_vs_z.png" in f),
    "overlap_heatmap":        lambda f, d: any(x.startswith("overlap_heatmap") for x in f),
    "momentum_before_after":  lambda f, d: ("wp_momentum_distribution_before_after.png" in f) or ("momentum_before_after.png" in f),
    "planewave_decomposition":lambda f, d: "planewave_decomposition.png" in f,
    "momentum_2d_map":        lambda f, d: "momentum_difference_map_2d.png" in f,
    "density_fourier":        lambda f, d: any(re.match(r"n_q_m\d", x) for x in f),
    "density_gifs":           lambda f, d: any(x.endswith(".gif") for x in f),
}
# dir-based: True if a results-relative dir path ends with any of these
DIR_OBS = {
    "density_total_vti":  ["vti/density_total", "vti/density_rt_total"],
    "density_system_vti": ["vti/density_system", "vti/density_rt_system"],
    "density_wp_vti":     ["vti/density_rt_wp", "vti/density_wp"],
    "density_delta_vti":  ["vti/density_delta", "vti/density_delta_coarse"],
    "wp_wavefunction_vti":["vti/wavefunction_wp", "vti/wavefunction_wp_rt"],
    "overlap_wp":         ["raw/overlap"],
    "overlap_full":       ["overlap_full"],
    "overlap_proxies":    ["overlap_proxies"],
    "leed_screens":       ["raw/screens", "analysis/screens"],
}
OBS_COLUMNS = list(FILE_OBS) + list(DIR_OBS)

META_COLUMNS = ["system", "run_name", "sim_type", "run_completed", "energy_ev",
                "sigma_bohr", "L_bohr", "n_electrons", "r_s", "dt_au", "n_steps",
                "total_time_au", "write_every", "norm_after", "max_overlap",
                "date_finished", "wall_time_s", "run_path"]
ALL_COLUMNS = META_COLUMNS + ["n_observables"] + OBS_COLUMNS


def parse_summary(path: Path) -> dict:
    d = {}
    try:
        for line in path.read_text(errors="replace").splitlines():
            if "=" in line and not line.strip().startswith("="):
                k, _, v = line.partition("=")
                d[k.strip()] = v.strip()
    except OSError:
        pass
    return d


def first_float(s):
    if not s:
        return None
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    return float(m.group()) if m else None


def collect_fs(results: Path):
    """Return (set of basenames, set of results-relative dir paths). Skips VTI
    frame contents (records only the series directory names) and build/."""
    files, dirs = set(), set()
    if not results.is_dir():
        return files, dirs
    for base, subdirs, fnames in os.walk(results):
        subdirs[:] = [s for s in subdirs if s != "build"]
        rel = os.path.relpath(base, results).replace(os.sep, "/")
        dirs.add(rel)
        if os.path.basename(base) == "vti":      # record series names, don't descend
            for s in subdirs:
                dirs.add(f"{rel}/{s}")
            subdirs[:] = []
            continue
        for fn in fnames:
            files.add(fn)
    return files, dirs


def derive_sim_type(run_name: str, run_type: str) -> str:
    rt = (run_type or "").lower()
    rn = run_name.lower()
    if "coronene" in rt or "leed" in rt:
        return "coronene"
    if "classical" in rt or rn.startswith("run_classical") or "_classical" in rn:
        return "classical"
    if "free" in rn or "vacuum" in rt:
        return "free_wp"
    if "wave-packet" in rt or "wave packet" in rt or rn.startswith("run_wp") or "_wp_" in rn:
        return "wp"
    return "unknown"


def parse_L(cell: str):
    """From '50^3 (cubic...)' or '35x35x60' -> (L_for_rs, label)."""
    if not cell:
        return None
    m = re.match(r"\s*(\d+(?:\.\d+)?)\s*\^\s*3", cell)
    if m:
        return float(m.group(1))
    m = re.match(r"\s*(\d+(?:\.\d+)?)\s*x\s*(\d+(?:\.\d+)?)\s*x\s*(\d+(?:\.\d+)?)", cell)
    if m:
        a, b, c = map(float, m.groups())
        return a if (a == b == c) else None  # r_s only well-defined for cubic
    return None


def scan_run(run_dir: Path, system: str) -> dict | None:
    results = run_dir / "results"
    summ = parse_summary(results / "raw" / "run_summary.txt")
    if not summ:
        summ = parse_summary(results / "run_summary.txt")
    files, dirs = collect_fs(results)
    if not summ and not files:
        return None  # not a real run dir

    run_name = summ.get("run_name") or run_dir.name
    cell = summ.get("cell_bohr", "")
    L = parse_L(cell)
    N = first_float(summ.get("n_electrons"))
    energy = first_float(summ.get("wp_energy_ev"))
    if energy is None:
        m = re.search(r"_E(\d+(?:p\d+)?)", run_name)
        energy = float(m.group(1).replace("p", ".")) if m else None
    sigma = first_float(summ.get("wp_sigma_bohr"))
    if sigma is None:
        m = re.search(r"sigma(\d+(?:p\d+)?)", run_name)
        sigma = float(m.group(1).replace("p", ".")) if m else None

    r_s = None
    if system == "jellium" and L and N:
        vol = L ** 3
        r_s = round((3.0 * vol / (4.0 * math.pi * N)) ** (1.0 / 3.0), 3)

    row = {
        "system": system,
        "run_name": run_name,
        "sim_type": derive_sim_type(run_name, summ.get("run_type", "")),
        "run_completed": summ.get("run_completed", ""),
        "energy_ev": energy if energy is not None else "",
        "sigma_bohr": sigma if sigma is not None else "",
        "L_bohr": L if L is not None else cell,
        "n_electrons": int(N) if N else "",
        "r_s": r_s if r_s is not None else "",
        "dt_au": summ.get("dt_au", ""),
        "n_steps": summ.get("rt_num_steps", ""),
        "total_time_au": summ.get("total_time_au", ""),
        "write_every": summ.get("write_every", ""),
        "norm_after": summ.get("norm_after", ""),
        "max_overlap": summ.get("max_overlap", ""),
        "date_finished": summ.get("date_finished", ""),
        "wall_time_s": summ.get("wall_time_s", ""),
        "run_path": str(run_dir),
    }
    n_obs = 0
    for name, test in FILE_OBS.items():
        present = bool(test(files, dirs))
        row[name] = int(present)
        n_obs += present
    for name, suffixes in DIR_OBS.items():
        present = any(dp.endswith(suf) for dp in dirs for suf in suffixes)
        row[name] = int(present)
        n_obs += present
    row["n_observables"] = n_obs
    return row


def find_runs(root: Path):
    for sysdir in sorted(p for p in root.iterdir() if p.is_dir()):
        system = sysdir.name
        for run_dir in sorted(sysdir.glob("run_*")):
            if run_dir.is_dir():
                yield run_dir, system


def load_csv(path: Path) -> dict:
    rows = {}
    if path.exists():
        with path.open(newline="") as fh:
            for r in csv.DictReader(fh):
                rows[r["run_name"]] = r
    return rows


def write_csv(path: Path, rows: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ALL_COLUMNS)
        w.writeheader()
        for name in sorted(rows):
            w.writerow({c: rows[name].get(c, "") for c in ALL_COLUMNS})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default="ResearchProject/systems")
    ap.add_argument("--out", default="docs/runs_catalogue.csv")
    ap.add_argument("--all", action="store_true", help="rebuild whole catalogue")
    ap.add_argument("--run", help="upsert a single run directory")
    args = ap.parse_args()
    out = Path(args.out)

    if args.run:
        run_dir = Path(args.run).resolve()
        system = run_dir.parent.name
        row = scan_run(run_dir, system)
        if not row:
            sys.exit(f"No run data found under {run_dir}")
        rows = load_csv(out)
        rows[row["run_name"]] = row
        write_csv(out, rows)
        print(f"Upserted {row['run_name']} ({row['n_observables']} observables) -> {out}")
        return

    if args.all:
        root = Path(args.root)
        rows, n = {}, 0
        for run_dir, system in find_runs(root):
            row = scan_run(run_dir, system)
            if row:
                rows[row["run_name"]] = row
                n += 1
        write_csv(out, rows)
        print(f"Catalogued {n} runs -> {out}")
        return

    ap.error("specify --all or --run <dir>")


if __name__ == "__main__":
    main()
