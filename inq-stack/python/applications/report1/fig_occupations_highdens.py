"""fig_occupations_highdens — WP-GS overlap evolution for high-density jellium.

Shows |⟨ψ_WP(t)|φ_n^GS⟩|² vs time for the WP orbital projected onto
GS basis states. Identifies which GS orbitals gain/lose weight as the
WP scatters through the jellium.

Run:
    python -m applications.report1.fig_occupations_highdens
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from applications.report1._shared_style import (
    apply_style, palette_sweep5, column_widths_in, panel_label, TufteCritic,
)

RUN = Path("ResearchProject/systems/jellium/run_wp_n162_L30_E100_highdens_sigma1")
OUT = "docs/reports/report1/figures/fig_occupations_highdens.png"

AU_TO_FS = 0.02418884


def main() -> None:
    apply_style()

    overlap_dir = RUN / "results" / "raw" / "observables" / "overlap"
    index_csv = overlap_dir / "index.csv"

    if not index_csv.exists():
        print(f"ERROR: No overlap index at {index_csv}")
        return

    idx_df = pd.read_csv(index_csv)
    # Deduplicate (index.csv sometimes has duplicate rows)
    idx_df = idx_df.drop_duplicates(subset=["step", "file"])

    wp_idx = 100
    for line in (RUN / "results" / "run_summary.txt").read_text().splitlines():
        if "wp_state_index" in line:
            wp_idx = int(line.split("=")[1].strip())
            break

    # Load each overlap file: row = WP state, columns = GS state overlaps
    times = []
    overlaps = []
    for _, row in idx_df.iterrows():
        fpath = overlap_dir / row["file"]
        if not fpath.exists():
            continue
        # Read header to get time
        with open(fpath) as f:
            header = f.readline()
        t_match = header.split("time_au=")[1].split()[0] if "time_au=" in header else "0"
        t_au = float(t_match)

        # Load overlap matrix (may be single row for wp_only mode)
        data = np.loadtxt(fpath, delimiter=",", comments="#")
        if data.ndim == 1:
            # Single row: |⟨ψ_WP|φ_n^GS⟩|² for n=0..N-1
            overlap_row = data
        else:
            # Full matrix: row=evolved state, col=GS state
            # WP is the last row (state wp_idx)
            overlap_row = data[-1, :] if data.shape[0] > wp_idx else data[0, :]

        times.append(t_au)
        overlaps.append(overlap_row)

    if not times:
        print("ERROR: No overlap data loaded")
        return

    times = np.array(times)
    overlaps = np.array(overlaps)  # (n_times, n_gs_states)
    t_fs = times * AU_TO_FS
    n_gs = overlaps.shape[1]

    print(f"Loaded {len(times)} overlap snapshots, {n_gs} GS states")
    print(f"Time range: [{t_fs.min():.4f}, {t_fs.max():.4f}] fs")

    # Sort by time
    order = np.argsort(times)
    t_fs = t_fs[order]
    overlaps = overlaps[order]

    # Find states with largest overlap at t=0 and t=end
    initial = overlaps[0]
    final = overlaps[-1]
    delta = final - initial

    # Top states by final overlap (largest weight = most populated after scattering)
    top_final = np.argsort(final)[::-1][:10]
    print(f"\nTop 10 GS states by final overlap:")
    for j in top_final:
        print(f"  State {j}: initial={initial[j]:.6f}, final={final[j]:.6f}, Δ={delta[j]:.6f}")

    # Plot: top 5 by largest increase and top 5 by largest absolute overlap
    gaining = np.argsort(delta)[::-1][:5]
    top_abs = np.argsort(final)[::-1][:5]
    # Merge and deduplicate
    show_states = list(dict.fromkeys(list(gaining) + list(top_abs)))[:8]

    W = column_widths_in["single"]
    fig, ax = plt.subplots(figsize=(W, W * 0.65))

    colors = ["#881818", "#C03828", "#783898", "#2070A0", "#185070",
              "#188048", "#D05838", "#48A0D0"]

    for i, j in enumerate(show_states):
        lw = 1.2 if j == wp_idx else 0.8
        ls = "-" if delta[j] > 0 else "--"
        ax.plot(t_fs, overlaps[:, j], ls, color=colors[i % len(colors)],
                linewidth=lw, label=f"State {j}")

    ax.set_xlabel(r"Time (fs)")
    ax.set_ylabel(r"$|\langle\psi_{\mathrm{WP}}(t)|\phi_n^{\mathrm{GS}}\rangle|^2$")
    ax.set_yscale("log")
    ax.set_ylim(bottom=1e-8)
    ax.legend(fontsize=4.5, loc="best", frameon=False, ncol=2)

    critic = TufteCritic()
    for iss in critic.critique(fig):
        print(f"  TufteCritic: {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"\nSaved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
