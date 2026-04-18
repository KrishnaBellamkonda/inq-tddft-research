#!/usr/bin/env python3
"""
Visualise results from the free Gaussian wavepacket propagation (run_wp).

Outputs:
  results/heatmap_animation.mp4    — 2D density at z=L/2 vs time
  results/broadening_comparison.png — measured σ(t) vs analytical
  results/isosurface_animation.mp4  — 3D isosurface via PyVista

Run from 03_free_gaussian_wp_propagation/ with quantum-wave-packet pyenv active:
  source ~/.bashrc
  source /path/to/quantum-wave-packet/bin/activate
  python plot_propagation.py
"""

import os
import glob
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

RESULTS = "results"

# ─────────────────────────────────────────────────────────────────────────────
# 1. Broadening comparison
# ─────────────────────────────────────────────────────────────────────────────
def plot_broadening():
    csv_path = os.path.join(RESULTS, "width_vs_time.csv")
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    t      = data[:, 0]
    sig_x  = data[:, 1]
    sig_y  = data[:, 2]
    sig_z  = data[:, 3]
    sig_an = data[:, 4]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, sig_x,  "C0-",  lw=1.5, label=r"$\sigma_x(t)$ measured")
    ax.plot(t, sig_y,  "C1--", lw=1.5, label=r"$\sigma_y(t)$ measured")
    ax.plot(t, sig_z,  "C2:",  lw=1.5, label=r"$\sigma_z(t)$ measured")
    ax.plot(t, sig_an, "k-",   lw=2.5, label=r"$\sigma_0\sqrt{1+t^2/(4\sigma_0^4)}$ analytical")
    ax.set_xlabel("Time (a.u.)")
    ax.set_ylabel("Width (bohr)")
    ax.set_title("Free Gaussian wavepacket — quantum broadening")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    out = os.path.join(RESULTS, "broadening_comparison.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
# 2. 2D heatmap animation
# ─────────────────────────────────────────────────────────────────────────────
def plot_heatmap_animation():
    slice_files = sorted(glob.glob(os.path.join(RESULTS, "slice_t*.txt")))
    if not slice_files:
        print("  No slice files found — skipping heatmap animation.")
        return
    print(f"  Found {len(slice_files)} slice files.")

    # Read first file to get metadata and set colour scale
    def read_slice(path):
        return np.loadtxt(path, comments="#")

    first = read_slice(slice_files[0])
    N = first.shape[0]

    # Collect max for colour scale from first frame (t=0 peak is highest)
    vmax = first.max()

    # Time values from filenames  (slice_t000.txt → parse header line)
    def get_time(path):
        with open(path) as f:
            line = f.readline()
        m = re.search(r"t=([0-9.e+\-]+)", line)
        return float(m.group(1)) if m else 0.0

    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(first, origin="lower", vmin=0, vmax=vmax,
                   cmap="inferno", extent=[0, 40, 0, 40])
    plt.colorbar(im, ax=ax, label=r"$|\psi|^2$")
    ax.set_xlabel("x (bohr)")
    ax.set_ylabel("y (bohr)")
    title = ax.set_title(f"t = 0.00 a.u.")

    def update(i):
        rho = read_slice(slice_files[i])
        im.set_data(rho)
        t = get_time(slice_files[i])
        title.set_text(f"t = {t:.2f} a.u.")
        return im, title

    ani = animation.FuncAnimation(fig, update, frames=len(slice_files),
                                  interval=60, blit=False)
    out = os.path.join(RESULTS, "heatmap_animation.mp4")
    ani.save(out, writer="ffmpeg", fps=15, dpi=120)
    plt.close(fig)
    print(f"  Saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
# 3. 3D isosurface animation via PyVista
# ─────────────────────────────────────────────────────────────────────────────
def plot_isosurface_animation():
    try:
        import pyvista as pv
    except ImportError:
        print("  PyVista not available — skipping 3D animation.")
        return

    pv.start_xvfb()          # headless rendering
    pv.global_theme.background = "black"
    pv.global_theme.font.color = "white"

    d3_files = sorted(glob.glob(os.path.join(RESULTS, "density3d_t*.txt")))
    if not d3_files:
        print("  No 3D density files found — skipping isosurface animation.")
        return
    print(f"  Found {len(d3_files)} 3D density files.")

    def read_d3(path):
        with open(path) as f:
            header = f.readline()
        m_dx     = re.search(r"dx=([0-9.e+\-]+)", header)
        m_stride = re.search(r"stride=([0-9]+)", header)
        m_NC     = re.search(r"NC=([0-9]+)", header)
        m_t      = re.search(r"t=([0-9.e+\-]+)", header)
        dx_full  = float(m_dx.group(1))
        stride   = int(m_stride.group(1))
        NC       = int(m_NC.group(1))
        t_val    = float(m_t.group(1))
        data     = np.loadtxt(path, comments="#").reshape(NC, NC, NC)
        dx_c     = dx_full * stride
        return data, dx_c, NC, t_val

    png_paths = []
    pl = pv.Plotter(off_screen=True, window_size=[600, 600])

    for i, path in enumerate(d3_files):
        rho, dx_c, NC, t_val = read_d3(path)

        grid = pv.ImageData()
        grid.dimensions = (NC+1, NC+1, NC+1)   # cell-centred
        grid.spacing    = (dx_c, dx_c, dx_c)
        grid.origin     = (0, 0, 0)
        grid.cell_data["rho"] = rho.flatten(order="C")

        grid_pt = grid.cell_data_to_point_data()

        # Isosurface at 10% of maximum
        iso_val = 0.10 * rho.max()

        pl.clear()
        if iso_val > 0:
            surface = grid_pt.contour([iso_val], scalars="rho")
            pl.add_mesh(surface, color="cyan", opacity=0.6)
        pl.add_text(f"t = {t_val:.2f} a.u.", position="upper_left",
                    font_size=14)
        pl.camera_position = "iso"
        pl.camera.zoom(1.3)

        png = os.path.join(RESULTS, f"iso_frame_{i:03d}.png")
        pl.screenshot(png)
        png_paths.append(png)

    pl.close()
    print(f"  Rendered {len(png_paths)} isosurface frames.")

    # Assemble video with ffmpeg (numbered-input form; concat-list has cwd issues)
    if png_paths:
        out = os.path.join(RESULTS, "isosurface_animation.mp4")
        ret = os.system(
            f"ffmpeg -y -framerate 10 "
            f"-i {os.path.join(RESULTS, 'iso_frame_%03d.png')} "
            f"-vf 'format=yuv420p' {out} 2>/dev/null"
        )
        if ret == 0:
            print(f"  Saved {out}")
        else:
            print("  ffmpeg failed — individual PNGs are in results/")


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=== Broadening comparison ===")
    plot_broadening()

    print("\n=== 2D heatmap animation ===")
    plot_heatmap_animation()

    print("\n=== 3D isosurface animation ===")
    plot_isosurface_animation()

    print("\nAll done.")
