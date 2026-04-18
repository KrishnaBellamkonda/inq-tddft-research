from pathlib import Path
import re
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
import imageio.v2 as iio
import pyvista as pv

# ============================================================
# Configuration
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

DENSITY_DIR = PROJECT_ROOT / "results" / "laser-petubation-density"
IONS_DIR = PROJECT_ROOT / "results" / "laser-petubation-ion-positions"
VIS_DIR = PROJECT_ROOT / "results" / "laser-petubation-visualisations"
VIS_DIR.mkdir(parents=True, exist_ok=True)

FPS = 8
WINDOW_SIZE = (1280, 960)   # divisible by 16
BACKGROUND = "#0b1020"      # dark background so density is visible
TEXT_COLOR = "white"
AXIS_COLOR = "white"

# Covalent radii in angstrom, converted to bohr
ANGSTROM_TO_BOHR = 1.889726125
COVALENT_RADII_ANG = {
    "H": 0.31,
    "F": 0.57,
}
COVALENT_RADII_BOHR = {
    k: v * ANGSTROM_TO_BOHR for k, v in COVALENT_RADII_ANG.items()
}

ATOM_COLORS = {
    "H": "#4ea3ff",
    "F": "#41d17d",
}
DEFAULT_ATOM_COLOR = "#ff6b6b"

# Visual scale for nuclear spheres
NUCLEAR_RADIUS_SCALE = 1.0

# Bond / trail radii in bohr
BOND_TUBE_RADIUS = 0.07
TRAIL_TUBE_RADIUS = 0.035

# Global zoom padding
ZOOM_MARGIN_FRACTION = 0.10

# Camera
CAMERA_DIRECTION = np.array([1.45, -1.15, 0.90], dtype=float)
CAMERA_DIRECTION /= np.linalg.norm(CAMERA_DIRECTION)
CAMERA_DISTANCE_FACTOR = 1.65
CAMERA_VIEW_UP = (0.0, 0.0, 1.0)

# Density cloud rendering
# These are adaptive per frame, but these numbers define the bands.
DENSITY_Q_LOW = 0.70
DENSITY_Q_MID = 0.88
DENSITY_Q_HIGH = 0.975

# Number of sampled points in each band
NPTS_LOW = 18000
NPTS_MID = 9000
NPTS_HIGH = 3500

# Cloud appearance
CLOUD_LOW_COLOR = "#4cc9f0"
CLOUD_MID_COLOR = "#90e0ef"
CLOUD_HIGH_COLOR = "#ffffff"

CLOUD_LOW_SIZE = 4.0
CLOUD_MID_SIZE = 6.5
CLOUD_HIGH_SIZE = 10.0

CLOUD_LOW_OPACITY = 0.05
CLOUD_MID_OPACITY = 0.10
CLOUD_HIGH_OPACITY = 0.22

# ============================================================
# Headless / off-screen rendering
# ============================================================

pv.OFF_SCREEN = True
if os.name != "nt":
    try:
        pv.start_xvfb(wait=0.1)
    except Exception:
        pass

# ============================================================
# File discovery
# ============================================================

step_re = re.compile(r"(\d+)")

def extract_step(path: Path) -> int:
    m = step_re.findall(path.stem)
    if not m:
        raise ValueError(f"Could not extract step from {path}")
    return int(m[-1])

density_files = sorted(DENSITY_DIR.glob("density_step_*.dat"), key=extract_step)
ion_files = sorted(IONS_DIR.glob("positions_step_*.dat"), key=extract_step)

if not density_files:
    raise FileNotFoundError(f"No density files found in {DENSITY_DIR}")
if not ion_files:
    raise FileNotFoundError(f"No ion position files found in {IONS_DIR}")

common_steps = sorted(set(map(extract_step, density_files)) & set(map(extract_step, ion_files)))
if not common_steps:
    raise RuntimeError("No common steps between density and ion-position snapshots.")

density_map = {extract_step(p): p for p in density_files}
ion_map = {extract_step(p): p for p in ion_files}

# ============================================================
# Readers
# ============================================================

def read_positions(path: Path):
    """
    Returns:
        time: float or None
        labels: list[str]
        pos: ndarray shape (natoms, 3)
    """
    time = None
    labels = []
    rows = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#"):
                if line.lower().startswith("# time"):
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            time = float(parts[2])
                        except ValueError:
                            time = None
                continue

            parts = line.split()
            if len(parts) >= 4:
                labels.append(parts[0])
                rows.append([float(parts[1]), float(parts[2]), float(parts[3])])

    pos = np.array(rows, dtype=float)
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(f"Could not parse positions from {path}")

    return time, labels, pos


def read_density(path: Path):
    """
    Reads columns:
      ix iy iz x y z rho

    Returns:
      dict with:
        nx, ny, nz
        x_axis, y_axis, z_axis
        rho3

    Axes are canonicalised to be strictly ascending, and rho3 is
    permuted consistently.
    """
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 7:
        raise ValueError(f"Unexpected density format in {path}")

    ix = data[:, 0].astype(int)
    iy = data[:, 1].astype(int)
    iz = data[:, 2].astype(int)
    x = data[:, 3]
    y = data[:, 4]
    z = data[:, 5]
    rho = data[:, 6]

    nx = ix.max() + 1
    ny = iy.max() + 1
    nz = iz.max() + 1

    x_axis = np.zeros(nx)
    y_axis = np.zeros(ny)
    z_axis = np.zeros(nz)

    for i in range(nx):
        x_axis[i] = np.mean(x[ix == i])
    for j in range(ny):
        y_axis[j] = np.mean(y[iy == j])
    for k in range(nz):
        z_axis[k] = np.mean(z[iz == k])

    rho3 = np.empty((nx, ny, nz), dtype=float)
    for a, b, c, val in zip(ix, iy, iz, rho):
        rho3[a, b, c] = val

    def sort_axis(axis_vals, rho_arr, axis_index, axis_name):
        order = np.argsort(axis_vals)
        sorted_axis = axis_vals[order]
        rho_sorted = np.take(rho_arr, order, axis=axis_index)

        diffs = np.diff(sorted_axis)
        if np.any(diffs <= 0.0):
            raise ValueError(
                f"Axis {axis_name} in {path.name} is not strictly increasing "
                f"after sorting. Cannot build the density cloud safely."
            )
        return sorted_axis, rho_sorted

    x_axis, rho3 = sort_axis(x_axis, rho3, 0, "x")
    y_axis, rho3 = sort_axis(y_axis, rho3, 1, "y")
    z_axis, rho3 = sort_axis(z_axis, rho3, 2, "z")

    return {
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "x_axis": x_axis,
        "y_axis": y_axis,
        "z_axis": z_axis,
        "rho3": rho3,
    }

# ============================================================
# Helpers
# ============================================================

def atom_color(label: str) -> str:
    return ATOM_COLORS.get(label, DEFAULT_ATOM_COLOR)

def atom_radius_bohr(label: str) -> float:
    return NUCLEAR_RADIUS_SCALE * COVALENT_RADII_BOHR.get(label, 0.80)

def make_tube_from_points(points: np.ndarray, radius: float):
    if len(points) < 2:
        return None
    spline = pv.Spline(points, n_points=max(100, len(points) * 10))
    return spline.tube(radius=radius)

def compute_bounds_from_points(points: np.ndarray):
    pts = np.asarray(points)
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    return mins, maxs

def combine_bounds(min_a, max_a, min_b, max_b):
    mins = np.minimum(min_a, min_b)
    maxs = np.maximum(max_a, max_b)
    return mins, maxs

def expand_bounds(mins, maxs, frac=0.1, minpad=0.5):
    span = np.maximum(maxs - mins, 1e-12)
    pad = np.maximum(frac * span, minpad)
    return mins - pad, maxs + pad

def bounds_to_pyvista(mins, maxs):
    return (
        float(mins[0]), float(maxs[0]),
        float(mins[1]), float(maxs[1]),
        float(mins[2]), float(maxs[2]),
    )

def set_global_camera(plotter: pv.Plotter, mins: np.ndarray, maxs: np.ndarray):
    center = 0.5 * (mins + maxs)
    spans = np.maximum(maxs - mins, 1e-12)
    span = float(np.max(spans))
    distance = CAMERA_DISTANCE_FACTOR * span
    position = center + CAMERA_DIRECTION * distance

    plotter.camera_position = [
        tuple(position),
        tuple(center),
        CAMERA_VIEW_UP,
    ]

    near = max(0.01, 0.05 * distance)
    far = 5.0 * distance
    plotter.camera.clipping_range = (near, far)

def add_scene_text(plotter: pv.Plotter, lines):
    text = "\n".join(lines)
    plotter.add_text(
        text,
        position="upper_left",
        font_size=11,
        color=TEXT_COLOR,
        shadow=False,
    )

def should_write_output(path: Path) -> bool:
    """
    If output exists, ask whether to rewrite it.
    If user says no, skip that section.
    In non-interactive mode, existing files are skipped.
    """
    if not path.exists():
        return True

    if not sys.stdin.isatty():
        print(f"{path.name} already exists. Non-interactive session detected, so skipping it.")
        return False

    while True:
        ans = input(f"{path.name} already exists. Rewrite it? [y/N]: ").strip().lower()
        if ans in ("y", "yes"):
            return True
        if ans in ("", "n", "no"):
            print(f"Skipping {path.name}")
            return False
        print("Please answer y or n.")

def axis_half_widths(axis: np.ndarray):
    """
    Approximate half-width around each grid point for jittering cloud points.
    """
    axis = np.asarray(axis, dtype=float)
    n = len(axis)

    if n == 1:
        return np.array([0.5], dtype=float)

    mids = 0.5 * (axis[1:] + axis[:-1])
    lower = np.empty(n, dtype=float)
    upper = np.empty(n, dtype=float)

    lower[1:] = mids
    upper[:-1] = mids

    lower[0] = axis[0] - (mids[0] - axis[0])
    upper[-1] = axis[-1] + (axis[-1] - mids[-1])

    widths = 0.5 * (upper - lower)
    widths = np.maximum(widths, 1e-6)
    return widths

def sample_density_band_points(
    density_dict: dict,
    mask: np.ndarray,
    n_target: int,
    rng: np.random.Generator,
):
    """
    Sample cloud points from masked density voxels, weighted by rho.
    Points are jittered within the local voxel scale so the cloud looks fluid.
    """
    rho3 = np.clip(density_dict["rho3"], 0.0, None)
    x_axis = density_dict["x_axis"]
    y_axis = density_dict["y_axis"]
    z_axis = density_dict["z_axis"]

    idx = np.argwhere(mask)
    if idx.shape[0] == 0:
        return np.empty((0, 3), dtype=float)

    vals = rho3[mask].astype(float)
    total = vals.sum()
    if total <= 0.0:
        return np.empty((0, 3), dtype=float)

    probs = vals / total

    choice = rng.choice(idx.shape[0], size=n_target, replace=True, p=probs)
    sel = idx[choice]

    ix = sel[:, 0]
    iy = sel[:, 1]
    iz = sel[:, 2]

    pts = np.column_stack([
        x_axis[ix],
        y_axis[iy],
        z_axis[iz],
    ])

    hx = axis_half_widths(x_axis)[ix]
    hy = axis_half_widths(y_axis)[iy]
    hz = axis_half_widths(z_axis)[iz]

    jitter = np.column_stack([
        rng.uniform(-0.45, 0.45, size=n_target) * hx,
        rng.uniform(-0.45, 0.45, size=n_target) * hy,
        rng.uniform(-0.45, 0.45, size=n_target) * hz,
    ])

    pts = pts + jitter
    return pts

# ============================================================
# Load first frame and metadata
# ============================================================

first_time, atom_labels, first_pos = read_positions(ion_map[common_steps[0]])
first_density = read_density(density_map[common_steps[0]])

# ============================================================
# Load all positions / times / bond lengths
# ============================================================

all_pos = []
all_times = []

for step in common_steps:
    t, labels, pos = read_positions(ion_map[step])
    all_times.append(t)
    all_pos.append(pos)

all_pos = np.array(all_pos)  # (nframes, natoms, 3)
nframes = len(common_steps)
natoms = all_pos.shape[1]

if natoms >= 2:
    bond_lengths = np.linalg.norm(all_pos[:, 1, :] - all_pos[:, 0, :], axis=1)
else:
    bond_lengths = np.full(nframes, np.nan)

if all(t is not None for t in all_times):
    time_axis = np.array(all_times, dtype=float)
    time_label = "time (a.u.)"
else:
    time_axis = np.array(common_steps, dtype=float)
    time_label = "step"

# ============================================================
# Global zoom box
# Use ions + strong density region over all frames
# ============================================================

ion_mins, ion_maxs = compute_bounds_from_points(all_pos.reshape(-1, 3))

density_bounds_found = False
density_mins = np.array([np.inf, np.inf, np.inf], dtype=float)
density_maxs = np.array([-np.inf, -np.inf, -np.inf], dtype=float)

for step in common_steps:
    d = read_density(density_map[step])
    rho = np.clip(d["rho3"], 0.0, None)
    positive = rho[rho > 0.0]
    if positive.size == 0:
        continue

    q = np.quantile(positive, DENSITY_Q_HIGH)
    idx = np.argwhere(rho >= q)
    if idx.size == 0:
        continue

    ix = idx[:, 0]
    iy = idx[:, 1]
    iz = idx[:, 2]

    x_local = d["x_axis"]
    y_local = d["y_axis"]
    z_local = d["z_axis"]

    frame_mins = np.array([
        x_local[ix].min(),
        y_local[iy].min(),
        z_local[iz].min(),
    ], dtype=float)
    frame_maxs = np.array([
        x_local[ix].max(),
        y_local[iy].max(),
        z_local[iz].max(),
    ], dtype=float)

    if not density_bounds_found:
        density_mins = frame_mins.copy()
        density_maxs = frame_maxs.copy()
        density_bounds_found = True
    else:
        density_mins, density_maxs = combine_bounds(
            density_mins, density_maxs, frame_mins, frame_maxs
        )

if density_bounds_found:
    scene_mins, scene_maxs = combine_bounds(ion_mins, ion_maxs, density_mins, density_maxs)
else:
    scene_mins, scene_maxs = ion_mins.copy(), ion_maxs.copy()

scene_mins, scene_maxs = expand_bounds(
    scene_mins,
    scene_maxs,
    frac=ZOOM_MARGIN_FRACTION,
    minpad=0.25,
)
scene_bounds = bounds_to_pyvista(scene_mins, scene_maxs)

# ============================================================
# Bond plot
# ============================================================

def make_bond_plot(out_path: Path):
    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    ax.plot(time_axis, bond_lengths, lw=2.2)
    ax.set_xlabel(time_label)
    ax.set_ylabel("bond length (bohr)")
    ax.set_title("HF bond-length oscillation")
    ax.grid(True, alpha=0.3)

    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    print(f"Saved {out_path}")

# ============================================================
# PyVista scene builders
# ============================================================

def base_plotter():
    pl = pv.Plotter(off_screen=True, window_size=WINDOW_SIZE)
    pl.set_background(BACKGROUND)

    try:
        pl.enable_depth_peeling(number_of_peels=8)
    except Exception:
        pass

    try:
        pl.enable_anti_aliasing("ssaa")
    except Exception:
        pass

    return pl

def add_density_cloud(plotter: pv.Plotter, density_dict: dict, frame_idx: int):
    rho3 = np.clip(density_dict["rho3"], 0.0, None)
    positive = rho3[rho3 > 0.0]

    if positive.size == 0:
        return {"q_low": 0.0, "q_mid": 0.0, "q_high": 0.0, "frame_max": 0.0}

    q_low = float(np.quantile(positive, DENSITY_Q_LOW))
    q_mid = float(np.quantile(positive, DENSITY_Q_MID))
    q_high = float(np.quantile(positive, DENSITY_Q_HIGH))
    frame_max = float(positive.max())

    rng = np.random.default_rng(12345 + frame_idx)

    mask_low = rho3 >= q_low
    mask_mid = rho3 >= q_mid
    mask_high = rho3 >= q_high

    pts_low = sample_density_band_points(density_dict, mask_low, NPTS_LOW, rng)
    pts_mid = sample_density_band_points(density_dict, mask_mid, NPTS_MID, rng)
    pts_high = sample_density_band_points(density_dict, mask_high, NPTS_HIGH, rng)

    if len(pts_low) > 0:
        plotter.add_points(
            pts_low,
            color=CLOUD_LOW_COLOR,
            opacity=CLOUD_LOW_OPACITY,
            point_size=CLOUD_LOW_SIZE,
            render_points_as_spheres=True,
        )

    if len(pts_mid) > 0:
        plotter.add_points(
            pts_mid,
            color=CLOUD_MID_COLOR,
            opacity=CLOUD_MID_OPACITY,
            point_size=CLOUD_MID_SIZE,
            render_points_as_spheres=True,
        )

    if len(pts_high) > 0:
        plotter.add_points(
            pts_high,
            color=CLOUD_HIGH_COLOR,
            opacity=CLOUD_HIGH_OPACITY,
            point_size=CLOUD_HIGH_SIZE,
            render_points_as_spheres=True,
        )

    return {
        "q_low": q_low,
        "q_mid": q_mid,
        "q_high": q_high,
        "frame_max": frame_max,
    }

def add_nuclei(plotter: pv.Plotter, pos: np.ndarray):
    for i in range(natoms):
        label = atom_labels[i] if i < len(atom_labels) else f"atom{i}"
        radius = atom_radius_bohr(label)
        color = atom_color(label)

        sphere = pv.Sphere(
            radius=radius,
            center=tuple(pos[i]),
            theta_resolution=48,
            phi_resolution=48,
        )
        plotter.add_mesh(
            sphere,
            color=color,
            smooth_shading=True,
            specular=0.30,
            ambient=0.25,
        )

def add_nuclear_reference_positions(plotter: pv.Plotter):
    for i in range(natoms):
        label = atom_labels[i] if i < len(atom_labels) else f"atom{i}"
        radius = atom_radius_bohr(label)
        color = atom_color(label)

        sphere = pv.Sphere(
            radius=radius,
            center=tuple(all_pos[0, i]),
            theta_resolution=32,
            phi_resolution=32,
        )
        plotter.add_mesh(
            sphere,
            color=color,
            opacity=0.12,
            smooth_shading=True,
        )

def add_bond(plotter: pv.Plotter, pos: np.ndarray):
    if natoms < 2:
        return
    pts = np.vstack([pos[0], pos[1]])
    line = pv.Line(pts[0], pts[1], resolution=1)
    tube = line.tube(radius=BOND_TUBE_RADIUS)
    plotter.add_mesh(tube, color="white", opacity=0.95, smooth_shading=True)

def add_trails(plotter: pv.Plotter, frame_idx: int):
    if frame_idx < 1:
        return

    for i in range(natoms):
        label = atom_labels[i] if i < len(atom_labels) else f"atom{i}"
        color = atom_color(label)
        pts = all_pos[:frame_idx + 1, i, :]
        tube = make_tube_from_points(pts, radius=TRAIL_TUBE_RADIUS)
        if tube is not None:
            plotter.add_mesh(tube, color=color, opacity=0.45, smooth_shading=True)

def add_axes_box(plotter: pv.Plotter):
    plotter.show_bounds(
        bounds=scene_bounds,
        grid=None,
        location="outer",
        all_edges=True,
        color=AXIS_COLOR,
        xtitle="x (bohr)",
        ytitle="y (bohr)",
        ztitle="z (bohr)",
        font_size=10,
    )

def add_legend(plotter: pv.Plotter):
    entries = []
    seen = set()
    for label in atom_labels:
        if label in seen:
            continue
        seen.add(label)
        entries.append([label, atom_color(label)])
    if entries:
        try:
            plotter.add_legend(entries, bcolor=(0.05, 0.07, 0.14), face="circle", border=True)
        except Exception:
            pass

# ============================================================
# Frame renderers
# ============================================================

def render_ions_frame(frame_idx: int):
    step = common_steps[frame_idx]
    time, labels, pos = read_positions(ion_map[step])

    pl = base_plotter()
    add_nuclear_reference_positions(pl)
    add_trails(pl, frame_idx)
    add_bond(pl, pos)
    add_nuclei(pl, pos)
    add_axes_box(pl)
    add_legend(pl)
    set_global_camera(pl, scene_mins, scene_maxs)

    tstr = "unknown" if time is None else f"{time:.6f} a.u."
    bond = bond_lengths[frame_idx] if natoms >= 2 else np.nan

    add_scene_text(pl, [
        f"HF ion motion | step {step}",
        f"time = {tstr}",
        f"bond = {bond:.6f} bohr",
    ])

    img = pl.screenshot(return_img=True)
    pl.close()
    return img

def render_density_frame(frame_idx: int):
    step = common_steps[frame_idx]
    d = read_density(density_map[step])

    pl = base_plotter()
    info = add_density_cloud(pl, d, frame_idx)
    add_axes_box(pl)
    set_global_camera(pl, scene_mins, scene_maxs)

    add_scene_text(pl, [
        f"HF electron density | step {step}",
        f"q_low  = {info['q_low']:.3e}",
        f"q_mid  = {info['q_mid']:.3e}",
        f"q_high = {info['q_high']:.3e}",
        f"max    = {info['frame_max']:.3e}",
    ])

    img = pl.screenshot(return_img=True)
    pl.close()
    return img

def render_combined_scene_frame(frame_idx: int):
    step = common_steps[frame_idx]
    time, labels, pos = read_positions(ion_map[step])
    d = read_density(density_map[step])

    pl = base_plotter()
    add_density_cloud(pl, d, frame_idx)
    add_nuclear_reference_positions(pl)
    add_trails(pl, frame_idx)
    add_bond(pl, pos)
    add_nuclei(pl, pos)
    add_axes_box(pl)
    add_legend(pl)
    set_global_camera(pl, scene_mins, scene_maxs)

    tstr = "unknown" if time is None else f"{time:.6f} a.u."
    bond = bond_lengths[frame_idx] if natoms >= 2 else np.nan

    add_scene_text(pl, [
        f"HF ions + electron density | step {step}",
        f"time = {tstr}",
        f"bond = {bond:.6f} bohr",
    ])

    img = pl.screenshot(return_img=True)
    pl.close()
    return img

# ============================================================
# Writers
# ============================================================

def write_video(filename: Path, frame_generator, fps: int = FPS):
    with iio.get_writer(
        filename,
        fps=fps,
        codec="libx264",
        quality=8,
        macro_block_size=16,
    ) as writer:
        for frame_idx in range(nframes):
            img = frame_generator(frame_idx)
            writer.append_data(img)
            print(f"[{filename.name}] frame {frame_idx + 1}/{nframes}")

    print(f"Saved {filename}")

def write_combined_video(filename: Path, fps: int = FPS):
    fig, (ax_img, ax_plot) = plt.subplots(
        1, 2,
        figsize=(12.8, 5.76),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.9, 1.0]},
    )
    canvas = FigureCanvasAgg(fig)

    scene_img = render_combined_scene_frame(0)
    im = ax_img.imshow(scene_img)
    ax_img.set_axis_off()
    ax_img.set_title("3D density and ionic motion")

    ax_plot.plot(time_axis, bond_lengths, lw=2.0)
    ax_plot.scatter([time_axis[0]], [bond_lengths[0]], s=70, color="red", zorder=3)
    ax_plot.set_xlabel(time_label)
    ax_plot.set_ylabel("bond length (bohr)")
    ax_plot.set_title("Bond-length evolution")
    ax_plot.grid(True, alpha=0.3)

    bmin = np.nanmin(bond_lengths)
    bmax = np.nanmax(bond_lengths)
    bpad = max(0.03 * max(bmax - bmin, 1e-12), 0.02)
    ax_plot.set_ylim(bmin - bpad, bmax + bpad)

    with iio.get_writer(
        filename,
        fps=fps,
        codec="libx264",
        quality=8,
        macro_block_size=16,
    ) as writer:
        for frame_idx in range(nframes):
            scene_img = render_combined_scene_frame(frame_idx)
            im.set_data(scene_img)

            ax_plot.cla()
            ax_plot.plot(time_axis, bond_lengths, lw=2.0)
            ax_plot.scatter(
                [time_axis[frame_idx]],
                [bond_lengths[frame_idx]],
                s=70,
                color="red",
                zorder=3,
            )
            ax_plot.set_xlabel(time_label)
            ax_plot.set_ylabel("bond length (bohr)")
            ax_plot.set_title("Bond-length evolution")
            ax_plot.grid(True, alpha=0.3)
            ax_plot.set_ylim(bmin - bpad, bmax + bpad)

            canvas.draw()
            rgba = np.asarray(canvas.buffer_rgba())
            rgb = rgba[:, :, :3].copy()
            writer.append_data(rgb)

            print(f"[{filename.name}] frame {frame_idx + 1}/{nframes}")

    plt.close(fig)
    print(f"Saved {filename}")

# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Density dir:  {DENSITY_DIR}")
    print(f"Ion dir:      {IONS_DIR}")
    print(f"Output dir:   {VIS_DIR}")
    print(f"Frames:       {nframes}")
    print(f"Atom labels:  {atom_labels}")
    print("Using covalent radii (bohr):")
    for label in atom_labels:
        print(f"  {label}: {atom_radius_bohr(label):.6f}")

    print(f"scene mins = {scene_mins}")
    print(f"scene maxs = {scene_maxs}")

    bond_plot_path = VIS_DIR / "hf_bond_length.png"
    ions_video_path = VIS_DIR / "hf_ions_3d.mp4"
    density_video_path = VIS_DIR / "hf_density_3d.mp4"
    combined_video_path = VIS_DIR / "hf_combined_3d.mp4"

    if should_write_output(bond_plot_path):
        make_bond_plot(bond_plot_path)

    if should_write_output(ions_video_path):
        write_video(ions_video_path, frame_generator=render_ions_frame, fps=FPS)

    if should_write_output(density_video_path):
        write_video(density_video_path, frame_generator=render_density_frame, fps=FPS)

    if should_write_output(combined_video_path):
        write_combined_video(combined_video_path, fps=FPS)

    print("Done.")
