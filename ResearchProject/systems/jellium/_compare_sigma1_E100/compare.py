"""σ=1 Bohr concentrated-WP comparison — free vs jellium at E=100 eV.

Builds the deliverable bundle the user requested 2026-05-17:

  1. Side-by-side density GIF (xz-slice; free left, jellium right).
  2. Difference density GIF (jellium minus free, frame-by-frame).
  3. Per-direction σ_p²(t), σ_r²(t), KL(t), <z>(t) plots with the IFW
     window HIGHLIGHTED (the new ifw_highlight helper).
  4. Three-way E_kin(t) — free WP, jellium WP, classical electron — with
     a clear caption that the classical kinetic-energy DROP is the
     direct signal of electronic stopping (force from bath density
     gradient under ehrenfest dynamics), not a numerical artefact.

Reads:
  ../run_free_wp_L50_E100_sigma1/results/...
  ../run_wp_n162_L50_E100_sigma1/results/...
  ../run_classical_n162_L50_E100_sigma1/results/...

Writes everything into ./ next to this script.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from inqview.pipeline._common import ifw_highlight

HA = 27.21138625
THIS = Path(__file__).parent
ROOT = THIS.parent

FREE = ROOT / "run_free_wp_L50_E100_sigma1" / "results"
JELL = ROOT / "run_wp_n162_L50_E100_sigma1"  / "results"
CLAS = ROOT / "run_classical_n162_L50_E100_sigma1" / "results"

# IFW window for σ=1: self-spread limited (see Cfg header)
T_IFW   = 9.5     # a.u. — where 3σ_density(t) reaches the far box face
T_TOTAL = 9.5     # we cap N_STEPS at the same value (no post-IFW region)


# ──────────────────────────────────────────────────────────────────────────
# VTI helpers — read frame to (cube, meta)
# ──────────────────────────────────────────────────────────────────────────
def _load_vti(path: Path):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    img = reader.GetOutput()
    nx, ny, nz = img.GetDimensions()
    flat = vtk_to_numpy(img.GetPointData().GetArray(0)).astype(np.float64)
    cube = flat.reshape((nz, ny, nx)).transpose(2, 1, 0)
    return cube, {
        "nx": nx, "ny": ny, "nz": nz,
        "origin": tuple(img.GetOrigin()),
        "spacing": tuple(img.GetSpacing()),
    }


_TRE = re.compile(r"_t(\d{6})\.vti$")
def _step_of(p: Path) -> int:
    m = _TRE.search(p.name); return int(m.group(1)) if m else -1


def _frames_for(vti_dir: Path) -> list[tuple[int, np.ndarray, dict]]:
    fs = sorted(vti_dir.glob("*_t*.vti"), key=_step_of)
    out = []
    for f in fs:
        cube, meta = _load_vti(f)
        out.append((_step_of(f), cube, meta))
    return out


# ──────────────────────────────────────────────────────────────────────────
# Density GIFs
# ──────────────────────────────────────────────────────────────────────────
def make_side_by_side_density_gif(out_gif: Path, dt_au: float) -> None:
    free_frames = _frames_for(FREE / "raw/vti/density_rt_total")
    jell_frames = _frames_for(JELL / "raw/vti/density_rt_total")
    if not (free_frames and jell_frames):
        print(f"  skipping side-by-side: free={len(free_frames)} jell={len(jell_frames)}")
        return

    n = min(len(free_frames), len(jell_frames))
    print(f"  side-by-side: {n} aligned frames")
    free_frames = free_frames[:n]
    jell_frames = jell_frames[:n]

    # Slice xz-plane (axis y, take centre)
    def xz(cube): return cube.take(cube.shape[1] // 2, axis=1)

    free_slices = [xz(c) for _, c, _ in free_frames]
    jell_slices = [xz(c) for _, c, _ in jell_frames]

    # Fixed colour scales (per-run, 99th percentile)
    free_max = float(np.percentile(np.concatenate([s.ravel() for s in free_slices]), 99))
    jell_max = float(np.percentile(np.concatenate([s.ravel() for s in jell_slices]), 99))

    meta_f = free_frames[0][2]
    meta_j = jell_frames[0][2]
    ext_f = [meta_f["origin"][0], meta_f["origin"][0] + meta_f["nx"]*meta_f["spacing"][0],
             meta_f["origin"][2], meta_f["origin"][2] + meta_f["nz"]*meta_f["spacing"][2]]
    ext_j = [meta_j["origin"][0], meta_j["origin"][0] + meta_j["nx"]*meta_j["spacing"][0],
             meta_j["origin"][2], meta_j["origin"][2] + meta_j["nz"]*meta_j["spacing"][2]]

    fig, axs = plt.subplots(1, 2, figsize=(9, 5), dpi=100)
    im_f = axs[0].imshow(free_slices[0].T, origin="lower", cmap="viridis",
                         extent=ext_f, vmin=0, vmax=free_max, aspect="equal")
    im_j = axs[1].imshow(jell_slices[0].T, origin="lower", cmap="viridis",
                         extent=ext_j, vmin=0, vmax=jell_max, aspect="equal")
    for a in axs:
        a.set_xlabel("x (Bohr)"); a.set_ylabel("z (Bohr)")
    axs[0].set_title("Free WP density")
    axs[1].set_title("Jellium WP density")
    plt.colorbar(im_f, ax=axs[0], shrink=0.8, label="n (Bohr⁻³)")
    plt.colorbar(im_j, ax=axs[1], shrink=0.8, label="n (Bohr⁻³)")
    sup = fig.suptitle(f"σ=1 Bohr, E=100 eV — t = 0.00 a.u.")
    fig.tight_layout()

    def update(i):
        im_f.set_data(free_slices[i].T)
        im_j.set_data(jell_slices[i].T)
        t = free_frames[i][0] * dt_au
        sup.set_text(f"σ=1 Bohr, E=100 eV — t = {t:.2f} a.u.")
        return im_f, im_j, sup

    anim = animation.FuncAnimation(fig, update, frames=n, interval=80, blit=False)
    anim.save(out_gif, writer="pillow", dpi=100)
    plt.close(fig)
    print(f"  wrote {out_gif}")


def make_difference_density_gif(out_gif: Path, dt_au: float) -> None:
    free_frames = _frames_for(FREE / "raw/vti/density_rt_total")
    jell_frames = _frames_for(JELL / "raw/vti/density_rt_total")
    if not (free_frames and jell_frames):
        print(f"  skipping diff: free={len(free_frames)} jell={len(jell_frames)}")
        return
    n = min(len(free_frames), len(jell_frames))
    free_frames = free_frames[:n]; jell_frames = jell_frames[:n]

    def xz(cube): return cube.take(cube.shape[1] // 2, axis=1)

    # Pre-compute the difference slices.  Use ρ_jell − ρ_jell_t0 minus
    # (ρ_free − ρ_free_t0) so we subtract each run's constant background
    # and see the difference in the *evolution* — clean comparison.
    def delta_t0(frames):
        ref = frames[0][1]
        return [c - ref for _, c, _ in frames]
    free_d = [xz(c) for c in delta_t0(free_frames)]
    jell_d = [xz(c) for c in delta_t0(jell_frames)]
    diff = [j - f for j, f in zip(jell_d, free_d)]

    vmax = float(np.percentile(np.abs(np.concatenate([d.ravel() for d in diff])), 99))
    if vmax <= 0: vmax = 1e-6

    meta = free_frames[0][2]
    ext = [meta["origin"][0], meta["origin"][0] + meta["nx"]*meta["spacing"][0],
           meta["origin"][2], meta["origin"][2] + meta["nz"]*meta["spacing"][2]]

    fig, ax = plt.subplots(figsize=(7, 5), dpi=100)
    im = ax.imshow(diff[0].T, origin="lower", cmap="seismic",
                   extent=ext, vmin=-vmax, vmax=vmax, aspect="equal")
    ax.set_xlabel("x (Bohr)"); ax.set_ylabel("z (Bohr)")
    plt.colorbar(im, ax=ax, shrink=0.85,
                 label=r"$\Delta\rho_{\rm jell}(t) - \Delta\rho_{\rm free}(t)$ (Bohr$^{-3}$)")
    title = ax.set_title("Difference density (jellium − free) — t = 0.00 a.u.")
    fig.tight_layout()

    def update(i):
        im.set_data(diff[i].T)
        t = free_frames[i][0] * dt_au
        title.set_text(f"Difference density (jellium − free) — t = {t:.2f} a.u.")
        return im, title

    anim = animation.FuncAnimation(fig, update, frames=n, interval=80, blit=False)
    anim.save(out_gif, writer="pillow", dpi=100)
    plt.close(fig)
    print(f"  wrote {out_gif}")


# ──────────────────────────────────────────────────────────────────────────
# Metric plots (IFW-highlighted)
# ──────────────────────────────────────────────────────────────────────────
def make_metric_plots() -> None:
    free_rs  = pd.read_csv(FREE / "raw/observables/wp_real_space_stats.csv", comment="#")
    free_mom = pd.read_csv(FREE / "raw/observables/wp_momentum_stats.csv", comment="#")
    jell_rs  = pd.read_csv(JELL / "raw/observables/wp_real_space_stats.csv", comment="#")
    jell_mom = pd.read_csv(JELL / "raw/observables/wp_momentum_stats.csv", comment="#")
    jell_kl  = pd.read_csv(JELL / "analysis/observables/kl_divergence.csv")
    clas_track = pd.read_csv(CLAS / "raw/observables/electron_track.csv")

    def H(ax): ifw_highlight(ax, T_IFW)

    # 1. σ_p² per direction
    fig, ax = plt.subplots(figsize=(9, 5)); H(ax)
    for col, c, lab in [("sigma_px2","C0","x"),("sigma_py2","C1","y"),("sigma_pz2","C3","z")]:
        ax.plot(free_mom["time_au"], free_mom[col], c+":", lw=1.2, label=f"free $\\sigma_{{p_{lab}}}^2$")
        ax.plot(jell_mom["time_au"], jell_mom[col], c+"-", lw=1.4, label=f"jellium $\\sigma_{{p_{lab}}}^2$")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$\sigma_p^2$ (Bohr$^{-2}$)")
    ax.set_title(r"σ=1 E=100: $\sigma_p^2(t)$ per direction — free (dotted) vs jellium (solid)")
    ax.legend(loc="best", fontsize=8); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(THIS / "sigma_p2_per_direction.png", dpi=150); plt.close(fig)

    # 2. σ_r² per direction with analytic free overlay
    fig, ax = plt.subplots(figsize=(9, 5)); H(ax)
    ax.plot(free_rs["time_au"], free_rs["sigma_z2"], "C3:", lw=1.4, label="free $\\sigma_z^2$")
    ax.plot(free_rs["time_au"], free_rs["sigma_x2"], "C0:", lw=1.2, label="free $\\sigma_x^2$")
    ax.plot(jell_rs["time_au"], jell_rs["sigma_z2"], "C3-", lw=1.6, label="jellium $\\sigma_z^2$")
    ax.plot(jell_rs["time_au"], jell_rs["sigma_x2"], "C0-", lw=1.4, label="jellium $\\sigma_x^2$")
    t = free_rs["time_au"].to_numpy()
    sigma_w = 1.0; s0 = sigma_w / math.sqrt(2.0)
    ax.plot(t, (s0**2)*(1.0 + (t/(sigma_w**2))**2), "k--", lw=1.0, label="analytic free")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$\sigma_r^2$ (Bohr$^2$)")
    ax.set_title(r"σ=1 E=100: $\sigma_r^2(t)$ per direction (free dotted, jellium solid)")
    ax.legend(loc="best", fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(THIS / "sigma_r2_per_direction.png", dpi=150); plt.close(fig)

    # 3. KL divergence (jellium only — free ≡ 0)
    fig, ax = plt.subplots(figsize=(9, 5)); H(ax)
    ax.plot(jell_kl["time_au"], jell_kl["kl_div"], "C3-", lw=1.8, label="jellium $KL(P_t \\| P_0)$")
    ax.axhline(0.0, color="C2", linestyle=":", lw=1.3, label="free $\\equiv 0$ (Gaussian preserved)")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("KL divergence (nats)")
    ax.set_title(r"σ=1 E=100: KL$(P_t \,\|\, P_0)$ — momentum-distribution drift from launch")
    ax.legend(loc="best", fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(THIS / "kl_divergence.png", dpi=150); plt.close(fig)

    # 4. <z>(t)
    fig, ax = plt.subplots(figsize=(9, 5)); H(ax)
    ax.plot(free_rs["time_au"], free_rs["z_mean"], "C2:", lw=1.6, label="free WP $\\langle z\\rangle$")
    ax.plot(jell_rs["time_au"], jell_rs["z_mean"], "C3-", lw=1.8, label="jellium WP $\\langle z\\rangle$")
    tc = clas_track["time_au"] if "time_au" in clas_track.columns else clas_track["t_au"]
    zc = clas_track["z"] if "z" in clas_track.columns else clas_track["z_bohr"]
    ax.plot(tc, zc, "C1-.", lw=1.4, label="classical electron $z(t)$")
    v_au = math.sqrt(2.0 * 100.0 / HA)
    ax.plot(t, -21.0 + v_au * t, "k--", lw=0.8, label=r"analytic $-21 + v\,t$")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$\langle z \rangle$ (Bohr)")
    ax.set_title(r"σ=1 E=100: $\langle z\rangle(t)$ — free WP, jellium WP, classical electron")
    ax.legend(loc="best", fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(THIS / "z_mean_trajectory.png", dpi=150); plt.close(fig)

    # 5. E_kin three-way with electronic-stopping annotation
    fig, ax = plt.subplots(figsize=(9, 5)); H(ax)
    ax.plot(free_mom["time_au"], free_mom["e_kin_ha"] * HA, "C2:", lw=1.5,
            label="free WP (native)")
    ax.plot(jell_mom["time_au"], jell_mom["e_kin_ha"] * HA, "C3-", lw=1.8,
            label="jellium WP (native)")
    vzc = clas_track["vz"] if "vz" in clas_track.columns else clas_track["v_z"]
    ke_clas = 0.5 * vzc**2 * HA
    ax.plot(tc, ke_clas, "C1-.", lw=1.4,
            label=r"classical $\frac{1}{2}m v_z^2$ (Ehrenfest)")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel(r"$E_{\rm kin}$ (eV)")
    ax.set_title(r"σ=1 E=100: projectile $E_{\rm kin}(t)$ — three-way")
    ax.legend(loc="best", fontsize=9); ax.grid(True, alpha=0.3)
    # Electronic-stopping note
    ax.text(0.02, 0.02,
            "Classical $E_{\\rm kin}$ drops because the\n"
            "Ehrenfest projectile feels forces from the\n"
            "self-consistent bath density gradient,\n"
            "transferring energy to electronic excitations\n"
            "of the jellium bath. THIS IS the electronic\n"
            "stopping signal — not numerics.",
            transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round", facecolor="#fff4cc", alpha=0.85),
            verticalalignment="bottom")
    fig.tight_layout(); fig.savefig(THIS / "e_kin_three_way.png", dpi=150); plt.close(fig)

    print(f"  wrote 5 metric plots to {THIS}")


def main(skip_gifs: bool = False) -> None:
    print(f"σ=1 free vs jellium @ E=100 — building comparison bundle")
    if not skip_gifs:
        print("Side-by-side density GIF...")
        make_side_by_side_density_gif(THIS / "density_side_by_side.gif", dt_au=0.020)
        print("Difference density GIF...")
        make_difference_density_gif(THIS / "density_difference.gif", dt_au=0.020)
    print("Metric plots...")
    make_metric_plots()
    print("done.")


if __name__ == "__main__":
    import sys
    main(skip_gifs="--no-gifs" in sys.argv)
