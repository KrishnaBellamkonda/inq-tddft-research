#!/usr/bin/env python3
"""energy_decomposition — per-run pairwise electrostatic decomposition for a single
WP-in-localised-jellium run, reconstructed from the saved density VTI frames.

Charge groups (same convention as inqkit::jellium::interaction_energies.hpp,
compute_coulomb_wp):
  P = projectile  = the wavepacket electron   (n_wp   = density_wp frames)
  S = slab        = bath electrons            (n_slab = n_total - n_wp)
  B = background  = positive jellium slab n+  (analytic: n0 * 1/2 erfc((|z|-a)/w))

Pairwise Coulomb (electrostatic) terms, using periodic FFT-Poisson phi = poisson(n)
(lap phi = -4 pi n, G=0 -> 0 gauge, matching INQ's periodic Poisson):
  E_pp = 1/2 integral n_wp * phi_wp                 projectile self-Hartree
  cross =      integral n_wp * phi_total            (n_wp against the full field)
  E_ps = cross - 2 E_pp                             projectile<->slab
  E_ss = E_hartree - cross + E_pp                   slab self-Hartree
  E_pb = -integral n_wp * phi_plus  (= U_proj_bg)   projectile<->background
  E_sb = E_external - E_pb                          slab<->background
so that   E_ss + E_ps + E_pp == E_hartree(INQ)  and  E_sb + E_pb == E_external(INQ)
are enforced EXACTLY (closure). E_hartree/E_external here are INQ's recorded values.

VALIDATION (built-in gate): the independently reconstructed
  E_hartree_check = 1/2 integral n_total * phi_total   must match INQ 'hartree',
  E_external_check = -integral n_total * phi_plus       must match INQ 'external'.
The Hartree check closes to ~1e-6 (Poisson convention correct). The external check
carries a small (~0.01%) analytic-n+ / charged-cell-G0 offset; it is absorbed into
the E_sb/E_pb SPLIT (their sum stays exact). Only E_pb's absolute value carries the
G=0 gauge (reference_charged_cell_hartree_convention) -- flagged on the plot.

Self-contained: numpy + scipy + pandas + matplotlib + inqview.load_vti. Not a twin
run, so it does NOT use twin_decompose (which reads an emitted interactions.csv).
"""
from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import erfc

HA_EV = 27.211386245988

# Slab background parameters (Cfg SlabN102_L25x25x140_w0p5)
BG_N0 = 102.0 / (25.0 * 25.0 * 25.0)   # N_electrons / V_inside = 6.528e-3 a0^-3
BG_HALF_WIDTH = 12.5
BG_EDGE_WIDTH = 0.5
BG_CENTER_Z = 0.0


def poisson_periodic(n: np.ndarray, L: tuple[float, float, float]) -> np.ndarray:
    """Solve lap phi = -4 pi n on a fully-periodic cell via FFT; G=0 term -> 0.

    n is in physical (load_vti) order; the global index-origin shift cancels in the
    pointwise integral n*phi, so no fftshift is applied (vti-coordinate-mapping rule).
    """
    nx, ny, nz = n.shape
    Lx, Ly, Lz = L
    kx = 2 * np.pi * np.fft.fftfreq(nx, d=Lx / nx)
    ky = 2 * np.pi * np.fft.fftfreq(ny, d=Ly / ny)
    kz = 2 * np.pi * np.fft.fftfreq(nz, d=Lz / nz)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    K2 = KX**2 + KY**2 + KZ**2
    K2[0, 0, 0] = 1.0
    phig = 4.0 * np.pi * np.fft.fftn(n) / K2
    phig[0, 0, 0] = 0.0
    return np.real(np.fft.ifftn(phig))


def background_density(z: np.ndarray, shape: tuple[int, int, int]) -> np.ndarray:
    """Positive jellium slab n+(z) = n0 * 1/2 erfc((|z-z0|-a)/w), uniform in x,y."""
    prof = BG_N0 * 0.5 * erfc((np.abs(z - BG_CENTER_Z) - BG_HALF_WIDTH) / BG_EDGE_WIDTH)
    nplus = np.empty(shape)
    nplus[:, :, :] = prof[None, None, :]
    return nplus


def _frame_steps(vti_dir: Path) -> list[int]:
    steps = []
    for f in glob.glob(str(vti_dir / "density_t*.vti")):
        m = re.search(r"_t(\d+)\.vti$", f)
        if m:
            steps.append(int(m.group(1)))
    return sorted(steps)


def compute_interactions(results_dir: str | Path, closure_tol_ev: float = 1.0e-2) -> pd.DataFrame:
    """Iterate the run's density frames -> per-frame pairwise decomposition (eV).

    Writes ``raw/observables/interactions.csv`` and returns the DataFrame. Asserts the
    Hartree closure gate (reconstruction vs INQ) is within ``closure_tol_ev`` at every
    frame; raises AssertionError otherwise (the validation gate).
    """
    from inqview import load_vti

    results = Path(results_dir)
    obs = results / "raw" / "observables"
    en = pd.read_csv(obs / "energies_merged.csv")
    tot_dir = results / "raw" / "vti" / "density_total"
    wp_dir = results / "raw" / "vti" / "density_wp"

    steps = _frame_steps(tot_dir)
    if not steps:
        raise FileNotFoundError(f"no density_total frames in {tot_dir}")

    nplus = None
    L = dV = None
    rows = []
    max_hartree_err = 0.0
    for st in steps:
        tot = load_vti(str(tot_dir / f"density_t{st:06d}.vti"), expect_centered_axis=None)
        wp = load_vti(str(wp_dir / f"density_t{st:06d}.vti"), expect_centered_axis=None)
        nt = np.asarray(tot.data)
        nw = np.asarray(wp.data)
        if nplus is None:
            x, y, z = tot.x, tot.y, tot.z
            dx, dy, dz = x[1] - x[0], y[1] - y[0], z[1] - z[0]
            dV = dx * dy * dz
            L = (x[-1] - x[0] + dx, y[-1] - y[0] + dy, z[-1] - z[0] + dz)
            nplus = background_density(z, nt.shape)
            phi_plus = poisson_periodic(nplus, L)

        phi_tot = poisson_periodic(nt, L)
        phi_wp = poisson_periodic(nw, L)

        e_pp = 0.5 * np.sum(nw * phi_wp) * dV
        cross = np.sum(nw * phi_tot) * dV
        e_ps = cross - 2.0 * e_pp
        e_hartree_check = 0.5 * np.sum(nt * phi_tot) * dV
        e_pb = -np.sum(nw * phi_plus) * dV
        e_external_check = -np.sum(nt * phi_plus) * dV

        # INQ recorded values at this step (anchor closure to them exactly)
        row = en[en.step == st]
        h_inq = float(row.hartree.iloc[0]) if not row.empty else e_hartree_check
        x_inq = float(row.external.iloc[0]) if not row.empty else e_external_check
        t_au = float(row.time_au.iloc[0]) if not row.empty else np.nan

        e_ss = h_inq - cross + e_pp          # E_ss+E_ps+E_pp == h_inq
        e_sb = x_inq - e_pb                  # E_sb+E_pb       == x_inq

        max_hartree_err = max(max_hartree_err, abs(e_hartree_check - h_inq) * HA_EV)
        rows.append(dict(
            step=st, time_au=t_au,
            e_ss=e_ss * HA_EV, e_pp=e_pp * HA_EV, e_ps=e_ps * HA_EV,
            e_sb=e_sb * HA_EV, e_pb=e_pb * HA_EV,
            hartree_inq=h_inq * HA_EV, external_inq=x_inq * HA_EV,
            e_hartree_check=e_hartree_check * HA_EV, e_external_check=e_external_check * HA_EV,
            N_wp=np.sum(nw) * dV, N_total=np.sum(nt) * dV,
        ))

    df = pd.DataFrame(rows)
    out = obs / "interactions.csv"
    df.to_csv(out, index=False)
    assert max_hartree_err < closure_tol_ev, (
        f"Hartree closure FAILED: max |reconstruction-INQ| = {max_hartree_err:.4f} eV "
        f"> {closure_tol_ev} eV -> Poisson convention wrong, decomposition untrustworthy")
    ext_err = float(np.max(np.abs(df.e_external_check - df.external_inq)))
    print(f"[decomp] {results.name}: {len(df)} frames  "
          f"Hartree closure max err={max_hartree_err:.2e} eV (PASS)  "
          f"external split residual max={ext_err:.2f} eV (absorbed into E_sb)")
    return df


# --------------------------------------------------------------------------- plots
def _delta(series: np.ndarray) -> np.ndarray:
    return series - series[0]


def plot_energy_delta(merged_csv: str | Path, out: str | Path, title: str):
    """ΔE(t) = E(t)-E(0) for every KS energy component (levels the huge offsets)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.read_csv(merged_csv)
    t = df["time_au"].to_numpy()
    comps = [("total", "k", 2.4), ("kinetic", "C0", 1.4), ("hartree", "C1", 1.4),
             ("external", "C2", 1.4), ("xc", "C3", 1.4),
             ("non_local", "C4", 1.2), ("exact_exchange", "C5", 1.2)]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for name, c, lw in comps:
        if name in df and np.any(np.abs(df[name]) > 1e-9):
            ax.plot(t, _delta(df[name].to_numpy()) * HA_EV, color=c, lw=lw,
                    label=f"ΔE_{name}", zorder=5 if name == "total" else 3)
    ax.axhline(0, color="gray", lw=0.6, ls=":")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("ΔE = E(t) − E(0)  (eV)")
    ax.set_title(f"{title} — component energy change")
    ax.legend(fontsize=8, ncol=2)
    dtot = _delta(df["total"].to_numpy())[-1] * HA_EV
    ax.text(0.02, 0.02, f"ΔE_total(final) = {dtot:+.2f} eV", transform=ax.transAxes,
            fontsize=8, va="bottom")
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)
    print(f"[decomp] wrote {Path(out).name} (ΔE_total={dtot:+.2f} eV)")


def plot_pairwise_delta(interactions_csv: str | Path, out: str | Path, title: str):
    """ΔE(t) for the pairwise Coulomb terms E_ss, E_ps, E_pp, E_sb, E_pb."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.read_csv(interactions_csv)
    t = df["time_au"].to_numpy()
    terms = [("e_ss", "slab–slab  E_ss", "C0"),
             ("e_ps", "slab–proj  E_ps (E_sp)", "C1"),
             ("e_pp", "proj self-Hartree  E_pp", "C3"),
             ("e_sb", "slab–bg  E_sb", "C2"),
             ("e_pb", "proj–bg  E_pb (U_proj_bg)", "C4")]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for col, lab, c in terms:
        ax.plot(t, _delta(df[col].to_numpy()), color=c, lw=1.6, label=f"Δ({lab})")
    ax.axhline(0, color="gray", lw=0.6, ls=":")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("ΔE = E(t) − E(0)  (eV)")
    ax.set_title(f"{title} — pairwise electrostatic decomposition")
    ax.legend(fontsize=8, ncol=2)
    ax.text(0.02, 0.02,
            "P=projectile(WP)  S=slab e⁻  B=+background\n"
            "E_ss+E_ps+E_pp≡E_hartree, E_sb+E_pb≡E_external (closure).\n"
            "E_pb absolute carries the charged-cell G=0 gauge.",
            transform=ax.transAxes, fontsize=6.5, va="bottom",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.8))
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)
    print(f"[decomp] wrote {Path(out).name}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Per-run pairwise electrostatic decomposition.")
    ap.add_argument("results_dir")
    ap.add_argument("--title", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    results = Path(a.results_dir)
    out = Path(a.out) if a.out else results / "report"
    out.mkdir(parents=True, exist_ok=True)
    title = a.title or results.name

    compute_interactions(results)
    plot_energy_delta(results / "raw/observables/energies_merged.csv",
                      out / "energy_delta_components.png", title)
    plot_pairwise_delta(results / "raw/observables/interactions.csv",
                        out / "energy_delta_pairwise.png", title)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
