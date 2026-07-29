#!/usr/bin/env python3
"""
Comparison ("phase") notebook builder: wide WP vs matched classical projectile, full PBC.

Usage:
  build_wp_vs_classical_pbc.py <wp_results_dir> <cl_results_dir> <out.ipynb> <out.png>

Produces:
  - <out.png>  : email figure — projectile KE vs path (initial-drag S) for WP & classical,
                 + re-ledgered ΔE_total (WP with the absorbed-orbital kinetic removed).
  - <out.ipynb>: executed study notebook (canonical theme) with the S(300 eV) comparison.

S is the INITIAL DRAG (light projectiles decelerate — .claude/rules/light-projectile-stopping.md):
  classical  S = -dKE_ion/ds over the early v>=0.85 v0 window (electron_track.csv);
  WP         S = -d(0.5 pz^2)/ds over the same window (wp_momentum_stats + wp_real_space_stats).
WP energy_total is re-ledgered by dropping the absorbed WP-orbital kinetic (e_kin_ha), which
otherwise rings as a phantom (reference_phantom_absorbed_wp_orbital_energy).
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

HA = 27.211386245988


def _read_hash(p: Path) -> pd.DataFrame:
    """Read a CSV whose first physical line is a '# ...' comment (WP stats files)."""
    with open(p) as f:
        first = f.readline()
    return pd.read_csv(p, skiprows=1 if first.lstrip().startswith("#") else 0)


def initial_drag_S(z, ke, v, v0, floors=(0.85, 0.70, 0.50)):
    """S = -dKE/ds over the early near-constant-velocity window. Returns (S_eV_per_Bohr, meanv, npts, frac)."""
    z = np.asarray(z, float); ke = np.asarray(ke, float); v = np.asarray(v, float)
    for fr in floors:
        m = v >= fr * v0
        # take the leading contiguous block from the start
        if m.sum() >= 8:
            idx = np.where(m)[0]
            # contiguous from first True
            cut = idx[0]
            end = cut
            while end + 1 < len(v) and v[end + 1] >= fr * v0:
                end += 1
            sl = slice(cut, end + 1)
            if (end - cut + 1) >= 8 and (z[sl].max() - z[sl].min()) > 1e-6:
                A = np.polyfit(z[sl], ke[sl] * HA, 1)   # ke in Ha -> eV; slope eV/Bohr
                return -A[0], float(v[sl].mean()), int(end - cut + 1), fr
    return np.nan, np.nan, 0, np.nan


def load_classical(cl: Path):
    trk = pd.read_csv(cl / "raw/observables/electron_track.csv")
    return dict(z=trk["z"].values, ke=trk["ke_ion_ha"].values, vz=trk["vz"].values,
                t=trk["time_au"].values)


def load_wp(wp: Path):
    mom = _read_hash(wp / "raw/observables/wp_momentum_stats.csv")
    rsp = _read_hash(wp / "raw/observables/wp_real_space_stats.csv")
    m = mom.merge(rsp[["step", "z_mean", "norm_check"]], on="step",
                  how="inner", suffixes=("", "_r"))
    pz = m["pz_mean"].values
    return dict(z=m["z_mean"].values, pz=pz, ke=0.5 * pz**2, t=m["time_au"].values,
                norm=m["norm_check_r"].values if "norm_check_r" in m else m["norm_check"].values)


def wp_deledgered_dE(wp: Path):
    """ΔE_total(t) [eV] with the phantom absorbed-WP-orbital kinetic removed."""
    o = pd.read_csv(wp / "raw/observables/observables.csv")
    mom = _read_hash(wp / "raw/observables/wp_momentum_stats.csv")
    m = o.merge(mom[["step", "e_kin_ha"]], on="step", how="inner")
    t = m["time_au"].values
    raw = (m["energy_total"].values - m["energy_total"].values[0]) * HA
    fix = ((m["energy_total"].values - m["e_kin_ha"].values)
           - (m["energy_total"].values[0] - m["e_kin_ha"].values[0])) * HA
    return t, raw, fix


def make_figure(wp_res: Path, cl_res: Path, out_png: Path):
    sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
    import inqview.visualisation.style as S
    S.apply_theme(); plt = S.plt
    v0 = 4.6957
    C = load_classical(cl_res); W = load_wp(wp_res)
    Sc = initial_drag_S(C["z"], C["ke"], C["vz"], v0)
    Sw = initial_drag_S(W["z"], W["pz"], W["pz"], v0)
    tE, dEraw, dEfix = wp_deledgered_dE(wp_res)

    fig, ax = plt.subplots(1, 2, figsize=(S.TWO_COL_W_IN, 3.2))
    ax[0].plot(C["z"], C["ke"] * HA, lw=1.2, color="C1", label=f"classical  S={Sc[0]:.2g} eV/Bohr")
    ax[0].plot(W["z"], W["ke"] * HA, lw=1.4, color="C0", label=f"WP  S={Sw[0]:.2g} eV/Bohr")
    ax[0].set_xlabel("projectile z (Bohr)"); ax[0].set_ylabel("projectile KE (eV)")
    ax[0].set_title("Initial-drag stopping (KE vs path)")
    ax[0].legend(fontsize=6.5, frameon=False)
    ax[1].plot(tE, dEraw, lw=1.0, color="C3", alpha=0.8, label="WP ΔE_total as logged")
    ax[1].plot(tE, dEfix, lw=1.4, color="C0", label="WP ΔE_total (orbital removed)")
    ax[1].axhline(0, ls=":", lw=0.8, color="0.6")
    ax[1].set_xlabel("time (a.u.)"); ax[1].set_ylabel(r"$\Delta E_\mathrm{total}$ (eV)")
    ax[1].set_title("WP energy ledger (PBC)")
    ax[1].legend(fontsize=6.5, frameon=False)
    S.save_presentation(fig, out_png)
    return Sc, Sw


def build_notebook(wp_res, cl_res, out_nb, out_png, Sc, Sw):
    import nbformat as nbf
    nb = nbf.v4.new_notebook()
    md = nb.cells.append
    md(nbf.v4.new_markdown_cell(
        "# Wide WP vs matched classical projectile — full PBC comparison\n\n"
        "Localised jellium slab, box 50×50×111 Bohr, dx=0.40, **periodicity 3 (full PBC)**, "
        "r_s=5.67. Matched two-sided CAP (η=−1.0 Ha, 14 Bohr/side). Projectile E=300 eV, "
        "σ=3.5 (σ_WP; classical σ_pot=σ_WP/√2=2.475). S = **initial drag** over the early "
        "v≥0.85·v0 window (light projectiles decelerate).\n\n"
        f"- **S(300 eV)_classical ≈ {Sc[0]:.2g} eV/Bohr** (mean v {Sc[1]:.2g}, {Sc[2]} pts, window {Sc[3]}·v0)\n"
        f"- **S(300 eV)_WP ≈ {Sw[0]:.2g} eV/Bohr** (mean v {Sw[1]:.2g}, {Sw[2]} pts, window {Sw[3]}·v0)\n"
        f"- **Quantum difference ΔS ≈ {Sw[0]-Sc[0]:.2g} eV/Bohr** (WP − classical)\n\n"
        f"![comparison]({out_png.name})\n\n"
        "The WP ΔE_total ‘as logged’ still rings from the phantom absorbed-orbital kinetic "
        "(a bookkeeping artifact, not physics); the ‘orbital removed’ trace is the physical "
        "energy. Both S(300 eV) sit near the numerical floor, consistent with the Lindhard "
        "high-velocity tail at v/v_F≈14. Run notebooks: `wp_pbc_E300_run_notebook.ipynb`, "
        "`classical_pbc_E300_run_notebook.ipynb`."))
    md(nbf.v4.new_markdown_cell(
        f"**Provenance** — WP: `{wp_res}` · classical: `{cl_res}`. Rebuild: "
        f"`build_wp_vs_classical_pbc.py <wp> <cl> <out.ipynb> <out.png>`."))
    import nbformat
    nbformat.write(nb, str(out_nb))


def main():
    wp_res, cl_res, out_nb, out_png = (Path(sys.argv[1]), Path(sys.argv[2]),
                                       Path(sys.argv[3]), Path(sys.argv[4]))
    Sc, Sw = make_figure(wp_res, cl_res, out_png)
    build_notebook(wp_res, cl_res, out_nb, out_png, Sc, Sw)
    print(f"S_classical={Sc[0]:.3g}  S_WP={Sw[0]:.3g} eV/Bohr ; wrote {out_nb} + {out_png}")


if __name__ == "__main__":
    main()
