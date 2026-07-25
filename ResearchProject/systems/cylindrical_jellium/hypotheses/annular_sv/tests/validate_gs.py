#!/usr/bin/env python3
"""Phase-1b GS validation battery + slice rendering for the annular jellium tube.

Loads a converged GS density VTI (PHYSICAL order — via inqview.load_vti, NEVER
fftshift) and runs the numeric battery, then renders xz/yz/xy slices of BOTH the
converged electron density n and the analytic positive background n+ (the annulus
mask is reconstructed in numpy from R_in/R_out/n0 — no C++ VTI needed).

Battery (per density):
  1. Neutrality   : ∫n dV ≈ N            (<1%)
  2. Radial profile n(d), d=√(x²+y²): flat plateau ≈ n0 in the wall
     (R_in+2w<d<R_out−2w); small on-axis BORE density (d<R_in); decay outside.
  3. Cylindrical symmetry: angular std / mean at the wall mid-radius is small.
Friedel oscillations + spill-out are reported (visual, in the radial plot).

Usage:
  validate_gs.py --vti <density.vti> --rin 5 --rout 13 --lz 48 --n 24 \
                 --n0 0.001105 --rs 6 --tag rs6 --outdir <dir>
Exit 0 = battery PASS, 1 = FAIL.
"""
from __future__ import annotations
import argparse, sys, math
from pathlib import Path
import numpy as np

ROOT = Path("/local/data/public/skcb2/tddft")
sys.path.insert(0, str(ROOT / "inq-stack/python"))


def radial_profile(n, x, y, z, nbins=60, dmax=None):
    """Cylindrically + axially averaged n(d), d=√(x²+y²) about the axis (0,0)."""
    X, Y = np.meshgrid(x, y, indexing="ij")          # (nx,ny)
    D = np.sqrt(X**2 + Y**2)                          # radial distance per (x,y)
    n_xy = n.mean(axis=2)                             # average over z (uniform tube)
    if dmax is None:
        dmax = min(x.max(), y.max())
    edges = np.linspace(0, dmax, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    prof = np.full(nbins, np.nan)
    for i in range(nbins):
        m = (D >= edges[i]) & (D < edges[i + 1])
        if m.any():
            prof[i] = n_xy[m].mean()
    return centers, prof, D, n_xy


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vti", required=True)
    ap.add_argument("--rin", type=float, required=True)
    ap.add_argument("--rout", type=float, required=True)
    ap.add_argument("--lz", type=float, required=True)
    ap.add_argument("--n", type=float, required=True)
    ap.add_argument("--n0", type=float, required=True)
    ap.add_argument("--rs", type=float, required=True)
    ap.add_argument("--w", type=float, default=1.0)
    ap.add_argument("--tag", default="rs")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    from inqview import load_vti
    vf = load_vti(args.vti)                            # physical order, no fftshift
    n = np.asarray(vf.data, dtype=np.float64)
    x, y, z = np.asarray(vf.x), np.asarray(vf.y), np.asarray(vf.z)
    dV = (x[1]-x[0]) * (y[1]-y[0]) * (z[1]-z[0])

    # --- 1. neutrality -----------------------------------------------------
    integ = float(n.sum() * dV)
    neutrality_err = abs(integ - args.n) / args.n
    pass_neutral = neutrality_err < 0.01

    # --- 2. radial profile -------------------------------------------------
    d, prof, D, n_xy = radial_profile(n, x, y, z, nbins=70, dmax=min(x.max(), y.max()))
    w = args.w
    wall = (d > args.rin + 2*w) & (d < args.rout - 2*w)
    bore = d < args.rin - 2*w
    outside = d > args.rout + 2*w
    plateau = float(np.nanmean(prof[wall])) if wall.any() else float("nan")
    plateau_rel = abs(plateau - args.n0) / args.n0
    bore_frac = float(np.nanmean(prof[bore])) / args.n0 if bore.any() else float("nan")
    out_frac = float(np.nanmean(prof[outside])) / args.n0 if outside.any() else float("nan")
    pass_plateau = plateau_rel < 0.10          # wall plateau within 10% of n0
    pass_bore = bore_frac < 0.30               # bore density small (tails only)
    pass_outside = out_frac < 0.10             # decays outside the wall

    # --- 3. cylindrical symmetry (angular variation at wall mid-radius) -----
    rmid = 0.5 * (args.rin + args.rout)
    ring = (D > rmid - 1.0) & (D < rmid + 1.0)
    ang_rel = float(np.nanstd(n_xy[ring]) / np.nanmean(n_xy[ring])) if ring.any() else float("nan")
    pass_symmetry = ang_rel < 0.15             # cubic grid imprints little anisotropy

    battery = {
        "neutrality (∫n=N <1%)": (pass_neutral, f"∫n={integ:.3f} vs N={args.n:.0f} ({neutrality_err*100:.2f}%)"),
        "wall plateau ≈ n0 (<10%)": (pass_plateau, f"plateau={plateau:.3e} vs n0={args.n0:.3e} ({plateau_rel*100:.1f}%)"),
        "small bore density (<30% n0)": (pass_bore, f"bore/n0={bore_frac:.2f}"),
        "decays outside wall (<10% n0)": (pass_outside, f"out/n0={out_frac:.3f}"),
        "cylindrical symmetry (<15%)": (pass_symmetry, f"angular std/mean={ang_rel:.3f}"),
    }
    print(f"\n=== GS battery [{args.tag}] (r_s={args.rs:.2f}, N={args.n:.0f}, n0={args.n0:.3e}) ===")
    all_pass = True
    for name, (ok, detail) in battery.items():
        all_pass &= ok
        print(f"  [{'PASS' if ok else 'FAIL'}] {name:32s} {detail}")

    # --- render slices (n and analytic n+) + radial profile ----------------
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from inqview.visualisation import style as _style  # canonical theme
        try: _style.apply()
        except Exception: pass

        ix0 = np.argmin(np.abs(x)); iy0 = np.argmin(np.abs(y)); iz0 = len(z)//2
        # analytic n+ on the same grid (annulus mask, erfc edges)
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        Dxy = np.sqrt(X**2 + Y**2)
        from math import erfc
        verfc = np.vectorize(lambda t: 0.5*math.erfc(t))
        nplus = args.n0 * verfc((Dxy-args.rout)/w) * (1.0 - verfc((Dxy-args.rin)/w))

        def panel(ax, img, ttl, ext):
            im = ax.imshow(img.T, origin="lower", extent=ext, aspect="equal")
            ax.set_title(ttl); plt.colorbar(im, ax=ax, fraction=0.046)

        fig, axes = plt.subplots(2, 3, figsize=(13, 8))
        ext_xy = [x.min(), x.max(), y.min(), y.max()]
        ext_xz = [x.min(), x.max(), z.min(), z.max()]
        ext_yz = [y.min(), y.max(), z.min(), z.max()]
        panel(axes[0,0], n[:,:,iz0],   f"n  xy (z=mid)", ext_xy)
        panel(axes[0,1], n[:,iy0,:],   f"n  xz (y=0)",   ext_xz)
        panel(axes[0,2], n[ix0,:,:],   f"n  yz (x=0)",   ext_yz)
        panel(axes[1,0], nplus[:,:,iz0], f"n+ xy (z=mid)", ext_xy)
        panel(axes[1,1], nplus[:,iy0,:], f"n+ xz (y=0)",   ext_xz)
        panel(axes[1,2], nplus[ix0,:,:], f"n+ yz (x=0)",   ext_yz)
        fig.suptitle(f"Annular tube GS {args.tag}: electron density n (top) vs background n+ (bottom)")
        fig.tight_layout()
        slc_png = outdir / f"gs_slices_{args.tag}.png"
        fig.savefig(slc_png, dpi=130); plt.close(fig)
        print(f"  wrote {slc_png}")

        fig2, ax2 = plt.subplots(figsize=(7, 4.5))
        ax2.plot(d, prof, "-o", ms=3, label="n(d) electron")
        ax2.plot(d, args.n0*verfc((d-args.rout)/w)*(1-verfc((d-args.rin)/w)), "--", label="n+(d) background")
        ax2.axvline(args.rin, color="grey", ls=":"); ax2.axvline(args.rout, color="grey", ls=":")
        ax2.axhline(args.n0, color="k", lw=0.5)
        ax2.set_xlabel("d = √(x²+y²) (Bohr)"); ax2.set_ylabel("density (a0⁻³)")
        ax2.set_title(f"Radial profile {args.tag}: bore | wall | spill-out (Friedel)")
        ax2.legend(); fig2.tight_layout()
        rad_png = outdir / f"gs_radial_{args.tag}.png"
        fig2.savefig(rad_png, dpi=130); plt.close(fig2)
        print(f"  wrote {rad_png}")
    except Exception as e:
        print(f"  (plotting skipped: {e})")

    print(f"\nGS BATTERY [{args.tag}]: {'PASS' if all_pass else 'FAIL'}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
