#!/usr/bin/env python3
"""Norm-corrected total energy -> post-processed stopping power for jellium WP runs.

Applies the kinetic-only WP norm correction (docs/plans/norm-corrected-stopping-power.md):

    E_corr(t) = E_total(t) - occ_WP * e_kin_ha(t) * (1 - norm_WP(t))

where INQ's reported energy_kinetic sums each orbital's kinetic / its own norm
(energy.hpp:83,55). Only the WP orbital loses norm here (CAP absorbs ~1 electron
= the WP; bath keeps norm~1), so we correct just its contribution.

No source edit, no re-run: reads the saved observables. Emits E_total vs E_corr
with the N-guards and a stopping-power estimate (window stated). The headline S
channel is chosen per the stopping-power-extraction skill; raw dE_corr/L_z is a
lower-tier estimate because a CAP-absorbed WP also removes its own KE from E_total.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/local/data/public/skcb2/tddft")
import sys
sys.path.insert(0, str(ROOT / "inq-stack" / "python"))
from inqview.visualisation.style import apply_theme  # noqa: E402
apply_theme()
HA = 27.211386


def _csv(run: Path, name: str) -> pd.DataFrame | None:
    p = run / "raw" / "observables" / name
    return pd.read_csv(p, comment="#") if p.exists() else None


def norm_correct(run_dir: str, occ_wp: float = 1.0):
    run = Path(run_dir)
    obs = _csv(run, "observables.csv")
    mom = _csv(run, "wp_momentum_stats.csv")
    rsp = _csv(run, "wp_real_space_stats.csv")
    en_file = _csv(run, "energies.csv")
    # energy_total column: observables.csv (jellium) or energies.csv (vacuum)
    if obs is not None and "energy_total" in obs:
        t = obs["time_au"].to_numpy(float); Etot = obs["energy_total"].to_numpy(float)
    elif en_file is not None:
        t = en_file["time_au"].to_numpy(float); Etot = en_file["total"].to_numpy(float)
    else:
        raise SystemExit(f"no energy_total in {run}")
    if mom is None:
        raise SystemExit(f"no wp_momentum_stats in {run} (need e_kin_ha)")

    tm = mom["time_au"].to_numpy(float)
    ekin = np.interp(t, tm, mom["e_kin_ha"].to_numpy(float))         # WP per-particle kinetic (Ha)
    # physical WP norm: prefer real-space (=∫|ψ|²dV, starts at 1); else momentum ratio
    if rsp is not None and "norm_check" in rsp:
        tr = rsp["time_au"].to_numpy(float); nr = rsp["norm_check"].to_numpy(float)
        normwp = np.interp(t, tr, nr / nr[0])
    else:
        nc = mom["norm_check"].to_numpy(float); normwp = np.interp(t, tm, nc / nc[0])

    E_corr = Etot - occ_wp * ekin * (1.0 - normwp)                   # the correction (Ha)

    # N guards
    eln = _csv(run, "electron_number.csv")
    Ntot = np.interp(t, eln["time_au"].to_numpy(float), eln["N_total"].to_numpy(float)) \
        if eln is not None else np.full_like(t, np.nan)
    N_bath = Ntot - occ_wp * normwp                                 # bath count (should be ~const)

    return dict(t=t, Etot=Etot, E_corr=E_corr, ekin=ekin, normwp=normwp,
                Ntot=Ntot, N_bath=N_bath, run=run.name)


def figure(d: dict, out: Path, L_z: float, title: str):
    t, Etot, Ecorr = d["t"], d["Etot"] * HA, d["E_corr"] * HA
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    # (1) energies, delta from t=0
    ax[0].plot(t, Etot - Etot[0], "C3", lw=1.8, label="reported ΔE_total")
    ax[0].plot(t, Ecorr - Ecorr[0], "C0", lw=1.8, label="norm-corrected ΔE_corr")
    ax[0].axhline(0, color="k", lw=.6)
    ax[0].set(xlabel="t (a.u.)", ylabel="ΔE from t=0 (eV)", title="energy: reported vs corrected")
    ax[0].legend(fontsize=8)
    # (2) WP norm + bath N guard
    ax[1].plot(t, d["normwp"], "C1", lw=1.8, label="WP norm ⟨ψ|ψ⟩")
    ax[1].plot(t, d["N_bath"] - d["N_bath"][0], "C2", lw=1.6, label="bath ΔN (guard)")
    ax[1].axhline(0, color="k", lw=.6)
    ax[1].set(xlabel="t (a.u.)", ylabel="norm / ΔN", title="absorption + bath-conservation guard")
    ax[1].legend(fontsize=8)
    # (3) corrected deposit / L_z running estimate
    S_run = (Ecorr - Ecorr[0]) / L_z
    ax[2].plot(t, S_run, "C0", lw=1.8)
    ax[2].axhline(0, color="k", lw=.6)
    ax[2].set(xlabel="t (a.u.)", ylabel=f"ΔE_corr / L_z  (eV/Bohr, L_z={L_z})",
              title="running deposit/thickness (caveat: incl. absorbed-WP KE)")
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--Lz", type=float, default=25.0, help="slab thickness (Bohr)")
    ap.add_argument("--occ-wp", type=float, default=1.0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--title", default=None)
    a = ap.parse_args(argv)
    d = norm_correct(a.run_dir, a.occ_wp)
    out = Path(a.out) if a.out else Path(a.run_dir) / "norm_corrected_stopping.png"
    figure(d, out, a.Lz, a.title or f"norm-corrected stopping — {d['run']}")

    dEtot = (d["Etot"][-1] - d["Etot"][0]) * HA
    dEcorr = (d["E_corr"][-1] - d["E_corr"][0]) * HA
    print(f"[nc] {d['run']}")
    print(f"  WP norm 1 -> {d['normwp'][-1]:.4f}   bath ΔN = {d['N_bath'][-1]-d['N_bath'][0]:+.4f} "
          f"(guard: |ΔN_bath|<~0.1 ⇒ single-orbital correction valid)")
    print(f"  ΔE_total   (reported)  = {dEtot:+.2f} eV")
    print(f"  ΔE_corr    (corrected) = {dEcorr:+.2f} eV")
    print(f"  raw S = ΔE_corr/L_z = {dEcorr/a.Lz:+.3f} eV/Bohr  (L_z={a.Lz}; "
          f"CAVEAT: includes absorbed-WP KE — use slab-face window for headline)")
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
