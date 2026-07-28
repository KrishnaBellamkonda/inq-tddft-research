#!/usr/bin/env python3
"""energy_diagnostics — CAP energy-artifact plots for a vacuum WP run, appended to
its run_report.ipynb.

Two plots per run (both mark t_enter_CAP and t_wrap = the time a classical projectile
at v=k0 reaches the +z box wall = when periodic re-entry begins):
  1. decomposed_energy.png  : every KS energy component vs t (E_total, E_kinetic,
     E_hartree, E_external, E_xc, E_non_local). In vacuum non-interacting only
     E_kinetic==E_total is non-zero (the rest are shown =0 as a bookkeeping check).
  2. compounded_energy.png  : norm(t); E_total(t); E_expected = norm*E0 (what E_total
     SHOULD be if absorption were clean); residual = E_total - norm*E0 (the ANOMALOUS
     energy stuck in the sharp CAP-boundary residual); and E_mean = E_total/norm.

Reads energies_merged.csv (per-step) + wp_real_space_stats.csv (norm_check, per
WF_EVERY) + run_summary.txt (k0, launch_z, LZ). Usage:
  energy_diagnostics.py <results_dir> [--append]     # --append -> add cells to run_report.ipynb
"""
from __future__ import annotations
import argparse, base64, re
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HA_EV = 27.211386
COMPS = ["total", "kinetic", "hartree", "external", "non_local", "xc", "exact_exchange"]


def _summary(results: Path) -> dict:
    d = {}
    for k, v in re.findall(r"([A-Za-z_]\w*)\s*=\s*(-?\d+\.?\d*)", (results/"run_summary.txt").read_text()):
        d[k.lower()] = float(v)
    return d


def _times(results: Path):
    d = _summary(results)
    k0 = d.get("k0", 5.421); launch = d.get("launch_z", -7.5); lz = d.get("cell_bohr", None)
    # cell_bohr parsed oddly (30x30x45); read LZ from the text
    m = re.search(r"cell_bohr\s*=\s*[\d.]+x[\d.]+x([\d.]+)", (results/"run_summary.txt").read_text())
    LZ = float(m.group(1)) if m else 45.0
    zc = d.get("z_cap0", LZ/2 - 15)
    t_enter = (zc - launch)/k0
    t_wrap = (LZ/2 - launch)/k0
    return t_enter, t_wrap, k0, launch, LZ


def _load(results: Path):
    obs = results/"raw/observables"
    ep = obs/"energies_merged.csv"
    if not ep.exists():
        ep = obs/"energies.csv"
    en = pd.read_csv(ep)
    rs = pd.read_csv(obs/"wp_real_space_stats.csv", comment="#")
    return en, rs


def plot_decomposed(results: Path, out: Path):
    en, _ = _load(results); t_enter, t_wrap, *_ = _times(results)
    t = en.time_au.to_numpy()
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    for c in COMPS:
        if c in en and np.any(np.abs(en[c]) > 1e-9):
            ax.plot(t, en[c]*HA_EV, lw=2.2 if c == "total" else 1.3,
                    label=f"E_{c}", zorder=5 if c == "total" else 3)
        elif c in en:
            ax.plot(t, en[c]*HA_EV, lw=0.8, ls=":", alpha=0.5, label=f"E_{c} (=0)")
    ax.axvline(t_enter, color="green", ls="--", lw=1.2, label=f"WP enters CAP (t={t_enter:.2f})")
    ax.axvline(t_wrap, color="red", ls="--", lw=1.4, label=f"wrap / re-entry (t={t_wrap:.2f})")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("energy (eV)")
    ax.set_title(f"{results.name} — decomposed KS energy vs time")
    ax.legend(fontsize=7, ncol=2); fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)


def plot_compounded(results: Path, out: Path):
    """The RESOLUTION plot. INQ reports occ·<psi|H|psi>/<psi|psi> (energy.hpp:55) =
    the intensive (per-particle MEAN) energy. Under a CAP the norm decays so this
    stays ~E0 even as the packet leaves. The physically-captured energy needs the
    EXTENSIVE energy  E_ext = E_reported * norm, which decays with the norm."""
    en, rs = _load(results); t_enter, t_wrap, *_ = _times(results)
    m = pd.merge_asof(en.sort_values("step"), rs[["step", "norm_check"]].sort_values("step"),
                      on="step", direction="nearest")
    t = m.time_au.to_numpy(); E = m.total.to_numpy()*HA_EV; N = np.clip(m.norm_check.to_numpy(), 1e-12, None)
    E0 = E[0]; E_ext = E*N; captured = E0*N[0] - E_ext
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.6, 6.4), sharex=True)
    a1.plot(t, E, "C4-", lw=1.6, label="E_reported (INQ, per-particle = ⟨H⟩/⟨ψ|ψ⟩)")
    a1.plot(t, E_ext, "k-", lw=2.2, label="E_ext = E_reported·norm  (EXTENSIVE, correct)")
    a1.plot(t, captured, "C3--", lw=1.6, label="captured = E₀ − E_ext")
    a1.set_ylabel("energy (eV)"); a1.legend(fontsize=7)
    a1.set_title(f"{results.name} — reported (per-particle) vs extensive energy\n"
                 f"captured = {captured[-1]:.1f} eV of {E0*N[0]:.1f} eV  "
                 f"(norm absorbed {100*(1-N[-1]/N[0]):.1f}%)")
    a2.semilogy(t, N, "C2-", lw=1.8, label="norm")
    a2.semilogy(t, np.clip(E_ext, 1e-3, None), "k-", lw=1.6, label="E_ext (eV)")
    a2.set_ylabel("norm  /  E_ext (eV, log)"); a2.set_xlabel("time (a.u.)"); a2.legend(fontsize=7)
    for a in (a1, a2):
        a.axvline(t_enter, color="green", ls="--", lw=1.0)
        a.axvline(t_wrap, color="red", ls="--", lw=1.2)
    a2.text(t_wrap, a2.get_ylim()[1], " wrap", color="red", fontsize=7, va="top")
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)


def plot_wrap_conservation(results: Path, out: Path):
    """Proof that periodic WRAP/re-entry is energy-conserving (no-CAP control):
    E_total(t) dead-flat across multiple wraps, with the WP peak_z sawtooth showing
    the wraps actually happen. Marks each wrap time n*t_wrap."""
    import glob, re
    from inqview import load_vti
    en, _ = _load(results); t_enter, t_wrap, k0, launch, LZ = _times(results)
    t = en.time_au.to_numpy(); E = en.total.to_numpy()*HA_EV; dE = E - E[0]
    fs = sorted(glob.glob(str(results/"raw/vti/density_wp/density_wp_t*.vti")))
    tz, pz = [], []
    for f in fs:
        st = int(re.search(r'_t(\d+)', f).group(1))
        v = load_vti(f, expect_centered_axis=None)
        tz.append(st*float(en.time_au.iloc[1]-en.time_au.iloc[0]))
        pz.append(float(v.z[np.asarray(v.data).sum(axis=(0, 1)).argmax()]))
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.6, 6.0), sharex=True)
    a1.plot(t, dE*1e3, "k-", lw=1.6)
    a1.set_ylabel("ΔE_total (meV)")
    a1.set_title(f"{results.name} — periodic wrap is energy-conserving "
                 f"(max|ΔE|={np.max(np.abs(dE))*1e3:.2f} meV over {len(en)} steps)")
    a2.plot(tz, pz, "C0.-", lw=1.0, ms=4, label="WP peak z (sawtooth = wraps)")
    a2.axhline(LZ/2, color="gray", ls=":", lw=0.8); a2.axhline(-LZ/2, color="gray", ls=":", lw=0.8)
    a2.set_ylabel("peak z (Bohr)"); a2.set_xlabel("time (a.u.)"); a2.legend(fontsize=8)
    n = 1
    while n*t_wrap <= t.max():
        for a in (a1, a2):
            a.axvline(n*t_wrap, color="red", ls="--", lw=1.1)
        a1.text(n*t_wrap, a1.get_ylim()[1], f" wrap {n}", color="red", fontsize=7, va="top")
        n += 1
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)
    print(f"[diag] wrote {out.name}  (max|dE|={np.max(np.abs(dE))*1e3:.3f} meV)")


def append_to_notebook(results: Path, pngs: list[Path]):
    import nbformat as nbf
    from nbformat.v4 import new_markdown_cell, new_code_cell, new_output
    nbp = results/"report"/"run_report.ipynb"
    nb = nbf.read(str(nbp), as_version=4)
    # drop any previously-appended diagnostics (idempotent)
    nb.cells = [c for c in nb.cells if "energy-diagnostics" not in c.get("metadata", {}).get("tag", "")]
    cells = [new_markdown_cell("## CAP energy-artifact diagnostics\n\n"
             "Green line = WP enters the CAP; **red line = wrap/re-entry** (a classical "
             "projectile at v=k₀ reaching the +z wall). The *residual* and *E_total/norm* "
             "expose energy that stays behind after the norm is absorbed.")]
    cells[0].metadata["tag"] = "energy-diagnostics"
    for p in pngs:
        md = new_markdown_cell(f"### {p.stem}"); md.metadata["tag"] = "energy-diagnostics"
        c = new_code_cell(f'from IPython.display import Image\nImage(filename="{p.name}")')
        c.metadata["tag"] = "energy-diagnostics"
        c.outputs = [new_output("display_data",
                     data={"image/png": base64.b64encode(p.read_bytes()).decode()}, metadata={})]
        cells += [md, c]
    nb.cells += cells
    nbf.write(nb, str(nbp))
    print(f"[diag] appended {len(pngs)} plots to {nbp}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir")
    ap.add_argument("--append", action="store_true")
    a = ap.parse_args(argv)
    r = Path(a.results_dir); rep = r/"report"; rep.mkdir(parents=True, exist_ok=True)
    p1, p2 = rep/"decomposed_energy.png", rep/"compounded_energy.png"
    plot_decomposed(r, p1); plot_compounded(r, p2)
    print(f"[diag] wrote {p1.name}, {p2.name}")
    if a.append:
        append_to_notebook(r, [p1, p2])


if __name__ == "__main__":
    main()
