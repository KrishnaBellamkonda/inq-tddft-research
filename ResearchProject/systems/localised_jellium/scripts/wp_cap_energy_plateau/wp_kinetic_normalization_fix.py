#!/usr/bin/env python3
"""wp_kinetic_normalization_fix — correct the localised-jellium CAP plateau for INQ's
per-particle (norm-normalized) energy reporting.

INQ's energy (inq/src/hamiltonian/energy.hpp:50-55 occ_sum) reports each orbital's
kinetic as occ*<psi|T|psi>/<psi|psi>. Under the CAP the WP orbital's norm decays, so
its kinetic contribution stays at the per-particle MEAN instead of leaving the ledger
(extensive). Hartree/external/xc are density-based (extensive, already correct). So the
ONLY correction is the WP kinetic:
    E_total_corrected(t) = E_total_reported(t) - occ_WP * <T_WP>(t) * (1/norm_WP(t) - 1)
with occ_WP=1 (inject_into_last_extra_state occ=1). <T_WP> and norm_WP are computed
from the saved complex wavefunction_wp frames (validated at t0: mean KE = 120.4 eV).

Outputs corrected_plateau.png (reported vs corrected E_total, nocap vs cap) into the
comparison dir and prints the reported vs corrected gap. See
[[reference_inq_reports_normalized_energy]].
"""
from __future__ import annotations
import glob, re
from pathlib import Path
import numpy as np, pandas as pd, vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HA = 27.211386; H = 0.5; DV = H**3; OCC = 1.0
BASE = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
            "scripts/wp_cap_energy_plateau/wp/results")


def _psi(path):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(str(path)); r.Update(); img = r.GetOutput()
    nx, ny, nz = img.GetDimensions(); pd_ = img.GetPointData()
    re_ = vtk_to_numpy(pd_.GetArray("wavefunction_real")).reshape(nz, ny, nx)
    im_ = vtk_to_numpy(pd_.GetArray("wavefunction_imag")).reshape(nz, ny, nx)
    return re_ + 1j*im_


def _kin_norm(psi):
    n = psi.size
    norm = float(np.sum(np.abs(psi)**2)*DV)
    kx = 2*np.pi*np.fft.fftfreq(psi.shape[0], d=H)
    ky = 2*np.pi*np.fft.fftfreq(psi.shape[1], d=H)
    kz = 2*np.pi*np.fft.fftfreq(psi.shape[2], d=H)
    K2 = (kx[:, None, None]**2 + ky[None, :, None]**2 + kz[None, None, :]**2)
    pk = np.fft.fftn(psi)
    T = float(0.5*DV*np.sum(K2*np.abs(pk)**2)/n)
    return T, norm


def wp_kinetic_series(tag: str) -> pd.DataFrame:
    fs = sorted(glob.glob(str(BASE/tag/"raw/vti/wavefunction_wp/wavefunction_t*.vti")))
    rows = []
    for f in fs:
        st = int(re.search(r'_t(\d+)', f).group(1))
        T, N = _kin_norm(_psi(f))
        rows.append(dict(step=st, T_wp_ha=T, norm_wp=N))
    return pd.DataFrame(rows)


def corrected(tag: str) -> pd.DataFrame:
    en = pd.read_csv(BASE/tag/"raw/observables/energies_merged.csv")
    wp = wp_kinetic_series(tag)
    m = pd.merge_asof(wp.sort_values("step"), en[["step", "time_au", "total"]].sort_values("step"),
                      on="step", direction="nearest")
    m["E_reported_eV"] = m.total*HA
    # correction = occ*<T>*(1/norm - 1) = occ*meanKE*(1-norm); subtract it
    m["correction_eV"] = OCC*m.T_wp_ha*(1.0/m.norm_wp - 1.0)*HA
    m["E_corrected_eV"] = m.E_reported_eV - m.correction_eV
    return m


def main():
    out = BASE/"comparison"; out.mkdir(parents=True, exist_ok=True)
    print("computing WP kinetic series (501 frames x 2 runs)...")
    nc = corrected("nocap"); cp = corrected("cap")
    g_rep = nc.E_reported_eV.iloc[-1] - cp.E_reported_eV.iloc[-1]
    g_cor = nc.E_corrected_eV.iloc[-1] - cp.E_corrected_eV.iloc[-1]
    print(f"nocap WP norm(tF)={nc.norm_wp.iloc[-1]:.4f}  cap WP norm(tF)={cp.norm_wp.iloc[-1]:.4f}")
    print(f"cap correction(tF) = {cp.correction_eV.iloc[-1]:.1f} eV  (meanKE*(1-norm))")
    print(f"REPORTED plateau gap  nocap-cap = {g_rep:.1f} eV")
    print(f"CORRECTED plateau gap nocap-cap = {g_cor:.1f} eV")

    fig, ax = plt.subplots(figsize=(7.8, 4.6))
    for d, lab, c in ((nc, "nocap", "C0"), (cp, "cap", "C3")):
        ax.plot(d.time_au, d.E_reported_eV, c+"--", lw=1.3, label=f"{lab} reported (INQ)")
        ax.plot(d.time_au, d.E_corrected_eV, c+"-", lw=2.0, label=f"{lab} corrected (WP kin. extensive)")
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("E_total (eV)")
    ax.set_title(f"Jellium plateau — WP-kinetic normalization correction\n"
                 f"gap nocap−cap: reported {g_rep:.0f} eV  →  corrected {g_cor:.0f} eV")
    ax.legend(fontsize=8)
    ax.text(0.02, 0.03, "corrected = reported − ⟨T_WP⟩·(1/norm−1); only the CAP run shifts "
            "(nocap WP norm≈1)", transform=ax.transAxes, fontsize=7,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.8))
    p = out/"corrected_plateau.png"; fig.tight_layout(); fig.savefig(p, dpi=140); plt.close(fig)
    print("wrote", p)
    nc.to_csv(out/"nocap_wp_kinetic.csv", index=False); cp.to_csv(out/"cap_wp_kinetic.csv", index=False)
    return nc, cp


if __name__ == "__main__":
    main()
