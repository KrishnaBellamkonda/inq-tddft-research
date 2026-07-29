#!/usr/bin/env python3
"""Phase-3 analysis — big-box (50×50×90) WP + classical, two-sided CAP, τ=100 a.u.

Adapts the P2.1 analysis to the 90-box geometry and ADDS the new energy-method
diagnostics requested for the production run:
  - `energy_method()`  : retained-energy stopping S=[E_total(t_f)−E_GS]/L_z, the
    convergence TRIPLE (residual WP norm / late |dE/dt| / E_total plateau width),
    the t=0 cross-term sanity (E_total(0)−⟨T_WP⟩−SIE vs E_GS), and the
    CAP-absorbed-energy LEDGER (WP-carried vs bath/collective-carried, diagnostic).
Everything else (classical KE(z), centroid/spread, reflection, momentum, energetics,
excitation, norm/absorption, heuristics, density-GIF battery) carries over with the
90-box constants.

Run (production = default; smoke = P3_TAG=p3v for testing while production runs):
  P3_TAG=p3 PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
  /local/data/public/skcb2/tddft/venv/bin/python3 analyse_phase3.py
"""
import os, glob, re, json, traceback
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.visualisation import style, make_density_gif_battery
from inqview.analysis import compute_heuristics
style.apply_theme()

HA = 27.211386
ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
TAG = os.environ.get("P3_TAG", "p3")             # "p3" = production, "p3v" = smoke
WP_RUN = f"{ROOT}/scripts/qsp_phase3/wp/results/{TAG}_wp"
CL_RUN = f"{ROOT}/scripts/qsp_phase3/classical/results/{TAG}_classical"
WP = f"{WP_RUN}/raw"
CL = f"{CL_RUN}/raw"
HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = f"{HERE}/figs"; os.makedirs(FIGS, exist_ok=True)

DT = 0.04                                  # production timestep
E_GS = -70.22568216820937                  # 90-box GS total energy (Ha)
SLAB_KIN_GS = 2.778065777432               # 90-box GS kinetic (Ha) — for the WP-kinetic split
SIE_EV = 4.40                              # Phase-1 SIE floor (eV) for σ=0.5
SLAB_FACE, CAP_IN, BOX = 12.5, 35.0, 45.0  # slab face, CAP inner face, box edge (|z|)
L_SLAB = 25.0                              # traversal length for dE/dx
SIGMA_WP = 0.5
Z0, V0, RS = -23.75, 2.711, 5.666          # equidistant launch, projectile speed, density
CAP_ETA, CAP_W_BOHR = -0.7, 10.0
TAU_AU = 100.0
T_EXIT = (SLAB_FACE - Z0) / V0             # reach far slab face +12.5  ≈ 13.37 a.u.
T_ENTER = (-SLAB_FACE - Z0) / V0           # reach near slab face −12.5 ≈  4.15 a.u.

def mark_exit(ax):
    ax.axvline(T_ENTER, ls=":", lw=.9, color="0.55")
    ax.axvline(T_EXIT, ls="--", lw=1.1, color="C2")
    yl = ax.get_ylim()
    ax.text(T_EXIT, yl[1], f" slab end\n t={T_EXIT:.1f}", fontsize=6.5, color="C2", va="top", ha="left")

R = {"notes": [], "TAG": TAG}
def log(m): print(m); R["notes"].append(m)
def guard(name):
    def deco(fn):
        def wrap():
            try: fn()
            except Exception as e:
                log(f"[SKIP {name}] {e}"); traceback.print_exc()
        return wrap
    return deco

# ---- Lindhard point reference ----
try:
    from inqview.analysis import lindhard_elf as L
    kF = L.kF_from_rs(RS)
    R["S_lindhard_point_100eV_eVbohr"] = HA * L.stopping_power_point(V0, kF)
    log(f"point-Lindhard S(100eV, r_s5.67) = {R['S_lindhard_point_100eV_eVbohr']:.3f} eV/Bohr")
except Exception as e:
    log(f"[lindhard skip] {e}"); R["S_lindhard_point_100eV_eVbohr"] = None

# ---- complex-field helpers ----
import vtk
from vtk.util.numpy_support import vtk_to_numpy
def ftime(p): return int(re.search(r"_t(\d+)\.vti", p).group(1)) * DT
def load_psi(path):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(path); r.Update(); img = r.GetOutput()
    nx, ny, nz = img.GetDimensions(); pd_ = img.GetPointData()
    re_ = vtk_to_numpy(pd_.GetArray("wavefunction_real")).reshape((nz, ny, nx)).transpose(2, 1, 0)
    im_ = vtk_to_numpy(pd_.GetArray("wavefunction_imag")).reshape((nz, ny, nx)).transpose(2, 1, 0)
    oz = img.GetOrigin()[2]; sx, sy, sz = img.GetSpacing()
    z = oz + (np.arange(nz) + 0.5) * sz
    return (re_ + 1j * im_), z, (sx, sy, sz)
def dpsi_dz(psi, sz):
    nz = psi.shape[2]; kz = 2*np.pi*np.fft.fftfreq(nz, d=sz)
    return np.fft.ifft(1j*kz[None, None, :]*np.fft.fft(psi, axis=2), axis=2)

# ============================ 1. CONVERGENCE (WP) ============================
@guard("convergence")
def convergence():
    o = pd.read_csv(f"{WP}/observables/observables.csv")
    t, E = o["time_au"].values, o["energy_total"].values
    n = pd.read_csv(f"{WP}/observables/electron_number.csv")
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 6.4), sharex=True)
    a1.plot(t, (E - E[0]) * HA, "C0-", lw=1.6); a1.axhline(0, color="0.6", lw=.6)
    a1.set_ylabel("E_total(t) − E_total(0)  (eV)")
    a1.set_title("WP run — total-energy convergence", fontsize=9); a1.grid(alpha=.25)
    a2.plot(n["time_au"], n["N_total"], "C1-", lw=1.6)
    a2.set_xlabel("time (a.u.)"); a2.set_ylabel("N_total"); a2.grid(alpha=.25)
    m = t >= 0.85*t.max(); slope = np.polyfit(t[m], E[m], 1)[0]*HA if m.sum() > 1 else float("nan")
    a2.set_title(f"N_total (late dE/dt = {slope:+.2f} eV/a.u.; ~0 ⇒ converged)", fontsize=9)
    mark_exit(a1); mark_exit(a2)
    fig.tight_layout(); fig.savefig(f"{FIGS}/conv_wp.png", dpi=180); plt.close(fig)
    R["wp_dEtot_eV"] = float((E[-1]-E[0])*HA); R["wp_late_slope_eV_au"] = float(slope)
    R["wp_Ntot_final"] = float(n["N_total"].values[-1]); R["t_final"] = float(t[-1])
    log(f"WP: ΔE_total={R['wp_dEtot_eV']:.1f} eV, late slope={slope:+.2f} eV/au, N_final={R['wp_Ntot_final']:.3f}")

# ===== 2. ENERGY METHOD (NEW): S + convergence triple + cross-term + CAP ledger =====
@guard("energy_method")
def energy_method():
    o = pd.read_csv(f"{WP}/observables/observables.csv")
    t, E, Ekin = o["time_au"].values, o["energy_total"].values, o["energy_kinetic"].values
    s = pd.read_csv(f"{WP}/observables/wp_real_space_stats.csv", comment="#")
    tn, norm = s["time_au"].values, s["norm_check"].values

    # --- t=0 cross-term sanity: E_total(0) − ⟨T_WP⟩ − SIE  vs  E_GS ---
    T_WP_analytic = 0.5*V0**2 + 3.0/(4*SIGMA_WP**2)   # drift ½k₀² + zero-point 3/(4σ²) = 6.675 Ha
    T_WP_run = float(Ekin[0] - SLAB_KIN_GS)           # run check: total kinetic(0) − slab GS kinetic
    E0 = float(E[0])
    cross_plus_sie_eV = (E0 - T_WP_analytic - E_GS) * HA     # residual after removing kinetic
    E_system0_minus_EGS_eV = (E0 - T_WP_analytic - SIE_EV/HA - E_GS) * HA  # after also removing SIE
    R["E_total_0_Ha"] = E0
    R["T_WP_analytic_Ha"] = float(T_WP_analytic); R["T_WP_run_Ha"] = T_WP_run
    R["t0_cross_plus_sie_eV"] = float(cross_plus_sie_eV)
    R["E_system0_minus_EGS_eV"] = float(E_system0_minus_EGS_eV)

    # --- convergence triple ---
    Ef = float(E[-1]); deposited = (Ef - E_GS) * HA
    m = t >= 0.85*t.max(); slope = float(np.polyfit(t[m], E[m], 1)[0]*HA) if m.sum() > 1 else float("nan")
    norm_f = float(norm[-1])
    tol_Ha = 1.0 / HA                                  # plateau = within 1 eV of E_total(t_f)
    plateau_start = float(t[-1])
    for i in range(len(E)-1, -1, -1):
        if abs(E[i]-Ef) < tol_Ha: plateau_start = float(t[i])
        else: break
    plateau_w = float(t[-1] - plateau_start)
    converged = bool(norm_f < 0.02 and abs(slope) < 0.2)
    R["wp_deposited_EminusEGS_eV"] = float(deposited)
    R["wp_S_eVbohr"] = float(deposited / L_SLAB)
    R["wp_norm_final"] = norm_f; R["wp_plateau_width_au"] = plateau_w
    R["wp_converged"] = converged; R["wp_late_slope_eV_au"] = slope

    # --- CAP-absorbed-energy ledger (cumulative) ---
    E_removed = (E0 - E) * HA                          # total energy the CAP has removed (eV)
    norm_on_t = np.interp(t, tn, norm)
    wp_absorbed = norm_on_t[0] - norm_on_t             # WP norm absorbed so far (fraction)
    wp_carried = wp_absorbed * T_WP_analytic * HA      # rough WP-carried energy (eV) — see caveat
    bath_carried = E_removed - wp_carried              # remainder = bath/collective + secondaries
    R["cap_removed_total_eV"] = float(E_removed[-1])
    R["cap_wp_carried_eV"] = float(wp_carried[-1])
    R["cap_bath_carried_eV"] = float(bath_carried[-1])

    # --- 3-panel figure ---
    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(7.6, 9.2), sharex=True)
    a1.plot(t, (E - E_GS) * HA, "C0-", lw=1.7)
    a1.axhline(deposited, ls=":", color="0.5",
               label=f"E_total(t_f)−E_GS = {deposited:.1f} eV  ⇒  S={deposited/L_SLAB:.2f} eV/Bohr")
    a1.set_ylabel("E_total(t) − E_GS  (eV)")
    a1.set_title(f"Retained 'system' energy (= S·L_z when converged); converged={converged}", fontsize=9)
    a1.legend(fontsize=7, frameon=False); a1.grid(alpha=.25); mark_exit(a1)

    a2.plot(tn, norm, "C1-", lw=1.6, label="WP orbital norm")
    a2.axhline(0.02, ls="--", color="0.6", lw=.8, label="convergence gate 0.02")
    a2.set_ylabel("WP orbital norm")
    a2.set_title(f"Convergence triple: norm_f={norm_f:.3f}, late dE/dt={slope:+.2f} eV/au, "
                 f"plateau={plateau_w:.1f} a.u.", fontsize=8.5)
    a2.legend(fontsize=7, frameon=False); a2.grid(alpha=.25); mark_exit(a2)

    a3.plot(t, E_removed, "C3-", lw=1.7, label=f"total CAP-removed ({E_removed[-1]:.0f} eV)")
    a3.plot(t, wp_carried, "C0--", lw=1.3, label="WP-carried ≈ norm_abs·⟨T_WP⟩ (approx)")
    a3.plot(t, bath_carried, "C2-.", lw=1.3, label="bath/collective + secondaries (diagnostic)")
    a3.axhline(0, color="0.6", lw=.5)
    a3.set_xlabel("time (a.u.)"); a3.set_ylabel("cumulative energy removed (eV)")
    a3.set_title("CAP-absorbed-energy ledger (WP vs bath split — bath part is the leakage diagnostic)", fontsize=8)
    a3.legend(fontsize=7, frameon=False); a3.grid(alpha=.25); mark_exit(a3)
    fig.tight_layout(); fig.savefig(f"{FIGS}/energy_method.png", dpi=180); plt.close(fig)

    log(f"energy-method: deposited(E_total(t_f)−E_GS)={deposited:.1f} eV ⇒ S={deposited/L_SLAB:.2f} eV/Bohr "
        f"[{'CONVERGED' if converged else 'NOT converged — UPPER/LOWER bound'}]")
    log(f"  convergence triple: norm_f={norm_f:.3f} (gate<0.02), late dE/dt={slope:+.2f} eV/au, "
        f"plateau={plateau_w:.1f} a.u.")
    log(f"  t=0 sanity: ⟨T_WP⟩ analytic={T_WP_analytic*HA:.1f} eV (run {T_WP_run*HA:.1f}); "
        f"E_total(0)−⟨T_WP⟩−E_GS={cross_plus_sie_eV:.2f} eV (cross+SIE); "
        f"E_system(0)−E_GS={E_system0_minus_EGS_eV:+.2f} eV (≈0 ⇒ bookkeeping OK)")
    log(f"  CAP ledger: total removed={E_removed[-1]:.0f} eV, WP-carried≈{wp_carried[-1]:.0f}, "
        f"bath/collective≈{bath_carried[-1]:.0f} eV (leakage diagnostic; approximate split)")

# ===================== 3. CLASSICAL KE + STOPPING ==========================
@guard("classical")
def classical():
    trk = pd.read_csv(f"{CL}/observables/electron_track.csv").drop_duplicates("step")
    t, z, ke = trk["time_au"].values, trk["z"].values, trk["ke_ion_ha"].values
    fig, ax = plt.subplots(figsize=(7, 4.3))
    ax.plot(t, ke*HA, "C3-", lw=1.6); ax.set_xlabel("time (a.u.)"); ax.set_ylabel("projectile KE (eV)")
    ax.set_title("Classical projectile KE — steady state after slab exit?", fontsize=9)
    mark_exit(ax); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(f"{FIGS}/classical_ke.png", dpi=180); plt.close(fig)
    try:
        zin = np.interp(-SLAB_FACE, z, t); zout = np.interp(SLAB_FACE, z, t)
        ke_in = np.interp(zin, t, ke); ke_out = np.interp(zout, t, ke)
        S = (ke_in - ke_out)*HA / L_SLAB
        R["classical_S_eVbohr"] = float(S)
        log(f"classical: KE {ke_in*HA:.1f}→{ke_out*HA:.1f} eV over {L_SLAB} Bohr ⇒ S={S:.3f} eV/Bohr")
    except Exception as e:
        log(f"[classical stopping skip] {e}")
    R["classical_ke0_eV"] = float(ke[0]*HA); R["classical_kef_eV"] = float(ke[-1]*HA)

# ============ 3b. CLASSICAL TRANSPORT: z(t) + KE(z) dip-and-recovery ==========
@guard("classical_transport")
def classical_transport():
    trk = pd.read_csv(f"{CL}/observables/electron_track.csv").drop_duplicates("step")
    t, z, ke = trk["time_au"].values, trk["z"].values, trk["ke_ion_ha"].values
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.5))
    a1.plot(t, z, "C3-", lw=1.8)
    for zz, lab, c in [(SLAB_FACE,"slab face","0.4"),(-SLAB_FACE,None,"0.4"),(0,"slab centre","C0"),
                       (CAP_IN,"CAP","C2"),(-CAP_IN,None,"C2"),(BOX,"box edge","0.6"),(-BOX,None,"0.6")]:
        a1.axhline(zz, ls="--", lw=.8, color=c, label=lab)
    a1.axvline(T_EXIT, ls=":", lw=1.0, color="C2")
    a1.set_xlabel("time (a.u.)"); a1.set_ylabel("ion z (Bohr)")
    a1.set_title("Classical ion position z(t)  (wraps past +45 ⇒ periodic re-entry)", fontsize=8.5)
    a1.legend(fontsize=6.5, frameon=False, loc="upper left"); a1.grid(alpha=.25)
    iout = int(np.argmax(z >= BOX)) or len(z)
    a2.plot(z[:iout], ke[:iout]*HA, "C3.-", lw=1.4, ms=3)
    a2.axvspan(-SLAB_FACE, SLAB_FACE, color="C0", alpha=.10, label="slab")
    kmin_i = int(np.argmin(ke[:iout]))
    a2.plot(z[kmin_i], ke[kmin_i]*HA, "ko", ms=6)
    a2.annotate(f"min at z={z[kmin_i]:+.1f}\nKE={ke[kmin_i]*HA:.1f} eV", (z[kmin_i], ke[kmin_i]*HA),
                textcoords="offset points", xytext=(8, 18), fontsize=7, arrowprops=dict(arrowstyle="->", lw=.7))
    kin = np.interp(-SLAB_FACE, z, ke); kout = np.interp(SLAB_FACE, z, ke)
    for zz in (-SLAB_FACE, SLAB_FACE): a2.axvline(zz, ls="--", lw=.8, color="0.4")
    S_face = (kin - kout)*HA / L_SLAB
    a2.annotate(f"equal-potential window (±12.5):\nΔKE={(kin-kout)*HA:.1f} eV ⇒ S={S_face:.3f} eV/Bohr",
                (0, (kin+kout)*HA/2), fontsize=7, ha="center",
                bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=.85))
    a2.set_xlabel("ion z (Bohr)"); a2.set_ylabel("projectile KE (eV)")
    a2.set_title("KE(z): conservative dip-and-recovery — only equal-potential loss is stopping", fontsize=8)
    a2.legend(fontsize=7, frameon=False); a2.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(f"{FIGS}/classical_transport.png", dpi=180); plt.close(fig)
    R["classical_ke_min_eV"] = float(ke[kmin_i]*HA); R["classical_ke_min_z"] = float(z[kmin_i])
    R["classical_S_facewindow_eVbohr"] = float(S_face)
    log(f"classical transport: KE min {ke[kmin_i]*HA:.1f} eV at z={z[kmin_i]:+.1f}; face-window S={S_face:.3f} eV/Bohr")

# ===================== 4. CENTROID + SPREADING ==============================
@guard("centroid")
def centroid():
    s = pd.read_csv(f"{WP}/observables/wp_real_space_stats.csv", comment="#")
    t, zc, sz = s["time_au"].values, s["z_mean"].values, np.sqrt(s["sigma_z2"].values)
    R["anchor_table"] = anchors(t, zc, sz)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 6.2), sharex=True)
    a1.fill_between(t, zc-sz, zc+sz, color="C0", alpha=.18); a1.plot(t, zc, "C0-", lw=1.8)
    for zz, c in [(SLAB_FACE,"0.4"),(-SLAB_FACE,"0.4"),(CAP_IN,"C2"),(-CAP_IN,"C2")]:
        a1.axhline(zz, ls="--", lw=.8, color=c)
    a1.set_ylabel("⟨z⟩ ± σ_z (Bohr)"); a1.set_title("WP centroid + spread", fontsize=9); a1.grid(alpha=.25)
    a2.plot(t, sz, "C1-", lw=1.6); a2.set_xlabel("time (a.u.)"); a2.set_ylabel("σ_z (Bohr)")
    a2.set_title(f"σ_z: {sz[0]:.2f}→{sz[-1]:.2f} ({sz[-1]/sz[0]:.0f}×)", fontsize=9); a2.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(f"{FIGS}/centroid.png", dpi=180); plt.close(fig)
    R["sigma_z_growth"] = float(sz[-1]/sz[0]); R["zc_max"] = float(zc.max())

def anchors(t, zc, sz):
    def cross(sig, tgt):
        d = sig - tgt
        for i in range(1, len(d)):
            if d[i-1] <= 0 <= d[i] or d[i-1] >= 0 >= d[i]:
                f = -d[i-1]/(d[i]-d[i-1]) if d[i] != d[i-1] else 0
                return float(t[i-1]+f*(t[i]-t[i-1]))
        return None
    return {"A1_lead_near": cross(zc+3*sz, -SLAB_FACE), "A2_cen_near": cross(zc, -SLAB_FACE),
            "A4_lead_far": cross(zc+3*sz, SLAB_FACE), "A5_cen_max": float(t[np.argmax(zc)]),
            "t_final": float(t[-1])}

# ===================== 5. REFLECTION — signed J_z at CAP edges ===============
@guard("reflection")
def reflection():
    wff = sorted(glob.glob(f"{WP}/vti/wavefunction_wp/wavefunction_t*.vti"), key=ftime)
    if not wff: log("[reflection] no wavefunction frames"); return
    t = np.array([ftime(p) for p in wff]); fp = np.zeros(len(wff)); fm = np.zeros(len(wff))
    for i, p in enumerate(wff):
        psi, z, (sx, sy, sz) = load_psi(p); Jz = np.imag(np.conj(psi)*dpsi_dz(psi, sz))
        ip = int(np.argmin(np.abs(z-CAP_IN))); im = int(np.argmin(np.abs(z+CAP_IN)))
        fp[i] = Jz[:, :, ip].sum()*sx*sy; fm[i] = Jz[:, :, im].sum()*sx*sy
    from scipy.integrate import cumulative_trapezoid
    Tc = cumulative_trapezoid(fp, t, initial=0); Rc = cumulative_trapezoid(-fm, t, initial=0)
    fig, ax = plt.subplots(figsize=(7, 4.3))
    ax.plot(t, Tc, "C0-", lw=1.7, label=f"+z CAP (transmit): {Tc[-1]:+.3f}")
    ax.plot(t, Rc, "C3-", lw=1.7, label=f"−z CAP (reflect): {Rc[-1]:+.3f}")
    ax.axhline(0, color="0.6", lw=.6); ax.set_xlabel("time (a.u.)"); ax.set_ylabel("cumulative ∫J_z dt (electrons)")
    ax.set_title("Per-CAP flux — transmission vs reflection (new two-sided CAP geometry)", fontsize=9)
    ax.legend(fontsize=7, frameon=False); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(f"{FIGS}/reflection.png", dpi=180); plt.close(fig)
    R["transmit_cum"] = float(Tc[-1]); R["reflect_cum"] = float(Rc[-1])
    log(f"flux: transmit(+z)={Tc[-1]:+.3f}  reflect(−z)={Rc[-1]:+.3f}  (reflectivity for the SEAM-free two-sided CAP)")

# ===================== 6. MOMENTUM n(k_z) at anchors ========================
@guard("momentum")
def momentum():
    wff = sorted(glob.glob(f"{WP}/vti/wavefunction_wp/wavefunction_t*.vti"), key=ftime)
    if not wff: return
    t = np.array([ftime(p) for p in wff]); anc = R.get("anchor_table", {})
    want = {"t0": 0.0, "A2 slab entry": anc.get("A2_cen_near"), "A4 far face": anc.get("A4_lead_far"),
            "A5 max": anc.get("A5_cen_max"), "final": t[-1]}
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for lab, tv in want.items():
        if tv is None: continue
        p = wff[int(np.argmin(np.abs(t-tv)))]
        psi, z, (sx, sy, sz) = load_psi(p)
        psik = np.fft.fftshift(np.fft.fftn(psi)); nk = np.abs(psik)**2
        kz = np.fft.fftshift(np.fft.fftfreq(psi.shape[2], d=sz))*2*np.pi
        nkz = nk.sum(axis=(0, 1)); nkz /= nkz.sum()*(kz[1]-kz[0])
        ax.plot(kz, nkz, lw=1.5, label=f"{lab} (t={tv:.1f})")
    ax.axvline(0, color="0.6", lw=.6); ax.axvline(V0, ls=":", color="0.5", lw=.8)
    ax.set_xlim(-6, 8); ax.set_xlabel("k_z (a.u.)  — k_z<0 = reflected"); ax.set_ylabel("n(k_z)")
    ax.set_title("Sign-resolved WP momentum n(k_z) at scattering anchors", fontsize=9)
    ax.legend(fontsize=7, frameon=False); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(f"{FIGS}/momentum_nkz.png", dpi=180); plt.close(fig)
    log("momentum n(k_z) at anchors written")

# ===================== 7. ENERGETICS (both runs) ============================
@guard("energetics")
def energetics():
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.3))
    for ax, path, ttl in [(axs[0], WP, "WP"), (axs[1], CL, "classical")]:
        o = pd.read_csv(f"{path}/observables/observables.csv"); t = o["time_au"]
        for col, c in [("energy_total","C0"),("energy_kinetic","C1"),
                       ("energy_hartree","C2"),("energy_xc","C3")]:
            if col in o: ax.plot(t, (o[col]-o[col].iloc[0])*HA, c, lw=1.3, label=col.replace("energy_",""))
        ax.set_xlabel("time (a.u.)"); ax.set_ylabel("ΔE component (eV)")
        sub = "monotonic drain (CAP absorbs)" if ttl == "WP" else "rises: ion not absorbed → re-entry"
        ax.set_title(f"{ttl} energetics — {sub}", fontsize=8.5)
        ax.legend(fontsize=7, frameon=False); ax.grid(alpha=.25); mark_exit(ax)
    fig.suptitle(f"CAP: η={CAP_ETA} Ha, {CAP_W_BOHR:.0f} Bohr/side (two-sided), region [±{CAP_IN:.0f},±{BOX:.0f}]   "
                 f"(dashed = mean-v slab exit t={T_EXIT:.1f} a.u.)", fontsize=8, y=1.02)
    fig.tight_layout(); fig.savefig(f"{FIGS}/energetics.png", dpi=180, bbox_inches="tight"); plt.close(fig)
    log("energetics written")

# ===================== 8. KS excitation (overlap_full) ======================
@guard("excitation")
def excitation():
    idxf = sorted(glob.glob(f"{WP}/observables/overlap_full/overlap_*.csv"))
    if len(idxf) < 1: log("[excitation] no overlap_full"); return
    M0 = pd.read_csv(idxf[0], comment="#", header=None).values   # skip "# step=..." header line
    Mf = pd.read_csv(idxf[-1], comment="#", header=None).values
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, M, ttl in [(axs[0], M0, "t=0 (identity)"), (axs[1], Mf, "t_final (excited)")]:
        im = ax.imshow(np.log10(np.maximum(M, 1e-6)), origin="lower", aspect="auto", cmap=style.cmap_for("sequential"))
        ax.set_xlabel("evolved orbital j"); ax.set_ylabel("GS orbital i")
        ax.set_title(f"|⟨ψ_i^GS|ψ_j(t)⟩|²  {ttl}", fontsize=9); fig.colorbar(im, ax=ax, label="log₁₀")
    fig.tight_layout(); fig.savefig(f"{FIGS}/excitation.png", dpi=180); plt.close(fig)
    nd = min(Mf.shape); offdiag = Mf[:nd, :nd].sum() - np.trace(Mf[:nd, :nd])
    R["ks_offdiag_final"] = float(offdiag); log(f"KS excitation off-diagonal weight (final) = {offdiag:.3f}")

# ===================== 9. E-field (post-hoc) ===============================
@guard("efield")
def efield():
    from inqview import load_vti
    frames = sorted(glob.glob(f"{WP}/vti/density_delta/*.vti"))
    if not frames: log("[efield] no density_delta"); return
    try:
        from inqview.analysis import efield as EF
        d = load_vti(frames[len(frames)//2])
        ef = EF.electric_field(d.data, (d.x[1]-d.x[0], d.y[1]-d.y[0], d.z[1]-d.z[0]))
        Ez = ef.ez
    except Exception as e:
        log(f"[efield] skipped: {e}"); return
    iy = len(d.y)//2
    fig, ax = plt.subplots(figsize=(7, 4.6))
    im = ax.imshow(Ez[:, iy, :].T, origin="lower", aspect="auto",
                   extent=[d.x[0], d.x[-1], d.z[0], d.z[-1]], cmap="RdBu_r")
    ax.set_xlabel("x (Bohr)"); ax.set_ylabel("z (Bohr)"); ax.set_title("Induced E_z (mid-time)", fontsize=9)
    fig.colorbar(im, ax=ax, label="E_z"); fig.tight_layout()
    fig.savefig(f"{FIGS}/efield.png", dpi=180); plt.close(fig); log("E-field map written")

# ============= 10. norm / boundary-absorption vs time (both runs) ===========
@guard("norm_absorption")
def norm_absorption_fig():
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    wn = pd.read_csv(f"{WP}/observables/electron_number.csv")
    cn = pd.read_csv(f"{CL}/observables/electron_number.csv")
    ax.plot(wn.time_au, wn.N_total - wn.N_total.iloc[0], "C0-", lw=1.7,
            label=f"WP: ΔN_total (absorbed {wn.N_total.iloc[0]-wn.N_total.iloc[-1]:.3f} e)")
    ax.plot(cn.time_au, cn.N_total - cn.N_total.iloc[0], "C3-", lw=1.7,
            label=f"classical: ΔN_total (absorbed {cn.N_total.iloc[0]-cn.N_total.iloc[-1]:.3f} e)")
    s = pd.read_csv(f"{WP}/observables/wp_real_space_stats.csv", comment="#")
    axr = ax.twinx()
    axr.plot(s.time_au, s.norm_check, "C0--", lw=1.2, alpha=.7,
             label=f"WP orbital norm (1→{s.norm_check.iloc[-1]:.3f})")
    axr.set_ylabel("WP orbital norm", color="C0"); axr.set_ylim(0, 1.05)
    ax.axhline(0, color="0.6", lw=.6); mark_exit(ax)
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("ΔN_total = −(boundary-absorbed e)")
    ax.set_title("Total norm & boundary absorption vs time", fontsize=9)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = axr.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, fontsize=7, frameon=False, loc="lower left"); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(f"{FIGS}/norm_absorption.png", dpi=180); plt.close(fig)
    log("norm/absorption figure written")

# ===================== 11. density GIF battery ==============================
@guard("gif")
def gifs():
    geom = dict(dt=DT, slab_face=SLAB_FACE, cap_inner=CAP_IN)
    cl_gifs, dvmax = make_density_gif_battery(CL_RUN, FIGS, run_label="classical", **geom)
    wp_gifs, _ = make_density_gif_battery(WP_RUN, FIGS, run_label="wp", density_vmax=dvmax, **geom)
    R["gif_density_vmax"] = float(dvmax) if dvmax else None
    R["gifs"] = {"wp": [[c, k, os.path.basename(p)] for c, k, p, _ in wp_gifs],
                 "classical": [[c, k, os.path.basename(p)] for c, k, p, _ in cl_gifs]}
    log(f"density GIF battery: {len(wp_gifs)} WP + {len(cl_gifs)} classical")

# ===================== 12. heuristics (groups A–I) ==========================
@guard("heuristics")
def heuristics():
    hwp = compute_heuristics(WP_RUN, rs=RS, v0=V0, z0=Z0, slab_half=SLAB_FACE, box_half=BOX, sigma_wp=SIGMA_WP)
    hcl = compute_heuristics(CL_RUN, rs=RS, v0=V0, z0=Z0, slab_half=SLAB_FACE, box_half=BOX)
    R["heuristics"] = {"wp": hwp.flat(), "classical": hcl.flat()}
    R["T_exit_au"] = T_EXIT; R["T_enter_au"] = T_ENTER
    eg = hwp.eg_scales
    n_periods = TAU_AU / eg['T_plasmon_au']
    log(f"heuristics: kF={eg['kF']:.3f} vF={eg['vF']:.3f} ω_p={eg['omega_p_ev']:.2f} eV "
        f"T_plasmon={eg['T_plasmon_au']:.1f} au (τ={TAU_AU:.0f} ⇒ ~{n_periods:.1f} plasmon periods); "
        f"t_exit={T_EXIT:.2f} au")

_fns = [convergence, energy_method, classical, classical_transport, centroid,
        reflection, momentum, energetics, excitation, efield,
        norm_absorption_fig, heuristics, gifs]
if os.environ.get("SKIP_GIFS"):
    _fns = [f for f in _fns if f is not gifs]
    log("[SKIP_GIFS] density-GIF battery skipped")
for fn in _fns:
    fn()

for tag, path in [("wp", WP), ("classical", CL)]:
    try:
        s = open(glob.glob(f"{path}/../run_summary.txt")[0]).read()
        m = re.search(r"wall_time_s = ([\d.]+)", s)
        if m: R[f"{tag}_wall_s"] = float(m.group(1))
    except Exception: pass

json.dump(R, open(f"{HERE}/results.json", "w"), indent=2)
log(f"\nwrote {HERE}/results.json  ({len(glob.glob(f'{FIGS}/*'))} figures, TAG={TAG})")
