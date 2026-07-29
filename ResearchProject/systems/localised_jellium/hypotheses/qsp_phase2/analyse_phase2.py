#!/usr/bin/env python3
"""Phase-2 P2.1 analysis — WP + classical convergence/CAP test (autonomous).

Produces figures (figs/*.png, figs/*.gif) + results.json consumed by
build_phase2_notebook.py. Each block is guarded: a partial failure logs and
continues so the notebook still builds from whatever succeeded.

Deliverables (objectives A): E_total convergence (WP) + residual norm; classical
KE steady state; CAP absorption + reflection (signed J_z); centroid/spreading;
stopping (WP ΔE_total/L_z gated; classical ΔKE_ion) vs point-Lindhard, SIE-aware;
momentum n(k_z)/n(k_⊥) at scattering anchors; E-field; KS excitation (overlap_full);
per-run energetics + xz density GIF.
"""
import os, glob, re, json, traceback
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.visualisation import style, make_density_gif_battery
from inqview.analysis import compute_heuristics
style.apply_theme()

HA = 27.211386
ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
WP_RUN = f"{ROOT}/scripts/qsp_phase2/wp/results/p2_wp"          # run dir (for batteries)
CL_RUN = f"{ROOT}/scripts/qsp_phase2/classical/results/p2_classical"
WP = f"{WP_RUN}/raw"
CL = f"{CL_RUN}/raw"
HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = f"{HERE}/figs"; os.makedirs(FIGS, exist_ok=True)
DT = 0.02
E_GS = -45.75884855          # Phase-1 GS (Ha)
SIE_EV = 4.40                # Phase-1 SIE floor (eV)
SLAB_FACE, CAP_IN, BOX = 12.5, 25.0, 35.0
L_SLAB = 25.0                # traversal length for dE/dx (slab thickness)
# projectile kinematics + CAP provenance (compiled value, not the stale summary string)
Z0, V0, RS = -22.0, 2.711, 5.666
CAP_ETA, CAP_W_BOHR = -0.7, 10.0
# slab-exit time at mean velocity: reach the FAR slab face (+SLAB_FACE) from z0
T_EXIT = (SLAB_FACE - Z0) / V0          # ≈ 12.73 a.u.
T_ENTER = (-SLAB_FACE - Z0) / V0        # ≈ 3.50 a.u.
def mark_exit(ax):
    """Dashed marker at the mean-velocity slab-entry and slab-exit times."""
    ax.axvline(T_ENTER, ls=":", lw=.9, color="0.55")
    ax.axvline(T_EXIT, ls="--", lw=1.1, color="C2")
    yl = ax.get_ylim()
    ax.text(T_EXIT, yl[1], f" slab end\n t={T_EXIT:.1f}", fontsize=6.5,
            color="C2", va="top", ha="left")
R = {"notes": []}            # results dict for the notebook
def log(m): print(m); R["notes"].append(m)
def guard(name):
    def deco(fn):
        def wrap():
            try: fn()
            except Exception as e:
                log(f"[SKIP {name}] {e}"); traceback.print_exc()
        return wrap
    return deco

# ---- Lindhard point reference at r_s=5.67 ----
try:
    from inqview.analysis import lindhard_elf as L
    kF = L.kF_from_rs(5.666)
    R["S_lindhard_point_100eV_eVbohr"] = HA * L.stopping_power_point(2.711, kF)
    log(f"point-Lindhard S(100eV, r_s5.67) = {R['S_lindhard_point_100eV_eVbohr']:.3f} eV/Bohr")
except Exception as e:
    log(f"[lindhard skip] {e}"); R["S_lindhard_point_100eV_eVbohr"] = None

# ---- complex-field helpers (from qa_ii, validated) ----
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
    # late-time slope = convergence gate
    m = t >= 0.8*t.max(); slope = np.polyfit(t[m], E[m], 1)[0]*HA
    a2.set_title(f"N_total (late dE/dt = {slope:+.2f} eV/a.u.; ~0 ⇒ converged)", fontsize=9)
    mark_exit(a1); mark_exit(a2)
    fig.tight_layout(); fig.savefig(f"{FIGS}/conv_wp.png", dpi=180); plt.close(fig)
    R["wp_dEtot_eV"] = float((E[-1]-E[0])*HA); R["wp_late_slope_eV_au"] = float(slope)
    R["wp_Ntot_final"] = float(n["N_total"].values[-1]); R["t_final"] = float(t[-1])
    R["wp_converged"] = bool(abs(slope) < 1.0)
    log(f"WP: ΔE_total={R['wp_dEtot_eV']:.1f} eV, late slope={slope:+.2f} eV/au, "
        f"N_final={R['wp_Ntot_final']:.3f}, converged={R['wp_converged']}")

# ===================== 2. CLASSICAL STEADY STATE + STOPPING ==================
@guard("classical")
def classical():
    trk = pd.read_csv(f"{CL}/observables/electron_track.csv").drop_duplicates("step")
    t, z, ke = trk["time_au"].values, trk["z"].values, trk["ke_ion_ha"].values
    fig, ax = plt.subplots(figsize=(7, 4.3))
    ax.plot(t, ke*HA, "C3-", lw=1.6); ax.set_xlabel("time (a.u.)")
    ax.set_ylabel("projectile KE (eV)")
    ax.set_title("Classical projectile KE — steady state after slab exit?", fontsize=9)
    mark_exit(ax); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(f"{FIGS}/classical_ke.png", dpi=180); plt.close(fig)
    # stopping = KE loss across the slab traversal (z from -12.5 to +12.5)
    try:
        zin = np.interp(-SLAB_FACE, z, t); zout = np.interp(SLAB_FACE, z, t)
        ke_in = np.interp(zin, t, ke); ke_out = np.interp(zout, t, ke)
        S = (ke_in - ke_out)*HA / L_SLAB
        R["classical_S_eVbohr"] = float(S)
        log(f"classical: KE {ke_in*HA:.1f}→{ke_out*HA:.1f} eV over {L_SLAB} Bohr ⇒ S={S:.3f} eV/Bohr")
    except Exception as e:
        log(f"[classical stopping skip] {e}")
    R["classical_ke0_eV"] = float(ke[0]*HA); R["classical_kef_eV"] = float(ke[-1]*HA)

# ============ 2b. CLASSICAL TRANSPORT: z(t) + KE(z) dip-and-recovery ==========
@guard("classical_transport")
def classical_transport():
    """Position vs time + KE vs z — shows the projectile slows to a minimum at the
    slab CENTRE and recovers on exit (a conservative mean-field-potential effect:
    energy borrowed and returned). Only the net loss between two points at EQUAL
    background potential (the symmetric slab faces ±12.5) is true stopping."""
    trk = pd.read_csv(f"{CL}/observables/electron_track.csv").drop_duplicates("step")
    t, z, ke = trk["time_au"].values, trk["z"].values, trk["ke_ion_ha"].values
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.5))
    # (a) z(t)
    a1.plot(t, z, "C3-", lw=1.8)
    for zz, lab, c in [(SLAB_FACE,"slab face","0.4"),(-SLAB_FACE,None,"0.4"),(0,"slab centre","C0"),
                       (CAP_IN,"CAP","C2"),(-CAP_IN,None,"C2"),(BOX,"box edge","0.6"),(-BOX,None,"0.6")]:
        a1.axhline(zz, ls="--", lw=.8, color=c, label=lab)
    a1.axvline(T_EXIT, ls=":", lw=1.0, color="C2")
    a1.set_xlabel("time (a.u.)"); a1.set_ylabel("ion z (Bohr)")
    a1.set_title("Classical ion position z(t)  (wraps past +35 ⇒ periodic re-entry)", fontsize=8.5)
    a1.legend(fontsize=6.5, frameon=False, loc="upper left"); a1.grid(alpha=.25)
    # (b) KE(z) — outbound leg until first box-edge crossing (before wrap)
    iout = np.argmax(z >= BOX) or len(z)
    a2.plot(z[:iout]*1, ke[:iout]*HA, "C3.-", lw=1.4, ms=3)
    a2.axvspan(-SLAB_FACE, SLAB_FACE, color="C0", alpha=.10, label="slab")
    kmin_i = int(np.argmin(ke[:iout]))
    a2.plot(z[kmin_i], ke[kmin_i]*HA, "ko", ms=6)
    a2.annotate(f"min at z={z[kmin_i]:+.1f}\nKE={ke[kmin_i]*HA:.1f} eV", (z[kmin_i], ke[kmin_i]*HA),
                textcoords="offset points", xytext=(8, 18), fontsize=7,
                arrowprops=dict(arrowstyle="->", lw=.7))
    kin = np.interp(-SLAB_FACE, z, ke); kout = np.interp(SLAB_FACE, z, ke)
    for zz in (-SLAB_FACE, SLAB_FACE):
        a2.axvline(zz, ls="--", lw=.8, color="0.4")
    S_face = (kin - kout)*HA / L_SLAB
    a2.annotate(f"equal-potential window (±12.5):\nΔKE={ (kin-kout)*HA:.1f} eV ⇒ S={S_face:.3f} eV/Bohr",
                (0, (kin+kout)*HA/2), fontsize=7, ha="center",
                bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=.85))
    a2.set_xlabel("ion z (Bohr)"); a2.set_ylabel("projectile KE (eV)")
    a2.set_title("KE(z): conservative dip-and-recovery — only equal-potential loss is stopping", fontsize=8)
    a2.legend(fontsize=7, frameon=False); a2.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(f"{FIGS}/classical_transport.png", dpi=180); plt.close(fig)
    R["classical_ke_min_eV"] = float(ke[kmin_i]*HA); R["classical_ke_min_z"] = float(z[kmin_i])
    R["classical_S_facewindow_eVbohr"] = float(S_face)
    log(f"classical transport: KE min {ke[kmin_i]*HA:.1f} eV at z={z[kmin_i]:+.1f}; "
        f"recovers to {kout*HA:.1f} eV at +12.5; face-window S={S_face:.3f} eV/Bohr "
        f"(conservative well borrowed/returned)")

# ===================== 3. WP STOPPING (ΔE_total/L_z, gated) ==================
@guard("wp_stopping")
def wp_stopping():
    o = pd.read_csv(f"{WP}/observables/observables.csv")
    Ef = float(o["energy_total"].values[-1])
    deposited = (Ef - E_GS)*HA                  # exact only at full absorption
    R["wp_deposited_EminusEGS_eV"] = float(deposited)
    R["wp_S_eVbohr"] = float(deposited / L_SLAB)
    gate = R.get("wp_converged", False)
    log(f"WP deposited (E_total(final)−E_GS) = {deposited:.1f} eV ⇒ S={deposited/L_SLAB:.2f} eV/Bohr "
        f"[{'CONVERGED' if gate else 'NOT converged ⇒ UPPER BOUND only'}]; SIE floor {SIE_EV} eV")

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
    ax.axhline(0, color="0.6", lw=.6); ax.set_xlabel("time (a.u.)")
    ax.set_ylabel("cumulative ∫J_z dt (electrons)")
    ax.set_title("Per-CAP flux — transmission vs reflection", fontsize=9)
    ax.legend(fontsize=7, frameon=False); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(f"{FIGS}/reflection.png", dpi=180); plt.close(fig)
    R["transmit_cum"] = float(Tc[-1]); R["reflect_cum"] = float(Rc[-1])
    log(f"flux: transmit(+z)={Tc[-1]:+.3f}  reflect(−z)={Rc[-1]:+.3f}")

# ===================== 6. MOMENTUM n(k_z) at anchors ========================
@guard("momentum")
def momentum():
    wff = sorted(glob.glob(f"{WP}/vti/wavefunction_wp/wavefunction_t*.vti"), key=ftime)
    if not wff: return
    t = np.array([ftime(p) for p in wff])
    anc = R.get("anchor_table", {})
    want = {"t0": 0.0, "A2 slab entry": anc.get("A2_cen_near"), "A4 far face": anc.get("A4_lead_far"),
            "A5 max": anc.get("A5_cen_max"), "final": t[-1]}
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for lab, tv in want.items():
        if tv is None: continue
        p = wff[int(np.argmin(np.abs(t-tv)))]
        psi, z, (sx, sy, sz) = load_psi(p)
        psik = np.fft.fftshift(np.fft.fftn(psi))
        nk = (np.abs(psik)**2)
        kz = np.fft.fftshift(np.fft.fftfreq(psi.shape[2], d=sz))*2*np.pi
        nkz = nk.sum(axis=(0, 1)); nkz /= nkz.sum()*(kz[1]-kz[0])
        ax.plot(kz, nkz, lw=1.5, label=f"{lab} (t={tv:.1f})")
    ax.axvline(0, color="0.6", lw=.6); ax.axvline(2.711, ls=":", color="0.5", lw=.8)
    ax.set_xlim(-6, 8); ax.set_xlabel("k_z (a.u.)  — k_z<0 = reflected")
    ax.set_ylabel("n(k_z)"); ax.set_title("Sign-resolved WP momentum n(k_z) at scattering anchors", fontsize=9)
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
        sub = "monotonic drain (CAP absorbs)" if ttl == "WP" else "rises: ion not absorbed → periodic re-entry"
        ax.set_title(f"{ttl} energetics — {sub}", fontsize=8.5)
        ax.legend(fontsize=7, frameon=False); ax.grid(alpha=.25); mark_exit(ax)
    fig.suptitle(f"CAP: η={CAP_ETA} Ha, {CAP_W_BOHR:.0f} Bohr/side, region [±{CAP_IN:.0f},±{BOX:.0f}]   "
                 f"(dashed = mean-v slab exit t={T_EXIT:.1f} a.u.)", fontsize=8, y=1.02)
    fig.tight_layout(); fig.savefig(f"{FIGS}/energetics.png", dpi=180, bbox_inches="tight"); plt.close(fig)
    log("energetics written")

# ===================== 8. KS excitation (overlap_full) ======================
@guard("excitation")
def excitation():
    idxf = sorted(glob.glob(f"{WP}/observables/overlap_full/overlap_*.csv"))
    if len(idxf) < 1: log("[excitation] no overlap_full"); return
    M0 = pd.read_csv(idxf[0]).select_dtypes("number").values
    Mf = pd.read_csv(idxf[-1]).select_dtypes("number").values
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, M, ttl in [(axs[0], M0, "t=0 (identity)"), (axs[1], Mf, "t_final (excited)")]:
        im = ax.imshow(np.log10(np.maximum(M, 1e-6)), origin="lower", aspect="auto", cmap=style.cmap_for("sequential"))
        ax.set_xlabel("evolved orbital j"); ax.set_ylabel("GS orbital i")
        ax.set_title(f"|⟨ψ_i^GS|ψ_j(t)⟩|²  {ttl}", fontsize=9); fig.colorbar(im, ax=ax, label="log₁₀")
    fig.tight_layout(); fig.savefig(f"{FIGS}/excitation.png", dpi=180); plt.close(fig)
    # off-diagonal weight at final = excitation magnitude
    nd = min(Mf.shape); offdiag = Mf[:nd, :nd].sum() - np.trace(Mf[:nd, :nd])
    R["ks_offdiag_final"] = float(offdiag); log(f"KS excitation off-diagonal weight (final) = {offdiag:.3f}")

# ===================== 9. E-field (post-hoc) ===============================
@guard("efield")
def efield():
    from inqview.analysis import efield as EF
    from inqview import load_vti
    frames = sorted(glob.glob(f"{WP}/vti/density_delta/*.vti"))
    if not frames: log("[efield] no density_delta"); return
    d = load_vti(frames[len(frames)//2])
    try:
        Ex, Ey, Ez = EF.efield_from_density(d.data, (d.x[1]-d.x[0], d.y[1]-d.y[0], d.z[1]-d.z[0]))
    except Exception:
        # generic fallback: gradient of Poisson(density)
        log("[efield] efield API fallback skipped"); return
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
    # WP orbital norm (right axis)
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

# ===================== 11. density GIF battery (9 WP + 3 classical) ==========
@guard("gif")
def gifs():
    geom = dict(dt=DT, slab_face=SLAB_FACE, cap_inner=CAP_IN)
    # classical first → its slab-tuned density scale is reused for the WP totals
    cl_gifs, dvmax = make_density_gif_battery(
        CL_RUN, FIGS, run_label="classical", **geom)
    wp_gifs, _ = make_density_gif_battery(
        WP_RUN, FIGS, run_label="wp", density_vmax=dvmax, **geom)
    R["gif_density_vmax"] = float(dvmax) if dvmax else None
    R["gifs"] = {"wp": [[c, k, os.path.basename(p)] for c, k, p, _ in wp_gifs],
                 "classical": [[c, k, os.path.basename(p)] for c, k, p, _ in cl_gifs]}
    log(f"density GIF battery: {len(wp_gifs)} WP + {len(cl_gifs)} classical "
        f"(shared total/bath density vmax={dvmax:.2e})")

# ===================== 12. heuristics (groups A–I) ==========================
@guard("heuristics")
def heuristics():
    hwp = compute_heuristics(WP_RUN, rs=RS, v0=V0, z0=Z0, slab_half=SLAB_FACE,
                             box_half=BOX, sigma_wp=0.5)
    hcl = compute_heuristics(CL_RUN, rs=RS, v0=V0, z0=Z0, slab_half=SLAB_FACE,
                             box_half=BOX)
    R["heuristics"] = {"wp": hwp.flat(), "classical": hcl.flat()}
    R["T_exit_au"] = T_EXIT; R["T_enter_au"] = T_ENTER
    eg = hwp.eg_scales
    log(f"heuristics: kF={eg['kF']:.3f} vF={eg['vF']:.3f} ω_p={eg['omega_p_ev']:.2f} eV "
        f"T_plasmon={eg['T_plasmon_au']:.1f} au (>τ=40 ⇒ plasmon under-resolved); "
        f"t_exit={T_EXIT:.2f} au; spread×{hwp.spreading.get('spread_factor',0):.0f}; "
        f"total_absorbed WP={hwp.norms.get('total_absorbed',0):.3f} "
        f"classical={hcl.norms.get('total_absorbed',0):.3f}")

_fns = [convergence, classical, classical_transport, wp_stopping, centroid,
        reflection, momentum, energetics, excitation, efield,
        norm_absorption_fig, heuristics, gifs]
if os.environ.get("SKIP_GIFS"):           # fast re-render of PNGs/results only
    _fns = [f for f in _fns if f is not gifs]
    log("[SKIP_GIFS] density-GIF battery skipped (reusing existing figs/*.gif)")
for fn in _fns:
    fn()

# ---- wall times (cost) ----
for tag, path in [("wp", WP), ("classical", CL)]:
    try:
        s = open(glob.glob(f"{path}/../run_summary.txt")[0]).read()
        m = re.search(r"wall_time_s = ([\d.]+)", s)
        if m: R[f"{tag}_wall_s"] = float(m.group(1))
    except Exception: pass

json.dump(R, open(f"{HERE}/results.json", "w"), indent=2)
log(f"\nwrote {HERE}/results.json  ({len(glob.glob(f'{FIGS}/*'))} figures)")
