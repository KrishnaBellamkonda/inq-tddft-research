#!/usr/bin/env python3
"""qsp_phase5 — cumulative S(E) plot (rebuilt + emailed after every WP run).

Overlays, on a single log-x energy axis (E = ½k₀²·HA, drift energy):
  • the quantum WP points measured so far (se_state.csv) — filled = converged,
    open ▽ + ↓ = UPPER bound (WP not fully absorbed by τ);
  • the BULK classical σ_WP=0.5 reference = the σ_q=0.354 `sigma0p35` set
    (the √2 convention — NOT sigma0p5), extracted by the bulk slope method;
  • the BULK point-charge Lindhard line (r_s=5.69 ≈ slab 5.666);
  • the one localised park-method classical point (v=2.0, 0.249 eV/Bohr) as a
    geometry-matched check.

Classical + Lindhard are STATIC → computed once and cached. Run with --email to
send the latest plot (threaded under [lj-wp-se-sweep]).

Usage:  python3 build_se_plot.py [--email] [--note "v=3.0 done"]
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HA   = 27.211386
RS   = 5.69
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/local/data/public/skcb2/tddft"
STATE = os.path.join(HERE, "se_state.csv")
FIGS  = os.path.join(HERE, "figs"); os.makedirs(FIGS, exist_ok=True)
OUTPNG = os.path.join(FIGS, "se_quantum_stopping.png")
CL_CACHE = os.path.join(HERE, "classical_sigma0p5_bulk.csv")
LR_CACHE = os.path.join(HERE, "lindhard_ref.npz")
THREAD = os.path.join(HERE, "email_thread.txt")
PARK = dict(E_eV=54.422772, S=0.24922839012161602, v=2.0)   # phase-4 localised park method
SIG35_DIR = "run_classical_n162_L50_sv_sigma0p35"            # σ_q=0.354 ⇔ σ_WP=0.5
SIG35_VTAGS = [(0.2, "v0p2"), (0.6, "v0p6"), (0.8, "v0p8"),
               (1.3, "v1p3"), (2.0, "v2p0"), (3.0, "v3p0")]

sys.path.insert(0, os.path.join(ROOT, "inq-stack", "python"))
from inqview.analysis import lindhard_elf as LE         # noqa: E402
try:
    from inqview.visualisation import style as STYLE     # noqa: E402
    STYLE.apply_theme()
except Exception as exc:  # noqa: BLE001
    print(f"[build_se_plot] theme unavailable ({exc}); matplotlib defaults")

KF = LE.kF_from_rs(RS)


def _classical_bulk():
    """Bulk classical σ_WP=0.5 S(E) points (eV/Bohr), cached. Returns DataFrame."""
    if os.path.exists(CL_CACHE):
        return pd.read_csv(CL_CACHE)
    rows = []
    try:
        sd = os.path.join(ROOT, "ResearchProject/systems/jellium/hypotheses/06_sigma_convergence")
        sys.path.insert(0, sd)
        import sigma_sweep_report as SSR  # noqa: E402
        for v0, vtag in SIG35_VTAGS:
            if not SSR._completed(SIG35_DIR, vtag, "nested"):
                continue
            try:
                r = SSR.extract(SIG35_DIR, vtag, "nested")
            except Exception as exc:  # noqa: BLE001
                print(f"[build_se_plot] classical {vtag}: {exc}"); r = None
            if r is None:
                continue
            rows.append(dict(v0=v0, v=r["v"], E_eV=0.5 * v0 * v0 * HA,
                             S_eVbohr=r["S_A"] * HA, S_err=r["S_A_err"] * HA))
    except Exception as exc:  # noqa: BLE001
        print(f"[build_se_plot] WARN classical overlay unavailable ({exc})")
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values("E_eV").reset_index(drop=True)
        df.to_csv(CL_CACHE, index=False)
    return df


def _lindhard_curve(emax_eV):
    vmax = float(np.sqrt(2.0 * max(emax_eV, 130.0) / HA)) + 0.5
    if os.path.exists(LR_CACHE):
        d = np.load(LR_CACHE)
        if d["vmax"] >= vmax - 1e-6:
            return d["E"], d["S"]
    vg = np.linspace(0.4, vmax, 130)
    S = np.array([LE.stopping_power_point(float(v), KF) for v in vg]) * HA
    E = 0.5 * vg * vg * HA
    np.savez(LR_CACHE, E=E, S=S, vmax=np.array(vmax))
    return E, S


def build(email=False, note=""):
    if not os.path.exists(STATE):
        print("[build_se_plot] no se_state.csv yet — nothing to plot"); return None
    wp = pd.read_csv(STATE).sort_values("E_eV").reset_index(drop=True)
    emax = float(max(wp["E_eV"].max(), PARK["E_eV"], 130.0))
    cl = _classical_bulk()
    Elr, Slr = _lindhard_curve(emax * 1.15)

    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ax.plot(Elr, Slr, "-", color="0.45", lw=1.6, zorder=1,
            label=r"Lindhard point-charge (bulk, $r_s$=5.69)")
    if len(cl):
        ax.plot(cl["E_eV"], cl["S_eVbohr"], "--o", color="C0", ms=5, lw=1.0, mfc="white",
                mec="C0", zorder=2, label=r"classical (bulk, $\sigma_{WP}$=0.5)")
    ax.plot([PARK["E_eV"]], [PARK["S"]], "s", color="C2", ms=9, mec="k", mew=0.5,
            zorder=3, label="classical localised park (v=2.0)")

    # WP points: filled o = converged; ▽↓ = upper bound; gray ✗ = ALIASED (excluded).
    # bound=="lower" (positive late energy slope) is exactly the grid-aliased set.
    # NOTE: no on-plot text labels (per-point S / aliased tags) — user asked for
    # annotation-free plots; the values live in the notebook's Result table.
    phys = [float(r["S_eVbohr"]) for _, r in wp.iterrows() if str(r.get("bound", "")) != "lower"]
    ytop = (max(phys) * 1.3) if phys else float(wp["S_eVbohr"].max()) * 1.1
    for _, r in wp.iterrows():
        b = str(r.get("bound", "")); E, S = float(r["E_eV"]), float(r["S_eVbohr"])
        if b == "lower":                                   # aliased — exclude from physics
            ax.annotate("", xy=(E, ytop * 0.99), xytext=(E, ytop * 0.86),
                        arrowprops=dict(arrowstyle="-|>", color="0.6", lw=1.3, ls=":"), zorder=4)
            ax.plot([E], [min(S, ytop * 0.93)], "x", color="0.6", ms=9, mew=2, zorder=5)
            continue
        if b == "exact":
            ax.plot([E], [S], "o", color="C3", ms=9, mec="k", mew=0.5, zorder=5)
        else:
            ax.plot([E], [S], "v", color="C3", ms=10, mfc="white", mec="C3", mew=1.6, zorder=5)
            ax.annotate("", xy=(E, S * 0.82), xytext=(E, S),
                        arrowprops=dict(arrowstyle="-|>", color="C3", lw=1.4), zorder=4)
    ax.plot([], [], "o", color="C3", mec="k", label="WP quantum (converged)")
    ax.plot([], [], "v", color="C3", mfc="white", mec="C3", label="WP quantum (↓ upper bound)")
    if (wp["bound"].astype(str) == "lower").any():
        ax.plot([], [], "x", color="0.6", mew=2, label="WP aliased (excluded, off-scale)")

    ax.set_ylim(0, ytop)
    ax.set_xscale("log")
    ax.set_xlabel("projectile energy  E = ½k₀²  (eV)")
    ax.set_ylabel("electronic stopping power  S  (eV/Bohr)")
    ax.set_title("Quantum (wavepacket) stopping power S(E) — localised jellium slab, "
                 r"$\sigma_{WP}$=0.5", fontsize=10)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.5, frameon=False, loc="upper right")

    # (no on-plot per-point table or provenance footnote — annotation-free plot per
    #  user request; the per-point table lives in the notebook body.)
    fig.tight_layout()
    fig.savefig(OUTPNG, dpi=180)
    plt.close(fig)
    print(f"[build_se_plot] wrote {OUTPNG}  ({len(wp)} WP point(s), {len(cl)} classical)")

    if email:
        _email(wp, cl, note)
    return OUTPNG


def _email(wp, cl, note):
    try:
        from inqview.email import send_run_email
    except Exception as exc:  # noqa: BLE001
        print(f"[build_se_plot] EMAIL skipped (import): {exc}"); return
    n = len(wp)
    tbl = ["E (eV)   v     S (eV/Bohr)   convergence"]
    for _, r in wp.iterrows():
        b = str(r.get("bound", ""))
        flag = "converged" if b == "exact" else f"{b} bound (WP not fully absorbed)"
        tbl.append(f"{r['E_eV']:6.0f}  {r['v']:4.1f}    {r['S_eVbohr']:6.2f}     {flag}")
    subj = f"[lj-wp-se-sweep] S(E) update — {n} WP point(s)" + (f" — {note}" if note else "")
    body = (
        f"Quantum (wavepacket) electronic stopping power S(E), localised jellium slab, σ_WP=0.5.\n"
        f"{n} WP point(s) measured so far (energy method S=[E_total(t_f)−E_GS]/L_z, L_z=25 Bohr).\n\n"
        + "\n".join(tbl) +
        "\n\nOverlays (BULK references — geometry estimate, ADR 0010):\n"
        "  • classical σ_WP=0.5 (σ_q=0.354 bulk slope method)\n"
        "  • point-charge Lindhard (r_s=5.69)\n"
        f"  • localised park-method classical point: S(v=2.0)={PARK['S']:.3f} eV/Bohr\n\n"
        "Converged points are true values; ↓ upper-bound points have the WP not fully\n"
        "absorbed by τ (residual WP energy still draining; true S sits below the marker).\n"
    )
    kw = {}
    if os.path.exists(THREAD):
        root = open(THREAD).read().strip()
        if root:
            kw = dict(in_reply_to=root, references=[root])
    try:
        mid = send_run_email(subj, body, attachments=[OUTPNG],
                             to="chiddukanna@gmail.com", **kw)
        print(f"[build_se_plot] emailed: {mid}")
        if not os.path.exists(THREAD):
            open(THREAD, "w").write(mid)
    except Exception as exc:  # noqa: BLE001
        print(f"[build_se_plot] EMAIL FAILED: {exc}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", action="store_true")
    ap.add_argument("--note", default="")
    a = ap.parse_args()
    build(email=a.email, note=a.note)
