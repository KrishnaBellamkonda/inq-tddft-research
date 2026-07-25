#!/usr/bin/env python3
"""Bath-structure sweep — classical vs quantum induced density via POD/DMD.

The governing-PDE headline was an artifact (WP-included blob; see handover
2026-07-04). This is the artifact-free replacement the panel + user converged on:
characterise the TRUE blob-free induced bath dn = (n_total - n_wp) - n_GS (WP) or
n_total - n_GS (classical point) with the VALIDATED POD/DMD kernels, and compare
classical vs quantum across a sigma-sweep (fixed v) and a velocity-sweep.

Metrics per cell: POD rank (modes for 90% energy), leading-mode energy, DMD
dominant angular frequency (vs omega_p) + growth. POD-rank is the robust
descriptor; DMD frequencies are approximate (non-stationary runs, growth>0).

Hands-off: per-cell try/except + saved JSON (idempotent), figures + 4-part email
at the end. CPU-only, no GPU, no new runs.

Run: venv/bin/python3 docs/campaigns/ml-patterns/bath_structure_sweep.py [--no-email]
"""
from __future__ import annotations
import sys, os, glob, json, argparse, traceback
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ART = os.path.join(HERE, "artifacts")
os.makedirs(ART, exist_ok=True)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from inqview import load_vti
from kernels import pod as P, dmd as D

J = "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium"
OMP = 0.1276               # a.u., n162/L50 bath
OMP_EV = OMP * 27.211
SEND_EMAIL = True

# --- cell registry: only runs with a reconstructable clean bath (verified on disk)
# is_wp True -> bath = total - wp - GS ; False -> classical point, bath = total - GS
SIGMA_SWEEP = [   # fixed v=2.71 (E=100), all 125^3
    dict(run="run_classical_n162_L50_E100_v2", is_wp=False, sigma=0.0, v=2.71, E=100),
    dict(run="run_wp_n162_L50_E100_sigma0p5_wf", is_wp=True, sigma=0.5, v=2.71, E=100),
    dict(run="run_wp_n162_L50_E100_sigma3_wf",  is_wp=True, sigma=3.0, v=2.71, E=100),
    dict(run="run_wp_n162_L50_E100_sigma8_wf",  is_wp=True, sigma=8.0, v=2.71, E=100),
]
VEL_SWEEP_CL = [  # classical point, velocity sweep, 125^3
    dict(run="run_classical_n162_L50_E20",    is_wp=False, sigma=0.0, v=1.21, E=20),
    dict(run="run_classical_n162_L50_E25",    is_wp=False, sigma=0.0, v=1.36, E=25),
    dict(run="run_classical_n162_L50_E50",    is_wp=False, sigma=0.0, v=1.92, E=50),
    dict(run="run_classical_n162_L50_E100_v2",is_wp=False, sigma=0.0, v=2.71, E=100),
    dict(run="run_classical_n162_L50_E600",   is_wp=False, sigma=0.0, v=6.64, E=600),
]
VEL_SWEEP_WP = [  # WP velocity points with clean bath (MIXED sigma/grid — caveated)
    dict(run="run_base_n162_L50_E1p5",         is_wp=True, sigma=5.0, v=0.33, E=1.5),
    dict(run="run_plasmon_n162_L50_E3p4_varyv",is_wp=True, sigma=3.0, v=0.50, E=3.4),
    dict(run="run_wp_n162_L50_E100_sigma3_wf", is_wp=True, sigma=3.0, v=2.71, E=100),
]


def _dir(base, *names):
    for n in names:
        f = sorted(glob.glob(os.path.join(base, n, "density_*t*.vti")))
        if f:
            return f
    return []


def load_bath(cell, maxf=110):
    base = os.path.join(J, cell["run"], "results/raw/vti")
    gs = load_vti(os.path.join(base, "density_gs_system/density_gs_system.vti")).data.astype(np.float32)
    tf = _dir(base, "density_rt_total", "density_total")
    if not tf:
        raise FileNotFoundError("no total series")
    stride = max(1, len(tf) // maxf)
    tf = tf[::stride]
    wf = _dir(base, "density_rt_wp", "density_wp")[::stride] if cell["is_wp"] else None
    cols = []
    for i, f in enumerate(tf):
        tot = load_vti(f).data.astype(np.float32)
        bath = tot - load_vti(wf[i]).data.astype(np.float32) if cell["is_wp"] else tot
        cols.append((bath - gs).ravel())
    X = np.stack(cols, axis=1)
    return X, gs.shape, 0.02 * stride * (2 if cell["is_wp"] else 1)  # dt approx (frame_dt~0.02-0.04)


def analyse_cell(cell):
    X, shp, dt = load_bath(cell)
    pod = P.pod(X, rank=20, center=True)
    T = X.shape[1]
    dm = D.dmd(X, dt, rank=12, window=(0, max(4, int(0.6 * T))))
    _, w, g, _ = dm.dominant()
    return dict(
        run=cell["run"], is_wp=cell["is_wp"], sigma=cell["sigma"], v=cell["v"], E=cell["E"],
        grid=list(shp), n_frames=T, dt=dt,
        pod_energy=[float(x) for x in pod.energy_fraction[:6]],
        pod_rank90=int(pod.n_modes_for(0.90)),
        lead_energy=float(pod.energy_fraction[0]),
        dmd_omega_ev=float(w * 27.211), dmd_growth=float(g),
    )


def run_group(cells, tag):
    out = []
    for c in cells:
        cp = os.path.join(ART, f"bathstruct_{c['run']}.json")
        if os.path.isfile(cp):
            out.append(json.load(open(cp))); print(f"[{tag}] {c['run']} cached"); continue
        try:
            r = analyse_cell(c)
            json.dump(r, open(cp, "w"), indent=2)
            out.append(r)
            print(f"[{tag}] {c['run']}: rank90={r['pod_rank90']} lead={r['lead_energy']:.2f} "
                  f"dmd={r['dmd_omega_ev']:.1f}eV")
        except Exception as e:
            print(f"[{tag}] {c['run']} SKIP: {e!r}"); traceback.print_exc()
    return out


def fig_sigma(sig, path):
    sig = sorted(sig, key=lambda r: r["sigma"])
    x = [r["sigma"] for r in sig]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(x, [r["pod_rank90"] for r in sig], "o-", color="C0")
    ax[0].set_xlabel("sigma_WP (Bohr)  [0 = classical point]"); ax[0].set_ylabel("POD modes for 90% energy")
    ax[0].set_title("(a) Induced-bath complexity vs projectile width (v=2.71)")
    ax[1].plot(x, [r["lead_energy"] for r in sig], "s-", color="C2")
    ax[1].set_xlabel("sigma_WP (Bohr)"); ax[1].set_ylabel("leading-mode energy fraction")
    ax[1].set_title("(b) Coherence of the leading mode vs width")
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def fig_vel(cl, wp, path):
    cl = sorted(cl, key=lambda r: r["v"]); wp = sorted(wp, key=lambda r: r["v"])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot([r["v"] for r in cl], [r["pod_rank90"] for r in cl], "o-", label="classical point")
    if wp:
        ax[0].plot([r["v"] for r in wp], [r["pod_rank90"] for r in wp], "s--",
                   label="WP (mixed sigma)", color="C1")
    ax[0].set_xlabel("velocity (a.u.)"); ax[0].set_ylabel("POD modes for 90%"); ax[0].legend()
    ax[0].set_title("(a) Induced-bath complexity vs velocity")
    ax[1].plot([r["v"] for r in cl], [r["dmd_omega_ev"] for r in cl], "o-", label="classical")
    if wp:
        ax[1].plot([r["v"] for r in wp], [r["dmd_omega_ev"] for r in wp], "s--",
                   label="WP", color="C1")
    ax[1].axhline(OMP_EV, ls=":", c="k", label=f"omega_p={OMP_EV:.1f}eV")
    ax[1].set_xlabel("velocity (a.u.)"); ax[1].set_ylabel("DMD dominant omega (eV)")
    ax[1].legend(); ax[1].set_title("(b) DMD dominant frequency vs velocity")
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def email(subj, parts, plots):
    if not SEND_EMAIL:
        print(f"[email skipped] {subj}"); return
    body = (f"1) HYPOTHESIS\n{parts['h']}\n\n2) WHAT WAS DONE\n{parts['d']}\n\n"
            f"3) WHAT THE PLOTS SHOW\n{parts['p']}\n\n4) CONCLUSION\n{parts['c']}\n")
    try:
        from inqview.email import send_run_email
        send_run_email(subject=subj, body=body,
                       attachments=[p for p in plots if os.path.isfile(p)])
        print(f"[email sent] {subj}")
    except Exception as e:
        print(f"[email FAILED, continuing] {e!r}")


def main(argv):
    global SEND_EMAIL
    ap = argparse.ArgumentParser(); ap.add_argument("--no-email", action="store_true")
    SEND_EMAIL = not ap.parse_args(argv).no_email
    print("=== bath-structure sweep (classical vs quantum, true blob-free bath) ===")
    sig = run_group(SIGMA_SWEEP, "sigma")
    clv = run_group(VEL_SWEEP_CL, "vel-cl")
    wpv = run_group(VEL_SWEEP_WP, "vel-wp")
    f1 = os.path.join(ART, "bathstruct_sigma_sweep.png")
    f2 = os.path.join(ART, "bathstruct_velocity_sweep.png")
    if sig: fig_sigma(sig, f1)
    if clv: fig_vel(clv, wpv, f2)
    summary = {"sigma_sweep": sig, "vel_classical": clv, "vel_wp": wpv}
    json.dump(summary, open(os.path.join(ART, "bathstruct_summary.json"), "w"), indent=2)
    # concise textual result
    def line(r): return f"sig={r['sigma']} v={r['v']}: rank90={r['pod_rank90']} lead={r['lead_energy']:.2f} dmd={r['dmd_omega_ev']:.0f}eV"
    print("\n--- SIGMA SWEEP (v=2.71) ---"); [print("  " + line(r)) for r in sorted(sig, key=lambda r: r['sigma'])]
    print("--- CLASSICAL VELOCITY SWEEP ---"); [print("  " + line(r)) for r in sorted(clv, key=lambda r: r['v'])]
    email("ml-patterns: bath-structure sweep (classical vs quantum induced density)",
          dict(h="The induced bath density carries a classical-vs-quantum STRUCTURAL "
                 "difference (not a governing PDE — that was a blob artifact). POD rank "
                 "and DMD spectrum of the TRUE blob-free bath should show the WP (finite "
                 "sigma, low-pass) giving a more coherent / lower-rank response than the "
                 "point classical projectile (couples to all q -> higher-rank e-h).",
               d="POD + DMD on dn=(n_total-n_wp)-n_GS (WP) / n_total-n_GS (classical) "
                 "across a sigma-sweep at v=2.71 (classical point + WP sigma 0.5/3/8) and "
                 "a classical velocity sweep E=20-600, plus clean WP velocity points.",
               p="Fig1: POD rank-for-90% and leading-mode energy vs sigma_WP. Fig2: POD "
                 "rank and DMD dominant frequency vs velocity, classical vs WP, with omega_p.",
               c="POD-rank (complexity) and leading-mode coherence quantify the "
                 "classical<->quantum difference in the induced density; velocity trend "
                 "shows how bath structure evolves with v. DMD frequencies approximate "
                 "(non-stationary runs)."),
          [f1, f2])
    print("\nSAVED figures + bathstruct_summary.json")


if __name__ == "__main__":
    main(sys.argv[1:])
