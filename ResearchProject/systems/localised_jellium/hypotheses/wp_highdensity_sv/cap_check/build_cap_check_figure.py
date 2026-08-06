#!/usr/bin/env python3
"""
CAP validation replica — summary figure.

Reads the three cap_check variants (free WP, sigma_WP=0.5, in an EMPTY
35x35x85 periodicity(2) box with the production CAP) and produces the one-page
dashboard the user reviews before the six production runs are launched.

Run from the repo root:
    venv/bin/python3 ResearchProject/systems/localised_jellium/hypotheses/\
wp_highdensity_sv/cap_check/build_cap_check_figure.py

Plan: docs/plans/wavepacket-highdensity-sv-twin.md (section 6b)
"""
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = "/rds/user/skcb2/hpc-work/tddft/inq-tddft-research"
SRC = f"{REPO}/ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/cap_check/results"
OUT = f"{REPO}/ResearchProject/systems/localised_jellium/hypotheses/wp_highdensity_sv/cap_check"

LZ, LXY = 85.0, 35.0
CAP_L = 12.5
Z_CAP_IN = LZ / 2 - CAP_L          # +30.0
SIGMA = 0.5
T_TRANSVERSE = 4.12                # 6 sigma_d = L_xy


def load(path):
    lines = [l for l in open(path) if not l.startswith("#")]
    return [{k: float(v) for k, v in d.items()} for d in csv.DictReader(lines)]


def main():
    os.makedirs(OUT, exist_ok=True)
    runs = {}
    for name in ("cap_v2p0", "nocap_v2p0"):
        rs = f"{SRC}/{name}/raw/observables/wp_real_space_stats.csv"
        ms = f"{SRC}/{name}/raw/observables/wp_momentum_stats.csv"
        if os.path.exists(rs) and os.path.exists(ms):
            runs[name] = (load(rs), load(ms))
    if not runs:
        raise SystemExit(f"no cap_check results under {SRC}")

    col = {"cap_v2p0": "C3", "nocap_v2p0": "C0"}
    lab = {"cap_v2p0": r"CAP on ($|\eta|$=1 Ha, 12.5 Bohr/face)",
           "nocap_v2p0": "CAP off (control)"}

    fig, ax = plt.subplots(2, 2, figsize=(12.5, 8.5))

    for n, (rs, ms) in runs.items():
        t = np.array([r["time_au"] for r in rs])
        nn = np.array([r["norm_check"] for r in rs])
        zc = np.array([r["z_mean_circ"] for r in rs])
        sz = np.sqrt([r["sigma_z2"] for r in rs])
        pz = np.array([m["pz_mean"] for m in ms])
        ax[0, 0].plot(t, nn / nn[0], col[n], label=lab[n])
        ax[0, 1].plot(t, zc, col[n], label=lab[n])
        ax[1, 0].plot(t, pz, col[n], label=lab[n])
        ax[1, 1].plot(t, sz, col[n], label=lab[n])

    ax[0, 0].set(xlabel="t (a.u.)", ylabel="norm / norm(0)", yscale="log",
                 title="A. CAP absorbs; the control stays unitary")
    ax[0, 0].axhline(1, ls=":", c="grey")

    ax[0, 1].set(xlabel="t (a.u.)", ylabel=r"$\langle z\rangle_{\rm circ}$ (Bohr)",
                 title="B. Without a CAP the ORBITAL WRAPS the +z face")
    for y, c in ((LZ / 2, "k"), (Z_CAP_IN, "C3"), (-LZ / 2, "k"), (-Z_CAP_IN, "C3")):
        ax[0, 1].axhline(y, ls="--", c=c, lw=0.8)
    ax[0, 1].text(1, LZ / 2 + 1.5, "box face +42.5", fontsize=7)
    ax[0, 1].text(1, Z_CAP_IN + 1.5, "CAP inner edge +30", fontsize=7, color="C3")
    ax[0, 1].annotate("wrap: +41.4 -> -28.2", xy=(40, -28), xytext=(20, -38),
                      arrowprops=dict(arrowstyle="->", color="C0"), color="C0", fontsize=9)

    ax[1, 0].set(xlabel="t (a.u.)", ylabel=r"$\langle p_z\rangle$ (Bohr$^{-1}$)",
                 title="C. WARNING — the CAP alone decelerates the packet\n"
                       "(vacuum: no bath, no forces, so the true value is flat)")
    ax[1, 0].axhline(2.0, ls=":", c="grey")
    ax[1, 0].text(22, 1.90, "true free-particle value = 2.0", fontsize=8, color="grey")

    tt = np.linspace(0, 20, 200)
    ax[1, 1].plot(tt, np.sqrt(SIGMA**2 / 2 + tt**2 / (2 * SIGMA**2)), "k:",
                  label=r"free law $\sqrt{\sigma^2/2+t^2/2\sigma^2}$")
    ax[1, 1].set(xlabel="t (a.u.)", ylabel=r"$\sigma_z$ (Bohr)",
                 title=r"D. $\sigma_{WP}$=0.5 spreads at 1.41 Bohr/a.u.")
    ax[1, 1].axhline(LXY / 6, ls="--", c="green", lw=0.8)
    ax[1, 1].text(25, LXY / 6 + 0.6, r"$6\sigma_d=L_{xy}$ (image overlap)",
                  fontsize=8, color="green")

    for a in (ax[0, 0], ax[1, 0], ax[1, 1]):
        a.axvline(T_TRANSVERSE, color="green", ls="-.", lw=1.2)
    ax[0, 0].text(T_TRANSVERSE + 0.8, 0.32, "fit window ends\nt = 4.12 a.u.",
                  fontsize=8, color="green")

    for a in ax.flat:
        a.grid(alpha=0.3)
        a.legend(fontsize=8, loc="best")

    fig.suptitle(r"CAP validation replica — free WP ($\sigma_{WP}$=0.5, $k_0$=2, "
                 r"launch $z=-24$) in an EMPTY 35$\times$35$\times$85 periodicity(2) box",
                 fontsize=12)
    fig.tight_layout()
    path = f"{OUT}/cap_validation.png"
    fig.savefig(path, dpi=140)
    print("wrote", path)

    # numeric summary the notebook/handover quotes
    for n, (rs, ms) in runs.items():
        nn0 = rs[0]["norm_check"]
        print(f"\n{n}:")
        print(f"  final norm/norm0 = {rs[-1]['norm_check']/nn0:.4e}")
        print(f"  min <p_z>        = {min(m['pz_mean'] for m in ms):.4f} "
              f"(negative => reflection)")
        for tw in (4.0, 8.0, 12.0, 16.0):
            i = min(range(len(rs)), key=lambda j: abs(rs[j]["time_au"] - tw))
            print(f"  t={tw:5.1f}  norm={rs[i]['norm_check']/nn0:.4f}  "
                  f"<p_z>={ms[i]['pz_mean']:.4f}  z_circ={rs[i]['z_mean_circ']:8.2f}")


if __name__ == "__main__":
    main()
