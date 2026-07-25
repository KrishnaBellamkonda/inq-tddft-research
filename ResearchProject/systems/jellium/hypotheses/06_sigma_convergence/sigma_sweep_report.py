#!/usr/bin/env python3
"""Cumulative S(v) convergence report + optional email — σ-convergence sweep.

Scans every (σ, v) run directory, extracts the stopping power by the energy-gain
regression (Method A, primary) and the kinetic/speed-integration cross-check
(Method B), and plots ALL completed points — one colour per σ — over the SINGLE
point-charge Lindhard reference. Robust to partial sweeps: only completed runs
(run_completed = true) are plotted.

Usage (venv):
  python3 sigma_sweep_report.py                 # just (re)build the figure
  python3 sigma_sweep_report.py --email "0.15"  # build + email, tagged "σ=0.15 done"

The figure is written to figures/sv_convergence.png. With --email it is sent to
chiddukanna@gmail.com via inqview.email.send_run_email (Gmail SMTP App Password).
"""
from __future__ import annotations

import os
import sys
import csv
import argparse

import numpy as np

sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.analysis import lindhard_elf as E
from inqview.analysis.stopping_extract import load_track
from inqview.visualisation import style as S

JELLIUM = "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium"
HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")

RS = 5.69
KF = E.kF_from_rs(RS)
MASS = 1.0
TRANSIENT_FRAC = 0.20

# σ → (run dir, on-disk layout). Two contracts coexist:
#   "nested" — this sweep's runs: results/<vtag>/raw/observables/, stamped
#              "run_completed   = true" on clean exit.
#   "flat"   — the pre-existing σ=0.5 anchor set (run_sv_sigma0p5): observables
#              live directly under results/<vtag>/, and the summaries carry
#              "run_completed  = false" because those runs were halted at the
#              boundary-stop rule. Their data is nonetheless usable (the original
#              extraction notebook used it), so we gate on file presence instead.
# NB: these σ are projectile CHARGE standard deviations (charge std), not the
# wavepacket σ. The two differ by √2: a wavepacket of σ_WP has density/charge std
# σ_WP/√2 (validated 2026-06-23, density_std = σ_WP/√2; see CONTEXT "σ-convention").
# So the matched WP label is σ_WP = √2·σ_q (e.g. σ_q=0.354 ⇔ σ_WP=0.5).
SIGMAS = [
    (0.15, "run_classical_n162_L50_sv_sigma0p15", "nested"),
    (0.25, "run_classical_n162_L50_sv_sigma0p25", "nested"),
    (0.35, "run_classical_n162_L50_sv_sigma0p35", "nested"),
    (0.50, "run_sv_sigma0p5",                     "flat"),
    (3.00, "run_classical_n162_L50_sv_sigma3p0",  "nested"),
]
# v0 → result subdir tag. σ=0.5 has no v0p2; only σ=3 has the v1p0 peak-refinement
# point — runs that lack a given vtag are simply skipped (file-presence gated).
VTAGS = [(3.0, "v3p0"), (2.0, "v2p0"), (1.3, "v1p3"), (1.0, "v1p0"),
         (0.8, "v0p8"), (0.6, "v0p6"), (0.2, "v0p2")]

HA_TO_EV = 27.211386245988  # CODATA Hartree → eV


def _obs_path(run_dir, vtag, layout):
    base = os.path.join(JELLIUM, run_dir, "results", vtag)
    return base if layout == "flat" else os.path.join(base, "raw", "observables")


def _completed(run_dir, vtag, layout):
    obsdir = _obs_path(run_dir, vtag, layout)
    if layout == "flat":
        # anchor set: usable iff both data files exist (completion flag is false)
        return (os.path.exists(os.path.join(obsdir, "observables.csv"))
                and os.path.exists(os.path.join(obsdir, "electron_track.csv")))
    rs = os.path.join(JELLIUM, run_dir, "results", vtag, "run_summary.txt")
    return os.path.exists(rs) and "run_completed   = true" in open(rs).read()


def load_observables(obsdir):
    """(t_au, E_total) from observables.csv, deduped by step, time-sorted."""
    path = os.path.join(obsdir, "observables.csv")
    t, Et, seen = [], [], set()
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if row.get("energy_total") in (None, ""):
                continue
            k = row["step"]
            if k in seen:
                continue
            seen.add(k)
            t.append(float(row["time_au"]))
            Et.append(float(row["energy_total"]))
    t = np.array(t); Et = np.array(Et)
    o = np.argsort(t)
    return t[o], Et[o]


def extract(run_dir, vtag, layout):
    """Return dict(S_A, S_A_err, S_B, v, vlo, vhi, r2, n) or None."""
    from scipy.stats import linregress
    obsdir = _obs_path(run_dir, vtag, layout)
    tr = load_track(os.path.join(obsdir, "electron_track.csv"), mass=MASS, axis="z")
    t_o, Et = load_observables(obsdir)
    if t_o.size < 4 or tr.t.size < 4:
        return None
    t_cut = tr.t.min() + TRANSIENT_FRAC * (tr.t.max() - tr.t.min())
    # Method A
    dE = Et - Et[0]
    m = t_o >= t_cut
    if m.sum() < 3:
        return None
    s_at = np.interp(t_o[m], tr.t, tr.s)
    fitA = linregress(s_at, dE[m])
    # Method B (kinetic / speed-integration)
    mk = tr.t >= t_cut
    path = np.trapezoid(tr.v[mk], tr.t[mk])
    dKE = tr.ke[mk][0] - tr.ke[mk][-1]
    S_B = dKE / path if path > 0 else float("nan")
    vv = tr.v[mk]
    return dict(S_A=fitA.slope, S_A_err=fitA.stderr, S_B=S_B,
                v=float(vv.mean()), vlo=float(vv.min()), vhi=float(vv.max()),
                r2=float(fitA.rvalue ** 2), n=int(m.sum()))


def collect():
    rows = []  # (sigma, run_dir, vtag, v0, result|None)
    for sigma, rd, layout in SIGMAS:
        for v0, vtag in VTAGS:
            if _completed(rd, vtag, layout):
                try:
                    r = extract(rd, vtag, layout)
                except Exception as exc:  # noqa: BLE001
                    r = None
                    print(f"  [warn] {rd}/{vtag}: {exc}")
                rows.append((sigma, rd, vtag, v0, r))
            else:
                rows.append((sigma, rd, vtag, v0, None))
    return rows


def _energy_ev(v):
    """Projectile kinetic energy E = ½ m v² in eV (electron, m=1 a.u.)."""
    return 0.5 * np.asarray(v, dtype=float) ** 2 * HA_TO_EV


def _reference_curve():
    """(vg, slr, v_pk) for the point-charge Lindhard reference. 90 pts: smooth
    curve, peak to ~0.03 a.u.; each point is a 4k×4k integral so keep it modest.
    Computed once and reused across the velocity and energy figures."""
    vg = np.linspace(0.12, 3.25, 90)
    slr = np.array([E.stopping_power_point(float(v), KF) for v in vg])
    return vg, slr, float(vg[int(np.argmax(slr))])


def _draw(ax, rows, ref, *, mode):
    """Render the S(v) convergence onto ax. mode='v' (velocity, a.u.) or
    'E' (projectile kinetic energy, eV, log x). Points sit at NOMINAL v0;
    the only error bar is the vertical linregress stderr on S. Returns
    (n_done, peak_x, kF_x) for axis/annotation use."""
    xf = (lambda v: v) if mode == "v" else _energy_ev
    vg, slr, v_pk = ref
    ax.plot(xf(vg), slr, "-", color="k", lw=1.6, zorder=2,
            label="Lindhard (point charge)")

    cmap = S.plt.get_cmap(S.cmap_for("sequential"))
    n_done = 0
    for i, (sigma, _rd, _layout) in enumerate(SIGMAS):
        # invert: smallest σ darkest (nearest the black point-charge curve)
        col = cmap(0.82 - 0.64 * i / max(len(SIGMAS) - 1, 1))
        # x = NOMINAL v0; yerr = linregress stderr; no horizontal error bar
        pts = [(v0, r["S_A"], r["S_A_err"])
               for (s, _rd, _vt, v0, r) in rows if s == sigma and r is not None]
        if not pts:
            continue
        pts.sort()
        for v0, sA, sErr in pts:
            n_done += 1
            ax.errorbar(xf(v0), sA, yerr=sErr,
                        fmt="o", color=col, ms=5.5, capsize=2.0, elinewidth=0.8,
                        ecolor=col, mec="k", mew=0.4, zorder=5)
        ax.plot([], [], "o", color=col, mec="k", mew=0.4, label=rf"$\sigma_q$={sigma}")
    return n_done, xf(v_pk), xf(KF)


def _annotate_lines(ax, peak_x, kF_x):
    """Two annotated vertical dashed lines: Lindhard peak and k_F."""
    for x, txt in ((peak_x, "Lindhard peak"), (kF_x, "k$_F$")):
        ax.axvline(x, ls="--", color="gray", lw=0.9, zorder=1)
        ax.annotate(txt, xy=(x, 0.97), xycoords=("data", "axes fraction"),
                    rotation=90, va="top", ha="right", fontsize=6,
                    color="gray", xytext=(-2, 0), textcoords="offset points")


def build_figure(rows):
    """Velocity-axis (primary) and energy-axis (companion) S(v) figures."""
    S.apply_theme()
    os.makedirs(FIGDIR, exist_ok=True)
    ref = _reference_curve()  # computed once, shared by both figures

    # --- velocity figure (primary) ---
    fig, ax = S.figure_one_col()
    n_done, peak_x, kF_x = _draw(ax, rows, ref, mode="v")
    _annotate_lines(ax, peak_x, kF_x)
    ax.set_xlabel("v  (a.u.)")
    ax.set_ylabel("S(v)  (Ha/Bohr)")
    ax.set_ylim(bottom=0)
    ax.set_title(f"σ-convergence of S(v) → point-charge Lindhard  "
                 f"(r$_s$={RS}, k$_F$={KF:.3f})", fontsize=7.5)
    ax.legend(fontsize=6, loc="upper right", framealpha=0.9,
              title=r"$\sigma_q$ = charge std  ($\sigma_\mathrm{WP}{=}\sqrt{2}\,\sigma_q$)",
              title_fontsize=5.5)
    out = os.path.join(FIGDIR, "sv_convergence.png")
    fig.savefig(out, dpi=200)
    S.plt.close(fig)
    print(f"wrote {out}  ({n_done} completed points)")

    # --- energy figure (companion, eV, log x) ---
    fig2, ax2 = S.figure_one_col()
    _, peak_e, kF_e = _draw(ax2, rows, ref, mode="E")
    _annotate_lines(ax2, peak_e, kF_e)
    ax2.set_xscale("log")
    ax2.set_xlabel("projectile kinetic energy  E = ½v²  (eV)")
    ax2.set_ylabel("S  (Ha/Bohr)")
    ax2.set_ylim(bottom=0)
    ax2.set_title(f"σ-convergence of S(E) → point-charge Lindhard  "
                  f"(r$_s$={RS}, k$_F$={KF:.3f})", fontsize=7.5)
    ax2.legend(fontsize=6, loc="upper left", framealpha=0.9,
               title=r"$\sigma_q$ = charge std  ($\sigma_\mathrm{WP}{=}\sqrt{2}\,\sigma_q$)",
               title_fontsize=5.5)
    out_e = os.path.join(FIGDIR, "sv_convergence_energy.png")
    fig2.savefig(out_e, dpi=200)
    S.plt.close(fig2)
    print(f"wrote {out_e}  ({n_done} completed points)")
    return out, out_e, n_done


def text_table(rows):
    lines = [f"{'sigma':>6} {'v0':>5} {'v_win':>6} {'S_A':>9} {'+/-':>8} "
             f"{'S_B':>9} {'A/B':>6} {'S_LR_pt':>9} {'sim/LR':>7} {'R2':>6}"]
    for sigma, _rd, _vt, v0, r in rows:
        if r is None:
            lines.append(f"{sigma:>6} {v0:>5.1f}  (pending)")
            continue
        slr = E.stopping_power_point(r["v"], KF)
        ab = r["S_A"] / r["S_B"] if r["S_B"] else float("nan")
        lines.append(f"{sigma:>6} {v0:>5.1f} {r['v']:>6.2f} {r['S_A']:>9.4f} "
                     f"{r['S_A_err']:>8.4f} {r['S_B']:>9.4f} {ab:>6.2f} "
                     f"{slr:>9.4f} {r['S_A']/slr:>7.2f} {r['r2']:>6.3f}")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", default=None,
                    help="σ label just completed (e.g. '0.15'); sends the figure")
    args = ap.parse_args()

    rows = collect()
    out, out_e, n_done = build_figure(rows)
    table = text_table(rows)
    print("\n" + table)

    if args.email is not None:
        from inqview.email import send_run_email
        n_total = len(SIGMAS) * len(VTAGS)
        subj = (f"[sigma-convergence] σ={args.email} ladder complete "
                f"({n_done} runs) — S(v) update")
        body = (
            f"σ = {args.email} ladder finished.\n\n"
            f"Cumulative S(v) convergence toward the single point-charge Lindhard "
            f"reference (r_s={RS}, k_F={KF:.4f}). {n_done} completed runs across "
            f"all σ.\n\n"
            f"Points sit at the NOMINAL launch velocity v0; the vertical bar is the "
            f"linear-regression stderr on S. Method A = energy-gain regression "
            f"(primary); Method B = kinetic/speed-integration cross-check; A/B "
            f"within 10% is a good sign.\n\n"
            f"{table}\n\n"
            f"Figures attached:\n"
            f"  sv_convergence.png         — S vs velocity (a.u.)\n"
            f"  sv_convergence_energy.png  — S vs projectile KE (eV, log x)\n"
        )
        try:
            msg_id = send_run_email(subj, body, attachments=[out, out_e],
                                    to="chiddukanna@gmail.com")
            print(f"emailed: {msg_id}")
        except Exception as exc:  # noqa: BLE001
            print(f"[EMAIL FAILED] {exc}")
            sys.exit(3)


if __name__ == "__main__":
    main()
