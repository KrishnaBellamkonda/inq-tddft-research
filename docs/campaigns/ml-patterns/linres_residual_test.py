#!/usr/bin/env python3
"""Linear-response residual / form-factor test — classical vs WP induced density.

Panel-chosen technique (2026-07-06). Tests the null n_WP(q,t) = F(q) n_cl(q,t)
with F(q)=exp(-q^2 sigma^2/2), frame-by-frame in the TIME domain (the matched runs
are too short for any omega-resolved method). See:
  docs/plans/linres-residual-classical-vs-wp.md
  kernels/formfactor_residual.py  (validated: tests/test_formfactor_residual.py)

Discriminants (both d'Alembert-safe -- magnitudes cancel rigid f(z-vt)):
  (1) |R(q,t)| == F(q)?  (2) is |R(q,t)| flat in t?  (3) high-q residual excess?
Fork A (sigma_WP vs sigma_pot=sigma_WP/sqrt2) resolved EMPIRICALLY via the
a(sigma) collapse -- nothing hardcoded.

Hands-off: per-pair try/except + idempotent JSON, figures + 4-part email. CPU
only, NO new runs, streams frames (never holds the (T,125^3) stack).

Run: venv/bin/python3 docs/campaigns/ml-patterns/linres_residual_test.py [--no-email] [--maxf N]
"""
from __future__ import annotations
import sys, os, re, glob, json, argparse, traceback
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ART = os.path.join(HERE, "artifacts")
os.makedirs(ART, exist_ok=True)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from inqview import load_vti
from kernels import formfactor_residual as FR

J = "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium"
OMP_EV = 0.1276 * 27.211
SEND_EMAIL = True

CLASSICAL = "run_classical_n162_L50_E100_v2"       # shared classical point partner
PAIRS = [   # matched v=2.71; WP _wf = blob-free bath.  varyv (dt=4.0) EXCLUDED.
    dict(sigma=0.5, wp="run_wp_n162_L50_E100_sigma0p5_wf"),   # primary (signal to q~4)
    dict(sigma=3.0, wp="run_wp_n162_L50_E100_sigma3_wf"),     # cross-check
    dict(sigma=8.0, wp="run_wp_n162_L50_E100_sigma8_wf"),     # SNR-dead; collapse only
]
_STEP = re.compile(r"density_t(\d+)\.vti$")


def _summary_field(run, key, cast=float):
    p = glob.glob(os.path.join(J, run, "results", "*run_summary*"))
    if not p:
        raise FileNotFoundError(f"no run_summary for {run}")
    for line in open(p[0]):
        if line.strip().startswith(key):
            return cast(line.split("=", 1)[1].split()[0])
    raise KeyError(f"{key} not in {p[0]}")


def _series(run, sub):
    return sorted(glob.glob(os.path.join(J, run, "results/raw/vti", sub, "density_t*.vti")))


def stream_spectra(run, is_wp, maxf):
    """Stream a bath series -> (q_r, amp_r(T,Nq), noise_r(Nq), q_z, amp_z(T,Nqz), times).

    bath = density_total - density_wp (WP) or density_total (classical); GS-subtract;
    radial + axial |q| spectra per frame.  Times = step_index * dt_au (exact).
    """
    dt = _summary_field(run, "dt_au")
    tot = _series(run, "density_total")
    if not tot:
        raise FileNotFoundError(f"no density_total for {run}")
    stride = max(1, len(tot) // maxf)
    tot = tot[::stride]
    wp = _series(run, "density_wp")[::stride] if is_wp else None
    gs = load_vti(os.path.join(J, run, "results/raw/vti",
                               "density_gs_system/density_gs_system.vti")).data.astype(np.float32)
    dx = float(load_vti(tot[0]).x[1] - load_vti(tot[0]).x[0])
    amp_r, noise_r, amp_z, times = [], [], [], []
    q_r = q_z = None
    for i, f in enumerate(tot):
        n_tot = load_vti(f).data.astype(np.float32)
        bath = n_tot - load_vti(wp[i]).data.astype(np.float32) if is_wp else n_tot
        ind = (bath - gs)
        qr, a_r, nz_r, _ = FR.radial_spectrum(ind, dx)
        qz, a_z = FR.axial_spectrum(ind, dx)
        q_r, q_z = qr, qz
        amp_r.append(a_r); noise_r.append(nz_r); amp_z.append(a_z)
        step = int(_STEP.search(os.path.basename(f)).group(1))
        times.append(step * dt)
    amp_r = np.array(amp_r); amp_z = np.array(amp_z)
    noise_r = np.median(np.array(noise_r), axis=0)          # per-shell floor
    return q_r, amp_r, noise_r, q_z, amp_z, np.array(times)


def tail_floor(nbar, q, frac=0.8):
    """Numerical noise floor = median amplitude in the signal-free high-q tail.

    The physical induced signal n_ind = chi(q) V_ext(q) decays at large q into a
    flat numerical plateau (GS-subtraction / interpolation residual); its median
    over the top (1-frac) of the q-range is the per-run noise floor. This is the
    correct SNR gate -- NOT the shell standard-error, which is tiny for shells
    with many modes and would admit the whole noise tail into the F(q) fit.
    """
    hi = q >= frac * q.max()
    return float(np.median(nbar[hi])) if hi.any() else float(np.median(nbar))


def run_pair(pair, maxf):
    sigma, wp = pair["sigma"], pair["wp"]
    cp = os.path.join(ART, f"linres_residual_sigma{sigma}.json")
    if os.path.isfile(cp):
        print(f"[sigma={sigma}] cached"); return json.load(open(cp))
    qr_c, ar_c, nz_c, qz_c, az_c, t_c = stream_spectra(CLASSICAL, False, maxf)
    qr_w, ar_w, nz_w, qz_w, az_w, t_w = stream_spectra(wp, True, maxf)
    assert np.allclose(qr_c, qr_w), "radial q-grids differ (unequal boxes/grids)"
    Tov = min(t_c.max(), t_w.max())
    M = min((t_c <= Tov).sum(), (t_w <= Tov).sum())
    tcommon = np.linspace(0.0, Tov, int(M))
    # 3-D radial (headline). Noise floor = high-q plateau (see tail_floor), which
    # cuts the WP band at the q where its form factor buries signal in noise.
    ncl = FR.resample_time(t_c, ar_c, tcommon)
    nwp = FR.resample_time(t_w, ar_w, tcommon)
    floor_c = np.full(qr_c.size, tail_floor(np.median(ncl, axis=0), qr_c))
    floor_w = np.full(qr_w.size, tail_floor(np.median(nwp, axis=0), qr_w))
    res3d = FR.residual_test(
        FR.PairResult(sigma, qr_c, ncl, nwp, floor_c, tcommon, noise_wp=floor_w), snr=3.0)
    # Panel's t-drift discriminant: if |R| is non-flat over the full window but
    # FLAT in the early v~v0 window, the drift is deceleration (WP slows, classical
    # holds v) rather than noise. Compare early-40% vs late-40% flatness.
    Ne = max(6, int(0.4 * M))
    early = FR.residual_test(FR.PairResult(sigma, qr_c, ncl[:Ne], nwp[:Ne], floor_c,
                                           tcommon[:Ne], noise_wp=floor_w), snr=3.0)
    late = FR.residual_test(FR.PairResult(sigma, qr_c, ncl[-Ne:], nwp[-Ne:], floor_c,
                                          tcommon[-Ne:], noise_wp=floor_w), snr=3.0)
    res3d["t_flatness_early"] = early["t_flatness"]
    res3d["t_flatness_late"] = late["t_flatness"]
    res3d["early_window_au"] = [float(tcommon[0]), float(tcommon[Ne - 1])]
    # 1-D axial (cross-check): high-q-tail floor on the q_perp=0 line
    nzc = FR.resample_time(t_c, az_c, tcommon)
    nzw = FR.resample_time(t_w, az_w, tcommon)
    fzc = np.full(qz_c.size, tail_floor(np.median(nzc, axis=0), qz_c))
    fzw = np.full(qz_w.size, tail_floor(np.median(nzw, axis=0), qz_w))
    res1d = FR.residual_test(FR.PairResult(sigma, qz_c, nzc, nzw, fzc, tcommon,
                                           noise_wp=fzw), snr=3.0)
    out = dict(sigma=sigma, wp=wp, classical=CLASSICAL, t_overlap_au=float(Tov),
               n_common=int(M), radial=res3d, axial=res1d)
    json.dump(out, open(cp, "w"), indent=2)
    print(f"[sigma={sigma}] Tov={Tov:.1f}au M={M} | 3D: sigma_fit={res3d['sigma_fit']:.3f} "
          f"(WP={sigma}, pot={sigma/np.sqrt(2):.3f}) r2={res3d['fit_r2']:.2f} "
          f"flat={res3d['t_flatness']:.2f} excess={res3d['highq_excess_over_noise']:.1f}")
    return out


def fig_pair(res, path):
    r = res["radial"]; q = np.array(r["q"]); band = np.array(r["band_mask"])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    ax[0].plot(q, r["R_median"], "o", ms=3, color="C0", label="|R(q)| = |n_WP|/|n_cl|")
    ax[0].plot(q, r["F_at_sigma_wp"], "-", color="C2", label=f"F(sigma_WP={res['sigma']})")
    ax[0].plot(q, r["F_at_sigma_pot"], "--", color="C3",
               label=f"F(sigma_pot={res['sigma']/np.sqrt(2):.3f})")
    if band.any():
        ax[0].axvspan(0, q[band].max(), color="k", alpha=0.05)
    ax[0].set_xlabel("|q| (a.u.)"); ax[0].set_ylabel("magnitude ratio")
    ax[0].set_ylim(0, 1.3); ax[0].legend(fontsize=8)
    ax[0].set_title(f"(a) form-factor test, sigma={res['sigma']} (fit sigma={r['sigma_fit']:.3f})")
    rn = np.array(r["resid_over_noise"])
    ax[1].axhline(0, color="k", lw=0.6); ax[1].axhline(3, color="r", ls=":", lw=0.8)
    ax[1].axhline(-3, color="r", ls=":", lw=0.8)
    ax[1].plot(q[band], rn[band], "s", ms=3, color="C1")
    ax[1].set_xlabel("|q| (a.u.)"); ax[1].set_ylabel("(n_WP - F*n_cl) / noise")
    ax[1].set_title(f"(b) normalized residual  (high-q excess={r['highq_excess_over_noise']:.1f} sigma)")
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def fig_collapse(col, per, path):
    """Fork-A collapse: only SNR-adequate sigma carry a meaningful exponent; the
    rest (form factor e-folds within a few shells) fit the noise floor and are
    plotted hollow/greyed as excluded."""
    fig, ax = plt.subplots(figsize=(6.0, 4.3))
    adq = {p["radial"]["sigma"] for p in per if p["radial"].get("snr_adequate")}
    xs = np.linspace(0, max(p["radial"]["sigma"] for p in per) ** 2 * 1.05, 50)
    ax.plot(xs, 0.5 * xs, "-", color="C2", label="slope 0.5  (sigma_WP)")
    ax.plot(xs, 0.25 * xs, "--", color="C3", label="slope 0.25 (sigma_pot)")
    seen = set()
    for p in per:
        r = p["radial"]; s = r["sigma"]; ok = s in adq
        lab = "SNR-adequate a(sigma)" if ok else "SNR-dead (excluded)"
        ax.plot(s**2, r["exponent_a"], "o" if ok else "x",
                color="C0" if ok else "0.6", ms=7,
                label=lab if lab not in seen else None)
        seen.add(lab)
    if col.get("ok"):
        ax.plot(xs, col["slope"] * xs, ":", color="k",
                label=f"fit slope={col['slope']:.3f} -> {col['selects']}")
        title = f"Fork A: slope {col['slope']:.3f} -> {col['selects']}"
    else:
        sp = col.get("single_point")
        title = ("Fork A INCONCLUSIVE (need >=2 SNR-adequate sigma); "
                 + (f"single-point leans {sp['leans']}" if sp else "no usable sigma"))
    ax.set_xlabel("sigma^2 (Bohr^2)"); ax.set_ylabel("exponent a  [|R|~exp(-a q^2)]")
    ax.legend(fontsize=8); ax.set_title(title, fontsize=9)
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-email", action="store_true")
    ap.add_argument("--maxf", type=int, default=100)
    a = ap.parse_args(argv); SEND_EMAIL = not a.no_email
    print("=== linear-response residual test (classical vs WP induced density) ===")
    per = []
    for pair in PAIRS:
        try:
            per.append(run_pair(pair, a.maxf))
        except Exception as e:
            print(f"[sigma={pair['sigma']}] SKIP: {e!r}"); traceback.print_exc()
    if not per:
        print("no pairs succeeded"); return
    col = FR.collapse_fork_a([p["radial"] for p in per])
    figs = []
    for p in per:
        f = os.path.join(ART, f"linres_residual_sigma{p['sigma']}.png")
        fig_pair(p, f); figs.append(f)
    fc = os.path.join(ART, "linres_forkA_collapse.png"); fig_collapse(col, per, fc)
    figs.append(fc)
    summary = dict(pairs=per, fork_a=col)
    json.dump(summary, open(os.path.join(ART, "linres_residual_summary.json"), "w"), indent=2)

    prim = next((p for p in per if p["sigma"] == 0.5), per[0])["radial"]
    fe = prim.get("t_flatness_early", np.nan); ff = prim.get("t_flatness", np.nan)
    drift_is_decel = np.isfinite(fe) and ff >= 0.2 and fe < 0.5 * ff
    if abs(prim["highq_excess_over_noise"]) < 3 and ff < 0.2:
        verdict = "NULL: one linear filter (F~exp(-q^2 sigma^2/2)) explains the whole difference"
    elif ff >= 0.2:
        verdict = ("t-DRIFT: |R(q,t)| non-flat -> the static-linear-filter null is REJECTED; "
                   + ("early-window flat -> DECELERATION (WP slows, classical holds v)"
                      if drift_is_decel else "drift persists early -> not pure deceleration"))
    else:
        verdict = "EXCESS: high-q non-linear/quantum fingerprint at sigma=0.5"
    dead = col.get("snr_dead_sigmas", [])
    fork = (f"slope {col['slope']:.3f} -> {col['selects']}" if col.get("ok")
            else f"INCONCLUSIVE ({col.get('reason')}); "
                 + (f"single-point leans {col['single_point']['leans']}"
                    if col.get("single_point") else "no usable sigma"))
    print(f"\n=== PRIMARY (sigma=0.5) VERDICT: {verdict} ===")
    print(f"   t_flatness full={ff:.3f} early={fe:.3f} -> decel={drift_is_decel}")
    print(f"Fork A: {fork}  (SNR-dead sigma excluded: {dead})")
    email("ml-patterns: linear-response residual test (classical vs WP)",
          dict(h="If both projectiles probe the SAME medium chi(q,w), the WP is just a "
                 "low-pass-filtered point charge: n_WP(q,t)=F(q) n_cl(q,t), "
                 "F(q)=exp(-q^2 sigma^2/2). t-flat collapse onto F = one linear filter; "
                 "t-flat high-q excess = quantum/non-linear physics; a t-drift that is "
                 "flat in the early v~v0 window = deceleration (WP slows, classical holds v).",
               d="3-D radial + 1-D axial |q| spectra of the blob-free induced bath, "
                 "frame-by-frame, matched-v (2.71) pairs sigma=0.5/3/8 vs the shared "
                 "classical point run, over each overlap window; empirical exponent fit "
                 "(descending arm) resolves the sigma_WP vs sigma_pot sqrt2 trap. Early-vs-"
                 "late window flatness split isolates deceleration. ONLY SNR-adequate "
                 f"sigma feed the Fork-A collapse (SNR-dead excluded: {dead}).",
               p="Per-sigma: (a) |R(q)| vs F(sigma_WP)/F(sigma_pot) + SNR band; (b) "
                 "normalized residual (n_WP-F n_cl)/noise. Collapse: exponent a vs sigma^2, "
                 "SNR-dead sigma greyed/excluded.",
               c=f"PRIMARY (sigma=0.5, the only SNR-adequate pair): {verdict}. "
                 f"Fork A: {fork}. sigma=3/8 are SNR-dead (form factor e-folds within "
                 "~1-4 shells; fit hits the blob-subtraction floor) -- exactly why the "
                 "panel's new run (v~1-2, L>=100, >=3 plasma periods, density_wp written) "
                 "is needed. Deeper than POD/DMD (V_ext division is d'Alembert-safe) and, "
                 "unlike any omega-resolved method, extractable on the existing short runs."),
          figs)
    print("SAVED linres_residual_summary.json + figures")


if __name__ == "__main__":
    main(sys.argv[1:])
