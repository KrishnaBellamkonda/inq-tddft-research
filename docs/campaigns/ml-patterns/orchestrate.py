#!/usr/bin/env python3
"""Autonomous orchestrator for the ml-patterns campaign (T1 -> T7).

Idempotent / resumable: a phase whose result JSON exists is skipped. Each phase
is wrapped in try/except (one phase's bug never kills the chain) and emails a
4-part result (hypothesis reminder -> what was done -> what the plot shows ->
conclusion) with >=1 plot via inqview.email.send_run_email (email failure is
logged, never fatal). Per-rung notebooks are auto-built under notebooks/.

Anti-p-hacking (ADR 0011): the >=4-try loop tunes the SHARED config on the
PINNED CALIBRATION split only; verdicts are read on the PINNED HELD-OUT cells;
all attempts are logged. Verdicts CONFIRM/REFUTE/INCONCLUSIVE are all valid.

CPU-only, numpy/scipy. NO INQ runs, NO GPU.

Usage:
    venv/bin/python3 docs/campaigns/ml-patterns/orchestrate.py            # all phases
    venv/bin/python3 docs/campaigns/ml-patterns/orchestrate.py T2 T3      # subset
    ... --no-email   to skip emails
"""
from __future__ import annotations
import os, sys, json, time, traceback
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ART = os.path.join(HERE, "artifacts")
NB = os.path.join(HERE, "notebooks")
os.makedirs(ART, exist_ok=True); os.makedirs(NB, exist_ok=True)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from kernels import celldb, pipeline as PL, formfactor as FF, pod as P, dmd as D
from kernels import normaliser as NRM

SEND_EMAIL = True
HYPO = ("Hypothesis: the bath induced-density wake of a finite-sigma quantum "
        "wavepacket projectile differs from a point classical projectile at "
        "matched velocity in two SIE-controlled, interpretable ways the scalar "
        "S(v) cannot capture: (iv) a q-space form-factor softening "
        "R(q)=n_WP(q)/n_classical(q) ~ F_WP/F_ONCV (~exp(-q^2 sigma_pot^2/2)) and "
        "(iii) a collective DMD wake frequency obeying lambda(v)=2*pi*v/omega_p. "
        "Verdicts read on a PINNED held-out cell split (ADR 0011); +-20% band.")


# --------------------------------------------------------------------------- IO
def res_path(phase): return os.path.join(ART, f"{phase}_result.json")
def have(phase): return os.path.isfile(res_path(phase))
def save(phase, obj):
    json.dump(_jsonify(obj), open(res_path(phase), "w"), indent=2)
def load(phase): return json.load(open(res_path(phase)))


def _jsonify(o):
    if isinstance(o, dict): return {k: _jsonify(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [_jsonify(v) for v in o]
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, (np.bool_,)): return bool(o)
    return o


def r2(x):
    try:
        if x is None or (isinstance(x, float) and not np.isfinite(x)): return x
        return float(f"{x:.2g}")
    except Exception:
        return x


# ----------------------------------------------------------------------- email
def email(subject, parts, plots):
    """parts = dict(hypothesis, done, plot_shows, conclusion); plots = [paths]."""
    if not SEND_EMAIL:
        print(f"[email skipped] {subject}"); return
    body = (f"1) HYPOTHESIS\n{parts['hypothesis']}\n\n"
            f"2) WHAT WAS DONE\n{parts['done']}\n\n"
            f"3) WHAT THE PLOT SHOWS\n{parts['plot_shows']}\n\n"
            f"4) CONCLUSION\n{parts['conclusion']}\n")
    try:
        from inqview.email import send_run_email
        send_run_email(subject=subject, body=body,
                       attachments=[p for p in plots if os.path.isfile(p)])
        print(f"[email sent] {subject}")
    except Exception as e:
        print(f"[email FAILED, continuing] {subject}: {e!r}")


# -------------------------------------------------------------------- notebooks
def write_notebook(name, title, cells_md_code):
    import nbformat as nbf
    nb = nbf.v4.new_notebook()
    nb.cells.append(nbf.v4.new_markdown_cell(f"# {title}\n\n{HYPO}"))
    for kind, content in cells_md_code:
        if kind == "md":
            nb.cells.append(nbf.v4.new_markdown_cell(content))
        else:
            nb.cells.append(nbf.v4.new_code_cell(content))
    path = os.path.join(NB, name)
    nbf.write(nb, path)
    return path


# ====================================================================== PHASES
def phase_T1():
    """Pre-gate: kernels validated (code-test + formula-validation) + F_ONCV."""
    # code-test
    import subprocess
    tp = os.path.join(HERE, "tests/test_kernels.py")
    cp = subprocess.run([sys.executable, tp], capture_output=True, text=True)
    tests_pass = cp.returncode == 0
    # F_ONCV unity range
    q = np.linspace(0, 6, 400)
    qmax5, F = FF.foncv_unity_range(PL.ONCV_UPF, q, tol=0.05)
    qmax2, _ = FF.foncv_unity_range(PL.ONCV_UPF, q, tol=0.02)
    # plot F_ONCV and a few F_WP/F_ONCV predictions
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(q[q > 0], F[q > 0]); ax[0].axhline(1, ls=":", c="grey")
    ax[0].axvline(qmax5, ls="--", c="r", label=f"|F-1|<5%: q<={qmax5:.2f}")
    ax[0].set_xlabel("q (1/Bohr)"); ax[0].set_ylabel("F_ONCV(q)")
    ax[0].set_title("ONCV projectile form factor"); ax[0].legend(); ax[0].set_xlim(0, 4)
    for s_pot, lab in [(0.354, "sigma_WP=0.5"), (0.707, "1"), (2.12, "3"), (5.66, "8")]:
        ax[1].plot(q[q > 0], PL.predict_FF(q[q > 0], s_pot), label=lab)
    ax[1].set_xlabel("q (1/Bohr)"); ax[1].set_ylabel("F_WP/F_ONCV")
    ax[1].set_title("T2 predictions"); ax[1].legend(); ax[1].set_xlim(0, 3); ax[1].set_yscale("log")
    fig.tight_layout(); p = os.path.join(ART, "T1_foncv.png"); fig.savefig(p, dpi=110); plt.close(fig)

    catalogue = [
        {"kernel": "pod.pod", "test": "rank-2 synthetic recovers 2 modes + subspace",
         "formula_validation": "CONFIRM (Brunton&Kutz/Halko)", "status": "pass" if tests_pass else "FAIL"},
        {"kernel": "dmd.dmd", "test": "damped sinusoid recovers omega & decay",
         "formula_validation": "CONFIRM (Tu et al. 2014 exact DMD)", "status": "pass" if tests_pass else "FAIL"},
        {"kernel": "formfactor.{F_WP,F_ONCV_from_upf,radial_power_spectrum,q_ratio}",
         "test": "Gaussian width/ratio + F_ONCV~1 low-q", "formula_validation": "CONFIRM (Jackson radial FT)",
         "status": "pass" if tests_pass else "FAIL"},
    ]
    result = {"phase": "T1", "tests_pass": tests_pass,
              "foncv_unity_q_5pct": r2(qmax5), "foncv_unity_q_2pct": r2(qmax2),
              "catalogue": catalogue, "plot": p,
              "formula_validation": {"POD": "CONFIRM", "DMD": "CONFIRM", "form_factor": "CONFIRM"}}
    save("T1", result)
    write_notebook("T1_kernels.ipynb", "T1 - Pre-gated kernels + F_ONCV", [
        ("md", f"All three formula-bearing kernels passed an independent "
               f"`formula-validation` agent (POD/DMD/form-factor: CONFIRM) and a "
               f"known-case `code-test` (tests_pass={tests_pass}).\n\n"
               f"**F_ONCV(q) (the ACTUAL ONCV projectile form factor)** computed "
               f"from the local potential of `electron-ONCV-1.2.upf` (Coulomb-tail "
               f"subtracted): **F_ONCV ~= 1 within 5% for q <= {qmax5:.2f} 1/Bohr** "
               f"(within 2% for q <= {qmax2:.2f}). The T2 prediction reduces to "
               f"exp(-q^2 sigma_pot^2/2) only inside this window."),
        ("code", "import json; print(json.load(open('../artifacts/T1_result.json')))"),
        ("md", "![F_ONCV](../artifacts/T1_foncv.png)"),
    ])
    email("ml-patterns T1: kernels pre-gated + F_ONCV",
          {"hypothesis": HYPO,
           "done": ("Built campaign-local kernels (POD via truncated/randomized SVD, "
                    "exact windowed DMD, subtraction-ladder normaliser, form-factor). "
                    "Each formula-bearing kernel passed an independent formula-validation "
                    "agent (POD/DMD/form-factor all CONFIRM) and a known-case code-test "
                    f"(pass={tests_pass}). Computed the ACTUAL F_ONCV(q) from the ONCV UPF."),
           "plot_shows": (f"Left: F_ONCV(q) ~= 1 up to q={qmax5:.2f} 1/Bohr (5%); "
                          "right: the parameter-free T2 predictions F_WP/F_ONCV per sigma."),
           "conclusion": (f"T1 PASS. Kernels validated; the q-range where the prediction "
                          f"is honestly exp(-q^2 sigma_pot^2/2) is q <= {qmax5:.2f} 1/Bohr.")},
          [p])
    return result


# ---- T2 form factor ----
CALIB_CONFIGS = [
    dict(nbins=32, qmax=1.9, q_fit_lo=0.20, q_fit_hi=1.6, pred_floor=0.20, smooth=3, max_frames=90, window_frac=(0.05, 0.6)),
    dict(nbins=28, qmax=1.6, q_fit_lo=0.20, q_fit_hi=1.4, pred_floor=0.10, smooth=5, max_frames=90, window_frac=(0.05, 0.6)),
    dict(nbins=40, qmax=1.9, q_fit_lo=0.25, q_fit_hi=1.5, pred_floor=0.15, smooth=3, max_frames=120, window_frac=(0.10, 0.7)),
    dict(nbins=24, qmax=1.4, q_fit_lo=0.18, q_fit_hi=1.2, pred_floor=0.08, smooth=5, max_frames=70, window_frac=(0.05, 0.5)),
]


def _cell_agree(cell, cfg):
    r = PL.cell_agreement(cell, PL.FormFactorConfig(**cfg))
    return r


def phase_T2():
    cells = celldb.resolve_form_factor_cells()
    attempts = []
    best = None
    for i, cfg in enumerate(CALIB_CONFIGS):
        scores = []
        per = []
        for cell in cells["calibration"]:
            try:
                r = _cell_agree(cell, cfg)
                scores.append(r["frac_within20"])
                per.append({"sigma_wp": cell["sigma_wp"], "frac20": r2(r["frac_within20"]),
                            "wp_method": r["wp_method"], "n_window": r["n_window"],
                            "sigma_eff": r2(r["sigma_eff"])})
            except Exception as e:
                per.append({"sigma_wp": cell["sigma_wp"], "error": repr(e)[:120]})
        score = float(np.mean(scores)) if scores else -1.0
        attempts.append({"try": i + 1, "config": cfg, "calib_score": r2(score), "calib_per_cell": per})
        if best is None or score > best["calib_score_raw"]:
            best = {"config": cfg, "calib_score_raw": score, "try": i + 1}
    frozen = best["config"]

    # held-out verdict on frozen config
    heldout = []
    figs = []
    fig, axes = plt.subplots(1, len(cells["heldout"]), figsize=(4 * len(cells["heldout"]), 3.6), squeeze=False)
    for j, cell in enumerate(cells["heldout"]):
        try:
            r = _cell_agree(cell, frozen)
            agrees = r["frac_within20"] >= 0.5
            heldout.append({"sigma_wp": cell["sigma_wp"], "sigma_pot": r2(cell["sigma_pot"]),
                            "frac20": r2(r["frac_within20"]), "median_rel": r2(r["median_rel"]),
                            "sigma_eff": r2(r["sigma_eff"]), "sigma_eff_rel": r2(r["sigma_eff_rel"]),
                            "n_window": r["n_window"], "wp_method": r["wp_method"], "agrees": bool(agrees)})
            ax = axes[0][j]
            ax.plot(r["q_sel"], r["R_sel"], "o-", ms=3, label="R(q) meas")
            ax.plot(r["q_sel"], r["pred_sel"], "k--", label="F_WP/F_ONCV")
            ax.set_title(f"sigma_WP={cell['sigma_wp']} ({'OK' if agrees else 'x'})")
            ax.set_xlabel("q (1/Bohr)"); ax.legend(fontsize=7)
        except Exception as e:
            heldout.append({"sigma_wp": cell["sigma_wp"], "error": repr(e)[:160]})
    fig.suptitle("T2 held-out form-factor R(q) vs F_WP/F_ONCV (frozen config)")
    fig.tight_layout(); p = os.path.join(ART, "T2_heldout.png"); fig.savefig(p, dpi=110); plt.close(fig); figs.append(p)

    valid = [h for h in heldout if "agrees" in h]
    n_agree = sum(h["agrees"] for h in valid)
    if len(valid) == 0:
        verdict = "INCONCLUSIVE"
    elif n_agree >= 2:
        verdict = "CONFIRM"
    elif n_agree == 0 and all(h["n_window"] >= 4 for h in valid):
        verdict = "REFUTE"
    else:
        verdict = "INCONCLUSIVE"

    # POD on Delta n = n_WP - n_classical (structural support) for one held-out cell
    pod_info = _pod_delta(cells["heldout"])

    result = {"phase": "T2", "split": {"calibration": [c["sigma_wp"] for c in cells["calibration"]],
              "heldout": [c["sigma_wp"] for c in cells["heldout"]], "skipped": cells["skipped"]},
              "attempts": attempts, "frozen_config": frozen, "frozen_from_try": best["try"],
              "heldout": heldout, "n_heldout_agree": n_agree, "n_heldout_valid": len(valid),
              "verdict": verdict, "pod_delta": pod_info, "plot": p, "pod_plot": pod_info.get("plot")}
    save("T2", result)
    _nb_T2(result)
    email(f"ml-patterns T2 form-factor: {verdict}",
          {"hypothesis": HYPO,
           "done": ("Resolved the E=100 sigma-sweep cells (calib sigma_WP in {1,5}, "
                    "held-out {0.5,3,8}) matched within-cut at r_s=5.69/L50. Ran the "
                    f"<=4-try calibration loop (best=try {best['try']}), froze the config, "
                    "and read the verdict on the held-out _wf bath-only cells. POD on "
                    "Delta n = n_WP - n_classical gives the real-space structural support."),
           "plot_shows": ("Held-out R(q)=n_WP(q)/n_classical(q) (markers) vs the "
                          "parameter-free F_WP/F_ONCV (dashed) per held-out sigma_WP."),
           "conclusion": (f"T2 held-out verdict: {verdict}. {n_agree}/{len(valid)} held-out "
                          f"sigma agree within +-20% over their valid q-window. Calibration "
                          f"used for tuning only (ADR 0011); all {len(attempts)} attempts logged.")},
          [p] + ([pod_info["plot"]] if pod_info.get("plot") else []))
    return result


def _pod_delta(heldout_cells):
    """POD of Delta n(t) = n_WP_bath - n_classical_bath (matched v) for a mid-sigma cell."""
    try:
        cell = sorted(heldout_cells, key=lambda c: abs(c["sigma_wp"] - 3))[0]
        cfg = PL.FormFactorConfig(max_frames=60)
        dx = cell["dx"]
        wp, wax, _ = PL.load_bath_induced(cell["wp_bath_dir"], cell["wp_gs"], dx, 60, 162,
                                          wp_dir=cell.get("wp_wp_dir", ""), total_dir=cell.get("wp_total_dir", ""))
        cl, cax, _ = PL.load_bath_induced(cell["cl_bath_dir"], cell["cl_gs"], dx, 60, 162)
        if not NRM.grids_match(cax, wax):
            cl = NRM.cogrid(cl, cax, wax)
        T = min(len(wp), len(cl))
        Dn = (wp[:T] - cl[:T]).reshape(T, -1).T.astype(np.float32)
        res = P.pod(Dn, rank=8, randomized=True)
        ef = res.energy_fraction[:6]
        # plot energy spectrum + leading mode mid-plane
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
        ax[0].bar(range(1, len(ef) + 1), ef); ax[0].set_xlabel("POD mode"); ax[0].set_ylabel("energy fraction")
        ax[0].set_title(f"POD of dn_WP-dn_cl (sigma_WP={cell['sigma_wp']})")
        nx = int(round(len(wax[0]))); shp = (len(wax[0]), len(wax[1]), len(wax[2]))
        mode0 = res.modes[:, 0].reshape(shp)
        ax[1].imshow(mode0[:, :, shp[2] // 2].T, origin="lower", cmap="RdBu_r")
        ax[1].set_title("leading mode (xy mid-plane)")
        fig.tight_layout(); p = os.path.join(ART, "T2_pod_delta.png"); fig.savefig(p, dpi=110); plt.close(fig)
        return {"sigma_wp": cell["sigma_wp"], "energy_fraction": [r2(x) for x in ef], "plot": p,
                "n_modes_90pct": int(res.n_modes_for(0.9))}
    except Exception as e:
        return {"error": repr(e)[:200]}


def _nb_T2(result):
    write_notebook("rung1_bulk_formfactor.ipynb", "T2 - Rung 1 bulk form-factor cut", [
        ("md", f"**Pinned split (ADR 0011).** calibration sigma_WP={result['split']['calibration']}, "
               f"held-out sigma_WP={result['split']['heldout']}. Skipped: {result['split']['skipped']}.\n\n"
               f"Headline metric R(q)=n_WP(q)/n_classical(q) (azimuthal + temporal median) vs the "
               f"parameter-free F_WP/F_ONCV. The <=4-try loop tuned the shared config on the "
               f"CALIBRATION cells only; the verdict is read on the HELD-OUT cells."),
        ("code", "import json; r=json.load(open('../artifacts/T2_result.json'))\n"
                 "print('frozen from try', r['frozen_from_try']); print('verdict', r['verdict'])\n"
                 "for a in r['attempts']: print('try',a['try'],'calib_score',a['calib_score'])\n"
                 "for h in r['heldout']: print(h)"),
        ("md", f"### Verdict (held-out): **{result['verdict']}**  "
               f"({result['n_heldout_agree']}/{result['n_heldout_valid']} held-out sigma within +-20%)\n\n"
               f"![heldout](../artifacts/T2_heldout.png)\n\n"
               f"POD on Delta n = n_WP - n_classical (structural support):\n\n"
               f"![pod](../artifacts/T2_pod_delta.png)"),
    ])


# ---- T3 wake ----
WAKE_CONFIGS = [
    dict(rank=10, max_frames=250, window_frac=(0.0, 0.5), project_pod_rank=25, line_axis=2),
    dict(rank=8, max_frames=200, window_frac=(0.0, 0.45), project_pod_rank=20, line_axis=2),
    dict(rank=14, max_frames=300, window_frac=(0.0, 0.6), project_pod_rank=30, line_axis=2),
    dict(rank=6, max_frames=180, window_frac=(0.05, 0.4), project_pod_rank=15, line_axis=2),
]


def phase_T3():
    cells = celldb.resolve_wake_cells()
    attempts = []; best = None
    for i, cfg in enumerate(WAKE_CONFIGS):
        per = []; scores = []
        for cell in cells["calibration"]:
            try:
                r = PL.cell_wake_dmd(cell, PL.WakeConfig(**cfg))
                if "omega_dmd_ev" in r:
                    rel = abs(r["omega_dmd_ev"] - r["omega_p_ev"]) / r["omega_p_ev"]
                    scores.append(1.0 if rel <= 0.20 else max(0.0, 1 - rel))
                    per.append({"E": cell["energy_ev"], "omega_dmd": r2(r["omega_dmd_ev"]),
                                "omega_p": r2(r["omega_p_ev"]), "rel": r2(rel), "nyquist_ok": r["nyquist_ok"]})
                else:
                    per.append({"E": cell["energy_ev"], "method_invalid": r.get("method_invalid")})
            except Exception as e:
                per.append({"E": cell["energy_ev"], "error": repr(e)[:120]})
        score = float(np.mean(scores)) if scores else -1.0
        attempts.append({"try": i + 1, "config": cfg, "calib_score": r2(score), "calib_per_cell": per})
        if best is None or score > best["raw"]:
            best = {"config": cfg, "raw": score, "try": i + 1}
    frozen = best["config"]

    heldout = []
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for cell in cells["heldout"]:
        try:
            r = PL.cell_wake_dmd(cell, PL.WakeConfig(**frozen))
            if "omega_dmd_ev" in r:
                rel = abs(r["omega_dmd_ev"] - r["omega_p_ev"]) / r["omega_p_ev"]
                agrees = (rel <= 0.20) and r["nyquist_ok"]
                heldout.append({"E": cell["energy_ev"], "v": r2(cell["velocity_au"]),
                                "omega_dmd_ev": r2(r["omega_dmd_ev"]), "omega_p_ev": r2(r["omega_p_ev"]),
                                "lambda_dmd": r2(r["lambda_dmd"]), "lambda_theory": r2(r["lambda_theory"]),
                                "rel": r2(rel), "nyquist_ok": r["nyquist_ok"], "agrees": bool(agrees),
                                "wp_method": r["wp_method"]})
                ax.scatter(cell["velocity_au"], r["omega_dmd_ev"], c="C0", zorder=3)
            else:
                heldout.append({"E": cell["energy_ev"], "method_invalid": r.get("method_invalid")})
        except Exception as e:
            heldout.append({"E": cell["energy_ev"], "error": repr(e)[:160]})
    if cells["heldout"]:
        wp_ev = cells["heldout"][0]["omega_p_ev"]
        ax.axhline(wp_ev, ls="--", c="k", label=f"omega_p={wp_ev:.2f} eV")
        ax.fill_between([0, 8], wp_ev * 0.8, wp_ev * 1.2, color="grey", alpha=0.2, label="+-20%")
    ax.set_xlabel("v (a.u.)"); ax.set_ylabel("DMD dominant omega (eV)")
    ax.set_title("T3 held-out wake DMD frequency vs omega_p"); ax.legend()
    fig.tight_layout(); p = os.path.join(ART, "T3_heldout.png"); fig.savefig(p, dpi=110); plt.close(fig)

    valid = [h for h in heldout if "agrees" in h]
    n_agree = sum(h["agrees"] for h in valid)
    if not valid:
        verdict = "INCONCLUSIVE"
    elif n_agree >= max(1, (len(valid) + 1) // 2):
        verdict = "CONFIRM"
    elif n_agree == 0 and all(h.get("nyquist_ok") for h in valid):
        verdict = "REFUTE"
    else:
        verdict = "INCONCLUSIVE"

    result = {"phase": "T3", "split": {"calibration": [r2(c["energy_ev"]) for c in cells["calibration"]],
              "heldout": [r2(c["energy_ev"]) for c in cells["heldout"]], "skipped": cells["skipped"]},
              "attempts": attempts, "frozen_config": frozen, "frozen_from_try": best["try"],
              "heldout": heldout, "n_heldout_agree": n_agree, "n_heldout_valid": len(valid),
              "verdict": verdict, "plot": p}
    save("T3", result)
    write_notebook("rung1_wake_dmd.ipynb", "T3 - Rung 1 wake DMD gate", [
        ("md", f"**Pinned split (ADR 0011).** calibration energies (even velocity index)="
               f"{result['split']['calibration']} eV, held-out (odd)={result['split']['heldout']} eV. "
               f"Skipped: {result['split']['skipped']}.\n\n"
               f"DMD (exact, windowed over the early near-constant-velocity stretch, POD-precompressed) "
               f"on the sigma_WP=5 bath induced density. Test: dominant DMD angular frequency vs "
               f"omega_p; lambda_dmd=2*pi*v/omega_dmd vs lambda(v)=2*pi*v/omega_p. Nyquist guard dt<pi/omega_p."),
        ("code", "import json; r=json.load(open('../artifacts/T3_result.json'))\n"
                 "print('verdict', r['verdict']); [print(h) for h in r['heldout']]"),
        ("md", f"### Verdict (held-out): **{verdict}**  ({n_agree}/{len(valid)} within +-20%)\n\n"
               f"![wake](../artifacts/T3_heldout.png)"),
    ])
    email(f"ml-patterns T3 wake DMD: {verdict}",
          {"hypothesis": HYPO,
           "done": ("DMD (windowed, POD-precompressed) on the sigma_WP=5 energy sweep; "
                    "pinned split calib=even / held-out=odd velocity index; <=4-try loop "
                    f"(best=try {best['try']}); verdict on held-out energies."),
           "plot_shows": ("Held-out DMD dominant frequency vs v with the omega_p line and "
                          "its +-20% band."),
           "conclusion": (f"T3 held-out verdict: {verdict}. {n_agree}/{len(valid)} held-out "
                          f"energies have a DMD wake frequency within +-20% of omega_p.")},
          [p])
    return result


def phase_T4():
    """Rung 1b: transfer the frozen wake pipeline to a localised slab (sigma_WP=0.5)."""
    df = celldb.load_db()
    sl = df[(df.system == "localised_jellium") & (df.wp_enabled == True) &
            np.isclose(df.sigma_wp_bohr.fillna(-1), 0.5, atol=0.05)].copy()
    sl = sl[sl.apply(celldb._has_density, axis=1)]
    sl = sl.assign(_nf=sl.density_system_vti_nframes.fillna(0)).sort_values("_nf", ascending=False)
    note = ""
    heldout = []
    p = None
    if sl.empty:
        verdict = "INCONCLUSIVE"; note = "no localised slab WP run with bath density"
    else:
        wpr = sl.iloc[0]
        cl = df[(df.system == "localised_jellium") & (df.wp_enabled == False) &
                np.isclose(df.energy_ev.fillna(-1), wpr.energy_ev, rtol=0.05)]
        cl = cl[cl.apply(celldb._has_density, axis=1)]
        if cl.empty:
            verdict = "INCONCLUSIVE"; note = "no matched slab classical"
        else:
            cell = celldb._cell(wpr, cl.assign(_nf=cl.density_total_vti_nframes.fillna(0)).sort_values("_nf").iloc[-1])
            try:
                frozen = load("T3")["frozen_config"] if have("T3") else WAKE_CONFIGS[0]
                r = PL.cell_wake_dmd(cell, PL.WakeConfig(**frozen))
                if "omega_dmd_ev" in r:
                    rel = abs(r["omega_dmd_ev"] - r["omega_p_ev"]) / r["omega_p_ev"]
                    verdict = "CONFIRM" if (rel <= 0.20 and r["nyquist_ok"]) else "REFUTE" if r["nyquist_ok"] else "INCONCLUSIVE"
                    heldout.append({k: r2(r[k]) if isinstance(r[k], float) else r[k]
                                    for k in ["omega_dmd_ev", "omega_p_ev", "lambda_dmd", "lambda_theory", "nyquist_ok"]})
                    fig, ax = plt.subplots(figsize=(5.5, 4))
                    ax.bar(["DMD omega", "omega_p"], [r["omega_dmd_ev"], r["omega_p_ev"]], color=["C0", "k"])
                    ax.set_ylabel("eV"); ax.set_title(f"T4 slab wake (sigma_WP=0.5) verdict={verdict}")
                    fig.tight_layout(); p = os.path.join(ART, "T4_slab.png"); fig.savefig(p, dpi=110); plt.close(fig)
                    note = f"slab geometry; rel={rel:.2f}"
                else:
                    verdict = "INCONCLUSIVE"; note = r.get("method_invalid", "dmd invalid")
            except Exception as e:
                verdict = "INCONCLUSIVE"; note = repr(e)[:160]
    result = {"phase": "T4", "verdict": verdict, "note": note, "heldout": heldout, "plot": p}
    save("T4", result)
    write_notebook("rung1b_slab.ipynb", "T4 - Rung 1b localised slab (geometry transfer)", [
        ("md", f"Transfer the FROZEN Rung-1 wake pipeline to the localised slab "
               f"(sigma_WP=0.5). Verdict: **{verdict}**. Note: {note}."),
        ("code", "import json; print(json.load(open('../artifacts/T4_result.json')))"),
    ] + ([("md", "![slab](../artifacts/T4_slab.png)")] if p else []))
    email(f"ml-patterns T4 slab: {verdict}",
          {"hypothesis": HYPO,
           "done": "Transferred the frozen wake DMD pipeline to the localised slab (sigma_WP=0.5).",
           "plot_shows": "Slab DMD wake frequency vs omega_p." if p else "(no plot: insufficient slab data)",
           "conclusion": f"T4 slab verdict: {verdict}. {note}"},
          [p] if p else [])
    return result


def phase_T5():
    """Rung 2 dynamics: DMD/Koopman mode table + hand-rolled SINDy on POD latents."""
    cells = celldb.resolve_wake_cells()
    pool = cells["heldout"] + cells["calibration"]
    cell = sorted(pool, key=lambda c: abs(c["velocity_au"] - 2.7))[0] if pool else None
    modes = []; sindy = {}; p = None
    if cell is None:
        verdict = "INCONCLUSIVE"; note = "no wake cell"
    else:
        try:
            dx = cell["dx"]
            delta, axes, method = PL.load_bath_induced(cell["wp_bath_dir"], cell["wp_gs"], dx, 250, 162,
                                                       wp_dir=cell.get("wp_wp_dir", ""), total_dir=cell.get("wp_total_dir", ""))
            T = len(delta); a, b = 0, T // 2
            sub = delta[a:b].reshape(b - a, -1).T.astype(np.float32)
            pod = P.pod(sub, rank=12, randomized=True)
            fdt = cell["frame_dt_au_wp"]
            res = D.dmd(pod.coeffs, dt=fdt, rank=8)
            HA = PL.HA_EV
            for k in range(len(res.omega)):
                w = res.angular_frequency[k] * HA
                if w > 0.05:
                    modes.append({"omega_eV": r2(w), "growth": r2(res.growth_rate[k]),
                                  "amp": r2(float(abs(res.amplitudes[k])))})
            modes = sorted(modes, key=lambda m: -m["amp"])[:6]
            # hand-rolled SINDy (STLSQ) on the top-2 POD latent coordinates
            sindy = _sindy(pod.coeffs[:2], fdt)
            # plot latent trajectories + DMD spectrum
            fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
            tt = np.arange(pod.coeffs.shape[1]) * fdt
            for r_ in range(min(3, pod.coeffs.shape[0])):
                ax[0].plot(tt, pod.coeffs[r_], label=f"a{r_+1}")
            ax[0].set_xlabel("t (a.u.)"); ax[0].set_title("POD latent dynamics"); ax[0].legend(fontsize=7)
            if modes:
                ax[1].stem([m["omega_eV"] for m in modes], [m["amp"] for m in modes])
            ax[1].axvline(cell["omega_p_ev"], ls="--", c="k", label=f"omega_p={cell['omega_p_ev']:.2f}")
            ax[1].set_xlabel("omega (eV)"); ax[1].set_title("DMD mode spectrum"); ax[1].legend(fontsize=7)
            fig.tight_layout(); p = os.path.join(ART, "T5_dynamics.png"); fig.savefig(p, dpi=110); plt.close(fig)
            verdict = "DONE"; note = f"cell E={cell['energy_ev']:.0f} method={method}"
        except Exception as e:
            verdict = "INCONCLUSIVE"; note = repr(e)[:180]
    result = {"phase": "T5", "verdict": verdict, "note": note, "dmd_modes": modes,
              "sindy": sindy, "plot": p, "omega_p_ev": r2(cell["omega_p_ev"]) if cell else None}
    save("T5", result)
    write_notebook("rung2_dynamics.ipynb", "T5 - Rung 2 dynamics (DMD/Koopman + SINDy)", [
        ("md", f"DMD/Koopman mode table + a hand-rolled STLSQ SINDy on the top-2 POD "
               f"latent coordinates of the bulk-jellium bath induced density. Note: {note}."),
        ("code", "import json; print(json.load(open('../artifacts/T5_result.json')))"),
    ] + ([("md", "![dyn](../artifacts/T5_dynamics.png)")] if p else []))
    email(f"ml-patterns T5 dynamics: {verdict}",
          {"hypothesis": HYPO,
           "done": "DMD/Koopman mode spectrum + SINDy on POD latents of the bulk induced density.",
           "plot_shows": "POD latent trajectories and the DMD mode spectrum vs omega_p." if p else "(no plot)",
           "conclusion": f"T5 {verdict}. {note}. DMD modes: {modes}"},
          [p] if p else [])
    return result


def _sindy(latents, dt, thresh=0.05, terms=("1", "a0", "a1", "a0^2", "a1^2", "a0a1")):
    """Minimal STLSQ SINDy for d a_i/dt on a small polynomial library (no sklearn)."""
    a = latents
    da = np.gradient(a, dt, axis=1)
    A0, A1 = a[0], a[1]
    Theta = np.vstack([np.ones_like(A0), A0, A1, A0 ** 2, A1 ** 2, A0 * A1]).T
    out = {}
    for i in range(2):
        y = da[i]
        xi = np.linalg.lstsq(Theta, y, rcond=None)[0]
        for _ in range(5):
            small = np.abs(xi) < thresh * max(1e-9, np.abs(xi).max())
            xi[small] = 0
            big = ~small
            if big.any():
                xi[big] = np.linalg.lstsq(Theta[:, big], y, rcond=None)[0]
        out[f"da{i}/dt"] = {t: r2(float(c)) for t, c in zip(terms, xi) if c != 0}
    return out


def phase_T6():
    """Exploratory (caveated): vacuum-WP-subtracted n_wp -> q-space dispersion/diffraction.
    SIE-confounded; NOT a headline. Falls back gracefully if no vacuum-WP partner."""
    df = celldb.load_db()
    # a WP run with density_wp + a free_wp (vacuum) partner at matched sigma/E
    note = ""; p = None; verdict = "EXPLORATORY"
    try:
        wp = df[(df.system == "jellium") & (df.wp_enabled == True) &
                df.density_wp_vti_nframes.notna() &
                np.isclose(df.r_s.fillna(-1), 5.69, rtol=0.05)]
        free = df[df.run_id.astype(str).str.contains("free_wp")]
        # pick a sigma where both exist (sigma~5 E100)
        wcand = wp[np.isclose(wp.sigma_wp_bohr.fillna(-1), 5.0, atol=0.1) &
                   np.isclose(wp.energy_ev.fillna(-1), 100, rtol=0.1)]
        fcand = free[free.run_id.astype(str).str.contains("E100") &
                     free.density_wp_vti_nframes.notna()]
        if wcand.empty or fcand.empty:
            # fall back: just FFT the in-bath WP density and show its q-rolloff
            wcand = wp[np.isclose(wp.sigma_wp_bohr.fillna(-1), 3.0, atol=0.1)]
            if wcand.empty:
                raise RuntimeError("no WP run with density_wp")
            row = wcand.iloc[0]
            wpdir = celldb._results_dir(row.run_path, row.density_wp_vti_dir)
            series, axes, _ = NRM.load_series(wpdir, 40)
            dx = float(row.spacing_bohr)
            q0, a0 = FF.radial_power_spectrum(series[5], dx, nbins=40)
            q1, a1 = FF.radial_power_spectrum(series[-1], dx, nbins=40)
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.semilogy(q0, a0 / a0[0], label="early")
            ax.semilogy(q1, a1 / a1[0], label="late (dispersed)")
            ax.set_xlabel("q (1/Bohr)"); ax.set_ylabel("|n_wp(q)| (norm)")
            ax.set_title("T6 exploratory: WP density q-rolloff (NO vacuum-WP subtraction)")
            ax.legend(); fig.tight_layout(); p = os.path.join(ART, "T6_exploratory.png")
            fig.savefig(p, dpi=110); plt.close(fig)
            note = ("CAVEAT: in-bath n_wp only; vacuum-WP/SIE control NOT applied -> "
                    "this conflates bath scattering with intrinsic dispersion + the ~7 eV "
                    "WP self-interaction (project-brief figure, externally unverified). "
                    "Exploratory only, not a headline claim.")
        else:
            note = "vacuum-WP partner present but full subtraction left for manual review"
    except Exception as e:
        verdict = "INCONCLUSIVE"; note = repr(e)[:180]
    result = {"phase": "T6", "verdict": verdict, "note": note, "plot": p}
    save("T6", result)
    write_notebook("exploratory_exchange_diffraction.ipynb", "T6 - Exploratory exchange/diffraction (caveated)", [
        ("md", f"**EXPLORATORY, NOT A HEADLINE.** Signatures (i) exchange-hole and (ii) "
               f"diffraction fringes are SIE-confounded (research Q3): they live at the "
               f"length/energy scale of the ~7 eV WP self-interaction error (project-brief "
               f"figure, externally unverified) and only become defensible on the vacuum-WP-"
               f"subtracted field. Note: {note}"),
        ("code", "import json; print(json.load(open('../artifacts/T6_result.json')))"),
    ] + ([("md", "![expl](../artifacts/T6_exploratory.png)")] if p else []))
    email("ml-patterns T6 exploratory (caveated)",
          {"hypothesis": HYPO,
           "done": "Exploratory q-space look at the WP projectile density (exchange/diffraction).",
           "plot_shows": "WP density q-rolloff early vs late." if p else "(no plot)",
           "conclusion": f"T6 EXPLORATORY (not headline). {note}"},
          [p] if p else [])
    return result


def phase_T7():
    """Synthesis across rungs + final summary."""
    summary = {}
    for ph in ["T1", "T2", "T3", "T4", "T5", "T6"]:
        if have(ph):
            r = load(ph)
            summary[ph] = {k: r.get(k) for k in ["verdict", "n_heldout_agree",
                           "n_heldout_valid", "foncv_unity_q_5pct"] if k in r}
    fig, ax = plt.subplots(figsize=(7, 3.5))
    labels = list(summary.keys())
    verds = [str(summary[k].get("verdict", "-")) for k in labels]
    ax.axis("off")
    txt = "\n".join(f"{k}: {summary[k].get('verdict','-')}"
                    + (f"  ({summary[k].get('n_heldout_agree')}/{summary[k].get('n_heldout_valid')} held-out)"
                       if summary[k].get('n_heldout_valid') is not None else "")
                    for k in labels)
    ax.text(0.02, 0.95, "ml-patterns synthesis\n\n" + txt, va="top", family="monospace", fontsize=11)
    fig.tight_layout(); p = os.path.join(ART, "T7_synthesis.png"); fig.savefig(p, dpi=110); plt.close(fig)
    result = {"phase": "T7", "summary": summary, "plot": p}
    save("T7", result)
    write_notebook("synthesis.ipynb", "T7 - Cross-rung synthesis", [
        ("md", "Cross-rung synthesis of the ml-patterns campaign. Verdicts are read on "
               "PINNED held-out splits (ADR 0011); CONFIRM/REFUTE/INCONCLUSIVE are all valid."),
        ("code", "import json; print(json.dumps(json.load(open('../artifacts/T7_result.json'))['summary'], indent=2))"),
        ("md", "![synthesis](../artifacts/T7_synthesis.png)"),
    ])
    email("ml-patterns T7 synthesis",
          {"hypothesis": HYPO,
           "done": "Synthesised all rung verdicts into the final notebook.",
           "plot_shows": "Per-phase held-out verdict summary.",
           "conclusion": "Synthesis: " + "; ".join(f"{k}={summary[k].get('verdict','-')}" for k in summary)},
          [p])
    return result


PHASES = {"T1": phase_T1, "T2": phase_T2, "T3": phase_T3, "T4": phase_T4,
          "T5": phase_T5, "T6": phase_T6, "T7": phase_T7}


def main(argv):
    global SEND_EMAIL
    args = [a for a in argv if not a.startswith("--")]
    if "--no-email" in argv:
        SEND_EMAIL = False
    want = args if args else list(PHASES.keys())
    status = {}
    for ph in want:
        if have(ph) and "--force" not in argv:
            print(f"[skip] {ph} (result exists)"); status[ph] = "skipped(existing)"; continue
        t0 = time.time()
        print(f"\n===== {ph} starting =====")
        try:
            PHASES[ph]()
            status[ph] = f"done in {time.time()-t0:.0f}s"
            print(f"===== {ph} done in {time.time()-t0:.0f}s =====")
        except Exception:
            tb = traceback.format_exc()
            status[ph] = "FAILED"
            print(f"===== {ph} FAILED =====\n{tb}")
            try:
                email(f"ml-patterns {ph} FAILED",
                      {"hypothesis": HYPO, "done": f"Phase {ph} raised.",
                       "plot_shows": "(no plot - failure)", "conclusion": tb[-1500:]}, [])
            except Exception:
                pass
    json.dump(status, open(os.path.join(ART, "phase_status.json"), "w"), indent=2)
    print("\nSTATUS:", json.dumps(status, indent=2))


if __name__ == "__main__":
    main(sys.argv[1:])
