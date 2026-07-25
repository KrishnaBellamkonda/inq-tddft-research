#!/usr/bin/env python3
"""Per-run quick analysis + validation + threaded email for a graphene CAP run.

Usage: post_and_email.py <run_outdir> <tag>
  <run_outdir> contains results/ (run_summary.txt, raw/observables/...).
Threads under the [graphene-cap] family via the stored message-id.
"""
import sys, os, csv
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

OUT = sys.argv[1]; TAG = sys.argv[2]
RES = os.path.join(OUT, "results")
THREAD = "/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/.graphene_cap_email_thread.txt"

def read_summary(p):
    d = {}
    if os.path.exists(p):
        for ln in open(p):
            if "=" in ln:
                k, v = ln.split("=", 1); d[k.strip()] = v.strip()
    return d

summ = read_summary(os.path.join(RES, "run_summary.txt"))

# --- figure: survival/absorbed(t) + energies(t) ---
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
inn = os.path.join(RES, "raw/observables/inner_norm_vs_time.csv")
if os.path.exists(inn):
    r = list(csv.DictReader(open(inn)))
    t = np.array([float(x["time_au"]) for x in r])
    surv = np.array([float(x["survival_inner_over_N0"]) for x in r])
    absb = np.array([float(x["absorbed_fraction"]) for x in r])
    ax[0].plot(t, surv, label="survival (inner/N0)")
    ax[0].plot(t, absb, label="absorbed fraction")
    ax[0].set_xlabel("time (a.u.)"); ax[0].set_ylabel("fraction"); ax[0].set_ylim(-0.02, 1.05)
    ax[0].legend(); ax[0].set_title("WP survival / absorption")
obs = os.path.join(RES, "raw/observables/observables.csv")
if os.path.exists(obs):
    r = list(csv.DictReader(open(obs)))
    def col(c): return np.array([float(x[c]) for x in r]) if r and c in r[0] else None
    t = col("time_au")
    for c, lab in [("energy_total", "E_total"), ("energy_kinetic", "E_kin")]:
        y = col(c)
        if y is not None: ax[1].plot(t, y - y[0], label=f"Δ{lab}")
    ax[1].set_xlabel("time (a.u.)"); ax[1].set_ylabel("ΔE (Ha)")
    ax[1].legend(); ax[1].set_title("Energy drift / exchange")
cap = summ.get("cap", "?"); eps = summ.get("epsilon_survival", "?"); absf = summ.get("absorbed_fraction", "?")
fig.suptitle(f"graphene CAP run: {TAG}  (CAP={cap}, E={summ.get('E_eV','?')} eV)  [feasibility replica, PROVISIONAL]")
fig.tight_layout()
figpath = os.path.join(RES, "quicklook.png"); fig.savefig(figpath, dpi=120)

# --- validation (best-effort) ---
val = "n/a"
try:
    from inqview.validation import validate_run
    rep = validate_run(RES); val = "PASS" if getattr(rep, "passed", None) else str(rep)
except Exception as e:
    val = f"(validator skipped: {e})"

# --- threaded email ---
prior = []
if os.path.exists(THREAD):
    prior = [l.strip() for l in open(THREAD) if l.strip()]
in_reply_to = prior[-1] if prior else None

body = f"""Graphene CAP run complete: {TAG}

run_completed   = {summ.get('run_completed','?')}
CAP             = {cap}   (eta={summ.get('eta_Ha','?')} Ha, L={summ.get('L_cap','?')} Bohr)
E               = {summ.get('E_eV','?')} eV,  k0={summ.get('k0','?')}
impact (cx,cy)  = ({summ.get('cx','?')}, {summ.get('cy','?')}) Bohr
survival eps    = {eps}
absorbed frac   = {absf}
WP norm@inject  = {summ.get('norm_after_inject','?')}, max_overlap={summ.get('max_overlap','?')}
N_steps/dt      = {summ.get('N_STEPS','?')} / {summ.get('dt','?')}
wall_s          = {summ.get('wall_s','?')}
validation      = {val}
outputs         = {RES}/raw/ (density VTIs total/system/wp, WP wavefunction, LEED screens, momentum, overlap)

Feasibility replica (24 C, 3x2). CAP results PROVISIONAL until inq-study Task #7.
Plan/handover: docs/plans/graphene-cap.md, docs/handovers/graphene-cap.md.
"""
try:
    from inqview.email import send_run_email
    mid = send_run_email(
        subject=f"[graphene-cap] run done: {TAG} (eps={eps})",
        body=body, attachments=[figpath] if os.path.exists(figpath) else None,
        to="chiddukanna@gmail.com", in_reply_to=in_reply_to,
        references=prior or None)
    with open(THREAD, "a") as f: f.write(mid + "\n")
    print(f"emailed {TAG}: {mid}")
except Exception as e:
    print(f"email failed for {TAG}: {e}")
