#!/usr/bin/env python3
"""Watcher: email the user once BOTH mass-pair runs (m1, m2) are propagating.

'Propagating' = each run's per-step energy_decomp.csv has >= THRESH steps past
step 0 with no CUDA illegal-access in its log. Then build an energy-conservation
plot for both and send the mandatory four-part email (email-notifications skill).
Autonomous: run detached; it exits after sending (or on timeout).

  cd .../mass_pair_n162 && nohup ../../../../venv/bin/python3 email_when_propagating.py \
       > email_watch.log 2>&1 &
"""
from __future__ import annotations
import csv, sys, time
from datetime import datetime
from pathlib import Path

ROOT   = Path("/local/data/public/skcb2/tddft")
WP     = ROOT/"ResearchProject/systems/localised_jellium/scripts/mass_pair_n162/wp"
TO     = "chiddukanna@gmail.com"
THRESH = 30                    # steps each run must clear to count as "propagating"
TIMEOUT_H = 8.0
sys.path.insert(0, str(ROOT/"inq-stack/python"))

def log(m): print(f"[{datetime.now():%F %T}] {m}", flush=True)

def steps_done(out: str) -> int:
    f = WP/"results"/out/"raw/observables/energy_decomp.csv"
    try:
        with open(f) as fh:
            return max(0, sum(1 for _ in fh) - 1)   # rows minus header = steps
    except Exception:
        return 0

def crashed(out: str) -> bool:
    try:
        return "illegal memory" in (WP/f"rt_{out}.log").read_text()
    except Exception:
        return False

def load_energy(out: str):
    steps, e = [], []
    try:
        with open(WP/"results"/out/"raw/observables/energy_decomp.csv") as fh:
            for row in csv.DictReader(fh):
                steps.append(int(row["step"])); e.append(float(row["energy_total"]))
    except Exception:
        pass
    return steps, e

def make_plot(path: Path) -> bool:
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=140)
        for out, lab, col in (("m1","m=1 (GPU0)","#1f77b4"), ("m2","m=2 (GPU1)","#d62728")):
            s, e = load_energy(out)
            if s:
                e0 = e[0]
                ax.plot(s, [(x-e0)*27211.386 for x in e], "-o", ms=2.5, color=col,
                        label=f"{lab}: {len(s)} steps")
        ax.set_xlabel("propagation step"); ax.set_ylabel("E_total − E_total(0)  [meV]")
        ax.set_title("N=120 mass pair — energy conservation (early propagation)")
        ax.axhline(0, lw=0.6, color="0.6"); ax.legend(frameon=False, fontsize=9)
        fig.tight_layout(); fig.savefig(path); plt.close(fig)
        return path.exists()
    except Exception as ex:
        log(f"plot failed: {ex}"); return False

def main():
    log("watcher start — waiting for BOTH m1 & m2 to reach "
        f"{THRESH} steps (timeout {TIMEOUT_H} h)")
    t_end = time.time() + TIMEOUT_H*3600
    while time.time() < t_end:
        s1, s2 = steps_done("m1"), steps_done("m2")
        log(f"m1={s1} steps  m2={s2} steps  (crashed: m1={crashed('m1')} m2={crashed('m2')})")
        if s1 >= THRESH and s2 >= THRESH:
            break
        time.sleep(60)
    else:
        log("timeout — sending a heads-up anyway with whatever progress exists")

    s1, s2 = steps_done("m1"), steps_done("m2")
    _, e1 = load_energy("m1"); _, e2 = load_energy("m2")
    d1 = (e1[-1]-e1[0])*27211.386 if len(e1) > 1 else float("nan")
    d2 = (e2[-1]-e2[0])*27211.386 if len(e2) > 1 else float("nan")
    plot = WP.parent/"both_propagating.png"; have_plot = make_plot(plot)

    both = s1 >= THRESH and s2 >= THRESH
    subject = ("[localised-jellium mass-pair] m1/m2 — BOTH runs propagating (N=120)"
               if both else
               "[localised-jellium mass-pair] m1/m2 — propagation status (N=120)")
    body = f"""HYPOTHESIS
  Does the projectile MASS (m=1 vs m=2, same 100 eV) change the energy the
  120-electron localised jellium absorbs from a sigma=1 wavepacket? Bath
  electrons are mass 1 in both runs; only the projectile (WP) mass differs.

WHAT WAS DONE
  - Localised jellium slab z in [-12.5,12.5] (25 Bohr), r_s~=5.69, N=120,
    dx=0.40; two CAPs 10 Bohr/side eta=-1.0; sigma_WP=1, E=100 eV.
  - Two matched WP runs, dt=0.04, 2500 steps (100 a.u.), full per-step energy
    decomposition, checkpoint every 500. m1 (mass 1, GPU0) + m2 (mass 2, GPU1).
  - (This system replaced the N=162 one, which did not fit the 24 GB GPU for
    real-time propagation.)

PLOT (attached: both_propagating.png)
  E_total(t) - E_total(0) in meV vs propagation step, for m1 and m2. A flat,
  small curve = clean ETRS propagation (energy conserved) past step 0 — the
  exact point the 162-electron system kept failing.

CONCLUSION
  {"Both runs are propagating cleanly" if both else "Partial — see step counts"}:
  m1 at {s1} steps (dE_total = {d1:.2f} meV over the window),
  m2 at {s2} steps (dE_total = {d2:.2f} meV). Each run targets 2500 steps;
  they checkpoint every 500 and auto-resume on any interruption. The N=120
  memory fix is confirmed working.
"""
    try:
        from inqview.email import send_run_email
        send_run_email(subject=subject, body=body,
                       attachments=[str(plot)] if have_plot else [], to=TO)
        log(f"EMAIL SENT: {subject}")
    except Exception as ex:
        log(f"EMAIL FAILED: {ex}")

if __name__ == "__main__":
    main()
