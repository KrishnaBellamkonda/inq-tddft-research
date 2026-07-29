#!/usr/bin/env python3
"""twosided_cap_vs_mask dispatcher — 3 phases, emailed on completion.

Plan: docs/plans/twosided-cap-vs-mask.md. Two-sided absorbers (L split ½ each end),
energy-scaled quasi-monochromatic packet (σ=4√2/k0). Reuses ONE env-driven binary
(scripts/twosided_cap_vs_mask/run) in two modes (cap | mask). inq/ untouched.

Phases (each emails its key figure on completion, threaded):
  1. CAP η-sweep at L=20  → ε(E) | η
  2. CAP L-sweep at η=-0.5 → ε(E) | L (CAP)
  3. mask L-sweep          → ε(E) | L (mask) + side-by-side
At the end: auto-build the study notebook (notebook-making skill) + final email.

    PYTHONPATH=.../inq-stack/python python3 dispatch.py     # ~120 runs, GPUs 0,1

ε PROVISIONAL until the inq-study engine regression (Task #7).
"""
import os, sys, math, time, subprocess
from pathlib import Path

ROOT = Path("/local/data/public/skcb2/tddft")
INQ_STUDY = ROOT / "inq-study"
SHARE = ROOT / "inq/install/share"
STACK = ROOT / "inq-stack/python"
SYS = ROOT / "ResearchProject/systems/vacuum"
BINARY = SYS / "scripts/twosided_cap_vs_mask/run"
SWEEP = SYS / "twosided_cap_vs_mask"
HYP = SYS / "hypotheses/twosided_cap_vs_mask"
GPUS = [0, 1]
TO = "chiddukanna@gmail.com"
HA = 27.211386245988

# grids (plan-locked); energies descending so cheap high-E runs finish first.
# 1 eV dropped: σ=4√2/k0 makes its box ~250 Bohr / ~1e5 steps (~50 min/run) —
# disproportionate cost; 2 eV (~10^0.3) anchors the low decade. (User flagged
# 1 eV droppable earlier.)
E_FULL = [1000, 300, 100, 64, 32, 16, 10, 4, 2]
E_ETA  = [300, 100, 32, 10, 2]          # reduced set for the η-sweep trend
LS     = [10, 16, 20, 26, 30]
ETAS   = [-0.30, -0.50, -0.70, -1.00]
ETA_STAR = -0.50
ANCHOR_L = 20
ANCHOR_E = 10

def k0_of(E):  return math.sqrt(2.0 * E / HA)

def env_for(gpu):
    e = dict(os.environ)
    e["INQ_SOURCE"] = str(INQ_STUDY)
    e["INQ_SHARE_PATH"] = str(SHARE)
    e["PSEUDOPOD_SHARE_PATH"] = str(SHARE / "pseudopod")
    e["CUDA_VISIBLE_DEVICES"] = str(gpu)
    e["PYTHONPATH"] = str(STACK)
    return e

def run_name(mode, E, L, eta):
    if mode == "cap":
        return f"run_cap_E{E}_L{L}_eta{abs(eta):.2f}"
    return f"run_mask_E{E}_L{L}"

def jobs_phase1():   # CAP η-sweep at anchor L
    return [dict(mode="cap", E=E, L=ANCHOR_L, eta=eta) for eta in ETAS for E in E_ETA]
def jobs_phase2():   # CAP L-sweep at η*
    return [dict(mode="cap", E=E, L=L, eta=ETA_STAR) for L in LS for E in E_FULL]
def jobs_phase3():   # mask L-sweep
    return [dict(mode="mask", E=E, L=L, eta=0.0) for L in LS for E in E_FULL]

def launch(j, gpu):
    d = SWEEP / run_name(j["mode"], j["E"], j["L"], j["eta"])
    d.mkdir(parents=True, exist_ok=True)
    if (d / "results/epsilon.txt").exists():
        return None                      # resumable: skip completed
    e = env_for(gpu)
    e.update(CAP_MODE=j["mode"], CAP_K0=f"{k0_of(j['E']):.6f}", CAP_L=str(j["L"]),
             CAP_ETA=str(j["eta"]), CAP_OUTDIR=str(d / "results"))
    log = open(d / "run.log", "w")
    return subprocess.Popen([str(BINARY)], cwd=d, env=e, stdout=log, stderr=subprocess.STDOUT)

def run_jobs(jobs, label):
    pending = list(jobs)
    running, t0 = {}, time.time()
    print(f"==> phase '{label}': {len(pending)} jobs", flush=True)
    while pending or running:
        for gpu in GPUS:
            if gpu not in running and pending:
                j = pending.pop(0)
                p = launch(j, gpu)
                if p is None:
                    print(f"  skip (done) {run_name(j['mode'],j['E'],j['L'],j['eta'])}", flush=True)
                    continue
                print(f"  launch {run_name(j['mode'],j['E'],j['L'],j['eta'])} GPU{gpu}", flush=True)
                running[gpu] = (p, j)
        for gpu in list(running):
            p, j = running[gpu]
            if p.poll() is not None:
                print(f"  done   {run_name(j['mode'],j['E'],j['L'],j['eta'])} rc={p.returncode}", flush=True)
                del running[gpu]
        time.sleep(2)
    print(f"==> phase '{label}' complete in {(time.time()-t0)/60:.1f} min", flush=True)

# ---- result parsing + phase figures + email -------------------------------
def parse(p):
    out = {}
    for ln in Path(p).read_text().splitlines():
        k, _, v = ln.partition(' ')
        try: out[k] = float(v)
        except ValueError: out[k] = v
    return out

def collect(mode):
    recs = []
    for d in sorted(SWEEP.glob(f"run_{mode}_*")):
        f = d / "results/epsilon.txt"
        if f.exists():
            try:
                r = parse(f); r["name"] = d.name; recs.append(r)
            except Exception: pass
    return recs

def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(STACK))
    try:
        from inqview.visualisation import style; style.apply_theme()
    except Exception: pass
    return plt

def fig_eps_vs_E_by(recs, key, fixed_desc, out, title):
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 4.6))
    groups = {}
    for r in recs:
        groups.setdefault(round(r.get(key, 0), 3), []).append(r)
    for g in sorted(groups):
        pts = sorted(groups[g], key=lambda r: r.get("E_eV", 0))
        ax.plot([p["E_eV"] for p in pts], [p["epsilon"]*100 for p in pts], "o-",
                label=f"{key}={g}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("E (eV)"); ax.set_ylabel(r"$\varepsilon$ (%)")
    ax.set_title(title); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)
    return out

def email_phase(subject, body, fig, prev):
    sys.path.insert(0, str(STACK))
    from inqview.email import send_run_email
    try:
        mid = send_run_email(subject=subject, body=body, attachments=[str(fig)] if fig else None,
                             to=TO, in_reply_to=prev[0], references=prev[1])
        prev[0] = mid; prev[1] = (prev[1] or []) + [mid]
        print(f"  emailed: {subject}", flush=True)
    except Exception as ex:
        print(f"  EMAIL FAILED ({subject}): {ex}", flush=True)

def best(recs):
    anchor = [r for r in recs if abs(r.get("E_eV", 0) - ANCHOR_E) < 0.5]
    pool = anchor or recs
    return min(pool, key=lambda r: r.get("epsilon", 9))

def main():
    if not BINARY.exists():
        print(f"BINARY missing: {BINARY}", file=sys.stderr); sys.exit(1)
    HYP.mkdir(parents=True, exist_ok=True)
    prev = [None, []]   # [in_reply_to, references]

    # Phase 1 — CAP η-sweep
    run_jobs(jobs_phase1(), "1: CAP η-sweep @ L=20")
    recs = [r for r in collect("cap") if abs(r.get("L_total",0)-ANCHOR_L) < 0.5]
    fig = fig_eps_vs_E_by(recs, "eta_Ha", "L=20", HYP/"phase1_eps_vs_E_by_eta.png",
                          "Two-sided CAP: ε(E) by depth η (L=20)")
    b = best(recs)
    email_phase("[twosided] Phase 1/3 done — CAP ε(E) by depth η (L=20)",
                f"CAP η-sweep complete ({len(recs)} runs). Best @~10eV: η={b.get('eta_Ha')} "
                f"ε={b.get('epsilon',0)*100:.2f}%. ε PROVISIONAL (Task #7).", fig, prev)

    # Phase 2 — CAP L-sweep at η*
    run_jobs(jobs_phase2(), "2: CAP L-sweep @ η=-0.5")
    recs = [r for r in collect("cap") if abs(r.get("eta_Ha",0)-ETA_STAR) < 0.01]
    fig = fig_eps_vs_E_by(recs, "L_total", "η=-0.5", HYP/"phase2_eps_vs_E_by_L_cap.png",
                          "Two-sided CAP (η=-0.5): ε(E) by width L")
    b = best(recs)
    email_phase("[twosided] Phase 2/3 done — CAP ε(E) by width L (η=-0.5)",
                f"CAP L-sweep complete ({len(recs)} runs). Best @~10eV: L={b.get('L_total')} "
                f"ε={b.get('epsilon',0)*100:.2f}%. ε PROVISIONAL (Task #7).", fig, prev)

    # Phase 3 — mask L-sweep
    run_jobs(jobs_phase3(), "3: mask L-sweep")
    mrecs = collect("mask")
    fig = fig_eps_vs_E_by(mrecs, "L_total", "mask", HYP/"phase3_eps_vs_E_by_L_mask.png",
                          "Two-sided mask: ε(E) by width L")
    bm = best(mrecs)
    email_phase("[twosided] Phase 3/3 done — mask ε(E) by width L",
                f"Mask L-sweep complete ({len(mrecs)} runs). Best @~10eV: L={bm.get('L_total')} "
                f"ε={bm.get('epsilon',0)*100:.2f}%. ε PROVISIONAL (Task #7).", fig, prev)

    # auto-build notebook (notebook-making skill) — once at end of batch
    builder = HYP / "build_twosided_report.py"
    if builder.exists():
        print("==> auto-build notebook", flush=True)
        subprocess.run([sys.executable, str(builder)],
                       env={**os.environ, "PYTHONPATH": str(STACK)}, check=False)
    nb = HYP / "twosided_cap_vs_mask_study.ipynb"
    email_phase("[twosided] ALL phases complete — study notebook built",
                f"All 3 phases done; notebook {'built' if nb.exists() else 'BUILD FAILED'} at {nb}. "
                f"CAP best @10eV L={best([r for r in collect('cap') if abs(r.get('eta_Ha',0)-ETA_STAR)<0.01]).get('L_total')}, "
                f"mask best @10eV L={bm.get('L_total')}. ε PROVISIONAL until Task #7.",
                str(nb) if nb.exists() else None, prev)
    print("==> dispatch complete", flush=True)

if __name__ == "__main__":
    main()
