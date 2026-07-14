#!/usr/bin/env python3
"""pbc_open_z orchestrator — Arm B: electrostatic z-periodicity vs the energy oscillation.

Campaign: docs/campaigns/localised_jellium/pbc-open-z-oscillation.md
Serial chain on ONE GPU (user requirement), idempotent resume (steps whose
run_summary shows run_completed=true are skipped), per-step try/except with
full traceback into the log; a failed gate BLOCKS the dependent RT runs
(correctness gate, allowed to block; cost never blocks).

Chain:
  0. p2 GS for slab_n52 (build+run via inq-run)  -> shared_gs/..._per2
  1. GS sanity gate: n(z) slab-shaped, tail below threshold at the CAP footprint
  2. rebuild cap_fix binary + 5-step p2 smoke (EM_PERIODICITY knob added)
  3. p2 twins of the three key p3 witnesses (same locked cap_fix harness):
       p2_two_eta0p2_700       (twin of run01: the 23.5 eV riser)
       p2_two_eta1p0_950       (twin of run06: the +31 eV above-zero riser)
       p2_wrap_eta2p0_w40_950  (twin of run11: the monotone winner)
  4. comparison table p2 vs p3 -> hypotheses/pbc_open_z/comparison.md
"""
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent                    # scripts/pbc_open_z
SYS = HERE.parent.parent                                  # systems/localised_jellium
CAPFIX = SYS / "scripts" / "cap_fix"
HYP = SYS / "hypotheses" / "pbc_open_z"
GS_DIR = SYS / "shared_gs" / "slab_n52_L40x40x80_dx0p333_per2"
VENVPY = "/local/data/public/skcb2/tddft/venv/bin/python3"
INQ_SOURCE = "/local/data/public/skcb2/tddft/inq-study"
GPU = os.environ.get("PBC_GPU", "0")
HA_EV = 27.211386
LOG = HERE / "orchestrate.log"

# p2 twin runs: (name, env overrides) — witnesses chosen to bracket the ladder
RUNS = [
    ("p2_two_eta0p2_700",
     {"EM_CAP_MODE": "two", "EM_CAP_ETA": "-0.2", "EM_N_STEPS": "700"}),
    ("p2_two_eta1p0_950",
     {"EM_CAP_MODE": "two", "EM_CAP_ETA": "-1.0", "EM_N_STEPS": "950"}),
    ("p2_wrap_eta2p0_w40_950",
     {"EM_CAP_MODE": "wrap", "EM_CAP_ETA": "-2.0", "EM_WRAP_WIDTH_BOHR": "40",
      "EM_N_STEPS": "950"}),
]
P3_TWINS = {  # existing cap_fix results for the comparison table
    "p2_two_eta0p2_700": "run01_baseline_two_eta0p2",
    "p2_two_eta1p0_950": "run06_poscontrol_eta1p0_950",
    "p2_wrap_eta2p0_w40_950": "run11_wrap_eta2p0_w40_950",
}


def log(msg):
    line = f"[{time.strftime('%F %T')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def sh(cmd, cwd, extra_env=None, timeout=14400):
    env = dict(os.environ, INQ_SOURCE=INQ_SOURCE, CUDA_VISIBLE_DEVICES=GPU)
    if extra_env:
        env.update(extra_env)
    log(f"RUN ({cwd}): {cmd}")
    r = subprocess.run(cmd, shell=True, cwd=cwd, env=env, timeout=timeout,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    tail = "\n".join(r.stdout.splitlines()[-12:])
    log(f"exit={r.returncode}\n{tail}")
    if r.returncode != 0:
        raise RuntimeError(f"command failed (exit {r.returncode}): {cmd}")
    return r.stdout


def completed(summary: Path) -> bool:
    try:
        return "run_completed = true" in summary.read_text()
    except OSError:
        return False


def step0_gs():
    if completed(HERE / "gs" / "results" / "run_summary.txt") and GS_DIR.exists():
        log("step0: p2 GS already complete — skip")
        return
    sh("inq-run", cwd=HERE / "gs", extra_env={"EM_PERIODICITY": "2"})


def step1_gate():
    code = f"""
import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
import numpy as np
from inqview import load_vti
f = load_vti("{HERE}/gs/results/density_gs/density_gs.vti", expect_centered_axis="z")
nz = f.data.mean(axis=(0, 1))          # planar-averaged n(z), physical order
z = f.z
interior = nz[np.abs(z) < 10].mean()
tail = nz[np.abs(z) > 32.5].max()
n0 = 52.0 / 40000.0
print(f"interior mean n(z) = {{interior:.4e}} (n0 = {{n0:.4e}})")
print(f"max tail n(|z|>32.5) = {{tail:.4e}}")
assert 0.7 * n0 < interior < 1.3 * n0, "interior density off n0"
assert tail < 1e-4, "GS tail reaches the CAP footprint"
print("GATE PASS")
"""
    out = sh(f"{VENVPY} - <<'EOF'\n{code}\nEOF", cwd=HERE)
    if "GATE PASS" not in out:
        raise RuntimeError("GS gate did not pass")


def step2_rebuild_smoke():
    smoke = CAPFIX / "results" / "smoke_p2" / "run_summary.txt"
    binary, src = CAPFIX / "run", CAPFIX / "run.cpp"
    if completed(smoke) and binary.exists() and binary.stat().st_mtime > src.stat().st_mtime:
        log("step2: binary fresh + p2 smoke done — skip")
        return
    sh("inq-run", cwd=CAPFIX, extra_env={
        "EM_PERIODICITY": "2", "EM_GS_DIR": str(GS_DIR), "EM_N_STEPS": "5",
        "EM_WRITE_EVERY": "1", "EM_OUT": "smoke_p2", "EM_CAP_ETA": "-0.2"})
    if not completed(smoke):
        raise RuntimeError("p2 smoke did not complete")


def step3_runs():
    for name, envx in RUNS:
        if completed(CAPFIX / "results" / name / "run_summary.txt"):
            log(f"step3: {name} already complete — skip")
            continue
        env = {"EM_PERIODICITY": "2", "EM_GS_DIR": str(GS_DIR), **envx}
        sh(f"./autoresearch.sh {name} {GPU} > results/{name}.log 2>&1",
           cwd=CAPFIX, extra_env=env)


def step4_compare():
    import pandas as pd  # noqa: F401 (venv python does the work below)
    code = f"""
import pandas as pd
HA = {HA_EV}
BASE = "{CAPFIX}/results"
pairs = {P3_TWINS!r}
rows = []
for p2, p3 in pairs.items():
    for tag, run in (("p2", p2), ("p3", p3)):
        df = pd.read_csv(f"{{BASE}}/{{run}}/raw/observables/observables.csv")
        dE = (df.energy_total - df.energy_total.iloc[0]) * HA
        im = int(dE.idxmin())
        rows.append((p2, tag, run, df.time_au[im], dE.min(), dE.iloc[-1] - dE[im],
                     max(0.0, dE.max())))
out = ["| pair | conv | run | t_min | drain (eV) | rise (eV) | excursion (eV) |",
       "|---|---|---|---|---|---|---|"]
for r in rows:
    out.append(f"| {{r[0]}} | {{r[1]}} | {{r[2]}} | {{r[3]:.1f}} | {{r[4]:+.2f}} "
               f"| {{r[5]:+.3f}} | {{r[6]:+.3f}} |")
text = "\\n".join(out)
print(text)
open("{HYP}/comparison.md", "w").write(
    "# Arm B: p2 (open-z electrostatics) vs p3 (fully periodic)\\n\\n" + text + "\\n")
"""
    sh(f"{VENVPY} - <<'EOF'\n{code}\nEOF", cwd=HERE)


def main():
    HYP.mkdir(parents=True, exist_ok=True)
    steps = [("step0_gs", step0_gs), ("step1_gate", step1_gate),
             ("step2_rebuild_smoke", step2_rebuild_smoke),
             ("step3_runs", step3_runs), ("step4_compare", step4_compare)]
    for name, fn in steps:
        log(f"=== {name} ===")
        try:
            fn()
        except Exception:
            log(f"FAILED in {name}:\n{traceback.format_exc()}")
            log("Chain blocked (dependent steps need this one). Fix and re-run — "
                "completed steps resume idempotently.")
            sys.exit(1)
    log("=== ALL STEPS COMPLETE ===")


if __name__ == "__main__":
    main()
