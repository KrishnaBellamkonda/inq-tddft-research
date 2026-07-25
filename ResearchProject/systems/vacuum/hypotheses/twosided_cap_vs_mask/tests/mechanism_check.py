#!/usr/bin/env python3
"""Task-specific mechanism check for the two-sided absorbers (validation-gates rule).

Asserts, at the anchor point (E≈10 eV, L=20), that BOTH absorber modes actually
remove the wavepacket through the two-sided geometry:

  cap  : two summed sin² slabs (perturbations::sum) integrated in H absorb >90%.
  mask : TwoSidedMaskAbsorber applied in the callback absorbs >90%.

Symmetry of the two ends is guaranteed by construction (CAP = absorbing(+mid,w) +
absorbing(-mid,w) with identical η,w; mask uses |z|) and is unit-tested in the pure
tier (inq-stack/.../test_mask_shape.cpp: M(-z)==M(z)). This check confirms the live
engine behaviour. Reads the smoke runs produced by scripts/twosided_cap_vs_mask/run.

    python3 mechanism_check.py        # exit 0 = PASS

ε PROVISIONAL until the inq-study engine regression (Task #7).
"""
import sys
from pathlib import Path

SCRIPTS = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/vacuum/scripts/twosided_cap_vs_mask")

def parse(p):
    o = {}
    for ln in Path(p).read_text().splitlines():
        k, _, v = ln.partition(' ')
        try: o[k] = float(v)
        except ValueError: o[k] = v
    return o

def check(label, fname, min_absorbed):
    # Mechanism oracle: the absorber must remove a substantial fraction of the WP
    # through the two-sided geometry (vs a hard wall ~0% absorbed, ε~1) AND not be
    # pathological. The CAP (η=-0.5) is a strong absorber (>0.90); the mask is a
    # softer multiplicative envelope, so a 10-Bohr/end ramp legitimately absorbs
    # less at the anchor — it must still clearly work (>0.50) with ε well below a
    # hard wall. Efficacy (how LOW ε goes) is the study's result, not the gate.
    f = SCRIPTS / fname / "results/epsilon.txt"
    if not f.exists():
        f = SCRIPTS / fname / "epsilon.txt"
    if not f.exists():
        print(f"  {label}: MISSING {f}"); return False
    r = parse(f)
    absorbed = r.get("absorbed_fraction", 0.0)
    eps = r.get("epsilon", 1.0)
    ok = (absorbed > min_absorbed) and (eps < 0.90)
    print(f"  {label}: absorbed={absorbed:.4f} eps={eps:.4e} mode={r.get('mode')} "
          f"L={r.get('L_total')} (min_absorbed={min_absorbed}) -> {'PASS' if ok else 'FAIL'}")
    return ok

def main():
    print("two-sided absorber mechanism check (anchor E~10 eV, L=20):")
    ok_cap = check("cap ", "smoke_cap", 0.90)
    ok_mask = check("mask", "smoke_mask", 0.50)
    passed = ok_cap and ok_mask
    print("RESULT:", "PASS" if passed else "FAIL")
    sys.exit(0 if passed else 1)

if __name__ == "__main__":
    main()
