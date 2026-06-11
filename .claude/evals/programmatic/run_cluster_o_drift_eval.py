#!/usr/bin/env python3
"""Runner for .claude/evals/clusters/cluster-o-min-observable-set.md (drift).

Canonical = minimum_observable_set.hpp. This checks that:
  (1) every canonical observable name in the .hpp is COVERED in the spec doc
      (docs/observables/minimum-set-spec.md) — caught directly if the spec uses
      the same key, else via the key->prose bridge below;
  (2) the tddft-simulations skill Phase 3 carries the canonical-reference
      sentinel (it defers to the .hpp, not a second source of truth).

The bridge documents the CURRENT key<->prose mapping (the spec uses prose names
where the .hpp uses keys). A canonical observable with neither a verbatim hit
nor a bridge entry is real drift (code added an observable the spec never
documented). Full cleanup (align the spec's keys to the .hpp) would remove the
bridge — tracked as the remaining Cluster-O item.

Pure stdlib:  python3 .claude/evals/programmatic/run_cluster_o_drift_eval.py
"""
from __future__ import annotations

import os
import re

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
HPP = os.path.join(REPO, "inq-stack/include/inqkit/observables/minimum_observable_set.hpp")
SPEC = os.path.join(REPO, "docs/observables/minimum-set-spec.md")
SKILL = os.path.join(REPO, ".claude/skills/tddft-simulations/SKILL.md")

# canonical key -> regex that matches its representation in the spec doc, for
# the keys the spec currently names in prose rather than by canonical key.
BRIDGE = {
    "gs_eigenvalues": r"GS eigenvalues",
    "gs_occupations": r"GS occupations",
    "gs_system_density": r"GS (orbital )?densit",
    "density_total_rt": r"density .*total|RT density",
    "density_wp_rt": r"density_wp|density \{system,wp,total\}",
    "leed_screen_config": r"LEED|screen_config",
}


def canonical_names():
    txt = open(HPP).read()
    return sorted(set(re.findall(r'(?:csv|vti|text)\("([^"]+)"', txt)))


def main() -> int:
    spec = open(SPEC).read()
    fails = []
    uncovered = []
    for name in canonical_names():
        if name in spec:
            continue                       # verbatim key in the spec
        pat = BRIDGE.get(name)
        if pat and re.search(pat, spec, re.IGNORECASE):
            continue                       # covered via the documented bridge
        uncovered.append(name)

    if uncovered:
        fails.append("canonical observables not covered in the spec doc "
                     "(real drift — add to spec or bridge): " + ", ".join(uncovered))

    skill = open(SKILL).read()
    if "min-obs-set: canonical" not in skill:
        fails.append("tddft-simulations Phase 3 missing the canonical-reference sentinel")

    if fails:
        print(f"FAIL: {len(fails)} drift checks failed")
        for f in fails:
            print("  -", f)
        return 1
    n = len(canonical_names())
    print(f"PASS: Cluster-O drift ({n} canonical observables all covered in spec; "
          f"Phase 3 sentinel present)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
