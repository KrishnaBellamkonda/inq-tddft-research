#!/usr/bin/env python3
"""compare.py <nocap_results> <cap_results> <out_png> — the campaign's headline
figure: total energy(t) for no-CAP vs CAP overlaid, so the PLATEAU GAP (energy
radiated to the boundaries and drained by the CAP) is read directly.
"""
import sys, glob
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HA_EV = 27.211386245988

def load_total(results: Path):
    files = sorted(glob.glob(str(results / "raw" / "observables" / "energies*.csv")))
    df = pd.concat([pd.read_csv(f) for f in files]).drop_duplicates("step").sort_values("step")
    return df["time_au"].to_numpy(), df["total"].to_numpy() * HA_EV

def main() -> int:
    nocap, cap, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    tn, en = load_total(nocap)
    tc, ec = load_total(cap)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(tn, en, "b-", lw=2.0, label="no-CAP (closed box)")
    ax.plot(tc, ec, "r-", lw=2.0, label="CAP (drains escaping flux)")
    e0 = en[0]
    ax.axhline(e0, color="grey", ls=":", lw=1, label="initial total")
    plateau_nocap = np.mean(en[int(0.8*len(en)):])
    plateau_cap = np.mean(ec[int(0.8*len(ec)):])
    gap = plateau_nocap - plateau_cap
    ax.set_xlabel("time (a.u.)"); ax.set_ylabel("total energy (eV)")
    ax.set_title("Energy-total plateau: no-CAP vs CAP")
    ax.legend(fontsize=9)
    ax.text(0.02, 0.03,
            f"plateau(no-CAP)={plateau_nocap:.1f} eV\nplateau(CAP)={plateau_cap:.1f} eV\n"
            f"gap (radiated) = {gap:.1f} eV",
            transform=ax.transAxes, fontsize=9, va="bottom",
            bbox=dict(boxstyle="round", fc="white", alpha=0.8))
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    print(f"plateau no-CAP={plateau_nocap:.2f} eV  CAP={plateau_cap:.2f} eV  gap={gap:.2f} eV -> {out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
