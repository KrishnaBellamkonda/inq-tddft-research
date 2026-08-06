"""Build the frozen (blinded) experimental-reference pack for the
bilayer-graphene stopping campaign.

Parses the raw NIST ESTAR HTML tables (graphite, matno 906) into a clean
CSV, records the Geelen 2019 quoted anchor values, and renders the
proof-of-extraction reference figure.

Provenance:
- ESTAR: physics.nist.gov/cgi-bin/Star/e_table-t.pl, matno=906
  (CARBON, GRAPHITE, rho=1.7 g/cm3, I=78.0 eV), fetched 2026-08-05.
- Geelen anchors: Geelen et al., PRL 123, 086802 (2019) (arXiv:1904.13152),
  values quoted verbatim from the text.
- CXRO f1/f2: henke.lbl.gov/optical_constants/sf/c.nff (10 eV-30 keV).

BLINDING: outputs live in this folder and are NOT to be consulted during
run design or analysis — Phase 4 comparison only.
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).parent
RHO = 1.7  # g/cm3, ESTAR graphite

def parse_estar(raw_html: str):
    """Rows: E[MeV] S_col S_rad S_tot [MeV cm2/g] CSDA[g/cm2] (+2 more)."""
    text = raw_html.replace("<br>", "\n")
    rows = []
    for line in text.splitlines():
        m = re.findall(r"\d\.\d{3}E[+-]\d{2}", line)
        if len(m) >= 5:
            rows.append([float(x) for x in m[:5]])
    return rows

rows = []
for name in ["estar_graphite_lowE_raw.html", "estar_graphite_default_raw.html"]:
    rows += parse_estar((HERE / name).read_text())
rows = sorted({r[0]: r for r in rows}.values())  # dedupe on energy, sort

with open(HERE / "estar_graphite_electrons.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["E_keV", "S_col_MeVcm2_per_g", "S_rad_MeVcm2_per_g",
                "S_tot_MeVcm2_per_g", "csda_range_g_per_cm2", "S_col_eV_per_A"])
    for e_mev, s_col, s_rad, s_tot, csda in rows:
        # eV/Angstrom = MeV cm2/g * g/cm3 * 1e6 eV/MeV * 1e-8 cm/A
        w.writerow([e_mev * 1e3, s_col, s_rad, s_tot, csda, s_col * RHO * 0.01])

# Geelen 2019 anchors — verbatim-quoted values only (never digitised curves).
geelen = [
    ("lambda_inel", "~0 eV", "approx 3", "graphene layers",
     "'a monotonic decrease from lambda_inel ~ 3 layers near 0 eV'"),
    ("lambda_inel", "25 eV", "approx 1", "graphene layers",
     "'to lambda_inel ~ 1 layer at 25 eV'"),
    ("lambda_el max (2LG)", "0-5 eV", "up to approx 80", "graphene layers",
     "'up to lambda_el ~ 80 layers for 2LG'"),
    ("lambda_el max (3LG/4LG)", "0-5 eV", "approx 20-30", "graphene layers",
     "'lambda_el ~ 20-30 layers for 3LG and 4LG'"),
    ("R max band (2LG+)", "5-15 eV", "reflectivity max / transmission min",
     "-", "above-vacuum interlayer bandgap; ABSENT for 1LG"),
    ("reflection minima count", "0-5 eV", "n-1 for nLG", "-",
     "interlayer transmission resonances"),
]
with open(HERE / "geelen2019_anchor_values.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["quantity", "energy", "value", "units", "verbatim_source_quote"])
    w.writerows(geelen)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    from inqview.visualisation import style  # canonical theme (ADR 0004)
    style.apply()
except Exception:
    pass

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.6))
e_kev = [r[0] * 1e3 for r in rows]
s_evA = [r[1] * RHO * 0.01 for r in rows]
ax1.loglog([e * 1e3 for e in e_kev], s_evA, "-", lw=1.5, label="ESTAR graphite $S_{col}$")
ax1.axvspan(15, 300, alpha=0.15, label="our runs (15–300 eV)")
ax1.set_xlabel("electron energy (eV)")
ax1.set_ylabel(r"$S_{col}$ (eV/$\mathrm{\AA}$)")
ax1.set_title("high-energy anchor: ESTAR (valid ≥ 1 keV)")
ax1.legend(frameon=False, fontsize=8)

en = [0.5, 25.0]
lam = [3.0, 1.0]
ax2.plot(en, lam, "o", ms=7, label=r"Geelen 2019 $\lambda_{inel}$ (quoted)")
ax2.axvspan(15, 25, alpha=0.15, label="overlap with our grid")
ax2.set_xlabel("electron energy above vacuum (eV)")
ax2.set_ylabel(r"$\lambda_{inel}$ (graphene layers)")
ax2.set_ylim(0, 4)
ax2.set_title("low-energy anchor: eV-TEM transmission")
ax2.legend(frameon=False, fontsize=8)
fig.suptitle("FROZEN reference data — do not consult before Phase 4", fontsize=9)
fig.tight_layout()
fig.savefig(HERE / "reference_stopping_figure.png", dpi=180)
print(f"rows: {len(rows)}, E range {min(e_kev):.3g}-{max(e_kev):.3g} keV, "
      f"S_col(1 keV) = {s_evA[0]:.2g} eV/A")
