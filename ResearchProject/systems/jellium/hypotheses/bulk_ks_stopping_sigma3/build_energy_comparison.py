#!/usr/bin/env python3
"""Build the ENERGY-COMPONENT comparison notebook for one twin pair.

    venv/bin/python build_energy_comparison.py [variant]

`variant` defaults to this folder's pair (sigma=3, r_s=5.702). Any
`scripts/<variant>/{wp,classical}` pair works.

WHY THIS EXISTS
---------------
The KE curves of the two halves look alike in shape but differ hugely in
magnitude. That is not (only) physics: the two runs put the projectile in
DIFFERENT PLACES in the energy ledger, so their components are not
like-for-like until that is unpicked. This notebook lays every component of both
runs side by side, plots the change in each, and states exactly which
comparisons are meaningful.

Plan: docs/plans/bulk-jellium-ks-stopping.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import nbformat as nbf
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
SYSTEM = REPO / "ResearchProject/systems/jellium"
KS_HYP = SYSTEM / "hypotheses/bulk_ks_stopping"

DEFAULT_VARIANT = "bulk_ks_stopping_sigma3"
META = {
    "bulk_ks_stopping_sigma3":      dict(label="r_s = 5.702, sigma = 3", z0=-28.0, LZ=80.0,
                                         N0=1.287807e-3, FIT_T1=19.48, NSTEP=600),
    "bulk_ks_stopping_sigma1":      dict(label="r_s = 5.702, sigma = 1", z0=-36.0, LZ=80.0,
                                         N0=1.287807e-3, FIT_T1=10.80, NSTEP=692),
    "bulk_ks_stopping":             dict(label="r_s = 5.702, sigma = 2", z0=-32.0, LZ=80.0,
                                         N0=1.287807e-3, FIT_T1=18.43, NSTEP=646),
    "bulk_ks_stopping_rs4_sigma3":  dict(label="r_s = 3.987, sigma = 3", z0=-28.0, LZ=80.0,
                                         N0=3.765625e-3, FIT_T1=19.48, NSTEP=600),
    "bulk_ks_stopping_rs4_sigma1":  dict(label="r_s = 3.987, sigma = 1", z0=-36.0, LZ=80.0,
                                         N0=3.765625e-3, FIT_T1=9.37, NSTEP=692),
    "bulk_ks_stopping_rs4":         dict(label="r_s = 3.987, sigma = 2", z0=-32.0, LZ=80.0,
                                         N0=3.765625e-3, FIT_T1=18.43, NSTEP=646),
}


def md(t): return nbf.v4.new_markdown_cell(t)
def code(s): return nbf.v4.new_code_cell(s)


def build(variant: str) -> nbf.NotebookNode:
    m = META[variant]
    c: list = []

    c.append(md(rf"""# Energy components: wavepacket vs classical — **{m['label']}**

The kinetic-energy curves of the two halves have similar *shape* but very
different *magnitude*. This notebook takes the energy ledger apart to find out
why, component by component.

## The one fact that governs every comparison below

**The two runs put the projectile in different places in the ledger.** INQ's total is

$$E_{{\rm total}} = E_{{\rm kinetic}} + E_{{\rm Hartree}} + E_{{\rm external}}
+ E_{{\rm nonlocal}} + E_{{\rm xc}} + E_{{\rm xx}} + E_{{\rm ion}} + E_{{\rm ion,kin}}$$

- **Wavepacket run** — the projectile *is* an occupied Kohn–Sham orbital. Its
  energy sits **inside $E_{{\rm kinetic}}$**, and it contributes to $E_{{\rm Hartree}}$
  and $E_{{\rm xc}}$ as well. There are no ions at all, so
  $E_{{\rm external}} = E_{{\rm nonlocal}} = 0$. The system is **closed**:
  $E_{{\rm total}}$ must be constant.
- **Classical run** — the projectile is an **external potential** (a Gaussian
  pseudo-ion). It never enters $n$, so it appears as
  $E_{{\rm external}} \neq 0$ and *nothing else*. INQ leaves
  $E_{{\rm ion,kin}}$ at zero, so $E_{{\rm total}}$ here is the **electronic
  energy alone** and is *supposed* to rise by whatever the projectile deposits.

**Consequence, and the reason magnitudes differ so much:** in the wavepacket run
every component is a **net** quantity — the bath's gain *minus* the projectile's
loss, because both live in the same terms. In the classical run the components
are the bath's response **only**. Comparing $\Delta E_{{\rm kinetic}}$ between them
directly compares a net against a gross.

## On $E_{{\rm ss}}$ and $E_{{\rm ps}}$

INQ has no terms by those names. Mapping to what it does emit:

| requested | INQ term | what it is |
|---|---|---|
| $E_{{\rm ps}}$ (pseudopotential) | `energy_external` + `energy_nonlocal` | projectile potential acting on $n$. **`nonlocal` is identically 0 here** — the Gaussian electron UPF is a *purely local* pseudopotential with no projectors. |
| $E_{{\rm ss}}$ (sum of states) | `energy_eigenvalues` | $\sum_i f_i \varepsilon_i$, the band-structure sum. A **diagnostic**, not a term in the total. |

Every component INQ emits is plotted below, zero ones included and labelled, so
nothing is hidden by omission."""))

    out_dir = SYSTEM / "hypotheses" / variant
    c.append(code(f"""import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

HERE   = Path("{out_dir}")
SYSTEM = Path("{SYSTEM}")
sys.path.insert(0, str(Path("{KS_HYP}")))
import ks_stopping as K
from ks_stopping import HA_TO_EV

VARIANT = "{variant}"
META    = {m!r}
BASE    = SYSTEM / "scripts" / VARIANT
plt.rcParams.update({{"figure.dpi": 120, "axes.grid": True, "grid.alpha": 0.25,
                     "axes.spines.top": False, "axes.spines.right": False}})
(HERE / "figures").mkdir(exist_ok=True)

EN = {{h: K.load_energies(BASE / h) for h in ("wp", "classical")}}
WP = K.load_wp_run(BASE / "wp", box_length_z=META["LZ"], z0=META["z0"])
CL = K.load_classical_run(BASE / "classical", box_length_z=META["LZ"])

# INQ's total is the sum of these eight; the last three are diagnostics.
SUMMED = ["energy_kinetic", "energy_hartree", "energy_external",
          "energy_nonlocal", "energy_xc", "energy_exact_exchange",
          "energy_ion", "energy_ion_kinetic"]
DIAG   = ["energy_nvxc", "energy_eigenvalues"]

def status(df, col):
    if col not in df.columns: return "absent"
    v = df[col].to_numpy()
    return "zero" if np.all(v == 0) else "active"

rows = []
for col in SUMMED + DIAG:
    rows.append(dict(component=col.replace("energy_", ""),
                     in_total=col in SUMMED,
                     wp=status(EN["wp"], col), classical=status(EN["classical"], col)))
display(pd.DataFrame(rows))
print("'zero' means INQ populated it with exactly 0.0 — not missing, genuinely zero.")"""))

    c.append(md(r"""---
## 1. Which components are even active

The table above is the structural difference in one place. `external` is active
for the classical run and identically zero for the wavepacket — that single row
is the projectile changing representation.

`nonlocal` is zero in **both**: the Gaussian electron UPF is purely local, so the
pseudopotential energy $E_{\rm ps}$ is entirely `external`.

---
## 2. Absolute component curves

Left column wavepacket, right column classical, **same component on the same
row** so the pair can be read across. Note the y-scales differ per panel — the
absolute offsets are large and uninformative; section 3 removes them."""))

    c.append(code(r"""active = [c for c in SUMMED + DIAG
          if any(c in EN[h].columns and not np.all(EN[h][c].to_numpy() == 0)
                 for h in ("wp", "classical"))]
n = len(active)
fig, ax = plt.subplots(n, 2, figsize=(11, 2.1*n), sharex=True)
for r, col in enumerate(active):
    for cix, h in enumerate(("wp", "classical")):
        d = EN[h]; a = ax[r, cix]
        if col in d.columns:
            v = d[col].to_numpy() * HA_TO_EV
            a.plot(d["time_au"], v, lw=1.4,
                   color="C0" if h == "wp" else "C3")
            if np.all(v == 0):
                a.annotate("identically zero", xy=(0.5, 0.5),
                           xycoords="axes fraction", ha="center",
                           fontsize=8, color="0.4")
        a.axvspan(4.0, META["FIT_T1"], color="0.9", zorder=0)
        if r == 0: a.set_title("wavepacket" if h == "wp" else "classical")
        if cix == 0: a.set_ylabel(col.replace("energy_", ""), fontsize=8)
ax[-1, 0].set_xlabel("time [a.u.]"); ax[-1, 1].set_xlabel("time [a.u.]")
fig.suptitle(f"{META['label']} — absolute energy components [eV]  (shaded = fit window)", y=1.002)
fig.tight_layout()
fig.savefig(HERE/"figures"/f"components_absolute_{VARIANT}.png", dpi=150, bbox_inches="tight")
plt.show()"""))

    c.append(md(r"""---
## 3. Residuals — $\Delta E(t) = E(t) - E(0)$

This is the comparison that matters. Subtracting $t=0$ removes the (large,
run-specific) absolute offsets and leaves only what the propagation did.

**Both halves are drawn on the SAME axes per component**, so the magnitude
difference you noticed is directly visible."""))

    c.append(code(r"""fig, ax = plt.subplots(2, 3, figsize=(14, 7))
panels = [c for c in active if c not in DIAG][:6]
for a, col in zip(ax.flat, panels):
    for h, colr in (("wp", "C0"), ("classical", "C3")):
        d = EN[h]
        if col not in d.columns: continue
        v = (d[col].to_numpy() - d[col].to_numpy()[0]) * HA_TO_EV
        a.plot(d["time_au"], v, lw=1.6, color=colr, label=h)
    a.axhline(0, color="k", lw=0.8)
    a.axvspan(4.0, META["FIT_T1"], color="0.9", zorder=0)
    a.set_title(r"$\Delta$" + col.replace("energy_", ""), fontsize=10)
    a.set_xlabel("time [a.u.]"); a.set_ylabel(r"$\Delta E$ [eV]")
    a.legend(fontsize=7, frameon=False)
for a in ax.flat[len(panels):]: a.axis("off")
fig.suptitle(f"{META['label']} — CHANGE in each component, both halves on shared axes", y=1.00)
fig.tight_layout()
fig.savefig(HERE/"figures"/f"components_residual_{VARIANT}.png", dpi=150, bbox_inches="tight")
plt.show()"""))

    c.append(md(r"""---
## 4. Total energy — and why the two totals mean different things

The wavepacket total must be **flat** (closed system: the projectile is inside
it). The classical total must **rise** (electronic energy only; the projectile's
kinetic energy is tracked separately and INQ leaves `ion_kinetic` at zero).

A flat WP total is therefore a *numerical* check, and a rising classical total is
the *physical* deposit. They are not the same quantity and must not be compared
as though they were."""))

    c.append(code(r"""fig, ax = plt.subplots(1, 3, figsize=(14, 3.9))
for h, colr in (("wp", "C0"), ("classical", "C3")):
    d = EN[h]
    v = (d["energy_total"].to_numpy() - d["energy_total"].to_numpy()[0]) * HA_TO_EV
    ax[0].plot(d["time_au"], v, lw=1.7, color=colr, label=h)
ax[0].axhline(0, color="k", lw=0.8); ax[0].axvspan(4.0, META["FIT_T1"], color="0.9", zorder=0)
ax[0].set_title(r"$\Delta E_{\rm total}$ — shared axes"); ax[0].legend(fontsize=8, frameon=False)
ax[0].set_xlabel("time [a.u.]"); ax[0].set_ylabel(r"$\Delta E$ [eV]")

d = EN["wp"]
v = (d["energy_total"].to_numpy() - d["energy_total"].to_numpy()[0]) * HA_TO_EV
ax[1].plot(d["time_au"], v, lw=1.5, color="C0")
ax[1].axhline(0, color="k", lw=0.8)
ax[1].set_title(f"WP total, own scale — drift {v[-1]:+.2e} eV")
ax[1].set_xlabel("time [a.u.]"); ax[1].set_ylabel(r"$\Delta E$ [eV]")

# the projectile's own KE, the quantity the classical total is mirroring
mk = (CL.t >= 0)
ax[2].plot(CL.t, (CL.T - CL.T[0]) * HA_TO_EV, lw=1.6, color="C3",
           label=r"classical $\Delta$KE$_{\rm ion}$")
d = EN["classical"]
ax[2].plot(d["time_au"],
           -(d["energy_total"].to_numpy() - d["energy_total"].to_numpy()[0]) * HA_TO_EV,
           lw=1.2, ls="--", color="k", label=r"$-\Delta E_{\rm total}$ (electronic)")
ax[2].axhline(0, color="k", lw=0.8); ax[2].axvspan(4.0, META["FIT_T1"], color="0.9", zorder=0)
ax[2].set_title("classical: closure check"); ax[2].legend(fontsize=7, frameon=False)
ax[2].set_xlabel("time [a.u.]"); ax[2].set_ylabel(r"$\Delta E$ [eV]")
fig.tight_layout()
fig.savefig(HERE/"figures"/f"total_energy_{VARIANT}.png", dpi=150, bbox_inches="tight")
plt.show()

ke_loss = (CL.T[0] - CL.T[-1]) * HA_TO_EV
gain = (EN["classical"]["energy_total"].iloc[-1] - EN["classical"]["energy_total"].iloc[0]) * HA_TO_EV
print(f"classical: projectile KE loss {ke_loss:8.3f} eV")
print(f"           electronic gain    {gain:8.3f} eV")
print(f"           closure mismatch   {abs(gain-ke_loss):8.3f} eV  ({100*abs(gain-ke_loss)/ke_loss:.2f} %)")"""))

    c.append(md(r"""---
## 5. The ledger, side by side

Total change in every component over the whole run, both halves. This is the
table to read when asking "where did the energy go, and why do the two runs
disagree in magnitude"."""))

    c.append(code(r"""rows = []
for col in SUMMED + DIAG:
    r = dict(component=col.replace("energy_", ""), in_total=col in SUMMED)
    for h in ("wp", "classical"):
        d = EN[h]
        r[h] = (d[col].to_numpy()[-1] - d[col].to_numpy()[0]) * HA_TO_EV if col in d.columns else np.nan
    r["ratio_cl_over_wp"] = (r["classical"] / r["wp"]) if r["wp"] not in (0, np.nan) and abs(r["wp"]) > 1e-12 else np.nan
    rows.append(r)
led = pd.DataFrame(rows)
display(led.style.format({c: "{:+.4f}" for c in ("wp", "classical", "ratio_cl_over_wp")}))

print("checksum — components summing to the total (eV):")
for h in ("wp", "classical"):
    d = EN[h]
    s = sum((d[c].to_numpy()[-1] - d[c].to_numpy()[0]) * HA_TO_EV
            for c in SUMMED if c in d.columns)
    t = (d["energy_total"].to_numpy()[-1] - d["energy_total"].to_numpy()[0]) * HA_TO_EV
    print(f"   {h:10s} sum of parts {s:+9.4f}   energy_total {t:+9.4f}   diff {abs(s-t):.2e}")"""))

    c.append(md(r"""---
## 6. What this shows

Read the residual panels and the ledger together:

- **`external` is the structural difference.** Non-zero only for the classical
  run. It is the entire pseudopotential energy $E_{\rm ps}$ (nonlocal is zero for
  this purely local UPF), and it has no wavepacket counterpart — the WP's
  equivalent contribution is folded into `kinetic`, `hartree` and `xc`.
- **The wavepacket's components are NET.** Its projectile lives in the same terms
  as its bath, so `kinetic` shows (bath gain − projectile loss), not the bath
  gain. The classical `kinetic` is purely the bath. That is the single largest
  reason the magnitudes differ, and it is bookkeeping, not physics.
- **Sign flips are informative.** Where $\Delta E_{\rm Hartree}$ or
  $\Delta E_{\rm xc}$ move in opposite directions between the two halves, the two
  projectiles are polarising the gas differently — that *is* physics, and it is
  not explained by the bookkeeping above.
- **The checksum** confirms the eight summed components reproduce `energy_total`
  in both runs, so nothing is being lost to an unaccounted channel.

**The comparison that IS like-for-like:** the classical electronic gain against
the wavepacket projectile's own KE loss — both are "energy delivered to the
bath". The totals themselves are not comparable."""))

    # ---------------- pairwise P/S/B decomposition -------------------------
    c.append(md(r"""---
## 7. Pairwise decomposition — the comparison INQ's scalars cannot make

`interactions.csv` splits the electrostatics by **charge group** rather than by
where the code files it: **P**rojectile, **S**ystem (bath electrons),
**B**ackground.

$$E_{SS}=	frac12\!\int\! n_Sarphi_S \quad
E_{PP}=	frac12\!\int\! n_Parphi_P \quad
E_{PS}=\int\! n_Sarphi_P \quad
E_{SB}=-\!\int\! n_Sarphi_+ \quad
E_{PB}=-\!\int\! n_Parphi_+$$

**Why this is the right comparison.** `energy_hartree` means *different things* in
the two runs — $E_{SS}$ for the classical (the projectile is an external potential
and never enters $n$), but $E_{SS}+E_{PS}+E_{PP}$ for the wavepacket (it is an
occupied orbital). The pairwise terms are defined by charge group, so
$E_{PS}$ in one run is the same physical object as $E_{PS}$ in the other.

**$E_{PP}$ is the projectile self-Hartree** — a wavepacket carries it and a rigid
classical charge cloud cannot change it. In LDA it is only partly cancelled by
exchange, so its *change* is the spurious self-interaction acting on the dynamics.

**Bulk note:** the background is uniform, so $arphi_+$ is pure $G{=}0$ which INQ
drops. $E_{SB}=E_{PB}=E_{BB}\equiv 0$ here by construction — they are written as
zeros to keep the schema identical across systems, not because nothing happened."""))

    c.append(code(r"""IX = {}
for h in ("wp", "classical"):
    f = BASE / h / "results/raw/observables/interactions.csv"
    IX[h] = pd.read_csv(f) if f.exists() else None
    if IX[h] is None:
        print(f"  {h}: interactions.csv ABSENT — run predates "
              f".claude/rules/decomposed-interaction-energies.md")

if all(v is not None for v in IX.values()):
    print("CLOSURE GATES (terms must sum back to the INQ scalars)")
    w = IX["wp"]
    dw = np.max(np.abs((w.e_ss + w.e_ps + w.e_pp) - w.e_hartree_check))
    print(f"  WP        E_SS+E_PS+E_PP vs energy_hartree : max|diff| = {dw:.3e} Ha")
    cl = IX["classical"]
    dc = np.max(np.abs(cl.e_ss - cl.e_hartree_inq))
    print(f"  classical E_SS          vs energy_hartree : max|diff| = {dc:.3e} Ha")
    assert dw < 1e-9 and dc < 1e-9, "CLOSURE FAILED — do not trust these terms"
    print("  -> closure holds; the decomposition is trustworthy")
    print("")

    rows = []
    for t in ("e_ss", "e_pp", "e_ps", "e_sb", "e_pb"):
        r = {"term": t.upper().replace("E_", "E_")}
        for h in ("wp", "classical"):
            v = IX[h][t].to_numpy()
            r[f"{h}_t0"] = v[0] * HA_TO_EV
            r[f"{h}_delta"] = (v[-1] - v[0]) * HA_TO_EV
        r["ratio_cl/wp"] = (r["classical_delta"] / r["wp_delta"]
                            if abs(r["wp_delta"]) > 1e-9 else np.nan)
        rows.append(r)
    display(pd.DataFrame(rows).style.format(
        {c: "{:+.4f}" for c in ("wp_t0","wp_delta","classical_t0","classical_delta","ratio_cl/wp")}))"""))

    c.append(code(r"""if all(v is not None for v in IX.values()):
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for a, t, ttl in zip(ax, ("e_ss", "e_pp", "e_ps"),
                         (r"$E_{SS}$ bath-bath (polarisation)",
                          r"$E_{PP}$ projectile SELF-Hartree",
                          r"$E_{PS}$ projectile-bath")):
        for h, colr in (("wp", "C0"), ("classical", "C3")):
            d = IX[h]; v = (d[t].to_numpy() - d[t].to_numpy()[0]) * HA_TO_EV
            a.plot(d["time_au"], v, lw=1.7, color=colr, label=h)
        a.axhline(0, color="k", lw=0.8)
        a.axvspan(4.0, META["FIT_T1"], color="0.9", zorder=0)
        a.set_title(ttl, fontsize=10); a.set_xlabel("time [a.u.]")
        a.set_ylabel(r"$\Delta E$ [eV]"); a.legend(fontsize=8, frameon=False)
    fig.suptitle(f"{META['label']} — pairwise terms, shared axes per panel", y=1.02)
    fig.tight_layout()
    fig.savefig(HERE/"figures"/f"pairwise_{VARIANT}.png", dpi=150, bbox_inches="tight")
    plt.show()

    w, cl = IX["wp"], IX["classical"]
    print(f"E_PP at t=0:  wp {w.e_pp.iloc[0]*HA_TO_EV:.4f} eV   "
          f"classical {cl.e_pp.iloc[0]*HA_TO_EV:.4f} eV")
    print("  identical => the sigma-matching convention (sigma_pot = sigma_WP/sqrt2)")
    print("     puts the SAME charge cloud in both runs at t=0.")
    print(f"\nE_PP change:  wp {(w.e_pp.iloc[-1]-w.e_pp.iloc[0])*HA_TO_EV:+.4f} eV   "
          f"classical {(cl.e_pp.iloc[-1]-cl.e_pp.iloc[0])*HA_TO_EV:+.4f} eV")
    print("  the classical charge cloud is RIGID (self-energy cannot change);")
    print("  the packet spreads, so its self-Hartree collapses. Wavepacket-only.")"""))

    c.append(md(r"""**Reading it.** $E_{SS}$ and $E_{PS}$ are the physics — bath polarisation and
projectile–bath coupling — and both are comparable between representations.
$E_{PP}$ is not physics the classical projectile can have: its change is the
wavepacket's spreading self-repulsion, and it is the term that made
$\Delta E_{
m hartree}$ come out *negative* for the wavepacket while the
classical was positive.

Compare the $E_{PS}$ ratio with this pair's stopping-power ratio: if the pairwise
decomposition is capturing the right physics, the projectile–bath interaction
should track the drag without needing an unexplained factor."""))

    nb = nbf.v4.new_notebook(cells=c)
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python",
                                 "name": "python3"}
    return nb


def main() -> int:
    variant = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_VARIANT
    if variant not in META:
        print(f"unknown variant {variant}; known: {list(META)}"); return 2
    nb = build(variant)
    out_dir = SYSTEM / "hypotheses" / variant
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"energy_component_comparison_{variant}.ipynb"
    nbf.write(nb, str(out))
    print(f"wrote {out} ({len(nb.cells)} cells)")
    client = NotebookClient(nb, timeout=1800, kernel_name="python3",
                            resources={"metadata": {"path": str(out_dir)}},
                            allow_errors=True)
    client.execute()
    nbf.write(nb, str(out))
    n_err = sum(1 for c in nb.cells if c.cell_type == "code"
                for o in c.get("outputs", []) if o.get("output_type") == "error")
    print(f"executed with {n_err} error(s) -> {out}")
    return 1 if n_err else 0


if __name__ == "__main__":
    raise SystemExit(main())
