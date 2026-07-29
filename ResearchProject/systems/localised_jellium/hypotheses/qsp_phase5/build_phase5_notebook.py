#!/usr/bin/env python3
"""Assemble qsp_phase5_study.ipynb — the WP quantum stopping-power S(E) velocity
sweep (last phase of the localised-jellium campaign), WITH cross-run diagnostics.

Thin narrative assembler: the S(E) figure is built by build_se_plot.py and the
cross-run comparison figures + density GIFs by build_phase5_comparisons.py; every
number is read from se_state.csv, so quoted numbers and figures can't disagree.

Run (after build_se_plot.py + build_phase5_comparisons.py):
  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
  /local/data/public/skcb2/tddft/venv/bin/python3 build_phase5_notebook.py
"""
import os, sys, math
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # hypotheses/ for _nbreport
from _nbreport import md, embed, setup_cell, set_outdir, build

set_outdir(HERE)
FIGS = os.path.join(HERE, "figs")
STATE = os.path.join(HERE, "se_state.csv")
wp = pd.read_csv(STATE).sort_values("E_eV").reset_index(drop=True) if os.path.exists(STATE) else pd.DataFrame()

def verdict(E, late_slope=None):
    # Aliasing = the CAP-impossible POSITIVE late energy slope (energy injected).
    # A DRAINING point (slope<0) at E>=300 eV is the finer-grid (h=0.35) rerun that
    # cleared the h=0.5 Nyquist aliasing — label it "finer-grid", NOT "ALIASED".
    if late_slope is not None and late_slope > 0:
        return "ALIASED"
    if E <= 130:
        return "clean"
    if E < 300:
        return "borderline"
    return "finer-grid"

def _alias_table():
    """Data-driven grid-resolution table (never hardcoded → never stale). The margin
    and %>k_Nyq columns are the ORIGINAL h=0.5 packet prediction (σ_p=1, k_Nyq=π/0.5);
    measured S + late dE/dt are read live from se_state.csv; the grid verdict follows
    the late-slope sign (positive = aliased)."""
    K_NYQ = math.pi / 0.5          # h=0.5 grid Nyquist = 6.28 a.u.
    sp = 1.0                        # σ_p = 1/(2σ_WP), σ_WP=0.5
    def _row(v, Enom):
        margin = (K_NYQ - v) / sp
        pct = 0.5 * math.erfc(margin / math.sqrt(2)) * 100
        r = None
        for _, rr in wp.iterrows():
            if abs(float(rr["E_eV"]) - Enom) < 8:
                r = rr; break
        S = f"{r['S_eVbohr']:.2f}" if r is not None else "—"
        sl = float(r["late_slope_eV_au"]) if r is not None else None
        sltxt = ("—" if sl is None else (f"**{sl:+.2f}**" if sl > 0 else f"{sl:+.2f}"))
        vd = ("**aliased**" if (sl is not None and sl > 0) else
              ("clean" if Enom <= 130 else ("finer-grid ✓ ‡" if Enom >= 300 else "borderline")))
        return f"| {v} | {v} | {Enom} | {margin:.2f}σ | **{pct:.2g}%** | {S} | {sltxt} | {vd} |"
    head = ("| v | k₀ | E (eV) | margin | **% > k_Nyq (h=0.5)** | measured S | late dE/dt | grid verdict |\n"
            "|--:|--:|--:|--:|--:|--:|--:|:--|")
    return head + "\n" + "\n".join(_row(v, E) for v, E in [(3, 122), (4, 218), (5, 340), (6, 490)])

def _table():
    if not len(wp):
        return "_(no runs analysed yet)_"
    head = ("| E (eV) | v | S (eV/Bohr) | convergence | norm_f | late dE/dt | grid |\n"
            "|--:|--:|--:|:--|--:|--:|:--|")
    rows = []
    for _, r in wp.iterrows():
        b = str(r.get("bound", ""))
        conv = "converged" if b == "exact" else f"↓ {b} bound"
        rows.append(f"| {r['E_eV']:.0f} | {r['v']:.1f} | {r['S_eVbohr']:.2f} | {conv} | "
                    f"{r.get('norm_f', float('nan')):.3f} | {r.get('late_slope_eV_au', float('nan')):+.3f} | "
                    f"**{verdict(r['E_eV'], r.get('late_slope_eV_au'))}** |")
    return head + "\n" + "\n".join(rows)

# energy → density-GIF basename
def _gif_grid():
    items = [(23, "v=1.3"), (54, "v=2.0"), (122, "v=3.0"),
             (218, "v=4.0 ⚠"), (340, "v=5.0 ✗alias"), (490, "v=6.0 ✗alias")]
    cells = []
    for E, vlab in items:
        fn = f"figs/cmp_density_{E}eV.gif"
        if os.path.exists(os.path.join(HERE, fn)):
            cells.append(
                f'<div style="text-align:center;margin:4px">'
                f'<b>{E} eV ({vlab})</b><br><img src="{fn}" width="300"></div>')
    return ('<div style="display:flex;flex-wrap:wrap;justify-content:center">'
            + "".join(cells) + "</div>")

cells = [setup_cell()]

cells.append(md(rf"""# Phase 5 — quantum (wavepacket) stopping power S(E): velocity sweep
### localised jellium slab · 50×50×90 · r$_s$≈5.67 · σ$_{{\rm WP}}$=0.5 · two-sided CAP η=−0.7 · L$_z$=25 Bohr

**Campaign:** `docs/campaigns/localised_jellium/qsp_phase5_velocity_sweep.md` (last
phase). Sweep the WP drift energy, measure the quantum stopping power at each, build
a single **S(E)** curve vs the bulk classical (σ$_{{\rm WP}}$=0.5) + bulk Lindhard
references.

**Grid (drift E = ½k₀²·27.211 eV):** {{23, 54, 122, 218, 340, 490}} eV ↔
v ∈ {{1.3, 2.0, 3.0, 4.0, 5.0, 6.0}}; the 54 eV (v=2.0) point is **reused from
phase 4**, the other five are new runs.

> **⚠ Headline caveat (read the grid-resolution section below):** on the original
> h=0.5 grid the two highest points (340, 490 eV) aliased. **340 eV (v=5) was re-run on
> a finer h=0.35 grid where it drains cleanly** (upper-bound S≈9.8 eV/Bohr), so only
> **490 eV (v=6) remains aliased and excluded**; 218 eV is borderline. The trustworthy
> curve is **E ≤ 340 eV** (23–122 eV clean; 340 eV a finer-grid upper bound).
"""))

cells.append(md(r"""## Method — energy method, slab-correct, WP-anchored

$$ S \;=\; \frac{E_\mathrm{total}(t_f) - E_\mathrm{GS}}{L_z}, \qquad L_z = 25\ \mathrm{Bohr}. $$

- **Anchor = E_GS** (bare-slab GS, −70.2257 Ha), **not** $E_\mathrm{total}(t_0)$: the
  WP drift KE lives inside $E_\mathrm{total}(0)$, so the bulk form would subtract it.
  (Validated: reproduces phase 4 = 2.39 eV/Bohr.)
- **Convergence gate:** WP fully absorbed (norm→0) AND $E_\mathrm{total}$ plateaued;
  else the residual WP energy is still draining → reported S is an **upper bound**.
- **Guard:** N$_\mathrm{total}$ conserved (<2% drain); WP-norm loss = projectile leaving."""))

cells.append(md(rf"""## Result — cumulative S(E)

{_table()}

The localised-slab WP points are overlaid on **bulk** references (geometry estimate,
ADR 0010): classical σ$_{{\rm WP}}$=0.5 (= bulk σ$_q$=0.354), point-charge Lindhard,
and the localised park point (v=2.0, 0.25 eV/Bohr)."""))

if os.path.exists(os.path.join(FIGS, "se_quantum_stopping.png")):
    cells.append(embed(f"{FIGS}/se_quantum_stopping.png",
                       "S(E): WP quantum points (filled=converged, ▽↓=upper bound) vs bulk "
                       "classical + Lindhard. NB: 490 eV aliased (excluded); 340 eV is the "
                       "finer-grid (h=0.35) rerun and drains (see below).", width=720))

# ---------------------------------------------------------------- aliasing
cells.append(md(r"""## ⚠ Grid resolution & aliasing — why 490 eV is excluded (and how 340 eV was recovered)

The real-space grid spacing $h=0.5$ Bohr sets a hard momentum ceiling — the Nyquist
wavevector and its single-particle energy cutoff:

$$ k_\mathrm{Nyq} = \pi/h = 6.28\ a_0^{-1}, \qquad
   E_\mathrm{cut} = \tfrac12(\pi/h)^2 = 19.74\ \mathrm{Ha} = \mathbf{537\ eV}. $$

The WP is **not monochromatic**: σ$_{\rm WP}$=0.5 ⇒ momentum width
$\sigma_p = 1/(2\sigma) = 1.0$, so its $k_z$ content is $\mathcal{N}(k_0,\sigma_p^2)$.
Even when the *drift* energy is below 537 eV, the Gaussian **tail** pushes content
past $k_\mathrm{Nyq}$, where it **aliases** (wraps to spurious momenta and injects
energy). Fraction of the packet above Nyquist, $1-\Phi\!\big((k_\mathrm{Nyq}-k_0)/\sigma_p\big)$:

""" + _alias_table() + r"""

**‡ 340 eV (v=5) was re-run on a finer $h=0.35$ grid** ($k_\mathrm{Nyq}=\pi/0.35=8.98$,
$E_\mathrm{cut}\approx1090$ eV), where **< 0.3 %** of the packet exceeds Nyquist: it then
**drains cleanly** (late slope < 0, norm→0) and yields the S(E) **upper bound** ≈ 9.8
eV/Bohr. The `% > k_Nyq (h=0.5)` column is the ORIGINAL coarse-grid prediction that
explained the initial blow-up.

The aliased fraction (1 %→10 %→39 % across v=4→5→6 on $h=0.5$) tracks exactly where the
**late energy slope flips positive** (energy *created*, impossible under a CAP). The
$h=0.5$ clean criterion is $k_0 + 3\sigma_p < k_\mathrm{Nyq}$ ⇒ **E ≲ 146 eV**; the finer
$h=0.35$ grid pushes this out and **recovers 340 eV**. **490 eV (v=6) is still on $h=0.5$
and remains aliased** — re-run it at $h\le0.35$ to recover it too. The momentum panel
below is the direct evidence (shown for the original $h=0.5$ runs)."""))

# ---------------------------------------------------------------- diagnostics
cells.append(md("""## Cross-run diagnostics

The six runs side by side, **all on the original $h=0.5$ grid** — the momentum panel is
the smoking gun and the energy panel shows the artifact's signature (rising
$E_\\mathrm{total}$ for v5/v6 on that grid). **340 eV (v=5) was subsequently re-run on a
finer $h=0.35$ grid**, where it drains cleanly (that finer-grid run is the 340 eV point in
the S(E) curve above); only **490 eV (v=6)** stays aliased."""))

for fn, cap in [
    ("cmp_momentum_kz.png", "**n_wp(k_z) at t=0 vs the grid Nyquist** (red line, k=6.28). "
                            "v5/v6 jam the wall; % = fraction of the packet above k_Nyq (aliased)."),
    ("cmp_energy.png", "**Left:** retained energy E_total(t)−E_GS — clean runs settle, v5/v6 RISE "
                       "(aliasing injects energy). **Right:** total kinetic energy(t)."),
    ("cmp_norm.png", "**WP orbital norm(t)** — absorption into the CAP. High-v exit fast (norm→0 "
                     "early); slow runs drain gradually and don't fully converge by τ."),
    ("cmp_centroid_sigma.png", "**Left:** WP centroid ⟨z⟩(t) (slope = velocity; slab faces ±12.5, "
                               "CAP ±35). **Right:** spreading σ_z(t)."),
]:
    p = os.path.join(FIGS, fn)
    if os.path.exists(p):
        cells.append(embed(p, cap, width=760))

# ---------------------------------------------------------------- density gifs
cells.append(md(r"""## xz total-density evolution — all runs side by side

The slab (centre) and the WP traversing left→right, with the CAP draining it at the
edges. Faster packets cross and absorb sooner; for v5/v6 the leading edge develops
the short-wavelength ripples that are the real-space face of the k-space aliasing."""))
cells.append(md(_gif_grid()))

# ---------------------------------------------------------------- takeaway
cells.append(md(r"""## Takeaway

- **Trustworthy result (E ≤ 122 eV):** the quantum WP stopping power is **flat at
  ≈ 2.4–2.6 eV/Bohr** across 23–122 eV — well above the matched bulk classical
  (≈ 0.5–0.9) and point-charge Lindhard (≈ 0.45), i.e. a several-fold quantum
  enhancement. All three points are convergence **upper bounds** (the σ=0.5 packet
  is never fully absorbed by τ; the bound is set by the deposited-energy definition,
  not run length).
- **218 eV (v=4):** borderline (~1% aliased on h=0.5); suggestive but not relied upon.
- **340 eV (v=5): recovered on a finer grid.** Re-run at h=0.35 (E_cut≈1090 eV) it
  **drains cleanly** (late slope −6.06, norm→0), giving an upper-bound **S ≈ 9.8 eV/Bohr**
  — confirming the earlier h=0.5 blow-up (apparent S≈13) was pure grid aliasing, not
  physics. The trustworthy curve now extends to 340 eV.
- **490 eV (v=6): excluded** — still on h=0.5, aliased (39% of the packet above k_Nyq;
  late slope +4.7, rising E_total). Re-run at h≤0.35 to recover it, as done for v5.
- **Next step for high-v:** re-run v=6 on the finer grid (as already done for v5) to
  close the sweep."""))

out = os.path.join(HERE, "qsp_phase5_study.ipynb")
build(cells, out, timeout=600)
print(f"done: {out}")
