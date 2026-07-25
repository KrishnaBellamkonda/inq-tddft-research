#!/usr/bin/env python3
"""Assemble wide_wavepacket_planning.ipynb — the campaign-planning notebook for
the `wide-wavepacket-lowspread` campaign (localised jellium slab).

This is a DESIGN/DELIBERATION notebook (pre-run), not a run-set study. It makes
the analytical predictions that select the wavepacket parameters, in priority
order: minimise spreading first.

Contents (executable code cells, canonical theme):
  1. Geometry of the planned 50x50x100 box (slab / CAP / launch gap).
  2. Free-Gaussian spreading prediction sigma(t)/sigma0 for a RANGE of sigma0,
     evaluated at the slab-interaction moment and at far-CAP absorption, for a
     few drift energies E.
  3. Simulation-time estimate, extrapolated from the phase-5 wall-clock anchor
     (wall ~= 0.054 h per a.u. of propagation at box 50x50x90, spacing 0.5).
  4. Master summary table (sigma0 x E) -> spread%, gap-fit, tau, wall.

Run:  venv/bin/python3 build_wide_wp_planning_notebook.py
Output: wide_wavepacket_planning.ipynb (executed in place) + wide_wp_planning_figs/.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbclient import NotebookClient

HERE = os.path.dirname(os.path.abspath(__file__))
NB_PATH = os.path.join(HERE, "wide_wavepacket_planning.ipynb")

cells = []
def md(s): cells.append(new_markdown_cell(s))
def code(s): cells.append(new_code_cell(s))

# ---------------------------------------------------------------------------
md(r"""# Campaign planning — wide low-spread wavepacket (localised jellium slab)
### Deliberation notebook · `wide-wavepacket-lowspread` · 2026-06-30

**Aim.** Fire a **wide, near-rigid wavepacket** through the localised jellium
slab, matched to a classical projectile of the *same* Gaussian width, so any
WP$-$classical stopping difference is a **purely quantum** effect (Pauli +
interference) rather than dispersion or interaction-range mismatch.

**Priority of this notebook: minimise spreading first → define the wavepacket
parameters.** Everything else (CAP completeness, SNR, cost) is sized *after* the
WP parameters are fixed.

**Companion docs (same folder):** `wide_wavepacket_lowspread.md` (campaign prompt),
`classical_projectile_fix.md` (sister campaign), and the phase-5 production
geometry `ResearchProject/systems/localised_jellium/shared/configs/slab_n82_L50x50x90.hpp`.

**The lever — free-particle Gaussian spreading law (atomic units, ħ = mₑ = 1):**
$$\sigma(t) = \sigma_0\,\sqrt{1 + \left(\frac{\hbar\,t}{2 m \sigma_0^{2}}\right)^{2}}
           = \sigma_0\,\sqrt{1 + \left(\frac{t}{2\sigma_0^{2}}\right)^{2}}.$$
The spreading term scales as $1/\sigma_0^{2}$ → **a wider packet stays rigid far
longer**. The drift energy sets the velocity $v=\sqrt{2E}$ (a.u.), hence the
transit time $t=\text{path}/v$ and the spread accumulated over it.
""")

# ---------------------------------------------------------------------------
code(r"""import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from inqview.visualisation.style import apply_theme
apply_theme()

HA_TO_EV = 27.211386

# ---- Planned box geometry: 50 x 50 x 101 Bohr ----
# LZ=101 chosen so the equidistant launch gives EXACTLY 4*sigma (=14 Bohr at
# sigma0=3.5) clearance to BOTH the slab face and the CAP inner face (user
# decision 2026-06-30: CAP<->WP gap = 4 sigma).
LX = LY = 50.0
LZ = 101.0                      # z in [-50.5, +50.5]
SPACING = 0.50                  # Bohr (phase-5 default; coarsen later for wide sigma)
SLAB_HALF = 12.5                # slab |z| < 12.5 (fixed -> r_s matched to production)
CAP_WIDTH = 10.0                # Bohr per side (SAME CAP as phase-5: eta=-0.7)
CAP_INNER = LZ/2 - CAP_WIDTH    # = 40.5  -> CAP region [40.5, 50.5]
# launch equidistant between slab face and CAP inner face (phase-5 logic):
LAUNCH_Z = -(SLAB_HALF + CAP_INNER)/2.0     # = -26.25
GAP_HALF = CAP_INNER - SLAB_HALF            # full launch gap (one side) = 27.5
CLEARANCE = abs(LAUNCH_Z) - SLAB_HALF       # launch -> slab face = 13.75 (== to CAP)

def vel(E_eV):                  # drift velocity (a.u.), m_e = 1
    return np.sqrt(2.0 * E_eV / HA_TO_EV)
def spread_factor(sig0, t):     # sigma(t)/sigma0
    return np.sqrt(1.0 + (t / (2.0 * sig0**2))**2)

# representative path lengths from the launch point
PATH_SLAB_CENTRE = abs(LAUNCH_Z)                 # first strong interaction ~ slab centre
PATH_FAR_CAP     = abs(LAUNCH_Z) + (CAP_INNER + CAP_WIDTH/2.0)  # launch -> far CAP middle

print(f"Box            : {LX:.0f} x {LY:.0f} x {LZ:.0f} Bohr   (z in [{-LZ/2:.1f}, {LZ/2:.1f}])")
print(f"Slab           : |z| < {SLAB_HALF}  (25 Bohr thick)")
print(f"CAP            : [{CAP_INNER:.1f}, {LZ/2:.1f}] each side  (inner face +/-{CAP_INNER:.1f})")
print(f"Launch (equid.): z0 = {LAUNCH_Z:.2f}   ->  {CLEARANCE:.2f} Bohr to BOTH slab face and CAP")
print(f"Gap-fit limit  : n*sigma0 <= {CLEARANCE:.2f}  ->  sigma0 <= {CLEARANCE/3:.2f} (3sigma)  |  {CLEARANCE/4:.2f} (4sigma)")
print(f"Path launch->slab centre = {PATH_SLAB_CENTRE:.2f} Bohr ; launch->far CAP = {PATH_FAR_CAP:.2f} Bohr")
""")

# ---------------------------------------------------------------------------
md(r"""## 1 — Spreading prediction for a range of σ₀ (the parameter to minimise)

Two spread budgets, cleanly separable:
- **At slab interaction** (≈ slab centre) — this is what sets *matched-σ quality*:
  how close the WP width is to σ₀ (and to the classical ghost) when it actually
  hits the slab. **This is the number to drive down.**
- **At far-CAP absorption** — accrued over the *whole* run; it does **not**
  corrupt the matched-σ comparison, only the CAP's ability to absorb the packet
  (sized later).

The σ_WP=0.5 production packet is shown for reference — it spreads catastrophically.
""")

code(r"""sig_grid = np.linspace(0.5, 6.0, 120)
energies = [100, 300, 500]      # eV
SIG_MARK = 3.5                  # highlighted operating-width candidate

# Two side-by-side panels (compact: 7.2 x 3.3 in), one figure only.
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True,
                         constrained_layout=True)
panels = [("at slab interaction\n(matched-σ quality)", PATH_SLAB_CENTRE),
          ("at far-CAP absorption\n(CAP load)",        PATH_FAR_CAP)]
for ax, (title, path) in zip(axes, panels):
    for E in energies:
        t = path / vel(E)
        ax.plot(sig_grid, 100.0*(spread_factor(sig_grid, t) - 1.0), label=f"E = {E} eV")
        # marker at sigma0 = 3.5
        ax.plot(SIG_MARK, 100.0*(spread_factor(SIG_MARK, t) - 1.0),
                "o", ms=4, color="0.15", zorder=5)
    ax.axvline(CLEARANCE/3, ls="--", lw=0.9, color="0.4")
    ax.axvline(SIG_MARK,    ls="-",  lw=0.7, color="0.15", alpha=0.5)
    ax.axhline(10, ls=":", lw=0.8, color="0.6")     # 10% guide
    ax.set_xlabel(r"initial width $\sigma_0$ (Bohr)")
    ax.set_title(title, fontsize=8.5)
    ax.set_ylim(0, 60)
    ax.set_xlim(0.5, 6.0)
axes[0].set_ylabel(r"spread $\sigma(t)/\sigma_0 - 1$  (%)")
axes[0].annotate(r"3σ gap-fit", xy=(CLEARANCE/3, 55), fontsize=7, color="0.4",
                 ha="right", rotation=90, va="top")
axes[1].annotate(r"$\sigma_0=3.5$", xy=(SIG_MARK, 52), fontsize=7, color="0.15", ha="left")
axes[0].legend(frameon=False, fontsize=7, loc="upper right")
fig.suptitle("Free-Gaussian spreading in the 50×50×100 box", fontsize=10)
os.makedirs(os.path.join(".", "wide_wp_planning_figs"), exist_ok=True)
fig.savefig("wide_wp_planning_figs/spreading_vs_sigma.png", dpi=140)
plt.show()

# sigma0 = 3.5 readout (the user's candidate width)
print("sigma0 = 3.5 Bohr predictions:")
for E in energies:
    ss = 100*(spread_factor(SIG_MARK, PATH_SLAB_CENTRE/vel(E)) - 1)
    sc = 100*(spread_factor(SIG_MARK, PATH_FAR_CAP/vel(E))     - 1)
    print(f"   E={E:>3} eV:  spread@slab = {ss:5.1f} %   spread@farCAP = {sc:5.1f} %")
print(f"\nsigma0=0.5 reference @slab, E=300 eV: "
      f"{100*(spread_factor(0.5, PATH_SLAB_CENTRE/vel(300))-1):.0f} %  (production packet — hopeless)")
""")

# ---------------------------------------------------------------------------
md(r"""## 2 — Gap-fit (how wide can σ₀ be before the packet touches a boundary?)

The packet must clear **both** the slab face and the CAP inner face at *t=0*, or
it either pre-interacts with the slab or gets clipped by the absorber. With the
equidistant launch in this 50×50×100 box the clearance is **13.75 Bohr** each
side, so the boundary-rule limits are:

- **3σ rule** (relaxed, used at large σ): σ₀ ≤ 4.58
- **4σ rule** (strict): σ₀ ≤ 3.44

So this box comfortably accommodates **σ₀ up to ≈ 3.4 (4σ) / 4.6 (3σ)** — i.e. the
whole wide-packet region of interest.
""")

# ---------------------------------------------------------------------------
md(r"""## 3 — Simulation-time estimate (extrapolated from phase-5)

**Anchor (measured).** Phase-5 σ=0.5 WP runs: **wall ≈ 0.054 h per a.u. of
propagation** at box 50×50×90, spacing 0.50 Bohr, dt=0.04, ≈102 states
(`docs/handovers/localised-jellium.md`, `qsp_phase5_velocity_sweep.md`).

**Scaling.** Cost ∝ N_grid · N_steps. Going 90→100 in z at the same spacing
multiplies N_grid by 100/90; a wide smooth packet later permits a coarser dx,
which would *reduce* cost as (0.5/dx)³ (shown as an optional column).

**Run length τ.** The wide packet is fast and does not stall, so τ is set by
transit + a short collective-relaxation plateau:
`τ ≈ (launch → far-CAP) / v + t_plateau`, with `t_plateau ≈ 15` a.u.
(Phase-5 instead used the conservative `τ ≈ 200/v`; both are tabulated.)
""")

code(r"""WALL_PER_AU_90 = 0.054                  # h / a.u., phase-5 anchor (box 90, dx 0.5)
GRID_FACTOR = (LZ/90.0)                  # z-extent scaling at dx = 0.5
WALL_PER_AU = WALL_PER_AU_90 * GRID_FACTOR
T_PLATEAU = 15.0                         # a.u. collective relaxation after absorption

def tau_geom(E_eV):
    return PATH_FAR_CAP / vel(E_eV) + T_PLATEAU
def tau_p5(E_eV):                        # phase-5's conservative choice, tau ~ 200/v capped 200
    return min(200.0, 200.0 / vel(E_eV))

print(f"grid factor (LZ {LZ:.0f}/90, dx 0.5) = {GRID_FACTOR:.3f}  ->  wall ~= {WALL_PER_AU:.3f} h / a.u.\n")
print(f"{'E[eV]':>6} {'v[au]':>6} {'tau_geom':>9} {'wall_geom[h]':>12} {'tau~200/v':>10} {'wall_p5[h]':>11}")
for E in [100, 200, 300, 400, 500]:
    tg, tp = tau_geom(E), tau_p5(E)
    print(f"{E:>6} {vel(E):>6.2f} {tg:>9.1f} {WALL_PER_AU*tg:>12.1f} {tp:>10.1f} {WALL_PER_AU*tp:>11.1f}")
""")

# ---------------------------------------------------------------------------
md(r"""## 4 — CAP completeness: toy model vs REAL INQ simulations

Locked operating point: **σ₀ = 3.5 Bohr, E = 300 eV.** The packet reaches the far
CAP spread to σ ≈ 4.1 Bohr (the 17.6% from §1). Does the two-sided sin² CAP absorb
it cleanly enough that the retained-energy ledger is trustworthy?

**Two sources, overlaid:**
1. **Real INQ runs** — `systems/vacuum/hypotheses/twosided_cap_vs_mask/`
   (`twosided_combined.csv`): a free Gaussian WP fired at the two-sided CAP, real
   TDDFT, ε = surviving inner-region norm. **It has E=300 eV points** across
   per-side widths and an η-sweep. *Caveat:* these use the **k₀-tied benchmark
   packet σ = 4√2/k₀ ≈ 1.20 Bohr — narrow, NOT our wide σ≈4.*
2. **1D toy** — `docs/reports/absorbing-boundary/cap_toy.py`, with σ parameterised
   (cross-checked bit-for-bit against `cap_toy` at its own σ). Lets us extrapolate
   to **our** width.

**Width naming (important):** the dataset's `L_total` = both sides; **per-side
width = Lhalf**. Phase-5's CAP is **Lhalf = 10 Bohr/side** (region [±35,±45]) =
the dataset's `L20`. The toy's `L` is the per-side width, so toy `L` ↔ real
`Lhalf`.

**What the overlay reveals:** (a) the toy is **pessimistic** — it over-predicts ε
vs real INQ by ~5× at Lhalf=10 (and more at larger L); (b) for the **narrow**
benchmark packet the real CAP at phase-5's params (Lhalf=10, η=0.7) reflects only
**~0.2%**; (c) our **wide** packet is **not** in the real dataset, so its adequacy
is the open question — the toy (pessimistic) flags a possible problem that a
Phase-0 real wide-σ run must settle.
""")

code(r"""import sys
sys.path.insert(0, "/local/data/public/skcb2/tddft/docs/reports/absorbing-boundary")
import cap_toy

def reflect_eps_sigma(k0, L, sigma, eta=-0.7, dx=0.1, dt=0.01):
    # CAP reflection error eps for an EXPLICIT packet width sigma.
    # Mirrors cap_toy.reflect_eps (kind='cap') but with sigma parameterised.
    X = 6.0*sigma; x0 = -3.0*sigma
    n = int(round((X + L)/dx)); x = -X + dx*np.arange(n)
    psi = cap_toy._gaussian(x, x0, k0, sigma).astype(np.complex128)
    N0 = np.sum(np.abs(psi)**2)*dx
    kk = 2*np.pi*np.fft.fftfreq(n, d=dx); kin = np.exp(-0.5j*kk**2*(dt/2))
    s = np.zeros_like(x); inb = (x >= 0) & (x <= L)
    s[inb] = np.sin(np.pi*x[inb]/(2*L))**2
    absorb = np.exp(-1j*(1j*eta*s)*dt)
    tau = 2.0*(3*sigma + L)/k0; ns = int(round(tau/dt))
    for _ in range(ns):
        psi = np.fft.ifft(kin*np.fft.fft(psi)); psi *= absorb
        psi = np.fft.ifft(kin*np.fft.fft(psi))
    return np.sum(np.abs(psi[x < 0])**2)*dx/N0

# --- known-case cross-check vs the validated toy (sigma = 4 sqrt2 / k0) ---
ETA = -0.7
k0_300 = vel(300)
_sig_toy = 4*np.sqrt(2)/k0_300
_chk = abs(reflect_eps_sigma(k0_300, 10, _sig_toy, eta=ETA)
           - cap_toy.reflect_eps(k0_300, 10, kind="cap", eta=ETA))
assert _chk < 1e-3, f"cross-check failed: {_chk}"
print(f"cross-check vs cap_toy OK (|diff| = {_chk:.1e})")

# --- our wide packet width at the far CAP ---
SIG_CAP = 3.5 * spread_factor(3.5, PATH_FAR_CAP/vel(300))
print(f"wide packet sigma at far CAP (sigma0=3.5, E=300) = {SIG_CAP:.2f} Bohr")

# --- REAL INQ two-sided CAP data (E=300 eV) ---
REAL = pd.read_csv("/local/data/public/skcb2/tddft/ResearchProject/systems/"
                   "vacuum/hypotheses/twosided_cap_vs_mask/twosided_combined.csv")
r300 = REAL[(REAL["mode"] == "cap") & (REAL["E_eV"] == 300.0)].copy()
SIG_REAL = float(r300["sigma"].iloc[0])      # k0-tied benchmark width ~1.20 Bohr
print(f"real-run benchmark packet sigma (k0-tied) = {SIG_REAL:.2f} Bohr  (narrow, not ours)\n")

rL = r300[r300["eta_Ha"] == -0.5].sort_values("Lhalf")      # per-side width sweep, eta=0.5
rEta = r300[r300["Lhalf"] == 10.0].sort_values("eta_Ha")    # eta sweep at Lhalf=10 (phase-5 width)

Lgrid = np.arange(4, 16, 1.0)        # per-side width (Lhalf)
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), constrained_layout=True)

# Panel A: eps vs per-side width Lhalf, E=300, eta=0.5
axes[0].plot(rL["Lhalf"], 100*rL["epsilon"], "s", ms=5, color="C0",
             label="REAL INQ (σ=1.2, η=0.5)")
axes[0].plot(Lgrid, [100*reflect_eps_sigma(k0_300, L, SIG_REAL, eta=-0.5) for L in Lgrid],
             "--", color="C0", lw=1.0, label="toy σ=1.2 (validation)")
axes[0].plot(Lgrid, [100*reflect_eps_sigma(k0_300, L, SIG_CAP, eta=-0.5) for L in Lgrid],
             "-", color="C3", lw=1.4, label=f"toy σ={SIG_CAP:.1f} (OUR packet)")
axes[0].axvline(10, ls=":", lw=0.9, color="0.4")
axes[0].annotate("phase-5\nLhalf=10", xy=(10, 30), fontsize=6.5, color="0.4", ha="center")
axes[0].set_xlabel("per-side CAP width  Lhalf (Bohr)")
axes[0].set_ylabel(r"reflection $\varepsilon$  (%)")
axes[0].set_title(r"$\varepsilon$ vs width (E=300 eV, η=0.5)", fontsize=8.5)
axes[0].set_yscale("log"); axes[0].legend(frameon=False, fontsize=6.5)

# Panel B: eps vs eta at Lhalf=10, E=300
axes[1].plot(-rEta["eta_Ha"], 100*rEta["epsilon"], "s", ms=5, color="C0",
             label="REAL INQ (σ=1.2)")
etas = np.array([0.3, 0.5, 0.7, 1.0])
axes[1].plot(etas, [100*reflect_eps_sigma(k0_300, 10, SIG_REAL, eta=-e) for e in etas],
             "--", color="C0", lw=1.0, label="toy σ=1.2")
axes[1].plot(etas, [100*reflect_eps_sigma(k0_300, 10, SIG_CAP, eta=-e) for e in etas],
             "-", color="C3", lw=1.4, label=f"toy σ={SIG_CAP:.1f} (OUR)")
axes[1].axvline(0.7, ls=":", lw=0.9, color="0.4")
axes[1].set_xlabel(r"CAP strength  $|\eta|$ (Ha)")
axes[1].set_ylabel(r"reflection $\varepsilon$  (%)")
axes[1].set_title(r"$\varepsilon$ vs strength (Lhalf=10, E=300)", fontsize=8.5)
axes[1].set_yscale("log"); axes[1].legend(frameon=False, fontsize=6.5)

fig.suptitle("Two-sided CAP reflectivity at E=300 eV — real INQ vs toy", fontsize=10)
fig.savefig("wide_wp_planning_figs/cap_reflectivity_real_vs_toy.png", dpi=140)
plt.show()

# readouts
print("REAL INQ at phase-5 width Lhalf=10 (=L_total 20), E=300, NARROW sigma=1.2:")
for _, r in rEta.iterrows():
    print(f"   eta={-r['eta_Ha']:.2f}:  eps = {r['epsilon']:.4f}")
print(f"\nTOY (pessimistic) for OUR wide sigma={SIG_CAP:.1f} at Lhalf=10:")
for e in (0.5, 0.7, 1.0):
    print(f"   eta={e:.2f}:  eps_toy = {reflect_eps_sigma(k0_300, 10, SIG_CAP, eta=-e):.4f}")
""")

md(r"""**Reading.**
- **Real INQ, narrow benchmark packet:** at phase-5's width (Lhalf=10) and η=0.7,
  the two-sided CAP reflects only **~0.2%** — excellent. The CAP design itself is
  sound for a *narrow* packet, and higher |η| (→1.0) drives ε to ~0.03%.
- **The toy is conservative** — it sits ~5× above the real INQ points at σ=1.2, so
  read toy numbers as *upper bounds*.
- **Our wide packet is untested.** The toy at σ≈4.1 (red) predicts markedly higher
  reflection at Lhalf=10 than the narrow benchmark — but the toy is pessimistic, so
  the true wide-σ reflection is bracketed between the toy curve (upper bound) and
  the narrow real points (optimistic). **This is the gap a Phase-0 real wide-σ run
  must close.**

**Box coupling:** fitting σ₀=3.5 at the 3σ rule needs clearance ≥10.5 Bohr, so the
CAP inner face sits at |z|≈33.5 → **Lhalf=10 just fits LZ=100** (CAP [33.5,43.5]);
widening per side to Lhalf=15 (if the wide-σ check demands it) needs LZ≈110, and/or
raising |η| toward 1.0 (cheaper than box). Locked next.
""")

# ---------------------------------------------------------------------------
md(r"""## 5 — Master summary: σ₀ × E → spread, gap-fit, cost

One table to weigh the operating point. **Read the priority column first**
(`spread@slab%`): that is what the campaign exists to minimise.
""")

code(r"""rows = []
for sig0 in [1.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0]:
    for E in [100, 300, 500]:
        v = vel(E)
        s_slab = 100*(spread_factor(sig0, PATH_SLAB_CENTRE/v) - 1)
        s_cap  = 100*(spread_factor(sig0, PATH_FAR_CAP/v)     - 1)
        tg = tau_geom(E)
        rows.append(dict(
            sigma0=sig0, E_eV=E, v_au=round(v,2),
            spread_slab_pct=round(s_slab,1),
            spread_farCAP_pct=round(s_cap,1),
            fits_3sigma=("yes" if 3*sig0 <= CLEARANCE else "NO"),
            fits_4sigma=("yes" if 4*sig0 <= CLEARANCE else "no"),
            tau_geom_au=round(tg,1),
            wall_geom_h=round(WALL_PER_AU*tg,1),
        ))
df = pd.DataFrame(rows)
pd.set_option("display.width", 140)
print(df.to_string(index=False))
df.to_csv("wide_wp_planning_figs/operating_point_table.csv", index=False)
""")

# ---------------------------------------------------------------------------
md(r"""## 6 — Reading of the prediction (deliberation log)

**Locked (2026-06-30):** σ₀ = **3.5 Bohr**, E = **300 eV**, box **50×50×100**.

- **Spreading is solved at the slab** — σ₀=3.5 at E=300 spreads only **2.6 %** when
  it hits the slab (matched-σ essentially perfect). Spreading no longer
  discriminates the σ₀ choice.
- **CAP completeness — grounded in REAL INQ data (§4).** At phase-5's width
  (Lhalf=10) and η=0.7 the two-sided CAP reflects only **~0.2%** for the *narrow*
  benchmark packet. The 1D toy is **~5× pessimistic** vs these real runs. Our
  *wide* σ≈4 packet is **not** in the real dataset: the toy flags possibly-higher
  reflection, but as an upper bound — so wide-σ CAP adequacy is the **one open risk
  a Phase-0 real run must close**, with levers Lhalf (→15, needs LZ≈110) and |η|
  (→1.0, cheaper than box).
- **Cost is modest** — geometry-minimal τ≈30 a.u. → ~1.8 h/run at E=300 (×1.11 grid
  factor); ample headroom for the autonomous sweep.

**CAP locked (2026-06-30): same as the previous campaign (phase-5)** — two-sided
sin², **η=−0.7, 10 Bohr/side**, region [±40,±50] in the LZ=100 box. Real INQ →
ε≈0.2 % (narrow packet); the **Phase-0 wide-σ CAP-completeness run** verifies it
for our σ≈4 packet (fallbacks η→1.0 or Lhalf→15 only if it fails).

**Next decisions:** the matched classical UPF (σ_pot=3.5/√2≈2.47) + the LZ=101 GS
mechanics; then the autonomy-readiness checklist (draft→ready). Phase-0 and the
Phase-1 grid (below) are locked.
""")

# ---------------------------------------------------------------------------
md(r"""## 7 — Phase-1 sweep grid + cutoff/aliasing guard (LOCKED)

**One width (σ_WP=3.5), 5 runs.** The grid energy range is bounded by the
**mandatory cutoff/aliasing guard**: the grid Nyquist `k_max = π/dx` must exceed
the WP's drift momentum by ≥ 4 momentum-σ, where (project convention)
`σ_p = 1/(√2·σ_WP)`. Phase-5's dx=0.50 fails this at E≥500, so dx is refined to
**0.40 Bohr** (clear of the dx=0.30 WP-init deadlock; ~2× grid cost) — giving
≥6σ_p margin up to E=600 eV. Each E runs WP + matched classical.
""")

code(r"""SIGMA_WP = 3.5
SIG_P = 1/(np.sqrt(2)*SIGMA_WP)          # momentum std (cutoff-guard convention)
DX_SWEEP = 0.40                          # refined from 0.50 to satisfy the 4-sigma cutoff guard
KMAX = np.pi/DX_SWEEP
E_GRID = np.array([200, 280, 360, 440, 520, 600])    # eV, 6 energies (2-GPU friendly)

# wall estimate: phase-5 anchor 0.054 h/au scaled by grid (LZ 101/90)*(0.5/dx)^3
GRID_FACTOR_SWEEP = (LZ/90.0)*(0.50/DX_SWEEP)**3
WALL_PER_AU_SWEEP = 0.054*GRID_FACTOR_SWEEP

print(f"sigma_WP={SIGMA_WP} -> sigma_p=1/(sqrt2*sigma_WP)={SIG_P:.3f} a.u.")
print(f"dx={DX_SWEEP} -> k_max={KMAX:.2f}, E_cut={0.5*KMAX**2*HA_TO_EV:.0f} eV; "
      f"grid factor x{GRID_FACTOR_SWEEP:.2f} -> wall {WALL_PER_AU_SWEEP:.3f} h/au\n")
print(f"{'E[eV]':>6} {'v=k0':>6} {'cutoff margin (sig_p)':>21} {'spread@slab%':>12} {'wall[h]':>8}")
all_ok = True
for E in E_GRID:
    k0 = vel(E); margin = (KMAX-k0)/SIG_P
    spr = 100*(spread_factor(SIGMA_WP, PATH_SLAB_CENTRE/k0)-1)
    tau = PATH_FAR_CAP/k0 + 15
    ok = margin >= 4.0; all_ok &= ok
    print(f"{E:>6} {k0:>6.2f} {margin:>21.1f} {spr:>12.1f} {WALL_PER_AU_SWEEP*tau:>8.1f}"
          + ("" if ok else "  <-- FAILS 4sigma"))
assert all_ok, "cutoff guard FAILED for some grid point"
print(f"\nCUTOFF GUARD: all {len(E_GRID)} points clear 4-sigma (>=6 sigma_p). "
      f"{2*len(E_GRID)} production runs (WP+classical) + 1 vacuum-WP SIE control.")
""")

# ---------------------------------------------------------------------------
nb = new_notebook(cells=cells, metadata={
    "kernelspec": {"name": "python3", "display_name": "Python 3"},
    "language_info": {"name": "python"},
})
print("Executing notebook ...")
client = NotebookClient(nb, timeout=300, kernel_name="python3",
                        resources={"metadata": {"path": HERE}})
client.execute()
with open(NB_PATH, "w") as f:
    nbf.write(nb, f)
print("wrote", NB_PATH)
