#!/usr/bin/env python3
"""Build (and execute) momentum_and_effective_width.ipynb.

Cross-run analysis of the localised-jellium (r_s=4.18) high-density WP S(v) sweep:
  1. WP momentum distribution start-vs-end (Gaussian reconstructed from saved
     moments) + exact moment time-series, for sigma_WP = 0.5 / 2.0 / 3.0.
  2. Effective time-averaged width <sigma_r>, showing the launch-sigma ordering
     inverts under dispersion.
  3. A new S(v) figure relabelled by <sigma_r>, and an S-vs-effective-width
     cross-plot reconciling the WP trend with the classical one.

Pure post-processing of exported sweep_data/ CSVs. No new simulation.
Run:  venv/bin/python3 build_momentum_width_notebook.py
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_IPYNB = os.path.join(HERE, "momentum_and_effective_width.ipynb")

nb = new_notebook()
cells = []

# ---------------------------------------------------------------- intro
cells.append(new_markdown_cell(r"""# WP momentum distributions and the effective (time-averaged) width

**System:** localised-jellium slab, $r_s = 4.18$ ("density A"), E-absorbed stopping
$S_B = (E_\text{tot}(t_f) - E_\text{GS})/L_\text{slab}$ under a two-sided CAP.

## The puzzle

Comparing classical and wavepacket (WP) stopping across widths:

- **Wavepacket:** at fixed velocity, $S_B$ *rises* with the launch width
  $\sigma_\text{WP}$  ($0.5 < 2.0 < 3.0$).
- **Classical:** at fixed velocity, $S_B$ *falls* with the potential width
  $\sigma_\text{pot}$.

The two trends have **opposite sign in $\sigma$**. This notebook tests one
explanation: a wavepacket does **not** keep its launch width — it *disperses* —
so labelling a WP run by $\sigma_\text{WP}$ mis-states the width the slab actually
sees. A narrow packet ($\sigma_\text{WP}=0.5$, spreading time $\tau=\sigma_0^2=0.25$
a.u.) spreads almost immediately; a wide one ($\sigma_\text{WP}=3$, $\tau=9$ a.u.)
stays compact. So the *time-averaged* width $\langle\sigma_r\rangle$ may **reorder**
the WP cases and put them on the same footing as the classical trend.

We (1) look at the momentum spectrum start-vs-end, (2) build $\langle\sigma_r\rangle$
per run, and (3) redraw $S(v)$ against this effective width.

> **Convention.** $\sigma$ always denotes the wavepacket width $\sigma_\text{WP}$;
> the classical Gaussian *potential* width is the derived $\sigma_\text{pot}=
> \sigma_\text{WP}/\sqrt2$ (methods only)."""))

# ---------------------------------------------------------------- setup
cells.append(new_code_cell(r"""import os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from inqview.visualisation import style

style.apply_theme()

HYP = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses"
WP  = f"{HYP}/wp_highdensity_sv"
SW  = f"{WP}/sweep_data"
CL  = f"{HYP}/classical_highdensity_sv/dyn_direct"

# sigma_WP -> (sweep_data prefix, colour, S-summary file)
SIG = {
    0.5: ("",      "#1f5fb4", f"{WP}/wp_S_summary.csv"),
    2.0: ("s2p0_", "#d1600a", f"{WP}/wp_S_summary_s2p0.csv"),
    3.0: ("s3p0_", "#2ca02c", f"{WP}/wp_S_summary_s3p0.csv"),
}
VTAG = {2.0: "v2p0", 2.5: "v2p5", 3.0: "v3p0", 3.5: "v3p5"}
VS   = [2.0, 2.5, 3.0, 3.5]


def load_concat(run_dir, base):
    "Concatenate base.csv + segment-suffixed base.fromNNNN.csv in step order."
    files = sorted(glob.glob(f"{run_dir}/{base}.csv")
                   + glob.glob(f"{run_dir}/{base}.from*.csv"))
    df = pd.concat([pd.read_csv(f, comment="#") for f in files])
    return df.drop_duplicates("step").sort_values("step").reset_index(drop=True)


def run_dir(sig, v):
    return f"{SW}/{SIG[sig][0]}{VTAG[v]}"


# slab is centred at z=0 with slab_half_width = 12.5 Bohr; the WP launches at
# z = -24 moving +z, so the NEAR (closer) face it reaches first is at z = -12.5.
SLAB_NEAR = -12.5


def widths(sig, v):
    "Per-step sigma_z(t), sigma_r(t), z_mean(t) and time (a.u.) for one run."
    rs = load_concat(run_dir(sig, v), "wp_real_space_stats")
    sig_z = np.sqrt(rs["sigma_z2"].to_numpy())
    sig_r = np.sqrt((rs["sigma_x2"] + rs["sigma_y2"] + rs["sigma_z2"]).to_numpy())
    return (rs["time_au"].to_numpy(), sig_z, sig_r, rs["z_mean"].to_numpy())


def entry_width(sig, v):
    "Effective width at slab entry: sigma_z, sigma_r and time when the centroid"
    " z_mean first reaches the near face (z=-12.5). Evaluated BEFORE the packet"
    " overlaps the slab, so it is free of interference-window contamination."
    t, sz, sr, zc = widths(sig, v)
    idx = np.argmax(zc >= SLAB_NEAR)          # first step past the near face
    if not (zc >= SLAB_NEAR).any():
        return np.nan, np.nan, np.nan
    return sz[idx], sr[idx], t[idx]


def moments(sig, v):
    "Per-step momentum moments for one run."
    m = load_concat(run_dir(sig, v), "wp_momentum_stats")
    return (m["time_au"].to_numpy(), m["pz_mean"].to_numpy(),
            np.sqrt(m["sigma_pz2"].to_numpy()), m["e_kin_ha"].to_numpy())


print("runs found:",
      sum(os.path.isdir(run_dir(s, v)) for s in SIG for v in VS), "/ 12")"""))

# ---------------------------------------------------------------- Section A
cells.append(new_markdown_cell(r"""## 1. Momentum distribution: start vs end

The sweep runs saved momentum **moments** per step (`pz_mean`, `sigma_pz2`,
`e_kin_ha`), not the full spectrum $P(p_z)$. We therefore reconstruct the
**Gaussian** $P(p_z)=\mathcal N(\langle p_z\rangle,\ \sigma_{p_z}^2)$ at the first
and last step — exact for the launch packet (a Gaussian), and a two-moment summary
thereafter.

> **Caveat (inference).** After strong scattering the true $P(p_z)$ can be
> non-Gaussian (a bimodal transmitted/reflected split); the Gaussian below
> captures only its mean and spread. The *exact* moment time-series underneath
> shows the real evolution.

Representative velocity $v=3.0$ a.u. (momentum aliasing $\le 1.2\%$ here for all
three widths; it reaches $\sim5\%$ only at $v=3.5$, $\sigma_\text{WP}=0.5$)."""))

cells.append(new_code_cell(r"""VREP = 3.0
fig, ax = style.figure_one_col()
pz = np.linspace(-1.5, 4.5, 600)

def gauss(x, mu, s):
    return np.exp(-0.5 * ((x - mu) / s) ** 2) / (s * np.sqrt(2 * np.pi))

for sig in SIG:
    t, mu, spz, ek = moments(sig, VREP)
    c = SIG[sig][1]
    ax.plot(pz, gauss(pz, mu[0],  spz[0]),  ls="--", color=c, lw=1.3,
            label=f"$\\sigma_{{WP}}$={sig}  start")
    ax.plot(pz, gauss(pz, mu[-1], spz[-1]), ls="-",  color=c, lw=1.6,
            label=f"$\\sigma_{{WP}}$={sig}  end")

ax.axvline(0, color="0.6", lw=0.7, zorder=0)
ax.set_xlabel(style.axis_label("momentum", symbol="$p_z$"))
ax.set_ylabel("$P(p_z)$  (reconstructed)")
ax.set_title(f"WP momentum spectrum, $v={VREP}$ a.u. (dashed = start, solid = end)",
             fontsize=9)
ax.legend(fontsize=6.5, frameon=False, ncol=1, loc="upper left")
fig.savefig(f"{WP}/momentum_start_end.png", dpi=300, bbox_inches="tight")
print("wrote momentum_start_end.png")
plt.show()"""))

cells.append(new_markdown_cell(r"""**Reading it.** Every packet starts centred at $p_z \approx v = 3$. By the end the
mean momentum has collapsed toward zero (deceleration) and the spread has
*narrowed* — the fast components are scattered and removed by the CAP. Momentum is
nearly invariant under *free* dispersion, so this collapse is a genuine
scattering / absorption signature, the exact complement of the real-space spreading
in Section 2. The narrow launch packet ($\sigma_\text{WP}=0.5$) starts with the
*broadest* momentum spread (position–momentum uncertainty), consistent with it
dispersing fastest in real space."""))

cells.append(new_code_cell(r"""# Exact moment time-series (not a reconstruction): p_z, sigma_pz, e_kin.
fig, axs = style.figure_two_col(height_in=2.6)
fig.clf()
axs = fig.subplots(1, 3)
for sig in SIG:
    t, mu, spz, ek = moments(sig, VREP)
    c = SIG[sig][1]
    axs[0].plot(t, mu,  color=c, lw=1.4, label=f"$\\sigma_{{WP}}$={sig}")
    axs[1].plot(t, spz, color=c, lw=1.4)
    axs[2].plot(t, ek,  color=c, lw=1.4)
for a in axs:
    a.set_xlabel("time (a.u.)")
axs[0].set_ylabel("$\\langle p_z\\rangle$ (a.u.)")
axs[1].set_ylabel("$\\sigma_{p_z}$ (a.u.)")
axs[2].set_ylabel("WP $E_\\mathrm{kin}$ (Ha)")
axs[0].legend(fontsize=6.5, frameon=False)
fig.suptitle(f"WP momentum moments over time, $v={VREP}$ a.u.", fontsize=9)
fig.tight_layout()
fig.savefig(f"{WP}/momentum_moments_time.png", dpi=300, bbox_inches="tight")
print("wrote momentum_moments_time.png")
plt.show()"""))

# ---------------------------------------------------------------- Section B
cells.append(new_markdown_cell(r"""## 2. Effective (time-averaged) width

Real-space spread $\sigma_r(t)=\sqrt{\sigma_x^2+\sigma_y^2+\sigma_z^2}$ and the
propagation-axis $\sigma_z(t)$. The time average over the run,
$\langle\sigma_r\rangle$, is the "width the slab sees" — the quantity the user
proposed as more representative than the launch label."""))

cells.append(new_code_cell(r"""fig, axs = style.figure_two_col(height_in=2.8)
fig.clf(); axs = fig.subplots(1, 2)
for sig in SIG:
    t, sz, sr, zc = widths(sig, VREP)
    c = SIG[sig][1]
    axs[0].plot(t, sz, color=c, lw=1.5, label=f"$\\sigma_{{WP}}$={sig}")
    axs[1].plot(t, sr, color=c, lw=1.5)
    # mark slab-entry (centroid at near face): the "effective width at entry"
    ez, er, et = entry_width(sig, VREP)
    axs[0].plot(et, ez, "o", color=c, ms=7, mec="k", mew=0.6, zorder=5)
    axs[1].plot(et, er, "o", color=c, ms=7, mec="k", mew=0.6, zorder=5)
for a, lab in zip(axs, ["$\\sigma_z(t)$", "$\\sigma_r(t)$"]):
    a.set_xlabel("time (a.u.)"); a.set_ylabel(f"{lab}  (Bohr)")
axs[0].legend(fontsize=7, frameon=False, title="launch width")
fig.suptitle(f"Real-space spreading, $v={VREP}$ a.u. (narrow packet disperses "
             "fastest; ● = slab entry)", fontsize=9)
fig.tight_layout()
fig.savefig(f"{WP}/width_time.png", dpi=300, bbox_inches="tight")
print("wrote width_time.png")
plt.show()"""))

cells.append(new_code_cell(r"""# Per-run width table: launch, sigma_r(0), entry widths, <sigma_z>, <sigma_r>.
rows = []
for sig in SIG:
    for v in VS:
        if not os.path.isdir(run_dir(sig, v)):
            continue
        t, sz, sr, zc = widths(sig, v)
        ez, er, et = entry_width(sig, v)
        rows.append(dict(sigma_WP=sig, v=v, sigma_r0=sr[0],
                         t_entry=et, sigma_z_entry=ez, sigma_r_entry=er,
                         mean_sigma_z=sz.mean(), mean_sigma_r=sr.mean()))
wtab = pd.DataFrame(rows)
# per-sigma_WP means across velocity (the headline effective widths)
summ = (wtab.groupby("sigma_WP")[["sigma_z_entry", "sigma_r_entry",
                                  "mean_sigma_z", "mean_sigma_r"]]
        .mean().round(2))
print("Per-run widths (Bohr):")
print(wtab.round(3).to_string(index=False))
print("\nEffective width per launch sigma_WP (mean over v):")
print(summ.to_string())
EFFR    = summ["mean_sigma_r"].to_dict()    # sigma_WP -> <sigma_r> (full-run avg)
ENTRY_R = summ["sigma_r_entry"].to_dict()   # sigma_WP -> sigma_r at slab entry
wtab"""))

cells.append(new_markdown_cell(r"""**The inversion.** Averaged over the run, the ordering of the three cases
*flips* relative to the launch label:

| launch $\sigma_\text{WP}$ | $\langle\sigma_z\rangle$ | $\langle\sigma_r\rangle$ |
|---|---|---|
| 0.5 | **largest** | **largest** |
| 2.0 | middle | middle |
| 3.0 | **smallest** | **smallest** |

The $\sigma_\text{WP}=0.5$ packet, launched narrowest, is the **widest** on average
because it disperses fastest; the $\sigma_\text{WP}=3$ packet stays the most compact.
Both metrics ($\sigma_z$ along the motion — least affected by the periodic transverse
box — and the full $\sigma_r$) give the same ordering, so the conclusion is robust.

### Effective width *at entry* (interference-free)

The full-run average blends in the post-collision phase, where slab interference
inflates $\sigma$. A cleaner "width the slab first sees" is $\sigma_r$ evaluated the
moment the **centroid crosses the near slab face** ($z=-12.5$ Bohr) — *before* the
packet overlaps the slab, so it carries no interference contamination. These are the
● markers on the spreading plot above and the `sigma_{z,r}_entry` columns in the
table.

What the entry width shows (mean over $v$, Bohr):

| $\sigma_\text{WP}$ | entry $\sigma_r$ | full-run $\langle\sigma_r\rangle$ |
|---|---|---|
| 0.5 | **10.8** | 21.8 |
| 2.0 | 3.8 | 18.4 |
| 3.0 | 4.2 | 17.0 |

Two things stand out, and they are **not** the same story:

- **$\sigma_\text{WP}=0.5$ is unambiguously the widest at entry** — it has already
  dispersed to $\sim$11 Bohr before touching the slab. This is robust and
  interference-free, and by itself explains why the "narrowest-launched" packet
  gives the *least* stopping.
- **$\sigma_\text{WP}=2.0$ and $3.0$ arrive nearly tied** ($\sim$3.8 vs 4.2 Bohr):
  at entry the width barely distinguishes them, yet their $S_B$ differs by $\sim$40%.
  Only the *full-run* $\langle\sigma_r\rangle$ separates them (2.0 disperses faster
  *after* entry, ending up wider on average). So the fine 2.0-vs-3.0 ordering leans
  on later-time width — precisely the interference-contaminated window this entry
  check isolates.

**Bottom line of the sanity check:** the dominant anomaly (the narrow packet) is a
genuine on-the-way-in dispersion effect; the 2.0-vs-3.0 refinement is weaker and
depends on the late-time average, so it should be read with caution."""))

# ---------------------------------------------------------------- Section C
cells.append(new_markdown_cell(r"""## 3. $S(v)$ redrawn against the effective width

If width drives stopping the *same way* for both models, then relabelling the WP
curves by $\langle\sigma_r\rangle$ (instead of $\sigma_\text{WP}$) should line the
WP trend up with the classical one: **narrower effective width $\to$ more
stopping.**"""))

cells.append(new_code_cell(r"""# (a) S(v): same curves, WP legend relabelled by <sigma_r>.
def load_S(path, scol, sigfilter=None):
    df = pd.read_csv(path)
    if sigfilter is not None:
        df = df[df["sigma_WP"] == sigfilter]
    return df.sort_values("v")

fig, ax = style.figure_one_col()
# classical (hollow, dashed) for reference, labelled by sigma_pot
CLS = [(0.5, f"{CL}/S_of_v_cap.csv", None),
       (2.0, f"{CL}/S_of_v_cap_sigma.csv", 2.0),
       (3.0, f"{CL}/S_of_v_cap_sigma.csv", 3.0)]
for sig, path, filt in CLS:
    d = load_S(path, "S_B_Eabs", filt)
    c = SIG[sig][1]
    ax.plot(d["v"], d["S_B_Eabs"], ls="--", lw=1.0, color=c, marker="o",
            ms=6, mfc="none", mec=c, mew=1.2,
            label=f"cl. $\\sigma_\\mathrm{{pot}}$={sig/np.sqrt(2):.2f}")
# WP (filled, solid), labelled by <sigma_r>
for sig in SIG:
    d = pd.read_csv(SIG[sig][2]).sort_values("v")
    c = SIG[sig][1]
    ax.plot(d["v"], d["S_deposit_corrected"], ls="-", lw=1.0, color=c, marker="o",
            ms=6, mfc=c, mec=c, mew=1.2,
            label=f"WP $\\langle\\sigma_r\\rangle$$\\approx${EFFR[sig]:.0f}")
ax.set_xlabel("projectile velocity $v$ (a.u.)")
ax.set_ylabel(style.axis_label("stopping_power", symbol="$S_B$"))
ax.set_ylim(bottom=0)
ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="upper right")
ax.set_title("E-absorbed stopping, WP relabelled by effective width", fontsize=9)
fig.savefig(f"{WP}/S_of_v_effective_width.png", dpi=300, bbox_inches="tight")
print("wrote S_of_v_effective_width.png")
plt.show()"""))

cells.append(new_code_cell(r"""# (b) The reconciliation: S vs effective width (log x), coloured by velocity.
#   classical at fixed sigma_pot ; WP at its per-run <sigma_r>.
vcol = {2.0: "#4477aa", 2.5: "#66ccee", 3.0: "#ee6677", 3.5: "#aa3377"}
fig, ax = style.figure_one_col()

# WP: per (sigma_WP, v) point at x = <sigma_r>(run)
wp_S = {sig: pd.read_csv(SIG[sig][2]).set_index("v")["S_deposit_corrected"]
        for sig in SIG}
for v in VS:
    xs, ys = [], []
    for sig in SIG:
        row = wtab[(wtab.sigma_WP == sig) & (wtab.v == v)]
        if row.empty:
            continue
        xs.append(float(row["mean_sigma_r"]))
        ys.append(float(wp_S[sig].loc[v]))
    order = np.argsort(xs)
    xs, ys = np.array(xs)[order], np.array(ys)[order]
    ax.plot(xs, ys, "-o", color=vcol[v], ms=6, mfc=vcol[v], mec=vcol[v],
            lw=1.0, label=f"WP  $v$={v}")
    # same S, but at the interference-free ENTRY width (hollow diamonds)
    xe = [float(wtab[(wtab.sigma_WP == s) & (wtab.v == v)]["sigma_r_entry"])
          for s in SIG if not wtab[(wtab.sigma_WP == s) & (wtab.v == v)].empty]
    ye = [float(wp_S[s].loc[v]) for s in SIG]
    oe = np.argsort(xe)
    ax.plot(np.array(xe)[oe], np.array(ye)[oe], ":D", color=vcol[v], ms=5,
            mfc="none", mec=vcol[v], mew=1.0, lw=0.8)
ax.plot([], [], marker="D", ms=5, mfc="none", mec="0.3", ls=":",
        label="WP at entry $\\sigma_r$")

# classical: at x = sigma_pot = sigma_WP/sqrt2 (fixed; potential does not disperse)
for sig, path, filt in CLS:
    d = load_S(path, "S_B_Eabs", filt).set_index("v")["S_B_Eabs"]
    xpot = sig / np.sqrt(2)
    for v in VS:
        if v in d.index:
            ax.plot(xpot, d.loc[v], marker="s", ms=6, mfc="none",
                    mec=vcol[v], mew=1.3, ls="none")
ax.plot([], [], marker="s", ms=6, mfc="none", mec="0.3", ls="none",
        label="classical (at $\\sigma_\\mathrm{pot}$)")

ax.set_xscale("log")
ax.set_xlabel("effective width  (Bohr) —  WP $\\langle\\sigma_r\\rangle$ / cl. $\\sigma_\\mathrm{pot}$")
ax.set_ylabel(style.axis_label("stopping_power", symbol="$S_B$"))
ax.set_ylim(bottom=0)
ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="upper right")
ax.set_title("Stopping vs effective width: one decreasing trend", fontsize=9)
fig.savefig(f"{WP}/S_vs_effective_width.png", dpi=300, bbox_inches="tight")
print("wrote S_vs_effective_width.png")
plt.show()"""))

cells.append(new_markdown_cell(r"""## Conclusion

- The WP momentum spectrum shifts down and narrows during a run — real
  scattering/CAP filtering, distinct from free dispersion.
- Averaged over the run the effective width **inverts** the launch-$\sigma$
  ordering: the $\sigma_\text{WP}=0.5$ packet is actually the *widest*
  ($\langle\sigma_r\rangle$ largest), the $\sigma_\text{WP}=3$ packet the
  *narrowest*, because narrow packets disperse fastest.
- The interference-free **entry width** (centroid at the near slab face) confirms
  the dominant effect: $\sigma_\text{WP}=0.5$ has already dispersed to $\sim$11 Bohr
  on arrival, so its low stopping is a genuine on-the-way-in dispersion effect, not a
  late-time interference artefact.
- The finer $\sigma_\text{WP}=2.0$-vs-$3.0$ ordering is weaker: they arrive nearly
  tied ($\sim$4 Bohr), so it separates only under the full-run $\langle\sigma_r\rangle$
  (which includes the interference window) — read with caution.
- Re-expressed against $\langle\sigma_r\rangle$, the WP stopping trend runs the
  **same direction as the classical one** — narrower effective width $\to$ more
  stopping. The apparent sign paradox is largely an artefact of labelling a
  dispersing packet by its launch width.

**Caveats.** $P(p_z)$ is a two-moment (Gaussian) reconstruction — the full spectrum
was not exported for these runs. $\langle\sigma_r\rangle$ includes a transverse
contribution inflated by the periodic box once the packet delocalises;
$\langle\sigma_z\rangle$ (propagation axis) gives the same ordering with less
ambiguity. Both point to the same reconciliation."""))

nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
    "language_info": {"name": "python"},
}

ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
ep.preprocess(nb, {"metadata": {"path": HERE}})
with open(OUT_IPYNB, "w") as f:
    nbf.write(nb, f)
print(f"\nwrote+executed {OUT_IPYNB}")
