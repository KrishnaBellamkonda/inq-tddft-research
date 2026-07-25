#!/usr/bin/env python3
"""Phase-5 synthesis: build + EXECUTE the cylindrical-jellium analysis notebook.

Produces hypotheses/annular_sv/annular_sv_report.ipynb from the production sweep
(9 classical runs) + the WP rung. Panels:
  1. S(v) + β(r_s)  — the headline friction-slope-vs-wall-density result.
  2. Induced wall current current_z(t) — the flow→induced-current (hydrovoltaic) signature.
  3. Wake — induced density Δn radial+axial cut (shared colorbar; linear + log).
  4. WP-vs-classical — quantum wavepacket vs matched classical ghost (r_s=6, v=0.30).
Plus the hydrovoltaic / quantum-friction framing. Canonical inqview theme; VTIs
loaded in physical order (load_vti, NEVER fftshift). Defensive: panels whose data
is missing render a "pending" note instead of failing.

Run (AFTER production + WP complete):
  venv/bin/python3 build_report.py            # build + execute
  venv/bin/python3 build_report.py --no-exec  # build only
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

SYS = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/cylindrical_jellium")
HYP = SYS / "hypotheses/annular_sv"
OUT = HYP / "annular_sv_report.ipynb"
PER_RUN_MANIFEST = HYP / "per_run_manifest.json"

HEADER = r"""
import sys, math, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
try:
    from inqview.visualisation import style as _st
    _st.apply()
except Exception as e:
    print("theme:", e)

SYS = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/cylindrical_jellium")
SWEEP = SYS / "annular_sv"
DENS = {6: dict(L_z=48, N=24), 4: dict(L_z=28, N=48), 2: dict(L_z=10, N=136)}
VELS = [0.15, 0.30, 0.45]
DT = 0.020

def run_label(rs, v): return f"rs{rs}_v{v:.2f}".replace(".", "p")
def find(rundir, name):
    rundir = Path(rundir)
    return next(rundir.glob(f"**/{name}"), None)
"""

S_EXTRACT = r"""
def extract_S(rs, v, vfrac=0.85):
    # S(v0) = INITIAL drag at the launch velocity: the light electron decelerates,
    # so S = -d(KE_ion)/ds over the early near-constant-velocity window vz>=vfrac*v0.
    rd = SWEEP / run_label(rs, v)
    trk = find(rd, "electron_track.csv")
    if not trk: return None
    T = pd.read_csv(trk)
    z = T["z"].to_numpy(); vz = T["vz"].to_numpy(); ke = T["ke_ion_ha"].to_numpy()
    if len(z) < 10 or not np.all(np.isfinite(ke)): return None
    s = np.abs(z - z[0])
    for vf in (vfrac, 0.70, 0.50):
        sel = vz >= vf * v; sel[:max(2,int(0.03*len(s)))] = False
        if sel.sum() >= 20: break
    if sel.sum() < 8: return None
    ss, kk = s[sel], ke[sel]
    coef = np.polyfit(ss, kk, 1); S = -coef[0]
    yhat = np.polyval(coef, ss); dof = max(1, len(ss)-2); sxx = np.sum((ss-ss.mean())**2)
    serr = math.sqrt(np.sum((kk-yhat)**2)/dof / sxx) if sxx > 0 else float("nan")
    return dict(rs=rs, v=v, S=float(S), S_err=float(serr),
                v_mean=float(vz[sel].mean()), npts=int(sel.sum()))

rows = [r for rs in (6,4,2) for v in VELS if (r:=extract_S(rs, v))]
df = pd.DataFrame(rows)
print("S(v) table (initial-drag, early-window):"); print(df if len(df) else "(no production data yet)")
"""

SV_PLOT = r"""
betas = {}
fig, ax = plt.subplots(figsize=(7,5))
if len(df):
    for rs in (6,4,2):
        sub = df[df.rs == rs].sort_values("v")
        if len(sub):
            ax.errorbar(sub.v, sub.S, yerr=sub.S_err, marker="o", capsize=3, label=f"r_s={rs}")
            if len(sub) >= 2:
                b = np.polyfit(sub.v, sub.S, 1)[0]; betas[rs] = b
    ax.set_xlabel("v (a.u.)"); ax.set_ylabel("S (Ha/Bohr)")
    ax.set_title("Annular-tube electronic stopping S(v) vs wall r_s"); ax.legend()
    # annotate beta
    txt = "\n".join(f"β(r_s={rs}) = {b:.4f}" for rs, b in betas.items())
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="w", alpha=0.7))
else:
    ax.text(0.5, 0.5, "production data pending", ha="center")
plt.show()
print("beta(r_s) = dS/dv:", {k: round(v,5) for k,v in betas.items()})
mono = (len(betas) == 3 and (betas[6] < betas[4] < betas[2] or betas[6] > betas[4] > betas[2]))
print("beta(r_s) monotonic across r_s={6,4,2}:", mono)
"""

CURRENT_PLOT = r"""
# Panel 2: induced wall current current_z(t) — the flow -> induced-current signature.
fig, ax = plt.subplots(figsize=(7.5,4.5))
any_c = False
for rs in (6,4,2):
    for v in VELS:
        obs = find(SWEEP / run_label(rs, v), "observables.csv")
        if not obs: continue
        O = pd.read_csv(obs)
        if "current_z" in O:
            ax.plot(O["time_au"], O["current_z"], lw=1, label=f"r_s={rs}, v={v}")
            any_c = True
ax.set_xlabel("t (a.u.)"); ax.set_ylabel("wall current  I_z (a.u.)")
ax.set_title("Induced axial current as the projectile glides (hydrovoltaic flow→current)")
if any_c: ax.legend(fontsize=7, ncol=3)
else: ax.text(0.5,0.5,"current data pending",ha="center")
plt.show()
"""

WAKE_PLOT = r"""
# Panel 3: wake — induced density Δn from a representative run's density_delta VTIs.
from inqview import load_vti
def latest_delta(rs, v):
    rd = SWEEP / run_label(rs, v)
    dd = sorted(rd.glob("**/density_delta/*.vti"))
    return dd[len(dd)//2] if dd else None   # a mid-time frame (wake established)
pick = None
for rs, v in [(2,0.45),(4,0.45),(6,0.45),(2,0.30)]:
    p = latest_delta(rs, v)
    if p: pick = (rs, v, p); break
if pick:
    rs, v, p = pick
    vf = load_vti(str(p))                    # physical order, NO fftshift
    d = np.asarray(vf.data); x, y, z = np.asarray(vf.x), np.asarray(vf.y), np.asarray(vf.z)
    iy0 = np.argmin(np.abs(y))
    sl = d[:, iy0, :]                        # xz plane (y=0): the wake along the tube
    clim = np.nanpercentile(np.abs(sl), 99)
    fig, axes = plt.subplots(1, 2, figsize=(13,4.5))
    ext = [z.min(), z.max(), x.min(), x.max()]
    im0 = axes[0].imshow(sl, origin="lower", extent=ext, aspect="auto",
                         cmap="RdBu_r", vmin=-clim, vmax=clim)
    axes[0].set_title(f"Δn wake xz (linear), r_s={rs} v={v}"); axes[0].set_xlabel("z"); axes[0].set_ylabel("x")
    plt.colorbar(im0, ax=axes[0], fraction=0.046)
    with np.errstate(divide="ignore"):
        logd = np.sign(sl)*np.log10(np.abs(sl)+1e-12)
    im1 = axes[1].imshow(logd, origin="lower", extent=ext, aspect="auto", cmap="RdBu_r")
    axes[1].set_title("Δn wake (signed log10)"); axes[1].set_xlabel("z"); axes[1].set_ylabel("x")
    plt.colorbar(im1, ax=axes[1], fraction=0.046)
    plt.show()
    print(f"wake from {p.relative_to(SWEEP)}")
else:
    print("wake VTIs pending")
"""

WP_PLOT = r"""
# Panel 4: WP (quantum) vs matched classical ghost at r_s=6, v=0.30.
fig, ax = plt.subplots(figsize=(7.5,4.5))
def energy_vs_path_classical(rs, v):
    rd = SWEEP / run_label(rs, v)
    obs, trk = find(rd, "observables.csv"), find(rd, "electron_track.csv")
    if not (obs and trk): return None
    O = pd.read_csv(obs); T = pd.read_csv(trk)
    z_at = np.interp(O["time_au"], T["time_au"], T["z"])
    return np.abs(z_at - T["z"].to_numpy()[0]), (O["energy_total"]-O["energy_total"].iloc[0]).to_numpy()
cl = energy_vs_path_classical(6, 0.30)
if cl: ax.plot(cl[0], cl[1], label="classical ghost (r_s=6, v=0.30)")
# WP: energy vs WP centroid path (wp_real_space_stats if present, else time proxy)
wp_dir = SWEEP / "wp_rs6_v0p30"
obsw = find(wp_dir, "observables.csv")
rsw = find(wp_dir, "wp_real_space_stats.csv")
if obsw:
    Ow = pd.read_csv(obsw); de = (Ow["energy_total"]-Ow["energy_total"].iloc[0]).to_numpy()
    if rsw is not None:
        try:
            R = pd.read_csv(rsw); zcol = [c for c in R.columns if c.lower() in ("z","mean_z","z_mean","r_z")]
            zc = R[zcol[0]].to_numpy() if zcol else None
        except Exception: zc = None
    else: zc = None
    if zc is not None and len(zc) == len(de):
        s = np.abs(zc - zc[0]); ax.plot(s, de, "--", label="WP quantum (r_s=6, v=0.30)")
    else:
        ax.plot(Ow["time_au"]*0.30, de, "--", label="WP quantum (path≈v·t)")
ax.set_xlabel("projectile path s (Bohr)"); ax.set_ylabel("ΔE_system (Ha)")
ax.set_title("Quantum WP vs classical ghost: system energy gain vs path (slope = S)")
ax.legend(); plt.show()
if not obsw: print("WP run pending (Phase 4)")
"""

def _img(rel, cap, width):
    """Path-referenced <img> markdown (keeps the .ipynb small; figs travel beside it)."""
    if not rel:
        return f"_{cap}: figure missing_"
    return f"*{cap}*\n\n<img src=\"{rel}\" width=\"{int(width)}\" alt=\"{cap}\">"


def per_run_cells():
    """One deep-dive subsection per projectile run from per_run_manifest.json:
    density matrix + z-t carpets + high-value pipeline observables + initial-drag
    stopping. Path-referenced (figures pre-generated by build_per_run_figs.py)."""
    md = new_markdown_cell
    if not PER_RUN_MANIFEST.exists():
        return [md("## Per-run deep dives\n\n_Manifest `per_run_manifest.json` not "
                   "found — run `build_per_run_figs.py` first to generate the per-run "
                   "density matrices + observables._")]
    man = json.loads(PER_RUN_MANIFEST.read_text())
    W_GIF, W_PNG = 360, 560
    CAT_TTL = {"total": "Total system", "wp": "Wavepacket |ψ|²",
               "bath": "Bath / wall (total − WP)"}
    KIND_ORDER = ("density", "delta0", "dstep")
    cells = [md(
        "# Per-run deep dives\n\n"
        "One section per projectile run (9 classical + 1 wavepacket). Each shows the "
        "**matrix of density visualisations** — {density `n`, induced `Δn=n(t)−n(0)`, "
        "instantaneous-flux `Δn=n(t+dt)−n(t)`} × {total[, wavepacket, bath]} — as xz "
        "mid-plane GIFs (vertical dashed lines = tube wall radii |x|=5,13 Bohr; the "
        "projectile glides up the *z* axis through the hollow bore), with the **moving "
        "projectile overlaid** (cyan marker + trail, from that run's own track); the "
        "**z–t carpets** (projectile z(t) overlaid); the **high-value observables** "
        "(energy decomposition, induced current, momentum, KL drift); the "
        "**stopping-power extraction**; and the **FFT** of the collective-response "
        "observables (every FFT-driven observable goes through the `fourier-analysis` "
        "skill's audited 6-stage panel — the pipeline's raw FFT plots are not used).\n\n"
        "> **Stopping method (`stopping-power-extraction` skill; Correa 2018).** Geometry "
        "is **continuous traversal** (periodic tube, medium fills all z) → **Method A**: "
        "the headline is the free-intercept slope of the **electronic deposit** "
        "`ΔE_total(x)` (`energy_total` rises as the ion loses KE). The fit window is the "
        "**early `v ≥ 0.85·v0` segment** — the light free-Ehrenfest electron decelerates "
        "and stops, so S is read *at* v0, NOT over the skill's default 20%-time remainder "
        "(which would average S across the whole deceleration; "
        "`.claude/rules/light-projectile-stopping.md`). **Guards:** `N(t)≈const` (no CAP) "
        "and `ΔE_total ≈ −ΔKE_ion` (energy conservation; the `−dKE_ion/dx` kinetic channel "
        "is the independent cross-check). **Divergent channels (>10%) / poor r² are FLAGGED, "
        "not averaged — the verdict on a flagged number is yours.**\n\n"
        "> **Not stored by these runs (noted, not fabricated):** per-state KS energies "
        "(`state_energies.csv`) → no KS eigen-energy bar-GIFs; GS `eigenvalues.csv` not "
        "retrofitted → no GS KS-excitation decomposition; no E-field pipeline phase; only "
        "1D |k| momentum is saved → no 2D (k_z,k_⊥) scattering map.")]

    order = [k for k in man if man[k]["rtype"] == "classical"] + \
            [k for k in man if man[k]["rtype"] == "wp"]
    for label in order:
        d = man[label]
        rtype = d["rtype"]
        s = d.get("stopping") or {}
        lz = d.get("launch_z")
        lz_txt = f"z = {lz:g} Bohr" if lz is not None else "z0 (see run_summary)"
        if s:
            fl = (" ⚠ **FLAGS:** " + "; ".join(s["flags"])) if s.get("flags") else \
                 " (channels agree, clean fit)"
            s_txt = (f" Stopping **S(v0) = {s['S']:.4f} ± {s.get('S_err',0):.4f} Ha/Bohr** "
                     f"(PRIMARY ΔE_total slope, r² = {s.get('r2',float('nan')):.2f}); "
                     f"kinetic cross-check S = {s.get('S_kinetic',float('nan')):.4f} "
                     f"(ratio {s.get('ratio',float('nan')):.2f}); N drained "
                     f"{s.get('N_drained',float('nan'))*100:.1f}%; mean v = "
                     f"{s.get('v_mean',float('nan')):.3f}, {s.get('npts','?')} pts.{fl}")
        else:
            s_txt = ""
        cells.append(md(
            f"## Run `{label}` — r_s = {d['rs']}, v0 = {d['v0']} ({rtype})\n\n"
            f"{'Quantum electron wavepacket (σ_WP=0.5, drift k0=0.30)' if rtype=='wp' else 'Classical Gaussian-charge electron (σ_pot=0.354, m_e, Ehrenfest)'} "
            f"launched on-axis at {lz_txt}, gliding +z down the bore.{s_txt}"))

        # Density matrix, grouped by category
        cells.append(md("### Density matrix (xz mid-plane GIFs)"))
        by_cat = {}
        for cat, kind, rel, cap in d["matrix"]:
            by_cat.setdefault(cat, {})[kind] = (rel, cap)
        for cat in ("total", "wp", "bath"):
            if cat not in by_cat:
                continue
            cells.append(md(f"**{CAT_TTL.get(cat, cat)}**"))
            for kind in KIND_ORDER:
                if kind in by_cat[cat]:
                    rel, cap = by_cat[cat][kind]
                    cells.append(md(_img(rel, cap, W_GIF)))

        # z-t carpets
        if d["carpets"]:
            cells.append(md("### z–t carpets"))
            for cap, rel in d["carpets"]:
                cells.append(md(_img(rel, cap, W_PNG)))

        # Stopping power (classical) — skill-compliant: primary + cross-check + flags
        if s.get("path"):
            cells.append(md("### Stopping power (Method A — electronic deposit slope)"))
            cap = ("Left: v_z(t) decelerating, early v≥0.85·v0 window shaded. "
                   "Middle: ΔE_total(s) electronic deposit with the free-intercept "
                   "fit (PRIMARY S). Right: the two channels (ΔE_total vs −ΔKE_ion) "
                   "overlaid — agreement = energy conservation.")
            if s.get("flags"):
                cap += "  ⚠ " + "; ".join(s["flags"])
            cells.append(md(_img(s["path"], cap, 760)))

        # High-value pipeline observables (time domain)
        for group, figs in d["pipeline"].items():
            if not figs:
                continue
            cells.append(md(f"### {group}"))
            for name, rel in figs:
                cells.append(md(_img(rel, name, W_PNG)))

        # FFT — via the fourier-analysis skill (audited 6-stage panel)
        if d.get("fft"):
            cells.append(md(
                "### Collective response — FFT (fourier-analysis skill)\n"
                "Audited 6-stage panel (raw → mean-baseline → Hann → ×4 zero-pad → "
                "|FFT| linear+log, detrend overlaid). Replaces the pipeline's raw "
                "FFT plots; peak located inside the physical plasmon band, never by "
                "global argmax."))
            for cap, rel in d["fft"]:
                cells.append(md(_img(rel, cap, 760)))
    return cells


def build(execute=True):
    md = new_markdown_cell
    cells = [
        md("# Cylindrical (annular) jellium tube — projectile-down-bore stopping power\n\n"
           "*Executed synthesis of the `cylindrical-jellium-projectile` campaign. A charged\n"
           "projectile glides on-axis down the hollow bore of a PERIODIC annular jellium tube;\n"
           "we measure its electronic stopping power S(v) and the low-velocity friction slope\n"
           "β(r_s)=dS/dv as a function of the wall density r_s.*"),
        md("## North star — nanotube hydrovoltaics / quantum friction\n"
           "In solid–liquid systems (water in carbon nanotubes), flow induces an electronic\n"
           "current via wall-electron coupling. We model the wall as **jellium of variable r_s**\n"
           "(the TDDFT-PENN idea: different materials = different free-electron gases; Penn 1987,\n"
           "Matias 2025). This campaign is the first rung: a single annular tube + on-axis\n"
           "projectile, S(v) vs wall r_s. The induced axial current (Panel 2) is the TDDFT analogue\n"
           "of the flow→current hydrovoltaic signature (Kavokine 2022; Coquinot/Lizée PRX 2023).\n\n"
           "**Geometry (locked):** R_in=5, R_out=13 Bohr (8 Bohr wall), L_xy=40, dx=0.5, erfc w=1,\n"
           "periodic tube; per-density L_z={48,28,10} for r_s={6,4,2} (N={24,48,136}, exact\n"
           "neutrality). Projectile = classical electron (Gaussian σ_pot=0.354 UPF, m_e, Ehrenfest),\n"
           "v={0.15,0.30,0.45}; dt=0.020, LDA."),
        new_code_cell(HEADER),
        md("## 1. Stopping power S(v) and friction slope β(r_s)\n"
           "S = d(ΔE_system)/ds from a linear regression of the electronic total energy vs the\n"
           "projectile path s=|z−z0|, discarding the first 20% transient. β(r_s)=dS/dv is the\n"
           "low-velocity friction coefficient — the hypothesis quantity."),
        new_code_cell(S_EXTRACT),
        new_code_cell(SV_PLOT),
        md("## 2. Induced wall current (flow → induced current)\n"
           "The axial current current_z(t) carried by the wall electron gas as the projectile\n"
           "passes — the TDDFT signature underlying nanotube hydrovoltaics."),
        new_code_cell(CURRENT_PLOT),
        md("## 3. Wake structure\n"
           "The induced density Δn(r,t) = n(t) − n(0): the lagging wall-charge wake whose backward\n"
           "field retards the projectile (cylindrical Echenique–Ritchie image stopping). Shown as an\n"
           "xz cut (along the tube) in linear and signed-log scale."),
        new_code_cell(WAKE_PLOT),
        md("## 4. Quantum rung — wavepacket vs classical ghost (r_s=6)\n"
           "An electron wavepacket (σ_WP=0.5, drift k0=0.30) gliding down the bore vs its matched\n"
           "classical Gaussian-electron ghost (the v=0.30 production run). Quantifies the quantum\n"
           "effect on S. Caveat: the electron-as-cation proxy rests on charge-even S at leading order\n"
           "(Barkas = the charge-odd correction; Lindhard 1976)."),
        new_code_cell(WP_PLOT),
        *per_run_cells(),
        md("## Conclusions (PROVISIONAL)\n"
           "- S(v) and β(r_s) read off above; success = β(r_s) monotonic and resolved across r_s.\n"
           "- r_s=6 is a small gas (~24 e) → finite-size/shell effects; cross-check the β trend\n"
           "  against r_s=4 (per the campaign guard rails).\n"
           "- Results provisional until the wake plateau / L_z adequacy is confirmed per run and the\n"
           "  WP-vs-classical comparison is complete.\n"),
    ]
    nb = new_notebook(); nb.cells = cells
    nb.metadata.kernelspec = {"name": "python3", "display_name": "Python 3"}
    if execute:
        from nbclient import NotebookClient
        client = NotebookClient(nb, timeout=1200, kernel_name="python3",
                                resources={"metadata": {"path": str(HYP)}})
        try:
            client.execute(); print("notebook executed")
        except Exception as e:
            print(f"execution error (saving partial): {e}")
    nbf.write(nb, str(OUT)); print("wrote", OUT)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--no-exec", action="store_true")
    a = ap.parse_args(); build(execute=not a.no_exec)
