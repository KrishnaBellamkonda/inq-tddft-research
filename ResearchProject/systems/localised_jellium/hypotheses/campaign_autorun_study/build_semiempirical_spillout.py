#!/usr/bin/env python3
"""Build semiempirical_spillout.ipynb — quantify electron spill-out and locate the source of
the non-zero far-field plateau in the semi-empirical model.

Campaign: docs/campaigns/localised_jellium_parameter_study_2/ ("Semi empirical model for the
total system"). Data: scripts/semiempirical_spillout/runs/<tag>/results/density_gs_system.vti
(p2 GS SCF matrix). Reference = lz160 (Lz=160, N=82, w=0).

  A box/Lz : lz90 lz120 lz160 lz240   B soft-w : w1 w2 w4
  C confine: N164 N328                 D solver : es60

FINDING the metrics are built around (seen in the raw data): the interior vacuum density DECAYS
(physical tail), but sharp-edge (w=0) runs show a NEAR-EDGE PILE-UP (~4e-5 e/Bohr^3) right at the
open-z box boundary, which softening (w>0) removes. So two DISTINCT things are measured:
  - interior tail n_e at FIXED distances (20, 30 Bohr) — box-independent if physical;
  - near-edge density (within 4 Bohr of the box edge) — the boundary artifact.

Execute: python3 -m nbconvert --to notebook --execute --inplace semiempirical_spillout.ipynb (venv)
"""
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

OUT = Path(__file__).resolve().parent

PRE = r"""import sys, glob
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.special import erfc
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview import load_vti
try:
    from inqview.visualisation import style as _st; _st.apply()
except Exception: pass

HA = 27.211386
AREA = 50.0 * 50.0            # Bohr^2
HALF = 12.5                   # slab half-width (Bohr)
RUNBASE = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/semiempirical_spillout/runs"

RUNS = {  # tag -> (Lz, N, w, extra_states, periodicity)
    "lz90":  dict(Lz=90,  N=82,  w=0, es=20, per=2), "lz120": dict(Lz=120, N=82,  w=0, es=20, per=2),
    "lz160": dict(Lz=160, N=82,  w=0, es=20, per=2), "lz240": dict(Lz=240, N=82,  w=0, es=20, per=2),
    "w1":    dict(Lz=160, N=82,  w=1, es=20, per=2), "w2":    dict(Lz=160, N=82,  w=2, es=20, per=2),
    "w4":    dict(Lz=160, N=82,  w=4, es=20, per=2),
    "N164":  dict(Lz=160, N=164, w=0, es=20, per=2), "N328":  dict(Lz=160, N=328, w=0, es=20, per=2),
    "es60":  dict(Lz=160, N=82,  w=0, es=60, per=2),
    # p3 (fully periodic / PBC) counterparts — the open-z boundary test (Q3):
    "p3_lz90":  dict(Lz=90,  N=82, w=0, es=20, per=3),
    "p3_lz160": dict(Lz=160, N=82, w=0, es=20, per=3),
    "p3_lz240": dict(Lz=240, N=82, w=0, es=20, per=3),
}
REF = "lz160"
def n0_of(N):  return N / (AREA * 2.0 * HALF)

# Resilient to a partially-complete run matrix: keep only runs whose density VTI
# exists (so the notebook builds/executes before all 10 GS runs finish, and simply
# fills in as they complete on re-execution).
RUNS = {t: m for t, m in RUNS.items()
        if glob.glob(f"{RUNBASE}/{t}/results/density_gs_system/*.vti")}
def AV(tags):  return [t for t in tags if t in RUNS]     # available subset, order-preserving
if REF not in RUNS: REF = next(iter(RUNS))

def load_ne(tag, symmetrise=True):
    f = glob.glob(f"{RUNBASE}/{tag}/results/density_gs_system/*.vti")[0]
    v = load_vti(f, expect_centered_axis="z")
    ne = v.data.mean(axis=(0, 1)); z = v.z
    if symmetrise: ne = 0.5 * (ne + ne[::-1])
    return z, ne

def nplus_of(z, N, w):
    n0 = n0_of(N); d = np.abs(z)
    mask = np.where(d < HALF, 1.0, 0.0) if w <= 0 else 0.5 * erfc((d - HALF) / w)
    return n0 * mask

def phi_sheet(rho, z):
    dz = z[1] - z[0]
    return np.array([-2*np.pi * np.sum(rho * np.abs(zi - z)) * dz for zi in z])

def at(z, ne, Z):  return ne[np.argmin(np.abs(np.abs(z) - Z))]

def analyse(tag):
    m = RUNS[tag]; z, ne = load_ne(tag)
    npz = nplus_of(z, m["N"], m["w"]); dz = z[1] - z[0]
    rho = npz - ne
    phi = phi_sheet(rho, z); E = -np.gradient(phi, z) * HA          # eV/Bohr
    zmax = z.max()
    # interior tail at FIXED distances (box-independent if physical)
    tail20, tail30 = at(z, ne, 20.0), at(z, ne, 30.0)
    # near-edge pile-up: mean density within 4 Bohr of the open-z box edge
    edge_dens = ne[np.abs(z) > (zmax - 4.0)].mean()
    # spill: electrons beyond a fixed |z| cutoff
    Ne = np.sum(ne) * AREA * dz; Np = np.sum(npz) * AREA * dz
    spill = {Z: np.sum(ne[np.abs(z) > Z]) * AREA * dz for Z in (HALF, 15, 20, 25, 30)}
    # enclosed net-charge deficit at a FIXED Z0=25 (Q_eff), and the field it implies
    Z0 = 25.0
    Q25 = np.sum(rho[np.abs(z) <= Z0]) * AREA * dz         # = electrons beyond 25 (deficit)
    E_at25 = 2*np.pi*Q25/AREA * HA                          # eV/Bohr, semi-empirical field at Z0
    Emax_vac = np.max(np.abs(E[np.abs(z) > HALF]))          # peak vacuum field (near-edge)
    return dict(tag=tag, z=z, ne=ne, npz=npz, rho=rho, phi=phi, E=E, Ne=Ne, Np=Np,
                tail20=tail20, tail30=tail30, edge_dens=edge_dens, spill=spill,
                Q25=Q25, E_at25=E_at25, Emax_vac=Emax_vac, **m)

A = {t: analyse(t) for t in RUNS}
print("analysed:", ", ".join(A))
r = A[REF]
print(f"reference {REF}: N_e={r['Ne']:.2f} N_+={r['Np']:.2f}  tail@20={r['tail20']:.1e} "
      f"tail@30={r['tail30']:.1e}  edge={r['edge_dens']:.1e}  spill>25={r['spill'][25]:.3f}e "
      f"Q25={r['Q25']:.3f}e  E@25={r['E_at25']:.4f} eV/Bohr")"""

Q1 = r"""# Q1 -- how much charge sits outside the slab, and where (reference lz160).
ref = A[REF]; z, ne = ref["z"], ref["ne"]; dz = z[1]-z[0]
rows = []
for Z in [12.5, 15, 20, 25, 27, 30, 40, 50]:
    inside = np.sum(ne[np.abs(z) <= Z]) * AREA * dz
    rows.append(dict(Z_Bohr=Z, e_inside=round(inside,4), pct_inside=round(100*inside/ref["Ne"],3),
                     e_outside=round(ref["Ne"]-inside,4), pct_outside=round(100*(ref["Ne"]-inside)/ref["Ne"],4)))
q1 = pd.DataFrame(rows); print(q1.to_string(index=False))
Zg = np.linspace(0, z.max(), 400)
frac_in = np.array([np.sum(ne[np.abs(z) <= Z])*AREA*dz for Z in Zg]) / ref["Ne"]
Qnet   = np.array([np.sum(ref["rho"][np.abs(z) <= Z])*AREA*dz for Z in Zg])
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
ax[0].plot(Zg, 100*frac_in, lw=2); ax[0].axvline(HALF, ls="--", color="0.5", label="slab face 12.5")
ax[0].set_xlabel("|z| cutoff (Bohr)"); ax[0].set_ylabel("% electrons inside |z|<Z")
ax[0].set_title("Cumulative electron fraction (lz160)"); ax[0].legend(frameon=False, fontsize=8)
ax[1].plot(Zg, Qnet, lw=2, color="#c0392b"); ax[1].axhline(0, color="0.6", lw=0.6)
ax[1].axvline(HALF, ls="--", color="0.5")
ax[1].set_xlabel("|z| cutoff (Bohr)"); ax[1].set_ylabel("enclosed NET charge Q(|z|<Z) (e)")
ax[1].set_title("Net charge deficit vs Z")
fig.tight_layout(); plt.show()
print(f"\\n-> {q1.iloc[4]['pct_outside']:.3f}% of electrons beyond |z|=27; deficit "
      f"+{Qnet[np.argmin(abs(Zg-30))]:.3f} e at Z=30, closing to ~0 only at the box edge.")"""

PROF_LZ = r"""# PROFILE OVERLAY -- Lz sweep (all w=0). Interior tail vs near-edge pile-up.
_lz = AV(["lz90","lz120","lz160","lz240"])
fig, ax = plt.subplots(figsize=(8.4, 4.8))
for t in _lz:
    z, ne = A[t]["z"], A[t]["ne"]
    ax.semilogy(z, ne, lw=1.6, label=f"{t} (Lz={A[t]['Lz']})")
ax.plot(A[REF]["z"], A[REF]["npz"], "--", color="0.4", lw=1, label="n_+ background")
ax.axvspan(-HALF, HALF, color="0.92")
ax.set_xlabel("z (Bohr)"); ax.set_ylabel("n_e(z)  (e/Bohr³, log)")
ax.set_title("Lz sweep (w=0): interior tail collapses; each box piles up at its own edge")
ax.set_ylim(1e-16, 1e-2); ax.legend(frameon=False, fontsize=8, ncol=2); fig.tight_layout(); plt.show()
print("interior tail (should be ~Lz-independent if physical):")
print("  " + "   ".join(f"{t}: n_e(20)={A[t]['tail20']:.2e} n_e(30)={A[t]['tail30']:.2e}" for t in _lz))"""

PROF_W = r"""# PROFILE OVERLAY -- w sweep (all Lz=160). Does softening remove the edge pile-up?
_ws = AV(["lz160","w1","w2","w4"])
fig, ax = plt.subplots(figsize=(8.4, 4.8))
for t in _ws:
    z, ne = A[t]["z"], A[t]["ne"]
    ax.semilogy(z, ne, lw=1.6, label=f"{t} (w={A[t]['w']})")
ax.axvspan(-HALF, HALF, color="0.92")
ax.set_xlabel("z (Bohr)"); ax.set_ylabel("n_e(z)  (e/Bohr³, log)")
ax.set_title("w sweep (Lz=160): edge softness vs the near-edge pile-up")
ax.set_ylim(1e-16, 1e-2); ax.legend(frameon=False, fontsize=8); fig.tight_layout(); plt.show()
print("near-edge density (within 4 Bohr of the box edge):")
print("  " + "   ".join(f"{t}(w={A[t]['w']}): {A[t]['edge_dens']:.2e}" for t in _ws))"""

Q2 = r"""# Q2 -- the semi-empirical field for the reference, and the enclosed-charge identity.
ref = A[REF]; z = ref["z"]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
ax[0].plot(z, ref["phi"]*HA, lw=2); ax[0].axvspan(-HALF, HALF, color="0.9")
ax[0].set_xlabel("z (Bohr)"); ax[0].set_ylabel("φ (eV)"); ax[0].set_title("Semi-empirical potential φ(z)")
ax[1].plot(z, ref["E"], lw=2, color="#1b6ca8"); ax[1].axvspan(-HALF, HALF, color="0.9")
ax[1].axhline(0, color="0.6", lw=0.5)
ax[1].set_xlabel("z (Bohr)"); ax[1].set_ylabel("E_z (eV/Bohr)")
ax[1].set_title("Semi-empirical field E(z)")
fig.tight_layout(); plt.show()
print(f"deficit at Z0=25: Q25 = {ref['Q25']:.3f} e -> field 2πQ/A = {ref['E_at25']:.4f} eV/Bohr; "
      f"peak vacuum |E| = {ref['Emax_vac']:.4f} eV/Bohr (near the box edge).")"""

METRICS = r"""# MASTER METRIC TABLE -- all runs.
cols = []
for t in RUNS:
    a = A[t]
    cols.append(dict(run=t, Lz=a["Lz"], N=a["N"], w=a["w"], es=a["es"],
                     N_e=round(a["Ne"],2),
                     tail_z20=a["tail20"], tail_z30=a["tail30"], edge_dens=a["edge_dens"],
                     spill_gt25_pct=round(100*a["spill"][25]/a["Ne"],3),
                     Q25_e=round(a["Q25"],3), E_at25=round(a["E_at25"],4),
                     Emax_vac=round(a["Emax_vac"],4)))
pd.set_option("display.float_format", lambda v: f"{v:.2e}" if (abs(v)<1e-2 and v!=0) else f"{v:.3f}")
mt = pd.DataFrame(cols); print(mt.to_string(index=False)); pd.reset_option("display.float_format")"""

def sweep(tags, xkey, xlabel, title):
    return f"""# {title}
tags = AV({tags!r})
tab = pd.DataFrame([dict(run=t, {xkey}=A[t][{xkey!r}],
    tail_z20=A[t]['tail20'], edge_dens=A[t]['edge_dens'],
    spill_gt25_pct=round(100*A[t]['spill'][25]/A[t]['Ne'],3),
    Q25_e=round(A[t]['Q25'],3), Emax_vac=round(A[t]['Emax_vac'],4)) for t in tags])
pd.set_option("display.float_format", lambda v: f"{{v:.2e}}" if (abs(v)<1e-2 and v!=0) else f"{{v:.3f}}")
print(tab.to_string(index=False) if tags else "(no runs available for this sweep yet)")
pd.reset_option("display.float_format")
if tags:
    xs = [A[t][{xkey!r}] for t in tags]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    ax[0].plot(xs, [A[t]['edge_dens'] for t in tags], "o-", lw=2, color="#c0392b")
    ax[0].set_xlabel("{xlabel}"); ax[0].set_ylabel("near-edge density (e/Bohr³)")
    ax[0].set_title("boundary pile-up"); ax[0].set_yscale("log")
    ax[1].plot(xs, [A[t]['Q25'] for t in tags], "s-", lw=2)
    ax[1].set_xlabel("{xlabel}"); ax[1].set_ylabel("deficit Q(|z|<25) (e)")
    ax[1].set_title("spill charge beyond 25 Bohr")
    fig.tight_layout(); plt.show()"""

BOX  = sweep(["lz90","lz120","lz160","lz240"], "Lz", "Lz (Bohr)",
    "SWEEP A -- box size. Interior tail (n_e@20,@30) should be Lz-independent (physical). Watch "
    "whether the near-edge pile-up and the deficit Q(<25) grow with Lz (uniform floor) or not.")
SOFTW = sweep(["lz160","w1","w2","w4"], "w", "edge width w (Bohr)",
    "SWEEP B -- softer background edge. Does softening remove the near-edge pile-up and the deficit? "
    "(lz160 = w=0 baseline.)")
CONF  = sweep(["lz160","N164","N328"], "N", "N (electrons; n0 ∝ N)",
    "SWEEP C -- confinement / n0. Higher N deepens the well. Does the deficit/spill shrink?")
SOLVER = r"""# SWEEP D -- solver check: extra_states 20 -> 60, same geometry.
for t in AV(["lz160","es60"]):
    a = A[t]
    print(f"{t:6s} es={a['es']:3d}: tail@20={a['tail20']:.2e} edge={a['edge_dens']:.2e} "
          f"spill>25={a['spill'][25]:.4f}e Q25={a['Q25']:.3f}e")
if "es60" in RUNS and "lz160" in RUNS:
    print(f"\\nΔQ25 (es60 − es20) = {A['es60']['Q25'] - A['lz160']['Q25']:+.4f} e "
          f"-> ~0 means the deficit is not set by the empty-state count.")
else:
    print("(es60 not finished yet)")"""

W_FIELD = r"""# Q4 -- semi-empirical E(z) and phi(z) for the w SWEEP (empirical n_e + matched n_+ + sheet stack).
# phi and E use each run's OWN background n_+(w); this is the verification that softening the
# background edge removes the non-zero far field (not just the density pile-up).
_ws = AV(["lz160","w1","w2","w4"])
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
for t in _ws:
    z = A[t]["z"]
    ax[0].plot(z, A[t]["phi"]*HA, lw=1.7, label=f"{t} (w={A[t]['w']})")
    ax[1].plot(z, A[t]["E"],       lw=1.7, label=f"{t} (w={A[t]['w']})")
for a in ax: a.axvspan(-HALF, HALF, color="0.92"); a.set_xlabel("z (Bohr)")
ax[0].set_ylabel("φ (eV)"); ax[0].set_title("semi-empirical potential φ(z) — w sweep")
ax[1].set_ylabel("E_z (eV/Bohr)"); ax[1].axhline(0, color="0.6", lw=0.5)
ax[1].set_title("semi-empirical field E(z) — w sweep")
ax[0].legend(frameon=False, fontsize=8); ax[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); plt.show()
print("far-field |E| (mean over |z|∈[30,40]) and peak vacuum |E| per w:")
for t in _ws:
    z=A[t]["z"]; E=A[t]["E"]; win=(np.abs(z)>=30)&(np.abs(z)<=40)
    print(f"  {t} (w={A[t]['w']}): far |E|={np.mean(np.abs(E[win])):.4f}  peak vac |E|={A[t]['Emax_vac']:.4f} eV/Bohr")"""

PBC = r"""# Q3 -- open-z (p2) vs PBC (p3): does the near-edge pile-up disappear with periodic boundaries?
# In p3 there is no open boundary in z (z is periodic); the slab has image copies at +/-Lz, so the
# midpoint |z|=Lz/2 is the inter-slab region, NOT an open edge. Overlay n_e(z) at matched Lz.
pairs = [("lz90","p3_lz90"), ("lz160","p3_lz160"), ("lz240","p3_lz240")]
pairs = [(a,b) for a,b in pairs if a in RUNS and b in RUNS]
if not pairs:
    print("(p3 runs not finished yet)")
else:
    fig, ax = plt.subplots(1, len(pairs), figsize=(5.2*len(pairs), 4.4), squeeze=False)
    for k,(p2,p3) in enumerate(pairs):
        for t,ls in ((p2,"-"),(p3,"--")):
            z,ne = A[t]["z"], A[t]["ne"]
            ax[0][k].semilogy(z, ne, ls, lw=1.6, label=f"{t} (per={A[t]['per']})")
        ax[0][k].axvspan(-HALF, HALF, color="0.92"); ax[0][k].set_ylim(1e-16,1e-2)
        ax[0][k].set_xlabel("z (Bohr)"); ax[0][k].set_ylabel("n_e(z) (e/Bohr³, log)")
        ax[0][k].set_title(f"Lz={A[p2]['Lz']}: open-z vs PBC"); ax[0][k].legend(frameon=False, fontsize=8)
    fig.tight_layout(); plt.show()
    print("near-edge / inter-slab density (within 4 Bohr of |z|=Lz/2):")
    for p2,p3 in pairs:
        print(f"  Lz={A[p2]['Lz']}: p2(open-z)={A[p2]['edge_dens']:.2e}   p3(PBC)={A[p3]['edge_dens']:.2e}")"""

def build():
    cells = [
        new_markdown_cell(
            "# Semi-empirical spill-out — where does the non-zero far field come from?\n\n"
            "*Campaign `localised_jellium_parameter_study_2`. The semi-empirical field of "
            "jellium+background does not vanish far from the slab; A3 tied this to a +0.39 e "
            "enclosed-charge deficit and **inferred** (unverified) a numerical density floor. This "
            "notebook tests that with a p2 GS matrix (`scripts/semiempirical_spillout/`, both GPUs) "
            "and quantifies the spill. Built by `build_semiempirical_spillout.py`.*\n\n"
            "**Two distinct things are measured, because the raw data shows two distinct features:** "
            "(1) an **interior vacuum tail** that decays with distance (n_e at fixed 20/30 Bohr — "
            "Lz-independent if physical); (2) a **near-edge pile-up** right at the open-z box boundary "
            "(density within 4 Bohr of the edge). The discriminator: a uniform numerical floor is "
            "Lz-scaling and w/n₀-insensitive; a boundary artifact lives only at the edge and may be "
            "w-sensitive; a physical tail is bounded and confinement-sensitive. *Verdict is yours.*"),
        new_code_cell(PRE),
        new_markdown_cell("## Q1 — how much charge is outside the slab, and where\n"
            "Cumulative electron fraction and enclosed NET charge for the reference run; read off any region."),
        new_code_cell(Q1),
        new_markdown_cell("## Profiles — Lz sweep (w=0)\n"
            "Overlaid `n_e(z)` (log). If the interior tail collapses across Lz but each curve turns up "
            "at its own box edge, the far density is a decaying tail plus a boundary pile-up — not a "
            "single uniform floor."),
        new_code_cell(PROF_LZ),
        new_markdown_cell("## Profiles — w sweep (Lz=160)\n"
            "Does softening the background edge suppress the near-edge pile-up? (Raw data suggested yes.)"),
        new_code_cell(PROF_W),
        new_markdown_cell("## Q4 — semi-empirical E(z) and φ(z) for the w sweep\n"
            "The verification you asked for: build the field and potential from each w-run's empirical "
            "density (+ its own background) via the sheet stack. If a small `w` removes the non-zero "
            "far field (not just the density pile-up), it confirms softening is the fix."),
        new_code_cell(W_FIELD),
        new_markdown_cell("## Q2 — the semi-empirical field and the enclosed-charge identity\n"
            "The field at Z₀=25 equals `2π·Q(|z|<25)/A`, so the non-zero field IS the spill charge beyond 25 Bohr."),
        new_code_cell(Q2),
        new_markdown_cell("## Master metric table (all 10 runs)"),
        new_code_cell(METRICS),
        new_markdown_cell("## Sweep A — box size (Lz)"),
        new_code_cell(BOX),
        new_markdown_cell("## Q3 — open-z (p2) vs PBC (p3): is the pile-up a boundary artifact?\n"
            "Matched-Lz overlays of `n_e(z)` for open-z vs fully-periodic cells. If the near-edge "
            "pile-up is an open-z artifact it should vanish (or change character) under PBC, where z "
            "has no free boundary."),
        new_code_cell(PBC),
        new_markdown_cell("## Sweep B — softer background edge (w)"),
        new_code_cell(SOFTW),
        new_markdown_cell("## Sweep C — confinement / background density (n₀ ∝ N)"),
        new_code_cell(CONF),
        new_markdown_cell("## Sweep D — solver check (empty-state count)"),
        new_code_cell(SOLVER),
        new_markdown_cell(
            "## What the evidence implies (verdict is yours)\n"
            "- **Interior tail (n_e@20, @30) Lz-independent** ⇒ the vacuum density in the *interior* is a "
            "genuine, bounded, decaying surface tail — small, and not the plateau's cause.\n"
            "- **Near-edge pile-up present for w=0 and suppressed by w>0** ⇒ the plateau is driven by a "
            "density feature pinned to the open-z box boundary — an artifact the sharp background seeds "
            "and softening removes. Then your `w` instinct was right: `w` is the lever.\n"
            "- **Deficit Q(<25) grows ∝ Lz with a constant floor, unmoved by w/n₀/states** ⇒ a uniform "
            "numerical/thermal floor filling the vacuum instead.\n"
            "- **Deficit shrinks with higher n₀** ⇒ a physical, confinement-controlled tail.\n"
            "Read which pattern the tables/plots above actually show and fill in the conclusion."),
    ]
    n = new_notebook(); n.cells = cells
    n.metadata.kernelspec = {"name": "python3", "display_name": "Python 3"}
    p = OUT / "semiempirical_spillout.ipynb"
    nbf.write(n, str(p)); print("wrote", p.name); return p

if __name__ == "__main__":
    build()
    print("execute: python3 -m nbconvert --to notebook --execute --inplace semiempirical_spillout.ipynb (venv)")
