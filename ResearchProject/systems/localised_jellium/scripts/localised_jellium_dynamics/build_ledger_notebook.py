#!/usr/bin/env python3
"""Phase 1/2 analysis for campaign localised-jellium-dynamics-analysis.

Phase 1: completed classical-vs-WP energy ledger with the new U_proj_bg columns
(reuse existing wp_r*_p2; new cl_r* carry energy_proj_bg_{ideal,impl}). Phase 2:
E_external + proj_bg vs r_cut at r=20. Robust to a partially-complete run set.
Writes ledger.png, rcut.png, and ledger_rcut.ipynb into hypotheses/localised_jellium_dynamics/.
"""
import sys, glob, csv
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
try:
    from inqview.visualisation import style as _st; _st.apply()
except Exception: pass

HA = 27.211386
LJ = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium")
DYN = LJ/"scripts/localised_jellium_dynamics"
WP_DIR = LJ/"scripts/campaign_autorun/runs/h0_p2"          # existing wp_r*_p2
OUT = LJ/"hypotheses/localised_jellium_dynamics"; OUT.mkdir(parents=True, exist_ok=True)
E_GS_P2 = 60.38307052445239   # Ha, verified anchor
RADII = (4, 12, 20, 28, 36, 40)

def row0(pattern):
    fs = glob.glob(pattern, recursive=True)
    if not fs: return None
    r = list(csv.reader(open(fs[0]))); h, d = r[0], r[1]
    return {k: float(v) for k, v in zip(h, d)}

def E(d, k): return d.get(k, 0.0) if d else 0.0

# ---------------- Phase 1 ledger ----------------
rows = []
for r in RADII:
    wp = row0(str(WP_DIR/f"wp_r{r}_p2/**/observables.csv"))
    cl = row0(str(DYN/f"runs/p1/cl_r{r}/**/observables.csv"))
    if wp is None or cl is None: continue
    dE_WP = (E(wp,"energy_total")-E_GS_P2)*HA
    dE_CL = (E(cl,"energy_total")-E_GS_P2)*HA
    dKin  = (E(wp,"energy_kinetic")-E(cl,"energy_kinetic"))*HA
    dXC   = (E(wp,"energy_xc")-E(cl,"energy_xc"))*HA
    dHE   = ((E(wp,"energy_hartree")+E(wp,"energy_external"))
             -(E(cl,"energy_hartree")+E(cl,"energy_external")))*HA
    upb_i = E(cl,"energy_proj_bg_ideal")*HA
    upb_m = E(cl,"energy_proj_bg_impl")*HA
    wmcl  = dE_WP - dE_CL
    # closure: the part of WP−CL not in kin/xc/(H+E) — what U_proj_bg should explain
    resid = wmcl - dKin - dXC - dHE
    rows.append(dict(r=r, dE_WP=round(dE_WP,2), dE_CL=round(dE_CL,2), WPmCL=round(wmcl,2),
                     dKin=round(dKin,2), dXC=round(dXC,2), dHE=round(dHE,2),
                     U_proj_bg_ideal=round(upb_i,2), U_proj_bg_impl=round(upb_m,2),
                     resid=round(resid,2), resid_plus_upb=round(resid+upb_i,2)))
led = pd.DataFrame(rows)
led.to_csv(OUT/"ledger.csv", index=False)
print("Phase 1 ledger:\n", led.to_string(index=False) if len(led) else "(no P1 runs yet)")

if len(led):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    ax[0].plot(led.r, led.WPmCL, "o-", label="WP−CL")
    ax[0].plot(led.r, led.dKin, "s--", label="dKin")
    ax[0].plot(led.r, led.dHE, "^--", label="d(U_H+U_ext)")
    ax[0].plot(led.r, led.U_proj_bg_ideal, "d-", label="U_proj_bg (ideal)")
    ax[0].set_xlabel("r (Bohr)"); ax[0].set_ylabel("energy (eV)")
    ax[0].set_title("Ledger components vs r"); ax[0].legend(frameon=False, fontsize=8)
    ax[1].axhline(0, color="0.6", lw=0.6)
    ax[1].plot(led.r, led.resid, "o-", color="#c0392b", label="WP−CL − dKin − dXC − d(H+E)")
    ax[1].plot(led.r, led.resid_plus_upb, "s-", color="#2ca02c", label="residual + U_proj_bg(ideal)")
    ax[1].set_xlabel("r (Bohr)"); ax[1].set_ylabel("residual (eV)")
    ax[1].set_title("Does U_proj_bg close the residual?"); ax[1].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(OUT/"ledger.png", dpi=130); plt.close(fig)
    print("wrote", OUT/"ledger.png")

# ---------------- r_cut=50 vs r_cut=120 (=Lz) comparison ----------------
cmp_rows = []
for r in RADII:
    wp = row0(str(WP_DIR/f"wp_r{r}_p2/**/observables.csv"))
    c50 = row0(str(DYN/f"runs/p1/cl_r{r}/**/observables.csv"))
    c120 = row0(str(DYN/f"runs/p1_rc120/cl_r{r}/**/observables.csv"))
    if not (wp and c50 and c120): continue
    def dHE(c): return ((E(wp,"energy_hartree")+E(wp,"energy_external"))
                        -(E(c,"energy_hartree")+E(c,"energy_external")))*HA
    def iminusi(c): return (E(c,"energy_proj_bg_ideal")-E(c,"energy_proj_bg_impl"))*HA
    def completed(c): return dHE(c)-E(c,"energy_proj_bg_impl")*HA
    dexc = (E(wp,"energy_xc")-E(c120,"energy_xc"))*HA          # added WP XC (E_xc(WP)-E_xc(GS))
    cmp_rows.append(dict(r=r, dHE_rc50=round(dHE(c50),2), dHE_rc120=round(dHE(c120),2),
                         ideal_minus_impl_rc50=round(iminusi(c50),2), ideal_minus_impl_rc120=round(iminusi(c120),2),
                         completed_rc50=round(completed(c50),2), completed_rc120=round(completed(c120),2),
                         dE_xc=round(dexc,2), completed_plus_xc=round(completed(c120)+dexc,2)))
cmp = pd.DataFrame(cmp_rows)
if len(cmp):
    cmp.to_csv(OUT/"rcut50_vs_120.csv", index=False)
    print("\nr_cut 50 vs 120 comparison:\n", cmp.to_string(index=False))

# ---------------- Phase 2 r_cut ----------------
# WP reference at r=20 (its full-reach electrostatic channel E_ext+E_H). The classical
# E_ext+E_H rises with r_cut and crosses this line at the WP's effective radial cutoff.
wp20 = row0(str(WP_DIR/"wp_r20_p2/**/observables.csv"))
WP_EH = (E(wp20,"energy_external")+E(wp20,"energy_hartree"))*HA if wp20 else float("nan")
rc_rows = []
for rc in (10, 20, 30, 40, 50):
    cl = row0(str(DYN/f"runs/p2/cl_r20_rc{rc}/**/observables.csv"))
    if cl is None: continue
    eh = (E(cl,"energy_external")+E(cl,"energy_hartree"))*HA
    rc_rows.append(dict(r_cut=rc, E_external=round(E(cl,"energy_external")*HA,2),
                        E_hartree=round(E(cl,"energy_hartree")*HA,2),
                        E_ext_plus_H=round(eh,2), E_total=round(E(cl,"energy_total")*HA,2),
                        proj_bg_ideal=round(E(cl,"energy_proj_bg_ideal")*HA,2),
                        proj_bg_impl=round(E(cl,"energy_proj_bg_impl")*HA,2),
                        WP_ref_E_ext_plus_H=round(WP_EH,2)))
rcut = pd.DataFrame(rc_rows)
rcut.to_csv(OUT/"rcut.csv", index=False)
print("\nPhase 2 r_cut:\n", rcut.to_string(index=False) if len(rcut) else "(no P2 runs yet)")

if len(rcut):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    ax[0].plot(rcut.r_cut, rcut.E_external, "o-", label="E_external")
    ax[0].plot(rcut.r_cut, rcut.E_total, "s--", label="E_total")
    ax[0].set_xlabel("r_cut (Bohr)"); ax[0].set_ylabel("energy (eV)")
    ax[0].set_title("Energy vs projectile r_cut (r=20)"); ax[0].legend(frameon=False, fontsize=8)
    ax[1].plot(rcut.r_cut, rcut.proj_bg_ideal, "d-", label="proj_bg ideal (should be flat)")
    ax[1].plot(rcut.r_cut, rcut.proj_bg_impl, "^-", label="proj_bg impl (r_cut-dependent)")
    ax[1].set_xlabel("r_cut (Bohr)"); ax[1].set_ylabel("E_proj_bg (eV)")
    ax[1].set_title("ideal r_cut-invariant vs impl"); ax[1].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(OUT/"rcut.png", dpi=130); plt.close(fig)
    print("wrote", OUT/"rcut.png")

# ---------------- assemble a lightweight notebook ----------------
try:
    import nbformat as nbf
    from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
    cells = [
        new_markdown_cell("# Localised-jellium-dynamics — Phase 1 ledger + Phase 2 r_cut\n\n"
            "*Completed classical-vs-WP energy ledger with the `U_proj_bg` diagnostic (Phase 1) and "
            "the projectile radial-cutoff sweep (Phase 2). Built by `build_ledger_notebook.py`; runs "
            "under `scripts/localised_jellium_dynamics/runs/`.*"),
        new_markdown_cell("## Phase 1 — completed ledger (eV)\n"
            "Per radius r (Bohr from the slab face): `dE_WP`/`dE_CL` vs the bare GS; `WPmCL` their "
            "difference; the per-component WP−CL splits `dKin` (WP localisation), `dXC`, `dHE = "
            "d(U_H+U_ext)` (electrostatic); and the two new projectile↔background columns "
            "`U_proj_bg_ideal` (∫n_proj·v_bg, true charge) and `U_proj_bg_impl` (−∫n₊·v_ion, as-implemented)."),
        new_code_cell("import pandas as pd, numpy as np\nimport matplotlib.pyplot as plt\n"
            "d = pd.read_csv('ledger.csv')\nd"),
        new_markdown_cell("### The two `U_proj_bg` estimates vs r\n"
            "Both grow ~linearly with r (the projectile–background interaction changes with distance). "
            "`ideal` is positive (open-z gauge on `v_bg`); `impl` is negative (UPF reference). Absolute "
            "signs are gauge-dependent — only differences are physical."),
        new_code_cell(
            "fig, ax = plt.subplots(figsize=(6.6,4.2))\n"
            "ax.plot(d.r, d.U_proj_bg_ideal, 'o-', label='U_proj_bg ideal  (∫ n_proj·v_bg)')\n"
            "ax.plot(d.r, d.U_proj_bg_impl,  's-', label='U_proj_bg impl   (−∫ n₊·v_ion)')\n"
            "ax.axhline(0, color='0.6', lw=0.6)\n"
            "ax.set_xlabel('r (Bohr from slab face)'); ax.set_ylabel('U_proj_bg (eV)')\n"
            "ax.set_title('Projectile–background Coulomb energy vs r (both methods)')\n"
            "ax.legend(frameon=False, fontsize=8); fig.tight_layout(); plt.show()"),
        new_markdown_cell("### Does U_proj_bg close the WP−classical electrostatic gap?\n"
            "The classical run is MISSING the ghost↔background term, so its electrostatic channel is "
            "incomplete. Completing it means subtracting `U_proj_bg` from the raw gap `dHE`: if the "
            "completed gap `dHE − U_proj_bg` becomes **r-independent (flat)**, the missing term accounts "
            "for the gap's r-dependence. The `impl` form is gauge-matched to INQ's classical `E_external` "
            "(both use `v_ion`); `ideal` uses `v_bg` (a different gauge)."),
        new_code_cell(
            "d['dHE_minus_ideal'] = d.dHE - d.U_proj_bg_ideal\n"
            "d['dHE_minus_impl']  = d.dHE - d.U_proj_bg_impl\n"
            "sp = lambda c: d[c].max()-d[c].min()\n"
            "print(f'spread across r (max-min):')\n"
            "print(f'  raw dHE                : {sp(\"dHE\"):7.1f} eV')\n"
            "print(f'  dHE - U_proj_bg_impl   : {sp(\"dHE_minus_impl\"):7.1f} eV  (flat = closed)')\n"
            "print(f'  dHE - U_proj_bg_ideal  : {sp(\"dHE_minus_ideal\"):7.1f} eV')\n"
            "fig, ax = plt.subplots(1,2, figsize=(12,4.3))\n"
            "ax[0].plot(d.r, d.dHE, 'o-', color='#c0392b', label='raw dHE = d(U_H+U_ext)')\n"
            "ax[0].plot(d.r, d.U_proj_bg_impl, 's--', label='U_proj_bg impl')\n"
            "ax[0].set_xlabel('r (Bohr)'); ax[0].set_ylabel('energy (eV)')\n"
            "ax[0].set_title('Raw electrostatic gap vs the missing term'); ax[0].legend(frameon=False, fontsize=8)\n"
            "ax[1].axhline(0, color='0.6', lw=0.6)\n"
            "ax[1].plot(d.r, d.dHE_minus_impl,  'o-', color='#2ca02c', label=f'dHE − impl  (spread {sp(\"dHE_minus_impl\"):.1f} eV)')\n"
            "ax[1].plot(d.r, d.dHE_minus_ideal, '^-', color='#1f77b4', label=f'dHE − ideal (spread {sp(\"dHE_minus_ideal\"):.1f} eV)')\n"
            "ax[1].axvline(25.0, ls=':', color='0.5'); ax[1].axvline(37.5, ls=':', color='0.7')\n"
            "ax[1].text(25.2, ax[1].get_ylim()[1], 'r=25\\n(ghost stops\\nspanning slab)', fontsize=7, va='top', color='0.4')\n"
            "ax[1].set_xlabel('r (Bohr)'); ax[1].set_ylabel('completed gap (eV)')\n"
            "ax[1].set_title('Completed electrostatic gap: flat ⇒ closed'); ax[1].legend(frameon=False, fontsize=8)\n"
            "fig.tight_layout(); plt.show()"),
        new_markdown_cell("**Reading it:** raw `dHE` spans ~171 eV across r; `dHE − U_proj_bg_impl` is flat "
            "to ~2 eV — the as-implemented ghost↔background term accounts for essentially the entire "
            "r-dependence of the electrostatic WP−CL difference. The residual constant is the small "
            "non-r-dependent offset (WP quantum electrostatic self-energy + charged-cell gauge: the WP "
            "cell is net −1 e, the classical cell is neutral). *Verdict is yours.*"),
        new_markdown_cell("### The `ideal` drift past r≈25 is the pseudopotential-cutoff geometry\n"
            "The Phase-1 ghost uses the FULL UPF (**r_cut = 50 Bohr**). With slab thickness 25, the ghost "
            "spans the whole slab only while `r + 25 < 50`, i.e. **r < 25**; it stops reaching the slab "
            "centre at r = 37.5. Past r≈25 the truncated potential no longer covers the far slab, so "
            "`E_external` (and `impl`, which uses the *same* cut `v_ion`) drop the far-slab interaction — "
            "but `ideal` (true un-cut charge) keeps counting it, so `dHE − ideal` drifts. The slope of "
            "raw `dHE` breaks at exactly this radius:"),
        new_code_cell(
            "sl = np.diff(d.dHE.values)/np.diff(d.r.values)\n"
            "print('slope of dHE between consecutive radii (eV/Bohr):')\n"
            "for i in range(len(sl)): print(f'  r={int(d.r.values[i]):>2}->{int(d.r.values[i+1]):>2}: {sl[i]:5.2f}   '\n"
            "      + ('(ghost spans slab)' if d.r.values[i+1] <= 25 else '(ghost NO LONGER spans slab)'))\n"
            "print('\\n-> constant ~5.4 eV/Bohr for r<=28, then breaks (3.9, 2.6) past the r=25 geometric threshold.')"),
        new_markdown_cell("## r_cut = 50 vs r_cut = 120 (= Lz): removing the truncation\n"
            "The Phase-1 ghost was truncated at 50 Bohr, so past r=25 it stopped spanning the slab. "
            "Re-running with an extended UPF (**r_cut = 120 = Lz**, pure Coulomb tail) makes the ghost "
            "span the whole box at every radius. Two tests: (1) does the `ideal − impl` drift flatten? "
            "(2) does the completed WP−CL gap `dHE − impl` become **r-independent**?"),
        new_code_cell("c = pd.read_csv('rcut50_vs_120.csv')\nc"),
        new_code_cell(
            "fig, ax = plt.subplots(1,2, figsize=(12,4.3))\n"
            "ax[0].plot(c.r, c.ideal_minus_impl_rc50,  'o-', label='r_cut=50  (drifts past r=25)')\n"
            "ax[0].plot(c.r, c.ideal_minus_impl_rc120, 's-', label='r_cut=120 (flat)')\n"
            "ax[0].axvline(25, ls=':', color='0.6'); ax[0].set_xlabel('r (Bohr)'); ax[0].set_ylabel('ideal − impl (eV)')\n"
            "ax[0].set_title('ideal − impl: truncation drift removed by r_cut=120'); ax[0].legend(frameon=False, fontsize=8)\n"
            "s50=c.completed_rc50.max()-c.completed_rc50.min(); s120=c.completed_rc120.max()-c.completed_rc120.min()\n"
            "ax[1].plot(c.r, c.completed_rc50,  'o-', label=f'r_cut=50  (spread {s50:.2f} eV)')\n"
            "ax[1].plot(c.r, c.completed_rc120, 's-', label=f'r_cut=120 (spread {s120:.2f} eV)')\n"
            "ax[1].set_xlabel('r (Bohr)'); ax[1].set_ylabel('completed WP−CL gap  dHE − impl (eV)')\n"
            "ax[1].set_title('Completed gap: r_cut=120 is r-INDEPENDENT'); ax[1].legend(frameon=False, fontsize=8)\n"
            "fig.tight_layout(); plt.show()\n"
            "print(f'completed-gap spread: r_cut=50 -> {s50:.2f} eV ;  r_cut=120 -> {s120:.2f} eV (flat)')\n"
            "print(f'ideal-impl: r_cut=50 drifts {c.ideal_minus_impl_rc50.min():.0f}->{c.ideal_minus_impl_rc50.max():.0f}; '\n"
            "      f'r_cut=120 constant ~{c.ideal_minus_impl_rc120.mean():.0f} eV')"),
        new_markdown_cell("**What this shows.** With `r_cut = 120` the `ideal − impl` offset is **dead-constant** "
            "across all radii (the truncation drift is gone), and the completed WP−CL gap `dHE − impl` is "
            "**flat to ~0.1 eV** (vs 1.8 eV at r_cut=50) — proving the r-dependence past r=25 was purely the "
            "pseudopotential truncation. **Caveat (the tell):** the absolute constants *shifted* (ideal−impl "
            "274→659 eV; gap 13.5→7.4 eV) because r_cut=120 ≫ Lx=Ly=50 also over-wraps the lateral periodic "
            "images and extends the Coulomb tail — a roughly uniform additive change, not an r-dependent one. "
            "So r_cut=Lz cleanly removes the *r-dependent* artifact (shape), but the *absolute* residual still "
            "carries the charged-cell gauge + lateral over-wrap and is not yet a pure quantum number. *Verdict is yours.*"),
        new_markdown_cell("## Phase 2 — r_cut sweep at r=20 (eV)\n"
            "`proj_bg_ideal` is flat vs r_cut (true charge, r_cut-invariant); the r_cut effect lives in "
            "`E_external`; `proj_bg_impl` varies (truncated potential)."),
        new_code_cell("rc = pd.read_csv('rcut.csv')\nrc"),
        new_code_cell(
            "fig, ax = plt.subplots(1,2, figsize=(12,4.3))\n"
            "ax[0].plot(rc.r_cut, rc.E_external, 'o-', label='E_external'); ax[0].plot(rc.r_cut, rc.E_total, 's--', label='E_total')\n"
            "ax[0].set_xlabel('r_cut (Bohr)'); ax[0].set_ylabel('energy (eV)'); ax[0].set_title('Energy vs r_cut (r=20)')\n"
            "ax[0].legend(frameon=False, fontsize=8)\n"
            "ax[1].plot(rc.r_cut, rc.proj_bg_ideal, 'd-', label='proj_bg ideal (flat)')\n"
            "ax[1].plot(rc.r_cut, rc.proj_bg_impl, '^-', label='proj_bg impl (r_cut-dependent)')\n"
            "ax[1].set_xlabel('r_cut (Bohr)'); ax[1].set_ylabel('E_proj_bg (eV)'); ax[1].set_title('ideal r_cut-invariant vs impl')\n"
            "ax[1].legend(frameon=False, fontsize=8); fig.tight_layout(); plt.show()"),
        new_markdown_cell("### Effective radial cutoff of the wavepacket\n"
            "The WP interacts with the system at full reach; the classical `E_ext+E_H` *rises* with r_cut "
            "as the ghost sees more of the slab. Where the classical curve crosses the **WP reference** "
            "(dashed line = WP `E_ext+E_H` at r=20) is the r_cut a classical ghost needs to match the WP "
            "— its **effective radial cutoff**."),
        new_code_cell(
            "wp_ref = rc.WP_ref_E_ext_plus_H.iloc[0]\n"
            "x = rc.r_cut.values.astype(float); y = rc.E_ext_plus_H.values.astype(float)\n"
            "reff = float(np.interp(wp_ref, y, x)) if (y.min() <= wp_ref <= y.max()) else float('nan')\n"
            "fig, ax = plt.subplots(figsize=(7.0,4.6))\n"
            "ax.plot(x, y, 'o-', color='#1b6ca8', label='classical E_ext+E_H (r=20)')\n"
            "ax.axhline(wp_ref, ls='--', color='#c0392b', label=f'WP reference E_ext+E_H = {wp_ref:.1f} eV')\n"
            "if np.isfinite(reff):\n"
            "    ax.axvline(reff, ls=':', color='0.5')\n"
            "    ax.plot([reff],[wp_ref],'k*',ms=12)\n"
            "    ax.annotate(f'effective r_cut ≈ {reff:.0f} Bohr', (reff, wp_ref), textcoords='offset points',\n"
            "                xytext=(8,-14), fontsize=9)\n"
            "ax.set_xlabel('classical projectile r_cut (Bohr)'); ax.set_ylabel('E_ext + E_H (eV)')\n"
            "ax.set_title('NAIVE effective r_cut (raw, incomplete classical — CONFOUNDED)')\n"
            "ax.legend(frameon=False, fontsize=8); fig.tight_layout(); plt.show()\n"
            "print(f'NAIVE effective radial cutoff ≈ {reff:.1f} Bohr  — but this compares an INCOMPLETE classical')"),
        new_markdown_cell("### Corrected comparison — complete the classical first (COMPLETED)\n"
            "The raw classical `E_ext+E_H` is MISSING the ghost↔background term, so crossing it with the WP "
            "conflates the missing term with a reach effect. The like-for-like comparison adds `U_proj_bg_impl` "
            "to the classical: **(E_ext_WP+E_H_WP) − (E_ext_CL+E_H_CL+U_proj_bg_impl)**. If this is flat vs "
            "r_cut, the WP−CL difference is a **quantum constant**, not a geometric effective-cutoff."),
        new_code_cell(
            "comp = rc.E_ext_plus_H + rc.proj_bg_impl            # completed classical electrostatic\n"
            "diff = rc.WP_ref_E_ext_plus_H - comp                # WP - completed classical\n"
            "print('WP - (E_ext_CL+E_H_CL+U_proj_bg_impl) vs r_cut (eV):')\n"
            "print('  ' + '  '.join(f'{rc.r_cut[i]:.0f}:{diff[i]:.2f}' for i in range(len(rc))))\n"
            "print(f'  spread = {diff.max()-diff.min():.2f} eV  (~flat = QUANTUM constant, not a reach)')\n"
            "fig, ax = plt.subplots(1,2, figsize=(12,4.3))\n"
            "ax[0].plot(rc.r_cut, comp, 'o-', color='#2ca02c', label='completed classical E_ext+E_H+impl')\n"
            "ax[0].axhline(rc.WP_ref_E_ext_plus_H.iloc[0], ls='--', color='#c0392b', label='WP reference')\n"
            "ax[0].set_xlabel('r_cut (Bohr)'); ax[0].set_ylabel('E_ext+E_H (eV)')\n"
            "ax[0].set_title('Completed classical is ~flat (no crossing)'); ax[0].legend(frameon=False, fontsize=8)\n"
            "ax[1].plot(rc.r_cut, diff, 'o-', color='#1b6ca8')\n"
            "ax[1].set_xlabel('r_cut (Bohr)'); ax[1].set_ylabel('WP − completed classical (eV)')\n"
            "ax[1].set_title(f'WP−CL electrostatic diff ≈ {diff.mean():.1f} eV (quantum, r_cut-independent)')\n"
            "fig.tight_layout(); plt.show()"),
        new_markdown_cell("**Corrected reading:** once the classical side is completed with `U_proj_bg_impl`, the "
            "WP−classical electrostatic difference is a nearly-constant **~14.5 eV**, independent of r_cut. "
            "That constant is the WP's genuine **quantum electrostatic self-energy** (a spread −1 e charge vs a "
            "point-like ghost) — NOT a geometric reach. So there is no clean 'effective radial cutoff' for the "
            "WP; the ~31 Bohr crossing above was an artifact of comparing an incomplete classical. The small "
            "residual r_cut-slope (~1.6 eV) is the only genuinely cutoff-sensitive piece. *Verdict is yours.*"),
        new_markdown_cell("## Gauge analysis — how much of the residual is charged-cell artifact?\n"
            "The completed WP−classical residual (~7.4 eV) survives `d(E_H+E_ext)` because the two cells "
            "have DIFFERENT net charge (WP cell net −1 e, classical cell neutral), so the 2D-periodic "
            "charged-cell (G=0) gauge does not cancel. Three independent estimates of that gauge:\n"
            "1. **self-Hartree probe** — WP charge's self-energy, p2 kernel vs free space;\n"
            "2. **neutralized scaling** — self-energy of (Gaussian −1 + uniform +1) vs box size Lx, Lz;\n"
            "3. **Makov–Payne** — analytic monopole finite-size term q²α/(2L)."),
        new_code_cell(
            "import numpy as np\n"
            "HA=27.211386; sig=0.5/np.sqrt(2.0)   # sigma_rho=0.354, |q|=1\n"
            "def poisson_p2(rho,x,z):\n"
            "    nx=len(x); dx=x[1]-x[0]; dz=z[1]-z[0]\n"
            "    rk=np.fft.fftn(rho,axes=(0,1)); k=2*np.pi*np.fft.fftfreq(nx,dx)\n"
            "    G=np.sqrt(k[:,None]**2+k[None,:]**2); absdz=np.abs(z[:,None]-z[None,:])\n"
            "    phik=np.empty_like(rk); Gq=np.round(G,10)\n"
            "    for g in np.unique(Gq):\n"
            "        m=Gq==g; cols=rk[m,:]\n"
            "        K=(-2*np.pi*absdz*dz) if g==0.0 else ((2*np.pi/g)*np.exp(-g*absdz)*dz)\n"
            "        phik[m,:]=cols@K.T\n"
            "    return np.real(np.fft.ifftn(phik,axes=(0,1)))\n"
            "def gaussian(LX,LZ,sp=0.5):\n"
            "    nx=int(round(LX/sp)); nz=int(round(LZ/sp))\n"
            "    x=(np.arange(nx)-nx//2)*sp; z=(np.arange(nz)-nz//2)*sp\n"
            "    X,Y,Z=x[:,None,None],x[None,:,None],z[None,None,:]\n"
            "    nG=np.exp(-(X**2+Y**2+Z**2)/(2*sig*sig)); nG/=(nG.sum()*sp**3)\n"
            "    return nG,x,z,sp\n"
            "# (1) self-Hartree: free space (analytic) vs p2 (bare charge)\n"
            "E_free=1.0/(2*sig*np.sqrt(np.pi))*HA\n"
            "nG,x,z,sp=gaussian(50,120); E_p2_bare=0.5*np.sum(nG*poisson_p2(nG,x,z))*sp**3*HA\n"
            "# (2) neutralized self-energy vs box\n"
            "def neut(LX,LZ):\n"
            "    nG,x,z,sp=gaussian(LX,LZ); rho=nG-1.0/(LX*LX*LZ)\n"
            "    return 0.5*np.sum(rho*poisson_p2(rho,x,z))*sp**3*HA\n"
            "Lx_scan=[40,50,75,100]; Ex=[neut(L,120) for L in Lx_scan]\n"
            "Lz_scan=[90,120,160,240]; Ez=[neut(50,L) for L in Lz_scan]\n"
            "slope,icpt=np.polyfit(1/np.array(Lx_scan), Ex, 1)\n"
            "gauge_lat=Ex[1]-icpt\n"
            "# (3) Makov-Payne\n"
            "E_MP=lambda L: 2.837/(2*L)*HA\n"
            "print(f'(1) self-Hartree: free={E_free:.2f} eV, p2(bare)={E_p2_bare:.2f} eV -> shift ~{E_p2_bare-E_free:+.2f} eV')\n"
            "print(f'(2) neutralized: lateral gauge (Lx=50 - L->inf intercept) = {gauge_lat:+.2f} eV;'\n"
            "      f' z-growth {Ez[0]:.2f}->{Ez[-1]:.2f} eV over Lz 90->240')\n"
            "print(f'(3) Makov-Payne @L=50: {E_MP(50):.2f} eV')\n"
            "print(f'--> all three ~1 eV; residual (~7.4 eV) is mostly PHYSICAL, not gauge.')"),
        new_code_cell(
            "fig, ax = plt.subplots(1,3, figsize=(15,4.2))\n"
            "ax[0].plot(1/np.array(Lx_scan), Ex, 'o-')\n"
            "ax[0].plot(1/np.array(Lx_scan), slope*(1/np.array(Lx_scan))+icpt, '--', color='0.6', label=f'intercept {icpt:.2f} eV')\n"
            "ax[0].set_xlabel('1/Lx (1/Bohr)'); ax[0].set_ylabel('neutralized self-energy (eV)')\n"
            "ax[0].set_title('Method 2: lateral scaling (flat -> tiny gauge)'); ax[0].legend(frameon=False, fontsize=8)\n"
            "ax[1].plot(Lz_scan, Ez, 's-', color='#c0392b')\n"
            "ax[1].set_xlabel('Lz (Bohr)'); ax[1].set_ylabel('neutralized self-energy (eV)')\n"
            "ax[1].set_title('Method 2: open-z scaling (~1 eV growth)')\n"
            "labels=['(1) self-H\\np2-free','(2) lateral','(2) z-growth','(3) Makov-\\nPayne']\n"
            "vals=[E_p2_bare-E_free, abs(gauge_lat), Ez[-1]-Ez[0], E_MP(50)]\n"
            "ax[2].bar(labels, vals, color=['#1f77b4','#2ca02c','#2ca02c','#7f7f7f'])\n"
            "ax[2].axhline(1.0, ls=':', color='0.5'); ax[2].set_ylabel('gauge estimate (eV)')\n"
            "ax[2].set_title('All three estimates: gauge ~1 eV vs 7.4 eV residual')\n"
            "for i,v in enumerate(vals): ax[2].text(i, v+0.03, f'{v:.2f}', ha='center', fontsize=8)\n"
            "fig.tight_layout(); plt.show()"),
        new_markdown_cell("## ★ MAIN RESULT — the classical-vs-quantum residual (~7.4 eV)\n"
            "With the pseudopotential extended to span the box (r_cut=120), the completed WP−classical "
            "electrostatic difference `d(U_H+U_ext) − U_proj_bg_impl` is **r-independent** — a flat "
            "constant across all projectile radii. This constant IS the quantum-vs-classical difference, "
            "and it is the plot the rest of the analysis builds on. Of it, ~1 eV is charged-cell gauge "
            "(from the section above) and the remainder is genuine WP electrostatic physics."),
        new_code_cell(
            "cc = pd.read_csv('rcut50_vs_120.csv')\n"
            "res = cc.completed_rc120.values; rvals = cc.r.values\n"
            "mean = res.mean(); spread = res.max()-res.min()\n"
            "gauge = 1.0                      # ~1 eV charged-cell gauge (from the three-method analysis)\n"
            "phys = mean - gauge\n"
            "fig, ax = plt.subplots(figsize=(8.0,4.8))\n"
            "ax.plot(rvals, res, 'o-', color='#1b6ca8', lw=2, ms=7, label=f'completed WP−CL residual (r_cut=120)')\n"
            "ax.axhline(mean, ls='--', color='#c0392b', lw=1, label=f'mean = {mean:.2f} eV (flat, spread {spread:.2f} eV)')\n"
            "ax.axhspan(mean-gauge, mean, color='0.75', alpha=0.5, label=f'~{gauge:.0f} eV charged-cell gauge')\n"
            "ax.axhspan(0, mean-gauge, color='#2ca02c', alpha=0.10, label=f'~{phys:.1f} eV genuine WP electrostatic physics')\n"
            "ax.set_ylim(0, mean*1.5)\n"
            "ax.set_xlabel('projectile radius r (Bohr from slab face)')\n"
            "ax.set_ylabel('WP − classical electrostatic residual (eV)')\n"
            "ax.set_title('MAIN RESULT: r-independent classical−quantum residual (r_cut=120)')\n"
            "ax.legend(frameon=False, fontsize=8, loc='upper right'); fig.tight_layout()\n"
            "fig.savefig('main_residual.png', dpi=140); plt.show()\n"
            "print(f'Completed WP−CL residual (r_cut=120): flat at {mean:.2f} eV across r={list(rvals)} (spread {spread:.2f} eV)')\n"
            "print(f'  decomposition:  ~{gauge:.0f} eV charged-cell gauge  +  ~{phys:.1f} eV genuine quantum electrostatic physics')"),
        new_markdown_cell("## Self-interaction test — does adding E_xc remove it?\n"
            "The 7.4 eV electrostatic residual contains the WP's self-Hartree (+21.7 eV free-space), which "
            "in exact DFT is cancelled by self-exchange. Test: add the WP's XC change "
            "`dE_xc = E_xc(WP) − E_xc(GS)` to the residual and see how much cancels. (E_xc(GS) is taken from "
            "the classical run, whose density = the bare GS.)"),
        new_code_cell(
            "cx = pd.read_csv('rcut50_vs_120.csv')\n"
            "print(cx[['r','completed_rc120','dE_xc','completed_plus_xc']].to_string(index=False))\n"
            "R=cx.completed_rc120; RX=cx.completed_plus_xc\n"
            "fig, ax = plt.subplots(figsize=(7.4,4.6))\n"
            "ax.axhline(0, color='0.6', lw=0.6)\n"
            "ax.plot(cx.r, R,  'o-', color='#1b6ca8', label=f'electrostatic residual R = dHE−impl  ({R.mean():.1f} eV)')\n"
            "ax.plot(cx.r, RX, 's-', color='#c0392b', label=f'R + dE_xc  ({RX.mean():.1f} eV, over-corrected)')\n"
            "ax.axhspan(-0.2, 0.2, color='0.85')\n"
            "ax.set_xlabel('r (Bohr)'); ax.set_ylabel('energy (eV)')\n"
            "ax.set_title('Adding full dE_xc over-corrects past zero (7.4 -> -9.1 eV)')\n"
            "ax.legend(frameon=False, fontsize=8); fig.tight_layout(); plt.show()\n"
            "print(f'\\ndE_xc = {cx.dE_xc.mean():.2f} eV (r-independent); R+dE_xc = {RX.mean():.2f} eV (flat).')\n"
            "print('Over-shoot past 0 => full dXC = self-exchange + WP-bath exchange, NOT pure self-interaction.')"),
        new_markdown_cell("**Interpretation.** `dE_xc = −16.5 eV` (r-independent) and `R + dE_xc = −9.1 eV` — the "
            "residual flips sign, i.e. the full XC channel *over-cancels*. The electrostatic residual R=7.4 eV "
            "is already `+21.7 (self-Hartree) − 14.3 (physical WP↔slab attraction)`, so it does not hold a bare "
            "self-interaction to cancel; subtracting the entire `dXC` (self-exchange **and** the physical "
            "WP↔bath exchange) removes too much. A *clean* SIE would be `WP self-Hartree + WP self-exchange` "
            "on the WP orbital alone (a few-eV residual in LDA). So E_xc removes most of the self-Hartree, but "
            "not as a simple subtraction from the electrostatic residual. *Verdict is yours.*"),
        new_markdown_cell("**Gauge verdict.** All three methods put the 2D-periodic charged-cell gauge at "
            "**order ~1 eV** (lateral ~0.1 eV, open-z ~1 eV, Makov–Payne ~0.77 eV). Against the ~7.4 eV "
            "completed residual, the gauge is a **small (~15%) contamination** — the residual is dominated by "
            "**genuine physics** (the WP electron's real electrostatic self-energy: 21.7 eV self-Hartree, "
            "largely cancelled by its attraction to the slab, plus its quantum spread vs a point ghost). "
            "*Caveat:* the B1-approximate p2 kernel (not INQ's exact Rozzi) is used for (1)/(2); three "
            "independent routes agreeing at ~1 eV makes the conclusion robust, but the exact number could "
            "move ~1 eV. The open-z gauge grows with Lz (a net charge in an aperiodic direction). *Verdict is yours.*"),

        # ── RESOLUTION (2026-07-13): the 7.4 eV is a pseudopotential representation artifact ──
        new_markdown_cell(
            "## ★★ RESOLVED — the 7.4 eV is a pseudopotential representation artifact, **not** a gauge\n\n"
            "The gauge speculation above is **superseded** by a direct verification in INQ's *own* p2 "
            "open-z convention (single-point eval `scripts/localised_jellium_dynamics/eval_projpot/run.cpp`; "
            "densities from the p5 at-rest run, r=12). The r-independent electrostatic residual is **not** a "
            "physical quantity to explain — it is the wavepacket's Hartree self-energy **minus a numerical "
            "artifact** in how the classical projectile's pseudopotential is represented on the grid.\n\n"
            "The classical projectile is an *external potential* with no self-energy in the ledger; the quantum "
            "WP is a real −1 charge cloud that repels itself. That self-repulsion is the whole residual — cleanly "
            "**+21.5 eV** (free-space) / **20.8 eV** (INQ p2). The reading of 7.4 eV is that self-Hartree minus a "
            "representation error δv = v_ion − V_proj_ideal (the erf(r/0.5)/r ≈ 1/r ghost tail wrapping and "
            "aliasing at r_cut=120) coupling to the slab electron spillout."),
        new_markdown_cell(
            "### Exact decomposition (reproduced to ±0.05 eV vs INQ)\n"
            "$$\\text{residual} \\;=\\; \\tfrac12 J[n_{WP},n_{WP}] \\;-\\; "
            "\\int (n_{\\text{slab}}-n_{+})\\,(v_{\\text{ion}}-V_{\\text{proj}}^{\\text{ideal}})\\,d^3r"
            "\\;=\\; E_{\\text{self-Hartree}} \\;-\\; \\text{(pseudopotential error}\\times\\text{spillout)}$$\n"
            "Terms: $\\tfrac12 J[n_{WP},n_{WP}]$ = WP Hartree self-energy; $n_{\\text{slab}}-n_{+}$ = slab net "
            "charge (electron spillout past the +background edges); $v_{\\text{ion}}$ = INQ's as-implemented "
            "projectile potential; $V_{\\text{proj}}^{\\text{ideal}}$ = Poisson potential of the true Gaussian charge. "
            "Provenance: `projpot_decomp.csv`, `projpot_wrap_findings.md`."),
        new_code_cell(
            "dc = pd.read_csv('projpot_decomp.csv')\n"
            "dc['resid_recon'] = dc.self_hartree_eV - dc.err_term_eV\n"
            "print(dc[['r_cut','ideal_eV','impl_eV','gap_eV','err_term_eV',"
            "'self_hartree_eV','resid_recon','residual_inq_eV']].to_string(index=False))\n"
            "print('\\nresidual = self_Hartree - INT(n_slab-n+).(v_ion-V_proj_ideal); reproduces INQ to +-0.05 eV.')\n"
            "print('  r_cut=50 : 20.81 - 6.82  = 14.0 eV  (INQ 13.99)')\n"
            "print('  r_cut=120: 20.81 - 13.45 = 7.4 eV   (INQ 7.36)')\n"
            "print('WP self-Hartree: 21.71 free-space / 21.49 periodic-FFT / 20.81 INQ p2  (gauge <1 eV).')\n"
            "print('The error term DOUBLES 6.8->13.5 eV as the erf/r tail wraps further at larger r_cut.')"),
        new_markdown_cell(
            "### The impl term is grid-pathological — the smoking gun\n"
            "Refining the grid (`dx`: 0.5→0.25 Bohr, no GS needed — the term depends only on the background + "
            "pseudopotential). A physical energy **converges** with grid; the r_cut=120 `impl` term **swings "
            "sign** (−524→−269→+224→+34 eV), while the `ideal` (true-Gaussian) term is grid-stable (~135 eV) and "
            "the r_cut=50 `impl` is stable (~−140 eV). Sign-changing non-convergence = numerical aliasing, not "
            "physics. Provenance: `projpot_gridsweep.csv`."),
        new_code_cell(
            "g = pd.read_csv('projpot_gridsweep.csv')\n"
            "print(g.to_string(index=False))\n"
            "fig, ax = plt.subplots(figsize=(7.4,4.6))\n"
            "for rc,mk,col in [(50,'o','#1b6ca8'),(120,'s','#c0392b')]:\n"
            "    s = g[g.r_cut==rc].sort_values('dx')\n"
            "    ax.plot(s.dx, s.impl_eV, mk+'-', color=col, label=f'impl  (-INT n+.v_ion),  r_cut={rc}')\n"
            "s0 = g[g.r_cut==50].sort_values('dx')\n"
            "ax.plot(s0.dx, s0.ideal_eV, '^--', color='#2e8b57', label='ideal (INT n_proj.v_bg) - grid-stable')\n"
            "ax.axhline(0, color='0.6', lw=0.6)\n"
            "ax.set_xlabel('grid spacing dx (Bohr)'); ax.set_ylabel('E_proj_bg (eV)')\n"
            "ax.set_title('r_cut=120 impl swings sign with dx (aliasing); r_cut=50 & ideal converge')\n"
            "ax.legend(frameon=False, fontsize=8); ax.invert_xaxis()\n"
            "fig.tight_layout(); fig.savefig('projpot_gridsweep.png', dpi=140); plt.show()"),
        new_markdown_cell(
            "**Takeaway.**\n"
            "- The r-independent electrostatic residual **is** the WP Hartree self-energy ≈ **20.8 eV** (INQ p2) / "
            "21.5 eV (free-space) — a real physical quantity a point ghost cannot carry.\n"
            "- The reading of **7.4 eV** is that self-Hartree minus a **pseudopotential representation artifact** "
            "(the long-range erf/r ghost tail aliasing at r_cut=120), proven three ways: exact decomposition, "
            "r_cut-doubling of the error term, and sign-changing grid non-convergence of `impl`.\n"
            "- The **charged-cell gauge is <1 eV** (self-Hartree channel only) — the ~1 eV figure above is right "
            "in magnitude but is *not* the origin of the 7.4 eV.\n"
            "- **r_cut = Lz is the worst choice, not the necessary one**: the ghost potential is long-range, so a "
            "bigger cutoff drags the tail across the box and aliases. Use the analytic `ideal` term, never a "
            "larger cutoff. See `reference_ghost_upf_tail_aliasing`."),
    ]
    nb = new_notebook(); nb.cells = cells
    nb.metadata.kernelspec = {"name":"python3","display_name":"Python 3"}
    nbf.write(nb, str(OUT/"ledger_rcut.ipynb"))
    print("wrote", OUT/"ledger_rcut.ipynb")
except Exception as e:
    print("notebook assembly skipped:", e)
