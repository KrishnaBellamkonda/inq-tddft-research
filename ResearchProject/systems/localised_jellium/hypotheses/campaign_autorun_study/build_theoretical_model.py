#!/usr/bin/env python3
"""build_theoretical_model.py — electrostatic (infinite sheet / slab) model of the
localised jellium, in SI, per the user's theoretical modelling.

Implements the three cases the user pasted, verbatim (φ=0 reference at the plane,
z_q measured from the slab centre):
  (1) single infinite sheet, areal σ at z=0:   φ = -σ/(2ε0)|z_q|,  U = qφ
  (2) uniform slab ρ0, thickness L, centred:    outside φ = -ρ0 L/(2ε0) z_q
      (plus the exact parabolic interior, shown for completeness)
  (3) slab collapsed to one sheet σ_tot = ρ0 L at z=0 (identical field outside)

The model gives a UNIFORM field outside → φ (and U) LINEAR in distance — i.e. the
interaction magnitude GROWS with distance. Graphs show that expectation. A final,
clearly-labelled EXTENSION superposes the positive-background sheet and the electron
sheet (the neutral slab) — within the same electrostatic model — so the two linear
terms can be compared with the KS-simulation curves. No verdict is drawn.

Run (venv): PYTHONPATH=.../inq-stack/python .../venv/bin/python3 build_theoretical_model.py
            python3 -m nbconvert --to notebook --execute --inplace theoretical_slab_model.ipynb
"""
from __future__ import annotations
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

REPO = "/local/data/public/skcb2/tddft"
LJ = f"{REPO}/ResearchProject/systems/localised_jellium"
OUT = Path(f"{LJ}/hypotheses/campaign_autorun_study")

CONST = """import numpy as np
import matplotlib.pyplot as plt   # kernel inline backend captures plt.show()
try:
    import sys; sys.path.insert(0, "%s/inq-stack/python")
    from inqview.visualisation import style as _st; _st.apply()
except Exception: pass""" % REPO + """

REPO_P = "/local/data/public/skcb2/tddft"

# ---- SI constants ---------------------------------------------------------
eps0 = 8.8541878128e-12     # F/m
e    = 1.602176634e-19      # C  (elementary charge, positive)
a0   = 5.29177210903e-11    # m  (Bohr)
J_per_eV = e                # 1 eV in J

# ---- slab parameters (SAME slab as the campaign_autorun runs) -------------
n0   = 1.312e-3             # e / Bohr^3  (positive-background number density)
half = 12.5                 # Bohr        (slab half-width)
L    = 2*half*a0            # m           (slab thickness = 25 Bohr)
rho0 = n0*e/a0**3           # C/m^3       (positive background CHARGE density, >0)
sigma_tot = rho0*L          # C/m^2       (collapsed areal density of the + sheet)
q_proj = -e                 # C           (projectile = one electron, negative)

print(f'L = {L*1e9:.3f} nm ({2*half:.0f} Bohr)')
print(f'rho0 = {rho0:.3e} C/m^3,  sigma_tot = rho0*L = {sigma_tot:.3f} C/m^2')

# axis: r = distance from the NEAR slab face (Bohr), matching the simulation.
# z_q (from centre) = r + half.  Model is valid outside the slab (r > 0).
r_bohr = np.linspace(0, 80, 400)
zq = (r_bohr + half)*a0     # m, distance from slab centre
def to_eV(U): return U/J_per_eV
def to_eVBohr(slope): return slope*a0/J_per_eV"""

SHEET = """# CASE 1 — single infinite sheet of areal density sigma at z=0.
#   phi(z_q) = -sigma/(2 eps0) |z_q|      U = q phi
# Uniform field E = sigma/(2 eps0): phi grows LINEARLY with |z_q| (no 1/r decay).
def phi_sheet(zq_m, sigma):   return -sigma/(2*eps0)*np.abs(zq_m)
def U_sheet(zq_m, sigma, q):  return q*phi_sheet(zq_m, sigma)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(r_bohr, phi_sheet(zq, sigma_tot), color='#8e44ad')
ax[0].set_title('phi(z_q): single + sheet (sigma = sigma_tot)')
ax[0].set_ylabel('phi (V)')
ax[1].plot(r_bohr, to_eV(U_sheet(zq, sigma_tot, q_proj)), color='#8e44ad')
ax[1].set_title('U(z_q) for q = -e (electron)  vs a + sheet')
ax[1].set_ylabel('U (eV)')
for a in ax: a.set_xlabel('r = z_q - L/2  (Bohr from near face)'); a.axhline(0, color='.6', lw=.7)
fig.suptitle('Infinite sheet: uniform field -> potential LINEAR in distance (magnitude grows)')
fig.tight_layout(); plt.show()
print('slope dU/d|z| =', f'{to_eVBohr(-sigma_tot*q_proj/(2*eps0)):.2f}', 'eV/Bohr (constant, does not decay)')"""

SLAB = """# CASE 2 & 3 — uniform slab (rho0, thickness L, centred) vs its collapsed sheet.
# Exact slab (phi=0 at centre):
#   inside  |z|<=L/2 :  phi = -rho0 z^2 /(2 eps0)          (parabolic)
#   outside |z| >L/2 :  phi = -rho0 L/(2 eps0) z + rho0 L^2/(8 eps0)
# User's slab-outside / collapsed sheet:  phi = -rho0 L/(2 eps0) z   (linear, no offset)
# -> identical FIELD (slope) outside; they differ by the constant rho0 L^2/(8 eps0).
def phi_slab_exact(z_m):
    z = np.abs(z_m); inside = z <= L/2
    out = -rho0*L/(2*eps0)*z + rho0*L**2/(8*eps0)
    ins = -rho0*z**2/(2*eps0)
    return np.where(inside, ins, out)
def phi_collapsed(z_m):  return -rho0*L/(2*eps0)*np.abs(z_m)

zc = np.linspace(-40, 80, 600)*a0
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(zc/a0, phi_slab_exact(zc), color='#c0392b', label='exact slab (parabolic inside)')
ax[0].plot(zc/a0, phi_collapsed(zc), '--', color='#1b6ca8', label='collapsed sheet (user formula)')
ax[0].axvspan(-half, half, color='.92', zorder=0)
ax[0].set_title('phi(z_q): slab vs collapsed sheet'); ax[0].set_ylabel('phi (V)')
ax[0].set_xlabel('z_q (Bohr from centre)'); ax[0].legend(frameon=False, fontsize=8)
# U for the projectile, outside only (r>0)
ax[1].plot(r_bohr, to_eV(q_proj*phi_slab_exact(zq)), color='#c0392b', label='exact slab')
ax[1].plot(r_bohr, to_eV(q_proj*phi_collapsed(zq)), '--', color='#1b6ca8', label='collapsed sheet')
ax[1].set_title('U(z_q) for q=-e outside the slab'); ax[1].set_ylabel('U (eV)')
ax[1].set_xlabel('r (Bohr from near face)'); ax[1].legend(frameon=False, fontsize=8)
for a in ax: a.axhline(0, color='.6', lw=.7)
fig.suptitle('Slab vs collapsed sheet: same slope outside (differ by constant rho0 L^2/(8 eps0))')
fig.tight_layout(); plt.show()
phi_off = rho0*L**2/(8*eps0)                       # potential offset (Volts)
print('constant potential offset rho0 L^2/(8 eps0) =', f'{phi_off:.1f}', 'V'
      '  -> energy offset q*offset =', f'{to_eV(q_proj*phi_off):.1f}', 'eV (constant, no r-dependence)')"""

EXPECT = """# WHAT THE MODEL SAYS WE SHOULD EXPECT — U vs distance over the simulation's range.
# For the projectile q=-e against the POSITIVE background sheet (sigma = +sigma_tot):
U_bg = to_eV(U_sheet(zq, +sigma_tot, q_proj))     # attraction: grows (more positive) with r
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.plot(r_bohr, U_bg, color='#8e44ad', lw=2,
        label='model: U(projectile vs + background sheet)')
ax.axhline(0, color='.6', lw=.7)
ax.set_xlabel('r (Bohr from near face)'); ax.set_ylabel('U (eV)')
ax.set_title('Model expectation: interaction LINEAR in distance (no decay)')
ax.legend(frameon=False)
fig.tight_layout(); plt.show()
print('MODEL: |U| grows ~linearly, slope', f'{to_eVBohr(sigma_tot*e/(2*eps0)):.2f}', 'eV/Bohr,'
      ' unbounded — it does NOT flatten or decay with distance.')"""

SUPERPOSE = """# EXTENSION (still the same electrostatic model) — the NEUTRAL slab = two sheets.
# The jellium slab is neutral: a POSITIVE background sheet (+sigma_tot) AND an
# electron sheet (-sigma_tot), both ~centred. Superpose their sheet potentials:
U_pos = to_eV(U_sheet(zq, +sigma_tot, q_proj))   # projectile vs + background
U_neg = to_eV(U_sheet(zq, -sigma_tot, q_proj))   # projectile vs - electron sheet
U_net = U_pos + U_neg                            # neutral slab (two coincident sheets)
fig, ax = plt.subplots(figsize=(7.6, 4.6))
ax.plot(r_bohr, U_pos, color='#8e44ad', label='+background sheet only (grows)')
ax.plot(r_bohr, U_neg, color='#c0392b', label='-electron sheet only (falls)')
ax.plot(r_bohr, U_net, color='k', lw=2.4, label='neutral slab = sum (two coincident sheets)')
ax.axhline(0, color='.6', lw=.7)
ax.set_xlabel('r (Bohr from near face)'); ax.set_ylabel('U (eV)')
ax.set_title('Superposition: two opposite coincident sheets cancel exactly outside')
ax.legend(frameon=False, fontsize=8); fig.tight_layout(); plt.show()
print('Two coincident opposite sheets: U_pos + U_neg =', f'{np.abs(U_net).max():.1e}', 'eV (exact cancellation).')
print('If the two sheets are instead SEPARATED (a dipole/finite slab), a residual remains —')
print('this is where the model connects to the finite, screened KS slab. Comparison is yours.')"""

VBG = """# NUMERICAL v_bg FROM THE POISSON SOLVER vs the infinite-plate prediction.
# v_bg = -poisson(n_+) dumped from INQ (the p2 Rozzi slab-truncated kernel that
# the localised_background_perturbation actually uses) via
# scripts/campaign_autorun/dump_vbg/run.cpp -- z-lineout at x=y=0. Atomic units:
# v_bg is the potential ENERGY an electron feels from the +background.
#   Infinite-plate prediction (atomic units, e=1): d^2 v_bg/dz^2 = +4*pi*n0 inside,
#   so v_bg is PARABOLIC inside (min at centre) and LINEAR outside with slope
#   4*pi*n0*a  (a = slab half-width); the FIELD E_z is a CONSTANT +/-4*pi*n0*a
#   outside -- the infinite-sheet hallmark (it does NOT decay to zero).
# The absolute offset of v_bg is a p2 G=0 gauge constant, so the potential panel
# references both curves to v_bg(0)=0; the field panel is gauge-free outright.
CSV = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
       "scripts/campaign_autorun/dump_vbg/results/vbg/vbg_lineout.csv")
HA = 27.211386
d = np.genfromtxt(CSV, delimiter=",", names=True)
z = d["z_bohr"]; vbg = d["v_bg_ha"]; nplus = d["n_plus"]
a_nom, n0v = 12.5, 1.312e-3
a_eff = z[nplus > 0].max()                 # discretised sharp-slab half-width
def vbg_plate(zz, a):                       # analytic infinite plate (atomic units), centred
    return np.where(np.abs(zz) <= a, 2*np.pi*n0v*zz**2,
                    2*np.pi*n0v*a**2 + 4*np.pi*n0v*a*(np.abs(zz)-a))
vbg_c  = (vbg - vbg[np.argmin(np.abs(z))]) * HA     # numeric, centred, eV
vplate = vbg_plate(z, a_eff) * HA                    # analytic, eV
Ez     = -np.gradient(vbg*HA, z)                     # field, eV/Bohr
Ez_out = 4*np.pi*n0v*a_eff*HA                        # analytic constant outside
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
ax[0].axvspan(-a_eff, a_eff, color="0.85", label="slab (+background)")
ax[0].plot(z, vbg_c, lw=2, label="INQ v_bg (p2 Poisson), centred")
ax[0].plot(z, vplate, "--", lw=1.6, label=f"infinite plate, slope 4πn₀a (a={a_eff:.2f})")
ax[0].set_xlabel("z (Bohr)"); ax[0].set_ylabel("v_bg − v_bg(0)  (eV)")
ax[0].set_title("Background well: numeric vs infinite plate")
ax[0].legend(frameon=False, fontsize=8)
ax[1].axvspan(-a_eff, a_eff, color="0.85")
ax[1].plot(z, Ez, lw=2, label="E_z = −dv_bg/dz (numeric)")
ax[1].axhline(+Ez_out, ls="--", color="k", lw=1, label=f"±4πn₀a = {Ez_out:.3f} eV/Bohr")
ax[1].axhline(-Ez_out, ls="--", color="k", lw=1); ax[1].axhline(0, color="0.6", lw=0.6)
ax[1].set_xlabel("z (Bohr)"); ax[1].set_ylabel("E_z (eV/Bohr)")
ax[1].set_title("Field: constant (≠0) outside = infinite-plate hallmark")
ax[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); plt.show()
m = (z > 20) & (z < 40)
slope_num = np.polyfit(z[m], (vbg*HA)[m], 1)[0]
print(f"a_nominal = {a_nom} Bohr,  a_effective (discretised sharp slab) = {a_eff:.2f} Bohr")
print(f"outside slope: numeric = {slope_num:.4f} eV/Bohr, analytic 4*pi*n0*a = {Ez_out:.4f} eV/Bohr "
      f"({100*slope_num/Ez_out:.1f}%)")
far = np.abs(Ez[(z > 45) & (z < 55)]).mean()
print(f"far-field |E_z| plateau (z=45..55) = {far:.4f} eV/Bohr  (NON-zero: a charged plate, as expected)")"""

EMPIRICAL = """# EMPIRICAL-DENSITY PLATE MODEL — use the REAL GS n(z) as a stack of infinite sheets.
# Instead of a sharp slab, take the planar-mean electron density n_-(z) measured in the
# GS, keep the analytic background n_+(z), and build phi(z_q) by summing infinite sheets:
#   phi(z_q) = -1/(2 eps0) integral rho(z') |z_q - z'| dz'     (1D Poisson, sheet stack)
from inqview import load_vti
GS_VTI = REPO_P + ("/ResearchProject/systems/localised_jellium/scripts/campaign_autorun/"
                   "runs/h2/gs_p2_lz120/results/density_gs_system/density_gs_system.vti")
v = load_vti(GS_VTI, expect_centered_axis='z')
ne_z = v.data.mean(axis=(0,1))               # planar-mean electron density n_-(z), e/Bohr^3
zb   = v.z                                    # Bohr
npz  = np.where(np.abs(zb) < half, n0, 0.0)   # background n_+(z), sharp slab
z_si = zb*a0; dz = z_si[1]-z_si[0]            # m
ne_si, np_si = ne_z/a0**3, npz/a0**3          # number density 1/m^3
def phi_stack(rho_charge, zq_m):              # rho_charge: C/m^3 array on z_si
    return np.array([-1/(2*eps0)*np.trapezoid(rho_charge*np.abs(zqi - z_si), z_si) for zqi in zq_m])
zq = (r_bohr + half)*a0
U_e  = to_eV(q_proj*phi_stack(-e*ne_si, zq))          # projectile vs electron sheet-stack
U_p  = to_eV(q_proj*phi_stack(+e*np_si, zq))          # projectile vs background sheet-stack
U_net= to_eV(q_proj*phi_stack(e*(np_si-ne_si), zq))   # vs NET (real, non-cancelling) density
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
ax[0].plot(zb, ne_z, color='#1b6ca8', label='n_-(z) electrons (measured)')
ax[0].plot(zb, npz, '--', color='#c0392b', label='n_+(z) background')
ax[0].axvspan(-half, half, color='.92', zorder=0); ax[0].set_xlim(-30, 30)
ax[0].set_xlabel('z (Bohr)'); ax[0].set_ylabel('planar-mean density (e/Bohr^3)')
ax[0].set_title('Empirical GS profile'); ax[0].legend(frameon=False, fontsize=8)
ax[1].plot(r_bohr, U_e, color='#1b6ca8', label='vs electrons only (grows)')
ax[1].plot(r_bohr, U_p, color='#c0392b', label='vs background only (grows, opposite)')
ax[1].plot(r_bohr, U_net, 'k', lw=2.4, label='vs NET real density (residual)')
ax[1].axhline(0, color='.6', lw=.7); ax[1].set_xlabel('r (Bohr from near face)')
ax[1].set_ylabel('U (eV)'); ax[1].set_title('U(r) from the empirical sheet stack')
ax[1].legend(frameon=False, fontsize=8)
fig.suptitle('Empirical-density plate model: real n(z) as a sheet stack')
fig.tight_layout(); plt.show()
print('NET residual U(r) span:', f'{U_net.min():.1f}..{U_net.max():.1f} eV',
      ' (surface-dipole residual — the two sheets do NOT cancel exactly)')

# NET-ONLY: the physically relevant curve for the (neutral) WP system — the
# projectile's Coulomb interaction with the REAL net density rho = e(n+ - n-).
fig, ax = plt.subplots(figsize=(7.4, 4.6))
ax.plot(r_bohr, U_net, 'k', lw=2.4)
ax.axhline(0, color='#c0392b', lw=.9, ls='--')
ax.fill_between(r_bohr, U_net, 0, color='.85')
ax.set_xlabel('r (Bohr from near face)'); ax.set_ylabel('U_net (eV)')
ax.set_title('Net real-density interaction only: ~0 at every r (neutral slab -> WP total flat)')
ax.set_ylim(-10, 10)
fig.tight_layout(); plt.show()
print(f'NET-only |U| <= {np.abs(U_net).max():.1f} eV across r=0..80 — the Coulomb attraction'
      ' from the real (non-cancelling) density is ~0 at ALL distances. This is what makes the'
      ' WP total energy flat vs r (validates the WP run).')"""

CUTOFF = """# PROJECTILE-CUTOFF TEST — the classical projectile UPF has z_valence=0 and r_max=50 Bohr.
# A POINT charge above an infinite sheet, with its Coulomb TRUNCATED at r_cut, sees only
# a disk of radius sqrt(r_cut^2 - dz^2); the 1D interaction kernel becomes
#   K(dz) = (r_cut - |dz|)  for |dz| < r_cut, else 0     (vs the uncut sheet's -|dz|)
# so U DECAYS to 0 as the slab passes beyond r_cut — no screening needed to get a decay.
rcut = 50.0*a0                                   # UPF radial cutoff (Bohr -> m)
def U_classical(zq_m, rc):
    out = []
    for zqi in zq_m:
        dzs = np.abs(zqi - z_si); m = dzs < rc
        out.append(e**2/(2*eps0)*np.trapezoid(np.where(m, ne_si*(rc - dzs), 0.0), z_si))
    return to_eV(np.array(out))                   # projectile(-e) vs electrons(-e): repulsive, +
rr = np.linspace(0, 90, 300); zqf = (rr + half)*a0
U_cut = U_classical(zqf, rcut)                    # with the true 50-Bohr cutoff
U_big = U_classical(zqf, 200*a0)                  # a far larger cutoff (what "no cutoff" trends to)
# overlay the ACTUAL simulated classical excess (extended sweep, Lz=200)
import glob as _g, csv as _csv, sys as _sys
from pathlib import Path as _Path
FAR = REPO_P + "/ResearchProject/systems/localised_jellium/scripts/campaign_autorun/runs/h0_p2_far"
_sys.path.insert(0, REPO_P + "/ResearchProject/systems/localised_jellium/scripts/campaign_autorun")
from analyse_phase import gs_energy
def _tot(run):
    f=next(iter(_g.glob(FAR+'/'+run+'/**/observables.csv', recursive=True)))
    rows=list(_csv.reader(open(f))); h,d=rows[0],rows[1]; return float(d[h.index('energy_total')])
EGSf = gs_energy(_Path(FAR)/'gs_p2_lz200/results'); HA=27.211386
r_sim = [4,12,20,28,36,44,52,60,68,76]
cl_sim = [ (_tot('cl_r'+str(r)+'_p2')-EGSf)*HA for r in r_sim ]
fig, ax = plt.subplots(figsize=(7.6, 4.8))
ax.plot(rr, U_cut, color='#c0392b', lw=2, label='point charge, cutoff r_cut=50 Bohr (UPF)')
ax.plot(rr, U_big, '--', color='.5', label='point charge, r_cut=200 Bohr (grows)')
ax.plot(r_sim, cl_sim, 'ks', ms=6, label='SIMULATED classical excess (KS)')
ax.axvline(50, ls=':', color='#c0392b'); ax.text(50.5, ax.get_ylim()[1]*0.6, 'UPF r_max', color='#c0392b', fontsize=8)
ax.axhline(0, color='.6', lw=.7); ax.set_xlabel('r (Bohr from near face)'); ax.set_ylabel('U / excess (eV)')
ax.set_title('Cutoff test: does the classical decay track the UPF r_max = 50 Bohr?')
ax.legend(frameon=False, fontsize=8); fig.tight_layout(); plt.show()
print(f'cutoff model U(r=4)={U_cut[0]:.0f} eV, ->0 near r~{50-half:.0f} Bohr; '
      f'simulated cl excess r=4 -> r=52: {cl_sim[0]:.0f} -> {cl_sim[6]:.1f} eV')"""

CUTOFFRUNS = """# EMPIRICAL CUTOFF SWEEP (KS runs) — 4 truncated projectile UPFs, r_cut in {10,20,30,40}
# Bohr. If the classical decay is set by the projectile potential's finite range, each
# dE_total(r) curve should drop to 0 near its OWN cutoff. Confirms/denies the artifact.
import glob as _g2, csv as _csv2, sys as _s2
from pathlib import Path as _P2
CT   = REPO_P + "/ResearchProject/systems/localised_jellium/scripts/campaign_autorun/runs/cutoff_test"
GSP2 = REPO_P + "/ResearchProject/systems/localised_jellium/scripts/campaign_autorun/runs/h2/gs_p2_lz120/results"
_s2.path.insert(0, REPO_P + "/ResearchProject/systems/localised_jellium/scripts/campaign_autorun")
from analyse_phase import gs_energy as _gse
EGS2 = _gse(_P2(GSP2)); HA2 = 27.211386
def _cl(rc, r):
    fs = _g2.glob(CT + ('/rc%d/cl_r%d_p2/**/observables.csv' % (rc, r)), recursive=True)
    if not fs: return None
    rows = list(_csv2.reader(open(fs[0]))); h,d = rows[0], rows[1]
    return {k: float(d[h.index('energy_'+k)]) for k in ('total','external','ion','ion_kinetic')}
CUTS = [10,20,30,40]; RAD = [2,4,8,12,16,20,24,28,32,36,40]
COLS = {10:'#1b6ca8',20:'#16a085',30:'#e67e22',40:'#c0392b'}
fig, ax = plt.subplots(figsize=(8.2, 5.0)); eion_max = 0.0; nrun = 0
for rc in CUTS:
    xs, dE = [], []
    for r in RAD:
        c = _cl(rc, r)
        if c is None: continue
        nrun += 1; xs.append(r); dE.append((c['total']-EGS2)*HA2)
        eion_max = max(eion_max, abs(c['ion']), abs(c['ion_kinetic']))
    if xs:
        ax.plot(xs, dE, 'o-', color=COLS[rc], label=f'UPF r_cut = {rc} Bohr')
        ax.axvline(rc, ls=':', color=COLS[rc], lw=1.2)
ax.axhline(0, color='.6', lw=.7); ax.set_xlabel('r (Bohr from near face)')
ax.set_ylabel('dE_total = E_tot(0) - E_GS  (eV)')
ax.set_title('Classical dE_total(r) vs projectile UPF cutoff (dotted = each cutoff radius)')
ax.legend(frameon=False); fig.tight_layout(); plt.show()
print(f'{nrun} cutoff runs read;  max |E_ion| = {eion_max:.1e} Ha across all of them')
print('E_ion (and E_ion_kinetic) ARE part of total() and ARE streamed — identically 0 here')
print('(z_valence=0 projectile; the jellium background is an external potential, not an Ewald ion).')"""

IMAGE = """# ANALYTICAL IMPROVEMENT — charge near a METAL surface: the image potential.
# Physical (unscreened-model-beating) expectation once the cutoff is removed: the jellium
# is a metal, so an external charge induces a screening image charge, giving
#   U_image(z) = - q^2 / (4 z)   (atomic units; z = distance from the image plane)
# which DECAYS as ~1/z (attractive), unlike the rigid sheet's linear growth.
z_im = np.linspace(2, 90, 300)            # Bohr, distance from image plane
U_image_eV = -(1.0)/(4.0*z_im)*27.211386  # q=1 (magnitude), atomic -> eV
fig, ax = plt.subplots(figsize=(7.2, 4.4))
ax.plot(z_im, U_image_eV, color='#16a085', lw=2, label='image potential  -q^2/(4z)  (metal)')
ax.axhline(0, color='.6', lw=.7); ax.set_xlabel('z (Bohr from surface)'); ax.set_ylabel('U (eV)')
ax.set_title('Analytical image-charge model (screened metal surface) — decays as ~1/z')
ax.legend(frameon=False); fig.tight_layout(); plt.show()
print('Image model: attractive, ~1/z decay. NOTE: in the CURRENT runs the classical decay is')
print('dominated by the UPF r_max=50 cutoff (above), so the image law is the EXPECTED physical')
print('behaviour once the projectile potential reaches far enough — a hypothesis to test, not a fit.')"""

ASSUME = """## Model assumptions (so you can decide what differs from the KS runs)

Stated as assumptions, **not** as an explanation — the comparison is yours to make.

| Model assumes | The KS simulation has |
|---|---|
| **Infinite** sheet/slab (uniform in x,y) | periodicity-2 slab (periodic in x,y ⇒ laterally infinite) — matches |
| **Rigid** charge (no response) | a **metallic** electron gas that **screens** the projectile |
| **Single** sheet (net charge) | a **neutral** slab (+background and −electrons together) |
| **Point** test charge q | a finite Gaussian projectile (σ_WP=0.5) / a real WP electron |
| φ, U grow **linearly**, unbounded | E_ext / U_H each move ~linearly near the slab, but the **total** is flat (WP) or decays (classical) |

The single-sheet model predicts a magnitude that **grows** with distance (uniform
field). The two-sheet (neutral) superposition above shows how, within the *same*
electrostatic model, two opposite sheets cancel — which is the bridge to the
neutral-slab behaviour seen in the runs. What this means for your interpretation
is yours to conclude."""

# ============================================================================
# PART II — classical vs WP from the KS runs: energy deconstruction + the
# screening / WP-potential test (new density-saving runs). Neutral throughout.
# ============================================================================

LEDGER = """# PART II - CLASSICAL vs WP, full energy deconstruction (KS insertion runs, periodicity 2).
# The h0_p2 runs insert either a WP electron or a classical Gaussian projectile at matched radii r,
# off the SAME Lz=120 open-z GS, streaming every energy component. Sum(parts)==total to ~1e-13,
# so the deconstruction is exact (checked below). E_ion = E_ion_kinetic = 0 in all of them.
import glob as _g, csv as _c, sys as _s
from pathlib import Path as _P
_CA = REPO_P + "/ResearchProject/systems/localised_jellium/scripts/campaign_autorun"
_s.path.insert(0, _CA); from analyse_phase import gs_energy as _gse
HA = 27.211386
EGS = _gse(_P(_CA + "/runs/h2/gs_p2_lz120/results"))
RADL = [4,12,20,28,36,40]; KE = ['total','kinetic','hartree','xc','external']
def _cmp(tag, r):
    f = _g.glob(_CA + f"/runs/h0_p2/{tag}_r{r}_p2/**/observables.csv", recursive=True)[0]
    rows = list(_c.reader(open(f))); h, d = rows[0], rows[1]
    return {k: float(d[h.index('energy_'+k)]) for k in KE}
W = {r: _cmp('wp', r) for r in RADL}; C = {r: _cmp('cl', r) for r in RADL}

# headline: dE_total(r) = E_tot(0) - E_GS, WP vs classical
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
ax[0].plot(RADL, [(W[r]['total']-EGS)*HA for r in RADL], 'o-', color='#1b6ca8', lw=2, label='WP electron')
ax[0].plot(RADL, [(C[r]['total']-EGS)*HA for r in RADL], 's--', color='#c0392b', lw=2, label='classical projectile')
ax[0].axhline(0, color='.6', lw=.7); ax[0].set_xlabel('r (Bohr from near face)')
ax[0].set_ylabel('dE_total = E_tot(0) - E_GS  (eV)')
ax[0].set_title('Insertion excess: WP flat (~self-energy) vs classical decaying'); ax[0].legend(frameon=False)
# WP - CL component differences (physical split): quantum self-energy vs electrostatic
dkin = [(W[r]['kinetic']-C[r]['kinetic'])*HA for r in RADL]
dxc  = [(W[r]['xc']-C[r]['xc'])*HA for r in RADL]
dele = [((W[r]['hartree']+W[r]['external'])-(C[r]['hartree']+C[r]['external']))*HA for r in RADL]
dtot = [(W[r]['total']-C[r]['total'])*HA for r in RADL]
ax[1].plot(RADL, dkin, 'o-', color='#16a085', label='d(kinetic) = WP zero-point (const)')
ax[1].plot(RADL, dxc,  'o-', color='#e67e22', label='d(xc) = WP self-XC (const)')
ax[1].plot(RADL, dele, 'o-', color='#8e44ad', label='d(Hartree+external) = electrostatic')
ax[1].plot(RADL, dtot, 'k-', lw=2.4, label='d(total) = WP - CL')
ax[1].axhline(0, color='.6', lw=.7); ax[1].set_xlabel('r (Bohr from near face)')
ax[1].set_ylabel('WP - classical  (eV)'); ax[1].set_title('WP - CL by component (Hartree+external summed: G=0-robust)')
ax[1].legend(frameon=False, fontsize=8)
fig.suptitle('Classical vs WP energy deconstruction (periodicity 2)'); fig.tight_layout(); plt.show()

# arithmetic ledger (eV), rounded
print(f"E_GS = {EGS*HA:.1f} eV   (WP zero-point 3/(4*0.5^2) = {3/(4*0.25)*HA:.1f} eV)")
print(f"{'r':>4}{'dE_WP':>8}{'dE_CL':>8}{'WP-CL':>8}{'  ||':>5}{'dKin':>8}{'dXC':>8}{'d(H+E)':>9}")
for r in RADL:
    print(f"{r:>4}{(W[r]['total']-EGS)*HA:>8.1f}{(C[r]['total']-EGS)*HA:>8.1f}"
          f"{(W[r]['total']-C[r]['total'])*HA:>8.1f}{'  ||':>5}"
          f"{(W[r]['kinetic']-C[r]['kinetic'])*HA:>8.1f}{(W[r]['xc']-C[r]['xc'])*HA:>8.1f}"
          f"{((W[r]['hartree']+W[r]['external'])-(C[r]['hartree']+C[r]['external']))*HA:>9.1f}")
# exactness check
_bad = max(abs(sum(W[r][k] for k in KE if k!='total') - W[r]['total']) for r in RADL)
print(f"max |sum(parts) - total| = {_bad:.1e} Ha (deconstruction exact)")
print("dKin ~ +82 eV and dXC ~ -16 eV are r-INDEPENDENT (the WP's quantum self-energy: zero-point")
print("KE + LDA self-XC of the added electron); d(Hartree+external) carries all r-dependence.")"""

XCDIFF = """# THE E_xc DIFFERENCE and the holistic total difference (Learning #1).
# E_xc(WP) - E_xc(classical): the classical projectile is a bare external potential (no electron ->
# no exchange-correlation); the WP is a REAL electron, so it adds LDA XC. The difference isolates
# that. d(total) = [WP quantum self-energy (const)] + [electrostatic difference (r-dependent)].
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
ax[0].plot(RADL, dxc, 'o-', color='#e67e22', lw=2)
ax[0].axhline(np.mean(dxc), color='.5', ls='--', lw=.9)
ax[0].set_xlabel('r (Bohr from near face)'); ax[0].set_ylabel('E_xc(WP) - E_xc(CL)  (eV)')
ax[0].set_title(f'XC difference ~ {np.mean(dxc):.0f} eV, r-independent')
# total difference decomposed into the constant quantum part + the electrostatic part
quantum = [dkin[i]+dxc[i] for i in range(len(RADL))]
ax[1].plot(RADL, dtot, 'k-', lw=2.4, label='WP - CL total')
ax[1].plot(RADL, quantum, '--', color='#16a085', label='quantum self-energy (dKin+dXC, const)')
ax[1].plot(RADL, dele, '--', color='#8e44ad', label='electrostatic d(H+E)')
ax[1].axhline(0, color='.6', lw=.7); ax[1].set_xlabel('r (Bohr from near face)')
ax[1].set_ylabel('WP - CL  (eV)'); ax[1].set_title('Total difference = constant quantum + r-dependent electrostatic')
ax[1].legend(frameon=False, fontsize=8)
fig.suptitle('E_xc difference and the WP-classical total-energy difference'); fig.tight_layout(); plt.show()
print(f'E_xc(WP)-E_xc(CL) = {np.mean(dxc):.1f} +/- {np.std(dxc):.2f} eV (essentially r-independent ->')
print(' at these separations it is the WP self-XC, NOT a slab-proximity screening signature).')
print(f'Electrostatic d(H+E): {dele[0]:.0f} eV (r=4)  ->  {dele[-1]:.0f} eV (r=40): the WP and classical')
print(' converge electrostatically far away; the residual d(total) -> dKin+dXC ~ '
      f'{quantum[-1]:.0f} eV is the pure cost of quantising the projectile.')
# comparison to the analytical model (Part I): the neutral-slab NET electrostatic ~ 0.
print('MODEL LINK: the empirical sheet-stack (Part I) gives NET projectile-slab electrostatics ~0')
print(' (neutral slab). The WP insertion excess is FLAT (~80 eV self-energy) -> its electrostatic')
print(' interaction with the slab is ~0 at all r, consistent with that cancellation. The classical')
print(' excess instead DECAYS (185 -> 12 eV) = the UPF-cutoff term (Part I cutoff sweep). Yours to weigh.')"""

WPPOT = """# *** SCREENING / WP-POTENTIAL TEST (Learning #2) *** -- the centrepiece.
# Does the wavepacket produce the SAME Coulomb potential as the classical Gaussian projectile?
# Potential equivalence reduces to SOURCE equivalence: if n_WP == the classical Gaussian charge
# (std sigma_rho = sigma_WP/sqrt2 = 0.354, since |psi|^2 halves the variance) then the two
# projectiles are electrostatically identical under ANY solver. The only distortion is the WP's
# orthogonalisation against occupied slab states. New runs saved n_WP = |psi_WP|^2 at t=0.
from inqview import load_vti as _lv
from scipy.special import erf as _erf
SCR = _CA + "/runs/screening_wp"
s_rho = 0.5/np.sqrt(2)                       # sigma_rho (Bohr); sigma_WP = 0.5 (label)
def _radial(vti, z0):
    v = _lv(vti); n = v.data
    ix = int(np.argmin(np.abs(v.x))); iy = int(np.argmin(np.abs(v.y)))
    iz = int(np.argmin(np.abs(v.z - z0)))
    X,Y,Z = np.meshgrid(v.x-v.x[ix], v.y-v.y[iy], v.z-v.z[iz], indexing='ij')
    R = np.sqrt(X**2+Y**2+Z**2); dV = (v.x[1]-v.x[0])*(v.y[1]-v.y[0])*(v.z[1]-v.z[0])
    return v, n, R, dV, (ix,iy,iz)
def _poisson_box(v, n, ctr, halfw=8.0):       # FFT-Poisson on a cube around the WP (validated cell)
    ix,iy,iz = ctr; dx=v.x[1]-v.x[0]; hw=int(halfw/dx)
    sub = n[ix-hw:ix+hw, iy-hw:iy+hw, iz-hw:iz+hw]
    m = sub.shape[0]; k = 2*np.pi*np.fft.fftfreq(m, d=dx)
    KX,KY,KZ = np.meshgrid(k,k,k,indexing='ij'); K2 = KX**2+KY**2+KZ**2; K2[0,0,0]=1
    sub_f = np.fft.ifftshift(sub)
    V = np.fft.fftshift(np.fft.ifftn(4*np.pi*np.fft.fftn(sub_f)/K2).real)
    c = m//2; ax_r = (np.arange(c, m)-c)*dx
    return ax_r, V[c, c, c:] - V.mean()

fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
Q = {}
for r, z0, col in [(12, -24.5, '#1b6ca8'), (4, -16.5, '#16a085')]:
    vti = SCR + f"/wp_r{r}_p2/results/wp_r{r}_p2/density_wp/density_wp.vti"
    v, n, R, dV, ctr = _radial(vti, z0)
    Q[r] = n.sum()*dV
    # radial scatter of n_WP vs the ideal Gaussian charge (both integrate to 1)
    sel = R < 3.0
    ax[0].plot(R[sel].ravel()[::37], n[sel].ravel()[::37], '.', ms=2, color=col, alpha=.35, label=f'n_WP r={r} (sigma=0.5)')
    rr, Vwp = _poisson_box(v, n, ctr)
    m = (rr > 0.3) & (rr < 6); Vana = _erf(rr/(np.sqrt(2)*s_rho))/np.maximum(rr, 1e-9)
    off = np.mean(Vwp[m] - Vana[m])
    ax[1].plot(rr[rr>0], (Vwp-off)[rr>0], 'o', ms=3, color=col, label=f'poisson(n_WP) r={r}')
rg = np.linspace(0.02, 3.0, 200)
ng = (2*np.pi*s_rho**2)**-1.5*np.exp(-rg**2/(2*s_rho**2))
ax[0].plot(rg, ng, 'k-', lw=1.6, label='ideal Gaussian charge (sigma_rho=0.354)')
ax[0].set_yscale('log'); ax[0].set_ylim(1e-4, 3); ax[0].set_xlabel('distance from WP centre (Bohr)')
ax[0].set_ylabel('n_WP (e/Bohr^3)'); ax[0].set_title('WP source charge vs ideal Gaussian'); ax[0].legend(frameon=False, fontsize=8)
rg2 = np.linspace(0.3, 6, 200); Vg = _erf(rg2/(np.sqrt(2)*s_rho))/rg2
ax[1].plot(rg2, Vg, 'k-', lw=1.6, label='classical projectile potential  erf(r/(sqrt2 sigma_rho))/r')
ax[1].set_xlabel('distance from centre (Bohr)'); ax[1].set_ylabel('V (Ha)')
ax[1].set_title('WP Coulomb potential vs classical Gaussian'); ax[1].legend(frameon=False, fontsize=8)
fig.suptitle('Screening / WP-potential test: is the WP the same Coulomb source as the classical projectile?')
fig.tight_layout(); plt.show()
print(f'int n_WP = {Q.get(12,float("nan")):.4f} (r=12), {Q.get(4,float("nan")):.4f} (r=4)  [1 electron]')
print('The WP source charge tracks the ideal Gaussian(sigma_rho=0.354) at the core; poisson(n_WP)')
print('overlays the analytic classical-projectile potential erf(r/(sqrt2 sigma_rho))/r. Small tail')
print('deviations come from orthogonalising the WP against occupied slab states. The classical UPF')
print('potential was separately verified == this Gaussian potential to RMS 0.000 Ha. So at t=0 the two')
print('projectiles are (nearly) the same electrostatic source -- any energy gap is quantum, not a')
print('different potential. Whether the residual distortion matters is yours to judge.')"""

BATH = """# SCREENING BASELINE at t=0, and why dynamical screening is out of CPU reach here.
# The slab density at t=0 is density_total (the 82 slab electrons; the WP is a separate saved
# orbital |psi_WP|^2). At insertion the slab orbitals are STILL the GS orbitals (the WP was
# orthogonalised against them but nothing has evolved), so n_slab(t=0) - n_GS must be ~0. It is
# EXACTLY 0 here (bit-identical). Screening -- the slab polarising around the WP -- is a DYNAMICAL
# response that builds over ~ a plasmon period; that is far beyond these 1-step CPU runs.
from inqview import load_vti as _lv2
GS = _CA + "/runs/h2/gs_p2_lz120/results/density_gs_system/density_gs_system.vti"
vg = _lv2(GS)
vt = _lv2(SCR + "/wp_r12_p2/results/wp_r12_p2/density_total/density_total.vti")   # slab @ t=0
vw = _lv2(SCR + "/wp_r12_p2/results/wp_r12_p2/density_wp/density_wp.vti")         # |psi_WP|^2
dV = (vg.x[1]-vg.x[0])*(vg.y[1]-vg.y[0])*(vg.z[1]-vg.z[0])
dn = vt.data - vg.data                                   # slab response (should be exactly 0)
# planar-mean profiles along z: GS slab, slab-with-WP-present (overlaid -> identical), and the WP blob
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
ax[0].plot(vg.z, vg.data.mean(axis=(0,1)), color='#c0392b', lw=2.4, label='n_slab (GS)')
ax[0].plot(vt.z, vt.data.mean(axis=(0,1)), '--', color='#1b6ca8', lw=1.4, label='n_slab (WP present, t=0)')
ax[0].plot(vw.z, vw.data.mean(axis=(0,1)), color='#16a085', lw=1.2, label='n_WP (separate orbital)')
ax[0].axvspan(-half, half, color='.92', zorder=0); ax[0].axvline(-24.5, color='#16a085', ls=':')
ax[0].set_xlabel('z (Bohr)'); ax[0].set_ylabel('planar-mean density (e/Bohr^3)'); ax[0].set_xlim(-40, 20)
ax[0].set_title('Slab unchanged by WP insertion (curves coincide); WP is a separate blob')
ax[0].legend(frameon=False, fontsize=8)
ax[1].plot(vg.z, dn.mean(axis=(0,1)), color='#8e44ad')
ax[1].axhline(0, color='.6', lw=.7); ax[1].set_xlabel('z (Bohr)')
ax[1].set_ylabel('planar-mean [n_slab(t=0) - n_GS]  (e/Bohr^3)')
ax[1].set_title(f'Screening baseline = 0 exactly (max |d| = {np.abs(dn).max():.0e})'); ax[1].set_xlim(-40, 20)
fig.suptitle('Screening baseline at t=0: no instantaneous slab response'); fig.tight_layout(); plt.show()
n0v = 82/(50*50*25); wp_pl = np.sqrt(4*np.pi*n0v); Tp = 2*np.pi/wp_pl
print(f'max |n_slab(t=0) - n_GS| = {np.abs(dn).max():.1e} e/Bohr^3 (bit-identical: slab is the GS slab).')
print(f'omega_p = sqrt(4 pi n0) = {wp_pl:.3f} Ha  ->  plasmon period T_p = {Tp:.0f} a.u. = {Tp/0.01:.0f} steps at dt=0.01.')
print('NOTE: inqkit density::total returns the 82 slab electrons here; the WP is captured separately')
print('as the saved orbital density (int n_WP = 1). Static screening of a LOCALISED charge would need')
print('re-converging the SCF with a FIXED external Gaussian (a ground-state run); dynamical screening')
print('needs ~T_p of propagation (~5000 CPU steps). Both are GPU follow-ups -- recorded, not fabricated.')"""

PART2_NOTE = """# Part II - classical vs WP from the KS runs

New work: the full **energy deconstruction** of matched classical-projectile vs WP-electron
insertion runs, and the **screening / WP-potential test** (new density-saving runs).
All neutral - the reader draws the verdict.

- **Deconstruction (exact).** Every energy component streamed; sum(parts) == total to ~1e-13;
  E_ion = 0. The WP - classical difference splits cleanly into an r-**independent** quantum
  self-energy (zero-point KE ~ +82 eV, self-XC ~ -16 eV) and an r-**dependent** electrostatic
  term (Hartree+external, summed to be robust to the periodicity-2 G=0 offset).
- **E_xc difference** ~ -16 eV, r-independent - the WP's own LDA self-XC, not a slab-proximity
  signature at these separations.
- **WP-potential test.** n_WP = |psi_WP|^2 saved at t=0; it tracks the ideal Gaussian charge
  (sigma_rho = sigma_WP/sqrt2 = 0.354) at the core, and poisson(n_WP) overlays the analytic
  classical-projectile potential - so the two projectiles are (nearly) the same electrostatic source.
- **Screening baseline** n_slab(t=0) - n_GS = 0 *exactly* (bit-identical); genuine screening is
  dynamical (~T_plasmon, ~5000 CPU steps) and is flagged as a GPU follow-up, not fabricated.
"""

# ============================================================================
# PART III — single-run component decomposition at PERIODICITY 3 (offset-free
# G=0), and the quantum self-energy ledger. Full-component p3 mirror of h0_p2.
# ============================================================================

PART3_NOTE = """# Part III - component decomposition & the quantum self-energy (periodicity 3)

Part II summed Hartree+external because **periodicity-2** inflates each of them by the
Poisson `G=0` term `0.5*rc^2` (poisson.hpp:49) - a large, box-dependent, opposite-sign
constant on a net-charged cell, so only their SUM is physical. To read each component
individually we remade the SAME radius sweep at **periodicity 3** (`G=0 -> 0`,
poisson.hpp:31), streaming every energy component: `runs/h0_p3/{wp,cl}_r{4..40}_p3`,
off the p3 GS `shared_gs/slab_n82_L50x50x120` (dispatcher `rerun_h0_p3.py`). Both
projectiles share the same GS/box/grid; sum(parts)==total to ~1e-13.

**What is clean, and what is not (ground-up):**
- **Kinetic** difference `= +3.00 Ha` at every r - exactly the WP zero-point `3/(4 sigma^2)`.
  The classical ghost has no kinetic energy; the WP is a real particle. Convention-free.
- **XC** difference `= -0.61 Ha` at every r - the WP's LDA self-XC. Convention-free.
- **The projectile self-Hartree** (its Coulomb self-energy) is physical and equals
  `1/(2 sigma_rho sqrt(pi)) = 0.80 Ha (22 eV)` in closed form, matched numerically by
  FFT-Poisson on the saved `n_WP` (`0.77 Ha`). BUT it does **not** appear as the raw
  `E_hartree(WP) - E_hartree(classical)`: inserting the WP makes the cell net -1 charged,
  so the charged-cell Poisson convention injects a `~N_slab * mean(V[n_WP])` term into the
  Hartree/external split. Evidence: raw `dHartree(r=40)` is `-1.1 Ha` at p3 but `-10.1 Ha`
  at p2 - same physics, different convention - while the physical self-energy is `+0.80 Ha`.
  So the self-energy is recovered from `n_WP` directly, not from a component subtraction.

Headline (yours to weigh): the classical ghost is missing exactly the WP's **quantum
self-energy** - zero-point KE `+3.0 Ha` (exact) and self-Hartree `+0.8 Ha`
(analytic == numeric) - net of its self-XC `-0.6 Ha`.
"""

P3DATA = """# Load the periodicity-3 full-component mirror (h0_p3). E_nonlocal=E_ion=E_xx=0 (LDA jellium).
import glob as _g3, csv as _c3
_CA3 = REPO_P + "/ResearchProject/systems/localised_jellium/scripts/campaign_autorun"
HA = 27.211386
RAD3 = [4,12,20,28,36,40]; COMP = ['total','kinetic','hartree','xc','external','nonlocal','ion','ion_kinetic','exact_exchange']
def _row3(tag, r):
    f = _g3.glob(_CA3 + f"/runs/h0_p3/{tag}_r{r}_p3/**/observables.csv", recursive=True)[0]
    rows = list(_c3.reader(open(f))); h, d = rows[0], rows[1]
    return {k: float(d[h.index('energy_'+k)]) for k in COMP}
W3 = {r: _row3('wp', r) for r in RAD3}; C3 = {r: _row3('cl', r) for r in RAD3}
R0 = 28                                              # the single representative run for the panels
# exactness of the decomposition
_parts = lambda d: d['kinetic']+d['hartree']+d['xc']+d['external']+d['nonlocal']+d['ion']+d['ion_kinetic']+d['exact_exchange']
_bad3 = max(max(abs(_parts(W3[r])-W3[r]['total']), abs(_parts(C3[r])-C3[r]['total'])) for r in RAD3)
print(f"periodicity-3 mirror loaded: {len(RAD3)} radii x (wp,cl); sum(parts)-total max = {_bad3:.1e} Ha (exact)")
print(f"single-run panels use r = {R0} Bohr")"""

P3TOTAL = """# PANEL 1 - the difference in TOTAL energy magnitude, WP vs classical.
dtot3 = {r: (W3[r]['total']-C3[r]['total'])*HA for r in RAD3}
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
# (a) the two totals at r=R0 as bars, with the difference annotated
b = ax[0].bar(['WP electron','classical ghost'], [W3[R0]['total']*HA, C3[R0]['total']*HA],
              color=['#1b6ca8','#c0392b'], width=.6)
ax[0].set_ylabel('E_total  (eV)'); ax[0].set_title(f'Total energy at r = {R0} Bohr')
lo = min(W3[R0]['total'], C3[R0]['total'])*HA; ax[0].set_ylim(lo-60, 0)
for rect,val in zip(b,[W3[R0]['total']*HA, C3[R0]['total']*HA]):
    ax[0].text(rect.get_x()+rect.get_width()/2, val, f'{val:.0f}', ha='center', va='top', fontsize=9)
ax[0].annotate(f'DE = E_WP - E_CL = +{dtot3[R0]:.1f} eV', xy=(0.5, 0.12), xycoords='axes fraction',
               ha='center', fontsize=10, bbox=dict(boxstyle='round', fc='#fff3cd', ec='.6'))
# (b) that difference across r (context for the single run)
ax[1].plot(RAD3, [dtot3[r] for r in RAD3], 'ko-', lw=2)
ax[1].axvline(R0, color='.6', ls=':'); ax[1].axhline(0, color='.6', lw=.7)
ax[1].set_xlabel('r (Bohr from near face)'); ax[1].set_ylabel('E_total(WP) - E_total(CL)  (eV)')
ax[1].set_title('Total-energy difference vs r (r-dependent electrostatics)')
fig.suptitle('Panel 1 - magnitude of the WP-vs-classical total-energy difference')
fig.tight_layout(); plt.show()
print(f'At r={R0}: E_tot(WP) = {W3[R0]["total"]*HA:.0f} eV, E_tot(CL) = {C3[R0]["total"]*HA:.0f} eV,'
      f' DE = +{dtot3[R0]:.1f} eV.')
print(f'DE(r) ranges {dtot3[RAD3[0]]:.0f} eV (r=4) -> +{dtot3[RAD3[-1]]:.0f} eV (r=40): the sign change is')
print('the r-dependent electrostatic part; the constant floor is the quantum self-energy (Panel 4).')"""

P3WATER = """# PANEL 2 - individual component decomposition of EACH run (r=R0), periodicity 3.
# p3 gives offset-free, physically-signed components: E_hartree>0 (e-e repulsion),
# E_external<0 (electrons in the +background well). Each waterfall sums to its E_total.
order = ['kinetic','external','hartree','xc']                 # nonlocal/ion/xx are identically 0
labels = ['kinetic','external','hartree','xc','= total']
def waterfall(axw, d, title, col):
    run = 0.0; xs = range(len(order)+1)
    for i,k in enumerate(order):
        val = d[k]*HA
        axw.bar(i, val, bottom=run, color=col, edgecolor='k', lw=.5, width=.7)
        axw.text(i, run+val+(6 if val>=0 else -6), f'{val:.0f}', ha='center',
                 va='bottom' if val>=0 else 'top', fontsize=8)
        run += val
    axw.bar(len(order), run, color='.35', edgecolor='k', lw=.5, width=.7)
    axw.text(len(order), run-6, f'{run:.0f}', ha='center', va='top', fontsize=8, color='w')
    axw.axhline(0, color='.6', lw=.7); axw.set_xticks(list(xs)); axw.set_xticklabels(labels, fontsize=8)
    axw.set_ylabel('energy (eV)'); axw.set_title(title)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
waterfall(ax[0], W3[R0], f'WP electron  (r={R0})', '#1b6ca8')
waterfall(ax[1], C3[R0], f'classical ghost  (r={R0})', '#c0392b')
fig.suptitle('Panel 2 - per-run energy decomposition (periodicity 3, offset-free G=0); bars sum to total')
fig.tight_layout(); plt.show()
for tag,d in [('WP',W3[R0]),('CL',C3[R0])]:
    print(f'{tag} r={R0} (eV): ' + '  '.join(f'{k}={d[k]*HA:.0f}' for k in order) + f'  total={d["total"]*HA:.0f}')
print('Classical kinetic/hartree/xc are the FROZEN GS values (the ghost is a pure external potential,')
print('so at step 0 the electron density is the GS density); ALL its r-dependence is in E_external.')"""

P3DIFF = """# PANEL 3 - component-wise difference WP - classical vs r (periodicity 3).
dkin3 = [(W3[r]['kinetic']-C3[r]['kinetic'])*HA for r in RAD3]
dxc3  = [(W3[r]['xc']-C3[r]['xc'])*HA for r in RAD3]
dH3   = [(W3[r]['hartree']-C3[r]['hartree'])*HA for r in RAD3]
dext3 = [(W3[r]['external']-C3[r]['external'])*HA for r in RAD3]
dtot3l= [(W3[r]['total']-C3[r]['total'])*HA for r in RAD3]
zp = 3/(4*0.25)*HA                                            # zero-point 3/(4 sigma^2), sigma=0.5
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
ax[0].plot(RAD3, dkin3, 'o-', color='#16a085', label='d(kinetic)')
ax[0].axhline(zp, color='#16a085', ls='--', lw=1, label=f'zero-point 3/(4s^2) = {zp:.0f} eV')
ax[0].plot(RAD3, dxc3, 's-', color='#e67e22', label='d(xc) = self-XC')
ax[0].axhline(0, color='.6', lw=.7); ax[0].set_xlabel('r (Bohr from near face)')
ax[0].set_ylabel('WP - classical  (eV)'); ax[0].set_title('Convention-free differences: kinetic & xc (flat)')
ax[0].legend(frameon=False, fontsize=8)
ax[1].plot(RAD3, dH3,   'o-', color='#8e44ad', label='d(hartree)  [charged-cell convention]')
ax[1].plot(RAD3, dext3, 'o-', color='#c0392b', label='d(external) [charged-cell convention]')
ax[1].plot(RAD3, [dH3[i]+dext3[i] for i in range(len(RAD3))], 'k-', lw=2, label='d(hartree+external)')
ax[1].plot(RAD3, dtot3l, '--', color='.4', lw=1.4, label='d(total)')
ax[1].axhline(0, color='.6', lw=.7); ax[1].set_xlabel('r (Bohr from near face)')
ax[1].set_ylabel('WP - classical  (eV)'); ax[1].set_title('Electrostatic components (net-charge convention-dependent)')
ax[1].legend(frameon=False, fontsize=8)
fig.suptitle('Panel 3 - component-wise WP-classical difference vs r'); fig.tight_layout(); plt.show()
print(f'd(kinetic) = {np.mean(dkin3):.1f} +/- {np.std(dkin3):.2f} eV (flat) == zero-point 3/(4s^2) = {zp:.1f} eV.')
print(f'd(xc)      = {np.mean(dxc3):.1f} +/- {np.std(dxc3):.2f} eV (flat) = WP self-XC.')
print(f'd(hartree) swings {dH3[0]:.0f} -> {dH3[-1]:.0f} eV and d(external) {dext3[0]:.0f} -> {dext3[-1]:.0f} eV:')
print(' individually r-dependent AND net-charge-convention-dependent (see Panel 4) - not directly the self-energy.')"""

P3SELF = """# PANEL 4 - the quantum self-energy ledger: what the classical ghost is MISSING.
# Three pieces the ghost lacks, each cross-checked against an independent calculation.
from inqview import load_vti as _lv4
s_rho4 = 0.5/np.sqrt(2)                                       # sigma_rho = sigma_WP/sqrt2
E_self_ana = 1.0/(2*s_rho4*np.sqrt(np.pi))                    # closed-form Gaussian self-Hartree (Ha)
# numeric self-Hartree: isolated FFT-Poisson on the saved n_WP (padded cube -> non-periodic)
_vti = _CA3 + "/runs/screening_wp/wp_r12_p2/results/wp_r12_p2/density_wp/density_wp.vti"
_v = _lv4(_vti); _n = _v.data; _dx=_v.x[1]-_v.x[0]; _dy=_v.y[1]-_v.y[0]; _dz=_v.z[1]-_v.z[0]
_i0 = np.unravel_index(np.argmax(_n), _n.shape); _hw = 8.0
_s = tuple(slice(max(0,_i0[a]-int(_hw/d)), _i0[a]+int(_hw/d)) for a,d in enumerate((_dx,_dy,_dz)))
_sub = np.pad(_n[_s], [(s,s) for s in _n[_s].shape])         # zero-pad x2 for isolation
_nx,_ny,_nz = _sub.shape
_kx=2*np.pi*np.fft.fftfreq(_nx,_dx);_ky=2*np.pi*np.fft.fftfreq(_ny,_dy);_kz=2*np.pi*np.fft.fftfreq(_nz,_dz)
_KX,_KY,_KZ=np.meshgrid(_kx,_ky,_kz,indexing='ij'); _K2=_KX**2+_KY**2+_KZ**2; _K2[0,0,0]=1
_Vk=4*np.pi*np.fft.fftn(_sub)/_K2; _Vk[0,0,0]=0
E_self_num = 0.5*np.sum(_sub*np.real(np.fft.ifftn(_Vk)))*_dx*_dy*_dz
# measured convention-free pieces (mean over r, they are flat)
dkin_meas = np.mean([(W3[r]['kinetic']-C3[r]['kinetic'])*HA for r in RAD3])/HA   # Ha
dxc_meas  = np.mean([(W3[r]['xc']-C3[r]['xc'])*HA for r in RAD3])/HA             # Ha
zp_ana = 3/(4*0.25)                                          # Ha
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
# (a) grouped bars: measured vs analytic for each self-energy piece
comps = ['zero-point KE','self-Hartree','self-XC']
meas  = [dkin_meas*HA, np.nan, dxc_meas*HA]                  # self-Hartree not cleanly measurable (see text)
ana   = [zp_ana*HA, E_self_ana*HA, np.nan]
num   = [np.nan, E_self_num*HA, np.nan]
x = np.arange(3); w = .27
ax[0].bar(x-w, meas, w, color='#1b6ca8', label='measured (dcomponent)')
ax[0].bar(x,   ana,  w, color='#16a085', label='analytic')
ax[0].bar(x+w, num,  w, color='#e67e22', label='numeric on n_WP')
ax[0].axhline(0, color='.6', lw=.7); ax[0].set_xticks(x); ax[0].set_xticklabels(comps, fontsize=9)
ax[0].set_ylabel('energy (eV)'); ax[0].set_title('Quantum self-energy the ghost lacks: cross-checks')
ax[0].legend(frameon=False, fontsize=8)
for xi,v in zip(x-w,meas):
    if not np.isnan(v): ax[0].text(xi, v, f'{v:.0f}', ha='center', va='bottom' if v>=0 else 'top', fontsize=8)
for xi,v in zip(x,ana):
    if not np.isnan(v): ax[0].text(xi, v, f'{v:.0f}', ha='center', va='bottom', fontsize=8)
for xi,v in zip(x+w,num):
    if not np.isnan(v): ax[0].text(xi, v, f'{v:.0f}', ha='center', va='bottom', fontsize=8)
# (b) the charged-cell convention caveat: raw dHartree(r=40) at p2 vs p3 vs the physical self-Hartree
p2 = {'wp': -90.7419, 'cl': -80.6643}                        # E_hartree (Ha), r=40, from h0_p2
raw = [ (p2['wp']-p2['cl'])*HA, (W3[40]['hartree']-C3[40]['hartree'])*HA, E_self_ana*HA ]
bars = ax[1].bar(['raw dHartree\\n(p2, r=40)','raw dHartree\\n(p3, r=40)','physical\\nself-Hartree'],
                 raw, color=['#c0392b','#8e44ad','#16a085'], width=.6)
ax[1].axhline(0, color='.6', lw=.7); ax[1].set_ylabel('energy (eV)')
ax[1].set_title('Why the self-Hartree is NOT the raw component difference')
for rect,val in zip(bars,raw): ax[1].text(rect.get_x()+rect.get_width()/2, val,
        f'{val:.0f}', ha='center', va='bottom' if val>=0 else 'top', fontsize=9)
fig.suptitle('Panel 4 - the "missing" quantum self-energy, and the charged-cell caveat')
fig.tight_layout(); plt.show()
print(f'Zero-point KE : measured d(kinetic) = {dkin_meas:.3f} Ha ({dkin_meas*HA:.1f} eV)  vs  3/(4s^2) = {zp_ana:.3f} Ha ({zp_ana*HA:.1f} eV)  [0.1%]')
print(f'Self-Hartree  : analytic 1/(2 s_rho sqrt(pi)) = {E_self_ana:.3f} Ha ({E_self_ana*HA:.1f} eV)  vs  numeric on n_WP = {E_self_num:.3f} Ha ({E_self_num*HA:.1f} eV)  [~3%]')
print(f'Self-XC       : measured d(xc) = {dxc_meas:.3f} Ha ({dxc_meas*HA:.1f} eV)  (WP LDA self-XC)')
print(f'CAVEAT: raw dHartree(r=40) = {(W3[40]["hartree"]-C3[40]["hartree"])*HA:.0f} eV (p3) vs {(p2["wp"]-p2["cl"])*HA:.0f} eV (p2) -- convention-dependent, matches NEITHER the')
print(f'physical self-Hartree (+{E_self_ana*HA:.0f} eV): inserting the WP makes the cell net -1 charged, so the Poisson')
print('G=0 convention injects ~N_slab*mean(V[n_WP]) into the Hartree/external split. The self-energy is')
print('recovered from n_WP directly (analytic == numeric), not from a component subtraction. Verdict yours.')"""


def nb(cells):
    n = new_notebook(); n.cells = cells
    n.metadata.kernelspec = {"name": "python3", "display_name": "Python 3"}
    return n

def build():
    cells = [
        new_markdown_cell(
            "# Electrostatic sheet / slab model — the linear-field expectation\n\n"
            "*Theoretical model of the localised jellium (user's modelling), in **SI** "
            "units. φ = 0 reference at the plane; z_q measured from the slab centre; the "
            "slab uses the SAME parameters as the `campaign_autorun` runs (n₀ = 1.312e-3 "
            "e/Bohr³, thickness 25 Bohr). Built by `build_theoretical_model.py`.*"),
        new_markdown_cell(
            "## The model\n\n"
            "**Single sheet** — areal charge density $\\sigma$ (C/m²) at $z=0$:\n\n"
            "$$\\phi(z_q) = -\\frac{\\sigma}{2\\varepsilon_0}\\,|z_q|,\\qquad "
            "U = q\\phi = -\\frac{\\sigma q}{2\\varepsilon_0}\\,|z_q|.$$\n\n"
            "**Slab** — volume density $\\rho_0$ (C/m³), thickness $L$, centred; for "
            "$z_q>L/2$ (distance $r=z_q-L/2$ from the near face):\n\n"
            "$$\\phi(z_q) = -\\frac{\\rho_0 L}{2\\varepsilon_0}\\,z_q,\\qquad "
            "U = -\\frac{\\rho_0 L\\,q}{2\\varepsilon_0}\\,z_q.$$\n\n"
            "**Slab collapsed to one sheet** at its centre, $\\sigma_\\text{tot}=\\rho_0 L$:\n\n"
            "$$\\phi(z_q) = -\\frac{\\sigma_\\text{tot}}{2\\varepsilon_0}\\,|z_q|,\\qquad "
            "U = -\\frac{\\sigma_\\text{tot}\\,q}{2\\varepsilon_0}\\,|z_q|.$$\n\n"
            "Slab and collapsed sheet are identical for the **field** outside "
            "($|z_q|>L/2$, since $\\sigma_\\text{tot}=\\rho_0 L$); they differ only inside, "
            "where the true slab gives a **parabolic** $\\phi$. The hallmark: an infinite "
            "sheet has a **uniform** field, so $\\phi$ and $U$ are **linear in distance** — "
            "the magnitude *grows*, with no $1/r$ decay."),
        new_code_cell(CONST),
        new_markdown_cell("## Case 1 — single infinite sheet"),
        new_code_cell(SHEET),
        new_markdown_cell("## Cases 2 & 3 — uniform slab vs its collapsed sheet"),
        new_code_cell(SLAB),
        new_markdown_cell("## What the model says we should expect"),
        new_code_cell(EXPECT),
        new_markdown_cell(
            "## Extension — the neutral slab as two sheets (same model)\n"
            "The jellium slab is **neutral**: a positive background sheet and an electron "
            "sheet. Superposing their single-sheet potentials (within this electrostatic "
            "model) shows how the two linear terms combine — the bridge to the neutral-slab "
            "behaviour. Presented as model output, not a conclusion."),
        new_code_cell(SUPERPOSE),
        new_markdown_cell(
            "## v_bg from the Poisson solver vs the infinite plate (sanity check)\n"
            "The previous sections are the *analytic* expectation. Here we test it against the "
            "**actual** background potential the simulation uses: `v_bg = −poisson(n₊)`, dumped "
            "straight from INQ's **p2 Rozzi slab-truncated** Poisson solver "
            "(`scripts/campaign_autorun/dump_vbg/run.cpp`) as a z-lineout. If the localised "
            "background behaves like an infinite uniform plate with all the charge in the slab, "
            "`v_bg` must be **parabolic inside** and **linear outside** with slope `4πn₀a`, and "
            "the field `E_z = −dv_bg/dz` must be a **constant ≠ 0** outside (the infinite-sheet "
            "hallmark — it does *not* decay to zero). The absolute offset of `v_bg` is a p2 G=0 "
            "gauge constant, so the potential panel centres both curves at `v_bg(0)=0`; the field "
            "panel is gauge-free. **Note:** the non-zero far field here is the *background alone* "
            "(a charged plate) — the zero-field-far-away expectation only applies to the **net** "
            "(background + electrons) neutral density, which is a separate check."),
        new_code_cell(VBG),
        new_markdown_cell(
            "## Empirical-density plate model (real n(z) as a sheet stack)\n"
            "Self-contained. Takes the **measured** planar-mean electron profile n₋(z) from "
            "the periodicity-2 GS, keeps the analytic background n₊(z), and builds "
            "φ(z_q) = −1/(2ε₀)∫ρ(z')|z_q−z'|dz' by summing infinite sheets — so the surface "
            "spill-out / Friedel structure (not a sharp slab) sets the result. The NET curve "
            "uses the real, non-cancelling ρ = e(n₊−n₋)."),
        new_code_cell(EMPIRICAL),
        new_markdown_cell(
            "## Projectile-cutoff test — does the classical decay track the UPF r_max?\n"
            "The classical projectile UPF has **z_valence = 0** and radial grid **r_max = 50 Bohr** "
            "(a `1/r` tail that simply ends there). A point charge above an infinite sheet with "
            "its Coulomb truncated at r_cut sees only a finite disk, giving a 1D kernel "
            "`K(dz)=(r_cut−|dz|)` that **decays to 0** as the slab passes beyond r_cut — no "
            "screening needed. This overlays that model (r_cut=50) on the **simulated** classical "
            "excess to test whether the observed decay is the cutoff."),
        new_code_cell(CUTOFF),
        new_markdown_cell(
            "## Empirical cutoff sweep — 4 projectile UPFs (KS runs)\n"
            "The decisive test. Four projectile UPFs truncated at r_cut = 10 / 20 / 30 / 40 Bohr, "
            "each run through the classical radius sweep (periodicity 2, off the Lz=120 GS, "
            "full energy decomposition). If the classical `E_total(r)` decay is an artifact of "
            "the projectile potential's finite range, **each curve should fall to ~0 near its "
            "own cutoff** (dotted lines). This also confirms `E_ion` (part of `total()`, streamed) "
            "stays 0. Whether this is an *appreciable* effect at the simulated distances is the "
            "reader's call."),
        new_code_cell(CUTOFFRUNS),
        new_markdown_cell(
            "## Analytical improvement — image potential (metal surface)\n"
            "Beyond the rigid sheet: the jellium is a **metal**, so an external charge induces a "
            "screening image → `U_image(z) = −q²/(4z)`, decaying as ~1/z. This is the expected "
            "*physical* law once the projectile potential reaches far enough — but note the "
            "current runs' decay is dominated by the UPF cutoff above, so this is a hypothesis to "
            "test (e.g. with a longer-range projectile UPF), not a fit to the present data."),
        new_code_cell(IMAGE),
        new_markdown_cell(ASSUME),
        # ---- PART II: classical vs WP KS runs (deconstruction + screening test) ----
        new_markdown_cell(PART2_NOTE),
        new_markdown_cell(
            "## Full energy deconstruction — classical projectile vs WP electron\n"
            "Matched insertion runs (periodicity 2, same Lz=120 GS), every energy component "
            "streamed. The headline `dE_total(r)` shows the WP **flat** (~self-energy) while the "
            "classical **decays** (the cutoff term of Part I). The WP−classical difference is split "
            "into the r-independent **quantum self-energy** (zero-point KE + self-XC) and the "
            "r-dependent **electrostatic** part (Hartree+external summed — robust to the p2 G=0 offset)."),
        new_code_cell(LEDGER),
        new_markdown_cell(
            "## The E_xc difference and the total-energy difference (Learning #1)\n"
            "`E_xc(WP) − E_xc(classical)` isolates the real WP electron's LDA exchange-correlation "
            "(the classical projectile has none). The total difference resolves into a **constant** quantum "
            "self-energy plus an **r-dependent** electrostatic term that vanishes far away — linking "
            "back to the neutral-slab NET≈0 electrostatics of the Part I model."),
        new_code_cell(XCDIFF),
        new_markdown_cell(
            "## ★ Screening / WP-potential test (Learning #2) — the centrepiece\n"
            "**Does the wavepacket produce the same Coulomb potential as the classical Gaussian?** "
            "Equivalent potentials ⇔ equivalent **source charge**. New runs saved n_WP = |ψ_WP|² at "
            "t=0. We compare it to the ideal Gaussian charge (σ_ρ = σ_WP/√2 = 0.354) and overlay "
            "`poisson(n_WP)` on the analytic classical-projectile potential `erf(r/(√2 σ_ρ))/r` (the FFT-"
            "Poisson was validated against that analytic form to RMS 3e-4 Ha). σ = 0.5 is the WP "
            "label; σ_ρ appears only here, in the potential-generation footnote."),
        new_code_cell(WPPOT),
        new_markdown_cell(
            "## Screening baseline and the dynamical-screening caveat\n"
            "`n_slab(t=0) − n_GS = 0` **exactly** (bit-identical): at insertion the slab orbitals are "
            "still the GS orbitals, so there is no *instantaneous* screening. Genuine screening (the "
            "slab polarising around the WP) is a **dynamical** response over ~a plasmon period "
            "(~5000 CPU steps here) — flagged as a GPU follow-up, recorded honestly, not fabricated."),
        new_code_cell(BATH),
        # ---- PART III: single-run component decomposition at periodicity 3 ----
        new_markdown_cell(PART3_NOTE),
        new_markdown_cell(
            "## Loading the periodicity-3 mirror\n"
            "Same radius sweep as Part II, remade at **periodicity 3** with the full energy "
            "decomposition (`rerun_h0_p3.py`) so each component is individually offset-free "
            "(`G=0 -> 0`). Both projectiles share the p3 GS/box/grid; the decomposition is exact."),
        new_code_cell(P3DATA),
        new_markdown_cell(
            "## Panel 1 — the total-energy magnitude difference\n"
            "The single-run headline: `E_total(WP)` vs `E_total(classical)` at r = 28, and how that "
            "difference behaves across r. The sign change with r is the electrostatic part; the "
            "constant floor underneath is the quantum self-energy (Panel 4)."),
        new_code_cell(P3TOTAL),
        new_markdown_cell(
            "## Panel 2 — per-run component decomposition (each run individually)\n"
            "Waterfall of the individual components for the WP run and the classical run at r = 28. "
            "Periodicity 3 makes each physically signed — `E_hartree > 0` (e–e repulsion), "
            "`E_external < 0` (electrons in the +background well) — and each bar-stack sums to its "
            "`E_total`. Note the classical kinetic/Hartree/xc are the **frozen GS values**; only its "
            "`E_external` moves with r."),
        new_code_cell(P3WATER),
        new_markdown_cell(
            "## Panel 3 — component-wise difference (WP − classical) vs r\n"
            "Which components carry the difference. `d(kinetic)` and `d(xc)` are **flat and "
            "convention-free** (the WP's zero-point KE and self-XC). `d(hartree)` and `d(external)` "
            "are r-dependent **and** carry the net-charge Poisson-convention term — so they are shown "
            "but not over-interpreted (Panel 4 makes the caveat explicit)."),
        new_code_cell(P3DIFF),
        new_markdown_cell(
            "## Panel 4 — the quantum self-energy ledger (and the charged-cell caveat)\n"
            "What the classical ghost is **missing** relative to the WP, each cross-checked "
            "independently: **zero-point KE** (measured `d(kinetic)` vs `3/(4σ²)`), **self-Hartree** "
            "(closed form `1/(2σ_ρ√π)` vs numeric FFT-Poisson on `n_WP`), and **self-XC** (measured "
            "`d(xc)`). The second panel shows *why* the self-Hartree is **not** the raw "
            "`E_hartree(WP) − E_hartree(classical)`: that subtraction is convention-dependent (p2 vs "
            "p3 disagree) because the WP makes the cell net-charged; the physical self-energy is "
            "recovered from `n_WP` directly. Verdict is yours."),
        new_code_cell(P3SELF),
    ]
    p = OUT / "theoretical_slab_model.ipynb"
    nbf.write(nb(cells), str(p)); print("wrote", p.name); return p

if __name__ == "__main__":
    build()
    print("execute: python3 -m nbconvert --to notebook --execute --inplace theoretical_slab_model.ipynb (venv)")
