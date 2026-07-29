"""Consistent-ideal energy-balance recomputation from saved p5 densities.

Tests whether the localised-jellium WP is a rigid Gaussian, and recomputes the
electrostatic residual d(E_H+E_ext)-E_proj_bg with the projectile represented by
the IDEAL analytic Gaussian (no pseudopotential cutoff/wrap), so the only
r-independent survivor is the WP self-Hartree + (slab-net-field x WP-distortion).

Key identity (both n_WP and n_proj integrate to 1 -> their difference is neutral
-> convention-free, no G=0 ambiguity):
   residual_ideal = 1/2 J[n_WP,n_WP]  +  J[n_slab - n_bg, n_WP - n_proj]
                  = (WP self-Hartree)  +  (slab net field . WP distortion)
"""
import numpy as np
from inqview.visualisation.field_io import load_vti

HA2EV = 27.211386
BASE = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "scripts/localised_jellium_dynamics/runs/p5/wp/results/wp/frames")

wp    = load_vti(f"{BASE}/wp/density_t000000.vti")
tot   = load_vti(f"{BASE}/total/density_t000000.vti")
x, y, z = wp.x, wp.y, wp.z
dx, dy, dz = x[1]-x[0], y[1]-y[0], z[1]-z[0]
dV = dx*dy*dz; V = (x[-1]-x[0]+dx)*(y[-1]-y[0]+dy)*(z[-1]-z[0]+dz)
n_wp   = wp.data
n_tot  = tot.data
n_slab = n_tot - n_wp

# --- reciprocal grid, G=0-dropped Coulomb kernel ---
kx = 2*np.pi*np.fft.fftfreq(len(x), d=dx)
ky = 2*np.pi*np.fft.fftfreq(len(y), d=dy)
kz = 2*np.pi*np.fft.fftfreq(len(z), d=dz)
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
G2 = KX*KX + KY*KY + KZ*KZ; G2[0,0,0] = np.inf
KER = 4*np.pi/G2

def ft(f):  return np.fft.fftn(f)*dV                 # \tilde f(G) = int f e^{-iGr}
def J(f, g):                                          # int int f(r)g(r')/|r-r'|
    return float((1.0/V)*np.sum(KER*np.real(np.conj(ft(f))*ft(g))))

# --- WP geometry: norm, centroid, width ---
norm = n_wp.sum()*dV
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
cz = (Z*n_wp).sum()*dV/norm; cx = (X*n_wp).sum()*dV/norm; cy = (Y*n_wp).sum()*dV/norm
var = (((X-cx)**2+(Y-cy)**2+(Z-cz)**2)*n_wp).sum()*dV/norm   # <r^2> (3D)
sig_wp = np.sqrt(var/3.0)                                    # per-axis density std
print("=== WP geometry (actual saved density) ===")
print(f"  norm            = {norm:.6f}  (expect 1)")
print(f"  centroid (x,y,z)= ({cx:.3f}, {cy:.3f}, {cz:.3f})  (expect ~0,0,-24.5)")
print(f"  density std/axis= {sig_wp:.4f} Bohr  (ideal sigma_pot = {0.5/np.sqrt(2):.4f})")

# --- ideal Gaussian projectile: charge 1, density std 0.354, at WP centroid ---
sig = 0.5/np.sqrt(2.0)
n_proj = np.exp(-(((X-cx)**2+(Y-cy)**2+(Z-cz)**2))/(2*sig*sig))
n_proj *= 1.0/(n_proj.sum()*dV)

# --- analytic +82 background slab (n0 for |z|<12.5, edge_width=0) ---
n0 = 82.0/(50.0*50.0*25.0)
n_bg = np.where(np.abs(Z) < 12.5, n0, 0.0)
print(f"  background charge= {n_bg.sum()*dV:.3f} (expect 82); slab electrons = {n_slab.sum()*dV:.3f}")

# --- energies ---
SH_wp    = 0.5*J(n_wp,   n_wp)   *HA2EV     # actual WP self-Hartree
SH_ideal = 0.5*J(n_proj, n_proj) *HA2EV     # ideal Gaussian self-Hartree (machinery check ~21.5)
AB       = J(n_slab - n_bg, n_wp - n_proj)*HA2EV   # slab-net-field . WP-distortion (convention-free)
res_ideal = SH_wp + AB
dE_xc = -16.47
combined = res_ideal + dE_xc

print("\n=== Self-Hartree ===")
print(f"  1/2 J[n_WP , n_WP ] (ACTUAL wp) = {SH_wp:8.3f} eV")
print(f"  1/2 J[n_prj, n_prj] (IDEAL   )  = {SH_ideal:8.3f} eV   (analytic 21.71; FFT check)")
print(f"  -> WP vs ideal self-Hartree     = {SH_wp-SH_ideal:+.3f} eV  (0 => WP is the rigid Gaussian)")

print("\n=== Consistent-IDEAL residual  d(E_H+E_ext) - E_proj_bg  (ideal projectile) ===")
print(f"  self-Hartree  1/2 J[n_WP,n_WP]            = {SH_wp:8.3f} eV")
print(f"  distortion    J[n_slab-n_bg, n_WP-n_proj] = {AB:8.3f} eV")
print(f"  ------------------------------------------------------")
print(f"  residual_ideal                           = {res_ideal:8.3f} eV")
print(f"  (measured impl residual: rc50 ~14.0, rc120 ~7.4 eV)")

print("\n=== Combined quantity  (dE_xc + residual) ===")
print(f"  dE_xc (WP self-XC)                       = {dE_xc:8.3f} eV")
print(f"  combined_ideal = residual_ideal + dE_xc  = {combined:8.3f} eV")
print(f"  free-space SIE = 21.71 + (-16.47)        = {21.712-16.47:8.3f} eV")
