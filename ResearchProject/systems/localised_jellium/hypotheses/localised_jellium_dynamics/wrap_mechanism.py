"""Decisive test: is the impl residual (7.4 vs clean 21.5 eV) the periodic WRAP
of the projectile pseudopotential? Reconstruct v_ion = erf(d/0.5)/d summed over
periodic images out to r_cut (INQ real-space UPF placement) vs minimum-image
(no wrap), and assemble the residual for each. If minimum-image -> ~21.5 and
real-space rcut=120 -> ~7.4, the wrap is the cause.
"""
import numpy as np
from scipy.special import erf
from inqview.visualisation.field_io import load_vti
HA2EV = 27.211386
BASE = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "scripts/localised_jellium_dynamics/runs/p5/wp/results/wp/frames")
wp = load_vti(f"{BASE}/wp/density_t000000.vti"); tot = load_vti(f"{BASE}/total/density_t000000.vti")
x, y, z = wp.x, wp.y, wp.z; dx, dy, dz = x[1]-x[0], y[1]-y[0], z[1]-z[0]; dV = dx*dy*dz
Lx, Ly, Lz = 50.0, 50.0, 120.0
n_wp = wp.data; n_slab = tot.data - n_wp
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
n0 = 82.0/(50.0*50.0*25.0); n_bg = np.where(np.abs(Z) < 12.5, n0, 0.0)
zc = -24.5; sig = 0.5/np.sqrt(2.0)                  # sigma_pot; UPF erf uses 0.5 = sqrt2*sig
n_proj = np.exp(-((X**2+Y**2+(Z-zc)**2))/(2*sig*sig)); n_proj *= 1.0/(n_proj.sum()*dV)

kx = 2*np.pi*np.fft.fftfreq(len(x), d=dx); ky = 2*np.pi*np.fft.fftfreq(len(y), d=dy); kz = 2*np.pi*np.fft.fftfreq(len(z), d=dz)
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij'); G2 = KX*KX+KY*KY+KZ*KZ; G2[0,0,0] = np.inf
def poisson(n): return np.real(np.fft.ifftn(np.fft.fftn(n)*(4*np.pi/G2)))
phi_plus = poisson(n_bg); v_bg = -phi_plus
ideal = float((n_proj*v_bg).sum()*dV)*HA2EV
print(f"INQ ideal (reported r=12) = 134.7 eV ;  recon ideal = {ideal:.2f} eV")

def v_ion_realspace(rcut, images):
    Vv = np.zeros_like(X)
    for ix in range(-images, images+1):
        for iy in range(-images, images+1):
            for iz in range(-1, 2):
                dxx = X-(ix*Lx); dyy = Y-(iy*Ly); dzz = Z-(zc+iz*Lz)
                d = np.sqrt(dxx*dxx+dyy*dyy+dzz*dzz)
                m = (d < rcut) & (d > 1e-9)
                Vv[m] += erf(d[m]/0.5)/d[m]
    return Vv

def v_ion_minimage():
    dxx = X-0.0; dxx -= Lx*np.round(dxx/Lx)
    dyy = Y-0.0; dyy -= Ly*np.round(dyy/Ly)
    dzz = Z-zc;  dzz -= Lz*np.round(dzz/Lz)
    d = np.sqrt(dxx*dxx+dyy*dyy+dzz*dzz); d = np.maximum(d, 1e-9)
    return erf(d/0.5)/d

def eproj_bg(v): return -float((n_bg*v).sum()*dV)*HA2EV
def e_e_proj(v): return  float((n_slab*v).sum()*dV)*HA2EV

SH = 0.5*float((poisson(n_wp)*n_wp).sum()*dV)*HA2EV
J_slab_wp = float((poisson(n_slab)*n_wp).sum()*dV)*HA2EV
int_nwp_vbg = float((n_wp*v_bg).sum()*dV)*HA2EV
print(f"\nself-Hartree 1/2 J[nWP,nWP]      = {SH:.2f} eV")
print(f"J[n_slab,n_WP] (WP-slab repuls)  = {J_slab_wp:.2f} eV")
print(f"int n_WP.v_bg (WP-background)     = {int_nwp_vbg:.2f} eV   (INQ ~134.6)")

for tag, v in [("minimage(no wrap)", v_ion_minimage()),
               ("realspace rcut=50", v_ion_realspace(50, 1)),
               ("realspace rcut=120", v_ion_realspace(120, 2))]:
    epb = eproj_bg(v); eep = e_e_proj(v)
    resid = SH + (J_slab_wp-eep) + (int_nwp_vbg-epb)
    print(f"\n[{tag}]")
    print(f"  E_proj_bg = -int n+.v_ion             = {epb:9.2f} eV   (INQ impl rc120 r12 = -524.55)")
    print(f"  e-proj    =  int n_slab.v_ion          = {eep:9.2f} eV   (INQ E_ext_CL-E_ext_GS rc120 = 532.8)")
    print(f"  residual  = SH +(Jsw-eep)+(nwpvbg-epb) = {resid:8.2f} eV   (INQ measured: rc50~14, rc120=7.4)")
