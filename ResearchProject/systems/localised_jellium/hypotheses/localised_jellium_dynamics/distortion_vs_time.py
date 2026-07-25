import numpy as np
from inqview.visualisation.field_io import load_vti
HA2EV = 27.211386
BASE = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "scripts/localised_jellium_dynamics/runs/p5/wp/results/wp/frames")
ref = load_vti(f"{BASE}/wp/density_t000000.vti")
x,y,z = ref.x,ref.y,ref.z
dx,dy,dz = x[1]-x[0],y[1]-y[0],z[1]-z[0]; dV=dx*dy*dz
V=(x[-1]-x[0]+dx)*(y[-1]-y[0]+dy)*(z[-1]-z[0]+dz)
kx=2*np.pi*np.fft.fftfreq(len(x),d=dx); ky=2*np.pi*np.fft.fftfreq(len(y),d=dy); kz=2*np.pi*np.fft.fftfreq(len(z),d=dz)
KX,KY,KZ=np.meshgrid(kx,ky,kz,indexing='ij'); G2=KX*KX+KY*KY+KZ*KZ; G2[0,0,0]=np.inf; KER=4*np.pi/G2
X,Y,Z=np.meshgrid(x,y,z,indexing='ij')
n0=82.0/(50.0*50.0*25.0); n_bg=np.where(np.abs(Z)<12.5,n0,0.0)
def ft(f): return np.fft.fftn(f)*dV
def J(f,g): return float((1.0/V)*np.sum(KER*np.real(np.conj(ft(f))*ft(g))))
sig=0.5/np.sqrt(2.0)
tot0=load_vti(f"{BASE}/total/density_t000000.vti").data
print(f"{'frame':>6} {'cz':>8} {'sig_wp':>8} {'selfH[eV]':>10} {'distort[eV]':>12}")
for t in (0,100,250,500):
    n_wp=load_vti(f"{BASE}/wp/density_t{t:06d}.vti").data
    nrm=n_wp.sum()*dV
    cx=(X*n_wp).sum()*dV/nrm; cy=(Y*n_wp).sum()*dV/nrm; cz=(Z*n_wp).sum()*dV/nrm
    var=(((X-cx)**2+(Y-cy)**2+(Z-cz)**2)*n_wp).sum()*dV/nrm; sig_wp=np.sqrt(var/3.0)
    n_proj=np.exp(-(((X-cx)**2+(Y-cy)**2+(Z-cz)**2))/(2*sig*sig)); n_proj*=1.0/(n_proj.sum()*dV)
    n_slab=load_vti(f"{BASE}/total/density_t{t:06d}.vti").data - n_wp
    SH=0.5*J(n_wp,n_wp)*HA2EV
    dist=J(n_slab-n_bg,n_wp-n_proj)*HA2EV
    print(f"{t:>6} {cz:8.3f} {sig_wp:8.4f} {SH:10.3f} {dist:12.3f}")
