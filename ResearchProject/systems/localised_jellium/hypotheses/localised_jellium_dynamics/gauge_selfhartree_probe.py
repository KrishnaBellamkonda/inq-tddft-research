import numpy as np

HA2EV = 27.211386

def self_hartree_periodic(L, dx, sigma_rho):
    """WP Hartree self-energy on a periodic grid, INQ convention (G=0 dropped).
    Charge magnitude 1 (electron), Gaussian density std = sigma_rho, centered."""
    Lx, Ly, Lz = L
    nx, ny, nz = int(round(Lx/dx)), int(round(Ly/dx)), int(round(Lz/dx))
    dV = (Lx/nx)*(Ly/ny)*(Lz/nz)
    x = (np.arange(nx) - nx//2)*(Lx/nx)
    y = (np.arange(ny) - ny//2)*(Ly/ny)
    z = (np.arange(nz) - nz//2)*(Lz/nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    r2 = X*X + Y*Y + Z*Z
    n = np.exp(-r2/(2.0*sigma_rho**2))
    n /= (n.sum()*dV)                      # normalize to charge 1
    ntil = np.fft.fftn(n)*dV               # n(G) = integral n e^{-iGr} dr ; n(0)=1
    kx = 2*np.pi*np.fft.fftfreq(nx, d=Lx/nx)
    ky = 2*np.pi*np.fft.fftfreq(ny, d=Ly/ny)
    kz = 2*np.pi*np.fft.fftfreq(nz, d=Lz/nz)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    G2 = KX*KX + KY*KY + KZ*KZ
    G2[0,0,0] = np.inf                     # drop G=0 (INQ / jellium convention)
    V = Lx*Ly*Lz
    EH = 0.5/V * np.sum(4*np.pi*np.abs(ntil)**2 / G2)
    return EH, (nx, ny, nz)

sigma_rho = 0.5/np.sqrt(2.0)               # = 0.3536 Bohr (density std; sigma_WP=0.5)
analytic = 1.0/(2.0*sigma_rho*np.sqrt(np.pi))
print(f"sigma_rho (density std)      = {sigma_rho:.4f} Bohr")
print(f"ANALYTIC free-space E_self   = {analytic:.4f} Ha = {analytic*HA2EV:.3f} eV")
print()

# 1) On the actual run grid (Lx=Ly=50, Lz=120, dx=0.5)
EH_run, dims = self_hartree_periodic((50.,50.,120.), 0.5, sigma_rho)
print(f"RUN GRID (50,50,120) dx=0.5  dims={dims}")
print(f"  periodic E_self (G=0 drop) = {EH_run:.4f} Ha = {EH_run*HA2EV:.3f} eV")
print(f"  shift vs free-space        = {(EH_run-analytic)*HA2EV:+.3f} eV")
print()

# 2) Box-size scaling: cubic boxes, show convergence to analytic as L grows
print("BOX-SIZE SCALING (cubic L, dx=0.5):")
print(f"  {'L':>5}  {'E_self[eV]':>11}  {'shift[eV]':>10}  {'-a/(2L)*Ha->eV predict':>22}")
for Lc in (20., 30., 40., 50., 70., 100., 150.):
    EH, _ = self_hartree_periodic((Lc,Lc,Lc), 0.5, sigma_rho)
    shift = (EH-analytic)*HA2EV
    mp = -2.837/(2*Lc)*HA2EV               # Makov-Payne monopole prediction (eV)
    print(f"  {Lc:5.0f}  {EH*HA2EV:11.3f}  {shift:10.3f}  {mp:22.3f}")
