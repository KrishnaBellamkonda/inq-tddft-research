#!/usr/bin/env python3
"""Extend the classical-ghost UPF radial cutoff to r_cut = 120 Bohr (= Lz).

The source `electron_gaussian_wpsigma0p5.upf` carries V(r)=+erf(r/0.5)/r Ha on a
uniform grid dr=0.01 out to r_max=50 Bohr; beyond ~50 Bohr this is a pure Coulomb
tail (erf saturated). We APPEND points 50.01..120.00 with the exact continuation
PP_LOCAL(r)=2/r Ry (= 2*erf/r, erf≈1), PP_RAB=dr, PP_RHOATOM=0, and update the
mesh_size + size= attributes. Result: the ghost potential spans the whole box, so
it never truncates the slab for any projectile position (campaign
localised-jellium-dynamics-analysis, r_cut=Lz test).
"""
import re
from pathlib import Path

SRC = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/"
           "pseudopotentials/electron_gaussian_wpsigma0p5.upf")
R_CUT = 120.0
OUT = SRC.parent / "electron_gaussian_wpsigma0p5_rc120.upf"

def floats(block): return [float(x) for x in block.split()]
def wrap(vals, ncol, fmt):
    return "\n".join("".join(fmt % v for v in vals[i:i+ncol]) for i in range(0, len(vals), ncol))

t = SRC.read_text()
R = floats(re.search(r'<PP_R\b[^>]*>(.*?)</PP_R>', t, re.S).group(1))
dr = R[1] - R[0]; r_last = R[-1]
n_new = int(round((R_CUT - r_last) / dr))
r_ext = [round(r_last + (i + 1) * dr, 6) for i in range(n_new)]   # 50.01 .. 120.00
N = len(R) + len(r_ext)
print(f"source: {len(R)} pts to r={r_last}; appending {len(r_ext)} pts to r={r_ext[-1]}; total N={N}")

ext = {
    "PP_R":       (r_ext,                       8, "%10.4f"),
    "PP_RAB":     ([dr] * len(r_ext),           8, "%10.4f"),
    "PP_LOCAL":   ([2.0 / r for r in r_ext],    4, "%21.10E"),   # Ry, pure Coulomb 2/r
    "PP_RHOATOM": ([0.0] * len(r_ext),          4, "%21.10E"),
}
for tag, (vals, ncol, fmt) in ext.items():
    m = re.search(rf'(<{tag}\b[^>]*>)(.*?)(</{tag}>)', t, re.S)
    head, data, tail = m.group(1), m.group(2), m.group(3)
    old = floats(data)
    head = re.sub(r'size="\s*\d+"', f'size="{N}"', head)
    newdata = "\n" + wrap(old + vals, ncol, fmt) + "\n"
    t = t[:m.start()] + head + newdata + tail + t[m.end():]
t = re.sub(r'mesh_size="\s*\d+"', f'mesh_size="{N}"', t)
OUT.write_text(t)
print("wrote", OUT)

# verify
t2 = OUT.read_text()
R2 = floats(re.search(r'<PP_R\b[^>]*>(.*?)</PP_R>', t2, re.S).group(1))
V2 = floats(re.search(r'<PP_LOCAL\b[^>]*>(.*?)</PP_LOCAL>', t2, re.S).group(1))
print(f"verify: PP_R n={len(R2)} r_max={R2[-1]:.2f}; PP_LOCAL n={len(V2)} "
      f"V[r=1]={V2[100]:.5f} V[r=100]={V2[9999]:.6f} (expect {2/100:.6f})")
assert len(R2) == N and abs(R2[-1] - R_CUT) < 1e-6 and len(V2) == N, "extension mismatch"
print("OK")
