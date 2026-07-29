#!/usr/bin/env python3
"""make_cutoff_upfs.py — generate classical-ghost UPFs truncated at several radial
cutoffs, to test whether the classical E_total(r) decay is set by the projectile
potential's finite range (the UPF r_max).

The source UPF `electron_gaussian_wpsigma0p5.upf` has a uniform radial grid
(dr = 0.01 Bohr, 5001 points -> r_max = 50 Bohr) carrying a 1/r tail. We truncate
the grid to r_cut in {10, 20, 30, 40} Bohr (keeping the first N = round(r_cut/dr)+1
points), re-wrapping PP_R / PP_RAB / PP_LOCAL / PP_RHOATOM and updating the size=
attributes + mesh_size. Beyond r_cut the potential is simply absent (z_valence=0 has
no long-range part), so the ghost's field ends at r_cut — the exact effect we test.

Output: cutoff_test/upfs/electron_gaussian_wpsigma0p5_rc{10,20,30,40}.upf
"""
from __future__ import annotations
import re
from pathlib import Path

SRC = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
           "shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf")
OUT = Path(__file__).resolve().parent / "upfs"; OUT.mkdir(parents=True, exist_ok=True)
CUTOFFS = [10, 20, 30, 40]   # Bohr

def _floats(block: str):
    return [float(x) for x in block.split()]

def _wrap(vals, ncol, fmt):
    lines = []
    for i in range(0, len(vals), ncol):
        lines.append("".join(fmt % v for v in vals[i:i+ncol]))
    return "\n".join(lines)

def truncate(text: str, N: int) -> str:
    # helper: replace a <PP_X ...>DATA</PP_X> block, keeping first N values, re-wrapping
    def repl(tag, ncol, fmt):
        nonlocal text
        m = re.search(rf'(<{tag}\b[^>]*>)(.*?)(</{tag}>)', text, re.S)
        assert m, f"{tag} not found"
        head, data, tail = m.group(1), m.group(2), m.group(3)
        vals = _floats(data)[:N]
        head = re.sub(r'size="\s*\d+"', f'size="{N}"', head)
        text = text[:m.start()] + head + "\n" + _wrap(vals, ncol, fmt) + "\n" + tail + text[m.end():]
    repl("PP_R",       8, "%10.4f")
    repl("PP_RAB",     8, "%10.4f")
    repl("PP_LOCAL",   4, "%21.10E")
    repl("PP_RHOATOM", 4, "%21.10E")
    # mesh_size in PP_HEADER
    text = re.sub(r'mesh_size="\s*\d+"', f'mesh_size="{N}"', text)
    return text

def main():
    src = SRC.read_text()
    # confirm uniform grid dr
    rblock = re.search(r'<PP_R\b[^>]*>(.*?)</PP_R>', src, re.S).group(1)
    R = _floats(rblock); dr = R[1] - R[0]
    print(f"source: {len(R)} points, dr={dr:.4f}, r_max={R[-1]:.2f} Bohr")
    for rc in CUTOFFS:
        N = int(round(rc / dr)) + 1
        out = truncate(src, N)
        p = OUT / f"electron_gaussian_wpsigma0p5_rc{rc}.upf"
        p.write_text(out)
        # verify
        R2 = _floats(re.search(r'<PP_R\b[^>]*>(.*?)</PP_R>', out, re.S).group(1))
        V2 = _floats(re.search(r'<PP_LOCAL\b[^>]*>(.*?)</PP_LOCAL>', out, re.S).group(1))
        print(f"  rc={rc}: N={N}, r_max={R2[-1]:.2f} Bohr, |PP_R|={len(R2)}, |PP_LOCAL|={len(V2)}  -> {p.name}")

if __name__ == "__main__":
    main()
