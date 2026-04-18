#!/usr/bin/env python3
"""
Generate coronene (C24H12) geometry with D6h symmetry.

Derivation
----------
Coronene consists of 7 fused hexagonal rings.  All C-C bonds are idealised
to d_CC = 1.421 Å (aromatic graphene value), C-H = 1.086 Å.

Carbon atoms fall on three concentric shells:
  Shell 1 (inner ring)  : 6 C at r = d_CC,   angles 0°, 60°, …, 300°
  Shell 2 (junctions)   : 6 C at r = 2·d_CC, angles 0°, 60°, …, 300°
  Shell 3 (perimeter)   : 12 C at r = d_CC·√7, alternating angle pairs
                          θ = arctan(√3/5) and θ = arctan(√3/2) in each 60° sector

The 12 H atoms are collinear with their bonded perimeter C and the origin,
at r = d_CC·√7 + d_CH.

Reference bond lengths: C-C 1.421 Å (graphene), C-H 1.086 Å (aromatic).
Both are standard values widely used in DFT benchmark sets.
"""

import numpy as np
import os

ANG_TO_BOHR = 1.8897259886

def generate_coronene(d_CC=1.421, d_CH=1.086):
    """
    Return (carbon_positions, hydrogen_positions) in Angstrom,
    with the molecule centred at the origin in the xy-plane (z=0).
    """
    sqrt3 = np.sqrt(3.0)
    sqrt7 = np.sqrt(7.0)

    C = []  # carbon positions
    H = []  # hydrogen positions

    # ── Shell 1 and Shell 2 (same angles, different radii) ──────────────────
    for k in range(6):
        theta = np.radians(k * 60.0)
        C.append((d_CC * np.cos(theta),   d_CC * np.sin(theta),   0.0))  # shell 1
        C.append((2*d_CC * np.cos(theta), 2*d_CC * np.sin(theta), 0.0))  # shell 2

    # ── Shell 3 (perimeter C) and H atoms ───────────────────────────────────
    r3 = d_CC * sqrt7
    # Two distinct angles per 60° sector, derived from geometry:
    #   atom at (5d/2, d√3/2)  → θ_A = arctan(√3/5)  ≈ 19.106°
    #   atom at (2d,   d√3  )  → θ_B = arctan(√3/2)  ≈ 40.894°
    theta_A = np.degrees(np.arctan2(d_CC * sqrt3 / 2.0, 5.0 * d_CC / 2.0))
    theta_B = np.degrees(np.arctan2(d_CC * sqrt3,       2.0 * d_CC))

    for k in range(6):
        for theta_deg in [theta_A + k * 60.0, theta_B + k * 60.0]:
            theta = np.radians(theta_deg)
            cx = r3 * np.cos(theta)
            cy = r3 * np.sin(theta)
            C.append((cx, cy, 0.0))
            # H directly outward from center through this C atom
            rH = r3 + d_CH
            H.append((rH * np.cos(theta), rH * np.sin(theta), 0.0))

    return np.array(C), np.array(H)


def verify_bonds(C, H, d_CC, d_CH, tol=1e-3):
    """Check that all expected C-C and C-H bonds are present.

    Coronene bond counts (Euler's formula, V=24, F=8):
      C-C bonds: 30  (inner C: 3 each, junction C: 3 each, perimeter C: 2 each)
      C-H bonds: 12  (perimeter C only)
    """
    print("Bond verification:")

    # Count C-C bonds per atom
    n_CC_bonds = 0
    ok = True
    for i in range(len(C)):
        cc_count = 0
        for j in range(len(C)):
            if i == j:
                continue
            d = np.linalg.norm(C[i] - C[j])
            if abs(d - d_CC) < tol:
                cc_count += 1
                n_CC_bonds += 1
        # Count C-H bonds for this C atom
        ch_count = sum(1 for h in H if abs(np.linalg.norm(h - C[i]) - d_CH) < tol)
        total = cc_count + ch_count
        if total != 3:
            print(f"  WARNING: C atom {i} r={np.linalg.norm(C[i]):.3f} Å"
                  f" has {cc_count} C-C + {ch_count} C-H = {total} bonds (expected 3)")
            ok = False

    n_CC_bonds //= 2  # each bond counted twice

    n_CH_bonds = sum(
        1 for h in H for c in C if abs(np.linalg.norm(h - c) - d_CH) < tol
    )

    print(f"  C-C bonds found: {n_CC_bonds}  (expected 30)")
    print(f"  C-H bonds found: {n_CH_bonds}  (expected 12)")
    if n_CC_bonds == 30 and n_CH_bonds == 12 and ok:
        print("  All bonds correct ✓")
    else:
        print("  Bond count mismatch!")


def main():
    d_CC = 1.421  # Angstrom
    d_CH = 1.086  # Angstrom

    C_ang, H_ang = generate_coronene(d_CC=d_CC, d_CH=d_CH)

    print(f"Coronene geometry: {len(C_ang)} C + {len(H_ang)} H = {len(C_ang)+len(H_ang)} atoms")
    print(f"C-C bond length  : {d_CC} Å")
    print(f"C-H bond length  : {d_CH} Å\n")

    verify_bonds(C_ang, H_ang, d_CC, d_CH)

    # ── Print radii summary ──────────────────────────────────────────────────
    r_C = np.linalg.norm(C_ang, axis=1)
    r_H = np.linalg.norm(H_ang, axis=1)
    print(f"\nC radii (Å): {np.sort(np.unique(np.round(r_C, 4)))}")
    print(f"H radius  (Å): {np.unique(np.round(r_H, 4))}")

    # ── Convert to bohr and write XYZ file ──────────────────────────────────
    out_dir = os.path.dirname(os.path.abspath(__file__))
    xyz_path = os.path.join(out_dir, "coronene.xyz")
    with open(xyz_path, "w") as f:
        n_atoms = len(C_ang) + len(H_ang)
        f.write(f"{n_atoms}\n")
        f.write(f"Coronene C24H12, D6h symmetry, C-C={d_CC} Ang, C-H={d_CH} Ang\n")
        for pos in C_ang:
            f.write(f"C  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}\n")
        for pos in H_ang:
            f.write(f"H  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}\n")
    print(f"\nWrote: {xyz_path}")

    # ── Write INQ-compatible C++ insert block (positions in bohr) ───────────
    cpp_path = os.path.join(out_dir, "coronene_inq_insert.txt")
    with open(cpp_path, "w") as f:
        f.write("// Coronene C24H12 — generated by gen_geometry.py\n")
        f.write("// D6h symmetry, C-C=1.421 Å, C-H=1.086 Å\n")
        f.write("// Molecule in xy-plane (z=0), centred at box midpoint (Lx/2, Ly/2, Lz/2)\n")
        f.write("// Positions in bohr, shifted to box centre\n")
        f.write("const double Lx = 34.76, Ly = 34.76, Lz = 59.91;  // bohr\n")
        f.write("const double cx = Lx/2, cy = Ly/2, cz = Lz/2;\n\n")
        for pos in C_ang:
            bx = pos[0] * ANG_TO_BOHR
            by = pos[1] * ANG_TO_BOHR
            f.write(f'ions.insert("C", {{(cx+{bx:.6f})_b, (cy+{by:.6f})_b, cz_b}});\n')
        for pos in H_ang:
            bx = pos[0] * ANG_TO_BOHR
            by = pos[1] * ANG_TO_BOHR
            f.write(f'ions.insert("H", {{(cx+{bx:.6f})_b, (cy+{by:.6f})_b, cz_b}});\n')
    print(f"Wrote: {cpp_path}")

    # ── Print bohr coordinates (for quick inspection) ──────────────────────
    print("\nC atoms in bohr (relative to centre):")
    for i, pos in enumerate(C_ang):
        bx, by = pos[0]*ANG_TO_BOHR, pos[1]*ANG_TO_BOHR
        r_b = np.sqrt(bx**2 + by**2)
        print(f"  C{i+1:2d}: ({bx:8.4f}, {by:8.4f})  r={r_b:.4f} bohr")
    print("\nH atoms in bohr (relative to centre):")
    for i, pos in enumerate(H_ang):
        bx, by = pos[0]*ANG_TO_BOHR, pos[1]*ANG_TO_BOHR
        r_b = np.sqrt(bx**2 + by**2)
        print(f"  H{i+1:2d}: ({bx:8.4f}, {by:8.4f})  r={r_b:.4f} bohr")


if __name__ == "__main__":
    main()
