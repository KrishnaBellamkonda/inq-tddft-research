#!/usr/bin/env python3
"""Q&A viii: the WP z-centroid stalls — the packet is NON-RIGID.

Plots the z-axis centroid <z>(t) for both runs to show that the wavepacket does
not translate rigidly: its (survival-weighted) centroid decelerates and stalls at
<z> ~ +5 Bohr, never reaching the far slab face (+12.5), while sigma_z balloons
~32x. The classical projectile (overlaid, for comparison) marches ballistically
through the slab to the +z box edge, decelerating only 2.71 -> 2.35 (the real
stopping signal).

Important: the WP <z> is the centroid of the *un-absorbed* (surviving) WP,
<z> = int z|psi|^2 / int |psi|^2. The +z CAP removes the fast forward components
(pz_mean 2.63 -> 1.71, e_kin 6.65 -> 4.31 Ha), so the surviving remainder is slow /
back-scattered and the centroid is dragged back. The stall is therefore a
*survival-weighted absorption signature*, NOT the packet physically stopping.

Inputs : wp_real_space_stats.csv (WP <z>, sigma_z) ; electron_track.csv (classical z, vz).
Outputs: qa_viii_centroid_trajectory.png + .csv in this directory.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from inqview.visualisation import style

style.apply_theme()

ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
WP = f"{ROOT}/scripts/fullsuite_wp/results/p5_wp/raw"
CL = f"{ROOT}/scripts/fullsuite_classical/results/p5_classical/raw"
OUT = f"{ROOT}/hypotheses/03_cap_stopping"

Z_FACE = 12.5      # slab faces +/-
Z_CAP = 17.5       # CAP inner edges +/-
Z_BOX = 25.0       # box edges +/-

# --- WP centroid + spread (survival-weighted) ---
wp = pd.read_csv(f"{WP}/observables/wp_real_space_stats.csv", comment="#")
tw = wp.time_au.values
zc = wp.z_mean.values
sz = np.sqrt(wp.sigma_z2.values)
vzc = np.gradient(zc, tw)            # d<z>/dt — centroid velocity

# --- classical projectile (unwrapped track) ---
cl = pd.read_csv(f"{CL}/observables/electron_track.csv").drop_duplicates("step")
tc = cl.time_au.values
zion = cl.z.values                   # unwrapped (reaches +30 by t=18)
vz = cl.vz.values

# --- summary numbers (for the notebook / handover) ---
z_stall = float(np.median(zc[tw >= 9.0]))
z_max = float(zc.max())
print(f"WP : sigma_z {sz[0]:.2f} -> {sz[-1]:.2f}  ({sz[-1]/sz[0]:.0f}x spread)")
print(f"WP : <z> launch {zc[0]:.2f} -> stall ~{z_stall:.2f} (max {z_max:.2f}); "
      f"never reaches far face +{Z_FACE}")
print(f"WP : centroid velocity d<z>/dt {vzc[0]:.2f} -> {vzc[-1]:.2f} (collapses)")
print(f"classical: z {zion[0]:.1f} -> {zion[-1]:.1f}; vz {vz[0]:.3f} -> {vz[-1]:.3f}")

# --- plot: (a) centroid trajectory, (b) centroid velocity ---
fig, (axA, axB) = plt.subplots(2, 1, figsize=(7.0, 7.2), sharex=True)

# slab band + region edges
for ax in (axA,):
    ax.axhspan(-Z_FACE, Z_FACE, color="0.85", alpha=0.6, zorder=0, label="jellium slab")
    for zz in (Z_CAP, -Z_CAP):
        ax.axhline(zz, ls="--", lw=1.0, color="C2", alpha=0.7)
    for zz in (Z_BOX, -Z_BOX):
        ax.axhline(zz, ls=":", lw=1.0, color="0.5")
    ax.axhline(Z_FACE, ls="-", lw=0.8, color="0.4")
    ax.axhline(-Z_FACE, ls="-", lw=0.8, color="0.4")

# (a) centroid trajectory
axA.fill_between(tw, zc - sz, zc + sz, color="C0", alpha=0.18,
                 label=r"WP $\langle z\rangle \pm \sigma_z$ (spread)")
axA.plot(tw, zc, "C0-", lw=2.0, label=r"WP $\langle z\rangle$ (survival-weighted)")
axA.plot(tc, zion, "C3-", lw=1.8, label="classical projectile $z$ (ballistic)")
axA.axhline(z_stall, ls="-.", lw=0.9, color="C0", alpha=0.6)
axA.annotate(rf"WP centroid stalls at $\langle z\rangle\!\approx\!{z_stall:.1f}$"
             "\n(never reaches far face +12.5)",
             xy=(12.5, z_stall), xytext=(8.5, -8.0), fontsize=8, color="C0",
             arrowprops=dict(arrowstyle="->", color="C0", lw=1.0))
axA.text(0.3, Z_FACE + 0.6, "far face +12.5", fontsize=7, color="0.3")
axA.text(0.3, Z_CAP + 0.6, "+z CAP edge +17.5", fontsize=7, color="C2")
axA.set_ylabel("z-centroid / position (Bohr)")
axA.set_ylim(-26, 32)
axA.set_title("z-centroid vs time — WP stalls (non-rigid), classical marches ballistically",
              fontsize=9)
axA.legend(fontsize=7, frameon=False, loc="upper left", ncol=1)
axA.grid(alpha=0.25)

# (b) centroid velocity
axB.axhline(0.0, ls="-", lw=0.8, color="0.5")
axB.plot(tw, vzc, "C0-", lw=2.0, label=r"WP $d\langle z\rangle/dt$ (centroid velocity)")
axB.plot(tc, vz, "C3-", lw=1.8, label=r"classical $v_z$ (2.71$\to$2.35)")
axB.set_xlabel("time (a.u.)")
axB.set_ylabel("velocity (a.u.)")
axB.set_title(r"Centroid velocity: WP collapses to ~0 (absorption-weighted), "
              "classical stays ~2.4", fontsize=9)
axB.legend(fontsize=7, frameon=False, loc="upper right")
axB.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(f"{OUT}/qa_viii_centroid_trajectory.png", dpi=200)
plt.close(fig)

pd.DataFrame({"time_au": tw, "wp_zc": zc, "wp_sigma_z": sz,
              "wp_centroid_vz": vzc}).to_csv(
    f"{OUT}/qa_viii_centroid_trajectory_wp.csv", index=False)
pd.DataFrame({"time_au": tc, "cl_z": zion, "cl_vz": vz}).to_csv(
    f"{OUT}/qa_viii_centroid_trajectory_cl.csv", index=False)
print("wrote qa_viii_centroid_trajectory.png + .csv")
