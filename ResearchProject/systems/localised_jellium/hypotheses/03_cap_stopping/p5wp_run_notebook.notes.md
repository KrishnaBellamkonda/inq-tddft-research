_Your notes/TODOs for this run. Edit freely — the builder pins this at the top of
`p5wp_run_notebook.ipynb` on every rebuild and never overwrites it._

- Need to plot the approximate exit time of the wave packet (given its travelling at
  the mean momentum, the time taken for it to reach the boundary on the other end of
  the slab.) — **[addressed: "WP transit & ballistic exit time" section, t_exit ≈ 10.33 a.u.]**
- Need to compare the stopping power estimate with the analytical value. Use the
  analytical value to guide where we might need the cutoff to be. It has to [be] made
  for the required density of jellium slab. — **[addressed: "Stopping power vs
  analytical Lindhard" section, r_s=3.995; curves bracket the measured S]**

## New notes (2026-06-23) — open, not yet done

### 1. CAP reflection in the wavepacket case
- Potential **reflection** may be happening off the CAP in the wavepacket case →
  need to **check the parameters of the CAP carefully**.
- Need to understand and **quantify the % absorption vs % reflection off EACH of the
  boundaries** (the −z and +z CAPs separately).

### 2. Quantum-case stopping power (equilibrated-system energy difference)
- Define the **"system"** as the state that **equilibrates at high time**: by then all
  energy transfer has happened, and — assuming any outgoing packets from
  **back-scattering and forward-scattering are absorbed** by the CAP — `E_total =
  E_system`. We are then left with only the system → this is **`E_final`**.
- **`E_initial`** is known (important: subtract the projectile's injected **−100 eV**
  to get the total energy of the system).
- `dE = E_final − E_initial`; `dx` = total width of the jellium slab the projectile
  went through. Then **`S = dE/dx`** is the quantum stopping power.
- Compare this to the classical case and **plot it on the SAME analytical Lindhard
  stopping** curve.

### 3. Future run: non-spreading wavepacket
- Design a future run where the wavepacket **does not spread appreciably**. This
  isolates the effect of the **quantum interactions** alone and excludes the
  interactions due to WP **spreading**. Need to design such a run.
