# A2 — launched-pair 100 eV kinetic-energy audit

Campaign: `docs/campaigns/localised_jellium_parameter_study_2/` (Energy book-keeping analysis), task A2.
Generated 2026-07-11. Verdict: **user's** (gate open).

## The pair (matched launched twins, qsp phase 3)

| | WP | classical |
|---|---|---|
| Run | `scripts/qsp_phase3/wp/results/p3_wp` | `scripts/qsp_phase3/classical/results/p3_classical` |
| Projectile | wavepacket σ_WP = 0.5, E = 100 eV, k₀ = 2.7110633401 | Gaussian-e ion σ_pot = 0.35 (= 0.5/√2 ✓), mass m_e, v₀ = 2.7110633401 |
| Cell / slab | 50×50×90 Bohr, spacing 0.5, half-width 12.5, N = 82, edge sharp, **3D-periodic** | idem |
| CAP | two-sided sin², η = −0.7 Ha, ±35..±45 | idem |
| Launch | z = −23.75 (r = 11.25 from face), dt = 0.04, 2500 steps | idem |
| GS | `scripts/qsp_phase3/gs` → `shared_gs/slab_n82_L50x50x90`, E_GS = −70.22568216820937 Ha | idem |

Note: "p3" in these run names = qsp **phase 3**, not periodicity 3 — but the box IS
3D-periodic. The A1 result (dKin, dXC periodicity-independent to < 0.05 eV) makes the
kinetic audit transferable to the locked-p2 context.

## Sub-claim (a): the classical projectile's 100 eV does NOT appear in E_total

1. `electron_track.csv` t=0: KE_ion = 3.674932 Ha = **100.00 eV** (= ½v₀², exact).
2. dE_CL(0) = E_total(0) − E_GS = 5.3987 Ha = **146.9 eV** — the same magnitude as the
   *at-rest* insertion energy at this r (h0_p3 ledger: 142.7 eV at r = 12, Lz = 120 box).
   No 100 eV on top.
3. Conservation (decisive): E_total(electronic) + KE_ion, step-matched:

| step | t (a.u.) | E_total (Ha) | KE_ion (Ha) | sum (Ha) |
|---:|---:|---:|---:|---:|
| 0 | 0.00 | −64.827 | 3.675 | −61.152 |
| 624 | 24.96 | −61.241 | 0.088 | −61.153 |
| 1248 | 49.92 | −66.327 | 5.148 | −61.179 |
| 2496 | 99.84 | −61.966 | 0.707 | −61.259 |

   The sum is constant to 2.9 eV over the run while each term individually swings by
   > 100 eV. If E_total contained KE_ion the sum could not be conserved.
   (The residual −2.9 eV drift is CAP absorption / numerics — untested split.)

## Sub-claim (b): the WP's 100 eV IS inside the electronic energies

- E_kin(WP, 0) − E_kin(CL, 0) = 9.42248 − 2.77807 = 6.64442 Ha = **180.8 eV**.
  (The classical E_kin(0) is the pure GS bath kinetic — the classical projectile adds
  no electrons, so the difference is the WP state's kinetic energy.)
- Independent cross-check (from the qsp phase-2 ledger work): ½Σ⟨p²⟩ of the WP at t=0
  = 6.644 Ha — identical.

## Sub-claim (c): dKin_WP−CL = KE_projectile + localisation energy

| quantity | Ha | eV |
|---|---:|---:|
| measured dKin(0) | 6.64442 | 180.8 |
| predicted: drift ½v₀² | 3.67493 | 100.0 |
| predicted: zero-point 3/(4σ²), σ = 0.5 | 3.00000 | 81.6 |
| predicted total | 6.67493 | 181.6 |
| **measured − predicted** | **−0.03051** | **−0.83** |

Agreement to 0.5%. (The −0.8 eV deficit is unattributed: candidate causes are
orthogonalisation of the WP state against the 82 bath states and finite-grid σ;
untested.)

## Incidental findings (flagged, not adjudicated)

1. **Mislabelled summary field**: `p3_classical/run_summary.txt` says
   `ke_ion_initial_ha = 0.669354436899`, but that value is ½·(final_vz)² — the FINAL
   ion KE. The true initial is 3.674932 Ha (track row 0). The summary writer records
   the post-run value under an "initial" name.
2. **The classical m_e projectile never entered the slab**: z went −23.75 → −14.7
   (face at −12.5); its KE oscillates (100 → 2.4 → 140 → 19 eV at the sampled steps)
   — it is bouncing in the attractive well outside the surface, consistent with the
   light-projectile deceleration rule.
3. WP−CL total at t=0 = 38.7 eV in this Lz=90 3D-periodic box (for context; the
   book-keeping decomposition of that number is the campaign's main thread).
