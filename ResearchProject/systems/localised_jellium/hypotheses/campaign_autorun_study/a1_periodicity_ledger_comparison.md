# A1 — periodicity 2 vs 3 energy-ledger comparison (insertion runs, t=0)

Campaign: `docs/campaigns/localised_jellium_parameter_study_2/` (Energy book-keeping analysis), task A1.
Generated 2026-07-11.

**Verdict (user, 2026-07-11): "use p2 for now"** — periodicity 2 locked for all
downstream tasks. Pragmatic lock: the p3−p2 WP-channel offset (4–7.5 eV, d(H+E)
only; dKin/dXC identical) stands as an open observation, not adjudicated.

## Provenance

| | periodicity 2 | periodicity 3 |
|---|---|---|
| Runs | `scripts/campaign_autorun/runs/h0_p2/{wp,cl}_r{4..40}_p2` | `scripts/campaign_autorun/runs/h0_p3/{wp,cl}_r{4..40}_p3` |
| GS | `runs/h2/gs_p2_lz120` (E_GS = 60.38307052445239 Ha) | `scripts/h0_base_difference/gs` == `runs/h2/gs_lz120` (E_GS = −108.5336851082701 Ha) |
| Slab | identical: 50×50×120 Bohr, half-width 12.5, edge_width 0 (sharp), N = 82, spacing 0.5, LDA | idem |
| Data row | `observables.csv` row 1 (t = 0 insertion, projectile at rest, σ_WP = 0.5, k0 = 0) | idem |

All energies in eV (1 Ha = 27.211 eV). Columns: dE_X = E_total(X, t=0) − E_GS; WP−CL = total difference;
dKin/dXC/d(H+E) = WP−CL per component, Hartree+external summed (G=0-robust per notebook cell 22).
Raw dHartree/dexternal individually are charged-cell-convention-poisoned (−274 eV p2 vs −29 eV p3 at r=40,
notebook cell 39) and are deliberately excluded.

## Periodicity 2 (E_GS = 1643.1 eV)

|   r | dE_WP | dE_CL | WP−CL |  dKin |   dXC | d(H+E) |
|----:|------:|------:|------:|------:|------:|-------:|
|   4 |  81.2 | 185.4 | −104.2 | 81.7 | −16.5 | −169.4 |
|  12 |  80.9 | 141.5 |  −60.6 | 81.7 | −16.5 | −125.8 |
|  20 |  80.6 |  97.5 |  −16.9 | 81.7 | −16.5 |  −82.2 |
|  28 |  80.2 |  54.6 |   25.6 | 81.7 | −16.5 |  −39.7 |
|  36 |  79.8 |  22.8 |   57.0 | 81.7 | −16.5 |   −8.3 |
|  40 |  79.6 |  12.3 |   67.3 | 81.7 | −16.5 |    2.0 |

## Periodicity 3 (E_GS = −2953.4 eV)

|   r | dE_WP | dE_CL | WP−CL |  dKin |   dXC | d(H+E) |
|----:|------:|------:|------:|------:|------:|-------:|
|   4 |  87.6 | 187.5 |  −99.9 | 81.8 | −16.5 | −165.2 |
|  12 |  87.0 | 142.7 |  −55.7 | 81.7 | −16.5 | −120.9 |
|  20 |  86.5 |  97.8 |  −11.3 | 81.7 | −16.5 |  −76.5 |
|  28 |  86.2 |  54.2 |   32.0 | 81.7 | −16.5 |  −33.3 |
|  36 |  86.0 |  21.9 |   64.1 | 81.7 | −16.5 |   −1.2 |
|  40 |  85.9 |  11.2 |   74.7 | 81.7 | −16.5 |    9.5 |

## p3 − p2 (row-wise, eV)

|   r | dE_WP | dE_CL | WP−CL | dKin |  dXC | d(H+E) |
|----:|------:|------:|------:|-----:|-----:|-------:|
|   4 |   6.4 |   2.1 |   4.2 |  0.0 |  0.0 |    4.2 |
|  12 |   6.1 |   1.2 |   4.9 |  0.0 |  0.0 |    4.9 |
|  20 |   6.0 |   0.3 |   5.7 |  0.0 |  0.0 |    5.7 |
|  28 |   6.0 |  −0.4 |   6.4 |  0.0 |  0.0 |    6.4 |
|  36 |   6.2 |  −0.9 |   7.1 |  0.0 |  0.0 |    7.1 |
|  40 |   6.3 |  −1.2 |   7.5 |  0.0 |  0.0 |    7.5 |

Max |p3 − p2|: 7.5 eV (WP−CL and d(H+E) at r = 40).

## Neutral observations (not a verdict)

- dKin and dXC — the WP quantum self-energy pieces — are **identical across periodicities to < 0.05 eV**.
- dE_CL agrees to ≤ 2.1 eV (within the campaign's 3 eV row tolerance).
- dE_WP differs by a near-constant **+6.0 to +6.4 eV** (p3 higher); the entire offset sits in the
  electrostatic d(H+E) channel and grows mildly with r (4.2 → 7.5 eV).
- Inference (unverified): the WP insertion makes the cell net −1 charged, and the 2D-periodic (p2)
  vs 3D-periodic (p3) Poisson G=0 conventions treat that net charge differently — consistent with the
  charged-cell caveat of notebook cell 39. The classical run stores the projectile in the *external*
  potential (cell stays 82 e⁻), which would explain why dE_CL is nearly convention-free.
