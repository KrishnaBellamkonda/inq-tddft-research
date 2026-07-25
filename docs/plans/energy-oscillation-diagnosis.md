# Plan: Localised-jellium ΔE_total energy-oscillation diagnosis

**Status:** design locked (grill-with-docs, 2026-07-13). Ready to author campaign +
scaffolds, then launch a background agent.

**Phenomenon** (`docs/notes/localised-jellium-energy-oscillation-investigation.md`,
glossary term in `CONTEXT.md`): in many localised-jellium RT runs `ΔE_total(t) =
E_total(t) − E_ref` does **not** decay monotonically once the CAP begins absorbing —
it *oscillates*, and in several runs rises **above 0**, which is unphysical (a closed
system has no energy influx; a CAP can only *remove* energy). Seen across WP,
effective-mass, heavier-electron, and some truncated-classical runs; contrast `p3_wp`
where ΔE_total decays to a stable plateau as expected.

**Goal (locked):** *diagnose + document only* — isolate the cause to **one confirmed
mechanism** via a decisive control experiment, documented in a notebook. Diagnostic
instrumentation (turning on extra energy-component output) is allowed; **no physics
fix is committed** in this campaign.

---

## Architecture (agent + advisor autonomous loop)

Locked via grill-with-docs. Vehicle = a **campaign** (ADR 0009) executed by a fresh
**background agent** that runs an adaptive loop:

```
                ┌─────────────────────────────────────────────┐
                │  INVESTIGATOR  (executing background agent)  │
                │  - runs ONE tiny probe (mine existing, or    │
                │    build+run an ablation run.cpp variant)    │
                │  - extracts energy components, plots          │
                │  - reports RAW results faithfully            │
                └───────────────┬─────────────────────────────┘
                                │ result.json + plots
                                ▼
                ┌─────────────────────────────────────────────┐
                │  ADVISOR  (single persistent methodologist,  │
                │  spawned via Agent tool each iteration)      │
                │  - reviews the probe against the ledger      │
                │  - rules candidate causes in/out             │
                │  - names the NEXT decisive probe  (BINDING)  │
                │  - or declares a cause CONFIRMED → stop       │
                └───────────────┬─────────────────────────────┘
                                │ ledger update + next probe
                                ▼
             loop until: cause confirmed  OR  ~8 experiments  OR  budget
```

- **Investigator** = the background agent (main loop). Runs the probe, reports raw
  numbers, never adjudicates.
- **Advisor** = one persistent subagent with a **TDDFT / stopping-power
  methodologist** lens (the E1 charter from `scientific-panel`, reused as a single
  standing critic). Verdict is **binding**: the investigator runs whatever probe the
  advisor names next. The advisor owns `hypothesis_ledger.md`.
- **Stop:** advisor declares one mechanism confirmed by a decisive control, OR the
  hard cap of ~8 tiny experiments, OR a wall-clock/token budget — whichever first.

## Candidate mechanisms (the ledger discriminates these)

| # | Hypothesis | Decisive control |
|---|---|---|
| a | CAP is a non-Hermitian energy **source** in the reported ledger (E_total not the conserved quantity under the absorber) | CAP-off run: does the ΔE>0 rise vanish? η-sweep: does amplitude scale with η? |
| b | Static background `v_bg` contribution is **absent from the reported energy functional** (E_ext/E_ion mis-accounts the perturbation) | `+v_bg` only (no projectile, no CAP): is E_total conserved? Compare to plain whole-cell jellium |
| c | Wrong subtracted **`E_ref`** (ΔE vs E_GS instead of E_total(0) of the RT run; charged-cell G=0 convention) | Recompute ΔE against E_total(0); inspect component baselines |
| d | **Propagator / grid numerics** (ETRS vs CN energy drift, dt too large, cutoff aliasing) | Pure-GS propagation conservation floor; dt-halving; cutoff guard |
| e | **Density-dependent (time-dependent) KS Hamiltonian** double-counting | Component decomposition: which term drifts; `energy_eigenvalues` (Σε_i) vs `energy_total`; `energy_nvxc` |

## Experiment strategy (locked: mine existing → ablate new)

**Phase 0 — mine existing (free):** read `observables.csv` of the named runs
(`hypotheses/muon_mass_fork/effmass_sigma1_*`, `p3_wp`, a truncated classical). Check
which energy-component columns exist. Localise the drifting component if columns are
present. **Expected gap:** defaults record only `energy_total`+`energy_kinetic`, so
components are likely absent → triggers Phase 1.

**Phase 1 — ablation ladder (new tiny probes).** Smallest cell that reproduces the
ΔE>0 rise (L≈30–50 Bohr, dx≈0.4 respecting the cutoff guard, ~100–300 steps, single
WP electron, **all energy-component ObservableSelection flags ON**). Strip one
subsystem at a time; each isolates a culprit:

1. **pure-GS propagation** — no projectile, no CAP, no kick. Conservation floor.
2. **+ `v_bg` background only** — localised jellium, no projectile, no CAP.
3. **+ CAP, no projectile** — CAP absorbs GS tail only.
4. **+ projectile (WP), no CAP** — projectile-induced density change only.
5. **single wrap-around CAP vs double CAP** (note Q1) — geometry effect.
6. **component decomposition** on whichever run drifts (note Q4): which of
   E_kin/E_H/E_ext/E_xc/E_ion carries the rise; `energy_eigenvalues`, `energy_nvxc`.

Advisor reorders adaptively and may add: occupation-number tracking (note Q3),
dt-halving, or a "does absorbed density's energy linger" probe (note Q2).

## Output contract (locked)

Folder `ResearchProject/systems/localised_jellium/hypotheses/energy_oscillation_diagnosis/`:

- **`<probe>_run_notebook.ipynb`** — a full standalone run-notebook per probe (via the
  `run-notebook` skill).
- **`energy_oscillation_diagnosis.ipynb`** — the master study notebook: one section
  **per experiment** with **Aim → Method → What was plotted → Results → advisor
  verdict**, plus intro (phenomenon + ledger mirror) and final synthesis (confirmed
  mechanism).
- **`hypothesis_ledger.md`** — living doc; advisor updates every iteration.
- **`probes/<probe>/`** — raw outputs (observables.csv, result.json, PNGs) per probe.

Run machinery in `scripts/energy_oscillation_diagnosis/`: `run_probe.py`
(build+run+extract one probe), `build_master_notebook.py` (assemble the master
notebook), ablation `run.cpp` variants (generated adaptively).

## Safety (always-on rules, restated)

- Correctness-only hard gates: cutoff/aliasing guard pre-launch, abort on NaN /
  complex energy, GS-present check. **No cost gate** — `checkpoint-don't-block`.
- GPU is default (NVML mismatch is not a blocker; verify with `cudaMemGetInfo`
  probe; warn if a GPU is occupied by another user).
- Per-iteration Gmail (hypothesis reminder → what was done → what the plot shows →
  conclusion; ≥1 plot attached).
- Every new run.cpp variant that adds logic ships its validation status; the
  observables writer flags are config, not new numerics (no new formula gate needed).

## Deliverables of this planning pass

1. This plan.
2. `docs/campaigns/localised_jellium/energy-oscillation-diagnosis.md` — self-contained
   campaign with embedded advisor charter + `<preflight>`.
3. `hypotheses/energy_oscillation_diagnosis/hypothesis_ledger.md` — seeded.
4. `scripts/energy_oscillation_diagnosis/run_probe.py` + `build_master_notebook.py`
   — working skeletons.
5. Handover `docs/handovers/energy-oscillation-diagnosis.md`.

Then: validate scaffolds, final user confirm, launch background agent.
