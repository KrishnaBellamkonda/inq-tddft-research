# Handover: energy-oscillation debugging (post-Arm-B thread)

Rolling handover for the thread that follows the cap_fix + pbc_open_z campaigns:
identify the CLOCK of the drain-then-rise E_total oscillation. Predecessors (both
done): `/local/data/public/skcb2/tddft/docs/handovers/cap-fix-experimentation.md`,
`/local/data/public/skcb2/tddft/docs/handovers/pbc-open-z-oscillation.md`.
Advisor spec-diff analysis (documented, committed f6a85a6):
`/local/data/public/skcb2/tddft/docs/notes/oscillating-vs-clean-run-spec-comparison.md`.

## 2026-07-14 (19:45) — m=1 engine-drift rerun of the clean p3_wp run LAUNCHED

**User request:** recreate the qsp_phase3 clean run (E_total decays to a fixed
value) exactly, m=1, "using the new library (inq-study) instead of inq".

**Premise correction (verified):** the original p3_wp run was ALREADY built
against inq-study — `CMakeCache.txt` in
`/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/qsp_phase3/wp/build/`
has `inq_SOURCE_DIR=/local/data/public/skcb2/tddft/inq-study`; `run_summary.txt`
says `engine = inq-study`. So the rerun is NOT an inq-vs-inq-study twin. It IS an
**engine-drift regression test**: six inq-study headers changed after the original
binary was built 2026-06-25 (`ks_hamiltonian.hpp`, `propagate.hpp`,
`electrons.hpp`, `laplacian.hpp`, `calculator.hpp`, `initial_guess.hpp` — the
mass-fork surgery). Question: does TODAY's inq-study at m=1 still reproduce the
clean decay?

**Config (exact original recipe):** `scripts/qsp_phase3/wp/run.cpp` unmodified;
`LJ_CAP=1 LJ_DT=0.04 LJ_N_STEPS=2500 LJ_LAUNCH_Z=-23.75 LJ_WRITE_EVERY=8
LJ_WF_EVERY=8`; GS `shared_gs/slab_n82_L50x50x90` (exists, untouched). Only
operational deltas: GPU 1 (GPU 0 occupied by another user), output
`results/p3_wp_m1_rerun` (original `results/p3_wp` preserved), and
`TMPDIR=$PWD/build/tmp` because the 9.8 GB `/tmp` partition is 100% full with
other users' files (killed the first build attempt with "No space left on
device"; do NOT delete others' /tmp files — always set TMPDIR to the big disk).

**Done + verified:**
- Rebuild against current inq-study: OK (full dep rebuild — libxc/catch2/spdlog).
- 50-step smoke on GPU 1: exit 0, all energies finite, propagation ended normally
  (`wp/smoke_m1_rerun.log`, results in `wp/results/smoke_m1_rerun/`).
- **KEY EARLY RESULT — bit-for-bit reproduction:** smoke energies match the
  original `prod_wp.log` to all 12 printed decimals at steps 1, 2, 10, 25, 50
  (e.g. step 50: e = -63.414464508084 both). The post-June-25 mass-fork engine
  edits are numerically inert at m=1. Expect the full run to reproduce the clean
  decay exactly.

**In flight:** production rerun, 2500 steps, launched ~19:45 on GPU 1,
log `wp/prod_m1_rerun.log`, output `wp/results/p3_wp_m1_rerun/`. ~4.2 s/step +
VTI-write spikes → ETA ~3.5–4.5 h (original wall was 5.7 h on GPU 0). No
checkpointing in this run.cpp (exactness prioritised over the checkpoint rule —
a kill loses the run; user aware).

**Not done:** final verdict (full-trace comparison vs original observables.csv),
notebook/doc update, commit of this milestone. S1/S2/S3 discriminators from the
advisor note remain designed-but-not-launched; S3 (absorbed-fraction analysis)
is zero-GPU and can run any time.

**Verify on completion:** exit 0 + `run_completed = true` in
`results/p3_wp_m1_rerun/run_summary.txt`; diff
`results/p3_wp_m1_rerun/raw/observables/observables.csv` against
`results/p3_wp/raw/observables/observables.csv` (expect identical or
noise-level); confirm E_total(t) decays to a fixed value with no late rise.

## 2026-07-14 (23:45) — RERUN COMPLETE, VERDICT: engine drift EXONERATED at m=1

Run finished clean: exit 0, `run_completed = true`, wall 14372 s (4.0 h, GPU 1).

**Primary result — BIT-IDENTICAL reproduction:** all 313 rows × all columns of
`observables.csv` match the original `p3_wp` exactly (max |dE_total| = 0.0).
Today's inq-study (post-mass-fork headers) at m=1 reproduces the 2026-06-25
clean run bit-for-bit over the full 2500 steps. The mass-fork engine edits are
numerically inert at m=1 — engine drift is ELIMINATED as a cause of the
oscillation. The oscillating campaign runs differ from the clean run by CONFIG
ONLY (advisor note suspects S1 CAP standoff / S2 perturbation strength).

**Secondary observation (present in BOTH runs, since bit-identical):** the
"clean" run's E_total minimum sits at t = 97.9 (of 100) with a terminal rise of
0.11 eV — marginally above the 0.1 eV noise floor. Under the period-lengthening
reading this hints even the clean box may turn, with a very long period; a
longer clean run would decide. Treat as a watch-item, not a claim.

**Shape numbers (2 s.f.):** drain 130 eV over the run; E_min at t = 98;
rise 0.11 eV; last-quarter slope −8.6e-2 eV/a.u. (still draining).

**Next:** S1/S2 discriminators (designed in the advisor note) and/or the
zero-GPU S3 absorbed-fraction analysis. Rerun outputs kept at
`wp/results/p3_wp_m1_rerun/`; smoke at `wp/results/smoke_m1_rerun/`.
