# Handover: bulk-jellium KS-orbital stopping power (classical + wavepacket twin)

**Rolling file. Latest milestone at top.**
**Repo:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research` (branch `quantum-stopping-power`)
**Machine:** CSD3, `ampere` partition (A100), account `mphil-nikiforakis-skcb2-sl2-gpu`
**Plan:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/bulk-jellium-ks-stopping.md`

---

## 2026-08-04 (latest) — RUN RECORD written; report-2 bulk figures deleted and rebuilt (4 figures)

### 1. Complete run record for the write-up

`ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/BULK_JELLIUM_RUN_RECORD.txt`
— a plain-text, paper-facing account of the whole campaign: the system and both
baths (r_s = 3.987 / 5.702 with their derived E_F, ω_p, v/v_F), numerics
(box, dx, cutoff, dt, states, SCF tol), both ground states, the two projectile
representations incl. the UPF provenance, the 12-run table, the analysis-window
rule, a file-by-file schema of everything each run wrote, the four S definitions,
the extracted results, the validation gates, and the limitations to state in a
paper. Every number in it was read back from the run artefacts; nothing retyped.

**Two facts it pins down that had been drifting:** (a) this campaign is **ONE
velocity** (100 eV) at both densities — the axes swept are density and σ, never
velocity; (b) the older r_s = 5.69 / 3.41 energy sweeps are *different systems*
(different box, N, σ, S-extraction) and must not be merged with it.

**Verified from the UPF meshes while writing it** (new, not previously recorded):
all three classical projectile files are exactly `V(r) = 2·erf(r/(√2 σ_pot))/r`
Ry to ≤ 2.2e-5, i.e. a Gaussian charge of std σ_pot = σ_WP/√2 — the σ-matching
convention is now confirmed from the potential itself, not only from E_PP(0).

### 2. Report-2 bulk figures: all deleted, four rebuilt

User instruction: remove every bulk plot in
`docs/reports/report2/drafts/draft1/figures/bulk_jellium/` and remake carefully.
**24 PNGs deleted** (that folder is gitignored, so this was not git-recoverable —
a full copy is in the session scratchpad `bulk_jellium_archive_2026-08-04/`, and
the 12 old generators are in `figures/bulk_jellium/_superseded/`).

New: one shared loader `_bulk_data.py` + four scripts. **S is now fitted over a
PATH window s = 10–50 Bohr** (user decision) instead of a time window, so the
decelerating classical projectile and the packet are compared over the same
40 Bohr of travel; each script prints the run-window value too.

| figure | headline numbers |
|---|---|
| `bulk_total_energy.png` | classical ledger closes to 0.52 eV = 0.96 %; WP total drifts 3.0e-4 eV |
| `bulk_dE_vs_position.png` | S = **0.906** (classical) / **0.156** (T_drift) / **−0.043** (T_int) / **0.114** (T_orb) eV/Bohr; ratio 5.8; additivity exact to 3e-17 |
| `bulk_interaction_energies.png` + `.csv` | ΔE_SS +26.13 cl vs +5.33 WP; ΔE_PS −5.47 vs −0.39; ΔE_PP 0.000 vs −4.62 |
| `bulk_S_vs_sigma.png` + `.csv` | classical S falls 1.083→0.819 with σ while WP S rises 0.094→0.175; ratio 11.6→4.7 |

All of these reproduce the repo's own `*_stopping_summary.json` when refitted on
the run's own window (0.885 / 0.157 / 0.115), so the path window is a change of
window, not of method.

**NEW FINDING — a closure gate in the rules does NOT hold.**
`.claude/rules/decomposed-interaction-energies.md` asserts
`E_SB + E_PS == energy_external` for a classical projectile. In
`bulk_ks_stopping_rs4/classical` it fails: `energy_external` = +58.97 Ha at t = 0
against E_PS = +0.0021 Ha, and the CHANGES differ too (−10.19 eV vs −5.47 eV,
ratio drifting 1.46 → 1.86 along the trajectory). The **Hartree** gates hold to
1.4e-11 eV in both halves, so the offline pairwise terms are sound; it is
`energy_external` that is unusable here. *Inference (labelled, untested):* INQ
evaluates the bare erf-Coulomb UPF in real space and retains a G = 0 /
minimum-image contribution that the offline Poisson solve drops, and it is
position-dependent because the UPF's unscreened tail (r_max = 50 Bohr) exceeds
half the 40 Bohr transverse cell. **To do:** either amend the rule with this
exception or test it (e.g. re-measure with a shorter-ranged UPF / larger cell).

Also corrected on the record: **E_PP is not "missing" in the classical half** — it
is present and identical at t = 0 (4.8163 eV both halves) and simply cannot
evolve, the cloud being rigid. And the electrostatic totals of the two halves do
NOT tally (+25.54 vs +5.21 eV); that 20.3 eV gap is the result, not an error.

Full per-figure detail, including the annotation-placement trap hit again:
`docs/reports/report2/drafts/draft1/plots_draft1_log.md` (top entry).

---

## 2026-08-02 — 100 eV HIGH-DENSITY CASE STUDY (report-ready figure set)

### What was asked

A case study on the **high-density bulk run at 100 eV**: classical KE(t); the
wavepacket's total T₁, the ⟨p⟩²/2m term and the var(p) term; position vs time
(WP centroid from **integrating ⟨p_z⟩**, i.e. Ehrenfest, not the density
centroid); T₁/T₂ vs position; the stopping power over the transient-excluded
window shown **as overlays on duplicated versions** of those plots; a text file
with S and its uncertainty; interaction-energy plots per projectile type plus a
difference plot. All in a downloadable folder. Mid-task: **"the plots must be
report ready."**

### Family selected — state this if it is revisited

All four pairs are 100 eV, so "100 eV" did not disambiguate. "High density" =
**r_s = 3.99** (`bulk_ks_stopping_rs4`), and **σ_WP = 2** was taken as the
reference σ. The builder is parameterised: `--family bulk_ks_stopping_rs4_sigma3`
produces the σ = 3 twin with no code change.

### Deliverables — all complete and verified

`ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping_rs4/case_study_100eV/`

| file | contents |
|---|---|
| `make_case_study.py` | the builder; every constant parsed from the run's own `run_summary.txt` / `wp_config.txt`, never retyped |
| `01`–`16` `*.png` | 16 figures, 600 dpi, canonical theme (ADR 0004), fixed canvas |
| `stopping_power.txt` | S ± uncertainty, energy budget, interaction-energy table, validation block, provenance |
| `case_study.md` | the written case study + production log (report-figures skill) |

### Results

| quantity | S (eV/Bohr) | r² |
|---|---|---|
| classical ½mv² | **0.88 ± 0.13** | 0.991 |
| wavepacket T₂ = ⟨p⟩²/2m (drift) | **0.16 ± 0.01** | 0.997 |
| wavepacket T₁ = ⟨p²⟩/2m (total) | 0.12 ± 0.02 | 0.993 |
| wavepacket T_var = var(p)/2m | −0.04 ± 0.02 | 0.902 |

Ratio classical/WP(T₂) = **5.6 ± 0.9**. Window t = 4.0–18.4 a.u. (from the run,
not chosen here). The **window systematic dominates**; the statistical error
rounds to 0.00 in every row.

Full-run budget: classical −54.4 eV, WP drift −8.7 eV, of which **+3.0 eV goes
into the packet's own momentum spread** and never reaches the bath.
Interaction-energy change: E_SS +26.1 (cl) vs +5.3 eV (WP); E_PP 0.000 (cl,
exactly constant) vs −4.6 eV (WP).

### Two findings worth carrying forward

1. **E_PP does not by itself explain the ~5.6× ratio.** The classical-minus-WP
   drift-loss gap is 45.7 eV; the packet's E_PP releases 4.6 eV, ~10% of it. The
   dominant single term is the **wake**: ΔE_SS 26.1 vs 5.3 eV, a 4.9× ratio that
   tracks the 5.6× stopping ratio. This *partially deflates* the hypothesis in
   `.claude/rules/decomposed-interaction-energies.md` that E_PP is the leading
   candidate for the bulk factor. Labelled inference; no ledger closure was
   attempted (bath kinetic + XC are outside this decomposition), so the 20.8 eV
   and 45.7 eV gaps must not be differenced naively.
2. **A number in the 2026-08-01 entry was wrong and is corrected here.** The WP
   zero-point momentum channel at t = 0 is **T₁ − T₂ = 5.10 eV**, not 2.55 eV —
   it is 3/(4σ_ψ²) with σ_ψ = 2 (the ψ-width convention these runs use), and the
   data confirms it to 6 digits. The earlier value used the density-width
   convention. Now asserted in the test suite so it cannot drift again.

### Validation

Closure vs INQ Hartree 4.98e-13 / 4.99e-13 Ha · Ehrenfest residual
(⟨z⟩ vs ∫⟨p⟩dt) **0.079 Bohr over 68.5 Bohr** · WP norm 0.999989595–1.0 ·
classical cloud never clips (window ends 18.4 a.u.) · S additivity
S(T₁) = S(T₂) + S(T_var) to 1e-10 · **E_PP(0) identical across halves**
(0.176996 Ha), which validates the σ_pot = σ_WP/√2 convention from the data ·
margin check: all 16 figures ≥ 8 px clear.

**Tests: 21 new in `hypotheses/bulk_ks_stopping/tests/test_case_study.py`;
suite now 123 passing (was 102).**

### Gotchas that cost time — read before touching the figure code

- **`hypotheses/<family>/tests` and `<family>/ks_stopping.py` are SYMLINKS** to
  `bulk_ks_stopping/`. `Path(__file__).resolve()` therefore always lands in the
  canonical suite, so a test must anchor on the hypotheses root and name its
  family explicitly — walking up from `__file__` silently targets the wrong
  family. (Cost one confusing collection error; the two "duplicate" test copies
  are literally the same file.)
- **Fixed-canvas figures crop silently.** No `bbox_inches="tight"` means an
  overrunning title is simply absent from the PNG with the build reporting
  success. `verify_margins()` now measures the ink bbox of every saved PNG and
  fails the build; `_shrink_overrunning_titles()` handles horizontal overrun.
- **numpy 2 `repr(np.float64(x))` is `"np.float64(0.3)"`** — writing synthetic
  test CSVs with f-string `!r` produces string columns and a TypeError. Use
  `DataFrame.to_csv`.
- Mathtext spans are parser-checked in the suite with a negative self-test; the
  `$\frac12$` class of bug is valid Python *and* valid TeX.

### 2026-08-02 follow-up — all six interaction terms

User asked for **all** the interaction energies plotted, not just the three that
move. Figures 14/15/16 now carry E_SS, E_PP, E_PS, E_SB, E_PB, E_BB; the results
file lists all six. E_SB/E_PB/E_BB are **bitwise exactly 0.0 at every step in
both halves** (verified, now asserted in the suite) because bulk's uniform
background makes phi_+ identically zero. 14/15 gained a lower panel magnified
50x showing them still flat -- that is what separates "structurally zero" from
"not computed". The three figures moved to TWO-COLUMN width: a six-entry legend
at one-column width is either half the panel tall or wider than the canvas, and
the margin check rejected both. The `legend()` helper now MEASURES the legend and
sizes headroom from it, rather than the fixed 1.30x that was silently calibrated
for three entries. Suite: **125 passing**.

### Not done (needs the user's call)

σ = 1 families still lack `interactions.csv` (GPU time) · the same case study for
the other three pairs (one command each) · `energy_component_comparison_*.ipynb`
still carry no density GIF (`.claude/rules/notebook-density-gif.md`).

---

## 2026-08-01 — NOTEBOOKS REBUILT + NEW PHASE-SPACE PAIR NOTEBOOKS

### What was asked

Rebuild the run notebooks so they carry the new interaction-energy graphs, and
add a **phase-space notebook per classical/wavepacket pair** to compare the two
representations directly. (User chose "phase-space" from three readings of
"phase notebook"; the repo's other meaning — `qsp_phase1..5`, an
investigation-STAGE notebook spanning a run-set — was explicitly not what was
wanted.)

### Status

`shared/bin/run-bulk-notebooks.slurm` — 12 tasks, 4 concurrent, 4 h wall each.
Tasks 0–7 rebuild the run notebooks (4 families x wp/classical); tasks 8–11
build the new phase notebooks.

| Submission | Tasks | Outcome |
|---|---|---|
| `32574476` | 8–11 (phase) | **COMPLETE** — 0 errors, 14 inline images each (9 density GIFs + 5 figures), 58–67 MB |
| `32574476` | 0–7 (run) | **FAILED**, 1 error each — the eaten-`\n` bug below |
| `32575606` | 0–7 (run) | **COMPLETE** after the fix — 0 errors, 35–52 MB, WP 33 cells/16 figures, classical 22 cells/7 figures |

**All 12 notebooks verified from their stored outputs (not just exit codes):
0 execution errors, and every closure gate PASSED.**

| run | closure vs INQ E_H | E_PP constant | clipping onset vs fit window |
|---|---|---|---|
| `bulk_ks_stopping` wp / cl | 4.86e-13 / 4.99e-13 Ha | 5.1e-11 Ha (305/324 rows) | t=24.56 vs 18.97 — CLEAR |
| `bulk_ks_stopping_sigma3` wp / cl | 4.98e-13 / 4.99e-13 Ha | 1.8e-11 Ha (260/301) | t=21.04 vs 19.48 — CLEAR |
| `bulk_ks_stopping_rs4` wp / cl | 4.98e-13 / 4.99e-13 Ha | 5.7e-11 Ha (323/324) | never clips |
| `bulk_ks_stopping_rs4_sigma3` wp / cl | 5.00e-13 / 4.96e-13 Ha | 1.6e-11 Ha (276/301) | t=22.32 vs 19.48 — CLEAR |

`norm_wp` stays in 0.999989 .. 1.000000 on all four WP halves (no CAP in this
study, so it must).

Products, per family (`ResearchProject/systems/jellium/hypotheses/<family>/`):
`<family>_wp.ipynb`, `<family>_classical.ipynb`, `run_pair_phase_<family>.ipynb`,
plus `phase_gifs_<family>/` and `phase_{portrait,velocity,kinetic,local_S,epp}_<family>.png`.

### 1. New kernels in `ks_stopping.py` (+12 tests, 23 total, all passing)

`/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/ks_stopping.py`

| Added | Does |
|---|---|
| `Interactions` + `load_interactions(run_dir, half)` | loads `interactions.csv` (+ resume segments), computes the closure residual against INQ's own `energy_hartree` |
| `Interactions.clip_index` / `.clip_time` / `.in_window()` | projectile-cloud clipping onset as a **trailing contiguous run**, and window truncation at it |
| `local_stopping(s, T, half_width)` | `S(s) = -dT/ds` as a **centred rolling OLS slope**, eV/Bohr |
| `PairPhase` + `load_pair(...)` | both halves on one footing; `.divergence(frac)`, `.epp_on_z(half)` |

Two things deliberately encoded, because getting either wrong silently corrupts
the comparison:

* **`T2`, not `T1`, is the classical analogue.** The WP's total orbital KE
  `T1 = <p^2>/2m` includes its own momentum spread; `T2 = <p>^2/2m` is the drift
  energy. Comparing `T_cl` against `T1` charges the packet ~2.6 eV
  (`3/(8 sigma^2)` at sigma=2) before it has moved. `T1-T2` is reported
  separately as the quantum-only channel.
* **`local_stopping` is a rolling OLS, not a finite difference.** Per-step dT is
  ~1e-5 Ha against ds ~0.1 Bohr, so `np.gradient` is pure noise. Edges are
  FILLED with the nearest interior value (documented + tested) — extrapolating
  them would manufacture a fake "S spikes at impact" feature.

### 2. Interaction-energy section added to all four run-notebook builders

The four `hypotheses/bulk_ks_stopping*/build_run_notebook.py` are clones
differing only in docstring, `CFG` literals and `run_dir`, so ONE insertion was
applied to all four by
`<scratchpad>/patch_ie_section.py` (idempotent; refuses a file already carrying
the marker). Verified: the inter-family diff is still exactly 41 lines, i.e. the
inserted block is byte-identical everywhere.

New **section 7 — Pairwise interaction energies**; the old takeaway became
section 8. Contains the closure gate vs INQ `energy_hartree`, the E_SB/E_PB/E_BB
== 0 bulk schema check, E_SS/E_PP/E_PS relative curves, the norm/clipping
diagnostic, and (classical) the E_PP-constancy gate with the clipping onset
compared against `FIT_T1`.

### 3. NEW phase-space pair notebooks

Builder: `hypotheses/bulk_ks_stopping/build_phase_notebook.py`
(`--all`, or one family; `--no-execute` to write without running).
Product: `hypotheses/<family>/run_pair_phase_<family>.ipynb`, 19 cells.

Sections: density-matrix GIF (rule-mandated, `make_twin_density_matrix`, 3x3,
30 frames) → **(z, v) phase portrait** with equal-time markers and the
divergence annotation → velocity histories + `dv` → KE(z) and the internal
`T1-T2` channel → **local S(z)** → **E_PP on the same z axis** + gauge-clean
`dE_PP` → takeaway.

`CFG` is **imported** from each family's own `build_run_notebook.py` rather than
retyped, so the run and phase notebooks for a family cannot disagree about box,
launch point or fit window.

### 4. Verified before submission (all four families, real data)

`<scratchpad>/validate_phase_cells.py` ran every analysis cell headless:

| family | sigma | r_s | v/v0 @ win-end (cl / wp) | diverge | mean S cl / wp | dE_PP |
|---|---|---|---|---|---|---|
| `bulk_ks_stopping` | 2 | 5.702 | 0.9195 / 0.9881 | t=15.4, z=+8.7 | 0.363 / 0.054 | −4.42 eV |
| `bulk_ks_stopping_sigma3` | 3 | 5.702 | 0.9250 / 0.9871 | t=16.9, z=+16.7 | 0.330 / 0.059 | −2.36 eV |
| `bulk_ks_stopping_rs4` | 2 | 3.987 | 0.8100 / 0.9671 | t=9.2, z=−7.6 | 0.875 / 0.151 | −4.41 eV |
| `bulk_ks_stopping_rs4_sigma3` | 3 | 3.987 | 0.8140 / 0.9630 | t=10.1, z=−1.4 | 0.808 / 0.165 | −2.38 eV |

All gates PASS (sigma-matching `dE_PP(0)`, rigid-cloud constancy) and every fit
window is CLEAR of its clipping onset. `dE_PP` reproduces
`summarise_epp_across_pairs.py` to ~0.02 eV (different window grid).

**The local estimator is cross-checked against the established windowed fits**,
now printed inside the notebook: classical local 0.363 vs global `S_cl_shared`
0.377; WP local 0.054 vs global `S_24` 0.057 — a few percent. `S_24` is the
right partner because it is the fit built on the DRIFT KE `T2`.

### 5. Three failed attempts worth recording (all now guarded by tests)

**THE HEADLINE LESSON.** The same bug class — *a builder that writes Python
source as strings, with a non-raw string containing a backslash* — bit **three
times in one day across two systems**. The escape fires when the BUILDER is
parsed, so the emitted cell is already corrupt:

| in the builder | becomes in the cell | symptom |
|---|---|---|
| `\r` (`\rangle`, `\rm`) | carriage return | SyntaxError (CR is a line terminator) |
| `\n` (in an f-string annotation) | real newline | unterminated string literal |
| `\a` (`\approx`), `\f` (`\frac`), `\v` (`\varphi`), `\b` (`\beta`) | BEL / FF / VT / BS | **compiles fine, renders wrong** |

Any new notebook builder in this repo should ship
`tests/test_notebook_cells.py` from the start, and use `code(r"""...""")` by
default.

**(a) Non-raw builder strings ate the LaTeX.** The first executed phase notebook
came back with **5 errors**, all one root cause: cell sources are emitted from
triple-quoted strings inside `build_phase_notebook.py`, and in a **non-raw**
builder string LaTeX like `\rangle` / `\rm` has its `\r` eaten as a carriage
return — so the emitted cell was a `SyntaxError` before it ever ran. (The 5th
error was just the cascade: `i1` is defined in one of the broken cells.) Fixed by
making those cells raw. **If you add a cell to this builder, use
`code(r"""...""")`.**

This is the SAME bug that hit
`cylindrical_jellium/hypotheses/channeling_twin/build_comparison_notebook.py` on
the same day — it is a property of the write-source-as-strings pattern, not of
either builder, so any new notebook builder in this repo should ship the guard
test below from the start.

**(b) `\frac12` is valid TeX but not valid mathtext.** `$\frac12mv^2$` compiles
as Python, passes a raw-string check, passes a blocklist of known-bad commands,
and still dies at RENDER time — matplotlib's mathtext requires
`\frac{num}{den}`. Caught **while the array job was already queued**, before the
phase tasks started, by the render test described next. Fixed to
`$\frac{1}{2}mv^2$`.

**(c) The same bug in the section I had just inserted — and in the ONE builder
the first guard did not cover.** All 8 run-notebook array tasks failed with 1
error each: the interactions figure cell's
`annotate(f"cloud clips the +z face\nt = …")` was emitted from a non-raw builder
string, so `\n` became a real newline and the cell was an unterminated literal.
Fixed by making that cell raw in all four builders (anchored on the 4-space
indent — three PRE-EXISTING cells share its first line at 8-space indent).

**Guard test — now covers EVERY builder, 12 notebooks (8 run + 4 phase):**
`hypotheses/bulk_ks_stopping/tests/test_notebook_cells.py`
(supersedes the phase-only `test_phase_notebook_cells.py`, deleted).

| Case | Catches |
|---|---|
| `test_every_code_cell_compiles` | (a) and (c) — CR and LF both terminate a line, so both land here |
| `test_no_control_characters_in_string_literals` | the SILENT variants — `\a \b \f \v` inside a cell literal, via the AST (runtime values) |
| `test_every_math_span_actually_renders` | (b) — every `$...$` span parsed by **matplotlib's own MathTextParser** |
| `test_no_unsupported_mathtext_commands` / `test_math_spans_are_balanced` | the enumerated `\le`/`\ge`/`\text{` habits; unterminated `$` |
| `test_run_notebook_has_the_interaction_energy_section` | the decomposed-energies rule — load, heading, closure gate, all three terms plotted |
| `test_phase_notebook_displays_the_density_gif` | the density-GIF rule — produced AND `display()`ed |
| 6 negative self-tests + 1 positive | each guard **provably fires** on the real defect, and correct cells do not trip |

**102 tests pass** (23 engine + 79 notebook-guard/parametrised).

**Two lessons, both about guard quality.**

1. **Prefer a real oracle to a list of remembered mistakes.** `\frac12` is valid
   Python AND valid TeX and only fails inside matplotlib's subset, so neither
   `compile()` nor a blocklist could catch it. Asking matplotlib to parse every
   span did.
2. **Test the emitted artefact at the right level, and prove the guard fires.**
   The first version flagged `"$T_1=\\langle …$"` (a correctly-escaped non-raw
   literal) and a genuine `"\n"` — 20 false positives — because it read raw
   SOURCE text. The cell-level signature of the bug is a **control character in
   a literal**, read from the AST. The negative self-tests exist because the
   first version of the CR self-test asserted the WRONG guard and passed
   vacuously until it was made to fail.

### Not done / not approved

* `interactions.csv` retrofit to the sigma=1 families (excluded by the user).
* The `energy_component_comparison_*.ipynb` notebooks do **not** carry a density
  GIF — a pre-existing gap against `.claude/rules/notebook-density-gif.md`, not
  introduced here and not in scope of this request.
* `interactions.csv` is buffered and only materialises at clean exit, so a
  KILLED run loses its rows. Harmless for completed runs; add `ix.flush()` in
  the callback if those `run.cpp` are touched again.

---

## 2026-08-01 — ALL REMAINING TWIN PAIRS RE-RUN WITH INTERACTION ENERGIES

### What was asked

Make classical twins of the recent bulk jellium sigma runs so classical and WP
can be compared, **ignoring sigma = 1**; run sigma = 2 and 3 where a WP analogue
exists; ensure interaction/decomposed energies AND the WP energetics are saved;
catalogue first, then build and submit autonomously.

### Catalogue (measured from disk, not assumed)

Every family ALREADY had both halves on disk — what was missing was
`interactions.csv`. Only `bulk_ks_stopping_sigma3` had it.

| family | sigma | r_s | N | box | dx | steps | status |
|---|---|---|---|---|---|---|---|
| `bulk_ks_stopping_sigma3` | 3 | 5.702 | 218 | 46x46x80 | 0.40 | 600 | DONE (302 rows both halves) |
| `bulk_ks_stopping` | 2 | 5.702 | 218 | 46x46x80 | 0.40 | 646 | **re-run** (tasks 0,1) |
| `bulk_ks_stopping_rs4` | 2 | 3.987 | 482 | 40x40x80 | 0.50 | 646 | **re-run** (tasks 2,3) |
| `bulk_ks_stopping_rs4_sigma3` | 3 | 3.987 | 482 | 40x40x80 | 0.50 | 600 | **re-run** (tasks 4,5) |
| `bulk_ks_stopping_sigma1` | 1 | 5.702 | — | — | — | — | EXCLUDED (user) |
| `bulk_ks_stopping_rs4_sigma1` | 1 | 3.987 | — | — | — | — | EXCLUDED (user); WP half was also incomplete, 88/348 rows |

**BOTH halves are re-run, not just the classical twins.** The quantity that
answers the question is the gauge-clean `dE_PP = E_PP(WP) - E_PP(classical)` in
the SAME cell; absolute E_PP carries the charged-cell G=0 gauge. A pair is
useless unless both halves carry `interactions.csv`.

### Wiring

6 x `run.cpp` patched by cloning the reference wiring from
`scripts/bulk_ks_stopping_sigma3/{wp,classical}` — 3 additive hunks each
(include, setup, callback). Verified by diffing each patched file against the
reference: the ONLY residual differences are header comments, `RUN_NAME`, the GS
path, and the config struct name. All executable code is identical.

Patch script (throwaway, not repo state):
`<scratchpad>/wire_ie.py`. Pre-re-run observables CSVs archived to
`<scratchpad>/preIE_csv_backup` (41 files, 61 KB) as insurance.

**Physics is unchanged by the wiring** — `compute_coulomb*` only reads the
density and solves Poisson, never touching the Hamiltonian. The re-runs
reproduce the existing trajectories and strictly ADD columns.

### NEW: closure verifier (validated before use)

`/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/ResearchProject/systems/jellium/scripts/verify_interactions_closure.py`

Implements the mandatory closure gates from
`.claude/rules/decomposed-interaction-energies.md`. Run per half; exits non-zero
on failure; invoked automatically by every array task.

Validated against the known-good sigma=3 pair BEFORE being wired in — and it
**caught a real effect on the first run**:

- **E_PP is not constant over the whole classical run** (spread 6.0e-4 Ha).
  Cause is NOT egg-box error (correlation with sub-grid phase = 0.03). It is the
  projectile's Gaussian cloud being **clipped by the +z box face**: `norm_proj`
  falls 1.000000 -> 0.994247 over the final 32 of 301 rows.
- Restricted to rows where the cloud is fully on the grid, E_PP is **bit-exactly
  constant** (spread 0.000e+00 at `norm_proj == 1.0`, 1.8e-11 Ha at
  `norm_proj >= 1-1e-9`) — the rigid-cloud expectation holds perfectly.
- **Clipping onset t = 21.04 a.u.** for sigma=3 r_s=5.702. The run's fit window
  ends at 19.48 (common cross-sigma window 9.37), so **the physics is untouched**.
  The verifier now REPORTS this onset as the hard upper bound on any fit window.

Two reporting bugs found and fixed while validating: a tolerance that assumed
exact constancy everywhere, and a clipping-onset that used a global `min()` and
so mislabelled the launch point (t=0, z=-28) as the onset. It must be the start
of the trailing clipped run; the first 3 rows sit ~4e-9 under norm 1.0 purely
from discretising the Gaussian, which is not clipping.

### Submitted (autonomous)

```
preflight build : 32526619   (builds all 6, fails the job if ANY fails)
production array: 32526620   (tasks 0-5, --dependency=afterok:32526619)
```

Scripts: `shared/bin/run-bulk-ie-build.slurm`, `shared/bin/run-bulk-ie-rerun.slurm`.

The preflight exists so a compile error stops ALL six before consuming GPU hours
— and prevents a PARTIAL set of twin pairs, which would be useless for the
comparison. Each array task runs its closure gate after the run and exits
non-zero if the terms do not sum back to INQ's scalars.

Cost: ~9.2 GPU-h total, ~2.2 h wall in parallel. Calibrated on the measured
sigma=3 pair (1:28 wp / 1:56 classical). `--time=08:00:00` (~4x headroom; the
classical half is not resumable). No new GS needed — GS depends only on box +
density. Disk growth ~0 (VTI filenames deterministic, overwritten in place).

### STATUS: COMPLETE — all 6 runs finished, every closure gate passed

| task | run | elapsed | rows | closure |
|---|---|---|---|---|
| 0 | `bulk_ks_stopping/wp` | 1:33:21 | 325 | ALL PASSED |
| 1 | `bulk_ks_stopping/classical` | 2:01:40 | 325 | ALL PASSED |
| 2 | `bulk_ks_stopping_rs4/wp` | 1:07:29 | 325 | ALL PASSED |
| 3 | `bulk_ks_stopping_rs4/classical` | 1:25:38 | 325 | ALL PASSED |
| 4 | `bulk_ks_stopping_rs4_sigma3/wp` | 1:02:06 | 302 | ALL PASSED |
| 5 | `bulk_ks_stopping_rs4_sigma3/classical` | 1:20:01 | 302 | ALL PASSED |

8.5 GPU-h. **All four bulk twin pairs now carry `interactions.csv` on both halves.**

### THE RESULT — self-Hartree collapse tracks WIDTH, not DENSITY

`hypotheses/bulk_ks_stopping_sigma3/summarise_epp_across_pairs.py` (new),
windowed to `[4.0, min(FIT_T1, clipping onset)]`:

| pair | E_PP(0) eV | dE_PP eV |
|---|---|---|
| r_s 5.702, sigma 2 | 4.8149 | **-4.3986** |
| r_s 5.702, sigma 3 | 3.0082 | **-2.3607** |
| r_s 3.987, sigma 2 | 4.8162 | **-4.4063** |
| r_s 3.987, sigma 3 | 3.0103 | **-2.3832** |

- **width effect** (sigma3 - sigma2): **+2.04 eV** at r_s 5.702, **+2.02 eV** at r_s 3.987
- **density effect** (r_s 3.987 - 5.702): **-0.008 eV** at sigma 2, **-0.023 eV** at sigma 3
- non-additivity: -0.015 eV

The width lever is ~100-250x stronger than the density lever, and the two
factorise cleanly.

**CONSEQUENCE — E_PP is DECOUPLED from the stopping ratio.** The classical/WP
stopping ratio moves 6.49 -> 5.65 (13 %) across the 2.92x density lever, but
dE_PP moves 0.2 % across that same lever. A mechanism cannot drive a 13 % effect
while itself flat to 0.2 %. **Self-Hartree collapse does not explain the
density dependence of the stopping ratio.**

This is the SECOND independent line of evidence against E_PP as the cause of the
classical/WP gap. The first: the within-slab gauge-clean twin comparison
(same entry, below) found slab dE_PP = -5.7 to -6.0 eV against bulk's -2.4,
i.e. MORE collapse where the user reports a SMALLER stopping discrepancy —
anti-correlated. Both point the same way: **look elsewhere for the ~2.2 residual.**

### INDEPENDENT VALIDATION of the E_PP values (1 degree of freedom)

For a Gaussian charge of std `sigma_pot = sigma_WP/sqrt2` the isolated
self-energy is `A/sigma_WP` with `A = 1/sqrt(2 pi) Ha.Bohr = 10.8556 eV.Bohr`.
INQ drops G=0 in a periodic cell, subtracting a sigma-INDEPENDENT constant C, so
`E_PP_measured = A/sigma - C`. Solving from the two sigmas:

| | A fitted | vs analytic | C (G=0 gauge) |
|---|---|---|---|
| r_s 5.702 | 10.8402 | **-0.14 %** | 0.6052 eV |
| r_s 3.987 | 10.8354 | **-0.19 %** | 0.6015 eV |

The measured coefficient reproduces the analytic Gaussian self-energy to <0.2 %
at BOTH densities, and the gauge offset agrees (0.605 vs 0.602 eV) across two
DIFFERENT box sizes. The decomposition is physically correct, not merely
internally consistent.

Two further gates pass on all four pairs: `dE_PP(0) = 2.6e-8 eV` (confirms the
sigma_pot = sigma_WP/sqrt2 matching NUMERICALLY) and classical `E_PP` drift
`5.6e-9 eV` (rigid cloud).

### Notebooks regenerated

`hypotheses/<variant>/energy_component_comparison_<variant>.ipynb` rebuilt for
the 3 re-run variants, 15 cells each, 0 errors, section 7 now populated. In-notebook
closure: WP `E_SS+E_PS+E_PP` vs `energy_hartree` = 8.6-9.0e-13 Ha; classical
`E_SS` vs `energy_hartree` = 0.00e+00 Ha.

**WINDOW CAVEAT — two different dE_PP numbers are in circulation, both correct.**
The notebooks report the END-OF-RUN change (e.g. -4.6286 eV for r_s 5.702
sigma 2 at t = 25.84); `summarise_epp_across_pairs.py` reports the
END-OF-FIT-WINDOW change (-4.3986 eV at t = 18.43). The windowed value is the
one to quote for physics — it excludes the region where the classical
projectile cloud is clipped by the box face.

### Known non-issue: interactions.csv appears empty mid-run

`std::ofstream ix` has no explicit flush, so the file shows 0 bytes until the
destructor runs at clean exit (verified mid-run on task 0 at step 32). This
matches the sigma=3 reference, which produced its full 302 rows. **A KILLED run
would lose its buffered rows**, unlike `observables.csv` whose writer flushes.
Hardening candidate (add `ix.flush()` in the callback) if these files are ever
touched again; not changed mid-flight.

---

## 2026-08-01 (later) — PAIRWISE INTERACTION ENERGIES. E_PP measured at last.

### Why this happened

User asked why the KE magnitudes differ so much between the two halves, and for
per-component energy plots. Building those exposed that **INQ's scalars cannot
answer it**: the two runs put the projectile in different ledger terms, so
`energy_hartree` means `E_SS` in the classical run and `E_SS + E_PS + E_PP` in
the WP run. Comparing them compares a NET against a GROSS.

### Delivered

**6 energy-component comparison notebooks**, one per (density, sigma) pair, at
`hypotheses/<variant>/energy_component_comparison_<variant>.ipynb` (11 cells each,
0 errors). Component activity table, absolute curves, **residuals on shared
axes**, total energy + closure check, full ledger with cl/wp ratios, checksum
(8 summed components reproduce `energy_total` to 1e-12 eV).

Builder: `hypotheses/bulk_ks_stopping_sigma3/build_energy_comparison.py`, takes a
variant argument, writes into that variant's OWN hypotheses folder.

**sigma=3 LOW-DENSITY pair (r_s = 5.702) RE-RUN with the decomposition wired in**
— jobs 32512952 (wp, 1:27:59) / 32512953 (classical, 1:55:48), both exit 0. That
notebook is now **15 cells** and has §7 with the pairwise terms.

### THE RESULT — the sign flip was a self-interaction artefact

Closure gates PASS (asserted in the notebook, not assumed):
`WP E_SS+E_PS+E_PP vs energy_hartree: 9.46e-13 Ha`;
`classical E_SS vs energy_hartree: 0.00e+00 Ha`.

| term | WP | classical | ratio |
|---|---|---|---|
| ΔE_SS (bath polarisation) | +1.278 eV | +4.738 eV | 3.71x |
| **ΔE_PP (projectile self-Hartree)** | **−2.568 eV** | **−0.016 eV** | — |
| ΔE_PS (projectile-bath) | +0.151 eV | +0.702 eV | **4.65x** |
| E_SB, E_PB, E_BB | 0 | 0 | (bulk: phi+ ≡ 0) |

1. **E_PP is IDENTICAL at t=0 in both runs: 3.0082 eV.** Independent confirmation
   that `sigma_pot = sigma_WP/sqrt2` puts genuinely the same charge cloud in both.
   Nothing forced this — it validates the sigma-matching convention.
2. **E_PP is constant for the classical projectile (−0.016 eV) and collapses for
   the WP (−2.568 eV, 85 % of it).** The classical Gaussian is RIGID — it cannot
   change its own self-energy. The packet spreads, so its self-Hartree falls.
   Wavepacket-only, by construction.
3. **This resolves the sign flip the user spotted.** WP Δ`hartree` = −1.14 eV vs
   classical +4.74 eV decomposes as −2.57 (self-Hartree collapse) + 1.28 (ordinary
   bath polarisation) + 0.15, summing exactly to the measured value. The two
   projectiles do NOT polarise the gas in opposite directions; a spurious self-term
   was swamping a normal one.
4. **ΔE_PS gives 4.65x against this pair's measured stopping ratio of 4.80x.**
   First quantity to track the drag WITHOUT the unexplained ~2.2 factor.
   ONE PAIR ONLY — suggestive, not established.

*Caution to keep:* E_PP is the self-**Hartree**, not the self-interaction ERROR.
LDA exchange cancels part of it; the SIE is the residual. E_PP BOUNDS the
contamination rather than measuring it.

### Standing requirement added (user: "all future runs, does not matter which system")

- **NEW RULE** `.claude/rules/decomposed-interaction-energies.md` — every
  projectile run, every system, writes `interactions.csv`. Closure gates, the
  bulk-vs-slab phi+ distinction, the gauge caveat, the `sigma_pot = sigma_WP/sqrt2`
  trap for classical `n_P`.
- **`tddft-simulations` skill §3a** — `interactions.csv` promoted to **Tier 1
  always-on**, with wiring snippet and both reference implementations named.

### IMPORTANT CORRECTION — the slab runs ALREADY have this data

A claim that "no simulation tracked these" is WRONG. Counts on disk:

| system | interactions*.csv |
|---|---|
| **localised_jellium** | **56** |
| jellium (bulk) | 2 (the sigma=3 low-density pair) |
| cylindrical_jellium / vacuum / coronene | 0 |

13 `run.cpp` under `localised_jellium` include the header (`slab_ks_wrap` cl+wp
across every density/velocity, `wp_highdensity_sv`, `classical_highdensity_sv`,
`localised_jellium_dynamics`, `sigma1_masspair`). Files carry the full schema
incl. closure columns, ~2762 rows each.

**CONSEQUENCE — a FREE test exists.** The user's original observation was that the
slab runs show a MUCH SMALLER WP/classical discrepancy than bulk. Those runs
already have E_PP on disk. If slab ΔE_PP is smaller than bulk's −2.57 eV, that
ties the smaller slab gap directly to less self-Hartree collapse — **no new
compute**. Slab packets are sigma_WP = 0.5-3 and cross a finite slab rather than
dispersing over 80 Bohr of bulk, so there is a real reason to expect a difference.
**DO THIS BEFORE spending ~10 GPU-h retrofitting the other 10 bulk runs.**

### Wiring notes for whoever extends this

- Bulk: `phi_plus` must be an explicitly ZERO-filled field (uniform background ->
  poisson is pure G=0, which INQ drops). Slab: `bg_pert.background_density(basis)`.
- Classical `n_P` is rebuilt EVERY step at the ion's current position via
  `gaussian_density(basis, center, sigma_pot)` — needs
  `inqkit/jellium/projectile_background_energy.hpp`, NOT just the interaction header.
- `StepContext` exposes `ions` and `electrons` as pointers; `energy_hartree` is a
  field on it, so the classical closure can be written per row without a re-read.

### Still open

1. **Slab E_PP comparison** — free, see above. Highest value per unit effort.
2. Retrofit `interactions.csv` to the other 10 bulk runs (~10 GPU-h).
3. classical-at-4-Bohr experiment (2 runs) — still unsubmitted.
4. `test_wp_circular_centroid_engine.cpp` — STILL never compiled or run.

---

## 2026-08-01 — CAMPAIGN COMPLETE. The gap factorises; the residual is unexplained.

**12/12 runs complete, 12/12 run notebooks built (0 errors), microscopy notebook
covers all six (density, sigma) pairs.**

### THE RESULT — final, all six pairs, time-matched sampling

| pair | S_ratio | lag asymmetry | residual |
|---|---|---|---|
| r_s=5.702 sigma=1 | 8.73 | 3.81 | 2.29 |
| r_s=5.702 sigma=2 | 5.21 | 2.46 | 2.12 |
| r_s=5.702 sigma=3 | 4.80 | 2.15 | 2.23 |
| r_s=3.987 sigma=1 | 7.24 | 3.10 | 2.34 |
| r_s=3.987 sigma=2 | 4.00 | 1.83 | 2.19 |
| r_s=3.987 sigma=3 | 3.60 | 1.45 | 2.49 |

    S_ratio  =  [lag asymmetry]  x  [2.2 +/- 0.2]

- **Lag asymmetry carries EVERYTHING about the medium and the geometry.** It is
  MONOTONIC in in-flight width in both baths (2.15 -> 2.46 -> 3.81 and
  1.45 -> 1.83 -> 3.10) and scales with density. It tracks S_ratio across the
  whole 8.73 -> 3.60 range.
- **The residual is FLAT: 2.12-2.49, spread factor 1.17**, across a 2.92x density
  change AND a 2.1x in-flight width change (3.45 -> 7.35 Bohr). Not monotonic in
  width in either bath.

### TWO HYPOTHESES PROPOSED AND BOTH REFUTED BY THEIR OWN PREDICTIONS

| hypothesis | its prediction | outcome |
|---|---|---|
| geometric back-action form factor | residual varies (monotonically) with packet width | **REFUTED** — flat, non-monotonic |
| charge / occupation mismatch | WP presents != 1 electron | **REFUTED** — measured directly: bath states occ 2.0, WP state occ **1.0**, sum 219 = 218 + 1; classical UPF `z_valence 1.00`, Coulomb tail coefficient 2.000000 Ry = 1 Ha. **Both projectiles carry exactly charge -1.** |

**The residual ~2.2 is an unexplained constant.** It is not the medium
(density-invariant), not geometry (width-invariant), and not charge
normalisation (verified matched). **Do not attach a third guess without a test
that predicts something.**

Untested candidate, named but NOT claimed: **self-interaction**. The WP is an
occupied KS orbital, so its own charge enters the Hartree potential acting on it,
uncancelled in LDA. Quantum-only, no classical counterpart. A test would be to
remove the orbital's self-interaction, or vary its occupation and rescale.

### THE ANSWER TO "why does the classical projectile slow faster", in order

1. **It stays compact.** S ~ w^(-0.244) (dilute) / w^(-0.286) (dense), |r| >= 0.999,
   and the SAME law governs the wavepacket. Not a quantum/classical distinction —
   a width difference that correlates with representation because only the
   quantum object disperses.
2. **It sits deeper in the non-linear regime** — 68 % local depletion of n_0 vs
   37 % for the WP. NEITHER is in linear response; say so anywhere Lindhard is
   mentioned.
3. **Its polarisation cloud lags more**, and that lag carries all the density and
   width dependence of the gap.
4. **A flat ~2.2 residual remains**, unexplained.

Energy channels (classical): deposited energy goes overwhelmingly into
**electron kinetic energy** (+17.5 of +22.6 eV dilute; +40.3 of +54.9 dense),
with Hartree second. Closure vs projectile KE loss: 0.22 % / 0.96 %.

### ANALYSIS BUG FOUND AND FIXED — read this before trusting any cross-run number

The six-pair build FIRST reported residual 3.32 for `r_s=3.987 sigma=1` and its
verdict cell flipped to "the form factor survives". **That was an artefact I
created.** Reducing `WRITE_EVERY` 2 -> 8 and re-running that one job left it with
87 density frames against 301-347 for every other run; a blanket `stride=8` then
sampled it **4x more sparsely**. Time-matched sampling gives **2.34**, and it
RESTORES monotonicity in the asymmetry ratio that the artefact had broken.

Hardenings, both now in the code:
- **`microscopy.stride_for(run_dir, kind, target_frames)`** — strides to a target
  COUNT, so runs with different cadences are sampled comparably. Its docstring
  names this incident. **Never hard-code a frame stride when comparing runs.**
- **The verdict cell now requires MONOTONICITY**, not just spread. The old form
  would call any spread > 1.5x "consistent with a form factor" including pure
  scatter.

**General lesson:** when a config change alters data layout, every analysis
parameter expressed in that layout's units becomes suspect — frame strides,
cadence-derived indices, window step-counts. Nothing errored; it produced a
plausible number that happened to match a hypothesis.

### SECOND DEFECT FOUND AND FIXED — energy cadence was coupled to VTI cadence

`RealTimeSession` takes `WRITE_EVERY` as its CALLBACK stride, and
`obs_writer.append()` lives inside that callback — so `observables.csv` (the
ENERGY series) followed the VTI cadence. Raising it 2 -> 8 thinned the classical
stopping fit from **68 points to ~17** in the common window. It would not have
failed; it would have silently produced noisier S_classical. The proposed
classical-at-4-Bohr experiment is *entirely* an energy-deposit measurement and
would have been the first casualty.

Fixed in all 12 run.cpp: callback now fires at **`OBS_EVERY = 2`** (energies
dense) with the VTI writes gated inside on `WRITE_EVERY = 8`. Storage saving
intact, measurement cadence restored. WP runs were never affected — their
stopping data is `wp_momentum_stats.csv` at `STATS_EVERY = 1` (693 rows, verified).
NOTE the classical runs needed a different patch (they build two sessions,
`rt_dens` + a per-step `rt_track`); the first sweep silently patched only the 6 WP
files. Caught by verifying the count, not by trusting the "12 patched" message.

### Artefacts

- 12 run notebooks, one per run, in `hypotheses/<variant>/` (ADR 0007), each with
  `ks_stopping.py` + `tests/` SYMLINKED to `hypotheses/bulk_ks_stopping/` so all
  twelve compute from ONE implementation and cannot disagree.
- `hypotheses/bulk_ks_stopping_microscopy/classical_vs_wavepacket_microscopy.ipynb`
  — 30 cells, per-density sections + the sigma section + verdict.
- Notebooks total ~550 MB against ~270 GB of raw output: **the figures a reader
  needs are already embedded**, so raw density frames are more prunable than they
  look if disk pressure returns.

### STILL OPEN

1. **classical-at-4-Bohr experiment** — 2 classical runs at sigma_charge ~3.5-4
   Bohr, inside the WP's in-flight width range. Deletes the extrapolation behind
   the 3-5x width-law residual and measures it instead of bounding it. Needs one
   new UPF (sigma_WP ~5-5.7). **PROPOSED, NOT APPROVED, NOT SUBMITTED.**
2. Self-interaction test for the flat ~2.2 residual (see above).
3. `test_wp_circular_centroid_engine.cpp` — still never compiled or run.

---

## 2026-07-31 (night) — DISK CRISIS, storage policy, campaign audit

### The filesystem hit 100 % and killed a run

`rs4_sigma1-wp` (32483981) aborted at step 688/692 with
`VTIImageDataWriter: failed while writing ... density_t000690.vti`. **Not a code
bug** — quota was at 1099.4 GB of a 1099.5 GB soft limit, 77 MB free. Physics was
healthy to the last step (energy drifting in the 12th decimal).

### Deletions (USER-APPROVED, 270 GB freed)

Quota 1099.4 → **829.7 GB (75 %)**.

1. **`results/raw/vti/density_delta`** — verified redundant BEFORE deleting:
   `inqkit/observables/density_delta.hpp` defines it "relative to a
   user-supplied reference density n(r, t₀)", i.e. exactly n(t) − n(0), and every
   run's `density_total` was confirmed to contain the t=0 frame (checked all 12).
   Zero information lost.
2. **`results/checkpoint`** (RT checkpoints, all 12). User approved this
   explicitly after I flagged that it contradicts
   `.claude/rules/final-timestep-checkpoint.md`.
   **CONSEQUENCE, now realised:** no run can be EXTENDED — any additional
   timesteps require a FULL re-run. This bit immediately: the aborted run above
   could have been finished by a 4-step resume; instead it needed a complete
   692-step re-run (32499918).
   **GS checkpoints were NOT touched** and remain intact (4.0 + 5.6 GB) — every
   run is still reproducible from scratch.
3. `wavefunction_wp` (25 GB) was offered and **declined** — still on disk.

**A `df -BG` figure I reported mid-way (31 GB) was wrong** — it read the raw 13 PB
Lustre pool, not the quota. Use `quota -s`, not `df`, on this filesystem.

### STORAGE POLICY (user instruction: "decrease the cadence of VTI production")

Applied to the two BASE config headers (the 4 sigma variants inherit):

| | before | after |
|---|---|---|
| `WRITE_EVERY` (density) | 2 | **8** |
| `WF_WRITE_EVERY` (complex psi) | 6 | **32** |

Plus `emit_raw_vti = false` in **all 12 run.cpp** — stops writing the redundant
full-resolution delta while keeping the coarse preview and the `density_l2`
series. Projected: WP run 36.5 → **6.3 GB** (−83 %), classical 19.7 → **2.5 GB**
(−87 %); a dilute 4-run set drops 112 → 18 GB. All six headers still compile.

**Costs nothing analytically:** the microscopy notebook already read with
`stride=8` on every-2-step output (= every 16th step). `WRITE_EVERY = 8` yields
87 frames natively — MORE than the 41 the wake analysis actually used.

**Perspective worth keeping:** every stopping number in this study came from
`observables/` CSVs, **under 1 GB across 14 runs**. The 270 GB was density frames
for pictures.

`microscopy.py` gained **`induced_series()`**, which derives Δn by subtracting the
t=0 frame of `density_total` — reproducing the deleted field exactly. The notebook
builder and its prose now use it, so the notebook rebuilds against runs that have
no `density_delta` directory. **Do not reintroduce reads of `density_delta`.**

### CAMPAIGN AUDIT — 11/12 complete, 1 re-running

| bath | sigma=1 | sigma=2 | sigma=3 |
|---|---|---|---|
| r_s = 5.702 | wp ✓ cl ✓ | wp ✓ cl ✓ | wp ✓ cl ✓ |
| r_s = 3.987 | wp **re-run 32499918** cl ✓ | wp ✓ cl ✓ | wp ✓ cl ✓ |

All reached full step count with `run_summary.txt`, except the re-run.
**All 12 have complete, usable science data** — the aborted run covered its fit
window [4, 9.37] nearly 3× over (reached t = 27.52), so no result changes. What
was missing was provenance (`run_summary.txt`), not physics. Note the re-run
writes at the NEW cadence, so it will have ~87 density frames against its
classical twin's 347 — harmless, the analysis interpolates in time.

### Housekeeping done

- **Stale gate label corrected** in 3 files (`bulk_ks_stopping{,_sigma1,_sigma3}/wp`):
  printed `sigma_pz^2 = 1/(4 sigma^2)` while checking the CORRECT `1/(2 sigma^2)`
  value. Label only — no computed value or gate outcome changed. **The existing
  run logs on disk still show the old text**; source and logs differ by this
  string alone.
- **Pauli-aware sigma_z gate ported** to `bulk_ks_stopping_sigma1/wp` (was only in
  `rs4_sigma1`). sigma=2 and sigma=3 keep the tight ±0.05 gate — their
  `max_overlap` is ~1e-5 or less, so no broadening occurs and the tight gate is
  still the right check there.

### Still outstanding

1. Notebooks for the **8 sigma-sweep runs** (only the four sigma=2 runs have them).
2. **Extend the microscopy notebook to sigma** — its per-density structure takes
   sigma sections directly. THIS IS THE ONE THAT TESTS THE OPEN PHYSICS: does the
   ~2.15 back-action residual move with sigma?
3. **classical-at-4-Bohr experiment** (2 runs) — proposed, NOT approved, NOT
   submitted. Turns the 3–5× residual from a bound into a measurement.
4. `test_wp_circular_centroid_engine.cpp` — still never compiled or run.

---

## 2026-07-31 (late) — SIGMA SWEEP RESULTS: width is real but does NOT close the gap

All fits on the COMMON window **[4.0, 9.37] a.u.** Two runs were still executing
when these were extracted, but **both had already passed t = 9.37** (classical at
t = 26, rs4_sigma1-wp at t = 19.1) and observables are written incrementally, so
every number below is over complete in-window data and will not change.

### Classical: a clean width law, and it factorises from density

| sigma_charge (Bohr) | S_cl r_s=5.702 | S_cl r_s=3.987 |
|---|---|---|
| 0.707 (sigma_WP=1) | **0.3237** | **0.8071** |
| 1.414 (sigma_WP=2) | 0.2765 | 0.6619 |
| 2.121 (sigma_WP=3) | 0.2471 | 0.5894 |

r² ≥ 0.990 throughout. Power law **S ~ w^(−0.244)** (dilute), **w^(−0.286)**
(dense), |r| ≥ 0.999. Relative width dependence is IDENTICAL at both densities
(1.17/1.00/0.89 vs 1.22/1.00/0.89) — **width and density factorise**.

A compact charge couples to short-wavelength density response; a broad one
cannot. Shrinking the classical charge 2.12 → 0.71 Bohr buys 24 % more drag.

### Wavepacket: ordered by IN-FLIGHT width, not sigma(0)

| bath | sigma_WP | sigma_d at t=9.37 | S₂ | r² |
|---|---|---|---|---|
| r_s=5.702 | 3 | 3.451 | 0.0515 | 0.992 |
| r_s=5.702 | 2 | 4.069 | 0.0531 | 0.997 |
| r_s=5.702 | **1** | **7.217** | **0.0371** | 1.000 |
| r_s=3.987 | 3 | 3.556 | 0.1637 | 0.996 |
| r_s=3.987 | 2 | 4.211 | 0.1653 | 0.999 |

**sigma = 1 gives the LEAST stopping of the three** — the crossover working
exactly as predicted before submission. It starts narrowest, disperses fastest,
and is the WIDEST projectile in the study by mid-window. Δ⟨p_z⟩ tracks S
independently (0.0100/0.0104/0.0073), so this is real deceleration, not a fit
artefact.

**The width response has a knee.** 3.45 → 4.07 Bohr (18 % wider): S moves 3 %,
i.e. fit noise. 4.07 → 7.22 Bohr (77 % wider): S drops 30 %. Flat over 3–4 Bohr,
falling away beyond — which is why the sigma=3 pair alone looked like saturation
and sigma=1 does not contradict it. Same curve, different parts.

### Ratios S_cl / S_wp

| | sigma=1 | sigma=2 | sigma=3 |
|---|---|---|---|
| r_s=5.702 | **8.73** | 5.21 | 4.80 |
| r_s=3.987 | *(pending 32483981)* | 4.00 | 3.60 |

sigma=1 swings the ratio to 8.7 because both effects compound: the classical
charge gets more compact AND the packet gets wider.

### THE HEADLINE: width explains a lot, and is NOT enough

Extrapolating the classical width law to the wavepacket's in-flight width:

| bath | WP in-flight w | S measured | width-law prediction | residual |
|---|---|---|---|---|
| r_s=5.702 | 4.069 | 0.0531 | 0.2119 | **4.0×** |
| r_s=5.702 | 7.217 | 0.0371 | 0.1843 | **5.0×** |
| r_s=3.987 | 4.211 | 0.1653 | 0.4844 | **2.9×** |

**A wavepacket sitting at 4.1 Bohr does NOT stop like a classical charge of
4.1 Bohr — it stops 3–5× less.** The residual is in the same range as the
independent ~2.1–2.2 back-action factor from the microscopy notebook (different
data entirely: induced-density asymmetry, not stopping fits). Two unrelated
routes both landing at "a few ×" for the non-geometric part.

**WEAK POINT, stated plainly:** the classical law is fitted over w = 0.71–2.12
and used at w = 3.5–7.2 — 2–3× beyond range, on three points. Excellent within
range, but a true form factor should STEEPEN once w exceeds the field scale, and
a steeper curve predicts less stopping at 4 Bohr and shrinks the residual. So
**3–5× is an UPPER BOUND on the non-geometric part, not a measurement of it.**

### THE DECISIVE FOLLOW-UP (proposed, awaiting user go-ahead)

**Classical-only runs at sigma_charge ≈ 3.5–4 Bohr**, i.e. directly inside the
wavepacket's in-flight width range. This DELETES the extrapolation and measures
the residual instead of bounding it. 2 runs, ~1.5 h each, needs one new UPF
(sigma_WP ≈ 5–5.7; the repo's widest is wpsigma3p5 = sigma_charge 2.47).

- If the classical at 4 Bohr still stops ~3× harder than a WP at 4 Bohr, the
  residual is real and genuinely quantum (back-action averaging / Pauli / SIE).
- If it drops to meet the WP, the whole gap was width and the power law simply
  understated the falloff.

### Job accounting

| Job | Run | State | Elapsed |
|---|---|---|---|
| 32482667 | sigma1-wp | COMPLETED | 1:44:14 |
| 32482668 | sigma1-classical | running (past window) | 2:10+ |
| 32482669 | sigma3-wp | COMPLETED | 1:36:25 |
| 32482670 | sigma3-classical | COMPLETED | 1:58:48 |
| 32482671 | rs4_sigma1-wp | **FAILED 4:0** (sigma_z gate) | 0:05:57 |
| 32482672 | rs4_sigma1-classical | COMPLETED | 1:37:43 |
| 32482673 | rs4_sigma3-wp | COMPLETED | 1:06:20 |
| 32482674 | rs4_sigma3-classical | COMPLETED | 1:23:20 |
| 32483981 | rs4_sigma1-wp (resubmit, patched gate) | running (past window) | 1:43+ |

### TODO carried forward

1. Write sigma-sweep results into the notebooks (extend
   `bulk_ks_stopping_microscopy` — its per-density structure takes sigma
   sections directly; do NOT rebuild).
2. Fix the stale `sigma_pz^2 = 1/(4 sigma^2)` label in
   `scripts/bulk_ks_stopping_{sigma1,sigma3}/wp/run.cpp` (value is correct, label
   is not) and port the Pauli-aware sigma_z gate from `rs4_sigma1`.
3. `test_wp_circular_centroid_engine.cpp` STILL never compiled or run (GPU node).

---

## 2026-07-31 (evening) — SIGMA SWEEP submitted (8 runs) + microscopy notebook

### Why: the width axis is what the density sweep left open

Density moved the S_classical/S_WP ratio only 13 %. Width is the remaining
candidate. Sweep: sigma_WP ∈ {1, 3} at BOTH densities, 100 eV, both projectile
representations = 8 runs. **No new GS needed** — the GS depends only on box +
density, not the projectile — so all eight launched at once with no afterok chain.

| Jobs | Variant | Bath | sigma | z0 | N_STEPS |
|---|---|---|---|---|---|
| 32482667/8 | `sigma1` wp/cl | r_s = 5.702 | 1 | −36 | 692 |
| 32482669/70 | `sigma3` wp/cl | r_s = 5.702 | 3 | −28 | 600 |
| 32482671/2 | `rs4_sigma1` wp/cl | r_s = 3.987 | 1 | −36 | 692 |
| 32482673/4 | `rs4_sigma3` wp/cl | r_s = 3.987 | 3 | −28 | 600 |

### THE COUNTER-INTUITIVE FACT THAT DRIVES THE DESIGN

**A narrower packet spreads FASTER.** sigma_d(t) = sqrt(sigma²/2 + t²/(2 sigma²)),
so small sigma has the larger t² coefficient. sigma=1 and sigma=2 **cross at
t = 2.0 a.u.**; beyond that sigma=1 is the *broader* projectile (13.05 vs 6.67
Bohr at t = 18.4). The width minimising ARRIVAL width is sigma_opt = sqrt(t) ≈ 3–4.

**In-flight width ordering is therefore sigma=1 (widest) > sigma=2 > sigma=3
(narrowest) — the REVERSE of the t=0 ordering.** Any reading of this sweep must
use in-flight sigma_d, not sigma(0). The user proposed sigma=1 as "more
concentrated"; it is, for the first 2 a.u. only (~5 Bohr of a 38 Bohr path).

sigma=4 was offered as the genuinely-narrowest-in-flight point; user judged it
too high and asked what else exists. Project has used sigma_WP ∈ {0.5, 0.53, 1,
2, 3, 3.5, 5, 8}; **sigma=3 chosen** — UPF already on disk, project precedent,
sigma_d(18.4) = 4.83 vs sigma=4's 4.31 (nearly as narrow, far less extreme).

### Launch position: BOUNDARY RULE (user reversed an earlier choice mid-turn)

z0 = −L_z/2 + 4 sigma → **−36 (s=1), −32 (s=2), −28 (s=3)**. Constant 4-sigma
clearance. **Consequence: path lengths differ per sigma, so cross-sigma
comparison MUST use a common TIME window, not raw path.**

**COMMON CROSS-SIGMA WINDOW = [4.0, 9.37] a.u.**, set by sigma=1 in the 40-Bohr
box (fastest spreader, smallest transverse cell). Pinned in every new header as
`COMMON_FIT_T1_AU` with a static_assert that it lies inside each run's own
interference limit. Validated by re-fitting the COMPLETED sigma=2 runs on it:
S₂ shifts only 5–8 %, systematics ~3× larger but still ~10 % of the value.

**Side finding from that check, important for interpretation:** the CLASSICAL S
drops **26 %** on the short window at both densities (0.375→0.277, 0.890→0.662).
Not noise — the classical projectile decelerates, and S ∝ ln v / v² rises as v
falls, so the full window averages over a velocity range. By
`.claude/rules/light-projectile-stopping.md` the short-window value is arguably
the *more* correct S(v₀). Ratios on the short window: **5.21× and 4.00×**.

### Verifications done before submitting

- **sigma=1 UPF GENERATED** (`electron_gaussian_wpsigma1p0.upf`, sigma_charge =
  1/√2). Validated independently vs the analytic erf form: max deviation
  **5.3e-11**, Coulomb tail coefficient exactly 2.000000, core depth exactly
  **2.000×** the sigma=2 file. sigma=3's UPF already existed.
- **Nyquist**: sigma=1 has the tightest margin in the study (1.30× at dx=0.5).
  Numerically verified on both grids — exact to machine precision at dx=0.4;
  **1.3e-3 %** error on T₁−T₂ at dx=0.5 (= 2.6e-4 eV of 20.4 eV, negligible).
- All four generated config headers **compile**, so every static_assert passes.

### CAVEAT discovered at run time (sigma=1 only)

`sigma1/wp` t=0 gates all PASS, but two numbers deserve care:

- **`max_overlap` = 1.3e-3, up 400× from sigma=2's 3.1e-6.** Real physics: a
  narrower packet has sigma_p = 0.707, which reaches into the occupied Fermi
  sphere (k_F = 0.337). **Gram–Schmidt now actually removes something** and Pauli
  blocking is no longer negligible. This is a quantum-only mechanism with no
  classical counterpart — relevant to any width conclusion.
- **sigma_z measured 3.2 % high** (0.7297 vs 0.7071) — the largest t=0 deviation
  in the study. At sigma=1 the density std is 0.707 Bohr against dx = 0.4, i.e.
  only **1.8 grid points per sigma**; `rs4_sigma1` at dx = 0.5 is coarser still
  (~1.4). Does NOT contaminate S (momentum moments are fine, and the circular
  centroid is exact at −36.0000), but **label any sigma_z(t) curve at sigma=1 as
  carrying a few-% real-space discretisation bias.**

### Known cosmetic defect, NOT fixed mid-flight

`scripts/bulk_ks_stopping_{sigma1,sigma3}/wp/run.cpp` print the gate label
`sigma_pz^2 = 1/(4 sigma^2)` while checking the CORRECT `1/(2 sigma^2)` value —
the stale label from the 2026-07-30 momentum-width error. Fixed in the `rs4`
copies; these two were copied from the original before that fix. Editing
`run.cpp` mid-run would desync source from the running binaries. **TODO: fix in
both dirs once the runs complete.**

### MICROSCOPY RESULT — the gap FACTORISES

`classical_vs_wavepacket_microscopy.ipynb` built and executed, **26 cells, 0
errors**. Both densities, both projectiles, common window [4, 9.37].

**Energy channels (classical), whole run:**

| channel | r_s = 5.702 | r_s = 3.987 |
|---|---|---|
| electronic kinetic | **+17.54** | **+40.29** |
| Hartree | +5.44 | +26.13 |
| xc | −1.07 | −1.31 |
| external | +0.72 | −10.19 |
| total | +22.63 | +54.92 |
| projectile KE loss | +22.58 | +54.39 |

Closure 0.22 % / 0.96 %. Deposited energy goes overwhelmingly into **electron
kinetic energy**, not field energy.

**Induced density — BOTH projectiles are strongly non-linear:**

| | classical peak \|Δn\| | WP peak \|Δn\| |
|---|---|---|
| r_s = 5.702 | **68.1 % of n₀** | 37.0 % of n₀ |
| r_s = 3.987 | **50.1 % of n₀** | 29.4 % of n₀ |

Nothing in this study is in the Lindhard/linear-response regime. Any comparison
with linear-response stopping theory must say so.

**THE SYNTHESIS — the headline of this notebook:**

```
              S ratio   asymmetry ratio   depletion ratio
r_s = 5.702    5.19x         2.46x            1.76x
r_s = 3.987    3.97x         1.83x            1.42x
```

S ratio falls ×0.765 with density; asymmetry ratio falls ×0.744. **They track to
3 %.** So the *entire density dependence* of the classical/WP gap lives in the
front/back lag of the induced cloud.

What remains is a near-constant multiplier:

    S_ratio / asymmetry_ratio  =  5.19/2.46 = 2.11
                               =  3.97/1.83 = 2.17

**2.11 and 2.17 — density-independent to 3 % across a 2.92× density change.**

*Inference (mine, NOT directly computed — the force integral has not been
evaluated):* this is the **back-action form factor**. Drag is force × charge; the
classical point charge samples the induced field at ONE point, while the
wavepacket samples it AVERAGED over its own width. A field varying on the
packet's scale partially cancels under that average — a purely geometric
suppression, which is exactly what should be density-blind. So:

    S_ratio = [lag asymmetry: density-dependent] x [extended-object averaging: ~2.15, geometric]

**This also resolves the dilution puzzle.** Simple charge dilution predicted
curvature in T₂(s) as σ grows; the data showed r² = 0.997 with σ quadrupling. A
back-action FORM FACTOR does not predict curvature — it saturates with the ratio
of field scale to packet scale rather than growing with σ.

**THE SIGMA SWEEP IS A DIRECT TEST OF THIS.** If the ~2.15 residual is geometric
it MUST move with sigma — larger for the wider in-flight packet (sigma=1),
smaller for the narrower one (sigma=3) — while the asymmetry ratio keeps tracking
density. If ~2.15 survives unchanged across all sigma, the form-factor
explanation is WRONG and something else is responsible. **Compute this residual
first when the sweep lands.**

### Microscopy notebook (build details)

`hypotheses/bulk_ks_stopping_microscopy/` — `microscopy.py` (kernels) +
`build_microscopy_notebook.py`. One notebook, one velocity, a section per
density: decomposed energies + component changes, projectile state (KS orbital
width vs fixed classical charge), kinematics, and the induced-density **wake in
the projectile frame**, then a cross-density synthesis of three ratios
(S, peak depletion, front/back asymmetry).

**Essential decomposition:** `density_delta` IS the bath response for the
classical run (projectile is an external potential, never in n) but is NOT for
the WP run (the packet is an occupied orbital and dominates its own delta). WP
panels therefore use n_bath = n_total − n_wp, minus its t=0 value. Comparing raw
deltas would show the WP inducing a far larger response — pure artefact.

**Prototype result (classical, r_s = 3.987):** Δn at the projectile is NEGATIVE
(it sits in a hole it made); asymmetry is negative throughout (more depletion
BEHIND than ahead — the lag that produces drag) and deepens as the projectile
slows, consistent with the 26 % window effect above. **Peak depletion −1.9e-3
against n₀ = 3.77e-3 — the gas is locally emptied by ~50 %, so the classical
projectile is a STRONGLY NON-LINEAR perturbation**, outside the regime where
Lindhard-type linear response applies. A smeared wavepacket of the same charge
sits much closer to linear. *Inference: the two projectiles may be in different
response regimes entirely, not merely different coupling strengths.*

---

## 2026-07-31 — r_s = 3.99 pair COMPLETE. Clearing hypothesis NOT supported.

Chain ran end-to-end unattended, all four jobs exit 0:

| Job | Stage | Elapsed |
|---|---|---|
| 32439807 | GS (8/8 analytic gates PASSED) | 12:26 |
| 32439808 | wp (9/9 t=0 gates PASSED) | 1:09:10 |
| 32439810 | classical | 1:28:56 |
| 32439811 | both notebooks, 0 errors | 7:44 |

Faster than the r_s = 5.702 pair despite 2× the states — dx = 0.5 halved the grid
and more than paid for the density. 3.95 s/step (wp) vs 8.28 s/step before.

### HEADLINE — the density lever barely moved the gap

Both pairs fitted on the SAME window [4, 18.427] (the r_s = 5.702 run was
**re-fitted**, not quoted from its own wider window, so this is like-for-like):

| | S_cl (deposit) | S_cl (kinetic) | S₂ (WP drift) | S₁ (WP full) | **S_cl/S₂** |
|---|---|---|---|---|---|
| r_s = 5.702 | 0.3746 | 0.3740 | 0.0577 | 0.0145 | **6.49×** |
| r_s = 3.987 | 0.8899 | 0.8846 | 0.1574 | 0.1151 | **5.65×** |

**A 2.92× density increase moved the WP/classical ratio by 13 %** (6.49 → 5.65).
Both sides scaled together: classical 2.38×, WP 2.73×. The gap is close to
density-independent.

**Verdict on the clearing hypothesis (2026-07-31 entry below): NOT SUPPORTED as
the dominant cause.** Clearing is a screening-length effect; a 2.9× density rise
should have shrunk the ratio substantially. It did not. The ratio did shrink
slightly, and the WP scaled marginally more strongly than the classical, so a
small clearing contribution is not excluded — but it cannot be the main term.

Energy conservation holds throughout: deposit vs kinetic channels agree to 0.17 %
(old) and 0.59 % (new); WP `energy_total` drift 3.0e-4 eV; WP norm drift 1.0e-5;
Ehrenfest residual |s₃−s₄| ≤ 0.015 Bohr over a 38.3 Bohr path.

### SOLID SECONDARY RESULT — S₁ vs S₂ is a fixed offset, not bath scattering

```
r_s=5.702:  S2=0.0577   d(T1-T2)/ds=+0.0431 eV/Bohr   S2 - that = 0.0145 = S1
r_s=3.987:  S2=0.1574   d(T1-T2)/ds=+0.0423 eV/Bohr   S2 - that = 0.1151 = S1
```

**S₂ tripled; the momentum-width term moved 1.9 %.** If T₁−T₂ growth were
scattering off the electron gas it would scale with n. It does not — it is a
property of the packet alone.

*Inference (mine, NOT independently verified):* most likely the **self-interaction
error** — the WP is an occupied KS orbital, its own charge enters the Hartree
potential acting on it, and LDA has no exact exchange to cancel it. The packet
pushes itself apart at a rate set by its own density (σ = 2, identical in both
runs), with the bath playing no role. Plain free dispersion cannot explain it: a
free Gaussian keeps σ_p exactly constant.

**Consequence, and it is firmer than the r_s = 5.702 run alone could support:
S₂ is the defensible KS-orbital stopping power. S₁ should be reported as
"S₂ minus a spreading term", never as an independent measurement.**

This also explains the earlier r² collapse. The offset is a fixed ~0.043 eV/Bohr:
at r_s = 5.702 that was 75 % of S₂ (so S₁ ≈ 0, r² = 0.62, error bar > value); at
r_s = 3.987 it is 27 % of a tripled S₂, so S₁ emerges clean at r² = 0.993. The
definition did not improve — the signal outgrew a fixed artefact.

### S₂ behaves like a real stopping power

S₂ scaling 2.73× against a 2.92× density ratio. Linear-in-n predicts 0.169,
measured 0.157 — 7 % below, the direction and size linear response expects for a
fast projectile as the Lindhard logarithm weakens with rising k_F.

### OPEN — and the dilution explanation has a problem

The natural fallback for the residual 5.65× is charge dilution: σ_z spreads
2.16 → 8.33 Bohr (3.9×) while the classical projectile holds a fixed 1.41 Bohr.
That is density-independent, which fits.

**But it predicts curvature the data does not show.** If drag fell as the packet
spread, T₂(s) would be steep early and flatten. Instead S₂ fits a straight line at
**r² = 0.997** across a window where σ nearly quadruples. Constant drag despite a
3.9× width change. Simple dilution does not survive that.

**NEXT EXPERIMENT — σ-sweep at fixed energy and fixed density** (σ ∈ {0.5, 1, 2, 4}).
If the ratio collapses as σ → small, width is the cause and the linearity needs
explaining. If the ratio is σ-independent too, the difference is intrinsic to
representing the projectile as a KS orbital, which is the more interesting result.
The density sweep cannot separate these; do NOT repeat it.

### Artefacts

- `hypotheses/bulk_ks_stopping_rs4/bulk_ks_stopping_rs4_{wp,classical}.ipynb`
  (52 MB / 34 MB, 0 errors, GIFs + figures embedded)
- `hypotheses/bulk_ks_stopping_rs4/figures/four_stopping_powers.png`
- `hypotheses/bulk_ks_stopping_rs4/{wp,classical}_stopping_summary.json`

---

## 2026-07-31 — r_s = 3.99 DENSITY-REPLICA pair SUBMITTED (autonomous chain)

Tests the density-clearing hypothesis below. **Everything about the projectile is
held fixed; density is the only variable that moves.**

### Jobs (submitted 2026-07-31 ~03:5x, queued behind a running wp-hd array)

| Job | Stage | Depends on | Wall |
|---|---|---|---|
| 32439807 | GS `gs_L40x40x80_orth_N482_dx0p50` | — | 4 h |
| 32439808 | wp | afterok:32439807 | 12 h |
| 32439810 | classical | afterok:32439807 | 12 h |
| 32439811 | both notebooks | afterok:32439808:32439810 | 3 h |

`afterok` throughout: the GS exits 3 on any failed analytic gate, so a bad
ground state stops the chain rather than seeding two production runs.

### Configuration and why

| | r_s = 5.702 pair (done) | r_s = 3.99 pair (this) |
|---|---|---|
| Cell | 46 × 46 × 80 | **40 × 40 × 80** (L_z unchanged) |
| N | 218 | **482** (even, closed shell) |
| n | 1.2878e-3 | **3.7656e-3** (2.92× denser) |
| dx | 0.40 | **0.50** |
| grid | 115×115×200 = 2.65 M | 80×80×160 = 1.024 M |
| states | 130 | 262 |
| σ_WP, E, z₀, dt, N_STEPS | 2, 100 eV, −32, 0.04, 646 | **identical** |
| ħω_p | 3.46 eV | 5.92 eV |
| plasma periods in run | 0.52 | **0.89** |
| fit window | [4, 18.97] | **[4, 18.43]** |

**r_s = 3 was rejected as infeasible, not as undesirable.** In the original
46×46×80 box it needs N = 1514 → 778 states → **~99 GB** of orbitals for ETRS's
three copies, against an 80 GB A100. Even 42×42×80 lands at ~74 GB (8 % headroom).
r_s ≈ 4 in a 40×40×80 box needs ~13 GB — comfortable. Recorded so nobody
re-proposes r_s = 3 at this box size.

**dx = 0.50 was VERIFIED, not assumed.** WP momentum moments computed on the exact
80×80×160 grid reproduce analytic values to **machine precision** (⟨p_z⟩, T₁ and
T₁−T₂ all zero error): the packet's spectral content reaches k = 3.461 against
Nyquist 6.283 (1.8× margin), bath k_F = 0.481 (13×), and the classical UPF form
factor is ~1e-17 at the Nyquist edge. **Deviation to declare:** every other jellium
run used dx = 0.40. The WP/classical ratio here is internally consistent (both
halves share the grid); a cross-pair comparison of *ratios* carries a second-order
grid difference.

**The transverse constraint now BINDS (inverted vs the L=46 pair).** T_transverse
= 18.43 < T_IFW = 18.97, so periodic images in x/y end the clean window, not the
+z face. FIT_T1 takes the min either way, so the analysis is correct, but the
window is 2.8 % shorter. To compare the two pairs strictly, re-fit the r_s = 5.702
run on [4, 18.43]. A `static_assert` pins the inversion so widening L_xy trips it.

### Analytic GS gates (exit 3 if any fail)

E_kinetic = **33.593666 Ha** (exact discrete plane-wave sum, ±2.0e-2);
E_hartree = 0 (±1.0e-2); E_total in [−38.30, −36.06] (PW92 prediction −37.179,
and the same model was 0.11 % low on the r_s = 5.702 run); n_occupied = 241;
occupations 2.0/0.0; gap 0.252 eV.

### Files

- `shared/configs/bulk_ks_stopping_L40x40x80_rs4.hpp`
- `save_gs/gs_L40x40x80_orth_N482_dx0p50/run.cpp`
- `scripts/bulk_ks_stopping_rs4/{wp,classical}/run.cpp`
- `hypotheses/bulk_ks_stopping_rs4/build_run_notebook.py`
  (`ks_stopping.py` and `tests/` are **symlinks** to the r_s = 5.702 folder — one
  implementation, so both pairs get bit-identical arithmetic)
- `shared/bin/run-bulk-{gs-rs4,ks-stopping-rs4,notebooks-rs4}.slurm`

Fixed in the rs4 copy only (the original is provenance for a completed run): the
t=0 gate label read `sigma_pz^2 = 1/(4 sigma^2)` while checking the correct
`1/(2 sigma^2)` value — a leftover from the 2026-07-30 momentum-width error.

### What to look at when it lands

The headline is the **ratio** S_classical / S_WP(drift), which was **6.6×** at
r_s = 5.702. If clearing drives the gap, higher density should shrink it markedly.
If it is roughly unchanged, the cause is the KE definitions or free-packet
dispersion instead — and a σ-sweep, not a density sweep, is the next lever.

---

## 2026-07-31 — OPEN HYPOTHESIS: initialisation density-clearing may explain the WP/classical gap

**Status (updated 2026-07-31, after the r_s = 3.99 pair): TESTED via the density
lever — NOT SUPPORTED as the dominant cause. A 2.92× density increase moved the
WP/classical ratio only 6.49× → 5.65×. See the r_s = 3.99 COMPLETE entry above
for the numbers and for what to try next (σ-sweep, not another density sweep).
The original hypothesis is kept below verbatim as the record of what was asked
and why.**

### The observation

| | S (eV/Bohr), same window |
|---|---|
| classical, bulk r_s = 5.70 | 0.377 |
| WP S₂ (drift), bulk r_s = 5.70 | 0.057 |
| ratio | **6.6×** |

The classical projectile shows a clear, clean energy decay. The wavepacket decays
too, but far less.

### The user's hypothesis (2026-07-31, verbatim intent)

> "This might be due to the clearing of density close to the wavepacket location.
> From preliminary analysis, the same comparison of the wavepacket and classical
> stopping power in the jellium slab cases is much closer. This suggests to me
> that the clearing of local density due to the initialisation step causes the
> problem. Need to check this more carefully."

So: **the WP injection step is suspected of evacuating bath density from the
projectile's neighbourhood**, leaving the packet travelling through a locally
rarefied medium and therefore under-reporting the drag. The **slab** runs
(`localised_jellium`) reportedly show a much smaller WP/classical discrepancy,
which is the comparative evidence pointing at initialisation rather than at the
stopping definitions themselves.

### Relevant data already in hand (context, not a verdict)

- `max_overlap = 3.12e-06` for this WP injection. The packet sits at k₀ = 2.711
  while the occupied manifold is confined to |k| ≤ k_F = 0.337, so **Gram–Schmidt
  orthogonalisation had almost nothing to remove**. *Inference (mine, unverified):
  if clearing is real here, orthogonalisation of the WP against the occupied
  manifold is unlikely to be the mechanism — look instead at the Hartree/ALDA
  self-interaction of the added orbital, which is not cancelled in LDA and pushes
  bath density away from the packet from the first step.*
- The bulk WP and classical halves both start from the **same bare jellium GS**
  (`gs_L46x46x80_orth_N218_dx0p40`), so neither begins with a pre-formed screening
  hole. The asymmetry is therefore not "one was pre-screened and the other wasn't".
- r_s = 5.70 is a **very dilute** bath. Screening length is long and the density is
  low, so any evacuated hole is a large fraction of the local density — the effect,
  if real, should be strongly density-dependent. This is what the r_s ≈ 3 twin pair
  below is designed to test.

### The concrete test (defined, NOT run)

1. From the saved density VTIs of the WP run, form the **bath** density
   n_bath = n_total − n_WP and plot n_bath in a shell around the packet centroid
   vs time. A depression at t → 0 that persists as a co-moving hole is the
   signature. Baseline: the uniform GS value n₀ = 1.2878e-3.
2. Compare against the classical run's induced density at the same path length —
   there the hole is physical (screening) and should be *deeper*, not shallower.
3. **Density lever (the r_s ≈ 3 pair, below):** if clearing drives the gap, raising
   the density shortens the screening length and should shrink the WP/classical
   ratio markedly. If the ratio is roughly density-independent, the cause lies in
   the KE definitions or in dispersion instead.
4. Cross-system: quantify the slab WP/classical ratio properly
   (`localised_jellium/hypotheses/wp_highdensity_sv`) rather than relying on the
   preliminary read.

**Competing explanation not yet excluded:** free-packet dispersion. σ_d grows
1.41 → 7.2 Bohr over the flight, so the WP's charge is smeared over a volume
growing as σ_d³ — a diluted perturbation couples more weakly regardless of any
initialisation artefact. A σ-sweep at fixed energy separates this from clearing;
the density sweep does not.

---

## 2026-07-31 — WP notebook: four-S comparison plot + skill cross-check

**Trigger.** User asked for the four stopping powers to be plotted together in the
WP notebook, extracted via the `stopping-power-extraction` skill, alongside the
per-definition KE(t) and s(t) curves.

**State before.** `bulk_ks_stopping_wp.ipynb` already had the per-definition curves
(cell 9: T₁, T₂, T₁−T₂ vs t; cell 12: s₃, s₃_naive, s₄, residual vs t), the 2×2
per-fit panel and the numeric table. **Missing:** a single figure comparing the four
S values, and any external reference to measure them against.

**Added** (in `build_run_notebook.py`, the reproducible source — the `.ipynb` is
generated, never hand-edited):

- New §5b, WP branch only: markdown + 2 code cells + reading notes.
- Notebook now imports the skill's kernels directly
  (`.claude/skills/stopping-power-extraction/stopping_power.py`) rather than
  re-implementing the fit, so this run is measured with the same code as every
  other run in the project.
- Figure `hypotheses/bulk_ks_stopping/figures/four_stopping_powers.png`
  (2 panels: ΔT(s) slopes overlaid | S values with error bars + classical band).
- `wp_stopping_summary.json` gains a `classical_reference` block.

### The methodological point (worth not re-deriving)

The skill's locked default for continuous traversal is **Method A** — slope of
ΔE_total(x) — with −dT/dx demoted to a cross-check. **That inverts for the WP half.**
The wavepacket is an occupied KS orbital, so it sits *inside* `energy_total`; the
system is closed and E_total is constant to **2.6e-4 eV** across the whole run.
Method A's fit target is identically zero — the method is **undefined**, not
imprecise. The −dT_i/ds_j slopes are the measurement; the classical twin supplies
the external reference. Recorded in plan §4.1.

### Classical reference (skill kernels, computed 2026-07-31)

| Channel | S (eV/Bohr) | r² |
|---|---|---|
| Method A, skill default (fixed 20% time cut) | 0.41 | 0.998 |
| Method A on the WP window [4, 18.97] | **0.38** | 0.994 |
| Sanity −dKE_ion/dx, same window | 0.38 | 0.994 |

Deposit vs kinetic channels agree to **0.18%** (independent CSVs) → conservation
confirmed. Ratios against the WP: classical/S₂ = **6.6×**, classical/S₁ = **24.5×**.

### Bug found and fixed — `norm_check` collision

`ks_stopping.py` merged the momentum and real-space stats files with suffixes
`_p`/`_r`; **both** carry a `norm_check` column meaning different things, and
`load_wp_run` took the momentum one. The notebook then printed
`norm = 48998109.640800 (expect 1)` — a healthy run displayed as catastrophically
broken. The momentum column is an unnormalised **Parseval sum over the FFT grid**
(~4.9e7), not a probability norm.

Fixed: `WPRun.norm` now takes `norm_check_r` (real space, = 1); the Parseval sum is
kept as a separate `WPRun.parseval` field and reported as a *constancy* diagnostic.
This is the second time this exact constant has misled — the first was the wrong
`wp_momentum_stats.hpp` docstring (see 2026-07-30 entry). Both are now documented
at the point of use.

11/11 tests in `hypotheses/bulk_ks_stopping/tests/` still pass. Verified in the
rebuilt notebook: `norm (real) = 1.000000`, drift 1.1e-5 over 646 steps; Parseval
sum 4.8998e7, varying by 1.1e-3 % (constant, as it should be).

### Rebuild verification (2026-07-31, 02:40)

`bulk_ks_stopping_wp.ipynb` — **28 cells, 0 errors**, 48 MB, 9 density GIFs
embedded, 5 figures. New §5b is cells 20–23. `wp_stopping_summary.json` now
carries the `classical_reference` block. **The classical notebook was NOT
rebuilt** — §5b is WP-only, and nothing in the classical branch changed.

Final numbers, fit window t ∈ [4.0, 18.97] a.u., n = 375:

| | KE | position | S (eV/Bohr) | r² |
|---|---|---|---|---|
| S₁₃ | ⟨p²⟩/2m | centroid | 0.015 ± 0.016 | 0.66 |
| S₁₄ | ⟨p²⟩/2m | ∫⟨p⟩dt | 0.015 ± 0.016 | 0.66 |
| S₂₃ | ⟨p⟩²/2m | centroid | **0.057 ± 0.004** | 0.999 |
| S₂₄ | ⟨p⟩²/2m | ∫⟨p⟩dt | **0.057 ± 0.004** | 0.999 |

**Read this as TWO numbers, not four.** S₁₃≈S₁₄ and S₂₃≈S₂₄ to 3 d.p. — the
Ehrenfest identity carried through to the answer, so the position axis is a
consistency check that passed, not a degree of freedom.

**S₁ is not distinguishable from zero** (uncertainty 0.016 > value 0.015) and
ΔT₁(s) is visibly curved — it RISES for the first ~20 Bohr before turning over,
hence r²=0.66. Quote S₂ = 0.057 as the headline drift stopping; report S₁ only
with the curvature caveat. The rise means momentum-space broadening is outrunning
drift loss early on: energy moving from the drift channel into the spread
channel, part real scattering and part free-particle dispersion.

### Note on the user's pasted configuration

The request included a config block (35×35×85 cell, two sin² CAP bands η=−1.0 Ha,
25 Bohr slab, σ_WP=0.5, launch z=−24, v∈{2.0,2.5,3.0,3.5}, `inq-study` engine).
That is **`localised_jellium/scripts/wp_highdensity_sv`**, not this bulk run — this
run is periodic, CAP-free, σ_WP=2, single 100 eV energy. Treated as a mis-paste and
**not** applied. Flagged to the user; if those slab runs are meant to receive the
same four-definition treatment that is a separate piece of work.

---

## 2026-07-30 (evening) — BOTH RUNS COMPLETE; stopping powers extracted

### Runs

| Job | Half | Wall | Steps | s/step | Exit |
|---|---|---|---|---|---|
| 32401711 | wp | 1h35m | 646 | 8.28 (5.67 compute / 10.88 write) | 0 |
| 32401322 | classical | 2h09m | 646 | 11.19 (10.22 compute / 12.15 write) | 0 |

Ran CONCURRENTLY on separate A100s (gpu-q-28, gpu-q-45). **The classical half is
the slower one** — it inserts a pseudo-ion, so it pays non-local projector
application plus an Ehrenfest `forces_stress` evaluation every step, neither of
which the ion-free wavepacket run does. Size future twin pairs off the classical
half. Outputs: wp 40 GB, classical 25 GB.

Job 32401321 (first wp attempt) FAILED exit 4 on a wrong t=0 gate — see below.

### Results — stopping power (eV/Bohr, fit window t in [4.0, 18.97] a.u.)

Uncertainties are stat (OLS) and syst (window edges moved +/- 3 a.u.) in
quadrature; the systematic dominates throughout.

| Definition | S | r^2 | significance |
|---|---|---|---|
| S_cl classical, shared window | **0.377 +/- 0.045** | 0.994 | 8.4 sigma |
| S_cl initial drag (v >= 0.85 v0) | 0.365 +/- 0.007 | 0.985 | |
| S_23 / S_24  = -d(<p>^2/2m)/ds | **0.057 +/- 0.004** | 0.999 | **14.6 sigma** |
| S_13 / S_14  = -d(<p^2>/2m)/ds | **0.015 +/- 0.016** | 0.658 | **consistent with ZERO** |

Classical projectile: z -32 -> +34.55 Bohr, v 2.7111 -> 2.3854 (88 % of v0),
KE lost 22.58 eV over 66.5 Bohr.

**PHYSICAL READING (stated carefully — an earlier draft of this file over-claimed
a "2:1 split", which the uncertainties do not support).** S_1 has a systematic
(0.016) LARGER than its central value (0.015), so the wavepacket's TOTAL kinetic
energy is CONSERVED within uncertainty: this run resolves no net energy transfer
to the bath from the <p^2> channel at all. Meanwhile the drift-momentum loss S_2
is significant at 14.6 sigma. The defensible statement is therefore:

> **the packet scatters almost elastically** — momentum is redistributed out of
> the drift component into the momentum spread, with little or no energy actually
> delivered to the bath that this run can resolve —

against the classical twin's unambiguous 0.377 eV/Bohr. The differing r^2 is part
of the same story: T2(s) is a clean line (drift decays smoothly, r^2 = 0.999)
while T1(s) is flat-and-noisy (r^2 = 0.658) because there is no trend to fit.

**Candidate cause of the quantum/classical gap** — the packet disperses (density
width 1.41 -> ~7 Bohr across the fit window, ~100x the volume), diluting the
charge that drives the induced response, while the classical Gaussian stays rigid
at 1.41 Bohr. **HYPOTHESIS, not demonstrated.** A sigma-sweep at fixed energy
would settle it. The notebooks plot measured sigma_z against the free-dispersion
law so the confinement question can at least be inspected.

**The position definitions agree to SIX significant figures**: S_13 vs S_14 differ
by 1.4e-6 eV/Bohr, S_23 vs S_24 by 6.9e-6. Two completely different routes (a
real-space phase-operator expectation vs a time-integral of an FFT moment) landing
on the same number. So the four combinations are really TWO: the position axis
collapses exactly as the Ehrenfest argument predicted before the runs, and all the
physical content is in T1 vs T2.

### Validation outcomes

| Check | Result |
|---|---|
| WP energy conservation (closed system) | drift **-0.0003 eV** over 25.8 a.u. — excellent |
| Classical electron/ion energy closure | electrons +22.631 eV vs projectile -22.582 eV -> **0.049 eV (0.22 %)** |
| Ehrenfest identity s3 == s4 | max abs diff **0.0055 Bohr** over a 40 Bohr path |
| WP orbital norm drift | -1.1e-5 relative — no leakage into the bath |
| circular vs naive centroid | **differ by up to 33.8 Bohr** |

**The circular centroid was decisive, not cosmetic.** Had definition 3 used the
naive integral of z|psi|^2, it would have been wrong by up to 33.8 Bohr and
S_13/S_23 would have been noise. Its correctness is independently confirmed by
the Ehrenfest identity holding to 0.0055 Bohr against s4, which is computed by a
completely different route (integrating <p_z>).

**`energy_ion_kinetic` is NOT populated by INQ** (all zeros; verified). So in the
classical run `energy_total` is the ELECTRONIC energy only and is SUPPOSED to
rise. Treating its increase as numerical drift would be a mistake — the correct
test is closure against the projectile's KE loss. `ks_stopping.conservation_check`
now takes `projectile_ke_loss_ev` and reports the closure; the notebook branches
on the run half.

### SECOND wrong-expectation incident (same class as the Rydberg one)

Job 32401321 aborted on 4 of its own t=0 gates. **The gates were wrong; the
wavepacket was correct.** I had used sigma_p = 1/(2 sigma), giving
sigma_p^2 = 1/(4 sigma^2) = 0.0625. The correct value is sigma_p^2 =
1/(2 sigma^2) = 0.125.

Decisive check: the real-space density std is sigma/sqrt2 = 1.41421 (which the
run reported exactly), so my value implied sigma_d * sigma_p = 0.354 < 1/2 —
a **violation of the Heisenberg bound**. A Gaussian is minimum-uncertainty and
must give exactly 1/2. With the correct sigma_p^2 every measured quantity matches
analytically to 5 decimals.

**Root cause: the "Known-case validation" docstring in
`inq-stack/include/inqkit/observables/wp_momentum_stats.hpp` was WRONG** and had
been since it was written (it claimed sigma_p = 1/(2 sigma_r), sigma^2_pd =
1/(4 sigma_r^2)). The CODE was always right. I copied the docstring into a gate
instead of deriving it. The header is now corrected, with the Heisenberg argument
included so the error is self-evident to the next reader, plus a note that the
momentum-space Parseval constant N is ~5e7 and must NOT be gated on as ~1.

Consequence for the physics: the localisation energy separating the two KE
definitions is **3/(4 sigma^2) = 5.102 eV** at sigma = 2, not the 2.551 eV
quoted earlier — twice the signal. No design decision changes.

**Standing lesson (now twice): never take an analytic expectation from a
docstring or from memory. Cross-check it against a physical invariant** — a
conservation law, the uncertainty bound, or a functional-independent sum — before
it becomes a gate. Both incidents were caught by gates before wasting GPU time,
which is the system working, but both cost a queue round-trip.

### Deliverables

- `hypotheses/bulk_ks_stopping/bulk_ks_stopping_wp.ipynb` — 24 cells, **0 errors**,
  9 embedded GIFs + 5 figures, 48 MB
- `hypotheses/bulk_ks_stopping/bulk_ks_stopping_classical.ipynb` — 17 cells,
  **0 errors**, 3 embedded GIFs + 3 figures, 35 MB
  (first attempt had 3 identical errors: `\frac12` is invalid matplotlib
  mathtext — needs `\frac{1}{2}`. Only affected an axis LABEL; `\tfrac12` in
  markdown cells is fine, MathJax accepts it.)
- `{wp,classical}_stopping_summary.json` — machine-readable results
- Density-matrix GIFs in `scripts/bulk_ks_stopping/<half>/results/report/`
  (9 for wp: {total,wp,bath} x {density,delta0,dstep}; 3 for classical)

### Next steps

1. Cross-run comparison notebook in `hypotheses/bulk_ks_stopping/`.
2. Run the W1 circular-centroid engine test on a GPU node (still NOT run — the
   feature is nonetheless validated in production by the Ehrenfest identity).
3. Open question for the user: is the 19x quantum/classical gap physical
   (dispersion diluting the perturbation) or an artefact of the packet width? A
   sigma-sweep at fixed energy would settle it.

---

## 2026-07-30 — ground state VALIDATED; both production runs launched

### Where it stands

| Item | Status |
|---|---|
| Design locked with user | DONE |
| W1 circular centroid in `WPRealSpaceStats` + test | Code DONE; **test NOT yet compiled/run** |
| W2 dispersion-aware IFW helper | DONE, `static_assert`s compile and pass |
| W3 Gaussian UPF sigma_pot = 1.4142 | DONE, validated to 5e-11 against analytic erf |
| W4 config header | DONE, all geometry `static_assert`s pass |
| W5 ground state | **DONE and PHYSICALLY VALIDATED** (job 32400615) |
| W6 both `run.cpp` | Written; **NOT yet compiled** — first compile happens inside jobs 32401321/32401322 |
| W7 SLURM scripts | DONE |
| W8 `ks_stopping.py` engine + 11 known-case tests | DONE, all pass |
| W9 notebook builder | Written; **not yet exercised on real data** |
| Production runs | **RUNNING** (submitted 2026-07-30 ~14:55) |

### Locked configuration

46 x 46 x 80 Bohr orthorhombic, periodic in x,y,z. N = 218 electrons
(r_s = 5.702, hbar*w_p = 3.462 eV, v_F = 0.3366). dx = 0.40 requested; **INQ
actually chose a 120 x 120 x 200 = 2.88 M grid** (dx_xy = 0.383, dx_z = 0.400).
129 states = 109 occupied + 20 extra. sigma_WP = 2 Bohr, E = 100 eV,
k0 = v = 2.7111, launch z0 = -32. dt = 0.04, N_STEPS = 646 (t = 25.84 a.u.).
Fit window t in [4.0, 18.97] a.u. = steps 100-474 = 40.6 Bohr of path.

### Ground state — validated against the analytic plane-wave result

Uniform jellium has no ions, so the KS potential is constant, the orbitals are
exactly plane waves and every energy term is predictable in closed form. Job
32400615, 17m48s, 84 SCF iterations at ~7.8 s/iter:

| Quantity | SCF | Analytic | Agreement |
|---|---|---|---|
| E_kinetic | 7.396581 Ha | 7.395907 Ha (exact discrete sum) | **6.7e-4 Ha** |
| E_total | -15.848701 Ha | -15.830 Ha (PW92) | **0.12 %** |
| E_hartree | 2.74e-4 Ha | 0 | 7 meV residual |
| E_xc | -23.245555 Ha | -23.226 Ha | 0.08 % |
| HOMO-LUMO gap | 0.25182 eV | 0.24354 eV | 0.008 eV |
| min occ (filled) | 1.9999988 | 2.0 | closed shell CONFIRMED |
| max occ (empty) | 2.26e-6 | 0.0 | |
| n_electrons / n_occupied | 218 / 109 | 218 / 109 | exact |

Checkpoint (5.6 GB):
`/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/ResearchProject/systems/jellium/checkpoints/gs_L46x46x80_orth_N218_dx0p40`

**Inverting the SCF identifies the functional**: (E_total - E_kinetic)/N gives
eps_xc = -0.106627 Ha, whose implied eps_c = -0.026274 Ha matches PW92's
-0.026189 to 0.3 %. So INQ's `options::theory{}.lda()` correlation is PW92-like.
This was not known before and is worth reusing.

### MISTAKE MADE AND CORRECTED — read this before writing another jellium gate

The GS job **exited 3 and its `afterok` dependants (jobs 32400779, 32400780)
were CANCELLED**, purely because of two bad gate constants I wrote:

1. **Units error (the important one).** I predicted E_total = -33.37 Ha using
   eps_x = -0.9163/r_s. **That expression is in RYDBERG.** In Hartree it is
   eps_x = -0.4582/r_s. The same trap applies to the PZ81 correlation constants
   (gamma = -0.1423 **Ry** = -0.07115 Ha). The factor of 2 in exchange made the
   predicted total roughly twice too negative. Corrected, the model reproduces
   the SCF to 0.12 %.
2. **Tolerance too tight.** The E_hartree gate demanded 0 +/- 1e-4 Ha; a
   converged-but-not-exact SCF leaves a 2.74e-4 Ha ripple (1.7e-5 of |E_total|).
   Genuine symmetry breaking would be 0.1-1 Ha. Tolerance widened to 1e-2 Ha.

Both constants are now fixed in
`save_gs/gs_L46x46x80_orth_N218_dx0p40/run.cpp`, with the units warning written
into the header block so it cannot be lost. **The ground state itself was never
wrong and was NOT recomputed** — the checkpoint from job 32400615 is the one in
use.

Lesson worth generalising: gate hard on quantities that are *unit-safe and
structural* (the discrete kinetic sum, occupations, electron count, the gap) and
treat functional-dependent totals as soft/reported. The E_kinetic gate is the
one that actually proved the ground state was right.

### Running now

| Job | Half | Submitted | Notes |
|---|---|---|---|
| 32401321 | `wp` | 2026-07-30 ~14:55 | no dependency — GS already validated |
| 32401322 | `classical` | 2026-07-30 ~14:55 | |

Both submitted WITHOUT an `afterok` dependency, deliberately: the GS checkpoint
exists and has been validated by hand, so re-running the GS just to satisfy a
now-fixed gate would waste 18 min of A100 time.

**Neither `run.cpp` has ever been compiled.** The first compile happens inside
these jobs (~10 min of the wall time). If they fail, it will be a compile error
in `scripts/bulk_ks_stopping/{wp,classical}/run.cpp`, not a physics problem —
check the job log for `error:` before assuming anything else.

Watch: `bulk-ks-32401321.out`, `bulk-ks-32401322.out` in the repo root.

### What each run must produce

WP (`scripts/bulk_ks_stopping/wp/results/`):
`raw/observables/wp_momentum_stats.csv` and `wp_real_space_stats.csv` (EVERY
step — these are the measurement), `observables.csv` (full energy decomposition),
`raw/vti/{density_total,density_wp,density_delta,wavefunction_wp}`,
`checkpoint/`, `rt_state.txt`, `run_summary.txt`.

Classical (`scripts/bulk_ks_stopping/classical/results/`): `electron_track.csv`
(z, v, KE every step), same energy decomposition, `density_total`/`density_delta`.

Both abort with exit 4 if their t=0 analytic gates fail (WP: `<p_z>` = k0,
sigma_pz^2 = 0.0625, T1 = 3.76864 Ha, T1-T2 = 2.551 eV, centroid = -32;
classical: mass 1.0, KE = 3.674905 Ha, launch z = -32).

### Next steps (autonomous — user asked for no further approval gates)

1. Wait for 32401321 / 32401322. A background watcher is polling them.
2. On success: `venv/bin/python ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/build_run_notebook.py wp` and `... classical`.
3. Run the W1 circular-centroid engine test (needs a GPU node; not yet done).
4. Cross-run comparison in `hypotheses/bulk_ks_stopping/`.

### Known limitations to carry into every downstream claim

- **2*pi/w_p = 49.4 a.u. exceeds the 25.8 a.u. run.** The bath never completes one
  plasma oscillation, so no steady wake forms. S here is an **initial-drag**
  stopping power, NOT a converged steady-state S(v), and must not be compared
  with Lindhard/Bethe as though it were. Forced by geometry: a light 100 eV
  electron crosses this box in ~26 a.u.
- s3 and s4 are related by an exact Ehrenfest identity (no ions, no CAP), so
  their agreement is a validation, not a second physics channel. The real
  contrast is T1 vs T2.
- `boundary_rule.hpp`'s legacy static `ifw_end_z` over-estimates the clean window
  by 22 % for this run (24.34 vs 18.97 a.u.). Use `ifw_end_t_dispersive`.

### Files added this session

- `docs/plans/bulk-jellium-ks-stopping.md`
- `inq-stack/include/inqkit/observables/wp_real_space_stats.hpp` (MODIFIED: +9 circular columns)
- `inq-stack/tests/include/inqkit/observables/test_wp_circular_centroid_engine.cpp` (+ CMake entry)
- `inq-stack/tests/python/inqview/io/test_gaussian_psp.py` (FIXED: hard-coded old-device path made all 8 tests silently skip since the migration)
- `ResearchProject/systems/jellium/shared/configs/boundary_rule.hpp` (MODIFIED: dispersive IFW helpers)
- `ResearchProject/systems/jellium/shared/configs/bulk_ks_stopping_L46x46x80.hpp`
- `ResearchProject/systems/jellium/shared/cpp/rt_state.hpp`
- `ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma2p0.upf`
- `ResearchProject/systems/jellium/save_gs/gs_L46x46x80_orth_N218_dx0p40/run.cpp`
- `ResearchProject/systems/jellium/scripts/bulk_ks_stopping/{wp,classical}/run.cpp`
- `ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/{ks_stopping.py,build_run_notebook.py,tests/test_ks_stopping.py}`
- `shared/bin/run-bulk-gs.slurm`, `shared/bin/run-bulk-ks-stopping.slurm`
- venv: installed `nbformat`, `nbclient`, `ipykernel` (were missing; login node has network, compute nodes do not)

### Validation status

| Check | Status |
|---|---|
| `ks_stopping.py` 11 known-case tests | **PASS** |
| `test_gaussian_psp.py` 8 tests | **PASS** (were silently skipping) |
| Gaussian UPF vs analytic erf | **PASS**, 5e-11 |
| `boundary_rule.hpp` static_asserts | **PASS** |
| config header static_asserts | **PASS** |
| Ground state vs analytic plane-wave | **PASS** (see table above) |
| W1 circular-centroid engine test | **NOT RUN** — needs a GPU node |
| Both production `run.cpp` compile | **NOT VERIFIED** — first compile is in the running jobs |

## 2026-08-02 — disk cleanup: raw VTI frames purged (four GIF'd sets)

- Deleted `results/raw/vti/` for `bulk_ks_stopping`, `bulk_ks_stopping_rs4`,
  `bulk_ks_stopping_sigma3`, `bulk_ks_stopping_rs4_sigma3` ({wp,classical} each,
  ~98 GB total). User-approved 2026-08-02. The sigma3 deletion is post-re-run:
  its phase GIFs were rebuilt 2026-08-01 19:45, after the interactions re-run's
  last VTI (05:59), so the banked analysis reflects the re-run data.
- NOT deleted: `bulk_ks_stopping_sigma1` and `bulk_ks_stopping_rs4_sigma1` VTIs
  (~45 GB) — their density-GIF battery was never built; deleting would make that
  permanent. Decision pending with the user.
- Final `checkpoint/` dirs untouched everywhere; runs remain extendable.
- Consequence: notebooks/GIFs in `hypotheses/bulk_ks_stopping*/` keep embedded
  outputs but cannot be rebuilt from raw fields without a re-run.

## 2026-08-03 — sigma1 VTIs purged (pending decision resolved)

The user approved immediate deletion of the `results/raw/vti/` trees of
`bulk_ks_stopping_sigma1/{wp,classical}` and
`bulk_ks_stopping_rs4_sigma1/{wp,classical}` (~47 GB) WITHOUT first building
their density-GIF battery — that battery is now permanently forfeited (the
2026-08-02 caveat no longer pending). Stopping results, observables CSVs,
`*_stopping_summary.json`, notebooks, and all final checkpoints are untouched.
Zero `.vti` files remain under either set (verified).
