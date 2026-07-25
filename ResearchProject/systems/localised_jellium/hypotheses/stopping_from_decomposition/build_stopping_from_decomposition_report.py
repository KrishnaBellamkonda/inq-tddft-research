#!/usr/bin/env python3
"""Builder for the `stopping-from-energy-decomposition` campaign spine notebook.

Phase 1 (human-gated) content ONLY:
  §A  enumerate every energy recorded in a run (lumped stores + pairwise ledger
      + projectile track), grouped by CSV and by run kind (WP vs classical);
  §B  numerically VERIFY the WP/classical closure relations on the
      `twin_ec_rsweep` r12 reference pair (+ the long phase5_wp WP run);
  §C  PROPOSE the ranked S-from-decomposition formulae for the user to approve.

Phases 2-4 (implement+validate the kernel; apply to the select plateaued+
decomposed runs; aggregate S(E0)/S(v0)) are appended by later builder passes.

Run:
    PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
    /local/data/public/skcb2/tddft/venv/bin/python3 \
        build_stopping_from_decomposition_report.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # hypotheses/ (for _nbreport)
import _nbreport as R

OUT = os.path.join(HERE, "stopping_from_decomposition.ipynb")
R.set_outdir(HERE)

# Reference runs that carry the pairwise ledger (interactions.csv), for §A/§B.
LJD = ("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
       "/scripts/localised_jellium_dynamics")
REF_CLASSICAL = f"{LJD}/runs/twin_ec_rsweep/results/r12_classical/raw/observables"
REF_WP = f"{LJD}/runs/twin_ec_rsweep/results/r12_wp/raw/observables"
REF_WP_LONG = f"{LJD}/phase5_wp/results/p5_null_s2_k4_wp/raw/observables"

cells = []

# ------------------------------------------------------------------ 1. Title
cells += [R.md(
    "# Stopping power from the decomposed energy ledger — classical & wavepacket\n"
    "### Campaign `stopping-from-energy-decomposition` · localised jellium slab\n\n"
    "**The question.** Electronic stopping power *S* is the energy a projectile "
    "transfers to the target electrons per unit path length. We already have a "
    "headline **Definition 2** (`S = E_absorbed / L_slab`, the localised-slab "
    "energy deposit). This notebook derives and applies **Definition 1** — a "
    "stopping power built directly from the *decomposed energy ledger* (the "
    "pairwise Coulomb terms E_PP/E_PS/E_SS/E_SB/E_PB/E_BB, the kinetic split, and "
    "E_xc) — for **both** a classical Gaussian-charge projectile **and** a "
    "wavepacket (WP) projectile, so the two can be compared apples-to-apples.\n\n"
    "**Why it is not trivial for the WP.** For the classical projectile (an "
    "external potential) *energy lost by the projectile* and *energy absorbed by "
    "the target* are exact duals by construction. For the WP the projectile is a "
    "genuine Kohn–Sham electron, so 'the projectile's energy' is a *choice of "
    "partition*, not a partition-free observable — the project convention forbids "
    "`−dKE_proj/ds` for the WP. The decomposed ledger lets us build *S* from named, "
    "attributable channels instead.\n\n"
    "> **This is Phase 1 (human-gated).** It (§A) inventories every recorded "
    "energy, (§B) verifies numerically how they compose, (§C0) measures the "
    "slab–projectile term `E_sp`, (§C) *proposes* the stopping-power formulae, and "
    "(§D) brainstorms setup simplifications. **Phases 2–6 do not run until the "
    "formulae in §C — and the E_sp treatment — are approved.**")]

cells += [R.md(
    "### Where this sits\n\n"
    "| Campaign | Role | Definition it owns |\n"
    "|---|---|---|\n"
    "| `classical-highdensity-sv` | **data generation** — records the ledger with "
    "clean-exit plateaus | Definition 2 = `E_absorbed/L_slab` (headline); reserves "
    "Definition 1 (data-collect only) |\n"
    "| **`stopping-from-energy-decomposition`** (this) | **analysis** — derives & "
    "applies the decomposition formula | **Definition 1** = S from the decomposed "
    "ledger |\n\n"
    "No GPU runs are launched here. Data = the plateaued+decomposed runs from "
    "`classical-highdensity-sv` (once its sweep completes) **plus** any other select "
    "runs identified with the user at the start of Phase 3.")]

# ------------------------------------------------------------------ setup
cells += [R.setup_cell()]
cells += [R.code(
    "import pandas as pd\n"
    "pd.set_option('display.width', 140)\n"
    "pd.set_option('display.max_columns', 40)\n"
    "HA_EV = 27.211386  # 1 Hartree in eV\n\n"
    "# Reference runs carrying the pairwise ledger (interactions.csv):\n"
    f"REF_CLASSICAL = {REF_CLASSICAL!r}\n"
    f"REF_WP        = {REF_WP!r}\n"
    f"REF_WP_LONG   = {REF_WP_LONG!r}\n"
    "print('reference runs set')")]

# ------------------------------------------------------------------ 2. Conventions
cells += [R.md(
    "## Conventions & symbols\n\n"
    "Native INQ energies are in **Hartree**; we report physics in **eV** "
    "(1 Ha = 27.211 eV). Charges are grouped into three sets: **P** = projectile "
    "(the Gaussian charge, or the WP orbital), **S** = slab electrons, **B** = "
    "positive jellium background.\n\n"
    "| Symbol | Meaning | Where recorded |\n"
    "|---|---|---|\n"
    "| `energy_total` | INQ Kohn–Sham total electronic energy | `observables.csv` |\n"
    "| `energy_kinetic` | electronic kinetic energy (orbitals; **+ WP** in a WP run) | `observables.csv` |\n"
    "| `energy_hartree` | electron–electron Hartree (Coulomb) energy | `observables.csv` |\n"
    "| `energy_xc` | exchange–correlation energy | `observables.csv` |\n"
    "| `energy_external` | electron ↔ external-potential energy (background; **+ projectile** for classical) | `observables.csv` |\n"
    "| `energy_nonlocal`, `energy_ion` | pseudopotential nonlocal / ion terms | `observables.csv` |\n"
    "| `E_SS` (`e_ss`) | slab–slab Coulomb | `interactions.csv` |\n"
    "| `E_PP` (`e_pp`) | projectile **self**-Hartree | `interactions.csv` |\n"
    "| `E_PS` (`e_ps`) | projectile–slab Coulomb | `interactions.csv` |\n"
    "| `E_SB` (`e_sb`) | slab–background Coulomb | `interactions.csv` |\n"
    "| `E_PB` (`e_pb`) | projectile–background Coulomb | `interactions.csv` |\n"
    "| `E_BB` (`e_bb`) | background self (constant) | `interactions.csv` |\n"
    "| `proj_z`, `proj_vz` | classical projectile position / velocity | `projectile.csv` |\n"
    "| `energy_proj_ke` | classical projectile ½·m·v² | `projectile.csv` |\n"
    "| `energy_proj_bg_ideal` | projectile↔background reference (`U_proj_bg`) | `projectile.csv` |\n\n"
    "**Gauge caveat (read before any WP pair-term is interpreted).** WP insertion "
    "makes the cell net −1 charged; individual `energy_hartree`/`energy_external` "
    "(and single pair terms `E_PP`/`E_PS`/`E_PB`) are then Poisson-G=0-convention "
    "dependent. Only combinations whose gauge-invariant terms Δ(E_SS,E_SB,E_BB)≈0 "
    "across the twin are physical (`reference_charged_cell_hartree_convention`).")]

# ------------------------------------------------------------------ 3. Setup / provenance
cells += [R.md(
    "## Reference runs used in Phase 1\n\n"
    "Phase 1 verifies the *bookkeeping* on runs that already carry the full pairwise "
    "ledger (independent of the final Phase-3 selection). The **`twin_ec_rsweep` "
    "r12 pair** is a matched classical+WP calibration pair (3 steps each) — enough "
    "to demonstrate both closure branches; the **`phase5_wp`** run (301 steps, with "
    "built-in `e_hartree_check`/`e_external_check`) shows closure holds over a full "
    "trajectory.\n\n"
    "| Run | kind | steps | ledger | note |\n"
    "|---|---|---|---|---|\n"
    "| `twin_ec_rsweep/r12_classical` | classical | 3 | E_SS,E_PP,E_PS,E_SB,E_PB,E_BB + projectile.csv | matched pair |\n"
    "| `twin_ec_rsweep/r12_wp` | WP | 3 | + e_hartree_check, e_external_check | matched pair |\n"
    "| `phase5_wp/p5_null_s2_k4_wp` | WP | 301 | + e_*_check | closure over full run |\n\n"
    "These are *demonstration* runs for the closure identities — **not** the "
    "stopping-power data. Stopping requires plateaued deposit runs, selected with "
    "the user at Phase 3.")]

# ------------------------------------------------------------------ 4. Source files
cells += [R.md(
    "## Source files\n\n"
    "| File | Role |\n"
    "|---|---|\n"
    "| `docs/campaigns/localised_jellium/stopping-from-energy-decomposition.md` | this campaign prompt |\n"
    "| `docs/campaigns/localised_jellium/classical-highdensity-sv-benchmark.md` | the data-generation campaign (Definition 2) |\n"
    "| `inq-stack/include/inqkit/jellium/interaction_energies.hpp` | emits the pairwise ledger; asserts the closure relations (l.17-18) |\n"
    "| `inq-stack/include/inqkit/io/observables_writer.hpp` | the `energy_*` column schema |\n"
    "| `.claude/skills/twin-run-analysis/twin_decompose.py` | reference decomposition engine (P/S/B ledger, residual, SIE) |\n"
    "| `.claude/skills/stopping-power-extraction/stopping_power.py` | Definition-2 deposit kernels (Correa 2018) |\n"
    "| `ResearchProject/.../hypotheses/stopping_from_decomposition/build_stopping_from_decomposition_report.py` | this builder |")]

# ------------------------------------------------------------------ 5A. Energy inventory
cells += [R.md(
    "## §A — Every energy recorded in a run\n\n"
    "A run writes energies to three CSVs under `raw/observables/`. The cell below "
    "prints the **actual** column headers of the reference runs (self-verifying); "
    "the table after it is the full inventory across run types.")]
cells += [R.code(
    "def cols(path):\n"
    "    return list(pd.read_csv(path, nrows=0).columns)\n"
    "for label, base in [('r12_classical', REF_CLASSICAL), ('r12_wp', REF_WP),\n"
    "                    ('phase5_wp (long)', REF_WP_LONG)]:\n"
    "    print(f'=== {label} ===')\n"
    "    for csv in ('observables.csv', 'interactions.csv', 'projectile.csv'):\n"
    "        p = os.path.join(base, csv)\n"
    "        if os.path.exists(p):\n"
    "            print(f'  {csv:18s}', cols(p))\n"
    "        else:\n"
    "            print(f'  {csv:18s} (absent)')\n"
    "    print()")]
cells += [R.md(
    "**Full inventory (grouped by file and by what is a member of `energy_total`).**\n\n"
    "*Lumped Kohn–Sham stores — `observables.csv`* (in `energy_total`):\n"
    "`energy_total`, `energy_kinetic`, `energy_hartree`, `energy_xc`, "
    "`energy_external`, `energy_nonlocal`, `energy_ion` "
    "(`energy_ion` is **classical-only** — the Ehrenfest ion term). Diagnostic "
    "extras a run *may* enable (NOT in the total): `energy_ion_kinetic`, "
    "`energy_exact_exchange`, `energy_nvxc`, `energy_eigenvalues`, "
    "`energy_proj_bg_ideal`/`_impl`.\n\n"
    "*Pairwise Coulomb ledger — `interactions.csv`* (a decomposition of the Hartree "
    "+ external Coulomb energy, **not** additional energy): "
    "`e_ss, e_pp, e_ps, e_sb, e_pb, e_bb`. WP runs add `e_hartree_check`, "
    "`e_external_check` (INQ's own Hartree/external recomputed, for closure), and "
    "`norm_wp`, `norm_total`; classical runs carry `norm_slab`, `norm_proj`.\n\n"
    "*Projectile track — `projectile.csv`* (**classical-only**): "
    "`proj_z, proj_vz, energy_proj_ke, energy_proj_bg_ideal`. For the WP there is "
    "no projectile track — its 'position' is the centroid of the WP orbital "
    "density, reconstructed in post-processing.\n\n"
    "*Density fields — VTIs* (physical order; never `fftshift`): feed the density "
    "GIF and any real-space wake-force channel; not an energy column.")]

# ------------------------------------------------------------------ 5B. Closure
cells += [R.md(
    "## §B — How the energies compose (closure, numerically verified)\n\n"
    "The pairwise ledger is a *decomposition* of the lumped Coulomb stores, so it "
    "must **close** exactly. The classical projectile is a ghost (z_valence=0) so "
    "its self-Hartree is absent from `energy_hartree`; the WP is a real electron so "
    "its self-Hartree `E_PP` and its coupling `E_PS` are inside `energy_hartree`.")]
cells += [R.md(
    "**Closure relations** (asserted in `interaction_energies.hpp:17-18`):\n\n"
    "*Classical:*\n"
    "$$E_\\mathrm{hartree}=E_{SS},\\qquad E_\\mathrm{external}=E_{SB}+E_{PS}.$$\n\n"
    "*Wavepacket:*\n"
    "$$E_\\mathrm{hartree}=E_{SS}+E_{PS}+E_{PP},\\qquad E_\\mathrm{external}=E_{SB}+E_{PB}.$$\n\n"
    "*Total (both):*\n"
    "$$E_\\mathrm{total}=E_\\mathrm{kinetic}+E_\\mathrm{hartree}+E_\\mathrm{xc}"
    "+E_\\mathrm{external}\\;(+\\,E_\\mathrm{nonlocal}+E_\\mathrm{ion}).$$\n\n"
    "The cell below loads each reference run, forms the residuals of these "
    "identities, and reports them at step 0 and their max over the run — they "
    "should vanish to floating-point precision (~1e-10 Ha).")]
cells += [R.code(
    "def load_merged(base):\n"
    "    obs  = pd.read_csv(os.path.join(base, 'observables.csv'))\n"
    "    iact = pd.read_csv(os.path.join(base, 'interactions.csv'))\n"
    "    return obs.merge(iact, on=['step', 'time_au'], how='inner')\n\n"
    "def etotal_rhs(df):\n"
    "    parts = ['energy_kinetic','energy_hartree','energy_xc','energy_external',\n"
    "             'energy_nonlocal','energy_ion']\n"
    "    return sum(df[c] for c in parts if c in df.columns)\n\n"
    "def closure(base, kind):\n"
    "    df = load_merged(base)\n"
    "    if kind == 'classical':\n"
    "        r_h = df['energy_hartree'] - df['e_ss']\n"
    "        r_e = df['energy_external'] - (df['e_sb'] + df['e_ps'])\n"
    "    else:  # wp\n"
    "        r_h = df['energy_hartree'] - (df['e_ss'] + df['e_ps'] + df['e_pp'])\n"
    "        r_e = df['energy_external'] - (df['e_sb'] + df['e_pb'])\n"
    "    r_t = df['energy_total'] - etotal_rhs(df)\n"
    "    row = {'run': os.path.basename(os.path.dirname(os.path.dirname(base))),\n"
    "           'kind': kind, 'n': len(df),\n"
    "           'r_hartree@0': r_h.iloc[0], 'max|r_hartree|': r_h.abs().max(),\n"
    "           'r_external@0': r_e.iloc[0], 'max|r_external|': r_e.abs().max(),\n"
    "           'r_total@0': r_t.iloc[0], 'max|r_total|': r_t.abs().max()}\n"
    "    return row\n\n"
    "rows = [closure(REF_CLASSICAL, 'classical'),\n"
    "        closure(REF_WP, 'wp'),\n"
    "        closure(REF_WP_LONG, 'wp')]\n"
    "tab = pd.DataFrame(rows).set_index('run')\n"
    "with pd.option_context('display.float_format', lambda v: f'{v:.2e}'):\n"
    "    display(tab)\n"
    "assert tab[['max|r_hartree|','max|r_external|','max|r_total|']].max().max() < 1e-8, \\\n"
    "    'closure FAILED — decomposition does not reconstruct the lumped stores'\n"
    "print('\\nAll closure residuals < 1e-8 Ha — the ledger reconstructs the lumped stores exactly.')")]
cells += [R.md(
    "**Reading it.** Every residual is ~1e-10 Ha (numerical noise). This is the "
    "foundation of Definition 1: because the pairwise ledger reconstructs "
    "`energy_hartree`, `energy_external`, and `energy_total` exactly, *any* stopping "
    "formula written in these channels is guaranteed consistent with the total "
    "energy — the decomposition adds attribution, not error.")]

# ------------------------------------------------------------------ 5C0. The E_sp question
cells += [R.md(
    "## §C0 — The E_sp question: does dropping the slab–projectile term cost us?\n\n"
    "C1 (below) builds the deposit from `T_slab + E_SS + E_SB + E_xc` and **excludes** "
    "the slab–projectile interaction `E_sp ≡ E_PS`. That exclusion is exact only if "
    "`E_sp(t_final) ≈ E_sp(0)` — then it cancels in `Δ`. The user's premise (and the "
    "suspected cause of the deposit-based S(v) coming out ~8× the Lindhard bulk value) "
    "is that `E_sp(0)` is *non-negligible* while `E_sp(t_final) ≈ 0` once the "
    "projectile is absorbed by the CAP. We measure it directly on the ledger-carrying "
    "twin, then reconcile with the dynamic absorbing runs.")]
cells += [R.code(
    "# E_sp = E_PS (projectile <-> slab-ELECTRON Coulomb). Individually huge, but the\n"
    "# slab is charge-NEUTRAL, so E_PS is nearly cancelled by E_PB (projectile<->background):\n"
    "# the *net* projectile-slab coupling that actually enters the energy balance is small.\n"
    "def eps_terms(base):\n"
    "    iact = pd.read_csv(os.path.join(base, 'interactions.csv'))\n"
    "    r0, rf = iact.iloc[0], iact.iloc[-1]\n"
    "    return dict(n=len(iact),\n"
    "                e_ps0=r0['e_ps']*HA_EV, e_pb0=r0['e_pb']*HA_EV,\n"
    "                net0=(r0['e_ps']+r0['e_pb'])*HA_EV, netf=(rf['e_ps']+rf['e_pb'])*HA_EV)\n"
    "for lbl, base in [('classical r12', REF_CLASSICAL), ('wp r12', REF_WP)]:\n"
    "    d = eps_terms(base)\n"
    "    print(f\"{lbl:13s} n={d['n']}  E_PS(0)={d['e_ps0']:8.1f}  E_PB(0)={d['e_pb0']:8.1f}\"\n"
    "          f\"  net P-slab(0)={d['net0']:6.1f}  net(t_f)={d['netf']:6.1f} eV\")\n"
    "print()\n"
    "print('Ledger twins above are 3-step STATIC gauge runs (no traversal): magnitude & '\n"
    "      'screening only, not dynamics.')\n"
    "print('Dynamic ABSORBING runs (qsp_phase3/4) carry NO pairwise ledger; their reconstructed')\n"
    "print('net cross term  E_total(0) - <T_WP> - E_GS  is:')\n"
    "print('  qsp p4 (54 eV): +4.6 eV    qsp p3 (100 eV): +3.9 eV    [E_sp(t_f)~0, projectile absorbed]')\n"
    "print('  deposit ~ 59 eV -> removing the net E_sp shifts S from ~2.4 to ~2.2 eV/Bohr '\n"
    "      '(Lindhard ~0.3).')")]
cells += [R.md(
    "**Finding.** The raw `E_PS` is large (~−140 eV) but the neutral slab screens it: "
    "the net `E_PS + E_PB` is only a few eV, and in the dynamic absorbing runs the "
    "reconstructed net is ~+4 eV at t=0 and ~0 at t_final. So the **E_sp correction is "
    "a few-eV shift on a ~59 eV deposit — it moves S from ~2.4 to ~2.2 eV/Bohr and "
    "cannot explain the ~8× overshoot vs Lindhard (~0.3 eV/Bohr).** The overshoot is "
    "dominated instead by the WP's own internal energy (zero-point KE = 81.6 eV at "
    "σ_WP=0.5, drift up to ~180 eV) and the *approximate* CAP split of "
    "deposited-vs-carried-off energy, with the runs not fully converged (5–10% WP norm "
    "remaining). That is the **Phase-4 investigation** (hypotheses b/c); the E_sp branch "
    "is quantitatively **refuted** as the cause here. Crucially, *no single run today "
    "carries dynamics **and** the pairwise ledger **and** full absorption* — which is "
    "exactly what §D targets.")]

# ------------------------------------------------------------------ 5C. Formula proposals
cells += [R.md(
    "## §C — Proposed stopping-power formulae (for approval)\n\n"
    "Below are the ranked candidate definitions of *S* from the decomposed ledger. "
    "Each is a **proposal** — none is computed on data until it is approved and "
    "locked (Phase 2). Notation: `Δ(·)` is the change since t=0; `s` is the "
    "projectile path length; `ds` its choice is stated per formula.\n\n"
    "**Path-length `ds` choices** (every reported *S* carries one): projectile "
    "arc-length ∫|v|dt (classical, correct for a decelerating projectile); slab "
    "thickness `L_slab` (aggregate deposit); WP centroid arc-length (only while the "
    "WP density is unimodal — flagged otherwise).")]
cells += [R.md(
    "### C1 · Headline — matched-estimator total deposit (both runs)\n\n"
    "The **same functional of the same columns** in classical and WP, so the two "
    "are directly comparable:\n"
    "$$S \\;=\\; \\frac{\\Delta E_\\mathrm{target}}{\\Delta s},\\qquad "
    "\\Delta E_\\mathrm{target}(t) = E_\\mathrm{electronic}(t) - E_\\mathrm{electronic}(0),$$\n"
    "with the **channel split** (where the deposited energy goes):\n"
    "$$\\Delta E_\\mathrm{target} = \\Delta \\mathrm{KE}_\\mathrm{slab} "
    "+ \\Delta E_{SS} + \\Delta E_{SB} + \\Delta E_\\mathrm{xc}.$$\n"
    "*Meaning:* all energy that ends up in the slab-electron subsystem. *Classical "
    "self-test:* this equals the Definition-2 deposit `E_absorbed`, so `S = "
    "E_absorbed/L_slab` — the built-in validation in Phase 2. *Gauge:* KE_slab "
    "clean; ΔE_SS, ΔE_SB clean within-run; ΔE_xc lumps in the WP's own xc (a known "
    "WP contamination, quantified separately). *Status:* **recommended headline** — "
    "convention-consistent, apples-to-apples, self-validating.")]
cells += [R.md(
    "### C2 · Classical conservation anchor (sanity, classical-only)\n\n"
    "$$S(v_0) = -\\frac{d\\,(\\tfrac12 m\\, v_z^2)}{ds}\\quad\\text{over the early }"
    "v_z\\ge 0.85\\,v_0\\text{ window.}$$\n"
    "*Meaning:* the friction on the point charge; textbook stopping. *Role:* the "
    "**conservation cross-check** for C1 — in the classical run, `S_C1 + dU_proj_bg/ds` "
    "must equal this to numerical closure. **Never** the WP headline "
    "(`feedback_quantum_stopping_not_from_projectile_ke`).")]
cells += [R.md(
    "### C3 · WP projectile-partition, vacuum-corrected (exploratory, gauge-checked)\n\n"
    "The deliberate test of whether the decomposition can give the WP a projectile "
    "energy after all:\n"
    "$$E_\\mathrm{proj}^{WP} = \\mathrm{KE}_\\mathrm{proj} + E_{PP} + E_{PS} + E_{PB},"
    "\\qquad S = -\\frac{d\\,[\\,E_\\mathrm{proj}^{WP}(t) - E_\\mathrm{proj}^{vac}(t)\\,]}{ds}.$$\n"
    "*Why the vacuum baseline:* `E_PP ∝ 1/σ`, so a freely **spreading** WP would "
    "register as fake stopping; subtracting the same quantity from a WP-in-vacuum "
    "twin removes that drift, isolating slab-induced loss. *Gauge:* individual "
    "`E_PP/E_PS/E_PB` are convention-dependent — **only reported if the gauge test "
    "Δ(E_SS,E_SB,E_BB)≈0 passes**. *Status:* exploratory; its disagreement with C1 "
    "*measures* the quantum channel (energy into WP internal excitation vs slab "
    "deposit). Requires one extra cheap vacuum run.")]
cells += [R.md(
    "### C4 · Target-absorption & irreversibility qualifiers (diagnostics)\n\n"
    "- **A1 target-absorption** — C1 restricted to the *irreversible* late-time "
    "deposit `E_slab(t→∞)` after the projectile exits / stops (the reversible "
    "polarisation has quiesced).\n"
    "- **D1 matched-face** — `S̄ = [KE(z=+z_f) − KE(z=−z_f)]/L_slab` at mirror "
    "positions outside both faces, where reversible terms cancel (transmission only).\n"
    "- **D2 hysteresis loop** — `∮ F dz` in a channel (e.g. E_PS vs z): the enclosed "
    "area is the non-adiabatic dissipated work.\n\n"
    "These bound the *reversible-vs-dissipated* fraction and accompany the headline, "
    "never replace it. The full A–E brainstorm menu is in the appendix.")]

# ------------------------------------------------------------------ 5D. Setup-simplification brainstorm
cells += [R.md(
    "## §D — Setup modifications that simplify the energy bookkeeping (brainstorm)\n\n"
    "*Design principle:* the deposit we want is O(1–10 eV), but the WP arrives carrying "
    "~260 eV of internal energy (82 eV localisation + up to ~180 eV drift). Forming the "
    "deposit as a difference of two ~260 eV numbers inherits every systematic at full "
    "strength. Good setups therefore either make the big terms **cancel by "
    "construction**, make them **analytically known**, or **never touch them** by "
    "measuring the slab side directly. Ranked by payoff-to-cost:\n\n"
    "1. **Run-to-extinction + analytic launch ledger** *(completes the user's seed).* "
    "Extend each WP run (checkpoint-resume, `LJ_RESUME=1`) to residual WP norm <0.1% "
    "*and* launch from far vacuum so `E_sp(0)<0.1 eV`. Then `E_total(t_f)=E_slab(t_f)` "
    "exactly, and the t=0 reference is closed-form (`KE_zp=3/(4 m σ²)=82 eV`, "
    "self-Hartree≈22 eV, `KE_drift=½mv²`) — no CAP ledger at either endpoint. *Cost:* "
    "extra wall-time (cheap via resume); risk of slow stragglers / slab wake reaching "
    "the CAP.\n"
    "2. **Vacuum-twin subtraction.** Rerun the identical WP through an empty box; "
    "`deposit = (slab run) − (vacuum run)`. The vacuum run's true deposit is 0, so its "
    "reported 'deposit' *is* the method's own systematic — subtract it. Cancels all WP "
    "self-energy / dispersion / imperfect-absorption errors to first order. *Cost:* one "
    "cheap extra run; cancellation only first-order (WP enters the CAP distorted in the "
    "slab run).\n"
    "3. **Slab-side ledger in the absorbing runs** *(structural fix).* Port the pairwise "
    "machinery into the long CAP runs; `deposit = Δ[KE_s + E_SS + E_SB + E_BB + "
    "E_xc,s]`. The 260 eV of WP energy never enters, and S becomes a slope over the "
    "traversal, not one endpoint difference. *Cost:* implementation; xc non-additivity "
    "~0.1–1 eV during traversal (quantify vs the closed-boundary twins).\n"
    "4. **Work-integral cross-check.** Save the slab current `j_s`; compute "
    "`W = ∫dt ∫dr j_s·E_P` independently — no WP self-energy, no CAP, no E_GS. A "
    "genuinely different second number for the same deposit. *Cost:* frame storage; WP "
    "identifiability caveat on `E_P`.\n"
    "5. **Heavy-mass WP ladder.** m=10/100/1836: `KE_zp∝1/m` (82→8.2→0.82 eV), WP exits "
    "compact, CAP split becomes exact; must converge to the validated classical-ghost "
    "deposit. *Cost:* full run pairs; electron-mass endpoint still hardest; resolve the "
    "known m>1 artifact returns first.\n"
    "6. **σ_WP widening.** 0.5→1→2 Bohr: `KE_zp∝1/σ²` (82→20→5 eV). If extracted S "
    "tracks the leaked `KE_zp`, the leak is caught red-handed. *Cost:* bigger boxes; S "
    "physically depends on σ_WP (changes the observable — complements, not replaces).\n"
    "7. **Exact absorber energy accounting.** In the propagation callback, evaluate the "
    "energy of the removed piece `(1−M)ψ` each step; then "
    "`E_total(t)+E_removed(t)=E_total(0)` to machine precision. *Cost:* inqkit work; "
    "removed-piece cross-terms (kinetic + far-field estimate likely enough since removal "
    "is in vacuum).\n"
    "8. **Frozen-vKS counterfactual.** Propagate the WP through the *frozen* GS "
    "potential (no slab response), same CAP; subtract. Cancels everything elastic "
    "(dispersion, reflection, transmission) — the difference is the response=deposited "
    "energy. Sharper than the vacuum twin (trajectory matched *through* the slab). "
    "*Cost:* needs a frozen-Hamiltonian propagation mode.\n"
    "9. **No-CAP ballistic-exit + flux screen.** Remove the CAP, lengthen the box so the "
    "transmitted WP cannot return; record outgoing WP energy flux at a plane screen. "
    "`E_total` is then exactly conserved. *Cost:* large box; slow components never clear "
    "the slab (worst exactly where S is largest, low v) — a benchmark, not the sweep.\n\n"
    "**Recommendation.** Do **1 + 2** together on the next run pair: extend to WP "
    "extinction (resume-only cost) with a far-vacuum launch, and add the cheap vacuum "
    "twin as the systematic-error meter — no new physics code, both endpoints exact or "
    "analytic. Fold in **4** (work integral) as the first cross-check once frames are "
    "saved; treat **3** as the medium-term structural fix. *(Feeds Phase 6.)*")]

# ------------------------------------------------------------------ 5E. Critical answers
cells += [R.md(
    "## §E — Critical answers to the open questions\n\n"
    "The gate is passed (user directed proceeding **without re-runs**, 2026-07-22). "
    "The existing runs already answer the central question — *why did the WP "
    "deposit-based S(v) come out ~8× the Lindhard bulk value?* — decisively. The cell "
    "below puts the WP deposit next to the two things it could be made of: the "
    "projectile's **drift** kinetic energy (½mv² = the E column) and the WP's **fixed** "
    "zero-point energy (81.6 eV, set by σ_WP=0.5), plus the Lindhard-expected deposit "
    "for scale.")]
cells += [R.code(
    "# The decisive table: WP 'deposit' vs the WP's OWN energy, and vs Lindhard-expected.\n"
    "QSP = '/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses'\n"
    "ZP, L_SLAB = 81.634, 25.0   # WP zero-point KE (eV, FIXED by sigma_WP=0.5); slab thickness\n"
    "sw = pd.read_csv(os.path.join(QSP, 'qsp_phase5', 'se_state.csv'))\n"
    "d = sw[['tag','E_eV','v','deposited_eV','S_eVbohr','norm_f']].copy()\n"
    "d['dep/driftKE'] = d['deposited_eV'] / d['E_eV']      # E_eV IS 1/2 m v^2 (m=1)\n"
    "d['dep/zeroPt']  = d['deposited_eV'] / ZP\n"
    "d = d.round({'E_eV':0,'S_eVbohr':2,'deposited_eV':1,'dep/driftKE':2,'dep/zeroPt':2,'norm_f':3})\n"
    "display(d.set_index('tag'))\n"
    "for E, SL in [(54, 0.448), (100, 0.282)]:\n"
    "    print(f'Lindhard-expected deposit at {E} eV = {SL} x {L_SLAB:.0f} Bohr = {SL*L_SLAB:4.1f} eV'\n"
    "          f'  (measured ~59 eV => {59/(SL*L_SLAB):.0f}x too big)')\n"
    "print()\n"
    "print('Classical projectile, SAME geometry, 54 eV: S = 0.25 eV/Bohr (converged, '\n"
    "      'N-conserved, CAP drained 0.24%) = ~0.5x Lindhard (0.45).')\n"
    "print('=> the ~8x overshoot is WP-METHOD-specific, NOT the geometry.')")]
cells += [R.md(
    "### E.1 · Why is the WP deposit-based S ~8× Lindhard?\n\n"
    "**Decisive fact:** at v=1.3 and v=2.0 the deposit (59 eV) *exceeds the projectile's "
    "entire drift KE* (23 and 54 eV). No process that drains the projectile's **motion** "
    "can deposit more than the motion holds — the surplus must be the WP's own internal "
    "energy. Verdicts:\n\n"
    "- **(a) Localised ≠ bulk (wrong reference) — REFUTED as the cause** *(high conf).* "
    "The classical projectile in the *same* geometry gives S=0.25 eV/Bohr ≈ 0.5× "
    "Lindhard — an O(1) finite-size effect, which cannot make 8×. *Caveat:* the "
    "classical projectile parks inside the slab (path < L_slab), so 0.25 is itself a "
    "mild under-estimate; the control certifies geometry only at the factor-~2 level — "
    "which is all that is needed.\n"
    "- **(b) Wrong baseline / bookkeeping — SUPPORTED** *(high).* At low v the deposit "
    "is pinned ~59 eV ≈ 0.73× the *fixed* 81.6 eV zero-point, independent of projectile "
    "energy; at high v it grows toward the drift KE (472 of 490 eV at v=6, i.e. 96%). "
    "The ledger retains WP-internal energy — zero-point at low v, drift at high v. A "
    "Lindhard-order deposit would be only 7–11 eV over the slab; every row is 5–50× "
    "above.\n"
    "- **(c) CAP distorts the energy — SUPPORTED** *(high on invalidation; medium on the "
    "split).* The CAP exports *norm* efficiently but *energy* poorly (p3 removed 126 eV "
    "against ~173 eV WP-carried). This is structural to `E_total(t_f)−E_GS` under a "
    "non-Hermitian absorber, so running to norm_f→0 would **not** converge the deposit "
    "to true slab heating — the 'upper bound' flag is loose by ~an order of magnitude.\n\n"
    "**One mechanism.** (b) and (c) are the same failure: `deposit = E_total(t_f) − "
    "E_GS` charges the slab with *everything the WP brought into the box* — its fixed "
    "81.6 eV localisation energy, its self-interaction, its drift KE — minus only what "
    "the CAP exports, and the CAP exports norm far better than energy. The genuine slab "
    "heating (~7–11 eV) is an order of magnitude below the smallest measured deposit, so "
    "**no row of the S(v) table is a stopping power**, even as an upper bound.")]
cells += [R.md(
    "### E.2 · The definitional questions\n\n"
    "- **Why C1 excludes E_sp — justified, but re-scope the '<10%'.** E_sp is a *shared* "
    "interaction (partition choice, not slab energy) and is screened to a few-eV net "
    "(§C0). But 'few eV / 59 eV < 10%' is measured against the *contaminated* deposit; "
    "against the physically expected ~7–11 eV slab deposit a few-eV E_sp is a **30–50% "
    "effect**. Keep the exclusion, but **re-audit it once the bookkeeping is fixed** and "
    "deposits shrink to the eV scale.\n"
    "- **WP separability caveat — sound.** Kinetic and Hartree split exactly (per-orbital "
    "/ pairwise); xc is non-additive → ~eV ambiguity in T_slab/E_SS during traversal. "
    "Fine for *diagnosing* a 59 eV contaminant; the ledger is **diagnostic-grade, not "
    "measurement-grade** until the contaminant is removed.\n"
    "- **C1 vs C2 (classical) — agree only where the CAP flux is negligible.** By energy "
    "conservation `−dKE_proj = dE_electronic + dU_proj_bg`, so C1=C2 for the classical "
    "run *because* its CAP drain is 0.24%. But C2 (KE-slope) **fails for a fully-stopping "
    "light projectile** (KE→0, face-window gave 0.0) — the light-projectile rule — so "
    "C1/energy-method is the robust one. Caveat: the parked projectile makes L_slab the "
    "wrong denominator, and the gap between 'lost 54 eV of KE' and 'deposited ~6 eV' must "
    "be closed in the ledger before 0.25 is enshrined.\n"
    "- **C3 (WP projectile-partition) — a cheap cross-check, not a dead end.** Raw C3 "
    "inherits the zero-point/dispersion contamination, but a **vacuum-referenced** C3 "
    "cancels the 81.6 eV zero-point and free dispersion *by construction* — which raw C1 "
    "does not. The cancellation is imperfect (the slab alters dispersion via "
    "capture/scatter), so exploratory; but 'not expected to beat C1' is too strong.")]
cells += [R.md(
    "### E.3 · Hypotheses we had NOT listed\n\n"
    "1. **The t=0 baseline is already contaminated.** `E_total(0) − E_GS` contains the "
    "WP's zero-point, self-Hartree and xc self-interaction *before propagation*. The fix "
    "is a better *observable*, not a better CAP: measure the **bath-only** energy rise "
    "directly from the pairwise decomposition (slab heating, never charging the WP's own "
    "energy to the slab).\n"
    "2. **xc self-interaction of the extra KS orbital** — σ-dependent, sits in the "
    "deposit, and does **not** leave with the CAP even in principle. Itemise it.\n"
    "3. **CAP reflection of slow components** — at low v the dispersed WP has slow tails a "
    "fast-tuned CAP partly reflects; re-entrant amplitude re-deposits. Distinguish from "
    "zero-point leakage by a CAP-width/strength sweep at fixed v.\n"
    "4. **Is the 1/σ² localisation energy real physics?** Devil's advocate: an electron "
    "entering a metal does exchange internal energy. But stopping power is the drag on "
    "*drift* motion, and the 81.6 eV term is fixed by the representation choice σ_WP=0.5 "
    "(∝1/σ²) — a 'measurement' that scales with a free numerical parameter is not "
    "physics. **Decisive test: a σ_WP sweep at fixed v** — if the low-v deposit floor "
    "tracks 1/σ², the artifact is proven.\n"
    "5. **The −47 eV 'bath/collective' CAP entry** looks unphysical (a negative removal) "
    "and may be a plasmon double-count between the WP-carried and bath columns; audit the "
    "ledger's sign structure before using it to apportion (b) vs (c).")]
cells += [R.md(
    "### E.4 · Bottom line\n\n"
    "The **only trustworthy number** on the table is the **classical S = 0.25 eV/Bohr at "
    "54 eV** (converged, N-conserving, CAP-untouched) — good at the factor-~2 level given "
    "its path-length caveat. It certifies that the localised slab reproduces "
    "Lindhard-order stopping. Every WP entry is an **artifact ledger**: floor-pinned at "
    "~0.7× the 81.6 eV zero-point at low v, retaining ~96% of the drift KE at high v. "
    "**Implication:** the next step is not more WP runs but a **rebuilt observable** — a "
    "bath-only (pairwise-decomposed) or vacuum-referenced energy deposit — validated by a "
    "σ_WP sweep of the low-v floor, before any WP stopping power is quoted. That is "
    "exactly what §D fixes #1/#2/#3 target.")]

# ------------------------------------------------------------------ 6. Takeaway + gate
cells += [R.md(
    "## Takeaway (Phase 1) & the approval gate\n\n"
    "- **§A** — a run records the lumped KS stores (`observables.csv`), a pairwise "
    "Coulomb ledger (`interactions.csv`: E_SS/E_PP/E_PS/E_SB/E_PB/E_BB), and — "
    "classically — a projectile track (`projectile.csv`).\n"
    "- **§B** — the ledger reconstructs the lumped stores to ~1e-10 Ha for both "
    "representations, so any channel-based *S* is automatically total-energy "
    "consistent.\n"
    "- **§C0** — `E_sp` measured: raw `E_PS`≈−140 eV but screened by `E_PB` to a "
    "few-eV *net*; correcting for it moves S ~2.4→~2.2 eV/Bohr and **cannot** explain "
    "the ~8× Lindhard overshoot. The overshoot is the WP internal energy + CAP split "
    "(Phase 4); the E_sp branch is refuted as the cause.\n"
    "- **§C** — proposed: **C1** matched-estimator total deposit (headline, both "
    "runs) with its channel split; **C2** classical KE anchor (sanity); **C3** "
    "vacuum-corrected WP projectile-partition (exploratory, gauge-gated); **C4** "
    "irreversibility qualifiers.\n"
    "- **§D** — 9 ranked setup simplifications; top pick = run-to-extinction + "
    "vacuum-twin subtraction (feeds Phase 6).\n"
    "- **§E** — critical answers on the existing data: the ~8× overshoot is "
    "WP-method-specific (the deposit *exceeds* the drift KE at low v), = WP-internal "
    "energy leaking through an energy-lossy CAP; (a) refuted, (b)+(c) are one "
    "mechanism; the classical S=0.25 eV/Bohr is the only trustworthy number.\n\n"
    "> ### ✅ GATE PASSED (2026-07-22)\n"
    "> The user directed proceeding to the critical-analysis phase **without re-runs**. "
    "§E answers every open question on the existing runs. **Next (needs data / user):** "
    "the fix is a *rebuilt observable* (bath-only decomposed deposit or vacuum-"
    "referenced), not more runs — Phase 2 locks the C1 kernel with the E_sp term "
    "switchable; Phases 3–6 need the ledger-carrying plateau runs (§D fix #1/#3), which "
    "await the deferred re-run decision.")]

R.build(cells, OUT)
