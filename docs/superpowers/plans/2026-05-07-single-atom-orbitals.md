# Single-atom orbitals — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:executing-plans (inline) — chosen by user. Steps use checkbox `- [ ]` syntax for tracking.

**Goal:** Generate ground-state Kohn–Sham orbitals (occupied + ≥30 extras) for an isolated H, Li, and Al atom in a 30 bohr cubic vacuum box, written as VTI series for ParaView visualisation.

**Architecture:** Three self-contained INQ tutorial programs under `Tutorial/single-atom-orbitals/{h,li,al}/run.cpp`. Each runs a stand-alone LDA SCF at 60 Ry, then loops over every Kohn–Sham state and writes both `|ψ_i|²` (real) and `ψ_i` (complex) via `inqkit::io::RealField3DWriter` / `ComplexField3DWriter`.

**Tech stack:** INQ (header-only C++17, GPU via CUDA), inqkit field/IO library, `inq-run` build wrapper.

---

## File structure

| Path | Responsibility | Created/Modified |
|---|---|---|
| `Tutorial/single-atom-orbitals/h/run.cpp`  | H atom GS + per-state writes  | Create |
| `Tutorial/single-atom-orbitals/li/run.cpp` | Li atom GS + per-state writes | Create |
| `Tutorial/single-atom-orbitals/al/run.cpp` | Al atom GS + per-state writes | Create |
| `docs/handovers/single-atom-orbitals.md`   | Rolling handover               | Create |

No shared header is introduced. The three `run.cpp` files are intentionally near-duplicates (differ only in atom symbol and a one-line comment about valence). The DRY violation is acceptable because each is a *tutorial*: a learner should be able to read one file end-to-end without chasing imports.

---

### Task 1: Create `h/run.cpp`

**Files:**
- Create: `Tutorial/single-atom-orbitals/h/run.cpp`

- [ ] **Step 1: Write the run program**

```cpp
// ============================================================================
// single-atom-orbitals/h: ground state of an isolated H atom in vacuum
//
// System: 1 H atom centred in a 30 bohr cubic finite cell.
//         LDA, 60 Ry cutoff, gamma-only, extra_states(30).
//         Restricted KS with Fermi smearing (T = 0.001 Ha) to handle the
//         single valence electron's fractional occupancy gracefully.
//
// Outputs (under results/):
//   density/                  total electron density (one VTI)
//   orbital_density/          |psi_i|^2 for every KS state i
//   orbitals/                 complex psi_i for every KS state i
//
// Purpose: visualise the s/p/d orbital ladder of an isolated atom in
//          ParaView. High-lying states (i >~ 15) are box modes, not
//          true Rydberg states — see spec for caveat.
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>

#include <cstdio>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;

int main() {
    auto const L = 30.0_bohr;
    auto cell = systems::cell::cubic(L).finite();

    systems::ions ions(cell);
    auto centre = L / 2.0;
    ions.insert("H", {centre, centre, centre});

    std::cout << "\n=== single-atom-orbitals/h ===\n";
    std::cout << "Cell: cubic finite, L = 30 bohr\n";

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .cutoff(60.0_Ry)
            .extra_states(30)
            .temperature(0.001_Ha),
        input::kpoints::gamma());

    ground_state::initial_guess(ions, electrons);

    auto gs = ground_state::calculate(
        ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(1e-6_Ha)
            .max_steps(1000)
            .broyden_mixing());

    std::cout << "\nSCF done.\n";
    std::cout << "  GS total energy = " << gs.energy.total() << " Ha\n";

    // ── Sanity check: total density integrates to N_valence ────────────────
    auto rho_total = inqkit::fields::density::total(electrons);
    auto n_int = operations::integral(rho_total);
    std::cout << "  integral of rho = " << n_int
              << " (expected ~ 1.0 valence electron)\n";

    // ── Eigenvalue ladder ──────────────────────────────────────────────────
    auto const & phi = electrons.kpin()[0];
    int const n_states = phi.spinor_set_size();
    std::cout << "  n_states = " << n_states << "\n";

    std::cout << "\nEigenvalue ladder (Ha):\n";
    for (int i = 0; i < n_states; ++i) {
        auto const eig = electrons.eigenvalues()[0][i];
        auto const occ = electrons.occupations()[0][i];
        std::printf("    [%3d]  eps = %+10.5f Ha   occ = %7.4f\n",
                    i, eig, occ);
    }

    // ── Write total density ────────────────────────────────────────────────
    std::cout << "\nWriting total density...\n";
    inqkit::io::RealField3DWriter(
        "results/density",
        {.field_name = "total_density", .include_meta = true},
        {.overwrite = true})
        .write(rho_total, "density_total");

    // ── Write every orbital (density and complex wavefunction) ─────────────
    std::cout << "Writing " << n_states << " orbital densities...\n";
    inqkit::io::RealField3DWriter rho_writer(
        "results/orbital_density",
        {.field_name = "orbital_density", .include_meta = true},
        {.overwrite = true});

    inqkit::io::ComplexField3DWriter psi_writer(
        "results/orbitals",
        {.field_name = "wavefunction", .include_meta = true},
        {.overwrite = true});

    for (int i = 0; i < n_states; ++i) {
        char tag[64];
        std::snprintf(tag, sizeof(tag), "orbital_%04d_density", i);
        rho_writer.write(inqkit::fields::density::orbital(electrons, i), tag);

        std::snprintf(tag, sizeof(tag), "orbital_%04d", i);
        psi_writer.write(inqkit::fields::orbital::wavefunction(electrons, i), tag);
    }

    std::cout << "\nDone. Output in results/\n";
    return 0;
}
```

- [ ] **Step 2: Build and run**

Run: `cd Tutorial/single-atom-orbitals/h && inq-run`
Expected:
- Compile succeeds (header-only INQ + inqkit).
- SCF converges within `max_steps = 1000`; final line `GS total energy = -0.4… Ha` (LDA H total energy is approximately `-0.46 Ha`; pseudopotential values may differ).
- `integral of rho` prints close to `1.0`.
- `n_states = 31` printed.
- 31 entries in the eigenvalue ladder, lowest with `occ ≈ 1.0`.
- `results/density/`, `results/orbital_density/`, `results/orbitals/` populated.

- [ ] **Step 3: Spot-check outputs**

Run: `ls Tutorial/single-atom-orbitals/h/results/orbital_density/ | wc -l`
Expected: ≥ 31 (one per orbital) plus a manifest file.

Run: `ls Tutorial/single-atom-orbitals/h/results/orbitals/ | wc -l`
Expected: ≥ 31 plus manifest.

If either count is wrong, do not proceed — debug the loop.

---

### Task 2: Create `li/run.cpp`

**Files:**
- Create: `Tutorial/single-atom-orbitals/li/run.cpp`

- [ ] **Step 1: Write the run program (identical to Task 1 except atom symbol and header comment)**

Differences from `h/run.cpp`:
- Header line `single-atom-orbitals/li: ground state of an isolated Li atom in vacuum`.
- Comment block notes `1 valence electron (1s² in core)`.
- `ions.insert("Li", {centre, centre, centre});`
- Sanity-check expected integral `~ 1.0` (Li with pseudopotential carries 1 valence electron).
- Reference LDA Li atom total energy with pseudopotential is approximately `-0.20 Ha` — vary per pseudopotential.

The rest of the file is byte-for-byte identical to `h/run.cpp`.

- [ ] **Step 2: Build and run**

Run: `cd Tutorial/single-atom-orbitals/li && inq-run`
Expected:
- SCF converges.
- `integral of rho ≈ 1.0`.
- `n_states = 31`.
- Output folders populated.

- [ ] **Step 3: Spot-check outputs**

Run: `ls Tutorial/single-atom-orbitals/li/results/orbital_density/ | wc -l`
Expected: ≥ 31 + manifest.

---

### Task 3: Create `al/run.cpp`

**Files:**
- Create: `Tutorial/single-atom-orbitals/al/run.cpp`

- [ ] **Step 1: Write the run program**

Differences from `h/run.cpp`:
- Header line `single-atom-orbitals/al: ground state of an isolated Al atom in vacuum`.
- Comment block notes `3 valence electrons (Ne core in pseudopotential)`.
- `ions.insert("Al", {centre, centre, centre});`
- Sanity-check expected integral `~ 3.0`.
- `n_states = 32` (ceil(3/2)=2 + 30) is the expected count.
- Reference LDA Al atom total energy with pseudopotential is approximately `-2.0 Ha` — vary per pseudopotential.

- [ ] **Step 2: Build and run**

Run: `cd Tutorial/single-atom-orbitals/al && inq-run`
Expected:
- SCF converges (may take more iterations than H/Li due to 3p partial occupation; smearing helps).
- `integral of rho ≈ 3.0`.
- `n_states = 32`.

- [ ] **Step 3: Spot-check outputs**

Run: `ls Tutorial/single-atom-orbitals/al/results/orbital_density/ | wc -l`
Expected: ≥ 32 + manifest.

---

### Task 4: Commit production code

- [ ] **Step 1: Stage and commit**

```bash
git add Tutorial/single-atom-orbitals/ \
        docs/superpowers/specs/2026-05-07-single-atom-orbitals-design.md \
        docs/superpowers/plans/2026-05-07-single-atom-orbitals.md
git commit -m "$(cat <<'EOF'
add single-atom-orbitals tutorial (H, Li, Al ground states)

Three isolated-atom GS runs in 30 bohr cubic finite cells, LDA at
60 Ry with extra_states(30). Each writes the total density plus the
complex wavefunction and |psi|^2 for every Kohn-Sham state, intended
for ParaView orbital visualisation.
EOF
)"
```

Note: the commit-message rule forbids the words `claude`, `anthropic`, `ai` (case-insensitive); the message above complies.

---

### Task 5: Write rolling handover

**Files:**
- Create: `docs/handovers/single-atom-orbitals.md`

- [ ] **Step 1: Write the handover with the standard sections**

Sections required (per `.claude/rules/handovers.md`):
- `## Current status` — what is done; what is verified.
- `## What changed` — list of new files.
- `## Files touched` — absolute paths of created files.
- `## Commands run` — `inq-run` invocations and observed outputs (energies, n_states).
- `## Tests and validation` — SCF tolerance, ∫ρ check, ladder printed.
- `## Trusted sources used` — INQ tutorial file `Tutorial/hf-gs-with-inqkit/run.cpp` (template); spec at `docs/superpowers/specs/2026-05-07-single-atom-orbitals-design.md`.
- `## Attribution notes` — adapted from existing `hf-gs-with-inqkit/run.cpp` pattern.
- `## Known issues / blockers` — any pseudopotential warnings; any SCF instability.
- `## Assumptions still in play` — high-i orbitals are box modes, not true Rydberg.
- `## Exact next steps` — open `results/orbital_density/*.vti` in ParaView; visual inspection.

- [ ] **Step 2: Commit**

```bash
git add docs/handovers/single-atom-orbitals.md
git commit -m "add handover for single-atom-orbitals tutorial"
```
