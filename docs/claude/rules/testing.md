# Rule: Testing and Validation

Apply to: `inq/src/`, `ResearchProject/`, `Tutorial/`, `docs/plans/`, `docs/handovers/`

## Rules

1. Every substantive code change must include a validation update before it is considered complete.

0. Use GPU execution whenever a GPU is available. For INQ, this means running `inq-run` (not `inq-run --cpu`). For Python scripts, prefer CUDA-enabled backends where relevant. Only fall back to CPU if GPU is unavailable or explicitly requested by the user.

2. Before implementing, define:
   - **Component tests**: one or more unit tests per core function (inputs, edge cases, expected outputs, units).
   - **Integration tests**: tests that exercise the assembled workflow end-to-end.
   - **Scientific validation tests**: tests against known analytic solutions, conservation laws, benchmark datasets, or published reference values.

3. Propose the full validation menu to the user before running expensive simulations, long test suites, or destructive workflows. The user decides which ones to authorise.

4. After implementation, every handover must state:
   - which tests were proposed
   - which were approved by the user
   - which were actually run
   - observed outcomes (pass/fail)
   - what remains unverified

5. Never claim a calculation is correct without evidence from at least one completed validation test.

## INQ-specific validation checklist

When adding or modifying INQ calculations:
- [ ] Energy conservation during real-time propagation (check drift rate)
- [ ] Ground state SCF convergence to tolerance
- [ ] Forces sum to zero (Newton's 3rd law) for isolated systems
- [ ] Known system benchmark (H2 total energy, N2 dipole response, etc.)
- [ ] GPU/CPU consistency check (run same input with and without CUDA, compare outputs)
- [ ] Restart check (load saved ground state, re-run few steps, compare to original)
