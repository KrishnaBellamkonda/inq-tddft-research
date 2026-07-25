# Rule: VTI coordinate mapping (index ↔ physical)

Apply to: any code that reads/plots INQ field data — `inq-stack/python/inqview/`,
`ResearchProject/systems/**/hypotheses/`, `**/make_*_postproc.py`, run
`analyse.py`, and notebook builders (`build_*_report.py`). Always on.

## The one rule

**inqkit VTIs are written in PHYSICAL order.** `inqkit::io::RealField3DWriter`
applies `fft_shift_index()` at write time and stamps `Origin = −L/2`, so array
index 0 already maps to the left-edge coordinate `−L/2`.

- **Never `np.fft.fftshift` a VTI** (or any field loaded from a `.vti`). Doing so
  swaps centre↔edge and silently produces flipped pictures — the recurring
  "slab/cluster appears at the edges, vacuum in the middle" bug. It looks
  plausible, so it is NOT caught by eye.
- **Load every VTI through `inqview.load_vti`** (canonical loader,
  `visualisation/field_io.py`). It returns data in physical order **plus** the
  `(x, y, z)` cell-centred coordinate axes and hard-asserts the axis/dim
  invariants. Use the returned axes for `extent=`; do not reconstruct coordinates
  by hand. Pass `expect_centered_axis="z"` for runs whose feature (slab/cluster)
  sits at the box centre — it fails loudly on a centre↔edge swap.
- **The ONLY data that needs `np.fft.fftshift` is LEED screen `.dat`** files —
  those are genuinely FFT-natural (peak lands at a corner otherwise). Load via
  `inqview.io.load_leed_pattern`, which already shifts.

## Why this keeps biting

The fact lived buried at line ~790 of the 825-line `tddft-simulations` skill, far
from where a model writing a loader re-derives the (wrong) assumption "INQ is
FFT-natural ⇒ fftshift". Hoisted here as an always-on rule so it triggers wherever
field-reading code is being written. See also `CONTEXT.md` → "Density
decomposition" and ADR on the canonical loader.
