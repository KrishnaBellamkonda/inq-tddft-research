#!/usr/bin/env python3
"""Email a CAP-in-jellium baseline run summary on completion (autonomous).

Usage: email_on_done.py <subdir> <label> <exit_status>
Reads results/<subdir>/run_summary.txt and sends it to the user, threaded under
the [cap-jellium-baseline] family.
"""
import sys
from pathlib import Path

sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.email import send_run_email  # noqa: E402

HERE = Path(__file__).resolve().parent
subdir, label, status = sys.argv[1], sys.argv[2], sys.argv[3]

summary_path = HERE / "results" / subdir / "run_summary.txt"
summary = summary_path.read_text() if summary_path.exists() else "(run_summary.txt missing)"

body = f"""Baseline 1 ({label}) finished. Process exit status: {status}.

This is an AUTONOMOUS update from the CAP-in-jellium baselines task.

--- run_summary.txt -------------------------------------------------
{summary}
---------------------------------------------------------------------

Written to results/{subdir}/:
  - raw/vti/density_system/   : ~300 density frames -> region-N (free vs slab),
                                 n(z,t), E-field (post), drainage profile.
  - raw/observables/electron_number.csv : N(t) every step (drainage curve).
  - raw/observables/state_energies.csv  : per-orbital <psi|H|psi>(t) + variance.
  - raw/observables/occupations_vs_time.csv, observables.csv, eigenvalues/.

NOTE: ALL absorption numbers PROVISIONAL until the inq-study engine regression
(Task #7) passes. This is Baseline 1 (CAP on, NO projectile) = the bath-drainage
reference. B2 (classical) and B3 (WP) are the next runs.
"""
mid = send_run_email(
    f"[cap-jellium-baseline] B1 {label} done",
    body,
    to="chiddukanna@gmail.com",
)
print("emailed:", mid)
