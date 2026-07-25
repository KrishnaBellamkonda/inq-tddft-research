#!/usr/bin/env python3
"""Generic autonomous run-completion emailer.

Usage: email_run.py <subdir> <label> <exit_status> <family> [extra_note]
Reads results/<subdir>/run_summary.txt and emails it, threaded under [<family>].
"""
import sys
from pathlib import Path

sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.email import send_run_email  # noqa: E402

HERE = Path(__file__).resolve().parent
subdir, label, status, family = sys.argv[1:5]
extra = sys.argv[5] if len(sys.argv) > 5 else ""

sp = HERE / "results" / subdir / "run_summary.txt"
summary = sp.read_text() if sp.exists() else "(run_summary.txt missing)"
body = f"""{label} finished. Process exit status: {status}.

AUTONOMOUS update.{(' ' + extra) if extra else ''}

--- run_summary.txt -------------------------------------------------
{summary}
---------------------------------------------------------------------
Output: results/{subdir}/raw/  (observables.csv, momentum_distribution.csv [WP],
  electron_track.csv [classical], state_energies, occupations, density VTI).

ALL absorption numbers PROVISIONAL until the inq-study engine regression (Task #7).
"""
mid = send_run_email(f"[{family}] {label} done", body, to="chiddukanna@gmail.com")
print("emailed:", mid)
