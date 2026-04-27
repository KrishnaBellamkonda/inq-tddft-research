"""Phase: ``summary`` — verify and (re)build ``results/run_summary.txt``.

The C++ run already writes a complete ``run_summary.txt`` (see
``coronene/shared/cpp/run_template.hpp``). This phase only:

  1. Confirms the file exists.
  2. Appends a ``post-processing`` block listing what artefacts the
     postprocess produced (run_name, host, timestamp, sizes).

Future expansion: validate that the file has every section required by the
spec (§3.2).
"""

from __future__ import annotations

import os
import platform
import socket
import time
from pathlib import Path

from . import pipeline as _pipeline


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    summary = results_dir / "run_summary.txt"
    if not summary.exists():
        _pipeline.skip(
            f"run_summary.txt missing at {summary} (the C++ run did not "
            "complete a final write; postprocess will continue but downstream "
            "tools may need to look at the stub)."
        )

    block = (
        "\n10. Post-processing\n"
        "-------------------\n"
        f"postprocess_host       = {socket.gethostname()}\n"
        f"postprocess_user       = {os.environ.get('USER', '?')}\n"
        f"postprocess_platform   = {platform.platform()}\n"
        f"postprocess_timestamp  = {time.strftime('%Y-%m-%dT%H:%M:%S')}\n"
        f"postprocess_run_name   = {run_name}\n"
    )

    text = summary.read_text()
    if "10. Post-processing" in text and not rebuild:
        return {"summary_path": str(summary), "appended": False}

    if "10. Post-processing" in text and rebuild:
        head, _, _ = text.partition("\n10. Post-processing")
        text = head.rstrip() + "\n"
    summary.write_text(text + block)
    return {"summary_path": str(summary), "appended": True}
