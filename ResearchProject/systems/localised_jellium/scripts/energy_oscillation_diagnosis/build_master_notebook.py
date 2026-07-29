#!/usr/bin/env python3
"""build_master_notebook.py — assemble the master study notebook for the
localised-jellium ΔE_total energy-oscillation diagnosis from per-probe result.json.

House narrative (per notebook-making): intro (phenomenon + ledger mirror) → one
section PER experiment (Aim → Method → What was plotted → Results → advisor
verdict) → final synthesis (confirmed mechanism). Rebuilt after each iteration so
the notebook is always current while the loop runs unattended.

Reads:  hypotheses/energy_oscillation_diagnosis/probes/*/result.json
        hypotheses/energy_oscillation_diagnosis/hypothesis_ledger.md
Writes: hypotheses/energy_oscillation_diagnosis/energy_oscillation_diagnosis.ipynb
        (executed via `jupyter nbconvert --execute` by the caller / analyse tail)

Per-probe advisor verdicts are read from result.json["advisor_verdict"] if the
agent has written it back (the loop appends the Advisor's JSON to each probe).
"""
from __future__ import annotations
import json
from pathlib import Path

import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HYP = Path(__file__).resolve().parents[2] / "hypotheses" / "energy_oscillation_diagnosis"
PROBES = HYP / "probes"
LEDGER = HYP / "hypothesis_ledger.md"
OUT = HYP / "energy_oscillation_diagnosis.ipynb"


def _load_probes() -> list[dict]:
    out = []
    for rj in sorted(PROBES.glob("*/result.json")):
        try:
            out.append(json.loads(rj.read_text()))
        except Exception as e:  # keep building even if one probe is mid-write
            out.append({"name": rj.parent.name, "_error": str(e)})
    return out


def _probe_section(p: dict) -> list:
    name = p.get("name", "?")
    s = p.get("summary", {})
    v = p.get("advisor_verdict") or {}
    md = [f"## Experiment: `{name}`", ""]
    md += [f"**Aim.** {p.get('aim','—')}", ""]
    md += [f"**Method.** {p.get('method','—')}", ""]
    if p.get("component_gap"):
        md += [f"> Component gap (requested but not recorded): `{p['component_gap']}`", ""]
    md += ["**What was plotted.** ΔE_total(t) in both reference conventions "
           "(vs E_GS and vs E_total(0) of the RT run) over the top panel, and the "
           "per-component drift Δ(E_kin/E_H/E_xc/E_ext/E_ion) over the bottom.", ""]
    md += ["**Results (raw).**", "",
           f"- ΔE vs RT-t0: final `{s.get('dE_total_vs_rt0_final')}` eV, "
           f"max `{s.get('dE_total_vs_rt0_max')}` eV",
           f"- ΔE vs E_GS: final `{s.get('dE_total_vs_gs_final')}` eV, "
           f"max `{s.get('dE_total_vs_gs_max')}` eV",
           f"- crosses 0 from below (unphysical rise): `{s.get('crosses_zero_above')}`",
           f"- N: `{s.get('N_initial')}` → `{s.get('N_final')}`",
           f"- Σε_i tracks E_total: `{s.get('eig_tracks_total')}`", ""]
    if v:
        md += ["**Advisor verdict.**", "",
               f"- reading: {v.get('probe_reading','—')}",
               f"- verdict: `{v.get('verdict','—')}`  confidence `{v.get('confidence','—')}`",
               f"- next: {(v.get('next_probe') or {}).get('name','—')}", ""]
    cells = [new_markdown_cell("\n".join(md))]
    png = p.get("plot_png")
    if png and Path(png).exists():
        cells.append(new_code_cell(
            "from IPython.display import Image\n"
            f"Image(filename={json.dumps(png)})"))
    return cells


def build() -> Path:
    probes = _load_probes()
    ledger = LEDGER.read_text() if LEDGER.exists() else "_(ledger not yet seeded)_"

    nb = new_notebook()
    nb.cells.append(new_markdown_cell(
        "# Localised-jellium ΔE_total energy-oscillation diagnosis\n\n"
        "Autonomous agent+advisor investigation. **Goal: diagnose + document only** "
        "— isolate the cause of the unphysical ΔE_total>0 rise to one confirmed "
        "mechanism; no physics fix.\n\n"
        "Phenomenon: `ΔE_total(t)=E_total(t)−E_ref` oscillates and rises above 0 once "
        "the CAP absorbs — unphysical (no energy influx; a CAP can only remove "
        "energy). Source note: "
        "`docs/notes/localised-jellium-energy-oscillation-investigation.md`."))
    nb.cells.append(new_markdown_cell("## Hypothesis ledger (live mirror)\n\n" + ledger))

    if not probes:
        nb.cells.append(new_markdown_cell("_No probes run yet._"))
    for p in probes:
        if p.get("_error"):
            nb.cells.append(new_markdown_cell(
                f"## Experiment: `{p['name']}`\n\n_result.json unreadable: {p['_error']}_"))
            continue
        nb.cells.extend(_probe_section(p))

    # Synthesis: pull the "Confirmed cause" section verbatim from the ledger if present.
    synth = ("_(advisor's final verdict — the confirmed mechanism and the decisive "
             "control that proves it — is written here at loop end.)_")
    if LEDGER.exists():
        txt = LEDGER.read_text()
        marker = "## Confirmed cause"
        if marker in txt:
            tail = txt.split(marker, 1)[1].strip()
            # stop at the next top-level heading if any
            for stop in ("\n## ",):
                if stop in tail:
                    tail = tail.split(stop, 1)[0].strip()
            if tail and "advisor fills in" not in tail:
                synth = tail
    nb.cells.append(new_markdown_cell("## Synthesis (confirmed cause)\n\n" + synth))

    OUT.write_text(nbf.writes(nb))
    print(f"wrote {OUT} ({len(nb.cells)} cells, {len(probes)} probes)")
    return OUT


if __name__ == "__main__":
    build()
