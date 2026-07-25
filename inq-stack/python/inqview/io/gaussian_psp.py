"""Generate an erf-smoothed Gaussian-charge electron pseudopotential (UPF).

A classical −1 (electron) projectile of Gaussian charge width ``sigma`` produces,
on the *bath* electrons, the repulsive local potential

    V(r) = C * erf(r / (sigma * sqrt(2))) / r          (C / r as r -> inf)

i.e. the potential of a Gaussian charge distribution. In reciprocal space this is
the bare Coulomb (4*pi/q^2) times the form factor exp(-q^2 sigma^2 / 2).

This module takes an existing **bare local Coulomb electron** UPF as a template
(z_valence = 0, number_of_proj = 0, sign-flipped so the tail is *positive* =
repulsive = electron) and replaces only its ``PP_LOCAL`` block with the
erf-smoothed form for a chosen sigma. Projectors and core charge are already
absent in the template and are left untouched.

The Coulomb coefficient ``C`` (and therefore the unit convention, Rydberg vs
Hartree) is read *from the template's own large-r tail*, so the generated file
matches the template's units exactly — no hard-coded factor of 2.

Reference for the erf-smoothed (Gaussian-charge) Coulomb form: standard result,
e.g. soft/regularised Coulomb pseudopotentials (the Gaussian charge -> erf
potential pair); see `ResearchProject/literature/.../soft-scf-pseudopotentials`.

SIGMA CONVENTION (UNIFIED 2026-06-21 — wavepacket is the single source of truth).
``generate_gaussian_psp`` takes ``sigma_wp``, the *wavepacket* sigma shared with
``inqkit::WavePacket`` (psi ~ exp(-r^2 / 2 sigma_wp^2), so the charge/density std is
sigma_charge = sigma_wp / sqrt(2)). The erf potential is built from sigma_charge, so a
CLASSICAL projectile and a WAVEPACKET given the SAME sigma present the IDENTICAL
charge cloud exp(-r^2 / sigma_wp^2) to the bath. The low-level ``v_erf_hartree`` still
takes the *charge std* (the honest math); only this high-level entry point converts.

WARNING: pre-2026-06-21 UPFs on disk (``electron_gaussian_sigma*.upf``) were generated
with the OLD convention where the filename sigma was the CHARGE STD directly
(= unified sigma_wp / sqrt(2)). See CONTEXT.md "sigma-convention unification" for the
full registry mapping every legacy file/run to its unified sigma_wp.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = ["GaussianPspResult", "v_erf_hartree", "generate_gaussian_psp"]


def v_erf_hartree(r: np.ndarray, sigma: float) -> np.ndarray:
    """erf-smoothed repulsive electron potential in **Hartree** (C=1).

    V(r) = erf(r / (sigma*sqrt2)) / r, with the finite r->0 limit
    V(0) = sqrt(2/pi) / sigma.

    Here ``sigma`` is the **charge-density std** (the honest math of the erf form).
    In the unified convention this equals sigma_wp / sqrt(2); callers wanting the
    wavepacket sigma should use ``generate_gaussian_psp(sigma_wp=...)``.
    """
    r = np.asarray(r, dtype=float)
    out = np.empty_like(r)
    small = r < 1e-12
    out[~small] = (
        np.vectorize(math.erf)(r[~small] / (sigma * math.sqrt(2.0))) / r[~small]
    )
    out[small] = math.sqrt(2.0 / math.pi) / sigma
    return out


@dataclass(frozen=True)
class GaussianPspResult:
    path: Path
    sigma_wp: float            # unified (wavepacket) sigma: psi ~ exp(-r^2/2 sigma_wp^2)
    sigma_charge: float        # charge-density std = sigma_wp/sqrt(2) (what the erf uses)
    coulomb_coeff: float       # C in V = C*erf(.)/r, in the template's units
    v0_template_units: float   # V(0) in template units
    v0_hartree: float          # V(0) in Hartree = sqrt(2/pi)/sigma_charge
    n_mesh: int


_FLOAT = re.compile(r"[-+]?\d+\.\d+(?:[EeDd][-+]?\d+)?")


def _read_block(text: str, tag: str) -> tuple[int, int, list[float]]:
    """Return (start_idx, end_idx, values) for a <tag ...> ... </tag> block.

    start_idx/end_idx are character offsets of the numeric payload (between the
    opening tag's closing '>' and the closing '</tag>').
    """
    open_m = re.search(rf"<{tag}\b[^>]*>", text)
    close_m = re.search(rf"</{tag}>", text)
    if open_m is None or close_m is None:
        raise ValueError(f"UPF template missing <{tag}> block")
    payload_start = open_m.end()
    payload_end = close_m.start()
    values = [float(x) for x in _FLOAT.findall(text[payload_start:payload_end])]
    return payload_start, payload_end, values


def _format_block(values: np.ndarray, columns: int = 4) -> str:
    """Format values as fixed-width columns matching QE/ONCV UPF style."""
    lines = []
    for i in range(0, len(values), columns):
        row = values[i : i + columns]
        lines.append("".join(f"{v:21.10E}" for v in row))
    return "\n" + "\n".join(lines) + "\n"


def generate_gaussian_psp(
    template_path: str | Path,
    sigma_wp: float,
    out_path: str | Path,
    *,
    tail_fit_range: tuple[float, float] = (8.0, 16.0),
) -> GaussianPspResult:
    """Write an erf-smoothed UPF for the given **wavepacket** sigma.

    ``sigma_wp`` is the UNIFIED (wavepacket) sigma (see module docstring). The
    charge-density std actually used in the erf is ``sigma_charge = sigma_wp/sqrt(2)``,
    so the classical projectile presents the same cloud exp(-r^2/sigma_wp^2) as a
    wavepacket of the same sigma. The Coulomb coefficient C and the unit convention
    are inferred from the template's large-r tail (V*r -> C). The generated PP_LOCAL
    is ``C * erf(r/(sigma_charge*sqrt2)) / r`` on the template's own radial mesh.
    """
    sigma_charge = sigma_wp / math.sqrt(2.0)
    template_path = Path(template_path)
    out_path = Path(out_path)
    text = template_path.read_text()

    _, _, r_vals = _read_block(text, "PP_R")
    loc_start, loc_end, loc_vals = _read_block(text, "PP_LOCAL")
    r = np.asarray(r_vals, dtype=float)
    v_template = np.asarray(loc_vals, dtype=float)
    if r.size != v_template.size:
        raise ValueError(
            f"mesh mismatch: PP_R has {r.size}, PP_LOCAL has {v_template.size}"
        )

    # Infer C (and units) from the template tail: V_template(r) ~ C/r.
    lo, hi = tail_fit_range
    mask = (r >= lo) & (r <= hi)
    if not mask.any():
        raise ValueError("no template mesh points in tail_fit_range")
    coulomb_coeff = float(np.mean(v_template[mask] * r[mask]))

    # New PP_LOCAL: C * erf(r/(sigma_charge*sqrt2)) / r, finite at r=0.
    v_new = coulomb_coeff * v_erf_hartree(r, sigma_charge)  # charge-std form

    new_block = _format_block(v_new, columns=4)
    new_text = text[:loc_start] + new_block + text[loc_end:]
    out_path.write_text(new_text)

    # Hartree V(0): template-units V(0)/ (C / C_hartree); C_hartree=1 form gives
    # V(0)_Ha = sqrt(2/pi)/sigma_charge regardless of template units.
    v0_template = coulomb_coeff * math.sqrt(2.0 / math.pi) / sigma_charge
    v0_hartree = math.sqrt(2.0 / math.pi) / sigma_charge
    return GaussianPspResult(
        path=out_path,
        sigma_wp=sigma_wp,
        sigma_charge=sigma_charge,
        coulomb_coeff=coulomb_coeff,
        v0_template_units=v0_template,
        v0_hartree=v0_hartree,
        n_mesh=r.size,
    )
