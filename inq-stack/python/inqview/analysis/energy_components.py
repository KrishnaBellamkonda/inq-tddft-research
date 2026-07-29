"""Functional energy-component flow (IV-M07).

Tracks where energy goes during a propagation by decomposing the total energy
into its KS functional components — kinetic, Hartree, exchange-correlation, and
the external/electron-ion term recovered as the residual. This sidesteps the
band-sum double-counting problem (the old energy_balance ledger) because it uses
the actual energy functional INQ writes to ``observables.csv``.

    E_ext(t) = E_total(t) - [ E_kin(t) + E_H(t) + E_xc(t) ]   (residual)

so by construction the four components sum to E_total at every step — that is
the kernel's exact test invariant. The companion renderers (initial-vs-final
bars, ΔE(t) lines, breakdown GIF) live in ``inqview.visualisation``.

Pure: numpy + pandas only (no matplotlib / VTK).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd

PathOrFrame = Union[str, Path, pd.DataFrame]

HA_TO_EV = 27.211386245988

_REQUIRED = ("time_au", "energy_total", "energy_kinetic", "energy_hartree", "energy_xc")


@dataclass(frozen=True)
class EnergyComponents:
    """Per-step energy decomposition (Hartree units unless noted).

    ``E_ext`` is the external/electron-ion residual ``E_total − (kin+H+xc)``.
    The ``d*`` properties are changes from the first recorded step.
    """

    time_au: np.ndarray
    E_total: np.ndarray
    E_kin: np.ndarray
    E_hartree: np.ndarray
    E_xc: np.ndarray
    E_ext: np.ndarray

    # --- changes from t0 (the "energy flow") ---
    @property
    def dE_total(self) -> np.ndarray:
        return self.E_total - self.E_total[0]

    @property
    def dE_kin(self) -> np.ndarray:
        return self.E_kin - self.E_kin[0]

    @property
    def dE_hartree(self) -> np.ndarray:
        return self.E_hartree - self.E_hartree[0]

    @property
    def dE_xc(self) -> np.ndarray:
        return self.E_xc - self.E_xc[0]

    @property
    def dE_ext(self) -> np.ndarray:
        return self.E_ext - self.E_ext[0]

    def component_sum(self) -> np.ndarray:
        """E_kin + E_H + E_xc + E_ext (== E_total by construction)."""
        return self.E_kin + self.E_hartree + self.E_xc + self.E_ext

    def breakdown(self, when: str = "final") -> dict[str, float]:
        """Component values (Ha) at the first/last step — for the bar chart."""
        i = 0 if when == "initial" else -1
        return {"kinetic": float(self.E_kin[i]), "hartree": float(self.E_hartree[i]),
                "xc": float(self.E_xc[i]), "external": float(self.E_ext[i]),
                "total": float(self.E_total[i])}

    def redistribution_ev(self) -> dict[str, float]:
        """Net ΔE per component from t0 to the final step, in eV (where it went)."""
        return {"kinetic": float(self.dE_kin[-1] * HA_TO_EV),
                "hartree": float(self.dE_hartree[-1] * HA_TO_EV),
                "xc": float(self.dE_xc[-1] * HA_TO_EV),
                "external": float(self.dE_ext[-1] * HA_TO_EV),
                "total": float(self.dE_total[-1] * HA_TO_EV)}


def compute(observables: PathOrFrame) -> EnergyComponents:
    """Decompose ``observables.csv`` (or a DataFrame) into energy components."""
    df = (observables if isinstance(observables, pd.DataFrame)
          else pd.read_csv(observables, comment="#"))
    missing = [c for c in _REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"observables missing columns {missing}; have {list(df.columns)}")

    t = df["time_au"].to_numpy(dtype=float)
    E_total = df["energy_total"].to_numpy(dtype=float)
    E_kin = df["energy_kinetic"].to_numpy(dtype=float)
    E_H = df["energy_hartree"].to_numpy(dtype=float)
    E_xc = df["energy_xc"].to_numpy(dtype=float)
    E_ext = E_total - (E_kin + E_H + E_xc)          # electron-ion residual
    return EnergyComponents(t, E_total, E_kin, E_H, E_xc, E_ext)
