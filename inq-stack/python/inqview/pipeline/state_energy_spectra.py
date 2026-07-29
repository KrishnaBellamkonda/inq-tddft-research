"""Phase: ``state_energy_spectra`` — per-state KS energy εN(t) FFT and the
opposite-phase pair diagnostic.

This phase is the workhorse of the **electron-hole vs plasmon** discriminator
described in the BCN:1719P quantum-kick paper (and Santervas-Arranz, Stengel,
Artacho, *Phys. Rev. Research* 7, 033292, 2025).

Idea
----
A single-particle electron-hole transition n -> n' at energy ω satisfies

    εN(t)  ≈ ε_n_0  + A * cos(ω t + φ)
    εN'(t) ≈ ε_n'_0 - A * cos(ω t + φ)        (opposite phase)

so the cross-spectrum |FFT(εN) · conj(FFT(εN'))| has a sharp peak at ω with
phase ≈ π (180°). A plasmon (collective mode) has no such single (n, n')
pair signature.

Inputs
------
- ``raw/observables/state_energies.csv`` (from
  ``inqkit::observables::StateEnergyWriter``); columns:
      step, time_au, kpoint_index, state_index, weight, occupation,
      E_expect_ha, E_variance_ha2

Outputs (under ``analysis/observables/state_energy_spectra/``)
--------------------------------------------------------------
- ``per_state_peak.csv``        — per-(kpoint, state) dominant FFT peak
                                  (omega_ev, amplitude) on plateau-detrend
                                  variant of εN(t).
- ``per_state_peak_scatter.png`` — scatter (state_index, ω_peak) coloured by
                                  occupation; eye-balled clustering shows
                                  which energies many states share.
- ``cross_spectrum_pairs.csv``  — top-N candidate (n, n') pairs ranked by
                                  cross-spectrum amplitude with phase ≈ π.
- ``cross_spectrum_pairs.png``  — bar chart of top pairs at ω ≈ FFT_peak_ev
                                  alongside the gamma_transitions histogram
                                  (when available).

Scope
-----
Single-rank only (kpoint loop is sequential). The per-state FFT cost is
O(n_states · N log N) — for the Li 54-atom run with n_states = 808 and
N ≈ 1550 (state_energies sampled every 10 propagation steps) this is
~10⁷ FLOPs, sub-second wall.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import _common

_HA_TO_EV = 27.21138625
_AU2FS = 0.02418884


def _plateau_detrend(s: np.ndarray, plateau_frac: float = 0.5) -> np.ndarray:
    """Match observables._build_variants['plateau_detrend']: subtract the
    mean of the last (1 - plateau_frac) fraction. Default frac=0.5 ⇒
    second-half mean (QBall recipe)."""
    N = s.size
    p0 = int(N * plateau_frac)
    return s - (s[p0:].mean() if p0 < N else s.mean())


def _hann_fft(s: np.ndarray, dt_au: float, pad: int = 8):
    """Plateau-detrend ⇒ Hann ⇒ pad ⇒ rfft. Returns (energy_ev, amplitude,
    phase_rad). Amplitude normalised by N (not N_padded)."""
    N = s.size
    if N < 8:
        return None
    win = np.hanning(N)
    sig = np.zeros(N * pad, dtype=np.float64)
    sig[:N] = s * win
    spec = np.fft.rfft(sig)
    freq_au = np.fft.rfftfreq(N * pad, d=dt_au)
    omega_au = 2 * np.pi * freq_au
    energy_ev = _HA_TO_EV * omega_au
    amp = np.abs(spec) / N
    phase = np.angle(spec)
    return energy_ev, amp, phase, spec / N


def _peak_in_band(energy_ev, amplitude, lo=0.5, hi=20.0):
    """Return (omega_peak_ev, amplitude_peak) with the peak masked to
    [lo, hi] eV (kills the DC tail and the high-ω noise)."""
    mask = (energy_ev > lo) & (energy_ev < hi)
    if not mask.any():
        return float("nan"), 0.0
    idx = np.argmax(amplitude[mask])
    return float(energy_ev[mask][idx]), float(amplitude[mask][idx])


def run(results_dir: Path, *, run_name: str, rebuild: bool, **opts) -> dict:
    csv = results_dir / "raw" / "observables" / "state_energies.csv"
    if not csv.exists():
        return {"skipped": f"missing {csv} (run uses old template?)"}

    df = pd.read_csv(csv)
    if df.empty:
        return {"skipped": "state_energies.csv is empty"}

    # Time step from the data (decoupled from propagation dt because the
    # writer is invoked every STATE_ENERGY_EVERY propagation steps).
    times = np.sort(df["time_au"].unique())
    if times.size < 8:
        return {"skipped": f"too few samples ({times.size}) for FFT"}
    dt_au = float(np.median(np.diff(times)))

    out_dir = _common.ensure_dir(
        results_dir / "analysis" / "observables" / "state_energy_spectra")

    # ── per-(k, state) FFT peak ───────────────────────────────────────────
    rows = []
    spectra = {}    # (k, state) -> (energy_ev, amp, phase, spec_normalised)
    for (k, st), grp in df.groupby(["kpoint_index", "state_index"]):
        grp = grp.sort_values("time_au")
        e_t = grp["E_expect_ha"].to_numpy() * _HA_TO_EV
        if e_t.size < 8:
            continue
        # Use only the rows aligned to the unique time grid.
        e_t_detrend = _plateau_detrend(e_t)
        result = _hann_fft(e_t_detrend, dt_au=dt_au, pad=8)
        if result is None:
            continue
        energy_ev, amp, phase, spec = result
        peak_ev, peak_amp = _peak_in_band(energy_ev, amp)
        occ = float(grp["occupation"].iloc[0])
        weight = float(grp["weight"].iloc[0])
        rows.append({
            "kpoint_index": int(k),
            "state_index":  int(st),
            "occupation":   occ,
            "weight":       weight,
            "epsilon_t0_ev": float(e_t[0]),
            "peak_omega_ev":  peak_ev,
            "peak_amplitude": peak_amp,
        })
        spectra[(int(k), int(st))] = (energy_ev, amp, phase, spec)

    if not rows:
        return {"skipped": "no per-state spectra computable"}

    per_state = pd.DataFrame(rows).sort_values(
        ["kpoint_index", "state_index"]).reset_index(drop=True)
    per_state_csv = out_dir / "per_state_peak.csv"
    if _common.need_rebuild(per_state_csv, rebuild):
        per_state.to_csv(per_state_csv, index=False)

    # ── per-state peak scatter ────────────────────────────────────────────
    scatter_png = out_dir / "per_state_peak_scatter.png"
    if _common.need_rebuild(scatter_png, rebuild):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
        sc = ax.scatter(
            per_state["state_index"],
            per_state["peak_omega_ev"],
            c=per_state["occupation"],
            cmap="viridis",
            s=12 * (per_state["peak_amplitude"]
                    / max(per_state["peak_amplitude"].max(), 1e-30)) + 1.0,
            alpha=0.7)
        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label("occupation")
        ax.set_xlabel("state index")
        ax.set_ylabel(r"per-state $\hbar\omega_{\rm peak}$ (eV)")
        ax.set_title(f"{run_name}: per-state KS-energy FFT peak\n"
                     "(point size ∝ peak amplitude; "
                     "clusters ⇒ shared excitation)")
        ax.set_ylim(0, 20)
        fig.tight_layout()
        fig.savefig(scatter_png)
        plt.close(fig)

    # ── opposite-phase pair diagnostic (e-h fingerprint) ──────────────────
    # Loop over each kpoint independently (vertical transitions stay
    # within one k by construction). Thresholds are normalised against
    # the per-kpoint full-occupation value (which in INQ equals the
    # kpoint weight times spin pairing — typically 0.25 for an 8-kpoint
    # spin-paired metal). occupied = occ > 50% of full; unoccupied =
    # occ < 10% of full.
    pairs_rows = []
    occ_full = float(per_state["occupation"].max())
    occ_thresh_high = 0.50 * occ_full
    occ_thresh_low  = 0.10 * occ_full
    for k_target in sorted(per_state["kpoint_index"].unique()):
        g = per_state[per_state["kpoint_index"] == int(k_target)]
        occ_states = g[g["occupation"] > occ_thresh_high]["state_index"].tolist()
        unocc_states = g[g["occupation"] < occ_thresh_low]["state_index"].tolist()
        # Prevent O(N²) blowup at k=0 by capping search to
        # |ε_n' - ε_n| < 10 eV (covers the paper's regimes 0–7 eV).
        eps0_by_state = dict(zip(g["state_index"], g["epsilon_t0_ev"]))
        for n in occ_states:
            spec_n = spectra.get((int(k_target), n))
            if spec_n is None:
                continue
            energy_ev_n, _, _, S_n = spec_n
            for np_ in unocc_states:
                de_ev = eps0_by_state[np_] - eps0_by_state[n]
                if de_ev <= 0.5 or de_ev > 10.0:
                    continue
                spec_np = spectra.get((int(k_target), np_))
                if spec_np is None:
                    continue
                _, _, _, S_np = spec_np
                # Cross-spectrum amplitude and phase at the candidate energy
                # bin. We sample at the FFT bin closest to de_ev.
                idx = int(np.argmin(np.abs(energy_ev_n - de_ev)))
                cs = S_n[idx] * np.conj(S_np[idx])
                amp = float(np.abs(cs))
                phase_deg = float(np.degrees(np.angle(cs)))
                # Opposite-phase metric: cos((180° - |phase_deg|) deg);
                # = 1.0 for perfect anti-phase, ≤ 0 for in-phase.
                opp_phase = float(np.cos(np.radians(180.0 - abs(phase_deg))))
                pairs_rows.append({
                    "kpoint_index":   int(k_target),
                    "n_occ":          int(n),
                    "n_unocc":        int(np_),
                    "delta_epsilon_ev": float(de_ev),
                    "cross_amp":      amp,
                    "phase_deg":      phase_deg,
                    "opposite_phase_metric": opp_phase,
                })

    pairs_csv = out_dir / "cross_spectrum_pairs.csv"
    if pairs_rows:
        pairs_df = pd.DataFrame(pairs_rows).sort_values(
            "cross_amp", ascending=False).head(50)
        if _common.need_rebuild(pairs_csv, rebuild):
            pairs_df.to_csv(pairs_csv, index=False)

        # Visualise top pairs
        pairs_png = out_dir / "cross_spectrum_pairs.png"
        if _common.need_rebuild(pairs_png, rebuild):
            import matplotlib.pyplot as plt
            top = pairs_df.head(15)
            fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
            colours = ["#cc4444" if m > 0.7 else "#888888"
                       for m in top["opposite_phase_metric"]]
            ax.barh(range(len(top)), top["cross_amp"], color=colours,
                    alpha=0.85)
            labels = [f"({n},{nu})  Δε={de:.2f} eV  φ={ph:+.0f}°"
                      for n, nu, de, ph in zip(
                          top["n_occ"], top["n_unocc"],
                          top["delta_epsilon_ev"], top["phase_deg"])]
            ax.set_yticks(range(len(top)))
            ax.set_yticklabels(labels, fontsize="x-small")
            ax.set_xlabel("cross-spectrum amplitude  "
                          "|FFT(εN) · conj(FFT(εN′))|")
            ax.set_title(f"{run_name}: e-h pair candidates "
                          "(red ⇒ opposite-phase: e-h fingerprint)")
            ax.invert_yaxis()
            fig.tight_layout()
            fig.savefig(pairs_png)
            plt.close(fig)

    return {
        "out_dir": str(out_dir),
        "n_states_with_spectrum": len(per_state),
        "n_pairs_evaluated": len(pairs_rows),
        "dt_au": dt_au,
    }
