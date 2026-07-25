"""fourier_panel.py — the project FFT-pipeline diagnostic panel (figure standard).

For EVERY FFT of a TDDFT time-series, emit a 3×2 panel that shows the processing
stages, so a spectrum is never a black box:

    row 1:  raw signal              |  baseline-removed signal (canonical 'mean')
    row 2:  windowed signal         |  zero-padded signal
    row 3:  |FFT| linear scale      |  |FFT| log scale

The canonical baseline is **mean** removal (user verdict 2026-06-25) and the FFT
panels overlay the **detrend** result as a dashed comparison line, so the choice
is always visible.

The stages are computed from `inqview.analysis.fourier.FourierTransform` itself
(its own `_apply_subtract`, `window.build`, `zero_pad`, coherent-gain norm), so the
panel is guaranteed to depict exactly what the kernel does — no re-derivation. The
final FFT panels reuse `FourierTransform.transform()` verbatim; a regression test
asserts the panel's stages and `transform()` stay consistent.

Use via the run-notebook / notebook-making skills whenever a signal is FFT'd.
"""
from __future__ import annotations

import numpy as np

from ..analysis.fourier import FourierTransform
from . import style as ST

HA_TO_EV = 27.211386245988


def fft_stages(time_au, values, ft: FourierTransform | None = None):
    """Return the intermediate arrays of `FourierTransform`'s pipeline as a dict.

    Mirrors `FourierTransform.transform` step-for-step (transient skip → subtract →
    window → zero-pad), reusing the kernel's own `_apply_subtract` and
    `window.build`, then calls `transform()` for the authoritative FFT. Keys:
    t, raw, detrended, window, windowed, t_pad, padded, freq_au, amplitude.
    """
    ft = ft or FourierTransform()
    t = np.asarray(time_au, float)
    v = np.asarray(values, float)

    # transient skip (identical rule to FourierTransform.transform §13.6)
    if ft.t_start_au > 0.0:
        mask = t >= (float(t[0]) + ft.t_start_au)
        t, v = t[mask], v[mask]

    n = len(t)
    dt = float(t[1] - t[0])
    detr = ft._apply_subtract(v)            # kernel's own baseline removal
    win = ft.window.build(n)                # kernel's own window
    windowed = detr * win
    if ft.zero_pad > 1:
        padded = np.zeros(n * ft.zero_pad, dtype=float)
        padded[:n] = windowed
    else:
        padded = windowed
    t_pad = t[0] + dt * np.arange(len(padded))

    res = ft.transform(time_au, values)     # authoritative FFT (re-applies skip)
    return dict(t=t, raw=v, detrended=detr, window=win, windowed=windowed,
                t_pad=t_pad, padded=padded,
                freq_au=res.frequency_au, amplitude=res.amplitude)


def fft_pipeline_panel(time_au, values, ft: FourierTransform | None = None, *,
                       label: str = "signal", freq_unit: str = "eV",
                       peak_band=None, fmax=None, title=None):
    """3×2 FFT-pipeline diagnostic panel. Returns a matplotlib Figure.

    Parameters
    ----------
    time_au, values : 1-D arrays (time in a.u., real signal).
    ft              : FourierTransform config (window/detrend/zero_pad/...). Default Hann.
    label           : signal name for axis labels (e.g. "dipole_z").
    freq_unit       : 'eV' (default) or 'au' for the frequency axis.
    peak_band       : optional (lo, hi) in `freq_unit` — search the dominant peak in
                      this physical band (avoids the DC bin) and annotate it.
    fmax            : optional upper frequency limit for the FFT panels (`freq_unit`).
    title           : suptitle.
    """
    ft = ft or FourierTransform()
    st = fft_stages(time_au, values, ft)
    ST.apply_theme()
    import matplotlib.pyplot as plt

    # fourier.py returns freq_au = rfftfreq (ordinary cycles/a.u.); the project
    # energy convention is the ANGULAR energy ħω = 2πf·E_h (see fourier_analysis.ipynb
    # `omega_eV` and pipeline/spectral_weight.py). Match it exactly.
    omega_au = st["freq_au"] * 2.0 * np.pi
    fx = omega_au * (HA_TO_EV if freq_unit == "eV" else 1.0)
    amp = st["amplitude"]
    funit = "eV" if freq_unit == "eV" else "a.u."

    # Comparison curve: the same pipeline (window/zero-pad/transient) with the
    # 'detrend' baseline, overlaid dashed on the FFT axes so the canonical
    # 'mean' choice is auditable (user verdict 2026-06-25). Skip if the primary
    # IS detrend (the curves would coincide).
    amp_cmp = None
    if ft.subtract != "detrend":
        ft_cmp = FourierTransform(
            window=ft.window, zero_pad=ft.zero_pad,
            smooth_sigma_bins=ft.smooth_sigma_bins,
            t_start_au=ft.t_start_au, subtract="detrend")
        amp_cmp = ft_cmp.transform(time_au, values).amplitude

    # dominant peak (in band if given, else skip the DC bin)
    if peak_band is not None:
        sel = (fx >= peak_band[0]) & (fx <= peak_band[1])
    else:
        sel = fx > (fx[1] if len(fx) > 1 else 0.0)
    ipk = None
    if sel.any():
        idx = np.where(sel)[0]
        ipk = idx[int(np.argmax(amp[idx]))]
    fpk = fx[ipk] if ipk is not None else None

    fig, axs = plt.subplots(3, 2, figsize=(7.0, 8.2), constrained_layout=True)

    # --- row 1: raw | detrended ---
    axs[0, 0].plot(st["t"], st["raw"], lw=0.8)
    axs[0, 0].set_title("1. raw signal")
    axs[0, 1].plot(st["t"], st["detrended"], lw=0.8, color="C1")
    axs[0, 1].set_title(f"2. baseline-removed  (subtract='{ft.subtract}')")
    for ax in axs[0]:
        ax.set_xlabel("t (a.u.)"); ax.set_ylabel(label)

    # --- row 2: windowed | padded ---
    env = np.max(np.abs(st["detrended"])) or 1.0
    axs[1, 0].plot(st["t"], st["windowed"], lw=0.8, color="C2")
    axs[1, 0].plot(st["t"], st["window"] * env, ls="--", lw=0.7, color="0.5",
                   label=f"{ft.window.name} envelope")
    axs[1, 0].set_title("3. windowed signal"); axs[1, 0].legend(fontsize=6, loc="upper right")
    axs[1, 1].plot(st["t_pad"], st["padded"], lw=0.6, color="C3")
    axs[1, 1].axvline(st["t"][-1], ls=":", lw=0.8, color="0.5",
                      label=f"zero-pad ×{ft.zero_pad}")
    axs[1, 1].set_title("4. zero-padded signal"); axs[1, 1].legend(fontsize=6, loc="upper right")
    for ax in axs[1]:
        ax.set_xlabel("t (a.u.)"); ax.set_ylabel(label)

    # --- row 3: FFT linear | log ---
    for col, (ax, logy) in enumerate(zip(axs[2], (False, True))):
        ax.plot(fx, amp, lw=0.9, color="C0", label=f"subtract='{ft.subtract}'")
        if amp_cmp is not None:
            ax.plot(fx, amp_cmp, lw=0.8, ls="--", color="0.45",
                    label="detrend (comparison)")
        if logy:
            ax.set_yscale("log"); ax.set_title("6. |FFT| — log scale")
        else:
            ax.set_title("5. |FFT| — linear scale")
        if peak_band is not None:
            ax.axvspan(peak_band[0], peak_band[1], color="C2", alpha=0.10)
        if fpk is not None:
            ax.axvline(fpk, ls="--", lw=0.8, color="C3")
            if not logy:
                ax.annotate(f"peak {fpk:.2f} {funit}", xy=(fpk, amp[ipk]),
                            xytext=(0.5, 0.9), textcoords="axes fraction",
                            fontsize=6, color="C3")
        ax.set_xlabel(rf"$\hbar\omega$ ({funit})"); ax.set_ylabel("|FFT| (coherent-gain norm.)")
        ax.set_xlim(0, fmax if fmax is not None else (fx.max()))
        ax.legend(fontsize=6, loc="upper right")

    fig.suptitle(title or f"FFT pipeline — {label}", fontsize=10)
    return fig
