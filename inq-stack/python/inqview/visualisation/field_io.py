"""inqview.visualisation.field_io — the ONE canonical VTI loader.

Why this exists
---------------
inqkit VTIs are written in **physical order**: `inqkit::io::RealField3DWriter`
applies `fft_shift_index()` at write time and stamps `Origin = -L/2`, so array
index 0 already maps to the left-edge coordinate `-L/2` (NOT the FFT-natural
centre). Therefore **VTI data must never be `np.fft.fftshift`-ed** — doing so
swaps centre↔edge and silently produces flipped pictures (the recurring
"slab-at-the-edges" bug). Only LEED screen `.dat` files are FFT-natural and need
a shift; those have their own loader (`inqview.io.load_leed_pattern`).

Every GIF/slice/profile — notebooks and `make_*_postproc.py` — must load through
`load_vti` and use the returned coordinate axes. No hand-rolled VTK reads, no
hand-rolled fftshift.

Layer note (ADR 0003): VTI reading needs VTK, and `inqview.io` is contractually
numpy-only, so this canonical loader lives in the VTK-allowed `visualisation`
layer and is re-exported lazily as `inqview.load_vti`.
"""
from __future__ import annotations

from pathlib import Path
from typing import NamedTuple, Optional, Tuple

import numpy as np


class VtiField(NamedTuple):
    """A VTI loaded in physical order.

    data : np.ndarray, shape (nx, ny, nz), indexed [ix, iy, iz]
    x, y, z : 1-D cell-centred coordinate axes (Bohr), monotonically increasing
    origin, spacing : the VTI ImageData origin/spacing tuples (x, y, z)
    """
    data: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    origin: Tuple[float, float, float]
    spacing: Tuple[float, float, float]

    def xz_slice(self, y: float = 0.0) -> np.ndarray:
        """Return the [z, x] density slice nearest plane y=`y`, ready for
        imshow(origin='lower', extent=[x0,x1,z0,z1]). Transposed so rows are z."""
        iy = int(np.argmin(np.abs(self.y - y)))
        return self.data[:, iy, :].T  # (nz, nx): rows=z, cols=x


def load_vti(
    path: str | Path,
    *,
    array: Optional[str] = None,
    expect_centered_axis: Optional[str] = None,
    expect_tol_bohr: Optional[float] = None,
) -> VtiField:
    """Load an inqkit VTI in physical order (no fftshift) with a hard self-check.

    Parameters
    ----------
    path : VTI file.
    array : optional DataArray name. Default (None) takes the FIRST array, which
        is what every single-field inqkit VTI contains. Complex fields written by
        `inqkit::io::ComplexField3DWriter` carry TWO arrays
        (`<name>_real`, `<name>_imag`) and the default would silently return only
        the real part — pass the name explicitly, or use `load_complex_vti`.
    expect_centered_axis : optional, one of {'x','y','z'}. If given, assert that
        the planar-summed density profile along that axis peaks near coordinate 0
        (the box centre) — a loud failure if the index→coordinate mapping is
        wrong. Use for runs whose feature (slab/cluster) sits at the centre.
    expect_tol_bohr : tolerance for the centred-feature check (default: 8 grid
        spacings, generous — it only needs to catch a centre↔edge swap).

    Returns
    -------
    VtiField (physical order).
    """
    try:
        import vtk
        from vtk.util.numpy_support import vtk_to_numpy
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(f"VTK is required to read {path}") from exc

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    img = reader.GetOutput()
    nx, ny, nz = img.GetDimensions()
    origin = tuple(float(v) for v in img.GetOrigin())
    spacing = tuple(float(v) for v in img.GetSpacing())

    point_data = img.GetPointData()
    if array is None:
        vtk_array = point_data.GetArray(0)
    else:
        vtk_array = point_data.GetArray(array)
        if vtk_array is None:
            available = [point_data.GetArrayName(i)
                         for i in range(point_data.GetNumberOfArrays())]
            raise KeyError(
                f"{path}: no DataArray named {array!r}; available: {available}")
    flat = vtk_to_numpy(vtk_array).astype(np.float64, copy=False)
    # VTK ImageData is x-fastest: reshape (nz,ny,nx) then transpose to (nx,ny,nz).
    data = flat.reshape((nz, ny, nx)).transpose(2, 1, 0)

    ox, oy, oz = origin
    sx, sy, sz = spacing
    # Cell-centred sample coordinates (matches the C++ writer / pipeline/density).
    x = ox + (np.arange(nx) + 0.5) * sx
    y = oy + (np.arange(ny) + 0.5) * sy
    z = oz + (np.arange(nz) + 0.5) * sz

    # ---- HARD self-check: the axis/dim invariants ------------------------
    assert data.shape == (nx, ny, nz), (
        f"reshape/transpose mismatch: data {data.shape} vs dims {(nx, ny, nz)}")
    for name, ax, s in (("x", x, sx), ("y", y, sy), ("z", z, sz)):
        assert s > 0.0, f"non-positive spacing on {name}: {s}"
        assert ax[1] > ax[0], f"{name} axis not increasing"
    # Physical order ⇒ first sample sits at the left edge, not the centre.
    assert abs(x[0] - (ox + 0.5 * sx)) < 1e-9, "x[0] is not the physical left edge"

    if expect_centered_axis is not None:
        # Discriminate a CENTRED feature (mass in the inner half of the box) from
        # an edge-split one (the centre↔edge swap a wrong fftshift produces). Both
        # are symmetric, so argmax/centre-of-mass cannot tell them apart — compare
        # inner-half vs outer-half |n| mass instead. Robust to Friedel peaks.
        axis_index = {"x": 0, "y": 1, "z": 2}[expect_centered_axis]
        coord = (x, y, z)[axis_index]
        other = tuple(a for a in (0, 1, 2) if a != axis_index)
        profile = np.abs(data).sum(axis=other)  # planar-summed |n| along the axis
        L_axis = nx * sx if axis_index == 0 else (ny * sy if axis_index == 1 else nz * sz)
        inner = np.abs(coord) < 0.25 * L_axis
        inner_mass = float(profile[inner].sum())
        outer_mass = float(profile[~inner].sum())
        assert inner_mass > outer_mass, (
            f"index→coordinate mapping looks WRONG on {expect_centered_axis}: "
            f"inner-half |n| mass {inner_mass:.3g} ≤ outer-half {outer_mass:.3g} — "
            f"the feature sits at the EDGES, not the centre. Did something "
            f"np.fft.fftshift a physical-order VTI?")

    return VtiField(data=data, x=x, y=y, z=z, origin=origin, spacing=spacing)


class ComplexVtiField(NamedTuple):
    """A complex orbital VTI loaded in physical order.

    data : complex np.ndarray, shape (nx, ny, nz), indexed [ix, iy, iz]
    x, y, z : 1-D cell-centred coordinate axes (Bohr), monotonically increasing
    """
    data: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    origin: Tuple[float, float, float]
    spacing: Tuple[float, float, float]


def load_complex_vti(
    path: str | Path,
    *,
    field_name: str = "wavefunction",
    expect_centered_axis: Optional[str] = None,
) -> ComplexVtiField:
    """Load a complex field written by `inqkit::io::ComplexField3DWriter`.

    That writer emits ONE .vti carrying two DataArrays, `<field_name>_real` and
    `<field_name>_imag`. Loading it with the plain `load_vti` silently returns
    only the real part (it takes array 0), which for a drifting wavepacket
    ψ ∝ e^{i k₀ z} is a cosine-modulated Gaussian — it *looks* like a plausible
    orbital, so the mistake is not caught by eye. Use this instead.

    Physical order is preserved, exactly as for `load_vti` — do NOT fftshift the
    result. See `.claude/rules/vti-coordinate-mapping.md`.
    """
    re = load_vti(path, array=f"{field_name}_real",
                  expect_centered_axis=expect_centered_axis)
    im = load_vti(path, array=f"{field_name}_imag")
    return ComplexVtiField(
        data=re.data + 1j * im.data,
        x=re.x, y=re.y, z=re.z, origin=re.origin, spacing=re.spacing)


def kz_marginal(field: ComplexVtiField) -> Tuple[np.ndarray, np.ndarray]:
    """Return (k_z, P(k_z)) — the normalised z-momentum marginal of an orbital.

    P(k_z) = Σ_{k_x,k_y} |ψ̃(k)|², normalised to ∫P dk_z = 1.

    For a Gaussian wavepacket ψ ∝ exp(-|r-b|²/2σ²) e^{i k₀·r} this marginal is
    EXACTLY Gaussian, N(k₀, σ_p²) with σ_p = 1/(√2 σ) — which is what makes it
    the right observable for a "is the momentum distribution still Gaussian?"
    check. The *radial* distribution n(|k|) is not Gaussian for a drifting
    packet, so `inqkit`'s radial MomentumDistribution cannot answer this.

    THE ORDERING TRAP. inqkit VTIs are stored in PHYSICAL order (index 0 ↔ −L/2),
    but `np.fft.fftn` expects FFT-natural order (index 0 ↔ the origin). So the
    field must be `ifftshift`-ed BEFORE transforming. This is not a violation of
    the never-fftshift-a-VTI rule — that rule governs *display* of real-space
    data; here we are converting to the FFT's own index convention in order to
    transform, and the real-space array itself is never reordered for plotting.
    Getting this wrong shifts the recovered ⟨k_z⟩ by a half-cell phase ramp and
    scrambles the profile.
    """
    psi = np.asarray(field.data)
    nx, ny, nz = psi.shape
    dz = float(field.spacing[2])

    # physical order -> FFT-natural order, transform, then order the k axes.
    psi_fft_order = np.fft.ifftshift(psi)
    psi_k = np.fft.fftn(psi_fft_order)

    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=dz)
    weight = (np.abs(psi_k) ** 2).sum(axis=(0, 1))

    order = np.argsort(kz)
    kz, weight = kz[order], weight[order]

    dk = kz[1] - kz[0]
    total = weight.sum() * dk
    if not total > 0.0:
        raise ValueError("kz_marginal: zero total momentum weight")
    return kz, weight / total


def gaussian_fit_quality(
    kz: np.ndarray,
    prob: np.ndarray,
    *,
    k0: float,
    sigma_p: float,
) -> dict:
    """Compare a measured k_z marginal against the ANALYTIC N(k0, sigma_p²).

    Deliberately compares against the *expected* Gaussian rather than a
    best-fit one: a packet that has been deformed into a narrower or shifted
    Gaussian would still fit a free 2-parameter Gaussian beautifully while being
    the wrong packet. Both are reported — `r2_analytic` is the one that matters,
    `r2_bestfit` says whether what remains is Gaussian-SHAPED at all.

    Returns mean/std/skewness/excess-kurtosis of the measured marginal (a true
    Gaussian has skew = 0, excess kurtosis = 0) plus the two R² values and the
    L1 residual against the analytic curve.
    """
    kz = np.asarray(kz, dtype=float)
    p = np.asarray(prob, dtype=float)
    dk = kz[1] - kz[0]

    mean = float((kz * p).sum() * dk)
    var = float(((kz - mean) ** 2 * p).sum() * dk)
    std = float(np.sqrt(var))
    skew = float(((kz - mean) ** 3 * p).sum() * dk / std ** 3) if std > 0 else float("nan")
    kurt = float(((kz - mean) ** 4 * p).sum() * dk / std ** 4 - 3.0) if std > 0 else float("nan")

    def _gauss(mu, sd):
        return np.exp(-0.5 * ((kz - mu) / sd) ** 2) / (sd * np.sqrt(2.0 * np.pi))

    def _r2(model):
        ss_res = float(((p - model) ** 2).sum())
        ss_tot = float(((p - p.mean()) ** 2).sum())
        return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    analytic = _gauss(k0, sigma_p)
    bestfit = _gauss(mean, std)

    return {
        "mean_kz": mean,
        "std_kz": std,
        "skewness": skew,
        "excess_kurtosis": kurt,
        "r2_analytic": _r2(analytic),
        "r2_bestfit": _r2(bestfit),
        "l1_residual_analytic": float(np.abs(p - analytic).sum() * dk),
        "expected_mean_kz": float(k0),
        "expected_std_kz": float(sigma_p),
    }


def kz_kperp_map(
    field: ComplexVtiField,
    *,
    n_kperp_bins: Optional[int] = None,
    kperp_max: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (k_z, k_perp, P) — the 2-D momentum density of one orbital.

    ``P[i, j]`` is the norm carried by grid points with longitudinal momentum
    ``k_z[i]`` and transverse momentum ``|k_perp|`` in bin ``j``. Normalised so
    ``P.sum() == 1``.

    WHY THIS AND NOT THE RADIAL ``n(|k|)`` HISTOGRAM. Binning by ``|k|`` alone
    folds the drift direction into the same coordinate as the transverse spread,
    so a packet that is decelerating (mean ``k_z`` falling) and one that is
    heating sideways (``k_perp`` growing) both show up as "the peak moved left".
    They are different physics and this function separates them: deceleration
    moves weight along the ``k_z`` axis, transverse heating moves it up the
    ``k_perp`` axis.

    NO BINNING ALONG z. ``k_z`` is a native FFT grid axis, so it is used exactly
    as-is (all ``nz`` values, sorted ascending). Only ``k_perp`` is binned,
    because it is a derived radial coordinate. This matters: a Gaussian packet
    of width ``sigma_p`` is carried by only ~10^2 grid points out of ~10^6, so
    any binning finer than the grid spacing produces a spiky comb that looks
    like structure and is not. The default bin width is one transverse grid
    spacing ``2*pi/L_xy``, which is the finest width the data supports.

    THE JACOBIAN IS INCLUDED, deliberately. ``P`` is a probability over
    ``(k_z, k_perp)``, i.e. the shell sum, so for an isotropic-in-plane Gaussian
    the ``k_perp`` marginal is a RAYLEIGH distribution peaking at
    ``k_perp = sigma_p``, not a Gaussian peaking at 0. That is the honest
    "how much norm sits at this transverse momentum" and is what a difference
    map should be built from. Divide by ``k_perp`` if the underlying 3-D density
    is wanted instead.

    DO NOT TAKE TRANSVERSE MOMENTS FROM THIS MAP. Along ``k_z`` the axis is the
    exact FFT grid, so ``sum(k_z^2 P)`` is exact. Along ``k_perp`` every point in
    a bin is assigned that bin's CENTRE, and because a Rayleigh tail falls
    steeply across a bin the mass really sits below the centre — so
    ``sum(k_perp^2 P)`` carries a systematic POSITIVE bias of order a few per
    cent (measured: +6.3 % at ``sigma_p / dk_perp = 3.0``, +9.4 % at 2.4, and it
    does not vanish with finer bins because finer bins are exactly what the grid
    cannot support). The unbinned moments are exact and are already written
    every step to ``wp_momentum_stats.csv`` — use those for numbers, and this map
    for SHAPE and for DIFFERENCES, where the bias largely cancels.

    Ordering: as in ``kz_marginal``, the field is ``ifftshift``-ed from physical
    order into FFT-natural order before transforming. See
    ``.claude/rules/vti-coordinate-mapping.md``.
    """
    psi = np.asarray(field.data)
    nx, ny, nz = psi.shape
    dx, dy, dz = (float(s) for s in field.spacing)

    psi_k = np.fft.fftn(np.fft.ifftshift(psi))
    w = np.abs(psi_k) ** 2
    total = w.sum()
    if not total > 0.0:
        raise ValueError("kz_kperp_map: zero total momentum weight")
    w = w / total

    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=dy)
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=dz)
    kperp = np.sqrt(kx[:, None] ** 2 + ky[None, :] ** 2)      # (nx, ny)

    if kperp_max is None:
        kperp_max = float(kperp.max())
    if n_kperp_bins is None:
        # one transverse grid spacing per bin — the finest the data supports.
        dk_xy = float(abs(kx[1] - kx[0]))
        n_kperp_bins = max(1, int(np.ceil(kperp_max / dk_xy)))

    edges = np.linspace(0.0, kperp_max, n_kperp_bins + 1)
    centers = 0.5 * (edges[1:] + edges[:-1])
    # np.digitize gives 1..n_kperp_bins inside range; clip the top edge in.
    idx = np.clip(np.digitize(kperp.ravel(), edges) - 1, 0, n_kperp_bins - 1)

    flat = w.reshape(nx * ny, nz)
    out = np.zeros((nz, n_kperp_bins), dtype=float)
    for j in range(n_kperp_bins):
        sel = idx == j
        if sel.any():
            out[:, j] = flat[sel, :].sum(axis=0)

    order = np.argsort(kz)
    return kz[order], centers, out[order, :]
