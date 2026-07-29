"""Known-case tests for inqview.analysis.heuristics (electron-gas scales,
timescales, zero-point, norm/absorption). Analytic HEG identities are the
oracle — no run data needed (deps-clean, numpy only)."""
import numpy as np

from inqview.analysis import heuristics as H

HA_EV = 27.211386245988


def test_electron_gas_scales_analytic():
    # rs=4: kF = (9π/4)^(1/3)/rs = 1.919158/4 = 0.479790 a0^-1
    s = H.electron_gas_scales(4.0)
    assert np.isclose(s["kF"], 1.9191583 / 4.0, rtol=1e-6)
    assert np.isclose(s["vF"], s["kF"])                       # vF = kF (a.u.)
    assert np.isclose(s["EF_ha"], 0.5 * s["kF"] ** 2)
    # ω_p = sqrt(3/rs^3) for the HEG (Ha)
    assert np.isclose(s["omega_p_ha"], np.sqrt(3.0 / 4.0 ** 3), rtol=1e-6)
    # k_TF = sqrt(4 kF/π)
    assert np.isclose(s["k_TF"], np.sqrt(4.0 * s["kF"] / np.pi), rtol=1e-6)
    # Friedel wavelength π/kF
    assert np.isclose(s["lambda_F_friedel"], np.pi / s["kF"], rtol=1e-6)
    # n0 = 3/(4π rs^3)
    assert np.isclose(s["n0"], 3.0 / (4.0 * np.pi * 4.0 ** 3), rtol=1e-6)


def test_timescales_constant_velocity():
    t = H.projectile_timescales(z0=-22.0, v=2.711, slab_half=12.5, box_half=35.0)
    assert np.isclose(t["t_enter_slab_au"], 9.5 / 2.711, rtol=1e-9)
    assert np.isclose(t["t_exit_slab_au"], 34.5 / 2.711, rtol=1e-9)   # reach far face
    assert np.isclose(t["t_cross_au"], 25.0 / 2.711, rtol=1e-9)
    assert np.isclose(t["t_reach_box_edge_au"], 57.0 / 2.711, rtol=1e-9)


def test_zero_point_sigma_half():
    z = H.wp_zero_point(0.5)
    assert np.isclose(z["zero_point_ke_ha"], 3.0)                     # 3/(4·0.25)
    assert np.isclose(z["zero_point_ke_ev"], 3.0 * HA_EV)            # ≈ 81.6 eV
    assert np.isclose(z["sigma_charge"], 0.5 / np.sqrt(2.0))


def test_norm_absorption_split():
    Nt = np.array([83.0, 82.5, 82.13])
    Nwp = np.array([1.0, 0.6, 0.136])
    a = H.norm_absorption(Nt, Nwp)
    assert np.isclose(a["total_absorbed"], 0.87, atol=1e-9)
    assert np.isclose(a["wp_orbital_absorbed"], 0.864, atol=1e-9)
    assert np.isclose(a["bath_overflow_absorbed"], 0.87 - 0.864, atol=1e-9)
    assert np.isclose(a["wp_fraction_absorbed"], 0.864, atol=1e-9)


def test_spreading_factor():
    sz = np.array([0.354, 5.0, 14.6])
    sp = H.spreading(sz)
    assert np.isclose(sp["spread_factor"], 14.6 / 0.354, rtol=1e-9)
    assert np.isclose(sp["sigma_z_max"], 14.6)
