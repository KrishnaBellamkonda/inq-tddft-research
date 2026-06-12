"""Known-case tests for inqview.analysis.stopping_extract on synthetic tracks.

Build a synthetic electron_track.csv with a KNOWN stopping law, integrate the
decelerating trajectory, and confirm the extractor recovers S(v).
"""
from __future__ import annotations

import numpy as np

from inqview.analysis import stopping_extract as SE


def _write_track(tmp_path, t, z, v):
    p = tmp_path / "electron_track.csv"
    with open(p, "w") as fh:
        fh.write("step,time_au,x,y,z,vx,vy,vz\n")
        for i in range(len(t)):
            fh.write(f"{i},{t[i]},0,0,{z[i]},0,0,{v[i]}\n")
    return str(p)


def test_constant_stopping_recovered(tmp_path):
    """Constant S0: KE decreases linearly in s, so S(v)=S0 everywhere."""
    S0 = 0.05
    dt = 0.02
    # integrate dv/dt = -S0/v * (ds/dt)/v ... simpler: impose KE(s)=KE0 - S0*s
    # with m=1, v=sqrt(2 KE). March in s.
    v = 1.0
    z = 0.0
    ts, zs, vs = [], [], []
    for i in range(4000):
        ts.append(i * dt)
        zs.append(z)
        vs.append(v)
        ds = v * dt
        ke = 0.5 * v * v - S0 * ds
        if ke <= 1e-4:
            break
        v = np.sqrt(2 * ke)
        z += ds
    path = _write_track(tmp_path, np.array(ts), np.array(zs), np.array(vs))
    tr = SE.load_track(path)
    vv, SS = SE.stopping_vs_v(tr, transient_bohr=3.0, window=21)
    assert SS.size > 10
    # recovered S should be ~S0 across the trajectory
    assert np.allclose(SS, S0, rtol=0.05), (SS.min(), SS.max())


def test_linear_friction_recovered(tmp_path):
    """Friction S(v)=Q v: recover Q from the low-v slope of S(v)."""
    Q = 0.1
    dt = 0.02
    v = 1.0
    z = 0.0
    ts, zs, vs = [], [], []
    for i in range(8000):
        ts.append(i * dt)
        zs.append(z)
        vs.append(v)
        ds = v * dt
        ke = 0.5 * v * v - (Q * v) * ds
        if ke <= 1e-4:
            break
        v = np.sqrt(2 * ke)
        z += ds
    path = _write_track(tmp_path, np.array(ts), np.array(zs), np.array(vs))
    tr = SE.load_track(path)
    vv, SS = SE.stopping_vs_v(tr, transient_bohr=3.0, window=21)
    # S/v should be ~Q
    ratio = SS / vv
    assert np.allclose(ratio, Q, rtol=0.08), (ratio.min(), ratio.max())
