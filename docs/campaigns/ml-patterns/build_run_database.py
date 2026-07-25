#!/usr/bin/env python3
"""
build_run_database.py -- Reproducibility-grade TDDFT run database builder.

Scans ResearchProject/systems/{jellium,localised_jellium,coronene,
cylindrical_jellium,graphene,vacuum} for run_summary.txt files, dedups to one
row per real run, parses the (3+) heterogeneous summary formats with per-system
field maps + sub-parsers for compound lines, applies config/run.cpp fallbacks,
computes derived electron-gas physics, detects observable channels by PATH, and
links classical<->wp twins.

Outputs (all under docs/):
  run_database.csv          -- canonical wide table, one row per run
  run_database.json         -- nested mirror (params + observable file lists)
  run_database_columns.json -- machine-readable data dictionary

Self-contained: stdlib + numpy only. Idempotent.

Run with the project venv:
  /local/data/public/skcb2/tddft/venv/bin/python3 \
      docs/campaigns/ml-patterns/build_run_database.py

NEVER invents values: absent parameters are written as the literal token NULL.

Derived-physics formulas use the standard homogeneous-electron-gas relations
(Ashcroft & Mermin, "Solid State Physics", Ch. 1-2; atomic/Hartree units,
hbar = m_e = e = 1):
    k_F   = (3 pi^2 n)^(1/3)            [A&M eq. 2.21]
    E_F   = k_F^2 / 2     (Hartree)     [A&M eq. 2.26, with hbar=m=1]
    v_F   = k_F           (a.u.)        [A&M eq. 2.24, p_F = hbar k_F, v=p/m]
    omega_p = sqrt(4 pi n) (Hartree)    [A&M eq. 1.38 in Gaussian/atomic units]
    r_s   = (3/(4 pi n))^(1/3)          [A&M eq. 1.2]
sigma_pot = sigma_WP / sqrt(2)  (classical Gaussian-potential width; project
    convention -- classical charge std equals WP density std).
"""

import csv
import json
import math
import os
import re
import subprocess
import sys

import numpy as np

HA_TO_EV = 27.211386  # 1 Hartree in eV (CODATA)
SQRT2 = math.sqrt(2.0)

ROOT = "/local/data/public/skcb2/tddft"
SYSTEMS_DIR = os.path.join(ROOT, "ResearchProject", "systems")
SYSTEMS = ["jellium", "localised_jellium", "coronene",
           "cylindrical_jellium", "graphene", "vacuum"]
OUT_DIR = os.path.join(ROOT, "docs")

NULL = "NULL"

# --------------------------------------------------------------------------
# Low-level parsing
# --------------------------------------------------------------------------

# Match one "key = value" assignment, value running until 2+ spaces precede the
# next "word =" token, or end of line.  Handles both the aligned structured
# format (one assignment/line) and the terse multi-assignment lines.
KV_RE = re.compile(r'([A-Za-z_][A-Za-z0-9_\[\]]*)\s*=\s*(.*?)(?=\s{2,}[A-Za-z_][A-Za-z0-9_\[\]]*\s*=|$)')


def parse_summary_text(text, want_raw=False):
    """Return dict of key -> value (last assignment wins) from a summary file.

    Section header lines ('1. Run identity', '----'), blank lines, and the
    'RUN SUMMARY' banner contain no '=' and are skipped automatically.

    If want_raw, also return d_raw mapping each line's FIRST key -> the entire
    remainder of that line (un-tokenised). Compound lines whose assignments are
    separated by a SINGLE space (e.g. cylindrical 'geometry = annular_tube
    R_in=5 R_out=13 L_z=10') need the raw remainder so the dedicated
    sub-parsers can extract every field.
    """
    d = {}
    d_raw = {}
    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if "=" not in line:
            continue
        # Skip pure separator lines.
        if set(line.strip()) <= set("=-"):
            continue
        first = None
        for m in KV_RE.finditer(line):
            key = m.group(1).strip()
            val = m.group(2).strip()
            if val == "":
                continue
            d[key] = val  # last wins (handles run_completed false->true)
            if first is None:
                first = key
        if first is not None:
            # remainder = everything after the first 'key ='
            mm = re.match(r'\s*[A-Za-z_][A-Za-z0-9_\[\]]*\s*=\s*(.*)$', line)
            if mm:
                d_raw[first] = mm.group(1).strip()
    if want_raw:
        return d, d_raw
    return d


def first_float(s):
    """Extract the first floating-point number from a string, or None."""
    if s is None:
        return None
    m = re.search(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', str(s))
    return float(m.group(0)) if m else None


def all_floats(s):
    return [float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', str(s))]


def parse_bool(s):
    if s is None:
        return None
    t = str(s).strip().lower()
    if t in ("true", "yes", "1", "on"):
        return True
    if t in ("false", "no", "0", "off"):
        return False
    return None


# --------------------------------------------------------------------------
# Compound-line sub-parsers
# --------------------------------------------------------------------------

def parse_cell(s):
    """Parse a cell spec into (x, y, z) in Bohr.

    Handles: '50^3 (cubic, periodic)', '50^3', '40 x 40 x 10', '50x50x90',
    '34.77 34.77 59.90', '50 x 50 x 70 (orthorhombic, periodic)',
    '20 22 60'.  Returns (None, None, None) if unparseable.
    """
    if s is None:
        return (None, None, None)
    s = str(s).strip()
    m = re.match(r'\s*([-+]?\d*\.?\d+)\s*\^\s*3', s)
    if m:
        v = float(m.group(1))
        return (v, v, v)
    # split on 'x' (with optional surrounding spaces) or whitespace
    body = re.split(r'\(', s)[0]  # drop trailing "(cubic, periodic)"
    parts = re.split(r'\s*[xX]\s*|\s+', body.strip())
    nums = []
    for p in parts:
        f = first_float(p)
        if f is not None:
            nums.append(f)
    if len(nums) >= 3:
        return (nums[0], nums[1], nums[2])
    if len(nums) == 1:
        return (nums[0], nums[0], nums[0])
    return (None, None, None)


def parse_projectile(s):
    """Sub-parse a 'projectile =' compound description.

    Returns dict with possible keys: projectile_kind, sigma_pot, sigma_wp,
    mass_au, ion_dynamics, wp_k0, energy_ev, projectile_charge.
    Examples:
      'classical electron (custom UPF + mass override)'
      'classical Gaussian-e ion (sigma_pot 0.354, mass m_e, ehrenfest)'
      'classical electron, erf Gaussian, m_e, ehrenfest'
      'electron WAVEPACKET sigma_WP 0.5 k0 0.3 (quantum)'
      'wavepacket sigma 0.5 E 340.142328125 eV k0 5'
    """
    out = {}
    if s is None:
        return out
    low = s.lower()
    if "wavepacket" in low or low.strip().startswith("wave"):
        out["projectile_kind"] = "electron wavepacket"
        out["sim_type"] = "wp"
    elif "classical" in low:
        out["projectile_kind"] = s.strip()
        out["sim_type"] = "classical"
    else:
        out["projectile_kind"] = s.strip()
    # sigma_pot (classical) -- 'sigma_pot 0.354'
    m = re.search(r'sigma_pot\s+([-+]?\d*\.?\d+)', s, re.I)
    if m:
        out["sigma_pot"] = float(m.group(1))
    # sigma_WP or sigma (wavepacket) -- 'sigma_WP 0.5' or 'sigma 0.5'
    m = re.search(r'sigma_?wp\s+([-+]?\d*\.?\d+)', s, re.I)
    if not m:
        m = re.search(r'(?<!_)\bsigma\s+([-+]?\d*\.?\d+)', s, re.I)
    if m and "sigma_pot" not in out:
        out["sigma_wp"] = float(m.group(1))
    # k0
    m = re.search(r'\bk0\s+([-+]?\d*\.?\d+)', s, re.I)
    if m:
        out["wp_k0"] = float(m.group(1))
    # E ... eV
    m = re.search(r'\bE\s+([-+]?\d*\.?\d+)\s*ev', s, re.I)
    if m:
        out["energy_ev"] = float(m.group(1))
    # mass
    if re.search(r'm_e\b', s):
        out["mass_au"] = 1.0
    # dynamics
    if "ehrenfest" in low:
        out["ion_dynamics"] = "ehrenfest"
    return out


def parse_geometry(s):
    """Sub-parse a 'geometry =' compound line (cylindrical).
    'annular_tube  R_in=5 R_out=13 L_z=10 (periodic)'
    Returns dict: geometry_kind, R_in, R_out, L_z.
    """
    out = {}
    if s is None:
        return out
    head = re.split(r'\s+', s.strip())[0]
    out["geometry_kind"] = head
    for key in ("R_in", "R_out", "L_z"):
        m = re.search(key + r'\s*=?\s*([-+]?\d*\.?\d+)', s, re.I)
        if m:
            out[key] = float(m.group(1))
    return out


def parse_cap(s):
    """Sub-parse a 'cap =' compound line. Returns cap_form, cap_eta,
    cap_width_frac (frac only when a fractional width is given).
    Examples:
      'on (two-sided sin2, eta -0.7 Ha, 10 Bohr/side, region +/-35..+/-45 ...)'
      'on (sin2 eta -0.5 mid +/-0.425 width 0.15)'
      '1'  (graphene uses cap = 1 with separate eta_Ha key)
    """
    out = {}
    if s is None:
        return out
    low = s.lower()
    for form in ("sin2", "monomial", "mono", "mask"):
        if form in low:
            out["cap_form"] = "sin2" if form == "sin2" else form
            break
    m = re.search(r'eta\s+([-+]?\d*\.?\d+)', s, re.I)
    if m:
        out["cap_eta"] = float(m.group(1))
    # fractional width: 'width 0.15' or 'width_frac 0.15' (a fraction < 1)
    m = re.search(r'width(?:_frac)?\s+([-+]?\d*\.?\d+)', s, re.I)
    if m:
        w = float(m.group(1))
        if w < 1.0:  # a fraction, not a Bohr length
            out["cap_width_frac"] = w
    return out


# --------------------------------------------------------------------------
# UPF PP_LOCAL inspection -- classical projectile potential form & width
# --------------------------------------------------------------------------
#
# UPF PP_LOCAL stores the local radial potential V(r) on the PP_R mesh, in
# RYDBERG, so for a Gaussian charge of charge std sigma_pot and |charge| = 1:
#     V(r) = (2/r) * erf( r / (sqrt(2) * sigma_pot) )      (Ry, magnitude)
#     V(0) = 2 * sqrt(2/pi) / sigma_pot                    (finite at the origin)
# A *bare/ONCV* projectile instead has a -Z/r long-range tail and a real
# pseudo core; it does NOT fit the erf form to machine precision. We therefore
# classify by FITTING the erf-Gaussian model (sigma_pot inferred from V(0)) and
# checking the residual over an intermediate-r window.  This is data-driven and
# never trusts the filename/header (per CONTEXT.md: verify Gaussian-ness by the
# PP_LOCAL V(r) DATA; the electron_gaussian_*.upf carry STALE Coulomb headers).
SQRT_2_OVER_PI = math.sqrt(2.0 / math.pi)
_UPF_CACHE = {}


def parse_upf_vr(path):
    """Parse (PP_R, PP_LOCAL) value arrays from a UPF. Returns (R, V) lists of
    equal length, or (None, None) if unreadable / malformed. Cached per path."""
    if path in _UPF_CACHE:
        return _UPF_CACHE[path]
    try:
        with open(path, "r", errors="ignore") as f:
            txt = f.read()
    except OSError:
        _UPF_CACHE[path] = (None, None)
        return (None, None)

    def block(tag):
        m = re.search(r'<%s\b[^>]*>(.*?)</%s>' % (tag, tag), txt, re.S)
        if not m:
            return None
        return [float(x) for x in re.findall(
            r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', m.group(1))]

    R = block("PP_R")
    V = block("PP_LOCAL")
    res = (R, V) if (R and V and len(R) == len(V)) else (None, None)
    _UPF_CACHE[path] = res
    return res


def classify_upf_potential(path):
    """Inspect PP_LOCAL V(r). Returns (form, sigma_pot, V0_stored).

    form = 'gaussian' if V(r) is the erf-Coulomb of a Gaussian charge (erf model
    fits to ~machine precision over 0.3<=r<=4 Bohr); 'coulombic' if a bare/ONCV
    pseudo (poor erf fit, -Z/r tail).  sigma_pot is the inferred charge std (Bohr)
    for gaussian, None for coulombic.  (None, None, None) if unreadable.
    """
    R, V = parse_upf_vr(path)
    if R is None:
        return (None, None, None)
    V0 = abs(V[0])
    if V0 <= 0:
        return (None, None, None)
    sigma_pot = 2.0 * SQRT_2_OVER_PI / V0
    res = []
    for r, v in zip(R, V):
        if 0.3 <= r <= 4.0:
            model = (2.0 / r) * math.erf(r / (SQRT2 * sigma_pot))
            res.append(abs(abs(v) - model) / max(abs(model), 1e-9))
    maxres = max(res) if res else 1.0
    # true Gaussian-charge UPFs (generated by the erf formula) fit to 0.0000;
    # the ONCV pseudo gives ~0.025. 5e-3 cleanly separates the two families.
    if maxres < 5e-3:
        return ("gaussian", sigma_pot, V0)
    return ("coulombic", None, V0)


_UPF_BASENAME_INDEX = None


def _upf_basename_index():
    """Map every electron*.upf basename found under any system's tree to its
    absolute path (built once).  Lets a run that references a sibling system's
    shared pseudopotential (e.g. localised run.cpp -> jellium's sigma0p35.upf)
    still resolve to the real file for V(r) inspection."""
    global _UPF_BASENAME_INDEX
    if _UPF_BASENAME_INDEX is not None:
        return _UPF_BASENAME_INDEX
    idx = {}
    for dirpath, _, filenames in os.walk(SYSTEMS_DIR):
        for fn in filenames:
            if fn.endswith(".upf"):
                idx.setdefault(fn, os.path.join(dirpath, fn))
    _UPF_BASENAME_INDEX = idx
    return idx


def _resolve_upf_path(p, system, anchor_dirs):
    """Resolve a (possibly relative) .upf path string to an existing file.
    Tries: as-is; joined to each anchor dir; <system>/shared/pseudopotentials/
    <basename>; then a global basename index across all systems. Returns an
    absolute path or None."""
    if not p:
        return None
    cands = [p]
    for a in anchor_dirs:
        if a:
            cands.append(os.path.normpath(os.path.join(a, p)))
    cands.append(os.path.join(SYSTEMS_DIR, system, p))
    cands.append(os.path.join(SYSTEMS_DIR, system, "shared",
                              "pseudopotentials", os.path.basename(p)))
    for c in cands:
        if c and os.path.isfile(c):
            return os.path.abspath(c)
    # global fallback: same basename anywhere under SYSTEMS_DIR
    g = _upf_basename_index().get(os.path.basename(p))
    return os.path.abspath(g) if g else None


def _grep_upf(src):
    """Return all '*.upf' tokens referenced in a source/config string."""
    return re.findall(r'([^\s"\'()<>]+\.upf)', src)


def locate_classical_upf(system, run_dir, run_cpp, cfg_hdr, d):
    """Best-effort projectile UPF path for a CLASSICAL run, by descending
    priority: summary pseudopotential field -> config header -> run.cpp ->
    system scripts/*.cpp scan (token-scored against the run name).  Returns an
    absolute path to an existing .upf, or None."""
    # 1. summary
    for key in ("projectile_upf", "pseudopotential", "psp"):
        v = d.get(key)
        if v:
            mm = re.search(r'(\S+\.upf)', str(v))
            if mm:
                r = _resolve_upf_path(mm.group(1), system, [run_dir])
                if r:
                    return r
    # 2. config header  3. run.cpp
    for srcfile in (cfg_hdr, run_cpp):
        if srcfile and os.path.isfile(srcfile):
            try:
                src = open(srcfile, "r", errors="ignore").read()
            except OSError:
                continue
            for p in _grep_upf(src):
                if "electron" not in os.path.basename(p).lower():
                    continue
                r = _resolve_upf_path(p, system, [os.path.dirname(srcfile)])
                if r:
                    return r
    # 4. scan the system's scripts/*.cpp for electron*.upf refs, token-scored
    search_root = os.path.join(SYSTEMS_DIR, system, "scripts")
    if not os.path.isdir(search_root):
        search_root = os.path.join(SYSTEMS_DIR, system)
    run_toks = set(re.split(r'[/_\-]', make_run_name(system, run_dir).lower()))
    best, best_score = None, -1
    for dirpath, _, filenames in os.walk(search_root):
        comps = dirpath.split(os.sep)
        if "build" in comps or "_deps" in comps or "CMakeFiles" in comps:
            continue
        for fn in filenames:
            if not (fn.endswith(".cpp") or fn.endswith(".hpp")):
                continue
            try:
                src = open(os.path.join(dirpath, fn), "r", errors="ignore").read()
            except OSError:
                continue
            for p in _grep_upf(src):
                if "electron" not in os.path.basename(p).lower():
                    continue
                r = _resolve_upf_path(p, system, [dirpath])
                if not r:
                    continue
                cdir_toks = set(re.split(r'[/_\-]',
                                         os.path.relpath(dirpath, search_root).lower()))
                score = len(cdir_toks & run_toks)
                if score > best_score:
                    best, best_score = r, score
    return best


# --------------------------------------------------------------------------
# Run discovery & dedup
# --------------------------------------------------------------------------

def find_summaries(system):
    base = os.path.join(SYSTEMS_DIR, system)
    out = []
    for dirpath, dirnames, filenames in os.walk(base):
        # skip CMake build trees (a path component exactly 'build', or _deps/
        # CMakeFiles) -- never contain real runs. NB 'build_smoke'/'buildsmoke'
        # ARE real smoke runs and must NOT be excluded.
        comps = dirpath.split(os.sep)
        if "build" in comps or "_deps" in comps or "CMakeFiles" in comps:
            continue
        if "run_summary.txt" in filenames:
            out.append(os.path.join(dirpath, "run_summary.txt"))
    return out


def run_key_for(summary_path):
    """Map a summary path to its canonical run directory (dedup key).

    Collapses a trailing 'raw' wrapper and a generic immediate 'results'
    wrapper; a NAMED results-subdir (e.g. results/p5_wp_v5p0) is itself the
    run identity.
    """
    sdir = os.path.dirname(summary_path)
    if os.path.basename(sdir) == "raw":
        sdir = os.path.dirname(sdir)
    if os.path.basename(sdir) == "results":
        run_dir = os.path.dirname(sdir)
    else:
        run_dir = sdir
    return run_dir


def make_run_name(system, run_dir):
    rel = os.path.relpath(run_dir, os.path.join(SYSTEMS_DIR, system))
    # tidy: drop '/results' segments and collapse consecutive duplicate parts
    parts = [p for p in rel.split(os.sep) if p != "results"]
    cleaned = []
    for p in parts:
        if not cleaned or cleaned[-1] != p:
            cleaned.append(p)
    return "/".join(cleaned)


# --------------------------------------------------------------------------
# Fallback: config header & run.cpp discovery, propagator extraction
# --------------------------------------------------------------------------

def find_run_cpp(system, run_dir, sim_type):
    """Best-effort location of the run.cpp that produced a run."""
    # 1. directly in the run dir (jellium convention)
    cand = os.path.join(run_dir, "run.cpp")
    if os.path.isfile(cand):
        return cand
    # 2. walk up looking for a sibling run.cpp (localised scripts/<x>/<role>/)
    d = run_dir
    for _ in range(5):
        c = os.path.join(d, "run.cpp")
        if os.path.isfile(c):
            return c
        d = os.path.dirname(d)
        if not d.startswith(SYSTEMS_DIR):
            break
    # 3. scripts tree matched by role (classical/wp/gs) for cylindrical etc.
    role = {"classical": "classical", "wp": "wp", "free_wp": "wp",
            "baseline": "gs"}.get(sim_type)
    if role:
        for dirpath, _, filenames in os.walk(os.path.join(SYSTEMS_DIR, system)):
            if "build" in dirpath or "_deps" in dirpath:
                continue
            if "run.cpp" in filenames and os.sep + role in dirpath + os.sep:
                return os.path.join(dirpath, "run.cpp")
    return None


def extract_propagator(run_cpp_path, summary_dict):
    """Determine propagator. run.cpp is authoritative; graphene summaries state
    it directly. INQ default (no explicit setter) is ETRS -- confirmed by the
    in-repo comments '// explicit; default' (run_positive_ion) and
    '// else default ETRS' (vacuum cap_probe)."""
    if "propagator" in summary_dict:
        v = summary_dict["propagator"].lower()
        if "crank" in v or v.strip() == "cn":
            return "crank_nicolson"
        if "etrs" in v:
            return "etrs"
    if run_cpp_path and os.path.isfile(run_cpp_path):
        try:
            with open(run_cpp_path, "r", errors="ignore") as f:
                src = f.read()
        except OSError:
            src = ""
        if re.search(r'\.crank_nicolson\s*\(', src):
            return "crank_nicolson"
        if re.search(r'\.etrs\s*\(', src):
            return "etrs"
        if re.search(r'real_time::propagate\s*\(', src):
            return "etrs"  # INQ default propagator (see docstring)
    return None


def find_config_header(system, run_cpp_path):
    """Resolve a shared/configs header #included by the run.cpp (provenance)."""
    if not run_cpp_path or not os.path.isfile(run_cpp_path):
        return None
    try:
        with open(run_cpp_path, "r", errors="ignore") as f:
            src = f.read()
    except OSError:
        return None
    for m in re.finditer(r'#include\s+"([^"]+\.hpp)"', src):
        inc = m.group(1)
        base = os.path.basename(inc)
        cand = os.path.join(SYSTEMS_DIR, system, "shared", "configs", base)
        if os.path.isfile(cand):
            return cand
        cand2 = os.path.normpath(os.path.join(os.path.dirname(run_cpp_path), inc))
        if "configs" in cand2 and os.path.isfile(cand2):
            return cand2
    return None


# --------------------------------------------------------------------------
# Observable channel detection (by PATH)
# --------------------------------------------------------------------------

# channel -> list of candidate basenames (dirs hold .vti frames; files are csv)
VTI_CHANNELS = {
    "density_total_vti": ["density_total", "density_rt_total"],
    "density_system_vti": ["density_system", "density_rt_system", "density_rt_target"],
    "density_wp_vti": ["density_wp", "density_rt_wp"],
    "density_delta_vti": ["density_delta", "density_rt_delta"],
    "wp_wavefunction_vti": ["wavefunction_wp", "density_rt_wp_wf"],
}
FILE_CHANNELS = {
    "observables_csv": ["observables.csv"],
    "wp_momentum_stats": ["wp_momentum_stats.csv"],
    "wp_realspace_stats": ["wp_real_space_stats.csv"],
    "state_energies": ["state_energies.csv"],
    "occupations": ["occupations_vs_time.csv", "occupations.csv"],
    "momentum_distribution": ["momentum_distribution.csv"],
    "gamma_transitions": ["gamma_transitions.csv"],
    "electron_track": ["electron_track.csv"],
    "report_md": ["REPORT.md"],
    "loss_function": ["loss_function.csv", "L_q_omega.csv"],
}
# dir-style channels (presence of a directory)
DIR_CHANNELS = {
    "eigenvalues": ["eigenvalues"],
    "overlap_wp": ["overlap"],
    "overlap_full": ["overlap_full"],
    "leed_screens": ["screens", "screens_leed_window", "leed"],
    "energy_decomposition": ["energy"],          # spectra/energy
    "density_fourier": ["density_fourier", "fourier"],
}


def _list_vti_frames(d):
    try:
        files = sorted(f for f in os.listdir(d) if f.endswith(".vti"))
    except OSError:
        return []
    return files


def detect_observables(summary_dir):
    """Walk summary_dir (bounded) and detect observable channels by path.
    Returns dict: channel -> {"dir": relpath, "nframes": int, "files": [..]}.
    Paths are relative to summary_dir.
    """
    found = {}
    # Build an index of every dir and file under summary_dir (skip deep junk).
    all_dirs = {}
    all_files = {}
    for dirpath, dirnames, filenames in os.walk(summary_dir):
        # don't descend into per-frame screen snapshots explosion beyond need
        rel = os.path.relpath(dirpath, summary_dir)
        all_dirs.setdefault(os.path.basename(dirpath), []).append(dirpath)
        for fn in filenames:
            all_files.setdefault(fn, []).append(os.path.join(dirpath, fn))

    def relp(p):
        return os.path.relpath(p, summary_dir)

    # VTI channels
    for chan, names in VTI_CHANNELS.items():
        for nm in names:
            if nm in all_dirs:
                d = all_dirs[nm][0]
                frames = _list_vti_frames(d)
                if frames:
                    found[chan] = {
                        "dir": relp(d),
                        "nframes": len(frames),
                        "files": [relp(os.path.join(d, f)) for f in frames],
                    }
                    break
    # gifs channel: any .gif under summary_dir
    gifs = [p for fn, ps in all_files.items() if fn.endswith(".gif") for p in ps]
    if gifs:
        gifs = sorted(gifs)
        found["density_gifs"] = {
            "dir": relp(os.path.dirname(gifs[0])),
            "nframes": len(gifs),
            "files": [relp(p) for p in gifs],
        }
    # FILE channels
    for chan, names in FILE_CHANNELS.items():
        for nm in names:
            if nm in all_files:
                p = all_files[nm][0]
                found[chan] = {"dir": relp(p), "nframes": NULL, "files": [relp(p)]}
                break
    # DIR channels
    for chan, names in DIR_CHANNELS.items():
        for nm in names:
            if nm in all_dirs:
                d = all_dirs[nm][0]
                try:
                    files = sorted(os.listdir(d))
                except OSError:
                    files = []
                found[chan] = {
                    "dir": relp(d),
                    "nframes": NULL,
                    "files": [relp(os.path.join(d, f)) for f in files],
                }
                break
    return found


# --------------------------------------------------------------------------
# Per-run extraction
# --------------------------------------------------------------------------

def pick_alias(d, *keys):
    for k in keys:
        if k in d and str(d[k]).strip() != "":
            return d[k]
    return None


def classify_sim_type(system, d, run_name, proj):
    rt = (pick_alias(d, "run_type") or "").lower()
    rn = run_name.lower()
    proj_str = (pick_alias(d, "projectile") or "").lower()
    # tokenise every path segment (split on _ and -) for keyword detection
    toks = set()
    for seg in re.split(r'[/_\-]', rn + "/" + rt):
        toks.add(seg)
    base = os.path.basename(rn)
    is_gs = (("ground_state_energy_ha" in d) or base == "gs" or
             rn.endswith("/gs") or "/gs/" in rn)
    has_proj = ("projectile" in d or "wp_enabled" in d or
                "wp_state_index" in d or "launch_z" in d or "v0_au" in d or
                "velocity_atu" in d or "proj_r0" in d)
    # vacuum absorber-tuning tests are baselines
    if system == "vacuum":
        return "baseline"
    if is_gs and not has_proj:
        return "baseline"
    is_wp = (proj.get("sim_type") == "wp" or "wave-packet" in rt or
             "wavepacket" in rt or "wavepacket" in proj_str or "wp" in toks or
             "wp_enabled" in d or "wp_state_index" in d)
    if is_wp and ("free" in toks):
        return "free_wp"
    if is_wp:
        return "coronene_wp" if system == "coronene" else "wp"
    if system == "coronene":
        # coronene non-wp with a checkpoint/geometry but no projectile => baseline
        return "coronene_wp" if has_proj else "baseline"
    is_classical = ("classical" in rt or "classical" in proj_str or
                    proj.get("sim_type") == "classical" or "cl" in toks or
                    "classical" in toks)
    if is_classical:
        return "classical"
    # CAP/absorber baselines and ground states with no projectile
    if not has_proj:
        return "baseline"
    return NULL


def geometry_kind_for(system, d):
    bg = (pick_alias(d, "background", "geometry") or "").lower()
    if system == "cylindrical_jellium":
        return "annular_tube"
    if system == "vacuum":
        return "vacuum"
    if system == "coronene":
        return "molecule"
    if system == "graphene":
        return "periodic_slab"
    if system == "localised_jellium":
        if "sphere" in bg:
            return "localised_sphere"
        if "slab" in bg:
            return "localised_slab"
        return "localised_slab"
    if system == "jellium":
        return "cubic_periodic"
    return NULL


def build_row(system, run_dir, summary_path, summary_text, all_summary_paths):
    d, d_raw = parse_summary_text(summary_text, want_raw=True)
    run_name = make_run_name(system, run_dir)
    run_id = system + "/" + run_name
    summary_dir = os.path.dirname(summary_path)

    proj = parse_projectile(d_raw.get("projectile", pick_alias(d, "projectile")))
    geom = parse_geometry(d_raw.get("geometry", pick_alias(d, "geometry")))
    cap = parse_cap(d_raw.get("cap", pick_alias(d, "cap")))

    row = {c: NULL for c in COLUMNS}

    # ---- A Identity / provenance ----
    row["run_id"] = run_id
    row["system"] = system
    row["run_name"] = run_name
    row["run_path"] = run_dir
    row["run_type"] = pick_alias(d, "run_type", "mode") or NULL
    sim_type = classify_sim_type(system, d, run_name, proj)
    row["sim_type"] = sim_type
    row["engine"] = pick_alias(d, "engine") or "inq"
    rc = parse_bool(pick_alias(d, "run_completed"))
    row["run_completed"] = rc if rc is not None else NULL
    row["date_finished"] = pick_alias(d, "date_finished") or NULL
    wt = first_float(pick_alias(d, "wall_time_s", "wall_s"))
    row["wall_time_s"] = wt if wt is not None else NULL

    # ---- run.cpp / config provenance + propagator ----
    run_cpp = find_run_cpp(system, run_dir, sim_type)
    row["run_cpp_path"] = run_cpp or NULL
    cfg_hdr = find_config_header(system, run_cpp)
    row["config_header_path"] = cfg_hdr or NULL
    prop = extract_propagator(run_cpp, d)
    row["propagator"] = prop or NULL

    # caveats
    caveats = []
    for k in ("PROVISIONAL", "ghost_background_term_omitted"):
        if k in d:
            caveats.append(f"{k}={d[k]}")
    row["known_caveats"] = "; ".join(caveats) if caveats else NULL
    row["is_pilot"] = any(t in run_name.lower() for t in
                          ("pilot", "smoke", "dryrun", "validate", "syntax",
                           "buildsmoke", "build_smoke"))

    # ---- rs compound line (jellium sv-sweep): 'rs = 5.69 (N=162, L=50, dx=0.40)'
    rs_raw = d_raw.get("rs") or d_raw.get("r_s")
    rs_L = rs_dx = rs_N = None
    if rs_raw:
        mL = re.search(r'\bL\s*=\s*([-+]?\d*\.?\d+)', rs_raw)
        mdx = re.search(r'\bdx\s*=\s*([-+]?\d*\.?\d+)', rs_raw)
        mN = re.search(r'\bN\s*=\s*([-+]?\d*\.?\d+)', rs_raw)
        rs_L = float(mL.group(1)) if mL else None
        rs_dx = float(mdx.group(1)) if mdx else None
        rs_N = float(mN.group(1)) if mN else None

    # ---- B Geometry / target ----
    cx, cy, cz = parse_cell(d_raw.get("cell_bohr", pick_alias(d, "cell_bohr")))
    if cx is None and rs_L is not None:
        cx = cy = cz = rs_L  # sv-sweep cubic box from rs-line L
    row["geometry_kind"] = geometry_kind_for(system, d)
    row["cell_x"] = cx if cx is not None else NULL
    row["cell_y"] = cy if cy is not None else NULL
    row["cell_z"] = cz if cz is not None else NULL
    row["boundary"] = pick_alias(d, "boundary", "periodicity") or NULL
    sp = first_float(pick_alias(d, "spacing_bohr", "spacing"))
    if sp is None and rs_dx is not None:
        sp = rs_dx
    row["spacing_bohr"] = sp if sp is not None else NULL
    cut = first_float(pick_alias(d, "cutoff_ha"))
    row["cutoff_ha"] = cut if cut is not None else NULL
    nel = first_float(pick_alias(d, "n_electrons", "num_electrons",
                                 "n_electrons_bath", "n_electrons_requested",
                                 "extra_electrons"))
    if nel is None and rs_N is not None:
        nel = rs_N
    row["n_electrons"] = nel if nel is not None else NULL
    nocc = first_float(pick_alias(d, "n_occupied"))
    row["n_occupied"] = nocc if nocc is not None else NULL
    ext = first_float(pick_alias(d, "extra_states"))
    row["extra_states"] = ext if ext is not None else NULL
    row["xc_functional"] = pick_alias(d, "xc_functional", "xc") or NULL
    row["spin"] = pick_alias(d, "spin") or NULL
    # r_s and n0 (explicit, else derived later)
    rs = first_float(pick_alias(d, "r_s", "r_s_eff", "rs"))
    n0 = first_float(pick_alias(d, "n0", "n0_a0m3"))
    row["r_s"] = rs if rs is not None else NULL
    row["n0"] = n0 if n0 is not None else NULL
    row["R_cl"] = first_float(pick_alias(d, "R_cl")) or NULL
    shw = first_float(pick_alias(d, "slab_half_width", "half_width"))
    if shw is None and "background" in d:
        shw = first_float(re.search(r'half_width\s+([-+]?\d*\.?\d+)',
                                    d["background"]).group(1)) \
            if re.search(r'half_width\s+([-+]?\d*\.?\d+)', d["background"]) else None
    row["slab_halfwidth"] = shw if shw is not None else NULL
    row["R_in"] = geom.get("R_in", first_float(pick_alias(d, "R_in"))) or NULL
    row["R_out"] = geom.get("R_out", first_float(pick_alias(d, "R_out"))) or NULL
    Lz = geom.get("L_z", first_float(pick_alias(d, "L_z")))
    row["L_z"] = Lz if Lz is not None else NULL
    row["n_ions"] = first_float(pick_alias(d, "n_ions", "n_carbon",
                                           "num_atoms")) or NULL
    row["geometry_file"] = pick_alias(d, "geometry_file") or NULL
    row["cap_form"] = cap.get("cap_form", pick_alias(d, "cap_form")) or NULL
    ce = cap.get("cap_eta")
    if ce is None:
        ce = first_float(pick_alias(d, "eta_Ha", "cap_eta", "eta"))
    row["cap_eta"] = ce if ce is not None else NULL
    row["cap_width_frac"] = cap.get("cap_width_frac",
                                    first_float(pick_alias(d, "cap_width_frac"))) or NULL
    row["kpoints"] = pick_alias(d, "kpoints") or NULL
    row["smearing_type"] = pick_alias(d, "smearing_type") or NULL
    st = first_float(pick_alias(d, "temperature_ev", "smearing_temp_ev"))
    row["smearing_temp_ev"] = st if st is not None else NULL
    row["net_charge"] = first_float(pick_alias(d, "net_charge")) or NULL
    row["target_pseudopotential"] = pick_alias(d, "pseudopotential", "psp",
                                               "target_pseudopotential") or NULL
    row["pseudo_family"] = pick_alias(d, "pseudo_family") or NULL
    sct = first_float(pick_alias(d, "scf_tol_ha", "scf_energy_tol"))
    row["scf_energy_tol"] = sct if sct is not None else NULL
    row["scf_mixing"] = pick_alias(d, "scf_mixing") or NULL
    row["engine_commit"] = pick_alias(d, "engine_commit") or NULL

    # ---- C Projectile / WP ----
    row["projectile_kind"] = proj.get("projectile_kind",
                                      pick_alias(d, "projectile_kind",
                                                 "projectile")) or NULL
    upf = pick_alias(d, "projectile_upf", "psp", "pseudopotential")
    # only treat as projectile UPF for projectile runs (not the GS target psp)
    if sim_type in ("classical", "wp", "free_wp", "coronene_wp") and upf:
        row["projectile_upf"] = upf
    else:
        row["projectile_upf"] = pick_alias(d, "projectile_upf") or NULL

    # --- classical projectile potential FORM + width, from the UPF V(r) ---
    # Locate the UPF (summary path / config header / run.cpp / scripts scan),
    # then classify by inspecting PP_LOCAL.  gaussian -> sigma_pot from V(0),
    # sigma_WP = sigma_pot*sqrt(2); coulombic (ONCV/point) -> width undefined.
    classical_form = NULL
    classical_sigpot = None  # gaussian charge std (None if coulombic/point)
    if sim_type == "classical":
        cl_upf = locate_classical_upf(system, run_dir, run_cpp, cfg_hdr, d)
        upf_form = upf_sigpot = None
        if cl_upf:
            upf_form, upf_sigpot, _v0 = classify_upf_potential(cl_upf)
            if row["projectile_upf"] == NULL:
                row["projectile_upf"] = cl_upf
        # FORM: prefer UPF V(r); else infer from the summary's projectile line
        # (an explicit 'sigma_pot N' => a Gaussian-charge projectile).
        if upf_form is not None:
            classical_form = upf_form
        elif proj.get("sigma_pot") is not None:
            classical_form = "gaussian"
        # WIDTH (gaussian only): prefer the run's OWN reported sigma_pot (the
        # per-run truth, e.g. localised 'sigma_pot 0.35'); else the UPF-fitted
        # sigma_pot.  coulombic/point => undefined (None).
        if classical_form == "gaussian":
            summary_sp = proj.get("sigma_pot")
            classical_sigpot = summary_sp if summary_sp is not None else upf_sigpot
    row["classical_potential_form"] = classical_form
    row["projectile_charge"] = first_float(pick_alias(d, "projectile_charge")) or NULL
    mass = proj.get("mass_au", first_float(pick_alias(d, "mass_au")))
    row["mass_au"] = mass if mass is not None else NULL

    # launch position
    lx = ly = lz = None
    lb = pick_alias(d, "launch_bohr", "wp_center_bohr", "proj_r0")
    if lb:
        v = all_floats(lb)
        if len(v) >= 3:
            lx, ly, lz = v[0], v[1], v[2]
    if lz is None:
        lz = first_float(pick_alias(d, "launch_z", "wp_offset_bohr"))
    row["launch_x"] = lx if lx is not None else NULL
    row["launch_y"] = ly if ly is not None else NULL
    row["launch_z"] = lz if lz is not None else NULL

    # velocity & energy
    vel = None
    vraw = pick_alias(d, "velocity_atu", "proj_v0")
    if vraw:
        vv = all_floats(vraw)
        if len(vv) >= 3:
            vel = math.sqrt(sum(x * x for x in vv))
        elif vv:
            vel = vv[0]
    if vel is None:
        vel = first_float(pick_alias(d, "v0_au", "v0", "velocity_au"))
    k0 = proj.get("wp_k0", first_float(pick_alias(d, "wp_k0_bohr_inv",
                                                  "wp_k0", "k0")))
    if k0 is not None and re.search(r'\s', str(pick_alias(d, "wp_k0_bohr_inv") or "")):
        kv = all_floats(pick_alias(d, "wp_k0_bohr_inv"))
        if len(kv) >= 3:
            k0 = math.sqrt(sum(x * x for x in kv))
    if vel is None and k0 is not None and sim_type in ("wp", "free_wp", "coronene_wp"):
        vel = k0  # de Broglie group velocity v = k0/m, m_e = 1 a.u.
    en = proj.get("energy_ev", first_float(pick_alias(
        d, "wp_energy_ev", "wp_ekin_ev", "KE_eV", "E_eV", "wp_E_drift_eV",
        "energy_ev", "KE_initial_eV", "projectile_KE_eV", "projectile_ke_ev")))
    # ALWAYS fill velocity<->energy when derivable (m_e = 1 a.u.):
    #   v = sqrt(2 * E_ha),  E_ha = energy_ev / HA_TO_EV  ;  E_ev = 0.5*v^2*HA_TO_EV
    if vel is None and en is not None:
        vel = math.sqrt(2.0 * (en / HA_TO_EV))
    if en is None and vel is not None:
        en = 0.5 * vel * vel * HA_TO_EV
    row["velocity_au"] = vel if vel is not None else NULL
    row["energy_ev"] = en if en is not None else NULL
    row["ion_dynamics"] = proj.get("ion_dynamics",
                                   pick_alias(d, "ion_dynamics")) or NULL

    wpen = parse_bool(pick_alias(d, "wp_enabled"))
    if wpen is None:
        wpen = sim_type in ("wp", "free_wp", "coronene_wp")
    row["wp_enabled"] = wpen
    # sigma_wp (headline).  CLASSICAL: derive ONLY from the UPF V(r) -- gaussian
    # -> sigma_WP = sigma_pot*sqrt(2); coulombic/point -> NULL (undefined).  Never
    # trust a 'sigma'/filename design label for a classical run (that is the
    # charge std sigma_pot and was the graphene sqrt(2) bug).  WP/free_wp: use the
    # explicit/projectile-line sigma (the 'sigma' key IS sigma_WP for a WP).
    sigwp = first_float(pick_alias(d, "wp_sigma_bohr", "sigma_wp", "sigma_WP"))
    if sigwp is None:
        sigwp = proj.get("sigma_wp")
    sigpot = proj.get("sigma_pot")
    if sim_type == "classical":
        sigpot = classical_sigpot  # from UPF V(r); None if coulombic/point
        sigwp = classical_sigpot * SQRT2 if classical_sigpot is not None else None
    else:
        if sigwp is None and sigpot is not None:
            sigwp = sigpot * SQRT2
        # WP 'sigma' key is the projectile cloud sigma_WP (e.g. graphene)
        if sigwp is None:
            sigwp = first_float(pick_alias(d, "sigma"))
    row["sigma_wp_bohr"] = sigwp if sigwp is not None else NULL
    row["wp_k0_bohr_inv"] = k0 if k0 is not None else NULL
    row["wp_direction"] = pick_alias(d, "wp_direction") or NULL
    row["wp_occupation"] = first_float(pick_alias(d, "wp_occupation")) or NULL
    na = first_float(pick_alias(d, "norm_after", "wp_norm_after"))
    row["wp_norm_after"] = na if na is not None else NULL
    mo = first_float(pick_alias(d, "max_overlap", "wp_max_overlap"))
    row["max_overlap"] = mo if mo is not None else NULL
    row["orthogonalised"] = parse_bool(pick_alias(d, "orthogonalised"))
    if row["orthogonalised"] is None:
        row["orthogonalised"] = NULL
    row["passed_tol"] = parse_bool(pick_alias(d, "passed_tol"))
    if row["passed_tol"] is None:
        row["passed_tol"] = NULL

    # ---- D RT / solver ----
    dt = first_float(pick_alias(d, "dt_au", "dt"))
    row["dt_au"] = dt if dt is not None else NULL
    ns = first_float(pick_alias(d, "rt_num_steps", "n_steps", "N_STEPS"))
    row["n_steps"] = ns if ns is not None else NULL
    tt = first_float(pick_alias(d, "total_time_au"))
    if tt is None and dt is not None and ns is not None:
        tt = dt * ns
    row["total_time_au"] = tt if tt is not None else NULL
    we = first_float(pick_alias(d, "write_every", "vti_every"))
    row["write_every"] = we if we is not None else NULL
    sse = first_float(pick_alias(d, "screen_snap_every"))
    row["screen_snap_every"] = sse if sse is not None else NULL

    # ---- E GS source ----
    row["gs_source"] = pick_alias(d, "gs_dir", "checkpoint_dir") or NULL

    # ---- derived n0 / r_s ----
    vol = None
    if all(isinstance(row[k], float) for k in ("cell_x", "cell_y", "cell_z")):
        vol = row["cell_x"] * row["cell_y"] * row["cell_z"]
    if isinstance(row["n0"], str) and isinstance(row["n_electrons"], float) \
            and vol and system in ("jellium",):
        # uniform jellium fills the cell: n = N / V  (A&M, homogeneous gas)
        row["n0"] = row["n_electrons"] / vol
    n_for_phys = row["n0"] if isinstance(row["n0"], float) else None
    if n_for_phys is None and isinstance(row["r_s"], float) and row["r_s"] > 0:
        # n = 3 / (4 pi r_s^3)   (A&M eq. 1.2)
        n_for_phys = 3.0 / (4.0 * math.pi * row["r_s"] ** 3)
    if isinstance(row["r_s"], str) and n_for_phys:
        row["r_s"] = (3.0 / (4.0 * math.pi * n_for_phys)) ** (1.0 / 3.0)

    # ---- Family 2 derived physics ----
    if n_for_phys and n_for_phys > 0:
        kF = (3.0 * math.pi ** 2 * n_for_phys) ** (1.0 / 3.0)  # A&M 2.21
        row["kF"] = kF
        row["E_F_ev"] = (kF ** 2 / 2.0) * HA_TO_EV               # A&M 2.26
        row["v_F"] = kF                                          # A&M 2.24 (a.u.)
        row["omega_p_ev"] = math.sqrt(4.0 * math.pi * n_for_phys) * HA_TO_EV  # A&M 1.38
        if isinstance(row["velocity_au"], float):
            row["v_over_vF"] = row["velocity_au"] / kF
    if isinstance(sigwp, float):
        row["sigma_pot_bohr"] = sigwp / SQRT2  # project convention
    # spreading ratio = 2 sigma^2 / t_transit, t_transit = L_along_v / v0
    if isinstance(sigwp, float) and isinstance(row["velocity_au"], float) \
            and row["velocity_au"] != 0 and isinstance(row["cell_z"], float):
        t_transit = row["cell_z"] / row["velocity_au"]
        if t_transit != 0:
            row["spreading_ratio"] = 2.0 * sigwp ** 2 / t_transit
    # grid dims
    if isinstance(sp, float) and sp > 0:
        for axis in ("x", "y", "z"):
            cval = row["cell_" + axis]
            if isinstance(cval, float):
                row["grid_n" + axis] = int(round(cval / sp))
    if isinstance(dt, float) and isinstance(we, float):
        row["frame_dt_au"] = dt * we

    # run data size
    try:
        out = subprocess.run(["du", "-sb", run_dir], capture_output=True,
                             text=True, timeout=120)
        if out.returncode == 0:
            nbytes = int(out.stdout.split()[0])
            row["run_data_gb"] = round(nbytes / 1e9, 4)
    except (subprocess.SubprocessError, ValueError, OSError):
        row["run_data_gb"] = NULL

    # ---- F observables ----
    obs = detect_observables(summary_dir)
    n_obs = 0
    for chan in OBS_CHANNELS:
        if chan in obs:
            row[chan + "_dir"] = obs[chan]["dir"]
            row[chan + "_nframes"] = obs[chan]["nframes"]
            n_obs += 1
        else:
            row[chan + "_dir"] = NULL
            row[chan + "_nframes"] = NULL
    row["n_observables"] = n_obs

    # relevant_to_induced_density
    has_density = any(c in obs for c in
                      ("density_delta_vti", "density_wp_vti",
                       "density_system_vti", "density_total_vti"))
    row["relevant_to_induced_density"] = bool(
        sim_type in ("classical", "wp", "free_wp", "coronene_wp") and has_density)

    return row, obs, {"sigma_pot": sigpot}


# --------------------------------------------------------------------------
# Column schema / data dictionary
# --------------------------------------------------------------------------

OBS_CHANNELS = [
    "density_total_vti", "density_system_vti", "density_wp_vti",
    "density_delta_vti", "wp_wavefunction_vti", "observables_csv",
    "eigenvalues", "wp_momentum_stats", "wp_realspace_stats", "state_energies",
    "occupations", "momentum_distribution", "gamma_transitions",
    "electron_track", "leed_screens", "overlap_wp", "overlap_full",
    "report_md", "loss_function", "energy_decomposition", "density_fourier",
    "density_gifs",
]

# (column, group, dtype, unit, source, description)
SCHEMA = [
    ("run_id", "A_identity", "str", "", "derived", "system/run_name unique key"),
    ("system", "A_identity", "str", "", "path", "owning system folder"),
    ("run_name", "A_identity", "str", "", "path", "cleaned run path under system"),
    ("run_path", "A_identity", "str", "", "path", "absolute run directory"),
    ("run_type", "A_identity", "str", "", "summary", "raw run_type/mode string"),
    ("sim_type", "A_identity", "str", "", "derived", "classical|wp|free_wp|coronene_wp|baseline"),
    ("engine", "A_identity", "str", "", "summary", "inq|inq-study"),
    ("run_completed", "A_identity", "bool", "", "summary", "run finished flag"),
    ("date_finished", "A_identity", "str", "ISO8601", "summary", "completion timestamp"),
    ("wall_time_s", "A_identity", "float", "s", "summary", "wall-clock runtime"),
    ("relevant_to_induced_density", "A_identity", "bool", "", "derived", "projectile run with delta-density output"),
    ("config_header_path", "A_identity", "str", "", "provenance", "shared/configs header #included by run.cpp"),
    ("run_cpp_path", "A_identity", "str", "", "provenance", "run.cpp that produced the run"),
    ("known_caveats", "A_identity", "str", "", "summary", "PROVISIONAL / omitted-term notes"),
    ("is_pilot", "A_identity", "bool", "", "derived", "name contains pilot/smoke/dryrun/validate/syntax"),
    ("geometry_kind", "B_geometry", "str", "", "derived", "cubic_periodic|localised_slab|localised_sphere|annular_tube|molecule|periodic_slab|vacuum"),
    ("cell_x", "B_geometry", "float", "Bohr", "summary", "cell length x"),
    ("cell_y", "B_geometry", "float", "Bohr", "summary", "cell length y"),
    ("cell_z", "B_geometry", "float", "Bohr", "summary", "cell length z"),
    ("boundary", "B_geometry", "str", "", "summary", "boundary/periodicity"),
    ("spacing_bohr", "B_geometry", "float", "Bohr", "summary", "real-space grid spacing"),
    ("cutoff_ha", "B_geometry", "float", "Ha", "summary", "plane-wave energy cutoff"),
    ("n_electrons", "B_geometry", "float", "", "summary", "electron count (bath/total)"),
    ("n_occupied", "B_geometry", "float", "", "summary", "occupied states"),
    ("extra_states", "B_geometry", "float", "", "summary", "extra/empty states"),
    ("xc_functional", "B_geometry", "str", "", "summary", "XC functional"),
    ("spin", "B_geometry", "str", "", "summary", "spin treatment"),
    ("r_s", "B_geometry", "float", "Bohr", "summary/derived", "Wigner-Seitz radius (A&M 1.2)"),
    ("n0", "B_geometry", "float", "Bohr^-3", "summary/derived", "electron density (N/V or 3/(4 pi r_s^3))"),
    ("R_cl", "B_geometry", "float", "Bohr", "summary", "classical/cluster radius"),
    ("slab_halfwidth", "B_geometry", "float", "Bohr", "summary", "slab half-width"),
    ("R_in", "B_geometry", "float", "Bohr", "summary", "annular inner radius"),
    ("R_out", "B_geometry", "float", "Bohr", "summary", "annular outer radius"),
    ("L_z", "B_geometry", "float", "Bohr", "summary", "tube/slab axial length"),
    ("n_ions", "B_geometry", "float", "", "summary", "number of nuclei"),
    ("geometry_file", "B_geometry", "str", "", "summary", "xyz geometry path"),
    ("cap_form", "B_geometry", "str", "", "summary", "absorber functional form"),
    ("cap_eta", "B_geometry", "float", "Ha", "summary", "CAP strength eta"),
    ("cap_width_frac", "B_geometry", "float", "", "summary", "CAP width fraction"),
    ("kpoints", "B_geometry", "str", "", "summary", "k-point grid"),
    ("smearing_type", "B_geometry", "str", "", "summary", "occupation smearing"),
    ("smearing_temp_ev", "B_geometry", "float", "eV", "summary", "electronic temperature"),
    ("net_charge", "B_geometry", "float", "e", "summary", "net charge"),
    ("target_pseudopotential", "B_geometry", "str", "", "summary", "target/background PSP"),
    ("pseudo_family", "B_geometry", "str", "", "summary", "PSP family"),
    ("scf_energy_tol", "B_geometry", "float", "Ha", "summary", "SCF tolerance"),
    ("scf_mixing", "B_geometry", "str", "", "summary", "SCF mixing scheme"),
    ("engine_commit", "B_geometry", "str", "", "summary", "engine git commit"),
    ("projectile_kind", "C_projectile", "str", "", "summary", "projectile description"),
    ("projectile_upf", "C_projectile", "str", "", "summary/located", "projectile UPF path (resolved from config/run.cpp if absent in summary)"),
    ("classical_potential_form", "C_projectile", "str", "", "computed (UPF V(r))", "gaussian|coulombic by PP_LOCAL erf-fit; NULL for non-classical"),
    ("projectile_charge", "C_projectile", "float", "e", "summary", "projectile charge"),
    ("mass_au", "C_projectile", "float", "m_e", "summary", "projectile mass"),
    ("launch_x", "C_projectile", "float", "Bohr", "summary", "launch x"),
    ("launch_y", "C_projectile", "float", "Bohr", "summary", "launch y"),
    ("launch_z", "C_projectile", "float", "Bohr", "summary", "launch z"),
    ("velocity_au", "C_projectile", "float", "a.u.", "summary/derived", "speed (|v| or k0/m_e)"),
    ("energy_ev", "C_projectile", "float", "eV", "summary", "projectile kinetic energy"),
    ("ion_dynamics", "C_projectile", "str", "", "summary", "ehrenfest/fixed"),
    ("wp_enabled", "C_projectile", "bool", "", "summary/derived", "wavepacket projectile present"),
    ("sigma_wp_bohr", "C_projectile", "float", "Bohr", "summary/derived", "wavepacket sigma_WP (headline)"),
    ("wp_k0_bohr_inv", "C_projectile", "float", "Bohr^-1", "summary", "wavepacket central wavevector"),
    ("wp_direction", "C_projectile", "str", "", "summary", "wp propagation axis"),
    ("wp_occupation", "C_projectile", "float", "", "summary", "wp state occupation"),
    ("wp_norm_after", "C_projectile", "float", "", "summary", "norm after injection"),
    ("max_overlap", "C_projectile", "float", "", "summary", "max overlap with occupied states"),
    ("orthogonalised", "C_projectile", "bool", "", "summary", "wp orthogonalised flag"),
    ("passed_tol", "C_projectile", "bool", "", "summary", "injection tolerance passed"),
    ("propagator", "D_solver", "str", "", "run.cpp/summary", "etrs|crank_nicolson (INQ default etrs)"),
    ("dt_au", "D_solver", "float", "a.u.", "summary", "RT time step"),
    ("n_steps", "D_solver", "float", "", "summary", "RT step count"),
    ("total_time_au", "D_solver", "float", "a.u.", "summary/derived", "total propagation time"),
    ("write_every", "D_solver", "float", "", "summary", "VTI write stride"),
    ("screen_snap_every", "D_solver", "float", "", "summary", "screen snapshot stride"),
    ("gs_source", "E_gs", "str", "", "summary", "ground-state checkpoint/gs_dir"),
    ("kF", "F2_derived", "float", "a.u.", "computed (A&M 2.21)", "Fermi wavevector (3 pi^2 n)^(1/3)"),
    ("E_F_ev", "F2_derived", "float", "eV", "computed (A&M 2.26)", "Fermi energy kF^2/2"),
    ("omega_p_ev", "F2_derived", "float", "eV", "computed (A&M 1.38)", "plasma frequency sqrt(4 pi n)"),
    ("v_F", "F2_derived", "float", "a.u.", "computed (A&M 2.24)", "Fermi velocity = kF"),
    ("v_over_vF", "F2_derived", "float", "", "computed", "velocity / v_F"),
    ("sigma_pot_bohr", "F2_derived", "float", "Bohr", "computed", "classical potential width sigma_WP/sqrt(2)"),
    ("spreading_ratio", "F2_derived", "float", "", "computed", "2 sigma^2 / t_transit"),
    ("grid_nx", "F2_derived", "int", "", "computed", "round(cell_x/spacing)"),
    ("grid_ny", "F2_derived", "int", "", "computed", "round(cell_y/spacing)"),
    ("grid_nz", "F2_derived", "int", "", "computed", "round(cell_z/spacing)"),
    ("frame_dt_au", "F2_derived", "float", "a.u.", "computed", "dt_au * write_every"),
    ("run_data_gb", "F2_derived", "float", "GB", "computed (du)", "run directory size"),
    ("n_observables", "F_observables", "int", "", "computed", "count of present observable channels"),
    ("twin_run_ids", "F4_twin", "str", "", "computed", "semicolon-joined matched classical<->wp twin run_ids (NULL if none), nearest-velocity first"),
    ("twin_count", "F4_twin", "int", "", "computed", "number of matched twins"),
    ("best_twin_id", "F4_twin", "str", "", "computed", "closest twin by |delta_v|"),
    ("match_type", "F4_twin", "str", "", "computed", "match_type of the best twin: point_vs_wp|sigma_matched_gauss|exact"),
    ("pair_width_matched", "F4_twin", "bool", "", "computed", "True iff an EXACT width-matched twin exists (classical sigma_pot == wp sigma_WP/sqrt2 within 1%)"),
]

COLUMNS = [c[0] for c in SCHEMA]
# insert per-channel observable columns just before n_observables
_obs_cols = []
for chan in OBS_CHANNELS:
    _obs_cols.append(chan + "_dir")
    _obs_cols.append(chan + "_nframes")
_idx = COLUMNS.index("n_observables")
COLUMNS = COLUMNS[:_idx] + _obs_cols + COLUMNS[_idx:]


def build_column_dict():
    out = {}
    for col, group, dtype, unit, source, desc in SCHEMA:
        out[col] = {"group": group, "dtype": dtype, "unit": unit,
                    "source": source, "description": desc}
    for chan in OBS_CHANNELS:
        out[chan + "_dir"] = {"group": "F_observables", "dtype": "str",
                              "unit": "", "source": "filesystem",
                              "description": f"{chan} directory/file (relative to results)"}
        out[chan + "_nframes"] = {"group": "F_observables", "dtype": "int",
                                  "unit": "", "source": "filesystem",
                                  "description": f"{chan} frame count (VTI series; NULL for scalar files)"}
    return out


# --------------------------------------------------------------------------
# Twin linkage
# --------------------------------------------------------------------------

def _num(v):
    return v if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def classify_match(a, b, extra):
    """Classify a velocity-matched classical<->wp pair.  Returns one of
    'point_vs_wp' | 'sigma_matched_gauss' | 'exact', or None if the pair fits
    none of the three twin categories (then it is NOT a twin).

    - point_vs_wp        : classical is coulombic/point (sigma_WP undefined),
                           matched to a finite-sigma WP.
    - exact              : both gaussian AND classical sigma_pot == wp sigma_WP/sqrt2
                           within 1% (strict width-matched subset).
    - sigma_matched_gauss: both gaussian AND sigma_WP within ~10% (incl. the
                           documented ~1% legacy 0.35-UPF near-matches).
    """
    cl = a if a["sim_type"] == "classical" else b
    wp = b if a["sim_type"] == "classical" else a
    cl_form = cl.get("classical_potential_form")
    cl_sigwp = _num(cl["sigma_wp_bohr"])
    wp_sigwp = _num(wp["sigma_wp_bohr"])
    # coulombic/point classical (sigma_WP genuinely undefined)
    if cl_form == "coulombic" or cl_sigwp is None:
        return "point_vs_wp"
    # classical is gaussian: need a finite WP width to compare
    if wp_sigwp is None:
        return None
    cl_sigpot = extra.get(cl["run_id"], {}).get("sigma_pot")
    if cl_sigpot is None:
        cl_sigpot = cl_sigwp / SQRT2
    target = wp_sigwp / SQRT2  # WP density std == matched classical charge std
    if abs(cl_sigpot - target) / max(abs(target), 1e-9) <= 0.01:
        return "exact"
    if abs(cl_sigwp - wp_sigwp) / max(abs(wp_sigwp), 1e-9) <= 0.10:
        return "sigma_matched_gauss"
    return None


def link_twins(rows, extra):
    """Symmetric velocity-based twin linker.  Each classical<->wp pair in the
    SAME system+geometry_kind with |delta_v|<=8% AND a valid match_type is a
    twin, linked in BOTH directions.  Stores twin_run_ids/twin_count/
    best_twin_id/match_type (+ _twin_matches for the JSON mirror)."""
    VTOL = 0.08
    for r in rows:
        r["twin_run_ids"] = NULL
        r["twin_count"] = 0
        r["best_twin_id"] = NULL
        r["match_type"] = NULL
        r["pair_width_matched"] = NULL
        r["_twin_matches"] = []

    matches = {r["run_id"]: {} for r in rows}  # run_id -> {other_id: (mtype, dv)}
    n = len(rows)
    for i in range(n):
        r = rows[i]
        if r["sim_type"] not in ("classical", "wp"):
            continue
        rv = _num(r["velocity_au"])
        if rv is None:
            continue
        want = "wp" if r["sim_type"] == "classical" else "classical"
        for j in range(n):
            if i == j:
                continue
            o = rows[j]
            if o["sim_type"] != want:
                continue
            if o["system"] != r["system"] or o["geometry_kind"] != r["geometry_kind"]:
                continue
            ov = _num(o["velocity_au"])
            if ov is None:
                continue
            dv = abs(ov - rv) / max(abs(rv), 1e-9)
            if dv > VTOL:
                continue
            mt = classify_match(r, o, extra)
            if mt is None:
                continue
            # symmetric: record in both directions (same mtype + dv)
            matches[r["run_id"]][o["run_id"]] = (mt, dv)
            matches[o["run_id"]][r["run_id"]] = (mt, dv)

    for r in rows:
        m = matches[r["run_id"]]
        if not m:
            continue
        ordered = sorted(m.items(), key=lambda kv: kv[1][1])  # by delta_v
        ids = [oid for oid, _ in ordered]
        r["twin_run_ids"] = ";".join(ids)
        r["twin_count"] = len(ids)
        r["best_twin_id"] = ids[0]
        r["match_type"] = ordered[0][1][0]
        r["pair_width_matched"] = any(mt == "exact" for _, (mt, _dv) in ordered)
        r["_twin_matches"] = [
            {"id": oid, "match_type": mt, "delta_v": round(dv, 4)}
            for oid, (mt, dv) in ordered]


# --------------------------------------------------------------------------
# Output writers
# --------------------------------------------------------------------------

def to_cell(v):
    if v is None:
        return NULL
    if isinstance(v, bool):
        return "true" if v else "false"
    return v


def write_csv(rows, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([to_cell(r.get(c, NULL)) for c in COLUMNS])


def write_json(rows, obs_map, path):
    objs = []
    for r in rows:
        params = {c: (None if r.get(c, NULL) == NULL else r.get(c))
                  for c in COLUMNS
                  if not (c.endswith("_dir") or c.endswith("_nframes"))}
        observables = {}
        for chan, info in obs_map[r["run_id"]].items():
            observables[chan] = {
                "dir": info["dir"],
                "nframes": (None if info["nframes"] == NULL else info["nframes"]),
                "files": info["files"],
            }
        objs.append({"run_id": r["run_id"], "params": params,
                     "twins": r.get("_twin_matches", []),
                     "observables": observables})
    with open(path, "w") as f:
        json.dump(objs, f, indent=1)


def write_columns(path):
    with open(path, "w") as f:
        json.dump(build_column_dict(), f, indent=2)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def richness(text):
    return len(parse_summary_text(text))


def main():
    runs = {}  # run_dir -> (system, summary_path, text)
    for system in SYSTEMS:
        for sp in find_summaries(system):
            rk = run_key_for(sp)
            try:
                with open(sp, "r", errors="ignore") as f:
                    text = f.read()
            except OSError:
                continue
            if rk not in runs or richness(text) > richness(runs[rk][2]):
                runs[rk] = (system, sp, text)

    rows = []
    obs_map = {}
    extra = {}
    for rk, (system, sp, text) in sorted(runs.items()):
        row, obs, ex = build_row(system, rk, sp, text, None)
        rows.append(row)
        obs_map[row["run_id"]] = obs
        extra[row["run_id"]] = ex

    link_twins(rows, extra)

    rows.sort(key=lambda r: r["run_id"])

    write_csv(rows, os.path.join(OUT_DIR, "run_database.csv"))
    write_json(rows, obs_map, os.path.join(OUT_DIR, "run_database.json"))
    write_columns(os.path.join(OUT_DIR, "run_database_columns.json"))

    # ---- console report ----
    print(f"Unique runs: {len(rows)}")
    from collections import Counter
    bd = Counter((r["system"], r["sim_type"]) for r in rows)
    for (s, t), n in sorted(bd.items()):
        print(f"  {s:20s} {t:12s} {n}")
    print(f"Columns: {len(COLUMNS)}")
    # NULL fraction per column
    nrun = len(rows)
    high_null = []
    for c in COLUMNS:
        nn = sum(1 for r in rows if to_cell(r.get(c, NULL)) == NULL)
        if nn / nrun >= 0.5:
            high_null.append((c, nn / nrun))
    print(f"Columns >=50% NULL: {len(high_null)}")
    for c, frac in high_null:
        print(f"  {c:30s} {frac*100:.0f}%")


if __name__ == "__main__":
    main()
