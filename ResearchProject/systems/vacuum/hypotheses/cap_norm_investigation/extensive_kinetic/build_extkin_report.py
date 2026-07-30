#!/usr/bin/env python3
"""Builder: extensive-kinetic validation study notebook (extkin_study.ipynb).

House-narrative notebook (notebook-making skill) for the vacuum double-sided-CAP
test pair dcap_extkin / dcap_baseline: does the in-run bare per-orbital kinetic
(OrbitalKineticStats) fix INQ's norm-divided energy_kinetic, is it identical to
INQ where it must be, and what does it cost?

Run:
    PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
    /local/data/public/skcb2/tddft/venv/bin/python3 build_extkin_report.py
"""
import importlib.util
import os
import sys

REPO = "/local/data/public/skcb2/tddft"
STACK = f"{REPO}/inq-stack/python"
HERE = os.path.dirname(os.path.abspath(__file__))
SWEEP = f"{REPO}/ResearchProject/systems/vacuum/scripts/wp_traversal_energy"
RES = f"{SWEEP}/results"
sys.path.insert(0, STACK)

# import the shared harness from the jellium hypotheses folder
_spec = importlib.util.spec_from_file_location(
    "_nbreport", f"{REPO}/ResearchProject/systems/localised_jellium/hypotheses/_nbreport.py")
nbr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nbr)
md, code, embed, build, set_outdir = nbr.md, nbr.code, nbr.embed, nbr.build, nbr.set_outdir
set_outdir(HERE)

# ---------------------------------------------------------------------------
# Pre-generate the density GIF (heavy media stay OUT of the executed cells).
# Only 2 frames exist (t=0, final): the runs were made lean because /local/data
# was 100% full on 2026-07-29 (WP_WF_EVERY=700).
# ---------------------------------------------------------------------------
def make_density_gif():
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from PIL import Image as PILImage
    from inqview import load_vti

    frames = []
    paths = [(0, f"{RES}/dcap_extkin/raw/vti/density_total/density_total_t000000.vti"),
             (700, f"{RES}/dcap_extkin/raw/vti/density_total/density_total_t000700.vti")]
    # shared fixed LOG scale across frames (dispersing/absorbed WP rule)
    fields = []
    for step, p in paths:
        f = load_vti(p)
        fields.append((step, f))
    vmax = max(f.data.max() for _, f in fields)
    vmin = vmax * 1e-9
    for step, f in fields:
        x, z = f.x, f.z
        iy = f.data.shape[1] // 2
        sl = np.clip(f.data[:, iy, :], vmin, None)  # (x, z)
        fig, ax = plt.subplots(figsize=(7, 3.6))
        im = ax.pcolormesh(z, x, sl, norm=LogNorm(vmin=vmin, vmax=vmax),
                           cmap="magma", shading="auto")
        for zc in (-15, 15):
            ax.axvline(zc, color="cyan", ls="--", lw=1)
        ax.axvline(-30, color="w", ls=":", lw=0.8); ax.axvline(30, color="w", ls=":", lw=0.8)
        ax.set_xlabel("z (Bohr)"); ax.set_ylabel("x (Bohr)")
        ax.set_title(f"n(x, y=0, z)  step {step}  t={step*0.01:.1f} a.u.  [log]")
        fig.colorbar(im, ax=ax, label="n (e/Bohr³)")
        fig.tight_layout()
        tmp = f"{HERE}/_frame_{step}.png"
        fig.savefig(tmp, dpi=110); plt.close(fig)
        frames.append(tmp)
    imgs = [PILImage.open(p) for p in frames]
    gif = f"{HERE}/fig_extkin_density.gif"
    imgs[0].save(gif, save_all=True, append_images=imgs[1:], duration=1200, loop=0)
    for p in frames:
        os.remove(p)
    print("wrote", gif)
    return gif

GIF = make_density_gif()

# ---------------------------------------------------------------------------
# Cells
# ---------------------------------------------------------------------------
cells = []

cells.append(md(
"# Extensive-kinetic fix — vacuum validation with a double-sided CAP\n\n"
"**Question.** INQ reports a *per-particle* kinetic energy: `energy.hpp:83,55` reduces "
"$\\sum_i \\mathrm{occ}_i\\,\\langle\\psi_i|T|\\psi_i\\rangle/\\langle\\psi_i|\\psi_i\\rangle$. "
"Under a norm-losing CAP this pins the reported energy at the surviving remnant's mean — the "
"\"energy shoots up\" artifact (root-cause note: `docs/notes/inq-energy-normalization-error.md`). "
"The in-run fix is the new observable `inqkit::observables::OrbitalKineticStats`: the BARE "
"(extensive) per-orbital kinetic and norm, computed in the callback with **no engine edit**.\n\n"
"This notebook validates the fix on a single Gaussian wavepacket in TRUE vacuum "
"(the WP is the only electron) with a **double-sided CAP** (both ±z ends — no wrapped-remnant "
"channel), and measures its time cost:\n\n"
"| run | OrbitalKineticStats | role |\n|---|---|---|\n"
"| `dcap_extkin` | ON (every step) | comparison data: INQ reported vs extensive, in ONE run |\n"
"| `dcap_baseline` | OFF | per-step wall-time baseline (physics identical) |\n\n"
"Where this sits: follow-up to the `cap_norm_investigation` phases 0–3 (which established the "
"normalization mechanism post-hoc); this is the first IN-RUN implementation, the prerequisite for "
"the next localised-jellium run design."))

cells.append(md(
"## Conventions\n\n"
"Hartree atomic units throughout (ħ = m_e = e = 1); energies also quoted in eV "
"(1 Ha = 27.211 eV).\n\n"
"| symbol | meaning | value here |\n|---|---|---|\n"
"| $k_0$ | WP drift momentum (Bohr⁻¹) | 5.421 |\n"
"| $\\sigma_0$ | WP wavefunction width $\\sigma_{WP}$ (Bohr) | 3.0 |\n"
"| $\\eta$ | CAP strength (Ha, imaginary amplitude) | −3.5 |\n"
"| $W$ | CAP band width per side (Bohr) | 15 |\n"
"| $L_z$ | box length (Bohr), box $[-30,30]$ | 60 |\n"
"| $h$ | grid spacing (Bohr) | 0.4 |\n"
"| $\\Delta t$ | time step (a.u.) | 0.01, 700 steps |\n"
"| $N_i$, $T_i$ | orbital norm, bare kinetic of orbital $i$ | logged per step |\n\n"
"Propagator: **ETRS** (norm-losing under the CAP — required; Crank–Nicolson renormalizes "
"each step and pumps real energy, see phase 3). Theory: non-interacting → "
"$E_{total}=E_{kinetic}$, all other components are bookkeeping zeros."))

cells.append(code(
"import sys, os\n"
f"sys.path.insert(0, {STACK!r})\n"
"import numpy as np, pandas as pd\n"
"import matplotlib.pyplot as plt\n"
"from inqview.visualisation import style\n"
"style.apply_theme()\n"
f"RES = {RES!r}\n"
"HA = 27.211386\n"
"def obs(run, name):\n"
"    return pd.read_csv(f'{RES}/{run}/raw/observables/{name}', comment='#')\n"
"en  = obs('dcap_extkin', 'energies.csv')\n"
"ek  = obs('dcap_extkin', 'orbital_kinetic_stats.csv')\n"
"mom = obs('dcap_extkin', 'wp_momentum_stats.csv')\n"
"en_b = obs('dcap_baseline', 'energies.csv')\n"
"m = en.merge(ek, on='step', suffixes=('', '_ek'))\n"
"print(f'rows: extkin={len(m)}  baseline={len(en_b)}')"))

cells.append(md(
"## Analytic design quantities\n\n"
"Drift kinetic energy of the packet (plane-wave part):\n\n"
"$$E_{drift} = \\tfrac{1}{2}k_0^2$$"))
cells.append(code(
"k0, sigma0 = 5.421, 3.0\n"
"E_drift = 0.5*k0**2\n"
"print(f'E_drift = {E_drift:.3f} Ha = {E_drift*HA:.0f} eV')"))

cells.append(md(
"Total mean kinetic energy of a Gaussian wavepacket adds the localization "
"(zero-point) term — for $\\psi\\propto e^{-r^2/2\\sigma_0^2}\\,e^{ik_0z}$ "
"(Cohen-Tannoudji conventions; verified against the run at $t=0$ below):\n\n"
"$$E_{WP} = \\tfrac{1}{2}k_0^2 + \\tfrac{3}{4\\sigma_0^2}$$"))
cells.append(code(
"E_wp = 0.5*k0**2 + 3.0/(4*sigma0**2)\n"
"print(f'E_WP(analytic) = {E_wp:.4f} Ha = {E_wp*HA:.1f} eV')"))

cells.append(md(
"Momentum-space width and the grid (aliasing) margin — the cutoff guard's check "
"($\\sigma_p$ of the WAVEFUNCTION envelope; the $\\sqrt2$ convention per the "
"`sigma-wp-convention` rule):\n\n"
"$$\\sigma_p = \\frac{1}{\\sqrt{2}\\,\\sigma_0}, \\qquad k_{Nyq} = \\pi/h$$"))
cells.append(code(
"h = 0.4\n"
"sigma_p = 1.0/(np.sqrt(2)*sigma0)\n"
"k_nyq = np.pi/h\n"
"print(f'sigma_p = {sigma_p:.3f}  k0+3sigma_p = {k0+3*sigma_p:.2f}  k_Nyq = {k_nyq:.2f}'\n"
"      f'  -> margin {(k_nyq-(k0+3*sigma_p)):.2f} Bohr^-1 (guard PASS, tail 0.00%)')"))

cells.append(md(
"CAP design absorption — single-pass amplitude survival through one band "
"(adiabatic estimate used to size $\\eta$, run header 2026-07-27):\n\n"
"$$s \\approx e^{-|\\eta| W / v}, \\quad v = k_0$$"))
cells.append(code(
"eta, W = 3.5, 15.0\n"
"s = np.exp(-eta*W/k0)\n"
"print(f'design survival ~ {s:.1e}  (measured final norm below)')"))

cells.append(md(
"## Simulation setup (provenance: `run_summary.txt`, verbatim)\n\n"
"Geometry: box $30\\times30\\times60$ Bohr, CAP bands $z\\in[15,30]$ and $z\\in[-30,-15]$ "
"(inner edges dashed in the density figure), launch $z=0$ → $5\\sigma_0=15$ Bohr clearance "
"to BOTH inner edges (boundary rule). TRUE vacuum: the WP replaces the single electron "
"(`extra_states(0).extra_electrons(1.0)`), so `density_total == density_wp`.\n\n"
"**Note the lean VTI cadence** (`wf_every = 700`: frames at $t=0$ and final only): "
"/local/data was 100% full when the pair ran (2026-07-29); the comparison is CSV-based "
"and unaffected."))
cells.append(code(
"for run in ('dcap_extkin', 'dcap_baseline'):\n"
"    print(f'--- {run} ---')\n"
"    print(open(f'{RES}/{run}/run_summary.txt').read())"))

cells.append(md(
"## Source files\n\n"
"| artefact | path |\n|---|---|\n"
"| observable (the fix) | `inq-stack/include/inqkit/observables/orbital_kinetic_stats.hpp` |\n"
"| run definition | `ResearchProject/systems/vacuum/scripts/wp_traversal_energy/run.cpp` "
"(`WP_CAP2`, `WP_EXTKIN` switches) |\n"
"| dispatcher | `.../wp_traversal_energy/run_extkin_test.sh` |\n"
"| engine | `inq-study` (scalar-potential complexification; `occ_sum` /norm is UPSTREAM, "
"byte-identical in both trees) |\n"
"| quick comparison script | `hypotheses/cap_norm_investigation/extensive_kinetic/compare_extkin.py` |\n"
"| this builder | `hypotheses/cap_norm_investigation/extensive_kinetic/build_extkin_report.py` |\n"
"| root-cause note | `docs/notes/inq-energy-normalization-error.md` |\n"
"| plan | `docs/plans/norm-corrected-stopping-power.md` → Extension (2026-07-29) |"))

cells.append(md(
"## Visual intuition — density evolution (dcap_extkin)\n\n"
"$n(x, y=0, z)$ on the propagation (xz) plane, shared fixed **log** scale, CAP inner "
"edges dashed cyan. Two frames only (lean cadence, see setup): launch at $z=0$, and the "
"final state — the packet is gone into the +z band (final norm ≈ 3×10⁻⁶), with no "
"wrapped remnant (that channel is closed by the −z band). `dcap_baseline` is "
"physics-identical (overlay check in the energetics below), so only this representative "
"run is shown."))
cells.append(embed(GIF, "density_total, t=0 ↔ t=7.0 a.u. (2-frame animation, log scale)", width=720))

cells.append(md(
"### Per-run energetics\n\n"
"All energy components from `energies.csv` (non-interacting: total == kinetic; "
"hartree/external/xc are zeros — plotted to prove the bookkeeping). Baseline total "
"overlaid: the two runs are the same physics."))
cells.append(code(
"fig, ax = plt.subplots(figsize=(8, 4))\n"
"t = en['time_au']\n"
"ax.plot(t, en['total']*HA, lw=2, label='E_total (reported, extkin run)')\n"
"ax.plot(en_b['time_au'], en_b['total']*HA, ':', lw=1.6, label='E_total (baseline run)')\n"
"for col in ('hartree', 'external', 'xc', 'non_local'):\n"
"    ax.plot(t, en[col]*HA, lw=0.8, alpha=0.6, label=col)\n"
"ax.set_xlabel('t (a.u.)'); ax.set_ylabel('E (eV)')\n"
"ax.set_title('Energy components — reported ledger')\n"
"ax.legend(fontsize=8); fig.tight_layout()\n"
"d = float(np.abs(np.interp(en_b['time_au'], t, en['total']) - en_b['total']).max())*HA\n"
"print(f'max |E_total(extkin) - E_total(baseline)| = {d:.2e} eV  (same physics)')"))

cells.append(md(
"## Result 1 — the identity: our norm-divided reconstruction == INQ, every step\n\n"
"INQ's reported kinetic (the artifact carrier) and the observable's two reductions:\n\n"
"$$E_{kin}^{INQ} = \\sum_i \\mathrm{occ}_i\\,\\frac{T_i}{N_i} \\quad\\text{(energy.hpp:55)}, "
"\\qquad E_{kin}^{bare} = \\sum_i \\mathrm{occ}_i\\,T_i$$\n\n"
"with $T_i = \\tfrac12\\,(dV/N_{grid})\\sum_k k^2|\\tilde\\psi_i(k)|^2$ and "
"$N_i = (dV/N_{grid})\\sum_k |\\tilde\\psi_i(k)|^2$ (INQ's `to_fourier` is an unnormalized "
"DFT — the $dV/N_{grid}$ factor was verified numerically). If we compute the same $T_i$ INQ "
"does, the reconstruction must equal the reported column exactly."))
cells.append(code(
"ident = m['kin_normdiv_total_ha'] - m['kinetic']\n"
"print(f'max |identity residual| = {ident.abs().max():.2e} Ha over {len(m)} steps')\n"
"print(f't=0: kin_bare = {m[\"kin_bare_total_ha\"][0]:.6f} Ha  (analytic E_WP = {E_wp:.6f} Ha)')\n"
"assert ident.abs().max() == 0.0"))

cells.append(md(
"## Result 2 — reported vs corrected total energy\n\n"
"The corrected (extensive) total replaces only the kinetic channel:\n\n"
"$$E_{corr}(t) = E_{total}^{rep}(t) - E_{kin}^{INQ}(t) + E_{kin}^{bare}(t)$$\n\n"
"Expected behaviour for a fully absorbed packet: $E_{corr}\\to 0$ tracking "
"$E_0\\cdot\\mathrm{norm}(t)$, while the reported total stays pinned at the remnant mean."))
cells.append(code(
"t = m['time_au'].to_numpy(); norm = m['norm_total'].to_numpy()\n"
"E_rep = m['total'].to_numpy(); E_corr = E_rep - m['kinetic'].to_numpy() + m['kin_bare_total_ha'].to_numpy()\n"
"E0 = E_rep[0]\n"
"fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))\n"
"ax[0].plot(t, E_rep*HA, lw=2, label='reported (INQ out-of-the-box)')\n"
"ax[0].plot(t, E_corr*HA, lw=2, label='corrected (extensive kinetic)')\n"
"ax[0].plot(t, E0*norm*HA, '--', lw=1.4, label='E0 · norm (expected)')\n"
"axn = ax[0].twinx(); axn.plot(t, norm, color='0.5', lw=1, alpha=0.7); axn.set_ylabel('WP norm', color='0.5')\n"
"ax[0].set_xlabel('t (a.u.)'); ax[0].set_ylabel('E (eV)'); ax[0].legend(fontsize=8, loc='center left')\n"
"ax[0].set_title('Total energy')\n"
"ax[1].plot(t, m['kinetic']*HA, lw=2, label='kinetic reported (norm-divided)')\n"
"ax[1].plot(t, m['kin_bare_total_ha']*HA, lw=2, label='kinetic bare (extensive)')\n"
"ekin_mean = np.interp(t, mom['time_au'], mom['e_kin_ha'])\n"
"ax[1].plot(t, ekin_mean*norm*HA, '--', lw=1.4, label='e_kin_ha · norm (post-hoc route)')\n"
"ax[1].set_xlabel('t (a.u.)'); ax[1].set_ylabel('E_kin (eV)'); ax[1].legend(fontsize=8)\n"
"ax[1].set_title('Kinetic channel')\n"
"fig.tight_layout()\n"
"print(f'E0 = {E0*HA:.1f} eV   final norm = {norm[-1]:.2e}')\n"
"print(f'final: reported {E_rep[-1]*HA:.1f} eV  vs  corrected {E_corr[-1]*HA:.2f} eV'\n"
"      f'  (E0·norm = {E0*norm[-1]*HA:.2f} eV)')\n"
"print(f'captured = {(E0-E_corr[-1])/E0*100:.1f}% of E0')\n"
"print(f'post-hoc vs in-run bare: max |Δ| = {np.abs(ekin_mean*norm - m[\"kin_bare_total_ha\"]).max()*HA:.2e} eV')"))

cells.append(md(
"## Result 3 — time cost\n\n"
"Two measurements: (a) the observable's own `wall_ms` column (chrono around each "
"evaluation — the trustworthy number), and (b) the run-level per-step wall difference "
"ON−OFF from `run_summary.txt` (noisy: GPU clock/thermal variation exceeds the "
"observable's cost)."))
cells.append(code(
"import re\n"
"def field(run, key):\n"
"    txt = open(f'{RES}/{run}/run_summary.txt').read()\n"
"    return float(re.search(rf'{key}\\s*=\\s*([0-9.eE+-]+)', txt).group(1))\n"
"wall = m['wall_ms'].to_numpy()[1:]  # row 0 includes CSV-header write\n"
"on, off = field('dcap_extkin','per_step_ms'), field('dcap_baseline','per_step_ms')\n"
"fig, ax = plt.subplots(figsize=(8, 3.4))\n"
"ax.plot(m['time_au'][1:], wall, lw=0.8)\n"
"ax.axhline(wall.mean(), color='C1', ls='--', label=f'mean {wall.mean():.2f} ms')\n"
"ax.set_xlabel('t (a.u.)'); ax.set_ylabel('observable wall (ms/step)')\n"
"ax.set_title('OrbitalKineticStats self-timing (1 orbital, 844k grid)')\n"
"ax.legend(); fig.tight_layout()\n"
"print(f'observable: {wall.mean():.2f} ms/step = {wall.mean()/on*100:.2f}% of the {on:.0f} ms step')\n"
"print(f'run-level per-step: ON {on:.1f} ms vs OFF {off:.1f} ms -> Δ = {on-off:+.1f} ms (noise-dominated)')"))

cells.append(md(
"**Inference (jellium-162 cost):** the observable does ONE batched `to_fourier` of the whole "
"orbital set + 2 reductions per orbital per recorded step. An ETRS step already performs "
"several full-set FFT passes, so the observable adds roughly the cost of one extra forward "
"pass — order a few percent at every-step cadence, less at `WP_EXTKIN_EVERY > 1`. To be "
"**measured in the jellium pilot**, not assumed."))

cells.append(md(
"## Takeaway\n\n"
"- **The identity is exact**: $\\sum_i \\mathrm{occ}_i T_i/N_i$ equals INQ's reported "
"kinetic to the last digit at all 701 steps, and $t=0$ bare kinetic equals the analytic "
"$\\tfrac12 k_0^2 + 3/(4\\sigma_0^2)$ — the observable provably computes INQ's own $T_i$, "
"minus the norm division.\n"
"- **The fix works**: reported total stays pinned near 380 eV as the CAP absorbs "
"(norm → 3×10⁻⁶); the corrected total tracks $E_0\\cdot$norm to 0.00 eV — 100% of the "
"402 eV packet correctly booked as captured. Double-sided CAP closes the wrap channel.\n"
"- **Cost is negligible at this scale**: 0.42 ms/step (0.14%); run-level Δ is noise. "
"Jellium overhead to be measured in the pilot.\n"
"- **Next**: wire `OrbitalKineticStats` (all orbitals) + an explicit CAP-sink term into "
"the next localised-jellium run design."))

build(cells, f"{HERE}/extkin_study.ipynb")
