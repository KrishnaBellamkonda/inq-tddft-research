# Runbook — High-L_z slab stopping-power runs on a bigger server

**Audience:** a fresh Claude instance (or operator) on a **different, larger-memory
machine**, tasked with extending the `slab_sE_sweep` S(E) campaign with **high slab
thickness (L_z) points**. Self-contained: everything needed to reproduce the method
is here. Origin machine = `lovelace` (this repo), origin repo root
`/local/data/public/skcb2/tddft` (call it `$ROOT`).

> Written 2026-10-19 from the `slab_sE_sweep` machinery + `slab-stopping-power`
> skill as they stand on the origin machine. Nothing here is invented — every
> constant is copied from `orchestrate_sE.py` / `launch_2gpu.sh` / the `run.cpp`
> env contracts.

### Resolved parameters (confirmed by the project owner, 2026-10-19)

- **Server is a FULL working-tree copy** of the origin repo — `inq-study/`, the
  pseudopotential share, `inqkit` headers, the venv, and the `shared_gs/slab_sE_sweep`
  ground states are already present. `git pull origin quantum-stopping-power` only
  refreshes the scripts/CSVs; no heavy rsync is needed. **Verify presence** (§4) and
  rebuild the binaries for the server's GPU arch (§5) — then run.
- **Scope: L_z = 60, 75, 90; CLASSICAL HALF FIRST** (all 3 thicknesses × 5 energies
  = 15 classical runs), validate the classical S(L_z) line, THEN run the WP half
  (another 15). Do **not** re-run L_z = 15/30/45.
- **Scheduler: Slurm** — launch each run as an `sbatch` job (§6 Slurm path). The
  setsid/plain-GPU recipe is retained only as a fallback.

---

## 1. TL;DR — what to do and why

The S(E) result is built by a **finite-size extrapolation**: at each projectile
energy we simulate a jellium slab of thickness L_z, measure the electronic
stopping power `S(L_z) = E_absorbed / L_z`, and extrapolate `1/L_z → 0` (quadratic
fit) to the bulk limit. We do this for a **classical** point-charge projectile and
a **quantum wavepacket (WP)** projectile; the gap between the two bulk limits is
the *quantum correction*.

**The limitation we are fixing:** each bulk point currently rests on only **3 slab
lengths, L_z = 15, 30, 45 Bohr**. We wanted L_z = 60, 75 too, but a 24 GB GPU
UVM-thrashes (~100× slowdown) on those boxes, so they were dropped
(`orchestrate_sE.py:53`). **This is a pure memory wall, not physics.** A
bigger-memory GPU removes it.

**Your job:** run the **high-L_z points** — **L_z = 60, 75, 90** — for all 5
energies, **classical half first (15 runs), then WP (15 runs)**, measure `S(L_z)`
for each with the locked method, and report the new rows so they can be merged into
`sweep_results.csv` and the `1/L_z` extrapolation re-fit on more points.

Nothing about geometry, engine, projectile, CAP, or analysis changes — **only L_z
grows**. Do **not** re-run L_z = 15/30/45 (that data already exists and is good).

---

## 2. The physical setup (locked — do not change)

Localised-jellium **slab** centred at z=0, projectile launched from the vacuum on
one side, traverses the slab, exits into a complex absorbing potential (CAP) on the
far side. r_s ≈ 4.18 jellium (n0 = 0.001288 a0⁻³). Two projectile representations,
**identical in every parameter except the representation**:

| Quantity | Value | Notes |
|---|---|---|
| Transverse box L_x = L_y | 50 Bohr | periodic in x,y (periodicity = 2, open z) |
| Box length along z | **L_z + 60** Bohr | 30 Bohr vacuum/CAP head-room each side |
| Slab half-width | **L_z / 2** | slab fills the central L_z of the box |
| Electron count N | **round(0.001288 · 50 · 50 · L_z)** | keeps n0 fixed as L_z grows |
| Jellium n0 | 0.001288 a0⁻³ (r_s ≈ 4.18) | LDA XC |
| Grid spacing dx | 0.5 Bohr | |
| Edge width | 1 Bohr | soft slab edge |
| σ_WP (wavepacket width) | **2.0 Bohr** | σ is ALWAYS the WP width (see note) |
| Launch z | **−L_z/2 − 10** Bohr | 10 Bohr run-up from the near face |
| CAP | two-sided sin², η = −1.0 Ha, width 10 Bohr/face | absorbs the transmitted projectile |
| Time step dt | 0.04 a.u. | ETRS propagator (CAP needs ETRS, not CN) |
| Classical mass / FD δ | 1.0 mₑ / 0.1 | moving Gaussian charge, velocity-Verlet |

**σ note:** reported σ is the **wavepacket** width σ_WP. The classical projectile's
Gaussian *potential* uses σ_pot = σ_WP/√2 internally — never relabel a curve by
σ_pot. (See `.claude/rules/sigma-wp-convention.md` on origin.)

### Energy → initial velocity map (KE = ½ m v₀²)

| E (eV) | v₀ (a.u.) |
|---|---|
| 31 | 1.5 |
| 54 | 2.0 |
| 100 | 2.71 |
| 167 | 3.5 |
| 276 | 4.5 |

### Derived quantities for the high-L_z boxes (verified)

| L_z | N | box L_z | half | launch z | nsteps classical | nsteps WP (base) |
|---|---|---|---|---|---|---|
| 60 | 193 | 120 | 30.0 | −40.0 | `ceil(1.6·(60+45)/v₀/0.04)` | `ceil((1.6·(60+45)/v₀+60)/0.04)` |
| 75 | 241 | 135 | 37.5 | −47.5 | `ceil(1.6·(75+45)/v₀/0.04)` | `ceil((1.6·(75+45)/v₀+60)/0.04)` |
| 90 | 290 | 150 | 45.0 | −55.0 | `ceil(1.6·(90+45)/v₀/0.04)` | `ceil((1.6·(90+45)/v₀+60)/0.04)` |

Step-count rules (from `orchestrate_sE.py`, SIZE_FAC = 1.6, WP_DRAIN = 60 a.u.):
- classical: `n_steps = ceil( 1.6·(L_z + 45)/v₀ / 0.04 )`
- WP: `n_steps = ceil( (1.6·(L_z + 45)/v₀ + 60) / 0.04 )`  — a drain margin so the
  dispersing packet fully CAP-absorbs. **Auto-extend** if it doesn't (see §6).

---

## 3. How S is measured (the method that makes the numbers trustworthy)

For **one completed run** → one `S(L_z)` value, via the `slab-stopping-power`
skill (`$ROOT/.claude/skills/slab-stopping-power/slab_stopping.py`,
`measure(obs_dir, half, e_gs_ha, L_z)`):

```
S = E_absorbed / L_z           (L_z = slab thickness, NOT the box length)
```

`E_absorbed` is anchored to the **projectile-free ground-state energy E_GS** (the
Ha number saved in the GS run_summary), never to E_total(t=0):

- **classical:** `E_absorbed = [ E_total_INQ(t_f) + E_PB(t_f) ] − E_GS`.
  INQ omits E_PB (projectile↔background); it must be added back from
  `interactions.csv` (`e_pb` column). A closure gate checks
  `energy_external == E_SB + E_PS` before trusting it.
- **wavepacket:** `E_absorbed = [ E_total_INQ(t_f) − T1·(1−norm) ] − E_GS`.
  E_PB is already inside energy_total (the WP is a real KS orbital). The norm term
  removes INQ's per-particle-mean inflation as the CAP drains the packet. **Valid
  ONLY if the WP is fully absorbed** — gate on `norm_wp_final ≤ 1e-3`.

The bulk limit is then `S_bulk = ` the `1/L_z = 0` intercept of a quadratic fit
`S = a + b/L_z + c/L_z²` through all L_z at that energy (`build_S_vs_invLz.py`).

---

## 4. Prerequisites — verify the full working-tree copy

The server **is a full working-tree copy** (owner-confirmed), so the heavy artefacts
are already present. Refresh the scripts and verify, don't rsync:

```bash
cd $ROOT && git fetch origin && git checkout quantum-stopping-power \
  && git pull origin quantum-stopping-power     # brings run.cpp, orchestrator, fit scripts, CSVs, this runbook

# verify the heavy artefacts are present (all must exist):
ls -d $ROOT/inq-study/ && ls $INQ_SHARE_PATH/pseudopod >/dev/null && echo "engine+share OK"
ls -d $ROOT/inq-stack/include/inqkit/ && echo "inqkit headers OK"
$PY -c "import inqview; print('inqview OK')"
ls -d $LJ/shared_gs/slab_sE_sweep/Lz60/checkpoint $LJ/shared_gs/slab_sE_sweep/Lz75/checkpoint \
  && echo "GS Lz60/Lz75 present (Lz90 to be generated)"
```

> **Path gotcha:** `orchestrate_sE.py` and `build_S_vs_invLz.py` import
> `slab_stopping` from a **hardcoded** path
> `/home/raid/skcb2/skcb2/tddft/.claude/skills/slab-stopping-power`. If the server's
> repo root differs, edit those two `sys.path.insert(...)` / `SKILL =` lines to point
> at `$ROOT/.claude/skills/slab-stopping-power/`.

Reference of what each piece is (should all already be on the full copy):

| Path (under `$ROOT`) | Role | Must be present |
|---|---|---|
| `inq/` | upstream INQ (header-only) + `install/share` pseudopotentials | yes |
| `inq-study/` | INQ replica that compiles the CAP absorber (stock `inq/` CANNOT) | **yes — the engine for every run** |
| `inq-stack/include/inqkit/` | our C++ headers (background, wavepacket, interaction energies, CAP, orbital stats) | yes |
| `inq-stack/python/inqview/` + the venv | post-processing (`pip install -e inq-stack/` into a venv) | yes |
| `ResearchProject/systems/localised_jellium/scripts/slab_sE_sweep/` | `classical/run.cpp`, `wp/run.cpp`, `orchestrate_sE.py`, `launch_2gpu.sh`, `gpu_rank_wrapper.sh`, `sweep_results.csv`, `sE_results.csv` | yes |
| `ResearchProject/systems/localised_jellium/scripts/campaign_autorun/gs/run.cpp` | the env-parametrised GS binary source | yes |
| `ResearchProject/systems/localised_jellium/hypotheses/slab_sE_sweep/` | `build_S_vs_invLz.py`, `build_S_vs_E.py`, `build_S_vs_Lz.py` | for re-fitting |
| `.claude/skills/slab-stopping-power/` | `slab_stopping.py` (the S measurement) | yes |
| `ResearchProject/systems/localised_jellium/shared_gs/slab_sE_sweep/{Lz60,Lz75}/` | **single-partition** GS for L_z 60/75 already exist (checkpoint + gs.log) | reuse if transferred; else regenerate (§5) |
| `ResearchProject/systems/.../pseudopotentials/` or `inq/install/share/pseudopod` | pseudopotential share | yes |

**Carry the existing `sweep_results.csv`** (it holds the L_z = 15/30/45 rows) so the
re-fit in §7 has the low-L_z points to combine with your new high-L_z ones.

Set at the top of your session (adjust to the target's actual repo root):

```bash
export ROOT=/path/to/tddft            # the repo root on THIS machine
export INQ_SHARE_PATH=$ROOT/inq/install/share
export PSEUDOPOD_SHARE_PATH=$ROOT/inq/install/share/pseudopod
export INQ_SOURCE=$ROOT/inq-study     # CAP absorber only compiles here
LJ=$ROOT/ResearchProject/systems/localised_jellium
BASE=$LJ/scripts/slab_sE_sweep
PY=$ROOT/venv/bin/python3             # the venv python (has inqview + VTK)
```

---

## 5. Build the engine + binaries

The CAP absorber does `double += complex` which **stock `inq/` refuses to
compile** — you MUST build against `inq-study`. The project wrapper `inq-run`
auto-detects the `.cpp`, builds, and runs; set `INQ_SOURCE=inq-study` first.

```bash
cd $BASE/classical && INQ_SOURCE=$ROOT/inq-study inq-run --build-only   # produces ./run
cd $BASE/wp        && INQ_SOURCE=$ROOT/inq-study inq-run --build-only   # produces ./run
cd $LJ/scripts/campaign_autorun/gs && INQ_SOURCE=$ROOT/inq-study inq-run --build-only  # GS binary
```

(If `inq-run` lacks `--build-only` on the target, run it once on a 1-step smoke;
it builds then runs. Check `build_*.log` / that `./run` appears.) The binaries must
match the **target GPU architecture** — rebuild on the target, do **not** copy the
origin's `run` ELF. Expect a long first compile (INQ is heavy templated C++17).

**Smoke it** (5 steps) on a small box first to confirm the engine links and energies
are finite before committing GPU-hours — see §8 env block with `LJ_N_STEPS=5`.

---

## 6. Run the high-L_z grid

### Launch model — one Slurm job per run, single GPU

The server uses **Slurm**, and each run goes on **one GPU** (the single-GPU design
sidesteps the 2-GPU GS-partition hazard in the fallback below). Order is
**classical first (Phase 1, 15 jobs) → WP (Phase 2, 15 jobs)**.

**Ground states first.** Lz60/Lz75 single-partition GS are already on the copy —
verify (`r_s ≈ 4.18`, a `checkpoint/` dir). Generate **only Lz90**, as a short GPU
job (`srun` interactively, or an sbatch wrapper of this block):

```bash
LZ=90; N=$($PY -c "print(round(0.001288*2500*$LZ))")
OUT=$LJ/shared_gs/slab_sE_sweep/Lz$LZ; mkdir -p $OUT && cd $OUT
LJ_LX=50 LJ_LY=50 LJ_LZ=$((LZ+60)) LJ_HALF=$($PY -c "print($LZ/2)") \
LJ_N=$N LJ_EDGE_W=1 LJ_PERIODICITY=2 LJ_SPACING=0.5 \
LJ_GS_DIR=$OUT/checkpoint LJ_TAG=Lz$LZ \
INQ_SHARE_PATH=$INQ_SHARE_PATH PSEUDOPOD_SHARE_PATH=$PSEUDOPOD_SHARE_PATH \
srun --gres=gpu:1 $LJ/scripts/campaign_autorun/gs/run 2>&1 | tee gs.log
# PASS: "GS energy = …" printed; run_summary.txt has r_s ≈ 4.18 and a checkpoint/ dir.
```

**Reusable run script** — save as `$BASE/slurm_run.sh`. Fill the cluster-specific
resource lines (partition/account/QoS/time) from the server's own conventions — I
cannot know those; everything below the `set -euo pipefail` is verified:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=sE_slab
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
### >>> FILL FOR THIS CLUSTER: --partition=... --account=... --qos=... <<<
set -euo pipefail
ROOT=/path/to/tddft                       # <-- server repo root
LJ=$ROOT/ResearchProject/systems/localised_jellium
BASE=$LJ/scripts/slab_sE_sweep
PY=$ROOT/venv/bin/python3
export INQ_SHARE_PATH=$ROOT/inq/install/share
export PSEUDOPOD_SHARE_PATH=$ROOT/inq/install/share/pseudopod

EV=$1; V0=$2; LZ=$3; HALF=$4              # e.g.  sbatch slurm_run.sh 276 4.5 60 classical
N=$($PY -c "print(round(0.001288*2500*$LZ))")
if [ "$HALF" = classical ]; then
  NS=$($PY -c "import math;print(math.ceil(1.6*($LZ+45)/$V0/0.04))")
else                                      # WP: + 60 a.u. CAP-drain margin
  NS=$($PY -c "import math;print(math.ceil((1.6*($LZ+45)/$V0+60)/0.04))")
fi
RUNDIR=$BASE/runs/E${EV}_Lz${LZ}; OUT=$RUNDIR/results/$HALF
mkdir -p "$RUNDIR"; cd "$RUNDIR"
RES=0                                      # resume if a checkpoint is short of target
if [ -f "$OUT/rt_state.txt" ]; then
  LAST=$(sed -n 's/^last_step=//p' "$OUT/rt_state.txt" | head -1)
  if [ -n "$LAST" ] && [ "$LAST" -ge "$NS" ]; then echo "done ($LAST>=$NS)"; exit 0; fi
  [ -n "$LAST" ] && RES=1
fi
export LJ_LX=50 LJ_LY=50 LJ_LZ=$((LZ+60)) LJ_HALF=$($PY -c "print($LZ/2)") LJ_N=$N \
  LJ_EDGE_W=1 LJ_PERIODICITY=2 LJ_SPACING=0.5 LJ_SIGMA=2.0 \
  LJ_LAUNCH_Z=$($PY -c "print(-$LZ/2-10)") LJ_K0=$V0 \
  LJ_N_STEPS=$NS LJ_DT=0.04 LJ_CAP_ETA=-1.0 LJ_CAP_L=10.0 \
  LJ_SAVE_EVERY=0 LJ_CKPT_EVERY=200 \
  LJ_GS_DIR=$LJ/shared_gs/slab_sE_sweep/Lz${LZ}/checkpoint LJ_OUT=$HALF LJ_RESUME=$RES
[ "$HALF" = classical ] && export LJ_MASS=1.0 LJ_DELTA=0.1
srun $BASE/$HALF/run
```

`srun` with `--gres=gpu:1` masks the allocation to one card, which INQ addresses as
device 0 — no `CUDA_VISIBLE_DEVICES` juggling. `LJ_RESUME` auto-detects a short
checkpoint, so re-submitting a killed/timed-out job just continues it.

**Phase 1 — classical, all 15 jobs:**
```bash
for EV_V in "31 1.5" "54 2.0" "100 2.71" "167 3.5" "276 4.5"; do
  set -- $EV_V
  for LZ in 60 75 90; do sbatch $BASE/slurm_run.sh $1 $2 $LZ classical; done
done
```
Then **measure (§7) and confirm the classical `S(L_z)` sits on a smooth, monotone
`1/L_z` trend** (compare to the existing L_z=15/30/45 classical points in
`sweep_results.csv`) before spending the WP budget.

**Phase 2 — WP, all 15 jobs (only after classical validates):**
```bash
for EV_V in "31 1.5" "54 2.0" "100 2.71" "167 3.5" "276 4.5"; do
  set -- $EV_V
  for LZ in 60 75 90; do sbatch $BASE/slurm_run.sh $1 $2 $LZ wp; done
done
```

- `LJ_SAVE_EVERY=0` → **no density VTIs** (they are ~GBs/run and fill the disk — see
  §9). Keep it 0 — we need only energies, not frames.
- `LJ_CKPT_EVERY=200` → interior checkpoints; a killed/timed-out job resumes on
  re-submit (the script sets `LJ_RESUME` automatically). The binaries implement this.

### WP absorption gate + auto-extend (mandatory for a valid WP S)

After a WP run finishes, read the last `norm_wp` from
`results/wp/raw/observables/interactions*.csv`. If **norm_wp_final > 1e-3**, the
packet is not fully absorbed and the WP S is invalid — **extend**: raise this run's
step target by ~1250 steps (≈ 50 a.u.) and re-submit, up to ~4 rounds. With the
`slurm_run.sh` template, bump `NS` by hand for that (E, L_z) — e.g. temporarily
`NS=$((NS+1250))` — and re-`sbatch`; `LJ_RESUME` continues from the checkpoint.
(`orchestrate_sE.py` automates this on an interactive node.)

### Autonomous alternative — the orchestrator (NON-Slurm / interactive nodes only)

`orchestrate_sE.py` drives the whole pipeline (GS → runs → gate/extend → measure →
re-fit) but it launches runs **directly on local GPUs** (`subprocess` +
`CUDA_VISIBLE_DEVICES`, pool `GPUS=[0,1]`), **not** through Slurm. Use it only if you
have an interactive GPU node you hold for the duration; on a batch-scheduled cluster
prefer the `sbatch` path above. To point it at the high-L_z set:

```python
# orchestrate_sE.py:53
SWEEP = [60, 75, 90]     # 15/30/45 already done — the orchestrator skips completed ones
```
It skips any GS/run already `done` in `sweep_status.json`, appends new rows to
`sweep_results.csv`, and re-runs the extrapolation. It does NOT do classical-first
(it interleaves both halves); if you want strict classical-first, use the `sbatch`
two-phase path. Email creds may be absent (reports degrade to disk/log — fine).

### Fallback — 2 GPUs per run (if each card is still only 24 GB)

If the target's individual GPUs are *not* big enough, use `launch_2gpu.sh`
(MPI, 2 ranks, VRAM split via `gpu_rank_wrapper.sh`):

```bash
cd $BASE && ./launch_2gpu.sh 276 60 classical     # <E_eV> <L_z> <classical|wp>
```

> **HAZARD (verified on origin 2026-09-18):** a 2-GPU run must load a GS saved in
> the **same 2-domain partition**. Loading a 1-GPU GS under 2 GPUs **silently
> corrupts** the orbitals (garbage E_total ≈ 2785 vs 923 Ha). `launch_2gpu.sh`
> therefore points at `shared_gs/slab_sE_sweep/Lz{LZ}_2gpu/checkpoint` and will
> FATAL if it's missing. You must generate a **2-GPU GS** first (same `gs/run` but
> launched under `mpirun -np 2 … gpu_rank_wrapper.sh`). This is exactly the
> complication a single big GPU avoids — prefer single-GPU if you can.

---

## 7. Measure S and report / re-fit

For each completed run, compute `S(L_z)` with the skill (pass the **slab thickness**
L_z, not the box length, and the GS energy from that L_z's run_summary):

```python
import sys; sys.path.insert(0, f"{ROOT}/.claude/skills/slab-stopping-power")
import slab_stopping as SS
obs = f"{BASE}/runs/E276_Lz60/results/classical/raw/observables"
e_gs = 923.xxx     # ground_state_energy_ha from shared_gs/slab_sE_sweep/Lz60 run_summary/gs.log
print(SS.measure(obs, "classical", e_gs, 60.0))   # -> {'E_absorbed_eV':…, 'S_eV_per_Bohr':…, 'L_z':60}
# WP: SS.measure(obs_wp, "wp", e_gs, 60.0) -> also returns norm_final; require norm_final <= 1e-3
```

**Report back** (the deliverable to bring to the origin machine): one row per
(E, L_z, half) with `E_eV, v0, L_z, E_GS_ha, S_classical, S_wp, E_abs_classical,
E_abs_wp, norm_wp_final` — i.e. the same schema as `sweep_results.csv`. The simplest
is to **append the new rows to the carried `sweep_results.csv`** (it already has
L_z = 15/30/45). If you ran the orchestrator it does this for you.

**Re-fit the extrapolation on all points** (old 3 + your new 2–3):

```bash
cd $LJ/hypotheses/slab_sE_sweep
for E in 31 54 100 167 276; do $PY build_S_vs_invLz.py $E; done   # re-reads all L_z rows, re-fits, upserts sE_results.csv
$PY build_S_vs_E.py                                               # final S(E) classical vs WP
$PY build_S_vs_Lz.py                                              # direct S vs L_z convergence view
```

`build_S_vs_invLz.py` automatically uses a quadratic fit when ≥3 L_z are present
(it already does with 3; with 5–6 it is much better constrained) and writes the new
bulk `S_cl_bulk`/`S_wp_bulk` into `sE_results.csv`. The website/report figures
downstream read `sE_results.csv`, so updating it propagates the improved extrapolation.

---

## 8. Smoke test block (run this BEFORE the production grid)

Confirm the engine links and energies are sane on the target (cheap, ~minutes):

```bash
# 5-step classical smoke at the smallest high-L_z box, on the Lz60 GS (one GPU)
cd /tmp && mkdir -p smoke && cd smoke
LJ_LX=50 LJ_LY=50 LJ_LZ=120 LJ_HALF=30 LJ_N=193 LJ_EDGE_W=1 LJ_PERIODICITY=2 \
LJ_SPACING=0.5 LJ_SIGMA=2.0 LJ_LAUNCH_Z=-40 LJ_K0=4.5 \
LJ_N_STEPS=5 LJ_DT=0.04 LJ_CAP_ETA=-1.0 LJ_CAP_L=10.0 LJ_SAVE_EVERY=0 \
LJ_GS_DIR=$LJ/shared_gs/slab_sE_sweep/Lz60/checkpoint LJ_OUT=classical \
LJ_MASS=1.0 LJ_DELTA=0.1 \
INQ_SHARE_PATH=$INQ_SHARE_PATH PSEUDOPOD_SHARE_PATH=$PSEUDOPOD_SHARE_PATH \
srun --gres=gpu:1 $BASE/classical/run
# PASS: finite E_total ~ (same order as the GS energy, NOT thousands of Ha off),
#       5 steps complete, no NaN/complex-energy. A garbage E_total ~2785 vs ~923
#       means a corrupted GS load (wrong partition) — STOP and fix the GS.
```

---

## 9. Caveats register (read before launching)

1. **Engine:** CAP runs compile **only against `inq-study`** (`INQ_SOURCE=$ROOT/inq-study`).
   Stock `inq/` hard compile-errors on `double += complex`. Never edit `inq/`.
2. **GS partition must match GPU count** (single vs 2-GPU) or orbitals corrupt
   silently — see §6 hazard. Prefer single big GPU.
3. **Disk fill = crash.** Density frames (`LJ_SAVE_EVERY>0`) are ~GBs/run and can
   fill the disk; on origin a full disk SIGBUS-crashes the agent runtime. Keep
   `LJ_SAVE_EVERY=0` for these runs (we only need energies, not frames). If a run
   dies oddly, `df -h` the data partition FIRST.
4. **Detach long runs** with `setsid`/`nohup`, not the harness background runner
   (which is SIGKILLed when the session ends).
5. **WP S is only valid once absorbed** (`norm_wp_final ≤ 1e-3`); extend otherwise.
   Don't quote a WP bulk point built from an unabsorbed packet.
6. **This is the slab-DEPOSIT method** (S = E_absorbed/L_z). Do **not** try to
   extract quantum S from −d(E_proj_KE)/ds — the WP orbital isn't identifiable as
   "the projectile". Total quantum stopping = total electronic energy deposition,
   which is exactly what the CAP-absorbed E_absorbed measures. (The
   `light-projectile-stopping` rule about KE-slope windows does NOT apply here.)
7. **Resume, don't restart:** every run checkpoints its last step; to add steps
   re-invoke with `LJ_RESUME=1` and a larger `LJ_N_STEPS`. A run that crashed
   *beyond* its last checkpoint leaves orphan outputs — if a resume aborts on an
   output-file collision, delete outputs with step > last_step and relaunch.
8. **Round N with the formula** `round(0.001288·2500·L_z)` (→ 60:193, 75:241,
   90:290). n0 is then recomputed from the integer N inside the binary; r_s stays
   ≈ 4.18. Verify r_s in the GS log.

---

## 10. Success criteria

- GS exists for each new L_z (checkpoint + r_s ≈ 4.18 in the log).
- For each (E, L_z) both halves complete; WP runs reach `norm_wp_final ≤ 1e-3`.
- `SS.measure(...)` returns a finite `S_eV_per_Bohr` for every run (and passes its
  internal closure gate for classical / norm gate for WP).
- New rows appended to `sweep_results.csv` in the documented schema.
- `build_S_vs_invLz.py` re-fits each energy on the enlarged L_z set and updates
  `sE_results.csv`; `S_vs_E.png` / `S_vs_Lz.png` regenerate without error.
- **Deliverable to carry back:** the updated `sweep_results.csv` (+ the new run
  output dirs if feasible) and the per-run S values.

---

## 11. Origin-machine reference paths (source of truth to rsync/consult)

```
$ROOT = /local/data/public/skcb2/tddft                 (on origin: lovelace)
Campaign dir   : $ROOT/ResearchProject/systems/localised_jellium/scripts/slab_sE_sweep/
  classical/run.cpp, wp/run.cpp         # TDDFT run definitions (env contract in header comments)
  orchestrate_sE.py                     # autonomous driver (constants §2 live here)
  launch_2gpu.sh, gpu_rank_wrapper.sh   # 2-GPU MPI fallback
  sweep_results.csv                     # per-run S (L_z=15/30/45 now) — CARRY THIS
  sE_results.csv                        # bulk S(E) — regenerated by the re-fit
GS binary src  : $ROOT/.../localised_jellium/scripts/campaign_autorun/gs/run.cpp
GS checkpoints : $ROOT/.../localised_jellium/shared_gs/slab_sE_sweep/{Lz15,Lz30,Lz45,Lz60,Lz75}
Fit/plots      : $ROOT/.../localised_jellium/hypotheses/slab_sE_sweep/build_S_vs_{invLz,E,Lz}.py
S measurement  : $ROOT/.claude/skills/slab-stopping-power/slab_stopping.py  (+ SKILL.md)
Engine         : $ROOT/inq-study  (CAP)   |  headers: $ROOT/inq-stack/include/inqkit/
Post-process   : $ROOT/inq-stack/python/inqview  (pip install -e inq-stack/ into a venv)
```

Open the `run.cpp` header comments and `orchestrate_sE.py` top block on the target
to re-confirm the env contract before launching — they are the authoritative source
if anything here is ambiguous.
```
