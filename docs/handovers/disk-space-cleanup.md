# Handover — disk-space cleanup for two new campaigns

**Rolling file. Latest milestone at top.**
**Repo:** `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research` (branch `quantum-stopping-power`)
**Task:** free space on RDS for two upcoming campaigns by pruning heavy files
from properly-analysed runs. User rule: nothing deleted without explicit
approval; cylindrical_jellium (proximity ladder + annular_sv) untouched.

---

## 2026-08-05 — survey + approved cleanup executed

### Survey result (before cleanup)

`ResearchProject/systems` totalled **795G**:

| Location | Size | Status |
|---|---|---|
| `localised_jellium/scripts/sigma56_sv` | 513G | σ=6 analysed (twin+synthesis notebooks, S CSVs); σ=5 runs complete, notebooks built except `run_wp_s5_v2.0` |
| `localised_jellium/scripts/slab_ks_wrap` | 84G | E_absorbed S(v) computed for all runs; σ=0.5 runs not settled |
| `localised_jellium/scripts/wp_highdensity_sv` | 45G | synthesis complete, cross-validated vs slab_ks_wrap |
| `jellium/scripts` (bulk_ks_stopping* + bulk_t0_density) | 70G | campaign wrapped 2026-08-04 (RUN_RECORD) |
| `cylindrical_jellium` | 48G | EXCLUDED (user: analysis pending) |
| `vacuum/scripts/wp_selfinteraction` | 13G | EXCLUDED (active SIC campaign) |
| GS stores (`jellium/checkpoints`, `save_gs`, `localised_jellium/shared_gs`) | ~17G | KEPT (shared ground states) |

Composition finding: VTIs were only ~1/3 of the space; RT checkpoints dominated
(interior `ckpt_step*` snapshots 276G + final `checkpoint/` dirs ~210G).

### Approved and executed (user approval 2026-08-05)

1. **Interior `ckpt_step*` checkpoints deleted — ~276G freed.**
   All 115 dirs in `sigma56_sv` (232G, 89 dirs) and `slab_ks_wrap` (44G,
   26 dirs). Safety gate: each deleted dir verified to have a non-empty
   sibling FINAL `checkpoint/` (115/115 passed, 0 skipped), so every run
   remains extendable per `.claude/rules/final-timestep-checkpoint.md`.
   Verified post-delete: `find ... -name 'ckpt_step*' | wc -l` → 0.

2. **VTI thinning to quarter cadence — ~196G freed.**
   User decision: do NOT delete all VTIs; keep every nth so density remains
   reconstructable at lower frequency. Implemented: per `raw/vti/<kind>/`
   dir, frames sorted by step stamp `_tNNNNNN.vti`; kept first + last +
   every 4th; dirs with <8 frames and un-stamped statics (e.g.
   `density_gs_system.vti`) untouched. Applied to `sigma56_sv`
   (194.1G / 9,425 frames deleted, 3,273 kept) and `wp_highdensity_sv`
   (1.9G / 210 deleted, 159 kept). Script:
   session scratchpad `thin_vti.py` (KEEP_EVERY=4, MIN_FRAMES=8).

### NOT approved / not done

- Final `checkpoint/` dirs everywhere KEPT (sigma56_sv 55G, slab_ks_wrap 44G,
  wp_highdensity_sv 42G, jellium bulk sweeps 71G) — user did not approve;
  deleting would break the extend-don't-restart guarantee.
- `run_wp_s5_v2.0` notebook still missing in
  `localised_jellium/hypotheses/sigma56_sv/` although the run
  (`scripts/sigma56_sv/wp/results/s5p0_v2p0`) completed — build it before any
  further VTI pruning of that run.
- All CSVs (`observables*`, `interactions*`, `projectile*`), `rt_state.txt`,
  `run_summary.txt`, notebooks, figures, hypotheses folders: untouched.

### Net effect (verified post-delete with du/df)

- `sigma56_sv` 513G → **116G**; `slab_ks_wrap` 84G → **42G**;
  `wp_highdensity_sv` 45G → **43G** — ~440G freed.
- Filesystem after cleanup: **660G available** of the 1.0T RDS quota
  (`/rds/user/skcb2` at 36% used).
- Every run still resumable from its final checkpoint; every VTI series still
  present at ¼ time resolution.
