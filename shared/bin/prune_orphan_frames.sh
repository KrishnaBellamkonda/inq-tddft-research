#!/bin/bash
# prune_orphan_frames.sh — make a checkpoint resume actually restartable.
#
#   Usage:  prune_orphan_frames.sh <results/<run_name>>        # reads rt_state.txt
#
# THE PROBLEM IT SOLVES (2026-08-03). Checkpoints are written every CKPT steps
# (581-872 here) but VTI frames every SAVE steps (8-14). So a run killed between
# checkpoints leaves frames for steps the checkpoint does NOT cover:
#
#     ckpt_step001744        <- resume starts here
#     density_t001750.vti ... density_t001960.vti   <- 16 ORPHAN frames
#
# `inqkit::io::VTIImageDataWriter` is deliberately overwrite=false (segment
# outputs must never clobber earlier data, .claude/rules/final-timestep-checkpoint.md),
# so on resume the first frame write past the checkpoint throws:
#
#     what():  VTIImageDataWriter: file already exists and overwrite=false:
#              results/s5p0_v2p0/raw/vti/density_total/density_t001750.vti
#
# and the run aborts ~45 s in. EVERY resume of a mid-checkpoint kill hits this.
# It cost 5 of 5 resume jobs on 2026-08-03 before it was diagnosed.
#
# WHY DELETING IS SAFE. Propagation from a checkpoint is deterministic, so the
# resumed segment recomputes exactly these steps. The orphans are strictly
# superseded, not lost data. The CSV side needs no such treatment: `_concat`
# already resolves the overlap with drop_duplicates(subset="step", keep="last").
#
# ONLY frames with step > last_step are touched. Everything at or below the
# checkpoint — and every checkpoint directory — is left alone.
set -euo pipefail

RUN="${1:?usage: prune_orphan_frames.sh <run_dir>}"
STATE="$RUN/rt_state.txt"

if [ ! -f "$STATE" ]; then
  echo "  prune: no rt_state.txt in $RUN — nothing to prune (fresh start)"
  exit 0
fi

LAST=$(sed -n 's/^last_step=\([0-9]*\).*/\1/p' "$STATE" | head -1)
if [ -z "$LAST" ]; then
  echo "  prune: could not read last_step from $STATE — refusing to guess" >&2
  exit 1
fi

VTI="$RUN/raw/vti"
[ -d "$VTI" ] || { echo "  prune: no raw/vti in $RUN"; exit 0; }

n=0
while IFS= read -r f; do
  # frame files are <kind>_t<zero-padded step>.vti
  step=$(basename "$f" | sed -n 's/.*_t\([0-9]\{1,\}\)\.vti$/\1/p')
  [ -z "$step" ] && continue
  if [ "$((10#$step))" -gt "$LAST" ]; then
    rm -f "$f"
    n=$((n + 1))
  fi
done < <(find "$VTI" -type f -name '*.vti')

echo "  prune: removed $n orphan frame(s) with step > $LAST from $RUN/raw/vti"
