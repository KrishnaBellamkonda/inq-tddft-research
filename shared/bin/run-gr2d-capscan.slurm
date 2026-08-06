#!/bin/bash
#SBATCH --job-name=gr2d-capscan
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-gr2d-capscan.slurm — Phase 1 CAP (eta, W, v) scan, two-definitions campaign.
# Plan: docs/plans/real-material-stopping-comparison.md
#
# Grid: eta in {1,2,3} Ha x W in {12,16,20} Bohr at Lz=80 (W=20 also at Lz=100),
# each at the two extreme velocities first (k=1.05, 4.70) — a (eta,W) point is
# only run at the FULL velocity grid {1.05,1.36,1.92,2.71,3.83,4.70} if both
# extremes pass residue < 5e-3 (cheap pre-filter; final gate 1e-3 in analyzer).
# Vacuum runs are 1-electron: ~minutes each. Idempotent (skips completed tags).
set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
[[ -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]] || { echo "bad REPO_ROOT" >&2; exit 1; }
. "$REPO_ROOT/shared/bin/csd3-env.sh"
export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"
. "$REPO_ROOT/shared/config.sh"
export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"
mkdir -p "$TMPDIR"

RUN_DIR="$REPO_ROOT/ResearchProject/systems/graphene/scripts/twodef_sv/capscan"
cd "$RUN_DIR"

echo "== gr2d-capscan job ${SLURM_JOB_ID:-none} on $(hostname) $(date) =="

# build once (inq-run compiles then execs the default variant — harmless)
if [ ! -x ./run ] || [ run.cpp -nt ./run ]; then
  CS_ETA=-1.0 CS_CAP_L=16 CS_LZ=80 CS_K0=2.71 inq-run || {
    [ -x ./run ] || { echo "[capscan] build failed" >&2; exit 1; }
  }
fi

VGRID="1.05 1.36 1.92 2.71 3.83 4.70"
export CS_SIGMA="${CS_SIGMA:-2.0}"    # production sigma_WP (decision 2026-08-05:
export CS_BUFFER="${CS_BUFFER:-2.0}"  # sigma=0.5 is not beam-like below ~250 eV)

mktag () {  # eta W Lz k0 -> run.cpp's tag (C++ %g-style number formatting)
  python3 -c "print(f\"eta{abs($1):g}_W${2}_Lz${3}_k{$4:g}_s{float('$CS_SIGMA'):g}\")"
}

run_pt () {  # eta W Lz k0
  local eta="$1" W="$2" lz="$3" k="$4"
  local tag; tag="$(mktag "$eta" "$W" "$lz" "$k")"
  if [ -f "results/$tag/run_summary.txt" ] \
     && grep -q "run_completed = true" "results/$tag/run_summary.txt"; then
    echo " [skip] $tag"; return 0
  fi
  echo " [run ] $tag  $(date +%H:%M:%S)"
  env INQ_SHARE_PATH="$INQ_SHARE_PATH" PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
      CS_ETA="$eta" CS_CAP_L="$W" CS_LZ="$lz" CS_K0="$k" \
      CS_SIGMA="$CS_SIGMA" CS_BUFFER="$CS_BUFFER" ./run
}

residue () {  # eta W Lz k0 -> prints norm_final (or 1.0 if missing)
  local tag; tag="$(mktag "$1" "$2" "$3" "$4")"
  awk -F' = ' '/^norm_final/{print $2; f=1} END{if(!f) print 1.0}' \
      "results/$tag/run_summary.txt" 2>/dev/null || echo 1.0
}

overall=0
for lz in 80 100; do
  for W in 12 16 20; do
    [ "$lz" = "100" ] && [ "$W" != "20" ] && continue   # Lz=100 only for W=20
    for eta in -1.0 -2.0 -3.0; do
      run_pt "$eta" "$W" "$lz" 1.05 || { overall=$?; continue; }
      run_pt "$eta" "$W" "$lz" 4.70 || { overall=$?; continue; }
      r_lo=$(residue "$eta" "$W" "$lz" 1.05)
      r_hi=$(residue "$eta" "$W" "$lz" 4.70)
      pass=$(python3 -c "print(1 if $r_lo < 5e-3 and $r_hi < 5e-3 else 0)")
      echo " [pre ] eta=$eta W=$W Lz=$lz residue(lo)=$r_lo residue(hi)=$r_hi pass=$pass"
      if [ "$pass" = "1" ]; then
        for k in $VGRID; do run_pt "$eta" "$W" "$lz" "$k" || overall=$?; done
      fi
    done
  done
done

echo "== capscan sweep done ($(date)); summaries: =="
grep -H "norm_final" results/*/run_summary.txt 2>/dev/null | sort
exit $overall
