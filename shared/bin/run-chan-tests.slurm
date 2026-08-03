#!/bin/bash
#SBATCH --job-name=chan-tests
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-chan-tests.slurm — the CORRECTNESS GATE for the annular-tube channeling
# twin. Plan: docs/plans/cylindrical-channeling-ks-stopping.md
#
# Builds and runs the library kernels this study introduces or newly relies on:
#
#   test_radial_occupancy_engine    f_bore(t) — the fraction of the wavepacket
#                                   still inside the bore. This is the study's
#                                   central premise made falsifiable; if it is
#                                   wrong, "the packet is channeling" cannot be
#                                   checked at all and the fit window is guessed.
#   test_projectile_force_minimum_image_engine
#                                   the classical twin launches 0.71 sigma_pot
#                                   from the -z face. A clipped Gaussian loses
#                                   24 % of its charge ASYMMETRICALLY across the
#                                   +/-delta finite difference, so the error does
#                                   not cancel out of the force.
#   test_gaussian_minimum_image_engine   the charge kernel that force depends on.
#   test_slab_occupancy_engine      the sibling observable radial_occupancy was
#                                   modelled on; runs here so a shared-pattern
#                                   regression is caught.
#   test_projectile (pure)          set_R / wrap_into_cell / unwrapped path.
#
# PREREQUISITE: both build trees must already be CONFIGURED on the LOGIN node —
# compute nodes have no outbound network (docs/handovers/csd3-device-setup.md):
#   cmake -S inq-stack/tests/include/engine -B inq-stack/tests/include/engine/build
#
# Nothing downstream should run if this fails: chain with afterok.
set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "ERROR: REPO_ROOT='$REPO_ROOT' is not the repo root." >&2; exit 1
fi
cd "$REPO_ROOT" || exit 1
. "$REPO_ROOT/shared/bin/csd3-env.sh" || { echo "csd3-env.sh failed"; exit 1; }

PURE_DIR="$REPO_ROOT/inq-stack/tests/include/build"
ENG_DIR="$REPO_ROOT/inq-stack/tests/include/engine/build"
JOBS="${SLURM_CPUS_PER_TASK:-12}"

export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"

echo "=============================================="
echo " job     : ${SLURM_JOB_ID:-none}  chan-tests"
echo " node    : $(hostname)"
echo " gpu     : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo " started : $(date)"
echo "=============================================="

for d in "$PURE_DIR" "$ENG_DIR"; do
  [[ -f "$d/CMakeCache.txt" ]] || { echo "ERROR: $d not configured (login node)." >&2; exit 1; }
done

rc=0

echo
echo "########## PURE: test_projectile ##########"
cmake --build "$PURE_DIR" --target test_projectile -j"$JOBS" || rc=1
"$PURE_DIR/test_projectile" || rc=1

for t in test_radial_occupancy_engine \
         test_inq_field_export_engine \
         test_wp_minimum_image_engine \
         test_projectile_force_minimum_image_engine \
         test_gaussian_minimum_image_engine \
         test_slab_occupancy_engine \
         test_wp_sic_engine; do
  echo
  echo "########## ENGINE: $t ##########"
  if ! cmake --build "$ENG_DIR" --target "$t" -j"$JOBS" > "$ENG_DIR/build_$t.log" 2>&1; then
    echo "BUILD FAILED — tail of $ENG_DIR/build_$t.log:"
    tail -60 "$ENG_DIR/build_$t.log"
    rc=1
    continue
  fi
  echo "build OK"
  ctest --test-dir "$ENG_DIR" -R "^$t\$" --output-on-failure --timeout 1800 || rc=1
done

echo "=============================================="
echo " exit status : $rc  ($([ $rc -eq 0 ] && echo 'ALL GATES PASSED' || echo 'FAILURES — production must not run'))"
echo " finished    : $(date)"
echo "=============================================="
exit $rc
