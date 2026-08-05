#!/bin/bash
# lzb-params.sh — box presets + step tables for the lz_bulk_sweep campaign.
# Plan: docs/plans/jellium-slab-extend-Lz.md
#
# SOURCE this with LZB_CFG set (and optionally VIDX=0..3) — it exports the box
# geometry and, when VIDX is given, the per-velocity step counts and cadences.
# It is the bash mirror of shared/configs/lzb_boxes.hpp and of the tables in
# hypotheses/lz_bulk_sweep/lzb_stopping.py; the three MUST agree. Values are
# LITERAL (not recomputed) so `complete` checks can never drift from the
# dispatcher through rounding.
#
#   N_STEPS = round(4.36 * (|launch_z| + L_z/2) / (v * dt)), dt = 0.04
#   SAVE ~ N/100 frames (VTI cadence ~3x coarser than sigma56 — disk decision,
#   user 2026-08-05), WF ~ N/40, CKPT = N/5 (>= 4 retained checkpoints).
#
# Geometry is per sigma FAMILY (user 2026-08-05): sigma = 0.5 runs replicate the
# 85-box vacuum layout (standoff 11.5, face->CAP 17.5) to match their existing
# L_slab = 25 anchor; sigma = 5 runs replicate the 105-box layout (15 / 27.5)
# to match the sigma56_sv anchor.

case "${LZB_CFG:?lzb-params.sh: LZB_CFG must be set}" in
  s0p5_L15)
    LZ=75;  HALF=7.5;  N_E=60;  EXTRA=15; LAUNCH=-19.0; SIGMA=0.5
    GS_TAG="slab_n60_L35x35x75_dx0p4_per2"
    N_TAB=(3079 2463 2053 1760); SAVE_TAB=(31 25 21 18)
    WF_TAB=(77 62 51 44);        CKPT_TAB=(616 493 411 352) ;;
  s0p5_L35)
    LZ=95;  HALF=17.5; N_E=140; EXTRA=34; LAUNCH=-29.0; SIGMA=0.5
    GS_TAG="slab_n140_L35x35x95_dx0p4_per2"
    N_TAB=(4169 3335 2780 2382); SAVE_TAB=(42 33 28 24)
    WF_TAB=(104 83 70 60);       CKPT_TAB=(834 667 556 476) ;;
  s5p0_L15)
    LZ=95;  HALF=7.5;  N_E=60;  EXTRA=15; LAUNCH=-22.5; SIGMA=5.0
    GS_TAG="slab_n60_L35x35x95_dx0p4_per2"
    N_TAB=(3815 3052 2543 2180); SAVE_TAB=(38 31 25 22)
    WF_TAB=(95 76 64 55);        CKPT_TAB=(763 610 509 436) ;;
  s5p0_L35)
    LZ=115; HALF=17.5; N_E=140; EXTRA=34; LAUNCH=-32.5; SIGMA=5.0
    GS_TAG="slab_n140_L35x35x115_dx0p4_per2"
    N_TAB=(4905 3924 3270 2803); SAVE_TAB=(49 39 33 28)
    WF_TAB=(123 98 82 70);       CKPT_TAB=(981 785 654 561) ;;
  *)
    echo "FATAL: unknown LZB_CFG '${LZB_CFG}' (want s0p5_L15, s0p5_L35, s5p0_L15, s5p0_L35)" >&2
    exit 2 ;;
esac

V_TAB=(2.0 2.5 3.0 3.5)
export LZ HALF N_E EXTRA LAUNCH SIGMA GS_TAG

if [ -n "${VIDX:-}" ]; then
  case "$VIDX" in (*[!0-9]*|"") echo "FATAL: bad VIDX '$VIDX'" >&2; exit 2 ;; esac
  [ "$VIDX" -ge 0 ] && [ "$VIDX" -le 3 ] || { echo "FATAL: VIDX out of range: $VIDX" >&2; exit 2; }
  V="${V_TAB[$VIDX]}"; NSTEPS="${N_TAB[$VIDX]}"
  SAVE="${SAVE_TAB[$VIDX]}"; WF="${WF_TAB[$VIDX]}"; CKPT="${CKPT_TAB[$VIDX]}"
  export V NSTEPS SAVE WF CKPT
fi
