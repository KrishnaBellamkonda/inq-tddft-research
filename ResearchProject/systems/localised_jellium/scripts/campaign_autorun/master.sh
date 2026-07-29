#!/bin/bash
# ============================================================================
# HEADLESS master dispatcher — localised-jellium GS ladder H1->H5 + cumulative.
# Runs UNATTENDED on GPU $GPU (default 1; GPU0 is the user's qsp_phase5 run).
# Each phase: run sims -> analyse_phase.py (plots + EMAILS, 4-part structure) ->
# continue. Robust: a phase failure emails a failure notice and the chain
# continues (data stays on disk for later re-analysis). H0 already done+emailed.
# ============================================================================
ROOT=/local/data/public/skcb2/tddft
export INQ_SHARE_PATH=$ROOT/inq/install/share
export PSEUDOPOD_SHARE_PATH=$INQ_SHARE_PATH/pseudopod
export INQ_SOURCE=$ROOT/inq-study
GPU=${GPU:-1}
PY=$ROOT/venv/bin/python3
INQRUN=$ROOT/shared/bin/inq-run
CA=$ROOT/ResearchProject/systems/localised_jellium/scripts/campaign_autorun
GSBIN=$CA/gs/run; WPBIN=$CA/wp/run; CLBIN=$CA/classical/run
RUNS=$CA/runs; mkdir -p "$RUNS"
GS120_P3=$ROOT/ResearchProject/systems/localised_jellium/scripts/h0_base_difference/gs/results  # H0 periodicity-3 L_z=120 GS
GS120_P3_CKPT=$ROOT/ResearchProject/systems/localised_jellium/shared_gs/slab_n82_L50x50x120
ANL="$PY $CA/analyse_phase.py"
log(){ echo "[$(date '+%F %T')] $*"; }

run_gs(){ # $1=rundir $2..=env assignments
  local rd=$1; shift; mkdir -p "$rd"; ( cd "$rd" && env "$@" LJ_GS_DIR="$rd/checkpoint" CUDA_VISIBLE_DEVICES=$GPU "$GSBIN" ) ; }
run_proj(){ # $1=bin $2=rundir $3..=env
  local bin=$1 rd=$2; shift 2; mkdir -p "$rd"; ( cd "$rd" && env "$@" CUDA_VISIBLE_DEVICES=$GPU "$bin" ) ; }
fail_email(){ $PY - "$1" "$2" <<'PY'
import sys
from inqview.email import send_run_email
send_run_email(subject=f"[localised-jellium GS] {sys.argv[1]} — PHASE FAILED",
  body=f"Phase {sys.argv[1]} failed in the headless run.\n{sys.argv[2]}\nData (if any) is under the runs dir; chain continues.",
  attachments=[], to="chiddukanna@gmail.com")
PY
}

log "MASTER start (GPU=$GPU)"

# ---- build wp + classical binaries (gs already built from the smoke) -------
log "build wp binary"; ( cd "$CA/wp" && CUDA_VISIBLE_DEVICES=$GPU LJ_GS_DIR="$GS120_P3_CKPT" LJ_OUT=buildsmoke "$INQRUN" ) >"$CA/wp/build.log" 2>&1 || log "wp build warn (see build.log)"
log "build classical binary"; ( cd "$CA/classical" && CUDA_VISIBLE_DEVICES=$GPU LJ_GS_DIR="$GS120_P3_CKPT" LJ_OUT=buildsmoke "$INQRUN" ) >"$CA/classical/build.log" 2>&1 || log "cl build warn (see build.log)"

# ---- H1: edge-width sweep (GS-only, periodicity 3, L_z=90) ------------------
{ log "H1 edge-width sweep"
  for w in 0 0.5 1 1.5 2; do
    log "  H1 gs w=$w"; run_gs "$RUNS/h1/gs_w$w" LJ_LX=50 LJ_LY=50 LJ_LZ=90 LJ_HALF=12.5 LJ_N=82 LJ_EDGE_W=$w LJ_PERIODICITY=3 LJ_TAG=h1_w$w || true
  done
  $ANL --phase H1 --base "$RUNS/h1" || fail_email H1 "analysis error"
} || fail_email H1 "phase error"

# ---- H2: Lz sweep (periodicity 3) + open-z GS (periodicity 2, L_z=120) ------
{ log "H2 Lz sweep + open-z"
  for lz in 50 70 90 120; do
    log "  H2 gs lz=$lz"; run_gs "$RUNS/h2/gs_lz$lz" LJ_LX=50 LJ_LY=50 LJ_LZ=$lz LJ_HALF=12.5 LJ_N=82 LJ_EDGE_W=0 LJ_PERIODICITY=3 LJ_TAG=h2_lz$lz || true
  done
  log "  H2 gs periodicity-2 L_z=120 (reused by H4/H5)"
  run_gs "$RUNS/h2/gs_p2_lz120" LJ_LX=50 LJ_LY=50 LJ_LZ=120 LJ_HALF=12.5 LJ_N=82 LJ_EDGE_W=0 LJ_PERIODICITY=2 LJ_TAG=h2_p2_lz120 || true
  $ANL --phase H2 --base "$RUNS/h2" || fail_email H2 "analysis error"
} || fail_email H2 "phase error"
GS120_P2_CKPT=$RUNS/h2/gs_p2_lz120/checkpoint
GS120_P2_RES=$RUNS/h2/gs_p2_lz120/results

# ---- H3: thickness sweep (GS-only, periodicity 3, L_z=90, N scaled) ---------
{ log "H3 thickness sweep"
  for an in "7.5 50" "12.5 82" "17.5 114" "22.5 148"; do
    set -- $an; a=$1; N=$2
    log "  H3 gs a=$a N=$N"; run_gs "$RUNS/h3/gs_a${a}_N${N}" LJ_LX=50 LJ_LY=50 LJ_LZ=90 LJ_HALF=$a LJ_N=$N LJ_EDGE_W=0 LJ_PERIODICITY=3 LJ_TAG=h3_a${a} || true
  done
  $ANL --phase H3 --base "$RUNS/h3" || fail_email H3 "analysis error"
} || fail_email H3 "phase error"

# ---- H4: WP energetics, r x {periodicity 3,2}, L_z=120 ----------------------
{ log "H4 WP energetics"
  for r in 4 16 28 40; do z=$(echo "scale=2; -12.5 - $r" | bc)
    run_proj "$WPBIN" "$RUNS/h4/wp_r${r}_p3" LJ_OUT=wp_r${r}_p3 LJ_LZ=120 LJ_PERIODICITY=3 LJ_LAUNCH_Z=$z LJ_K0=0 LJ_SIGMA=0.5 LJ_GS_DIR="$GS120_P3_CKPT" || true
    if [ -d "$GS120_P2_CKPT" ]; then
      run_proj "$WPBIN" "$RUNS/h4/wp_r${r}_p2" LJ_OUT=wp_r${r}_p2 LJ_LZ=120 LJ_PERIODICITY=2 LJ_LAUNCH_Z=$z LJ_K0=0 LJ_SIGMA=0.5 LJ_GS_DIR="$GS120_P2_CKPT" || true
    fi
  done
  $ANL --phase H4 --base "$RUNS/h4" --gs120-p3 "$GS120_P3" --gs120-p2 "$GS120_P2_RES" || fail_email H4 "analysis error"
} || fail_email H4 "phase error"

# ---- H5: classical mirror, r x {periodicity 3,2}, L_z=120 -------------------
{ log "H5 classical mirror"
  for r in 4 16 28 40; do z=$(echo "scale=2; -12.5 - $r" | bc)
    run_proj "$CLBIN" "$RUNS/h5/cl_r${r}_p3" LJ_OUT=cl_r${r}_p3 LJ_LZ=120 LJ_PERIODICITY=3 LJ_LAUNCH_Z=$z LJ_GS_DIR="$GS120_P3_CKPT" || true
    if [ -d "$GS120_P2_CKPT" ]; then
      run_proj "$CLBIN" "$RUNS/h5/cl_r${r}_p2" LJ_OUT=cl_r${r}_p2 LJ_LZ=120 LJ_PERIODICITY=2 LJ_LAUNCH_Z=$z LJ_GS_DIR="$GS120_P2_CKPT" || true
    fi
  done
  $ANL --phase H5 --base "$RUNS/h5" --h4-base "$RUNS/h4" --gs120-p3 "$GS120_P3" --gs120-p2 "$GS120_P2_RES" || fail_email H5 "analysis error"
} || fail_email H5 "phase error"

# ---- cumulative email (all highlight plots) --------------------------------
log "cumulative email"
$PY - <<'PY' || true
from pathlib import Path
from inqview.email import send_run_email
ca = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium")
pngs = []
for p in ["hypotheses/h0_base_difference/H0_base_difference.png",
          "scripts/campaign_autorun/runs/h1/H1_edge_model.png",
          "scripts/campaign_autorun/runs/h2/H2_gs_convergence.png",
          "scripts/campaign_autorun/runs/h3/H3_surface_energetics.png",
          "scripts/campaign_autorun/runs/h4/H4_wp_energetics.png",
          "scripts/campaign_autorun/runs/h5/H5_classical_subtraction.png"]:
    f = ca / p
    if f.exists(): pngs.append(str(f))
send_run_email(
  subject="[localised-jellium GS] CAMPAIGN COMPLETE — all phases (cumulative)",
  body=("The localised-jellium GS parameter-study ladder has finished its headless run.\n"
        "Attached: the highlight plot from each completed phase (H0-H5).\n\n"
        "H0 base gap (artifact-dominated) -> H1 edge model -> H2 GS convergence/open-z\n"
        "-> H3 surface energetics (sigma_s RAW, E_self caveat) -> H4 E_SIE + BC verdict\n"
        "-> H5 classical image error (route-2 ghost-bg flagged).\n\n"
        "Review the per-phase emails for the hypothesis/method/plot/conclusion of each.\n"
        "Flagged for your judgement: H2 work-function (Phi extractor), H3 sigma_s E_self\n"
        "correction, H5 ghost-background term. Cumulative notebook to follow."),
  attachments=pngs, to="chiddukanna@gmail.com")
print("cumulative sent with", len(pngs), "plots")
PY
log "MASTER done"
