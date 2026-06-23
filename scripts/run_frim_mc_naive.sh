#!/usr/bin/env bash
# Run FRIM naive MC solver and evaluate the resulting xi.
#
# Usage:
#   bash scripts/run_frim_mc_naive.sh
#   GNAME=ca-GrQc_1 FRIM_MC=20 bash scripts/run_frim_mc_naive.sh
#   RUN_NAIVE=0 EVAL_MC=10000 bash scripts/run_frim_mc_naive.sh   # eval only

set -euo pipefail

_SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$_SCRIPT_DIR/.." || exit 1
# shellcheck source=scripts-topo/_run_log.sh
. "$_SCRIPT_DIR/../scripts-topo/_run_log.sh"
SCRIPT_TAG="[$(basename "$0" .sh)]"

# ---- configurable ----
GNAME=${GNAME:-"soc-facebook"}
DIR=${DIR:-"graphInfo"}
PDIST=${PDIST:-"wc"}
RAND_SEED=${RAND_SEED:-1}
XI_LO=${XI_LO:-0.5}
FRIM_MC=${FRIM_MC:-100}
FRIM_SWEEPS=${FRIM_SWEEPS:-4}
EVAL_MC=${EVAL_MC:-10000}
RUN_NAIVE=${RUN_NAIVE:-1}
RUN_EVAL=${RUN_EVAL:-1}
LOG_DIR=${LOG_DIR:-"./result/logs/frim_mc_naive"}

# ---- derived ----
OUT_NAME="${GNAME}_${RAND_SEED}_frim_mc_naive_${PDIST}_mc${FRIM_MC}"
XI_FILE="./result/xi/xi_${OUT_NAME}"
EVAL_FILE="./result/evaluation/eval_${OUT_NAME}"
LOG_FILE="${LOG_DIR}/${OUT_NAME}.log"

mkdir -p "$LOG_DIR"

run_log "$SCRIPT_TAG gname=$GNAME pdist=$PDIST rand_seed=$RAND_SEED frim_mc=$FRIM_MC sweeps=$FRIM_SWEEPS eval_mc=$EVAL_MC"

if [ ! -x ./frim ]; then
    run_log "$SCRIPT_TAG building ./frim"
    make
fi

if [ ! -f "${DIR}/${GNAME}.nodehyper.vec" ]; then
    run_log "$SCRIPT_TAG formatting graph (missing node hyper params)"
    ./frim -func=format -gname="$GNAME" -dir="$DIR" -pdist="$PDIST"
fi

{
    if [ "$RUN_NAIVE" -eq 1 ]; then
        echo "=== FRIM MC naive ==="
        echo "out: $OUT_NAME"
        ./frim \
            -func=frim_mc_naive \
            -gname="$GNAME" \
            -dir="$DIR" \
            -pdist="$PDIST" \
            -rand_seed="$RAND_SEED" \
            -xi_lo="$XI_LO" \
            -frim_mc="$FRIM_MC" \
            -frim_sweeps="$FRIM_SWEEPS"
        echo "xi saved: $XI_FILE"
    else
        echo "=== skip naive (RUN_NAIVE=0) ==="
    fi

    if [ "$RUN_EVAL" -eq 1 ]; then
        echo "=== eval ==="
        ./frim \
            -func=eval \
            -gname="$GNAME" \
            -dir="$DIR" \
            -pdist="$PDIST" \
            -rand_seed="$RAND_SEED" \
            -method=frim_mc_naive \
            -frim_mc="$FRIM_MC" \
            -eval_mc="$EVAL_MC"
        echo "eval saved: $EVAL_FILE"
        if [ -f "$EVAL_FILE" ]; then
            echo "--- eval result ---"
            cat "$EVAL_FILE"
        fi
    else
        echo "=== skip eval (RUN_EVAL=0) ==="
    fi
} 2>&1 | tee "$LOG_FILE"

run_log "$SCRIPT_TAG done log=$LOG_FILE"
