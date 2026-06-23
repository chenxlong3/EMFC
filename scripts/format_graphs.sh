#!/usr/bin/env bash
# Format graphs: generate q/tau/lam and forward/reverse graph files.
#
# Usage:
#   bash scripts/format_graphs.sh
#   PDIST=uniform PEDGE=0.1 bash scripts/format_graphs.sh

set -euo pipefail

_SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$_SCRIPT_DIR/.." || exit 1
# shellcheck source=scripts-topo/_run_log.sh
. "$_SCRIPT_DIR/../scripts-topo/_run_log.sh"
SCRIPT_TAG="[$(basename "$0" .sh)]"

# ---- graph list (edit here) ----
GRAPH_NAMES=(
    # ca-GrQc
    # soc-Epinions
    soc-facebook
)

# ---- configurable ----
DIR=${DIR:-"graphInfo"}
PDIST=${PDIST:-"wc"}
PEDGE=${PEDGE:-0.1}

if [ ! -x ./frim ]; then
    run_log "$SCRIPT_TAG building ./frim"
    make
fi

for name in "${GRAPH_NAMES[@]}"; do
    graph_file="${DIR}/${name}"
    if [ ! -f "$graph_file" ]; then
        run_log "$SCRIPT_TAG skip $name (not found: $graph_file)"
        continue
    fi

    run_log "$SCRIPT_TAG dataset=$name func=format pdist=$PDIST"
    if [ "$PDIST" = "uniform" ]; then
        ./frim -func=format -gname="$name" -dir="$DIR" -pdist=uniform -pedge="$PEDGE"
    else
        ./frim -func=format -gname="$name" -dir="$DIR" -pdist="$PDIST"
    fi
done

run_log "$SCRIPT_TAG done (${#GRAPH_NAMES[@]} graphs in list)"
