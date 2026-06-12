#!/bin/bash
# Scans model_invert slope from 0.44 to 0.56 (step 0.02).
# For each slope: patches PC_StepLadder_Correction.h, recompiles,
# runs all source runs, saves results under slope_scan/run<NNN>/slope_<val>.root

set -e
CORRECTION_H="../Armory/PC_StepLadder_Correction.h"
SCAN_DIR="slope_scan"

process_run() {
    local wrun=$(printf "%03d" "$1")
    local slope="$2"
    local outdir="${SCAN_DIR}/run${wrun}"
    mkdir -p "$outdir"
    local out="${outdir}/slope_${slope}.root"
    root -q -l -b -x "../../ANASEN_analysis/data/${DATASET}_Data/${PREFIX}${wrun}_mapped.root" \
         -e "tree->Process(\"TrackRecon.C+\", \"${out}\")" > /dev/null 2>&1
    [ -f "$out" ] && echo "  run $wrun slope $slope OK" || echo "  run $wrun slope $slope FAILED"
}
export -f process_run
export SCAN_DIR

for slope_x100 in 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 72 74 76 78 80 82 84 86 88 90; do
# for slope_x100 in 68 70 72 74 76 78 80 82 84 86 88 90; do
    slope=$(awk "BEGIN{printf \"%.2f\", $slope_x100/100}")
    echo "=== slope=${slope} ==="

    sed -i "s/double slope = [0-9.]*/double slope = ${slope}/" "$CORRECTION_H"

    echo "  Compiling..."
    root -q -l -b -e '.L TrackRecon.C++O' 2>/dev/null

    # 27Al alpha+gas runs (9, 12, 13)
    export DATASET="27Al" PREFIX="Run_"
    echo "  27Al runs 9 12 13..."
    export source_vertex=53.44; export timecut_low=400.0; process_run 12 "$slope"
    export source_vertex=-5.36; export timecut_low=12.0; export timecut_high=120.0; process_run 9 "$slope"
    unset timecut_low
    unset timecut_high

    # 17F alpha runs with per-run source vertices
    export DATASET="17F" PREFIX="SourceRun_"
    echo "  17F runs 18-21..."
    export source_vertex=53.44;  process_run 18 "$slope"
    export source_vertex=14.24;  process_run 19 "$slope"
    export source_vertex=-24.96; process_run 20 "$slope"
    export source_vertex=-73.96; process_run 21 "$slope"

    echo "  done"
done

sed -i "s/double slope = [0-9.]*/double slope = 0.60/" "$CORRECTION_H"
echo "Restored slope = 0.60"
echo ""
echo "All done. Now run:"
echo "  root -l 'scratch/plot_slope_scan.C'"