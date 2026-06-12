#!/bin/bash
# Scans dither_sigma from 0.0 to 10.0 mm.
# No recompilation needed per step; values are passed via environment variable!

set -e
SCAN_DIR="dither_scan"

echo "=== Compiling TrackRecon.C once... ==="
root -q -l -b -e '.L TrackRecon.C++O' 2>/dev/null

process_run() {
    local wrun=$(printf "%03d" "$1")
    local dither="$2"
    local outdir="${SCAN_DIR}/run${wrun}"
    
    mkdir -p "$outdir"
    local out="${outdir}/dither_${dither}.root"

    # Pass the variable to C++
    export DITHER_SIGMA="$dither"

    root -q -l -b -x "../../ANASEN_analysis/data/${DATASET}_Data/${PREFIX}${wrun}_mapped.root" \
         -e "tree->Process(\"TrackRecon.C+\", \"${out}\")" > /dev/null 2>&1
         
    [ -f "$out" ] && echo "  run $wrun dither $dither OK" || echo "  run $wrun dither $dither FAILED"
}

export -f process_run
export SCAN_DIR

# Loop through the dither amounts (in mm) you want to test
for dither in 0.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0; do
    echo "=== Scanning Dither Sigma = ${dither} mm ==="
    
    # 27Al alpha+gas runs (9, 12)
    export DATASET="27Al" PREFIX="Run_"
    echo "  27Al runs 9 12..."
    # FIX: Changed "$slope" to "$dither"
    export source_vertex=53.44; export timecut_low=400.0; process_run 12 "$dither"
    export source_vertex=-5.36; export timecut_low=12.0; export timecut_high=120.0; process_run 9 "$dither"
    unset timecut_low
    unset timecut_high

    # 17F alpha runs with per-run source vertices
    export DATASET="17F" PREFIX="SourceRun_"
    echo "  17F runs 18-21..."
    export source_vertex=53.44;  process_run 18 "$dither"
    export source_vertex=14.24;  process_run 19 "$dither"
    export source_vertex=-24.96; process_run 20 "$dither"
    export source_vertex=-73.96; process_run 21 "$dither"
done

echo "=== Dither scan complete! ==="