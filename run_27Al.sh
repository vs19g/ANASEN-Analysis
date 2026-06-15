#!/bin/bash

export DATASET="27Al"
export PREFIX="Run_"
export OUT_DIR="Output_27Al"
export Gain=2

# Clean up previous runs
rm -f ${OUT_DIR}/results_run*.root output_27Al.root

echo "Pre-compiling TrackRecon.C safely on a single core..."
root -q -l -b -e '.L TrackRecon.C++O'

process_run() {
    local wrun=$(printf "%03d" "$1")
    local prefix="${PREFIX:-Run_}"
    local outdir="${OUT_DIR:-Output_default}"
    local infile="../ANASEN_analysis/data/${DATASET}_Data/${prefix}${wrun}_mapped.root"
    local out="${outdir}/results_run${wrun}.root"

    # Ensure the directory exists so ROOT doesn't fail silently
    mkdir -p "$outdir"

    root -q -l -b -x "$infile" \
         -e "tree->Process(\"TrackRecon.C+\", \"${out}\")" > /dev/null 2>&1

    if [ -f "$out" ]; then
        echo "Run $wrun completed successfully in $outdir."
    else
        echo "ERROR: Run $wrun failed to generate $out"
    fi
}

export -f process_run

echo "Starting parallel processing..."
parallel --bar -j 4 process_run ::: {50..52}

echo "Merging files..."
hadd -k -j 4 ${OUT_DIR}/output_27Al.root ${OUT_DIR}/results_run*.root

unset DATASET