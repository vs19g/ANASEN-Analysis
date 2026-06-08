#!/bin/bash

export DATASET="27Al"

# Clean up previous runs
rm -f 27Al_output/results_run*.root output_27Al.root

echo "Pre-compiling TrackRecon.C safely on a single core..."
root -q -l -b -e '.L TrackRecon.C++O'

process_run() {
    local wrun=$(printf "%03d" $1)
    local out="Output_27Al/results_run${wrun}.root"

    root -q -l -b -x "../ANASEN_analysis/data/27Al_Data/Run_${wrun}_mapped.root" \
         -e "tree->Process(\"TrackRecon.C+\", \"${out}\")" > /dev/null 2>&1
         
    if [ -f "$out" ]; then
        echo "Run $wrun completed successfully."
    else
        echo "ERROR: Run $wrun failed to generate $out"
    fi
}

export -f process_run

echo "Starting parallel processing..."
parallel --bar -j 4 process_run ::: {50..52}

echo "Merging files..."
hadd -k -j 4 output_27Al.root Output_27Al/results_run*.root

unset DATASET