#!/bin/bash

export DATASET="27Al"
export PREFIX="Run_"
export OUT_DIR="Output_27Al"
export reactiondata=1
export CO2percent=3
export pressure_in_torr=250
export CATHODE_GAIN=1
export BEAM_AXIS_X=-15
export BEAM_AXIS_Y=-5

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
unset reactiondata
unset CO2percent
unset pressure_in_torr
unset CATHODE_GAIN
unset A1C1_LOWBAND_RFACTOR
unset A1C1_Z_SCALE_QQQ
unset A1C1_Z_OFF_QQQ
unset A1C1_Z_OFF_SX3
unset BEAM_AXIS_X
unset BEAM_AXIS_Y