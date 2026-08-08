#!/bin/bash

export DATASET="27Al"
export PREFIX="Run_"
export OUT_DIR="Output_27Al"
export reactiondata=1
export CO2percent=3
export pressure_in_torr=250
export CATHODE_GAIN=3.0
export source_vertex=-200.0
export BEAM_AXIS_X=-15
export BEAM_AXIS_Y=-5

# Clean up previous runs
rm -f Output_27Al/*.root

echo "Pre-compiling TrackRecon.C safely on a single core..."
root -q -l -b -e '.L TrackRecon.C++O'

process_run() {
    local wrun=$(printf "%03d" "$1")
    local prefix="${PREFIX:-Run_}"
    local infile="../ANASEN_analysis/data/${DATASET}_Data/${prefix}${wrun}_mapped.root"
    local out="Output_27Al/results_run${wrun}.root"

    mkdir -p Output_27Al

    root -q -l -b -x "$infile" \
         -e "tree->Process(\"TrackRecon.C+\", \"${out}\")" > /dev/null 2>&1

    if [ -f "$out" ]; then
        echo "Run $wrun completed successfully in Output_27Al."
    else
        echo "ERROR: Run $wrun failed to generate $out"
    fi
}

export -f process_run

echo "Starting parallel processing..."
time parallel --bar -j 6 process_run ::: {24..41} 
time parallel --bar -j 3 process_run ::: 44 45 46
# time parallel --bar -j 8 process_run ::: {50..59}
# time parallel --bar -j 4 process_run ::: 62 63 66 67 68
# time parallel --bar -j 1 process_run ::: 73
# time parallel --bar -j 1 process_run ::: 74
# time parallel --bar -j 4 process_run ::: {78..89}

echo "Merging files..."
hadd -k -j 4 Output_27Al/output_27Al.root Output_27Al/results_run*.root

# rootbrowse Output_27Al/output_27Al.root

unset DATASET
unset PREFIX
unset OUT_DIR
unset reactiondata
unset CO2percent
unset pressure_in_torr
unset CATHODE_GAIN
unset source_vertex
unset A1C1_LOWBAND_RFACTOR
unset A1C1_Z_SCALE_QQQ
unset A1C1_Z_OFF_QQQ
unset A1C1_Z_OFF_SX3
unset BEAM_AXIS_X
unset BEAM_AXIS_Y
echo "Script execution finished."