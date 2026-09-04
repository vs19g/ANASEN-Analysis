#!/bin/bash
export DATASET="17F"
export reactiondata=1
export CO2percent=4
export pressure_in_torr=250

run_once() {
    local wrun=$(printf "%03d" "$1")
    local prefix="${PREFIX:-Run_}"
    local outdir="${OUT_DIR:-Output_default}"
    local infile="../ANASEN_analysis/data/${DATASET}_Data/${prefix}${wrun}_mapped.root"
    local out="${outdir}/results_run${wrun}.root"

    # Skip if the input file doesn't exist (return, not continue — this is a function)
    if [ ! -f "$infile" ]; then
        echo "SKIP: input $infile not found"
        return
    fi

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

export -f run_once

export DATASET="17F"
export PREFIX="Run_"
export OUT_DIR="Output_17F"
export source_vertex=-200.0
export CATHODE_GAIN=1
# A1C1 cfrac low-band fold + z autocal (see TrackRecon.C Begin()). Defaults are the
# 17F values; override here to re-tune without recompiling.
export A1C1_LOWBAND_RFACTOR=7.0
export CUTLIST=cuts_list.txt
rm -f ${OUT_DIR}/*.root

# Pre-compile TrackRecon.C ONCE on a single core so parallel jobs don't race on ACLiC
echo "Pre-compiling TrackRecon.C..."
root -q -l -b -e '.L TrackRecon.C++O'

# 3% CO2
# parallel --bar -j 6 run_once ::: {325..400}
parallel --bar -j 11  run_once ::: {351..400}
# parallel --bar -j 7 run_once ::: 351 353 355 358 359 360 362 367
hadd -j 4 -k ${OUT_DIR}/Output_17F.root ${OUT_DIR}/results_run*.root

unset source_vertex
unset DATASET
unset reactiondata
unset CO2percent
unset pressure_in_torr
unset CATHODE_GAIN
unset A1C1_LOWBAND_RFACTOR
unset CUTLIST
unset A1C1_Z_SCALE_QQQ
unset A1C1_Z_OFF_QQQ
unset A1C1_Z_OFF_SX3