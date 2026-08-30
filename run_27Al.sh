#!/bin/bash

export DATASET="27Al"
export PREFIX="Run_"
export OUT_DIR="Output_27Al"
export reactiondata=1
export CO2percent=3
export pressure_in_torr=250
export CATHODE_GAIN=3.0
export source_vertex=-200.0
export DEDX_SCALE=0.89
export CUTLIST=cuts_list.txt
export BEAM_AXIS_X=0
export BEAM_AXIS_Y=0

echo "Pre-compiling TrackRecon.C safely on a single core..."
root -q -l -b -e '.L TrackRecon.C++O'

process_run() {
    local wrun=$(printf "%03d" "$1")
    local prefix="${PREFIX:-Run_}"
    local infile="../ANASEN_analysis/data/${DATASET}_Data/${prefix}${wrun}_mapped.root"
    
    # Dynamically point to the correct output directory for this X/Y iteration
    # local current_out_dir="Output_27Al_X${BEAM_AXIS_X}_Y${BEAM_AXIS_Y}"
    local current_out_dir="Output_27Al"
    local out="${current_out_dir}/results_run${wrun}.root"

    root -q -l -b -x "$infile" \
         -e "tree->Process(\"TrackRecon.C+\", \"${out}\")" > /dev/null 2>&1

    if [ -f "$out" ]; then
        echo "Run $wrun completed successfully in ${current_out_dir}."
    else
        echo "ERROR: Run $wrun failed to generate $out"
    fi
}

export -f process_run

# for x in -5 5
# do 
#     BEAM_AXIS_X=$x  
#     for y in -5 5 
#     do 
#         BEAM_AXIS_Y=$y  

        # Define and create a clean directory name BEFORE running parallel tasks
        # CURRENT_OUT_DIR="Output_27Al_X${BEAM_AXIS_X}_Y${BEAM_AXIS_Y}"
        CURRENT_OUT_DIR="Output_27Al"
        mkdir -p "$CURRENT_OUT_DIR"

        echo "Running Eloss.py with a scaling parameter of $DEDX_SCALE"
        echo "running with a beam offset of $BEAM_AXIS_X $BEAM_AXIS_Y"
        python3 eloss_calculations/Eloss.py
        
        echo "Starting parallel processing..."
        time parallel --bar -j 10 process_run ::: {24..41} 
        time parallel --bar -j 3 process_run ::: 44 45 46
        time parallel --bar -j 8 process_run ::: {50..59}
        # time parallel --bar -j 4 process_run ::: 62 63 66 67 68
        # time parallel --bar -j 1 process_run ::: 73
        # time parallel --bar -j 1 process_run ::: 74
        # time parallel --bar -j 4 process_run ::: {78..89}

        echo "Merging files..."
        # Fixed: Safely merge using the clean directory variable (added -f to overwrite if re-running)
        hadd -k -f -j 4 "${CURRENT_OUT_DIR}/output_27Al.root" "${CURRENT_OUT_DIR}/results_run"*.root
#     done
# done

# Cleanup
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
unset CUTLIST
unset DEDX_SCALE
unset CURRENT_OUT_DIR
echo "Script execution finished."