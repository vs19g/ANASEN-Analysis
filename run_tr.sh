#!/bin/bash

# ==========================================
# ANASEN Data Processing Script
# ==========================================

echo "Pre-compiling TrackRecon.C safely on a single core..."
root -q -l -b -e '.L TrackRecon.C++O'

# Function to process a single run
process_run() {
    local wrun=$(printf "%03d" "$1")
    # Use PREFIX if set, otherwise default to "Run_"
    local prefix="${PREFIX:-Run_}" 
    local out="source_run${wrun}.root"

    root -q -l -b -x "../ANASEN_analysis/data/${DATASET}/${prefix}${wrun}_mapped.root" \
         -e "tree->Process(\"TrackRecon.C+\", \"${out}\")" > /dev/null 2>&1
         
    if [ -f "$out" ]; then
        echo "Run $wrun completed successfully."
    else
        echo "ERROR: Run $wrun failed to generate $out"
    fi
}
export -f process_run

# --- Block 1: 27Al Source Runs No Gas (1-8) ---
if [[ 1 -eq 0 ]]; then
    export DATASET="27Al_Data"
    export PREFIX="Run_"
    export CO2percent=3.0
    
    echo "Starting parallel processing for 27Al runs..."
    rm -f p_output/all.root
    
    parallel --bar -j 6 process_run ::: {1..8}
fi

# --- Block 2: 27Al Alpha+Gas Runs (9, 12) ---
if [[ 1 -eq 0 ]]; then
    export DATASET="27Al_Data"
    export PREFIX="Run_"
    export timecut_low=400.0
    export source_vertex=53.44
    export CO2percent=3.0
    
    echo "Processing 27Al alpha+gas runs..."
    parallel --bar -j 6 process_run ::: 9 12
    
    unset timecut_low
fi

# --- Block 3: 27Al Protons+Gas Runs (15, 17-22) ---
if [[ 1 -eq 0 ]]; then
    export DATASET="27Al_Data"
    export PREFIX="Run_"
    export source_vertex=-200.0 # Source on the entrance window
    export CO2percent=3.0
    echo "Starting parallel processing for 27Al proton runs..."
    rm -f p_output/all.root
    
    parallel --bar -j 6 process_run ::: 15 {17..22}
    
    # Uncomment to merge via hadd when ready
    # hadd -j 4 -k p_output/all.root p_output/results_run*.root
    exit
fi

# --- Block 4: 17F Source Runs (5-14) ---
if [[ 1 -eq 0 ]]; then
    export DATASET="17F_Data"
    export PREFIX="Source_"
    echo "Starting parallel processing for 17F source runs..."
    rm -f p_output/all.root
    
    parallel --bar -j 6 process_run ::: {5..13}
fi

# --- Block 5: 17F Alpha Run with Gas (18-21) ---
if [[ 1 -eq 0 ]]; then
    export DATASET="17F_Data"
    export PREFIX="SourceRun_"
    export CO2percent=3.0
    echo "Processing 17F alpha runs with dynamic source vertices..."
    
    # Running sequentially since the source_vertex variable changes per run
    export source_vertex=53.44;  process_run 18
    export source_vertex=14.24;  process_run 19
    export source_vertex=-24.96; process_run 20
    export source_vertex=-73.96; process_run 21
fi

# --- Block 6: 17F Proton Data  ---
if [[ 1 -eq 1 ]]; then
    export DATASET="17F_Data"
    export PREFIX="ProtonRun_"
    
    export source_vertex=-57.28
    export CO2percent=4.0
    # parallel --bar -j 6 process_run ::: {38..42} #4% CO2
    export CO2percent=3.0
    parallel --bar -j 6 process_run ::: {44..48} #3CO2 
    # hadd -j 4 -k p_output/all.root p_output/results_run*.root
    exit

fi

# ==========================================
# Cleanup Environment Variables
# ==========================================
unset flipa
unset flipc
unset anode_offset
unset cathode_offset
unset source_vertex  # Fixed typo here
unset DATASET
unset PREFIX
unset flip180
unset timecut_low
unset timecut_high
unset CO2percent

echo "Script execution finished."