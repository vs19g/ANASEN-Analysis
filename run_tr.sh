#!/bin/bash

# ==========================================
# ANASEN Data Processing Script
# ==========================================

echo "Pre-compiling TrackRecon.C safely on a single core..."
root -q -l -b -e '.L TrackRecon.C++O'

# Function to process a single run
process_run() {
    local wrun=$(printf "%03d" "$1")
    local prefix="${PREFIX:-Run_}" 
    local outdir="${OUT_DIR:-Output_default}" 
    local out="${outdir}/results_run${wrun}.root"

    # Ensure the directory exists so ROOT doesn't fail silently
    mkdir -p "$outdir"

    root -q -l -b -x "../ANASEN_analysis/data/${DATASET}_Data/${prefix}${wrun}_mapped.root" \
         -e "tree->Process(\"TrackRecon.C+\", \"${out}\")" > /dev/null 2>&1
         
    if [ -f "$out" ]; then
        echo "Run $wrun completed successfully in $outdir."
    else
        echo "ERROR: Run $wrun failed to generate $out"
    fi
}
export -f process_run

export CO2percent=3

if [[ 1 -eq 0 ]]; then
    export DATASET=27Al Gain=2 a1c1_calibrate=1
    # root -q -l -b -e '
    #     TChain *t = new TChain("tree");
    #     t->Add("../ANASEN_analysis/data/27Al_Data/Run_015_mapped.root");
    #     t->Add("../ANASEN_analysis/data/27Al_Data/Run_017_mapped.root");
    #     t->Add("../ANASEN_analysis/data/27Al_Data/Run_018_mapped.root");
    #     t->Add("../ANASEN_analysis/data/27Al_Data/Run_019_mapped.root");
    #     t->Add("../ANASEN_analysis/data/27Al_Data/Run_020_mapped.root");
    #     t->Add("../ANASEN_analysis/data/27Al_Data/Run_021_mapped.root");
    #     t->Add("../ANASEN_analysis/data/27Al_Data/Run_022_mapped.root");
    #     t->Process("TrackRecon.C+", "Output_cal/calib_27Al.root");
        # '
    export DATASET=17F Gain=1 
    root -q -l -b -e '
        TChain *t = new TChain("tree");
        // proton runs, 3% CO2
        t->Add("../ANASEN_analysis/data/17F_Data/ProtonRun_044_mapped.root");
        t->Add("../ANASEN_analysis/data/17F_Data/ProtonRun_045_mapped.root");
        t->Add("../ANASEN_analysis/data/17F_Data/ProtonRun_046_mapped.root");
        t->Add("../ANASEN_analysis/data/17F_Data/ProtonRun_047_mapped.root");
        t->Add("../ANASEN_analysis/data/17F_Data/ProtonRun_048_mapped.root");
        // proton runs, 4% CO2
        t->Add("../ANASEN_analysis/data/17F_Data/ProtonRun_038_mapped.root");
        t->Add("../ANASEN_analysis/data/17F_Data/ProtonRun_039_mapped.root");
        t->Add("../ANASEN_analysis/data/17F_Data/ProtonRun_040_mapped.root");
        t->Add("../ANASEN_analysis/data/17F_Data/ProtonRun_041_mapped.root");
        t->Add("../ANASEN_analysis/data/17F_Data/ProtonRun_042_mapped.root");
        // source runs (extra stats; same Gain=1)
        t->Add("../ANASEN_analysis/data/17F_Data/SourceRun_018_mapped.root");
        t->Add("../ANASEN_analysis/data/17F_Data/SourceRun_019_mapped.root");
        t->Add("../ANASEN_analysis/data/17F_Data/SourceRun_020_mapped.root");
        t->Add("../ANASEN_analysis/data/17F_Data/SourceRun_021_mapped.root");
        t->Process("TrackRecon.C+", "Output_cal/calib_17F.root");
        '
    unset a1c1_calibrate
    unset DATASET
    unset Gain
    exit
fi

# --- Block 1: 27Al Source Runs No Gas (1-8) ---
if [[ 1 -eq 0 ]]; then
    export DATASET="27Al"
    export PREFIX="Run_"
    export OUT_DIR="Output_av"
    echo "Starting parallel processing for 27Al runs..."
    rm -f ${OUT_DIR}/all.root
    
    parallel --bar -j 6 process_run ::: {1..8}
fi

# --- Block 2: 27Al Alpha+Gas Runs (9, 12) ---
if [[ 1 -eq 1 ]]; then
    export DATASET="27Al"
    export PREFIX="Run_"
    export OUT_DIR="Output_a"
    export Gain=2
    echo "Processing 27Al alpha+gas runs..."
    export source_vertex=-5.36; export timecut_low=12.0; export timecut_high=119.0; process_run 9 "$slope"
    unset timecut_high
    export source_vertex=53.44; export timecut_low=400.0; process_run 12 "$slope"
    unset Gain
    unset timecut_low
fi

# --- Block 3: 27Al Protons+Gas Runs (15, 17-22) ---
if [[ 1 -eq 0 ]]; then
    export DATASET="27Al"
    export PREFIX="Run_"
    export OUT_DIR="Output_p"       
    export Gain=2
 
    rm -f ${OUT_DIR}/Al_protons.root
    export source_vertex=-200.0 # Source on the entrance window
    echo "Starting parallel processing for 27Al proton runs..."

    # process_run 17
    parallel --bar -j 6 process_run ::: 15 {17..22}
    hadd -j 4 -k ${OUT_DIR}/Al_protons.root ${OUT_DIR}/results_run0{15..22}.root
    unset Gain

fi

# --- Block 4: 17F Source Runs (5-14) ---
if [[ 1 -eq 0 ]]; then
    export DATASET="17F"
    export PREFIX="Source_"
    export OUT_DIR="Output_av"
    echo "Starting parallel processing for 17F source runs..."
    
    parallel --bar -j 6 process_run ::: {5..13}
fi

# --- Block 5: 17F Alpha Run with Gas (18-21) ---
if [[ 1 -eq 1 ]]; then
    export DATASET="17F"
    export PREFIX="SourceRun_"
    export OUT_DIR="Output_a"
    echo "Processing 17F alpha runs with dynamic source vertices..."
    
    # Running sequentially since the source_vertex variable changes per run
    export source_vertex=53.44;  process_run 18
    export source_vertex=14.24;  process_run 19
    export source_vertex=-24.96; process_run 20
    export source_vertex=-73.96; process_run 21
fi

# --- Block 6: 17F Proton Data  ---
if [[ 1 -eq 1 ]]; then
    export DATASET="17F"
    export PREFIX="ProtonRun_"
    export OUT_DIR="Output_p"
    rm -f ${OUT_DIR}/*pc*.root 
    export source_vertex=-200.0
    
    parallel --bar -j 6 process_run ::: {44..48} #3% CO2 
    hadd -j 4 -k ${OUT_DIR}/3pc.root ${OUT_DIR}/results_run0{44..48}.root
    
    export CO2percent=4
    parallel --bar -j 6 process_run ::: {38..42} #4% CO2
    hadd -j 4 -k ${OUT_DIR}/4pc.root ${OUT_DIR}/results_run0{38..42}.root
fi

# ==========================================
# Cleanup Environment Variables
# ==========================================
unset source_vertex  
unset DATASET
unset OUT_DIR
unset PREFIX
unset CO2percent
unset timecut_low
unset timecut_high

echo "Script execution finished."