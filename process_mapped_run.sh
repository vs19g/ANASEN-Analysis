#!/bin/bash

# ==========================================
# CONFIGURATION
# ==========================================
DATA_DIR="/mnt/d/Remapped_files/27Al_data/root_data"
MACRO="TrackRecon.C"

# SAFETY SETTINGS
JOBS=2          # Keep low (2-4) to prevent WSL crashes
MIN_MEM="1G"    # Wait if RAM is full
# ==========================================

# 1. CHECK ARGUMENTS
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <Start_Run> <End_Run>"
    echo "Example: $0 10 50"
    exit 1
fi

START_RUN=$1
END_RUN=$2

# 2. COMPILE MACRO
# Compiling once is mandatory for parallel execution
echo "Compiling ${MACRO}..."
root -l -b -q -e "gROOT->ProcessLine(\".L ${MACRO}+\");"
if [ $? -ne 0 ]; then
    echo "Error: Compilation failed."
    exit 1
fi

# 3. DEFINE WORKER FUNCTION
run_job() {
    file_path="$1"
    macro_name="$2"
    
    echo "Processing: $file_path"
    
    # Execute ROOT
    nice -n 15 root -l -x -b -q "$file_path" -e "tree->Process(\"${macro_name}+\");" > "${file_path}.log" 2>&1
}
export -f run_job

# 4. QUEUE BUILDER (The "Skip" Logic is here)
echo "Checking runs $START_RUN to $END_RUN..."

for (( i=$START_RUN; i<=$END_RUN; i++ ))
do
    # Construct the input filename
    # Logic: Run_0 + number -> Run_0115
    file="${DATA_DIR}/Run_0${i}_mapped.root"
    
    # ------------------------------------------------------------
    # SKIP LOGIC
    # We check if the log file exists. If so, we assume it's done.
    # ------------------------------------------------------------
    log_file="${file}.log"
    
    if [ -f "$log_file" ]; then
        # >&2 redirects to stderr so it doesn't get fed into 'parallel'
        echo "Skipping Run $i: Log file already exists." >&2
        continue
    fi

    # Only add to queue if the INPUT file actually exists
    if [ -f "$file" ]; then
        echo "$file"
    else
        echo "Warning: Input file for run $i not found." >&2
    fi

done | parallel --jobs $JOBS --memfree $MIN_MEM --retries 2 run_job {} "$MACRO"