#!/bin/bash
# run_misc_ede.sh -- run the standalone MISC E-dE telescope plots (MiscEdE.C)
# over raw, pre-mapping event-built files. Independent of the main TrackRecon
# pipeline: no mapping/calibration/reconstruction, works on beam-monitor runs.
#
# The macro reads the flat sn/ch/e event-built tree and produces per-channel 
# spectra + every pairwise E-dE 2D, so you can see all correlations and confirm 
# telescopes directly from the data.
#
# Point EVT_DIR / PREFIX / SUFFIX at wherever your event-built files live and
# how they're named, then list the run numbers. Each run is processed into its
# own MiscEdE_<stem>.root; set MERGE=1 to also hadd them together.
#
# Usage:
#   ./run_misc_ede.sh                 # uses the RUNS list below
#   ./run_misc_ede.sh 351 353 355     # override runs on the command line

set -u

# ---- configure these to your layout ----
EVT_DIR="${EVT_DIR:-../ANASEN_analysis/data/17F_Data}"   # where the event-built files are
PREFIX="${PREFIX:-Run_}"                                 # filename prefix before the run number
SUFFIX="${SUFFIX:-_2000.root}"                           # filename suffix after the run number
PAD="${PAD:-3}"                                          # zero-pad width for the run number (Run_053 -> 3)
OUT_DIR="${OUT_DIR:-Output_misc_ede}"
MERGE="${MERGE:-1}"                                      # 1 = also hadd all per-run outputs
MACRO="${MACRO:-MiscEdE.C}"
CUT_DIR="${CUT_DIR:-/home/vsitaraman/ANASEN_analysis}"   # Where 17FCut.root and 16OCut.root live
# ----------------------------------------

# Runs: command-line args win, else the lollipop/Si-monitor test runs below.
if [ "$#" -gt 0 ]; then
    RUNS=("$@")
else
    RUNS=( 52 53 54 55 63 65 68 69 72 74 75 76 78 92 168 170 226 262 324 327 337 334 369 370 371 372 373 374)
fi

mkdir -p "$OUT_DIR"
OUT_ABS="$(cd "$OUT_DIR" && pwd)"
MACRO_ABS="$(pwd)/$MACRO"
export OUT_ABS MACRO_ABS EVT_DIR PREFIX SUFFIX PAD CUT_DIR

built_path() {
    local run; run=$(printf "%0${PAD}d" "$1")
    echo "${EVT_DIR}/${PREFIX}${run}${SUFFIX}"
}
export -f built_path

# One run -> one MiscEdE_<stem>.root, runnable under GNU parallel.
misc_ede_one() {
    local r="$1"
    local infile; infile=$(built_path "$r")
    if [ ! -f "$infile" ]; then
        echo "SKIP: $infile not found"
        return
    fi
    local infile_abs; infile_abs="$(cd "$(dirname "$infile")" && pwd)/$(basename "$infile")"
    echo "=== MISC E-dE: run $r ($infile_abs) ==="
    
    # EXPORT METADATA FOR C++ SCRIPT
    export CURRENT_RUN="$r"
    
    # Macro names its own output MiscEdE_<stem>.root in the CWD, so run inside
    # OUT_DIR to keep outputs together. Absolute paths avoid relative-path
    # surprises from the cd.
    ( cd "$OUT_ABS" && root -l -b -q "${MACRO_ABS}(\"${infile_abs}\")" )
}
export -f misc_ede_one

rm -f "$OUT_DIR"/*.root

JOBS="${JOBS:-7}"   # match the -j 7 you use for run_once
if command -v parallel >/dev/null 2>&1; then
    parallel --bar -j "$JOBS" misc_ede_one ::: "${RUNS[@]}"
else
    echo "GNU parallel not found; running serially."
    for r in "${RUNS[@]}"; do misc_ede_one "$r"; done
fi

if [ "$MERGE" -eq 1 ]; then
    shopt -s nullglob 
    parts=()
    for f in "$OUT_DIR"/MiscEdE_*.root; do
        [ "$(basename "$f")" = "MiscEdE_ALL.root" ] && continue
        parts+=("$f")
    done
    if [ "${#parts[@]}" -gt 1 ]; then
        echo "=== merging ${#parts[@]} outputs -> ${OUT_DIR}/MiscEdE_ALL.root ==="
        hadd -f "${OUT_DIR}/MiscEdE_ALL.root" "${parts[@]}"
    fi
fi

echo "Done. Outputs in ${OUT_DIR}/"