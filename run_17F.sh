rm results_run*.root
export DATASET="17F"
export reactiondata=1

#17F reaction data
#declare -i run=231 #49
#while [[ $run -lt 258 ]]; do #392
#    wrun=$(printf "%03d" $run)
#    file_exists=$(test -f ../ANASEN_analysis/data/17F_Data/Run_"$wrun"_mapped.root)
#    if [[ $file_exists -ne 0 ]]; then continue; fi
#    root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Run_"$wrun"_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run$wrun.root;
#    run=run+1
#done

function run_once() {
    wrun=$(printf "%03d" $1)
    file_exists=$(test -f ../ANASEN_analysis/data/17F_Data/Run_"$wrun"_mapped.root)
    if [[ $file_exists -ne 0 ]]; then continue; fi
    root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Run_"$wrun"_mapped.root -e 'tree->Process("TrackRecon.C+O","Analyzer_17F.root")'; mv Analyzer_17F.root 17F_output/results_run$wrun.root;
}

export -f run_once
# run_once 350
# parallel -j 6 --ctag run_once {1} ::: {325..400}
parallel -j 6 --ctag run_once {1} ::: {351,353,355,358,359,360,362,367}
rm output_17F.root
hadd -j 4 -k output_17F.root 17F_output/results_run3*.root

unset souce_vertex
unset DATASET
unset flip180
unset flipa
unset anode_offset
unset reactiondata
