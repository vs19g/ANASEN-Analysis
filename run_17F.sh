rm results_run*.root
export DATASET="17F"
export flip180="0"
export flipa=0
export anode_offset=1
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Source_005_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run05.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Source_006_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run06.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Source_007_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run07.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Source_008_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run08.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Source_009_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run09.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Source_010_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run10.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Source_011_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run11.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Source_012_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run12.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Source_013_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run13.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Source_014_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run14.root;

#17F pulser runs
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/PulserRun_015_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run15.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/PulserRun_016_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run16.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/PulserRun_017_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run17.root;

#17F alpha run with gas
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/SourceRun_018_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run18.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/SourceRun_019_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run19.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/SourceRun_020_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run20.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/SourceRun_021_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run21.root;

#17F reaction data
export flip180="1"
declare -i run=231 #49
while [[ $run -lt 235 ]]; do #392
    wrun=$(printf "%03d" $run)
    file_exists=$(test -f ../ANASEN_analysis/data/17F_Data/Run_"$wrun"_mapped.root)
    if [[ $file_exists -ne 0 ]]; then continue; fi
    root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Run_"$wrun"_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run$wrun.root;
    run=run+1
done

rm output_17F.root
hadd -j 4 -k output_17F.root results_run2*.root

unset souce_vertex
unset DATASET
unset flip180
unset flipa
unset anode_offset
