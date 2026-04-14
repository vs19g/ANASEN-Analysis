#Alpha runs at different spacer positions
#rm results_run*.root
#root -q -l -x ../ANASEN_analysis/data/27Al_Data/Run_009_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run09.root;
#root -q -l -x ../ANASEN_analysis/data/27Al_Data/Run_001_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run01.root;
#root -q -l -x ../ANASEN_analysis/data/27Al_Data/Run_002_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run02.root;
#root -q -l -x ../ANASEN_analysis/data/27Al_Data/Run_003_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run03.root;
#root -q -l -x ../ANASEN_analysis/data/27Al_Data/Run_004_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run04.root;

#alpha+gas 27Al
export DATASET="26Al"
export flip180="0"
#root -q -b -x ../ANASEN_analysis/data/27Al_Data/Run_009_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run09.root;
#root -q -b -x ../ANASEN_analysis/data/27Al_Data/Run_012_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run12.root;
#root -q -b -x ../ANASEN_analysis/data/27Al_Data/Run_013_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run13.root;
#exit

#protons+gas, 27Al
export flip180="1"
#root -q -b -x ../ANASEN_analysis/data/27Al_Data/Run_015_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run15.root;
#root -q -b -x ../ANASEN_analysis/data/27Al_Data/Run_017_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run17.root;
#root -q -b -x ../ANASEN_analysis/data/27Al_Data/Run_018_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run18.root;
#root -q -b -x ../ANASEN_analysis/data/27Al_Data/Run_019_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run19.root;
#root -q -b -x ../ANASEN_analysis/data/27Al_Data/Run_020_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run20.root;
#root -q -b -x ../ANASEN_analysis/data/27Al_Data/Run_021_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run21.root;
#root -q -b -x ../ANASEN_analysis/data/27Al_Data/Run_022_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run22.root;
#exit

#27Al reaction data
#root -b -q -l -x ../ANASEN_analysis/data/27Al_Data/Run_051_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run51.root;
#root -b -q -l -x ../ANASEN_analysis/data/27Al_Data/Run_078_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run78.root;
#root -b -q -l -x ../ANASEN_analysis/data/27Al_Data/Run_081_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run81.root;

#root -l -x  results_run19.root results_run12.root -e "new TBrowser"
#exit
export DATASET="17F"
export flip180="0"
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
#export source_vertex=-57.28; root -q -l -b -x ../ANASEN_analysis/data/17F_Data/ProtonRun_035_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run35.root;
export source_vertex=-8.28; root -q -l -b -x ../ANASEN_analysis/data/17F_Data/ProtonRun_036_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root resulrs_run36.root;
export source_vertex=-27.88; root -q -l -b -x ../ANASEN_analysis/data/17F_Data/ProtonRun_037_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run37.root;
export source_vertex=11.32; root -q -l -b -x ../ANASEN_analysis/data/17F_Data/ProtonRun_038_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run38.root;
export source_vertex=30.92; root -q -l -b -x ../ANASEN_analysis/data/17F_Data/ProtonRun_039_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run39.root;
export source_vertex=50.52; root -q -l -b -x ../ANASEN_analysis/data/17F_Data/ProtonRun_041_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run41.root;
export source_vertex=70.12; root -q -l -b -x ../ANASEN_analysis/data/17F_Data/ProtonRun_042_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run42.root;
export source_vertex=109.32; root -q -l -b -x ../ANASEN_analysis/data/17F_Data/ProtonRun_043_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run43.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/ProtonRun_043_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run43.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Run_099_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run99.root;
#root -q -l -b -x ../ANASEN_analysis/data/17F_Data/Run_104_mapped.root -e 'tree->Process("MakeVertex.C+O")'; mv Analyzer_SX3.root results_run104.root;
#mv Analyzer_SX3.root results_run19.root;

unset souce_vertex
unset DATASET
unset flip180
