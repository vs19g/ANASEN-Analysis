rm results_run*.root
export DATASET="27Al"
export flip180="0"
export flipa=0
export anode_offset=2
#declare -i run=28
#while [[ $run -lt 34 ]]; do #runs 1 to 84
#    wrun=$(printf "%03d" $run)
#    root -q -l -b -x ../ANASEN_analysis/data/27Al_Data/Run_"$wrun"_mapped.root -e 'tree->Process("Make#Vertex.C+O")'; mv Analyzer_SX3.root 27Al_output/results_run$wrun.root;
#    run=run+1
#done

declare -i run=50
while [[ $run -lt 51 ]]; do #runs 1 to 84
    wrun=$(printf "%03d" $run)
    root -q -l -b -x ../ANASEN_analysis/data/27Al_Data/Run_"$wrun"_mapped.root -e 'tree->Process("MakeVertex.C+O","Analyzer_27Al.root")'; mv Analyzer_27Al.root 27Al_output/results_run$wrun.root;
    run=run+1
done

rm output.root
hadd -k -j 4 output.root 27Al_output/results_run*.root
mv output.root output_27Al.root

#root -q -l -b -x -e '.L MakeVertex.C+O';
#halfproc=3
#parallel --ctag --bar -j $halfproc ./run.sh ::: {028..034} ::: 27Al #color-tag, linebuffer, then run run.sh in parallel

unset souce_vertex
unset DATASET
unset flip180
