#!/bin/bash
#parallel -j 6 echo ./ProcessRun.sh {1} 2000 0 ::: {020..400}

#parallel --results log/log_{}.txt --ctag -j 6 ./ProcessRun.sh {1} 2000 0 ::: {020..400} # for 17F

parallel --results log/log_{}.txt --ctag -j 6 ./ProcessRun.sh {1} 2000 0 ::: {001..021}  
