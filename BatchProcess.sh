#!/bin/bash
#parallel -j 6 echo ./ProcessRun.sh {1} 2000 0 ::: {020..400}

parallel  --jobs 1 --results log/log_{}.txt --memfree 1G --ctag -j 1 ./ProcessRun.sh {1} 2000 0 ::: {109..400} # for 17F

# parallel --results log/log_{}.txt --ctag -j 6 ./ProcessRun.sh {1} 4000 0 ::: {5..21}  
