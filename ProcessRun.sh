#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 runID timeWindow_ns"
    echo "Exiting..."
    exit 1
fi

runID=$1
timeWindow=$2

rawFolder=/home/tandem/Desktop/analysis/data
rootFolder=/home/tandem/Desktop/analysis/data/root_data

rsync -a splitpole@128.186.111.223:/media/nvmeData/ANASEN27Alap/*.fsu /home/tandem/Desktop/analysis/data

fileList=`\ls -1 ${rawFolder}/Run_${runID}_*.fsu`

./EventBuilder ${timeWindow} 0 0 10000000 ${fileList}

outFile=${rawFolder}/*${runID}*${timeWindow}.root

mv -vf ${outFile} ${rootFolder}/.

./Mapper ${rootFolder}/*${runID}*${timeWindow}.root

root "processRun.C(\"${rootFolder}/Run_${runID}_mapped.root\")"
