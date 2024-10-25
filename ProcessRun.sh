#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 runID timeWindow_ns option"
    echo "Exiting..."
    exit 1
fi

runID=$1
timeWindow=$2

option=$3

rawFolder=/home/tandem/data1/2024_09_17Fap/data
rootFolder=/home/tandem/data1/2024_09_17Fap/data/root_data

if [ $option -eq 0 ]; then

    rsync -auh --info=progress2 splitpole@128.186.111.223:/media/nvmeData/2024_09_17Fap/*.fsu /home/tandem/data1/2024_09_17Fap/data

    fileList=`\ls -1 ${rawFolder}/*Run_${runID}_*.fsu`

    ./EventBuilder ${timeWindow} 0 0 10000000 ${fileList}

    outFile=${rawFolder}/*${runID}*${timeWindow}.root

    mv -vf ${outFile} ${rootFolder}/.

    ./Mapper ${rootFolder}/*${runID}*${timeWindow}.root
fi

root "processRun.C(\"${rootFolder}/Run_${runID}_mapped.root\")"
