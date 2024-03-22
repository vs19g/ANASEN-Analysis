#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 runID timeWindow_ns"
    echo "Exiting..."
    exit 1
fi

runID=$1
timeWindow=$2

rawFolder=/media/nvmeData/ANASEN_test/analysis/data
rootFolder=/media/nvmeData/ANASEN_test/root_data

fileList=`\ls -1 ${rawFolder}/*${runID}*.fsu`

./EventBuilder ${timeWindow} 0 0 500000 ${fileList}

outFile=${rawFolder}/*${runID}*${timeWindow}.root

mv -vf ${outFile} ${rootFolder}/.

./Mapper ${rootFolder}/*${runID}*${timeWindow}.root
