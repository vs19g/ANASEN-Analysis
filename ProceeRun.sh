#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 runID timeWindow_ns"
    echo "Exiting..."
    exit 1
fi

runID=$1
timeWindow=$2

rawFolder=/media/nvmeData/ANASEN_test
rootFolder=/media/nvmeData/ANASEN_test/root_data

fileList=`\ls -1 ${rawFolder}/*${runID}*.fsu`

./EventBuilderNoTrace ${timeWindow} 0 ${fileList}

outFile=${rawFolder}/*${runID}*${timeWindow}_noTrace.root

mv -vf ${outFile} ${rootFolder}/.

./Mapper ${rootFolder}/*${runID}*${timeWindow}_noTrace.root