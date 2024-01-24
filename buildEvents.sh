#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 runID timeWindow_ns"
    echo "Exiting..."
    exit 1
fi

runID=$1
timeWindow=$2

fileList=`\ls -1 ../*${runID}*.fsu`

./EventBuilderNoTrace ${timeWindow} 0 ${fileList}

mv -vf ../*${runID}*${timeWindow}_noTrace.root .