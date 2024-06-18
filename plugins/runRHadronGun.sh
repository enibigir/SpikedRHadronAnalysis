#!/bin/bash

for file in data/IndividualRHadronProcesses/*.txt; do
    name=${file##*/}
    base=${name%.txt}
    echo "Running RHadron $base"
    cmsRun test/RHadronGun_cfg.py inputFiles=data/IndividualRHadronProcesses/$base.txt outputFile=data/IndividualRHadronRoots/$base.root
done