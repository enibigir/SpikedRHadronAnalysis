#!/bin/bash

# CONTROL CENTRE ----------------------

mass=1800
events=1000
cmEnergy="13TeV"
prefix="Pre"
#genFilePath="test/EXO-RunIISummer20UL18GENSIM-00010_1_cfg_v3.py"

# -------------------------------------

# example dir_name: $dir_name
dir_name="M"$mass"_"$cmEnergy"_pythia8"
echo "All files will be saved in $dir_name"

if [! -f $genFilePath ]; then
    echo "Gen file not found"
    exit 0
fi

if [ ! -d "data/$dir_name" ]; then
    mkdir -p data/$dir_name
    echo "creating data/${dir_name}"
fi

# gen-sim output files
genSimRoot="data/"$prefix"2015_EXO-RunIISummer20UL18GENSIM-00010-v3_"$events"Events.root"
#genSimOut="gensimM"$mass".out"

# digi-L1-digi2ray output files
digiRawRoot=$prefix"2015digirawM"$mass"_"$events"Events.root"
digiRawOut=$prefix"2015digirawM"$mass"_"$events"Events.out"

# reco output files
recoRoot=$prefix"2015recoM"$mass"_"$events"Events.root"
recoOut=$prefix"2015recoM"$mass"_"$events"Events.out"

# miniAOD 
#miniAODRoot="miniAODM"$mass".root"
#miniAODOut="miniAODM"$mass".out"

#echo "Starting step 0: GEN-SIM"
#cmsDriver.py $genFilePath \
#    --fileout file:tests/$dir_name/$genSimRoot \
#    --mc \
#    --eventcontent RAWSIM \
#    --datatier GEN-SIM \
#    --conditions auto:run2_mc \
#    --step GEN,SIM \
#    --python_filename tests/$dir_name/step0_cfg.py \
#    --geometry DB:Extended \
#    --customise SimG4Core/CustomPhysics/Exotica_HSCP_SIM_cfi.customise,SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 \
#    --era Run2_2018 \
#    -n $events >& tests/$dir_name/$genSimOut
#echo "Step 0 completed"

echo "Starting step 1: DIGI-L1-DIGI2RAW"
cmsDriver.py --filein file:data/$genSimRoot \
    --fileout file:data/$dir_name/$digiRawRoot\
    --mc \
    --eventcontent RAWSIM \
    --datatier GEN-SIM-RAW \
    --conditions auto:run2_mc \
    --step DIGI,L1,DIGI2RAW \
    --python_filename data/$dir_name/step1_cfg.py \
    --geometry DB:Extended \
    --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 \
    --era Run2_2018 \
    -n -1 >& data/$dir_name/$digiRawOut
echo "Step 1 completed"

echo "Starting step 2: RAW2DIGI-L1Reco-RECO"
cmsDriver.py --filein file:data/$dir_name/$digiRawRoot \
    --fileout file:data/$dir_name/$recoRoot \
    --mc \
    --eventcontent FEVTDEBUGHLT \
    --datatier GEN-SIM-DIGI-AOD \
    --conditions auto:run2_mc \
    --step HLT:GRun,RAW2DIGI,L1Reco,RECO \
    --python_filename data/$dir_name/step2_cfg.py \
    --geometry DB:Extended \
    --era Run2_2018 \
    -n -1 >& data/$dir_name/$recoOut
echo "Step 2 completed"

#echo "Starting step 3: PAT"
#cmsDriver.py --filein file:tests/$dir_name/$recoRoot \
#   --fileout file:tests/$dir_name/$miniAODRoot \
#   --mc \
#   --eventcontent MINIAODSIM \
#   --datatier MINIAODSIM \
#   --conditions auto:run2_mc \
#   --step PAT \
#   --python_filename tests/$dir_name/step3_cfg.py \
#   --magField 38T_PostLS1 \
#   --geometry Extended2015 \
#   --era Run2_25ns \
#   --runUnscheduled \
#   -n -1 >& tests/$dir_name/$miniAODOut
#echo "Step 3 completed"