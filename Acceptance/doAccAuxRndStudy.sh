#!/bin/bash

analysisIs2D=$1
config=$2
runMode=$3

if [ ${#config} -eq 0 ] ; then
    echo "doAccAuxRndStudy.sh analysisIs2D configFile [runMode]"
    exit
fi

if [ ${#runMode} -eq 0 ] ; then
    runMode=DYTools::NORMAL_RUN
fi

runModeStr=
if [ "${runMode}" == "DYTools::DEBUG_RUN" ] || \
    [ "${runMode}" == "DEBUG_RUN" ] ; then
    runModeStr="debugRun-"
    runMode=DYTools::DEBUG_RUN
fi

seedMin=1001
seedMax=1020

systMode=DYTools::FSR_RND_STUDY
rndStudyStr="FSR_RND_STUDY"

iSeed=$seedMin

while [ ${iSeed} -le ${seedMax} ] ; do
    root -l -q -b plotDYAcceptance.C+\(${analysisIs2D},\"${config}\",${runMode},$systMode,\"${rndStudyStr}${iSeed}\"\) | tee log-accRndStudy-$((${analysisIs2D}+1))D-${runModeStr}-${rndStudyStr}${iSeed}.out
    iSeed=$(( $iSeed + 1 ))
done
