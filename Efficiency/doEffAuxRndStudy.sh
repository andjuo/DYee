#!/bin/bash

analysisIs2D=$1
config=$2
systMode=$3
runMode=$4

if [ ${#systMode} -eq 0 ] ; then
    echo "doEffAuxRndStudy.sh analysisIs2D configFile systMode [runMode]"
    exit
fi

seedMin=1001
seedMax=1020

rndStudyStr=
if [ "${systMode}" == "DYTools::FSR_RND_STUDY" ] || \
    [ "${systMode}" == "FSR_RND_STUDY" ] ; then
  rndStudyStr="FSR_RND_STUDY"
  systMode=DYTools::FSR_RND_STUDY
elif [ "${systMode}" == "DYTools::PU_RND_STUDY" ] || \
    [ "${systMode}" == "PU_RND_STUDY" ] ; then
  rndStudyStr="PU_RND_STUDY"
  systMode=DYTools::PU_RND_STUDY
  seedMin=$(($seedMin-1000))
  seedMax=$(($seedMax-1000))
else
  echo "systMode=${systMode} is not expected"
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

iSeed=$seedMin

while [ ${iSeed} -le ${seedMax} ] ; do
    root -l -q -b plotDYEfficiency.C+\(${analysisIs2D},\"${config}\",${runMode},$systMode,\"${rndStudyStr}${iSeed}\"\) | tee log-effRndStudy-$((${analysisIs2D}+1))D-${runModeStr}-${rndStudyStr}${iSeed}.out
    iSeed=$(( $iSeed + 1 ))
done
