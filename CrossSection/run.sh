#!/bin/bash

analysisIs2D=0
the_cases="ESFUnregEn ESFregEn ESFFSR5minus ESFFSR5plus ESFPU5minus ESFPU5plus"
count=0

for cs in ${the_cases} ; do
    dim=$(( $analysisIs2D + 1 ))
    rootb calcCrossSection.C+\(${analysisIs2D},\"defaultAdHocRemote\",DYTools::NORMAL_RUN,DYTools::NO_SYST,DYTools::_cs_None,${count}\) | tee log-xsec-${dim}D-${cs}.out
    count=$(( $count + 1 ))
done
