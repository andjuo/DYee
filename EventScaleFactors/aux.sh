#!/bin/bash

filename="../config_files/fall8TeV-vilnius.input"
filename="../config_files/data_vilnius8TeV_regSSD.conf.py"
triggerSet="Full2012_hltEffOld"

analysisIs2D=0
anTag="-$((${analysisIs2D}+1))D"

tnpFullRun_eval=1

tnpFullRun_recalc="data11111mc11111sf0"
tnpFullRun_recalc="data11000mc11000sf0"
#tnpFullRun_recalc=1

tnpFullRun_eval="data00000mc00000sf0"
tnpFullRun_eval="data11000mc11000sf0"
#tnpFullRun_eval=-1

timeStamp="-`date +%Y%m%d-%H%M`"

logDir="./"
logDir="../logs-${anTag}"
logDir="dir-logEtaMax24/"
#logDir="dir-logx/"

#debugMode="DYTools::DEBUG_RUN"
debugMode="DYTools::NORMAL_RUN"
systMode="DYTools::NO_SYST"
systMode="DYTools::UNREGRESSED_ENERGY"

source evaluateESF.sh ${filename} ${debugMode} ${tnpFullRun_eval} ${systMode} | tee ${logDir}/out${timeStamp}-12-evaluateESF-efficiencyScaleFactors${anTag}.log
source recalcESF.sh ${analysisIs2D} ${filename} ${debugMode} ${tnpFullRun_recalc} ${systMode} | tee ${logDir}/out${timeStamp}-12-recalcESF-efficiencyScaleFactors${anTag}.log
