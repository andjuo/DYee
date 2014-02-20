#!/bin/bash

anTag="2D"
filename="../config_files/fall8TeV-vilnius.input"
filename="../config_files/data_vilnius8TeV_egamma.conf.py"
triggerSet="Full2012_hltEffOld"

tnpFullRun_eval=1

tnpFullRun_recalc="data11111mc11111sf0"
tnpFullRun_recalc="data00000mc00000sf1"
#tnpFullRun_recalc=1
#tnpFullRun_recalc="data00000mc10000sf0"

tnpFullRun_eval="data11111mc11111sf0"
tnpFullRun_eval=-1
tnpFullRun_eval="data00000mc00000sf1"
tnpFullRun_recalc=-1

timeStamp="-`date +%Y%m%d-%H%M`"

logDir="./"
logDir="../logs-${anTag}"
#logDir="dir-logEtaMax24/"
#logDir="dir-logEtaMax24corr5/"
#logDir="dir-logx/"

#debugMode="DYTools::DEBUG_RUN"
debugMode="DYTools::NORMAL_RUN"

source evaluateESF.sh ${filename} ${debugMode} ${tnpFullRun_eval} | tee ${logDir}/out${timeStamp}-12-evaluateESF-efficiencyScaleFactors${anTag}.log
#source recalcESF.sh ${filename} ${debugMode} ${tnpFullRun_recalc} | tee ${logDir}/out${timeStamp}-12-recalcESF-efficiencyScaleFactors${anTag}.log
