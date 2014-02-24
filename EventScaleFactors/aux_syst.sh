#!/bin/bash

filename="../config_files/data_vilnius8TeV_regSSD-tagSyst.conf.py"
triggerSet="Full2012_hltEffOld"

tnpFullRun_eval=1
tnpFullRun_eval="data00111mc00111sf0"
tnpFullRun_eval="data00011mc00011sf0"
tnpFullRun_eval="data00000mc00001sf0"

tnpFullRun_recalc="data11111mc11111sf0"
tnpFullRun_recalc="data00000mc00000sf1"
tnpFullRun_recalc=1
tnpFullRun_recalc="data00000mc00011sf0"
#tnpFullRun_recalc=${tnpFullRun_eval}

#tnpFullRun_eval=-1

timeStamp="-`date +%Y%m%d-%H%M`"

logDir="./"
logDir="../logs-${anTag}"
logDir="./logs-tagSyst/"
#logDir="dir-logx/"

#debugMode="DYTools::DEBUG_RUN"
debugMode="DYTools::NORMAL_RUN"

no_syst=0

if [ ${no_syst} -eq 1 ] ; then
  systMode="DYTools::NO_SYST"; systName="noSyst";
  source recalcESF.sh ${filename} ${debugMode} ${tnpFullRun_recalc} ${systMode} | tee ${logDir}/out${timeStamp}-12-recalcESF-efficiencyScaleFactors${anTag}-${systName}.log
fi


systMode="DYTools::NO_SYST"; systName=""
systMode="DYTools::RESOLUTION_STUDY"; systName="-tagLowerPt"
systMode="DYTools::PILEUP_5plus"; systName="-pu5plus"
systMode="DYTools::PILEUP_5minus"; systName="-pu5minus"
#systMode="DYTools::FSR_STUDY"; systName="-tagMediumID"

source evaluateESF.sh ${filename} ${debugMode} ${tnpFullRun_eval} ${systMode} | tee ${logDir}/out${timeStamp}-12-evaluateESF-efficiencyScaleFactors${anTag}-${systName}.log
source recalcESF.sh ${filename} ${debugMode} ${tnpFullRun_recalc} ${systMode} | tee ${logDir}/out${timeStamp}-12-recalcESF-efficiencyScaleFactors${anTag}-${systName}.log


exit

tnpFullRun_eval="data00000mc00011sf0"

#systMode="DYTools::RESOLUTION_STUDY"; systName="-tagLowerPt"
systMode="DYTools::FSR_STUDY"; systName="-tagMediumID"
systMode="DYTools::PILEUP_5plus"; systName="-pu5plus"

source evaluateESF.sh ${filename} ${debugMode} ${tnpFullRun_eval} ${systMode} | tee ${logDir}/out${timeStamp}-12-evaluateESF-efficiencyScaleFactors${anTag}-${systName}.log
source recalcESF.sh ${filename} ${debugMode} ${tnpFullRun_recalc} ${systMode} | tee ${logDir}/out${timeStamp}-12-recalcESF-efficiencyScaleFactors${anTag}-${systName}.log
