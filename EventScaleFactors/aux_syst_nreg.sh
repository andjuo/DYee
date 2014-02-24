#!/bin/bash

filename="../config_files/data_vilnius8TeV_regSSD-tagSyst.conf.py"
triggerSet="Full2012_hltEffOld"

tnpFullRun_eval=1
tnpFullRun_eval="data00111mc00111sf0"
tnpFullRun_eval="data00011mc00011sf0"
tnpFullRun_eval="data00011mc00011sf0"

tnpFullRun_recalc="data11111mc11111sf0"
tnpFullRun_recalc="data00000mc00000sf1"
tnpFullRun_recalc=1
tnpFullRun_recalc="data00011mc00011sf0"
#tnpFullRun_recalc=${tnpFullRun_eval}

#tnpFullRun_eval=-1

timeStamp="-`date +%Y%m%d-%H%M`"

logDir="./"
logDir="../logs-${anTag}"
logDir="./logs-unregEn/"
#logDir="dir-logx/"

#debugMode="DYTools::DEBUG_RUN"
debugMode="DYTools::NORMAL_RUN"

no_syst=0

if [ ${no_syst} -eq 1 ] ; then
  systMode="DYTools::NO_SYST"; systName="noSyst";
  source recalcESF.sh ${filename} ${debugMode} ${tnpFullRun_recalc} ${systMode} | tee ${logDir}/out${timeStamp}-12-recalcESF-efficiencyScaleFactors${anTag}-${systName}.log
fi


#systMode="DYTools::NO_SYST"; systName=""
systMode="DYTools::UNREGRESSED_ENERGY"; systName="-nonReg"


source evaluateESF.sh ${filename} ${debugMode} ${tnpFullRun_eval} ${systMode} | tee ${logDir}/out${timeStamp}-12-evaluateESF-efficiencyScaleFactors${anTag}-${systName}.log
source recalcESF.sh ${filename} ${debugMode} ${tnpFullRun_recalc} ${systMode} | tee ${logDir}/out${timeStamp}-12-recalcESF-efficiencyScaleFactors${anTag}-${systName}.log


