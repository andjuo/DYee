#!/bin/bash

inpFile="../config_files/data_vilnius8TeV_regSSD.conf.py"

runPrepYields=1
runSubtractBkg=1
runPlotSyst=0

runMode=DYTools::NORMAL_RUN

systList="DYTools::ESCALE_DIFF_0000 DYTools::ESCALE_DIFF_0005 DYTools::ESCALE_DIFF_0010 DYTools::ESCALE_DIFF_0015 DYTools::ESCALE_DIFF_0020"

systList="DYTools::ESCALE_DIFF_0005 DYTools::ESCALE_DIFF_0010 DYTools::ESCALE_DIFF_0015 DYTools::ESCALE_DIFF_0020"

#systList="DYTools::ESCALE_DIFF_0000"
timeStamp="-`date +%Y%m%d-%H%M`"

runModeStr=${runMode/DYTools::/}

for systMode in ${systList} ; do

  systModeStr=${systMode/DYTools::/}
  echo ${systModeStr}

  if [ ${runPrepYields} -eq 1 ] ; then
    root -l -q -b prepareYieldsR9.C+\(\"${inpFile}\",${runMode},${systMode}\) \
	| tee log-prepYieldsR9-${runModeStr}-${systModeStr}-${timeStamp}.log
  fi
  if [ ${runSubtractBkg} -eq 1 ] ; then
    root -l -q -b subtractBackgroundR9.C+\(\"${inpFile}\",${runMode},${systMode}\) \
	| tee log-subtractBackgroundR9-${runModeStr}-${systModeStr}-${timeStamp}.log
  fi
done
