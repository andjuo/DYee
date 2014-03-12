#!/bin/bash

#
# Check if the environment variables are set. Assign values if they are empty
#

debugMode=
for arg in $@ ; do
  tmp=${arg/--debug/}
  if [ ${#arg} -ne ${#tmp} ] ; then
      if [ ${arg} == "$1" ] ; then
	  shift # shift argument list
      fi
    debugMode=${arg:7:2}
  fi
done
if [ ${#debugMode} -eq 0 ] ; then debugMode=0; fi
#debugMode=1


confInputFile=$1
fullRun=$2

if [ ${#confInputFile} -eq 0 ] || [ "${confInputFile}" == "default" ] ; then
    confInputFile="../config_files/data_vilnius8TeV_regSSD.conf.py"
fi

# check whether the full run was requested, overriding internal settings
if [ ${#fullRun} -eq 0 ] ; then
  fullRun=1
fi


# 
#  Individual flags to control the calculation
#

doFsrStudy=0
doPUStudy=0
doResolutionStudy=0
doShapeReweight=0
doCalcUnfoldingSyst=0

#
#  Modify flags if fullRun=1
#

if [ ${#fullRun} -eq 6 ] ; then

  echo "5 Flags decoding of '${fullRun}'"
  doFsrStudy=${fullRun:0:1}
  doPUStudy=${fullRun:1:1}
  doResolutionStudy=${fullRun:2:1}
  doShapeReweight=${fullRun:3:1}
  doCalcUnfoldingSyst=${fullRun:5:1}

elif [ ${#fullRun} -eq 1 ] && [ ${fullRun} -eq 1 ] ; then
    doFsrStudy=1; doPUStudy=1; doResolutionStudy=1; doShapeReweight=1;
    doCalcUnfoldingSyst=1
fi

runMode=DYTools::NORMAL_RUN
if [ ${debugMode} -eq 1 ] ; then
  runMode=DYTools::DEBUG_RUN
elif [ ${debugMode} -eq -1 ] ; then
  runMode=DYTools::LOAD_DATA
fi

echo
echo
echo "evaluateUnfoldingSyst.sh: ${runMode}"
echo "    confInputFile=${confInputFile}"
echo "    runFlags=${fullRun}"
echo "       - doFsrStudy=${doFsrStudy}"
echo "       - doPUStudy =${doPUStudy}"
echo "       - doResolutionStudy=${doResolutionStudy}"
echo "       - doShapeReweight=${doShapeReweight}"
echo "       + doCalcUnfoldingSyst=${doCalcUnfoldingSyst}"
echo 
echo



#
#  Flag of an error
#
noError=1

# some other variables
systMode=DYTools::NO_SYST
#FSRreweight=1.0
#FSRmassDiff=1.


# --------------------------------
#    Define functions to run
# --------------------------------

checkFile() { 
  scriptName="evaluateUnfoldingSyst.sh"
  if [ ${noError} -eq 0 ] ; then 
      noError=0
      echo "! ${scriptName}: error present before checking for the file(s) $@"
      return
  fi

  # double underscore not to mess a possible variable with the same name in the script
  for __fname in $@ ; do
# if "-f" does not work (check plain file), 
# one can use -e (name exists), but then directory will return 'true' as well
# 2. else... echo ".. ok" can be removed
  if [ ${#__fname} -gt 0 ] ; then 
      if [ ! -f ${__fname} ] ; then 
	  echo "! ${scriptName}: file ${__fname} is missing"
	  noError=0
      else
	  echo "! ${scriptName}: file <${__fname}> checked ok"
      fi
  fi
  done
}



runPlotDYUnfoldingMatrix() {
  root -b -q -l  plotUnfoldingMatrix.C+\(\"${confInputFile}\",${runMode},${systMode}\)
  #errCode=$?
  #echo "errCode=${errCode}"
  if [ $? != 7 ] ; then noError=0; fi

  echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
  echo 
  if [ ${noError} -eq 0 ] ; then echo "ERROR"; fi
  echo "DONE: plotUnfoldingMatrix.C(\"${confInputFile}\",\"${runMode}\",${systMode}\)"
  echo 
  echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

}

#info:
#int plotUnfoldingMatrix(const TString conf,
#			DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
#			DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST,
#			double FSRreweight=1.0, double FSRmassDiff=1.) {


runCalcUnfoldingSystematics() {
  root -b -q -l  calcUnfoldingSystematics.C+\(\"${confInputFile}\",${debugMode}\)
  if [ $? != 7 ] ; then noError=0;
  else 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: calcUnfoldingSystematics.C+(\"${confInputFile}\",debugMode=${debugMode})"
     echo 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
  fi
}


# --------------------------------
#    Main sequence
# --------------------------------

#
#  Compile header files
#
root -b -q -l rootlogon.C+
if [ $? != 0 ] ; then noError=0; fi

# 
#   Check that the codes compile
#

storeConfInputFile=${confInputFile}
confInputFile="_DebugRun_"
deriveSyst=$(( ${doFsrStudy} + ${doPUStudy} + ${doResolutionStudy} + ${doShapeReweight} ))
# do debug compilation only if both codes have to be run
if [ ${deriveSyst} -gt 0 ] && [ ${doCalcUnfoldingSyst} -gt 0 ] ; then
if [ ${noError} -eq 1 ] && [ ${deriveSyst} -gt 0 ] ; then 
    runPlotDYUnfoldingMatrix;
    checkFile plotUnfoldingMatrix_C.so
fi
if [ ${noError} -eq 1 ] && [ ${doCalcUnfoldingSyst} -gt 0 ] ; then 
    runCalcUnfoldingSystematics; 
    checkFile calcUnfoldingSystematics_C.so
fi
fi
confInputFile=${storeConfInputFile}


#
#   Calculations
#

if [ ${doFsrStudy} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  systMode="DYTools::FSR_STUDY"
  runPlotDYUnfoldingMatrix
fi

if [ ${doPUStudy} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  systMode="DYTools::PU_STUDY"
  runPlotDYUnfoldingMatrix
fi

if [ ${doShapeReweight} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  systMode="DYTools::ESCALE_RESIDUAL"
  runPlotDYUnfoldingMatrix
fi


if [ ${doResolutionStudy} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  systMode="DYTools::RESOLUTION_STUDY"
  runPlotDYUnfoldingMatrix
fi

if [ ${doCalcUnfoldingSyst} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  runCalcUnfoldingSystematics
fi


# return the error code
exit ${noError}
