#!/bin/bash

#
# use for debug, tuning and recalculation of the efficiencies or scale factors
# evaluateESF.sh is the main script
#

debugMode="DYTools::NORMAL_RUN"
fullRun=1

if [ ${#1} -gt 0 ] ; then confInputFile=$1; fi
if [ ${#2} -gt 0 ] ; then debugMode=$2; fi

if [ ${#3} -gt 0 ] ; then 
    fullRun=$3;
    echo "fullRun=<${fullRun}>"
    if [ ${#fullRun} -eq 2 ] && [ ${fullRun} -eq -1 ] ; then
	echo -e "\n\t fullRun=-1. Skipping recalcESF.sh\n\n"
	exit
    fi
fi

collectEvents=0  # it is recommended to have collectEvents=1 in evaluateESF!

# if you do not want to have the time stamp, comment the line away 
# or set timeStamp=
timeStamp="-`date +%Y%m%d-%H%M`"
#timeStamp=

#
# Check if the environment variables are set. Assign values if they are empty
#
if [ ${#confInputFile} -eq 0 ] ; then
    confInputFile="../config_files/data8TeV.input" # used in CalcEventEff.C
fi
if [ ${#debugMode} -eq 0 ] ; then
    debugMode="DYTools::NORMAL_RUN"
fi

echo
echo
echo "recalcESF.sh:"
echo "    confInputFile=${confInputFile}"
echo "    timeStamp=${timeStamp}"
echo "    debugMode=${debugMode}"
echo 
echo

# 
#  Individual flags to control the calculation
#

runMC_Reco=0
runMC_Id=0
runMC_Hlt=0
runMC_Hlt_leg1=0
runMC_Hlt_leg2=0
runData_Reco=0
runData_Id=0
runData_Hlt=0
runData_Hlt_leg1=0
runData_Hlt_leg2=0
runCalcEventEff=1

#
#  Modify flags if fullRun=1
#

longFullRunFlag=0
if [ ${#fullRun} -eq 7 ] ; then
    longFullRunFlag=1
  if [ ${fullRun:0:1} -eq 1 ] ; then runMC_Reco=1; else runMC_Reco=0; fi
  if [ ${fullRun:1:1} -eq 1 ] ; then runMC_Id=1; else runMC_Id=0; fi
  if [ ${fullRun:2:1} -eq 1 ] ; then runMC_Hlt=1; else runMC_Hlt=0; fi
  if [ ${fullRun:3:1} -eq 1 ] ; then runData_Reco=1; else runData_Reco=0; fi
  if [ ${fullRun:4:1} -eq 1 ] ; then runData_Id=1; else runData_Id=0; fi
  if [ ${fullRun:5:1} -eq 1 ] ; then runData_Hlt=1; else runData_Hlt=0; fi
  if [ ${fullRun:6:1} -eq 1 ] ; then runCalcEventEff=1; else runCalcEventEff=0; fi

elif [ ${#fullRun} -eq 11 ] ; then
    longFullRunFlag=1
  if [ ${fullRun:0:1} -eq 1 ] ; then runData_Reco=1; else runData_Reco=0; fi
  if [ ${fullRun:1:1} -eq 1 ] ; then runData_Id=1; else runData_Id=0; fi
  if [ ${fullRun:2:1} -eq 1 ] ; then runData_Hlt=1; else runData_Hlt=0; fi
  if [ ${fullRun:3:1} -eq 1 ] ; then runData_Hlt_leg1=1; else runData_Hlt_leg1=0; fi
  if [ ${fullRun:4:1} -eq 1 ] ; then runData_Hlt_leg2=1; else runData_Hlt_leg2=0; fi
  if [ ${fullRun:5:1} -eq 1 ] ; then runMC_Reco=1; else runMC_Reco=0; fi
  if [ ${fullRun:6:1} -eq 1 ] ; then runMC_Id=1; else runMC_Id=0; fi
  if [ ${fullRun:7:1} -eq 1 ] ; then runMC_Hlt=1; else runMC_Hlt=0; fi
  if [ ${fullRun:8:1} -eq 1 ] ; then runMC_Hlt_leg1=1; else runMC_Hlt_leg1=0; fi
  if [ ${fullRun:9:1} -eq 1 ] ; then runMC_Hlt_leg2=1; else runMC_Hlt_leg2=0; fi
  if [ ${fullRun:10:1} -eq 1 ] ; then runCalcEventEff=1; else runCalcEventEff=0; fi

elif [ ${#fullRun} -eq 19 ] ; then  # "dataXXXXXmcXXXXXsfX"
    longFullRunFlag=1
  if [ ! ${fullRun:0:4} == "data" ] || [ ! ${fullRun:9:2} == "mc" ] || [ ! ${fullRun:16:2} == "sf" ] ; then
      echo "recalcESF.sh the expected format for a string of length 19 is 'dataXXXXXmcXXXXXsfX'"
      exit 0
  fi
  if [ ${fullRun:4:1} -eq 1 ] ; then runData_Reco=1; else runData_Reco=0; fi
  if [ ${fullRun:5:1} -eq 1 ] ; then runData_Id=1; else runData_Id=0; fi
  if [ ${fullRun:6:1} -eq 1 ] ; then runData_Hlt=1; else runData_Hlt=0; fi
  if [ ${fullRun:7:1} -eq 1 ] ; then runData_Hlt_leg1=1; else runData_Hlt_leg1=0; fi
  if [ ${fullRun:8:1} -eq 1 ] ; then runData_Hlt_leg2=1; else runData_Hlt_leg2=0; fi
  if [ ${fullRun:11:1} -eq 1 ] ; then runMC_Reco=1; else runMC_Reco=0; fi
  if [ ${fullRun:12:1} -eq 1 ] ; then runMC_Id=1; else runMC_Id=0; fi
  if [ ${fullRun:13:1} -eq 1 ] ; then runMC_Hlt=1; else runMC_Hlt=0; fi
  if [ ${fullRun:14:1} -eq 1 ] ; then runMC_Hlt_leg1=1; else runMC_Hlt_leg1=0; fi
  if [ ${fullRun:15:1} -eq 1 ] ; then runMC_Hlt_leg2=1; else runMC_Hlt_leg2=0; fi
  if [ ${fullRun:18:1} -eq 1 ] ; then runCalcEventEff=1; else runCalcEventEff=0; fi

elif [ ${#fullRun} -eq 1 ] && [ ${fullRun} -eq 1 ] ; then
  runMC_Reco=1; runMC_Id=1; runMC_Hlt=1
  runMC_Hlt_leg1=1;
  runMC_Hlt_leg2=1;
  runData_Reco=1; runData_Id=1; runData_Hlt=1
  runData_Hlt_leg1=1;
  runData_Hlt_leg2=1;
  runCalcEventEff=1   # it prepares the skim for event efficiencies

else
  echo
  echo "recalcESF.sh failed to recognize the string fullRun[${#fullRun}]=\"${fullRun}\""
  echo
  exit 0
fi

if [ ${longFullRunFlag} -eq 1 ] ; then
  echo 
  echo " recalcESF.sh special fullRun string obtained (${fullRun}), decoded as:"
  echo "  runData_Reco=${runData_Reco}, runData_Id=${runData_Id}, runData_Hlt=${runData_Hlt}, runData_Hlt_leg1=${runData_Hlt_leg1}, runData_Hlt_leg2=${runData_Hlt_leg2}"
  echo "  runMC_Reco=${runMC_Reco}, runMC_Id=${runMC_Id}, runMC_Hlt=${runMC_Hlt} runMC_Hlt_leg1=${runMC_Hlt_leg1}, runMC_Hlt_leg2=${runMC_Hlt_leg2}"
  echo "  runCalcEventEff=${runCalcEventEff}"
  echo
fi


#
#  Flag of an error
#
noError=1


# --------------------------------
#    Define functions to run
# --------------------------------

checkFile() { 
  if [ ${noError} -eq 0 ] ; then 
      noError=0
      echo "! recalcESF.sh: error present before checking for the file(s) $@"
      return
  fi

  # double underscore not to mess a possible variable with the same name in the script
  for __fname in $@ ; do
# if "-f" does not work (check plain file), 
# one can use -e (name exists), but then directory will return 'true' as well
# 2. else... echo ".. ok" can be removed
  if [ ${#__fname} -gt 0 ] ; then 
      if [ ! -f ${__fname} ] ; then 
	  echo "! recalcESF.sh: file ${__fname} is missing"
	  noError=0
      else
	  echo "! recalcESF.sh: file <${__fname}> checked ok"
      fi
  fi
  done
}

runCalcEff() {
 effKind=$1
# calculate
 root -l -q -b  ${LXPLUS_CORRECTION} calcEff.C+\(\"${inpFile}\",\"${effKind}\",${onData}\) \
     | tee log${timeStamp}-calcEff-${dataKind}-${effKind}.out
  if [ $? != 0 ] ; then noError=0;
  else 
     checkFile calcEff_C.so
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: calcEff(\"$inpFile\",\"${effKind}\")"
     echo 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
  fi
}

runCalcEventEff() {
 _collectEvents=$1
 echo "_collectEvents=${_collectEvents}"
 if [ ${#_collectEvents} -eq 0 ] ; then _collectEvents=0; fi
 root -l -q -b  ${LXPLUS_CORRECTION} calcEventEff.C+\(\"${inpFile}\",${_collectEvents},${debugMode}\) \
     | tee log${timeStamp}-calcEventEff.out
  if [ $? != 0 ] ; then noError=0;
  else 
      checkFile calcEventEff_C.so
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: calcEventEff(\"${inpFile}\",collectEvents=${_collectEvents},debug=${debugMode})"
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
checkFile tnpSelectEvents_hh.so

# 
#   Check that the codes compile
#

inpFile="_check_"
onData=0
dataKind="check"
runAny=$(( ${runMC_Reco} + ${runMC_Id} + ${runMC_Hlt} + ${runData_Reco} + ${runData_Id} + ${runData_Hlt} ))
if [ ${runAny} -eq 0 ] ; then
    runAny=$(( ${runMC_Hlt_leg1} + ${runMC_Hlt_leg2} + ${runData_Hlt_leg1} + ${runData_Hlt_leg2} ))
fi
if [ ${runAny} -gt 0 ] && [ ${noError} -eq 1 ] ; then runCalcEff; fi
if [ ${runCalcEventEff} -eq 1 ] && [ ${noError} -eq 1 ] ; then runCalcEventEff; fi
if [ ${noError} -eq 1 ] ; then echo; echo "  -=- Resuming normal calculation -=-"; echo; fi

inpFile=${confInputFile}


# Process MC

onData=0
dataKind="mc"

if [ ${runMC_Reco} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  runCalcEff "RECO"
fi

if [ ${runMC_Id} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  runCalcEff "ID"
fi

if [ ${runMC_Hlt} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  runCalcEff "HLT"
fi

if [ ${runMC_Hlt_leg1} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  runCalcEff "HLTleg1"
fi

if [ ${runMC_Hlt_leg2} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  runCalcEff "HLTleg2"
fi


# Process data

onData=1
dataKind="data"

if [ ${runData_Reco} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    runCalcEff "RECO"
fi

if [ ${runData_Id} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    runCalcEff "ID"
fi

if [ ${runData_Hlt} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    runCalcEff "HLT"
fi

if [ ${runData_Hlt_leg1} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    runCalcEff "HLTleg1"
fi

if [ ${runData_Hlt_leg2} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    runCalcEff "HLTleg2"
fi


# Calculate efficiency scale factors

if [ ${runCalcEventEff} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    runCalcEventEff ${collectEvents}
fi


# return the error code
exit ${noError}

