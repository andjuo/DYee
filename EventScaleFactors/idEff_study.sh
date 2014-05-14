#!/bin/bash

params="dxyz invEminusInvP Aeff relPFIso dEta dPhi sigmaIEtaIEta HoverE"
count=0
for p in ${params} ; do count=$(( $count+1 )); done
echo -e "\tThere are $count parameters"

applyIMin=$1
applyIMax=$2

if [ ${#applyIMin} -eq 0 ] || [ ${#applyIMax} -eq 0 ] ; then
    echo "please provide indices to apply"
    exit
fi

factors="0.9 1.1"
isDataStr="0 1"
logDir="./dir-idSyst/"
timeStamp="-`date +%Y%m%d-%H%M`"

inpFile="default"
inpFile="cernRemote"

count=0
for p in ${params} ; do
  if [ ${count} -ge ${applyIMin} ] && [ ${count} -le ${applyIMax} ] ; then
     echo "working with param=${p}"
     for f in ${factors} ; do
       effIdString="ID_idSyst${p}_${f}"
       effIdStringNoDot=${effIdString/"."/}
       for isData in ${isDataStr} ; do
	 echo "work on \"${inpFile}\", effIdString=${effIdString}, effIdStringNoDot=${effIdStringNoDot}, isData=${isData}"
	 root -l -q -b eff_IdHlt.C+\(0,\"${inpFile}\",\"${effIdString}\",${isData},DYTools::NORMAL_RUN,DYTools::UNREGRESSED_ENERGY\) | tee ${logDir}log-effID-${timeStamp}-${isData}-${effIdStringNoDot}-sel.log
	 root -l -q -b calcEff.C+\(0,\"${inpFile}\",\"${effIdStringNoDot}\",${isData},0,DYTools::UNREGRESSED_ENERGY\) | tee ${logDir}log-effID-${timeStamp}-${isData}-${effIdStringNoDot}-calc.log
       done
     done
  fi
  count=$(( $count + 1 ))
done

