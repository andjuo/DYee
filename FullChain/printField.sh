#!/bin/bash

fname=$1
fieldType=$2
fieldName=$3

if [ ${#fname} -eq 0 ] ; then
    echo -e "printField.sh  fileName  [fieldType]  fieldName"
    echo -e "  fieldType= M(default) or V"
    exit
fi

if [ ${#fieldType} -gt 1 ] && [ ${#fieldName} -eq 0 ] ; then
    fieldName=${fieldType}
    fieldType=M
fi

root -l -q -b ../FullChain/printField.C+\(\"${fname}\",\"${fieldType}\",\"${fieldName}\"\)
