#!/bin/bash

x1=$1
x2=$2
x3=$3

if [ ${#x1} -eq 0 ] ; then
    echo "addInQuad.sh x1 x2 [x3]";
    exit
fi

if [ ${#x2} -eq 0 ] ; then x2=0; fi
if [ ${#x3} -eq 0 ] ; then x3=0; fi

root -l -q -b addInQuad.C+\($x1,$x2,$x3\)
