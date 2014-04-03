#!/bin/bash

binMin=0
binMax=41

binMin=24
binMax=156

bin=${binMin}

while [ ${bin} -lt ${binMax} ] ; do
  echo "bin=${bin}"
  root -l -q -b plotEtEtaPairs.C+\(${bin}\)
  bin=$(( $bin + 1 ))
done
