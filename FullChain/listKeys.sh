#!/bin/bash

fnames=$@

if [ ${#fnames} -eq 0 ] ; then
    echo -e "printField.sh  fileName1 [fileName2 ...]"
    exit
fi

for f in ${fnames} ; do
  root -l -q -b ../FullChain/printField.C+\(\"${f}\"\)
done
