#!/bin/bash

idx=1
if [ ${#1} -gt 0 ] ; then idx=$1; fi

line="\n-----------------------------------------\n"

echo -e ${line}
root -l -q -b printUnfMatrix.C+\(1,${idx}\) | tee tmp1.log
echo -e ${line}
root -l -q -b printUnfMatrix.C+\(2,${idx}\) | tee tmp2.log
echo -e ${line}
diff tmp1.log tmp2.log
echo -e ${line}

