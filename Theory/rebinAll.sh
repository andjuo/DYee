#!/bin/bash

files=`ls 2D*.root`
#files=`ls 2D*200to1500*.root`

processed=

for f in ${files} ; do
  ftest=${f/-mdf.root/}
  if [ ${#ftest} -eq ${#f} ] ; then
      echo "processing file $f"
      root -l -q -b doRebin.C+\(\"$f\"\)
      processed="$processed $f"
  else
      echo "skipping $f"
  fi
done

echo
echo "processed:"
echo

for p in $processed ; do
  echo "  $p"
done
