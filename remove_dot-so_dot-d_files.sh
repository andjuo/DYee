#!/bin/bash +x

if [ ${#2} -eq 0 ] ; then
  echo
  echo The script removes *.so *.d *~ files in ALL subdirectories
  echo remove_dot-so_dot-d_tilda_files.sh [dirName [\$PWD]]
  echo
fi

start="."
parentDir="${PWD}"
if [ ${#2} -gt 0 ] ; then parentDir=$2; fi
if [ ${#1} -gt 0 ] ; then start=$1; fi
#echo
#echo "start=<$start>, parentDir=<$parentDir>"
echo "currDir=$PWD"
#echo

cd ${start}
dirCands=`ls .`

#echo $dirCands

for d in ${dirCands} ; do
  if [ -d ${d} ] ; then
      #echo "directory=$d"
      cd $d
      source ${parentDir}/remove_dot-so_dot-d_tilda_files.sh . ${parentDir}
      cd ..
  fi
done

#ls *.so *.d *~ *.dll *.stackdump
#rm -f *.so *.d *~ *.dll *.stackdump
rm -f *.so *.d *.dll *.stackdump

