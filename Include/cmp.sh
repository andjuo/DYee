#!/bin/bash


dir1=$PWD
dir2=$1

if [ ${#dir2} -eq 0 ] ; then
  echo 
  echo "Error: remote file name to compare is not provided"
  echo 
  echo "cmp.sh: A script to compare two versions of the same file"
  echo "Use:"
  echo "  ./cmp.sh remote_file_including_full_path" 
  echo "This will compare the given file to the file in current directory"
  echo "To change the order of files supplied to diff, add '-r'"
  echo
  exit
fi

fname=${dir2##*/}
#echo fname=${fname}

if [ "$2" == "-r" ] ; then
  diff ${dir1}/${fname} ${dir2}
else
  diff ${dir2} ${dir1}/${fname}
fi

