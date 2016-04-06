#!/bin/bash
# 
# Name: pwdist.sh
# 
# Author: Taylor Neely
#
# Function: Creates pairwise 
# distance files for all .afa 
# files
# 
# Input: Folder path
#
# Notes:
#

# check if arg present
if [ "$1" == "" ]
	then echo "Working in current directory."
else
	# cd into folder path
	echo "Working in " $1 " directory."
	cd $1
fi

# select .afa files
files=($(ls *.afa))

# look at each file
for file in "${files[@]}"
do
	out="${file%.*}.dist"
	# generate pw distances
	raxmlHPC-PTHREADS-SSE3 -s $file -m PROTGAMMAJTT -f x -T 2 -n dist -p 1256
	# rename output
	mv RAxML_distances.dist $out
	# clean out scraps
	rm RAxML_*
done