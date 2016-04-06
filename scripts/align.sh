#!/bin/bash
# 
# Name: align.sh
# 
# Author: Taylor Neely
#
# Function: Aligns all .faa files
# in the given directory
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

# select .faa files
files=($(ls *.faa))

# look at each file
for file in "${files[@]}"
do
	out="${file%.*}.afa"
	muscle -in $file -out $out
done