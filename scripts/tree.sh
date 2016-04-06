#!/bin/bash
# 
# Name: tree.sh
# 
# Author: Taylor Neely
#
# Function: Creates phylogenetic  
# trees from .dist files
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
	out="${file%.*}.tree"
	# generate tree
	FastTree $file > $out
done