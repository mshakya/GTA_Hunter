#!/bin/bash
# 
# Name: addscore.sh
# 
# Author: Taylor Neely
#
# Function: Adds SVM score and class 
# to the tree file for tree edits
# 
# Input: score file path and change tree path
#
# Notes:
#

# check if arg present
if [ "$1" == "" ] || [ "$2" == "" ]
	then echo "Must specify score output file and tree destination file paths."
else
	scoreFile="$1"
	treeFile="$2"
	out="${treeFile%.*}_scored.tree"
	cat $treeFile > $out
	# read each line of file
	while read line
	do
		if echo "$line" | grep -q \> # avoid irrelevant lines
			then
				liner=$line
				set -- $liner
				liner=$*
				liner=`echo "$liner"`
				# get gi
				gi=`echo "$liner" | cut -d " " -f1 | cut -d ">" -f2`
				# get org
				org=`echo "$liner" | cut -d "[" -f2 | cut -d "]" -f1`
				# get class
				class=`echo "$liner" | rev | cut -d " " -f1 | rev`
				# get score
				score=`echo "$liner" | rev | cut -d " " -f2 | rev`
				sed -i -- "s#$gi#'$gi $org [$class ($score)]'#g" $out
		fi
	done < $scoreFile
fi