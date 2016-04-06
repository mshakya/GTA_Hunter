#!/bin/bash
# 
# Name: refaa.sh
# 
# Author: Taylor Neely
#
# Function: Renames genes in 
# all .faa files to gi and org
# name in the given directory
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
	#print status
	echo Starting $file
	# read each line of file
	while read line
	do
		if echo "$line" | grep -q \>gi # avoid irrelevant lines
			then
				gi=`echo "$line" | cut -d "|" -f2` # get gi
				name="["`echo "$line" | cut -d "[" -f2` # get name
				# change line
				sed -i "/$gi/c\>$gi $name" $file
		fi
	done < $file
done