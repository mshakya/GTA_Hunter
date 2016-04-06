#!/bin/bash
# 
# Name: process.sh
# 
# Author: Taylor Neely
#
# Function: Does all formating and file
# generation needed (except add score)
# 
# Input: score file path and change tree path
#
# Notes:
#

# check if arg present
if [ "$1" == "" ]
	then 
	echo "Working in current directory."
	path="./"
else
	# specify path
	path=$1
fi

# rename files
./refaa.sh $path
# align files
./align.sh $path
# create pwd files
./pwdist.sh $path
# make trees
./tree.sh $path