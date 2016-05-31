#!/bin/bash
# 
# Name: test.sh
# 
# Author: Taylor Neely
#
# Function: Batch testing
# 
# Input: None
#
# Notes:
#

# set gta_tool params
c=0.1
d=0.03
k=4

# select training set
gtas=($(ls ../../data/2016-2-23/gta/*gta.faa))
virals=($(ls ../../data/2016-2-23/viral/*viral.faa))

# select weights
wgtas=($(ls ../../data/2016-2-23/gta/*gta.dist))
wvirals=($(ls ../../data/2016-2-23/viral/*viral.dist))

# select test set
testers=($(ls ../../data/2016-2-23/test/phage/*test.faa))

# confirm sets align
if [ ${#gtas[@]} == ${#virals[@]} ] && [ ${#gtas[@]} == ${#testers[@]} ]
	then
	# test
	for i in "${!gtas[@]}"
	do
		# get files
		gta="${gtas[$i]}"
		wgta="${wgtas[$i]}"
		viral="${virals[$i]}"
		wviral="${wvirals[$i]}"
		tester="${testers[$i]}"
		out="../../data/results/2016-2-26/"$(basename "${tester%.*}.res")
		# make sure numbers align
		gn=`echo $(basename "$gta") | cut -d "_" -f1`
		wgn=`echo $(basename "$wgta") | cut -d "_" -f1`
		vn=`echo $(basename "$viral") | cut -d "_" -f1`
		wvn=`echo $(basename "$wviral") | cut -d "_" -f1`
		tn=`echo $(basename "$tester") | cut -d "_" -f1`
		if [ $gn == $vn ] && [ $gn == $tn ] && [ $gn == $wgn ] && [ $gn == $wvn ]
			then
			# run gta_tool
			python ../GTA_Hunter.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -q $tester > $out
			echo "We good on gene $gn."
		else
			echo "No good! Gene numbers off (g=$gn, wg=$wgn, v=$vn, wv=$wvn, t=$tn)"
		fi

	done
else
	echo "Batch run cannot work, file numbers do not line up."
fi