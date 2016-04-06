#!/bin/bash
# 
# Name: xval.sh
# 
# Author: Taylor Neely
#
# Function: Batch xval over a variety
# of conditions
# 
# Input: None
#
# Notes:
#

# xval params
folds=10

# select training set
gtas=($(ls ../../data/2016-2-23/gta/*gta.faa))
virals=($(ls ../../data/2016-2-23/viral/*viral.faa))

# select weights
wgtas=($(ls ../../data/2016-2-23/gta/*gta.dist))
wvirals=($(ls ../../data/2016-2-23/viral/*viral.dist))

# confirm sets align
if [ ${#gtas[@]} == ${#virals[@]} ] && [ ${#wgtas[@]} == ${#wvirals[@]} ]
	then
	# test
	for i in "${!gtas[@]}"
	do
		# get files
		gta="${gtas[$i]}"
		wgta="${wgtas[$i]}"
		viral="${virals[$i]}"
		wviral="${wvirals[$i]}"
		# make sure numbers align
		gn=`echo $(basename "$gta") | cut -d "_" -f1`
		wgn=`echo $(basename "$wgta") | cut -d "_" -f1`
		vn=`echo $(basename "$viral") | cut -d "_" -f1`
		wvn=`echo $(basename "$wviral") | cut -d "_" -f1`
		if [ $gn == $vn ] && [ $gn == $wgn ] && [ $gn == $wvn ]
			then
			# define out by number
			echo "Starting gene $gn"
			out="../../data/results/2016-2-25/"$gn"_xval.txt"
			echo > $out
			# set gta_tool params
			c=0.01
			d=0
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.1
			d=0
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=1
			d=0
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=100
			d=0
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=10000
			d=0
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.01
			d=0.01
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.1
			d=0.01
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=1
			d=0.01
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=100
			d=0.01
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=10000
			d=0.01
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.01
			d=0.02
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.1
			d=0.02
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=1
			d=0.02
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=100
			d=0.02
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=10000
			d=0.02
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.01
			d=0.03
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.1
			d=0.03
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=1
			d=0.03
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=100
			d=0.03
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=10000
			d=0.03
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.01
			d=0.04
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.1
			d=0.04
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=1
			d=0.04
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=100
			d=0.04
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=10000
			d=0.04
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.01
			d=0.05
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.1
			d=0.05
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=1
			d=0.05
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=100
			d=0.05
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=10000
			d=0.05
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.01
			d=0.1
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=0.1
			d=0.1
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=1
			d=0.1
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=100
			d=0.1
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			# set gta_tool params
			c=10000
			d=0.1
			k=3
			# xval
			echo "Params:c=$c, d=$d, k=$k, folds=$folds" >> $out
			python ../GTA_Tool.py -g $gta -v $viral -w $wgta $wviral -d $d -c $c -x $folds -m >> $out
			
			echo "We good on gene $gn."
		else
			echo "No good! Gene numbers off (g=$gn, wg=$wgn, v=$vn, wv=$wvn)"
		fi

	done
else
	echo "Batch run cannot work, file numbers do not line up."
fi