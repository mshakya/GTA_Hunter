#!/bin/bash
# 
# Name: xval_custom.sh
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

# directories
ddir="../../data/2016-3-28/"
outd="../../data/results/2016-5-8/prophage/"

#gta
gtas=$ddir"gta/"
#viral
virals=$ddir"viral/"
#test
tests=$ddir"test/prophage/"

### XVALS ###

# Gene 2
gene=2
echo "Starting gene $gene"
# set gta_tool params
k="-k 3"
lam="-p 25"
chem="-y"
c=0.1
d="-d 0.01"
q="-q "$tests$gene"_test.faa"
out=$outd$gene"_res.txt"
> $out
# xval
echo "Gene $gene. Params:c=$c, d=$d, k=$k, lam=$lam, chem=$chem" >> $out
python ../GTA_Hunter.py -g $gtas$gene"_gta.faa" -v $virals$gene"_viral.faa" -w $gtas$gene"_gta.dist" $virals$gene"_viral.dist" $q $d -c $c $k $lam $chem >> $out

# Gene 3
gene=3
echo "Starting gene $gene"
# set gta_tool params
k="-k 3"
lam="-p 25"
chem="-y"
c=1
d="-d 0"
q="-q "$tests$gene"_test.faa"
out=$outd$gene"_res.txt"
> $out
# xval
echo "Gene $gene. Params:c=$c, d=$d, k=$k, lam=$lam, chem=$chem" >> $out
python ../GTA_Hunter.py -g $gtas$gene"_gta.faa" -v $virals$gene"_viral.faa" -w $gtas$gene"_gta.dist" $virals$gene"_viral.dist" $q $d -c $c $k $lam $chem >> $out

# Gene 4
gene=4
echo "Starting gene $gene"
# set gta_tool params
k="-k 2"
lam="-p 5"
chem="-y"
c=1
d="-d 0.05"
q="-q "$tests$gene"_test.faa"
out=$outd$gene"_res.txt"
> $out
# xval
echo "Gene $gene. Params:c=$c, d=$d, k=$k, lam=$lam, chem=$chem" >> $out
python ../GTA_Hunter.py -g $gtas$gene"_gta.faa" -v $virals$gene"_viral.faa" -w $gtas$gene"_gta.dist" $virals$gene"_viral.dist" $q $d -c $c $k $lam $chem >> $out

# Gene 5
gene=5
echo "Starting gene $gene"
# set gta_tool params
k="-k 2"
lam="-p 15"
chem="-y"
c=1
d="-d 0.02"
q="-q "$tests$gene"_test.faa"
out=$outd$gene"_res.txt"
> $out
# xval
echo "Gene $gene. Params:c=$c, d=$d, k=$k, lam=$lam, chem=$chem" >> $out
python ../GTA_Hunter.py -g $gtas$gene"_gta.faa" -v $virals$gene"_viral.faa" -w $gtas$gene"_gta.dist" $virals$gene"_viral.dist" $q $d -c $c $k $lam $chem >> $out

# Gene 6
gene=6
echo "Starting gene $gene"
# set gta_tool params
k="-k 2"
lam="-p 5"
chem="-y"
c=0.01
d="-d 0.01"
q="-q "$tests$gene"_test.faa"
out=$outd$gene"_res.txt"
> $out
# xval
echo "Gene $gene. Params:c=$c, d=$d, k=$k, lam=$lam, chem=$chem" >> $out
python ../GTA_Hunter.py -g $gtas$gene"_gta.faa" -v $virals$gene"_viral.faa" -w $gtas$gene"_gta.dist" $virals$gene"_viral.dist" $q $d -c $c $k $lam $chem >> $out

# Gene 8
gene=8
echo "Starting gene $gene"
# set gta_tool params
k="-k 2"
lam="-p 3"
chem="-y"
c=0.1
d="-d 0.01"
q="-q "$tests$gene"_test.faa"
out=$outd$gene"_res.txt"
> $out
# xval
echo "Gene $gene. Params:c=$c, d=$d, k=$k, lam=$lam, chem=$chem" >> $out
python ../GTA_Hunter.py -g $gtas$gene"_gta.faa" -v $virals$gene"_viral.faa" -w $gtas$gene"_gta.dist" $virals$gene"_viral.dist" $q $d -c $c $k $lam $chem >> $out

# Gene 12
gene=12
echo "Starting gene $gene"
# set gta_tool params
k="-k 4"
lam="-p 3"
chem="-y"
c=0.1
d="-d 0.04"
q="-q "$tests$gene"_test.faa"
out=$outd$gene"_res.txt"
> $out
# xval
echo "Gene $gene. Params:c=$c, d=$d, k=$k, lam=$lam, chem=$chem" >> $out
python ../GTA_Hunter.py -g $gtas$gene"_gta.faa" -v $virals$gene"_viral.faa" -w $gtas$gene"_gta.dist" $virals$gene"_viral.dist" $q $d -c $c $k $lam $chem >> $out

# Gene 13
gene=13
echo "Starting gene $gene"
# set gta_tool params
k="-k 2"
lam="-p 15"
chem="-y"
c=0.1
d="-d 0"
q="-q "$tests$gene"_test.faa"
out=$outd$gene"_res.txt"
> $out
# xval
echo "Gene $gene. Params:c=$c, d=$d, k=$k, lam=$lam, chem=$chem" >> $out
python ../GTA_Hunter.py -g $gtas$gene"_gta.faa" -v $virals$gene"_viral.faa" -w $gtas$gene"_gta.dist" $virals$gene"_viral.dist" $q $d -c $c $k $lam $chem >> $out


# Gene 14
gene=14
echo "Starting gene $gene"
# set gta_tool params
k="-k 3"
lam="-p 15"
chem="-y"
c=0.01
d="-d 0.03"
q="-q "$tests$gene"_test.faa"
out=$outd$gene"_res.txt"
> $out
# xval
echo "Gene $gene. Params:c=$c, d=$d, k=$k, lam=$lam, chem=$chem" >> $out
python ../GTA_Hunter.py -g $gtas$gene"_gta.faa" -v $virals$gene"_viral.faa" -w $gtas$gene"_gta.dist" $virals$gene"_viral.dist" $q $d -c $c $k $lam $chem >> $out

# Gene 15
gene=15
echo "Starting gene $gene"
# set gta_tool params
k="-k 4"
lam="-p 15"
chem="-y"
c=100
d="-d 0.04"
q="-q "$tests$gene"_test.faa"
out=$outd$gene"_res.txt"
> $out
# xval
echo "Gene $gene. Params:c=$c, d=$d, k=$k, lam=$lam, chem=$chem" >> $out
python ../GTA_Hunter.py -g $gtas$gene"_gta.faa" -v $virals$gene"_viral.faa" -w $gtas$gene"_gta.dist" $virals$gene"_viral.dist" $q $d -c $c $k $lam $chem >> $out
