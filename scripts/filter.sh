#!/bin/bash
# 
# Name: filter.sh
# 
# Author: Taylor Neely
#
# Function: Batch filtering
# 
# Input: None
#
# Notes:
#

cd ..

number=2
python Loader.py -g ../data/2016-2-16/gta/"$number"_gta.faa -v ../data/2016-2-16/viral/"$number"_viral.faa -q ../data/2016-2-16/test/phage/"$number"_test.faa -o "$number"_test.faa

number=3
python Loader.py -g ../data/2016-2-16/gta/"$number"_gta.faa -v ../data/2016-2-16/viral/"$number"_viral.faa -q ../data/2016-2-16/test/phage/"$number"_test.faa -o "$number"_test.faa

number=4
python Loader.py -g ../data/2016-2-16/gta/"$number"_gta.faa -v ../data/2016-2-16/viral/"$number"_viral.faa -q ../data/2016-2-16/test/phage/"$number"_test.faa -o "$number"_test.faa

number=5
python Loader.py -g ../data/2016-2-16/gta/"$number"_gta.faa -v ../data/2016-2-16/viral/"$number"_viral.faa -q ../data/2016-2-16/test/phage/"$number"_test.faa -o "$number"_test.faa

number=6
python Loader.py -g ../data/2016-2-16/gta/"$number"_gta.faa -v ../data/2016-2-16/viral/"$number"_viral.faa -q ../data/2016-2-16/test/phage/"$number"_test.faa -o "$number"_test.faa

number=7
python Loader.py -g ../data/2016-2-16/gta/"$number"_gta.faa -v ../data/2016-2-16/viral/"$number"_viral.faa -q ../data/2016-2-16/test/phage/"$number"_test.faa -o "$number"_test.faa

number=8
python Loader.py -g ../data/2016-2-16/gta/"$number"_gta.faa -v ../data/2016-2-16/viral/"$number"_viral.faa -q ../data/2016-2-16/test/phage/"$number"_test.faa -o "$number"_test.faa

number=11
python Loader.py -g ../data/2016-2-16/gta/"$number"_gta.faa -v ../data/2016-2-16/viral/"$number"_viral.faa -q ../data/2016-2-16/test/phage/"$number"_test.faa -o "$number"_test.faa

number=12
python Loader.py -g ../data/2016-2-16/gta/"$number"_gta.faa -v ../data/2016-2-16/viral/"$number"_viral.faa -q ../data/2016-2-16/test/phage/"$number"_test.faa -o "$number"_test.faa

number=13
python Loader.py -g ../data/2016-2-16/gta/"$number"_gta.faa -v ../data/2016-2-16/viral/"$number"_viral.faa -q ../data/2016-2-16/test/phage/"$number"_test.faa -o "$number"_test.faa

number=14
python Loader.py -g ../data/2016-2-16/gta/"$number"_gta.faa -v ../data/2016-2-16/viral/"$number"_viral.faa -q ../data/2016-2-16/test/phage/"$number"_test.faa -o "$number"_test.faa

number=15
python Loader.py -g ../data/2016-2-16/gta/"$number"_gta.faa -v ../data/2016-2-16/viral/"$number"_viral.faa -q ../data/2016-2-16/test/phage/"$number"_test.faa -o "$number"_test.faa
