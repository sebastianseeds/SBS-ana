#!/bin/sh

# SSeeds - 10.23.23 - shell script to run all kinematics through full parsing

## Usage
#./run_all_parse.sh <int: e' momentum method> <int: cluster selection method> <int: pass>

nkine=6
kine=(4 7 11 14 8 9)

for ((i=0; i<$nkine; i++))
do

    root -l -b -q 'parse_barebones.C('${kine[i]}',2,4,false,true,true,false)'

wait

done
