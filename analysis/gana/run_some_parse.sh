#!/bin/sh

# SSeeds - 10.23.23 - shell script to run all kinematics through full parsing

## Usage
#./run_all_parse.sh <int: e' momentum method> <int: cluster selection method> <int: pass>

epm=$1
cluster_method=$2
pass=$3
verbose=$4

nkine=3
kine=(4 8 9)

for ((i=0; i<$nkine; i++))
do

    root -l -b -q 'parse.C('${kine[i]}','$epm','$cluster_method','$pass','$verbose')'

wait

done
